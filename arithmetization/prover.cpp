#include "poly.hpp"
#include "matrixFr.hpp"
#include <unordered_set>
#include <utility>
#include <tuple>
#include <unordered_map>
#include <algorithm>

using namespace std::chrono;
typedef vector<vector<vector<uint32_t>>> vvv_int;

void inttAndCreatePolyff(std::vector<Fr> &vec, polyff &result)
{
    Ntt<Fr> ntt;
    ntt.init(nextPowOf2(vec.size()));
    ntt.intt(&vec[0]);
    result.data = std::move(vec);
}

std::tuple<vector<polyff>, vector<vector<Fr>>, matrixFr>
permutation(const vector<Fr> &W, size_t row, size_t col, const vector<vector<vector<uint32_t>>> &positions)
{
    matrixFr index(row, col);
    uint64_t N = nextPowOf2(W.size());
    Ntt<Fr> ntt;
    ntt.init(N);

    vector<polyff> S(col);
    vector<vector<Fr>> S_ij(row, vector<Fr>(col));

    for (size_t i = 0; i < col; i++)
    {

        vector<Fr> scaledW = W;
        Fr scale = i + 1;
        mulveccons(scaledW, scale);
        index.replaceColumn(i, scaledW);
    }

    matrixFr index_origin = index;
    index.swapMultipleElements(positions);

    for (size_t i = 0; i < col; i++)
    {

        vector<Fr> column = index.getColumn(i);
        if (column.size() < N)
        {
            column.resize(N, Fr(0));
        }
        ntt.intt(&column[0]);
        S[i] = polyff(column.begin(), column.end());
    }

    S_ij = index.data;

    return {S, S_ij, index_origin};
}

polyff copyconstraints(const matrixFr &A, const vector<polyff> &S, const vector<vector<Fr>> &S_ij, const matrixFr &index_origin, Fr alpha, Fr beta, const vvv_int &positions)
{
    auto [row, col] = A.getSize();

    vector<Fr> accs;
    accs.push_back(Fr(1));

    matrixFr modifiedA = A;
    modifiedA.duplicateMultipleElements(positions);

    Fr acc = 1;
    for (size_t i = 0; i < row; ++i)
    {
        for (size_t j = 0; j < col; ++j)
        {

            acc *= (modifiedA[i][j] + alpha * index_origin[i][j] + beta) / (modifiedA[i][j] + alpha * S_ij[i][j] + beta);
        }

        if (i < (row - 1))
        {
            accs.push_back(acc);
        }
    }

    vector<Fr> accs_copy = accs;
    Ntt<Fr> ntt;
    ntt.init(nextPowOf2(accs_copy.size()));

    ntt.intt(&accs_copy[0]);

    polyff acc_x;
    acc_x.data = accs_copy;

    return acc_x;
}

string vecToStr(const vector<Fr> &vec)
{
    string str;
    for (const auto &v : vec)
    {
        str += v.getStr(10) + ",";
    }
    return str;
}

tuple<vector<Fr>, vector<Fr>, vector<Fr>>
lookup(const matrixFr &trace, const matrixFr &table, uint32_t selector, uint32_t numberofselector, uint32_t numberofwire, Fr &alpha)
{
    vector<vector<Fr>> f;
    vector<vector<Fr>> t;

    unordered_map<string, vector<Fr>> tableMap;

    for (size_t j = 0; j < table.getSize().first; ++j)
    {
        string key = vecToStr(vector<Fr>(table[j].begin(), table[j].begin() + numberofwire));
        tableMap[key] = table[j];
    }

    for (size_t i = 0; i < trace.getSize().first; ++i)
    {
        if (trace[i][selector].isOne())
        {
            vector<Fr> traceSegment(trace[i].begin() + numberofselector, trace[i].begin() + numberofselector + numberofwire);
            string key = vecToStr(traceSegment);
            auto it = tableMap.find(key);
            if (it != tableMap.end())
            {
                f.push_back(it->second);
            }
        }
        else
        {
            f.push_back(table[table.getSize().first - 1]);
        }
    }

    matrixFr matf;
    matf.data = f;

    vector<Fr> alphaPowers(matf.getSize().second);
    Fr alphaPow = 1;
    for (uint32_t i = 0; i < matf.getSize().second; i++)
    {
        alphaPowers[i] = alphaPow;
        alphaPow *= alpha;
    }

    vector<Fr> f_result(trace.getSize().first, Fr(0));
    vector<Fr> t_result(trace.getSize().first, Fr(0));

    for (size_t col = 0; col < matf.getSize().second; ++col)
    {

        vector<Fr> f_col = matf.getColumn(col);
        vector<Fr> t_col = table.getColumn(col);

        mulveccons(f_col, alphaPowers[col]);
        mulveccons(t_col, alphaPowers[col]);

        if (col == 0)
        {
            f_result = f_col;
            t_result = t_col;
        }
        else
        {
            f_result = addVectors(f_result, f_col);
            t_result = addVectors(t_result, t_col);
        }
    }

    return make_tuple(f_result, t_result, matf.getColumn(5));
}

tuple<polyff, polyff, polyff, polyff, polyff, polyff, polyff, polyff, polyff>
plookup_copyconstraints(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma)
{
    auto start = high_resolution_clock::now();
    auto checkpoint = start;

    size_t m = f.size();
    size_t n = t.size();

    size_t N = Max(m, n);
    N = nextPowOf2(N);

    if (m < N)
    {
        Fr lastElementOfT = t.back();
        f.resize(N, lastElementOfT);
    }

    if (n < N)
    {
        Fr lastElementOfT = t.back();
        t.resize(N, lastElementOfT);
    }

    checkpoint = high_resolution_clock::now();
    cout << "Resize time: " << duration_cast<microseconds>(checkpoint - start).count() << " microseconds\n";

    unordered_map<Fr, size_t> frequencyMap;
    for (auto &element : f)
    {
        frequencyMap[element]++;
    }

    auto map_end1 = high_resolution_clock::now();
    cout << "Frequency map creation time1: " << duration_cast<microseconds>(map_end1 - checkpoint).count() << " microseconds\n";
    checkpoint = map_end1;

    vector<Fr> s;
    s.resize(f.size() + t.size());

    size_t s_index = 0;
    for (auto &element : t)
    {
        s[s_index++] = element;
        if (frequencyMap.find(element) != frequencyMap.end())
        {
            size_t count = frequencyMap[element];
            frequencyMap[element] = 0;
            for (size_t i = 0; i < count; ++i)
            {
                s[s_index++] = element;
            }
        }
    }

    auto map_end2 = high_resolution_clock::now();
    cout << "Frequency map creation time2: " << duration_cast<microseconds>(map_end2 - checkpoint).count() << " microseconds\n";
    checkpoint = map_end2;

    t.push_back(t.back());
    s.push_back(t.back());

    vector<Fr> accs;
    Fr acc = 1;
    accs.push_back(acc);

    for (size_t i = 1; i < N + 1; ++i)
    {
        acc *= ((1 + beta) * (gamma + f[i - 1]) * (gamma * (1 + beta) + t[i - 1] + beta * t[i])) /
               ((gamma * (1 + beta) + s[i - 1] + beta * s[i]) * (gamma * (1 + beta) + s[i + N - 1] + beta * s[i + N]));

        if (i < (N))
        {
            accs.push_back(acc);
        }
    }

    auto acc_end = high_resolution_clock::now();
    cout << "Accumulator calculation time: " << duration_cast<microseconds>(acc_end - checkpoint).count() << " microseconds\n";
    checkpoint = acc_end;

    vector<Fr> accs_copy = accs;
    Ntt<Fr> ntt;
    ntt.init(nextPowOf2(accs_copy.size()));

    ntt.intt(&accs_copy[0]);

    polyff zx;
    zx.data = accs_copy;
    vector<Fr> zwx_copy(accs.begin() + 1, accs.end());
    zwx_copy.push_back(Fr(1));
    vector<Fr> fx_copy(f.begin(), f.end());
    vector<Fr> tx_copy(t.begin(), t.end() - 1);
    vector<Fr> twx_copy(t.begin() + 1, t.end());
    vector<Fr> h1x_copy(s.begin(), s.begin() + N);
    vector<Fr> h1wx_copy(s.begin() + 1, s.begin() + N + 1);
    vector<Fr> h2x_copy(s.begin() + N, s.begin() + 2 * N);
    vector<Fr> h2wx_copy(s.begin() + N + 1, s.begin() + (2 * N + 1));

    polyff zwx, fx, tx, twx, h1x, h1wx, h2x, h2wx;

    if (zwx_copy.size() > std::pow(2, 15))
    {

        std::vector<std::thread> threads;

        threads.emplace_back(inttAndCreatePolyff, std::ref(zwx_copy), std::ref(zwx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(fx_copy), std::ref(fx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(tx_copy), std::ref(tx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(twx_copy), std::ref(twx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h1x_copy), std::ref(h1x));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h1wx_copy), std::ref(h1wx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h2x_copy), std::ref(h2x));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h2wx_copy), std::ref(h2wx));

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {

        inttAndCreatePolyff(zwx_copy, zwx);
        inttAndCreatePolyff(fx_copy, fx);
        inttAndCreatePolyff(tx_copy, tx);
        inttAndCreatePolyff(twx_copy, twx);
        inttAndCreatePolyff(h1x_copy, h1x);
        inttAndCreatePolyff(h1wx_copy, h1wx);
        inttAndCreatePolyff(h2x_copy, h2x);
        inttAndCreatePolyff(h2wx_copy, h2wx);
    }

    auto end = high_resolution_clock::now();
    cout << "Total time: " << duration_cast<microseconds>(end - start).count() << " microseconds\n";

    return make_tuple(zx, zwx, fx, tx, twx, h1x, h1wx, h2x, h2wx);
}

tuple<polyff, polyff, polyff>
plookup_poly_ft(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma)
{
    size_t m = f.size();
    size_t n = t.size();

    size_t N = max(m, n);
    N = nextPowOf2(N);

    if (m < N)
    {
        Fr lastElementOfF = f.back();
        f.resize(N, lastElementOfF);
    }

    if (n < N)
    {
        Fr lastElementOfT = t.back();
        t.resize(N, lastElementOfT);
    }

    vector<Fr> fx_copy(f.begin(), f.end());
    vector<Fr> tx_copy(t.begin(), t.end());
    vector<Fr> twx_copy(t.begin() + 1, t.end());
    twx_copy.push_back(t.back());

    Ntt<Fr> ntt;
    ntt.init(N);

    ntt.intt(&fx_copy[0]);
    ntt.intt(&tx_copy[0]);
    ntt.intt(&twx_copy[0]);

    polyff fx;
    fx.data = fx_copy;
    polyff tx;
    tx.data = tx_copy;
    polyff twx;
    twx.data = twx_copy;

    return make_tuple(fx, tx, twx);
}

tuple<polyff, polyff, polyff, polyff, polyff, polyff, polyff>
plookup_copyconstraints_1(vector<Fr> f, vector<Fr> t)
{
    auto start = high_resolution_clock::now();
    auto checkpoint = start;

    size_t m = f.size();
    size_t n = t.size();

    size_t N = Max(m, n);
    N = nextPowOf2(N);

    if (m < N)
    {
        Fr lastElementOfT = t.back();
        f.resize(N, lastElementOfT);
    }

    if (n < N)
    {
        Fr lastElementOfT = t.back();
        t.resize(N, lastElementOfT);
    }

    checkpoint = high_resolution_clock::now();

    unordered_map<Fr, size_t> frequencyMap;
    for (auto &element : f)
    {
        frequencyMap[element]++;
    }

    auto map_end1 = high_resolution_clock::now();

    checkpoint = map_end1;

    vector<Fr> s;
    s.resize(f.size() + t.size());

    size_t s_index = 0;
    for (auto &element : t)
    {
        s[s_index++] = element;
        if (frequencyMap.find(element) != frequencyMap.end())
        {
            size_t count = frequencyMap[element];
            frequencyMap[element] = 0;
            for (size_t i = 0; i < count; ++i)
            {
                s[s_index++] = element;
            }
        }
    }

    auto map_end2 = high_resolution_clock::now();

    checkpoint = map_end2;

    t.push_back(t.back());
    s.push_back(t.back());

    auto acc_end = high_resolution_clock::now();

    checkpoint = acc_end;

    vector<Fr> fx_copy(f.begin(), f.end());
    vector<Fr> tx_copy(t.begin(), t.end() - 1);
    vector<Fr> twx_copy(t.begin() + 1, t.end());
    vector<Fr> h1x_copy(s.begin(), s.begin() + N);
    vector<Fr> h1wx_copy(s.begin() + 1, s.begin() + N + 1);
    vector<Fr> h2x_copy(s.begin() + N, s.begin() + 2 * N);
    vector<Fr> h2wx_copy(s.begin() + N + 1, s.begin() + (2 * N + 1));

    polyff fx, tx, twx, h1x, h1wx, h2x, h2wx;

    if (fx_copy.size() > std::pow(2, 15))
    {

        std::vector<std::thread> threads;

        threads.emplace_back(inttAndCreatePolyff, std::ref(fx_copy), std::ref(fx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(tx_copy), std::ref(tx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(twx_copy), std::ref(twx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h1x_copy), std::ref(h1x));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h1wx_copy), std::ref(h1wx));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h2x_copy), std::ref(h2x));
        threads.emplace_back(inttAndCreatePolyff, std::ref(h2wx_copy), std::ref(h2wx));

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {

        inttAndCreatePolyff(fx_copy, fx);
        inttAndCreatePolyff(tx_copy, tx);
        inttAndCreatePolyff(twx_copy, twx);
        inttAndCreatePolyff(h1x_copy, h1x);
        inttAndCreatePolyff(h1wx_copy, h1wx);
        inttAndCreatePolyff(h2x_copy, h2x);
        inttAndCreatePolyff(h2wx_copy, h2wx);
    }

    auto end = high_resolution_clock::now();

    return make_tuple(fx, tx, twx, h1x, h1wx, h2x, h2wx);
}

tuple<polyff, polyff>
Plookup_poly_z(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma)
{
    auto start = high_resolution_clock::now();
    auto checkpoint = start;

    size_t m = f.size();
    size_t n = t.size();

    size_t N = Max(m, n);
    N = nextPowOf2(N);

    if (m < N)
    {
        Fr lastElementOfT = t.back();
        f.resize(N, lastElementOfT);
    }

    if (n < N)
    {
        Fr lastElementOfT = t.back();
        t.resize(N, lastElementOfT);
    }

    checkpoint = high_resolution_clock::now();

    unordered_map<Fr, size_t> frequencyMap;
    for (auto &element : f)
    {
        frequencyMap[element]++;
    }

    auto map_end1 = high_resolution_clock::now();

    checkpoint = map_end1;

    vector<Fr> s;
    s.resize(f.size() + t.size());

    size_t s_index = 0;
    for (auto &element : t)
    {
        s[s_index++] = element;
        if (frequencyMap.find(element) != frequencyMap.end())
        {
            size_t count = frequencyMap[element];
            frequencyMap[element] = 0;
            for (size_t i = 0; i < count; ++i)
            {
                s[s_index++] = element;
            }
        }
    }

    auto map_end2 = high_resolution_clock::now();

    checkpoint = map_end2;

    t.push_back(t.back());
    s.push_back(t.back());

    vector<Fr> accs;

    accs.resize(N);
    Fr acc = 1;
    accs[0] = 1;

    for (size_t i = 1; i < N + 1; ++i)
    {
        acc *= ((1 + beta) * (gamma + f[i - 1]) * (gamma * (1 + beta) + t[i - 1] + beta * t[i])) /
               ((gamma * (1 + beta) + s[i - 1] + beta * s[i]) * (gamma * (1 + beta) + s[i + N - 1] + beta * s[i + N]));

        if (i < (N))
        {
            accs[i] = acc;
        }
    }

    auto acc_end = high_resolution_clock::now();

    checkpoint = acc_end;

    vector<Fr> accs_copy = accs;
    Ntt<Fr> ntt;
    ntt.init(nextPowOf2(accs_copy.size()));

    ntt.intt(&accs_copy[0]);

    polyff zx;
    zx.data = accs_copy;
    vector<Fr> zwx_copy(accs.begin() + 1, accs.end());
    zwx_copy.push_back(Fr(1));

    polyff zwx, fx, tx, twx, h1x, h1wx, h2x, h2wx;

    if (zwx_copy.size() > std::pow(2, 32))
    {

        std::vector<std::thread> threads;

        threads.emplace_back(inttAndCreatePolyff, std::ref(zwx_copy), std::ref(zwx));

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {

        inttAndCreatePolyff(zwx_copy, zwx);
    }

    auto end = high_resolution_clock::now();

    return make_tuple(zx, zwx);
}
