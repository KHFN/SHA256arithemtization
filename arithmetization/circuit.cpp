#include "prover.hpp"
#include "poly.hpp"
#include "matrixFr.hpp"
#include "prover.hpp"
#include "SHA256.hpp"

using namespace std;

class circuit
{

public:
    uint32_t numberofselector;
    uint32_t numberofwire;

    matrixFr excutiontrace;
    matrixFr table;

    vector<vector<vector<uint32_t>>> permutation_pos_interior;

    vector<polyff> poly_table;
    vector<polyff> poly_selector;
    vector<polyff> poly_permutation;
    vector<polyff> poly_wire;

    polyff zx, zwx, fx, tx, twx, h1x, h1wx, h2x, h2wx;
    polyff fx_t;
    polyff accx;
    polyff accwx;

    polyff tz;

    circuit() = default;
    circuit(uint32_t row, uint32_t numberofselector, uint32_t numberofwire)
    {

        matrixFr A(row, (numberofselector + numberofwire));
        excutiontrace = A;
    }

    void init_table(string tablefilename)
    {

        std::vector<std::vector<Fr>> loadedData;
        loadFromCSV(loadedData, tablefilename);
        table.data = loadedData;
    }

    void routing(vector<vector<uint32_t>> pos)
    {

        if (pos.size() != 2 || pos[0].size() != 2 || pos[1].size() != 2)
        {
            throw std::invalid_argument("Position vectors must have exactly two elements (row and column).");
        }

        excutiontrace.duplicateElement(pos[0], pos[1]);

        permutation_pos_interior.push_back(pos);
    }
};

Fr xor3Fr(Fr a, Fr b, Fr c)
{

    int64_t xor3 = (a.getInt64()) ^ (b.getInt64()) ^ (c.getInt64());

    return Fr(xor3);
}

Fr choiceFr(Fr a, Fr b, Fr c)
{

    int64_t choice = ((a.getInt64()) & (b.getInt64())) ^ (~(a.getInt64()) & (c.getInt64()));

    return Fr(choice);
}

Fr majorFr(Fr a, Fr b, Fr c)
{

    int64_t major = ((a.getInt64()) & (b.getInt64())) ^ ((a.getInt64()) & (c.getInt64())) ^ ((b.getInt64()) & (c.getInt64()));

    return Fr(major);
}

tuple<vector<Fr>, vector<Fr>>
splitBinaryString(unsigned int num, vector<int> lengths)
{

    reverse(lengths.begin(), lengths.end());
    string binary = toBinaryString(num);
    int binaryLength = binary.length();
    int totalLength = 0;
    for (int length : lengths)
    {
        totalLength += length;
    }

    if (totalLength < binaryLength)
    {
        throw runtime_error("Error: 분할 길이 합이 이진수의 길이보다 작습니다.");
    }

    if (totalLength > binaryLength)
    {
        binary = string(totalLength - binaryLength, '0') + binary;
    }

    reverse(binary.begin(), binary.end());

    vector<Fr> split;
    vector<int> starts;
    vector<Fr> coeff;

    int start = 0;
    for (int i = lengths.size() - 1; i >= 0; --i)
    {
        int length = lengths[i];
        string part = binary.substr(start, length);
        reverse(part.begin(), part.end());
        split.push_back(Fr(part, 2));
        coeff.push_back(Fr(1 << start));
        start += length; //
    }

    Fr acc = 0;
    for (int i = 0; i < split.size(); i++)
    {

        acc += split[i] * (coeff[i]);
    }

    split.push_back(acc);

    return {split, coeff};
}

tuple<vector<Fr>, vector<Fr>>
Gate_decomposition(uint32_t chunck, vector<int> lengths)
{

    auto [a, b] = splitBinaryString(chunck, lengths);

    vector<Fr> subcoeff1(b.begin(), b.begin() + 6);
    vector<Fr> subcoeff2(b.begin() + 6, b.end());
    subcoeff2.insert(subcoeff2.begin(), 1);
    subcoeff2.resize(6, 0);

    vector<Fr> suba1(a.begin(), a.begin() + 6);
    suba1.push_back(dotProduct(suba1, subcoeff1));

    vector<Fr> suba2(a.begin() + 6, a.end() - 1);
    suba2.insert(suba2.begin(), suba1.back());
    suba2.resize(6, 0);
    suba2.push_back(dotProduct(suba2, subcoeff2));

    vector<Fr> r1 = concatenateVectors(suba1, subcoeff1);
    vector<Fr> r2 = concatenateVectors(suba2, subcoeff2);

    return make_tuple(r1, r2);
}

tuple<vector<Fr>, vector<Fr>>
split_ss0_rotr7(uint32_t chunck)
{

    vector<int> lengths = {4, 3, 2, 4, 4, 4, 3, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_ss0_rotr18(uint32_t chunck)
{

    vector<int> lengths = {3, 4, 4, 4, 3, 2, 4, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_ss0_shr3(uint32_t chunck)
{

    vector<int> lengths = {3, 2, 4, 4, 4, 3, 4, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_ss1_rotr17(uint32_t chunck)
{

    vector<int> lengths = {3, 4, 2, 4, 4, 2, 3, 4, 4, 2};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_ss1_rotr19(uint32_t chunck)
{

    vector<int> lengths = {2, 3, 4, 2, 4, 4, 2, 3, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_ss1_shr10(uint32_t chunck)
{

    vector<int> lengths = {2, 4, 4, 2, 3, 4, 4, 2, 3, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs0_rotr6(uint32_t chunck)
{

    vector<int> lengths = {2, 4, 3, 4, 2, 4, 4, 4, 2, 3};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs0_rotr11(uint32_t chunck)
{

    vector<int> lengths = {2, 3, 2, 4, 3, 4, 2, 4, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs0_rotr25(uint32_t chunck)
{

    vector<int> lengths = {2, 4, 4, 4, 2, 3, 2, 4, 3, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs1_rotr2(uint32_t chunck)
{

    vector<int> lengths = {2, 2, 4, 4, 2, 3, 4, 3, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs1_rotr13(uint32_t chunck)
{

    vector<int> lengths = {3, 4, 4, 2, 2, 4, 4, 2, 3, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<Fr>, vector<Fr>>
split_bs1_rotr22(uint32_t chunck)
{

    vector<int> lengths = {2, 3, 4, 3, 4, 4, 2, 2, 4, 4};
    return Gate_decomposition(chunck, lengths);
}

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>>
bigsigma_1(uint32_t numberofselector, uint32_t numberofwire, uint32_t chunk)
{

    uint32_t row = 18;
    circuit A(row, numberofselector, numberofwire);

    auto [a, b] = split_bs0_rotr6(chunk);
    auto [c, d] = split_bs0_rotr11(chunk);
    auto [e, f] = split_bs0_rotr25(chunk);

    A.excutiontrace.replaceSubmatrix({0, numberofselector}, {a});
    A.excutiontrace.replaceSubmatrix({1, numberofselector}, {b});
    A.excutiontrace[0][0] = 1;
    A.excutiontrace[1][0] = 1;
    A.routing({{0, 9}, {1, numberofselector}});

    A.excutiontrace.replaceSubmatrix({2, numberofselector}, {c});
    A.excutiontrace.replaceSubmatrix({3, numberofselector}, {d});
    A.excutiontrace[2][0] = 1;
    A.excutiontrace[3][0] = 1;
    A.routing({{2, 9}, {3, numberofselector}});

    A.excutiontrace.replaceSubmatrix({4, numberofselector}, {e});
    A.excutiontrace.replaceSubmatrix({5, numberofselector}, {f});
    A.excutiontrace[4][0] = 1;
    A.excutiontrace[5][0] = 1;
    A.routing({{4, 9}, {5, numberofselector}});

    A.excutiontrace[6][2] = 1;
    A.excutiontrace[7][2] = 1;
    A.excutiontrace[8][2] = 1;
    A.excutiontrace[9][2] = 1;
    A.excutiontrace[10][2] = 1;
    A.excutiontrace[11][2] = 1;
    A.excutiontrace[12][2] = 1;
    A.excutiontrace[13][2] = 1;
    A.excutiontrace[14][2] = 1;
    A.excutiontrace[15][2] = 1;
    A.excutiontrace[16][0] = 1;
    A.excutiontrace[17][0] = 1;

    A.routing({{0, 2 + 3}, {6, numberofselector}});
    A.routing({{0, 3 + 3}, {7, numberofselector}});
    A.routing({{0, 4 + 3}, {8, numberofselector}});
    A.routing({{0, 5 + 3}, {9, numberofselector}});
    A.routing({{1, 1 + 3}, {10, numberofselector}});
    A.routing({{1, 2 + 3}, {11, numberofselector}});
    A.routing({{1, 3 + 3}, {12, numberofselector}});
    A.routing({{1, 4 + 3}, {13, numberofselector}});
    A.routing({{0, 0 + 3}, {14, numberofselector}});
    A.routing({{0, 1 + 3}, {15, numberofselector}});

    A.routing({{2, 4 + 3}, {6, numberofselector + 1}});
    A.routing({{2, 5 + 3}, {7, numberofselector + 1}});
    A.routing({{3, 1 + 3}, {8, numberofselector + 1}});
    A.routing({{3, 2 + 3}, {9, numberofselector + 1}});
    A.routing({{3, 3 + 3}, {10, numberofselector + 1}});
    A.routing({{3, 4 + 3}, {11, numberofselector + 1}});
    A.routing({{2, 0 + 3}, {12, numberofselector + 1}});
    A.routing({{2, 1 + 3}, {13, numberofselector + 1}});
    A.routing({{2, 2 + 3}, {14, numberofselector + 1}});
    A.routing({{2, 3 + 3}, {15, numberofselector + 1}});

    A.routing({{5, 3 + 3}, {6, numberofselector + 2}});
    A.routing({{5, 4 + 3}, {7, numberofselector + 2}});
    A.routing({{4, 0 + 3}, {8, numberofselector + 2}});
    A.routing({{4, 1 + 3}, {9, numberofselector + 2}});
    A.routing({{4, 2 + 3}, {10, numberofselector + 2}});
    A.routing({{4, 3 + 3}, {11, numberofselector + 2}});
    A.routing({{4, 4 + 3}, {12, numberofselector + 2}});
    A.routing({{4, 5 + 3}, {13, numberofselector + 2}});
    A.routing({{5, 1 + 3}, {14, numberofselector + 2}});
    A.routing({{5, 2 + 3}, {15, numberofselector + 2}});

    A.excutiontrace[6][6] = xor3Fr(A.excutiontrace[6][3], A.excutiontrace[6][4], A.excutiontrace[6][5]);
    A.excutiontrace[7][6] = xor3Fr(A.excutiontrace[7][3], A.excutiontrace[7][4], A.excutiontrace[7][5]);
    A.excutiontrace[8][6] = xor3Fr(A.excutiontrace[8][3], A.excutiontrace[8][4], A.excutiontrace[8][5]);
    A.excutiontrace[9][6] = xor3Fr(A.excutiontrace[9][3], A.excutiontrace[9][4], A.excutiontrace[9][5]);
    A.excutiontrace[10][6] = xor3Fr(A.excutiontrace[10][3], A.excutiontrace[10][4], A.excutiontrace[10][5]);
    A.excutiontrace[11][6] = xor3Fr(A.excutiontrace[11][3], A.excutiontrace[11][4], A.excutiontrace[11][5]);
    A.excutiontrace[12][6] = xor3Fr(A.excutiontrace[12][3], A.excutiontrace[12][4], A.excutiontrace[12][5]);
    A.excutiontrace[13][6] = xor3Fr(A.excutiontrace[13][3], A.excutiontrace[13][4], A.excutiontrace[13][5]);
    A.excutiontrace[14][6] = xor3Fr(A.excutiontrace[14][3], A.excutiontrace[14][4], A.excutiontrace[14][5]);
    A.excutiontrace[15][6] = xor3Fr(A.excutiontrace[15][3], A.excutiontrace[15][4], A.excutiontrace[15][5]);

    A.routing({{6, 6}, {16, numberofselector}});
    A.routing({{7, 6}, {16, numberofselector + 1}});
    A.routing({{8, 6}, {16, numberofselector + 2}});
    A.routing({{9, 6}, {16, numberofselector + 3}});
    A.routing({{10, 6}, {16, numberofselector + 4}});
    A.routing({{11, 6}, {16, numberofselector + 5}});
    A.routing({{12, 6}, {17, numberofselector + 1}});
    A.routing({{13, 6}, {17, numberofselector + 2}});
    A.routing({{14, 6}, {17, numberofselector + 3}});
    A.routing({{15, 6}, {17, numberofselector + 4}});

    A.excutiontrace.replaceSubmatrix({16, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 3), Fr(1 << 7), Fr(1 << 9), Fr(1 << 13), Fr(1 << 17)}});
    A.excutiontrace.replaceSubmatrix({17, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 21), Fr(1 << 23), Fr(1 << 26), Fr(1 << 28), 0}});

    vector<Fr> e5(A.excutiontrace[16].begin() + numberofselector, A.excutiontrace[16].begin() + numberofselector + 6);
    vector<Fr> e6(A.excutiontrace[16].begin() + numberofselector + 7, A.excutiontrace[16].begin() + numberofselector + 7 + 6);
    A.excutiontrace[16][numberofselector + 6] = dotProduct(e5, e6);
    A.routing({{16, numberofselector + 6}, {17, numberofselector}});

    vector<Fr> e7(A.excutiontrace[17].begin() + numberofselector, A.excutiontrace[17].begin() + numberofselector + 6);
    vector<Fr> e8(A.excutiontrace[17].begin() + numberofselector + 7, A.excutiontrace[17].begin() + numberofselector + 7 + 6);

    A.excutiontrace[17][numberofselector + 6] = dotProduct(e7, e8);

    return make_tuple(A.excutiontrace.data, A.permutation_pos_interior);
}

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>>
bigsigma_0(uint32_t numberofselector, uint32_t numberofwire, uint32_t chunk)
{

    uint32_t row = 18;
    circuit A(row, numberofselector, numberofwire);

    auto [a, b] = split_bs1_rotr2(chunk);
    auto [c, d] = split_bs1_rotr13(chunk);
    auto [e, f] = split_bs1_rotr22(chunk);

    A.excutiontrace.replaceSubmatrix({0, numberofselector}, {a});
    A.excutiontrace.replaceSubmatrix({1, numberofselector}, {b});
    A.excutiontrace[0][0] = 1;
    A.excutiontrace[1][0] = 1;
    A.routing({{0, 9}, {1, numberofselector}});

    A.excutiontrace.replaceSubmatrix({2, numberofselector}, {c});
    A.excutiontrace.replaceSubmatrix({3, numberofselector}, {d});
    A.excutiontrace[2][0] = 1;
    A.excutiontrace[3][0] = 1;
    A.routing({{2, 9}, {3, numberofselector}});

    A.excutiontrace.replaceSubmatrix({4, numberofselector}, {e});
    A.excutiontrace.replaceSubmatrix({5, numberofselector}, {f});
    A.excutiontrace[4][0] = 1;
    A.excutiontrace[5][0] = 1;
    A.routing({{4, 9}, {5, numberofselector}});

    A.excutiontrace[6][2] = 1;
    A.excutiontrace[7][2] = 1;
    A.excutiontrace[8][2] = 1;
    A.excutiontrace[9][2] = 1;
    A.excutiontrace[10][2] = 1;
    A.excutiontrace[11][2] = 1;
    A.excutiontrace[12][2] = 1;
    A.excutiontrace[13][2] = 1;
    A.excutiontrace[14][2] = 1;
    A.excutiontrace[15][2] = 1;
    A.excutiontrace[16][0] = 1;
    A.excutiontrace[17][0] = 1;

    A.routing({{0, 1 + 3}, {6, numberofselector}});
    A.routing({{0, 2 + 3}, {7, numberofselector}});
    A.routing({{0, 3 + 3}, {8, numberofselector}});
    A.routing({{0, 4 + 3}, {9, numberofselector}});
    A.routing({{0, 5 + 3}, {10, numberofselector}});
    A.routing({{1, 1 + 3}, {11, numberofselector}});
    A.routing({{1, 2 + 3}, {12, numberofselector}});
    A.routing({{1, 3 + 3}, {13, numberofselector}});
    A.routing({{1, 4 + 3}, {14, numberofselector}});
    A.routing({{0, 0 + 3}, {15, numberofselector}});

    A.routing({{2, 4 + 3}, {6, numberofselector + 1}});
    A.routing({{2, 5 + 3}, {7, numberofselector + 1}});
    A.routing({{3, 1 + 3}, {8, numberofselector + 1}});
    A.routing({{3, 2 + 3}, {9, numberofselector + 1}});
    A.routing({{3, 3 + 3}, {10, numberofselector + 1}});
    A.routing({{3, 4 + 3}, {11, numberofselector + 1}});
    A.routing({{2, 0 + 3}, {12, numberofselector + 1}});
    A.routing({{2, 1 + 3}, {13, numberofselector + 1}});
    A.routing({{2, 2 + 3}, {14, numberofselector + 1}});
    A.routing({{2, 3 + 3}, {15, numberofselector + 1}});

    A.routing({{5, 2 + 3}, {6, numberofselector + 2}});
    A.routing({{5, 3 + 3}, {7, numberofselector + 2}});
    A.routing({{5, 4 + 3}, {8, numberofselector + 2}});
    A.routing({{4, 0 + 3}, {9, numberofselector + 2}});
    A.routing({{4, 1 + 3}, {10, numberofselector + 2}});
    A.routing({{4, 2 + 3}, {11, numberofselector + 2}});
    A.routing({{4, 3 + 3}, {12, numberofselector + 2}});
    A.routing({{4, 4 + 3}, {13, numberofselector + 2}});
    A.routing({{4, 5 + 3}, {14, numberofselector + 2}});
    A.routing({{5, 1 + 3}, {15, numberofselector + 2}});

    A.excutiontrace[6][6] = xor3Fr(A.excutiontrace[6][3], A.excutiontrace[6][4], A.excutiontrace[6][5]);
    A.excutiontrace[7][6] = xor3Fr(A.excutiontrace[7][3], A.excutiontrace[7][4], A.excutiontrace[7][5]);
    A.excutiontrace[8][6] = xor3Fr(A.excutiontrace[8][3], A.excutiontrace[8][4], A.excutiontrace[8][5]);
    A.excutiontrace[9][6] = xor3Fr(A.excutiontrace[9][3], A.excutiontrace[9][4], A.excutiontrace[9][5]);
    A.excutiontrace[10][6] = xor3Fr(A.excutiontrace[10][3], A.excutiontrace[10][4], A.excutiontrace[10][5]);
    A.excutiontrace[11][6] = xor3Fr(A.excutiontrace[11][3], A.excutiontrace[11][4], A.excutiontrace[11][5]);
    A.excutiontrace[12][6] = xor3Fr(A.excutiontrace[12][3], A.excutiontrace[12][4], A.excutiontrace[12][5]);
    A.excutiontrace[13][6] = xor3Fr(A.excutiontrace[13][3], A.excutiontrace[13][4], A.excutiontrace[13][5]);
    A.excutiontrace[14][6] = xor3Fr(A.excutiontrace[14][3], A.excutiontrace[14][4], A.excutiontrace[14][5]);
    A.excutiontrace[15][6] = xor3Fr(A.excutiontrace[15][3], A.excutiontrace[15][4], A.excutiontrace[15][5]);

    A.routing({{6, 6}, {16, numberofselector}});
    A.routing({{7, 6}, {16, numberofselector + 1}});
    A.routing({{8, 6}, {16, numberofselector + 2}});
    A.routing({{9, 6}, {16, numberofselector + 3}});
    A.routing({{10, 6}, {16, numberofselector + 4}});
    A.routing({{11, 6}, {16, numberofselector + 5}});
    A.routing({{12, 6}, {17, numberofselector + 1}});
    A.routing({{13, 6}, {17, numberofselector + 2}});
    A.routing({{14, 6}, {17, numberofselector + 3}});
    A.routing({{15, 6}, {17, numberofselector + 4}});

    A.excutiontrace.replaceSubmatrix({16, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 2), Fr(1 << 6), Fr(1 << 10), Fr(1 << 12), Fr(1 << 15)}});
    A.excutiontrace.replaceSubmatrix({17, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 19), Fr(1 << 22), Fr(1 << 26), Fr(1 << 30), 0}});

    vector<Fr> e5(A.excutiontrace[16].begin() + numberofselector, A.excutiontrace[16].begin() + numberofselector + 6);
    vector<Fr> e6(A.excutiontrace[16].begin() + numberofselector + 7, A.excutiontrace[16].begin() + numberofselector + 7 + 6);
    A.excutiontrace[16][numberofselector + 6] = dotProduct(e5, e6);
    A.routing({{16, numberofselector + 6}, {17, numberofselector}});

    vector<Fr> e7(A.excutiontrace[17].begin() + numberofselector, A.excutiontrace[17].begin() + numberofselector + 6);
    vector<Fr> e8(A.excutiontrace[17].begin() + numberofselector + 7, A.excutiontrace[17].begin() + numberofselector + 7 + 6);

    A.excutiontrace[17][numberofselector + 6] = dotProduct(e7, e8);

    return make_tuple(A.excutiontrace.data, A.permutation_pos_interior);
}

tuple<vector<Fr>, vector<Fr>>
bitdecomposition_4(uint32_t chunk)
{

    vector<int> lengths = {4, 4, 4, 4, 4, 4, 4, 4};
    auto [a, b] = Gate_decomposition(chunk, lengths);

    return {a, b};
}

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>>
majorgate(uint32_t numberofselector, uint32_t numberofwire, uint32_t A, uint32_t B, uint32_t C)
{

    uint32_t row = 16;
    circuit M(row, numberofselector, numberofwire);
    M.excutiontrace[0][0] = 1;
    M.excutiontrace[1][0] = 1;
    M.excutiontrace[2][0] = 1;
    M.excutiontrace[3][0] = 1;
    M.excutiontrace[4][0] = 1;
    M.excutiontrace[5][0] = 1;
    M.excutiontrace[6][2] = 1;
    M.excutiontrace[7][2] = 1;
    M.excutiontrace[8][2] = 1;
    M.excutiontrace[9][2] = 1;
    M.excutiontrace[10][2] = 1;
    M.excutiontrace[11][2] = 1;
    M.excutiontrace[12][2] = 1;
    M.excutiontrace[13][2] = 1;
    M.excutiontrace[14][0] = 1;
    M.excutiontrace[15][0] = 1;

    M.excutiontrace[6][7] = 1;
    M.excutiontrace[7][7] = 1;
    M.excutiontrace[8][7] = 1;
    M.excutiontrace[9][7] = 1;
    M.excutiontrace[10][7] = 1;
    M.excutiontrace[11][7] = 1;
    M.excutiontrace[12][7] = 1;
    M.excutiontrace[13][7] = 1;

    auto [a, b] = bitdecomposition_4(A);
    M.excutiontrace.replaceSubmatrix({0, numberofselector}, {a, b});
    vector<Fr> m1(a.begin(), a.begin() + 6);
    vector<Fr> m2(b.begin() + 1, b.begin() + 3);
    vector<Fr> m = concatenateVectors(m1, m2);
    vector<vector<Fr>> tmp;
    tmp.push_back(m);
    vector<vector<Fr>> mt;

    auto [c, d] = bitdecomposition_4(B);
    M.excutiontrace.replaceSubmatrix({2, numberofselector}, {c, d});
    m.clear();
    m1.clear();
    m2.clear();
    m1.assign(c.begin(), c.begin() + 6);
    m2.assign(d.begin() + 1, d.begin() + 3);
    m = concatenateVectors(m1, m2);
    tmp.push_back(m);

    auto [e, f] = bitdecomposition_4(C);
    M.excutiontrace.replaceSubmatrix({4, numberofselector}, {e, f});
    m.clear();
    m1.clear();
    m2.clear();
    m1.assign(e.begin(), e.begin() + 6);
    m2.assign(f.begin() + 1, f.begin() + 3);
    m = concatenateVectors(m1, m2);
    tmp.push_back(m);
    mt = transpose(tmp);

    M.excutiontrace.replaceSubmatrix({6, numberofselector}, mt);

    for (int i = 6; i < 14; i++)
    {

        M.excutiontrace[i][6] = majorFr(M.excutiontrace[i][3], M.excutiontrace[i][4], M.excutiontrace[i][5]);
    }

    M.routing({{0 + 6, 6}, {14, numberofselector}});
    M.routing({{1 + 6, 6}, {14, numberofselector + 1}});
    M.routing({{2 + 6, 6}, {14, numberofselector + 2}});
    M.routing({{3 + 6, 6}, {14, numberofselector + 3}});
    M.routing({{4 + 6, 6}, {14, numberofselector + 4}});
    M.routing({{5 + 6, 6}, {14, numberofselector + 5}});
    M.routing({{6 + 6, 6}, {15, numberofselector + 1}});
    M.routing({{7 + 6, 6}, {15, numberofselector + 2}});

    M.excutiontrace.replaceSubmatrix({14, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 4), Fr(1 << 8), Fr(1 << 12), Fr(1 << 16), Fr(1 << 20)}});
    M.excutiontrace.replaceSubmatrix({15, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 24), Fr(1 << 28), 0, 0, 0}});

    vector<Fr> e5(M.excutiontrace[14].begin() + numberofselector, M.excutiontrace[14].begin() + numberofselector + 6);
    vector<Fr> e6(M.excutiontrace[14].begin() + numberofselector + 7, M.excutiontrace[14].begin() + numberofselector + 7 + 6);
    M.excutiontrace[14][numberofselector + 6] = dotProduct(e5, e6);
    M.routing({{14, numberofselector + 6}, {15, numberofselector}});

    vector<Fr> e7(M.excutiontrace[15].begin() + numberofselector, M.excutiontrace[15].begin() + numberofselector + 6);
    vector<Fr> e8(M.excutiontrace[15].begin() + numberofselector + 7, M.excutiontrace[15].begin() + numberofselector + 7 + 6);

    M.excutiontrace[15][numberofselector + 6] = dotProduct(e7, e8);

    return make_tuple(M.excutiontrace.data, M.permutation_pos_interior);
}

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>>
choicegate(uint32_t numberofselector, uint32_t numberofwire, uint32_t A, uint32_t B, uint32_t C)
{

    uint32_t row = 16;
    circuit M(row, numberofselector, numberofwire);
    M.excutiontrace[0][0] = 1;
    M.excutiontrace[1][0] = 1;
    M.excutiontrace[2][0] = 1;
    M.excutiontrace[3][0] = 1;
    M.excutiontrace[4][0] = 1;
    M.excutiontrace[5][0] = 1;
    M.excutiontrace[6][2] = 1;
    M.excutiontrace[7][2] = 1;
    M.excutiontrace[8][2] = 1;
    M.excutiontrace[9][2] = 1;
    M.excutiontrace[10][2] = 1;
    M.excutiontrace[11][2] = 1;
    M.excutiontrace[12][2] = 1;
    M.excutiontrace[13][2] = 1;
    M.excutiontrace[14][0] = 1;
    M.excutiontrace[15][0] = 1;

    M.excutiontrace[6][7] = 2;
    M.excutiontrace[7][7] = 2;
    M.excutiontrace[8][7] = 2;
    M.excutiontrace[9][7] = 2;
    M.excutiontrace[10][7] = 2;
    M.excutiontrace[11][7] = 2;
    M.excutiontrace[12][7] = 2;
    M.excutiontrace[13][7] = 2;

    auto [a, b] = bitdecomposition_4(A);
    M.excutiontrace.replaceSubmatrix({0, numberofselector}, {a, b});
    vector<Fr> m1(a.begin(), a.begin() + 6);
    vector<Fr> m2(b.begin() + 1, b.begin() + 3);
    vector<Fr> m = concatenateVectors(m1, m2);
    vector<vector<Fr>> tmp;
    tmp.push_back(m);
    vector<vector<Fr>> mt;

    auto [c, d] = bitdecomposition_4(B);
    M.excutiontrace.replaceSubmatrix({2, numberofselector}, {c, d});
    m.clear();
    m1.clear();
    m2.clear();
    m1.assign(c.begin(), c.begin() + 6);
    m2.assign(d.begin() + 1, d.begin() + 3);
    m = concatenateVectors(m1, m2);
    tmp.push_back(m);

    auto [e, f] = bitdecomposition_4(C);
    M.excutiontrace.replaceSubmatrix({4, numberofselector}, {e, f});
    m.clear();
    m1.clear();
    m2.clear();
    m1.assign(e.begin(), e.begin() + 6);
    m2.assign(f.begin() + 1, f.begin() + 3);
    m = concatenateVectors(m1, m2);
    tmp.push_back(m);
    mt = transpose(tmp);

    M.excutiontrace.replaceSubmatrix({6, numberofselector}, mt);

    for (int i = 6; i < 14; i++)
    {

        M.excutiontrace[i][6] = choiceFr(M.excutiontrace[i][3], M.excutiontrace[i][4], M.excutiontrace[i][5]);
    }

    M.routing({{0 + 6, 6}, {14, numberofselector}});
    M.routing({{1 + 6, 6}, {14, numberofselector + 1}});
    M.routing({{2 + 6, 6}, {14, numberofselector + 2}});
    M.routing({{3 + 6, 6}, {14, numberofselector + 3}});
    M.routing({{4 + 6, 6}, {14, numberofselector + 4}});
    M.routing({{5 + 6, 6}, {14, numberofselector + 5}});
    M.routing({{6 + 6, 6}, {15, numberofselector + 1}});
    M.routing({{7 + 6, 6}, {15, numberofselector + 2}});

    M.excutiontrace.replaceSubmatrix({14, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 4), Fr(1 << 8), Fr(1 << 12), Fr(1 << 16), Fr(1 << 20)}});
    M.excutiontrace.replaceSubmatrix({15, numberofselector + 7},
                                     {{Fr(1 << 0), Fr(1 << 24), Fr(1 << 28), 0, 0, 0}});

    vector<Fr> e5(M.excutiontrace[14].begin() + numberofselector, M.excutiontrace[14].begin() + numberofselector + 6);
    vector<Fr> e6(M.excutiontrace[14].begin() + numberofselector + 7, M.excutiontrace[14].begin() + numberofselector + 7 + 6);
    M.excutiontrace[14][numberofselector + 6] = dotProduct(e5, e6);
    M.routing({{14, numberofselector + 6}, {15, numberofselector}});

    vector<Fr> e7(M.excutiontrace[15].begin() + numberofselector, M.excutiontrace[15].begin() + numberofselector + 6);
    vector<Fr> e8(M.excutiontrace[15].begin() + numberofselector + 7, M.excutiontrace[15].begin() + numberofselector + 7 + 6);

    M.excutiontrace[15][numberofselector + 6] = dotProduct(e7, e8);

    return make_tuple(M.excutiontrace.data, M.permutation_pos_interior);
}

vector<vector<Fr>>
add5gate(uint32_t numberofselector, uint32_t numberofwire, vector<Fr> wire)
{

    if (wire.size() > 5)
    {
        throw std::invalid_argument("add5gate : 벡터 길이가 5이하여야 합니다.");
    }
    vector<vector<Fr>> result(1, vector<Fr>((numberofselector + numberofwire), 0));

    result[0][1] = 1;

    for (uint8_t i = 0; i < 5; i++)
    {

        result[0][numberofselector + i] = wire[i];
    }

    uint64_t W_new = (result[0][numberofselector] + result[0][numberofselector + 1] + result[0][numberofselector + 2] + result[0][numberofselector + 3] + result[0][numberofselector + 4]).getInt64();
    uint64_t W_new_copy = W_new;

    W_new = W_new & (4294967295);
    result[0][numberofselector + 6] = Fr(W_new);
    result[0][numberofselector + 5] = result[0][numberofselector + 6] - Fr(W_new_copy);

    result[0][1] = 1;
    return result;
}
