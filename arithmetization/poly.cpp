
#include <mcl/bls12_381.hpp>
#include <mcl/ntt.hpp>
#include <iostream>
#include <cmath>
#include <thread>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <mcl/lagrange.hpp>

using namespace mcl;
using namespace mcl::bn;
using namespace std;

typedef std::vector<Fr> FrVec;

template <class T>
inline void initPowSeq(T *sec, const T &u, size_t n)
{
    if (n == 0)
        return;
    sec[0] = 1;
    if (n == 1)
        return;
    sec[1] = u;
    for (size_t i = 2; i < n; i++)
    {
        T::mul(sec[i], sec[i - 1], u);
    }
}

bool IsPowerOfTwo(uint64_t v)
{
    return (v & (v - 1)) == 0;
}

uint32_t Max(uint32_t x, uint32_t y)
{
    return (x < y) ? y : x;
}

uint32_t Min(uint32_t x, uint32_t y)
{
    return (x < y) ? x : y;
}

uint64_t nextPowOf2(const uint64_t v)
{
    if (v == 0)
    {
        return 1;
    }
    return uint64_t(1) << static_cast<uint64_t>(log2(v - 1) + 1);
}

uint32_t nextPowOf2_32(const uint64_t v)
{
    if (v == 0)
    {
        return 1;
    }
    return uint32_t(1) << static_cast<uint64_t>(log2(v - 1) + 1);
}

void parallelComposeWithMonomial(Fr *resultData, const Fr *inputData, const Fr *powSeq, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        resultData[i] = inputData[i] * powSeq[i];
    }
}

class polyff
{
public:
    vector<Fr> data;

    polyff() = default;
    polyff(std::initializer_list<Fr> initList) : data(initList) {}
    polyff(int val) : data(1, Fr(val)) {}
    polyff(Fr val) : data(1, val) {}
    polyff(vector<Fr>::const_iterator begin, vector<Fr>::const_iterator end) : data(begin, end) {}
    polyff(size_t size, int initialValue)
    {
        data.resize(size, Fr(initialValue));
    }
    polyff(size_t size, Fr initialValue)
    {
        data.resize(size, initialValue);
    }

    size_t size() const
    {
        return data.size();
    }

    void resize(size_t newSize, const Fr &value = Fr(0))
    {
        data.resize(newSize, value);
    }

    Fr &operator[](size_t index)
    {
        if (index >= data.size())
        {
            throw std::out_of_range("Index out of range");
        }
        return data[index];
    }
    const Fr &operator[](size_t index) const
    {
        if (index >= data.size())
        {
            throw std::out_of_range("Index out of range");
        }
        return data[index];
    }

    Fr *getDataPtr()
    {
        if (data.empty())
        {
            throw std::runtime_error("Data is empty");
        }
        return &data[0];
    }

    const Fr *getDataPtr() const
    {
        if (data.empty())
        {
            throw std::runtime_error("Data is empty");
        }
        return &data[0];
    }

    polyff PolyCondense() const
    {
        uint32_t d = data.size();
        if (d == 0)
        {
            throw runtime_error("IsPolyZeroError");
        }

        uint32_t tmpd = d;
        for (uint32_t i = d - 1; i > 0; i--)
        {
            if (!data[i].isZero())
            {
                break;
            }
            tmpd--;
        }
        return polyff(data.begin(), data.begin() + tmpd);
    }

    bool operator==(const polyff &rhs) const
    {
        polyff p = this->PolyCondense();
        polyff q = rhs.PolyCondense();

        if (p.data.size() != q.data.size())
        {
            return false;
        }

        for (size_t i = 0; i < p.data.size(); i++)
        {
            if (p.data[i] != q.data[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!() const
    {
        if (data.empty())
        {
            throw runtime_error("IsPolyZeroError");
        }
        for (const auto &val : data)
        {
            if (!val.isZero())
            {
                return false;
            }
        }
        return true;
    }

    polyff operator-() const
    {
        polyff result;
        for (const auto &coeff : data)
        {
            result.data.push_back(-coeff);
        }
        return result;
    }

    polyff divideByXnMinusOne(int n) const
    {
        if (n <= 0)
        {
            throw std::runtime_error("n must be positive");
        }

        std::vector<Fr> divisor(n, Fr(0));
        divisor.push_back(Fr(1));
        divisor[0] = -Fr(1);

        polyff quotient, remainder = *this;
        while (remainder.size() >= divisor.size())
        {
            Fr factor = remainder.data.back();
            size_t index = remainder.size() - divisor.size();
            quotient.data.insert(quotient.data.begin(), factor);

            for (size_t i = 0; i < divisor.size(); ++i)
            {
                remainder[index + i] = remainder[index + i] - factor * divisor[i];
            }
            remainder.data.pop_back();
        }

        return quotient;
    }

    Fr evaluate(const Fr &x) const
    {
        Fr result(0);
        for (auto it = data.rbegin(); it != data.rend(); ++it)
        {
            result = result * x + *it;
        }
        return result;
    }

    polyff composeWithMonomial(const Fr &k) const
    {
        size_t n = this->size();
        FrVec powSeq(n);
        initPowSeq(powSeq.data(), k, n);

        polyff result = *this;
        result.resize(n);

        const size_t threshold = 131072;
        if (n >= threshold)
        {

            size_t hardwareConcurrency = std::thread::hardware_concurrency();
            size_t perThread = n / hardwareConcurrency;
            std::vector<std::thread> threads;

            for (size_t i = 0; i < hardwareConcurrency; ++i)
            {
                size_t start = i * perThread;
                size_t end = (i + 1 == hardwareConcurrency) ? n : (i + 1) * perThread;
                threads.emplace_back(parallelComposeWithMonomial, &result.data[0], &this->data[0], &powSeq[0], start, end);
            }

            for (auto &t : threads)
            {
                if (t.joinable())
                {
                    t.join();
                }
            }
        }
        else
        {

            for (size_t i = 0; i < n; ++i)
            {
                result.data[i] = this->data[i] * powSeq[i];
            }
        }

        return result;
    }

    polyff syntheticDivide(const Fr &root) const
    {
        if (data.empty())
        {
            throw std::runtime_error("Cannot divide an empty polynomial.");
        }

        std::vector<Fr> resultCoefficients(data.size() - 1);

        Fr carry = 0;
        for (size_t i = 0; i < data.size() - 1; ++i)
        {
            Fr updatedCoefficient = data[i] + carry;
            resultCoefficients[i] = updatedCoefficient;
            carry = updatedCoefficient * (-root);
        }

        return polyff(resultCoefficients.begin(), resultCoefficients.end()).PolyCondense();
    }

    polyff add(int value) const
    {
        polyff result = *this;
        if (result.size() == 0)
        {
            result.data.push_back(Fr(value));
        }
        else
        {
            result[0] = result[0] + value;
        }
        return result;
    }

    polyff add(uint32_t value) const
    {

        return this->add(static_cast<int>(value));
    }

    friend polyff operator+(const polyff &a, const polyff &b);
    friend polyff operator*(const polyff &a, const polyff &b);
    friend polyff operator/(const polyff &A, const polyff &B);
    friend polyff operator%(const polyff &A, const polyff &B);
    friend std::ostream &operator<<(std::ostream &os, const polyff &p);
};

void parallelAdd(polyff &result, const polyff &a, const polyff &b, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        result[i] = a[i] + b[i];
    }
}

polyff operator+(const polyff &a, const polyff &b)
{
    size_t ad = a.size();
    size_t bd = b.size();
    size_t n = std::max(ad, bd);

    polyff c(n, 0);

    if (n >= std::thread::hardware_concurrency() * 1500000)
    {
        std::vector<std::thread> threads;
        size_t perThread = n / std::thread::hardware_concurrency();
        size_t start = 0;

        for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i)
        {
            size_t end = std::min(start + perThread, n);
            threads.push_back(std::thread(parallelAdd, std::ref(c), std::ref(a), std::ref(b), start, end));
            start = end;
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {

        for (size_t i = 0; i < n; ++i)
        {
            if (i < ad)
                c[i] = c[i] + a[i];
            if (i < bd)
                c[i] = c[i] + b[i];
        }
    }

    return c.PolyCondense();
}

void parallelSub(polyff &result, const polyff &a, const polyff &b, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        if (i < a.size())
            result[i] = a[i];
        if (i < b.size())
            result[i] = result[i] - b[i];
    }
}

polyff operator-(const polyff &a, const polyff &b)
{
    size_t ad = a.size();
    size_t bd = b.size();
    size_t n = std::max(ad, bd);

    polyff c(n, 0);

    if (n >= std::thread::hardware_concurrency() * 150000)
    {
        std::vector<std::thread> threads;
        size_t perThread = n / std::thread::hardware_concurrency();
        size_t start = 0;

        for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i)
        {
            size_t end = std::min(start + perThread, n);
            threads.push_back(std::thread(parallelSub, std::ref(c), std::ref(a), std::ref(b), start, end));
            start = end;
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {
        for (size_t i = 0; i < n; ++i)
        {
            if (i < ad)
                c[i] = a[i];
            if (i < bd)
                c[i] = c[i] - b[i];
        }
    }

    return c.PolyCondense();
}

void parallelMultiply(Fr *result, const Fr *a, const Fr *b, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        result[i] = a[i] * b[i];
    }
}

polyff operator*(const polyff &a, const polyff &b)
{
    uint32_t N = a.size() + b.size() - 1;
    N = nextPowOf2(N);
    polyff a_copy = a.PolyCondense();
    polyff b_copy = b.PolyCondense();
    a_copy.resize(N, 0);
    b_copy.resize(N, 0);

    polyff res;
    res.resize(N, 0);

    Ntt<Fr> ntt;
    ntt.init(N);

    std::thread t1([&]
                   { ntt.ntt(&a_copy[0]); });
    std::thread t2([&]
                   { ntt.ntt(&b_copy[0]); });

    t1.join();
    t2.join();

    if (N >= std::thread::hardware_concurrency() * 1024)
    {
        std::vector<std::thread> threads;
        size_t perThread = N / std::thread::hardware_concurrency();
        size_t start = 0;

        for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i)
        {
            size_t end = std::min(start + perThread, static_cast<size_t>(N));
            threads.push_back(std::thread(parallelMultiply, &res[0], &a_copy[0], &b_copy[0], start, end));
            start = end;
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {
        for (size_t i = 0; i < N; ++i)
        {
            res[i] = a_copy[i] * b_copy[i];
        }
    }

    ntt.intt(&res[0]);
    return res.PolyCondense();
}

std::ostream &operator<<(std::ostream &os, const polyff &p)
{
    os << "[";
    for (size_t i = 0; i < p.data.size(); ++i)
    {
        if (i > 0)
        {
            os << ", ";
        }
        os << p.data[i];
    }
    os << "]";
    return os;
}

FrVec PolyCondense(FrVec a)
{

    uint32_t d = a.size();
    if (d == 0)
    {
        throw runtime_error("IsPolyZeroError");
    }

    uint32_t tmpd = d;
    for (uint32_t i = d - 1; i > 0; i--)
    {

        if (a[i].isZero() != true)
        {
            break;
        }
        tmpd--;
    }
    FrVec res(a.begin(), a.begin() + tmpd);
    return res;
}

FrVec PolyMul(FrVec a, FrVec b)
{

    uint32_t N = a.size() + b.size() - 1;
    N = nextPowOf2(N);
    a.resize(N, 0);
    b.resize(N, 0);

    FrVec res(N);

    Ntt<Fr> ntt;

    ntt.init(N);

    Fr *pt_a = &a[0];
    Fr *pt_b = &b[0];
    Fr *pt_res = &res[0];

    ntt.ntt(pt_a);
    ntt.ntt(pt_b);

    for (uint32_t i = 0; i < res.size(); i++)
    {
        res[i] = a[i] * b[i];
    }

    ntt.intt(pt_res);

    return PolyCondense(res);
}

polyff operator/(const polyff &A, const polyff &B)
{
    if (B.size() == 0)
    {
        throw std::runtime_error("Division by zero polynomial");
    }

    if (A.size() < B.size())
    {

        throw std::runtime_error("Divisor polynomial degree is greater than dividend");
    }

    polyff Q;
    polyff R = A;

    while (R.size() >= B.size())
    {
        size_t diff = R.size() - B.size();
        Fr scale = R.data.back() / B.data.back();
        polyff T(diff, 0);
        T.data.push_back(scale);

        Q = Q + T;
        R = R - T * B;
    }

    return Q;
}

polyff operator%(const polyff &A, const polyff &B)
{
    if (B.size() == 0)
    {
        throw std::runtime_error("Division by zero polynomial");
    }

    polyff R = A;

    while (R.size() >= B.size())
    {
        size_t diff = R.size() - B.size();
        Fr scale = R.data.back() / B.data.back();
        polyff T(diff, 0);
        T.data.push_back(scale);

        R = R - T * B;
    }

    return R;
}

polyff dividePolyByRoot(const polyff &numerator, const Fr &root)
{
    polyff quotient, remainder = numerator;
    Fr factor;

    while (remainder.size() > 1)
    {
        factor = remainder.data.back() / root;
        quotient.data.insert(quotient.data.begin(), factor);

        for (size_t i = remainder.size() - 1; i > 0; --i)
        {
            remainder[i] -= factor * root;
        }
        remainder.data.pop_back();
    }

    return quotient;
}

polyff dividePolyByRoots(const polyff &numerator, const vector<Fr> &roots)
{
    polyff result = numerator;
    for (const Fr &root : roots)
    {
        result = dividePolyByRoot(result, root);
    }
    return result;
}

std::vector<std::vector<Fr>> divideRootsIntoThreads(const std::vector<Fr> &roots, size_t threads)
{
    std::vector<std::vector<Fr>> dividedRoots(threads);
    size_t rootsPerThread = roots.size() / threads;
    size_t remainingRoots = roots.size() % threads;

    size_t start = 0;
    for (size_t i = 0; i < threads; ++i)
    {
        size_t end = start + rootsPerThread + (remainingRoots-- > 0 ? 1 : 0);
        dividedRoots[i] = std::vector<Fr>(roots.begin() + start, roots.begin() + end);
        start = end;
    }

    return dividedRoots;
}

polyff computePolynomialInThread(const std::vector<Fr> &threadRoots)
{
    polyff result = {1};
    for (const auto &root : threadRoots)
    {
        result = result * polyff({-root, 1});
    }
    return result;
}

polyff multiplyPolynomials(const std::vector<polyff> &polynomials)
{
    polyff result = {1};
    for (const auto &poly : polynomials)
    {
        result = result * poly;
    }
    return result;
}

polyff createPolynomialFromRoots(const std::vector<Fr> &roots)
{
    size_t threads = std::thread::hardware_concurrency();
    bool useParallel = roots.size() >= threads * 8;

    if (useParallel)
    {
        std::vector<std::vector<Fr>> dividedRoots = divideRootsIntoThreads(roots, threads);
        std::vector<std::thread> threadPool;
        std::vector<polyff> threadResults(threads);

        for (size_t i = 0; i < threads; ++i)
        {
            threadPool.emplace_back([&threadResults, &dividedRoots, i]()
                                    { threadResults[i] = computePolynomialInThread(dividedRoots[i]); });
        }

        for (auto &t : threadPool)
        {
            t.join();
        }

        return multiplyPolynomials(threadResults);
    }
    else
    {

        return computePolynomialInThread(roots);
    }
}

void parallelHadamardMultiply(Fr *result, const Fr *a, const Fr *b, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        result[i] = a[i] * b[i];
    }
}

polyff HadamardMultiply(const polyff &a, const std::vector<Fr> &b)
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument("Vectors must be of the same size for Hadamard multiply.");
    }

    polyff result(a.size(), 0);

    size_t N = a.size();

    if (N >= std::thread::hardware_concurrency() * 1024)
    {
        std::vector<std::thread> threads;
        size_t perThread = N / std::thread::hardware_concurrency();
        size_t start = 0;

        for (size_t i = 0; i < std::thread::hardware_concurrency(); ++i)
        {
            size_t end = std::min(start + perThread, N);
            threads.push_back(std::thread(parallelHadamardMultiply, &result.data[0], &a.data[0], &b[0], start, end));
            start = end;
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
    }
    else
    {
        for (size_t i = 0; i < N; ++i)
        {
            result[i] = a[i] * b[i];
        }
    }

    return result;
}

Fr hashToFr(const string &hash)
{
    Fr result;
    result.setHashOf(hash.c_str(), hash.size());
    return result;
}

template <typename T>
T hornerEvaluate(const std::vector<T> &coefficients, const T &x)
{
    T result = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it)
    {
        result = result * x + *it;
    }
    return result;
}

void evaluatePolynomial(Fr &result, const Fr *coefficients, size_t size, const Fr &x)
{
    result = hornerEvaluate(std::vector<Fr>(coefficients, coefficients + size), x);
}

class KZG
{
private:
    G1 g1;

    vector<G1> sPowersG1;

public:
    void init(size_t maxDegree)
    {

        Fr s;
        s.setByCSPRNG();

        G1 g1;
        hashAndMapToG1(g1, "g1");

        this->g1 = g1;

        sPowersG1.resize(maxDegree + 1);

        Fr sPow = 1;
        for (size_t i = 0; i <= maxDegree; ++i)
        {
            G1::mul(sPowersG1[i], g1, sPow);

            sPow *= s;
        }
    }

    G1 commit(polyff &polynomial)
    {
        G1 commitment;
        commitment.clear();

        Fr *pt = &polynomial.data[0];
        G1::mulVecMT(commitment, &sPowersG1[0], pt, polynomial.size(), 8);
        return commitment;
    }
};
