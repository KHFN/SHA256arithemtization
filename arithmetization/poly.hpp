#ifndef POLY_HPP
#define POLY_HPP

#include "poly.cpp"
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


inline void MCL_init() {
    initPairing(mcl::BLS12_381);
}


template<class T>
void initPowSeq(T* sec, const T& u, size_t n);

bool IsPowerOfTwo(uint64_t v);

uint32_t Max(uint32_t x, uint32_t y);

uint32_t Min(uint32_t x, uint32_t y);

uint64_t nextPowOf2(const uint64_t v);

void parallelComposeWithMonomial(Fr* resultData, const Fr* inputData, const Fr* powSeq, size_t start, size_t end);

class polyff; //{
// public:
//     vector<Fr> data;

//     polyff();
//     polyff(std::initializer_list<Fr> initList);
//     polyff(int val);
//     polyff(Fr val);
//     polyff(vector<Fr>::const_iterator begin, vector<Fr>::const_iterator end);
//     polyff(size_t size, int initialValue);
//     polyff(size_t size, Fr initialValue);

//     size_t size() const;
//     void resize(size_t newSize, const Fr& value = Fr(0));
//     Fr& operator[](size_t index);
//     const Fr& operator[](size_t index) const;
//     Fr* getDataPtr();
//     const Fr* getDataPtr() const;
//     polyff PolyCondense() const;
//     bool operator==(const polyff& rhs) const;
//     bool operator!() const;
//     polyff divideByXnMinusOne(int n) const;
//     Fr evaluate(const Fr& x) const;

//     friend polyff operator+(const polyff& a, const polyff& b);
//     friend polyff operator-(const polyff& a, const polyff& b);
//     friend polyff operator*(const polyff& a, const polyff& b);
//     friend polyff operator/(const polyff& A, const polyff& B);
//     friend polyff operator%(const polyff& A, const polyff& B);
//     friend std::ostream& operator<<(std::ostream& os, const polyff& p);
// };

polyff operator+(const polyff& a, const polyff& b);
polyff operator-(const polyff& a, const polyff& b);
polyff operator*(const polyff& a, const polyff& b);
polyff operator/(const polyff& A, const polyff& B);
polyff operator%(const polyff& A, const polyff& B);
std::ostream& operator<<(std::ostream& os, const polyff& p);


void parallelAdd(polyff& result, const polyff& a, const polyff& b, size_t start, size_t end);
void parallelSub(polyff& result, const polyff& a, const polyff& b, size_t start, size_t end);
void parallelMultiply(Fr* result, const Fr* a, const Fr* b, size_t start, size_t end);
polyff HadamardMultiply(const polyff& a, const std::vector<Fr>& b);
Fr hashToFr(const string& hash);
template<typename T>
T hornerEvaluate(const std::vector<T>& coefficients, const T& x);
void evaluatePolynomial(Fr& result, const Fr* coefficients, size_t size, const Fr& x);


class KZG;

#endif
