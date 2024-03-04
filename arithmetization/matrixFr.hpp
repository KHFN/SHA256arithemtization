// matrixFr.h
#ifndef MATRIXFR_HPP
#define MATRIXFR_HPP

#include "mcl/bls12_381.hpp"
#include "matrixFr.cpp"
#include <vector>
#include <iostream>
#include <fstream>

class matrixFr; 
// {
// public:
//     std::vector<std::vector<mcl::bn::Fr>> data;

//     matrixFr();
//     matrixFr(size_t outerSize, size_t innerSize);

//     std::vector<mcl::bn::Fr>& operator[](size_t index);
//     const std::vector<mcl::bn::Fr>& operator[](size_t index) const;

//     void insertColumn(size_t position, const std::vector<mcl::bn::Fr>& column);
//     std::vector<mcl::bn::Fr> getColumn(size_t position) const;
//     void replaceColumn(size_t position, const std::vector<mcl::bn::Fr>& newColumn);
//     void replaceRow(size_t position, const std::vector<mcl::bn::Fr>& newRow);

//     friend std::ostream& operator<<(std::ostream& os, const matrixFr& m);
// };
template<typename T>
void printVector(const vector<T>& vec);
class matrixFr;
void multiplyPartByConstant(vector<Fr>& vec, const Fr& constant, size_t start, size_t end);
void mulveccons(vector<Fr>& vec, const Fr& constant);
void addVectorSegment(const vector<Fr>& a, const vector<Fr>& b, vector<Fr>& result, size_t start, size_t end);
vector<Fr> addVectors(const vector<Fr>& a, const vector<Fr>& b);
void saveToCSV(const vector<vector<Fr>>& data, const string& filename);
void loadFromCSV(vector<vector<Fr>>& data, const string& filename);
template <typename T>
vector<T> concatenateVectors(const vector<T>& vec1, const vector<T>& vec2);
template<typename T>
T dotProduct(const vector<T>& v1, const vector<T>& v2);
void saveToCSV(const vector<vector<Fr>>& data, const string& filename);
void loadFromCSV(vector<vector<Fr>>& data, const string& filename);

#endif // MATRIXFR_H