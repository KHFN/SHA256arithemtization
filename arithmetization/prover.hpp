
#ifndef PROVER_HPP
#define PROVER_HPP

#include "prover.cpp"
#include "poly.hpp"
#include "matrixFr.hpp"
#include <utility>


typedef vector<vector<vector<uint32_t>>> vvv_int;

template<typename T>
void inttAndCreatePolyff(std::vector<Fr>& vec, polyff& result);


std::tuple<vector<polyff>, vector<vector<Fr>>, matrixFr> 
permutation(const vector<Fr>& W, size_t row, size_t col, const vvv_int& positions);
string vecToStr(const vector<Fr>& vec);
tuple<vector<Fr>, vector<Fr>, vector<Fr>> 
lookup(const matrixFr& trace, const matrixFr& table, uint32_t selector, uint32_t numberofselector, uint32_t numberofwire, Fr& alpha);
tuple<polyff, polyff, polyff, polyff, polyff, polyff, polyff, polyff, polyff>
plookup_copyconstraints(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma);
tuple<polyff, polyff, polyff> 
plookup_poly_ft(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma);
tuple<polyff, polyff, polyff, polyff, polyff, polyff, polyff> 
plookup_copyconstraints_1(vector<Fr> f, vector<Fr> t);
tuple<polyff, polyff> 
Plookup_poly_z(vector<Fr> f, vector<Fr> t, Fr beta, Fr gamma);



#endif 