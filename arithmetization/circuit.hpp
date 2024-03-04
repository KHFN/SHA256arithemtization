#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "prover.hpp"
#include "poly.hpp"
#include "matrixFr.hpp"
#include "prover.hpp"
#include "circuit.cpp"

class circuit;

Fr xor3Fr(Fr a, Fr b, Fr c);

Fr choiceFr(Fr a, Fr b, Fr c);

Fr majorFr(Fr a, Fr b, Fr c);

tuple<vector<Fr>, vector<Fr>>
splitBinaryString(unsigned int num, vector<int> lengths);

tuple<vector<Fr>,vector<Fr>>
Gate_decomposition(uint32_t chunck, vector<int> lengths);

tuple<vector<Fr>,vector<Fr>>
split_ss0_rotr7(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_ss0_rotr18(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_ss0_shr3(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_ss1_rotr17(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_ss1_rotr19(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_ss1_shr10(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs0_rotr6(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs0_rotr11(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs0_rotr25(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs1_rotr2(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs1_rotr13(uint32_t chunck);

tuple<vector<Fr>,vector<Fr>>
split_bs1_rotr22(uint32_t chunck);

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>> 
bigsigma_1(uint32_t numberofselector, uint32_t numberofwire, uint32_t chunk);

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>> 
bigsigma_0(uint32_t numberofselector, uint32_t numberofwire, uint32_t chunk);

tuple<vector<Fr>, vector<Fr>> 
bitdecomposition_4(uint32_t chunk);

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>> 
majorgate(uint32_t numberofselector, uint32_t numberofwire, uint32_t A, uint32_t B, uint32_t C);

tuple<vector<vector<Fr>>, vector<vector<vector<uint32_t>>>> 
choicegate(uint32_t numberofselector, uint32_t numberofwire, uint32_t A, uint32_t B, uint32_t C);

vector<vector<Fr>> 
add5gate(uint32_t numberofselector, uint32_t numberofwire, vector<Fr> wire);

#endif
