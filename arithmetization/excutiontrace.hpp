#ifndef EXCUTIONTRACE_HPP
#define EXCUTIONTRACE_HPP

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <unordered_set>
#include <utility>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <tuple>
#include <bitset>
#include <cstdint>
#include <stdexcept>
#include "matrixFr.hpp"
#include "poly.hpp"
#include "prover.hpp"
#include "circuit.hpp"
#include "SHA256.hpp"
#include "excutiontrace.cpp"

using namespace std;


bool compareHexStrings(string a, string b);
string generateRandomString(size_t length);
tuple<circuit, double> generateExcutiontrace(uint32_t bitSize);

#endif 