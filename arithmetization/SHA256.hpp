#ifndef SHA256_HPP
#define SHA256_HPP

#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <map>
#include "matrixFr.hpp"
#include "SHA256.cpp"

extern const std::vector<uint32_t> initialHashValues;

extern const uint32_t K[64];

std::string binToHex(const std::string &binStr);

std::string toBinaryString(uint32_t value);

std::string toBinaryString0(uint32_t value);

uint32_t rightRotate(uint32_t value, int bits);

void printPaddedMessage(const std::vector<std::vector<uint32_t>> &blocks);

std::vector<std::vector<uint32_t>> padMessage(const std::string &input);

std::vector<uint32_t> messageSchedule(const std::vector<uint32_t> &block);

std::tuple<std::vector<std::vector<uint32_t>>, std::vector<uint32_t>, std::string>
compressAndDigest(const std::vector<uint32_t> &W, bool isFirstBlock, const std::vector<uint32_t> &prevHashValues);

std::tuple<std::vector<std::vector<uint32_t>>, std::vector<std::vector<std::vector<uint32_t>>>, std::string>
calculateSHA256(const std::string &input);

#endif