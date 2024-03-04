#include <iostream>
#include <vector>
#include <bitset>
#include <tuple>
#include <string>
#include <chrono>
#include <map>
#include "matrixFr.hpp"

using namespace std;

const vector<uint32_t> initialHashValues = {
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
    0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};

const uint32_t K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

string binToHex(const string &binStr)
{
    string hexStr;
    map<string, char> binToHexMap = {
        {"0000", '0'}, {"0001", '1'}, {"0010", '2'}, {"0011", '3'}, {"0100", '4'}, {"0101", '5'}, {"0110", '6'}, {"0111", '7'}, {"1000", '8'}, {"1001", '9'}, {"1010", 'a'}, {"1011", 'b'}, {"1100", 'c'}, {"1101", 'd'}, {"1110", 'e'}, {"1111", 'f'}};

    for (size_t i = 0; i < binStr.length(); i += 4)
    {

        string bin = binStr.substr(i, 4);
        hexStr += binToHexMap[bin];
    }

    return hexStr;
}

string toBinaryString(uint32_t value)
{
    if (value == 0)
        return "0";
    bitset<32> bits(value);
    string binaryString = bits.to_string();

    return binaryString.substr(binaryString.find('1'));
}

string toBinaryString0(uint32_t value)
{
    if (value == 0)
        return "0";
    bitset<32> bits(value);
    string binaryString = bits.to_string();

    return binaryString;
}

uint32_t rightRotate(uint32_t value, int bits)
{
    return (value >> bits) | (value << (32 - bits));
}

void printPaddedMessage(const vector<vector<uint32_t>> &blocks)
{
    for (const auto &block : blocks)
    {
        for (uint32_t word : block)
        {
            cout << bitset<32>(word) << endl;
        }
        cout << endl;
    }
}

vector<vector<uint32_t>> padMessage(const std::string &input)
{
    std::vector<std::vector<uint32_t>> blocks;
    size_t inputBits = input.size() * 8;
    size_t totalBits = inputBits + 1;
    size_t paddingBits = 0;

    if (totalBits % 512 > 448)
    {
        paddingBits = 512 - (totalBits % 512) + 448;
    }
    else
    {
        paddingBits = 448 - (totalBits % 512);
    }
    totalBits += paddingBits + 64;

    blocks.resize(totalBits / 512, std::vector<uint32_t>(16, 0));

    for (size_t i = 0; i < input.size(); ++i)
    {
        blocks[i / 64][(i % 64) / 4] |= static_cast<uint32_t>(input[i]) << ((3 - (i % 4)) * 8);
    }

    size_t bitPosition = inputBits;
    blocks[bitPosition / 512][(bitPosition % 512) / 32] |= 1U << (31 - ((bitPosition % 32) % 32));

    blocks.back()[14] = static_cast<uint32_t>(inputBits >> 32);
    blocks.back()[15] = static_cast<uint32_t>(inputBits & 0xFFFFFFFF);

    return blocks;
}

vector<uint32_t> messageSchedule(const vector<uint32_t> &block)
{
    vector<uint32_t> W(64);
    int i;

    for (i = 0; i < 16; ++i)
    {
        W[i] = block[i];
    }

    for (i = 16; i < 64; ++i)
    {
        uint32_t s0 = rightRotate(W[i - 15], 7) ^ rightRotate(W[i - 15], 18) ^ (W[i - 15] >> 3);
        uint32_t s1 = rightRotate(W[i - 2], 17) ^ rightRotate(W[i - 2], 19) ^ (W[i - 2] >> 10);
        W[i] = W[i - 16] + s0 + W[i - 7] + s1;
    }

    return W;
}

tuple<vector<vector<uint32_t>>, vector<uint32_t>, string>
compressAndDigest(const vector<uint32_t> &W, bool isFirstBlock, const vector<uint32_t> &prevHashValues)
{
    vector<uint32_t> H = isFirstBlock ? initialHashValues : prevHashValues;
    vector<vector<uint32_t>> workingVars;

    for (int i = 0; i < 64; ++i)
    {
        uint32_t S1 = rightRotate(H[4], 6) ^ rightRotate(H[4], 11) ^ rightRotate(H[4], 25);
        uint32_t ch = (H[4] & H[5]) ^ (~H[4] & H[6]);
        uint32_t temp1 = H[7] + S1 + ch + K[i] + W[i];
        uint32_t S0 = rightRotate(H[0], 2) ^ rightRotate(H[0], 13) ^ rightRotate(H[0], 22);
        uint32_t maj = (H[0] & H[1]) ^ (H[0] & H[2]) ^ (H[1] & H[2]);
        uint32_t temp2 = S0 + maj;

        H[7] = H[6];
        H[6] = H[5];
        H[5] = H[4];
        H[4] = H[3] + temp1;
        H[3] = H[2];
        H[2] = H[1];
        H[1] = H[0];
        H[0] = temp1 + temp2;

        workingVars.push_back(H);
    }

    vector<uint32_t> finalHashValues = isFirstBlock ? initialHashValues : prevHashValues;
    for (int i = 0; i < 8; ++i)
    {
        finalHashValues[i] += H[i];
    }

    string finalHashBinaryString;
    for (uint32_t value : finalHashValues)
    {
        bitset<32> bits(value);
        finalHashBinaryString += bits.to_string();
    }

    return make_tuple(workingVars, finalHashValues, finalHashBinaryString);
}

tuple<vector<vector<uint32_t>>, vector<vector<vector<uint32_t>>>, string>
calculateSHA256(const string &input)
{

    auto blocks = padMessage(input);

    vector<vector<uint32_t>> allMessageSchedules;    
    vector<vector<vector<uint32_t>>> allWorkingVars; 
    vector<uint32_t> H = initialHashValues;          

    for (size_t i = 0; i < blocks.size(); ++i)
    {

        auto W = messageSchedule(blocks[i]);
        allMessageSchedules.push_back(W);

        bool isFirstBlock = (i == 0);
        auto [workingVars, finalHashValues, _] = compressAndDigest(W, isFirstBlock, H);

        allWorkingVars.push_back(workingVars);
        H = finalHashValues; 
    }

    string finalHashHex;
    for (uint32_t value : H)
    {
        finalHashHex += binToHex(toBinaryString0(value));
    }

    return make_tuple(allMessageSchedules, allWorkingVars, finalHashHex);
}
