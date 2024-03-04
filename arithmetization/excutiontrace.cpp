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

using namespace std;

bool compareHexStrings(string a, string b)
{

    transform(a.begin(), a.end(), a.begin(), ::tolower);
    transform(b.begin(), b.end(), b.begin(), ::tolower);

    return a == b;
}

string generateRandomString(size_t length)
{
    static const char characters[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    string result;
    result.reserve(length);

    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distribution(0, sizeof(characters) - 2);

    for (size_t i = 0; i < length; ++i)
    {
        result += characters[distribution(generator)];
    }

    return result;
}

tuple<circuit, double> generateExcutiontrace(uint32_t bitSize)
{

    uint32_t row;
    uint32_t numberofselector = 3;
    uint32_t numberofwire = 13;

    size_t length = bitSize / 8;
    string input = generateRandomString(length);

    auto [messageSchedules, workingVars, finalHashHex] = calculateSHA256(input);

    row = 6856 * (messageSchedules.size() - 1) + 1728 + 80 * 64 + 8;
    circuit A(row, numberofselector, numberofwire);

    cout << "해시 입력 블록 수: " << messageSchedules.size() << endl;
    cout << "Gate 수: " << row << endl;

    auto start1 = chrono::high_resolution_clock::now();
    vector<vector<uint32_t>> pos_permutation;
    vector<uint32_t> initialworkVar;

    for (uint32_t i = 0; i < messageSchedules.size(); i++)
    {

        for (uint32_t j = 6856 * i + 0; j < (6856 * i + 36 * 48); j += 36)
        {

            auto [a, b] = split_ss0_rotr7(messageSchedules[i][(j - 6856 * i) / 36 + 1]);

            vector<vector<Fr>> tmp;
            tmp.push_back(a);
            tmp.push_back(b);
            A.excutiontrace.replaceSubmatrix({j, numberofselector}, tmp);
            A.excutiontrace[j][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 1][numberofselector - numberofselector] = 1;

            auto [c, d] = split_ss0_rotr18(messageSchedules[i][(j - 6856 * i) / 36 + 1]);

            tmp.clear();
            tmp.push_back(c);
            tmp.push_back(d);
            A.excutiontrace.replaceSubmatrix({j + 2, numberofselector}, tmp);
            A.excutiontrace[j + 2][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 3][numberofselector - numberofselector] = 1;

            auto [e, f] = split_ss0_shr3(messageSchedules[i][(j - 6856 * i) / 36 + 1]);

            tmp.clear();
            tmp.push_back(e);
            tmp.push_back(f);
            A.excutiontrace.replaceSubmatrix({j + 4, numberofselector}, tmp);
            A.excutiontrace[j + 4][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 5][numberofselector - numberofselector] = 1;

            A.excutiontrace[j + 6][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 7][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 8][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 9][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 10][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 11][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 12][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 13][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 14][numberofselector - numberofselector + 2] = 1;

            A.routing({{j, numberofselector + 2}, {j + 6, numberofselector}});
            A.routing({{j, numberofselector + 3}, {j + 7, numberofselector}});
            A.routing({{j, numberofselector + 4}, {j + 8, numberofselector}});
            A.routing({{j, numberofselector + 5}, {j + 9, numberofselector}});
            A.routing({{j + 1, numberofselector + 1}, {j + 10, numberofselector}});
            A.routing({{j + 1, numberofselector + 2}, {j + 11, numberofselector}});
            A.routing({{j + 1, numberofselector + 3}, {j + 12, numberofselector}});
            A.routing({{j, numberofselector}, {j + 13, numberofselector}});
            A.routing({{j, numberofselector + 1}, {j + 14, numberofselector}});

            A.routing({{j + 2, numberofselector + 5}, {j + 6, numberofselector + 1}});
            A.routing({{j + 3, numberofselector + 1}, {j + 7, numberofselector + 1}});
            A.routing({{j + 3, numberofselector + 2}, {j + 8, numberofselector + 1}});
            A.routing({{j + 3, numberofselector + 3}, {j + 9, numberofselector + 1}});
            A.routing({{j + 2, numberofselector}, {j + 10, numberofselector + 1}});
            A.routing({{j + 2, numberofselector + 1}, {j + 11, numberofselector + 1}});
            A.routing({{j + 2, numberofselector + 2}, {j + 12, numberofselector + 1}});
            A.routing({{j + 2, numberofselector + 3}, {j + 13, numberofselector + 1}});
            A.routing({{j + 2, numberofselector + 4}, {j + 14, numberofselector + 1}});

            A.routing({{j + 4, numberofselector + 1}, {j + 6, numberofselector + 2}});
            A.routing({{j + 4, numberofselector + 2}, {j + 7, numberofselector + 2}});
            A.routing({{j + 4, numberofselector + 3}, {j + 8, numberofselector + 2}});
            A.routing({{j + 4, numberofselector + 4}, {j + 9, numberofselector + 2}});
            A.routing({{j + 4, numberofselector + 5}, {j + 10, numberofselector + 2}});
            A.routing({{j + 5, numberofselector + 1}, {j + 11, numberofselector + 2}});
            A.routing({{j + 5, numberofselector + 2}, {j + 12, numberofselector + 2}});
            A.routing({{j + 5, numberofselector + 3}, {j + 13, numberofselector + 2}});
            A.excutiontrace[j + 14][numberofselector + 2] = 0;

            A.excutiontrace[j + 6][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 6][numberofselector], A.excutiontrace[j + 6][numberofselector + 1], A.excutiontrace[j + 6][numberofselector + 2]);
            A.excutiontrace[j + 7][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 7][numberofselector], A.excutiontrace[j + 7][numberofselector + 1], A.excutiontrace[j + 7][numberofselector + 2]);
            A.excutiontrace[j + 8][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 8][numberofselector], A.excutiontrace[j + 8][numberofselector + 1], A.excutiontrace[j + 8][numberofselector + 2]);
            A.excutiontrace[j + 9][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 9][numberofselector], A.excutiontrace[j + 9][numberofselector + 1], A.excutiontrace[j + 9][numberofselector + 2]);
            A.excutiontrace[j + 10][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 10][numberofselector], A.excutiontrace[j + 10][numberofselector + 1], A.excutiontrace[j + 10][numberofselector + 2]);
            A.excutiontrace[j + 11][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 11][numberofselector], A.excutiontrace[j + 11][numberofselector + 1], A.excutiontrace[j + 11][numberofselector + 2]);
            A.excutiontrace[j + 12][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 12][numberofselector], A.excutiontrace[j + 12][numberofselector + 1], A.excutiontrace[j + 12][numberofselector + 2]);
            A.excutiontrace[j + 13][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 13][numberofselector], A.excutiontrace[j + 13][numberofselector + 1], A.excutiontrace[j + 13][numberofselector + 2]);
            A.excutiontrace[j + 14][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 14][numberofselector], A.excutiontrace[j + 14][numberofselector + 1], A.excutiontrace[j + 14][numberofselector + 2]);

            A.excutiontrace[j + 15][0] = 1;
            A.excutiontrace[j + 16][0] = 1;
            A.routing({{j + 6, numberofselector + 3}, {j + 15, numberofselector + 0}});
            A.routing({{j + 7, numberofselector + 3}, {j + 15, numberofselector + 1}});
            A.routing({{j + 8, numberofselector + 3}, {j + 15, numberofselector + 2}});
            A.routing({{j + 9, numberofselector + 3}, {j + 15, numberofselector + 3}});
            A.routing({{j + 10, numberofselector + 3}, {j + 15, numberofselector + 4}});
            A.routing({{j + 11, numberofselector + 3}, {j + 15, numberofselector + 5}});
            A.routing({{j + 12, numberofselector + 3}, {j + 16, numberofselector + 1}});
            A.routing({{j + 13, numberofselector + 3}, {j + 16, numberofselector + 2}});
            A.routing({{j + 14, numberofselector + 3}, {j + 16, numberofselector + 3}});

            A.excutiontrace.replaceSubmatrix({j + 15, numberofselector + 7},
                                             {{Fr(1 << 0), Fr(1 << 2), Fr(1 << 6), Fr(1 << 10), Fr(1 << 14), Fr(1 << 17)}});
            A.excutiontrace.replaceSubmatrix({j + 16, numberofselector + 7},
                                             {{Fr(1 << 0), Fr(1 << 21), Fr(1 << 25), Fr(1 << 29), 0, 0}});

            vector<Fr> e1(A.excutiontrace[j + 15].begin() + numberofselector, A.excutiontrace[j + 15].begin() + numberofselector + 6);
            vector<Fr> e2(A.excutiontrace[j + 15].begin() + numberofselector + 7, A.excutiontrace[j + 15].begin() + numberofselector + 7 + 6);
            A.excutiontrace[j + 15][numberofselector + 6] = dotProduct(e1, e2);
            A.routing({{j + 15, numberofselector + 6}, {j + 16, numberofselector}});

            vector<Fr> e3(A.excutiontrace[j + 16].begin() + numberofselector, A.excutiontrace[j + 16].begin() + numberofselector + 6);
            vector<Fr> e4(A.excutiontrace[j + 16].begin() + numberofselector + 7, A.excutiontrace[j + 16].begin() + numberofselector + 7 + 6);

            A.excutiontrace[j + 16][numberofselector + 6] = dotProduct(e3, e4);
            A.routing({{j + 16, numberofselector + 6}, {j + 35, numberofselector + 1}});
            A.excutiontrace[j + 35][1] = 1;

            auto [gg, ff] = split_ss1_rotr17(messageSchedules[i][(j - 6856 * i) / 36 + 14]);

            tmp.clear();
            tmp.push_back(gg);
            tmp.push_back(ff);
            A.excutiontrace.replaceSubmatrix({j + 17, numberofselector}, tmp);
            A.excutiontrace[j + 17][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 17 + 1][numberofselector - numberofselector] = 1;

            auto [hh, ii] = split_ss1_rotr19(messageSchedules[i][(j - 6856 * i) / 36 + 14]);

            tmp.clear();
            tmp.push_back(hh);
            tmp.push_back(ii);
            A.excutiontrace.replaceSubmatrix({j + 17 + 2, numberofselector}, tmp);
            A.excutiontrace[j + 17 + 2][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 17 + 3][numberofselector - numberofselector] = 1;

            auto [jj, kk] = split_ss1_shr10(messageSchedules[i][(j - 6856 * i) / 36 + 14]);

            tmp.clear();
            tmp.push_back(jj);
            tmp.push_back(kk);
            A.excutiontrace.replaceSubmatrix({j + 17 + 4, numberofselector}, tmp);
            A.excutiontrace[j + 17 + 4][numberofselector - numberofselector] = 1;
            A.excutiontrace[j + 17 + 5][numberofselector - numberofselector] = 1;

            A.excutiontrace[j + 17 + 6][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 7][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 8][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 9][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 10][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 11][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 12][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 13][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 14][numberofselector - numberofselector + 2] = 1;
            A.excutiontrace[j + 17 + 15][numberofselector - numberofselector + 2] = 1;

            A.routing({{j + 17, numberofselector + 5}, {j + 17 + 6, numberofselector}});
            A.routing({{j + 17 + 1, numberofselector + 1}, {j + 17 + 7, numberofselector}});
            A.routing({{j + 17 + 1, numberofselector + 2}, {j + 17 + 8, numberofselector}});
            A.routing({{j + 17 + 1, numberofselector + 3}, {j + 17 + 9, numberofselector}});
            A.routing({{j + 17 + 1, numberofselector + 4}, {j + 17 + 10, numberofselector}});
            A.routing({{j + 17, numberofselector}, {j + 17 + 11, numberofselector}});
            A.routing({{j + 17, numberofselector + 1}, {j + 17 + 12, numberofselector}});
            A.routing({{j + 17, numberofselector + 2}, {j + 17 + 13, numberofselector}});
            A.routing({{j + 17, numberofselector + 3}, {j + 17 + 14, numberofselector}});
            A.routing({{j + 17, numberofselector + 4}, {j + 17 + 15, numberofselector}});

            A.routing({{j + 17 + 3, numberofselector + 1}, {j + 17 + 6, numberofselector + 1}});
            A.routing({{j + 17 + 3, numberofselector + 2}, {j + 17 + 7, numberofselector + 1}});
            A.routing({{j + 17 + 3, numberofselector + 3}, {j + 17 + 8, numberofselector + 1}});
            A.routing({{j + 17 + 3, numberofselector + 4}, {j + 17 + 9, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector}, {j + 17 + 10, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector + 1}, {j + 17 + 11, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector + 2}, {j + 17 + 12, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector + 3}, {j + 17 + 13, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector + 4}, {j + 17 + 14, numberofselector + 1}});
            A.routing({{j + 17 + 2, numberofselector + 5}, {j + 17 + 15, numberofselector + 1}});

            A.routing({{j + 17 + 4, numberofselector + 3}, {j + 17 + 6, numberofselector + 2}});
            A.routing({{j + 17 + 4, numberofselector + 4}, {j + 17 + 7, numberofselector + 2}});
            A.routing({{j + 17 + 4, numberofselector + 5}, {j + 17 + 8, numberofselector + 2}});
            A.routing({{j + 17 + 5, numberofselector + 1}, {j + 17 + 9, numberofselector + 2}});
            A.routing({{j + 17 + 5, numberofselector + 2}, {j + 17 + 10, numberofselector + 2}});
            A.routing({{j + 17 + 5, numberofselector + 3}, {j + 17 + 11, numberofselector + 2}});
            A.routing({{j + 17 + 5, numberofselector + 4}, {j + 17 + 12, numberofselector + 2}});
            A.routing({{j + 17 + 4, numberofselector}, {j + 17 + 13, numberofselector + 2}});
            A.routing({{j + 17 + 4, numberofselector + 1}, {j + 17 + 14, numberofselector + 2}});
            A.routing({{j + 17 + 4, numberofselector + 2}, {j + 17 + 15, numberofselector + 2}});
            A.excutiontrace[j + 30][numberofselector + 2] = 0;
            A.excutiontrace[j + 31][numberofselector + 2] = 0;
            A.excutiontrace[j + 32][numberofselector + 2] = 0;

            A.excutiontrace[j + 17 + 6][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 6][numberofselector], A.excutiontrace[j + 17 + 6][numberofselector + 1], A.excutiontrace[j + 17 + 6][numberofselector + 2]);
            A.excutiontrace[j + 17 + 7][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 7][numberofselector], A.excutiontrace[j + 17 + 7][numberofselector + 1], A.excutiontrace[j + 17 + 7][numberofselector + 2]);
            A.excutiontrace[j + 17 + 8][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 8][numberofselector], A.excutiontrace[j + 17 + 8][numberofselector + 1], A.excutiontrace[j + 17 + 8][numberofselector + 2]);
            A.excutiontrace[j + 17 + 9][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 9][numberofselector], A.excutiontrace[j + 17 + 9][numberofselector + 1], A.excutiontrace[j + 17 + 9][numberofselector + 2]);
            A.excutiontrace[j + 17 + 10][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 10][numberofselector], A.excutiontrace[j + 17 + 10][numberofselector + 1], A.excutiontrace[j + 17 + 10][numberofselector + 2]);
            A.excutiontrace[j + 17 + 11][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 11][numberofselector], A.excutiontrace[j + 17 + 11][numberofselector + 1], A.excutiontrace[j + 17 + 11][numberofselector + 2]);
            A.excutiontrace[j + 17 + 12][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 12][numberofselector], A.excutiontrace[j + 17 + 12][numberofselector + 1], A.excutiontrace[j + 17 + 12][numberofselector + 2]);
            A.excutiontrace[j + 17 + 13][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 13][numberofselector], A.excutiontrace[j + 17 + 13][numberofselector + 1], A.excutiontrace[j + 17 + 13][numberofselector + 2]);
            A.excutiontrace[j + 17 + 14][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 14][numberofselector], A.excutiontrace[j + 17 + 14][numberofselector + 1], A.excutiontrace[j + 17 + 14][numberofselector + 2]);
            A.excutiontrace[j + 17 + 15][numberofselector + 3] = xor3Fr(A.excutiontrace[j + 17 + 15][numberofselector], A.excutiontrace[j + 17 + 15][numberofselector + 1], A.excutiontrace[j + 17 + 15][numberofselector + 2]);

            A.excutiontrace[j + 17 + 16][0] = 1;
            A.excutiontrace[j + 17 + 17][0] = 1;
            A.routing({{j + 17 + 6, numberofselector + 3}, {j + 17 + 16, numberofselector + 0}});
            A.routing({{j + 17 + 7, numberofselector + 3}, {j + 17 + 16, numberofselector + 1}});
            A.routing({{j + 17 + 8, numberofselector + 3}, {j + 17 + 16, numberofselector + 2}});
            A.routing({{j + 17 + 9, numberofselector + 3}, {j + 17 + 16, numberofselector + 3}});
            A.routing({{j + 17 + 10, numberofselector + 3}, {j + 17 + 16, numberofselector + 4}});
            A.routing({{j + 17 + 11, numberofselector + 3}, {j + 17 + 16, numberofselector + 5}});
            A.routing({{j + 17 + 12, numberofselector + 3}, {j + 17 + 17, numberofselector + 1}});
            A.routing({{j + 17 + 13, numberofselector + 3}, {j + 17 + 17, numberofselector + 2}});
            A.routing({{j + 17 + 14, numberofselector + 3}, {j + 17 + 17, numberofselector + 3}});
            A.routing({{j + 17 + 15, numberofselector + 3}, {j + 17 + 17, numberofselector + 4}});

            A.excutiontrace.replaceSubmatrix({j + 17 + 16, numberofselector + 7},
                                             {{Fr(1 << 0), Fr(1 << 2), Fr(1 << 5), Fr(1 << 9), Fr(1 << 13), Fr(1 << 15)}});
            A.excutiontrace.replaceSubmatrix({j + 17 + 17, numberofselector + 7},
                                             {{Fr(1 << 0), Fr(1 << 18), Fr(1 << 22), Fr(1 << 24), Fr(1 << 28), 0}});

            vector<Fr> e5(A.excutiontrace[j + 17 + 16].begin() + numberofselector, A.excutiontrace[j + 17 + 16].begin() + numberofselector + 6);
            vector<Fr> e6(A.excutiontrace[j + 17 + 16].begin() + numberofselector + 7, A.excutiontrace[j + 17 + 16].begin() + numberofselector + 7 + 6);
            A.excutiontrace[j + 17 + 16][numberofselector + 6] = dotProduct(e5, e6);
            A.routing({{j + 17 + 16, numberofselector + 6}, {j + 17 + 17, numberofselector}});

            vector<Fr> e7(A.excutiontrace[j + 17 + 17].begin() + numberofselector, A.excutiontrace[j + 17 + 17].begin() + numberofselector + 6);
            vector<Fr> e8(A.excutiontrace[j + 17 + 17].begin() + numberofselector + 7, A.excutiontrace[j + 17 + 17].begin() + numberofselector + 7 + 6);

            A.excutiontrace[j + 17 + 17][numberofselector + 6] = dotProduct(e7, e8);
            A.routing({{j + 17 + 17, numberofselector + 6}, {j + 35, numberofselector + 3}});
            A.excutiontrace[j + 35][1] = 1;

            A.excutiontrace[j + 35][numberofselector] = messageSchedules[i][(j - 6856 * i) / 36];
            A.excutiontrace[j + 35][numberofselector + 2] = messageSchedules[i][(j - 6856 * i) / 36 + 9];

            uint64_t W_new = (A.excutiontrace[j + 35][numberofselector] + A.excutiontrace[j + 35][numberofselector + 1] + A.excutiontrace[j + 35][numberofselector + 2] + A.excutiontrace[j + 35][numberofselector + 3] + A.excutiontrace[j + 35][numberofselector + 4] + A.excutiontrace[j + 35][numberofselector + 5]).getInt64();
            uint64_t W_new_copy = W_new;

            W_new = W_new & (4294967295);
            A.excutiontrace[j + 35][numberofselector + 6] = Fr(W_new);
            A.excutiontrace[j + 35][numberofselector + 5] = A.excutiontrace[j + 35][numberofselector + 6] - Fr(W_new_copy);
        }

        if (i == 0)
        {

            initialworkVar = initialHashValues;
        }

        const uint32_t startingrow = 1728;
        uint32_t presentrow;

        for (uint32_t k = 6856 * i; k < 6856 * i + 80 * 64; k += 80)
        {

            vector<uint32_t> workVar;
            uint32_t m = (k - 6856 * i) / 80;

            if (m == 0)
            {

                workVar = initialworkVar;
            }
            else
            {

                workVar = workingVars[i][m - 1];
            }

            presentrow = startingrow + k;
            auto [region1, b] = bigsigma_1(numberofselector, numberofwire, workVar[4]);
            for (int i = 0; i < b.size(); i++)
            {

                b[i][0][0] += presentrow;
                b[i][1][0] += presentrow;
            }
            A.permutation_pos_interior = concatenateVectors(A.permutation_pos_interior, b);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region1);

            presentrow = presentrow + region1.size();
            auto [region2, b2] = choicegate(numberofselector, numberofwire, workVar[4], workVar[5], workVar[6]);
            for (int i = 0; i < b2.size(); i++)
            {

                b2[i][0][0] += presentrow;
                b2[i][1][0] += presentrow;
            }
            A.permutation_pos_interior = concatenateVectors(A.permutation_pos_interior, b2);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region2);

            presentrow = presentrow + region2.size();
            Fr Wi;
            if (m < 16)
            {
                Wi = Fr(messageSchedules[i][m]);
                A.excutiontrace[presentrow][7] = Wi;
            }
            else
            {
                A.routing({{(6856 * i + (36 * (m - 15) - 1)), 9}, {presentrow, 7}});
                Wi = A.excutiontrace[(6856 * i + (36 * (m - 15) - 1))][9];
            }

            Fr h = workVar.back();
            A.excutiontrace[presentrow][numberofselector] = h;

            uint32_t temprow = (1728 + k + (region1.size() - 1));
            A.routing({{temprow, 9}, {presentrow, numberofselector + 1}});
            Fr bsig1 = A.excutiontrace[presentrow][numberofselector + 1];

            temprow = (1728 + k + (region1.size()) + (region2.size()) - 1);
            A.routing({{temprow, 9}, {presentrow, numberofselector + 2}});
            Fr ch = A.excutiontrace[presentrow][numberofselector + 2];

            Fr Ki = K[m];
            A.excutiontrace[presentrow][numberofselector + 3] = Ki;
            vector<Fr> wire(A.excutiontrace[presentrow].begin() + numberofselector,
                            A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region3 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, {region3});

            presentrow = presentrow + region3.size();
            auto [region4, b4] = bigsigma_0(numberofselector, numberofwire, workVar[0]);
            for (int i = 0; i < b4.size(); i++)
            {

                b4[i][0][0] += presentrow;
                b4[i][1][0] += presentrow;
            }
            A.permutation_pos_interior = concatenateVectors(A.permutation_pos_interior, b4);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region4);

            presentrow = presentrow + region4.size();
            auto [region5, b5] = majorgate(numberofselector, numberofwire,
                                           workVar[0], workVar[1], workVar[2]);
            for (int i = 0; i < b5.size(); i++)
            {

                b5[i][0][0] += presentrow;
                b5[i][1][0] += presentrow;
            }
            A.permutation_pos_interior = concatenateVectors(A.permutation_pos_interior, b5);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region5);

            presentrow = presentrow + region5.size();
            temprow = 1728 + k + region1.size() + region2.size() + region3.size() + region4.size() - 1;
            A.routing({{temprow, 9}, {presentrow, numberofselector + 1}});
            temprow += region5.size();
            A.routing({{temprow, 9}, {presentrow, numberofselector}});
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region6 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region6);

            presentrow = presentrow + region6.size();
            temprow = 1728 + k + region1.size() + region2.size() + region3.size() - 1;
            A.routing({{temprow, 9}, {presentrow, numberofselector}});
            temprow += (region4.size() + region5.size() + region6.size());
            A.routing({{temprow, 9}, {presentrow, numberofselector + 1}});
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region7 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region7);

            presentrow = presentrow + region7.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[3];
            temprow = 1728 + k + region1.size() + region2.size() + region3.size() - 1;
            A.routing({{temprow, 9}, {presentrow, numberofselector + 1}});
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region8 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region8);

            presentrow = presentrow + region8.size();
            temprow = presentrow - region8.size() - 1;
            A.routing({{temprow, 9}, {presentrow, numberofselector}});
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region9 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region9);

            presentrow = presentrow + region9.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[0];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region10 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region10);

            presentrow = presentrow + region10.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[1];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region11 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region11);

            presentrow = presentrow + region11.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[2];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region12 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region12);

            presentrow = presentrow + region12.size();
            temprow = presentrow - (region9.size() + region10.size() + region11.size() + region12.size()) - 1;
            A.routing({{temprow, 9}, {presentrow, numberofselector}});
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region13 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region13);

            presentrow = presentrow + region13.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[4];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region14 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region14);

            presentrow = presentrow + region14.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[5];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region15 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region15);

            presentrow = presentrow + region15.size();
            A.excutiontrace[presentrow][numberofselector] = workVar[6];
            wire.assign(A.excutiontrace[presentrow].begin() + numberofselector,
                        A.excutiontrace[presentrow].begin() + numberofselector + 5);
            auto region16 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({presentrow, 0}, region16);
        }

        uint32_t temp;
        vector<Fr> wire;
        for (uint32_t t = 6856 * (i + 1) - 8; t < 6856 * (i + 1); t++)
        {

            A.routing({{(t - 8), 9}, {t, numberofselector}});
            A.excutiontrace[t][numberofselector + 1] = initialworkVar[t - (6856 * i + 6848)];
            wire.assign(A.excutiontrace[t].begin() + numberofselector,
                        A.excutiontrace[t].begin() + numberofselector + 5);
            auto region17 = add5gate(numberofselector, numberofwire, wire);
            A.excutiontrace.replaceSubmatrix({t, 0}, region17);
            initialworkVar[t - (6856 * i + 6848)] = (bitset<32>(A.excutiontrace[t][9].getInt64())).to_ulong();
        }
    }
    auto end1 = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> duration1 = end1 - start1;
    cout << endl;
    cout << "생성된 문자열: " << input << endl;
    cout << "회로 최종 결과 : ";
    string tmpstr = "";
    for (uint32_t elem : initialworkVar)
    {

        tmpstr += bitset<32>(elem).to_string();
    }
    cout << binToHex(tmpstr) << endl;
    cout << "해시값(Hex): " << finalHashHex << endl;
    if (compareHexStrings(binToHex(tmpstr), finalHashHex))
    {
        cout << "검사: true" << endl;
    }
    else
    {
        cout << "검사: false" << endl;
    }

    return make_tuple(A, duration1.count());
}
