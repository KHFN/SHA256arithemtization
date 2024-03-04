#include "mcl/bls12_381.hpp"
#include <vector>
#include <iostream>
#include <thread>
#include <fstream>

using namespace mcl::bn;
using namespace std;

template <typename T>
void printVector(const vector<T> &vec)
{
    cout << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        cout << vec[i];
        if (i < vec.size() - 1)
        {
            cout << ", ";
        }
    }
    cout << "]" << endl;
}

template <typename Fr>
void printMatrix(const vector<vector<Fr>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (const auto &elem : row)
        {
            cout << elem << " ";
        }
        cout << "\n"; // 행이 끝날 때마다 줄바꿈
    }
}

class matrixFr
{

public:
    vector<vector<Fr>> data;

    matrixFr() = default;

    matrixFr(size_t outerSize, size_t innerSize)
    {
        data.resize(outerSize, vector<Fr>(innerSize, 0));
    }

    vector<Fr> &operator[](size_t index)
    {
        return data[index];
    }

    const vector<Fr> &operator[](size_t index) const
    {
        return data[index];
    }

    pair<size_t, size_t> getSize() const
    {
        size_t numRows = data.size();
        size_t numCols = numRows ? data[0].size() : 0;
        return {numRows, numCols};
    }

    void initializeWithSizeOf(const matrixFr &other)
    {
        auto [rows, cols] = other.getSize();
        data.resize(rows, vector<Fr>(cols));
    }

    void insertColumn(size_t position, const vector<Fr> &column)
    {
        if (data.empty())
        {
            throw std::invalid_argument("Matrix is empty.");
        }
        if (position > data[0].size())
        {
            throw std::out_of_range("Column insert position is out of range.");
        }
        if (column.size() != data.size())
        {
            throw std::invalid_argument("Column size does not match matrix row size.");
        }
        for (size_t i = 0; i < data.size(); ++i)
        {
            data[i].insert(data[i].begin() + position, column[i]);
        }
    }

    vector<Fr> getColumn(size_t position) const
    {

        if (data.empty() || position >= data[0].size())
        {
            throw std::out_of_range("Column position is out of range.");
        }

        vector<Fr> column(data.size());
        for (size_t i = 0; i < data.size(); ++i)
        {
            column[i] = data[i][position];
        }

        return column;
    }

    void replaceColumn(size_t position, const vector<Fr> &newColumn)
    {
        if (data.empty())
        {
            throw std::invalid_argument("Matrix is empty.");
        }
        if (position >= data[0].size())
        {
            throw std::out_of_range("Column position is out of range.");
        }
        if (newColumn.size() != data.size())
        {
            throw std::invalid_argument("New column size does not match matrix row size.");
        }
        for (size_t i = 0; i < data.size(); ++i)
        {
            data[i][position] = newColumn[i];
        }
    }

    void replaceRow(size_t position, const vector<Fr> &newRow)
    {
        if (position >= data.size())
        {
            throw std::out_of_range("Row position is out of range.");
        }
        if (!data.empty() && newRow.size() != data[0].size())
        {
            throw std::invalid_argument("New row size does not match matrix column size.");
        }
        data[position] = newRow;
    }

    void replaceSubmatrix(const vector<uint32_t> &topLeftIndex, const vector<vector<Fr>> &replacement)
    {

        if (topLeftIndex.size() != 2)
        {
            throw std::invalid_argument("topLeftIndex must have exactly two elements.");
        }

        size_t startRow = topLeftIndex[0];
        size_t startCol = topLeftIndex[1];

        if (startRow + replacement.size() > data.size() || (startCol + replacement[0].size() > data[0].size()))
        {
            throw std::out_of_range("Replacement data exceeds matrix bounds.");
        }

        for (size_t i = 0; i < replacement.size(); ++i)
        {
            for (size_t j = 0; j < replacement[i].size(); ++j)
            {
                if (startRow + i < data.size() && startCol + j < data[i].size())
                {
                    data[startRow + i][startCol + j] = replacement[i][j];
                }
            }
        }
    }

    template <typename T>
    void swapElements(const vector<T> &pos1, const vector<T> &pos2)
    {

        if (pos1.size() != 2 || pos2.size() != 2)
        {
            throw std::invalid_argument("Position vectors must have exactly two elements (row and column).");
        }

        size_t row1 = static_cast<size_t>(pos1[0]);
        size_t col1 = static_cast<size_t>(pos1[1]);
        size_t row2 = static_cast<size_t>(pos2[0]);
        size_t col2 = static_cast<size_t>(pos2[1]);

        if (row1 >= data.size() || row2 >= data.size() || col1 >= data[0].size() || col2 >= data[0].size())
        {
            throw std::out_of_range("Row or column position is out of range.");
        }

        Fr temp = data[row1][col1];
        data[row1][col1] = data[row2][col2];
        data[row2][col2] = temp;
    }

    template <typename T>
    void duplicateElement(const vector<T> &sourcePos, const vector<T> &destPos)
    {
        if (sourcePos.size() != 2 || destPos.size() != 2)
        {
            throw std::invalid_argument("Position vectors must have exactly two elements (row and column).");
        }

        size_t sourceRow = static_cast<size_t>(sourcePos[0]);
        size_t sourceCol = static_cast<size_t>(sourcePos[1]);
        size_t destRow = static_cast<size_t>(destPos[0]);
        size_t destCol = static_cast<size_t>(destPos[1]);

        if (sourceRow >= data.size() || destRow >= data.size() ||
            sourceCol >= data[0].size() || destCol >= data[0].size())
        {
            throw std::out_of_range("Row or column position is out of range.");
        }

        data[destRow][destCol] = data[sourceRow][sourceCol];
    }

    template <typename T>
    void swapMultipleElements(const vector<vector<T>> &positions)
    {
        for (const auto &posPair : positions)
        {
            if (posPair.size() != 2 || posPair[0].size() != 2 || posPair[1].size() != 2)
            {
                throw std::invalid_argument("Each position pair must contain exactly two elements (row and column).");
            }
            swapElements(posPair[0], posPair[1]);
        }
    }

    template <typename T>
    void duplicateMultipleElements(const vector<vector<T>> &positionPairs)
    {
        for (const auto &positionPair : positionPairs)
        {
            if (positionPair.size() != 2 || positionPair[0].size() != 2 || positionPair[1].size() != 2)
            {
                throw std::invalid_argument("Each position pair must contain exactly two elements (row and column).");
            }
            duplicateElement(positionPair[0], positionPair[1]);
        }
    }

    friend ostream &operator<<(ostream &os, const matrixFr &m)
    {
        for (const auto &row : m.data)
        {
            for (const auto &elem : row)
            {
                os << elem.getStr() << ' ';
            }
            os << '\n';
        }
        return os;
    }
};

void multiplyPartByConstant(vector<Fr> &vec, const Fr &constant, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        vec[i] *= constant;
    }
}

void mulveccons(vector<Fr> &vec, const Fr &constant)
{
    size_t length = vec.size();
    const size_t minSizeForParallel = (1 << 17);

    if (length < minSizeForParallel)
    {

        multiplyPartByConstant(vec, constant, 0, length);
    }
    else
    {

        unsigned numThreads = std::thread::hardware_concurrency();
        size_t partSize = length / numThreads;
        vector<thread> threads(numThreads);

        for (unsigned i = 0; i < numThreads; ++i)
        {
            size_t start = i * partSize;
            size_t end = (i == numThreads - 1) ? length : start + partSize;
            threads[i] = thread(multiplyPartByConstant, ref(vec), constant, start, end);
        }

        for (thread &t : threads)
        {
            t.join();
        }
    }
}

void addVectorSegment(const vector<Fr> &a, const vector<Fr> &b, vector<Fr> &result, size_t start, size_t end)
{
    for (size_t i = start; i < end; ++i)
    {
        result[i] = a[i] + b[i];
    }
}

vector<Fr> addVectors(const vector<Fr> &a, const vector<Fr> &b)
{
    const size_t threshold = 1 << 18;
    vector<Fr> result(a.size());

    if (a.size() < threshold || a.size() != b.size())
    {

        for (size_t i = 0; i < a.size(); ++i)
        {
            result[i] = a[i] + b[i];
        }
    }
    else
    {

        unsigned numThreads = thread::hardware_concurrency();
        vector<thread> threads(numThreads);
        size_t blockSize = a.size() / numThreads;

        for (unsigned i = 0; i < numThreads; ++i)
        {
            size_t start = i * blockSize;
            size_t end = (i + 1 == numThreads) ? a.size() : start + blockSize;
            threads[i] = thread(addVectorSegment, cref(a), cref(b), ref(result), start, end);
        }

        for (thread &t : threads)
        {
            t.join();
        }
    }

    return result;
}

template <typename T>
vector<T> concatenateVectors(const vector<T> &vec1, const vector<T> &vec2)
{
    vector<T> result;
    result.reserve(vec1.size() + vec2.size());
    result.insert(result.end(), vec1.begin(), vec1.end());
    result.insert(result.end(), vec2.begin(), vec2.end());
    return result;
}

template <typename T>
T dotProduct(const vector<T> &v1, const vector<T> &v2)
{
    if (v1.size() != v2.size())
    {
        throw std::invalid_argument("벡터의 크기가 서로 다릅니다.");
    }

    T result = 0;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

template <typename Fr>
vector<vector<Fr>> transpose(const vector<vector<Fr>> &matrix)
{
    if (matrix.empty())
        return {}; // 입력 벡터가 비어있는 경우 빈 벡터 반환

    size_t rows = matrix.size();    // 원본 행렬의 행 개수
    size_t cols = matrix[0].size(); // 원본 행렬의 열 개수

    vector<vector<Fr>> transposed(cols, vector<Fr>(rows));

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

void saveToCSV(const vector<vector<Fr>> &data, const string &filename)
{
    ofstream file(filename);

    for (const auto &row : data)
    {
        for (size_t i = 0; i < row.size(); ++i)
        {
            if (i > 0)
                file << ",";
            file << row[i].getStr(); // Fr 객체를 문자열로 변환하여 저장
        }
        file << "\n";
    }

    file.close();
}

void loadFromCSV(vector<vector<Fr>> &data, const string &filename)
{
    ifstream file(filename);
    string line;

    while (getline(file, line))
    {
        vector<Fr> row;
        stringstream ss(line);
        string value;

        while (getline(ss, value, ','))
        {
            Fr fr;
            fr.setStr(value);
            row.push_back(fr);
        }

        data.push_back(row);
    }

    file.close();
}
