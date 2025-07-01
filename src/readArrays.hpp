#ifndef READARRAYS_HPP_
#define READARRAYS_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

// Error codes:
// -1: File open error
//  0: Parsing error
//  1: Success
template <typename T>
int read4D(const std::string &filename, std::vector<T> &data, int dim1, int dim2, int dim3, int dim4) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return -1;
    }

    if (dim1 == 0 || dim2 == 0 || dim3 == 0 || dim4 == 0) {
        std::cerr << "Error: Dimensions must be greater than zero." << std::endl;
        return 0;
    }

    std::string line;
    std::getline(file, line);

    int file_dim2, file_dim3, file_dim4;
    std::istringstream iss_dim(line);
    if (!(iss_dim >> file_dim2 >> file_dim3 >> file_dim4)) {
        std::cerr << "Error: First line of " << filename << " does not contain three dimension integers." << std::endl;
        return 0;
    }
    if (file_dim2 != dim2 || file_dim3 != dim3 || file_dim4 != dim4) {
        std::cerr << "Error: Dimension mismatch in " << filename << ". Expected (" 
                  << dim2 << ", " << dim3 << ", " << dim4 << "), got (" 
                  << file_dim2 << ", " << file_dim3 << ", " << file_dim4 << ")." << std::endl;
        return 0;
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int i, j, k;
        if (!(iss >> i >> j >> k)) {
            std::cerr << "Error parsing indices in line: " << line << std::endl;
            return 0;
        }

        i--;
        j--;
        k--;

        int col = 0;
        T value;
        while (iss >> value) {
            size_t idx = ((col * dim2 + i) * dim3 + j) * dim4 + k;

            if (idx >= data.size()) {
                std::cerr << "Index out of bounds while reading " << filename << ": " << idx << "(" << col << ", " << i << ", " << j << ", " << k << \
                    ") for data size " << data.size() << std::endl;
                return 0;
            }

            data[idx] = value;
            col++;
        }
    }

    return 1;
}

#endif