#ifndef WRITEARRAYS_HPP_
#define WRITEARRAYS_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

template <typename T>
void write4D(const std::string& filename, std::vector<T> data,
    int dim1, int dim2, int dim3, int dim4,
    bool VTK = false) {

    std::cout << "Writing data to " << filename << std::endl;
    if (VTK == false) {
        std::ofstream outfile(filename, std::ios::out);
        if (!outfile) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }
        // Write the dimensions in the first line
        outfile << std::setw(6) << dim2 << " " << std::setw(5) << dim3 << " " << std::setw(5) << dim4 << "\n";

        // Write the data in the requested format
        for (int i = 0; i < dim2; i++) {
            for (int j = 0; j < dim3; j++) {
                for (int k = 0; k < dim4; k++) {
                    outfile << std::setw(6) << i+1 << " "
                        << std::setw(5) << j+1 << " "
                        << std::setw(5) << k+1 << " ";
                    for (int p = 0; p < dim1; p++) {
                        int index = p * dim2 * dim3 * dim4 + i * dim3 * dim4 + j * dim4 + k;
                        outfile << std::setw(14) << std::scientific << std::setprecision(7) << data[index] << " ";
                    }
                    outfile << "\n";
                }
            }
        }

        outfile.close();
    }
    else {
        std::cerr << "VTK format is not implemented yet." << std::endl;
    }
}
#endif