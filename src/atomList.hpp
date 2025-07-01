#ifndef ATOMLIST_HPP_
#define ATOMLIST_HPP_

#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

#include "constants.hpp"

struct atomList {
    int nPhase, nStruc;

    double q00[3];
    double aC[3][3];
    double** xAtom;
    std::complex<double>** fAtom;
    double** fAtomR, ** fAtomI;
    double*** su1, *** su2, *** su3;

    int nAtom;

    bool read(const std::string& filename, int nPhase, int nStruc);
    bool validate();
    bool setup();

    ~atomList() {
        for (int i = 0; i < nAtom; ++i) {
            delete[] xAtom[i];
        }
        delete[] xAtom;

        for (int i = 0; i < nPhase; ++i) {
            delete[] fAtom[i];
            delete[] fAtomR[i];
            delete[] fAtomI[i];
        }
        delete[] fAtom;
        delete[] fAtomR;
        delete[] fAtomI;

        for (int i = 0; i < nStruc; ++i) {
            for (int j = 0; j < nPhase; ++j) {
                delete[] su1[i][j];
                delete[] su2[i][j];
                delete[] su3[i][j];
            }
            delete[] su1[i];
            delete[] su2[i];
            delete[] su3[i];
        }
        delete[] su1;
        delete[] su2;
        delete[] su3;
    }
};

bool atomList::read(const std::string& filename, int nPhase, int nStruc) {
    this->nPhase = nPhase;
    this->nStruc = nStruc;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::vector<double> numbers;

    while (std::getline(infile, line)) {

        auto excl = line.find('!');
        if (excl != std::string::npos) line = line.substr(0, excl);

        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](int ch) { return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](int ch) { return !std::isspace(ch); }).base(), line.end());
        if (line.empty()) continue;

        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            numbers.push_back(val);
        }
    }

    int idx = 0;
    q00[0] = numbers[idx++];
    q00[1] = numbers[idx++];
    q00[2] = numbers[idx++];

    aC[0][0] = numbers[idx++];
    aC[0][1] = numbers[idx++];
    aC[0][2] = numbers[idx++];
    aC[1][0] = numbers[idx++];
    aC[1][1] = numbers[idx++];
    aC[1][2] = numbers[idx++];
    aC[2][0] = numbers[idx++];
    aC[2][1] = numbers[idx++];
    aC[2][2] = numbers[idx++];

    nAtom = static_cast<int>(numbers[idx++]);

    if (nAtom <= 0) {
        std::cerr << "Invalid number of atoms: " << nAtom << std::endl;
        return false;
    }

    std::cout << "Reading atom coordinates" << std::endl;
    xAtom = new double* [nAtom];
    for (int i = 0; i < nAtom; ++i) {
        xAtom[i] = new double[3];
        xAtom[i][0] = numbers[idx++];
        xAtom[i][1] = numbers[idx++];
        xAtom[i][2] = numbers[idx++];
    }

    std::cout << "Reading atom form factors" << std::endl;
    fAtom = new std::complex<double>*[nPhase];
    fAtomR = new double* [nPhase];
    fAtomI = new double* [nPhase];

    su1 = new double** [nStruc];
    su2 = new double** [nStruc];
    su3 = new double** [nStruc];
    for (int i = 0; i < nStruc; i++) {
        su1[i] = new double* [nPhase];
        su2[i] = new double* [nPhase];
        su3[i] = new double* [nPhase];
        for (int j = 0; j < nPhase; j++) {
            su1[i][j] = new double[nAtom];
            su2[i][j] = new double[nAtom];
            su3[i][j] = new double[nAtom];
        }
    }

    for (int i = 0; i < nPhase; i++) {
        fAtom[i] = new std::complex<double>[nAtom];
        fAtomR[i] = new double[nAtom];
        fAtomI[i] = new double[nAtom];

        int m;
        m = static_cast<int>(numbers[idx++]) - 1;
        std::cout << "Read input parameters for phase # " << m + 1 << std::endl;

        for (int j = 0; j < nAtom; j++) {
            fAtomR[m][j] = numbers[idx++];
            fAtomI[m][j] = numbers[idx++];
            fAtom[m][j] = std::complex<double>(fAtomR[m][j], fAtomI[m][j]);
        }

        for (int j = 0; j < nAtom; j++) {
            for (int k = 0; k < nStruc; k++) {
                su1[k][m][j] = numbers[idx++];
                su2[k][m][j] = numbers[idx++];
                su3[k][m][j] = numbers[idx++];
            }
        }
    }

    return true;
}

bool atomList::validate() {
    std::cout << "Validating read atom parameters" << std::endl;
    std::cout << "q00: ";
    for (int i = 0; i < 3; ++i) std::cout << q00[i] << " ";
    std::cout << std::endl;

    std::cout << "aC:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            std::cout << aC[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << "nAtom: " << nAtom << std::endl;

    std::cout << "xAtom:" << std::endl;
    for (int i = 0; i < nAtom; ++i) {
        std::cout << "  ";
        for (int j = 0; j < 3; ++j)
            std::cout << xAtom[i][j] << " ";
        std::cout << std::endl;
    }

    std::cout << "fAtom (per phase):" << std::endl;
    for (int i = 0; i < nPhase; ++i) {
        std::cout << "  Atom " << i << ": ";
        for (int p = 0; p < nAtom; ++p) {
            std::cout << "(" << fAtomR[i][p] << "," << fAtomI[i][p] << ") ";
        }
        std::cout << std::endl;
    }

    // su1, su2, su3
    std::cout << "su1, su2, su3 (per structure, phase, atom):" << std::endl;
    for (int i = 0; i < nPhase; ++i) {
        std::cout << "  Phase " << i << ":" << std::endl;
        for (int s = 0; s < nAtom; ++s) {
            std::cout << "    Atom " << s << ": ";
            for (int p = 0; p < nStruc; ++p) {
                std::cout << su1[p][i][s] << " " << su2[p][i][s] << " " << su3[p][i][s] << std::endl;
            }
            std::cout << std::endl;
        }
    }


    return true;
}

bool atomList::setup() {
    std::cout << "Setting up atom parameters" << std::endl;

    for (int i = 0; i < nStruc; i++) {
        for (int j = 0; j < nPhase; j++) {
            for (int k = 0; k < nAtom; k++) {
                su1[i][j][k] /= (l0/os0);
                su2[i][j][k] /= (l0/os0);
                su3[i][j][k] /= (l0/os0);
            }
        }
    }

    return true;
}

#endif