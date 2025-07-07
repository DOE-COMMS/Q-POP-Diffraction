#ifndef SIZECONTEXT_HPP_
#define SIZECONTEXT_HPP_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cufftXt.h>

struct sizeContext {
    double lx, ly, lz;
    int nx, ny, nz;
    int ns, nf;
    int nPhase, nStruc;
    double strainAvg[6];

    int n;
    double dx, dy, dz;
    int k1, k2;
    double kt;

    int ns1;
    int k0;

    double hf, hs, h1, h2, hd, h, hstep;

    int Rn1, Rn2, Rn3;
    int Cn1, Cn2, Cn3;

    bool read(const std::string& filename);
    bool setup();
};

bool sizeContext::read(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    std::vector<double> numbers;

    while (std::getline(infile, line)) {

        auto excl = line.find('!');
        if (excl !=  std::string::npos) line = line.substr(0, excl);

        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](int ch) { return !std::isspace(ch); }));
        line.erase(std::find_if(line.rbegin(), line.rend(), [](int ch) { return !std::isspace(ch); }).base(), line.end());
        if (line.empty()) continue;

        std::istringstream iss(line);
        double val;
        while (iss >> val) {
            numbers.push_back(val);
        }
    }

    if (numbers.size() < 16) {
        std::cerr << "Not enough parameters in input file." << std::endl;
        return false;
    }

    int idx = 0;
    lx = numbers[idx++];
    ly = numbers[idx++];
    lz = numbers[idx++];
    nx = static_cast<int>(numbers[idx++]);
    ny = static_cast<int>(numbers[idx++]);
    nz = static_cast<int>(numbers[idx++]);
    ns = static_cast<int>(numbers[idx++]);
    nf = static_cast<int>(numbers[idx++]);
    nPhase = static_cast<int>(numbers[idx++]);
    nStruc = static_cast<int>(numbers[idx++]);
    strainAvg[0] = numbers[idx++];
    strainAvg[1] = numbers[idx++];
    strainAvg[2] = numbers[idx++];
    strainAvg[3] = numbers[idx++];
    strainAvg[4] = numbers[idx++];
    strainAvg[5] = numbers[idx++];

    std::cout << "\nRead parameters:" << std::endl;
    std::cout << "  lx = " << lx << ", ly = " << ly << ", lz = " << lz << std::endl;
    std::cout << "  nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
    std::cout << "  ns = " << ns << ", nf = " << nf << std::endl;
    std::cout << "  nPhase = " << nPhase << ", nStruc = " << nStruc << std::endl;
    std::cout << "  strainAvg = [" << strainAvg[0] << ", " << strainAvg[1] << ", " 
              << strainAvg[2] << ", " << strainAvg[3] << ", " << strainAvg[4] 
              << ", " << strainAvg[5] << "]\n" << std::endl;

    if (nPhase < 1 || nPhase > 12) {
        std::cerr << "This program allows 1~12 phases only. You cannot claim " << nPhase << " phases." << std::endl;
        return false;
    }

    if (nStruc < 0 || nStruc > 12) {
        std::cerr << "This program allows 0~12 structural order parameters only. You cannot claim "
                  << nStruc << " structural order parameters." << std::endl;
        return false;
    }

    return true;
}

bool sizeContext::setup() {
    n = nx * ny * nz;

    Rn1 = nx;
    Rn2 = ny;
    Rn3 = nz;

    Cn1 = nx;
    Cn2 = ny;
    Cn3 = floor(nz/2) + 1;

    dx = lx / nx;
    dy = ly / ny;
    dz = lz / nz;

    k1 = ns + 1;
    k2 = ns + nf;

    if (nf == 0) {
        k1 = 1;
        k2 = nz;
    }

    ns1 = 0;
    k0 = ns1 + 1;

    if (nf != 0) {
        k1 = ns + 1; // film bottom
        k2 = ns + nf; // film surface

        if (k2 > nz) k2 = nz;
    }
    else {
        hf = nf * dz;
        hs = (ns - ns1 - 1) * dz;
        h1 = dz * (k1 - 1);
        h2 = dz * (k2 - 1);
        hd = h2 - h1;

        h = hs > hf ? hs : hf;
        hstep = dz;
    }

    // std::cout << "Size context setup complete:" << std::endl;
    // std::cout << "  lx = " << lx << ", ly = " << ly << ", lz = " << lz << std::endl;
    // std::cout << "  nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
    // std::cout << "  ns = " << ns << ", nf = " << nf << std::endl;
    // std::cout << "  nPhase = " << nPhase << ", nStruc = " << nStruc << std::endl;
    // std::cout << "  dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;
    // std::cout << "  k1 = " << k1 << ", k2 = " << k2 << std::endl;
    // std::cout << "  ns1 = " << ns1 << ", k0 = " << k0 << std::endl;
    // std::cout << "  hf = " << hf << ", hs = " << hs << ", h1 = " << h1
    //           << ", h2 = " << h2 << ", hd = " << hd << std::endl;

    return true;
}

#endif