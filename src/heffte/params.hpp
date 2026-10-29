#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <array>
#include <string>
#include <vector>

#include <hdf5.h>

#include "util.hpp"

// Inputs of qpop_diffraction_mgpu, written by `qpop-diffraction --mgpu` (io_utils.write_mgpu_params).
struct mgpuParams {
    std::array<long, 3> full;       // FerroDyn grid (nx, ny, nz)
    std::array<long, 6> crop;       // x0, x1, y0, y1, z0, z1, inclusive
    std::array<double, 3> spacing;  // dx, dy, dz in nm
    std::array<double, 3> G, qC;    // nm^-1
    double V;                       // unit-cell volume, nm^3
    double uScale;                  // FerroDyn displacement (m) -> nm

    int nAtom, nPhase, nStruc, nMaterial;
    std::vector<double> R;          // (nAtom, 3) atom positions, nm
    std::vector<double> fRe, fIm;   // (nPhase, nAtom) form factors
    std::vector<double> su;         // (3, nStruc, nPhase, nAtom), nm per C/m^2
    std::vector<int> phaseOf;       // phase of each material, -1 if it does not scatter
    std::vector<double> dc;         // c_m(T) / aC_z - 1 of each material

    std::string plattFile, plattPrefix, displacementFile, materialFile, phaseFile, regionFile, outputFile;

    long nx() const { return crop[1] - crop[0] + 1; }
    long ny() const { return crop[3] - crop[2] + 1; }
    long nz() const { return crop[5] - crop[4] + 1; }

    void read(const std::string& path);
};

template <typename T>
std::vector<T> readDataset(hid_t file, const char* name, hid_t memType, std::vector<hsize_t>& dims) {
    hid_t dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) fail(std::string("Missing dataset ") + name + " in the multi-GPU parameters.");
    hid_t space = H5Dget_space(dset);
    dims.resize(H5Sget_simple_extent_ndims(space));
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);

    hsize_t n = 1;
    for (auto d : dims) n *= d;
    std::vector<T> data(n);
    H5Dread(dset, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    H5Sclose(space);
    H5Dclose(dset);
    return data;
}

inline std::string readString(hid_t file, const char* name) {
    if (H5Aexists(file, name) <= 0) return "";
    hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
    hid_t type = H5Aget_type(attr);
    std::string s(H5Tget_size(type), '\0');
    H5Aread(attr, type, s.data());
    H5Tclose(type);
    H5Aclose(attr);
    return s.substr(0, s.find('\0'));
}

inline void mgpuParams::read(const std::string& path) {
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) fail("Cannot open multi-GPU parameters " + path);

    std::vector<hsize_t> dims;
    auto grid = readDataset<long>(file, "grid", H5T_NATIVE_LONG, dims);
    auto limits = readDataset<long>(file, "crop", H5T_NATIVE_LONG, dims);
    auto h = readDataset<double>(file, "spacing", H5T_NATIVE_DOUBLE, dims);
    auto g = readDataset<double>(file, "G", H5T_NATIVE_DOUBLE, dims);
    auto q = readDataset<double>(file, "qC", H5T_NATIVE_DOUBLE, dims);
    for (int d = 0; d < 3; d++) {
        full[d] = grid[d];
        spacing[d] = h[d];
        G[d] = g[d];
        qC[d] = q[d];
    }
    for (int d = 0; d < 6; d++) crop[d] = limits[d];
    V = readDataset<double>(file, "V", H5T_NATIVE_DOUBLE, dims)[0];
    uScale = readDataset<double>(file, "u_scale", H5T_NATIVE_DOUBLE, dims)[0];

    R = readDataset<double>(file, "R", H5T_NATIVE_DOUBLE, dims);
    nAtom = dims[0];
    fRe = readDataset<double>(file, "fAtom_re", H5T_NATIVE_DOUBLE, dims);
    nPhase = dims[0];
    fIm = readDataset<double>(file, "fAtom_im", H5T_NATIVE_DOUBLE, dims);
    su = readDataset<double>(file, "su", H5T_NATIVE_DOUBLE, dims);
    nStruc = dims[1];
    phaseOf = readDataset<int>(file, "phase_of", H5T_NATIVE_INT, dims);
    nMaterial = dims[0];
    dc = readDataset<double>(file, "dc", H5T_NATIVE_DOUBLE, dims);

    plattFile = readString(file, "platt_file");
    plattPrefix = readString(file, "platt_prefix");
    displacementFile = readString(file, "displacement_file");
    materialFile = readString(file, "material_file");
    phaseFile = readString(file, "phase_file");
    regionFile = readString(file, "region_file");
    outputFile = readString(file, "output_file");

    H5Fclose(file);

    if (nStruc != 3) fail("qpop_diffraction_mgpu expects nStruc = 3 (P_x, P_y, P_z).");
}

#endif
