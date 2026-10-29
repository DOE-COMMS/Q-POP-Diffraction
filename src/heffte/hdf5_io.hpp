#ifndef HDF5_IO_HPP_
#define HDF5_IO_HPP_

#include <array>
#include <string>
#include <vector>

#include <hdf5.h>
#include <mpi.h>

#include "util.hpp"

// FerroDyn writes one dataset per component, shaped (nz, ny, nx) with x fastest.
// Every rank reads and writes a block of it collectively; blocks keep that layout.

inline hid_t openParallel(const std::string& path, unsigned flags, MPI_Comm comm) {
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL);
    hid_t file = H5Fopen(path.c_str(), flags, fapl);
    H5Pclose(fapl);
    if (file < 0) fail("Cannot open " + path);
    return file;
}

inline hid_t collective() {
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
    return dxpl;
}

// Read the block start (z, y, x), count (z, y, x) of a FerroDyn dataset.
template <typename T>
void readBlock(hid_t file, const std::string& path, const std::string& name, hid_t memType,
               const std::array<long, 3>& full, std::array<hsize_t, 3> start, std::array<hsize_t, 3> count,
               std::vector<T>& out) {
    hid_t dset = H5Dopen2(file, name.c_str(), H5P_DEFAULT);
    if (dset < 0) fail("Missing dataset '" + name + "' in " + path);
    hid_t fspace = H5Dget_space(dset);

    int rank = H5Sget_simple_extent_ndims(fspace);
    std::vector<hsize_t> dims(rank);
    H5Sget_simple_extent_dims(fspace, dims.data(), nullptr);
    if (rank != 3 || dims[0] != hsize_t(full[2]) || dims[1] != hsize_t(full[1]) || dims[2] != hsize_t(full[0])) {
        std::string shape;
        for (auto d : dims) shape += (shape.empty() ? "" : ", ") + std::to_string(d);
        fail("Dataset '" + name + "' in " + path + " has shape (" + shape + "), expected (nz, ny, nx) = (" +
             std::to_string(full[2]) + ", " + std::to_string(full[1]) + ", " + std::to_string(full[0]) +
             ") from system_setting.json.  Files written before 2026-05-26 are (nx, ny, nz): "
             "convert them with qpop_monitor.migrate_hdf5_layout.");
    }

    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
    hid_t mspace = H5Screate_simple(3, count.data(), nullptr);
    hid_t dxpl = collective();
    out.resize(count[0] * count[1] * count[2]);
    if (H5Dread(dset, memType, mspace, fspace, dxpl, out.data()) < 0) fail("Cannot read '" + name + "' from " + path);

    H5Pclose(dxpl);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dset);
}

// Write /intensity (nz, ny, nx) into the file prepared by qpop-diffraction --mgpu.
// The rank holds rows y' = (yLow + j + ny/2) % ny for j < nyl, sorted by y' (see the intensity kernel),
// so its rows are [y0', ny) and, when the fftshift wraps, [0, nyl - (ny - y0')).
inline void writeIntensity(const std::string& path, MPI_Comm comm, long nx, long ny, long nz,
                           long yLow, long nyl, const std::vector<double>& data) {
    hid_t file = openParallel(path, H5F_ACC_RDWR, comm);
    if (H5Lexists(file, "intensity", H5P_DEFAULT) > 0) H5Ldelete(file, "intensity", H5P_DEFAULT);

    hsize_t dims[3] = {hsize_t(nz), hsize_t(ny), hsize_t(nx)};
    hid_t fspace = H5Screate_simple(3, dims, nullptr);
    hid_t dset = H5Dcreate2(file, "intensity", H5T_IEEE_F64LE, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    long y0 = (yLow + ny / 2) % ny;
    long high = std::min(nyl, ny - y0);
    hsize_t start[3] = {0, hsize_t(y0), 0};
    hsize_t count[3] = {hsize_t(nz), hsize_t(high), hsize_t(nx)};
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
    if (high < nyl) {
        start[1] = 0;
        count[1] = nyl - high;
        H5Sselect_hyperslab(fspace, H5S_SELECT_OR, start, nullptr, count, nullptr);
    }

    hsize_t mdims[3] = {hsize_t(nz), hsize_t(nyl), hsize_t(nx)};
    hid_t mspace = H5Screate_simple(3, mdims, nullptr);
    hid_t dxpl = collective();
    if (H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, dxpl, data.data()) < 0) fail("Cannot write /intensity to " + path);

    H5Pclose(dxpl);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dset);
    H5Fclose(file);
}

#endif
