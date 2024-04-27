# Q-POP-Diffraction

```sh
tar -cJf examples.tar.xz Examples/
tar -xJf examples.tar.xz
```

## Acknowledgement
This development of the software was supported as part of the Computational Materials Sciences Program funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences, under Award No. DE-SC0020145.

## Compiling with CMake
This repository uses CMake as the build system.
### Executable environment (more than one)
* cmake version 3.27.7
* ninja version 1.10.1
* Intel oneAPI 2023 or 2024
* build muprosdk and openmupro

### Building

**Note that we are using the `MUPRODEV` option here, which signifies that we are in developer mode and not install mupro::openmupro library, so it's necessary to specify the `compilation paths for mupro::openmupro` in the CMake file under the `src` folder:**
```cmake
if(MUPRODEV)
    message(STATUS the muprod dev is true)
    cmake_path(GET CMAKE_BINARY_DIR FILENAME CurrentDir)
    set(root "/home/tanwen/muprosdk")
    set(openmupro_DIR ${root}/openmupro/out/build/${CurrentDir}/library)
    set(license_DIR ${root}/license/out/build/${CurrentDir}/library)
endif()
```
Generate the build files using CMake:
```cmake
cmake --preset="linux-Debug-dev" .
cd out/build/debug
ninja
```
The generated executable file `Diffraction` are located in `out/build/debug/src`


