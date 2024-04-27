set(exename Diffraction)

if(MUPRODEV)
    message(STATUS the muprod dev is true)
    cmake_path(GET CMAKE_BINARY_DIR FILENAME CurrentDir)
    set(root "/home/tanwen/muprosdk")
    # set(mupro_DIR ${root}/library/out/build/${CurrentDir})
    set(openmupro_DIR ${root}/openmupro/out/build/${CurrentDir}/library)
    set(license_DIR ${root}/license/out/build/${CurrentDir}/library)
endif()

message(STATUS the open mupro dir ${openmupro_DIR})

find_package(openmupro)
# find_package(mupro REQUIRED)
find_package(MPI REQUIRED)
find_package(toml-f REQUIRED)

add_executable(${exename})

target_sources(${exename} PRIVATE diffraction_main.f90 diffraction_routines.f90)

target_link_libraries(${exename} PUBLIC mupro::openmupro)

# target_compile_options(${exename} PUBLIC -g -traceback)
target_compile_options(${exename} PUBLIC "-O3 -align rec16byte -align array64byte -align all -r8  -traceback")

set_property(TARGET ${exename} PROPERTY LINKER_LANGUAGE Fortran)

set_target_properties(
  ${exename}
  PROPERTIES POSITION_INDEPENDENT_CODE TRUE
             OUTPUT_NAME ${exename}
             Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
    