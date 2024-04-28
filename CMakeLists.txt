cmake_minimum_required(VERSION 3.20)

project(
  Q-POP-DIFFRACTION
  VERSION 0.1.0
  DESCRIPTION "The mupro metal model"
  LANGUAGES Fortran C CXX)

  set(exename Diffraction)

  # if(MUPRODEV)
      # message(STATUS "Using openmupro from source")
      cmake_path(GET CMAKE_BINARY_DIR FILENAME BuildType)
      set(openmupro_root "/home/tanwen/muprosdk/openmupro")
      set(openmupro_DIR ${openmupro_root}/out/build/${BuildType}/library CACHE STRING "The openmupro directory")
  # endif()

  
  message(STATUS "The open mupro dir ${openmupro_DIR}")
  
  find_package(openmupro)
  find_package(MPI REQUIRED)
  find_package(toml-f REQUIRED)
  
  add_executable(${exename})
  
  target_sources(${exename} PRIVATE src/diffraction_main.f90 src/diffraction_routines.f90)
  
  target_link_libraries(${exename} PUBLIC mupro::openmupro)
  
  set_property(TARGET ${exename} PROPERTY LINKER_LANGUAGE Fortran)
  
  set_target_properties(
    ${exename}
    PROPERTIES POSITION_INDEPENDENT_CODE TRUE
               OUTPUT_NAME ${exename}
               Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
  
  add_custom_target(
    examples ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Examples
    COMMENT "Extract the example files"
  )

  add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Examples
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND tar -xf ${CMAKE_CURRENT_SOURCE_DIR}/docs/examples.tar.xz -C ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Decompress"
  VERBATIM USES_TERMINAL)