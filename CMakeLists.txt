cmake_minimum_required(VERSION 2.8)
project(REFPROP Fortran)

SET (USE_OPENMP OFF CACHE BOOL "Use OpenMP" FORCE)
SET (CMAKE_BUILD_TYPE "TESTING" CACHE INTERNAL "Set test build" FORCE)
         
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

message(STATUS "REFPROP: Build enabled")
message(STATUS "REFPROP: Each REFPROP instance will use up to ${OMP_NUM_PROCS} threads with OpenMP")

FILE(GLOB refprop_files   "${CMAKE_CURRENT_SOURCE_DIR}/fortran/*.FOR")
LIST(SORT refprop_files)

LIST(REMOVE_ITEM refprop_files
     "${CMAKE_CURRENT_SOURCE_DIR}/fortran/COMMONS.FOR"
     "${CMAKE_CURRENT_SOURCE_DIR}/fortran/COMTRN.FOR"
)

if(UNIX OR MINGW)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-ffast-math -fno-common -ffloat-store")
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/fortran")
add_library(REFPRP64 SHARED ${refprop_files})
add_executable(RGPTableGen ${refprop_files} "${CMAKE_CURRENT_SOURCE_DIR}/RGP.f90")
