cmake_minimum_required(VERSION 3.3)
SET(CMAKE_BUILD_TYPE TESTING)

project(RGPTableGen Fortran)


if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()


# Build Environment
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

if (NOT REFPROP_FLAGS_ADDED)
  set(REFPROP_FLAGS_ADDED "true" CACHE INTERNAL "Has refprop flags been added?")  
  
  if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
  else()
    SET (CMAKE_BUILD_TYPE "TESTING" CACHE INTERNAL "Set test build")
    INCLUDE("${CMAKE_MODULE_PATH}/SetFortranFlags.cmake") 
    INCLUDE("${CMAKE_MODULE_PATH}/SetCompileFlag.cmake")
  endif()
  
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-fpic")
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-fno-common")
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-fall-intrinsics")
  
#  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-ffast-math")
#  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "-ffloat-store")
endif()


# Options
SET (USE_OPENMP OFF CACHE BOOL "Use OpenMP")

find_package(PythonInterp REQUIRED)

#####################
# CONFIGURE REFPROP #
#####################

set(REFPROP_BASE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Refprop code location")

# Set locations to store output
set(REFPROP_FLUID_FOLDER "${CMAKE_BINARY_DIR}")
set(REFPROP_DLL_FILE     "${CMAKE_BINARY_DIR}") 


# Copy REFPROP fluid files
message(STATUS "REFPROP: REFPROP fluids found at ${REFPROP_BASE_PATH}, copying fluid files")
file(COPY "${REFPROP_BASE_PATH}/fluids" DESTINATION "${CMAKE_BINARY_DIR}")
message(STATUS "REFPROP: Fluid files copied")

# Get fortran sources
FILE(GLOB refprop_files   "${REFPROP_BASE_PATH}/fortran/*.FOR")

# Remove non-required items
LIST(REMOVE_ITEM refprop_files
"${REFPROP_BASE_PATH}/fortran/COMMONS.FOR" # include file
"${REFPROP_BASE_PATH}/fortran/COMTRN.FOR"  # include file
"${REFPROP_BASE_PATH}/fortran/UTILITY.FOR"
)

# Make a copy of files (which are to be patched)
message(STATUS "REFPROP: Patching fortran files")
file(COPY "${REFPROP_BASE_PATH}/fortran/COMMONS.FOR" DESTINATION "${CMAKE_BINARY_DIR}/fortran")
file(COPY "${REFPROP_BASE_PATH}/fortran/COMTRN.FOR"  DESTINATION "${CMAKE_BINARY_DIR}/fortran")
file(RENAME "${CMAKE_BINARY_DIR}/fortran/COMMONS.FOR" "${CMAKE_BINARY_DIR}/fortran/commons.for")
file(RENAME "${CMAKE_BINARY_DIR}/fortran/COMTRN.FOR" "${CMAKE_BINARY_DIR}/fortran/comtrn.for")

execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/patchUtility.py"
                "${REFPROP_BASE_PATH}/fortran/UTILITY.FOR"
                "${CMAKE_BINARY_DIR}/UTILITY.FOR")

LIST(APPEND refprop_files "${CMAKE_BINARY_DIR}/UTILITY.FOR")

# Set include directory
include_directories("${CMAKE_BINARY_DIR}/fortran")

# Setup build
include_directories("${REFPROP_BASE_PATH}/fortran")

# Use common core method
add_library(REFPROPObjects OBJECT ${refprop_files})
add_executable(RGPTableGen $<TARGET_OBJECTS:REFPROPObjects> "${CMAKE_CURRENT_SOURCE_DIR}/RGP.for")

target_compile_options(RGPTableGen PRIVATE -std=f2008 -ffree-form)
