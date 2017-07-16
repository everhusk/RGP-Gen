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

#####################
# CONFIGURE REFPROP #
#####################

# See if we can find REFPROP
option(USE_INSTALLED_REFPROP "Use the fortran from REFPROP folder" ON)
if (WIN32)
  # Refprop (new version) is a 64 bit binary
  find_path(REFPROP_BASE_PATH "REFPRP64.dll" "c:/program files/REFPROP" "c:/program files (x86)/REFPROP" "c:/REFPROP")
  message(STATUS "REFPROP Path: ${REFPROP_BASE_PATH}")
endif()

# REFPROP sources must be located in the current build folder
if (NOT USE_INSTALLED_REFPROP OR NOT DEFINED REFPROP_BASE_PATH OR NOT REFPROP_BASE_PATH)
  set(REFPROP_BASE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Refprop code location")
  message(STATUS "REFPROP Path: Using REFPROP in source directory")
endif()

# Copy REFPROP fluid files
FILE(GLOB fluid_files "${REFPROP_BASE_PATH}/fluids/*.fld")
if (NOT fluid_files)
    message(FATAL_ERROR "REFPROP Path: Fluids were not found in ${REFPROP_BASE_PATH}/fluids/")
endif()

message(STATUS "REFPROP: REFPROP fluids found at ${REFPROP_BASE_PATH}, copying fluid files")
file(COPY "${REFPROP_BASE_PATH}/fluids" DESTINATION "${CMAKE_BINARY_DIR}")
message(STATUS "REFPROP: Fluid files copied")

# Setup build
include_directories("${REFPROP_BASE_PATH}/fortran")

# Use common core method
add_library(REFPROPObjects OBJECT 
"${REFPROP_BASE_PATH}/fortran/CMNS.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_ANC.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_BWR.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_CPP.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_DE.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_ECS.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_FEQ.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_MLT.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_PH0.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_PR.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_QUI.FOR"
"${REFPROP_BASE_PATH}/fortran/CORE_STN.FOR"
"${REFPROP_BASE_PATH}/fortran/FLASH2.FOR"
"${REFPROP_BASE_PATH}/fortran/FLSH_SUB.FOR"
"${REFPROP_BASE_PATH}/fortran/IDEALGAS.FOR"
"${REFPROP_BASE_PATH}/fortran/MIX_AGA8.FOR"
"${REFPROP_BASE_PATH}/fortran/MIX_HMX.FOR"
"${REFPROP_BASE_PATH}/fortran/PASS_FTN.FOR"
"${REFPROP_BASE_PATH}/fortran/PROP_SUB.FOR"
"${REFPROP_BASE_PATH}/fortran/REALGAS.FOR"
"${REFPROP_BASE_PATH}/fortran/SAT_SUB.FOR"
"${REFPROP_BASE_PATH}/fortran/SETUP.FOR"
"${REFPROP_BASE_PATH}/fortran/SETUP2.FOR"
"${REFPROP_BASE_PATH}/fortran/TRNS_ECS.FOR"
"${REFPROP_BASE_PATH}/fortran/TRNS_TCX.FOR"
"${REFPROP_BASE_PATH}/fortran/TRNS_VIS.FOR"
"${REFPROP_BASE_PATH}/fortran/TRNSP.FOR"
"${REFPROP_BASE_PATH}/fortran/UTILITY.FOR")

add_executable(RGPTableGen $<TARGET_OBJECTS:REFPROPObjects> "${CMAKE_CURRENT_SOURCE_DIR}/RGP.for")
target_compile_options(RGPTableGen PRIVATE -std=f2008 -ffree-form)

