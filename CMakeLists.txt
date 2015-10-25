cmake_minimum_required(VERSION 2.8)
project(REFPROP Fortran)

SET (USE_OPENMP ON CACHE BOOL "Use OpenMP")
SET (CMAKE_BUILD_TYPE "RELEASE" CACHE INTERNAL "Set test build")
         
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp)


# running on windows and refprop is installed? then copy what we need
if (WIN32)
  # Set locations to store output
  set(REFPROP_FLUID_FOLDER "${CMAKE_BINARY_DIR}")
  set(REFPROP_DLL_FILE     "${CMAKE_BINARY_DIR}")
  
  # Refprop (new version) is a 64 bit binary
  find_path(REFPROP_BASE_PATH "REFPRP64.dll" "c:/program files/REFPROP" "c:/program files (x86)/REFPROP")
endif()


if (USE_OPENMP)
  message(STATUS "REFPROP: Each REFPROP instance may use up to ${OMP_NUM_PROCS} threads with OpenMP")
endif()

# Copy fluid folder to output directory
if (REFPROP_BASE_PATH)
  message(STATUS "REFPROP: REFPROP found at ${REFPROP_BASE_PATH}, copying/using existing fluid/fortran files")
else()
  set(REFPROP_BASE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Code Path")
endif()  

# Remove non required files 
message("${REFPROP_BASE_PATH}")
file(COPY "${REFPROP_BASE_PATH}/fluids" DESTINATION "${CMAKE_BINARY_DIR}")
FILE(GLOB refprop_files   "${REFPROP_BASE_PATH}/fortran/*.FOR")
LIST(REMOVE_ITEM refprop_files
"${REFPROP_BASE_PATH}/fortran/COMMONS.FOR"
"${REFPROP_BASE_PATH}/fortran/COMTRN.FOR"
"${REFPROP_BASE_PATH}/fortran/UTILITY.FOR"
)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "patchUtility.py" "${REFPROP_BASE_PATH}/fortran/UTILITY.FOR" "${CMAKE_CURRENT_SOURCE_DIR}/UTILITY.FOR")
LIST(APPEND refprop_files "${CMAKE_CURRENT_SOURCE_DIR}/UTILITY.FOR")
message("${refprop_files}")


# Add some nice flags
if(UNIX OR MINGW)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}" Fortran "${CMAKE_Fortran_FLAGS} -ffast-math -fno-common -ffloat-store" FORCE)
endif()

# Setup build
include_directories("${REFPROP_BASE_PATH}/fortran")

# Use common core method
add_library(REFPROPObjects OBJECT ${refprop_files})
add_library(REFPRP64 SHARED $<TARGET_OBJECTS:REFPROPObjects>)
add_executable(RGPTableGen $<TARGET_OBJECTS:REFPROPObjects> "${CMAKE_CURRENT_SOURCE_DIR}/RGP.f90")
