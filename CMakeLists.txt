cmake_minimum_required(VERSION 3.18)
project(QuaDRiGa-NG LANGUAGES C CXX)

# Options for Quadriga-lib components
option(ENABLE_MATLAB "Enable MATLAB MEX API" ON)
option(ENABLE_OCTAVE "Enable Octave MEX API" ON)
option(ENABLE_PYTHON "Enable Python API" ON)

# If HDF5_QD is ON, HDF5 is build from sources and statically linked
option(HDF5_QD "Link HDF5 statically" OFF)

# If ARMA_QD is ON, the Armadillo library from Quadriga-lib is used
option(ARMA_QD "Use armadillo headers provided with quadriga-lib" ON)

# Versions for external libraries (adjust as needed)
set(quadriga_lib_version "0.10.0")

# Set the C++ standard and basic compile flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -fopenmp")

# Config options for HDF5
set(BUILD_HDF5 ${HDF5_QD})
if (NOT BUILD_HDF5)
    find_package(HDF5 COMPONENTS C HL)
    if(HDF5_FOUND)
        get_filename_component(HDF5_LIB_DIR "${HDF5_LIBRARIES}" DIRECTORY)
        message(STATUS "HDF5 Library dir: ${HDF5_LIB_DIR}")
    else()
        set(BUILD_HDF5 ON)
    endif()
endif()

# Config options for Quadriga-lib
include(ExternalProject)
set(QUADRIGA_LIB_ZIP "${CMAKE_SOURCE_DIR}/external/quadriga_lib-${quadriga_lib_version}.zip")
set(QUADRIGA_LIB_ARGS
    URL ${QUADRIGA_LIB_ZIP}
    BINARY_DIR "${CMAKE_BINARY_DIR}/quadriga_lib-${quadriga_lib_version}_build"
    INSTALL_DIR "${CMAKE_SOURCE_DIR}/quadriga_lib"
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        -DENABLE_MATLAB=${ENABLE_MATLAB}
        -DENABLE_OCTAVE=${ENABLE_OCTAVE}
        -DENABLE_PYTHON=${ENABLE_PYTHON}
        -DHDF5_STATIC=${BUILD_HDF5}
        -DARMA_EXT=${ARMA_QD}
)
if (MSVC)
    list(APPEND QUADRIGA_LIB_ARGS DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
endif()
set(QUADRIGA_LIB_A "${CMAKE_SOURCE_DIR}/quadriga_lib/lib/libquadriga.a")
set(QUADRIGA_LIB_HDF5 "${CMAKE_SOURCE_DIR}/quadriga_lib/lib/libhdf5.a")
ExternalProject_Add(quadriga_lib ${QUADRIGA_LIB_ARGS})
