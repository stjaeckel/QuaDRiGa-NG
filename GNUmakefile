# This Makefile is for Linux / GCC environments

all:   quadriga-lib   moxunit-lib

quadriga-lib:
	cmake -B build_linux -D CMAKE_INSTALL_PREFIX=.
	cmake --build build_linux -j32 
	cmake --install build_linux

quadriga-lib_hdf5:
	cmake -B build_hdf5 -D HDF5_QD=ON -D CMAKE_INSTALL_PREFIX=.
	cmake --build build_hdf5 -j32 
	cmake --install build_hdf5

moxunit-lib:
	- rm -rf external/MOxUnit-master
	unzip external/MOxUnit.zip -d external/

clean:
	- rm -rf build
	- rm -rf build*
	- rm -rf quadriga_lib
