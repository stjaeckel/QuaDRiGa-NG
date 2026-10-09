# This Makefile is for Linux / GCC environments

octave = ON
matlab = ON

OCTAVE_VERSION := $(shell mkoctfile -v 2>/dev/null)
MATLAB_BIN = matlab24

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

# Tests
test:   moxunit-lib
ifeq ($(octave),ON)
ifneq ($(OCTAVE_VERSION),)
	octave --eval "cd xunit_new; run_tests;"
endif
endif
ifeq ($(matlab),ON)
	$(MATLAB_BIN) -batch "run('xunit_new/run_tests.m');"
endif

clean:
	- rm -rf build
	- rm -rf build*
	- rm -rf quadriga_lib
