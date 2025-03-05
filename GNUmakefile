# This Makefile builds Qadriga-Lib for Linux / GCC environments

# Set path to your MATLAB installation (optional):
# If left empty, the build script trys to autodetect the MATLAB path. 
# If none is found, MATLAB targets are not compiled.
MATLAB_PATH = #/usr/local/MATLAB/R2023a

# External libraries
# External libraries are located in the 'external' folder. Set the version numbers here.
quadriga_lib_version = 0.3.2

# Autodetect MATLAB path
ifeq ($(MATLAB_PATH),)
	MATLAB_PATH := $(shell readlink -f $(shell readlink -f $(shell which matlab)) | sed 's/\/bin\/matlab//')
endif

# Check if Octave is installed by trying to run mkoctfile
OCTAVE_VERSION := $(shell mkoctfile -v 2>/dev/null)

all:  quadriga-lib   moxunit-lib

quadriga-lib:   clean
	unzip external/quadriga_lib-$(quadriga_lib_version).zip -d external/

	( cd external/quadriga_lib-$(quadriga_lib_version) && make dirs && make armadillo-lib && make pugixml-lib && make pybind11-lib )
	( cd external/quadriga_lib-$(quadriga_lib_version) && make python -j16 && make documentation )

ifneq ($(MATLAB_PATH),)
	( cd external/quadriga_lib-$(quadriga_lib_version) && make mex_matlab -j16 )
endif
ifneq ($(OCTAVE_VERSION),)
	( cd external/quadriga_lib-$(quadriga_lib_version) && make mex_octave -j16 )
endif

	- mkdir quadriga_lib
	- mkdir quadriga_lib/+quadriga_lib
	- mkdir quadriga_lib/python

	cp -R external/quadriga_lib-$(quadriga_lib_version)/+quadriga_lib/* quadriga_lib/+quadriga_lib/
	cp -R external/quadriga_lib-$(quadriga_lib_version)/html_docu quadriga_lib/
	cp -R external/quadriga_lib-$(quadriga_lib_version)/lib/*.so quadriga_lib/python

	cp external/quadriga_lib-$(quadriga_lib_version)/LICENSE quadriga_lib/

	rm -rf external/quadriga_lib-$(quadriga_lib_version)

moxunit-lib:
	- rm -rf external/MOxUnit-master
	unzip external/MOxUnit.zip -d external/

# Cleanup
clean:
	- rm -rf external/quadriga_lib-$(quadriga_lib_version)
	- rm -rf quadriga_lib/+quadriga_lib/*.m
	- rm -rf quadriga_lib/+quadriga_lib/*.mex
	- rm -rf quadriga_lib/+quadriga_lib/*.mexa64
	- rm -rf quadriga_lib/html_docu
	- rm -rf quadriga_lib/python

tidy:   clean
	- rm -rf quadriga_lib/
	- rm -rf external/MOxUnit-master
