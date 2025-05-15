# This Makefile builds Qadriga-Lib for Windows

all:   quadriga-lib   moxunit-lib

quadriga-lib:
	cmake -B build_windows -D HDF5_QD=ON -D CMAKE_INSTALL_PREFIX=.
	cmake --build build_windows --config Release -- /m
	cmake --install build_windows

moxunit-lib:
	cd external
	tar -xf MOxUnit.zip
	cd ..

clean:
	- rmdir /s /q build_windows
	- rmdir /s /q quadriga_lib
