# This Makefile builds Qadriga-Lib for Windows

# Set path to your MATLAB installation (optional):
MATLAB_PATH = C:\Program Files\MATLAB\R2022b

# External libraries
quadriga_lib_version = 0.3.0

all:  quadriga-lib

quadriga-lib:   clean
	tar -xf external/quadriga_lib-$(quadriga_lib_version).zip
	move quadriga_lib-$(quadriga_lib_version) external

	cd external\quadriga_lib-$(quadriga_lib_version)
	nmake armadillo-lib
	nmake pugixml-lib
	nmake hdf5-lib
	nmake
	cd ..
	cd ..

	- mkdir quadriga_lib
	- mkdir "quadriga_lib\+quadriga_lib"

	xcopy "external\quadriga_lib-$(quadriga_lib_version)\+quadriga_lib\*" "quadriga_lib\+quadriga_lib\" /i
	xcopy "external\quadriga_lib-$(quadriga_lib_version)\LICENSE" "quadriga_lib\+quadriga_lib\" /i

	- rmdir /s /q external\quadriga_lib-$(quadriga_lib_version)

clean:
	- rmdir /s /q external\quadriga_lib-$(quadriga_lib_version)
	- del "quadriga_lib\+quadriga_lib"\*.mexw64

