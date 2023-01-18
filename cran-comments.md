## Test environments

* local OS X install, R 4.2.2
* Rhub (Windows Server 2022, R-devel, 64 bit, R-devel, 64 bit; Ubuntu Linux 20.04.1 LTS, R-release, GCC; Fedora Linux, R-devel, clang, gfortran; Debian Linux, R-devel, GCC ASAN/UBSAN)
* win-builder (devel and release)

## R CMD check results

On all platforms the package passed checks successfully 

## Changes since last update

Minor bugs fixed. Changed all calls to sprintf() in the C code to snprintf(), as sprintf() is deprecated


