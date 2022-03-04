## Test environments

* local OS X install, R 4.1.2
* Rhub (Windows Server 2022, R-devel, 64 bit; Ubuntu Linux 20.04.1 LTS, R-release, GCC; Fedora Linux, R-devel, clang, gfortran; Debian Linux, R-devel, GCC ASAN/UBSAN)
* win-builder (devel and release)

## R CMD check results

Package fails checks on Rhub - Debian Linux, R-devel, GCC ASAN/UBSAN with the error message that the compilation of Rcpp fails.
This leads to further failures, in particular as a result of the Rcpp installation failure the packages 'shiny', 'shinyjs', 
'shinydashboard' and 'shinydashboardPlus' that are used by the deBif package can not be installed.


On all other platforms the package passed checks successfully with only a note that some subdirectories (in particular the
subdirectory containing the user manual) are 1Mb or more

## Changes since last update

First release
