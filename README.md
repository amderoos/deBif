# deBif

This directory contains the R package **deBif** which is aimed at bifurcation and phaseplane anaylysis of system of ordinary differential equations.

### Currently implemented methods:

* Phaseplane analysis of dynamical systems in terms of 1 or 2 ordinary differential equations
* Bifurcation analysis of dynamical systems in terms of ordinary differential equations.  
  Implemented methods:
    - Equilibrium continuation as function of one model parameter
    - Detection of saddle-node bifurcation points (limit points), transcritical bifurcation points (branching points) and Hopf bifurcation points
    - Continuation of limit cycles as function of one model parameter
    - Continuation of saddle-node bifurcation points (limit points), transcritical bifurcation points (branching points) and Hopf bifurcation points as function of 2 model parameters

The package **deBif** is primarily intended to be used in education. It is in no way meant to replace the ultimate program to properly analyse the bifurcations in dynamical systems, which is [Matcont](https://sourceforge.net/projects/matcont/). 

### Installation

You can install this package by issuing the following command from the R command line:
```
devtools::install_github("amderoos/deBif")
```

*****
