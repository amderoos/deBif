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

I developed the package **deBif** primarily to use in the modelling courses that I teach. It is in no way meant to replace the ultimate program to properly analyse the bifurcations in dynamical systems, which is [Matcont](https://sourceforge.net/projects/matcont/). If you are interested in a tutorial to learn more about bifurcation theory and for step-by-step exercises to use the **deBif** package, you can check out the [reader that I am using for teaching model analysis](https://staff.fnwi.uva.nl/a.m.deroos/projects/BifurcationTheory/).

### Installation

You can install this package by issuing the following command from the R command line:
```
devtools::install_github("amderoos/deBif")
```

*****
