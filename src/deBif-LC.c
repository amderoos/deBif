/*
  NAME
    deBif

  PURPOSE
    Interface to the routines that are used for locating points on steady state
    bifurcation curves of non-linear ODE systems and curves determining dynamic
    regimes of such systems in two-parameter domains. More generally, the routines
    allow for computing fixed points of a system of non-linear equations.

    Copyright (C) 2021, Andre M. de Roos, University of Amsterdam

    This file is part of the deBif software package.

    deBif is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    deBif is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with deBif. If not, see <http://www.gnu.org/licenses/>.

    Last modification: AMdR - Mar 01, 2021
*/

#include "globals.h"

#define FEMTO                     1.0E-15

/*
 *====================================================================================================================================
 *  Definition of global variables and parameters
 *====================================================================================================================================
 */

#if (defined(CMDLINEDEBUG))

const int                         pntDim = 84, sysOdeDim = 2, freeparsdim = 1, glorder = 4, ninterval = 10, CDfinemeshdim = 41, MaxIter = 20;
static int                        CurveType = LC;

static double                     ATol = 1.0E-9, RTol = 1.0E-7, CTol = 1.0E-8;

#define STARTPNT                  1
#if (STARTPNT == 1)
// This is the starting point from the HP in rosenzweig
static double                     point[] = {0.12390952,
                                             0.02867143,  0.10991672, 0.02867020, 0.10993486, 0.02866653, 0.10995255,
                                             0.02866053, 0.10996937, 0.02865233, 0.10998489, 0.02864214, 0.10999872,
                                             0.02863021, 0.11001054, 0.02861683, 0.11002005, 0.02860233, 0.11002702,
                                             0.02858707, 0.11003126, 0.02857143, 0.11003269, 0.02855579, 0.11003126,
                                             0.02854053, 0.11002702, 0.02852603, 0.11002005, 0.02851265, 0.11001054,
                                             0.02850072, 0.10999872, 0.02849053, 0.10998489, 0.02848233, 0.10996937,
                                             0.02847632, 0.10995255, 0.02847266, 0.10993486, 0.02847143, 0.10991672,
                                             0.02847266, 0.10989857, 0.02847632, 0.10988088, 0.02848233, 0.10986406,
                                             0.02849053, 0.10984855, 0.02850072, 0.10983471, 0.02851265, 0.10982289,
                                             0.02852603, 0.10981338, 0.02854053, 0.10980641, 0.02855579, 0.10980217,
                                             0.02857143, 0.10980074, 0.02858707, 0.10980217, 0.02860233, 0.10980641,
                                             0.02861683, 0.10981338, 0.02863021, 0.10982289, 0.02864214, 0.10983471,
                                             0.02865233, 0.10984855, 0.02866053, 0.10986406, 0.02866653, 0.10988088,
                                             0.02867020, 0.10989857, 0.02867143, 0.10991672,
                                             54.17611241};

static double                     tanvec[] = { 0.000000e+00,
                                               1.444861e-01,  2.414423e-10,  1.427073e-01,  2.621384e-02,  1.374145e-01,  5.178220e-02,
                                               1.287381e-01,  7.607552e-02,  1.168917e-01,  9.849560e-02,  1.021671e-01,  1.184904e-01,
                                               8.492681e-02,  1.355676e-01,  6.559533e-02,  1.493066e-01,  4.464867e-02,  1.593692e-01,
                                               2.260261e-02,  1.655077e-01, -2.800177e-10,  1.675707e-01, -2.260261e-02,  1.655077e-01,
                                              -4.464867e-02,  1.593692e-01, -6.559533e-02,  1.493066e-01, -8.492681e-02,  1.355676e-01,
                                              -1.021671e-01,  1.184904e-01, -1.168917e-01,  9.849560e-02, -1.287381e-01,  7.607552e-02,
                                              -1.374145e-01,  5.178220e-02, -1.427073e-01,  2.621384e-02, -1.444861e-01, -2.414423e-10,
                                              -1.427073e-01, -2.621384e-02, -1.374145e-01, -5.178220e-02, -1.287381e-01, -7.607552e-02,
                                              -1.168917e-01, -9.849560e-02, -1.021671e-01, -1.184904e-01, -8.492681e-02, -1.355676e-01,
                                              -6.559533e-02, -1.493066e-01, -4.464867e-02, -1.593692e-01, -2.260261e-02, -1.655077e-01,
                                               2.800177e-10, -1.675707e-01,  2.260261e-02, -1.655077e-01,  4.464867e-02, -1.593692e-01,
                                               6.559533e-02, -1.493066e-01,  8.492681e-02, -1.355676e-01,  1.021671e-01, -1.184904e-01,
                                               1.168917e-01, -9.849560e-02,  1.287381e-01, -7.607552e-02,  1.374145e-01, -5.178220e-02,
                                               1.427073e-01, -2.621384e-02,  1.444861e-01,  2.414423e-10,
                                               0.000000e+00};

static double                     CDupoldp[] = {-1.938025e-09,  1.159770e+00,
                                                -1.564345e-01,  1.145492e+00,
                                                -3.090170e-01,  1.103007e+00,
                                                -4.539905e-01,  1.033363e+00,
                                                -5.877853e-01,  9.382740e-01,
                                                -7.071068e-01,  8.200815e-01,
                                                -8.090170e-01,  6.816960e-01,
                                                -8.910065e-01,  5.265248e-01,
                                                -9.510565e-01,  3.583888e-01,
                                                -9.876883e-01,  1.814281e-01,
                                                -1.000000e+00, -1.671042e-09,
                                                -9.876883e-01, -1.814281e-01,
                                                -9.510565e-01, -3.583888e-01,
                                                -8.910065e-01, -5.265248e-01,
                                                -8.090170e-01, -6.816960e-01,
                                                -7.071068e-01, -8.200815e-01,
                                                -5.877853e-01, -9.382740e-01,
                                                -4.539905e-01, -1.033363e+00,
                                                -3.090170e-01, -1.103007e+00,
                                                -1.564345e-01, -1.145492e+00,
                                                 1.938025e-09, -1.159770e+00,
                                                 1.564345e-01, -1.145492e+00,
                                                 3.090170e-01, -1.103007e+00,
                                                 4.539905e-01, -1.033363e+00,
                                                 5.877853e-01, -9.382740e-01,
                                                 7.071068e-01, -8.200815e-01,
                                                 8.090170e-01, -6.816960e-01,
                                                 8.910065e-01, -5.265248e-01,
                                                 9.510565e-01, -3.583888e-01,
                                                 9.876883e-01, -1.814281e-01,
                                                 1.000000e+00,  1.671042e-09,
                                                 9.876883e-01,  1.814281e-01,
                                                 9.510565e-01,  3.583888e-01,
                                                 8.910065e-01,  5.265248e-01,
                                                 8.090170e-01,  6.816960e-01,
                                                 7.071068e-01,  8.200815e-01,
                                                 5.877853e-01,  9.382740e-01,
                                                 4.539905e-01,  1.033363e+00,
                                                 3.090170e-01,  1.103007e+00,
                                                 1.564345e-01,  1.145492e+00,
                                                -1.938025e-09,  1.159770e+00};

#elif (STARTPNT == 2)
// This is the starting point for the LC number 66 starting from the HP in rosenzweig
static double                     point[] = { 0.13186143,
                                              0.05269729, 0.09623171, 0.05459896, 0.09957497, 0.05563825, 0.10319565, 0.05582370, 0.10701612, 0.05520802, 0.11095244,
                                              0.05387226, 0.11491474, 0.05191043, 0.11880671, 0.04942278, 0.12252634, 0.04651215, 0.12596693, 0.04328328, 0.12901963,
                                              0.03984294, 0.13157730, 0.03629951, 0.13354174, 0.03276100, 0.13483060, 0.02933147, 0.13538493, 0.02610472, 0.13517799,
                                              0.02315737, 0.13421905, 0.02054403, 0.13255374, 0.01829543, 0.13025976, 0.01641908, 0.12743876, 0.01490553, 0.12420439,
                                              0.01373451, 0.12067351, 0.01288091, 0.11695871, 0.01232133, 0.11316217, 0.01203774, 0.10937423, 0.01202001, 0.10567329,
                                              0.01226715, 0.10212689, 0.01278846, 0.09879380, 0.01360326, 0.09572633, 0.01474045, 0.09297213, 0.01623671, 0.09057565,
                                              0.01813375, 0.08857932, 0.02047048, 0.08702356, 0.02327381, 0.08594554, 0.02654627, 0.08537791, 0.03024804, 0.08534577,
                                              0.03428362, 0.08586358, 0.03849582, 0.08693256, 0.04267469, 0.08853925, 0.04658020, 0.09065497, 0.04998315, 0.09323710,
                                              0.05269729, 0.09623171,
                                             56.05687471};

static double                     tanvec[] = { 0.0041731514,
                                               0.0079112818, -0.0038699427,  0.0084893444, -0.0030582381,  0.0087256846, -0.0021643494,
                                               0.0086526044, -0.0012026157,  0.0083180474, -0.0001886375,  0.0077730010,  0.0008599207,
                                               0.0070616052,  0.0019216841,  0.0062210097,  0.0029703525,  0.0052824604,  0.0039740028,
                                               0.0042743467,  0.0048951538,  0.0032261007,  0.0056915607,  0.0021704154,  0.0063199598,
                                               0.0011436418,  0.0067402310,  0.0001839419,  0.0069208225, -0.0006730296,  0.0068456886,
                                              -0.0014003249,  0.0065180197, -0.0019840155,  0.0059609693, -0.0024247646,  0.0052144258,
                                              -0.0027370140,  0.0043284633, -0.0029432930,  0.0033545527, -0.0030693343,  0.0023397827,
                                              -0.0031399510,  0.0013229137, -0.0031752747,  0.0003321937, -0.0031898289, -0.0006132053,
                                              -0.0031921242, -0.0015016426, -0.0031846404, -0.0023269750, -0.0031638662, -0.0030862033,
                                              -0.0031197608, -0.0037775897, -0.0030350980, -0.0043991819, -0.0028847890, -0.0049477566,
                                              -0.0026348831, -0.0054177347, -0.0022446391, -0.0058008628, -0.0016703991, -0.0060865848,
                                              -0.0008735477, -0.0062626946,  0.0001655520, -0.0063171188,  0.0014314375, -0.0062400275,
                                               0.0028627494, -0.0060259064,  0.0043543968, -0.0056744797,  0.0057725841, -0.0051910811,
                                               0.0069897273, -0.0045854854,  0.0079112818, -0.0038699427,
                                               0.9991182275};

static double                     CDupoldp[] = {  0.091738925,  0.126042827,
                                                  0.058325770,  0.138474570,
                                                  0.024137031,  0.147946308,
                                                 -0.008747712,  0.154208602,
                                                 -0.038935729,  0.157018478,
                                                 -0.065641196,  0.156135779,
                                                 -0.088481775,  0.151319226,
                                                 -0.107296668,  0.142357701,
                                                 -0.122004610,  0.129113446,
                                                 -0.132530838,  0.111583174,
                                                 -0.138794223,  0.089967827,
                                                 -0.140758543,  0.064740605,
                                                 -0.138504654,  0.036691293,
                                                 -0.132311764,  0.006927776,
                                                 -0.122717673, -0.023195487,
                                                 -0.110499456, -0.052224153,
                                                 -0.096585910, -0.078787514,
                                                 -0.081916728, -0.101775174,
                                                 -0.067295261, -0.120468878,
                                                 -0.053285841, -0.134555022,
                                                 -0.040184665, -0.144076471,
                                                 -0.028044899, -0.149336917,
                                                 -0.016730133, -0.150773661,
                                                 -0.005976690, -0.148864707,
                                                  0.004556817, -0.144061514,
                                                  0.015253578, -0.136751270,
                                                  0.026513063, -0.127236626,
                                                  0.038725097, -0.115739466,
                                                  0.052235875, -0.102411731,
                                                  0.067291800, -0.087356699,
                                                  0.083953231, -0.070652176,
                                                  0.101960971, -0.052399053,
                                                  0.120590091, -0.032753056,
                                                  0.138524995, -0.011953075,
                                                  0.153827178,  0.009640996,
                                                  0.164129642,  0.031554341,
                                                  0.167093911,  0.053224164,
                                                  0.161052655,  0.074052128,
                                                  0.145588275,  0.093452477,
                                                  0.121746728,  0.110921036,
                                                  0.091738925,  0.126042827};
#endif

static double                     CDwi[] = {0.07777778, 0.35555556, 0.13333333, 0.35555556, 0.07777778};
static double                     CDwt[] = { 0.52520833,  0.8078098, -0.508159,  0.2143279, -0.03918701,
                                            -0.04082282,  0.6735161,  0.475506, -0.1283069,  0.02010762,
                                             0.02010762, -0.1283069,  0.475506,  0.6735161, -0.04082282,
                                            -0.03918701,  0.2143279, -0.508159,  0.8078098,  0.52520833};

static double                     CDwpvec[] = {-5.4645619,  7.703378, -3.211852,  1.171819, -0.1987826,
                                               -0.1119474, -4.530089,  5.542103, -1.046150,  0.1460829,
                                               -0.1460829,  1.046150, -5.542103,  4.530089,  0.1119474,
                                                0.1987826, -1.171819,  3.211852, -7.703378,  5.4645619};

static double                     CDwpdt[] = {-5.4645619,  0.0000000, -0.1119474,  0.0000000, -0.1460829,  0.0000000,  0.1987826,  0.0000000,
                                               0.0000000, -5.4645619,  0.0000000, -0.1119474,  0.0000000, -0.1460829,  0.0000000,  0.1987826,
                                               7.7033780,  0.0000000, -4.5300890,  0.0000000,  1.0461500,  0.0000000, -1.1718190,  0.0000000,
                                               0.0000000,  7.7033780,  0.0000000, -4.5300890,  0.0000000,  1.0461500,  0.0000000, -1.1718190,
                                              -3.2118520,  0.0000000,  5.5421030,  0.0000000, -5.5421030,  0.0000000,  3.2118520,  0.0000000,
                                               0.0000000, -3.2118520,  0.0000000,  5.5421030,  0.0000000, -5.5421030,  0.0000000,  3.2118520,
                                               1.1718190,  0.0000000, -1.0461500,  0.0000000,  4.5300890,  0.0000000, -7.7033780,  0.0000000,
                                               0.0000000,  1.1718190,  0.0000000, -1.0461500,  0.0000000,  4.5300890,  0.0000000, -7.7033780,
                                              -0.1987826,  0.0000000,  0.1460829,  0.0000000,  0.1119474,  0.0000000,  5.4645619,  0.0000000,
                                               0.0000000, -0.1987826,  0.0000000,  0.1460829,  0.0000000,  0.1119474,  0.0000000,  5.4645619};

#else

static int                        pntDim, sysOdeDim = 0, freeparsdim;
static int                        CurveType;

static SEXP                       R_UserFunc;
static SEXP                       R_PntVals;
static SEXP                       R_FixedPars;

static long int                   N_Protected;

static int                        glorder, ninterval, CDfinemeshdim;
static double                     *CDupoldp = NULL, *CDwi = NULL, *CDwt = NULL, *CDwpvec = NULL, *CDwpdt = NULL;

#endif

static int                        dGlobalMemDim = 0;
static double                     *dGlobalMem = NULL;
static double                     *state0, *rhs, *ups, *ficdmat, *ficd, *xp, *t, *blockjac, *partjac, *ic;

static int                        state0Dim, rhsDim, upsRows, upsCols, ficdmatRows, ficdmatCols, ficdDim, xpRows, xpCols,
                                  tRows, tCols, blockjacRows, blockjacCols, partjacRows, partjacCols, icDim;

/*
 *====================================================================================================================================
 *  Implementation of problem specification routines
 *====================================================================================================================================
 */

#if (defined(CMDLINEDEBUG))

int EQsystem(double *argument, double *result)

{
  double  K = argument[0];
  double  R = argument[1];
  double  C = argument[2];

  double  r = 0.5, a = 5.0, h = 3.0, eps = 0.5, mu = 0.05;

  result[0] = r * R * (1 - R/K) - a * R * C/(1 + a * h *R);
  result[1] = eps * a * R * C/(1 + a * h * R) - mu * C;

  return SUCCES;
}

#else

int EQsystem(double *argument, double *result)

{
  SEXP                            Result;                                           // Protected
  SEXP                            R_pnt;                                            // Unprotected

  //******************************************************************************
  // Map current estimate of solution to R variables

  memcpy(REAL(R_PntVals), argument, (freeparsdim + sysOdeDim) * sizeof(double));

  R_pnt = LCONS(R_UserFunc, LCONS(R_NilValue, LCONS(R_PntVals, LCONS(R_FixedPars, R_NilValue))));
  N_Protected++;

  Result = PROTECT(R_forceAndCall(R_pnt, 3, R_GlobalEnv));

  if (!sysOdeDim) sysOdeDim = length(coerceVector(VECTOR_ELT(Result, 0), REALSXP));

  memcpy(result, REAL(coerceVector(VECTOR_ELT(Result, 0), REALSXP)), sysOdeDim * sizeof(double));
  UNPROTECT(1);
  N_Protected--;

  return SUCCES;
}

#endif

/*==================================================================================================================================*/

static inline double dot(int dim, double *x, double *y)
{
  double                          ss = 0;

  for (int ii = 0; ii < dim; ii++) ss += x[ii] * y[ii];

   return ss;
}

static inline void bialt2AI(int matdim, double *A, double *result)
{
  double                          *dblpnt = result;

  for (int p = 1; p < matdim; p++)
    for (int q = 0; q < p; q++)
      for (int r = 1; r < matdim; r++)
        for (int s = 0; s < r; s++)
          {
            if (r == q) *dblpnt = -A[p*matdim + s];
            else if ((r != p) && (s == q)) *dblpnt = A[p*matdim + r];
            else if ((r == p) && (s == q)) *dblpnt = A[p*matdim + p] + A[q*matdim + q];
            else if ((r == p) && (s != q)) *dblpnt = A[q*matdim + s];
            else if (s == p) *dblpnt = -A[q*matdim + r];
            else *dblpnt = 0.0;
            dblpnt++;
          }
  return;
}

// Computes the matrix product A x B  (A : (rowsA * cArB) matrix; B : (cArB * colsB) matrix)
static inline void matXmat(int rowsA, int cArB, int colsB, double *A, double *B, double *AxB)
{
  memset(AxB, 0, (rowsA * colsB) * sizeof(double));
  for (int ii = 0; ii < rowsA; ii++)
    for (int jj = 0; jj < colsB; jj++)
      for (int kk = 0; kk < cArB; kk++)
        AxB[ii * colsB + jj] += A[ii * cArB + kk] * B[kk * colsB + jj];

   return;
}

// Multiplies the vector V with the constant c
static inline void conXvec(int dim, double c, double *V)
{
  for (int ii = 0; ii < dim; ii++) V[ii] *= c;

   return;
}

/*==================================================================================================================================*/

int BPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), int method, double *res)

/*
 * BPcondition - Routine computes the additional conditions determining the location of
 *               a branching point, see the Matcont documentation (Branch point locator,
 *               page 36, eq. 41)
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that the
 *                        point consists of:
 *                          p:  the first bifurcation parameter
 *                          q:  the second bifurcation parameter
 *                          x:  the solution point
 *                          b:  an additional value
 *                          v:  the eigenvector (same dimension as x)
 *              y       : Pointer to an array containing the values of the unknowns
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              res     : Pointer to the result vector to adjust and/or compute. If NULL
 *                        the routine is used to initialize the additional variables for
 *                        the BP continuation
 */
{
  int                             resindx, retcode  = SUCCES;
  double                          *dBaseMem, *Jac;
  double                          eval, *evec;

  // Some of the following vectors are over-sized, but that's OK
  Jac = dBaseMem = calloc((pntdim * pntdim), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in BPcondition()");

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The resulting Jacobian equals the following matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2. Notice that all coefficients pertaining to yi are to be found
   * in ROW i + 2 (as opposed to column i). The Jacobian matrix is hence stored in column-wise
   * (fortran-style) form.
   */

  /*
   *  The extended system to solve for is:
   *
   *     F(x, p) + b*v   = 0
   *     (F_x(x, p))^T v = 0
   *     v^T F_p(x, p)   = 0
   *     v^T v - 1       = 0
   *
   *     with initial conditions b = 0 and v the eigenvector of the matrix (F_x(x, p))^T pertaining to
   *     the eigenvalue with the smallest norm. The unknowns are:
   *
   *     p:  the first bifurcation parameter
   *     q:  the second bifurcation parameter
   *     x:  the solution point
   *     b:  an additional value
   *     v:  the eigenvector (same dimension as x)
   *
   */
  // Adjust the base equations
  eval = y[freeparsdim + sysOdeDim];
  evec = y + freeparsdim + sysOdeDim + 1;

  //  F(x, p) + b*v   = 0
  for (resindx = 0; resindx < sysOdeDim; resindx++) res[resindx] += eval*evec[resindx];

  // (F_x(x, p))^T v = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, evec, Jac + (freeparsdim + ii) * sysOdeDim);

  // v^T F_p(x, p)   = 0 : Take the second parameter in the list of arguments
  res[resindx++] = dot(sysOdeDim, evec, Jac + sysOdeDim);

  // v^T v - 1       = 0
  res[resindx++] = dot(sysOdeDim, evec, evec) - 1;

  free(dBaseMem);

  return retcode;
}


/*==================================================================================================================================*/

int HPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), int method, double *res)

/*
 * HPcondition - Routine computes the additional conditions determining the location of
 *               a Hopf bifurcation point, see equations (10.40) and (10.41) on page 485-486
 *               in Kuznetsov (1996), Elements of applied bifurcation analysis.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the bifurcation parameters) plus the dimension of
 *                        the vector of state variables
 *              y       : Pointer to an array containing as first one or two element
 *                        the value of the free parameter(s) and as subsequent elements
 *                        the values of the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *              res     : Pointer to the result
 */
{
  const int                       bialtdim = sysOdeDim*(sysOdeDim - 1)/2;
  int                             retcode = SUCCES;
  double                          *dBaseMem, *Jac, *Jx, *JI;

  // Some of the following vectors are over-sized, but that's OK
  Jac = dBaseMem = calloc(pntdim*pntdim + sysOdeDim*sysOdeDim + bialtdim*bialtdim, sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in HPcondition()");

  Jx    = Jac + pntdim*pntdim;
  JI    = Jx  + sysOdeDim*sysOdeDim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The resulting Jacobian equals the following matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2. Notice that all coefficients pertaining to yi are to be found
   * in ROW i + 2 (as opposed to column i). The Jacobian matrix is hence stored in column-wise
   * (fortran-style) form. Below it is transposed into a row-wise (C-style) form.
   */

  // Extract the restricted Jacobian in transposed (row-wise, C-style) form
  for (int ii = 0; ii < sysOdeDim; ii++)
    for (int jj = 0; jj < sysOdeDim; jj++) Jx[ii*sysOdeDim + jj] = Jac[(freeparsdim + jj)*sysOdeDim + ii];

  // Now construct the bialternate matrix product of 2J and I
  bialt2AI(sysOdeDim, Jx, JI);

  retcode = Determinant(bialtdim, JI, res + sysOdeDim, NULL);
  if ((retcode != SUCCES) && (retcode != SINGULARITY))
    {
      ErrorMsg("Failed to compute determinant of bialternate matrix product in HPcondition()");
      free(dBaseMem);
      return retcode;
    }

  free(dBaseMem);

  return retcode;
}


/*==================================================================================================================================*/

int LPcondition(const int pntdim, double *y, int (*fnc)(double *, double *), const int method, double *res)

/*
 * LPcondition -  Routine computes the factor determining the location of a limit point, i.e.
 *                the parameter component of the tangent vector. This component always has
 *                index 0 in the vector of the solution point.
 *
 * Arguments -  pntdim  : The dimension of the argument vector 'y'. Notice that this
 *                        equals 2 (for the bifurcation parameters) plus the dimension of
 *                        the vector of state variables
 *              y       : Pointer to an array containing as first one or two element
 *                        the value of the free parameter(s) and as subsequent elements
 *                        the values of the state variables y.
 *              fnc     : Pointer to function specifying the system of
 *                        equations. The function must have a (double)
 *                        pointer as first argument, containing the point
 *                        in which to evaluate the system and a (double)
 *                        pointer as second argument, containing the
 *                        results after evaluation of the equations.
 *              method  : Method to use for differential computation: FORWARD or CENTRAL
 *
 *   Continuation of a limitpoint is carried out using the defining system
 *   (10.97) on page 515 of Kuznetsov (1998):
 *
 *    F(x, p)                 = 0
 *    F_x(x, p) q             = 0
 *    (F_x(x, p))^T p - eps*p = 0
 *    q^T q - 1               = 0
 *    p^T p - 1               = 0
 *
 *    with initial conditions b = 0 and v the eigenvector of the matrix
 *    (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
 *    The unknowns are:
 *
 *    p:   the bifurcation parameter
 *    x:   the solution point
 *    eps: an additional value
 *    q:   the eigenvector of F_x pertaining to the 0 eigenvalue
 *    p:   the eigenvector of (F_x)^T pertaining to the 0 eigenvalue
 *
 *   The advantage of this defining system is that detection of Bogdanov-Takens
 *   points and cusp point are straightforward
 */
{
  register int                    resindx;
  double                          *dBaseMem, *Jac, *Jx, *JxT, eps, *pv, *qv;

  Jac = dBaseMem = calloc((2 * pntdim * sysOdeDim), sizeof(double));
  if (!dBaseMem) return ErrorMsg("Memory allocation error in LPcondition()");

  Jx = Jac + pntdim * sysOdeDim;

  // Determine the Jacobian of the extended system (variable plus parameter dependence).
  Jacobian(pntdim, y, sysOdeDim, Jac, fnc, method);

  /*
   * The Jacobian equals the following (n+2)x(n) matrix of partial derivatives:
   *
   *           |dF1/dp1 ... dFn/dp1|
   *           |dF1/dp2 ... dFn/dp2|
   *           |dF1/dy1 ... dFn/dy1|
   *      Df = |   .           .   |
   *           |   .           .   |
   *           |   .           .   |
   *           |dF1/dyn ... dFn/dyn|
   *
   * In which n = pntdim-2 (i.e. equal to the number of state variables).
   * Notice that all coefficients pertaining to yi are to be found in ROW i + 2
   * (as opposed to column i). The matrix is hence stored in column-wise (fortran) style.
   */

  // Extract the restricted Jacobian in transposed (row-wise, C-style) form
  for (int ii = 0; ii < sysOdeDim; ii++)
    for (int jj = 0; jj < sysOdeDim; jj++) Jx[ii * sysOdeDim + jj] = Jac[(freeparsdim + jj)*sysOdeDim + ii];

  // The transposed and restricted Jacobian of the system is found at row 2..(sysOdeDim+2)
  JxT = Jac + freeparsdim * sysOdeDim;

  // Extract the additional variables from the state
  eps  = y[freeparsdim + sysOdeDim];
  qv   = y + freeparsdim +     sysOdeDim + 1;
  pv   = y + freeparsdim + 2 * sysOdeDim + 1;

  resindx = sysOdeDim;
  // Add the additional values
  // A qv = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, qv, Jx + ii * sysOdeDim);

  // A^T pv - eps*pv = 0
  for (int ii = 0; ii < sysOdeDim; ii++, resindx++) res[resindx] = dot(sysOdeDim, pv, JxT + ii * sysOdeDim) - eps * pv[ii];

  // <qv, qv) - 1 =
  res[resindx++] = dot(sysOdeDim, qv, qv) - 1;

  // <pv, pv) - 1 =
  res[resindx++] = dot(sysOdeDim, pv, pv) - 1;

  free(dBaseMem);

  return SUCCES;
}


/*==================================================================================================================================*/

static void ExtSystemLCblockjac(double *xp, double *state0, double lcperiod, double *partjacout)

// ExtSystemLCblockjac <- function(xp, state, parms, curveData, nopts = NULL) {

{
  int                             CDwpdtRows  = (sysOdeDim * (glorder + 1));
  int                             CDwpdtCols  = (sysOdeDim * (glorder));
  double                          *sysjac;

  sysjac = blockjac + sysOdeDim;

  /* The first term on the right-hand side of the blockjac assignment below amounts to a direct copy
   * of CDwpdt into rows 1..(1 + CDwpdtRows) of partjac
   *
   *   wploc = curveData$wp/dt;
   *   blockjac[range1,1+(1:(jaccol-2))] <- wploc[range1,] - state["LCperiod"]*fastkron(glorder, sysOdeDim, t(CDwt[,j]), sysjac)
   */
  memcpy(partjacout + partjacCols, CDwpdt, (CDwpdtRows * CDwpdtCols) * sizeof(double));

  // Evaluate function value on each collocation point
  for (int ii = 0; ii < glorder; ii++)
    {
      memcpy(state0 + 1, xp + ii * sysOdeDim, sysOdeDim * sizeof(double));
      memset(blockjac, 0, (blockjacRows * blockjacCols) * sizeof(double));
      Jacobian((freeparsdim + sysOdeDim), state0, sysOdeDim, blockjac, EQsystem, FORWARD);

      // Derivatives w.r.t. the bifurcation parameter in first row
      memcpy(partjacout + (ii * sysOdeDim), blockjac, sysOdeDim * sizeof(double));

      // blockjac[range1,1+(1:(jaccol-2))] <- wploc[range1,] - state["LCperiod"]*fastkron(glorder, sysOdeDim, t(CDwt[,j]), sysjac)
      for (int jj = 0; jj < (glorder + 1); jj++)
        for(int kk = 0; kk < sysOdeDim; kk++)
          for(int ll = 0; ll < sysOdeDim; ll++)
            {
              partjacout[(1 + jj * sysOdeDim + kk) * partjacCols + (ii * sysOdeDim) + ll]  -= lcperiod * CDwt[ii * (glorder + 1) + jj] * sysjac[kk * sysOdeDim + ll];
            }

      // Derivatives w.r.t. the bifurcation parameter in last row
      EQsystem(state0, partjacout + ((partjacRows - 1) * partjacCols) + (ii * sysOdeDim));
    }

  // Multiply the parameter derivatives with -state["LCperiod"]
  for (int jj = 0; jj < partjacCols; jj++) partjacout[jj] *= -lcperiod;

  // Multiply the last row with -1
  for (int jj = 0, kk = ((partjacRows - 1) * partjacCols); jj < partjacCols; jj++) partjacout[kk + jj] *= -1;

  return;
}



static int ExtSystemLCjac(const int pntdim, double *y, const int fncdim, double *fulljac, int (*fnc)(double *, double *), int method)

{
  double                          dt = 1.0 / ((double)ninterval);

  memset(dGlobalMem, 0, dGlobalMemDim * sizeof(double));
  memcpy(state0, y, (freeparsdim + sysOdeDim) * sizeof(double));                    // Copy bifurcation parameter and state at t=0

  // In LCfunc.R the limit cycle period is added, but that seems unnecessary
  // state0[(freeparsdim + sysOdeDim)] = y[pntdim -1];
  double LCperiod = y[pntdim -1];

  // Extract a matrix with state variable values.
  // In R:  ups: sysOdeDim rows; finemesh  columns (fortran column-wise storage)
  // In C:  ups: finemesh  rows; sysOdeDim columns (row-wise storage)
  memcpy(ups, y + freeparsdim, (upsRows * upsCols) * sizeof(double));

  /* The jacobian has pntdim rows and fncdim columns and has the following layout:
   *
   *           |dF1/dp1   ... dFn/dp1  |
   *           |dF1/dy1   ... dFn/dy1  |
   *           |   .             .     |
   *      Df = |   .             .     |
   *           |   .             .     |
   *           |dF1/dyn-1 ... dFn/dyn-1|
   *           |dF1/dyn   ... dFn/dyn  |
   *
   * Here p1 is the bifurcation parameter and y1..yn-1 are the nodal points
   * at which the LC is approximated (dimension: finemeshdim * sysOdeDim)
   * and yn is the LC period.
   *
   * F1 .. Fn are the values returned by LCcondition, which include:
   *
   *  - (ninterval * glorder * sysOdeDim) values determined by the
   *    right-hand side of the model at the nodal points at which the LC
   *    is approximated
   *  - sysOdeDim values determined by the circular boundary condition
   *  - 1 value (the last) corresponding to the integral condition
   *
   */

  memset(fulljac, 0, pntdim * fncdim * sizeof(double));

  // Evaluate the block jacobians
  for (int ii = 0; ii < ninterval; ii++)
    {
      // Value of polynomial on each collocation point
      //     xp <- ups[, range1] %*% curveData$wt
      // Notice that CDwt (curveData$wt) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(xpRows, (glorder + 1), xpCols, CDwt, ups + ii * glorder * sysOdeDim, xp);

      /* partjac dimensions:
       *
       * Columns:   partjacCols = (sysOdeDim * glorder) (= blockrow in R)
       * Rows   :   partjacRows = (sysOdeDim * (glorder + 1) + 2) rows (= blockcol in R)
       *
       * Notice that just like for fulljac, in partjac all the derivative values
       * w.r.t. a particular variable should be in the same ROW.
       */

      // partjac <- ExtSystemLCblockjac(xp, state0, parms, curveData, nopts)
      memset(partjac, 0, (partjacRows * partjacCols) * sizeof(double));
     ExtSystemLCblockjac(xp, state0, LCperiod, partjac);

      // The derivative w.r.t. the bifurcation parameter
      memcpy(fulljac + ii * partjacCols, partjac, partjacCols * sizeof(double));

      //  fulljac[rowrange, 1 + colrange] <- partjac[(1:blockrow), 1 + (1:blockcol)]
      for (int jj = 0; jj < sysOdeDim * (glorder+1); jj++)
        memcpy(fulljac + (1 + ii * (sysOdeDim * glorder) + jj) * fncdim + ii * partjacCols, partjac + (1 + jj) * partjacCols, partjacCols * sizeof(double));

      // Last row of jacobian
      memcpy(fulljac + (pntdim - 1) * fncdim + ii * partjacCols, partjac + (partjacRows - 1) * partjacCols, partjacCols * sizeof(double));

      // Derivative of the integral constraint
      // p <- dt*(CDupoldp[,range1]*curveData$pwi)                                  //pwi is sysOdeDim stacked rows of wi
      // ic[1 + colrange] <- ic[1 + colrange] + p[1:blockcol]
      for (int jj =0; jj < (glorder + 1); jj++)
        for (int kk = 0; kk < sysOdeDim; kk++)
          ic[1 + (ii * glorder + jj) * sysOdeDim + kk] += dt * (CDupoldp[(ii * glorder + jj) * sysOdeDim + kk] * CDwi[jj]);
    }

  // Derivative of the boundary condition into the one-to-last column of fulljac
  for (int ii = 0; ii < sysOdeDim; ii++)
    {
      fulljac[(ii + 1) * fncdim + (fncdim - (sysOdeDim + 1) + ii)] = 1;                            // Derivative w.r.t. ups[kk] in row 2 and further
      fulljac[(pntDim - (sysOdeDim + 1) + ii) * fncdim + (fncdim - (sysOdeDim + 1) + ii)] = -1;   // Derivative w.r.t. ups[(CDfinemeshdim - 1) * sysOdeDim + kk] in end rows (but not the last)
    }

  // Derivative of the integral constraint into the last column of fulljac
  for (int ii = 0; ii < pntdim; ii++) fulljac[ii * fncdim + (fncdim - 1)] = ic[ii];

  return SUCCES;
}

static int LCcondition(double *argument, double *result)

// ExtSystemLC <- function(t, state, parms, curveData, nopts = NULL) {

{

  memset(dGlobalMem, 0, dGlobalMemDim * sizeof(double));
  memcpy(state0, argument, (freeparsdim + sysOdeDim) * sizeof(double));                    // Copy bifurcation parameter and state at t=0

  // In LCfunc.R the limit cycle period is added, but that seems unnecessary
  // state0[(freeparsdim + sysOdeDim)] = argument[pntDim -1];
  double LCperiod = argument[pntDim -1];

  // Extract a matrix with state variable values.
  // In R:  ups: sysOdeDim rows; finemesh  columns (fortran column-wise storage)
  // In C:  ups: finemesh  rows; sysOdeDim columns (row-wise storage)
  memcpy(ups, argument + freeparsdim, (upsRows * upsCols) * sizeof(double));

  // Compute the values for the integral condition: element-wise multiplication
  // of ups and upoldp and summing over the columns
  memset(ficd, 0, CDfinemeshdim * sizeof(double));
  for (int ii = 0; ii < upsRows; ii++)
    for (int jj = 0; jj < upsCols; jj++)
      ficd[ii] += ups[ii * upsCols + jj] * CDupoldp[ii * upsCols + jj];

  // Evaluate the derivatives at all the nodal points of the mesh
  double dt = 1.0 / ((double)ninterval);
  for (int ii = 0; ii < ninterval; ii++)
    {
      // Value of polynomial on each collocation point
      //     xp <- ups[, range1] %*% curveData$wt
      // Notice that CDwt (curveData$wt) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(xpRows, (glorder + 1), xpCols, CDwt, ups + ii * glorder * sysOdeDim, xp);

      // Derivative of polynomial on each collocation point
      // t  <- (ups[, range1] %*% CDwpvec)/dt
      // Notice that CDwpvec (curveData$wpvec) should be stored as a matrix with glorder
      // rows and (glorder + 1) columns
      matXmat(tRows, (glorder + 1), tCols, CDwpvec, ups + ii * glorder * sysOdeDim, t);
      conXvec((tRows  * tCols), ninterval, t);

      // Evaluate function value on each collocation point
      for (int jj = 0; jj < glorder; jj++)
        {
          memcpy(state0 + 1, xp + jj * sysOdeDim, sysOdeDim * sizeof(double));
          EQsystem(state0, rhs);

          for (int kk = 0; kk < sysOdeDim; kk++)
            result[(ii * glorder + jj) * sysOdeDim + kk] = t[jj * tCols + kk] - LCperiod * rhs[kk];
        }

      // Put the appropriate values of the integral condition in the matrix
      memcpy(ficdmat + ii * ficdmatCols, ficd + ii * glorder, ficdmatCols * sizeof(double));
    }

  // Ciruclar boundary conditions
  for (int kk = 0; kk < sysOdeDim; kk++)
    result[(CDfinemeshdim - 1) * sysOdeDim + kk] = ups[kk] - ups[(CDfinemeshdim - 1) * sysOdeDim + kk];

  //  Integral constraint
  result[CDfinemeshdim * sysOdeDim] = 0.0;
  for (int ii = 0; ii < ninterval; ii++)
    result[CDfinemeshdim * sysOdeDim] += dt * dot(ficdmatCols, CDwi, ficdmat + ii * ficdmatCols);

  return SUCCES;
}


/*==================================================================================================================================*/

int AllEquations(double *argument, double *result)

{
  int                             retval = SUCCES;

  // Compute the basic system of equations
  EQsystem(argument, result);

  //==================================================================================================================================
  // Add the final value in case of BP, HP or LP continuation

  if (CurveType == BP)
    retval = BPcondition(pntDim, argument, EQsystem, CENTRAL, result);
  else if (CurveType == HP)
    retval = HPcondition(pntDim, argument, EQsystem, CENTRAL, result);
  else if (CurveType == LP)
    retval = LPcondition(pntDim, argument, EQsystem, CENTRAL, result);

  return retval;
}


/*==================================================================================================================================*/

#if (defined(CMDLINEDEBUG))

int main(int argc, char **argv)
{
  int matsize = ((sysOdeDim * (glorder + 1)) * (sysOdeDim * glorder));

  // Elements of curveData$wp are only used in ExtSystemLCblockjac and
  // all divided by dt = 1.0 / ninterval. Take care of that multiplication here
  for (int ii = 0; ii < matsize; ii++) CDwpdt[ii] *= ninterval;

  // Allocate the global memory
  state0Dim    = freeparsdim + sysOdeDim + 1;
  rhsDim       = sysOdeDim;
  upsRows      = CDfinemeshdim;
  upsCols      = sysOdeDim;
  ficdmatRows  = ninterval;
  ficdmatCols  = glorder+1;
  ficdDim      = CDfinemeshdim;
  xpRows       = glorder;
  xpCols       = sysOdeDim;
  tRows        = glorder;
  tCols        = sysOdeDim;
  blockjacRows = freeparsdim + sysOdeDim;
  blockjacCols = sysOdeDim;
  partjacRows  = sysOdeDim * (glorder+1) + freeparsdim + 1;
  partjacCols  = sysOdeDim * glorder;
  icDim        = pntDim;

  dGlobalMemDim  = 0;
  dGlobalMemDim += state0Dim;
  dGlobalMemDim += rhsDim;
  dGlobalMemDim += upsRows      * upsCols;
  dGlobalMemDim += ficdmatRows  * ficdmatCols;
  dGlobalMemDim += ficdDim;
  dGlobalMemDim += xpRows       * xpCols;
  dGlobalMemDim += tRows        * tCols;
  dGlobalMemDim += blockjacRows * blockjacCols;
  dGlobalMemDim += partjacRows  * partjacCols;
  dGlobalMemDim += icDim;

  state0 = dGlobalMem = calloc(dGlobalMemDim, sizeof(double));
  if (!dGlobalMem)
    {
      ErrorMsg("Memory allocation error in deBif()");
      return 1;
    }

  rhs       = state0    + state0Dim;
  ups       = rhs       + rhsDim;
  ficdmat   = ups       + upsRows      * upsCols;
  ficd      = ficdmat   + ficdmatRows  * ficdmatCols;
  xp        = ficd      + ficdDim;
  t         = xp        + xpRows       * xpCols;
  blockjac  = t         + tRows        * tCols;
  partjac   = blockjac  + blockjacRows * blockjacCols;
  ic        = partjac   + partjacRows  * partjacCols;

  Jacobian_Step = 1.0E-4;
  int nIter = MaxIter;

  int retcode = FindPoint(pntDim, freeparsdim, point, tanvec, RTol, ATol, CTol, MaxIter, &nIter, LCcondition, ExtSystemLCjac);

  if (retcode == SUCCES)
    ErrorMsg("FindPoint: Succes!");
  else
    ErrorMsg("FindPoint: Error!");

  double                          fulljac[pntDim * pntDim];

  retcode = TangentVec(pntDim, point, fulljac, tanvec, LCcondition, ExtSystemLCjac, NULL);

  if (retcode == SUCCES)
    ErrorMsg("TangentVec: Succes!");
  else
    ErrorMsg("TangentVec: Error!");

  for (int ii = 0; ii < pntDim; ii++)
    printf("%E, ", tanvec[ii]);
  printf("\n");

  ExtSystemLCjac(pntDim, point, pntDim - 1, fulljac, EQsystem, FORWARD);

  if (retcode == SUCCES)
    ErrorMsg("ExtSystemLCjac: Succes!");
  else
    ErrorMsg("ExtSystemLCjac: Error!");

  free(dGlobalMem);

  return 0;
}

#else

/*==================================================================================================================================*/

SEXP deBif(SEXP curveType, SEXP userFunc, SEXP initVals, SEXP fixedParVals, SEXP tanVec,
           SEXP rTol, SEXP aTol, SEXP cTol, SEXP jacStep, SEXP maxIter, SEXP glOrder, SEXP nInterval, SEXP cData)

{
  int                             MaxIter, retcode = FAILURE, listel = 0, nIter = 10;
  double                          RTol, ATol, CTol;
  double                          *dBaseMem, *point, *tanvec, *JacExport, *dblPnt;
  SEXP                            outputList = R_NilValue, nms, outputListEl[3], R_VarNames;
  SEXP                            R_VarNamesLC, R_cDataNames;

  N_Protected = 0L;
  glorder = ninterval = CDfinemeshdim = 0;

  //============================== Process the curve type argument ===================================================================

  if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "BP"))      CurveType = BP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "EQ")) CurveType = EQ;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "HP")) CurveType = HP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "LP")) CurveType = LP;
  else if (!strcmp(CHAR(STRING_ELT(curveType, 0)), "LC")) CurveType = LC;

  if ((CurveType == EQ) || (CurveType == LC)) freeparsdim = 1;
  else freeparsdim = 2;

  //============================== Process the user funcion argument =================================================================

  R_UserFunc = PROTECT(duplicate(userFunc));
  N_Protected++;

  //============================== Process the numerical option argument =============================================================

  RTol          = asReal(rTol);
  ATol          = asReal(aTol);
  CTol          = asReal(cTol);
  Jacobian_Step = asReal(jacStep);
  MaxIter       = asInteger(maxIter);
  glorder       = asInteger(glOrder);
  ninterval     = asInteger(nInterval);

  if (CurveType == LC)
    {
      char    optname[MAX_STR_LEN];
      int     ncols = length(cData);

      // Get the dimensions from curveData
      R_cDataNames = PROTECT(getAttrib(cData, R_NamesSymbol));
      N_Protected++;

      for (int ii = 0, nset = 0; ii < ncols; ii++)
        {
          strcpy(optname, CHAR(STRING_ELT(R_cDataNames, ii)));
          if (strcmp(optname, "finemeshdim") == 0)
            {
              CDfinemeshdim = asInteger(VECTOR_ELT(cData, ii));
              nset++;
            }
          else if (strcmp(optname, "statedim") == 0)
            {
              sysOdeDim = asInteger(VECTOR_ELT(cData, ii));
              nset++;
            }
          if (nset == 2) break;
        }
    }

  //============================== Process the initial point argument ================================================================

  pntDim = length(initVals);

  point = dBaseMem = calloc((2 * pntDim + (pntDim * pntDim)), sizeof(double));      // point, tanvec, JacExport
  if (!dBaseMem)
    {
      ErrorMsg("Memory allocation error in deBif()");
      UNPROTECT(N_Protected);

      return outputList;
    }

  tanvec    = point     + pntDim;
  JacExport = tanvec    + pntDim;

  memcpy(point, REAL(initVals), pntDim * sizeof(double));

  R_VarNames = PROTECT(getAttrib(initVals, R_NamesSymbol));
  N_Protected++;

  if (CurveType == LC)
    {
      char    optname[MAX_STR_LEN];
      int     ncols = length(cData);

      // Allocate the global memory
      state0Dim    = freeparsdim + sysOdeDim + 1;
      rhsDim       = sysOdeDim;
      upsRows      = CDfinemeshdim;
      upsCols      = sysOdeDim;
      ficdmatRows  = ninterval;
      ficdmatCols  = glorder+1;
      ficdDim      = CDfinemeshdim;
      xpRows       = glorder;
      xpCols       = sysOdeDim;
      tRows        = glorder;
      tCols        = sysOdeDim;
      blockjacRows = freeparsdim + sysOdeDim;
      blockjacCols = sysOdeDim;
      partjacRows  = sysOdeDim * (glorder+1) + freeparsdim + 1;
      partjacCols  = sysOdeDim * glorder;
      icDim        = pntDim;

      dGlobalMemDim  = 0;
      dGlobalMemDim += state0Dim;
      dGlobalMemDim += rhsDim;
      dGlobalMemDim += upsRows      * upsCols;
      dGlobalMemDim += ficdmatRows  * ficdmatCols;
      dGlobalMemDim += ficdDim;
      dGlobalMemDim += xpRows       * xpCols;
      dGlobalMemDim += tRows        * tCols;
      dGlobalMemDim += blockjacRows * blockjacCols;
      dGlobalMemDim += partjacRows  * partjacCols;
      dGlobalMemDim += icDim;

      int dGlobalMemDim2 = 0;
      dGlobalMemDim2 += (pntDim - 2);                                                // CDupoldp (exclude par & period)
      dGlobalMemDim2 += (glorder + 1);                                               // CDwi
      dGlobalMemDim2 += (glorder)                   * (glorder + 1);                 // CDwt
      dGlobalMemDim2 += (glorder)                   * (glorder + 1);                 // CDwpvec
      dGlobalMemDim2 += (sysOdeDim * (glorder + 1)) * (sysOdeDim * glorder);         // CDwp

      state0 = dGlobalMem = calloc(dGlobalMemDim + dGlobalMemDim2, sizeof(double));
      if (!dGlobalMem)
        {
          ErrorMsg("Memory allocation error in deBif()");
          free(dBaseMem);

          return outputList;
        }

      rhs       = state0    + state0Dim;
      ups       = rhs       + rhsDim;
      ficdmat   = ups       + upsRows      * upsCols;
      ficd      = ficdmat   + ficdmatRows  * ficdmatCols;
      xp        = ficd      + ficdDim;
      t         = xp        + xpRows       * xpCols;
      blockjac  = t         + tRows        * tCols;
      partjac   = blockjac  + blockjacRows * blockjacCols;
      ic        = partjac   + partjacRows  * partjacCols;

      CDupoldp  = ic        + icDim;
      CDwi      = CDupoldp  + (pntDim - 2);
      CDwt      = CDwi      + (glorder + 1);
      CDwpvec   = CDwt      + ((glorder) * (glorder + 1));
      CDwpdt    = CDwpvec   + ((glorder) * (glorder + 1));

      for (int ii = 0, nset = 0; ii < ncols; ii++)
        {
          strcpy(optname, CHAR(STRING_ELT(R_cDataNames, ii)));
          if (strcmp(optname, "upoldp") == 0)
            {
              memcpy(CDupoldp, REAL(VECTOR_ELT(cData, ii)), (pntDim - 2) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wi") == 0)
            {
              memcpy(CDwi, REAL(VECTOR_ELT(cData, ii)), (glorder + 1) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wt") == 0)
            {
              memcpy(CDwt, REAL(VECTOR_ELT(cData, ii)), ((glorder) * (glorder + 1)) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wpvec") == 0)
            {
              memcpy(CDwpvec, REAL(VECTOR_ELT(cData, ii)), ((glorder) * (glorder + 1)) * sizeof(double));
              nset++;
            }
          else if (strcmp(optname, "wp") == 0)
            {
              int matsize = ((sysOdeDim * (glorder + 1)) * (sysOdeDim * glorder));
              memcpy(CDwpdt, REAL(VECTOR_ELT(cData, ii)), matsize * sizeof(double));

              // Elements of curveData$wp are only used in ExtSystemLCblockjac and
              // all divided by dt = 1.0 / ninterval. Take care of that multiplication here
              for (int ii = 0; ii < matsize; ii++) CDwpdt[ii] *= ninterval;
              nset++;
            }
          if (nset == 5) break;
        }

      R_PntVals = PROTECT(allocVector(REALSXP, (freeparsdim + sysOdeDim + 1)));
      N_Protected++;

      // create names for single state vector
      R_VarNamesLC = PROTECT(allocVector(STRSXP, (freeparsdim + sysOdeDim + 1)));
      N_Protected++;

      for (int ii = 0; ii < (freeparsdim + sysOdeDim); ii++)
        SET_STRING_ELT(R_VarNamesLC, ii, STRING_ELT(R_VarNames, ii));

      SET_STRING_ELT(R_VarNamesLC, (freeparsdim + sysOdeDim), mkChar("LCperiod"));

      setAttrib(R_PntVals, R_NamesSymbol, R_VarNamesLC);

      UNPROTECT(1);
      N_Protected--;
    }
  else
    {
      R_PntVals = PROTECT(duplicate(initVals));
      N_Protected++;
    }

  //================================ Process the parameters argument =================================================================

  R_FixedPars = PROTECT(duplicate(fixedParVals));
  N_Protected++;

  //================================ Process the tangent argument ====================================================================

  memcpy(tanvec, REAL(tanVec), pntDim * sizeof(double));

  //============================= Compute the solution point =========================================================================

  retcode = FAILURE;
  if (CurveType == LC)
    {
        retcode = FindPoint(pntDim, freeparsdim, point, tanvec, RTol, ATol, CTol, MaxIter, &nIter, LCcondition, ExtSystemLCjac);
    }
  else
    {
      retcode = FindPoint(pntDim, freeparsdim, point, tanvec, RTol, ATol, CTol, MaxIter, &nIter, AllEquations, Jacobian);
    }

  //============================= Return the solution point =========================================================================

  if (retcode == SUCCES)
    {
      PROTECT(outputListEl[listel] = allocVector(REALSXP, pntDim));
      N_Protected++;

      dblPnt = REAL(outputListEl[listel]);
      memcpy(dblPnt, point, pntDim * sizeof(double));

      setAttrib(outputListEl[listel], R_NamesSymbol, R_VarNames);
      listel++;

      PROTECT(outputListEl[listel] = allocVector(INTSXP, 1));
      N_Protected++;

      // return the number of iterations
      INTEGER(outputListEl[listel])[0] = nIter;
      listel++;

      if (CurveType == LC)
        retcode = TangentVec(pntDim, point, JacExport, tanvec, LCcondition, ExtSystemLCjac, NULL);
      else
        retcode = TangentVec(pntDim, point, JacExport, tanvec, AllEquations, Jacobian, NULL);

      if (retcode == SUCCES)
        {
          PROTECT(outputListEl[listel] = allocVector(REALSXP, pntDim));
          N_Protected++;

          dblPnt = REAL(outputListEl[listel]);
          memcpy(dblPnt, tanvec, pntDim * sizeof(double));

          setAttrib(outputListEl[listel], R_NamesSymbol, R_VarNames);
          listel++;

          PROTECT(outputListEl[listel] = allocMatrix(REALSXP, pntDim, pntDim));
          N_Protected++;

          dblPnt = REAL(outputListEl[listel]);
          // memcpy(dblPnt, JacExport, pntDim * sysOdeDim * sizeof(double));
          for (int ii = 0; ii < pntDim; ii++)
            {
              memcpy((dblPnt + ii * pntDim), (JacExport + ii * (pntDim - 1)), (pntDim - 1) * sizeof(double));
              dblPnt[ii * pntDim + pntDim - 1] = tanvec[ii];
            }
          listel++;
        }

      outputList = PROTECT(allocVector(VECSXP, listel));
      N_Protected++;

      // create names
      nms = PROTECT(allocVector(STRSXP, listel));
      N_Protected++;

      switch (listel)
        {
          case 4:
            SET_STRING_ELT(nms, 2, mkChar("Jacobian"));
          case 3:
            SET_STRING_ELT(nms, 1, mkChar("tangent"));
          case 2:
            SET_STRING_ELT(nms, 1, mkChar("niter"));
          default:
            SET_STRING_ELT(nms, 0, mkChar("y"));
        }

      for (int ii = 0; ii < listel; ii++) SET_VECTOR_ELT(outputList, ii, outputListEl[ii]);
      setAttrib(outputList, R_NamesSymbol, nms);
    }

  UNPROTECT(N_Protected);

  if (dGlobalMem) free(dGlobalMem);
  free(dBaseMem);

  return outputList;
}

#endif
