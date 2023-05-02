// Included libraries
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>

/**************************************************************************************
// Adjoint Methods Programming Assignment Poisson Solver File Header
// Code for Poisson solver provided by Dr. Carl Olivier-Gooch, modified to work with 
// adjoint methods for error estimation.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef POISSON_H

#define POISSON_H

// Include user created header files here to avoid them being called multiple times 

// Define number of points to be used
#define IMAX 40
#define JMAX IMAX
#define MAX(a,b) ((a) > (b) ? (a) : (b)) // Return max of two numbers 

// Define problem parameters
const double dDX = 1./(double) IMAX; // dx
const double dDY = 1./(double) JMAX; // dy
const double dDX2 = 1./(double) (IMAX*IMAX); // dx^2
const double dDY2 = 1./(double) (JMAX*JMAX); // dy^2
const double dTol = 1.e-12; // Tolerance
const double dA = 1.; // Area for unit square

// Included functions
double dExact(const double dX, const double dY);
double dTInit(const int i, const int j);
void vApplyHomoBC(double a2dPrimal[IMAX+2][JMAX+2]);
double dPointGSUpdate(double a2dSoln[IMAX+2][JMAX+2], double a2dSource[IMAX+2][JMAX+2], void (*vApplyBC)(double [IMAX+2][JMAX+2]), const double dOmega);
double dPointJacUpdate(double a2dT[IMAX+2][JMAX+2], double a2dSource[IMAX+2][JMAX+2], void (*vApplyBC)(double [IMAX+2][JMAX+2]), const double dOmega);
double Primal_Source(const double dX, const double dY);
double Dual_Source(const double dX, const double dY);
void run_Poisson(void (*vApplyBC)(double [IMAX+2][JMAX+2]), double (*dSource)(const double dX, const double dY), double (&a2dSoln)[IMAX + 2][JMAX + 2], double (&a2dSource)[IMAX+2][JMAX+2]);
double calc_grad(double u_1, double u_2, double dDX);
double calc_L2Norm(double a2dVec[IMAX + 2][JMAX + 2]);
void calc_functional(double a2dSoln[IMAX+2][JMAX+2], double a2dSource[IMAX + 2][JMAX + 2], double (&a2dFunctional));
void storeVTKSolution(double a2dSoln[IMAX + 2][JMAX + 2], std::string fileName);
void computeLu(double Soln2nd[IMAX + 2][JMAX + 2], double Soln4th[IMAX + 4][JMAX + 4], double (&Lu)[IMAX + 2][JMAX + 2]);
void vApplyBC_uQ3(double a2dPrimal[IMAX+2][JMAX+2]);
void vApplyBC_vQ3(double a2dPrimal[IMAX+2][JMAX+2]);
double compute_h(const double dX);
double compute_e(const double dX);
// Included classes and objects 

#endif
