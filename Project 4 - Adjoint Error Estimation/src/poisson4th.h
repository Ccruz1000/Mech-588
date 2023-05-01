/**************************************************************************************
// Adjoint Methods Programming Assignment 4th Order Poisson Solver File Header
// Code for Poisson solver provided by Dr. Carl Olivier-Gooch, modified to work with 
// adjoint methods for error estimation.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

#ifndef POISSON4TH_H

#define POISSON4TH_H

// Include user created header files here to avoid them being called multiple times 
#include "poisson.h"

// Included Functions
void vApplyHomoBC_4thOrder(double a2dPrimal[IMAX+4][JMAX+4]);
void run_Poisson_4thOrder(void (*vApplyBC)(double [IMAX+4][JMAX+4]), double (*dSource)(const double dX, const double dY), double (&a4dSoln)[IMAX + 4][JMAX + 4], double (&a2dSource)[IMAX+4][JMAX+4]);
double dPointGSUpdate_4thOrder(double a2dSoln[IMAX+4][JMAX+4], double a2dSource[IMAX+4][JMAX+4], void (*vApplyBC)(double [IMAX+4][JMAX+4]), const double dOmega);
void storeVTKSolution_4thOrder(double a2dSoln[IMAX + 4][JMAX + 4], std::string fileName);
double calc_L2Norm_4thOrder(double a2dVec[IMAX + 4][JMAX + 4]);

#endif