#include "poisson.h"

/**************************************************************************************
// Adjoint Methods Programming Assignment Poisson Solver File Header
// Code for Poisson solver provided by Dr. Carl Olivier-Gooch, modified to work with 
// adjoint methods for error estimation.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Define exact solution
double dExact(const double dX, const double dY)
{
  return (sin(M_PI*dX) * sin(M_PI*dY));
}

// Define initial condition 
double dTInit(const int i, const int j)
{
  double dX = (i-0.5)*dDX;
  double dY = (j-0.5)*dDY;
  return 0;
}

// Apply boundary conditions
// 2 extra points for ghost cells at boundary
void vApplyBC(double a2dPrimal[IMAX+2][JMAX+2])
{
  int i, j;
  double dX;
  double dY = (j-0.5)*dDY;
  
  // Fixed temperature at y = 0 and y = 1
  // Top and bottom boundary
  for (i = 1; i <= IMAX; i++) 
  {
    dX = (i-0.5)*dDX;
    a2dPrimal[i][0] = 2 * dExact(dX, 0) - a2dPrimal[i][1]; // Bottom boundary
    a2dPrimal[i][JMAX+1] = 2 * dExact(dX, 1) - a2dPrimal[i][JMAX]; // Top boundary
  }
   
  // Fixed temperature at x = 0 and x = 1
  // Left and right boundary
  for (j = 1; j <= JMAX; j++) 
  {
    dY = (j-0.5)*dDY;
    a2dPrimal[0][j] = 2 * dExact(0, dY) - a2dPrimal[1][j]; // Left boundary
    a2dPrimal[IMAX+1][j] = 2 * dExact(1, dY) - a2dPrimal[IMAX][j]; // Right boundary
  }
}

/* Point Gauss-Seidel update with SOR. <- SOR = successive over-relaxation  */
// Used to solve system of linear equations. SOR just gives faster convergence
double dPointGSUpdate(double a2dSoln[IMAX+2][JMAX+2], double a2dSource[IMAX+2][JMAX+2], 
                      void (*vApplyBC)(double [IMAX+2][JMAX+2]), const double dOmega)
{
  double dMaxChange = 0;
  int i, j;
  
/*   for (j = JMAX; j >= 1; j --) { */
/*     for (i = IMAX; i >= 1; i --) { */
  for (i = 1; i <= IMAX; i ++) 
  {
    for (j = 1; j <= JMAX; j ++) 
    {
      double dDeltaT = (dDY2 * (a2dSoln[i+1][j] + a2dSoln[i-1][j]) +
			                 dDX2 * (a2dSoln[i][j+1] + a2dSoln[i][j-1]) -
		                   dDX2*dDY2 * a2dSource[i][j]) / 
	                     (2 * (dDX2 + dDY2)) - a2dSoln[i][j];

      a2dSoln[i][j] += dDeltaT * dOmega;
      // printf("%2d %2d %10f %10f\n", i, j, dDeltaT, dMaxChange);
      dMaxChange = MAX(dMaxChange, fabs(dDeltaT * dOmega));
    }
  }
  vApplyBC(a2dSoln);
  return (dMaxChange);
}      

/* Point Jacobi update with SOR. SOR = successive over-relaxation */
// Alternate method to solve system of linear equations? I dont think this is used.
double dPointJacUpdate(double a2dT[IMAX+2][JMAX+2], double a2dSource[IMAX+2][JMAX+2],
		                   void (*vApplyBC)(double [IMAX+2][JMAX+2]), const double dOmega)
{
  double dMaxChange = 0;
  double a2dTTemp[IMAX+2][JMAX+2];
  int i, j;

  for (j = JMAX+1; j >= 0; j --) 
  {
    for (i = IMAX+1; i >= 0; i --) 
    {
      a2dTTemp[i][j] = a2dT[i][j];
    }
  }

  for (j = JMAX; j >= 1; j --) 
  {
    for (i = IMAX; i >= 1; i --)
   {
      // for (i = 1; i <= IMAX; i ++) {
      // for (j = 1; j <= JMAX; j ++) {
      double dDeltaT = (dDY2 * (a2dTTemp[i+1][j] + a2dTTemp[i-1][j]) +
			                 dDX2 * (a2dTTemp[i][j+1] + a2dTTemp[i][j-1]) -
		                   dDX2*dDY2 * a2dSource[i][j]) / 
	                     (2 * (dDX2 + dDY2)) - a2dTTemp[i][j];
                       a2dT[i][j] += dDeltaT * dOmega;
      dMaxChange = MAX(dMaxChange, dDeltaT * dOmega);
    }
  }
  vApplyBC(a2dT);
  return (dMaxChange);
}      

double dSource(const double dX, const double dY)
{
  // static const double dVal = -M_PI*M_PI*M_PI/4; // Primal
  static const double dVal = (M_PI*M_PI*M_PI*M_PI*M_PI/2);
  //return  dVal * sin(M_PI*dX) * sin(M_PI*dY);
  return dVal * dX * (1 - dX) * dY * (1 - dY);
}


