// Included libraries
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define number of points in X and Y 
#define IMAX 40
#define JMAX IMAX

// Return max of two numbers 
#define MAX(a,b) ((a) > (b) ? (a) : (b))

const double dDX = 1./(double) IMAX; // dx
const double dDY = 1./(double) JMAX; // dy
const double dDX2 = 1./(double) (IMAX*IMAX); // dx^2
const double dDY2 = 1./(double) (JMAX*JMAX); // dy^2
const double dTol = 1.e-12; // Tolerance
const double dA = 1.; // Area for unit square

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
  for (i = 1; i <= IMAX; i++) 
  {
    dX = (i-0.5)*dDX;
    a2dPrimal[i][0] = 2 * dExact(dX, 0) - a2dPrimal[i][1];
    a2dPrimal[i][JMAX+1] = 2 * dExact(dX, 1) - a2dPrimal[i][JMAX];
  }
   
  // Fixed temperature at x = 0 and x = 1
  for (j = 1; j <= JMAX; j++) 
  {
    dY = (j-0.5)*dDY;
    a2dPrimal[0][j] = 2 * dExact(0, dY) - a2dPrimal[1][j];
    a2dPrimal[IMAX+1][j] = 2 * dExact(1, dY) - a2dPrimal[IMAX][j];
  }
}

/* Point Gauss-Seidel update with SOR. <- SOR = successive over-relaxation  */
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
double dPointJacUpdate(double a2dT[IMAX+2][JMAX+2], 
		                   double a2dSource[IMAX+2][JMAX+2],
		                   void (*vApplyBC)(double [IMAX+2][JMAX+2]),
		                   const double dOmega)
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
  static const double dVal = -M_PI*M_PI*M_PI/4;
  return  dVal * sin(M_PI*dX) * sin(M_PI*dY);
}

int main()
{
  double a2dSoln[IMAX+2][JMAX+2];
  double a2dSource[IMAX+2][JMAX+2];
  int i, j;
  for (i = 0; i <= IMAX+1; i++) 
  {
    for (j = 0; j <= JMAX+1; j++)
     {
      double dX = (i-0.5) * dDX;
      double dY = (j-0.5) * dDY;
      a2dSoln[i][j] = 0;
      /* Let's try a bit of integration here, four-point Gauss quadrature... */
      a2dSource[i][j] = 0.25* (dSource(dX + dDX/sqrt(12), dY + dDY/sqrt(12)) +
	                             dSource(dX + dDX/sqrt(12), dY - dDY/sqrt(12)) +
	                             dSource(dX - dDX/sqrt(12), dY + dDY/sqrt(12)) +
	                             dSource(dX - dDX/sqrt(12), dY - dDY/sqrt(12)));
    }
  }
  vApplyBC(a2dSoln);
  {
    int iIter = 0;
    double dMaxChange;
    do {
      dMaxChange = dPointGSUpdate
	(a2dSoln, a2dSource, 
	 vApplyBC, 1.7 /*1.4 + 0.4*sin(iIter*0.3)*/);
	 /*      dMaxChange = dPointJacUpdate(a2dSoln, 1.);*/
      iIter ++;
      /* fprintf(stderr, "%4d %12.5G\n", iIter, dMaxChange); */
    } while (iIter < 50000 && dMaxChange > dTol);
  }

  {
    double dL2 = 0;
    double dAve = 0;
    for (i = 1; i <= IMAX; i++) {
      for (j = 1; j <= JMAX; j++) {
	double dExactCV = a2dSource[i][j] / (-2*M_PI*M_PI);
	double dErr = dExactCV - a2dSoln[i][j];
	dL2 += dErr*dErr;
	dAve += a2dSoln[i][j];
      }
    }
    dL2 = sqrt(dL2 / (IMAX*JMAX));
    dAve /= (IMAX*JMAX);
    fprintf(stderr, "L2 error norm: %12.5G\n", dL2);
    fprintf(stderr, "Average temp: %12.10G\n", dAve);
    fprintf(stderr, "Temp at center: %12.10G\n",
	    0.25*(a2dSoln[IMAX/2][JMAX/2] + a2dSoln[IMAX/2][JMAX/2+1] +
		  a2dSoln[IMAX/2+1][JMAX/2] + a2dSoln[IMAX/2+1][JMAX/2+1]));
  }

  return 0;
}
