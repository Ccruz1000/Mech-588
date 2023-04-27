// Included header files
#include "poisson.h"

int main()
{
  double a2dSoln[IMAX+2][JMAX+2];
  double a2dSource[IMAX+2][JMAX+2];
  int i, j;
  // Loop through all cells
  for (i = 0; i <= IMAX+1; i++) 
  {
    for (j = 0; j <= JMAX+1; j++)
     {
      // dX and dY represent the current X and Y coordinate, at the 
      // left and bottom cell face. 
      double dX = (i-0.5) * dDX;
      double dY = (j-0.5) * dDY;
      a2dSoln[i][j] = 0; // Initialize solution to be 0
      /* Let's try a bit of integration here, four-point Gauss quadrature... */
      a2dSource[i][j] = 0.25* (dSource(dX + dDX/sqrt(12), dY + dDY/sqrt(12)) +
	                             dSource(dX + dDX/sqrt(12), dY - dDY/sqrt(12)) +
	                             dSource(dX - dDX/sqrt(12), dY + dDY/sqrt(12)) +
	                             dSource(dX - dDX/sqrt(12), dY - dDY/sqrt(12)));
    }
  }
  vApplyBC(a2dSoln); // Apply boundary conditions
  
  int iIter = 0;
  double dMaxChange; // Used to find maximum possible change in temperature, before solution is 
  // considered steady state.

  // Here we will calculate the max change of all the cells. We will then
  // increment our iteration. If we reach the maximum number of iterations, or
  // if the solution reaches steady state, then we will break from the loop.
  do 
  {
    dMaxChange = dPointGSUpdate(a2dSoln, a2dSource, vApplyBC, 1.7);
	 // vApplyBC, 1.7 /*1.4 + 0.4*sin(iIter*0.3)*/);
	 /*      dMaxChange = dPointJacUpdate(a2dSoln, 1.);*/
    iIter ++;
      /* fprintf(stderr, "%4d %12.5G\n", iIter, dMaxChange); */
  } 

  while (iIter < 50000 && dMaxChange > dTol);
  
  double dL2 = 0;
  double dAve = 0;

  // Loop over interior cells to calculate the exact solution, error and L2Norm
  for (i = 1; i <= IMAX; i++) 
  {
    for (j = 1; j <= JMAX; j++) 
    {
      double dExactCV = a2dSource[i][j] / (-2*M_PI*M_PI);
      double dErr = dExactCV - a2dSoln[i][j];
      	dL2 += dErr*dErr; // Used to compute the L2Norm of the vector
      	dAve += a2dSoln[i][j]; // Average solution temperature
    }
  }

  dL2 = sqrt(dL2 / (IMAX*JMAX));
  dAve /= (IMAX*JMAX);
  fprintf(stderr, "L2 error norm: %12.5G\n", dL2);
  fprintf(stderr, "Average temp: %12.10G\n", dAve);
  fprintf(stderr, "Temp at center: %12.10G\n",
          0.25*(a2dSoln[IMAX/2][JMAX/2] + a2dSoln[IMAX/2][JMAX/2+1] +
	        a2dSoln[IMAX/2+1][JMAX/2] + a2dSoln[IMAX/2+1][JMAX/2+1]));
  

  return 0;
}