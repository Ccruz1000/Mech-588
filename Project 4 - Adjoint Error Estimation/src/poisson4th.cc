// Included libraries and headers
#include "poisson4th.h"

/**************************************************************************************
// Adjoint Methods Programming Assignment 4th Order Poisson Solver File 
// Code for Poisson solver provided by Dr. Carl Olivier-Gooch, modified to work with 
// adjoint methods for error estimation.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Apply homogenous boundary conditions, used for 4th order accuracy 
// 4 extra points for ghost cells at boundary
void vApplyHomoBC_4thOrder(double a2dPrimal[IMAX+4][JMAX+4])
{
  int i, j;
  double dX;
  double dY = (j-0.5)*dDY;
  
  // Fixed temperature at y = 0 and y = 1
  // Top and bottom boundary
  for (i = 2; i <= IMAX + 2; i++) 
  {
    dX = (i-0.5)*dDX;
    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    a2dPrimal[i][0] = -a2dPrimal[i][3];
    a2dPrimal[i][1] = -a2dPrimal[i][2];
    a2dPrimal[i][JMAX + 2] = -a2dPrimal[i][JMAX + 1];
    a2dPrimal[i][JMAX + 3] = -a2dPrimal[i][JMAX];
  }
   
  // Fixed temperature at x = 0 and x = 1
  // Left and right boundary
  for (j = 2; j <= JMAX + 2; j++) 
  {
    dY = (j-0.5)*dDY;
    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    a2dPrimal[0][j] = -a2dPrimal[3][j];
    a2dPrimal[1][j] = -a2dPrimal[2][j];
    a2dPrimal[IMAX + 2][j] = -a2dPrimal[IMAX + 1][j];
    a2dPrimal[IMAX + 3][j] = -a2dPrimal[IMAX][j];
  }
}


/*
WHERE I AM FOR TOMORROW - 
Currently, we have question 1 working properly. I need to implement the 4th order solver that Carl showed us. As of right now, 
the 4th order functions work, but the 4th order scheme is not implemented. To ensure that the indices are correct, I will
run the 4th order functions with everything the same as the 2nd order. This way, I should get the same solution for both the 2nd order and the
4th order functions. If I compare the L2Norm for each of these, I should get the exact same answer. Once this works I can modify the GS solver 
for the 4th order function to use the scheme from the 4th order. I can verify that it is indeed 4th order by using grid resolution of the L2Error,
and ensuring a slope of 4. Once this works, I can solve the 2nd question and being the third :)
*/
void run_Poisson_4thOrder(void (*vApplyBC)(double [IMAX+4][JMAX+4]), double (*dSource)(const double dX, const double dY), double (&a2dSoln)[IMAX + 4][JMAX + 4], double (&a2dSource)[IMAX+4][JMAX+4])
{
  int i, j;
  // Loop through all cells
  for (i = 0; i <= IMAX+3; i++) 
  {
    for (j = 0; j <= JMAX+3; j++)
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
    dMaxChange = dPointGSUpdate_4thOrder(a2dSoln, a2dSource, vApplyBC, 1.7);
   // vApplyBC, 1.7 /*1.4 + 0.4*sin(iIter*0.3)*/);
   /*      dMaxChange = dPointJacUpdate(a2dSoln, 1.);*/
    iIter ++;
      /* fprintf(stderr, "%4d %12.5G\n", iIter, dMaxChange); */
  } 

  while (iIter < 50000 && dMaxChange > dTol);
  
  double dL2 = 0;
  double dAve = 0;

  // Loop over interior cells to calculate the exact solution, error and L2Norm
  for (i = 2; i <= IMAX + 1; i++) 
  {
    for (j = 2; j <= JMAX + 1; j++) 
    {
      double dExactCV = dExact((i - 0.5) * dDX, (j - 0.5) * dDY);
      double dErr = dExactCV - a2dSoln[i][j];
        dL2 += dErr*dErr; // Used to compute the L2Norm of the vector
        dAve += a2dSoln[i][j]; // Average solution temperature
    }
  }

  dL2 = sqrt(dL2 / (IMAX*JMAX));
  dAve /= (IMAX*JMAX);
  fprintf(stderr, "L2 error norm: %1.16e\n", dL2);
  // fprintf(stderr, "Average temp: %12.10G\n", dAve);
  // fprintf(stderr, "Temp at center: %12.10G\n",
  //         0.25*(a2dSoln[IMAX/2][JMAX/2] + a2dSoln[IMAX/2][JMAX/2+2] +
  //         a2dSoln[IMAX/2+2][JMAX/2] + a2dSoln[IMAX/2+2][JMAX/2+2])); 
}

double dPointGSUpdate_4thOrder(double a2dSoln[IMAX+4][JMAX+4], double a2dSource[IMAX+4][JMAX+4], void (*vApplyBC)(double [IMAX+4][JMAX+4]), const double dOmega)
{
  double dMaxChange = 0;
  int i, j;
  
  // Loop through interior cells
  for (i = 2; i <= IMAX + 1; i ++) 
  {
    for (j = 2; j <= JMAX + 1; j ++) 
    {
      double dDeltaT = (dDY2 * (a2dSoln[i+1][j] + a2dSoln[i-1][j]) +
			                 dDX2 * (a2dSoln[i][j+1] + a2dSoln[i][j-1]) -
		                   dDX2*dDY2 * a2dSource[i][j]) / 
	                     (2 * (dDX2 + dDY2)) - a2dSoln[i][j];

                       /* 
                       [(u_(i, j+1) + u_(i, j-1))/dx^2 + (u_(i+1,j) + u_(i-1,j))/dy^2 - u_(i, j) / (dx^2*dy^2)] / (2 * dx^2 * dy^2) - u_(i, j)
                       */ 

      a2dSoln[i][j] += dDeltaT * dOmega;
      // printf("%2d %2d %10f %10f\n", i, j, dDeltaT, dMaxChange);
      dMaxChange = MAX(dMaxChange, fabs(dDeltaT * dOmega));
    }
  }
  vApplyBC(a2dSoln);
  return (dMaxChange);
}


void storeVTKSolution_4thOrder(double a2dSoln[IMAX + 4][JMAX + 4], std::string fileName)
{
  const char * fileChar = fileName.c_str();
  FILE *vtkFile;
  vtkFile = fopen(fileChar, "w");
  unsigned long i, j;
  fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
  fprintf(vtkFile,"TITLE = \"Quad data\"\n");
  fprintf(vtkFile,"ASCII\n");
  fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
  fprintf(vtkFile,"DIMENSIONS %d %d %d\n", IMAX + 1, JMAX + 1, 1);
  fprintf(vtkFile,"POINTS %d DOUBLE\n", (IMAX + 1) * (JMAX + 1));
  // Write point data. Indexed so that i and j sit on the i - 0.5 face. 
  for(int i = 2; i <= IMAX + 2; i++)
    for(int j = 2; j <= JMAX + 2; j++)
    {
      double dX = (i - 2.0) * dDX;
      double dY = (j - 2.0) * dDY;
      fprintf(vtkFile, "%16.16e %1.16e %1.16e\n", dX, dY, 0.0);
    }
    // Cell data. indexed at current cell
  fprintf(vtkFile,"CELL_DATA %d\n", IMAX * JMAX);
  fprintf(vtkFile,"SCALARS temperature DOUBLE 1\n");
  fprintf(vtkFile,"LOOKUP_TABLE default\n");
  for(int i = 2; i <= IMAX + 1; i++)
    for(int j = 2; j <= JMAX + 1; j++)
    {
      fprintf(vtkFile, "%1.16e\n", a2dSoln[i][j]);
    }
  fclose(vtkFile);
}