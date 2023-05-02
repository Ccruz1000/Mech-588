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

  
  // Fixed temperature at y = 0 and y = 1
  // Top and bottom boundary
  for (i = 2; i <= IMAX + 1; i++) 
  {
    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    // 5th Order Interpolation
    a2dPrimal[i][0] = (1.0 / 12.0) * (-505.0 * a2dPrimal[i][2] + 335.0 * a2dPrimal[i][3] - 145.0 * a2dPrimal[i][4] + 27.0 * a2dPrimal[i][5]);
    a2dPrimal[i][1] = (1.0 / 12.0) * (-77.0 * a2dPrimal[i][2] + 43.0 * a2dPrimal[i][3] - 17.0 * a2dPrimal[i][4] + 3.0 * a2dPrimal[i][5]);
    a2dPrimal[i][JMAX + 2] = (1.0 / 12.0) * (-77.0 * a2dPrimal[i][JMAX + 1] + 43.0 * a2dPrimal[i][JMAX] - 17.0 * a2dPrimal[i][JMAX - 1] + 3.0 * a2dPrimal[i][JMAX - 2]);
    a2dPrimal[i][JMAX + 3] = (1.0 / 12.0) * (-505.0 * a2dPrimal[i][JMAX + 1] + 335.0 * a2dPrimal[i][JMAX] - 145.0 * a2dPrimal[i][JMAX - 1] + 27.0 * a2dPrimal[i][JMAX - 2]);



    // Fourth Order interpolation
    // a2dPrimal[i][0] = (1.0 / 3.0) * (-70.0 * a2dPrimal[i][2] + 32.0 * a2dPrimal[i][3] - 7.0 * a2dPrimal[i][4]);
    // a2dPrimal[i][1] = (1.0 / 3.0) * (-13.0 * a2dPrimal[i][2] + 5.0 * a2dPrimal[i][3] - a2dPrimal[i][4]);
    // a2dPrimal[i][JMAX + 2] = (1.0 / 3.0) * (-13.0 * a2dPrimal[i][JMAX + 1] + 5.0 * a2dPrimal[i][JMAX] - a2dPrimal[i][JMAX - 1]);
    // a2dPrimal[i][JMAX + 3] = (1.0 / 3.0) * (-70.0 * a2dPrimal[i][JMAX + 1] + 32.0 * a2dPrimal[i][JMAX] - 7.0 * a2dPrimal[i][JMAX - 1]);
  }
   
  // Fixed temperature at x = 0 and x = 1
  // Left and right boundary
  for (j = 2; j <= JMAX + 1; j++) 
  {

    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    
    // Fifth Order interpolation
    a2dPrimal[0][j] = (1.0 / 12.0) * (-505.0 * a2dPrimal[2][j] + 335.0 * a2dPrimal[3][j] - 145.0 * a2dPrimal[4][j] + 27.0 * a2dPrimal[5][j]);
    a2dPrimal[1][j] = (1.0 / 12.0) * (-77.0 * a2dPrimal[2][j] + 43.0 * a2dPrimal[3][j] - 17.0 * a2dPrimal[4][j] + 3.0 * a2dPrimal[5][j]);
    a2dPrimal[IMAX + 2][j] = (1.0 / 12.0) * (-77.0 * a2dPrimal[IMAX + 1][j] + 43.0 * a2dPrimal[IMAX][j] - 17.0 * a2dPrimal[IMAX - 1][j] + 3.0 * a2dPrimal[IMAX - 2][j]);
    a2dPrimal[IMAX + 3][j]  = (1.0 / 12.0) * (-505.0 * a2dPrimal[IMAX +  1][j] + 335.0 * a2dPrimal[IMAX][j] - 145.0 * a2dPrimal[IMAX - 1][j] + 27.0 * a2dPrimal[IMAX - 2][j]);
 
    // 4th Order Interpolation 
    // a2dPrimal[0][j] = (1.0 / 3.0) * (-70.0 * a2dPrimal[2][j] + 32.0 * a2dPrimal[3][j] - 7.0 * a2dPrimal[4][j]);
    // a2dPrimal[1][j] = (1.0 / 3.0) * (-13.0 * a2dPrimal[2][j] + 5.0 * a2dPrimal[3][j] - a2dPrimal[4][j]);
    // a2dPrimal[IMAX + 2][j] = (1.0 / 3.0) * (-13.0 * a2dPrimal[IMAX + 1][j] + 5.0 * a2dPrimal[IMAX][j] - a2dPrimal[IMAX - 1][j]);
    // a2dPrimal[IMAX + 3][j]  = (1.0 / 3.0) * (-70.0 * a2dPrimal[IMAX +  1][j] + 32.0 * a2dPrimal[IMAX][j] - 7.0 * a2dPrimal[IMAX - 1][j]);
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
void run_Poisson_4thOrder(void (*vApplyHomoBC_4thOrder)(double [IMAX+4][JMAX+4]), double (*dSource)(const double dX, const double dY), double (&a2dSoln)[IMAX + 4][JMAX + 4], double (&a2dSource)[IMAX+4][JMAX+4])
{
  int i, j;
  // Loop through all cells
  for (i = 0; i <= IMAX+3; i++) 
  {
    for (j = 0; j <= JMAX+3; j++)
    {
     // dX and dY represent the current X and Y coordinate, at the 
     // left and bottom cell face. 
     double dX = (i - 1.5) * dDX;
     double dY = (j - 1.5) * dDY; // We have a -1.5 here to account for the two ghost cells.
     a2dSoln[i][j] = 0; // Initialize solution to be 0
     /* Let's try a bit of integration here, four-point Gauss quadrature... */
     a2dSource[i][j] = 0.25 * (dSource(dX + dDX/sqrt(12), dY + dDY/sqrt(12)) +
                            dSource(dX + dDX/sqrt(12), dY - dDY/sqrt(12)) +
                            dSource(dX - dDX/sqrt(12), dY + dDY/sqrt(12)) +
                            dSource(dX - dDX/sqrt(12), dY - dDY/sqrt(12)));
    }
  }
  vApplyHomoBC_4thOrder(a2dSoln); // Apply boundary conditions
  
  int iIter = 0;
  double dMaxChange; // Used to find maximum possible change in temperature, before solution is 
  // considered steady state.

  // Here we will calculate the max change of all the cells. We will then
  // increment our iteration. If we reach the maximum number of iterations, or
  // if the solution reaches steady state, then we will break from the loop.
  do 
  {
    dMaxChange = dPointGSUpdate_4thOrder(a2dSoln, a2dSource, vApplyHomoBC_4thOrder, 1.7);
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
      double dExactCV = dExact((i - 1.5) * dDX, (j - 1.5) * dDY);
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

double dPointGSUpdate_4thOrder(double a2dSoln[IMAX+4][JMAX+4], double a2dSource[IMAX+4][JMAX+4], void (*vApplyHomoBC_4thOrder)(double [IMAX+4][JMAX+4]), const double dOmega)
{
  double dMaxChange = 0;
  int i, j;
  
  // Loop through interior cells
  for (i = 2; i <= IMAX + 1; i ++) 
  {
    for (j = 2; j <= JMAX + 1; j ++) 
    {
      // SECOND ORDER DISCRETIZATION
      // double dDeltaT = (dDY2 * (a2dSoln[i+1][j] + a2dSoln[i-1][j]) +
			   //               dDX2 * (a2dSoln[i][j+1] + a2dSoln[i][j-1]) -
		    //                dDX2*dDY2 * a2dSource[i][j]) / 
	     //                 (2 * (dDX2 + dDY2)) - a2dSoln[i][j];

      // BUG IS HERE, NEED TO DISCRETIZE SOURCE TERM ACCURATELY 
      double dDeltaT = (4.0 * dDX2 / 3.0) * (a2dSoln[i + 1][j] - 2.0 * a2dSoln[i][j] + a2dSoln[i - 1][j]) - (dDX2 / 12.0) * (a2dSoln[i + 2][j] - 2 * a2dSoln[i][j] - a2dSoln[i - 2][j]) + // X discretization
                       (4.0 * dDY2 / 3.0) * (a2dSoln[i][j + 1] - 2.0 * a2dSoln[i][j] + a2dSoln[i][j - 1]) - (dDY2 / 12.0) * (a2dSoln[i][j + 2] - 2 * a2dSoln[i][j] - a2dSoln[i][j - 2]) - // Y discretization
                       dDX2 * dDY2 * a2dSource[i][j] - // Source term
                       a2dSoln[i][j]; // Subtract solution to find deltaT

      // printf("Cell (%d, %d) - Source %f - DDX2 %f - DDY2 %f\n", i, j, a2dSource[i][j], dDX2, dDY2);
      // double dDeltaT = (dDX2 / 12.0) *  (-a2dSoln[i + 2][j] + 16.0 * a2dSoln[i + 1][j] - 30.0 * a2dSoln[i][j] + 16.0 * a2dSoln[i - 1][j] - a2dSoln[i - 2][j]) +
      //                  (dDY2 / 12.0) * (-a2dSoln[i][j + 2] + 16.0 * a2dSoln[i][j + 1] - 30.0 * a2dSoln[i][j] + 16.0 * a2dSoln[i][j - 1] - a2dSoln[i][j - 2]) -
      //                  dDX2 * dDY2 * a2dSource[i][j] - a2dSoln[i][j];

      a2dSoln[i][j] += dDeltaT * dOmega;
      // printf("%2d %2d %10f %10f\n", i, j, dDeltaT, dMaxChange);
      dMaxChange = MAX(dMaxChange, fabs(dDeltaT * dOmega));
    }
  }
  vApplyHomoBC_4thOrder(a2dSoln);
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

double calc_L2Norm_4thOrder(double a2dVec[IMAX + 4][JMAX + 4])
{
  double L2Norm = 0.0; // Store L2Norm of vector
  // Loop over interior points 
  for(int i = 2; i <= IMAX + 1; i++)
  {
    for(int j = 2; j <= JMAX + 1; j++)
    {
    L2Norm += a2dVec[i][j] * a2dVec[i][j];      
    }
  }
  return sqrt(L2Norm);
}