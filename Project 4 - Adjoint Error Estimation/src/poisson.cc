#include "poisson.h"

/**************************************************************************************
// Adjoint Methods Programming Assignment Poisson Solver File 
// Code for Poisson solver provided by Dr. Carl Olivier-Gooch, modified to work with 
// adjoint methods for error estimation.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Define exact solution <- NOT FOR THE RIGHT PROBLEM
double dExact(const double dX, const double dY)
{
  return (M_PI/8)*(sin(M_PI*dX) * sin(M_PI*dY)); // Would need to calculate this
}

// Define initial condition <- STEADY STATE SO IRRELEVANT
double dTInit(const int i, const int j)
{
  double dX = (i-0.5)*dDX;
  double dY = (j-0.5)*dDY;
  return 0;
}

// Apply homogenous boundary conditions
// 2 extra points for ghost cells at boundary
void vApplyHomoBC(double a2dPrimal[IMAX+2][JMAX+2])
{
  int i, j;
  double dX;
  double dY = (j-0.5)*dDY;
  
  // Fixed temperature at y = 0 and y = 1
  // Top and bottom boundary
  for (i = 1; i <= IMAX; i++) 
  {
    dX = (i-0.5)*dDX;
    // a2dPrimal[i][0] = 2 * dExact(dX, 0) - a2dPrimal[i][1]; // Bottom boundary
    // a2dPrimal[i][JMAX+1] = 2 * dExact(dX, 1) - a2dPrimal[i][JMAX]; // Top boundary
    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    a2dPrimal[i][0] = -a2dPrimal[i][1];
    a2dPrimal[i][JMAX + 1] = -a2dPrimal[i][JMAX];
  }
   
  // Fixed temperature at x = 0 and x = 1
  // Left and right boundary
  for (j = 1; j <= JMAX; j++) 
  {
    dY = (j-0.5)*dDY;
    // a2dPrimal[0][j] = 2 * dExact(0, dY) - a2dPrimal[1][j]; // Left boundary
    // a2dPrimal[IMAX+1][j] = 2 * dExact(1, dY) - a2dPrimal[IMAX][j]; // Right boundary
    // Below we will set the boundary conditions based on problem 1. Here we have 
    // all boundary conditions set to be 0, as they are homogenous on all sides
    a2dPrimal[0][j] = -a2dPrimal[1][j];
    a2dPrimal[IMAX + 1][j] = -a2dPrimal[IMAX][j];
  }
}

/* Point Gauss-Seidel update with SOR. <- SOR = successive over-relaxation  */
// Used to solve system of linear equations. SOR just gives faster convergence
double dPointGSUpdate(double a2dSoln[IMAX+2][JMAX+2], double a2dSource[IMAX+2][JMAX+2], 
                      void (*vApplyBC)(double [IMAX+2][JMAX+2]), const double dOmega)
{
  double dMaxChange = 0;
  int i, j;
  
  //for (j = JMAX; j >= 1; j --) {
  //  for (i = IMAX; i >= 1; i --) {
  for (i = 1; i <= IMAX; i ++) 
  {
    for (j = 1; j <= JMAX; j ++) 
    {
      double dDeltaT = (dDY2 * (a2dSoln[i+1][j] + a2dSoln[i-1][j]) +
			                 dDX2 * (a2dSoln[i][j+1] + a2dSoln[i][j-1]) -
		                   dDX2 * dDY2 * a2dSource[i][j]) / 
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

// source term for primal function
double Primal_Source(const double dX, const double dY)
{
  static const double dVal = -M_PI*M_PI*M_PI / 4.0; // Primal
  return  dVal * sin(M_PI*dX) * sin(M_PI*dY);
  // return 0.0;
}

double Dual_Source(const double dX, const double dY)
{
  static const double dVal = (M_PI*M_PI*M_PI*M_PI*M_PI / 2.0); // Dual
  return dVal * dX * (1 - dX) * dY * (1 - dY);
}

// Used to run the poisson solver given the different test cases
void run_Poisson(void (*vApplyBC)(double [IMAX+2][JMAX+2]), double (*dSource)(const double dX, const double dY), double (&a2dSoln)[IMAX + 2][JMAX + 2], double (&a2dSource)[IMAX+2][JMAX+2])
{
  // double a2dSoln[IMAX+2][JMAX+2];
  // double a2dSource[IMAX+2][JMAX+2];
  int i, j;
  // Loop through all cells
  for (i = 0; i <= IMAX+1; i++) 
  {
    for (j = 0; j <= JMAX+1; j++)
    {
     // Indexed so that i and j represent the current face i+1/2, j+1/2. i+1/2 - 0.5 would give the 
     // position at the cell center.

     double dX = (i-0.5) * dDX;
     double dY = (j-0.5) * dDY;
     a2dSoln[i][j] = 0.0; // Initialize solution to be 0
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
  //        0.25*(a2dSoln[IMAX/2][JMAX/2] + a2dSoln[IMAX/2][JMAX/2+1] +
  //        a2dSoln[IMAX/2+1][JMAX/2] + a2dSoln[IMAX/2+1][JMAX/2+1])); 
}

// This function will calculate the gradient of the solution at the cell center
double calc_grad(double u_1, double u_2, double dDx)
{
  // u_1 will be the solution at the +1 increment, while u_2 will be the 
  // solution at the -1 increment. For example, if we are calculating the gradient
  // with respect to x, u_1 will be the cell averaged solution for the cell i+1, and 
  // u_2 will be the cell averaged solution for the cell i-1. 
  return (u_1 - u_2) / (2 * dDx);
}

// This function will calculate the L2Norm of a 2D varray
double calc_L2Norm(double a2dVec[IMAX + 2][JMAX + 2])
{
  double L2Norm = 0.0; // Store L2Norm of vector
  // Loop over interior points 
  for(int i = 1; i <= IMAX; i++)
  {
    for(int j = 1; j <= JMAX; j++)
    {
    L2Norm += a2dVec[i][j] * a2dVec[i][j];      
    }
  }
  return sqrt(L2Norm);
}

// This function will calculate the functional, given a solution term and a source term
// In our case, the functional is J=(u, g)=(v, f)
void calc_functional(double a2dSoln[IMAX + 2][JMAX + 2], double a2dSource[IMAX + 2][JMAX + 2], double (&a2dFunctional))
{
  // We have already computed the Gauss quadrature integration for our solution. As such, we do not
  // need to compute this again.

  // Start by looping over the cells internal cells
  for(int i = 1; i <= IMAX; i++)
  {
    for(int j = 1; j <= JMAX; j++)
    {
      // To start we will compute the gradient of the solution within the cell
      double grad_x = calc_grad(a2dSoln[i + 1][j], a2dSoln[i - 1][j], dDX);
      double grad_y = calc_grad(a2dSoln[i][j + 1], a2dSoln[i][j - 1], dDY);

      /*
        Now we can compute the temperature at the Gauss quadrature points.
        ----------------
        |              |
        |   1     2    |
        |      c       |
        |   3     4    |
        |              |
        ----------------
        Above, we can see our Gauss-quadrature points. We can call the centroid XC and YC. Point 1
        would be at (XC-dx/sqrt(12), YC+dx/sqrt(12)), point 2 would be at (XC+dx/sqrt(12), YC+dx/sqrt(12))
        point 3 would be at (XC-dx/sqrt(12), YC-dx/sqrt(12)) and point 4 would be at
        (XC+dx/sqrt(12), YC-dx/sqrt(12)). Since it is all relative, we can just use dx/sqrt(12), 
        -dx/sqrt(12), dy/sqrt(12), -dy/sqrt(12).
      */
      double GCellSoln1 = a2dSoln[i][j] - grad_x * (dDX/sqrt(12)) + grad_y * (dDY/sqrt(12));
      double GCellSoln2 = a2dSoln[i][j] + grad_x * (dDX/sqrt(12)) + grad_y * (dDY/sqrt(12));
      double GCellSoln3 = a2dSoln[i][j] - grad_x * (dDX/sqrt(12)) - grad_y * (dDY/sqrt(12));
      double GCellSoln4 = a2dSoln[i][j] - grad_x * (dDX/sqrt(12)) - grad_y * (dDY/sqrt(12));
      double GCellSolnTot = 0.25 * (GCellSoln1 + GCellSoln2 + GCellSoln3 + GCellSoln4);

      // Now we can compute the functional at the cell by multiplying the 
      // gauss integral of the solution in the cell with the source term in the cell
      // Remember we have to multiply by the cell area as well. 
      a2dFunctional += a2dSource[i][j] * GCellSolnTot * dDX * dDY;
    }
  }
}

void storeVTKSolution(double a2dSoln[IMAX + 2][JMAX + 2], std::string fileName)
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
  for(int i = 1; i <= IMAX + 1; i++)
    for(int j = 1; j <= JMAX + 1; j++)
    {
      double dX = (i - 1.0) * dDX;
      double dY = (j - 1.0) * dDY;
      fprintf(vtkFile, "%16.16e %1.16e %1.16e\n", dX, dY, 0.0);
    }
    // Cell data. indexed at current cell
  fprintf(vtkFile,"CELL_DATA %d\n", IMAX * JMAX);
  fprintf(vtkFile,"SCALARS temperature DOUBLE 1\n");
  fprintf(vtkFile,"LOOKUP_TABLE default\n");
  for(int i = 1; i <= IMAX; i++)
    for(int j = 1; j <= JMAX; j++)
    {
      fprintf(vtkFile, "%1.16e\n", a2dSoln[i][j]);
    }
  fclose(vtkFile);
}