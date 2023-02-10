#ifndef TRI_THOMAS_H
#define TRI_THOMAS_H

/* Code for solving block-tridiagonal matrix problems, using the Thomas
   algorithm.  The only subroutine in here that you'll -need- to call is
   SolveThomas, although things like Add3x3 or AddVec might be useful,
   too. */

/* LHS array is sized as [*][3][3][3].  The last two indices identify
   the element within a block; the third index is the row and the fourth
   is the column of the Jacobian matrix.  The second index tells which
   block it is: 0 is below the main diagonal, 1 is on the main diagonal,
   2 is above the main diagonal.  The first index tells which block row
   you're looking at (the i or j index from the discretization). */

/* RHS array is [*][3].  The second index tells which element of the
   solution vector, and the first is the block row. */

#include <math.h>
#include <stdio.h>
#include "settings.h"

void SpewMatrix(double Source[3][3]);

void SpewVector(double Source[3]);

void CopyVec(const double Source[3],
             double Target[3]);

void Copy3x3(double Source[3][3],
             double Target[3][3]);

void Mult3x3(double A[3][3],
             double B[3][3],
             double C[3][3]);

void MultVec(double A[3][3],
             const double Vec[3],
             double Result[3]);

void Add3x3(double A[3][3],
            double B[3][3],
            const double Factor,
            double C[3][3]);

void AddVec(const double A[3],
            const double B[3],
            const double Factor,
            double C[3]);

void Invert3x3(double Block[3][3],
               double Inverse[3][3]);

void SolveBlockTri(double LHS[MAXSIZE][3][3][3],
                   double RHS[MAXSIZE][3],
                   int iNRows);

void InitLHS(double LHS[100][3][3][3], const int NRows);

#endif
