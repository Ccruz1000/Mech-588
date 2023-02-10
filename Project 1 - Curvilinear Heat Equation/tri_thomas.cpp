#include "tri_thomas.h"

void SpewMatrix(double Source[3][3])
{
  printf("%10.6f %10.6f %10.6f\n", Source[0][0], Source[1][0], Source[2][0]);
  printf("%10.6f %10.6f %10.6f\n", Source[0][1], Source[1][1], Source[2][1]);
  printf("%10.6f %10.6f %10.6f\n", Source[0][2], Source[1][2], Source[2][2]);
}

void SpewVector(double Source[3])
{
  printf("%10.6f %10.6f %10.6f\n", Source[0], Source[1], Source[2]);
}

void CopyVec(const double Source[3],
             double Target[3])
{
  Target[0] = Source[0];
  Target[1] = Source[1];
  Target[2] = Source[2];
}

void Copy3x3(double Source[3][3],
             double Target[3][3])
{
  Target[0][0] = Source[0][0];
  Target[0][1] = Source[0][1];
  Target[0][2] = Source[0][2];

  Target[1][0] = Source[1][0];
  Target[1][1] = Source[1][1];
  Target[1][2] = Source[1][2];

  Target[2][0] = Source[2][0];
  Target[2][1] = Source[2][1];
  Target[2][2] = Source[2][2];
}

void Mult3x3(double A[3][3],
             double B[3][3],
             double C[3][3])
{
  C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
  C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
  C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];

  C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
  C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
  C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];

  C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
  C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
  C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

void MultVec(double A[3][3],
             const double Vec[3],
             double Result[3])
{
  Result[0] = A[0][0] * Vec[0] + A[0][1] * Vec[1] + A[0][2] * Vec[2];
  Result[1] = A[1][0] * Vec[0] + A[1][1] * Vec[1] + A[1][2] * Vec[2];
  Result[2] = A[2][0] * Vec[0] + A[2][1] * Vec[1] + A[2][2] * Vec[2];
}

void Add3x3(double A[3][3],
            double B[3][3],
            const double Factor,
            double C[3][3])
{
  C[0][0] = A[0][0] + Factor * B[0][0];
  C[0][1] = A[0][1] + Factor * B[0][1];
  C[0][2] = A[0][2] + Factor * B[0][2];

  C[1][0] = A[1][0] + Factor * B[1][0];
  C[1][1] = A[1][1] + Factor * B[1][1];
  C[1][2] = A[1][2] + Factor * B[1][2];

  C[2][0] = A[2][0] + Factor * B[2][0];
  C[2][1] = A[2][1] + Factor * B[2][1];
  C[2][2] = A[2][2] + Factor * B[2][2];
}

void AddVec(const double A[3],
            const double B[3],
            const double Factor,
            double C[3])
{
  C[0] = A[0] + Factor * B[0];
  C[1] = A[1] + Factor * B[1];
  C[2] = A[2] + Factor * B[2];
}

void Invert3x3(double Block[3][3],
               double Inverse[3][3])
{
  double DetInv = 1. / (+Block[0][0] * Block[1][1] * Block[2][2] + Block[0][1] * Block[1][2] * Block[2][0] + Block[0][2] * Block[1][0] * Block[2][1] - Block[0][2] * Block[1][1] * Block[2][0] - Block[0][1] * Block[1][0] * Block[2][2] - Block[0][0] * Block[1][2] * Block[2][1]);

  /* Expand by minors to compute the inverse */
  Inverse[0][0] = +DetInv * (Block[1][1] * Block[2][2] -
                             Block[2][1] * Block[1][2]);
  Inverse[1][0] = -DetInv * (Block[1][0] * Block[2][2] -
                             Block[2][0] * Block[1][2]);
  Inverse[2][0] = +DetInv * (Block[1][0] * Block[2][1] -
                             Block[2][0] * Block[1][1]);
  Inverse[0][1] = -DetInv * (Block[0][1] * Block[2][2] -
                             Block[2][1] * Block[0][2]);
  Inverse[1][1] = +DetInv * (Block[0][0] * Block[2][2] -
                             Block[2][0] * Block[0][2]);
  Inverse[2][1] = -DetInv * (Block[0][0] * Block[2][1] -
                             Block[2][0] * Block[0][1]);
  Inverse[0][2] = +DetInv * (Block[0][1] * Block[1][2] -
                             Block[1][1] * Block[0][2]);
  Inverse[1][2] = -DetInv * (Block[0][0] * Block[1][2] -
                             Block[1][0] * Block[0][2]);
  Inverse[2][2] = +DetInv * (Block[0][0] * Block[1][1] -
                             Block[1][0] * Block[0][1]);
}

void SolveBlockTri(double LHS[MAXSIZE][3][3][3],
                   double RHS[MAXSIZE][3],
                   int iNRows)
{
  int j;
  double Inv[3][3];

  for (j = 0; j < iNRows - 1; j++)
  {
    /* Compute the inverse of the main block diagonal. */
    Invert3x3(LHS[j][1], Inv);
    /* Scale the right-most block diagonal by the inverse. */
    {
      double Temp[3][3];
      Mult3x3(Inv, LHS[j][2], Temp);
      Copy3x3(Temp, LHS[j][2]);
    }

    /* Scale the right-hand side by the inverse. */
    {
      double Temp[3];
      MultVec(Inv, RHS[j], Temp);
      CopyVec(Temp, RHS[j]);
    }

    /* Left-multiply the jth row by the sub-diagonal on the j+1st row
       and subtract from the j+1st row.  This involves the
       super-diagonal term and the RHS of the jth row. */
    {
      /* First the LHS manipulation */
#define A LHS[j + 1][0]
#define B LHS[j + 1][1]
#define C LHS[j][2]
      double Temp[3][3], Temp2[3][3];
      double TVec[3], TVec2[3];
      Mult3x3(A, C, Temp);
      Add3x3(B, Temp, -1., Temp2);
      Copy3x3(Temp2, B);

      /* Now the RHS manipulation */
      MultVec(A, RHS[j], TVec);
      AddVec(RHS[j + 1], TVec, -1., TVec2);
      CopyVec(TVec2, RHS[j + 1]);
#undef A
#undef B
#undef C
    }
  } /* Done with forward elimination loop */
  /* Compute the inverse of the last main block diagonal. */
  j = iNRows - 1;
  Invert3x3(LHS[j][1], Inv);

  /* Scale the right-hand side by the inverse. */
  {
    double Temp[3];
    MultVec(Inv, RHS[j], Temp);
    CopyVec(Temp, RHS[j]);
  }

  /* Now do the back-substitution. */
  for (j = iNRows - 2; j >= 0; j--)
  {
    /* Matrix-vector multiply and subtract. */
#define C LHS[j][2]
    RHS[j][0] -= (C[0][0] * RHS[j + 1][0] +
                  C[0][1] * RHS[j + 1][1] +
                  C[0][2] * RHS[j + 1][2]);
    RHS[j][1] -= (C[1][0] * RHS[j + 1][0] +
                  C[1][1] * RHS[j + 1][1] +
                  C[1][2] * RHS[j + 1][2]);
    RHS[j][2] -= (C[2][0] * RHS[j + 1][0] +
                  C[2][1] * RHS[j + 1][1] +
                  C[2][2] * RHS[j + 1][2]);
#undef C
  }
}

void InitLHS(double LHS[100][3][3][3], const int NRows)
{
  int i;
  for (i = 0; i < NRows; i++)
  {
    LHS[i][0][0][0] = 1. - 2.;
    LHS[i][0][0][1] = 2.;
    LHS[i][0][0][2] = 3.;
    LHS[i][0][1][0] = 4.;
    LHS[i][0][1][1] = 5. - 2.;
    LHS[i][0][1][2] = 6.;
    LHS[i][0][2][0] = 7.;
    LHS[i][0][2][1] = 8.;
    LHS[i][0][2][2] = 0. - 2.;

    LHS[i][1][0][0] = 1.;
    LHS[i][1][0][1] = 2.;
    LHS[i][1][0][2] = 3.;
    LHS[i][1][1][0] = 4.;
    LHS[i][1][1][1] = 5.;
    LHS[i][1][1][2] = 6.;
    LHS[i][1][2][0] = 7.;
    LHS[i][1][2][1] = 8.;
    LHS[i][1][2][2] = 0.;

    LHS[i][2][0][0] = 1. - 3.;
    LHS[i][2][0][1] = 2.;
    LHS[i][2][0][2] = 3.;
    LHS[i][2][1][0] = 4.;
    LHS[i][2][1][1] = 5. - 3.;
    LHS[i][2][1][2] = 6.;
    LHS[i][2][2][0] = 7.;
    LHS[i][2][2][1] = 8.;
    LHS[i][2][2][2] = 0. - 3.;
  }

  LHS[0][0][0][0] =
      LHS[0][0][0][1] =
          LHS[0][0][0][2] =
              LHS[0][0][1][0] =
                  LHS[0][0][1][1] =
                      LHS[0][0][1][2] =
                          LHS[0][0][2][0] =
                              LHS[0][0][2][1] =
                                  LHS[0][0][2][2] = 0.;

  LHS[NRows - 1][2][0][0] =
      LHS[NRows - 1][2][0][1] =
          LHS[NRows - 1][2][0][2] =
              LHS[NRows - 1][2][1][0] =
                  LHS[NRows - 1][2][1][1] =
                      LHS[NRows - 1][2][1][2] =
                          LHS[NRows - 1][2][2][0] =
                              LHS[NRows - 1][2][2][1] =
                                  LHS[NRows - 1][2][2][2] = 0.;
}
