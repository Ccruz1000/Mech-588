// Included header
#include "Steger.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Steger Warming Flux Vector Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// The following will compute the flux using the Steger-Warming scheme 
void Steger_warming(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx)
{	
	// Initialize all required variables, vectors, and arrays
	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	std::vector<double> c = std::vector<double>(nx); // Store speed of sound
	std::array<std::vector<double>, 3> F_plus = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::array<std::vector<double>, 3> F_minus = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	double lam1, lam2, lam3; // Store eigenvalues
	
	// calculate primitive variables from conserved variables
	calc_Prim(U, Prim, gamma, nx);

	// Calculate speed of sound
	for(int i = 0; i < nx; i++)
	{
		c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
	}

	// apply boundary conditions to the positive, negative, and final flux vectors
	calc_bound(F, Prim, nx);
	calc_bound(F_plus, Prim, nx);
	calc_bound(F_minus, Prim, nx);
 

	// Calculate flux at each cell
	for(int i = 1; i < nx - 1; i++)
	{
		// We will start by calculating F+
		double A_plus[3][3] = {0}; // Store Jacobian 
		double lambda_plus[3][3]; // Store matrix of eigenvalues
		double lam_plus1, lam_plus2, lam_plus3; // Store positive eigenvalues

		// Store solution vector at current cell in double array to allow for calculation of F
		double U_cell[3] = {U[0][i], U[1][i], U[2][i]};

		// Initialize and assemble matrices
		double du_dv[3][3], XR[3][3], lambda_plus[3][3], lambda_minus[3][3], XL[3][3], dv_du[3][3];
		assemble_matrices(Prim[0][i], Prim[1][i], c[i], gamma, du_dv, XR, XL, dv_du);

		// Calculate Eiegenvalues
		lam1 = (Prim[1][i] + Prim[1][i +1]) / 2;
		lam2 = (Prim[1][i] + c[i] + Prim[1][i +1] + c[i]) / 2;
		lam2 = (Prim[1][i] - c[i] + Prim[1][i +1] - c[i]) / 2;

		// Assemble Lambda+
		lam_plus1 = 0.5 * (lam1 + std::sqrt(lam1 * lam1 + epsilon));
		lam_plus2 = 0.5 * (lam2 + std::sqrt(lam2 * lam2 + epsilon));
		lam_plus2 = 0.5 * (lam3 + std::sqrt(lam3 * lam3 + epsilon));
		assemble_lambda(lam_plus1, lam_plus2, lam_plus3, lambda_plus);

		// Compute A+
		calc_Jacobian(du_dv, XR, lambda_plus, XL, dv_du, A_plus);

		// Compute F+
		double TempF_plus[3] = {0};
		matrix_vector(A_plus, U_cell, TempF_plus);
		F_plus[0][i] = TempF_plus[0];
		F_plus[1][i] = TempF_plus[1];
		F_plus[2][i] = TempF_plus[2];

		// Then we will calculate F-
		double A_minus[3][3] = {0}; // Store Jacobian 
		double lambda_minus[3][3]; // Store matrix of eigenvalues
		double lam_minus1, lam_minus2, lam_minus3; // Store positive eigenvalues

		// Assemble Lambda-
		lam_minus1 = 0.5 * (lam1 - std::sqrt(lam1 * lam1 + epsilon));
		lam_minus2 = 0.5 * (lam2 - std::sqrt(lam2 * lam2 + epsilon));
		lam_minus2 = 0.5 * (lam3 - std::sqrt(lam3 * lam3 + epsilon));
		assemble_lambda(lam_minus1, lam_minus2, lam_minus3, lambda_minus);

		// Compute A-
		calc_Jacobian(du_dv, XR, lambda_minus, XL, dv_du, A_minus);

		// Compute F+
		double TempF_minus[3] = {0};
		matrix_vector(A_minus, U_cell, TempF_minus);
		F_minus[0][i] = TempF_minus[0];
		F_minus[1][i] = TempF_minus[1];
		F_minus[2][i] = TempF_minus[2];
	}

	// Populate flux vector
	for(int i = 1; i < nx - 1; i++)
		for(int j = 0; j < 3; j++)
		{
			F[j][i] = F_plus[j][i] + F_minus[j][i];
		}

	// apply boundary conditions to the positive, negative, and final flux vectors
	calc_bound(F, Prim, nx);
	calc_bound(F_plus, Prim, nx);
	calc_bound(F_minus, Prim, nx);	
}
