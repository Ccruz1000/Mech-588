// Included header
#include "Steger.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Steger Warming Flux Vector Splitting File Header
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// This function will take the density, velocity, and speed of sound from the appropriate 
// side of the cell. It will then calculate the Steger flux values F+ and F-
void calc_Steger_flux(double rho, double u, double c, double gamma, double lam[3], double U[3], double (&F)[3])
{
	// Initialize and assemble matrices
	double du_dv[3][3], XR[3][3], lambda_plus[3][3], lambda_minus[3][3], XL[3][3], dv_du[3][3];
	assemble_matrices(rho, u, c, gamma, du_dv, XR, XL, dv_du);

	// Initialize and assemble lambda
	double lambda[3][3];
	assemble_lambda(lam[0], lam[1], lam[2], lambda);

	// Initialize and compute A
	double A[3][3] = {0};
	calc_Jacobian(du_dv, XR, lambda, XL, dv_du, A);

	// Compute flux
	matrix_vector(A, U, F);
}

// The following will compute the flux using the Steger-Warming scheme 
void Steger_warming(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx)
{	
	// Initialize all required variables, vectors, and arrays
	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	std::vector<double> c = std::vector<double>(nx); // Store speed of sound
	double lam1, lam2, lam3; // Store eigenvalues
	
	// calculate primitive variables from conserved variables
	calc_Prim(U, Prim, gamma, nx);

	// Calculate speed of sound
	for(int i = 0; i < nx; i++)
	{
		c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
	} 

	// Loop through internal cell faces
	for(int i = 1; i < nx; i++)
	{
		// We can start by assigning values based on the cell to the left, and the cell to the right
		double rho_left, rho_right;
		rho_left = Prim[0][i - 1]; 
		rho_right = Prim[0][i];

		double u_left, u_right;
		u_left = Prim[1][i - 1];
		u_right = Prim[1][i];

		double c_left, c_right;
		c_left = c[i - 1];
		c_right = c[i];

		// Now we can calculate the eigenvalues using the averaged data at the control volume boundary
		lam1 = (u_left + u_right) / 2.0;
		lam2 = (u_left + c_left + u_right + c_right) / 2.0;
		lam3 = (u_left - c_left + u_right - c_right) / 2.0;

		// We will start by calculating F+, using the solution values in the left cell
		// Store solution vector at current cell in double array based on left cell to allow for calculation of F+
		double U_plus[3] = {U[0][i - 1], U[1][i - 1], U[2][i - 1]};

		// Compute positive eigenvalues
		double lam_plus[3] = {0.5 * (lam1 + std::sqrt(lam1 * lam1 + epsilon)), 0.5 * (lam2 + std::sqrt(lam2 * lam2 + epsilon)), 0.5 * (lam3 + std::sqrt(lam3 * lam3 + epsilon))};
		
		// Initialize F+ using values from left cell
		double F_plus[3] = {0};

		// Calculate F+
		calc_Steger_flux(rho_left, u_left, c_left, gamma, lam_plus, U_plus, F_plus);

		// Now we can repeat the above process for F-
		// Store solution vector at current cell in double array based on right cell to allow for calculation of F-
		double U_minus[3] = {U[0][i], U[1][i], U[2][i]};

		// Compute negative eigenvalues
		double lam_minus[3] = {0.5 * (lam1 - std::sqrt(lam1 * lam1 + epsilon)), 0.5 * (lam2 - std::sqrt(lam2 * lam2 + epsilon)), 0.5 * (lam3 - std::sqrt(lam3 * lam3 + epsilon))};

		// Initialize F-
		double F_minus[3] = {0};

		// Compute F- using values from right cell
		calc_Steger_flux(rho_right, u_right, c_right, gamma, lam_minus, U_minus, F_minus);

		// Assign values to flux vector
		for(int j = 0; j < 3; j++)
		{
			F[j][i] = F_plus[j] + F_minus[j];
		}
	}

	// Apply boundary conditions to cell faces missed by above loop 
	calc_bound(F, Prim, nx);
}
