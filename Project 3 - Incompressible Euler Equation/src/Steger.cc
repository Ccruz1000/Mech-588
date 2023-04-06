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
	std::array<std::vector<double>, 3> Fplus = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)}; // Store positive flux vector
	std::array<std::vector<double>, 3> Fminus = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)}; // Store negative flux vector
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
	calc_bound(Fplus, Prim, nx);
	calc_bound(Fminus, Prim, nx);

	// Matrices to calculate Roes flux
	// double lambda_plus[3][3]; // Store positive eigen value matrix
	// double lamda_minus[3][3]; // Store negative eigenvalue matrix
	// double A_plus[3][3] = {0}; // Store Jacobian matrix based on positive eigenvalues
	double A_minus[3][3] = {0}; // Store Jacobian matrix based on negative eigenvalues

	// Calculate flux at each cell
	for(int i = 1; i < nx - 1; i++)
	{
		// Initialize matrices
		double du_dv[3][3], XR[3][3], lambda_plus[3][3], lambda_minus[3][3], XL[3][3], dv_du[3][3];

		assemble_matrices(Prim[0][i], Prim[1][i], c[i], gamma, du_dv, XR, XL, dv_du);

		//  // Values to store the 3 individual eigenvalues

		// lam1 = (0.5 * (Prim[1][i] + Prim[1][i + 1]));
		// lam2 = (0.5 * (Prim[1][i] + c[i] + Prim[1][i + 1] + c[i + 1]));
		// lam3 = (0.5 * (Prim[1][i] - c[i] + Prim[1][i + 1] - c[i + 1]));

		// // Populate eigenvector matrices
		// double lambda_plus[3][3] = {{0.5 * (lam1 + std::sqrt(lam1 * lam1 + epsilon)), 0, 0},
		// 			               {0, 0.5 * (lam2 + std::sqrt(lam2 * lam2 + epsilon)), 0},
		// 			               {0, 0, 0.5 * (lam3 + std::sqrt(lam3 * lam3 + epsilon))}};

		// double lambda_minus[3][3] = {{0.5 * (lam1 - std::sqrt(lam1 * lam1 + epsilon)), 0, 0},
		// 			               {0, 0.5 * (lam2 - std::sqrt(lam2 * lam2 + epsilon)), 0},
		// 			               {0, 0, 0.5 * (lam3 - std::sqrt(lam3 * lam3 + epsilon))}};


		// // Now that all the matrices are defined, we can calculate Aplus and Aminus
		// calc_Jacobian(du_dv, XR, lambda_plus, XL, dv_du, A_plus);
		// calc_Jacobian(du_dv, XR, lambda_minus, XL, dv_du, A_minus);

		// // Now, we know that the flux is equivalent to (1/2)(Fi+Fi+1) - 1/2A~(Ui+1-Ui)
		// // Store values in a double array, to be transferred to the vector later
		// double U_cell[3] = {U[0][i + 1] - U[0][i], U[1][i + 1] - U[1][i], U[2][i + 1] - U[2][i]};
		// double Roe[3] = {0};
		// matrix_vector(A_tild, U_dif, Roe);

		// F[0][i] = 0.5 * (Fnorm[0][i] + Fnorm[0][i + 1]) - 0.5 * Roe[0];
		// F[1][i] = 0.5 * (Fnorm[1][i] + Fnorm[1][i + 1]) - 0.5 * Roe[1];
		// F[2][i] = 0.5 * (Fnorm[2][i] + Fnorm[2][i + 1]) - 0.5 * Roe[2];
	}
}