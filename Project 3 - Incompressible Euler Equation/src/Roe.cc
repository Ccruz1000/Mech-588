// Included libraries found in header file
#include "Roe.h"


/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Goated video - https://www.youtube.com/watch?v=F_PsFHvt8IU&ab_channel=P.S.Volpiani

// This function will calculate the Roe averaged values
void calc_Roe_avg(const std::array<std::vector<double>, 4> Prim, double &rho_tild, double &u_tild, double &c_tild, double gamma, int nx)
{
	double h_left, h_right, h_tild, P_tild; // Store enthalpy and pressure 

	// Loop through cells to calculate Rho averaged values
	for(int i = 1; i < nx; i++)
	{
		rho_tild = std::sqrt(Prim[0][i - 1] * Prim[0][i]); // Populate Roe averaged density 
		u_tild = (Prim[0][i] * Prim[1][i] + rho_tild * Prim[1][i - 1]) / (Prim[0][i] * rho_tild); // Populate Roe averaged velocity

		// Calculate enthalpy
		h_left = (gamma / (gamma - 1)) * (Prim[2][i - 1] / Prim[0][i - 1]);
		h_right = (gamma / (gamma - 1)) * (Prim[2][i] / Prim[0][i]);
		h_tild = (Prim[0][i] * h_right + rho_tild * h_left) / (Prim[0][i] + rho_tild); 

		//Calculate Roe averaged pressure from Roe avereaged enthalpy
		P_tild = rho_tild * h_tild * (gamma - 1) / gamma;

		// Use Roe averaged pressure to find Roe averaged speed of sound
		c_tild = std::sqrt(gamma * P_tild / rho_tild);
	}
}


// The following will compute the flux using the Roe scheme
void Roe_flux(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx)
{
	// Initialize all required variables, vectors and arrays
	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	double lam1, lam2, lam3; // Store eigenvalues
	double abslam1, abslam2, abslam3; // Store absolute value of eigenvalues
	double rho_tild; // Store Roe averaged density
	double u_tild; // Store Roe averaged velocity 
	double c_tild; // Store Roe averaged speed of sound

	// Calculate primitive variables from conserved variables
	calc_Prim(U, Prim, gamma, nx);

	// Calculate Roe Averaged values
	calc_Roe_avg(Prim, rho_tild, u_tild, c_tild, gamma, nx);

	// Loop through cells
	for(int i = 1; i < nx; i++)
	{
		// Intialize arrays to calculate flux 
		double F_roe[3] = {0}; // Store Roe's flux
		double F_normR[3] = {0}; // Store analytical Euler flux based on values in cell to right of face
		double F_normL[3] = {0}; // Store analytical Euler flux based on values in cell to left of face

		// Initialize matrices for Roe flux calculation
		double du_dv[3][3], XR[3][3], lambda_abs[3][3], XL[3][3], dv_du[3][3];
		assemble_matrices(rho_tild, u_tild, c_tild, gamma, du_dv, XR, XL, dv_du);
			
		// Calculate eigenvalues
		lam1 = 0.5 * (u_tild);
		lam2 = 0.5 * (u_tild + c_tild);
		lam3 = 0.5 * (u_tild - c_tild);
		
		// Assemble lambda
		abslam1 = std::sqrt(lam1 * lam1 + epsilon);
		abslam2 = std::sqrt(lam2 * lam2 + epsilon);
		abslam3 = std::sqrt(lam3 * lam3 + epsilon);
		assemble_lambda(abslam1, abslam2, abslam3, lambda_abs);

		// Compute A~
		double Atild[3][3] = {0};
		calc_Jacobian(du_dv, XR, lambda_abs, XL, dv_du, Atild);

		// Compute Euler and Roe fluxes
		F_normR[0] = Prim[0][i] * Prim[1][i]; // rho*u
		F_normR[1] = Prim[0][i] * Prim[1][i] * Prim[1][i] + Prim[2][i]; // rho*u^2 + P
		F_normR[2] = Prim[1][i] * (Prim[2][i] + Prim[3][i]); // u*(P+Et)

		F_normL[0] = Prim[0][i - 1] * Prim[1][i - 1]; // rho*u
		F_normL[1] = Prim[0][i - 1] * Prim[1][i - 1] * Prim[1][i - 1] + Prim[2][i - 1]; // rho*u^2 + P
		F_normL[2] = Prim[1][i - 1] * (Prim[2][i - 1] + Prim[3][i - 1]); // u*(P+Et)

		double tempF_Roe[3] = {0}; // Temporary vector to store Roe flux
		double U_cell[3]; // Temporary solution vector to store solution at cell

		// Populate U_cell
		for(int j = 0; j < 3; j++)
		{
			U_cell[j] = U[j][i] - U[j][i - 1];
		}

		// Compute Roe flux
		matrix_vector(Atild, U_cell, F_roe);

		// Calculate flux at face from Euler and Roe flux
		for(int j = 0; j < 3; j++)
		{
			F[j][i] = 0.5 * (F_normL[j] + F_normR[j]) - 0.5 * F_roe[j]; 
		}

		calc_bound(F, Prim, nx);
	}

	// // Apply boundary conditions to Euler and Roe flux
	// calc_bound(F_norm, Prim, nx);
	// calc_bound(F_roe, Prim, nx);

	// // Calculate flux at interior points
	// for(int i = 1; i < nx - 1; i++)
	// 	for(int j = 0; j < 3; j++)
	// 	{
	// 		F[j][i] = 0.5 * (F_norm[j][i + 1] + F_norm[j][i]) - 0.5 * (F_roe[j][i + 1] - F_roe[j][i]);
	// 	}

	// // Apply boundary conditions
	// calc_bound(F, Prim, nx);
}
