// Included libraries found in header file
#include "Roe.h"


/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Goated video - https://www.youtube.com/watch?v=F_PsFHvt8IU&ab_channel=P.S.Volpiani

// The following will compute the flux using the Roe scheme
void Roe_flux(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, int nx)
{

	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	std::vector<double> c = std::vector<double>(nx); // Store speed of sound
	// Compute the primitive variables from the conserved variables
	for(int i = 1; i < nx - 1; i++)
	{
		Prim[0][i] = U[0][i]; // Density
		Prim[1][i] = U[1][i] / U[0][i]; // Velocity
		Prim[3][i] = U[2][i]; // Energy
		Prim[2][i] = 0.5 * (1 - gamma) * U[1][i] * U[1][i] / U[0][i] + (gamma - 1) * U[2][i];  // Pressure (P = 0.5 * (1 - gamma) * U2^2/U1 + (gamma - 1) * U3
		// c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
		// calculate enthalpy
		
		if(Prim[0][i] < 0)
		{
			printf("Warning! Negative density found.\n");
		}

		if(Prim[0][i] == 0)
		{
			printf("Warning! Zero density found\n");
		}

		if(Prim[2][i] < 0)
		{
			printf("Warning! Negative pressure found.\n");
		}
	}
	// Store normal flux
	std::array<std::vector<double>, 3> Fnorm = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)}; 

	for(int i = 1; i < nx - 1; i++)
	{
		Fnorm[0][i] = Prim[0][i] * Prim[1][i]; // rho*u
		Fnorm[1][i] = Prim[0][i] * Prim[1][i] * Prim[1][i] + Prim[2][i]; // rho*u^2 + P
		Fnorm[2][i] = Prim[1][i] * (Prim[2][i] + Prim[3][i]); // u*(P+Et)
	}

	// Boundary conditions
	Fnorm[0][0] = 0; // UPDATE THIS TO BE PROPER VALUE
	Fnorm[1][0] = Prim[2][1]; // UPDATE THIS TO PROPER VALUE
	Fnorm[2][0] = 0; // UPDATE THIS TO PROPER VALUE

	Fnorm[0][nx - 1] = 0; // UPDATE THIS TO BE PROPER VALUE
	Fnorm[1][nx - 1] = Prim[2][nx - 2]; // UPDATE THIS TO PROPER VALUE
	Fnorm[2][nx - 1] = 0; // UPDATE THIS TO PROPER VALUE

	// Store Roe averaged values
	std::array<std::vector<double>, 3> U_tild = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::vector<double> u_tild = std::vector<double>(nx);
	double rho_tild, h_tild, h_left, h_right, P_tild;
	for(int i = 1; i < nx - 1; i++)
	{	
		// Calculate rho averaged density and velocity
		rho_tild = std::sqrt(Prim[0][i + 1] * Prim[0][i]); 
		u_tild[i] = (Prim[0][i + 1] * Prim[1][i + 1] + rho_tild * Prim[1][i]) / (Prim[0][i + 1] + rho_tild);

		// Calculate enthalpy
		h_left = (gamma / (gamma - 1)) * (Prim[2][i] / Prim[0][i]);
		h_right = (gamma / (gamma - 1)) * (Prim[2][i + 1] / Prim[0][i + 1]);
		h_tild = (Prim[0][i + 1] * h_right + rho_tild * h_left) / (Prim[0][i + 1] + rho_tild);

		// Use Roe averaged enthalpy to find the Roe averaged pressure
		P_tild = rho_tild * h_tild * (gamma - 1) / gamma;

		// Use Roe averaged pressure to find Roe averaged speed of sound
		c[i] = std::sqrt(gamma * P_tild / rho_tild);
		
		// Populate Roe averaged solution vector
		U_tild[0][i] = rho_tild;
		U_tild[1][i] = rho_tild * u_tild[i];
		U_tild[2][i] = h_tild;
	} 

	// Matrices to calculate Roes flux
	double lam1, lam2, lam3;
	double A_tild[3][3] = {0};

	// Boundary condition
	F[0][0] = 0; 
	F[1][0] = Prim[2][1]; 
	F[2][0] = 0; 

	F[0][nx - 1] = 0; 
	F[1][nx - 1] = Prim[2][nx - 2]; 
	F[2][nx - 1] = 0;

	// Calculate flux at each cell
	for(int i = 1; i < nx - 1; i++)
	{
		// Assemble matrices to perform calculation of Jacobian matrix
		double XR[3][3] = { {1.0, U_tild[0][i] / c[i], -U_tild[0][i] / c[i]},
					        {0, 1.0, 1.0},
					        {0, U_tild[0][i] * c[i], -U_tild[0][i] * c[i]} };
		/*
				|u     rho/c       -rho/c    |
		XR =    |0     1           1         |
				|0     rho*c       -rho*c    |
		*/

		double XL[3][3] = { {1, 0, -1 / (c[i] * c[i])},
					        {0, 1.0 / 2.0, 1.0 / (2.0 * U_tild[0][i] * c[i])},
					        {0, 1.0 / 2.0, -1.0 / (2.0 * U_tild[0][i] * c[i])} };
		        
		/*
				|1     0             -1/c^2      |
		XL =    |0     1/2           1/(2(rho*c) |
				|0     1/2           -1/(2(rho*c)|
		*/ 

		double du_dv[3][3] = {{1, 0, 0},
		                      {u_tild[i], U_tild[0][i], 0}, 
		                      {0.5 * u_tild[i] * u_tild[i], U_tild[0][i] * u_tild[i], 1.0 / (gamma - 1)}};
		/*
				|1     0     0              |
		du_dv = |u     rho   0              |
				|u^2/2 rho*u 1 / (gamma - 1)|
		*/

		double dv_du[3][3] = {{1, 0, 0},
					          {-u_tild[i] / U_tild[0][i], 1 / U_tild[0][i], 0},
					          {(gamma - 1) * 0.5 * u_tild[i] * u_tild[i], -(gamma - 1) * u_tild[i], gamma - 1}};

		/*
				|1                0              0        |
		dv_du = |-u/rho           1/rho          0        |
				|(gamma-1)u^2/2   -(gamma-1)u    gamma - 1|
		*/

		 // Values to store the 3 individual eigenvalues

		lam1 = (0.5 * (u_tild[i] + u_tild[i + 1]));
		lam2 = (0.5 * (u_tild[i] + c[i] + u_tild[i + 1] + c[i + 1]));
		lam3 = (0.5 * (u_tild[i] - c[i] + u_tild[i + 1] - c[i + 1]));

		double lambda_abs[3][3] = {{std::sqrt(lam1 * lam1 + epsilon), 0, 0},
					               {0, std::sqrt(epsilon + lam2 * lam2), 0},
					               {0, 0, std::sqrt(epsilon + lam3 * lam3)}};

		/*
				    |(((ui+ui+1)/2)^2+epsilon)^(1/2) 0                                             0                                            |
		lambda =    |0                               ((((ui+c)+((ui+1)+(ci+1))/2)^2+epsilon)^(1/2) 0                                            |
				    |0                               0                                             ((((ui+c)+((ui+1)+(ci+1))/2)^2+epsilon)^(1/2)|
		*/ 

		// Now that all the matrices are defined, we can calculate A~
		// A~ = (du_dv) * XR * lambda_abs * XL * (dv_du)
		// Intermediate matrices to store values
		double inter1[3][3] = {0};
		double inter2[3][3] = {0};
		// A = du_dv*XR
		matrix_mult3x3(du_dv, XR, inter1);
		// B = A * lamda_abs
		matrix_mult3x3(inter1, lambda_abs, inter2);
		// Zero inter1, and then C = B * XL
		for(int k = 0; k < 3; k++)
			for(int j = 0; j < 3; j++)
			{
				inter1[k][j] = 0; 
			}
		matrix_mult3x3(inter2, XL, inter1);
		// Finally, A~ = C * dv_du
		matrix_mult3x3(inter1, dv_du, A_tild);

		// Now, we know that the flux is equivalent to (1/2)(Fi+Fi+1) - 1/2A~(Ui+1-Ui)
		// Store values in a double array, to be transferred to the vector later
		double U_dif[3] = {U[0][i + 1] - U[0][i], U[1][i + 1] - U[1][i], U[2][i + 1] - U[2][i]};
		double Roe[3] = {0};
		matrix_vector(A_tild, U_dif, Roe);

		F[0][i] = 0.5 * (Fnorm[0][i] + Fnorm[0][i + 1]) - 0.5 * Roe[0];
		F[1][i] = 0.5 * (Fnorm[1][i] + Fnorm[1][i + 1]) - 0.5 * Roe[1];
		F[2][i] = 0.5 * (Fnorm[2][i] + Fnorm[2][i + 1]) - 0.5 * Roe[2];
	}
}
