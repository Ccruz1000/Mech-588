// Included libraries found in header file
#include "Roe.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Goated video - https://www.youtube.com/watch?v=F_PsFHvt8IU&ab_channel=P.S.Volpiani

// Read input values from input file
void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond, double &t_final)
{
	/*
		This function will read the input file, which will be used to select the 
		problem parameters. If they are not defined, or the input file is not defined
		properly, then default parameters will be assigned.
	*/

	// Assign default parameters if input file is not input properly
	if(input_file == "Null")
	{
		// Assign the default values to provide a contact discontinuity
		CFL = 0.95;
		nx = 201;
		xbegin = 0;
		xend = 1;
		dt = 2.5e-3;
		gamma = 1.4;
		init_cond[0] = 1.01; // left density
		init_cond[1] = 0; // left velocity
		init_cond[2] = 1; // left pressure
		init_cond[3] = 1; // right density
		init_cond[4] = 0; // right velocity
		init_cond[5] = 1; // right pressure
		t_final = 0.15;
	}

	// If input file is properly provided, try to read values from it
	else
	{
		std::ifstream input(input_file); // Read input file

		// Check if mesh is opened, return error otherwise
		if(input.is_open())
		{
			std::string var; // Store the variable that is being populated
			std::string val; // Store the value as a string

			/*
			The following will be used to ensure that each variable was assigned properly.
			If there is an error in the input file (for example missing :), then one of the lines
			will be read, but the variable will not be assigned the value. This will inform the 
			user if this happens
			*/
			int numval = 0; 
			int numline = 0;

			// Loop through each line of the input file
			while(input >> var >> val)
			{
				// Check which variable is being populated and assign value appropriately
				// Definitely better way to do this, but if it ain't broke don't fix it
				if (var == "CFL:")
				{
					CFL = std::stod(val);
					numval += 1;
				}

				else if (var == "nx:")
				{
					nx = std::stoi(val) + 2; // Add 2 for ghost cells
					numval += 1;					
				}

				else if (var == "xbegin:")
				{
					xbegin = std::stod(val);
					numval += 1;					
				}

				else if (var == "xend:")
				{
					xend = std::stod(val);	
					numval += 1;				
				}

				else if(var == "dt:")
				{
					dt = std::stod(val);
					numval += 1;					
				}

				else if (var == "gamma:")
				{
					gamma = std::stod(val);	
					numval += 1;				
				}

				else if (var == "rhol:")
				{
					init_cond[0] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "ul:")
				{
					init_cond[1] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "Pl:")
				{
					init_cond[2] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "rhor:")
				{
					init_cond[3] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "ur:")
				{
					init_cond[4] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "Pr:")
				{
					init_cond[5] = std::stod(val);	
					numval += 1;				
				}

				else if (var == "tfinal:")
				{
					t_final = std::stod(val);
					numval += 1;
				}

				numline += 1;
			}

			// Check that all variables were assigned
			if(numval != numline)
			{
				printf("\n****************************************************\n");
				printf("Warning! All values may not have been assigned.\nPlease check input file.\n");
				printf("\n****************************************************\n");
			}

			// Close input file once we are done with it
			input.close();
		}
		// Give warning that input file was not read properly
		else
		{
			printf("Input file not read successfully, please try again.\n");
			exit(1);
		}
	}
}

// Assign values from input, based on whether input file is given
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond, double &t_final)
{
	// Check if input file was supplied
	if(argc < 2)
	{
		printf("\nNo input file was provided. Please check.\nDefault parameters will be assigned to all input parameters.\n");
		
		std::string input_file = "Null";
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma, init_cond, t_final);
	}

	// Let user know that input was successfully provided
	else 
	{
		// Convert the inputs file from an array of char to a string
		std::string input_file = argv[1];

		// Inform user that the input file was successfully provided, and let them know which one was provided.
		// std::cout << "\nInput file provided successfully" << std::endl;
		// std::cout << "File given: " << input_file << std::endl;
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma, init_cond, t_final);
	}
}

// Assign initial conditions based on input values
void initial_conditions(std::array<std::vector<double>, 4> &Prim, int nx, double *init_cond, double gamma)
{
	int mid = nx/2; // Point where membrane is contained

	// Assign values to each point
	for(int i = 0; i < nx; i++)
	{
		// Left hand side if before mid
		if(i < mid)
		{
			Prim[0][i] = init_cond[0];  // Store density
			Prim[1][i] = init_cond[1];  // Store velocity
			Prim[2][i] = init_cond[2];  // Store pressure
		}
		// Right hand side after mid
		else
		{
			Prim[0][i] = init_cond[3];  // Store density
			Prim[1][i] = init_cond[4];  // Store velocity
			Prim[2][i] = init_cond[5];  // Store pressure
		}
		// Calculate energy from the assigned initial conditions
		Prim[3][i] = Prim[2][i] / (gamma - 1) + 0.5 * Prim[0][i] * Prim[1][i] * Prim[1][i]; // Store energy
		// Et = P/(gamma - 1) + rho * u^2 / 2
	}
}

// Save solution to file, which can be plotted later using python
void save_results(std::string fileName, std::array<std::vector<double>, 4> &Prim, double dx, int nx)
{
	const char *fileChar = fileName.c_str();
	FILE *solutionFile;
	solutionFile = fopen(fileChar, "w");

	// Print file header
	fprintf(solutionFile, "X rho u P E\n");

	// Print values
	for(int i = 1; i < nx - 1; i++)
	{
		fprintf	(solutionFile, "%14.12e %14.12e %14.12e %14.12e %14.12e\n", i * dx, Prim[0][i], Prim[1][i], Prim[2][i], Prim[3][i]);
	}
}

// The following will compute the product of 2 3x3 matrices
void matrix_mult3x3(double A[3][3], double B[3][3], double C[3][3])
{
	// Empty matrix to store results in
	double D[3][3] = {{0}};
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			for(int k = 0; k < 3; k++)
			{

				C[i][j] += A[i][k] * B[k][j];
			}
			// C[i][j] = D[i][j];
		}
}

// The following will compute a matrix vector multiplication of a 3x3 matrix, and a 3x1 vector
void matrix_vector(double A[3][3], double B[3], double C[3])
{
	double D[3] = {0};
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			D[i] += A[i][j] * B[j];
		}
		C[i] = D[i];
	}
}

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
	F[0][0] = 0; // UPDATE THIS TO BE PROPER VALUE
	F[1][0] = Prim[2][1]; // UPDATE THIS TO PROPER VALUE
	F[2][0] = 0; // UPDATE THIS TO PROPER VALUE

	F[0][nx - 1] = 0; // UPDATE THIS TO BE PROPER VALUE
	F[1][nx - 1] = Prim[2][nx - 2]; // UPDATE THIS TO PROPER VALUE
	F[2][nx - 1] = 0; // UPDATE THIS TO PROPER VALUE

	// Calculate flux at each cell
	for(int i = 1; i < nx - 1; i++)
	{
		double du_dv[3][3] = {{1, 0, 0},
		                      {u_tild[i], U_tild[0][i], 0}, 
		                      {0.5 * u_tild[i] * u_tild[i], U_tild[0][i] * u_tild[i], 1.0 / (gamma - 1)}};
		/*
				|1     0     0              |
		du_dv = |u     rho   0              |
				|u^2/2 rho*u 1 / (gamma - 1)|
		*/

		double XR[3][3] = {{1.0, U_tild[0][i] / c[i], -U_tild[0][i] / c[i]},
					       {0, 1, 1},
					       {0, U_tild[0][i] * c[i], -U_tild[0][i] * c[i]}};
		/*
				|u     rho/c       -rho/c    |
		XR =    |0     1           1         |
				|0     rho*c       -rho*c    |
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

		double XL[3][3] = {{1, 0, -1/(c[i]*c[i])},
					       {0, 1/2, 1/(2*U_tild[0][i] * c[i])},
					       {0, 1/2, -1/(2*U_tild[0][i] * c[i])}};

		/*
				|1     0             -1/c^2      |
		XL =    |0     1/2           1/(2(rho*c) |
				|0     1/2           -1/(2(rho*c)|
		*/

		double dv_du[3][3] = {{1, 0, 0},
					          {-u_tild[i] / U_tild[0][i], 1 / U_tild[0][i], 0},
					          {(gamma - 1) * 0.5 * u_tild[i] *u_tild[i], -(gamma - 1) * u_tild[i], gamma - 1}};

		/*
				|1                0              0        |
		dv_du = |-u/rho           1/rho          0        |
				|(gamma-1)u^2/2   -(gamma-1)u    gamma - 1|
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

		F[0][i] = 0.5 * (Fnorm[0][i] + Fnorm[0][i + 1]) + 0.5 * Roe[0];
		F[1][i] = 0.5 * (Fnorm[1][i] + Fnorm[1][i + 1]) + 0.5 * Roe[1];
		F[2][i] = 0.5 * (Fnorm[2][i] + Fnorm[2][i + 1]) + 0.5 * Roe[2];
	}
}

// // The following will compute the flux using the Steger Warming scheme
// void Steger_Warming(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, std::vector<double> c, int nx, double *init_cond)
// {
// 	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
// 	std::vector<double> htot = std::vector<double>(nx); // Store enthalpy
// 	// Compute the primitive variables from the conserved variables
// 	for(int i = 1; i < nx - 1; i++)
// 	{
// 		Prim[0][i] = U[0][i]; // Density
// 		Prim[1][i] = U[1][i] / U[0][i]; // Velocity
// 		Prim[3][i] = U[2][i]; // Energy
// 		Prim[2][i] = 0.5 * (1 - gamma) * U[1][i] * U[1][i] / U[0][i] + (gamma - 1) * U[2][i];  // Pressure (P = 0.5 * (1 - gamma) * U2^2/U1 + (gamma - 1) * U3
// 		c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
// 		// calculate enthalpy
		
// 		if(Prim[0][i] < 0)
// 		{
// 			printf("Warning! Negative density found.\n");
// 		}

// 		if(Prim[2][i] < 0)
// 		{
// 			printf("Warning! Negative pressure found.\n");
// 		}
// 	}

// 	// Matrices to calculate Roes flux
// 	double lambda_plus[3][3] = {0};
// 	double lambda_minus[3][3] = {0};
// 	double lam1, lam2, lam3;
// 	double A_hat[3][3] = {0};

// 	// Boundary condition
// 	F[0][0] = 0; // UPDATE THIS TO BE PROPER VALUE
// 	F[1][0] = Prim[2][1]; // UPDATE THIS TO PROPER VALUE
// 	F[2][0] = 0; // UPDATE THIS TO PROPER VALUE

// 	F[0][nx - 1] = 0; // UPDATE THIS TO BE PROPER VALUE
// 	F[1][nx - 1] = Prim[2][nx - 2]; // UPDATE THIS TO PROPER VALUE
// 	F[2][nx - 1] = 0; // UPDATE THIS TO PROPER VALUE

// 	// Calculate flux at each cell
// 	for(int i = 1; i < nx - 1; i++)
// 	{
// 		double du_dv[3][3] = {{1, 0, 0},
// 		                      {Prim[1][i], Prim[0][i], 0}, 
// 		                      {0.5 * Prim[1][i] * Prim[1][i], Prim[0][i] * Prim[1][i], 1.0 / (gamma - 1)}};
// 		/*
// 				|1     0     0        |
// 		du_dv = |u     rho   0        |
// 				|u^2/2 rho*u gamma - 1|
// 		*/
// 		double XR[3][3] = {{1.0, Prim[0][i] / c[i], -Prim[0][i] / c[i]},
// 					       {0, 1, 1},
// 					       {0, Prim[0][i] * c[i], -Prim[0][i] * c[i]}};
// 		/*
// 				|u     c(u+c)/rho  c(c-u)/rho|
// 		XR =    |0     1           1         |
// 				|0     rho*c       -rho*c    |
// 		*/
// 		 // Values to store the 3 individual eigenvalues

// 		lam1 = Prim[1][i];
// 		lam2 = Prim[1][i] + c[i];
// 		lam3 = Prim[1][i] - c[i];

// 		double lambda_plus[3][3] = {{lam1 + std::sqrt(lam1 * lam1 + epsilon), 0, 0},
// 					               {0, lam2 + std::sqrt(lam2 * lam2 + epsilon), 0},
// 					               {0, 0, lam3 + std::sqrt(lam3 * lam3 + epsilon)}};

// 		double lambda_minus[3][3] = {{lam1 - std::sqrt(lam1 * lam1 + epsilon), 0, 0},
// 					               {0, lam2 - std::sqrt(lam2 * lam2 + epsilon), 0},
// 					               {0, 0, lam3 - std::sqrt(lam3 * lam3 + epsilon)}};

// 		/*
// 				         |lam1 +/- abs(lam1+eps)     0                   0                 |
// 		lambda(+/-) =    |0                      lam2 +/- abs(lam2+eps)  0                 |
// 				         |0                      0                   lam3 +/- abs(lam3+eps)|
// 		*/

// 		double XL[3][3] = {{1, 0, -1/(c[i]*c[i])},
// 					       {0, 1/2, 1/(2*Prim[0][i] * c[i])},
// 					       {0, 1/2, -1/(2*Prim[0][i] * c[i])}};

// 		/*
// 				|1     0             -1/c^2      |
// 		XL =    |0     1/2           1/(2(rho*c) |
// 				|0     1/2           -1/(2(rho*c)|
// 		*/

// 		double dv_du[3][3] = {{1, 0, 0},
// 					          {-u_tild[i] / U_tild[0][i], 1 / U_tild[0][i], 0},
// 					          {(gamma - 1) * 0.5 * u_tild[i] *u_tild[i], -(gamma - 1) * u_tild[i], gamma -1}};

// 		/*
// 				|1                0              0        |
// 		dv_du = |-u/rho           1/rho          0        |
// 				|(gamma-1)u^2/2   -(gamma-1)u    gamma - 1|
// 		*/       

// 		// Now that all the matrices are defined, we can calculate Ahat
// 		// Ahat = (du_dv) * XR * lambda(+/-) * XL * (dv_du)
// 		// Intermediate matrices to store values
// 		double inter1[3][3] = {0};
// 		double inter2[3][3] = {0};
// 		// We will do lamdaplus first
// 		// A = du_dv*XR
// 		matrix_mult3x3(du_dv, XR, inter1);
// 		// B = A * lamda_abs
// 		matrix_mult3x3(inter1, lambda_plus, inter2);
// 		// Zero inter1, and then C = B * XL
// 		for(int k = 0; k < 3; k++)
// 			for(int j = 0; j < 3; j++)
// 			{
// 				inter1[k][j] = 0; 
// 			}
// 		matrix_mult3x3(inter2, XL, inter1);
// 		// Finally, A~ = C * dv_du
// 		// matrix_mult3x3(inter1, dv_du, A_tild);

// 		// Now, we know that the flux is equivalent to (1/2)(Fi+Fi+1) - 1/2A~(Ui+1-Ui)
// 		// Store values in a double array, to be transferred to the vector later
// 		double U_dif[3] = {U[0][i + 1] - U[0][i], U[1][i + 1] - U[1][i], U[2][i + 1] - U[2][i]};
// 		double Roe[3] = {0};
// 		matrix_vector(A_tild, U_dif, Roe);

// 		// Calculate total flux
// 		// Boundary conditions
// 		// if(i == 0)
// 		// {
// 		// 	F[0][i] = 0;
// 		// 	F[1][i] = init_cond[2]; // UPDATE THIS TO BE PROPER VALUE
// 		// 	F[2][i] = 0;
// 		// }

// 		// else if (i == nx - 1)
// 		// {
// 		// 	F[0][i] = 0;
// 		// 	F[1][i] = init_cond[5]; // 
// 		// 	F[2][i] = 0;
// 		// }

// 		// else
// 		// {
// 		F[0][i] = 0.5 * (Fnorm[0][i] + Fnorm[0][i + 1]) - 0.5 * Roe[0];
// 		F[1][i] = 0.5 * (Fnorm[1][i] + Fnorm[1][i + 1]) - 0.5 * Roe[1];
// 		F[2][i] = 0.5 * (Fnorm[2][i] + Fnorm[2][i + 1]) - 0.5 * Roe[2];
// 		// }
// 		// F[0][i] = 0.5 * (Fnorm[0][i] + Fnorm[0][i + 1]) - 0.5 * Roe[0];
// 		// F[1][i] = 0.5 * (Fnorm[1][i] + Fnorm[1][i + 1]) - 0.5 * Roe[1];
// 		// F[2][i] = 0.5 * (Fnorm[2][i] + Fnorm[2][i + 1]) - 0.5 * Roe[2];
// 	}
// }

// Timestepping schemes
void forward_euler(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> &F, double dt, double dx, int nx)
{
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = U[0][i] + (dt / dx) * (F[0][i] - F[0][i - 1]);
		U[1][i] = U[1][i] + (dt / dx) * (F[1][i] - F[1][i - 1]);
		U[2][i] = U[2][i] + (dt / dx) * (F[2][i] - F[2][i - 1]);
	}
}

