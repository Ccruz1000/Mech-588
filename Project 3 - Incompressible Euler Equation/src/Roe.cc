// Included libraries found in header file
#include "Roe.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Goated video - https://www.youtube.com/watch?v=F_PsFHvt8IU&ab_channel=P.S.Volpiani

// Read input values from input file
void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond)
{
	/*
		This function will read the input file, which will be used to select the 
		problem parameters. If they are not defined, or the input file is not defined
		properly, then default parameters will be assigned.
	*/

	// Assign default parameters if input file is not input properly
	if(input_file == "Null")
	{
		CFL = 0.95;
		nx = 200;
		xbegin = 0;
		xend = 1;
		dt = 2.5e-2;
		gamma = 1.4;
		init_cond[0] = 6; // left density
		init_cond[1] = 0; // left velocity
		init_cond[2] = 12; // left pressure
		init_cond[3] = 1; // right density
		init_cond[4] = 0; // right velocity
		init_cond[5] = 1; // right pressure
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
					nx = std::stoi(val);
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
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma, double *init_cond)
{
	// Check if input file was supplied
	if(argc < 2)
	{
		printf("\nNo input file was provided. Please check.\nDefault parameters will be assigned to all input parameters.\n");
		
		std::string input_file = "Null";
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma, init_cond);
	}

	// Let user know that input was successfully provided
	else 
	{
		// Convert the inputs file from an array of char to a string
		std::string input_file = argv[1];

		// Inform user that the input file was successfully provided, and let them know which one was provided.
		std::cout << "\nInput file provided successfully" << std::endl;
		std::cout << "File given: " << input_file << std::endl;
		read_input(input_file, CFL, nx, xbegin, xend, dt, gamma, init_cond);
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
		Prim[3][i] = Prim[2][i] / ((gamma - 1) * Prim[0][i]) + (Prim[1][i] * Prim[1][i] / 2); // Store energy
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
	for(int i = 0; i < nx; i++)
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
void Roe_flux(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, double dx, double gamma, std::vector<double> c, int nx)
{
	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	std::vector<double> htot = std::vector<double>(nx);
	// Compute the primitive variables from the conserved variables
	for(int i = 0; i < nx; i++)
	{
		Prim[0][i] = U[0][i]; // Density
		Prim[1][i] = U[1][i] / U[0][i]; // Velocity
		Prim[3][i] = U[2][i] / U[0][i]; // Energy
		Prim[2][i] = (gamma - 1) * Prim[0][i] * (Prim[3][i] - 0.5 * Prim[2][i] * Prim[2][i]);  // Pressure (P = (gamma - 1) * rho * (E - 0.5 u^2))
		// calculate enthalpy
		htot[i] = gamma / (gamma - 1) * Prim[2][i] / Prim[0][i] + 0.5 * Prim[1][i] * Prim[1][i];
	}

	// Initialize Roe averaged values and Roe flux vector
	std::array<std::vector<double>, 3> Phi = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	double R;  // Quantity inside square root R_j+1/2
	double rho_hat; 
	double u_hat;
	double H_hat;
	double a_hat;

	// Initialize vectors to store required values for flux calculation
	double U_diff[3] = {0};
	double P[3][3] = {0};
	double lamda_hat[3][3] = {0};
	double Pinv[3][3] = {0};
	double A_hat[3][3] = {0};
	double alpha1, alpha2; // Variables used to compute P and Pinv

	// // Loop through cells
	// for(int i = 0; i < nx; i++)
	// {
	// 	// Compute Roe averages
	// 	R = sqrt(Prim[0][i + 1] / Prim[0][i]);
	// 	rho_hat = R*Prim[0][i];
	// 	u_hat = R * Prim[1][i + 1] + Prim[1][i] / (1 + R);
	// 	H_hat[i] = R*htot[i + 1] * htot / (1 + R);
	// 	a_hat[i] = sqrt((gamma - 1) * (H_hat - 0.5 * Prim[1][i] * Prim[1][i]));

	// 	// Compute alphas o find P and Pinv
	// 	alpha1 = (gamma - 1) * u_hat * u_hat / (2 * a_hat * a_hat);
	// 	alpha2 = (gamma - 1) / (a_hat * a_hat);

	// 	// Compute values for flux calculations 
	// 	U_diff[0] = = U[0][i + 1] - U[0][i];
	// 	U_diff[1] = U[1][i + 1] - U[1][i];
	// 	U_diff[2] = U[2][i + 1] - U[2][i];

	// 	P = {{1, 1, 1}, 
	// 	    {(u_hat - a_hat), u_hat, u_hat + a_hat}, 
	// 	    {H_hat - u_hat * a_hat, 0.5 * u_hat * u_hat, H_hat + a_hat * u_hat}};

	// 	Pinv = {{0.5 * (alpha1 + u_hat / a_hat), -0.5 * (alpha2 * u_hat + 1 / a_hat), alpha2 / 2}, 
	// 	        {1 - alpha1, alpha2 * u_hat, -alpha2},
	// 	        {0.5 * (alpha1 - u_hat / a_hat), - 0.5 * (alpha2 * u_hat - 1 / a_hat), alpha2 / 2}};

	// 	lambda_hat = {{std::fabs(u_hat - a_hat + epsilon), 0, 0}, 
	// 				  {0, std::fabs(u_hat), 0},
	// 				  {0, 0, std::fabs(u_hat + a_hat)}};

	// 	// Compute A_hat Jacobian
	// 	double intermediate[3][3] = {0}; // Intermediate matrix to store values for matrix multiplication
	// 	matrix_mult3x3(P, lamda_hat, intermediate);
	// 	matrix_mult3x3(intermediate, Pinv, A_hat);



	// }


}

