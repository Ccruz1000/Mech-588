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
void Roe_flux(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double dx, double gamma, std::vector<double> c, int nx)
{
	double epsilon = 1e-8; // Value to smooth out absolute values of eigenvectors near 0
	std::vector<double> htot = std::vector<double>(nx); // Store enthalpy
	// Compute the primitive variables from the conserved variables
	for(int i = 0; i < nx; i++)
	{
		Prim[0][i] = U[0][i]; // Density
		Prim[1][i] = U[1][i] / U[0][i]; // Velocity
		Prim[3][i] = U[2][i]; // Energy
		Prim[2][i] = 0.5 * (1 - gamma) * U[1] * U[1] / U[0] + (gamma - 1) * U[2];  // Pressure (P = 0.5 * (1 - gamma) * U2^2/U1 + (gamma - 1) * U3
		// calculate enthalpy
		htot[i] = gamma / (gamma - 1) * Prim[2][i] / Prim[0][i] + 0.5 * Prim[1][i] * Prim[1][i];
	}

	// Store normal flux
	std::array<std::vector<double>, 3> Fnorm = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)}; 

	for(int i = 0; i < nx; i++)
	{
		F[0][i] = Prim[0][i] * Prim[1][i]; // rho*u
		F[1][i] = Prim[0][i] * Prim[1][i] * Prim[1][i] + Prim[2][i]; // rho*u^2 + P
		F[2][i] = Prim[1][i] * (Prim[2][i] + Prim[3][i]); // u*(P+Et)
	}

	// Store Roe averaged values
	std::array<std::vector<double>, 3> U_tild = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	double rho_tild, u_tild, h_tild, h_left, h_right;
	for(int i = 0; i < nx; i++)
	{	
		// Calculate rho averaged density and velocity
		rho_tild = std::sqrt(Prim[0][i + 1] + Prim[0][i]); 
		u_tild = (Prim[0][i + 1] * Prim[1][i + 1] + rho_tild * Prim[1][i]) / (Prim[0][i + 1] + rho_tild);

		// Calculate enthalpy
		h_left = Prim[3][i] + Prim[2][i] / Prim[0][i];
		h_right = Prim[3][i + 1] + Prim[2][i + 1] / Prim[0][i + 1];
		h_tild = (Prim[0][i + 1] * h_right + rho_tild * h_left) / (Prim[0][i + 1] + rho_tild);
		
		// Populate Roe averaged solution vector
		U_tild[0][i] = rho_tild;
		U_tild[1][i] = rho_tild * u_tild;
		U_tild[2][i] = h_tild;
	} 

	// Matrices to calculate Roes flux
	double A_tild[3][3] = {0};
	double du_dv[3][3] = {0};
	double XR[3][3] = {0};
	double lambda_abs[3][3] = {0};
	double dv_du[3][3] = {0};

	// Calculate flux at each cell
	for(int i = 0; i < nx; i++)
	{
		du_dv = {{1, 0, 0},
		         {U_tild[1][i] / U_tild[0][i], U_tild[0][i], 0}, 
		         {0.5 * U_tild[1][i] * U_tild[1][i], U_tild[0][i] * U_tild[1][i], gamma - 1}};
		/*
				|1     0     0        |
		du_dv = |u     rho   0        |
				|u^2/2 rho*u gamma - 1|
		*/
		XR = {{U_tild[1][i] / U_tild[0][i], c[i] * ((U_tild[1][i] / U_tild[0][i]) + c[i]) / U_tild[0][i], c[i] * ((c[i] - U_tild[1][i] / U_tild[0][i])) / U_tild[0][i]},
		      {0, 1, 1},
		      {0, U_tild[0][i] * c[i], -U_tild[0][i] * c[i]}}
		/*
				|u     c(u+c)/rho  c(c-u)/rho|
		XR =    |0     1           1         |
				|0     rho*c       -rho*c    |
		*/

		lambda_abs = 0.5 * {{std::sqrt(epsilon + (0.5 * (U_tild[1][i] / U_tild[0][i])) * (0.5 * (U_tild[1][i + 1] / U_tild[0][i + 1]))), 0, 0},
		                    {0, std::sqrt(epsilon + (0.5 * ((U_tild[1][i] / U_tild[0][i]) + c[i])) * (0.5 * ((U_tild[1][i + 1] / U_tild[0][i + 1]) + c[i + 1]))), 0},
		                    {0, 0, std::sqrt(epsilon + (0.5 * ((U_tild[1][i] / U_tild[0][i]) - c[i])) * (0.5 * ((U_tild[1][i + 1] / U_tild[0][i + 1]) - c[i + 1])))}};    
	}

}

