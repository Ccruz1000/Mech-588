// Included libraries found in header file
#include "functions.h"

/**************************************************************************************
// Compressible Flow Programming Assignment General Function File
// This file will contain all general purpose functions and classes. This includes 
// things such as saving files, matrix multiplication, and reading input data.
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

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
	int mid = nx / 2; // Point where membrane is contained

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

				D[i][j] += A[i][k] * B[k][j];
			}
			C[i][j] = D[i][j];
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

// This function will compute the primitive variables using the conserved solution vector
void calc_Prim(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, double gamma, double nx)
{
		// Compute the primitive variables from the conserved variables
	for(int i = 1; i < nx - 1; i++)
	{
		Prim[0][i] = U[0][i]; // Density
		Prim[1][i] = U[1][i] / U[0][i]; // Velocity
		Prim[3][i] = U[2][i]; // Energy
		Prim[2][i] = 0.5 * (1 - gamma) * U[1][i] * U[1][i] / U[0][i] + (gamma - 1) * U[2][i];  // Pressure (P = 0.5 * (1 - gamma) * U2^2/U1 + (gamma - 1) * U3
		
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
}

// This function will set the boundary conditions of the problem
void calc_bound(std::array<std::vector<double>, 3> &F, const std::array<std::vector<double>, 4> Prim, double nx)
{
	// Boundary condition
	F[0][0] = 0; 
	F[1][0] = Prim[2][1]; 
	F[2][0] = 0; 

	F[0][nx - 1] = 0; 
	F[1][nx - 1] = Prim[2][nx - 2]; 
	F[2][nx - 1] = 0;
}

// This function will assemble the matrices required to calculate the Jacobian matrix 
void assemble_matrices(double rho, double u, double c, double gamma, double du_dv[3][3], double XR[3][3], double XL[3][3], double dv_du[3][3])
{
	// Temporary matrices to allow us to store the expected matrix. We can then assign the values to the 
	// given matrices to return the assembled matrix. This helps with headaches about matrix initialization
	// Assemble matrices to perform calculation of Jacobian matrix
		double xrtemp[3][3] = { {1.0, rho / c, -rho / c},
					        {0, 1.0, 1.0},
					        {0, rho * c, -rho * c} };
		/*
				|1     rho/c       -rho/c    |
		XR =    |0     1           1         |
				|0     rho*c       -rho*c    |
		*/

		double xltemp[3][3] = { {1, 0, -1 / (c * c)},
					        {0, 1.0 / 2.0, 1.0 / (2.0 * rho * c)},
					        {0, 1.0 / 2.0, -1.0 / (2.0 * rho * c)} };
		        
		/*
				|1     0             -1/c^2      |
		XL =    |0     1/2           1/(2(rho*c) |
				|0     1/2           -1/(2(rho*c)|
		*/ 

		double du_dvtemp[3][3] = {{1, 0, 0},
		                      {u, rho, 0}, 
		                      {0.5 * u * u, rho * u, 1.0 / (gamma - 1)}};
		/*
				|1     0     0              |
		du_dv = |u     rho   0              |
				|u^2/2 rho*u 1 / (gamma - 1)|
		*/

		double dv_dutemp[3][3] = {{1, 0, 0},
					          {-u / rho, 1 / rho, 0},
					          {(gamma - 1) * 0.5 * u * u, -(gamma - 1) * u, gamma - 1}};

		/*
				|1                0              0        |
		dv_du = |-u/rho           1/rho          0        |
				|(gamma-1)u^2/2   -(gamma-1)u    gamma - 1|
		*/

		// Now that all of the matrices are assembled, we will assign the values to the given matrices to be returned
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
			{
				XR[i][j] = xrtemp[i][j];
				XL[i][j] = xltemp[i][j];
				du_dv[i][j] = du_dvtemp[i][j];
				dv_du[i][j] = dv_dutemp[i][j];
			}
}

// This function will assemble lambda given the appropriate Eigenvalue functions (for example, in Roe it will need to be given the absolute value of each eigenvalue)
// but for Steger it will need to give iether the positive or negative eigenvalue 
void assemble_lambda(double lam1, double lam2, double lam3, double lambda[3][3])
{
	// We will start by overwriting each value of lambda to be 0
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			lambda[i][j] = 0.0;
		}

	// Now we will overwrite the diagnol to be the values provided
	lambda[0][0] = lam1;
	lambda[1][1] = lam2;
	lambda[2][2] = lam3;
}

// This function will take in the matrices for du_dv, XR, lambda, XL, and dv_du to calculate the A jacobian matrix
void calc_Jacobian(double du_dv[3][3], double XR[3][3], double lambda[3][3], double XL[3][3], double dv_du[3][3], double A[3][3])
{
	// Create intermediate matrices to perform matrix multiplication
	double inter1[3][3] = {0};
	double inter2[3][3] = {0};

	// Start by multiplying du_dv and XR
	matrix_mult3x3(du_dv, XR, inter1);

	// Next multiply (du_dv*XR) and lambda
	matrix_mult3x3(inter1, lambda, inter2);

	// Zero out inter1
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			inter1[i][j] = 0.0;
		}

	// Multiply (du_dv*XR*lambda) and XL
	matrix_mult3x3(inter2, XL, inter1);

	// Finally multiply (du_dv*XR*lambda*XL) and dv_du to get A
	matrix_mult3x3(inter1, dv_du, A);
}

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