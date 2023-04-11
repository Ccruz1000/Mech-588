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
void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &gamma, double *init_cond, double &t_final, int &num_save, int &timescheme, int &limitscheme)
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
		nx = 200; // number of cells
		xbegin = 0;
		xend = 1;
		gamma = 1.4;
		init_cond[0] = 1.01; // left density
		init_cond[1] = 0; // left velocity
		init_cond[2] = 1.0; // left pressure
		init_cond[3] = 1.0; // right density
		init_cond[4] = 0; // right velocity
		init_cond[5] = 1.0; // right pressure
		t_final = 0.15;
		num_save = 20;
		timescheme = 2; // Midpoint default
		limitscheme = 0; // superbee default
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

				if (var == "num_save:")
				{
					num_save = std::stoi(val);
					numval += 1;
				}

				if (var == "timescheme:")
				{
					timescheme = std::stoi(val);
					numval += 1;
				}

				if (var == "limitscheme:")
				{
					limitscheme = std::stoi(val);
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
void assign_input(int argc, char* argv[], double &CFL, int &nx, double &xbegin, double &xend, double &gamma, double *init_cond, double &t_final, int &num_save, int &timescheme, int &limitscheme)
{
	// Check if input file was supplied
	if(argc < 2)
	{
		printf("\nNo input file was provided. Please check.\nDefault parameters will be assigned to all input parameters.\n");
		
		std::string input_file = "Null";
		read_input(input_file, CFL, nx, xbegin, xend, gamma, init_cond, t_final, num_save, timescheme, limitscheme);
	}

	// Let user know that input was successfully provided
	else 
	{
		// Convert the inputs file from an array of char to a string
		std::string input_file = argv[1];

		// Inform user that the input file was successfully provided, and let them know which one was provided.
		// std::cout << "\nInput file provided successfully" << std::endl;
		// std::cout << "File given: " << input_file << std::endl;
		read_input(input_file, CFL, nx, xbegin, xend, gamma, init_cond, t_final, num_save, timescheme, limitscheme);
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

// The following function will run the simulation
void run_sim(int argc, char* argv[], int scheme)
{
	// Declare run parameters
	double CFL, xbegin, xend, dx, dt, gamma, t_final;
	int nx, num_save, timescheme, limitscheme;
	std::string title; // Save initial and final data
	std::string savepath; // Save timestep data

	// We will store the initial condition values in an array. The first 3 entries will be for the left side
	// and will represent rho, u, and P in that order. The next 3 entries will be for the right side, and will 
	// represent rho, u, and P in that order

	double init_cond[6] = {0};

	// Read input file, or provide default values if none given
	assign_input(argc, argv, CFL, nx, xbegin, xend, gamma, init_cond, t_final, num_save, timescheme, limitscheme);

	// Calculate the spatial step
	dx = (xend - xbegin) / nx;

	// Create array to store primitives. Prim[0][i] is rho, Prim[1][i] is velocity, Prim[2][i] is pressure, and Prim[3][i] is energy
	std::array<std::vector<double>, 4> Prim = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};

	// Assign initial condition values
	initial_conditions(Prim, nx, init_cond, gamma);

	// Calculate maximum initial timestep
	calc_timestep(Prim, dt, gamma, nx, CFL, dx);

	// Create and populate array to store conserved variables and flux. 
	std::array<std::vector<double>, 3> U = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};  // Solution vector
	std::array<std::vector<double>, 3> FH = {std::vector<double>(nx + 1), std::vector<double>(nx + 1), std::vector<double>(nx + 1)}; // Flux vector for high order scheme
	std::array<std::vector<double>, 3> FL = {std::vector<double>(nx + 1), std::vector<double>(nx + 1), std::vector<double>(nx + 1)}; // Flux vector for low order scheme
	std::array<std::vector<double>, 3> F = {std::vector<double>(nx + 1), std::vector<double>(nx + 1), std::vector<double>(nx + 1)}; // Limited flux vector

	// Initialize vector to store limiting function
	std::array<std::vector<double>, 3> phi = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};

	// Populate solution using primitives
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = Prim[0][i];  // rho
		U[1][i] = Prim[0][i] * Prim[1][i];  // rho * u
		U[2][i] = Prim[3][i]; // Et
	}

	// Select title to save files with
	switch(scheme)
	{
		case 0:
		{
			printf("\n*******************************************************\n");
			printf("Roe's flux difference splitting scheme slected.\n");
			title = "_Roe_";
			savepath = "RoeTimesteps/inter_";
			break;
		}

		case 1:
		{
			printf("\n*******************************************************\n");
			printf("Steger-Warming flux vector splitting scheme selected.\n");
			title = "_Steger_";
			savepath = "StegerTimesteps/inter_";
			break;
		}

		case 2:
		{
			printf("\n*******************************************************\n");
			printf("Lax-Friedrich flux scheme selected.\n");
			title = "_Lax_";
			savepath = "LaxTimesteps/inter_";
			break;
		}

		default:
		{
			printf("Invalid scheme selected.\nFor Roe enter 0.\nFor Steger-Warming enter 1.\nFor Lax-Friedrich enter 2.\nExiting.\n");
			exit(1);
		}
	}

	switch(timescheme)
	{
		case 0:
		{
			printf("Forward Euler timestepping scheme selected.\n");
			printf("*******************************************************\n\n");
			break;
		}

		case 1:
		{
			printf("RK2 timestepping scheme selected.\n");
			printf("*******************************************************\n\n");
			break;
		}

		case 2:	
		{
			printf("Midpoint timestepping scheme selected.\n");
			printf("*******************************************************\n\n");
			break;
		}

		default:
		{
			printf("Invalid timescheme selected.\nFor forward Euler, use 0.\nFor RK2, use 1.\nFor midpoint use 2.\nExiting\n");
			printf("*******************************************************\n");
			exit(1);
		}
	}


	switch(limitscheme)
	{
		case 0:
		{
			printf("Superbee flux limiting scheme selected.\n");
			printf("*******************************************************\n\n");
			break;
		}

		case 1:
		{
			printf("Van-Leer flux limiting scheme selected.\n");
			printf("*******************************************************\n\n");
			break;
		}

		default:
		{
			printf("Invalid flux limiting scheme selected.\nFor Superbee, use 0.\nFor Van-Leer, use 1.\nExiting\n");
			printf("*******************************************************\n");
			exit(1);
		}
	}

	// Calculate initial flux
	switch(scheme)
	{
		case 0:
		{	
			// Start by calculating flux with higher order scheme
			// lax_wendroff(U, Prim, FH, gamma, dt, dx, nx);
			// Next calculate flux with lower order scheme
			Roe_flux(U, Prim, FL, dx, gamma, nx);
			// // Finally apply the limiter function to the flux to pass to the timestepping
			// calc_limiter(U, phi, nx, limitscheme);
			// apply_limiter(F, FL, FH, phi, nx);
			break;
		}

		case 1:
		{
			// // Start by calculating flux with higher order scheme
			// lax_wendroff(U, Prim, FH, gamma, dt, dx, nx);
			// // Next calculate flux with lower order scheme
			Steger_warming(U, Prim, F, dx, gamma, nx);
			// Finally apply the limiter function to the flux to pass to the timestepping
			// calc_limiter(U, phi, nx, limitscheme);
			// apply_limiter(F, FL, FH, phi, nx);
			break;
		}

		case 2:
		{
			lax_wendroff(U, Prim, F, gamma, dt, dx, nx);
			break;
		}
	}

	// Save initial condition results
	std::string fileName = "Results/Init" + title + "cond.txt";
	save_results(fileName, Prim, dx, nx);
	
	// Begin timestepping
	double t_curr = 0;
	int itime = 0;
	while(t_curr < t_final)
	{
		switch(timescheme)
		{
			case 0:
			{
				forward_euler(U, F, dt, dx, nx);
				break;
			}

			case 1:
			{
				RK2(U, F, dt, dx, nx, gamma, scheme);
				break;
			}

			case 2:
			{
				midpoint(U, F, dt, dx, nx, gamma, scheme);
				break;
			}
		}

		switch(scheme)
		{
			case 0:
			{	
			// Start by calculating flux with higher order scheme
			// lax_wendroff(U, Prim, FH, gamma, dt, dx, nx);
			// Next calculate flux with lower order scheme
			Roe_flux(U, Prim, FL, dx, gamma, nx);
			// // Finally apply the limiter function to the flux to pass to the timestepping
			// calc_limiter(U, phi, nx, limitscheme);
			// apply_limiter(F, FL, FH, phi, nx);
			break;
			}

		case 1:
			{
			// // Start by calculating flux with higher order scheme
			// lax_wendroff(U, Prim, FH, gamma, dt, dx, nx);
			// // Next calculate flux with lower order scheme
			Steger_warming(U, Prim, F, dx, gamma, nx);
			// Finally apply the limiter function to the flux to pass to the timestepping
			// calc_limiter(U, phi, nx, limitscheme);
			// apply_limiter(F, FL, FH, phi, nx);
			break;
			}

		case 2:
			{
			lax_wendroff(U, Prim, F, gamma, dt, dx, nx);
			break;
			}
		}


		// Update primitives
		calc_Prim(U, Prim, gamma, nx);

		// Calculate maximum stable timestep
		calc_timestep(Prim, dt, gamma, nx, CFL, dx);


		itime += 1;
		t_curr += dt;
		std::string inter_results;
		inter_results = savepath + std::to_string(itime) + ".txt";
		if (itime % num_save == 0)
		{
			save_results(inter_results, Prim, dx, nx);
		} 
		printf("Timestep: %i Time: %f\n", itime, t_curr);
	}

	std::string final_results = "Results/Final" + title + "cond.txt";
	save_results(final_results, Prim, dx, nx);
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
void matrix_mult3x3(double A[3][3], double B[3][3], double (&C)[3][3])
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

// Format printing a matrix out
void print_matrix(double A[3][3])
{
	printf("\n");
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			printf("%16.16e ", A[i][j]); 
		}
		std::cout << std::endl;
	}
	printf("\n");
}


// The following will compute a matrix vector multiplication of a 3x3 matrix, and a 3x1 vector
void matrix_vector(double A[3][3], double B[3], double (&C)[3])
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
	for(int i = 0; i < nx; i++)
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
void calc_bound(std::array<std::vector<double>, 3> &F, const std::array<std::vector<double>, 4> &Prim, double nx)
{
	// Boundary condition
	F[0][0] = 0; 
	F[1][0] = Prim[2][0]; 
	F[2][0] = 0; 

	F[0][nx] = 0; 
	F[1][nx] = Prim[2][nx - 1]; 
	F[2][nx] = 0;
}

// This function will assemble the matrices required to calculate the Jacobian matrix 
void assemble_matrices(double rho, double u, double c, double gamma, double (&du_dv)[3][3], double (&XR)[3][3], double (&XL)[3][3], double (&dv_du)[3][3])
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
void assemble_lambda(double lam1, double lam2, double lam3, double (&lambda)[3][3])
{
	// We will start by overwriting each value of lambda to be 0
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			lambda[i][j] = 0.0;
		}

	// Now we will overwrite the diagonal to be the values provided
	lambda[0][0] = lam1;
	lambda[1][1] = lam2;
	lambda[2][2] = lam3;
}

// This function will take in the matrices for du_dv, XR, lambda, XL, and dv_du to calculate the A jacobian matrix
void calc_Jacobian(double du_dv[3][3], double XR[3][3], double lambda[3][3], double XL[3][3], double dv_du[3][3], double (&A)[3][3])
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

// This function will calculate the maximum stable timestep given the primitive variables
void calc_timestep(std::array<std::vector<double>, 4> Prim, double &dt, double gamma, int nx, double CFL, double dx)
{
	// Calculate maximum stable timestep
	// The below is used to help find the maximum of the absolute value of u + c, since we do
	// not know c at every point
	std::vector<double> tot_u = std::vector<double>(nx); 
	for(int i = 0; i < nx; i++)
	{
		tot_u[i] = std::fabs(Prim[1][i] + std::sqrt(gamma * Prim[2][i] / Prim[0][i]));
	}
	// find maximum velocity used to find maximum stable timestep
	double maxvel = *std::max_element(tot_u.begin(), tot_u.end());
	dt = CFL * dx / maxvel;
}


// The following function will compute the flux using the 2 step Lax-Friedrichs scheme
void lax_wendroff(std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 4> &Prim, std::array<std::vector<double>, 3> &F, double gamma, double dt, double dx, double nx)
{
	/*
	 This function will compute the flux at the cell faces using the second order Lax-Friedrichs scheme
	 The flux for this scheme is given as follows:
	 F_(i+1/2) = 0.5 * (F(U_i) + F(U_(i + 1)) - (dx / (2dt)) *(U_(i + 1) - U_i)
	*/

	// Calculate primitive variables from solution variables
	calc_Prim(U, Prim, gamma, nx);

	// Loop through internal faces. This loop is indexed so that i gives the cell to the right of the face
	for(int i = 1; i < nx; i++)
	{	
		// Initialize vectors to store flux and solution from right and left cell
		double F_left[3] = {0};
		double F_right[3] = {0};
		double U_left[3] = {0};
		double U_right[3] = {0};

		// Populate flux from left and right cell
		F_left[0] = Prim[0][i - 1] * Prim[1][i - 1];
		F_left[1] = Prim[0][i - 1] * Prim[1][i - 1] * Prim[1][i - 1] + Prim[2][i - 1];
		F_left[2] = Prim[1][i - 1] * (Prim[2][i - 1] + Prim[3][i - 1]);
		
		F_right[0] = Prim[0][i] * Prim[1][i];
		F_right[1] = Prim[0][i] * Prim[1][i] * Prim[1][i] + Prim[2][i];
		F_right[2] = Prim[1][i] * (Prim[2][i] + Prim[3][i]);

		// Populate solution from left and right cell
		U_left[0] = U[0][i - 1];
		U_left[1] = U[1][i - 1];
		U_left[2] = U[2][i - 1];

		U_right[0] = U[0][i];
		U_right[1] = U[1][i];
		U_right[2] = U[2][i];

		// Populate flux vector
		for(int j = 0; j < 3; j++)
		{
			F[j][i] = 0.5 * (F_left[j] + F_right[j]) - (dx / (2 * dt)) * (U_right[j] - U_left[j]);
		}
	}

	// Apply boundary conditions
	calc_bound(F, Prim, nx);
}


// Timestepping schemes
void forward_euler(std::array<std::vector<double>, 3> &U, const std::array<std::vector<double>, 3> &F, double dt, double dx, int nx)
{
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < 3; j++)
	{
		U[j][i] = U[j][i] - (dt / dx) * (F[j][i + 1] - F[j][i]);
	}
}


void RK2(std::array<std::vector<double>, 3> &U, const std::array<std::vector<double>, 3> &F, double dt, double dx, int nx, double gamma, int scheme)
{
	// Store intermediate step
	std::array<std::vector<double>, 3> U_int = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::array<std::vector<double>, 3> F_int = {std::vector<double>(nx + 1), std::vector<double>(nx + 1), std::vector<double>(nx + 1)};
	std::array<std::vector<double>, 4> Prim_int = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < 3; j++)
		{
			U_int[j][i] = U[j][i] - (dt / dx) * (F[j][i + 1] - F[j][i]);
		}

	// Calculate interediate flux based on above values
	switch(scheme)
		{
			case 0:
			{
				Roe_flux(U_int, Prim_int, F_int, dx, gamma, nx);
				break;
			}
			case 1:
			{	
				Steger_warming(U_int, Prim_int, F_int, dx, gamma, nx);
				break;
			}
		}

	// Calculate U at next timestep
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < 3; j++)
		{
			U[j][i] = U[j][i] - (F[j][i + 1] - F[j][i] + F_int[j][i + 1] - F_int[j][i]) * (dt / (2.0 * dx));
		}
}

void midpoint(std::array<std::vector<double>, 3> &U, const std::array<std::vector<double>, 3> &F, double dt, double dx, int nx, double gamma, int scheme)
{
	// Store intermediate step
	std::array<std::vector<double>, 3> U_int = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::array<std::vector<double>, 3> F_int = {std::vector<double>(nx + 2), std::vector<double>(nx + 2), std::vector<double>(nx + 2)};
	std::array<std::vector<double>, 4> Prim_int = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < 3; j++)
		{

			U_int[j][i] = U[j][i] - (dt / 2.0) * (F[j][i + 1] - F[j][i]) / dx;
		}

	// Calculate interediate flux based on above values
	switch(scheme)
		{
			case 0:
			{
				Roe_flux(U_int, Prim_int, F_int, dx, gamma, nx);
				break;
			}
			case 1:
			{	
				Steger_warming(U_int, Prim_int, F_int, dx, gamma, nx);
				break;
			}
		}

	// Calculate U at next timestep
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < 3; j++)
		{
			U[j][i] = U_int[j][i] - dt * (F_int[j][i + 1] - F_int[j][i]) / dx;
		}
}

// This function will apply the flux limiter to the flux 
void apply_limiter(std::array<std::vector<double>, 3> &F, std::array<std::vector<double>, 3> &FL, std::array<std::vector<double>, 3> &FH, std::array<std::vector<double>, 3> phi, int nx)
{
	for(int i = 0; i < nx + 1; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			F[j][i] = FL[j][i] + phi[j][i] * (FH[j][i] - FL[j][i]);
		}
	}
}

// This function will compute the flux limiter vector based on selected flux limiting scheme
void calc_limiter(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> &phi, int nx, int limitscheme)
{
	switch(limitscheme)
	{
		case 0:
			{
				superbee(U, phi, nx);
				break;
			}

		case 1:
			{
				vanleer(U, phi, nx);
				break;
			}
	}
}


void superbee(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> &phi, int nx)
{
	// This function will limit the flux based on the superbee method. 
	// This will give a sharper resolution than van_leer (I think), but may be less
	// stable.

	double r[3]; // Store the r value for the cell

	// Loop through cells. 
	for(int i = 0; i < nx; i++)
	{
		// for(int j = 0; j < 3; j++)
		// {
		// 	// Calculate R of the cell
		// 	if (i == 0)
		// 	{
		// 		r[j] = ((U[j][i] - (2.0 * U[j][i] - U[j][i + 1])) / ((U[j][i + 1]) - U[j][i]));
		// 	}

		// 	else if (i == (nx - 1))
		// 	{
		// 		r[j] = (U[j][i] - U[j][i - 1]) / ((2.0 * U[j][i] - U[j][i - 1]) - U[j][i]);
		// 		printf("cell: %d r[0]: %f r[1]: %f r[2]: %f\n", i, r[0], r[1], r[2]);
		// 	}

		// 	else
		// 	{
		// 		r[j] = (U[j][i] - U[j][i - 1]) / (U[j][i + 1] - U[j][i]);
		// 	}

		// 	// Populate Phi-vector
		// 	phi[j][i] = std::max({0.0, std::min({2.0 * r[j], 1.0}), std::min({r[j], 2.0})});
		// 	if(i == nx - 1)
		// 	{
		// 		printf("phi[0]: %f phi[1]: %f phi[2]: %f\n", phi[0][i], phi[1][i], phi[2][i]);
		// 	}

		// }
		for(int j = 0; j < 3; j++)
		{
			phi[j][i] = 0.0;
		}
	}
}

void vanleer(const std::array<std::vector<double>, 3> &U, std::array<std::vector<double>, 3> &phi, int nx)
{
	printf("vanleer\n");
}

double minmod(double a, double b)
{
	if ((fabs(a) < fabs(b)) && ((a * b) > 0.0))
	{
		return a;
	}

	else if ((fabs(b) < fabs(a)) && ((a * b) > 0.0))
	{
		return b;
	}

	else 
	{
		return 0.0;
	}
}

double maxmod(double a, double b)
{
	if ((fabs(a) > fabs(b))  && ((a * b) > 0.0))
	{
		return a;
	}

	else if ((fabs(a) > fabs(b))  && ((a * b) > 0.0))
	{
		return b;
	}

	else
	{
		return 0.0;
	}
}
