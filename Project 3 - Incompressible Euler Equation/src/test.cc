// Included headers
#include "functions.h"
#include "Steger.h"
#include "Roe.h"

/********************************************************************
// Compressible Flow Programming Assignment Test File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
*********************************************************************/
int tests_failed = 0; 

void print_vector(double A[3])
{
	for(int i = 0; i < 3; i++)
	{
		printf("%f\n", A[i]);
	}
}

void test_matrix_mult()
{
	double tol = 1e-14; // Tolerance for round off errors from copy and pasting values out of matlab

	// Test 1 - Array populated with integers
	double ans1[3][3] = {{92, 77, 15}, {84, 92, 49}, {14, 23, 112}};
	double A1[3][3] = {{0, 4, 7}, {1, 2, 8}, {3, 0, 1}};
	double B1[3][3] = {{2, 4, 37}, {9, 0, 2}, {8, 11, 1}};
	double C1[3][3] = {0};

	// Test 2 - Array populated with random positive decimals
	double ans2[3][3] = { {0.634769461616549,   0.757334032772758,   1.044774278904069},
						  {0.741916598978334,   1.055395097349698,   0.972094721338084},
						  {0.709618498718592,   0.894767870137602,   1.106059893597152} };

	double A2[3][3] = { {7.278408103560474e-02,   9.521485708286341e-01,   6.307239147279311e-01},
				        {8.428842861262987e-01,   8.122607403152180e-01,   4.182929741803995e-01},
				        {2.199665616099711e-01,   9.837728737665031e-01,   6.799655931937661e-01} };

	double B2[3][3] = { {0.307120683562613,   0.634773165368482,   0.138743584256598},
				        {0.425320486721343,   0.269551875844950,   0.934499317967610},
				        {0.328904083825703,   0.720567600661611,   0.229726130875622} };
	double C2[3][3] = {0};

	// Test 3 - Array populated with random decimals
	double ans3[3][3] = { {5.526416317742069e-01,  -5.937847883095778e-01,   6.525443493228130e-01},
						  {-3.157238466175765e-01,  -4.749014532523135e-01,   5.768061852106061e-02},
						  {5.389485672746748e-01,   1.129945910594600e-01,   2.625426036829367e-01} };

	double A3[3][3] = { {7.278408103560474e-02,   9.521485708286341e-01,   6.307239147279311e-01},
				        {-8.428842861262987e-01,   8.122607403152180e-01,   -4.182929741803995e-01},
				        {2.199665616099711e-01,   -9.837728737665031e-01,   6.799655931937661e-01} };

	double B3[3][3] = { {-1.437976888804915e-02,   4.730366299813864e-01,  -2.507484980777674e-01},
				        {2.726217641215412e-02,  -3.413505473226895e-01,   2.016985815622616e-01},
				        {8.367061964938479e-01,  -4.807141157234681e-01,   7.590449973033032e-01} };
	double C3[3][3] = {0};

	// Test 4 - Check that it can do multiple multiplications in a row - D = A * B * C
	double ans4[3][3] = { {1.234072316759171e+00, 2.843943072863482e-01, 2.717423998780961e-01},
				   {-1.365259524838104e+00,  -2.745957382351005e-01,  -2.178431637403847e-01},
				   {4.472213099753786e-01,   1.471741762896765e-01,   9.968185436187840e-02} };

	double A4[3][3] = { {3.724937551620764e-01,   7.519666831431802e-01,   7.734272699527877e-01},
				        {-6.489828655470700e-01,  -9.480388309613396e-01,  -1.412100765478300e-01},
				        {9.272194933728584e-02,   3.813228882969860e-01,   2.517556456848080e-01} };

	double B4[3][3] = { {0.449749054716817,  -0.942583463752558,   0.547062318340494},
				        {0.102986260240843,  -0.251420445924351,   0.676374786118009},
				        {0.267563303913484,  -0.245700271981016,   0.318850420028027} };

	double C4[3][3] = { {0.559342518121196,  -0.183462331105557,   0.761166069531486},
				        {-0.660921500415916,   0.498641662348991,   0.315121376277699},
				        {0.520018351366212,   0.762682886958945,   0.164613160491391} };

	double D4[3][3] = {0};
	double intermediate[3][3] = {0};

	// Test 5 - Check that I am doing F=A*B*C*D*E properly (required for calculating A)
	double ans5[3][3] = { {0.244968854310229, 0.111299728712215, 0.102381466641423},
				   		  {1.703828050883825, 0.788071231268292, 0.694344887199920},
				   	      {1.250532177000984, 0.568677515181263, 0.523004865156332} };

	double A5[3][3] = { {0.186247990375911, 0.186735956047455,  -0.309209610694437},
				        {-1.393950088689942,  -0.444610592984647,   3.709441728101316},
				        {0.588956295942359,   0.843078798690492,  -1.685259834637217} };

	double B5[3][3] = { {-0.374508125079890,  -0.143722481332039,   0.176611789300290},
				        {0.441575322441360,   0.649063404349929,   0.411791312186093},
				        {-0.107741278396972,   0.162291863363063,   0.171732939187435} };

	double C5[3][3] = { {-2.223355103336122e-01,   6.682175092774623e-02,  -4.793242237512246e-01},
				        {9.870013518250698e-02,  -1.606853582186608e+00,   5.412056732710324e-02},
				        {-3.073825877160580e-01,   6.349739882515631e-01,  -3.617916667619571e-01} };

	double D5[3][3] = { {-1.971860233147539,  -0.783243535772908,  -0.897352793774935},
				        {-1.341637347237744,  -0.819435367664623,  -0.677935685181097},
				        {-1.598616549967567,  -0.754605609438523,  -0.303529804533033} };

	double E5[3][3] = { {0.985930116573769,   0.391621767886454,  0.448676396887467},
				        {0.670818673618872,   0.409717683832312,   0.338967842590549},
				        {0.799308274983784,   0.377302804719262,   0.151764902266516} };

	double F5[3][3] = {0};
	double intermediate2[3][3] = {0};


	matrix_mult3x3(A1, B1, C1);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if (std::fabs(C1[i][j] - ans1[i][j]) > tol)
			{
				tests_failed += 1;
				printf("Matrix multiplication not properly working for test case 1\n");
				print_matrix(C1);
				return;
			}
		}

	matrix_mult3x3(A2, B2, C2);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if (std::fabs(C2[i][j] - ans2[i][j]) > tol)
			{
				tests_failed += 1;
				printf("Matrix multiplication not properly working for test case 2\n");
				print_matrix(C2);
				return;
			}
		}

	matrix_mult3x3(A3, B3, C3);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if (std::fabs(C3[i][j] - ans3[i][j]) > tol)
			{
				tests_failed += 1;
				printf("Matrix multiplication not properly working for test case 3\n");
				print_matrix(C3);
				return;
			}
		}

	// Test 4
	// Start by obtaining A4 * B4
	matrix_mult3x3(A4, B4, intermediate);
	// Next obtain the above product * C4
	matrix_mult3x3(intermediate, C4, D4);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if (std::fabs(D4[i][j] - ans4[i][j]) > tol)
			{
				tests_failed += 1;
				printf("Matrix multiplication not properly working for test case 4\n");
				print_matrix(D4);
				return;
			}
		}

	// Test 5
	// Start by assigning 0 values to intermediate, then find A5*B5
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			intermediate[i][j] = 0;
		}

	matrix_mult3x3(A5, B5, intermediate);
	// Now find A5*B5*C5 and 0 intermediate again
	matrix_mult3x3(intermediate, C5, intermediate2);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			intermediate[i][j] = 0;
		}
	// Now find A5*B5*C5*D5 and 0 intermediate2;
	matrix_mult3x3(intermediate2, D5, intermediate);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			intermediate2[i][j] = 0;
		}
	// Finally find F5 = A5*B5*C5*D5*E5
	matrix_mult3x3(intermediate, E5, F5);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if (std::fabs(F5[i][j] - ans5[i][j]) > tol)
			{
				printf("Error - %14.12e\n", std::fabs(F5[i][j] - ans5[i][j]));
				tests_failed += 1;
				printf("Matrix multiplication not properly working for test case 5\n");
				print_matrix(F5);
				return;
			}
		}

	printf("Matrix multiplication testing passed\n");
}

void test_mat_vec_mult()
{
	double tol = 1e-15; // Tolerance for roundoff errors
	// Test 1 - Simple matric and vector of integers
	double ans1[3] = {31, 23, 33};
	double A1[3] = {1, 2, 3};
	double B1[3][3] = {{2, 4, 7}, {1, 2, 6}, {9, 12, 0}};
	double C1[3] = {0};

	matrix_vector(B1, A1, C1);

	for(int i = 0; i < 3; i++)
	{
		if (std::fabs(C1[i] - ans1[i]) > tol)
			{
				tests_failed += 1;
				printf("Matrix-vector multiplication not properly working for test case 1\n");
				return;
			}
	}

	// Test 2 - Matrix and vector of positive decimals
	double ans2[3] = {1.374531889323082, 0.657613100144160, 0.775249935796417};
	double A2[3] = {0.918097952428487, 0.499734321941770, 0.145922859523648};
	double B2[3][3] = {{0.972339603800191,   0.736531596412903,   0.779581571100408}, 
					   {0.160971842164582,   0.736558241903146,   0.971347109896470},
					   {0.235097721171814,   0.960551209444794,   0.544032595611690}};
	double C2[3] = {0};

	matrix_vector(B2, A2, C2);

	for(int i = 0; i < 3; i++)
	{
		if (std::fabs(C2[i] - ans2[i]) > tol)
			{
				tests_failed += 1;
				printf("Matrix-vector multiplication not properly working for test case 2\n");
				return;
			}
	}

	// Test 3 - Matrix and vector of positive decimals
	double ans3[3] = {0.279317005401431, -0.120182844437230, -0.200914580479742};
	double A3[3] = {-0.442560055209259, 0.253928902616895, -0.182317723245457};
	double B3[3][3] = {{-1.777890172715650e-01,   1.598287298518278e-01,  -8.778606450055515e-01}, 
					   {1.767302293342071e-01,  -1.486301460967837e-01,   2.318817071942647e-02},
					   {5.898059933482745e-01,   7.879613941416386e-02,  -2.199532474050649e-01}};
	double C3[3] = {0};

	matrix_vector(B3, A3, C3);

	for(int i = 0; i < 3; i++)
	{
		if (std::fabs(C3[i] - ans3[i]) > tol)
			{
				tests_failed += 1;
				printf("Matrix-vector multiplication not properly working for test case 3\n");
				return;
			}
	}

	printf("Matrix-vector multiplication passed\n");
}

// The following will check if a matrix is actually an identity matrix
void test_identity(double matrix[3][3], int test, int cell)
{
	double tol = 1e-15;
	std::string message; // Store the print message based the test number
	switch(test)
	{
		case 0:
			message = "XR * XL calculation test failed at cell ";
			break;
		case 1:
			message = "du_dv * dv_du test failed at cell ";
			break;
		case 2:
			message = "Matrix multiplication to compute A failed at cell ";
			break;
	}

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{
			if(i == j)
			{
				if(std::fabs(matrix[i][j] - 1.0) > tol)
				{
					std::cout << message << cell <<std::endl;
					print_matrix(matrix);
					tests_failed += 1;
					return;
				}
			}

			else
			{
				if(std::fabs(matrix[i][j] - 0.0) > tol)
				{
					std::cout << message << cell <<std::endl;
					print_matrix(matrix);
					tests_failed += 1;
					return;
				}
			}
		}
}

// The following function will take place due to the code refactor. It will perform the above 4 tests using my newly defined
// functions. 
void test_total_jacobian(std::array<std::vector<double>, 3> U, std::array<std::vector<double>, 4> Prim, double gamma, double nx)
{
	double tol = 1e-15;
	// Initialize and calculate F(U), du, and F(U+du)
	double delta = 1e-6; // Multiplier for du
	std::array<std::vector<double>, 3> F = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::array<std::vector<double>, 3> Fdu = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	std::array<std::vector<double>, 3> du = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};
	for(int i = 0; i < nx; i++)
	{
		// Store du
		du[0][i] = delta * U[0][i];
		du[1][i] = delta * U[1][i];
		du[2][i] = delta * U[2][i];

		// Store F(U)
		F[0][i] = U[1][i]; // Density * velocity
		F[1][i] = 0.5 * (3 - gamma) * U[1][i] * U[1][i] / U[0][i] + (gamma - 1) * U[2][i]; // (rho*u^2 + P)
		F[2][i] = gamma * (U[1][i] * U[2][i] / U[0][i]) + 0.5 * (gamma - 1) * (U[1][i] * U[1][i] * U[1][i]) / (U[0][i] * U[0][i]); // Energy
 		
 		// Store F(U + du)
 		Fdu[0][i] = U[1][i] * (1 +delta); // Density * velocity
		Fdu[1][i] = (0.5 * (3 - gamma) * U[1][i] * U[1][i] / U[0][i] + (gamma - 1) * U[2][i]) * (1 + delta); // (rho*u^2 + P)
		Fdu[2][i] = (gamma * (U[1][i] * U[2][i] / U[0][i]) + 0.5 * (gamma - 1) * (U[1][i] * U[1][i] * U[1][i]) / (U[0][i] * U[0][i])) * (1 + delta); // Energy
	}

	// Initialize and calculate speed of sound
	std::vector<double> c = std::vector<double>(nx);
	for(int i = 0; i < nx; i++)
	{
		c[i] = std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
	}

	// Initialize all of the matrices
	double XR[3][3], du_dv[3][3], lambda[3][3], XL[3][3], dv_du[3][3], A[3][3];
	// Initialize identity matrix
	double I[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	// Initialize eigenvalues
	double lam1, lam2, lam3;

	double test_num; // Store test index

	// Assemble the matrices and lambda at each cell
	for(int i = 0; i < nx; i++)
	{
		// Calculate eigenvalues 
		lam1 = Prim[1][i];
		lam2 = Prim[1][i] + c[i];
		lam3 = Prim[1][i] - c[i];

		// Assemble matrices
		assemble_matrices(Prim[0][i], Prim[1][i], c[i], gamma, du_dv, XR, XL, dv_du);
		
		// Assemble Lambda matrix
		assemble_lambda(lam1, lam2, lam3, lambda);

		// Set test index to 0
		test_num = 0;

		// Test if XR * XL gives an identity matrix
		double X_test[3][3] = {0};
		matrix_mult3x3(XR, XL, X_test);
		test_identity(X_test, test_num, i);
		test_num += 1;

		// If test passes for every cell, inform user
		if(i == (nx - 1))
		{
			printf("XR * XL test passed\n");
		}

		// Test if du_dv * dv_du gives an identity matrix
		double du_test[3][3] = {0};
		matrix_mult3x3(du_dv, dv_du, du_test);
		test_identity(du_test, test_num, i);
		test_num += 1;

		// If test passes for every cell, inform user
		if(i == (nx - 1))
		{
			printf("du_dv * dv_du test passed\n");
		}

		// Test if du_dv * XR * I * XL * dv_du gives an identity matrix
		double inter1[3][3] = {0};
		double inter2[3][3] = {0};
		double result[3][3] = {0};

		matrix_mult3x3(du_dv, XR, inter1);
		matrix_mult3x3(inter1, I, inter2);
		// Zero inter1
		for(int k = 0; k < 3; k++)
			for(int j = 0; j < 3; j++)
			{
				inter1[k][j] = 0.0;
			}
		matrix_mult3x3(inter2, XL, inter1);
		matrix_mult3x3(inter1, dv_du, result);
		test_identity(result, test_num, i);

		// If test passes for every cell, inform user
		if(i == (nx - 1))
		{
			printf("A calculation test passed\n");
		}

		// Final test will check that F(U) - F(U + du) = Adu
		// Calculate A
		double inter3[3][3] = {0};
		double inter4[3][3] = {0};
		double A[3][3] = {0};

		matrix_mult3x3(du_dv, XR, inter3);
		matrix_mult3x3(inter3, lambda, inter4);
		// Zero inter3
		for(int k = 0; k < 3; k++)
			for(int j = 0; j < 3; j++)
			{
				inter3[k][j] = 0.0;
			}
		matrix_mult3x3(inter4, XL, inter3);
		matrix_mult3x3(inter3, dv_du, A);

		double store[3] = {du[0][i], du[1][i], du[2][i]}; // Store du as an array of doubles
		double check[3] = {0}; // Store Adu
		matrix_vector(A, store, check);

		// Check that they are equivalent everywhere
		for(int l = 0; l < 3; l++)
		{
			if(std::fabs(Fdu[l][i] - F[l][i] - check[l]) > tol)
			{
				printf("Jacobian matrix test failed for cell %d at F[%d]\n", i, l);
				tests_failed += 1;
				return;
			}
		}

		// If test passes for every cell, inform user
		if(i == (nx - 1))
		{
			printf("Jacobian calculation test passed\n");
		}	
	}
}

// This function will make sure the maximum stable timestep is being found
void test_stable_dt()
{
	double tol = 1e-15;
	double gamma = 1.0;
	double CFL = 1.0;
	double dx = 1.0; 
	// This function will work by first defining a primitive vector that has data for 4 cells
	// We will set rho to 1, and gamma to 1. We also set the (absolute) velocity to be the highest in the third cell
	// We need to find the maximum absolute value, and the index of this cell to calculate the local speed of sound
	// In this case, the highest velocity is 1.0 in cell 3. Here we have P = 4, rho = 1. c = sqrt(gamma*P/rho) = sqrt(4) = 2
	// Our cfl is 1, dx is 1, so we expect to find a value for the CFL of dt = CFL * dx / max(abs(u) + c) which gives
	// dt = 1 * 1 / (1. + 2) = 1/3
	std::vector<double> rho = {1, 1, 1, 1};
	std::vector<double> u = {0.9, 0.6, -1.0, 0.8};
	std::vector<double> P = {1, 1, 4, 1};
	std::vector<double> E = {1, 1, 1, 1};

	std::array<std::vector<double>, 4> Prim = {rho, u, P, E};

	std::vector<double> tot_u = std::vector<double>(4);
	for(int i = 0; i < 4; i++)
	{
		tot_u[i] = std::fabs(Prim[1][i]) + std::sqrt(gamma * Prim[2][i] / Prim[0][i]);
	}
	double maxvel = *std::max_element(tot_u.begin(), tot_u.end());
	double dt = CFL * dx / maxvel;
	
	if(dt != 1.0 / 3.0)
	{
		printf("Maximum stable time step test failed\n");
		tests_failed += 1;
		return;
	}

	printf("Maximum stable time step test passed\n");
}

// This function will test that the flux integral is 0 given uniform initial conditions
void test_uniform_init()
{
	double tol = 1e-15; // Tolerance
	// Initialize variables 
	double CFL, xbegin, xend, dx, dt, gamma, t_final;
	int nx, num_save, timescheme;
	double init_cond[6] = {0};

	// Assign values to provide uniform initial conditions
	CFL = 0.6;
	nx = 201;
	xbegin = 0;
	xend = 1;
	gamma = 1.4;
	init_cond[0] = 1.0; // left density
	init_cond[1] = 0; // left velocity
	init_cond[2] = 1.0; // left pressure
	init_cond[3] = 1.0; // right density
	init_cond[4] = 0; // right velocity
	init_cond[5] = 1.0; // right pressure
	t_final = 0.15;
	num_save = 20;
	timescheme = 0; // euler default

	// Define dx
	dx = (xend - xbegin) / (nx - 1);

	// Create array to store primitives. Prim[0][i] is rho, Prim[1][i] is velocity, Prim[2][i] is pressure, and Prim[3][i] is energy
	std::array<std::vector<double>, 4> Prim = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};

	// Create initial conditions based on uniform input
	initial_conditions(Prim, nx, init_cond, gamma);

	// Create and populate array to store conserved variables. Calculate speed of sound in each cell at the same time
	std::array<std::vector<double>, 3> U = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};  // Solution vector
	std::array<std::vector<double>, 3> F_roe = {std::vector<double>(nx + 2), std::vector<double>(nx + 2), std::vector<double>(nx + 2)}; // Flux vector
	std::array<std::vector<double>, 3> F_steger = {std::vector<double>(nx + 2), std::vector<double>(nx + 2), std::vector<double>(nx + 2)}; // Flux vector

	// Populate solution and using primitives
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = Prim[0][i];  // rho
		U[1][i] = Prim[0][i] * Prim[1][i];  // rho * u
		U[2][i] = Prim[3][i]; // Et
	}
	
	double F_cell[3]; // Store flux integral at cell

	// Test Roe Scheme
	Roe_flux(U, Prim, F_roe, dx, gamma, nx);

	// Calculate flux integral
	for(int i = 1; i < nx - 1; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			F_cell[j] = -(F_roe[j][i] - F_roe[j][i - 1]);
			if(fabs(F_cell[j] > tol))
			{
				printf("Roe flux integral failed at cell %d.\nFI[0]: %14.12e FI[1]: %14.12e FI[2]: %14.12e\n", i, F_cell[0], F_cell[1], F_cell[2]);
				tests_failed += 1;
				return;
			}
		}
	}

	// Let user know if Roe scheme flux integral passes
	printf("Roe flux integral test passed.\n");

	// Test Steger-Warming Scheme
	Steger_warming(U, Prim, F_steger, dx, gamma, nx);

	// Calculate flux integral
	for(int i = 1; i < nx - 1; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			F_cell[j] = -(F_steger[j][i] - F_steger[j][i - 1]);
			if(fabs(F_cell[j] > tol))
			{
				printf("Steger-Warming flux integral failed at cell %d.\nFI[0]: %14.12e FI[1]: %14.12e FI[2]: %14.12e\n", i, F_cell[0], F_cell[1], F_cell[2]);
				tests_failed += 1;
				return;
			}
		}
	}

	// Let user know if Steger-Warming scheme flux integral passes
	printf("Steger-Warming flux integral test passed.\n");
}


/*
Enter the values in Wolfram Alpha
wolfram alpha 
dudv
{{1, 0, 0},{u, r, 0},{u^2/2, r*u, 1/(g-1)}}
dvdu
{{1, 0, 0},{-u/r, 1/r, 0},{(g-1)*u^2/2,-(g-1)*u, g-1}}
XR
{{1, r/c, -r/c}, {0, 1, 1}, {0, r*c, -r*c}}
XL
{{1, 0, -1/c^2},  {0, 1/2, 1/(2*r*c)}, {0, 1/2, -1/(2*r*c)}}
lamda1
{{l, 0, 0}, {0, 0, 0}, {0, 0, 0}}
lamda2
{{0, 0, 0}, {0, l, 0}, {0, 0, 0}}
lamda3
{{0, 0, 0}, {0, 0, 0}, {0, 0, l}}
du_dv*XR
{{1, r/c, -r/c}, {u, ru/c+r, r-ru/c}, {u^2/2, ru^2/(2c)+ru, -rc/(g-1)-ru^2/(2c)+ru}}

XL*dv_du
{{1 - (g-1)u^2/(2c^2), (g-1)u/c, -(g-1)/c^2}, {(g-1)u^2/(4rc), 1/(2r)-(g-1)u/(2rc), (g-1)/(2rc)}, {-(g-1)u^2/(4rc)-u/(2r), (g-1)u/(2rc)+1/(2r), -(g-1)/(2rc)}}

*/

int main(int argc, char *argv[])
{
	// Calculate runtime
	clock_t tstart = clock();
	
	// Declare run parameters
	double CFL, xbegin, xend, dx, gamma, t_final;
	int nx, num_save, timescheme;

	// We will store the initial condition values in an array. The first 3 entries will be for the left side
	// and will represent rho, u, and P in that order. The next 3 entries will be for the right side, and will 
	// represent rho, u, and P in that order
	double init_cond[6] = {0};

	// Read input file, or provide default values if none given
	assign_input(argc, argv, CFL, nx, xbegin, xend, gamma, init_cond, t_final, num_save, timescheme);

	// printf("\nCFL: %f xbegin: %f xend: %f dt: %f gamma: %f nx: %d\n", CFL, xbegin, xend, dt, gamma, nx);
	// printf("\nInitial conditions:\nLeft: rho - %f u - %f P - %f\nRight: rho - %f u - %f P - %f\n", init_cond[0], init_cond[1], init_cond[2], init_cond[3], init_cond[4], init_cond[5]);

	// Calculate the spatial step
	dx = (xend - xbegin) / (nx - 1);

	// Create array to store primitives. Prim[0][i] is rho, Prim[1][i] is velocity, Prim[2][i] is pressure, and Prim[3][i] is energy
	std::array<std::vector<double>, 4> Prim = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};

	// Assign initial condition values
	initial_conditions(Prim, nx, init_cond, gamma);

	// Create and populate array to store conserved variables. Calculate speed of sound in each cell at the same time
	std::array<std::vector<double>, 3> U = {std::vector<double>(nx), std::vector<double>(nx), std::vector<double>(nx)};  // Solution vector
	// Populate solution and sound speed vector using primitives
	for(int i = 0; i < nx; i++)
	{
		U[0][i] = Prim[0][i];  // rho
		U[1][i] = Prim[0][i] * Prim[1][i];  // rho * u
		U[2][i] = Prim[3][i]; // Et
	}

	printf("\n\n*******************************************\n           TESTING\n*******************************************\n\n");

	// Run tests
	test_matrix_mult();
	test_mat_vec_mult();
	test_total_jacobian(U, Prim, gamma, nx);
	test_stable_dt();
	test_uniform_init();

	// Let user know if testing failed or passed
	if(tests_failed == 0)
	{
		printf("\n*******************************************\n           TESTING PASSED\n*******************************************\n\n");
	}

	else
	{
		if(tests_failed == 1)
		{
			printf("\n*******************************************\n        %i TEST FAILED\n*******************************************\n\n", tests_failed);
		}

		else
		{
			printf("\n*******************************************\n        %i TESTS FAILED\n*******************************************\n\n", tests_failed);
		}
	}

	// Let user know how long testing took
	printf("Tests succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}