// Included headers
#include "Roe.h"

/********************************************************************
// Compressible Flow Programming Assignment Test File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
*********************************************************************/
int tests_failed = 0; 

void print_matrix(double A[3][3])
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			printf("%14.12e ", A[i][j]); 
		}
		std::cout << std::endl;
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

	printf("Matric-vector multiplication passed\n");
}

int main()
{
	// Calculate runtime
	clock_t tstart = clock();
	printf("\n\n*******************************************\n           TESTING\n*******************************************\n");

	// Run tests
	test_matrix_mult();
	test_mat_vec_mult();

	// Let user know if testing failed or passed
	if(tests_failed == 0)
	{
		printf("*******************************************\n           TESTING PASSED\n*******************************************\n\n");
	}

	else
	{
		if(tests_failed == 1)
		{
			printf("*******************************************\n        %i TEST FAILED\n*******************************************\n\n", tests_failed);
		}

		else
		{
			printf("*******************************************\n        %i TESTS FAILED\n*******************************************\n\n", tests_failed);
		}
	}

	// Let user know how long testing took
	printf("Tests succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}