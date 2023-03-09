// Included headers
#include "Mesh_Read.h"
#include "field_functions.h"

/********************************************************************
// Unstructured Mesh Programming Assignment Test File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

int tests_failed = 0;

// Test reading of vertex coordinates
void test_vertex(Mesh mesh)
{
	double coord_tol = 1e-8; // Tolerance to account for rounding/truncation errors
	double x_test, y_test; // Store x and y values for comparison with those read from mesh
	double x_error, y_error; 
	// Open test file to obtain data
	std::string testfile = "Test_Files/vert.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 6; i++)
		{
			test >> x_test >> y_test; // Read coordinates from test file
			// Calculate x and y errors
			x_error = (fabs(mesh.Vert[0][i] - x_test));
			y_error = (fabs(mesh.Vert[1][i] - y_test));

			// Check if error has been exceeded
			if (x_error > coord_tol)
			{
				tests_failed += 1;
				printf("Error reading mesh X coordinate at vertex %i\n", i);
				return;
			}

			if (y_error > coord_tol)
			{
				tests_failed += 1;
				printf("Error reading mesh Y coordinate at vertex %i\n",  i);
				return;
			}
		}
	}
	else
	{
		printf("Vertex test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

void test_cell_area(Mesh mesh)
{
	double area_tol = 1e-6; // Tolerance to account for rounding errors
	double test_area; // Error read from file
	double area_error; // Calculated error

	std::string testfile = "Test_Files/area.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 4; i++)
		{
			test >> test_area;

			area_error = fabs(mesh.Cell_area[i] - test_area);

			if(mesh.Cell_area[i] < 0)
			{
				tests_failed += 1;
				printf("Error, negative cell area detected at cell %i\n", i);
				return;
			}

			if(test_area < 0)
			{
				tests_failed += 1;
				printf("Error, negative test area detected at cell %i\n", i);
				return;
			}

			if(area_error > area_tol)
			{
				tests_failed += 1;
				printf("Area error exceed tolerance for cell %i\n", i);
				return;
			}
		}
	}
	else
	{
		printf("Area test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

void test_edge_length(Mesh mesh)
{
	double tol = 1e-6; // Tolerance to account for rounding errors
	double test_value; // Error read from file
	double error; // Calculated error

	std::string testfile = "Test_Files/length.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 8; i++)
		{
			test >> test_value; // Read value from file

			error = fabs(mesh.Edge_length[i] - test_value); // Compute error

			// Ensure length isnt negative
			if(mesh.Edge_length[i] < 0)
			{
				tests_failed += 1;
				printf("Error, negative edge length detected at edge %i\n", i);
				return;
			}

			// Ensure test length isnt negative
			if(test_value < 0)
			{
				tests_failed += 1;
				printf("Error, negative test edge length detected at edge %i\n", i);
				return;
			}

			// Ensure tolerance is met
			if(error > tol)
			{
				tests_failed += 1;
				printf("Error exceed tolerance for cell %i\n", i);
				return;
			}
		}
	}

	else
	{
		printf("Edge length test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

void test_edge_norm(Mesh mesh)
{
	double tol = 1e-8; // Tolerance to account for rounding/truncation errors
	double x_test, y_test; // Store x and y values for comparison with those read from mesh
	double x_error, y_error; 
	// Open test file to obtain data
	std::string testfile = "Test_Files/norm.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 8; i++)
		{
			test >> x_test >> y_test; // Read coordinates from test file
			// Calculate x and y errors (abs of both because orientation could be different)
			x_error = (fabs(fabs(mesh.Edge_norm[0][i]) - fabs(x_test)));
			y_error = (fabs(fabs(mesh.Edge_norm[1][i]) - fabs(y_test)));

			// Check if error has been exceeded
			if (x_error > tol)
			{
				printf("Error calculating normal vector X component at edge %i\n", i);
				return;
			}

			if (y_error > tol)
			{
				printf("Error calculating normal vector Y component at edge %i\n",  i);
				return;
			}

			if (sqrt(pow(mesh.Edge_norm[0][i], 2) + pow(mesh.Edge_norm[1][i], 2)) - 1 > tol)
			{
				printf("Normal vector length is not 1 at edge %i\n", i);
				return;
			}
		}
	}

	else
	{
		printf("Norm test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

void test_edge_midpoint(Mesh mesh)
{
	double tol = 1e-8; // Tolerance to account for rounding/truncation errors
	double x_test, y_test; // Store x and y values for comparison with those read from mesh
	double x_error, y_error; 
	// Open test file to obtain data
	std::string testfile = "Test_Files/midpoint.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 8; i++)
		{
			test >> x_test >> y_test; // Read coordinates from test file
			// Calculate x and y errors 
			x_error = (fabs(mesh.Edge_centroid[0][i] - x_test));
			y_error = (fabs(mesh.Edge_centroid[1][i] - y_test));

			// Check if error has been exceeded
			if (x_error > tol)
			{
				tests_failed += 1;
				printf("Error calculating normal vector X component at edge %i\n", i);
				return;
			}

			if (y_error > tol)
			{
				tests_failed += 1;
				printf("Error calculating normal vector Y component at edge %i\n",  i);
				return;
			}
		}
	}

	else
	{
		printf("Midpoint test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

void test_cell_centroid(Mesh mesh)
{
	double tol = 1e-8; // Tolerance to account for rounding/truncation errors
	double x_test, y_test; // Store x and y values for comparison with those read from mesh
	double x_error, y_error; 
	// Open test file to obtain data
	std::string testfile = "Test_Files/centroid.txt";
	std::ifstream test(testfile);

	if(test.is_open())
	{
		for(int i = 0; i < 4; i++)
		{
			test >> x_test >> y_test; // Read coordinates from test file
			// Calculate x and y errors (abs of both because orientation could be different)
			x_error = (fabs(fabs(mesh.Cell_centroid[0][i]) - fabs(x_test)));
			y_error = (fabs(fabs(mesh.Cell_centroid[1][i]) - fabs(y_test)));

			// Check if error has been exceeded
			if (x_error > tol)
			{
				tests_failed += 1;
				printf("Error calculating normal vector X component at edge %i\n", i);
				return;
			}

			if (y_error > tol)
			{
				tests_failed += 1;
				printf("Error calculating normal vector Y component at edge %i\n",  i);
				return;
			}

		}
	}
	else
	{
		printf("Cell centroid test file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

int main()
{
	std::string analytical = "Face-Cell/analytical.mesh";
	Mesh mesh = read_mesh(analytical);

	printf("\n\n*******************************************\n           TESTING\n*******************************************\n");
	
	test_vertex(mesh);
	printf("Vertices passed\n");
	test_cell_area(mesh);
	printf("Cell area passed\n");
	test_edge_length(mesh);
	printf("Edge length passed\n");
	test_edge_norm(mesh);
	printf("Norm passed\n");
	test_edge_midpoint(mesh);
	printf("Midpoint passed\n");
	test_cell_centroid(mesh);
	printf("Centroid passed\n");

	// After running tests inform user if tests failed or not
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


return 0;
}
