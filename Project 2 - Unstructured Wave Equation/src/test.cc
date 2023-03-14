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
	printf("Vertices passed\n");
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
	printf("Cell area passed\n");
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
	printf("Edge length passed\n");
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
	printf("Norm passed\n");
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
	printf("Midpoint passed\n");
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
	printf("Centroid passed\n");
}

void test_grad_calc(Mesh mesh)
{
	double tol = 1e-8; // Acceptable error in calculating gradients 
	std::vector<double> temp_cent = std::vector<double>(mesh.iNCell); // Store temperatures at cell centroids
	std::array<std::vector<double>, 2> Exact_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)}; // Store gradient
	// Populate initial value for cell temperatures
	for(int i = 0; i < mesh.iNCell; i++)
	{
		temp_cent[i] = 5 * mesh.Cell_centroid[0][i] + 6 * mesh.Cell_centroid[1][i]; // 5*x + 6*y linear cell mesh
		Exact_Grad[0][i] = 5.0;
		Exact_Grad[1][i] = 6.0;  // Expected gradient of 5 w.r.t x, and 6 w.r.t y
	}
 
	std::array<std::vector<double>, 2> Cell_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)}; // Store gradient
	calc_grad(mesh, temp_cent, Cell_Grad); // Calculate gradient
	std::array<std::vector<double>, 2> Grad_error = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)}; // Store gradient
	
	for(int i = 0; i < mesh.iNCell; i++)
	{
		Grad_error[0][i] = Exact_Grad[0][i] - Cell_Grad[0][i];
		Grad_error[1][i] = Exact_Grad[1][i] - Cell_Grad[1][i];
		
		// Warn if error exceeds tolerance
		if(Grad_error[0][i] > tol)
		{
			tests_failed += 1;
			printf("Error in calculating cell gradient x-component at cell %i\n", i);
			return;
		}

		if(Grad_error[1][i] > tol)
		{
			tests_failed += 1;
			printf("Error in calculating cell gradient y-component at cell %i\n", i);
			return;
		}
	}
	printf("Gradient calculation passed \n");
}

void test_grad_indices(Mesh mesh)
{
	for(int cell = 0; cell < mesh.iNCell; cell++)
	{
		// Check how many cell edges sit on the boundary.
		int bound_edge = 0;
		for(int edge = 0; edge < 3; edge++)
		{
			if (mesh.Edge[0][mesh.Cell_edge[edge][cell]] == -1 || mesh.Edge[1][mesh.Cell_edge[edge][cell]] == -1)
			{
				bound_edge += 1;
			} 	
		}

		int neighbor[3] = {-1, -1, -1}; // Store the neighbors in this 
		int num_neighbor = 0; // Calculate the neighbors found, used for indexing

		// If there are 0 boundary edges, we will use the centroids of each neighbor
		if(bound_edge == 0)
		{
			for(int i = 0; i < 3; i++)
			{
				neighbor[num_neighbor] = mesh.Cell_neighbor[i][cell];
				num_neighbor += 1;
			}
			// printf("Cell %i has neighbors %i %i %i\n", cell, neighbor[0], neighbor[1], neighbor[2]);
		}

		// If there is 1 boundary edge, we will use the centroids of neighbors a and b, as well as a neighbor of a
		if(bound_edge == 1)
		{
			// a is the index of neighbor cell, a. b is the index of neighbor cell b
			int a;
			for(int i = 0; i < 3; i++)
			{
				if(mesh.Cell_neighbor[i][cell] != -1)
				{ 
					a = mesh.Cell_neighbor[i][cell];
					neighbor[num_neighbor] = a;
					num_neighbor += 1; // Increment index
					// Loop through list of neighbors for cell a
					for(int j = 0; j < 3; j++)
					{
						// Check if neighbor of cell a is a boundary, or the current cell i
						if(mesh.Cell_neighbor[j][a] != -1 && mesh.Cell_neighbor[j][a] != cell)
						{
							neighbor[num_neighbor] = mesh.Cell_neighbor[j][a];
							// If for some reason the loop continues after the list is full, stop
							if(num_neighbor == 2);
							{
								break;
							}
							num_neighbor += 1;
						}
					}
				}
			}
		}

		// If there are 2 boundary edges we will use the centroid of neighbor a, as well as two neighbors of a
		if(bound_edge == 2)
		{
			// a is the index of cell a
			int a;
			for(int i = 0; i < 3; i++)
			{
				if(mesh.Cell_neighbor[i][cell] != -1)
				{
					a = mesh.Cell_neighbor[i][cell];
					neighbor[num_neighbor] = a;
					num_neighbor += 1; // Increment index
					// Loop through list of neighbors for cell a
					for(int j = 0; j < 3; j++)
					{
						// Check if neighbor of cell a is a boundary, or the current cell i
						if(mesh.Cell_neighbor[j][a] != -1 && mesh.Cell_neighbor[j][a] != cell)
						{
							neighbor[num_neighbor] = mesh.Cell_neighbor[j][a];
							num_neighbor += 1;
						}
					}
				}
			}
		}

		// Check if any neighbor is still -1, or if any neighbors repeat
		for(int i = 0; i < 3; i++)
		{
			if(neighbor[i] == -1)
			{
				tests_failed += 1;
				printf("Error, cell %i has neighboring cells that do not exist (VALUE OF -1 FOR BOUNDARY)\n", cell);
				return;
			}

			for(int j = i + 1; j < 3; j++)
			{
				if(neighbor[i] == neighbor[j])
				{
					tests_failed += 1;
					printf("Error, cell %i has repeating indices for gradient calculation\n", cell);
					return;
				}
			}
		}
	}
	printf("Gradient indices passed\n");
}

void test_upwind_calc()
{
	// Read unit square mesh
	std::string meshfile = "Face-Cell/Unit_square.mesh";
	Mesh mesh = read_mesh(meshfile);
	std::string fileName = "Results/Test.vtk";
	// Write values of 0 to allow saving vtk file
	std::vector<double> temp_cent = std::vector<double>(mesh.iNCell);
	for(int i = 0; i < mesh.iNCell; i++)
	{
		temp_cent[i] = 0.0;
	}

	std::array<std::vector<double>, 2> vel = {std::vector<double>(mesh.iNEdge), std::vector<double>(mesh.iNEdge)};
	// Set value of velocity to [1; 0] ( 1 in x direction)
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		vel[0][i] = 1.0;
		vel[1][i] = 0.0;
	}
	// Calculate dot products, and upwind cells based on a horizontal velocity with magnitude 1
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		mesh.dot_product[i] = vel[0][i] * mesh.Edge_norm[0][i] + vel[1][i] * mesh.Edge_norm[1][i]; // compute dot product for cell edge
		// Add left cell as upwind if dot product is positive
		if(mesh.dot_product[i] > 0)
		{
			mesh.Edge_upwind[0][i] = mesh.Edge[0][i];
			mesh.Edge_upwind[1][i] = mesh.Edge[1][i];
		}
		else
		{
			mesh.Edge_upwind[0][i] = mesh.Edge[1][i];
			mesh.Edge_upwind[1][i] = mesh.Edge[0][i];
		}
	}
	// Check dot products for normal velocity
	if(mesh.dot_product[0] != -1.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for horizontal velocity at edge 0\n");
		return;
	}

	if(mesh.dot_product[1] != 0.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for horizontal velocity at edge 1\n");
		return;
	}

	if(mesh.dot_product[2] != 1.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for horizontal velocity at edge 2\n");
		return;
	}

	if(mesh.dot_product[3] != 0.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for horizontal velocity at edge 3\n");
		return;
	}

	if(mesh.Edge_upwind[0][4] != 0 && mesh.Edge_upwind[1][4] != 1)
	{
		tests_failed += 1;
		printf("Error calculating upwind and downwind cell for horizontal velocity at edge 4\n");
		return;
	}

	printf("Upwind cells succesfully tested for horizontal velocity test case\n");

	// Set value of velocity to [0; 1] ( 1 in y direction)
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		vel[0][i] = 0.0;
		vel[1][i] = 1.0;
	}

	// Calculate dot products, and upwind cells based on a vertical velocity with magnitude 1
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		mesh.dot_product[i] = vel[0][i] * mesh.Edge_norm[0][i] + vel[1][i] * mesh.Edge_norm[1][i]; // compute dot product for cell edge
		// Add left cell as upwind if dot product is positive
		if(mesh.dot_product[i] > 0)
		{
			mesh.Edge_upwind[0][i] = mesh.Edge[0][i];
			mesh.Edge_upwind[1][i] = mesh.Edge[1][i];
		}
		else
		{
			mesh.Edge_upwind[0][i] = mesh.Edge[1][i];
			mesh.Edge_upwind[1][i] = mesh.Edge[0][i];
		}
	}
	// Check dot products for normal velocity
	if(mesh.dot_product[0] != 0.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for vertical velocity at edge 0\n");
		return;
	}

	if(mesh.dot_product[1] != 1.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for vertical velocity at edge 1\n");
		return;
	}

	if(mesh.dot_product[2] != 0.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for vertical velocity at edge 2\n");
		return;
	}

	if(mesh.dot_product[3] != -1.0)
	{
		tests_failed += 1;
		printf("Error calculating dot product of edge norm and velocity for vertical velocity at edge 3\n");
		return;
	}

	if(mesh.Edge_upwind[0][4] != 1 && mesh.Edge_upwind[1][4] != 0)
	{
		tests_failed += 1;
		printf("Error calculating upwind and downwind cell for vertical velocity at edge 4\n");
		return;
	}

	printf("Upwind cells succesfully tested for vertical velocity test case\n");	
}

int main()
{
	// Calculate runtime
	time_t start, end;
	time(&start);
	// Read test meshes
	 
	std::string verycoarse = "Face-Cell/mech511-square-verycoarse.mesh";
	std::string coarse = "Face-Cell/mech511-square-coarse.mesh";
	std::string medium = "Face-Cell/mech511-square-medium.mesh";
	std::string fine = "Face-Cell/mech511-square-fine.mesh";
	std::string veryfine = "Face-Cell/mech511-square-veryfine.mesh";
	std::string analytical = "Face-Cell/analytical.mesh";
	Mesh mesh = read_mesh(analytical);
	// Mesh verycoarse_mesh = read_mesh(verycoarse);
	// Mesh coarse_mesh = read_mesh(coarse);
	Mesh medium_mesh = read_mesh(medium);
	// Mesh fine_mesh = read_mesh(fine);
	// Mesh veryfine_mesh = read_mesh(veryfine);

	printf("\n\n*******************************************\n           TESTING\n*******************************************\n");
	
	test_vertex(mesh);
	test_cell_area(mesh);
	test_edge_length(mesh);
	test_edge_norm(mesh);
	test_edge_midpoint(mesh);
	test_cell_centroid(mesh);
	// Run tests on all meshes
	// test_grad_calc(verycoarse_mesh);
	// test_grad_calc(coarse_mesh);
	test_grad_calc(medium_mesh);
	// test_grad_calc(fine_mesh);
	// test_grad_calc(veryfine_mesh);
	// test_grad_indices(verycoarse_mesh);
	// test_grad_indices(coarse_mesh);
	test_grad_indices(medium_mesh);
	// test_grad_indices(fine_mesh);
	// test_grad_indices(veryfine_mesh);
	test_upwind_calc();

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

time(&end);
double time_taken = double(end - start); 
printf("Tests succesfully ran in %.8f seconds\n", time_taken);
return 0;
}
