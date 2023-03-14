// Included header
#include "field_functions.h"
constexpr double pi() { return std::atan(1)*4; } // Define pi

/********************************************************************
// Unstructured Mesh Programming Assignment Field Functions File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

void initial_condition(const Mesh &mesh, std::vector<double> &temp_cent, std::array<std::vector<double>, 2> &vel)
{
	// This function will initialize both the velocities u, and v at the midpoints. 
	// It will also initialize the temperature T at the cell centroid

	// Loop through edge midpoints to calculate velocity at these points
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		vel[0][i] = mesh.Edge_centroid[1][i] * pi(); 
		vel[1][i] = -mesh.Edge_centroid[0][i] * pi();
	}

	// for(int i = 0; i < mesh.iNEdge; i++)
	// {
	// 	printf("Edge %i has velocity U: %f V: %f\n", i, vel[0][i], vel[1][i]);
	// }

	// Loop through cell centroids to calculate velocity at these points
	for(int i = 0; i < mesh.iNCell; i++)
	{
		temp_cent[i] = exp(-5 * (mesh.Cell_centroid[0][i] * mesh.Cell_centroid[0][i] + (mesh.Cell_centroid[1][i] - 1) * (mesh.Cell_centroid[1][i] - 1)));
	}
	// temp_cent = {100.0, 102.0, 101.0, 97.0}; // Used for testing flux based on analytical 
	// for(int i = 0; i < mesh.iNCell; i++)
	// {
	// 	printf("Cell %i has initial temp %14.12e\n", i, temp_cent[i]);
	// }
}

void save_VTK(std::string fileName, const Mesh mesh, const std::vector<double> temp_cent)
{
	const char * fileChar = fileName.c_str();
	FILE *vtkFile;
	vtkFile = fopen(fileChar, "w");

	// Print data to VTK file
	fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
	fprintf(vtkFile, "Mech 588 VTK File\n");
	fprintf(vtkFile, "ASCII\n");
	fprintf(vtkFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(vtkFile, "POINTS %i DOUBLE\n", mesh.iNVert);

	// Print vertex coordinates
	for(int i = 0; i < mesh.iNVert; i++)
	{
		fprintf(vtkFile, "%f %f %f \n", mesh.Vert[0][i], mesh.Vert[1][i], 0.0);
	}

	// Print cell connectivity
	fprintf(vtkFile, "CELLS %i %i\n", mesh.iNCell, 4*mesh.iNCell);
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%i %i %i %i\n", 3, mesh.Cell_vert[0][i], mesh.Cell_vert[1][i], mesh.Cell_vert[2][i]);
	}

	// Print cell type for each cell (triangle = 5)
	fprintf(vtkFile, "CELL_TYPES %i\n", mesh.iNCell);
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%i\n", 5);
	}

	// Create lookup table for temperature
	fprintf(vtkFile, "CELL_DATA %i\n", mesh.iNCell);
	fprintf(vtkFile, "SCALARS temperature DOUBLE 1\n");
	fprintf(vtkFile, "LOOKUP_TABLE temperature\n");
	for(int i = 0; i < mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%14.12e\n", temp_cent[i]);
	}

	fprintf(vtkFile, "VECTORS Velocity DOUBLE\n");
	for(int i = 0; i< mesh.iNCell; i++)
	{
		fprintf(vtkFile, "%14.12e %14.12e %14.12e\n", mesh.Cell_centroid[1][i]*pi(), -mesh.Cell_centroid[0][i]*pi(), 0.0);
	}

	fclose(vtkFile);
}

void calc_grad(const Mesh mesh, std::vector<double> temp_cent, std::array<std::vector<double>, 2> &Cell_Grad)
{
	/*
	This function will calculate the gradient in each cell. We will use the method derived in the analytical section of this project. 
	The overall order of the function is as follows:
	-> Loop through each cell
	-> For each cell, count how many edges of the cell are boundary edges. 
	-> Once we know how many boundary edges of the cell are boundary edges, we can determine which cells to use for the flux calculation. 
	For each cell that has a boundary edge, we will use another edge from a neighboring cell. For example if cell i has neighboring cells a b and c,
	we will use the a b and c to calculate the gradient in i. However, if cell i is on the boundary, and has neighbors a and b, where b is neighbors with 
	d, then we will use a b and d to calculate the gradient in i. Finally, if cell i is on the boundary, and only has neighbor a (2 edges are boundary edges),
	where a has neighbors e and f, we will use a e and f to calculate the gradient in i.
	*/
	// Loop through cells 
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
		// printf("Cell %i will use %i %i and %i in gradient calcs\n", cell, neighbor[0], neighbor[1], neighbor[2]);
		// Now that we have found the cells to use in the flux calculations, we can determine the gradient 
		// using the math outline in analytical. 
		double x0, x1, x2, x3; // Store x coordinates of centroid for each cell
		double y0, y1, y2, y3; // Store y coordinates of centroid for each cell
		double T0, T1, T2, T3; // Store cell averaged temperature for each cell
		double DX2, DY2, DXDY, DXT, DYT; // Store amounts used for solving 2x2 system
		
		x0 = mesh.Cell_centroid[0][cell];
		x1 = mesh.Cell_centroid[0][neighbor[0]];
		x2 = mesh.Cell_centroid[0][neighbor[1]];
		x3 = mesh.Cell_centroid[0][neighbor[2]];

		y0 = mesh.Cell_centroid[1][cell];
		y1 = mesh.Cell_centroid[1][neighbor[0]];
		y2 = mesh.Cell_centroid[1][neighbor[1]];
		y3 = mesh.Cell_centroid[1][neighbor[2]];

		T0 = temp_cent[cell];
		T1 = temp_cent[neighbor[0]];
		T2 = temp_cent[neighbor[1]];
		T3 = temp_cent[neighbor[2]];

		DX2 = pow((x1 - x0), 2) + pow((x2 - x0), 2) + pow((x3 - x0), 2);
		DY2 = pow((y1 - y0), 2) + pow((y2 - y0), 2) + pow((y3 - y0), 2);
		DXDY = ((x1 - x0) * (y1 - y0)) + ((x2 - x0) * (y2 - y0)) + ((x3 - x0) * (y3 - y0));
		//printf("DX2 %f\nDY2 %f\nDXDY %f\n", DX2, DY2, DXDY);
		DXT = ((x1 - x0) * (T1 - T0)) + ((x2 - x0) * (T2 - T0)) + ((x3 - x0) * (T3 - T0));
		DYT = ((y1 - y0) * (T1 - T0)) + ((y2 - y0) * (T2 - T0)) + ((y3 - y0) * (T3 - T0));
		//printf("DXT %f\nDYT %f\n", DXT, DYT);

		Cell_Grad[0][cell] = (1 / (DX2 * DY2 - pow(DXDY, 2))) * (DY2 * DXT - DXDY * DYT);
		Cell_Grad[1][cell] = (1 / (DX2 * DY2 - pow(DXDY, 2))) * (-DXDY * DXT + DX2 * DYT);
	}
}

void calc_upwind(Mesh &mesh, std::array<std::vector<double>, 2> vel)
{
	/* 
	This function will loop through each edge, and take the dot product of the edge normal and the
	velocity at the edge midpoint. Based on the sign of this dot proudct, we can determine the upwind
	and downwind cell. We can then asssign the indices to our vector based on the following convention:
	Our normal point into the right cell if you walk from edge origin to vertex. Therefore, if the dot product 
	is positive, the upwind cell is the left cell, and if it is negative the upwind cell is the right cell. 
	NOTE: IF IM WRONG WE SHOULD SEE ROTATION THAT IS CCW INSTEAD OF CW, AND WE CAN CHANGE INDICES.
	*/ 
	double dot_product; // Store calculate value for the dot product at each edge
	// Loop through edges 
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
		// printf("Edge %i has normal velocity %f with upwind cell %i and downwind cell %i\n", i, dot_product, mesh.Edge_upwind[0][i], mesh.Edge_upwind[1][i]);
	}
}

/*
Now that above works, we need to calculate the flux at each edge. Currently, below calculates
the solution at the edge midpoint based on the upwind cell properly. We just need to find the total 
flux across the edge by multiplying by edge length, and then loop through the cells to add the total flux integral.
After that we just have to timestep and we DONE BABY!!!!
*/

void calc_flux(Mesh &mesh, std::vector<double> temp_cent, std::array<std::vector<double>, 2> Cell_Grad, std::vector<double> &Edge_flux)
{
	/*
	This function will begin by calculating the solution at each edge
	based on the upwind cell.
	*/
	int upwind_cell; // Store upwind cell index for the edge
	double dxf, dyf; // Store difference in x and y between cell centroid, and edge midpoint coordinates
	double Ti; // Store solution at edge midpoint
	// Loop through edges to calculate solution at edge midpoint
	for(int i = 0; i < mesh.iNEdge; i++)
	{
		upwind_cell = mesh.Edge_upwind[0][i];
		dxf = mesh.Edge_centroid[0][i] - mesh.Cell_centroid[0][upwind_cell];
		dyf = mesh.Edge_centroid[1][i] - mesh.Cell_centroid[1][upwind_cell];
		Ti = temp_cent[upwind_cell] + Cell_Grad[0][upwind_cell] * dxf + Cell_Grad[1][upwind_cell] * dyf;
		Edge_flux[i] = mesh.dot_product[i] * Ti * mesh.Edge_length[i];
		// printf("Edge %i has solution %f\n", i, Ti);
	}
}
