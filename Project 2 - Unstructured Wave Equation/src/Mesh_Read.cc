// Included libraries found in header file
#include "Mesh_Read.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/
// Constructor to initialize object
Mesh::Mesh(int ncell, int nedge, int nbdry, int nvert)
           : iNCell(ncell), iNEdge(nedge), iNBdry(nbdry), iNVert(nvert)
{
	// Initialize cell edge and vert data to -1. They will be updated later
	for(int i = 0; i < iNCell; i++)
	{
		Cell_edge[0][i] = -1;
		Cell_edge[1][i] = -1;
		Cell_edge[2][i] = -1;

		Cell_neighbor[0][i] = -1;
		Cell_neighbor[1][i] = -1;
		Cell_neighbor[2][i] = -1;
	}
}

// Member functions
void Mesh::print_edges()
{
	for(int i = 0; i < iNEdge; i++)
	{
		printf("Edge %i: %i %i %i %i\n", i, Edge[0][i], Edge[1][i], Edge[2][i], Edge[3][i]);
	}

	for(int i = 0; i < iNBdry; i++)
	{
		printf("Bdry %i: %i %i %i\n", i, Bdry[0][i], Bdry[1][i], Bdry[2][i]);
	}
}

void Mesh::print_vert_coord()
{
	for(int i = 0; i < iNVert; i++)
	{
		printf("Vertex: %i, X: %f, Y: %f\n", i, Vert[0][i], Vert[1][i]);
	}	
}

void Mesh::calc_cell_edge()
{
	/*
	Below, we have cell_pos which is initialized to have the length of the number of cells
	and all values of 0. This will be used to select the index for assigning cells their edges.
	We will loop through the list of edges. Once an edge is assigned to the cell, we will increment 
	the value for that cell in Cell_pos, so that the next edge that contains this cell will be placed 
	into the next position in the Cell_edge array. This should be relatively efficient, having only 
	to loop through the list of edges once to assign all the edges.
	*/
	std::vector<int> Cell_pos = std::vector<int>(iNCell); 
	// Initialize Cell_pos to be all 0's.
	int cntr = 0; 
	for(int i = 0; i < iNCell; i++)
	{
		Cell_pos[i] = 0;
	}

	// Loop through each edge and assign indices to cell
	for(int edge = 0; edge < (iNEdge); edge++)
	{
		// Check left edges
		Cell_edge[Cell_pos[Edge[0][edge]]][Edge[0][edge]] = edge;
		Cell_pos[Edge[0][edge]] += 1; // Increment Cell_pos for the cell that was just populated
		// Check right edges, ignoring if its a boundary edge
		if(Edge[1][edge] != -1)
		{
			Cell_edge[Cell_pos[Edge[1][edge]]][Edge[1][edge]] = edge;
			Cell_pos[Edge[1][edge]] += 1; 
		}
	}

	// for(int i = 0; i < iNCell; i++)
	// {
	// 	printf("Cell %i Edges: %i %i %i\n", i, Cell_edge[0][i], Cell_edge[1][i], Cell_edge[2][i]);
	// }
}

void Mesh::calc_cell_neighbor()
{
	/*
	This loop will calculate the neighbors for each cell. To do this there are a few steps. 
	1. Count the number of neighbors. Here we can look at each edge the cell has, and check if it is a 
	boundary edge. If it is we move to the next edge, if it is not we will incremement the number of neighbors
	by 1. This way, when we are searching for neighbors, we can break from the loop, and move onto the next cell
	in order to improve code efficiency
	2. Check each subsequent cell to see if one of its edges matches the edge of the current cell. If it does then
	we can add the cell to our list of neighbors, both ways.
	*/
	int num_neighbor = 0; // Used to calculate the number of neighbors each cell has
	int neighbors_found = 0; // Determine how many neighbors have been found

	std::vector<int> Cell_pos = std::vector<int>(iNCell); // Used for indexing to add cell neighbors
	// Initialize Cell_pos to be all 0's.
	int cntr = 0; 
	for(int i = 0; i < iNCell; i++)
	{
		Cell_pos[i] = 0;
	}

	// Loop through each cell
	for(int cell_index = 0; cell_index < iNCell; cell_index++)
	{
		// Count number of neighbors that cell has
		for(int edge = 0; edge < 3; edge++)
		{
			if(Edge[1][Cell_edge[edge][cell_index]] != -1)
			{
				num_neighbor += 1;
			}
		}

		// Determine neighbors for each cell by looking at all of the subsequent cells
		for(int sub_cells = cell_index + 1; sub_cells < iNCell; sub_cells++)
		{
			/* Compare edge 0 of cell a (the one we're finding neighbors for) to edge 0 1 and 2 of cell
			 b (the one at index sub_cells). Repeat this for edge 1 and 2 of cell b. If a neighbor is found
			add it to the list of neighbors for both cell a and b, and incrememnt neighbors_found.
			If neighbors_found = num_neighbor, break from loop, and move to next cell for cell a. 
			*/
			for(int edge_a = 0; edge_a < 3; edge_a++)
			{
				//printf("Edge a for cell %i = %i\n", cell_index, edge_a);
				for(int edge_b = 0; edge_b < 3; edge_b++)
				{
					//printf("Edge b for cell %i = %i\n", sub_cells, edge_b);
					if(Cell_edge[edge_a][cell_index] == Cell_edge[edge_b][sub_cells])
					{
						Cell_neighbor[Cell_pos[cell_index]][cell_index] = sub_cells; // Add cell b to cell a neighbors
						Cell_neighbor[Cell_pos[sub_cells]][sub_cells] = cell_index; // Add cell a to cell b neighbors
						Cell_pos[sub_cells] += 1; // Increase increment of where to add cell b neighbor
						Cell_pos[cell_index] += 1; // Increase increment of where to add cell a neighbor
						neighbors_found += 1; // Increase neighbors found
					}
				}
				if(neighbors_found == num_neighbor)
				{
					num_neighbor = 0;
					neighbors_found = 0;
					break;
				}
			} 

		}
		// Reset number of neighbors
		num_neighbor = 0;
		neighbors_found = 0;
	}
	// for(int i = 0; i < iNCell; i++)
	// {
	// 	printf("%i %i %i - Neighbors for Cell %i\n", Cell_neighbor[0][i], Cell_neighbor[1][i], Cell_neighbor[2][i], i);
	// }
}

void Mesh::calc_edge_length()
{
	/*
	In this function, we will loop over each edge, and calculate its length. We will
	achieve this by using the pythagorean theorem, and finding dx and dy between the 
	points.
	*/

	int index_origin, index_dest; // Find indices of vertices at each end of edge
	double x_origin, y_origin; // Find coordinates at origin
	double x_dest, y_dest; // Find coordinates at destination
	double length; // Calculate length of edge
	double dx, dy; // Find distance in x and y of edge

	for(int i = 0; i < iNEdge; i++)
	{
		// Determine indices of edge origin and destination
		index_origin = Edge[2][i];
		index_dest = Edge[3][i];

		// Determine x and y position at origin
		x_origin = Vert[0][index_origin];
		y_origin = Vert[1][index_origin];

		// Determine x and y position at destination
		x_dest = Vert[0][index_dest];
		y_dest = Vert[1][index_dest];

		// Calculate length of edge
		dx = x_dest - x_origin;
		dy = y_dest - y_origin;
		length = sqrt((dx * dx) + (dy * dy));

		Edge_length[i] = length;
	}

	for(int i = 0; i < iNEdge; i++)
	{
		printf("Edge %i has length %f\n", i, Edge_length[i]);
	}
}

// Add comments, and refactor this
void Mesh::calc_cell_vert()
{
	std::vector<int> Cell_pos = std::vector<int>(iNCell);
	int cntr = 0;
	int max_ed,min_ed;
	for(int i = 0; i < iNCell; i++)
	{
	Cell_pos[i] = 0;
	}
	for( int i= 0; i < iNCell ; i++)
	{
	max_ed = std::max(Cell_edge[0][i], Cell_edge[1][i]);
	max_ed = std::max(max_ed,Cell_edge[2][i]);

	min_ed = std::min(Cell_edge[0][i], Cell_edge[1][i]);
	min_ed= std::min(min_ed, Cell_edge[2][i]);

	Cell_vert[0][i]= Edge[2][max_ed];
	Cell_vert[1][i]= Edge[3][max_ed];
	Cell_vert[2][i]= Edge[2][min_ed];
	if(Edge[2][min_ed]== Edge[3][max_ed] || Edge[2][min_ed]== Edge[2][max_ed] )
	        {
	            Cell_vert[2][i]= Edge[3][min_ed];
	        }

	printf("Cell %i verts: %i %i %i\n", i, Cell_vert[0][i], Cell_vert[1][i], Cell_vert[2][i]);

	}
}



// void Mesh::calc_cell_centroid()
// {
// 	int i = 0;
// }

Mesh read_mesh(std::string meshname)
{
	std::ifstream mesh(meshname); // Read mesh
	
	// Check if mesh is opened, return error otherwise
	if(mesh.is_open())
	{
		// Read number of cells, edges, boundary edges, and vertices
		int iNCell, iNEdge, iNBdry, iNVert;
		mesh >> iNCell >> iNEdge >> iNBdry >> iNVert;
		printf("Number of Cells - %i\nNumber of Edges - %i\nNumber of Boundary Edges - %i\nNumber of Vertices - %i\n", iNCell, iNEdge, iNBdry, iNVert);
		
		// Store X and Y coordinates for each vertex (X is Vert[0], Y is Vert[1])
		std::array<std::vector<double>, 2> Vert = {std::vector<double>(iNVert), std::vector<double>(iNVert)};
		for(int i = 0; i < iNVert; i++)
		{
			mesh >> Vert[0][i] >> Vert[1][i];
		}

		/*
		Store boundary edges and internal edges seperately
		For Bdry - Bdry[0]: Connected cell, Bdry[1]: Edge origin vertex, Bdry[2]: Edge destination vertex
		For Edge - Edge[0]: Left cell Edge[1]: Right cell, Edge[2]: Edge origin vertex, Edge[3]: Edge destination vertex
		*/
		std::array<std::vector<int>, 3> Bdry = {std::vector<int>(iNBdry), std::vector<int>(iNBdry), std::vector<int>(iNBdry)};
		std::array<std::vector<int>, 4> Edge = {std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge), std::vector<int>(iNEdge)};
		std::array<std::vector<int>, 4> Tot_Edge = {std::vector<int>(iNEdge + iNBdry), std::vector<int>(iNEdge + iNBdry), std::vector<int>(iNEdge + iNBdry), std::vector<int>(iNEdge + iNBdry)};
		// int total_edge = iNBdry + iNEdge; // Total number of edges found in the mesh
		int left, right, origin, destination; // Temporarily store values to allow for selecting which array they belong to
		// Counters to use for indices
		int Bdrycntr = 0;
		int Edgecntr = 0;
		for(int i = 0; i < iNEdge; i++)
		{
			mesh >> left >> right >> origin >> destination;
			Edge[0][i] = left;
			Edge[1][i] = right;
			Edge[2][i] = origin;
			Edge[3][i] = destination;
			if(right == -1)
			{
				Bdry[0][Bdrycntr] = left;
				Bdry[1][Bdrycntr] = origin;
				Bdry[2][Bdrycntr] = destination;
				Bdrycntr += 1;
			}
		}

		Mesh temp(iNCell, iNEdge, iNBdry, iNVert);
		// Populate array of vectors for vertex coordinate data
		for(int i = 0; i < iNVert; i++)
		{
			temp.Vert[0][i] = Vert[0][i];
			temp.Vert[1][i] = Vert[1][i];
		}

		// Populate array of vectors for boundary index data
		for(int i = 0; i < iNBdry; i++)
		{
			temp.Bdry[0][i] = Bdry[0][i];
			temp.Bdry[1][i] = Bdry[1][i];
			temp.Bdry[2][i] = Bdry[2][i];
		}

		// Populate array of vectors for edge index data
		for(int i = 0; i < iNEdge; i++)
		{
			temp.Edge[0][i] = Edge[0][i];
			temp.Edge[1][i] = Edge[1][i];
			temp.Edge[2][i] = Edge[2][i];
			temp.Edge[3][i] = Edge[3][i];
		}

		return temp;
	}

	else
	{
		printf("Mesh file not found. Please check that proper filename was used\n");
		exit(1);
	}
}

int main()
{
// Calculate runtime
time_t start, end;
time(&start);
// Read mesh 
std::string verycoarse = "Face-Cell/mech511-square-verycoarse.mesh";
std::string coarse = "Face-Cell/mech511-square-coarse.mesh";
std::string medium = "Face-Cell/mech511-square-medium.mesh";
std::string fine = "Face-Cell/mech511-square-fine.mesh";
std::string veryfine = "Face-Cell/mech511-square-veryfine.mesh";
std::string analytical = "Face-Cell/analytical.mesh";


Mesh mesh = read_mesh(analytical);
// mesh.print_edges();
mesh.calc_cell_edge();
mesh.calc_cell_neighbor();
mesh.calc_edge_length();
mesh.calc_cell_vert();
time(&end);
double time_taken = double(end - start);
printf("Code succesfully run in %.3f seconds\n", time_taken);
return 0;
}