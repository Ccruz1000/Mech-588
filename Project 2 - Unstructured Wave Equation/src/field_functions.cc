// Included header
#include "field_functions.h"
constexpr double pi() { return std::atan(1)*4; } // Define pi

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

	for(int i = 0; i < mesh.iNEdge; i++)
	{
		printf("Edge %i has velocity U: %f V: %f\n", i, vel[0][i], vel[1][i]);
	}

	// Loop through cell centroids to calculate velocity at these points
	for(int i = 0; i < mesh.iNCell; i++)
	{
		temp_cent[i] = exp(-5 * (mesh.Cell_centroid[0][i] * mesh.Cell_centroid[0][i] + (mesh.Cell_centroid[1][i] - 1) * (mesh.Cell_centroid[1][i] - 1)));
	}

	for(int i = 0; i < mesh.iNCell; i++)
	{
		printf("Cell %i has initial temp %14.12e\n", i, temp_cent[i]);
	}
}
