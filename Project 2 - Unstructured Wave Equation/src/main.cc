// Included headers
#include "Mesh_Read.h"
#include "field_functions.h"

/********************************************************************
// Unstructured Mesh Programming Assignment Main File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

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

std::vector<double> temp_cent = std::vector<double>(mesh.iNCell); // Store temperatures at cell centroids
std::array<std::vector<double>, 2> vel = {std::vector<double>(mesh.iNEdge), std::vector<double>(mesh.iNEdge)}; // Store velocity u and v at edge midpoints

initial_condition(mesh, temp_cent, vel);
std::string file = "Results/Test.vtk";
save_VTK(file, mesh, temp_cent);

std::array<std::vector<double>, 2> Cell_Grad = {std::vector<double>(mesh.iNCell), std::vector<double>(mesh.iNCell)};
calc_grad(mesh, temp_cent, Cell_Grad);

time(&end);
double time_taken = double(end - start); 
printf("Code succesfully run in %.3f seconds\n", time_taken);
return 0;
}