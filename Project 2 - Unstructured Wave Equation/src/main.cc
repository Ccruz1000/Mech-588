#include "Mesh_Read.h"

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
// mesh.print_vert_coord();
// mesh.calc_all_param();
time(&end);
double time_taken = double(end - start); 
printf("Code succesfully run in %.3f seconds\n", time_taken);
return 0;
}