// Included header files
#include "poisson.h"
#include "poisson4th.h"

int main()
{
	// Calculate runtime
	clock_t tstart = clock();

	double u[IMAX + 2][JMAX + 2]; // Solution to primal problem
	double v[IMAX + 2][JMAX + 2]; // Solution to dual problem

	double f[IMAX + 2][JMAX + 2]; // Source to primal problem
	double g[IMAX + 2][JMAX + 2]; // Source to dual problem

	double u_4thOrder[IMAX + 4][JMAX + 4]; // u for the 4th order problem
	double v_4thOrder[IMAX + 4][JMAX + 4]; // v for the 4th order problem

	double f_4thOrder[IMAX + 4][JMAX + 4]; // f for the 4th order problem
	double g_4thOrder[IMAX + 4][JMAX + 4]; // g for the 4th order problem

	printf("\n************************************************\n               PRIMAL FUNCTIONS\n************************************************\n");
	printf("Running Primal\n");

	run_Poisson(vApplyHomoBC, Primal_Source, u, f);
	printf("2nd Order L2 - %f\n", calc_L2Norm(u));
	printf("Running 4th Order Primal\n");
	
	run_Poisson_4thOrder(vApplyHomoBC_4thOrder, Primal_Source, u_4thOrder, f_4thOrder);
	printf("4th Order L2 - %f\n", calc_L2Norm_4thOrder(u_4thOrder));
	printf("************************************************\n");

	printf("\n************************************************\n               DUAL FUNCTIONS\n************************************************\n");
	printf("Running Dual\n");

	run_Poisson(vApplyHomoBC, Dual_Source, v, g);
	printf("Running 4th Order Dual\n");
	run_Poisson_4thOrder(vApplyHomoBC_4thOrder, Dual_Source, v_4thOrder, g_4thOrder);
	printf("************************************************\n\n");

	// Calculate functional
	double functional1 = 0.0; // J=(u,g)
	double functional2 = 0.0; // J=(v,f)

	calc_functional(u, g, functional1);
	calc_functional(v, f, functional2);

	std::string u_4 = "Results/u_4th.vtk";
	std::string u_save = "Results/u.vtk";
	std::string v_save = "Results/v.vtk";
	std::string v_4 = "Results/v_4th.vtk";
	storeVTKSolution(u, u_save);
	storeVTKSolution_4thOrder(u_4thOrder, u_4);
	storeVTKSolution(v, v_save);
	storeVTKSolution_4thOrder(v_4thOrder, v_4);


	// Print L2Norms
	printf("Using %d points\nJ=(u,g) functional: %1.16e\nJ=(v,f) functional:%1.16e\n", IMAX, functional1, functional2);
	printf("J=(u,g) error: %1.16e\nJ=(v,f) error: %1.16e\n", 1.0 - functional1, 1.0 - functional2);

 	// Let user know runtime
	printf("\nCode succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}