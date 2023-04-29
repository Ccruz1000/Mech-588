// Included header files
#include "poisson.h"

int main()
{
	// Calculate runtime
	clock_t tstart = clock();

	double u[IMAX + 2][JMAX + 2]; // Solution to primal problem
	double v[IMAX + 2][JMAX + 2]; // Solution to dual problem

	double f[IMAX + 2][JMAX + 2]; // Source to primal problem
	double g[IMAX + 2][JMAX + 2]; // Source to dual problem

	printf("Running Primal\n");
	run_Poisson(vApplyHomoBC, Primal_Source, u, f);
	printf("Running Dual\n");
	run_Poisson(vApplyHomoBC, Dual_Source, v, g);

	// Calculate functional
	double functional1 = 0.0; // J=(u,g)
	double functional2 = 0.0; // J=(v,f)

	calc_functional(u, g, functional1);
	calc_functional(v, f, functional2);

	// Compute L2Norm of functionals


	// Print L2Norms
	printf("J=(u,g) functional: %f\nJ=(v,f) functional:%f\n", functional1, functional2);

 	// Let user know runtime
	printf("\nCode succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}