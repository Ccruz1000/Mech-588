// Included header files
#include "poisson.h"

int main()
{
	double u[IMAX + 2][JMAX + 2]; // Solution to primal problem
	double v[IMAX + 2][JMAX + 2]; // Solution to dual problem

	double f[IMAX + 2][JMAX + 2]; // Source to primal problem
	double g[IMAX + 2][JMAX + 2]; // Source to dual problem

	printf("Running Primal\n");
	run_Poisson(vApplyBC, Primal_Source, u, f);
	printf("Runing Dual\n");
	run_Poisson(vApplyBC, Dual_Source, v, g);
 
	return 0;
}