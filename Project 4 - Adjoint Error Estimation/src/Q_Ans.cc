// Include Libraries
#include "Q_Ans.h"

void Q1()
{
	printf("\n*****************************************************\n");
	printf("*****************************************************\n");
	printf("\nRUNNING Q1\n");

	double u[IMAX + 2][JMAX + 2]; // Solution to primal problem
	double v[IMAX + 2][JMAX + 2]; // Solution to dual problem

	double f[IMAX + 2][JMAX + 2]; // Source to primal problem
	double g[IMAX + 2][JMAX + 2]; // Source to dual problem


	printf("\n************************************************\n               PRIMAL FUNCTIONS\n************************************************\n");
	run_Poisson(vApplyHomoBC, Primal_Source, u, f);
	printf("2nd Order L2 - %f\n", calc_L2Norm(u));
	printf("************************************************\n");
	printf("\n************************************************\n               DUAL FUNCTIONS\n************************************************\n");
	run_Poisson(vApplyHomoBC, Dual_Source, v, g);
	printf("2nd Order L2 - %f\n", calc_L2Norm(v));
	printf("************************************************\n");

	// Calculate functional
	double functional1 = 0.0; // J=(u,g)
	double functional2 = 0.0; // J=(v,f)

	calc_functional(u, g, functional1);
	calc_functional(v, f, functional2);

	std::string u_save = "Results/u.vtk";
	std::string v_save = "Results/v.vtk";

	storeVTKSolution(u, u_save);
	storeVTKSolution(v, v_save);

	printf("Using %d points\nJ=(u,g) functional: %1.16e\nJ=(v,f) functional:%1.16e\n", IMAX, functional1, functional2);
	printf("J=(u,g) error: %1.16e\nJ=(v,f) error: %1.16e\n", 1.0 - functional1, 1.0 - functional2);
	printf("*****************************************************\n");
	printf("*****************************************************\n\n");
}

void Q2()
{
	printf("\n*****************************************************\n");
	printf("*****************************************************\n");
	printf("\nRUNNING Q2\n");

	double u[IMAX + 2][JMAX + 2]; // Solution to primal problem
	double v[IMAX + 2][JMAX + 2]; // Solution to dual problem

	double f[IMAX + 2][JMAX + 2]; // Source to primal problem
	double g[IMAX + 2][JMAX + 2]; // Source to dual problem

	double u_4thOrder[IMAX + 4][JMAX + 4]; // u for the 4th order problem
	double v_4thOrder[IMAX + 4][JMAX + 4]; // v for the 4th order problem

	double f_4thOrder[IMAX + 4][JMAX + 4]; // f for the 4th order problem
	double g_4thOrder[IMAX + 4][JMAX + 4]; // g for the 4th order problem


	printf("\n************************************************\n               PRIMAL FUNCTIONS\n************************************************\n");
	printf("Running 2nd Order Primal\n");
	run_Poisson(vApplyHomoBC, Primal_Source, u, f);
	printf("2nd Order L2 - %f\n", calc_L2Norm(u));
	printf("Running 4th Order Primal\n");
	run_Poisson_4thOrder(vApplyHomoBC_4thOrder, Primal_Source, u_4thOrder, f_4thOrder);
	printf("4th Order L2 - %f\n", calc_L2Norm_4thOrder(u_4thOrder));
	printf("************************************************\n");

	printf("\n************************************************\n               DUAL FUNCTIONS\n************************************************\n");
	printf("Running 2nd Order Dual\n");
	run_Poisson(vApplyHomoBC, Dual_Source, v, g);
	printf("2nd Order L2 - %f\n", calc_L2Norm(v));
	printf("Running 4th Order Dual\n");
	run_Poisson_4thOrder(vApplyHomoBC_4thOrder, Dual_Source, v_4thOrder, g_4thOrder);
	printf("************************************************\n");

	std::string u_4 = "Results/u_4th.vtk";
	std::string v_4 = "Results/v_4th.vtk";
	storeVTKSolution_4thOrder(u_4thOrder, u_4);
	storeVTKSolution_4thOrder(v_4thOrder, v_4);

	printf("\nComputing Functional Correction Term\n");
	// Compute the functional error estimate
	double functional1 = 0.0;
	calc_functional(u, g, functional1);
	double Lu4thOrder[IMAX + 2][JMAX + 2];
	double deltaJ = 0.0;
	computeLu(u, u_4thOrder, Lu4thOrder);
	calc_functional(Lu4thOrder, v, deltaJ);

	printf("Corrected functional = J - deltaJ: %1.16e, Corrected functional Error - %1.16e\n\n", functional1 - deltaJ, 1 - (functional1 - deltaJ));
	printf("*****************************************************\n");
	printf("*****************************************************\n\n");
}

void Q3()
{
	printf("\n*****************************************************\n");
	printf("*****************************************************\n");
	printf("\nRUNNING 3\n");
	double u_q3[IMAX + 2][JMAX + 2];
	double f_q3[IMAX + 2][IMAX + 2];
	run_Poisson(vApplyBC_uQ3, Primal_Source, u_q3, f_q3);
	std::string u_q3_save = "Results/u_Q3.vtk";
	storeVTKSolution(u_q3, u_q3_save);

	double v_q3[IMAX + 2][JMAX + 2];
	double g_q3[IMAX + 2][JMAX + 2];
	run_Poisson(vApplyBC_vQ3, Dual_Source, v_q3, g_q3);
	std::string v_q3_save = "Results/v_Q3.vtk";
	storeVTKSolution(v_q3, v_q3_save);
	
	// Compute the functional from the primal function, J = (u,g)+(un,h)
	double functional3 = 0.0; // J=(u,g)+(un,h) 
	double functional4 = 0.0; // J=(v,f)+(vn,e)

	compute_q3_primal_functional(u_q3, g_q3, functional3);

	printf("J=(u,g)+(un,h) functional: %1.16e\nJ=(v,f)+(vn,e) functional: %1.16e\n", functional3, functional4);
	printf("*****************************************************\n");
	printf("*****************************************************\n\n");
}