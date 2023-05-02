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

	printf("\nComputing Functional Correction Term\n");
	// Compute the functional error estimate
	double Lu4thOrder[IMAX + 2][JMAX + 2];
	double deltaJ = 0.0;
	computeLu(u, u_4thOrder, Lu4thOrder);
	calc_functional(Lu4thOrder, v, deltaJ);

	// Q3
	printf("Running Q3\n");
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

	double h[IMAX + 2]; // Store h for contour integral
	double e[IMAX + 2]; // store store e for contour integral

	double dv_dn[IMAX + 2]; // Store dv_dn for functional calculation
	double du_dn[IMAX + 2]; // Store du_dn for functional calculation
	
	double contour_int1 = 0.0; // Store value of countour integral (un, h) to compute functional
	double contour_int2 = 0.0; // Store value of contour integral (vn, e) to compute functional 
	// For this we will loop through the top surface to calculate h and e
	// We can also compute dv_dn and du_dn
	for(int i = 1; i < IMAX; i++)
	{
		double dX = (i - 0.5) * dDX;
		double offset = dDX / sqrt(12);
		h[i] = 0.25 * (2.0 * compute_h(dX + offset) + 2.0 * compute_h(dX - offset));
		e[i] = 0.25 * (2.0 * compute_e(dX + offset) + 2.0 * compute_e(dX - offset));

		// Use 2nd order one sided finite difference to compute derivatives
		dv_dn[i] = 2.0 * dDX * (3.0 * v_q3[i][JMAX] - 4.0 * v_q3[i][JMAX - 1] + v_q3[i][JMAX - 2]); // Bottom boundary
		printf("X = %f, dv_dn = %f\n", dX, dv_dn[i]);
		// Above is positive because we have (vn, e) where e is at the top boundary so we get derivative 
		// of v w.r.t y
		du_dn[i] = 2.0 * dDX * (3.0 * u_q3[i][1] - 4.0 * u_q3[i][2] + u_q3[i][3]); // Top boundary  

		contour_int1 += -du_dn[i] * h[i] * dDX * dDY; // (un, h)
		contour_int2 += dv_dn[i] * e[i] * dDX * dDY; // (vn, e)
	}

	double functional3 = 0.0; // J=(u,g)+(un,h) 
	double functional4 = 0.0; // J=(v,f)+(vn,e)

	printf("\n\n\n1. %f\n2. %f\n\n\n", contour_int1, contour_int2);
	calc_functional(u_q3, g_q3, functional3);
	calc_functional(v_q3, f_q3, functional4);

	functional3 += contour_int1; 
	functional4 += contour_int2;

	// Print L2Norms
	printf("Using %d points\nJ=(u,g) functional: %1.16e\nJ=(v,f) functional:%1.16e\n", IMAX, functional1, functional2);
	printf("J=(u,g) error: %1.16e\nJ=(v,f) error: %1.16e\n", 1.0 - functional1, 1.0 - functional2);
	printf("Corrected functional = J - deltaJ: %1.16e, Corrected functional Error - %1.16e\n", functional1 - deltaJ, 1 - (functional1 - deltaJ));
	printf("Q3 Begins Here\n");
	printf("J=(u,g)+(un,h) functional: %1.16e\nJ=(v,f)+(vn,e) functional: %1.16e\n", functional3, functional4);

 	// Let user know runtime
	printf("\nCode succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}