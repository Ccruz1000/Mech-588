#include "Fcurv.h"
#include "Flux.h"

/********************************************************************
// Curvillinear Mesh Programming Assignment Flux Integral File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

void calc_flux(Vector &I, Solution &s)
{
	for(int i = 1; i < s.nx - 1; i++)
		for(int j = 1; j < s.ny - 1; j++)
		{
			I(i, j) = s.T(i + 1, j + 1) * 0.25 *  (-s.Si(i, j) - s.Sj(i, j)) +
					  s.T(i + 1, j) * (s.Di(i, j) + 0.25 * (-s.Si(i, j) + s.Sj(i, j))) +
					  s.T(i + 1, j - 1) * 0.25 * (s.Si(i, j) + s.Sj(i, j)) +
					  s.T(i, j + 1) * (s.Dj(i, j) + 0.25 * (-s.Si(i, j) + s.Si(i - 1, j))) +
					  s.T(i, j) * (-s.Di(i, j) - s.Dj(i, j) - s.Di(i - 1, j) - s.Dj(i, j - 1)) +
					  s.T(i, j - 1) * (s.Dj(i, j - 1) + 0.25 * (s.Si(i, j) - s.Si(i - 1, j))) + 
					  s.T(i - 1, j + 1) * 0.25 * (s.Si(i - 1, j) + s.Sj(i, j)) +
					  s.T(i - 1, j) * (s.Di(i - 1, j) + 0.25 * (s.Sj(i, j) - s.Sj(i, j - 1))) +
					  s.T(i - 1, j - 1) * 0.25 * (-s.Si(i - 1, j) - s.Sj(i, j - 1));

		}
}

int main()
{
	Solution test = readmesh("Flux_testmesh-17x17.vel");
	test.T = test.u;
	test.calc_mesh_metric();
	Vector I(test.nx, test.ny);
	calc_flux(I, test);
	for(int i = 0; i < test.nx; i++)
		for(int j = 0; j < test.ny; j++)
		{
			printf("I(i, j): %f\n", I(i, j));
		}
	//test.printmesh();
	return 0;
}