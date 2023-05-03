// Included header files
#include "poisson.h"
#include "poisson4th.h"
#include "Q_Ans.h"

int main()
{
	// Calculate runtime
	clock_t tstart = clock();
	Q1();
	Q2();
	Q3(); 

 	// Let user know runtime
	printf("\nCode succesfully ran in %6.3e seconds\n", (double)(clock() - tstart)/CLOCKS_PER_SEC);
	return 0;
}