// Included libraries found in header file
#include "Roe.h"

/**************************************************************************************
// Compressible Flow Programming Assignment Roe's Flux Difference Splitting Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
***************************************************************************************/

// Goated video - https://www.youtube.com/watch?v=F_PsFHvt8IU&ab_channel=P.S.Volpiani

void read_input(std::string input_file, double &CFL, int &nx, double &xbegin, double &xend, double &dt, double &gamma)
{
	/*
		This function will read the input file, which will be used to select the 
		problem parameters. If they are not defined, or the input file is not defined
		properly, then default parameters will be assigned.
	*/

	// Assign default parameters if input file is not input properly
	if(input_file == "Null")
	{
		CFL = 0.95;
		nx = 200;
		xbegin = 0;
		xend = 1;
		dt = 2.5e-2;
		gamma = 1.4;
	}

	// If input file is properly provided, try to read values from it
	else
	{
		std::ifstream input(input_file); // Read input file

		// Check if mesh is opened, return error otherwise
		if(input.is_open())
		{
			std::string var; // Store the variable that is being populated
			std::string val; // Store the value as a string
			// Loop through each line of the input file
			while(true)
			{
				input >> var >> val;

				// Check which variable is being populated and assign value appropriately
				// Definitely better way to do this, but if it ain't broke don't fix it
				if (var == "CFL:")
				{
					CFL = std::stod(val);
				}

				else if (var == "nx:")
				{
					nx = std::stoi(val);					
				}

				else if (var == "xbegin:")
				{
					xbegin = std::stod(val);					
				}

				else if (var == "xend:")
				{
					xend = std::stod(val);					
				}

				else if(var == "dt:")
				{
					dt = std::stod(val);					
				}

				else if (var == "gamma:")
				{
					gamma = std::stod(val);					
				}

				// Leave loop once end of line is reached
				if(input.eof())
				{
					break;
				}
			}
			// Close input file once we are done with it
			input.close();
		}
		// Give warning that input file was not read properly
		else
		{
			printf("Input file not read successfully, please try again.\n");
			exit(1);
		}
	}
}

