#include <iostream>
#include <fstream>
#include <string>

#define SIZE 10000


std::string coursemesh = "mesh-8x32.vel";
std::string medmesh = "mesh-16x64.vel";
std::string finemesh = "mesh-32x128.vel";



int main()
{
	// std::ifstream f(coursemesh);

	// if(f.is_open())
	// 	std::cout << f.rdbuf();
	std::string line; // Read file lines
	int i =0; // Variable to keep count of lines
	std::string arr[SIZE];  //Array to store lines 

	std::ifstream myfile (coursemesh); // Opens a new ifstream object (used for reading files) called myfile using 
								  // coursemesh filename
	if(myfile.is_open()) // check that myfile is open
	{
		while(!myfile.eof()) // ensure that myfile isnt at the last line
		{
			getline(myfile, line);
			arr[i] = line;
			i++;
		}
		myfile.close();
	}
	else
	{
		std::cout << "Couldnt open the file\n";
	}

	for(int j=0; j<i; j++)
	{
		std::cout << arr[j] << std::endl;
	}

	std::cout << "Testing to find each individual line\n";
	std::cout << arr[0][2] << std::endl;
	return 0;
}

