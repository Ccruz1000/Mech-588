#include <iostream>
#include <fstream>
#include <string>


std::string coursemesh = "mesh-8x32.vel";
std::string medmesh = "mesh-16x64.vel";
std::string finemesh = "mesh-32x128.vel";

int main()
{
	std::ifstream f(coursemesh);

	if(f.is_open())
		std::cout << f.rdbuf();
}

