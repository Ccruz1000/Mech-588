#include "Fcurv.h"
/********************************************************************
// Curvillinear Mesh Programming Assignment Function File
// Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
// Christian Rowsell (40131393)
********************************************************************/

/********************************************************************
// Operations from Mech 587 - Used in operator overloading
********************************************************************/
//Allocate memory for a 1D array
void Operations::alloc1D(double **V, const unsigned long &N, bool isInit, double def) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}

	*V = new double[N];

	if(isInit)
		for(unsigned long i = 0; i < N; i++)
			(*V)[i] = 0.0;
}

// Deallocate memory for a 1D array
void Operations::deAlloc1D(double **V, const unsigned long &N) const
{
	if(*V != nullptr){
		delete[] (*V);
		(*V) = nullptr;
	}
}

// Copy array to new space in memory
void Operations::copyArray(double *S, double *T, const unsigned long &N) const
{
	for(int i = 0; i < N; i++)
		T[i] = S[i];
}

// Addition used for vector + operator overloading
void Operations::addVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] + V2[i];
}

// Subtraction used for vector - operator overloading
void Operations::subtractVec(double *R, const double *V1, const double *V2, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = V1[i] - V2[i];
}

// Scale vector by scalar
void Operations::scaleVec(double *R, const double &a, const double *V1, const unsigned long &N) const
{
	//Create an assertion error here.
	for(unsigned long i = 0; i < N; i++)
		R[i] = a*V1[i];
}

// Take dot product of two vectors
void Operations::dotVec(double &R, const double *V1, const double *V2, const unsigned long &N) const
{
	R = 0.0;
	for(unsigned long i = 0; i < N; i++)
		R = R + (V1[i]*V2[i]);
}

// Find maximum value of the absolute value of the vector
void Operations::findAbsMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(fabs(V2[i]) > R){
			idx = i;
			R = fabs(V2[i]);
		}
}

// Find maximum value in vector
void Operations::findMax(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = -1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] > R){
			idx = i;
			R = V2[i];
		}
}

// Find minimum value in vector
void Operations::findMin(double &R, unsigned long &idx, const double *V2, const unsigned long &N) const
{
	R = 1e20; idx = 0;
	for(unsigned long i = 0; i < N; i++)
		if(V2[i] < R){
			idx = i;
			R = V2[i];
		}
}

/********************************************************************
// Vector class from Mech 587 - Wrapper for 1D array to 2D array
********************************************************************/
// Constructors
Vector::Vector() : Nx(0), Ny(0), N(0), b(nullptr), isInit(false){}

Vector::Vector(unsigned long iNx, unsigned long iNy, bool initflag, double initVal): Nx(iNx), Ny(iNy), N(iNx*iNy), isInit(false), b(nullptr)
{
	O.alloc1D(&b,N,initflag,initVal);
	isInit = true;
}

Vector::Vector(const Vector &V1) : Nx(V1.Nx), Ny(V1.Ny), N(V1.N), isInit(false), b(nullptr)
{
	O.alloc1D(&b,N,false);
	isInit = true;
	O.copyArray(V1.b, b, N);
}

// Destructor
Vector::~Vector()
{
	O.deAlloc1D(&b,N);
	isInit = false;
	Nx = 0, Ny = 0; N = 0;
}

// Set size of vector and allocate memory 
void Vector::setSize(unsigned long iNx, unsigned long iNy)
{
	Nx = iNx, Ny = iNy;
	N = Nx*Ny;
	b = nullptr;
	O.alloc1D(&b,N);
	isInit = true;
}

// Operator overloads
double& Vector::operator()(unsigned long i, unsigned long j)
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else{
		printf("attempted to access invalid memory location. Please check. Exiting\n");
		exit(0);
	}
}

double Vector::operator()(unsigned long i, unsigned long j) const
{
	if(!isInit){
		printf("Matrix has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}
	unsigned long idxNode = i*Ny+j;
	if(isInit && isValid(idxNode))
		return b[idxNode];
	else
		return 0.0;
}

double& Vector::operator()(unsigned long i)
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else{
		printf("attempted to access invalid memory location. Please check. Exiting\n");
		exit(0);
	}
}

double Vector::operator()(unsigned long i) const
{
	if(!isInit){
		printf("Vector has not been dimensioned. Please check. Exiting\n");
		exit(0);
	}

	if(isValid(i))
		return b[i];
	else
		return 0.0;
}

Vector Vector::operator = (const Vector &V1)
{
	Nx = V1.Nx;
	Ny = V1.Ny;

	if(N == V1.N)
		O.copyArray(V1.b, b, V1.N);
	else{
		O.deAlloc1D(&b, N);
		N = V1.N;
		O.alloc1D(&b, N);
		O.copyArray(V1.b, b, N);
	}

	return *this;
}

Vector Vector::operator - ()
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b, -1.0, this->b, this->N);

	return R;
}

Vector Vector::operator + (Vector const &V)
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V)
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a)
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}

Vector Vector::operator + (Vector const &V) const
{
	if(V.N != this->N){
		printf("Invalid vector addition operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny, false);
	O.addVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator - (const Vector &V) const
{
	if(V.N != this->N){
		printf("Invalid vector subtraction operation. Vector dimensions do not match. Please check. Exiting \n");
		exit(0);
	}

	Vector R(this->Nx, this->Ny,false);
	O.subtractVec(R.b,this->b,V.b,this->N);

	return R;
}

Vector Vector::operator * (const double &a) const
{
	Vector R(this->Nx, this->Ny, false);
	O.scaleVec(R.b,a,this->b,this->N);
	return R;
}

Vector operator * (const double &a, const Vector &V1)
{
	Operations O;
	Vector R(V1.Nx, V1.Ny, false);
	O.scaleVec(R.b,a,V1.b,V1.N);
	return R;
}

// Computer L2Norm of vector
const double Vector::L2Norm()
{
	double L2;
	O.dotVec(L2,this->b, this->b, this->N);
	L2 = sqrt(L2/(this->N));
	return L2;
}

// Computer L infinity norm of vector
const double Vector::LinfNorm(unsigned long &ix, unsigned long &iy)
{
	double Linf = -1e20;
	unsigned long idx;
	O.findAbsMax(Linf, idx, this->b, this->N);
	ix = idx/Ny;
	iy = idx%Ny;
	return Linf;
}

// Find maximum value in vector
const double Vector::GetMax()
{
	double max = -1e20;
	unsigned long idx;
	O.findMax(max,idx,this->b, this->N);
	return max;
}

// Find minimum value in vector
const double Vector::GetMin()
{
	double min = 1e20;
	unsigned long idx;
	O.findMin(min, idx, this->b, this->N);
	return min;
}

// Store vector in file
void Vector::storeV(char filename[50])
{
	FILE *f;
	f = fopen(filename,"w");
	for(unsigned long i = 0; i < Nx; i++){
		for(unsigned long j = 0; j < Ny; j++)
			fprintf(f,"%lu %lu %10e\n",i, j, (*this)(i,j));
		fprintf(f,"\n");
	}
}

/********************************************************************
// Solution class constructors and member functions
// Used to store solution, and perform any operations in solving 
********************************************************************/

// Null constructor, initializes everything empty
Solution::Solution() : nx(0), ny(0), N(0), X(), Y(), u(), v(), T()
{}

// Void constructor, initializes vectors with all 0's
Solution::Solution(unsigned long iNx, unsigned long iNy) : nx(iNx), ny(iNy), N(iNx * iNy), X(iNx, iNy), Y(iNx, iNy), u(iNx, iNy), v(iNx, iNy), T(iNx, iNy)
{}

// Constructor, create solution object with pre defined vectors
Solution::Solution(unsigned long iNx, unsigned long iNy, Vector iX, Vector iY, Vector iu, Vector iv, Vector iT)
: nx(iNx) , ny(iNy), N(nx * ny), X(iX), Y(iY), u(iu), v(iv), T(iT)
{}

// Member functions
void Solution::printmesh()
{
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
			{
				printf("X: %f, Y: %f, u: %f, v: %f, T: %f\n", X(i, j), Y(i, j), u(i, j), v(i, j), T(i, j));
				if(j == ny - 1)
				{
					printf("\n");
				}
			}
}


/********************************************************************
// Useful functions - Various functions to be used in the main file
********************************************************************/
Solution readmesh(std::string meshname)
{
	std::ifstream mesh(meshname); // Open mesh file for reading  
	// Read individual lines without whitespace
	long unsigned int iNx, iNy, numpoints; // Number of X and Y points in mesh 
	int ifile, jfile; // i and j to be read from the file
	double x_store, y_store, u_store, v_store; // Store values to input into array later
	// Read values for nx and ny from first line of meshfile
	mesh >> iNx >> iNy;
	numpoints = iNx * iNy;
	// Open mesh to return later
	Solution test(iNx, iNy);
	if(mesh.is_open()) // Check that file is opened
	{
		// Create arrays to store x, y, u and v from mesh
		double x[numpoints], y[numpoints], u[numpoints], v[numpoints], T[numpoints];
		int cntr = 0; // Counter to index position in the array
		// Read through other lines
		while(true)
		{
			mesh >> ifile >> jfile >> x_store >> y_store >> u_store >> v_store;
			// Store read values into mesh
			x[cntr] = x_store;
			y[cntr] = y_store;
			u[cntr] = u_store;
			v[cntr] = v_store;
			cntr += 1;
			// Exit loop once all lines are read
			if(mesh.eof())
			{
				break;
			}
		}
		// Create empty solution to assign values from 1D double array above
		int index = 0;
		for(int i = 0; i < iNx; i++)
			for(int j = 0; j < iNy; j++)
			{
				test.X(i, j) = x[i * iNy + j];
				test.Y(i, j) = y[i * iNy + j];
				test.u(i, j) = u[i * iNy + j];
				test.v(i, j) = v[i * iNy + j];
			}
		mesh.close();
		}
	else
	{
		std::cout << "Mesh not read properly, please try again. Exiting\n";
	}
	std::cout << "Mesh read successfully, exiting\n";
	return test;
} 

void storeVTKSolution(Solution &s, std::string fileName)
{
	const char * fileChar = fileName.c_str();
	FILE *vtkFile;
	vtkFile = fopen(fileChar, "w");
	unsigned long i, j;
	fprintf(vtkFile, "# vtk DataFile Version 2.0\n");
	fprintf(vtkFile,"TITLE = \"Quad data\"\n");
	fprintf(vtkFile,"ASCII\n");
	fprintf(vtkFile,"DATASET STRUCTURED_GRID\n");
	fprintf(vtkFile,"DIMENSIONS %lu %lu %d\n", s.ny, s.nx, 1);
	fprintf(vtkFile,"POINTS %lu FLOAT\n",s.N);
	for(int i = 0; i < s.nx; i++)
		for(int j = 0; j < s.ny; j++)
		{
			fprintf(vtkFile, "%f %f %f\n", s.X(i, j), s.Y(i, j), 0.0);
		}

	fprintf(vtkFile,"POINT_DATA %lu\n", s.nx * s.ny);
	fprintf(vtkFile,"SCALARS temperature FLOAT 1\n");
	fprintf(vtkFile,"LOOKUP_TABLE default\n");
	for(int i = 0; i < s.nx; i++)
		for(int j = 0; j < s.ny; j++)
		{
			fprintf(vtkFile, "%lf\n", s.T(i, j));
		}
	fprintf(vtkFile,"VECTORS velocity FLOAT\n");
	for(int i = 0; i < s.nx; i++)
		for(int j = 0; j < s.ny; j++)
		{
			fprintf(vtkFile, "%lf %lf %lf\n", s.u(i, j), s.v(i, j), 0.0);
		}
	fclose(vtkFile);
}

/********************************************************************
// Test functions - Used to validate creation of solution class
********************************************************************/

