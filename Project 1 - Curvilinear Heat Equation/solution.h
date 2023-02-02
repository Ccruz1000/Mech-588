#include <iostream>
#include <string>
#include <vector>

#ifndef SOLUTION_H

#define  SOLUTION_H

class Operations{
public:
	void alloc1D(double **V, const unsigned long &N, bool isInit = true, double def = 0.0) const;
	void deAlloc1D(double **V, const unsigned long &N) const;
};

class Solution
{
private:


public:

};

// Based on one of the answers to https://stackoverflow.com/questions/21943621/how-to-create-a-contiguous-2d-array-in-c?noredirect=1&lq=1
template <typename T>
class Array_2D
{
private: 
	T *data; // Pointer to data contained within the array 
public:
	// Values to check maximum i and j for indices later on
    unsigned long nx_max;
    unsigned long ny_max;
	Array_2D(unsigned long iNx, unsigned long iNy); // Constructor
	~Array_2D(); // Destructor
	// Member functions
	// Overload the () operator so that we can use (i, j) indices
    T operator()(int i, int j)
    {
    	return data[get_index(i, j)];
    };

    T operator=

    // Return the index so we can use it as (i, j) like normal
    unsigned long get_index(unsigned long i, unsigned long j)
    {
    	if(i >= 0 and i < nx_max and j >= 0 and j < ny_max)
   		{
      		return (i * nx_max + j);
   		}
    }; 
};

// template <typename T>
// class Array_2D
// {
// private:
//     T *data;
// public:
// 	// Values to check maximum i and j for indices later on
//     unsigned long nx_check;
//     unsigned long ny_check;
//     Array_2D(unsigned long row, int column);
//     ~Array_2D();

//     // Overload the () operator so that we can use (i, j) indices
//     T operator()(int iNx, int iNy)
//     {
//         return data_inside[get_index(iNx, iNy)];
//     }

//     int get_index(int index1, int index2){
//         if(index1>=0 and index1<size[0] and index2>=0 and index2<=size[1]){
//             return index1*size[0] + index2;
//         }else{
//             assert("wrong index for array!" == "True");
//         }
//     }

// };

// template <typename T>
// Array_2D<T>::Array_2D(int row, int column)
// {
//     size[0] = row;
//     size[1] = column;
//     data_inside = new T[row*column];
// }

// template <typename T>
// Array_2D<T>::~Array_2D()
// {
//     // Delete array
//     delete[] data_inside;
// }

#endif
