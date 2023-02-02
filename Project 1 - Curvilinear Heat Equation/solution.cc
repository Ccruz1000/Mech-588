#include "solution.h"

// Constructor to allocate memory for 2D array
template <typename T>
Array_2D<T>::Array_2D(unsigned long iNx, unsigned long iNy)
{
   // Set size array to allow checking of indices later
   nx_max = iNx;
   ny_max = iNy; 
   // Allocate memory for new array
   data = new T[iNx * iNy];
};

// Destructor to deallocate memory for 2D array
template <typename T>
Array_2D<T>::~Array_2D()
{
   // Deallocate memory for array
   delete[] data;
};

int main()
{
   unsigned long nx, ny;
   nx = 3;
   ny = 3; 
   double slow[nx][ny];
   Array_2D<double> fast(nx, ny);

   double cntr = 0.0;

   for(int i = 0; i < nx; i++)
      for(int j = 0; j < ny; j++)
      {
         slow[i][j] = cntr;
         fast(i, j) = cntr;
         cntr += 1.0; 
      }
   // Print slow
   // int numop = 0;
 
   // for(int i = 0; i < nx; i++)
   //    for(int j = 0; j < ny; j++)
   //    {
   //       std::cout << slow[i][j] << " ";
   //       if((numop + 1) % 3 == 0)
   //       {
   //          std::cout << std::endl;
   //       }
   //       numop += 1;
   //    }

   // // Print slow
   // int numop = 0;
 
   // for(int i = 0; i < nx; i++)
   //    for(int j = 0; j < ny; j++)
   //    {
   //       std::cout << fast(i, j) << " ";
   //       if((numop + 1) % 3 == 0)
   //       {
   //          std::cout << std::endl;
   //       }
   //       numop += 1;
   //    }


   return 0;
}

// int main()
// {
//    int i = 3;
//    int j = 3;

//    Array_2D<double> array(i, j);

//    double slow[3][3];
//    double cntr = 0.0;
//    int index = 0;

//    for(int i = 0; i < 3; i++)
//       for(int j = 0; j < 3; j++)
//       {
//          slow[i][j] = cntr;
//          array(i, j);
//          cntr += 1.0;
//          index += 1;
//       }




//    return 0;
// }
// int main()
// {
   // double slow[3][3];
   // double fast[9];
   // double cntr = 0.0;
   // int index = 0;

   // for(int i = 0; i < 3; i++)
   //    for(int j = 0; j < 3; j++)
   //    {
   //       slow[i][j] = cntr;
   //       fast[index] = cntr;
   //       cntr += 1.0;
   //       index += 1;
//       }
//    int numop = 0;
//    std::cout << "Slow\n";
//    for(int i = 0; i < 3; i++)
//       for(int j = 0; j < 3; j++)
//       {
//          std::cout << slow[i][j] << " ";
//          if((numop + 1) % 3 == 0)
//          {
//             std::cout << std::endl;
//          }
//          numop += 1;
//       }

//    std::cout << "\n\n\nFast using single loop\n";
//    for(int i = 0; i < 9; i++)
//    {
//       std::cout << fast[i] << " ";
//       if((i + 1) % 3 == 0)
//       {
//          std::cout << "\n";
//       }
//    }
 
//    std::cout << "\n\n\nFast using nested loop\n";
//    int fastcntr = 0;
//    int numop2 = 0;
//    for(int i = 0; i < 3; i++)
//       for(int j = 0; j < 3; j++)
//       {
//          fastcntr = (3 * i + j);
//          std::cout << fast[fastcntr] << " ";
//          if((numop2 + 1) % 3 == 0)
//          {
//             std::cout << std::endl;
//          }
//          numop2 += 1;
//       }
//    int icheck = 2;
//    int jcheck = 0;

//    std::cout << "Fast: " << fast[(3 * icheck + jcheck)] << std::endl;
//    std::cout << "Slow: " << slow[icheck][jcheck] << std::endl;

// }