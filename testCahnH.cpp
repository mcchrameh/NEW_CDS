#include "CahnH.h"
#include "VTK.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <time.h> 
#include "mpi.h"

int main(int argc, char **argv)
{
/*  
  clock_t begin, end;
  double time_spent;
  begin = clock();

   CahnHill2D SEQ(8,8);
   SEQ.Solver2();
   end = clock();
   time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  std::  cout<<""<<std::endl;
   std:: cout<<"Time spent in serial part="<<time_spent<<std::endl;
*/


 double t1=0.0, t2=0.0;
 //std::ifstream file("initialData.dat");

 MPI_Init(&argc, &argv);
 
  
 
   
 Parallel_CahnHill2D  ParVer( 2048,2048 );
 
 t1 = MPI_Wtime();
 //ParVer.ReadFile(file);


ParVer.parallel_solver();
 
 
 t2 = MPI_Wtime();
 printf("MPI_Wtime measured a 1 second sleep to be: %1.2f\n", t2-t1);fflush(stdout);

 MPI_Finalize();

 return 0;

}
