#ifndef _CAHNHILL_H
#define _CAHNHILL_H
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "VTK.h"
#include "mpi.h"
class CahnHill2D
{
  protected:
       int nx, ny;
       double delta_t, M, delta_x, delta_y,u,b,k, C1,C2;
       double tau,A,F,v,D,dxx,B,r;
       std::vector<std::vector<double> >phi;
       double **PHI,**PHI_old,**gamma,**Laplacian2;
       int *down1x, *down2x, *up1x, *up2x;
       //void setBC();
      // void computeBoundaryPoints();
       //void computeInnerPoints();
       void UpdateSolution();
       void initialCondition();
       void setLaplacianBase();
       void setIndex();
       double g(double phi);
       void SetSecondLaplacian2();
       void FiniteDifferenceScheme();
      // double gamma(int i, int j);


  public:
        CahnHill2D();
        CahnHill2D(int Nx,int Ny);
        virtual ~CahnHill2D();
       void  display();
       std::string make_output_filename(int index)
           {
              std::ostringstream ss;
              
              ss << "output_" << index << ".dat";
             // ss << "output_" << index << ".vtk";
              return ss.str();
           }

       void result(int count);
       void Solve();
       void Solver2();
       friend void save_vtk(double *U, int Nx, int Ny) ;
       friend double **alloc_2d_int(int rows,int cols);             

};


class Parallel_CahnHill2D: public  CahnHill2D
{
   protected:
            int  nlocalx,nlocaly, remainder, N, Nx,Ny;
            int  rank, size, tag,low, ROW,COL,my2drank, right_nbr, left_nbr, up_nbr,down_nbr;
            int  offset,startrow,endrow, startcol,endcol;
            double **PHI_p,**PHI_old_p,**gamma_p,**Laplacian2_p,*u_lastrow,*u_firstrow, *Boundary_bottom, *Boundary_top,*u_firstcol,*u_lastcol;
            double *Boundary_top1, *Boundary_bottom1, *array_gamma_top,*array_gamma_bottom, *u_store,*bc,*U;
            int  *right1x,*left1x;   
                       
           double *Boundary_col_left,*Boundary_col_right;         
   public: 
     
         Parallel_CahnHill2D();
         Parallel_CahnHill2D(int N1x, int N1y);
        ~Parallel_CahnHill2D();
         void  Initialize_parallel();
         void  ExchangeData(MPI_Comm new_comm, double **array);
      
         void  ComputeLaplacianBase_p(MPI_Comm new_comm);
         void  setIndex_p();
         void  setSecond_laplacian(MPI_Comm new_comm);
         void  FiniteDifferenceScheme_p(MPI_Comm new_comm);
         void  UpdateSolution_p(MPI_Comm new_comm);
         void  WriteToFile_MPI(MPI_Comm new_comm);
         void  SendToMaster(MPI_Comm new_comm,int count); //this  increases communication thus less efficient. 
         void  parallel_solver();
         void  WriteToFile(MPI_Comm new_comm);
         void  ReadFile(std::ifstream& infile);
         MPI_Comm  CreatCartesianTopology();
         void  ExchangeRowCol();
         void  FreeCommunicator(MPI_Comm new_comm);
         
         int* NextNearestNeighborProcess(int myrank,int *array_rank,int num_processes,int *dims,MPI_Comm comm);
         friend double **alloc_2d_int(int rows,int cols);
 
 
   

};










#endif
