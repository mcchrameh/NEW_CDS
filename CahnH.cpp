
#include "CahnH.h"
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <time.h>
#include "mpi.h"
#define Random_min  -0.05 //-0.25 // for mim random number generator
#define Random_max  0.05 //0.25
#define PHI_old(i,j)     PHI_old[i][j]
#define gamma(i,j)       gamma[i][j]
#define Laplacian2(i,j)   Laplacian2[i][j]
#define gamma_p(i,j)     gamma_p[i][j]
#define PHI_old_p(i,j)   PHI_old_p[i][j]
#define Laplacian2_p(i,j) Laplacian2_p[i][j]
using namespace std;




CahnHill2D::CahnHill2D()
      {
        nx =0;
        ny=0;
        //phi=0; 
       
       }
CahnHill2D::CahnHill2D(int Nx, int Ny)
      {
          
          nx=Nx;ny=Ny;delta_x=1.0;delta_y=1.0; delta_t =0.001;M =1.0;b=M;k=M, u=0.5;
          C1=1.0/6.0; C2=1.0/12.0; //C1=1.0/6, C2=1.0/12
          dxx = delta_x;
          B=0.005; D=0.5;A=1.3;v=1.5;tau=0.36;F=0.5;r=0.5;
        for(int i=0;i<nx;i++) 
            {
              phi.push_back(vector<double>(ny));           }
         for(int i=0;i<nx;i++)
            {
             for(int j=0;j<ny;j++)
                 phi[i][j]=0.0; 
            }
          
          PHI=new double*[nx];
          PHI_old = new double*[nx];
          gamma = new double*[nx];
          Laplacian2 = new double*[nx];
         for(int i=0;i<nx; i++)
            {
              PHI[i]=new double[ny];
              PHI_old[i]=new double[ny]; 
              gamma[i]=new double[ny];
              Laplacian2[i]=new double [ny];
            }
          
          for(int i=0;i<nx;i++)
          {
              for(int j=0;j<ny;j++)
              {
                  PHI[i][j]=0.0;
                  PHI_old[i][j]=0.0;
                  gamma[i][j]=0.0;
                  Laplacian2[i][j]=0.0;
                
              }
          }
          
          
          //use these for getting the right indices
          down1x = new int [nx]; memset(down1x, 0, nx*sizeof(int));
          down2x = new int [nx]; memset(down2x, 0, nx*sizeof(int));
          up1x   = new int [nx]; memset(up1x, 0, nx*sizeof(int));
          up2x   = new int [nx]; memset(up2x, 0, nx*sizeof(int));
          cout<<"Constructor called"<<endl;
       }
CahnHill2D::~CahnHill2D()
      {

          for(int i=0;i<nx;i++)
             {
               delete [] PHI[i];
             }
           delete [] PHI;
 
            for(int i=0;i<nx;i++)
             {
               delete [] PHI_old[i];
             }
          delete [] PHI_old;
          
          for(int i=0;i<nx;i++)
          {
              delete [] gamma[i];
          }
          for(int i=0;i<nx;i++)
          {
              delete [] Laplacian2[i];
          }
          
          delete [] gamma;
          delete [] down1x;
          delete [] down2x;
          delete [] up1x;
          delete [] up2x;
          delete [] Laplacian2;
        cout <<"Destructor called"<<endl;

      }
double **alloc_2d_int(int rows,int cols)
        {

  
      double  *data = (double *)malloc((rows+2)*(cols+2)*sizeof(double));
      double  **array= (double **)malloc((rows+2)*sizeof(double*));
        for (int i=0; i<rows+2; i++)
            array[i] = &(data[cols*i]);

    return array;
   }



        
void CahnHill2D:: setIndex()
    {
        for(int s=0; s<nx; s++)
           {
                up1x[s]=s+1;
                up2x[s]=s+2;
                down1x[s]=s-1;
                down2x[s]=s-2;
           }
        
        
        
            down1x[0]=nx-1;
            down2x[0]=nx-2;
            down2x[1]=nx-1;
            up1x[nx-1]=0;
            up2x[nx-1]=1;
            up2x[nx-2]=0;
        
        
        

    }

double CahnHill2D::g(double phi)
{
    double q=0.0;
   // q=(1.0 + tau - A*pow((1.0-2.0*F),2))*PHI_old(i,j)-v*(1.0-2.0*F)*pow(PHI_old(i,j),2)-u*pow(PHI_old(i,j),3);
    //q=(1.0 + tau - A*pow((1.0-2.0*F),2))*phi-v*(1.0-2.0*F)*pow(phi,2)-u*pow(phi,3);
    q = A*phi- (A/3.0)*pow(phi,3); // map from paper
   // cout<<"g(i,j)="<<q<<"\t";
    return q;
    
}
void CahnHill2D::setLaplacianBase()
     {
      //   cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";
         double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;
         setIndex();
         
         for(int i=0;i<nx;i++)
         {
             for(int j=0;j<ny;j++)
             {
                 
                 AP=C1*(PHI_old(up1x[i],j) + PHI_old(down1x[i],j)
                        + PHI_old(i,up1x[j]) + PHI_old(i,down1x[j]));
                 BP=C2*(PHI_old(down1x[i],up1x[j]) + PHI_old(down1x[i],down1x[j])
                        +PHI_old(up1x[i],up1x[j]) + PHI_old(up1x[i],down1x[j]));
                 ATP = AP + BP;
                // gtemp= (1.0 + tau - A*pow((1.0-2.0*F),2))*PHI_old(i,j)-v*(1.0-2.0*F)*pow(PHI_old(i,j),2)-u*pow(PHI_old(i,j),3);
                // cout<<"g(i,j)m="<<gtemp<<"\t";
                 //g_temp =g(i,j);
                 
                 gamma(i,j)=g(PHI_old(i,j))+ D*(ATP-PHI_old(i,j))-PHI_old(i,j);
        //         cout<<"gamma="<<gamma(i,j)<<"\t";
                 
                 
              }
             cout<<""<<endl;
         }
         
     }
void CahnHill2D::SetSecondLaplacian2()
{
    double AP=0.0, BP=0.0;
    setIndex();
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            AP=C1*(gamma(up1x[i],j) + gamma(down1x[i],j)
                   + gamma(i,up1x[j]) + gamma(i,down1x[j]));
            BP=C2*(gamma(down1x[i],up1x[j]) + gamma(down1x[i],down1x[j])
                   +gamma(up1x[i],up1x[j]) + gamma(up1x[i],down1x[j]));
            Laplacian2(i,j) = AP + BP;

            
            
        }
        
    }
    
    
}
void CahnHill2D::FiniteDifferenceScheme()
{
 //   double r=0.5;
    cout<<"in finite difference scheme"<<endl;
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
          // PHI[i][j]= PHI_old(i,j) + delta_t*(Laplacian2[i][j]-gamma(i,j) + B*(1.0)*PHI_old(i,j));
            PHI[i][j]= PHI_old(i,j) - B*(PHI_old(i,j)-1.0 + 2*r)+ gamma(i,j)-Laplacian2(i,j);
            
        }
        
    }
    
}

void CahnHill2D:: display()
        {

           for(int i=0; i<nx;i++)
              {
                for(int j=0;j<ny;j++)
                  {
                    cout<<""<<PHI_old[i][j]<<"\t";
                  }
                cout <<""<<endl;
              }


        }


void CahnHill2D::initialCondition()
     {
  /*     
         for(int i=0; i<nx;i++)
            {
              for(int j=0;j<ny;j++)
                 {
                    double range =Random_max-Random_min;
                    double div =RAND_MAX / range;
                    PHI_old[i][j]=Random_min + (rand()/div);
                     
                 }

            }
*/
 std::ifstream infile("initialData.dat");
   for(int i=0; i<nx;i++)
            {
              for(int j=0;j<ny;j++)
                 {
                    
                   infile>> PHI_old[i][j];

                 }

            }
       

 //       PHI_old[0][0]=0.1;PHI_old[0][1]=0.3;PHI_old[0][2]=-0.2;
 //       PHI_old[1][0]=0.4;PHI_old[1][1]=0.2;PHI_old[1][2]=-0.1;
 //       PHI_old[2][0]=-0.3;PHI_old[2][1]=-0.5;PHI_old[2][2]=0.5;







 

     }


void CahnHill2D::result(int count)
{
    
        FILE *file = fopen(make_output_filename(count).c_str(), "w");
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            
            fprintf(file,"%15.11f \t",PHI[i][j]);
            //fprintf(file,"%15.11f \t",PHI_old[i][j]);
            
        }
        fprintf(file,"\n");
    }
    fclose(file);
    
    
}


    


void CahnHill2D::UpdateSolution()
{
    
    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<ny;j++)
        {
            PHI_old[i][j]=PHI[i][j];
        }
        
    }
    
    
}

void CahnHill2D::Solver2()
{

   initialCondition();
    int count =0;
    double t=0.0;
    while(t<12)
    {
        std::cout<<"in while loop solver"<<std::endl;
        setLaplacianBase();
        SetSecondLaplacian2();
        FiniteDifferenceScheme();
        UpdateSolution();
        if (count==10|| count==0)
            result( count);
       
        t+=1.0;
        count++;
        cout<<"time="<<t<<"\t";
        
    }

  
  
    cout<<""<<endl;
  
    
    
}

void CahnHill2D::Solve()
{
    initialCondition();
    
    //setBC();
    int count =0;
    double t=0.0;
    while(t<13)
    {
        //display();
       // computeBoundaryPoints();

      //  computeInnerPoints();
      //  computeBoundaryPoints();
        UpdateSolution();
       if (count==10 || count==0)
          result( count);
        //t+=delta_t;
        t+=1.0;
        count++;
     }
}


Parallel_CahnHill2D::Parallel_CahnHill2D(int N1x,int N1y):CahnHill2D( N1x, N1y)
{
 //   size =numProc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  //  MPI_Status status;
   //the number of points (x-direction) in each process. Assuming is a square grid
    std::cout<<"size="<<size<<std::endl;
  N =N1y;
     Nx=N1x;
     Ny = N1y;
 
/*
  nlocalx = (int)Nx/size;

    startrow = rank*nlocalx;
    endrow = startrow + nlocalx-1;

    if(rank==size-1)
      {
        endrow=Nx-1;
        nlocalx=endrow-startrow + 1;
      }

*/

 
   
 
  if(Nx%size!=0)
   {
    if(rank<size-1)
     {
      nlocalx = (int)Nx/sqrt(size);
    //  nlocalx = (int)Nx/size;
       
      startrow = rank*nlocalx;
       endrow = startrow + nlocalx-1;
       printf("rank=%d, nlocal=%d\n",rank,nlocalx);

     }
     
    else 
      {
       nlocalx = (int)Nx/sqrt(size);
       // nlocalx = (int)Nx/size;
        startrow = rank*nlocalx;
        endrow=Nx-1;
        nlocalx=endrow-startrow + 1;
     //   printf("rank=%d, nlocal=%d\n",rank,nlocalx);
      }
    }
   else
     {
       nlocalx = (int)Nx/sqrt(size);
     // nlocalx = (int)Nx/size;
     }



 
   if(Ny%size!=0)
     {
      
       if(rank<size-1)
         { 
          nlocaly = (int)Ny/sqrt(size);
        // nlocaly = (int)Ny/size;
         startcol = rank*nlocalx;
         endcol = startcol + nlocaly-1;
         printf("rank=%d, nlocaly=%d\n",rank,nlocaly);
     
         }
        else
           { 
            nlocaly = (int)Ny/sqrt(size);
           // nlocaly = (int)Ny/size;
            startcol = rank*nlocalx;
            endcol=Ny-1;
            nlocaly=endcol-startcol + 1;
     //   printf("rank=%d, nlocal=%d\n",rank,nlocalx);
            }
    
      }
    else
     {
       nlocaly =  (int)Ny/sqrt(size);
      // nlocaly = (int)Ny/size;
     }
 

//     nlocaly=nlocalx;
     tag=201;
   /*--------------------------------------------------------------------
    create local array (matrices) on each process
  ---------------------------------------------------------------------*/
    
       
     PHI_p=new double*[nlocalx+2];
     PHI_old_p=new double*[nlocalx+2];
     gamma_p = new double*[nlocalx+2];
     Laplacian2_p= new double*[nlocalx+2];
    //for transfering boundary points for each process
    u_lastrow=new double[nlocaly]; // change Ny to nlocaly
    u_firstrow=new double[nlocaly];
    u_firstcol=new double[nlocalx];
    u_lastcol =new double[nlocalx];
    Boundary_top = new double[nlocaly];
    Boundary_bottom =new double[nlocaly];
    Boundary_top1 = new double[nlocaly];
    Boundary_bottom1=new double[nlocaly];
    Boundary_col_left =new double[nlocaly];
    Boundary_col_right =new double[nlocalx];

    u_store = new double[(nlocalx+2)*(nlocaly+2)]; // use by each process to store its solution and latter send it to the master for result viewing. 
    right1x = new int[nlocaly];
    left1x = new int [nlocaly];
    array_gamma_top = new double[nlocaly+2];
    array_gamma_bottom = new double[nlocaly+2];
    bc = new double[nlocaly];
    U = new double[Nx*Ny];
    std::memset(U, 0, Nx*Ny*sizeof(double));   
    std::memset(array_gamma_top, 0, (nlocaly+2)*sizeof(double));
    std::memset(array_gamma_bottom, 0, (nlocaly+2)*sizeof(double));
    std::memset(u_store, 0, (nlocalx+2)*(nlocaly+2)*sizeof(double));
    std::memset(bc, 0, nlocaly*sizeof(double));   
    std::memset(left1x, 0, nlocaly*sizeof(int));
    std::memset(u_firstrow, 0,( nlocaly)*sizeof(double));
    std::memset(u_firstcol, 0, (nlocalx)*sizeof(double));
    std::memset(u_lastrow, 0, (nlocalx)*sizeof(double));


    std::memset(u_lastrow, 0, (nlocaly)*sizeof(double));
    std::memset(Boundary_top, 0, (nlocaly)*sizeof(double));
    std::memset(Boundary_bottom, 0, (nlocaly)*sizeof(double));
    std::memset(Boundary_bottom1, 0, (nlocaly)*sizeof(double));
    std::memset(Boundary_top1, 0, (nlocaly)*sizeof(double));


         for(int i=0;i<nlocalx+2; i++)
            {
              PHI_p[i]=new double[nlocaly+2];
              PHI_old_p[i]=new double[nlocaly+2];
              gamma_p[i]=new double[nlocaly+2];
              Laplacian2_p[i]=new double [nlocaly+2];
            }

          for(int i=0;i<nlocalx+2  ;i++)
          { 
             for(int j=0;j<nlocaly+2;j++)
              {
                  PHI_p[i][j]=0.0;
                  PHI_old_p[i][j]=0.0;
                  gamma_p[i][j]=0.0;
                  Laplacian2_p[i][j]=0.0;

              }
          }
   
  std::cout<<"Parallel constructor called"<<std::endl;


}


Parallel_CahnHill2D::~Parallel_CahnHill2D()
{
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //MPI_Comm_size(MPI_COMM_WORLD, &size);

  for(int i =0; i<nlocalx+2; i++)
    {
      delete [] PHI_p[i];
      delete [] PHI_old_p[i];
      delete [] gamma_p[i];
      delete [] Laplacian2_p[i];
    }
     delete []PHI_p;
     delete []PHI_old_p;
     delete []gamma_p;
     delete []Laplacian2_p;
     delete []u_lastrow;
    delete []u_firstrow;
     delete []Boundary_top;
     delete []Boundary_bottom;
     delete []left1x;
     delete []right1x;
     delete []array_gamma_top;
     delete []array_gamma_bottom;
     delete []u_store;
     delete []bc;
     delete []u_firstcol;
     delete []u_lastcol;
     delete []Boundary_col_left;
     delete []Boundary_col_right;
     delete []U;
 
 
//    MPI_Comm_free(&new_comm);

    std::cout<<"destroyed parallel matrices"<<std::endl;
     
}


void Parallel_CahnHill2D::Initialize_parallel()
{

         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MPI_Comm_size(MPI_COMM_WORLD, &size);
 
         srand(time(NULL) + rank);
         printf("Nlocals=(%d,%d)\n",nlocalx,nlocaly);
         std::cout<<"I am rank="<<rank<<std::endl;
         for(int i=1; i<=nlocalx ;i++)
            {
              for(int j=1;j<=nlocaly;j++)
                 {
                    double range =Random_max-Random_min;
                    double div =RAND_MAX / range;
                    PHI_old_p[i][j]=Random_min + (rand()/div);
                            
                  //  std::cout<<PHI_old_p[i][j]<<"\t";

                 }
                  // std::cout<<""<<std::endl;

            }
   
 std::cout<<" intialize data structure"<<std::endl;   
 
}

void Parallel_CahnHill2D:: setIndex_p()
    {   
        for(int s=0; s<N; s++)
           {    
                right1x[s]=s+1;
                
                left1x[s]=s-1;
                
           }


            
            left1x[0]=N-1;
            right1x[N-1]=0;
            
            



    
    }

void Parallel_CahnHill2D::FiniteDifferenceScheme_p(MPI_Comm new_comm)
{
    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);
    //double r=0.5;
    cout<<"in finite difference scheme parallel"<<endl;
    for(int i=1;i<nlocalx+1;i++)
    {
        for(int j=1;j<nlocaly+1;j++)
        {
          // PHI[i][j]= PHI_old(i,j) + delta_t*(Laplacian2[i][j]-gamma(i,j) + B*(1.0)*PHI_old(i,j));
            PHI_p[i][j]= PHI_old_p(i,j) - B*(PHI_old_p(i,j)-1.0 + 2*r)+ gamma_p(i,j)-Laplacian2_p(i,j);
            
        }
        
    }
    
}

void Parallel_CahnHill2D::UpdateSolution_p(MPI_Comm new_comm)
    {
       MPI_Comm_rank(new_comm, &my2drank);
       MPI_Comm_size(new_comm, &size);      
       for(int i=1;i<nlocalx+1;i++)
          {
           for(int j=1;j<nlocaly+1;j++)
              {
                PHI_old_p[i][j]=PHI_p[i][j];
             }

          }

    }
void Parallel_CahnHill2D::FreeCommunicator(MPI_Comm new_comm)
     {
       MPI_Comm_free(&new_comm);
  
     }
    //

void Parallel_CahnHill2D::ExchangeRowCol()
    {
        }
void Parallel_CahnHill2D::ExchangeData(MPI_Comm new_comm, double **array)
{
       MPI_Status status;
       
       double value[2]; value[0]=0.0; value[1]=0.0; //used  for point to point comm to fill up ghost cells
       double value1[2]; value1[0]=0.0; value1[1]=0.0;
    //-----------Variables for creating Cartesian Topology----------/ 
   int  ROW=0,COL=1;

     //int periods[2];
     //periods[0]=1;periods[1]=1;//to enable wrap around
    int coords[2];
     int dims[2];
     dims[ROW] =sqrt(size);dims[COL]=sqrt(size); //for square topology
  //---------------------------------------/
 

    //  new_comm = CreatCartesianTopology();
     

         //Exchange left and right
       //start with even column number
       MPI_Comm_rank(new_comm,&my2drank);
       MPI_Cart_coords(new_comm,my2drank,2,coords);

       if(coords[1]%2==0)
         {

           for(int i=1;i<=nlocalx;i++)
             {
              
               u_lastcol[i]=array[i][nlocaly];
             }
            MPI_Comm_rank(new_comm,&my2drank);
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
          //  printf("my coord=(%d,%d) rank=%d, neighbors=(%d,%d)\n",coords[0],coords[1],my2drank,left_nbr,right_nbr);
            MPI_Send( u_lastcol, nlocalx+1, MPI_DOUBLE, right_nbr, tag, new_comm );
            MPI_Recv( u_firstcol, nlocalx+1,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
            for(int i=1;i<=nlocalx;i++)
                {
                  array[i][nlocaly+1] =u_firstcol[i];
                //  printf("last col: Phi_old_p[%d][%d]=%lf\n",i,nlocaly,PHI_old_p[i][nlocaly+1]);
                }
        }
        else
          {

           for(int i=1;i<=nlocalx;i++)
             {
               u_firstcol[i]=array[i][1];
               
             }
           MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
           MPI_Recv( u_lastcol, nlocalx+1,  MPI_DOUBLE, left_nbr, MPI_ANY_TAG, new_comm, &status );
           MPI_Send( u_firstcol, nlocalx+1, MPI_DOUBLE, left_nbr, tag, new_comm );
           for(int i=1;i<=nlocalx;i++)
                {
                  array[i][0] =u_lastcol[i];
                //  printf("first col: Phi_old_p[%d][0]=%lf\n",i,PHI_old_p[i][0]);
                }


          }
    
      if(coords[1]%2==0)
        {
           for(int i=1;i<=nlocalx;i++)
             {
               u_firstcol[i]=array[i][1];
               
             }
           MPI_Comm_rank(new_comm,&my2drank);
           MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
           MPI_Send( u_firstcol, nlocalx+1, MPI_DOUBLE, left_nbr, tag, new_comm );
           MPI_Recv( u_lastcol, nlocalx+1,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
           for(int i=1;i<=nlocalx;i++)
                {
                  array[i][0] =u_lastcol[i];
         //         printf("rank=%d::recv last col: Phi_old_p[%d][0]=%lf\n",my2drank,i,PHI_old_p[i][0]);
                }
         }
      else
        {

           for(int i=1;i<=nlocalx;i++)
             {
               
               u_lastcol[i]=array[i][nlocaly];
             }
           MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr );
           MPI_Recv( u_firstcol, nlocalx+1,  MPI_DOUBLE, right_nbr, MPI_ANY_TAG, new_comm, &status );
           MPI_Send( u_lastcol, nlocalx+1, MPI_DOUBLE, right_nbr, tag, new_comm );
           for(int i=1;i<=nlocalx;i++)
                {
                  array[i][nlocaly+1] =u_firstcol[i];
         //         printf("rank=%d::recv first col: Phi_old_p[%d][%d]=%lf\n",my2drank,i,nlocaly+1,PHI_old_p[i][nlocaly+1]);
                }

  

        }
     
        //Exchange rows  
      if (coords[0]%2==0) //even rows send data downwards
         {

             for(int j=1;j<=nlocaly;j++)
             {

              
               u_lastrow[j-1]=array[nlocalx][j];
              }


            MPI_Comm_rank(new_comm,&my2drank);
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Send( u_lastrow, nlocaly, MPI_DOUBLE, down_nbr, tag, new_comm );
            MPI_Recv( u_firstrow, nlocaly,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            for(int j=1;j<=nlocaly;j++)
                {
                  array[nlocalx+1][j] =u_firstrow[j-1];
        //         printf("rank=%d::recv first row: Phi_old_p[%d][%d]=%lf\n",my2drank,j,nlocalx+1,PHI_old_p[nlocalx+1][j]);
                }
          }
        else
         {


              for(int j=1;j<=nlocaly;j++)
             {

               u_firstrow[j-1]=array[1][j];
              
              }


             MPI_Comm_rank(new_comm,&my2drank);
             MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
             MPI_Recv( u_lastrow, nlocaly,  MPI_DOUBLE, up_nbr, MPI_ANY_TAG, new_comm, &status );
             MPI_Send( u_firstrow, nlocaly, MPI_DOUBLE, up_nbr, tag, new_comm );
            for(int j=1;j<=nlocaly;j++)
                {
                  array[0][j] =u_lastrow[j-1];
          //        printf("rank=%d::recv last row: Phi_old_p[0][%d]=%lf\n",my2drank,j,PHI_old_p[0][j]);
                }



         }
       if(coords[0]%2==0) //even rows send data upwards
         {
             for(int j=1;j<=nlocaly;j++)
             {

               u_firstrow[j-1]=array[1][j];
              
              }

             MPI_Comm_rank(new_comm,&my2drank);
             MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
             MPI_Send( u_firstrow, nlocaly, MPI_DOUBLE, up_nbr, tag, new_comm );
             MPI_Recv( u_lastrow, nlocaly,  MPI_DOUBLE,up_nbr, tag, new_comm, &status );
             for(int j=1;j<=nlocaly;j++)
                {
                  array[0][j] =u_lastrow[j-1];
                 // printf("rank=%d::recv boundary: from:: up_nbr=%d:: boundary[%d]=%lf\n",my2drank,up_nbr,j,Boundary_bottom[j]);
                }
         }
        else if(coords[0]%2!=0)
        {
            for(int j=1;j<=nlocaly;j++)
             {

               u_lastrow[j-1]=array[nlocalx][j];
        //       printf("rank=%d, Boundary_bottom[%d]=%lf\n",my2drank,j,Boundary_bottom[j]);
              }


            MPI_Comm_rank(new_comm,&my2drank);
            printf("I am rank=%d\n",my2drank);
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Recv( u_firstrow, nlocaly,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( u_lastrow,  nlocaly, MPI_DOUBLE, down_nbr, tag, new_comm );
            for(int j=1;j<=nlocaly;j++)
                {
                  array[nlocalx+1][j] =u_firstrow[j-1];
                 // printf("rank=%d::recv first row: Phi_old_p[%d][%d]=%lf\n",my2drank,nlocalx+1,j,PHI_old_p[nlocalx+1][j]);
                }

         }
     std::cout<<"-----check point"<<std::endl;
    MPI_Barrier(new_comm);
    MPI_Cart_coords(new_comm,my2drank,2,coords);


     //point to point communication to completely fill up ghost points
     if(coords[0]==0 && coords[1]==0) //processes at the edge with exchange data first, ie (0,0) with (0,dim[COL])
       {
       
         
           value[0]=array[0][1];
           MPI_Comm_rank(new_comm,&my2drank);
           MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
           MPI_Send( &value[0],  1, MPI_DOUBLE, left_nbr, tag, new_comm );
           MPI_Recv( &value1[0], 1,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
                  array[0][0]=value1[0]; 
         }
        else if (coords[0]==0 && coords[1]==dims[COL]-1)
          {
              value1[0]=array[0][nlocaly];
              MPI_Comm_rank(new_comm,&my2drank);
              MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
              MPI_Recv( &value[0], 1,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
              MPI_Send( &value1[0],  1, MPI_DOUBLE, right_nbr, tag, new_comm );
              array[0][nlocaly+1]=value[0];  
          }
       else
         {
           
         }
      
    //Exchange Bottom edge points 
    if ((coords[0]=dims[ROW]-1) && (coords[1]==0))
       {
            
           value[0]=array[nlocalx+1][1];
           MPI_Comm_rank(new_comm,&my2drank);
           MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
           MPI_Send( &value[0],  1, MPI_DOUBLE, left_nbr, tag, new_comm );
           MPI_Recv( &value1[0], 1,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
                  array[nlocalx+1][0]=value1[0];
        //   printf("value1=%lf, phi[%d][0]=%lf\n",value1[0],nlocalx+1,array[nlocalx+1][0]); 
       }
     else if ((coords[0]==dims[ROW]-1) && (coords[1]==dims[COL]-1))
          {
             
             value1[0]=array[nlocalx+1][nlocaly];
              MPI_Comm_rank(new_comm,&my2drank);
              MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
              MPI_Recv( &value[0], 1,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
              MPI_Send( &value1[0],  1, MPI_DOUBLE, right_nbr, tag, new_comm );
              array[nlocalx+1][nlocaly+1]=value[0];  
          }
       else
         {

         }

       //   printf("myrank=%d, coords=(%d,%d)\n",my2drank,coords[0],coords[1]);
     
MPI_Barrier(new_comm);

  //Exchange right most and left most edge points between adjacent processors
  // Start with even rows then with odd rows
  MPI_Cart_coords(new_comm,my2drank,2,coords);

  if(coords[0]%2==0)  //even rows exchange
    {
       MPI_Comm_rank(new_comm,&my2drank);
       if((my2drank%2==0) && (coords[1]!=dims[COL]-1))
         {
            value[0]=array[0][nlocaly];
            value[1]=array[nlocalx+1][nlocaly];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Send( value,  2, MPI_DOUBLE, right_nbr, tag, new_comm );
            MPI_Recv( value1, 2,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
            array[0][nlocaly+1]=value1[0];
            array[nlocalx+1][nlocaly+1]=value1[1];
         }
        else
         {
            value1[0]=array[0][1];
            value1[1]=array[nlocalx+1][1];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Recv( value, 2,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( value1,  2, MPI_DOUBLE, left_nbr, tag, new_comm );

            array[0][0]=value[0];
            array[nlocalx+1][0]=value[1];
        }
     // MPI_Comm_rank(new_comm,&my2drank);
      MPI_Cart_coords(new_comm,my2drank,2,coords);
      if((my2drank%2!=0)&& (coords[1]!=dims[COL]-1)) //for odd ranks
        {
            value[0]=array[0][nlocaly];
            value[1]=array[nlocalx+1][nlocaly];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Send( value,  2, MPI_DOUBLE, right_nbr, tag, new_comm );
            MPI_Recv( value1, 2,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
            array[0][nlocaly+1]=value1[0];
            array[nlocalx+1][nlocaly+1]=value1[1];

        }
     else if ((my2drank%2==0)&&(coords[1]!=0))
       {
           
            value1[0]=array[0][1];
           // value1[1]=array[nlocalx+1][nlocaly];
            value1[1]=array[nlocalx+1][1];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Recv( value, 2,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( value1,  2, MPI_DOUBLE, left_nbr, tag, new_comm );
            array[0][0]=value[0];
            array[0][nlocaly+1]=value[1];
      }

}
else //odd rows exchange
   {
     MPI_Comm_rank(new_comm,&my2drank);
     if(my2drank%2==0)
         {
            value[0]=array[0][nlocaly];
            value[1]=array[nlocalx+1][nlocaly];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Send( value,  2, MPI_DOUBLE, right_nbr, tag, new_comm );
            MPI_Recv( value1, 2,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
            array[0][nlocaly+1]=value1[0];
            array[nlocalx+1][nlocaly+1]=value1[1];
         }
      else
         {

            value1[0]=array[0][1];
            //value1[1]=array[nlocalx+1][nlocaly];
            value1[1]=array[nlocalx+1][1];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Recv( value, 2,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( value1,  2, MPI_DOUBLE, left_nbr, tag, new_comm );

            array[0][0]=value[0];
            array[nlocalx+1][0]=value[1];
        }
     // MPI_Comm_rank(new_comm,&my2drank);
      MPI_Cart_coords(new_comm,my2drank,2,coords);
      if((my2drank%2!=0)&& (coords[1]!=dims[COL]-1))
        {
            value[0]=array[0][nlocaly];
            value[1]=array[nlocalx+1][nlocaly];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Send( value,  2, MPI_DOUBLE, right_nbr, tag, new_comm );
            MPI_Recv( value1, 2,  MPI_DOUBLE,right_nbr, MPI_ANY_TAG, new_comm, &status );
            array[0][nlocaly+1]=value1[0];
            array[nlocalx+1][nlocaly+1]=value1[1];

        }
     else if ((my2drank%2==0)&&(coords[1]!=0))
       {

            value1[0]=array[0][1];
          //  value1[1]=array[nlocalx+1][nlocaly];
            value1[1]=array[nlocalx+1][1];
            MPI_Cart_shift( new_comm, 1, 1, &left_nbr, &right_nbr);
            MPI_Recv( value, 2,  MPI_DOUBLE,left_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( value1,  2, MPI_DOUBLE, left_nbr, tag, new_comm );
            array[0][0]=value[0];
            array[0][nlocaly+1]=value[1];
      }




   }

// implement odd ranks for more than 2x2 grid


MPI_Cart_coords(new_comm,my2drank,2,coords);
if(coords[1]==0)
  {
       MPI_Cart_coords(new_comm,my2drank,2,coords);

       if(coords[0]%2==0 && coords[0]!=dims[ROW]-1)
         {
            value[0]=array[nlocalx][0];
         
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Send( &value[0],  1, MPI_DOUBLE, down_nbr, tag, new_comm );
            MPI_Recv( &value1[0], 1,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            array[nlocalx+1][0]=value1[0];

         }
       else
         {
            value1[0]=array[1][0];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Recv( &value[0], 1,  MPI_DOUBLE,up_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( &value1[0],  1, MPI_DOUBLE, up_nbr, tag, new_comm );
              array[0][0]=value[0];
         }

      //processes on odd row number send downwards
     
     
       MPI_Cart_coords(new_comm,my2drank,2,coords);
       if((coords[0]%2!=0)&&(coords[0]!=dims[ROW]-1))
         {
             //xcoord = coords[0];
            value[0]=array[nlocalx][0];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Send( &value[0],  1, MPI_DOUBLE, down_nbr, tag, new_comm );
            MPI_Recv( &value1[0], 1,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            array[nlocalx+1][0]=value1[0];
         }
      else if(coords[0]%2==0 && coords[0]!=0 )
        {
            value1[0]=array[1][0];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Recv( &value[0], 1,  MPI_DOUBLE,up_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( &value1[0],  1, MPI_DOUBLE, up_nbr, tag, new_comm );
              array[0][0]=value[0];
        }
      else
        {

        }
} 

 

//Right boundary processes
MPI_Cart_coords(new_comm,my2drank,2,coords);
if(coords[1]==dims[COL]-1)
  {
    if((coords[0]%2==0)&& (coords[0]!=dims[ROW]-1))
      {
            value[0]=array[nlocalx][nlocaly+1];

            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Send( &value[0],  1, MPI_DOUBLE, down_nbr, tag, new_comm );
            MPI_Recv( &value1[0], 1,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            array[nlocalx+1][nlocaly+1]=value1[0];
      }
     else
      {
            value1[0]=array[1][nlocaly+1];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Recv( &value[0], 1,  MPI_DOUBLE,up_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( &value1[0],  1, MPI_DOUBLE, up_nbr, tag, new_comm );
            array[0][nlocaly+1]=value[0];
      }
     MPI_Cart_coords(new_comm,my2drank,2,coords);
       if((coords[0]%2!=0)&&(coords[0]!=dims[ROW]-1))
         {
            
            value[0]=array[nlocalx][nlocaly+1];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Send( &value[0],  1, MPI_DOUBLE, down_nbr, tag, new_comm );
            MPI_Recv( &value1[0], 1,  MPI_DOUBLE,down_nbr, MPI_ANY_TAG, new_comm, &status );
            array[0][nlocaly+1]=value1[0];
         }
      else if((coords[0]%2==0) && (coords[0]!=0) )
        {
            value1[0]=array[1][nlocaly+1];
            MPI_Cart_shift( new_comm, 0, 1, &up_nbr, &down_nbr);
            MPI_Recv( &value[0], 1,  MPI_DOUBLE,up_nbr, MPI_ANY_TAG, new_comm, &status );
            MPI_Send( &value1[0],  1, MPI_DOUBLE, up_nbr, tag, new_comm );
              array[0][nlocaly+1]=value[0];
        }
      else
        {

        }

  }
  
 MPI_Barrier(new_comm);

/*
       if(my2drank==0)
          {
        
           printf("I am rank=%d\n", my2drank);
           for(int i=0;i<nlocalx+2;i++)
              {
               for(int j=0;j<nlocaly+2;j++)
                 {

                   printf("Phi-old_p[%d][%d]=%lf\t",i,j,array[i][j]);
                  }
                printf("\n");
               }
             }*/



  

}
int* Parallel_CahnHill2D::NextNearestNeighborProcess(int myrank, int *array_rank,int num_processes,int *dims,MPI_Comm comm)
    {
     
      int new_rank=0;
      int mycoords[2], process_coords[2];
   //   std::memset(array_rank, 0, 4*sizeof(int));

      if (num_processes==4 )
        {
          
          
           int mycoords[2], process_coords[2]={0,0};

           MPI_Cart_coords(comm,myrank,2,mycoords);
           process_coords[0]=mycoords[0]+1;
           process_coords[1]=mycoords[1]+1;
           if (process_coords[0]>dims[0]|| process_coords[0]>dims[1])
              {
           //     array_rank=MPI_PROC_NULL;
              }
            else
              {
               MPI_Cart_rank(comm,process_coords,&new_rank);
               array_rank[0]=new_rank;
               new_rank=0;
               return array_rank;
              }
           
          
            
        }
     else
       {
        MPI_Cart_coords(comm,myrank,2,mycoords);
        process_coords[ROW]=mycoords[ROW]+1;
        process_coords[COL]=mycoords[COL]+1;
        MPI_Cart_rank(comm,process_coords,&new_rank);
        array_rank[0]=new_rank;
        new_rank=0;
        process_coords[ROW]=mycoords[ROW]-1;
        process_coords[COL]=mycoords[COL]+1;
        MPI_Cart_rank(comm,process_coords,&new_rank);
        array_rank[1]=new_rank;
        new_rank=0;
        process_coords[ROW]=mycoords[ROW]-1;
        process_coords[COL]=mycoords[COL]-1;
        MPI_Cart_rank(comm,process_coords,&new_rank);
        array_rank[2]=new_rank;
        new_rank=0;
        process_coords[ROW]=mycoords[ROW]+1;
        process_coords[COL]=mycoords[COL]-1;
        MPI_Cart_rank(comm,process_coords,&new_rank);
        array_rank[3]=new_rank;
       }
     
 

  
      return array_rank; 

    }
MPI_Comm Parallel_CahnHill2D:: CreatCartesianTopology()
    {
       //MPI_Status status;
       MPI_Comm   new_comm;
      int  ROW=0,COL=1;
    
       int periods[2];
        periods[0]=1;periods[1]=1;//to enable wrap around
      int coords[2];
      int dims[2];
      dims[ROW] =sqrt(size);dims[COL]=sqrt(size); //for square topology
     //---------------------------------------/
     MPI_Comm_size( MPI_COMM_WORLD, &size );
     MPI_Comm_rank(MPI_COMM_WORLD,&rank);
     MPI_Cart_create( MPI_COMM_WORLD, 2,dims, periods, 0, &new_comm ); //5th parameter set to 0 inorder to keep  the rank order in the  new communicator
     MPI_Comm_rank(new_comm,&my2drank);
     MPI_Cart_coords(new_comm,my2drank,2,coords);
   
     MPI_Comm_size( new_comm, &size );
    return new_comm;

    }

void Parallel_CahnHill2D::ComputeLaplacianBase_p( MPI_Comm new_comm)
    {
      // ExchangeData(PHI_old_p);
       MPI_Comm_rank(new_comm, &my2drank);
       MPI_Comm_size(new_comm, &size);
       //ExchangeData(new_comm, PHI_old_p);    
       //cout<<"gamma(0,0)="<<gamma(0,0)<<"\t";
        double AP=0.0, BP=0.0,  ATP=0.0;//gtemp=0.0;
        

         for(int i=1;i<=nlocalx;i++)
         {
             for(int j=1;j<=nlocaly;j++)
             {

                 AP=C1*(PHI_old_p(i+1,j) + PHI_old_p(i-1,j)
                        + PHI_old_p(i,j+1) + PHI_old_p(i,j-1));
                 BP=C2*(PHI_old_p(i-1,j+1) + PHI_old_p(i-1,j-1)
                        +PHI_old_p(i+1,j+1) + PHI_old_p(i+1,j-1));
                 ATP = AP + BP;
              
                 gamma_p(i,j)=CahnHill2D::g(PHI_old_p(i,j))+ D*(ATP-PHI_old_p(i,j))-PHI_old_p(i,j);
              
               }
          }
   
}
void Parallel_CahnHill2D::setSecond_laplacian(MPI_Comm new_comm)
{

    MPI_Comm_rank(new_comm, &my2drank);
    MPI_Comm_size(new_comm, &size);
   // ExchangeData(new_comm,gamma_p);

    double AP=0.0, BP=0.0;

    for(int i=1;i<=nlocalx;i++)
    {   
        for(int j=1;j<=nlocaly;j++)
        {   
            AP=C1*(gamma_p(i+1,j) + gamma_p(i-1,j)
                   + gamma_p(i,j+1) + gamma_p(i,j-1));
            BP=C2*(gamma_p(i-1,j+1) + gamma_p(i-1,j-1)
                   +gamma_p(i+1,j+1) + gamma_p(i+1,j-1));
            Laplacian2_p(i,j) = AP + BP;
         
         
         
        }
     
    }

 
}

void Parallel_CahnHill2D::WriteToFile_MPI(MPI_Comm new_comm)
    {
       
       char const *fmt ="%20.8f";
       char const *endfmt="%20.8f\n";
       MPI_Status status;
       MPI_File fileh;
       //MPI_Offset disp;
       MPI_Datatype num_as_string;
       MPI_Datatype localarray;
       const int charspernum=16;
       MPI_Comm_rank(new_comm, &rank);
       MPI_Comm_size(new_comm, &size);
      
       MPI_Type_contiguous(charspernum, MPI_CHAR,&num_as_string);
       MPI_Type_commit(&num_as_string);
       double **data_new;
       data_new= alloc_2d_int(nlocalx,nlocaly);
       for(int i=1;i<nlocalx+1;i++)
          {
          for(int j=1;i<nlocaly+1;j++)
             {

              data_new[i][j]=PHI_p[i][j];

             }


          }
       /*----convert data into txt-----*/
       char *data_as_txt=(char *)malloc((nlocalx+2)*(nlocaly+2)*charspernum*sizeof(char));
       int count =0;
       
       for(int i=1;i<nlocalx+1;i++)
          {
           for(int j=1;j<nlocaly;j++)
              {
                sprintf(&data_as_txt[count*charspernum],fmt,data_new[i][j]);
                count++;
              }
             sprintf(&data_as_txt[count*charspernum],endfmt,data_new[i][nlocaly]);
             count++;

          }
       
       char  mpifilename[]="mpifile";
       char  NATIVE[]="native"; 
       int globalsizes[2]={Nx, Ny};
       int localsizes[2]={nlocalx,nlocaly};
       int starts[2]={rank*nlocalx,rank*nlocaly};
       int order = MPI_ORDER_C;
       MPI_Type_create_subarray(2,globalsizes,localsizes,starts,order,MPI_DOUBLE,&localarray);
       MPI_Type_commit(&localarray);
 
       MPI_File_open(MPI_COMM_WORLD, mpifilename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fileh);
       MPI_File_set_view(fileh, 0, MPI_DOUBLE, localarray, NATIVE, MPI_INFO_NULL);
       MPI_File_write_all(fileh,PHI_p,(nlocalx+2)*(nlocaly+2),MPI_DOUBLE,&status);
      // MPI_File_set_view(fileh, rank*nlocalx*nlocaly*sizeof(double), MPI_DOUBLE, MPI_DOUBLE, NATIVE, MPI_INFO_NULL);
      // MPI_File_write_all(fileh,data_as_txt,nlocalx*nlocaly,num_as_string,&status);
       MPI_File_close(&fileh);
       MPI_Type_free(&localarray);
       MPI_Type_free(&num_as_string);
       free(data_new[0]);
       free(data_new);


     


   }

void Parallel_CahnHill2D::SendToMaster(MPI_Comm new_comm,int count)
    {
        //MPI_Status status;      
        MPI_Comm_rank(new_comm, &rank);
        MPI_Comm_size(new_comm, &size);
       

        for(int i=1;i<nlocalx+1;i++)
           {
                for(int j=1;j<nlocaly+1;j++)
                   {
                           u_store[i*nlocaly+j]=PHI_p[i][j];
                   }
                }

        if(rank==0)
          {
         //         U = new double[Nx*Ny];
           //       std::memset(U, 0, Nx*Ny*sizeof(double));
          }
       MPI_Gather( u_store, nlocalx*nlocaly, MPI_DOUBLE, U, nlocalx*nlocaly, MPI_DOUBLE, 0, new_comm);
     
        if (rank==0)
           {
                  FILE *file = fopen(make_output_filename(count).c_str(), "w");
                   for(int i=1;i<Nx;i++)
                      {
                            for(int j=1;j<Ny;j++)
                              {

                                 fprintf(file,"%15.11f \t",U[i*Ny+j]);

                              }
                              fprintf(file,"\n");
                     }
                fclose(file);

           //      delete []U;
           }
    // delete []U;

   //  save_vtk(U, Nx,Ny);
    }
void Parallel_CahnHill2D::parallel_solver()
    {
       Initialize_parallel();
     //  ReadFile(output_0.dat);

       MPI_Comm new_comm;
       new_comm = CreatCartesianTopology();
     


       int count =0;
       double t=0.0;
       while(t<1001)
           {
            // std::cout<<"in while loop solver"<<std::endl;
             ExchangeData(new_comm,PHI_old_p);
             ComputeLaplacianBase_p(new_comm);
             ExchangeData(new_comm,gamma_p);
             setSecond_laplacian(new_comm);
             ExchangeData(new_comm,Laplacian2_p);
             FiniteDifferenceScheme_p(new_comm);
             UpdateSolution_p(new_comm);
             t+=1.0;
             count++;
       
           }
          //SendToMaster(new_comm, count);
           WriteToFile( new_comm);
          // WriteToFile_MPI(new_comm);


            cout<<"count="<<count<<endl;
    
    MPI_Barrier(new_comm);
    
  /*
     if(my2drank==0)
          {
        
           printf("I am rank=%d\n", my2drank);
           for(int i=0;i<nlocalx+2;i++)
              {
               for(int j=0;j<nlocaly+2;j++)
                 {

                   printf("Laplacian2_p[%d][%d]=%lf\t",i,j,Laplacian2_p(i,j));
                 //  printf("gamma_p[%d][%d]=%lf\t",i,j,gamma_p(i,j));
                  //  printf("Phi-old_p[%d][%d]=%lf\t",i,j,PHI_old_p(i,j));

                  }
                printf("\n");
               }
             }
*/
 
   FreeCommunicator(new_comm);
    }


void Parallel_CahnHill2D::WriteToFile(MPI_Comm new_comm)
    {
     //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     //  MPI_Comm_size(MPI_COMM_WORLD, &size);
     MPI_Comm_rank(new_comm, &rank);
    MPI_Comm_size(new_comm, &size);

      char fname[200];
      sprintf(fname,"solution%d.dat",rank);
      FILE* fp;
      fp = fopen(fname,"w");
      for(int i=1;i<nlocalx+1;i++)
         {
           for(int j =1; j<nlocaly+1;j++)
             {
                fprintf(fp,"%20.16f", PHI_p[i][j]);
               
             }
           fprintf(fp,"\n");
   
        }
     fclose(fp);  


    }

void Parallel_CahnHill2D::ReadFile(std::ifstream& infile)
    {

 

            MPI_Comm new_comm;
            new_comm= CreatCartesianTopology();
         
            MPI_Status status;
            MPI_Comm_rank(new_comm, &rank);
            MPI_Comm_size(new_comm, &size);
            tag = 201;
            low=rank*nlocalx;
            int coords[2];
           //high=low+nlocalx;
          double **u1, *u2;
          u2 =new double[Nx*Ny];
          
           
           //u2 =new double[Nx*Ny];
           std::memset(u2, 0, Nx*Ny*sizeof(double));
           u1 = new double*[Nx+2];
           for(int i=0;i<Nx+2;i++)
              {       
                 u1[i]=new double[Ny+2];
              }
 
               for(int i=0;i<Nx+2;i++)
                  {
                   for(int j=0;j<Ny+2;j++)
                      {
                        u1[i][j]=0.0;
                       }
                   }
                   
         // u1 =alloc_2d_int(Nx,Ny);

       
      
       
      //if(rank==0)
         
          for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {   
                         
            //           infile>>u2[i*Ny+j];
                      }
          
                    }
     
  
          for(int i=0;i<Nx;i++)
                 {
                  for(int j=0;j<Ny;j++)
                     {
                         
                       infile>>u1[i][j];
                      }

                    }
   
  // MPI_Bcast(u2,Ny*Nx,MPI_DOUBLE,0,new_comm);
           
               MPI_Barrier(new_comm);
     //MPI_Comm_rank(new_comm,&rank);
     MPI_Cart_coords(new_comm,rank,2,coords);

             for(int i=1;i<=nlocalx;i++)
                {
                  
                   for(int j=1;j<=nlocaly;j++)
                       {
                        
                        
                           PHI_old_p[i][j]=u1[coords[0]*nlocalx + i-1][coords[1]*nlocaly + j-1];
                          


                      
                         printf("rank=%d, phi_old_p[%d][%d]=%lf\n",rank,i,j,PHI_old_p[i][j]);
                       }
            
               }
             for(int i =0;i<Nx+2;i++)  
            {     
              delete []u1[i];
            }    
               
            delete []u1;
         // free(u1[0]);
        //  free(u1);
            delete []u2;
            FreeCommunicator( new_comm);
                                                                                                                                                                                                                                                                                                                                  

     }
