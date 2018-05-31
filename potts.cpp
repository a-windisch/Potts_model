/////////////////////////////////////////////////
//        Program q-state Potts model          //
//           by Andreas Windisch               //
//         andreas.windisch@yahoo.com          //
//   Implemented for Computational Physics I   //
//   taught by Prof. Christof Gattringer       //
//    University of Graz, Austria,  Sep. 2009  //
/////////////////////////////////////////////////  

///////////////////////////////////////////////////////////////////////
// Feel free to use this code as you please. Andreas Windisch, 2018. //
///////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>

using namespace std;

const int L=64;//length along one dimension of the 2-d lattice
const int q=2;//there are q states in a q-state Potts model

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> function prototypes <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void fillspins(int S[L*L],int start);                                   //initializes spins on lattice
void update(int S[L*L],double J,double M,int neib[L*L][4],int nsweeps); //performs nsweep sweeps
void neibinit(int neib[L*L][4]);                                        //initializes neighbors
int  kronecker(int,int);                                                //Kronecker delta: 1 if i=j and 0 otherwise
void magnetization(int S[L*L],double mag[q]);                           //measures magnetization


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin main <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

int main(void)
{
   //variables
   const int start = 0;     //1...hotstart,0...coldstart
   int S[L*L];              //spin array
   int neib[L*L][4];        //in 2d, each point has four neighbors
   double J0 = 0.73;        //J0...starting value
   double JINC = 0.03;      //increment
   double JSTEPS = 10;      //number of steps
   double J = J0;           //J value
   double M = 0;            //ext. field
   double H = 0;            //energy
   int nequi = 500;         //number of equilibration steps
   int nmeas = 1000;        //number of measurements
   int nskip = 30;          //number of updates between measurements
   double aux = 0;          //aux variable used for sorting below
   double mag[q] = {0};     //array for magnetization
   ofstream outfile("potts.dat");  //output stream

   cout << " ///////////////////////////////////////////////// " << endl;
   cout << " //       Program q-state Potts model           // " << endl;
   cout << " //           by Andreas Windisch               // " << endl;
   cout << " //         andreas.windisch@yahoo.com          // " << endl;
   cout << " //   Implemented for Computational Physics I   // " << endl;
   cout << " //   taught by Prof. Christof Gattringer       // " << endl;
   cout << " //    University of Graz, Austria,  Sep. 2009  // " << endl;
   cout << " ///////////////////////////////////////////////// " << endl;
   cout << endl;
   cout << " ///////////////////////////////////////////////////////////////////////" << endl;
   cout << " // Feel free to use this code as you please. Andreas Windisch, 2018. //" << endl;
   cout << " ///////////////////////////////////////////////////////////////////////" << endl;
   cout << endl;
   cout << endl;
   
   
   cout << "Simulating q=" << q << "states Potts Model.\n";
   cout << endl;
   cout << "Parameters:\n";
   cout << "length L ................. " << L << "\n";
   cout << "states q ................. " << q << "\n";
   cout << "Jstart   ................. " << J0 << "\n";
   cout << "Jsteps   ................. " << JSTEPS << "\n";
   cout << "J-incr   ................. " << JINC << "\n";
   cout << "M(ext.)  ................. " << M << "\n";
   cout << "nmeas    ................. " << nmeas << "\n";
   cout << "nskip    ................. " << nskip << "\n";
   cout << "nequi    ................. " << nequi << "\n";
   
   outfile << "Simulating q = " << q << " state Potts model.\n";
   outfile << "Parameters:\n";
   outfile << "length L ................. " << L << "\n";
   outfile << "states q ................. " << q << "\n";
   outfile << "Jstart   ................. " << J0 << "\n";
   outfile << "Jsteps   ................. " << JSTEPS << "\n";
   outfile << "J-incr   ................. " << JINC << "\n";
   outfile << "M(ext.)  ................. " << M << "\n";
   outfile << "nmeas    ................. " << nmeas << "\n";
   outfile << "nskip    ................. " << nskip << "\n";
   outfile << "nequi    ................. " << nequi << "\n";

   outfile.width(18);//parameters for output
   outfile.precision(10);
   outfile.setf(ios_base::fixed);
   outfile << "\n";

   neibinit(neib);//get to know your neighbors
   cout << "Neighbors initialzed.\n";

   fillspins(S,start);//initialize spins
   cout << "Spins initialized.\n";

   for(int i=0; i<JSTEPS; i++) //start a measurement
   {
      J = J0 + i * JINC;
      update( S, J, M, neib, nequi); //equilibrate the system
      cout << "Equilibrated for j = " << J << ".\n";
      for(int imeas=0; imeas<nmeas; imeas++) //measurements
      {
         update( S, J, M, neib, nskip); //'empty' updates
         magnetization( S, mag);        //measure observable magnetization
         for(int i=0; i<q-1; i++)       //sort the result (allowed by symmetry)
         {
            for(int ii=i+1; ii<q; ii++)
            {
               if(mag[i] < mag[ii])
               {
                  aux     = mag[ii];
                  mag[ii] = mag[i];
                  mag[i]  = aux;
               }
            }
         }//end of sorting

         for(int ii=0; ii<q; ii++)
         {
            outfile << mag[ii] << " ";  //write to file
            mag[ii] = 0;
         }

         for(int iii=0; iii<L*L; iii++) //measure energy
         {
            H +=    kronecker(S[iii],S[neib[iii][0]])
                 +  kronecker(S[iii],S[neib[iii][1]]);
         }

         H = 2L*L-H;            //energy observable
         outfile << H << "\n"; //write to file
         H = 0;
      }
      cout << "j = " << J << " done.\n";
   }
   cout << "Done! \n";
   return 0;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end  main <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

void fillspins(int S[L*L],int start) //initialize spins
{
   if(start==1)
   {
      for(int i=0; i<L*L; i++)
      {
         S[i] = rand()%q + 1;//random spin orientation
      }
   }
   else
   {
      for(int ii=0; ii<L*L; ii++)
      {
         S[ii]=1; //all spins 1
      }
   }
}

//-----------------------------------------------------------------------------------------

void neibinit(int neib[L*L][4]) //initialize neighbors
{
   int i, n1p, n1m, n2p, n2m;
   for(int n1=0; n1<L; n1++)
   {
      for(int n2=0; n2<L; n2++)
      {
         n1p=n1+1;
         n1m=n1-1;
         n2p=n2+1;
         n2m=n2-1;
         if(n1p == L)
            n1p = 0;  //periodic boundary conditions
         if(n1m == -1)
            n1m = L-1;
         if(n2p == L)
            n2p = 0;
         if(n2m == -1)
            n2m = L-1;
         i=n1+L*n2; //use a one dimensional index for the two dimensional lattice

         neib[i][0]=n1p+L*n2;   //easten   neighbor
         neib[i][1]=n1+L*n2p;   //northern neighbor
         neib[i][2]=n1m+L*n2;   //western  neighbor
         neib[i][3]=n1+L*n2m;   //southern Nneighbor
      }
   }
}

//-----------------------------------------------------------------------------------------

void update(int S[L*L],double J,double M,int neib[L*L][4],int nsweeps) //update configuration
{
   int offer;//offer a new spin orientation
   double rho,r;//probabilities
   for(int n=0; n<nsweeps; n++)
   {
      for(int i=0; i<L*L; i++)
      {
         offer = rand()%q + 1; //choose a orientation
         rho =    exp(J*(kronecker(offer,S[neib[i][0]])
                + kronecker(offer,S[neib[i][1]])
                + kronecker(offer,S[neib[i][2]])
                + kronecker(offer,S[neib[i][3]])
                - kronecker(S[i],S[neib[i][0]])
                - kronecker(S[i],S[neib[i][1]])
                - kronecker(S[i],S[neib[i][2]])
                - kronecker(S[i],S[neib[i][3]]))
                + M*(kronecker(offer,1)
                - kronecker(S[i],1)));
         r = rand()/(RAND_MAX+1.); // uniform random variable 
         if( r <= rho)
         {
            S[i] = offer; //accept offer
         }
      }
   }
}

//-----------------------------------------------------------------------------------------

int kronecker(int a,int b)//Kronecker delta
{
   if(a==b)
   {
      return 1;
   }
   else
   {
      return 0;
   }
}

//-----------------------------------------------------------------------------------------

void magnetization(int S[L*L],double mag[q])
{
   for(int i=0;i<L*L;i++)//measure magnetization of a configuration
   {
      mag[S[i]-1]++;
   }
}

//-----------------------------------------------------------------------------------------
                                 // END OF FILE
//-----------------------------------------------------------------------------------------
