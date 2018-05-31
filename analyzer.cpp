/////////////////////////////////////////////////
//  Program analyzer for q-state Potts model   //
//           by Andreas Windisch               //
//         andreas.windisch@yahoo.com          //
//   Implemented for Computational Physics I   //
//   taught by Prof. Christof Gattringer       //
//    University of Graz, Austria,  Sep. 2009  //
/////////////////////////////////////////////////  

///////////////////////////////////////////////////////////////////////
// Feel free to use this code as you please. Andreas Windisch, 2018. //
///////////////////////////////////////////////////////////////////////


#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

istream &nxl(istream &is) //advance line by line
{
   char xcr[257];
   is.getline(xcr,256);
   return is;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> begin main <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


int main(void)
{
   //variables
   string dummy;                               //variable to take up irrelevant parts of input
   int L = 0;                                  //length along one dimension
   int q = 0;                                  //number of states in q state Potts model
   double J0 = 0;                              //start value of J
   double JSTEPS = 0;                          //number of steps
   double JINC = 0;                            //increment
   double M = 0;                               //external field
   double H = 0;                               //energy
   double epsH = 0;                            //error for H
   double C = 0;                               //heat capacity
   double aux = 0;                             //aux variable used for reading data below
   double j = 0;                               //value of j
   double mag[8] = {0};                        //magnetization array
   double epsm[8] = {0};                       //error for magnetization
   double chi[8] = {0};                        //susceptibility
   double h[2000] = {0};                       //means of h
   int nmeas = 0;                              //number of measurements
   int nskip = 0;                              //number of updates between measurements
   int nequi = 0;                              //number of equilibration steps

   cout << " ///////////////////////////////////////////////// " << endl;
   cout << " //  Program analyzer for q-state Potts model   // " << endl;
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


   cout << "Reading parameters...\n";          //read in raw data produced by potts.cpp
   cout << endl; 

   ifstream infile("potts.dat");               //open two streams
   ifstream infile2("potts.dat");
   infile >> nxl;                              //nextline
   infile2 >> nxl;                             //nextline for stream 2
   infile >> nxl;
   infile2 >> nxl;
   infile >> dummy >> dummy >> dummy >> L;
   infile2 >> nxl;
   infile >> dummy >> dummy >> dummy >> q;
   infile2 >> nxl;
   infile >> dummy >> dummy >> J0;
   infile2 >> nxl;
   infile >> dummy >> dummy >> JSTEPS;
   infile2 >> nxl;
   infile >> dummy >> dummy >> JINC;
   infile2 >> nxl;
   infile >> dummy >> dummy >> M;
   infile2 >> nxl;
   infile >> dummy >> dummy >> nmeas;
   infile2 >> nxl;
   infile >> dummy >> dummy >> nskip;
   infile2 >> nxl;
   infile >> dummy >> dummy >> nequi;
   infile2 >> nxl;
   infile >> nxl;
   infile2 >> nxl;


   cout << "length L ............... " << L << "\n";
   cout << "states q ............... " << q << "\n";
   cout << "Jstart   ............... " << J0 << "\n";
   cout << "Jsteps   ............... " << JSTEPS << "\n";
   cout << "Jincr    ............... " << JINC << "\n";
   cout << "M(ext.)  ............... " << M << "\n";
   cout << "nmeas    ............... " << nmeas << "\n";
   cout << "nskip    ............... " << nskip << "\n";
   cout << "nequi    ............... " << nequi << "\n";
   cout << "\n\n\n";
   cout << "Out pattern:\n";
   cout << "j\t [ Mr/V \tEpsilon Mr ] (for all states)\t[ Chir/(beta V) ] (for all states)\tC/V\tH/V\tEpsilon H\t\n\n";

   ofstream outfile("finaldata.dat");  //prepare output
   outfile.width(18);          
   outfile.precision(10);
   outfile.setf(ios_base::fixed);

   j = J0;                               //set initial value

   for(int i=0; i<JSTEPS; i++)
   {
      j += JINC;
      H = 0;  //reset H

      for(int k=0; k<q; k++)
      {
         mag[k]=0; //init mag array
         chi[k]=0; //init chi array
      }

      for(int ii=0; ii<nmeas; ii++) //observable <M>
      {
         for(int iii=0; iii<q; iii++)
         {
            infile >> aux;
            mag[iii] += aux;
            //cout << aux << " ";
            aux=0;
         }
         infile >> aux; //observable <H>
         H+=aux;
         aux=0;
      } //end for nmeasloop

      H = H/nmeas;   //mean energy
      h[i] = H/(L*L);//store mean of H

      for(int iv=0; iv<q; iv++)//mean magnetization
      {
         mag[iv] = mag[iv]/nmeas;
      }

      for(int ii=0; ii<nmeas; ii++)//observable Chi
      {
         for(int iii=0; iii<q; iii++)
         {
            infile2 >> aux;
            chi[iii] += (aux - mag[iii]) * (aux - mag[iii]);
            aux = 0;
         }

         infile2 >> aux; //observable C
         C += (aux - H) *  (aux - H);
         aux = 0;
      } //end for nmeas loop chi

      for(int i=0; i<q; i++)
      {
         mag[i]  = mag[i]/(L*L);
         chi[i]  = chi[i]/nmeas;
         epsm[i] = sqrt(chi[i]/(nmeas - 1)); //error of Mr
         epsm[i]/=(L*L);
         chi[i]  =chi[i]/(L*L);
      }

      C     = C/nmeas;
      epsH  = sqrt(C/(nmeas - 1)); //error of H
      epsH /= (L*L);
      C     = C/(L*L);
      cout    << j << "  "; //output to screen and file
      outfile << j << "  ";
      H=H/(L*L);

      for(int i=0; i<q; i++)//output magnetization and error
      {
         cout    << mag[i] << "  " << epsm[i] << "  ";
         outfile << mag[i] << "  " << epsm[i] << "  ";
      }

      for(int i=0; i<q; i++)//output chi
      {
         cout    << chi[i] << "  ";
         outfile << chi[i] << "  ";
      }

      cout    << C << "  ";//output C
      outfile << C << "  ";

      cout    << H << "  " << epsH << "\n";//output H
      outfile << H << "  " << epsH << "\n";
   } //end observables

   cout << endl;
   cout << "Output written to 'finaldata.dat'." << endl;
   cout << "Done. " << endl;
   return 0;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  end main  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

