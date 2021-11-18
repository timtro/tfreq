//______________________________________________________________________________
//
// Program: pos2vel.
// Author : Timothy A.V. Teatro <timothy.teatro@uoit.ca>
// Date   : Jan 2009
//
// Description:
//    This program is intended to differentiate the positions in an MD
//    trajectory and create a vel file. Here, we use a 'centered'
//    differentiation scheme. I.E.
//       f'(t) = [ f(t+dt)-f(t-dt) ] / 2dt
// Usage:
//   Simply compile:
//     g++ -O3 -o thisfile thisfile.cpp
//   (optionally, substitute your compiler of choice).
//   then run from command line:
//     ./thisfile case.pos > case.vel
//
// WARNINGS:
//   There is no error checking, so don't forget to specify the pos file and to
//   verify its integrity.
//
//BEGIN_________________________________________________________________________
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//
// Define constants.
//
const int       nat=96;            // Number of Atoms in System

const double    H_PLANCK_SI      = 6.6260693E-34;   // J s
const double    HARTREE_SI       = 4.35974417E-18;  // J
const double    K_BOLTZMANN_SI   = 1.3806505E-23;   // J K^-1
const double    BOHR_RADIUS_SI   = 0.5291772108E-10;// m
const double    AMU_SI           = 1.66053886E-27;  // Kg
const double    AU_SEC           = 2.41888432E-17;  // s

//
// User defined variables
//

double          dt               = 15*2*2.41888432E-17;  // A.U.
int           pwcp               = 1;   // Factor to differentiate PW from CP time step units.

//
// Define global variables.
//
int             ts=0;              // Stores the current MD time step.
double          mdTimeA=0.00;       // Stores the MD time in fs.
double          mdTimeB=0.00;       // Stores the MD time in fs.
double          mdTimeC=0.00;       // Stores the MD time in fs.
FILE           *posfile;           // The input stream.
double        **atomPosA;          // Matrix of atom velocities, t-dt.
double        **atomPosB;          // Matrix of atom velocities, t.
double        **atomPosC;          // Matrix of atom velocities, t+dt.
double        **atomVel;           // Atomic velocity vector at t.
//
// Prototype function(s).
//
inline int      getPos ();         // Function to read atom velocities from
                                   //   posfile.
inline int      copyto (double **A, double **B); // Copy A onto B

main (int argc, char **argv)
{
   int             i=0,              // Dummy.
                   j=0,              // Dummy.
                   TotalTimeSteps=0; // Count of timesteps read.
   //
   // Start by opening the file and dimensioning the velocity matrix.
   //
   if ( (posfile = fopen ( argv[1], "r" )) == NULL )
   {
      printf("ERROR: Cannot open file\n");
      exit(1);
   }

   atomPosA = (double**)malloc( nat*sizeof(double*) );
   atomPosB = (double**)malloc( nat*sizeof(double*) );
   atomPosC = (double**)malloc( nat*sizeof(double*) );
   atomVel = (double**)malloc( nat*sizeof(double*) );
   for(i=0; i<nat; i++) atomPosA[i] = (double*)malloc( 3*sizeof(double) );
   for(i=0; i<nat; i++) atomPosB[i] = (double*)malloc( 3*sizeof(double) );
   for(i=0; i<nat; i++) atomPosC[i] = (double*)malloc( 3*sizeof(double) );
   for(i=0; i<nat; i++) atomVel[i] = (double*)malloc( 3*sizeof(double) );

   //
   // We'll set out our initioal atomPosA/B/C, and after this, we only need to
   // read one position at a time and just shuffle them around in memory. Note,
   // if I get more time, look into doing this in a more cleaver way with
   // pointers to avoid the excessive memory operations.
   //

   getPos();
   copyto( atomPosC, atomPosA );
   mdTimeA = mdTimeC;
   getPos();
   copyto( atomPosC, atomPosB );
   mdTimeB = mdTimeC;
   TotalTimeSteps=2;
   while ( getPos() != EOF )
   {
      TotalTimeSteps++;

//      printf("ts: %6d     A: %13.9f,    B: %13.9f,      C: %13.9f\n", TotalTimeSteps, atomPosA[0][0], atomPosB[0][0], atomPosC[0][0]);
      printf("%7d %14.5f\n", ts-1, mdTimeB);
      for ( i = 0; i < nat; i++ )
      {
         for ( j = 0; j < 3; j++ )
         {
            atomVel[i][j] = ( atomPosC[i][j]-atomPosA[i][j] ) * BOHR_RADIUS_SI / (2*dt);
         }
         printf("%18.9E%18.9E%18.9E\n", atomVel[i][0], atomVel[i][1], atomVel[i][2]);
      }
      // Move around the data for a fresh entry to atomPosC at the start of the
      // loop.
      copyto( atomPosB, atomPosA );
      mdTimeA = mdTimeB;
      copyto( atomPosC, atomPosB );
      mdTimeB = mdTimeC;
   }
   fclose(posfile);
}
//
// getPos() is a simple function to read the velocity coordinates out of the
// specified input file, which contains velocities in atomic units (bohr/AU).
//
inline int
getPos ()
{
   int             i=0, a=0;
   //
   // Read out the time-step and MD time in the file, then enter a loop which 
   // will read off the coordinates of the 'nat' atoms and performs conversions.
   //
   a=fscanf(posfile, "%d%lg", &ts, &mdTimeC);
   //
   // Now read the velocities of each atom.
   //
   for (; i<nat; i++)
   {
      fscanf(posfile, "%lG%lG%lG", &atomPosC[i][0],
                                   &atomPosC[i][1],
                                   &atomPosC[i][2]);
   }
   return a;
}
inline int
copyto (double **A, double**B)
{
   int i,j;
   for ( i = 0; i < nat; i++)
   for ( j = 0; j < 3; j++)
   {
      B[i][j] = A[i][j];
   }
   return 0;
}
