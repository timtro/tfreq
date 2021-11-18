// ______________________________________________________________
//
// Program: tfreq
// Author : Timothy A.V. Teatro <timothy.teatro@uoit.ca>
// Date : July 2008
//
// Description:
// This code reads in a set of atomic trajectories in .pos
// format and computes the velocity autocorrelation function
// (VAC). The VAC algorithm was adapted from the formulas
// given in J. Kohanoff Comp. Mat. Sci. 2 221-232 (1994).
// The estimator formulation used here is unbiased and
// statistically consistent, but may yield a non-valid VAC.
// Using the FFTW package, the Fourier transform is computed
// and the resultant frequency spectrum is printed to stdout.
//
// NB: Alpha stage: FFT data is not yet abscissiated properly.
//
// Usage:
// Simply compile:
// :> make all
// (optionally, substitute your compiler of choice). Then
// run from the command line:
//    ./tfreq case.vel > case.freq
//
// WARNINGS:
// There is no error checking, so don't forget to specify the
// vel file and verify its integrity.
//
// BEGIN_________________________________________________________
//
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//
// User configurable variables
//
int padd = 15;
//
// Define constants.
//
const double H_PLANCK_SI = 6.6260693E-34;        // J s
const double HARTREE_SI = 4.35974417E-18;        // J
const double K_BOLTZMANN_SI = 1.3806505E-23;     // J K^-1
const double BOHR_RADIUS_SI = 0.5291772108E-10;  // m
const double AMU_SI = 1.66053886E-27;            // Kg
const double AU_SEC = 2.41888432E-17;            // s
//
// Define global variables.
//
int ts = 0;             // Stores the current MD time step.
double mdTime = 0.00f;  // Stores the MD time in fs.
int nat = 1;            // Number of atoms.
FILE *velfile;          // The input stream.
double ***v;            // Matrix of atom velocities.
//
// Prototype function(s).
//
void getnat();  // Get the number of atoms.

int main(int argc, char **argv) {
  int i = 0;  // Dummy.
  int j = 0;  // Dummy.
  int m = 0;  // Dummy.
  int n = 0;  // Dummy.
  int M = 0;  // Total number of time steps.
  char file[256] = "";
  fftw_complex *Z;                  // Velocity autocorrelation vector.
  double a;                         // Junk/Temporary storage.
  char eof;                         // EOF state of velfile.
  int tsstart = 0;                  // Starting time-step.
  int tsstop = 0;                   // Starting time-step.
  double tstart = 0.0f;             //
  double tstop = 0.0f;              //
  double sqtp = 2.506628274631000;  // sqrt of 2 Pi
  fftw_complex norm = 1.0f;
  double sigma = 1.0f;
  fftw_complex *dft_out;
  fftw_plan dft_plan;
  FILE *zfile;
  FILE *ftfile;
  //
  // Open file
  //
  printf("\n");
  printf("#        TFREQ\n");
  printf("#        VDOS Calculator\n");
  printf("#        Version 0.5.1 beta\n");
  printf("#        Timothy A.V. Teatro <timothy.teatro@uoit.ca>\n");
  printf("#\n# Starting Program.\n");

  if ((velfile = fopen(argv[1], "r")) == NULL) {
    printf("ERROR: Cannot open file\n");
    exit(1);
  }
  strncpy(file, "Z.", 2);
  strncat(file, argv[1], 256);
  if ((zfile = fopen(file, "w")) == NULL) {
    printf("ERROR: Cannot open Z file\n");
    exit(1);
  }
  strncpy(file, " ", 256);
  strncpy(file, "FT.", 3);
  strncat(file, argv[1], 256);
  if ((ftfile = fopen(file, "w")) == NULL) {
    printf("ERROR: Cannot open FT file\n");
    exit(1);
  }
  printf("# All files succesffully opened\n");
  //
  // Count each time step in the file so we can dimension our variables.
  // Obviously, this isn't ideal. I'd like to dynamically dimension the
  // memory and store it in the same loop, but it's more trouble than it's
  // worth for now.
  //
  fscanf(velfile, "%d%lg", &tsstart, &tstart);
  rewind(velfile);
  getnat();
  rewind(velfile);
  while (1) {
    fscanf(velfile, "%d%lg", &tsstop, &tstop);
    if (feof(velfile)) break;
    for (i = 0; i < nat; i++) {
      fscanf(velfile, "%lG%lG%lG", &a, &a, &a);
    }
    M++;
  }
  rewind(velfile);
  printf("# Number of atoms: %d\n", nat);
  printf("# MD Start time is %f ps, and end time is %f ps.\n", tstart, tstop);
  printf("# MD Duration is %f ps\n", tstop - tstart);
  printf(
      "# Data starts at time step %d and ends at step %d, totalling %d "
      "steps.\n#\n",
      tsstart, tsstop, M);
  //
  // Dimension the matrices.
  //
  printf("# Allocating memory...\n");
  v = (double ***) malloc(M * sizeof(double **));
  for (i = 0; i < M; i++) {
    v[i] = (double **) malloc(nat * sizeof(double *));
    for (j = 0; j < nat; j++) v[i][j] = (double *) malloc(3 * sizeof(double));
  }
  Z = (fftw_complex *) fftw_malloc(M * padd * sizeof(fftw_complex));
  dft_out = (fftw_complex *) fftw_malloc(M * padd * sizeof(fftw_complex));
  //
  // Rewind and start populating our velocity matrix. This is a 3D matrix
  // which contains the velocities of each atom at each time step,
  // indexed as v[t][a][x] where t is the timestep (0..M-1), a is
  // the number of atoms (0..nat-1) and x is the coordinate specifier (0..2).
  //
  printf("#   done.\n# Reading %d velocity array into primary memory...\n", M);
  for (i = 0; i < M; i++) {
    fscanf(velfile, "%d%lg", &ts, &mdTime);
    for (j = 0; j < nat; j++) {
      fscanf(velfile, "%lG%lG%lG", &v[i][j][0], &v[i][j][1], &v[i][j][2]);
    }
  }
  printf("# Successfully read %d time steps of velocity data.\n#\n", M);
  fclose(velfile);
  // Zero the Z vector -
  for (m = 0; m < M; m++) Z[m] = 0.000000f;
  //
  // Now compute the VAC sequence.
  //
  // This algorithm was adapted from the formulas given in
  // J. Kohanoff Comp. Mat. Sci. 2 221-232 (1994). The estimator
  // formulation used here is unbiased and statistically consistent,
  // but may yield a non-valid VAC.
  //
  //
  // Looping through all time origins to collect an average -
  //
  printf("# Starting autocorrelation...\n");
  for (m = 0; m < M; m++)
    //
    // Looping through each calculation that is part of Z[m] -
    //
    for (n = 0; n < M - m - 1; n++)
      //
      // Looping to take the average over all atoms -
      //
      for (i = 0; i < nat; i++)
        //
        // Dot product -
        //
        for (j = 0; j < 3; j++) Z[m] += v[n + m][i][j] * v[n][i][j];
  //
  // VAC calculation done. Applying scaling.
  //
  printf("#   done.\n# Scaling Z[m]\n");
  for (m = 0; m < M; m++) {
    Z[m] /= (M - m);
  }
  for (m = 1; m < M; m++) {
    Z[m] /= Z[0];
  }
  Z[0] = 1.0f;
  //
  // Zero our padding and apply gaussian smoothing, then rescale.
  //
  printf("#   done.\n# Padding with zeros and gaussian smoothing...\n");
  for (i = M; i < M * padd; i++) Z[i] = 0.000f;
  sigma = M / 2.50;
  for (i = 0; i < M * padd; i++)
    Z[i] *= exp(-i * i / (2 * sigma * sigma)) / (sigma * sqtp);
  for (m = 1; m < M * padd; m++) {
    Z[m] /= Z[0];
  }
  Z[0] = 1.0f;
  //
  // DFT
  //
  dft_plan =
      fftw_plan_dft_1d(M * padd, Z, dft_out, FFTW_FORWARD, FFTW_ESTIMATE);
  printf("#   done.\n# Executing FFT\n");
  fftw_execute(dft_plan);
  //
  // Normalize the data (two methods)...
  //
  printf("#   done.\n# Normalizing data...\n");
  //
  // Normalize the peak so that < F | F* > = 1
  //
  // norm = 0.000f;
  // for ( m = 0; m < M*padd; m++ )
  //    norm += dft_out[m] * conj(dft_out[m]);
  // norm *= (tstop - tstart)*1E-12 / M;
  // norm = 1/norm;
  //
  // Normalize max peak to 1
  //
  norm = 0.00f;
  for (m = 0; m <= (M * padd) / 2; m++)
    if (creal(dft_out[m] * conj(dft_out[m])) > creal(norm))
      norm = dft_out[m] * conj(dft_out[m]);
  norm = 1 / norm;
  //
  // Print Output:
  //
  printf("#   done.\n# Writing output...\n");
  for (m = 0; m < M * padd; m++) {
    fprintf(zfile, "%14.9E %14.9f\n", tstart + (tstop - tstart) * m / M, Z[m]);
  }
  for (m = 0; m <= ((M * padd) / 2); m++) {
    fprintf(ftfile, "%17.9E %17.9E %17.9E %17.9E\n",
            m / (2.99792458E10 * (tstop - tstart) * padd * 1E-12),
            creal(dft_out[m]), cimag(dft_out[m]),
            norm * dft_out[m] * conj(dft_out[m]));
  }
  printf("#   done.\n# !! Program Complete !!\n");
}

void getnat() {
  char item[20];
  char eof;

  fscanf(velfile, "%d%lg", &ts, &mdTime);
  nat = 0;

  do {
    eof = fscanf(velfile, " %s", item);
    nat++;
  } while (strchr(item, '.') != NULL && eof != EOF);
  nat /= 3;
}
