###############################################################################
# Sample Makefile for use of FFTW library
# Jemmy Hu <jemmyhu@sharcnet.ca>
# March, 2006
#
#compile MPI c++ code using FFTW-2.1.5 library and Pathscale compiler on Narwhal
#double precision is used
###############################################################################
#Compiler
CC = gcc

SRC = tfreq.c
PSRC = pos2vel.c

OBJ=$(SRC:.c=.o)
POBJ=$(PSRC:.c=.o)

#compiler flags
#MPIDIR1=/opt/hpmpi
#MPIDIR=/opt/hpmpi/lib/linux_amd64
FFTWDIR=/opt/sharcnet/fftw/3.1
#ACMLDIR=/opt/sharcnet/acml/acml-pathscale-64bit-current/pathscale64

IDIRS  =  -I $(FFTWDIR)/include
#IDIRS2  = -I $(MPIDIR1)/include
FLAGS  =  $(IDIRS)
#LAPACK  = -L $(ACMLDIR)/lib -lacml -lpathfortran
#MPI    = -L $(MPIDIR) -lmpi
LIBS   = -L $(FFTWDIR)/lib -lfftw3 -lm

#compilation 
tfreq: $(OBJ)
	$(CC) -O3 $(IDIRS) $(OBJ) $(LIBS) -o $@

pos2vel: $(POBJ)
	$(CC) -O3 $(IDIRS) $(OBJ) $(LIBS) -o $@

%.o: %.c
	$(CC) $(FLAGS) -c $<

clean : 
	\rm -f *.o *~ *.mod
