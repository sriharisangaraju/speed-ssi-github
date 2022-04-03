SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

FC_PC=mpif90


FC_PC_FLAGS=-O5 -g -c -cpp -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -fallow-argument-mismatch -DPETSC_AVOID_MPIF_H
LD_PC_FLAGS=-O5  -fopenmp /Users/srihari/apps/metis-4.0.3-clang/libmetis.a 

#/opt/homebrew/Cellar/metis/5.1.0/lib/libmetis.dylib 
#/Users/srihari/apps/metis-4.0.3/libmetis.a

PETSC_LD_FLAGS=
PETSC_FC_FLAGS=


#PETSC_LD_FLAGS=-L/cineca/prod/opt/libraries/petsc/3.7/openmpi--1-10.3--gnu--6.1.0/lib -lpetsc
#PETSC_FC_FLAGS=-I/cineca/prod/opt/libraries/petsc/3.7/openmpi--1-10.3--gnu--6.1.0/include 



EXEC=SPEED


.PHONY: all
all: $(OBJS) $(EXEC)

$(EXEC): $(OBJS)
	$(FC_PC) -o $@  $(OBJS) $(LD_PC_FLAGS) $(PETSC_LD_FLAGS)


$(OBJS):%.o: %.f90  MODULES.o kdtree2.o 
	$(FC_PC) $(FC_PC_FLAGS) $(PETSC_FC_FLAGS) $< -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod SPEED
