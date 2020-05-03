SRCS=$(wildcard *.f90)
OBJS=$(SRCS:.f90=.o)

FC_PC=mpif90

PKG_CONFIG_PATH=/usr/local/petsc/lib/pkgconfig



LD_PC_FLAGS=-O5 -fopenmp -lgomp /opt/metis-4.0.3/libmetis.a
FC_PC_FLAGS=-O5 -g -c -cpp -ffree-form -ffree-line-length-none -fopenmp -fbounds-check -DPETSC_AVOID_MPIF_H 



PETSC_LD_FLAGS=$(shell PKG_CONFIG_PATH=/usr/local/petsc/lib/pkgconfig pkg-config --libs PETSc)
PETSC_FC_FLAGS=$(shell PKG_CONFIG_PATH=/usr/local/petsc/lib/pkgconfig pkg-config --cflags PETSc)

EXEC=SPEED


.PHONY: all
all: $(OBJS) $(EXEC)


$(EXEC): $(OBJS)
	$(FC_PC) -o $@  $(OBJS) $(LD_PC_FLAGS) 

#$(PETSC_LD_FLAGS)


$(OBJS):%.o: %.f90 MODULES.o  
	$(FC_PC) $(FC_PC_FLAGS)  $< -o $@

#$(PETSC_FC_FLAGS)

.PHONY: clean
clean:
	-rm -f *.o *.mod SPEED
