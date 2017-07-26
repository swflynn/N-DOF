CXX = g++
FC = gfortran
CXXFLAGS = -O2 -Wall
FFLAGS = -O2 -Wall

LD = $(CXX)
LDFLAGS = $(CXXFLAGS)


all: scp_Lmon

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^


scp_Lmon: scp_Lmon.o TIP4P.o cg.o sobol.o sobol_stdnormal.o \
	mbpol/libmbpol.a
	$(FC) -o $@ $^ mbpol/libmbpol.a $(FFLAGS) -lstdc++


clean:
	rm -f libX2O.a *.o *.mod scp_lmon

