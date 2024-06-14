MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 		 = $(wildcard *.cpp src-fluid/*.cpp  src-fluid/mesh_interface/*.cpp)
FCOMPILER        = gfortran -O2
CXXFLAGS        += -w -DANSI_DECLARATORS

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

f: $(CSOURCES:.cpp=.o) 
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -lboost_system -std=c++11

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

clearall:
	@$ rm *.o src-fluid/*.o  *f results/*.vtu

clearinput:
	@$ rm Main.o results/*.vtu *f

run1:
	./f

run2:
	@$ mpirun -np 2 ./f

run4:
	@$ mpirun -np 4 ./f

run8:
	@$ mpirun -np 8 ./f


