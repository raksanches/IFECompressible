static char help[] = "Code for solving compressible flows";

#include "src-fluid/FluidDomain.h"

int main(int argc, char **args)
{
	// Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, help);
    int rank, size;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

	//#include "test.h"
    //#include "Example1.h"
    //#include "Example2.h"
    //#include "Example3.h"
    //#include "Example4.h"
    #include "Example5.h"    
     
	PetscFinalize();
	return 0;
}
