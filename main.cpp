#include "compute.h"

int main(int argc, char* argv[]) {
    // Initialize QUESO environment

#ifndef QUESO_HAS_MPI
    MPI_Init(&argc,&argv);
    FullEnvironment* env = new FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
    FullEnvironment* env = new FullEnvironment(argv[1],"",NULL);
#endif

    //call apllication
    compute(*env);

    //Finalize QUESO environment
    delete env;
#ifndef QUESO_HAS_MPI
    MPI_Finalize();
#endif
    
    return 0;
} 
