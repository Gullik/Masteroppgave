/* DiP3D */
/* Main function */
/* Author: Wojciech Jacek Miloch */
/* University of Oslo, Norway */
/* 2009 */

/* Simplified version to work only
 * on the MG part
 *
   Gullik Vetvik Killie, 2015*/

#include "const.h"


int main(int argc, char *argv[]) {

    /*Initialize MPI*/
    int numtasks=1;
    int rank=0;
#ifdef MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks); //no of processesors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //numbers of processes
#endif

    if(rank==0)
        printf("DiP3D, Running only MG-part \n");
    if(rank==0)
        convert();


//    if(rank == 0)
//        mglin_init();


    return(0);
}
