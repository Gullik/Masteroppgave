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



    //input.c
    convert();                  //Reading input parameters
    readdata(argc, argv);       //nGridPoints ->  ngx, ngy, ngz

    //grid.c
    memorygrid();               //Clearing memory for ***rho

    //mg/fmg_P.c
    mglin_init(ngx, ngy, ngz);  //

    /***************Now I think I have what I need ********************/
    /* Time to assign a potential field and try to use the MG solver  */
    printf("# of grid points in each direction: %d \n ", ngx);

    int nPoints = ngx*ngy*ngz;  //Number of grid points
    int i;


    rho = dvecmem(0, nPoints - 1 );

    for(i = 0; i < nPoints ; i++)
        rho[i] = 1;

    /********* Let's try to us it **********/
    mglin(rho, NCYCLES);

    return(0);
}
