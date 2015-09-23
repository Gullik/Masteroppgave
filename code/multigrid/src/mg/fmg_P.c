/* Based on Numerical Recipes in C. */
/* To solve eliptic PDE on 2D grid with the Full Multi Grid method*/
/* If not possible to use fmm, then use Gauss-Seidel method*/

#include <stdio.h>
#include<math.h>
#define NRANSI
#include "nrutil.h"
#include "../const.h"
//#include "funct.h"
#define NPRE 4
#define NPOST 4
//#define NGMAX 15
inline int ix(int off, int i, int j, int k);

// we have external
//int fmg_ng, fmg_nnx, fmg_nny, fmg_mingridx, fmg_mingridy;
//double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];

/*3D*/
void mglin_init(int nx, int ny, int nz)
{
    int nn,ng,ng1,nn1,nnx,nny,nnz,ngrid,nn1x,nn1y,nn1z, k,j;
    printf("using Poisson solver with periodic boundary conditions\n");

    ng=0;
    fmg_nnx=nx;
    fmg_nny=ny;
    fmg_nnz=nz;

    if(nx==ny)
    {
        if(nx==nz)
        {
            nn=nx;
            fmg_mingridx=fmg_mingridy=fmg_mingridz=3;
            while (nn >>= 1) ng++;
            fmg_ng=ng;
            printf("start ng %d nn %d comp %d\n", ng,nn, 1+(1L<<ng)); //getchar();
            if (nx != 1+(1L << ng))
            {
                printf("WARNING: cells number is not a power of 2 in mglin. Will use a non-optimal field solver!\n");
                //find minimal grid
                nn1=nx;
                ng1=1;
                while((nn1 % 2) == 1)
                {
                    ng1++; nn1=1+(nn1-1)/2;
                }
                fmg_mingridx=fmg_mingridy=fmg_mingridz=nn1;
                fmg_ng=ng1;
            }
        }
    }
    else
    {
        printf("WARNING: cells number is not a power of 2 in mglin. Will use a non-optimal field solver!\n");
        printf("Tutaj algorytm chyba nie dziala i zaraz padnie!!\n");
        nn1=nx;
        ng1=1;
        while((nn1 % 2) == 1)
        {
            ng1++; nn1=1+(nn1-1)/2;
        }
        ng=ng1;
        nn1=ny;
        ng1=1;
        while((nn1 % 2) == 1)
        {
            ng1++; nn1=1+(nn1-1)/2;
        }
        if(ng1 < ng)
            ng=ng1;

        ng1=1;
        nn1=nz;
        while((nn1 % 2) == 1)
        {
            ng1++; nn1=1+(nn1-1)/2;
        }
        if(ng1 < ng)
            ng=ng1;


        nn1x=nx;
        nn1y=ny;
        nn1z=nz;

        for(k=0; k<ng-1; k++)
        {
            nn1x=1+(nn1x-1)/2;
            nn1y=1+(nn1y-1)/2;
            nn1z=1+(nn1z-1)/2;
        }
        //account for ghost points
        fmg_mingridx=nn1x;
        fmg_mingridy=nn1y;
        fmg_mingridz=nn1z;

        fmg_nnx=nx;
        fmg_nny=ny;
        fmg_nnz=nz;

        nnx=nx;
        nny=ny;
        nnz=nz;

        fmg_ng=ng;
    }
    printf("mglin init (excluding ghost points): ng %d mingrid %d and %d and %d  nnx %d %d %d\n", fmg_ng, fmg_mingridx, fmg_mingridy, fmg_mingridz, fmg_nnx, fmg_nny, fmg_nnz);
    if (fmg_ng > NGMAX) nrerror("increase NGMAX in mglin.");

    // allocate memory for grid arrays
    //periodic ghost points are included
    if(fmg_ng>1)
    {
        nnx=(fmg_nnx/2+1);
        nny=(fmg_nny/2+1);
        nnz=(fmg_nnz/2+1);
        ngrid=fmg_ng-1;
        irho[ngrid]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
        fill0(irho[ngrid],  nnx,  nny, nnz);

        printf("rho ngrid %d : %E \n", ngrid, irho[ngrid][1][1][1]);
        while (nnx > fmg_mingridx) {
            nnx=nnx/2+1;
            nny=nny/2+1;
            nnz=nnz/2+1;
            irho[--ngrid]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            fill0(irho[ngrid],  nnx,  nny, nnz);

        }

        nnx=fmg_mingridx;
        nny=fmg_mingridy;
        nnz=fmg_mingridz;
        iu[1]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
        irhs[1]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
        fill0(iu[1],  nnx,  nny, nnz);
        fill0(irhs[1],  nnx,  nny, nnz);

        printf("rhs hmm %E \n", irhs[1][0][0][0]);

        for (j=2;j<=fmg_ng;j++) {
            nnx=(2*nnx-1);
            nny=(2*nny-1);
            nnz=(2*nnz-1);
            //	printf("%d %d %d for %d\n", nnx, nny, nnz, j);
            //	getchar();
            iu[j]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            irhs[j]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            ires[j]=d3tensor(0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);

            fill0(iu[j],  nnx,  nny, nnz);
            fill0(irhs[j],  nnx,  nny, nnz);
            fill0(ires[j],  nnx,  nny, nnz);

            //	printf("rhs %d \n", irhs[j][0][0][0]);
        }
    }

}

/*3D*/
void mglin_destroy()
{
    int nnx,nny,nnz,j;
    if(fmg_ng>1)
    {
        for(nnx=fmg_nnx,nny=fmg_nny,nnz=fmg_nnz,j=fmg_ng;j>=2;j--,nnx=(nnx/2+1),nny=(nny/2+1),nnz=(nnz/2+1))
        {
            free_d3tensor(ires[j],0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            free_d3tensor(irhs[j],0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            free_d3tensor(iu[j],0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
            if (j != fmg_ng) free_d3tensor(irho[j],0,(nnx-1)+2,0,(nny-1)+2,0,(nnz-1)+2);
        }
        free_d3tensor(irhs[1],0,fmg_mingridx-1+2,0,fmg_mingridy-1+2,0,fmg_mingridz-1+2);
        free_d3tensor(iu[1],0,fmg_mingridx-1+2,0,fmg_mingridy-1+2,0,fmg_mingridz-1+2);
        free_d3tensor(irho[1],0,fmg_mingridx-1+2,0,fmg_mingridy-1+2,0,fmg_mingridz-1+2);
    }
}


void mglin(double *u, int ncycle)
{
    void addint(double ***uf, double ***uc, double ***res, int nfx, int nfy, int nfz);
    void copy(double ***aout, double ***ain, int nx, int ny, int nz);
    void copy0(double ***aout, double *ain, int nx, int ny, int nz);
    void copyfinal(double *aout, double ***ain, int nx, int ny, int nz);
    void fill0(double ***u, int nx, int ny, int nz);
    void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz);
    void relax(double ***u, double ***rhs, int nx, int ny, int nz);
    void resid(double ***res, double ***u, double ***rhs, int nx, int ny, int nz);
    void rstrct(double ***uc, double ***uf, int ncx, int ncy, int ncz);
    void rstrct0(double ***uc, double *uf, int ncx, int ncy, int ncz);
    void slvsml(double ***u, double ***rhs);
    void slvsml2(double ***u, double ***rhs, int nx, int ny, int nz);

    unsigned int j,jcycle,jj,jpost,jpre,ng=0,ngrid;
    unsigned int  nnx, nny, nnz, nfx,nfy,nfz;

    printf("Solving field equations in mglin\n");

    /*initialize fmg*/
    //minus ghost points
    ng=fmg_ng;
    nnx=fmg_nnx;
    nny=fmg_nny;
    nnz=fmg_nnz;

    if(ng>1)
    {
        nnx=nnx/2+1;
        nny=nny/2+1;
        nnz=nnz/2+1;
        ngrid=ng-1;
        printf("parameters ngrid (excluding ghost points) %d %d %d %d\n", ngrid, nnx, nny, nnz);
        //	printf("parameters irho ngrid %E\n", irho[6][1][1][1]);
        rstrct0(irho[ngrid],u,nnx,nny,nnz);

        int zi,zj,zk;
        printf("original ngx %d ngy %d ngz %d\n", ngx, ngy, ngz);
        for(zi=0; zi<ngx; zi++){
            for(zj=0; zj<ngy; zj++)
            {
                for(zk=0; zk<ngz; zk++)
                    ;//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //printf("\n");
            }//printf("\n");
        }
        printf("it was original ngx %d ngy %d ngz %d\n", ngx, ngy, ngz);
        //getchar();
        printf("rescrtrcted0 ngx %d ngy %d ngz %d\n", nnx, nny, nnz);
        for(zi=1; zi<nnx+1; zi++){
            for(zj=1; zj<nny+1; zj++)
            {
                for(zk=1; zk<nnz+1; zk++)
                    ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //printf("%E\t", irho[ngrid][zi][zj][zk]);
                //printf("\n");
            }//printf("\n");
        }
        printf("it was rescricted0 ngx %d ngy %d ngz %d\n", nnx, nny, nnz);
        //getchar();


        /*TUTAJ*/
        while (nnx > fmg_mingridx) {
            nnx=nnx/2+1;
            nny=nny/2+1;
            nnz=nnz/2+1;
            --ngrid;
            rstrct(irho[ngrid],irho[ngrid+1],nnx,nny,nnz);

            {printf("restricted ngx %d ngy %d ngz %d\n", nnx, nny, nnz);
                for(zi=1; zi<nnx+1; zi++){
                    for(zj=1; zj<nny+1; zj++)
                    {
                        for(zk=1; zk<nnz+1; zk++)
                            ; 	  //printf("%E\t", u[ix(0,zi,zj,zk)]);
                        //  printf("%E\t", irho[ngrid][zi][zj][zk]);
                        //printf("\n");
                    }
                    //printf("\n");
                }
            }
        }
        nnx=fmg_mingridx;
        nny=fmg_mingridy;
        nnz=fmg_mingridz;
        printf("got here 1\n"); //getchar();
        if((fmg_mingridx==3) && (fmg_mingridy==3) && (fmg_mingridz==3))
        {slvsml(iu[1],irho[1]);  printf("I am using the full multigrid method\n"); //getchar();
        }
        else
        {
            printf("I am using the semi - multigrid method (not the best choice of grid size)\n");
            slvsml2(iu[1],irho[1],fmg_mingridx,fmg_mingridy,fmg_mingridz);
        }
        ngrid=ng;

        {printf("solved initial ngx %d ngy %d ngz %d\n", nnx, nny, nnz);
            printf("check here periodic \n");
            //getchar();
            for(zi=0; zi<nnx+2; zi++)
            {
                for(zj=0; zj<nny+2; zj++)
                {
                    for(zk=0; zk<nnz+2; zk++)
                        ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                    //printf("%E\t", iu[1][zi][zj][zk]);
                    //printf("\n");
                }
                //printf("\n");
            }
            //getchar();
        }
        for (j=2;j<=ngrid;j++) {

            printf("OBS check here before interp periodic 2 %d \n", nnx);
            //getchar();
            for(zi=0; zi<nnx+2; zi++)
            {
                for(zj=0; zj<nny+2; zj++)
                {
                    for(zk=0; zk<nnz+2; zk++)
                        ;  //printf("%E\t", u[ix(0,zi,zj,zk)]);
                    //  printf("%E\t", iu[j-1][zi][zj][zk]);
                    //printf("\n");
                }
                //printf("\n");
            }
            // getchar();

            nnx=2*nnx-1;
            nny=2*nny-1;
            nnz=2*nnz-1;
            printf("now interp?\n");
            interp(iu[j],iu[j-1],nnx,nny,nnz);

            printf("OBS check here periodic 2 %d \n", nnx);
            //  getchar();
            for(zi=0; zi<nnx+2; zi++)
            {
                for(zj=0; zj<nny+2; zj++)
                {
                    for(zk=0; zk<nnz+2; zk++)
                        ;	  //printf("%E\t", u[ix(0,zi,zj,zk)]);
                    //	  printf("%E\t", iu[j][zi][zj][zk]);
                    //printf("\n");
                }
                //			  printf("\n");
            }
            //		  getchar();


            //	printf("rhs %d \n", irhs[2][0][0][0]); getchar();
            printf("I got here j: %d and ng %d\n",j,ng);
            if(j!=ngrid)
            {
                //  printf("try to copy %E\n", irhs[j][0][0][0]);
                copy(irhs[j], irho[j], nnx,nny,nnz);
                //via periodic

            }
            else
            {
                // printf("j=%d\n",j);
                copy0(irhs[j], u,nnx,nny,nnz);

            }
            for (jcycle=1;jcycle<=ncycle;jcycle++) {
                nfx=nnx;
                nfy=nny;
                nfz=nnz;
                for (jj=j;jj>=2;jj--) {
                    for (jpre=1;jpre<=NPRE;jpre++)
                        relax(iu[jj],irhs[jj],nfx,nfy,nfz);


                    resid(ires[jj],iu[jj],irhs[jj],nfx,nfy,nfz);
                    nfx=nfx/2+1;
                    nfy=nfy/2+1;
                    nfz=nfz/2+1;
                    rstrct(irhs[jj-1],ires[jj],nfx,nfy,nfz);
                    fill0(iu[jj-1],nfx,nfy,nfz);
                }
                if((fmg_mingridx==3) && (fmg_mingridy==3) && (fmg_mingridz==3))
                    slvsml(iu[1],irhs[1]);
                else
                    slvsml2(iu[1],irhs[1],fmg_mingridx,fmg_mingridy,fmg_mingridz);

                {printf("solved initial 2 ngx %d ngy %d ngz %d\n", nnx, nny, nnz);
                    for(zi=0; zi<3+2; zi++)
                    {
                        for(zj=0; zj<3+2; zj++)
                        {
                            for(zk=0; zk<3+2; zk++)
                                ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                            //printf("%E\t", iu[1][zi][zj][zk]);
                            //printf("\n");
                        }
                        //printf("\n");
                    }
                }
                nfx=fmg_mingridx;
                nfy=fmg_mingridy;
                nfz=fmg_mingridz;
                for (jj=2;jj<=j;jj++) {
                    nfx=2*nfx-1;
                    nfy=2*nfy-1;
                    nfz=2*nfz-1;
                    addint(iu[jj],iu[jj-1],ires[jj],nfx,nfy,nfz);

                    for (jpost=1;jpost<=NPOST;jpost++)
                        relax(iu[jj],irhs[jj],nfx,nfy,nfz);
                }
            }
        }
        printf("finished the fmg solver and now copying final\n");
        copyfinal(phi,iu[ngrid],fmg_nnx,fmg_nny,fmg_nnz);
    }
    else
    {
        int i,ipass, isw,j,jsw, sweep, k;
        int iter=0;
        double h2,hx,hy,hz,error,errorcheck,toler, factor;
        double ***testphi;
        jsw=1;
        toler=0.00000001;
        printf("Wrong grid size to use the multigrid method\nI am solving the potential with the slow converging Gauss-Seidel method\n");

        testphi=d3tensor(0,(fmg_nnx-1)+2,0,(fmg_nny-1)+2,0,(fmg_nnz-1)+2);
        {
            int ii,jj,kk;
            for (ii=0;ii<fmg_nnx;ii++)
                for (jj=0;jj<fmg_nny;jj++)
                    for(kk=0;kk<fmg_nnz;kk++)
                        testphi[ii+1][jj+1][kk+1]=0.0;
            //u[ix(0,ii,jj,kk)];
        }
        //copyfinal(rho,testphi,fmg_nnx,fmg_nny,fmg_nnz);

        fill0(testphi, fmg_nnx, fmg_nny, fmg_nnz);
        hx=Lx/(fmg_nnx-1);
        hy=Ly/(fmg_nny-1);
        hz=Lz/(fmg_nnz-1);

        h2=hx*hy;

        factor=1.0/6;
        do
        {


            //periodic boundaries
            //bottom

            for(i=0; i<fmg_nny+1; i++)
                for(j=0; j<fmg_nny+1; j++) {
                    testphi[0][i][j]=testphi[fmg_nnx-1][i][j];
                    testphi[fmg_nnx+1][i][j]=testphi[2][i][j];
                }
            for(i=0; i<fmg_nnx+1; i++)
                for(j=0; j<fmg_nny+1; j++) {
                    testphi[i][j][0]=testphi[i][j][fmg_nnz-1];
                    testphi[i][j][fmg_nnz+1]=testphi[i][j][2];
                }
            for(i=0; i<fmg_nnx+1; i++)
                for(j=0; j<fmg_nnz+1; j++) {
                    testphi[i][0][j]=testphi[i][fmg_nny-1][j];
                    testphi[i][fmg_nny+1][j]=testphi[i][2][j];
                }

            iter++;
            errorcheck=0.5*toler;
            for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
                for(i=0; i<fmg_nny+1; i++)
                    for(j=0; j<fmg_nny+1; j++) {
                        testphi[0][i][j]=testphi[fmg_nnx-1][i][j];
                        testphi[fmg_nnx+1][i][j]=testphi[2][i][j];
                    }
                for(i=0; i<fmg_nnx+1; i++)
                    for(j=0; j<fmg_nny+1; j++) {
                        testphi[i][j][0]=testphi[i][j][fmg_nnz-1];
                        testphi[i][j][fmg_nnz+1]=testphi[i][j][2];
                    }
                for(i=0; i<fmg_nnx+1; i++)
                    for(j=0; j<fmg_nnz+1; j++) {
                        testphi[i][0][j]=testphi[i][fmg_nny-1][j];
                        testphi[i][fmg_nny+1][j]=testphi[i][2][j];
                    }





                sweep=jsw;
                for (i=1;i<fmg_nnx+1;i++)
                    for (j=1;j<fmg_nny+1;j+=1)
                        for(k=(2-(i+j+sweep)%2); k<fmg_nnz+1; k+=2)
                        {
                            error=testphi[i][j][k];
                            testphi[i][j][k]=factor*(testphi[i+1][j][k]+testphi[i-1][j][k]+testphi[i][j+1][k]+testphi[i][j-1][k]+testphi[i][j][k-1]+testphi[i][j][k+1]+h2*u[ix(0,i-1,j-1,k-1)]);
                            //we have plus sign here, because the rho is positive here
                            error=testphi[i][j][k]-error;
                            if(errorcheck < fabs(error))
                                errorcheck=fabs(error);
                        }

                //  printf("%E\n", testphi[2][2][2]);
                //	  getchar();

                //this is only temporary to avoid huge numbers that can arise due to loops!
                double shift=testphi[0][0][0];
                for (i=0;i<fmg_nnx+2;i++)
                    for (j=0;j<fmg_nny+2;j++)
                        for(k=0; k<fmg_nnz+2; k++)
                            testphi[i][j][k]-=shift;



            }
            if(iter % 100 == 0) printf("iteration %d errorcheck= %E\n", iter, errorcheck);

        }
        //  while(errorcheck>toler);
        while(iter < 1000);
        printf("END %d iterations\n", iter);
        copyfinal(phi,testphi,fmg_nnx,fmg_nny,fmg_nnz);
        free_d3tensor(testphi,0,(fmg_nnx-1)+2,0,(fmg_nny-1)+2,0,(fmg_nnz-1)+2);
    }
}
#undef NPRE
#undef NPOST
#undef NGMAX
#undef NRANSI

/*3D*/
//periodic = dirichlet
void rstrct(double ***uc, double ***uf, int ncx, int ncy, int ncz)
{
    int ic,iif,jc,jf,kc,kf,nccx, nccy, nccz;
    /*set ncc to be the dimension of the larger grid*/
    nccx = 2*ncx-2+1; /*previously*/
    nccy = 2*ncy-2+1; /*previously*/
    nccz = 2*ncz-2+1;
    double wf=1.0/12;
    //printf("in rstrct ncx %d\n", ncx); //getchar();

    int i,j;
    for(i=0; i<nccy+1; i++)
        for(j=0; j<nccz+1; j++) {
            uf[0][i][j]=uf[nccx-1][i][j];
            uf[nccx+1][i][j]=uf[2][i][j];
        }
    for(i=0; i<nccx+1; i++)
        for(j=0; j<nccy+1; j++) {
            uf[i][j][0]=uf[i][j][nccz-1];
            uf[i][j][nccz+1]=uf[i][j][2];
        }
    for(i=0; i<nccx+1; i++)
        for(j=0; j<nccz+1; j++) {
            uf[i][0][j]=uf[i][nccy-1][j];
            uf[i][nccy+1][j]=uf[i][2][j];
        }

    {
        //printf("!!! UF before in normal restriction weighted ngx %d ngy %d ngz %d\n", ncx, ncy, ncz);
        int ii,jj,kk;
        for(ii=1; ii<nccx+1; ii++)
        {
            for(jj=1; jj<nccy+1; jj++)
            {
                for(kk=1; kk<nccz+1; kk++)
                    ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //			printf("%E\t", uf[ii][jj][kk]);
                //		printf("\n");
            }
            //	printf("\n");
        }
    }//	getchar();

    for(kf=1,kc=1;kc<ncz+1;kc++,kf+=2)
        for (jf=1,jc=1;jc<ncy+1;jc++,jf+=2){
            for (iif=1,ic=1;ic<ncx+1;ic++,iif+=2){
                uc[ic][jc][kc]=0.5*uf[iif][jf][kf]+wf*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]+uf[iif][jf+1][kf]+uf[iif][jf-1][kf]+uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
                //printf("table %E\t%E\t%E\t%E\t%E\t%E\t%E\n", uf[iif][jf][kf], uf[iif+1][jf][kf], uf[iif-1][jf][kf], uf[iif][jf+1][kf], uf[iif][jf-1][kf], uf[iif][jf][kf+1], uf[iif][jf][kf-1]);
            }
        }


    {
        //printf("!!! UF after in normal restriction weighted ngx %d ngy %d ngz %d\n", ncx, ncy, ncz);
        int ii,jj,kk;
        for(ii=0; ii<ncx+2; ii++)
        {
            for(jj=0; jj<ncy+2; jj++)
            {
                for(kk=0; kk<ncz+2; kk++)
                    ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //printf("%E\t", uc[ii][jj][kk]);
                //printf("\n");
            }
            //printf("\n");
        }
    }	//getchar();

    /*boundary points*/
    /*planes*/
    /*do bottom and top plane*/
    /*
  for (jf=1,jc=1;jc<ncy+1;jc++,jf+=2){
    for (iif=1,ic=1;ic<ncx+1;ic++,iif+=2){
      uc[ic][jc][1]=uf[iif][jf][1];
      uc[ic][jc][ncz]=uf[iif][jf][nccz];
    }
  }
  //do face and back (depth -> y)
  for (kf=1,kc=1;kc<ncz+1;kc++,kf+=2){
    for (iif=1,ic=1;ic<ncx+1;ic++,iif+=2){
      uc[ic][1][kc]=uf[iif][1][kf];
      uc[ic][ncz][kc]=uf[iif][nccy][kf];
    }
  }
  //do left and right (depth -> y)
  for (kf=1,kc=1;kc<ncz+1;kc++,kf+=2){
    for (jf=1,jc=1;jc<ncy+1;jc++,jf+=2){
      uc[1][jc][kc]=uf[1][jf][kf];
      uc[ncx][jc][kc]=uf[nccx][jf][kf];
    }
  }

  //do edges
  //bottom and top along x
  for (jc=1,ic=1;ic<ncx+1;ic++,jc+=2) {
    uc[ic][1][1]=uf[jc][1][1];
    uc[ic][ncy][1]=uf[jc][nccy][1];
    uc[ic][1][ncz]=uf[jc][1][nccz];
    uc[ic][ncy][ncz]=uf[jc][nccy][nccz];
  }

  //bottom and top along y
  for (jc=1,ic=1;ic<ncy+1;ic++,jc+=2) {
    uc[1][ic][1]=uf[1][jc][1];
    uc[ncx][ic][1]=uf[nccx][jc][1];
    uc[1][ic][ncz]=uf[1][jc][nccz];
    uc[ncx][ic][ncz]=uf[nccx][jc][nccz];
  }

  //vertically
  for (jc=1,ic=1;ic<ncz+1;ic++,jc+=2) {
    uc[1][1][ic]=uf[1][1][jc];
    uc[ncx][1][ic]=uf[nccx][1][jc];
    uc[1][ncy][ic]=uf[1][nccy][jc];
    uc[ncx][ncy][ic]=uf[nccx][nccy][jc];
  }

    //test rstrct
    for (jc=1,ic=1;ic<ncz;ic++,jc+=2) {

    }
    printf("we will test restriction here\n");

    */

}

/*************************************************/

/*3D*/
//periodic /= dirichlet
void rstrct0(double ***uc, double *uf, int ncx, int ncy, int ncz)
{
    int ic,iif,jc,jf,kc,kf,nccx, nccy, nccz;
    /*set ncc to be the dimension of the larger grid*/
    nccx = 2*ncx-2+1; /*previously*/
    nccy = 2*ncy-2+1; /*previously*/
    nccz = 2*ncz-2+1;
    double wf=1.0/12;

    //prepare boundaries in rho
    double *** rhord;
    int ii,jj,kk;

    rhord=d3tensor(0,(nccx-1)+2,0,(nccy-1)+2,0,(nccz-1)+2);
    {
        for (ii=0;ii<nccx;ii++)
            for (jj=0;jj<nccy;jj++)
                for(kk=0;kk<nccz;kk++)
                    rhord[ii+1][jj+1][kk+1]=uf[ix(0,ii,jj,kk)];
        //u[ix(0,ii,jj,kk)];
    }


    {printf("in restriction1 rescricted ngx %d ngy %d ngz %d\n", nccx, nccy, nccz);
        for(ii=0; ii<nccx+1; ii++)
            for(jj=0; jj<nccy+1; jj++)
            {
                for(kk=0; kk<nccz+1; kk++)
                    ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //				printf("%E\t", rhord[ii][jj][kk]);
                //			printf("\n");
            }}
    //getchar();

    //periodic boundaries rhord
    //periodic boundaries
    //bottom
    int i,j;
    for(i=0; i<nccy+1; i++)
        for(j=0; j<nccz+1; j++) {
            rhord[0][i][j]=rhord[nccx-1][i][j];
            rhord[nccx+1][i][j]=rhord[2][i][j];
        }
    for(i=0; i<nccx+1; i++)
        for(j=0; j<nccy+1; j++) {
            rhord[i][j][0]=rhord[i][j][nccz-1];
            rhord[i][j][nccz+1]=rhord[i][j][2];
        }
    for(i=0; i<nccx+1; i++)
        for(j=0; j<nccz+1; j++) {
            rhord[i][0][j]=rhord[i][nccy-1][j];
            rhord[i][nccy+1][j]=rhord[i][2][j];
        }


    {printf("in restriction restricted ngx %d ngy %d ngz %d\n", nccx, nccy, nccz);
        for(ii=0; ii<nccx+1; ii++)
        {
            for(jj=0; jj<nccy+1; jj++)
            {
                for(kk=0; kk<nccz+1; kk++)
                    ;	//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //					printf("%E\t", rhord[ii][jj][kk]);
                //			printf("\n");
            }
            //		printf("\n");
        }
    }


    //getchar();
    //from 1 to ncz
    //for(kf=3,kc=2;kc<ncz;kc++,kf+=2)
    //	for (jf=3,jc=2;jc<ncy;jc++,jf+=2){
    //		for (iif=3,ic=2;ic<ncx;ic++,iif+=2){
    //		uc[ic][jc][kc]=0.5*rhord[iif][jf][kf]+wf*(rhord[iif+1][jf][kf]+rhord[iif-1][jf][kf]+rhord[iif][jf+1][kf]+rhord[iif][jf-1][kf]+rhord[iif][jf][kf+1]+rhord[iif][jf][kf-1]);
    //		}
    //	}

    for(kf=1,kc=1;kc<ncz+1;kc++,kf+=2)
        for (jf=1,jc=1;jc<ncy+1;jc++,jf+=2){
            for (iif=1,ic=1;ic<ncx+1;ic++,iif+=2){
                //	printf("%E\t%E\n",uf[ix(0,iif,jf,kf)], uc[ic][jc][kc]);
                //	  printf("%d %d %d %E\n", ic, jc, kc, uf[ix(0,iif,jf,kf)]);
                uc[ic][jc][kc]=0.5*rhord[iif][jf][kf]+wf*(rhord[iif+1][jf][kf]+rhord[iif-1][jf][kf]+rhord[iif][jf+1][kf]+rhord[iif][jf-1][kf]+rhord[iif][jf][kf+1]+rhord[iif][jf][kf-1]);
                //	printf("%E\t%E\n",uf[ix(0,iif,jf,kf)], uc[ic][jc][kc]);
            }
        }
    //  printf("in rstrct0 buss errrror byl here nncx %d \n", nccx);
    //getchar();

    {
        //printf("!!! in restriction weighted ngx %d ngy %d ngz %d\n", ncx, ncy, ncz);
        for(ii=0; ii<ncx+2; ii++)
        {
            for(jj=0; jj<ncy+2; jj++)
            {
                for(kk=0; kk<ncz+2; kk++)
                    ;//printf("%E\t", u[ix(0,zi,zj,zk)]);
                //				printf("%E\t", uc[ii][jj][kk]);
                //			printf("\n");
            }
            //		printf("\n");
        }
    }//	getchar();

    /*boundary points*/
    //accounted for ghost points 0 -> 1, ngy-1 -> ngy
    //but at uf we work in old coordinates
    /*planes*/
    /*do bottom and top plane*/

    /*
  for (jf=2,jc=1;jc<ncy;jc++,jf+=2){
    for (iif=2,ic=1;ic<ncx;ic++,iif+=2){
        uc[ic+1][jc+1][1]=uf[ix(0,iif,jf,0)]=2;
        uc[ic+1][jc+1][ncz]=uf[ix(0,iif,jf,nccz-1)]=2;
    }
  }

  //do face and back (depth -> y)
  for (kf=2,kc=1;kc<ncz;kc++,kf+=2){
    for (iif=2,ic=1;ic<ncx;ic++,iif+=2){
        uc[ic+1][1][kc+1]=uf[ix(0,iif,0,kf)];
        uc[ic+1][ncz][kc+1]=uf[ix(0,iif,nccy-1,kf)];
    }
  }
  //do left and right (depth -> y)
  for (kf=2,kc=1;kc<ncz;kc++,kf+=2){
    for (jf=2,jc=1;jc<ncy;jc++,jf+=2){
        uc[1][jc+1][kc+1]=uf[ix(0,0,jf,kf)];
        uc[ncx][jc+1][kc+1]=uf[ix(0,nccx-1,jf,kf)];
    }
  }

  //do edges
  //bottom and top along x
  for (jc=0,ic=0;ic<ncx;ic++,jc+=2) {
    uc[ic+1][1][1]=uf[ix(0,jc,0,0)];
    uc[ic+1][ncy][1]=uf[ix(0,jc,nccy-1,0)];
    uc[ic+1][1][ncz]=uf[ix(0,jc,0,nccz-1)];
    uc[ic+1][ncy][ncz]=uf[ix(0,jc,nccy-1,nccz-1)];
  }

  //bottom and top along y
  for (jc=0,ic=0;ic<ncy;ic++,jc+=2) {
    uc[1][ic+1][1]=uf[ix(0,0,jc,0)];
    uc[ncx][ic+1][1]=uf[ix(0,nccx-1,jc,0)];
    uc[1][ic+1][ncz]=uf[ix(0,0,jc,nccz-1)];
    uc[ncx][ic+1][ncz]=uf[ix(0,nccx-1,jc,nccz-1)];
  }

  //vertically
  for (jc=0,ic=0;ic<ncz;ic++,jc+=2) {
    uc[1][1][ic+1]=uf[ix(0,0,0,jc)];
    uc[ncx][1][ic+1]=uf[ix(0,nccx-1,0,jc)];
    uc[1][ncy][ic+1]=uf[ix(0,0,nccy-1,jc)];
    uc[ncx][ncy][ic+1]=uf[ix(0,nccx-1,nccy-1,jc)];
  }
*/

    free_d3tensor(rhord, 0,(nccx-1)+2,0,(nccy-1)+2,0,(nccz-1)+2);


}

/*************************************************/

/*3D ??? test indices here else ti should be OK*/
//periodic == dirichlet
void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz)
{
    int ic,iif,jc,jf,kf,kc,ncx,ncy,ncz,jk;
    ncx=nfx/2+1;
    ncy=nfy/2+1;
    ncz=nfz/2+1;
    printf("in interp ncx %d\n", ncx); //getchar();

    /*do copies*/
    for(kc=1, kf=1; kc<ncz+1; kc++, kf+=2) //kc kf
        for (jc=1,jf=1;jc<ncy+1;jc++,jf+=2) //jc jf
            for (ic=1;ic<ncx+1;ic++)  //ic
            {

                uf[2*ic-1][jf][kf]=uc[ic][jc][kc];
                //  printf("copies %d %d %d %d %d %d %E %E\n", 2*ic, jf, kf, ic, jc, kc, uf[2*ic][jf][kf],uc[ic][jc][kc]);
            }
    // printf("in interp OK\n");
    /*interp every second slice*/
    for(kf=1; kf<nfz+1; kf+=2)
    {
        for(jf=1;jf<nfy+1;jf+=2)
            for(iif=2;iif<nfx;iif+=2)
                uf[iif][jf][kf]=0.5*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]);
        for(jf=2;jf<nfy;jf+=2)
            for(iif=1;iif <nfx+1;iif++)
                uf[iif][jf][kf]=0.5*(uf[iif][jf+1][kf]+uf[iif][jf-1][kf]);
    }
    /*interpolate betwen slices in depth*/
    for(kf=2; kf<nfz; kf+=2)
    {
        for(jf=1;jf<nfy+1;jf++)
            for(iif=1; iif<nfx+1; iif++)
                uf[iif][jf][kf]=0.5*(uf[iif][jf][kf-1]+uf[iif][jf][kf+1]);
    }
    int i,j,k;
    for(i=1;i<nfx+1;i++)
        for(j=1;j<nfy+1;j++)
            for(k=1; k<nfz+1; k++)
            {
                //	printf("in loop\n");
                //	printf("%E\n", uf[i][j][k]);
                //	aout[i][j][k]=ain[i][j][k];
                //	printf("after %d %d\n", ain[i][j][k], aout[i][j][k]);
            }

    for(i=0; i<nfy+1; i++)
        for(j=0; j<nfz+1; j++) {
            uf[0][i][j]=uf[nfx-1][i][j];
            uf[nfx+1][i][j]=uf[2][i][j];
        }
    for(i=0; i<nfx+1; i++)
        for(j=0; j<nfy+1; j++) {
            uf[i][j][0]=uf[i][j][nfz-1];
            uf[i][j][nfz+1]=uf[i][j][2];
        }
    for(i=0; i<nfx+1; i++)
        for(j=0; j<nfz+1; j++) {
            uf[i][0][j]=uf[i][nfy-1][j];
            uf[i][nfy+1][j]=uf[i][2][j];
        }

    //	printf("END interp\n");
}
/**************************************************/
/*3D*/
//periodic == dirichlet
void addint(double ***uf, double ***uc, double ***res, int nfx, int nfy, int nfz)
{
    void interp(double ***uf, double ***uc, int nfx, int nfy, int nfz);
    int i,j,k;

    interp(res,uc,nfx,nfy,nfz);
    for(i=0;i<nfx+2;i++)
        for (j=0;j<nfy+2;j++)
            for(k=0;k<nfz+2;k++)
                uf[i][j][k] += res[i][j][k];


    for(i=0; i<nfy+1; i++)
        for(j=0; j<nfz+1; j++) {
            uf[0][i][j]=uf[nfx-1][i][j];
            uf[nfx+1][i][j]=uf[2][i][j];
        }
    for(i=0; i<nfx+1; i++)
        for(j=0; j<nfy+1; j++) {
            uf[i][j][0]=uf[i][j][nfz-1];
            uf[i][j][nfz+1]=uf[i][j][2];
        }
    for(i=0; i<nfx+1; i++)
        for(j=0; j<nfz+1; j++) {
            uf[i][0][j]=uf[i][nfy-1][j];
            uf[i][nfy+1][j]=uf[i][2][j];
        }
}
/**************************************************/
/*3D*/
//periodic work here
void slvsml(double ***u, double ***rhs)
{
    void fill0(double ***u, int nx, int ny, int nz);
    double h2=Lx*0.5;
    double hx,hy;

    //ghosts are accounted for
    fill0(u,3,3,3);
    int nx,ny,nz;


    //dirichlet was
    //  u[1][1][1] = h*h*rhs[1][1][1]/6.0;

    //do a few iterations
    int iter=0;
    int i,ipass,isw,j,jsw=1,k, sweep;
    double factor=1.0/6;



    //maybe should be inside the loop?
    //u[1][1][1]=0;
    //	u[1][1][1]=0;
    nx=ny=nz=3;

    hx=Lx/(nx-1);
    hy=Ly/(ny-1);
    h2=hx*hy;
    while(iter < 10000) {

        iter++;
        //periodic boundaries
        //bottom
        for(i=0; i<ny+1; i++)
            for(j=0; j<nz+1; j++) {
                u[0][i][j]=u[nx-1][i][j];
                u[nx+1][i][j]=u[2][i][j];
            }
        for(i=0; i<nx+1; i++)
            for(j=0; j<ny+1; j++) {
                u[i][j][0]=u[i][j][nz-1];
                u[i][j][nz+1]=u[i][j][2];
            }
        for(i=0; i<nx+1; i++)
            for(j=0; j<nz+1; j++) {
                u[i][0][j]=u[i][ny-1][j];
                u[i][ny+1][j]=u[i][2][j];
            }


        for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
            sweep=jsw;
            for (i=1;i<nx+1;i++)
                for (j=1;j<ny+1;j+=1)
                    for(k=(2-(i+j+sweep)%2); k<nz+1; k+=2)
                    {
                        u[i][j][k]=factor*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]+h2*rhs[i][j][k]);
                    }
            //shift to avoid large numbers in the loop...
            double shift=u[2][2][2];

            for (i=0;i<nx+2;i++)
                for (j=0;j<ny+2;j++)
                    for(k=0; k<nz+2; k++)
                    {
                        //u[i][j][k]-=u[2][2][2];
                        u[i][j][k]-=shift;

                        //printf("u[2][2][2] %E\n",u[i][j][k]);
                    }

        }
    }


    //exit(1);


    //fill0(u,3,3,3);

}

/**************************************************/
/*3D*/
void slvsml2(double ***u, double ***rhs, int nx, int ny, int nz)
{
    int i,ipass,isw,j,jsw=1,k, sweep;
    double hx,hy,h2;
    double toler, factor, error, errorcheck;
    void fill0(double ***u, int nx, int ny, int nz);
    int iter=0;
    fill0(u,nx,ny,nz);
    hx=Lx/(nx-1);
    hy=Ly/(ny-1);
    h2=hx*hy;
    toler=0.00001;
    factor=1.0/6;


    do
    {

        //periodic boundaries
        //bottom
        for(i=0; i<ny+1; i++)
            for(j=0; j<nz+1; j++) {
                u[0][i][j]=u[nx-1][i][j];
                u[nx+1][i][j]=u[2][i][j];
            }
        for(i=0; i<nx+1; i++)
            for(j=0; j<ny+1; j++) {
                u[i][j][0]=u[i][j][nz-1];
                u[i][j][nz+1]=u[i][j][2];
            }
        for(i=0; i<nx+1; i++)
            for(j=0; j<nz+1; j++) {
                u[i][0][j]=u[i][ny-1][j];
                u[i][ny+1][j]=u[i][2][j];
            }

        iter++;
        errorcheck=0.5*toler;
        for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {


            sweep=jsw;
            for (i=1;i<nx+1;i++)
                for (j=1;j<ny+1;j+=1)
                    for(k=(2-(i+j+sweep)%2); k<nz+1; k+=2)
                    {
                        error=u[i][j][k];
                        u[i][j][k]=factor*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]+h2*rhs[i][j][k]);
                        error=u[i][j][k]-error;
                        if(errorcheck <	fabs(error))
                            errorcheck=fabs(error);
                    }
            //shift to avoid large numbers in the loop...
            double shift=u[2][2][2];
            for (i=0;i<nx+2;i++)
                for (j=0;j<ny+2;j++)
                    for(k=0; k<nz+2; k++)
                    {
                        //u[i][j][k]-=u[2][2][2];
                        u[i][j][k]-=shift;

                    }
        }
    }
    while(iter < 10000);
    //while(errorcheck>toler);
}

/**************************************************/
/*3D*/
void relax(double ***u, double ***rhs, int nx, int ny, int nz)
{
    int i,ipass,isw,j,jsw=1,k, sweep;
    double hx,hy,h2,factor;
    hx=Lx/(nx-1);
    hy=Ly/(ny-1);
    h2=hx*hy;
    factor=1.0/6;

    //periodic
    for(i=0; i<ny+1; i++)
        for(j=0; j<nz+1; j++) {
            u[0][i][j]=u[nx-1][i][j];
            u[nx+1][i][j]=u[2][i][j];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<ny+1; j++) {
            u[i][j][0]=u[i][j][nz-1];
            u[i][j][nz+1]=u[i][j][2];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<nz+1; j++) {
            u[i][0][j]=u[i][ny-1][j];
            u[i][ny+1][j]=u[i][2][j];
        }


    //go red first
    for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
        sweep=jsw;
        for (i=1;i<nx+1;i++)
            for (j=1;j<ny+1;j+=1)
                for(k=(2-(i+j+sweep)%2); k<nz+1; k+=2)
                    u[i][j][k]=factor*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k-1]+u[i][j][k+1]+h2*rhs[i][j][k]);
        //periodic	sides

        //we could do it here, but does not seem to be necessary - however could turn out to be important in future
        double shift=u[2][2][2];

        for (i=0;i<nx+2;i++)
            for (j=0;j<ny+2;j++)
                for(k=0; k<nz+2; k++)
                    ;//u[i][j][k]-=shift;
    }
}

/**************************************************/
/*3D*/
void resid(double ***res, double ***u, double ***rhs, int nx, int ny, int nz)
{
    int i,j,k,jf,jc, iif, ic, kf, kc;
    double hx,hy,h2i;
    //PERIODIC
    hx=Lx/(nx-1);
    hy=Ly/(ny-1);
    h2i=1.0/(hx*hy);

    //periodic

    for(i=0; i<ny+1; i++)
        for(j=0; j<nz+1; j++) {
            u[0][i][j]=u[nx-1][i][j];
            u[nx+1][i][j]=u[2][i][j];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<ny+1; j++) {
            u[i][j][0]=u[i][j][nz-1];
            u[i][j][nz+1]=u[i][j][2];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<nz+1; j++) {
            u[i][0][j]=u[i][ny-1][j];
            u[i][ny+1][j]=u[i][2][j];
        }



    /*interior*/
    for (i=1; i<nx+1;i++)
        for (j=1;j<ny+1;j++)
            for(k=1;k<nz+1;k++)
                res[i][j][k] = h2i*(u[i+1][j][k]+u[i-1][j][k]+u[i][j+1][k]+u[i][j-1][k]+u[i][j][k+1]+u[i][j][k-1]-6.0*u[i][j][k])+rhs[i][j][k];

    /*boundary points*/
    // for (i=0;i<nx;i++)
    //   res[i][0]=res[i][ny-1]=0.0;
    //  for (i=0;i<ny;i++)
    //   res[0][i]=res[nx-1][i]=0.0;

    /*boundary points*/
    /*planes*/
    /*do bottom and top plane*/
    for (jc=0;jc<ny+1;jc++){
        for (ic=0;ic<nx+1;ic++){
            res[ic][jc][0]=res[ic][jc][nz-1];
            //res[ic][jc][nz-1]=res[ic][jc][nz-1];
        }
    }
    /*do face and back (depth -> y)*/
    for (kc=0;kc<nz+1;kc++){
        for (ic=0;ic<nx+1;ic++){
            res[ic][0][kc]=res[ic][nz-1][kc];
        }
    }
    /*do left and right (depth -> y)*/
    for (kc=0;kc<nz+1;kc++){
        for (jc=0;jc<ny+1;jc++){
            res[0][jc][kc]=res[nx-1][jc][kc];
        }
    }

    /*do edges*/
    //bottom and top along x
    for (ic=0;ic<nx+1;ic++) {
        res[ic][0][0]=res[ic][ny-1][nz-1];
        res[ic][ny+1][0]=res[ic][2][nz-1];
    }

    //bottom and top along y
    for (ic=0;ic<ny+1;ic++) {
        res[0][ic][0]=res[nx+1][ic][nz-1];
        res[nx+1][ic][0]=res[2][ic][nz-1];
    }

    //vertically
    for (ic=0;ic<nz+1;ic++) {
        res[0][0][ic]=res[nx-1][ny-1][ic];
        res[nx+1][0][ic]=res[2][ny-1][ic];
    }
}


/***************************************************/
/*3D*/
void copy(double ***aout, double ***ain, int nx, int ny, int nz)
{
    //  printf("in copy %E %E\n", ain[0][ny-1][nz-1], ain[0][ny-1][nz-1]);
    //aout[1][ny-1][nz-1]=ain[1][ny-1][nz-1];

    //  printf("in copy %E %E\n", aout[1][ny-1][nz-1], ain[1][ny-1][nz-1]);
    //accounting for ghosts

    int i,j,kk;
    for(i=0;i<nx+2;i++)
        for(j=0;j<ny+2;j++)
            for(kk=0; kk<nz+2; kk++)
            {
                //	  printf("in loop ijk %d %d %d\n", i,j,kk);
                //  printf("%E %E\n", ain[i][j][kk], aout[i][j][kk]);
                aout[i][j][kk]=ain[i][j][kk];
                //	printf("after %E %E\n", ain[i][j][k], aout[i][j][k]);
            }
    //  printf("finished\n");
}
/***************************************************/
/*3D*/
void copy0(double ***aout, double *ain, int nx, int ny, int nz)
{
    int i,j,k;
    //accounting for ghosts

    for (i=0;i<nx;i++)
        for (j=0;j<ny;j++)
            for(k=0;k<nz;k++){
                //	printf("in copy %E %E\n", ain[ix(0,i,j,k)], aout[i][j][k]);
                aout[i+1][j+1][k+1]=ain[ix(0,i,j,k)];}

    //boundaries - ghost points
    for(i=0; i<ny+1; i++)
        for(j=0; j<nz+1; j++) {
            aout[0][i][j]=aout[nx-1][i][j];
            aout[nx+1][i][j]=aout[2][i][j];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<ny+1; j++) {
            aout[i][j][0]=aout[i][j][nz-1];
            aout[i][j][nz+1]=aout[i][j][2];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<nz+1; j++) {
            aout[i][0][j]=aout[i][ny-1][j];
            aout[i][ny+1][j]=aout[i][2][j];
        }
    printf("copied 0\n");
    //minus sign is needed cause we have L*u=rhs, and our rhs should be -rho
}
/***************************************************/
/*3D*/
void copyfinal(double *aout, double ***ain, int nx, int ny, int nz)
{

    int i,j,k;

    for(i=0; i<ny+1; i++)
        for(j=0; j<nz+1; j++) {
            ain[0][i][j]=ain[nx-1][i][j];
            ain[nx+1][i][j]=ain[2][i][j];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<ny+1; j++) {
            ain[i][j][0]=ain[i][j][nz-1];
            ain[i][j][nz+1]=ain[i][j][2];
        }
    for(i=0; i<nx+1; i++)
        for(j=0; j<nz+1; j++) {
            ain[i][0][j]=ain[i][ny-1][j];
            ain[i][ny+1][j]=ain[i][2][j];
        }

    printf("copyfinal ain[1][1][1] %E\n", ain[1][1][1]);
    printf("copyfinal ain[ngx][ngy][ngz] %E\n", ain[nx][ny][nz]);

    //03.09.2013 modifications WJM -> getting the average potential as the reference.
    double shift=0.0;
    int counterr=0;
    //ain[1][1][1];
    //getchar();
    //take only interior points (exclude ghost points!!!) to avoid bias
    for (i=1;i<nx+1;i++)
        for (j=1;j<ny+1;j++)
            for(k=1; k<nz+1; k++)
            {
                shift+=ain[i][j][k];
                counterr++;
            }
    shift/=counterr;

    //ain[i][j][k]-=ain[1][1][1];
    // ain[i][j][k]-=shift;
    //return mean zero potential;

    //again only interior points
    for (i=1;i<nx+1;i++)
        for (j=1;j<ny+1;j++)
            for(k=1; k<nz+1; k++)
                ain[i][j][k]-=shift;
    //test
    shift=0.0;
    for (i=1;i<nx+1;i++)
        for (j=1;j<ny+1;j++)
            for(k=1; k<nz+1; k++)
                //ain[i][j][k]-=ain[1][1][1];
                shift+=ain[i][j][k];
    printf("this is mean potential %E\n", shift);
    //	getchar();


    printf("in copyfinal %E\n", ain[3][3][3]);
    printf("in copyfinal %E\n", ain[nx+1][ny+1][nz+1]);
    printf("in copyfinal 1 1 1 %E\n", ain[1][1][1]);

    //only inerior points are returned
    for (i=1;i<nx+1;i++)
        for (j=1;j<ny+1;j++)
            for(k=1;k<nz+1;k++)
            {
                aout[ix(0,i-1,j-1,k-1)]=ain[i][j][k];
                //		  printf()
            }
}
/***************************************************/

/*3D*/
void fill0(double ***u, int nx, int ny, int nz)
{
    int i,j,k;
    //accounting for ghosts
    for(i=0;i<nx+2;i++)
        for(j=0;j<ny+2;j++)
            for(k=0;k<nz+2;k++)
                u[i][j][k]=0.0;
}
/***************************************************/
