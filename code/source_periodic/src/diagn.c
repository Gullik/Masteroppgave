/* DiP3D program */
/* Diagnostics*/
/* Author: Wojciech Jacek Miloch */
/* University of Oslo, Norway */
/* 2009 */
#include"const.h"

/*
  POINTERS' LEGEND
  2 DIMENSIONAL DIAGNOSTICS  
  (pot2Dclr); //potential clr
  (pot2D); //potential
  (efx); //Efield x
  (efy); //E field y
  (frho); //charge density
  (idens); //ion density
  (edens); //electron density
  (epe); // ele potetnial energy
  (eke); // ele kinetic energy
  (ipe); // ion potetnial energy
  (ike); // ion kinetic energy
*/

void diagn_open()
{
  //OPEN FILES FOR DIAGNOSTICS
  if(rank==0)
    {
      //TWO DIMESIONAL DIAGNOSTICS open and print scale!
      pot2Dclr=my_file_open("./data/pot3Dclr.dat", "w");
      printscale(pot2Dclr);      
      pot2D=my_file_open("./data/pot3D.dat", "w");
      printscale(pot2D);
      pot2Dav=my_file_open("./data/pot3Dav.dat", "w");
      printscale(pot2Dav);
      efx=my_file_open("./data/efieldx.dat", "w");
      printscale(efx);      
      efy=my_file_open("./data/efieldy.dat", "w");
      printscale(efy);      
      efz=my_file_open("./data/efieldz.dat", "w");
      printscale(efz);   
      frho=my_file_open("./data/density.dat", "w");
      printscale(frho);      
      idens=my_file_open("./data/idensity.dat", "w");
      printscale(idens);      
      edens=my_file_open("./data/edensity.dat", "w");
      printscale(edens);	
      bdens=my_file_open("./data/bdensity.dat", "w");
      printscale(bdens);  
      iavvel=my_file_open("./data/iavvel.dat", "w");
      printscale(iavvel);	  
      eavvel=my_file_open("./data/eavvel.dat", "w");
      printscale(eavvel);
      epe=my_file_open("./data/PEelectron.dat", "w");
      printscale(epe); 
      ipe=my_file_open("./data/PEion.dat", "w");
      printscale(ipe);
      pe=my_file_open("./data/PEtotal.dat", "w");
      printscale(pe);      
      pe_time=my_file_open("./data/PEtime.dat", "w");
      epe_time=my_file_open("./data/PEelec_time.dat", "w");
      ipe_time=my_file_open("./data/PEion_time.dat", "w");	
    }
  //KINETIC ENERGIES
  if((numtasks>1)&& (rank==1))
    {
      eke=my_file_open("./data/KEelectron.dat", "w");
      printscale(eke);
      ike=my_file_open("./data/KEion.dat", "w");	
      printscale(ike); 
      eke_time=my_file_open("./data/KEelec_time.dat", "w");
      ike_time=my_file_open("./data/KEion_time.dat", "w");	
    }
  else 
    if(rank==0)
      {
	eke=my_file_open("./data/KEelectron.dat", "w");
	printscale(eke);
	ike=my_file_open("./data/KEion.dat", "w");	  
	printscale(ike); 
	eke_time=my_file_open("./data/KEelec_time.dat", "w");
	ike_time=my_file_open("./data/KEion_time.dat", "w");
      }
  
  if(rank==0)
    {
      dustshape=my_file_open("./data/dustshape.dat", "w");
      dustshapet=my_file_open("./data/dustshapetime.dat", "w");
      //initial!
      printdustshape();      
      dustcharge=my_file_open("./data/dustcharge.dat", "w");
      printscale(dustcharge);
	  newprobe=my_file_open("./data/newprobe.dat", "w");
	  printscale(newprobe);
		dust_time=my_file_open("./data/dustcharge_time.dat", "a");
      convergence=my_file_open("./data/convpot_time.dat", "w");	
      dhist=my_file_open("./data/dust_history.dat", "w");      
      numberofprints=0;
      dens_err=my_file_open("./data/dens_err.dat","w");
    
    }
  //phase-space -> open on all nodes
  // evxphs=my_file_open("./data/evxphs.dat", "a");
  // ivxphs=my_file_open("./data/ivxphs.dat", "a"); 
  //evxphs2=my_file_open("./data/evxphs2.dat", "a");
  // ivxphs2=my_file_open("./data/ivxphs2.dat", "a");
  // evxphs3=my_file_open("./data/evxphs3.dat", "a");
  //  ivxphs3=my_file_open("./data/ivxphs3.dat", "a");     
}

void diagn_close()
{  
  if(rank==0)
    {
      //NEW DIAGN 2D
      fclose(pot2Dclr); //potential clr
      fclose(pot2D); //potential
      fclose(pot2Dav); //potential av
      fclose(efx); //Efield x
      fclose(efy); //E field y
      fclose(efz);
      fclose(frho); //charge density
      fclose(idens); //ion density
      fclose(edens); //electron density
      fclose(bdens);
      fclose(epe); //potetnial energy electron
      fclose(ipe); //potetnial energy ion
      fclose(pe);
      fclose(pe_time);      
      fclose(epe_time);
      fclose(ipe_time);
      fclose(dustcharge);
      fclose(newprobe);
		fclose(dust_time);
      fclose(convergence);
      fclose(dhist);
      //but printalso number of diagnostics prints
      FILE *timesteps;
      timesteps=my_file_open("./data/timesteps.dat", "w");
      fprintf(timesteps, "%d\n", numberofprints);
      fclose(timesteps);
      fclose(dens_err);
      fclose(dustshape);  
      fclose(dustshapet);
      fclose(eavvel);
      fclose(iavvel);
      //fclose(evxphs);
      //fclose(ivxphs);
    }
  if((numtasks>1)&& (rank==1))
    {
      fclose(eke);     //kinetic energy electron
      fclose(ike);     //kinetic energy ion
      fclose(eke_time);
      fclose(ike_time);
    }
  else 
    if(rank==0)
      {
	fclose(eke);
	fclose(ike);
	fclose(eke_time);
	fclose(ike_time);
      }
  if(rank==0)
    {
 
      FILE *timer;  
      timer=my_file_open("./data/timer.dat", "w");
      fprintf(timer, "%d\n", timerprobes);
      fclose(timer);
    }
}


/***PRINT GRID QUANTITIES FOR A GIVEN TIME***/
void printgrid(int t)
{
  int hj;
  int i,l; //particles
  int j,k,ll; //grid
  int di,dj;
  double q;
  double x,y,z;
  double pdensity;
  double temp;
  int offset; 
  double newweight;
  //printf("average %d in time %d",average,t);
  average=1;
  for(hj=0; hj< diagint; hj++)
    if(t>= diagint_st[hj])
      {
	//printf("dddhj %d and %d st %d\n",hj, diagint_av[hj], diagint_st[hj]);
	average=diagint_av[hj];
      }
  printf("rank %d DIAGN: Average %d in time %d\n", rank, average,t);
  
  weight = 1/(1.0*average);
  //CLEAR GRID DIAGNOSTICS
  //  int counter;
  
  if(t==0) //initially clear vectors for diagnostics
    {
      /*clear vectors*/
      for(i=0;i<pdensMAX; i++)
	pdens[i]=0.0;
      for(i=0;i<KEMAX; i++)
	KE[i]=0.0;
      if(rank==0)
	{
	  for(i=0;i<PEMAX; i++)
	    PE[i]=0.0;
	  for(i=0;i<PEtotalMAX; i++)
	    PEtotal[i]=0.0;
	  for(i=0;i<qdensMAX; i++)
	    qdens[i]=0.0;
	  for(i=0;i<phiavMAX; i++)
	    phiav[i]=0.0;	
	}
    }


  int pphs=0;
  int spp;
  int delt=15;
  spp=20;
#ifdef MPI
  if(rank==1){
#endif
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs20.dat", "a");
      ivxphs=my_file_open("./data/ivxphs20.dat", "a"); 
      pphs=1;
    }
  spp=1000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs1000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs1000.dat", "a"); 
      pphs=1;
    }
  spp=2500;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs2500.dat", "a");
      ivxphs=my_file_open("./data/ivxphs2500.dat", "a"); 
      pphs=1;
    }

  spp=3000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs3000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs3000.dat", "a"); 
      pphs=1;
    }
  
 spp=4000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs4000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs4000.dat", "a"); 
      pphs=1;
    }

 spp=5000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs5000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs5000.dat", "a"); 
      pphs=1;
    }
 spp=7000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs7000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs7000.dat", "a"); 
      pphs=1;
    }
 spp=10000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs10000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs10000.dat", "a"); 
      pphs=1;
    }
 spp=13000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs13000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs13000.dat", "a"); 
      pphs=1;
    }

 spp=15000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs15000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs15000.dat", "a"); 
      pphs=1;
    }

 spp=17000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs17000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs17000.dat", "a"); 
      pphs=1;
    }

 spp=20000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs20000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs20000.dat", "a"); 
      pphs=1;
    }
 spp=25000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs25000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs25000.dat", "a"); 
      pphs=1;
    }

 spp=30000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs30000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs30000.dat", "a"); 
      pphs=1;
    }


 spp=35000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs35000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs35000.dat", "a"); 
      pphs=1;
    }


 spp=40000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs40000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs40000.dat", "a"); 
      pphs=1;
    }


 spp=45000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs45000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs45000.dat", "a"); 
      pphs=1;
    }


 spp=50000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs50000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs50000.dat", "a"); 
      pphs=1;
    }


 spp=55000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs55000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs55000.dat", "a"); 
      pphs=1;
    }


 spp=60000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs60000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs60000.dat", "a"); 
      pphs=1;
    }


 spp=65000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs65000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs65000.dat", "a"); 
      pphs=1;
    }


 spp=70000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs70000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs70000.dat", "a"); 
      pphs=1;
    }


 spp=75000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs75000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs75000.dat", "a"); 
      pphs=1;
    }


 spp=80000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs80000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs80000.dat", "a"); 
      pphs=1;
    }


 spp=85000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs85000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs85000.dat", "a"); 
      pphs=1;
    }


 spp=90000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs90000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs90000.dat", "a"); 
      pphs=1;
    }


 spp=95000;
  if(t>=spp && t < spp + delt)
    {
      evxphs=my_file_open("./data/evxphs95000.dat", "a");
      ivxphs=my_file_open("./data/ivxphs95000.dat", "a"); 
      pphs=1;
    }


#ifdef MPI
  }
#endif

  for(i=0;i<S;i++) //collect data
    {
      double vx,vy,vz;
      offset=pdens_off*i; 
      q=chargeandnorm[i]/(dV);
      pdensity=1/((dV)*(dV)); //one dV for weighting second for density: rho=dpart/dV
      //  printf("assigning i: %d  npart[i] %d\n", i, npart[i]);
      //  getchar();
      for(l=0;l<npart[i];l++)
	{                      	  
	  vx=spec[i].part[l].vx;
	  vy=spec[i].part[l].vy;
	  vz=spec[i].part[l].vz;
	  j=spec[i].part[l].x/dx; //take the integer
	  k=spec[i].part[l].y/dy;
	  ll=spec[i].part[l].z/dz;
	  x=spec[i].part[l].x-j*dx;
	  y=spec[i].part[l].y-k*dy;	
	  z=spec[i].part[l].z-ll*dz;
	 
	  pdens[ix(offset,j,k,ll)]+=pdensity*(dx-x)*(dy-y)*(dz-z);
	  pdens[ix(offset,j+1,k,ll)]+=pdensity*x*(dy-y)*(dz-z);
	  pdens[ix(offset,j+1,k+1,ll)]+=pdensity*x*y*(dz-z);
	  pdens[ix(offset,j,k+1,ll)]+=pdensity*(dx-x)*y*(dz-z);
	  pdens[ix(offset,j,k,ll+1)]+=pdensity*(dx-x)*(dy-y)*z;
	  pdens[ix(offset,j+1,k,ll+1)]+=pdensity*x*(dy-y)*z;
	  pdens[ix(offset,j+1,k+1,ll+1)]+=pdensity*x*y*z;
	  pdens[ix(offset,j,k+1,ll+1)]+=pdensity*(dx-x)*y*z; 
	
	  vxvec[ix(offset,j,k,ll)]+=vx*pdensity*(dx-x)*(dy-y)*(dz-z);
	  vxvec[ix(offset,j+1,k,ll)]+=vx*pdensity*x*(dy-y)*(dz-z);
	  vxvec[ix(offset,j+1,k+1,ll)]+=vx*pdensity*x*y*(dz-z);
	  vxvec[ix(offset,j,k+1,ll)]+=vx*pdensity*(dx-x)*y*(dz-z);
	  vxvec[ix(offset,j,k,ll+1)]+=vx*pdensity*(dx-x)*(dy-y)*z;
	  vxvec[ix(offset,j+1,k,ll+1)]+=vx*pdensity*x*(dy-y)*z;
	  vxvec[ix(offset,j+1,k+1,ll+1)]+=vx*pdensity*x*y*z;
	  vxvec[ix(offset,j,k+1,ll+1)]+=vx*pdensity*(dx-x)*y*z; 
	
	  vyvec[ix(offset,j,k,ll)]+=vy*pdensity*(dx-x)*(dy-y)*(dz-z);
	  vyvec[ix(offset,j+1,k,ll)]+=vy*pdensity*x*(dy-y)*(dz-z);
	  vyvec[ix(offset,j+1,k+1,ll)]+=vy*pdensity*x*y*(dz-z);
	  vyvec[ix(offset,j,k+1,ll)]+=vy*pdensity*(dx-x)*y*(dz-z);
	  vyvec[ix(offset,j,k,ll+1)]+=vy*pdensity*(dx-x)*(dy-y)*z;
	  vyvec[ix(offset,j+1,k,ll+1)]+=vy*pdensity*x*(dy-y)*z;
	  vyvec[ix(offset,j+1,k+1,ll+1)]+=vy*pdensity*x*y*z;
	  vyvec[ix(offset,j,k+1,ll+1)]+=vy*pdensity*(dx-x)*y*z; 
	  
	  vzvec[ix(offset,j,k,ll)]+=vz*pdensity*(dx-x)*(dy-y)*(dz-z);
	  vzvec[ix(offset,j+1,k,ll)]+=vz*pdensity*x*(dy-y)*(dz-z);
	  vzvec[ix(offset,j+1,k+1,ll)]+=vz*pdensity*x*y*(dz-z);
	  vzvec[ix(offset,j,k+1,ll)]+=vz*pdensity*(dx-x)*y*(dz-z);
	  vzvec[ix(offset,j,k,ll+1)]+=vz*pdensity*(dx-x)*(dy-y)*z;
	  vzvec[ix(offset,j+1,k,ll+1)]+=vz*pdensity*x*(dy-y)*z;
	  vzvec[ix(offset,j+1,k+1,ll+1)]+=vz*pdensity*x*y*z;
	  vzvec[ix(offset,j,k+1,ll+1)]+=vz*pdensity*(dx-x)*y*z; 	  
			
	  //  if(t==337) exit(1);		  		  
	  /*print v-r phase space*/
#ifdef MPI
	  if(pphs==1 && l<6000 && rank == 1)
#else
	  if(pphs==1 && l<6000)
#endif
	    {
	      //  if(j>=25  && j <28) //here I can choose other counters...
	      //	if(j>=15 && j<16 && k>=15 && k <16 && ll>=15 && ll<16)
	      // if(k>=15 && k <16)	
	      // if(ll>=15 && ll<16)
	      //		if(l<1000)
	      {
		printf("PRINTING VX PHSPACE no. %d !!\n\n\n", l);
		if(i==0 || i==3)
		  {
		    temp=spec[i].part[l].x*normx;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		    temp=spec[i].part[l].y*normx;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		    temp=spec[i].part[l].z*normx;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		    temp=spec[i].part[l].vx*normvel;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		    temp=spec[i].part[l].vy*normvel;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		    temp=spec[i].part[l].vz*normvel;
		    fwrite(&temp, sizeof(double), 1,evxphs);
		  } 
		if(i==1 || i==2)
		  {
		    temp=spec[i].part[l].x*normx;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		    temp=spec[i].part[l].y*normx;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		    temp=spec[i].part[l].z*normx;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		    temp=spec[i].part[l].vx*normvel;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		    temp=spec[i].part[l].vy*normvel;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		    temp=spec[i].part[l].vz*normvel;
		    fwrite(&temp, sizeof(double), 1,ivxphs);
		  } 
		printf("FINISHED PRINTING VX PHSPACE!!\n\n\n");
		
	      }

	    }	  
	}
    }
#ifdef MPI 
  if(rank==1){
#endif
  if(pphs==1)
    {
      fclose(evxphs);
      fclose(ivxphs);
    }

#ifdef MPI
    }
#endif


  //  printf("DIAGN.C: data collected\n");
      //collected data from all particles      
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(pdens,rpdens,pdensMAX,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(vxvec,rvxvec,pdensMAX,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(vyvec,rvyvec,pdensMAX,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
  MPI_Reduce(vzvec,rvzvec,pdensMAX,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);  
//  printf("rank %d Reduced part 1\n", rank);
  for(di=0; di<noofdusts; di++)
    for(dj=0; dj<dpartlast[di]; dj++)
      {	
		//printf("heeeh\n");
		//printf("rank %d reducing %E %E\n", rank, dpart[di][dj].q, dpart[di][dj].q);
	//MPI_Reduce(&dpart[di][dj].x,&rdpart[di][dj].x,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	//MPI_Reduce(&dpart[di][dj].y,&rdpart[di][dj].y,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	//MPI_Reduce(&dpart[di][dj].z,&rdpart[di][dj].z,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	   MPI_Reduce(&dpart[di][dj].q,&rdpart[di][dj].q,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	   MPI_Reduce(&dpart[di][dj].ecnt,&rdpart[di][dj].ecnt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	   MPI_Reduce(&dpart[di][dj].icnt,&rdpart[di][dj].icnt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	   MPI_Reduce(&dpart[di][dj].bcnt,&rdpart[di][dj].bcnt,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		  
      }
//printf("rank %d Reduced part 2\n", rank);
#endif
  
  // average the qdensity, potential energies  
  if(rank==0)
    {
      for(i=0; i<qdensMAX; i++)
	{
#ifdef MPI
	  qdens[i]+=rrho[i];       //  gp[i][j].qdens+=rrho[ix(0,i,j)]; //averaged charge density
#else
	  qdens[i]+=rho[i];	       //  gp[i][j].qdens+=rho[ix(0,i,j)];
#endif	      
#ifdef MPI
	  PE[i]+=rpdens[i]*phi[i];
	  PE[PEMAXhalf+i]+=rpdens[PEMAXhalf+i]*phi[i];
	  PEtotal[i]+=rrho[i]*phi[i];
	  // gp[i][j].PE[0]+=rpdens[ix(0,i,j)]*gp[i][j].phi; //to be multipl with charge also
	  // gp[i][j].PE[1]+=rpdens[ix(pdens_off,i,j)]*gp[i][j].phi;
	  //gp[i][j].PEtotal+=rrho[ix(0,i,j)]*gp[i][j].phi;
#else
	  PE[i]+=pdens[i]*phi[i];
	  PE[PEMAXhalf+i]+=pdens[PEMAXhalf+i]*phi[i];
	  PEtotal[i]+=rho[i]*phi[i];	      
	  //  gp[i][j].PE[0]+=pdens[ix(0,i,j)]*gp[i][j].phi; //to be multipl with charge also
	  //  gp[i][j].PE[1]+=pdens[ix(pdens_off,i,j)]*gp[i][j].phi;
	  //  gp[i][j].PEtotal+=rho[ix(0,i,j)]*gp[i][j].phi;
#endif
	  //no reducing here so only on rank 0	      
	  // gp[i][j].phiav+=gp[i][j].phi;
	  phiav[i]+=phi[i];
	}

//modification WJM 18.09.2013 - extra printing
      if((t % average) == 10)
	{
	  newweight=0.1;
          printf("diagnostics printing %d\n", t);// getchar();	      
	  //need to print before density
	  print_avpvel(eavvel, 0,t,  newweight);
	  print_avpvel(iavvel, 1,t,  newweight);
	  printf("in in diagnostics \n");
	  printdensity(edens, 0, t, newweight); //electron density
	  printdensity(idens, 1, t, newweight); //ion density
#ifdef BEAM
 	  printdensity(bdens, 2, t, newweight);
#endif
	  printqdensity(frho, t, newweight); //charge density
	  printefield(efx, 0, t, newweight); //electric field in x direction (0)
       	  printefield(efy, 1, t, newweight); //electric field in y direction (1)
	  printefield(efz, 2, t, newweight); //electric field in y direction (2)
	  printpotential(pot2D, t, newweight); //potential
	  printavpotential(pot2Dav, t, newweight);
	  printPE(epe, epe_time, 0, t, newweight); //Potential energy, electron 0
	  printPE(ipe, ipe_time, 1, t, newweight); //Potential energy, ions 1
	  printPEtotal(pe, pe_time, t, newweight); //Potential energy density+total vs. time
	  printdustcharge(dustcharge,t,newweight);
		printf("print new probe\n");
	  printnewprobe(newprobe, t, newweight);
		printf("printed new probe\n");
		numberofprints++;

	}

      //printf("DIAGN.C: before printing diagnostics t %d average %d\n", t, average);      
      if((t % average) == 0)
	{
	  //       	  printf("DIAGN.C: in printing diagnostics\n");      

	  printf("diagnostics printing %d\n", t);// getchar();	      
	  //need to print before density
	  print_avpvel(eavvel, 0,t,  weight);
	  print_avpvel(iavvel, 1,t,  weight);
	  printf("in in diagnostics \n");
	  printdensity(edens, 0, t, weight); //electron density
	  printdensity(idens, 1, t, weight); //ion density
#ifdef BEAM
 	  printdensity(bdens, 2, t, weight);
#endif
	  printqdensity(frho, t, weight); //charge density
	  printefield(efx, 0, t, weight); //electric field in x direction (0)
       	  printefield(efy, 1, t, weight); //electric field in y direction (1)
	  printefield(efz, 2, t, weight); //electric field in y direction (2)
	  printpotential(pot2D, t, weight); //potential
	  printavpotential(pot2Dav, t, weight);
	  printPE(epe, epe_time, 0, t, weight); //Potential energy, electron 0
	  printPE(ipe, ipe_time, 1, t, weight); //Potential energy, ions 1
	  printPEtotal(pe, pe_time, t, weight); //Potential energy density+total vs. time
	  printdustcharge(dustcharge,t,weight);
		printf("print new probe\n");
	  printnewprobe(newprobe, t, weight);
		printf("printed new probe\n");
		numberofprints++;
	  /*clear vectors on node 1*/	 
	  if(rank==0)
	    {
	      
	      //  for(i=0; i<ngx; i++)
	      //	for(j=0; j<ngy; j++)
	      //	  {	
	      //	    for(counter=0; counter<S; counter++)
	      //	      {
	      //		gp[i][j].PE[counter]=0.0;
	      //	      }
	      
	      //	    gp[i][j].qdens=0.0; //averaged charge density
	      //	    gp[i][j].PEtotal=0.0;
	      //	     gp[i][j].phiav=0.0;
	      //	  }
	      
	      for(i=0; i<phiavMAX; i++)
		{
		  phiav[i]=0.0;
		  PEtotal[i]=0.0;
		  qdens[i]=0.0;	
		}  
	      for(i=0; i<PEMAX; i++)
		{
		  PE[i]=0.0;
		}  
	    }
	  
	  //convergence of the potential unit!!!
	  if(t >= CONVTEST)
	    printconvpot(convergence,t,CONVTEST);
	  for(i=0; i<phiMAX; i++)
	    potconv[t % CONVTEST][i]=phi[i];
	}
    }
/*clear density vectors on all nodes if printed*/
if((t % average)==0)     
  for(i=0;i<pdensMAX; i++)
    pdens[i]=vxvec[i]=vyvec[i]=vzvec[i]=0.0;

// printf("DIAGN.C : after printing diagnostics\n");
	
//clear NEWPROBE
if((t % average)==0)     
	for(i=0; i<noofdusts; i++)
		for(j=0; j<dpartlast[i]; j++)
		{	
			//we are on all ranks
#ifdef MPI			
			rdpart[i][j].ecnt=rdpart[i][j].icnt=rdpart[i][j].bcnt=0;
#endif
			dpart[i][j].ecnt=dpart[i][j].icnt=dpart[i][j].bcnt=0;
		}	

}

    
/********PRINT TIME CONVERGENCE**********/
void printconvpot(FILE *fpointer, int t, int step)
{
  double totalerror=0.0;
  double partial;
  int i;

  for(i=0; i<phiMAX; i++) //phimax == npoints
  {
     partial=phi[i]-potconv[t % step][i];
     totalerror+=partial*partial;
   }
//  fprintf(fpointer,"%E\t%E\n", t*dt*normtime, sqrt(totalerror)*normpot/(ngx*ngy));
  fprintf(fpointer,"%E\t%E\n", t*dt*normtime, sqrt(totalerror)*normpot/(ngx*ngy*ngz));
//fwrite(fpointer);
//here we have error per grid point?

}


/******PRINT SCALE IN 2 DIMENSIONAL DIAGNOSTICS*******/
/*****must be called after normalization*************/
void printscale(FILE *fpointer)
{
 // int i,j;
  //printscale unformatted output - check if it works for one element!
  //number of grid points
  double tempd;
  fwrite(&ngx, sizeof(int), 1, fpointer);
  fwrite(&ngy, sizeof(int), 1, fpointer);
  fwrite(&ngz, sizeof(int), 1, fpointer);
  //Box size -> position of grid point to be calculated in the analysis software: Lx*i/ngx
  tempd=Lx*normx;
  fwrite(&tempd, sizeof(double), 1, fpointer);
  tempd=Lx*normx;
  fwrite(&tempd, sizeof(double), 1, fpointer);
  tempd=Lx*normx;
  fwrite(&tempd, sizeof(double), 1, fpointer);
  //fprintf(fpointer, "%d\n", ngx); //no of grid points
  //fprintf(fpointer, "%d\n", ngy);  

 //for(i=0; i<ngx; i++)
   // fprintf(fpointer, "%E\n", i*dx*normx);
 // for(j=0; j<ngy; j++)
   // fprintf(fpointer, "%E\n", j*dy*normx);
}

/************PRINT CHARGE DENSITY *********************/
void printqdensity(FILE *fpointer, int t, double weight)
{
  int i,j,k,mltx, mlty, mltz; 
  double *temp, time;
  //PRINT DENSITY
  //print time
  //  fprintf(fpointer, "%E\n", t*dt*normtime);
  time=t*dt*normtime;
  fwrite(&time,sizeof(double), 1, fpointer);
  //print grid quantities
  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<ngx; i++)
    { 
      mltx=1;	
      if((i==0) || (i==(ngx-1))) mltx=2;
      for(j=0; j<ngy; j++)
	{
	  mlty=1;
	  if((j==0) || (j==(ngy-1))) mlty=2;
	  for(k=0; k<ngz; k++)
	    {
	      mltz=1;
	      if((k==0) || (k==(ngz-1))) mltz=2;	
	      temp[ix(0,i,j,k)]=weight*qdens[ix(0,i,j,k)]*normqdens*mltx*mlty*mltz;
	    }
	}
    }
  fwrite(temp, sizeof(double), qdensMAX, fpointer);
  free_dvecmem(temp, 0, qdensMAX-1);
//  for(i=0; i<ngx; i++)
 //   for(j=0; j<ngy; j++)
  //    {	
//	mlt=1;
//	if((j==0) || (j==(ngy-1))) mlt*=2;
//	if((i==0) || (i==(ngx-1))) mlt*=2;	
//	fprintf(fpointer,"%E\n", weight*gp[i][j].qdens*normqdens*mlt);
  //    } 
}


/*************PRINT PARTICLE DENSITY******************/
void printdensity(FILE *fpointer, int kk, int t, double weight)
{

  int i,j,k,mltx, mlty, mltz, offset; 
 double *temp, time;
  double normxcube=normx*normx*normx;
  offset=pdens_off*kk;
  //PRINT DENSITY
  //print time
  //  fprintf(fpointer, "%E\n", t*dt*normtime);
  time=t*dt*normtime; 
 fwrite(&time, sizeof(double), 1, fpointer);
//fprintf(fpointer, "%E\n", time);  //print grid quantities
	//03.09.2013 Modifications - to have correct densities printed at the edges ...
#ifdef PERIODIC	
	double mlti;
	printf("DIAGN: periodic weighting, test it\n");
	//getchar();
	//sides (interior)
	int ij,ik, il;
#ifndef MPI	
	mlti=1.0;
	for(ij=1; ij<ngx-1; ij++)
		for(ik=1; ik<ngy-1; ik++)
		{						
			pdens[ix(offset,ij,ik,0)]=mlti*(pdens[ix(offset,ij,ik,0)]+pdens[ix(offset,ij,ik,ngz-1)]);
			pdens[ix(offset,ij,ik,ngz-1)]=pdens[ix(offset,ij,ik,0)];
		}
	for(ik=1; ik<ngy-1; ik++)
		for(il=1; il<ngz-1; il++)
		{
			pdens[ix(offset,0,ik,il)]=mlti*(pdens[ix(offset,0,ik,il)]+pdens[ix(offset,ngx-1,ik,il)]);
			pdens[ix(offset,ngx-1,ik,il)]=pdens[ix(offset,0,ik,il)];
		}
	for(ij=1; ij<ngx-1; ij++)
		for(il=1; il<ngz-1; il++)
		{
			pdens[ix(offset,ij,0,il)]=mlti*(pdens[ix(offset,ij,0,il)]+pdens[ix(offset,ij,ngy-1,il)]);
			pdens[ix(offset,ij,ngy-1,il)]=pdens[ix(offset,ij,0,il)];
		}
	
	//edges
	mlti=1.0;
	for(ij=1; ij<ngx-1; ij++)
	{
		pdens[ix(offset,ij,0,0)]=mlti*(pdens[ix(offset,ij,0,0)]+pdens[ix(offset,ij,ngy-1,ngz-1)]+pdens[ix(offset,ij,0,ngz-1)]+pdens[ix(offset,ij,ngy-1,0)]);				
		pdens[ix(offset,ij,0,ngz-1)]=pdens[ix(offset,ij,ngy-1,0)]=pdens[ix(offset,ij,ngy-1,ngz-1)]=pdens[ix(offset,ij,0,0)];
	}
	for(ik=1; ik<ngy-1; ik++)
	{
		pdens[ix(offset,0,ik,0)]=mlti*(pdens[ix(offset,0,ik,0)]+pdens[ix(offset,ngy-1,ik,ngz-1)]+pdens[ix(offset,0,ik,ngz-1)]+pdens[ix(offset,ngx-1,ik,0)]);				
		pdens[ix(offset,0,ik,ngz-1)]=pdens[ix(offset,ngx-1,ik,0)]=pdens[ix(offset,ngx-1,ik,ngz-1)]=pdens[ix(offset,0,ik,0)];
	}
	for(il=1; il<ngz-1; il++)
	{
		pdens[ix(offset,0,0,il)]=mlti*(pdens[ix(offset,0,0,il)]+pdens[ix(offset,ngx-1,ngy-1,il)]+pdens[ix(offset,ngx-1,0,il)]+pdens[ix(offset,0,ngy-1,il)]);				
		pdens[ix(offset,ngx-1,0,il)]=pdens[ix(offset,0,ngy-1,il)]=pdens[ix(offset,ngx-1,ngy-1,il)]=pdens[ix(offset,0,0,il)];
	}
	mlti=1.0;	//corners
	pdens[ix(offset,0,0,0)]=mlti*(pdens[ix(offset,ngx-1,ngy-1,ngz-1)]+pdens[ix(offset,0,ngy-1,ngz-1)]+pdens[ix(offset,ngx-1,0,ngz-1)]+pdens[ix(offset,ngx-1,ngy-1,0)]
								  +pdens[ix(offset,0,0,ngz-1)]+pdens[ix(offset,0,ngy-1,0)]+pdens[ix(offset,ngx-1,0,0)]+pdens[ix(offset,0,0,0)]);
	pdens[ix(offset,ngx-1,ngy-1,ngz-1)]=pdens[ix(offset,0,ngy-1,ngz-1)]=pdens[ix(offset,ngx-1,0,ngz-1)]=pdens[ix(offset,ngx-1,ngy-1,0)]=pdens[ix(offset,0,0,ngz-1)]=pdens[ix(offset,0,ngy-1,0)]=pdens[ix(offset,ngx-1,0,0)]=pdens[ix(offset,0,0,0)];
#endif
	
#ifdef MPI
	mlti=1.0;
	for(ij=1; ij<ngx-1; ij++)
		for(ik=1; ik<ngy-1; ik++)
		{						
			rpdens[ix(offset,ij,ik,0)]=mlti*(rpdens[ix(offset,ij,ik,0)]+rpdens[ix(offset,ij,ik,ngz-1)]);
			rpdens[ix(offset,ij,ik,ngz-1)]=rpdens[ix(offset,ij,ik,0)];
		}
	for(ik=1; ik<ngy-1; ik++)
		for(il=1; il<ngz-1; il++)
		{
			rpdens[ix(offset,0,ik,il)]=mlti*(rpdens[ix(offset,0,ik,il)]+rpdens[ix(offset,ngx-1,ik,il)]);
			rpdens[ix(offset,ngx-1,ik,il)]=rpdens[ix(offset,0,ik,il)];
		}
	for(ij=1; ij<ngx-1; ij++)
		for(il=1; il<ngz-1; il++)
		{
			rpdens[ix(offset,ij,0,il)]=mlti*(rpdens[ix(offset,ij,0,il)]+rpdens[ix(offset,ij,ngy-1,il)]);
			rpdens[ix(offset,ij,ngy-1,il)]=rpdens[ix(offset,ij,0,il)];
		}
	
	//edges
	mlti=1.0;
	for(ij=1; ij<ngx-1; ij++)
	{
		rpdens[ix(offset,ij,0,0)]=mlti*(rpdens[ix(offset,ij,0,0)]+rpdens[ix(offset,ij,ngy-1,ngz-1)]+rpdens[ix(offset,ij,0,ngz-1)]+rpdens[ix(offset,ij,ngy-1,0)]);				
		rpdens[ix(offset,ij,0,ngz-1)]=rpdens[ix(offset,ij,ngy-1,0)]=rpdens[ix(offset,ij,ngy-1,ngz-1)]=rpdens[ix(offset,ij,0,0)];
	}
	for(ik=1; ik<ngy-1; ik++)
	{
		rpdens[ix(offset,0,ik,0)]=mlti*(rpdens[ix(offset,0,ik,0)]+rpdens[ix(offset,ngy-1,ik,ngz-1)]+rpdens[ix(offset,0,ik,ngz-1)]+rpdens[ix(offset,ngx-1,ik,0)]);				
		rpdens[ix(offset,0,ik,ngz-1)]=rpdens[ix(offset,ngx-1,ik,0)]=rpdens[ix(offset,ngx-1,ik,ngz-1)]=rpdens[ix(offset,0,ik,0)];
	}
	for(il=1; il<ngz-1; il++)
	{
		rpdens[ix(offset,0,0,il)]=mlti*(rpdens[ix(offset,0,0,il)]+rpdens[ix(offset,ngx-1,ngy-1,il)]+rpdens[ix(offset,ngx-1,0,il)]+rpdens[ix(offset,0,ngy-1,il)]);				
		rpdens[ix(offset,ngx-1,0,il)]=rpdens[ix(offset,0,ngy-1,il)]=rpdens[ix(offset,ngx-1,ngy-1,il)]=rpdens[ix(offset,0,0,il)];
	}
	mlti=1.0;	//corners
	rpdens[ix(offset,0,0,0)]=mlti*(rpdens[ix(offset,ngx-1,ngy-1,ngz-1)]+rpdens[ix(offset,0,ngy-1,ngz-1)]+rpdens[ix(offset,ngx-1,0,ngz-1)]+rpdens[ix(offset,ngx-1,ngy-1,0)]
								  +rpdens[ix(offset,0,0,ngz-1)]+rpdens[ix(offset,0,ngy-1,0)]+rpdens[ix(offset,ngx-1,0,0)]+rpdens[ix(offset,0,0,0)]);
	rpdens[ix(offset,ngx-1,ngy-1,ngz-1)]=rpdens[ix(offset,0,ngy-1,ngz-1)]=rpdens[ix(offset,ngx-1,0,ngz-1)]=rpdens[ix(offset,ngx-1,ngy-1,0)]=rpdens[ix(offset,0,0,ngz-1)]=rpdens[ix(offset,0,ngy-1,0)]=rpdens[ix(offset,ngx-1,0,0)]=rpdens[ix(offset,0,0,0)];
	
	
#endif	
	
	
#endif	
	
	
	
	//03.09.2013 End of modifications		
	
  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<ngx; i++)
  { 
   mltx=1;	
   if((i==0) || (i==(ngx-1))) mltx=1;
   for(j=0; j<ngy; j++)
     {
       mlty=1;
       if((j==0) || (j==(ngy-1))) mlty=1;
       for(k=0; k<ngz; k++)
	 {
	   mltz=1;
	   if((k==0) || (k==(ngz-1))) mltz=1;	
#ifdef MPI	
	   //   temp[i]=weight*qdens[i]*normqdens*mlt;
	   temp[ix(0,i,j,k)]=mltx*mlty*mltz*weight*ratio*rpdens[ix(offset,i,j,k)]/(normxcube);
#else
	   temp[ix(0,i,j,k)]=mltx*mlty*mltz*weight*ratio*pdens[ix(offset,i,j,k)]/(normxcube);
	 //  fprintf(fpointer, "%E\n", temp[ix(0,i,j,k)]);
	 //  printf("%E\n", temp[ix(0,i,j,k)]);
#endif		
	   //normxcube because of normalized dV in particle collection
	 }
     }
  }
  fwrite(temp, sizeof(double), qdensMAX, fpointer);
  //here qdens because of offset
  printf("finished printdensity\n"); //getchar();
  free_dvecmem(temp, 0, qdensMAX-1); 
}

/*************PRINT AVERAGE PARTICLE VELOCITY******************/
void print_avpvel(FILE *fpointer, int k, int t, double weight)
{
 printf("started print_avpvel\n");
  int i,offset; 
  double *tempx, *tempy, *tempz, time;
  //PRINT DENSITY
  //print time
  int size=pdens_off-1;
  tempx=dvecmem(0,size);
  tempy=dvecmem(0,size);
  tempz=dvecmem(0,size);
  offset=pdens_off*k;
  time=t*dt*normtime;
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
  //  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<pdens_off; i++)
    {        
#ifdef MPI
      tempx[i]=weight*normvel*rvxvec[offset+i]/rpdens[offset+i];
      tempy[i]=weight*normvel*rvyvec[offset+i]/rpdens[offset+i];
      tempz[i]=weight*normvel*rvzvec[offset+i]/rpdens[offset+i];
#else
      tempx[i]=weight*normvel*vxvec[offset+i]/pdens[offset+i];
      tempy[i]=weight*normvel*vyvec[offset+i]/pdens[offset+i];
      tempz[i]=weight*normvel*vzvec[offset+i]/pdens[offset+i];	
#endif
    }  
	printf("in print_avpvel\n");  
  fwrite(tempx, sizeof(double), qdensMAX, fpointer);
  fwrite(tempy, sizeof(double), qdensMAX, fpointer);
  fwrite(tempz, sizeof(double), qdensMAX, fpointer);

  free_dvecmem(tempx, 0, size); 
  free_dvecmem(tempy, 0, size);
  free_dvecmem(tempz, 0, size);
  printf("finished print_avpvel\n");

}

/*************PRINT POTENTIAL ************************/
void printpotential(FILE *fpointer, int t, double weight)
{
  int i,j,k;
  double *temp, time;
  //print time  
  time=t*dt*normtime;
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
// for(i=0; i<ngx; i++)
// for(j=0; j<ngy; j++)
// for(k=0; k<ngz; k++)
 // phi[ix(0,i,j,k)]=i;

    temp=dvecmem(0,phiMAX-1);
  for(i=0; i<phiMAX; i++)
    temp[i]=phi[i]*normpot;
  fwrite(temp, sizeof(double), phiMAX, fpointer);
  free_dvecmem(temp,0,phiMAX-1); 
}

/*************PRINT POTENTIAL AV *********************/
void printavpotential(FILE *fpointer, int t, double weight)
{
  int i;
  double *temp, time;
  //print time  
  time=t*dt*normtime;
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
  temp=dvecmem(0,phiMAX-1);  
  for(i=0; i<phiMAX; i++)
    temp[i]=phiav[i]*normpot*weight;
  fwrite(temp, sizeof(double), phiMAX, fpointer);
  free_dvecmem(temp,0,phiMAX-1);
}
/*************PRINT ELECTRIC FIELD********************/
void printefield(FILE *fpointer, int help, int t, double weight)
{
  int i, offset;
  double *temp, time;
  //print time  
  time=t*dt*normtime;
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
  temp=dvecmem(0,phiMAX-1);
  if(help==0)
    offset=0;	
  if(help==1)
    offset=FsEy;
  if(help==2)
    offset=FsEz;	
  for(i=0; i<phiMAX; i++)
    temp[i]=Fs[offset+i]*normEfield*dVdt;
  //additional factor dVdt is due to the electric_field(), in which we multiply the E field by dt/dV
  fwrite(temp, sizeof(double), phiMAX, fpointer);
  free_dvecmem(temp,0,phiMAX-1);  
}

/***************FUNCTION called in main()*********************/
void printKEall(int t)
{
  int i;
  if((numtasks>1 && rank==1)||(numtasks==1 && rank==0))
    {
      //diagn.c - print KE
      printKE(eke, eke_time, 0, t, weight); 
      printKE(ike, ike_time, 1, t, weight); 
    }
  for(i=0;i<KEMAX; i++)
    KE[i]=0.0;
}

/**************PRINT KINETC ENERGY***********************/
/*KE must be callculated after acceleration, but before moving*/
void printKE(FILE *fpointer, FILE *fpointer2, int specie, int t, double weight)
{
  double totalke;
  int i,j,k,mltx, mlty, mltz,offset; 
  double *temp,time;
  double normxcube=normx*normx*normx;
  offset=KE_off*specie;
  totalke=0.0; 
  //PRINT DENSITY
  //print time
  time=t*dt*normtime;
  //  fprintf(fpointer, "%E\n", t*dt*normtime);
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<ngx; i++)
    { 
      mltx=1;	
      if((i==0) || (i==(ngx-1))) mltx=2;
      for(j=0; j<ngy; j++)
	{
	  mlty=1;
	  if((j==0) || (j==(ngy-1))) mlty=2;
	  for(k=0; k<ngz; k++)
	    {
	      mltz=1;
	      if((k==0) || (k==(ngz-1))) mltz=2;	
#ifdef MPI	
	      temp[ix(0,i,j,k)]=mltx*mlty*mltz*weight*0.25*mass[specie]*rKE[ix(offset,i,j,k)]*normvel*normvel/(dV)/cellvolume; 
#else 
	      temp[ix(0,i,j,k)]=mltx*mlty*mltz*weight*0.25*mass[specie]*KE[ix(offset,i,j,k)]*normvel*normvel/(dV)/cellvolume; 
#endif
	      //  temp[i]=weight*qdens[i]*normqdens*mlt;
	      //	temp[ix(0,i,j,k)]=mlt*weight*ratio*rpdens[ix(offset,i,j,k)]/(normxcube);	
	      //normxcube because of normalized dV in particle collection
	      totalke+=temp[i]; 
	    }
	}
    }
  fwrite(temp, sizeof(double), qdensMAX, fpointer);
  //here qdens because of offset
  free_dvecmem(temp,0,qdensMAX-1);
  fprintf(fpointer2,"%E\t%E\n",t*dt*normtime,totalke*cellvolume);
}

/**************PRINT POTENTIAL ENERGY ************************/
void printPE(FILE *fpointer,FILE *fpointer2, int specie, int t, double weight)
{   
  double totalpe, *temp, time;
  int i,j,k,mltx, mlty, mltz,offset; 
  double normxcube=normx*normx*normx;
  offset=PEMAXhalf*specie;
  totalpe=0.0; 
  //PRINT DENSITY
  //print time
  //  fprintf(fpointer, "%E\n", t*dt*normtime);
  time=t*dt*normtime;
  fwrite(&time, sizeof(double), 1, fpointer);
  //print grid quantities
  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<qdensMAX; i++)
    { 
      //I am not sure about multipliers here
      temp[i]=weight*charge[specie]*PE[offset+i]*normpot/(normxcube);
      //  temp[i]=weight*qdens[i]*normqdens*mlt;
      //	temp[ix(0,i,j,k)]=mlt*weight*ratio*rpdens[ix(offset,i,j,k)]/(normxcube);	
      //normxcube because of normalized dV in particle collection
      totalpe+=temp[i]; 
    }   
  fwrite(temp, sizeof(double), qdensMAX, fpointer);
  free_dvecmem(temp,0,qdensMAX-1);
  fprintf(fpointer2,"%E\t%E\n",t*dt*normtime,totalpe*cellvolume);
}  

/**************PRINT POTENTIAL ENERGY ************************/
void printPEtotal(FILE *fpointer, FILE *fpointer2, int t, double weight)
{
  double totalpe;
  double *temp, time; 
  int i;   
  totalpe=0.0; 
  //print time
  time=t*dt*normtime;
  //fprintf(fpointer, "%E\n", t*dt*normtime);
  fwrite(&time, sizeof(double), 1, fpointer);
  temp=dvecmem(0,qdensMAX-1);
  for(i=0; i<qdensMAX; i++)
   { 
   temp[i]=weight*PEtotal[i]*normPE;
   totalpe+=temp[i]; 
	}   
  fwrite(temp, sizeof(double), qdensMAX, fpointer);
  free_dvecmem(temp,0,qdensMAX-1);
  //it is charge density * potential, both normalized
  fprintf(fpointer2,"%E\t%E\n",t*dt*normtime,totalpe*cellvolume);
}


/*************PRINT PARTICLE DENSITY******************/
void printnewprobe(FILE *fpointer, int t, double weight)
{
	int i,j; 
	int alldata=0;
	for(i=0; i<noofdusts; i++)
		alldata+=dpartlast[i];
	//PRINT DENSITY
	//print time
	if(t==0)
		fprintf(fpointer,"%d\n", alldata);
	fprintf(fpointer, "%E\n", t*dt*normtime);
	//print grid density
#ifdef MPI  
	// MPI_Reduce(rdpart,dpart,alldata,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	for(i=0; i<noofdusts; i++)
		for(j=0; j<dpartlast[i]; j++)
		{	
			//we are on rank 0
			fprintf(fpointer,"%E\t%E\t%E\t%E\t%E\t%E\n", rdpart[i][j].x*normx, rdpart[i][j].y*normx, rdpart[i][j].z*normx, rdpart[i][j].ecnt, rdpart[i][j].icnt, rdpart[i][j].bcnt);
		}
#else  
	for(i=0; i<noofdusts; i++)
		for(j=0; j<dpartlast[i]; j++)
			fprintf(fpointer,"%E\t%E\t%E\t%E\t%E\t%E\n", dpart[i][j].x*normx, dpart[i][j].y*normx, dpart[i][j].z*normx, dpart[i][j].ecnt, dpart[i][j].icnt, dpart[i][j].bcnt);      
#endif
} 

/*************PRINT PARTICLE DENSITY******************/
void printdustcharge(FILE *fpointer, int t, double weight)
{
  int i,j; 
  int alldata=0;
  for(i=0; i<noofdusts; i++)
    alldata+=dpartlast[i];
  //PRINT DENSITY
  //print time
  if(t==0)
    fprintf(fpointer,"%d\n", alldata);
  fprintf(fpointer, "%E\n", t*dt*normtime);
  //print grid density
#ifdef MPI  
  // MPI_Reduce(rdpart,dpart,alldata,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for(i=0; i<noofdusts; i++)
    for(j=0; j<dpartlast[i]; j++)
      {	
	//we are on rank 0
	fprintf(fpointer,"%E\t%E\t%E\t%E\n", rdpart[i][j].x*normx, rdpart[i][j].y*normx, rdpart[i][j].z*normx, rdpart[i][j].q*normcharge);
      }
#else  
  for(i=0; i<noofdusts; i++)
    for(j=0; j<dpartlast[i]; j++)
      fprintf(fpointer,"%E\t%E\t%E\t%E\n", dpart[i][j].x*normx, dpart[i][j].y*normx, dpart[i][j].z*normx, dpart[i][j].q*normcharge);      
#endif
} 


/**********dust charge vs time ***************/
void printdustchargetime(FILE *fpointer, int t, double weight)
{
  int i,j; 
  double total=0.0;
  //PRINT DENSITY
  //print time
  fprintf(fpointer, "%E\t", t*dt*normtime);
  //print grid density
#ifdef MPI
  int alldata=0;
  for(i=0; i<noofdusts; i++)
    alldata+=dpartlast[i];
  //reduced already
  //  MPI_Reduce(rdpart,dpart,alldata,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for(i=0; i<noofdusts; i++)
    {
      total=0;
      for(j=0; j<dpartlast[i]; j++)
	total+=rdpart[i][j].q;
      fprintf(fpointer, "%E\t", total*normcharge);
    }
  fprintf(fpointer, "\n");
#else  
  for(i=0; i<noofdusts; i++)
    {
      total=0.0;
      for(j=0; j<dpartlast[i]; j++)
	total+=dpart[i][j].q;
      fprintf(fpointer, "%E\t", total*normcharge);
    }
  fprintf(fpointer, "\n");	      
  printf("DIAGN: printdustchargetime()\n");
  fclose(fpointer);
   fpointer=my_file_open("./data/dustcharge_time.dat", "a");

#endif      
}


//3DUPTOHERE
/**********dust potential history vs time ***************/
void printdth(FILE *fpointer, int t)
{
  int i,j; 
  double totalpoten=0.0;
  int noofgp=0;
  //PRINT DENSITY
  //  fprintf(fpointer,"%E\n", t*dt*normtime);
  
  for(i=0; i<phiMAX; i++)
    {	
#ifdef MPI
      //	if((j==0) || (j==(ngy-1))) rrho[ix(rhoMAXhalf,i,j)]*=2;
      //if((i==0) || (i==(ngx-1))) rrho[ix(rhoMAXhalf,i,j)]*=2;
      //	total+= rrho[ix(rhoMAXhalf,i,j)];
      //if(gp[i][j].marker>BOUND)
      /*{
	totalpoten+=phi[i];
	noofgp++;
      }
      */
      // fprintf(fpointer, "%d\t%E\t%d\t%E\t%E\t%E\n",i,dx*normx,j,dy*normx,gp[i][j].phi*normpot,rrho[ix(rhoMAXhalf,i,j)]*normqdens);
#else 
      //	if((j==0) || (j==(ngy-1))) rho[ix(rhoMAXhalf,i,j)]*=2;
      //if((i==0) || (i==(ngx-1))) rho[ix(rhoMAXhalf,i,j)]*=2;
      //	total+= rho[ix(rhoMAXhalf,i,j)];

      /* if(gp[i][j].marker>BOUND)
	{
	  totalpoten+=gp[i][j].phi;
	  noofgp++;
	  }*/

      //fprintf(fpointer, "%d\t%E\t%d\t%E\t%E\t%E\n",i,dx*normx,j,dy*normx,gp[i][j].phi*normpot,rho[ix(rhoMAXhalf,i,j)]*normqdens);
#endif
    } 
  /*
  fprintf(fpointer, "%E\t%E\n", t*dt*normtime, totalpoten*normpot/noofgp);
  printf("DIAGN: printdth()\n");*/
}

/********PRINT PROBESHAPE **************************/
void printdustshape()
{
  //static
  /*
  int i,j,allcorners;
  allcorners=0;
  for(i=0; i<noofdusts; i++)
      allcorners+=ncorners[i];

  fprintf(dustshape, "%d\n", allcorners+noofdusts-1);

  for(i=0; i<noofdusts; i++)
    {
	for(j=0; j<ncorners[i]; j++)
	fprintf(dustshape,"%E\t%E\n",dustx[i][j]*dx*normx, dusty[i][j]*dy*normx);
      fprintf(dustshape,"%E\t%E koniec\n",dustx[i][0]*dx*normx, dusty[i][0]*dy*normx);

    }
	*/
}

/********PRINT PROBESHAPE **************************/
void printdustshapetime(int t)
{
  if(noofdusts>0)
    {
      //moving
      int i,j,allcorners;
      allcorners=0;
      for(i=0; i<noofdusts; i++)
	allcorners+=ncorners[i];
      
      fprintf(dustshapet, "%E\n", t*normtime*dt);
      fprintf(dustshapet, "%d\n", allcorners+noofdusts-1);
      
      for(i=0; i<noofdusts; i++)
	{
	  for(j=0; j<ncorners[i]; j++)
	    fprintf(dustshapet,"%E\t%E\n",dustx[i][j]*dx*normx, dusty[i][j]*dy*normx);
	  fprintf(dustshapet,"%E\t%E\n",dustx[i][0]*dx*normx, dusty[i][0]*dy*normx);
	}
    }
}
/*********potential probes**************************/
void pot_probes(int t)
{
  int i1,i2,i3,j,k;
  double time;
  double *temp1, *temp2, *temp3;
  //, *temp4, *temp5, *temp6, *temp7, *temp8, *temp9;
  int size; 
  int centr=(int)((ngy-1)/2);
  int c11, c12, c13, c21, c22, c23;

  c12=c22=centr;
  c11=c21=(int)(centr/2);
  c13=c23=(int)(centr+centr/2);

  if(rank==0){
      probes11=my_file_open("./data/probes11.dat", "a");
      probes12=my_file_open("./data/probes12.dat", "a");	
      probes13=my_file_open("./data/probes13.dat", "a");
      probes14=my_file_open("./data/probes14.dat", "a");
      probes15=my_file_open("./data/probes15.dat", "a");	
      probes16=my_file_open("./data/probes16.dat", "a");   
      probes17=my_file_open("./data/probes17.dat", "a");
      probes18=my_file_open("./data/probes18.dat", "a");	
      probes19=my_file_open("./data/probes19.dat", "a");

      probes21=my_file_open("./data/probes21.dat", "a");
      probes22=my_file_open("./data/probes22.dat", "a");	
      probes23=my_file_open("./data/probes23.dat", "a");
      probes24=my_file_open("./data/probes24.dat", "a");
      probes25=my_file_open("./data/probes25.dat", "a");	
      probes26=my_file_open("./data/probes26.dat", "a");   
      probes27=my_file_open("./data/probes27.dat", "a");
      probes28=my_file_open("./data/probes28.dat", "a");	
      probes29=my_file_open("./data/probes29.dat", "a");

      probes31=my_file_open("./data/probes31.dat", "a");
      probes32=my_file_open("./data/probes32.dat", "a");	
      probes33=my_file_open("./data/probes33.dat", "a");
      probes34=my_file_open("./data/probes34.dat", "a");
      probes35=my_file_open("./data/probes35.dat", "a");	
      probes36=my_file_open("./data/probes36.dat", "a");   
      probes37=my_file_open("./data/probes37.dat", "a");
      probes38=my_file_open("./data/probes38.dat", "a");	
      probes39=my_file_open("./data/probes39.dat", "a");
  }
 if(rank==0){
  i1=centr;
  i2=centr;
  i3=centr;
  size=ngy;
  //print time
  timerprobes++;
  time=t*dt*normtime;

  fwrite(&time, sizeof(double), 1, probes11);
  fwrite(&time, sizeof(double), 1, probes12);
  fwrite(&time, sizeof(double), 1, probes13);     
  fwrite(&time, sizeof(double), 1, probes14);
  fwrite(&time, sizeof(double), 1, probes15);
  fwrite(&time, sizeof(double), 1, probes16);  
  fwrite(&time, sizeof(double), 1, probes17);
  fwrite(&time, sizeof(double), 1, probes18);
  fwrite(&time, sizeof(double), 1, probes19);  

  fwrite(&time, sizeof(double), 1, probes21);
  fwrite(&time, sizeof(double), 1, probes22);
  fwrite(&time, sizeof(double), 1, probes23);     
  fwrite(&time, sizeof(double), 1, probes24);
  fwrite(&time, sizeof(double), 1, probes25);
  fwrite(&time, sizeof(double), 1, probes26);  
  fwrite(&time, sizeof(double), 1, probes27);
  fwrite(&time, sizeof(double), 1, probes28);
  fwrite(&time, sizeof(double), 1, probes29);  
   
  fwrite(&time, sizeof(double), 1, probes31);
  fwrite(&time, sizeof(double), 1, probes32);
  fwrite(&time, sizeof(double), 1, probes33);     
  fwrite(&time, sizeof(double), 1, probes34);
  fwrite(&time, sizeof(double), 1, probes35);
  fwrite(&time, sizeof(double), 1, probes36);  
  fwrite(&time, sizeof(double), 1, probes37);
  fwrite(&time, sizeof(double), 1, probes38);
  fwrite(&time, sizeof(double), 1, probes39);  

  //print potential - grid quatnities

  //along x
  for(j=0; j<ngy; j++)
    //    for(k=0;k<ngz; k++)
      {
	//along x
	ptemp11[j]=phi[ix(0,j,c11,c21)]*normpot;     
	ptemp12[j]=phi[ix(0,j,c12,c21)]*normpot; 
	ptemp13[j]=phi[ix(0,j,c13,c21)]*normpot; 
	ptemp14[j]=phi[ix(0,j,c11,c22)]*normpot;     
	ptemp15[j]=phi[ix(0,j,c12,c22)]*normpot; 
	ptemp16[j]=phi[ix(0,j,c13,c22)]*normpot; 
	ptemp17[j]=phi[ix(0,j,c11,c23)]*normpot;     
	ptemp18[j]=phi[ix(0,j,c12,c23)]*normpot; 
	ptemp19[j]=phi[ix(0,j,c13,c23)]*normpot; 
	//along y
	ptemp21[j]=phi[ix(0,c11,j,c21)]*normpot;     
	ptemp22[j]=phi[ix(0,c12,j,c21)]*normpot; 
	ptemp23[j]=phi[ix(0,c13,j,c21)]*normpot; 
	ptemp24[j]=phi[ix(0,c11,j,c22)]*normpot;     
	ptemp25[j]=phi[ix(0,c12,j,c22)]*normpot; 
	ptemp26[j]=phi[ix(0,c13,j,c22)]*normpot; 
	ptemp27[j]=phi[ix(0,c11,j,c23)]*normpot;     
	ptemp28[j]=phi[ix(0,c12,j,c23)]*normpot; 
	ptemp29[j]=phi[ix(0,c13,j,c23)]*normpot; 
	//along z
	ptemp31[j]=phi[ix(0,c11,c21,j)]*normpot;     
	ptemp32[j]=phi[ix(0,c12,c21,j)]*normpot; 
	ptemp33[j]=phi[ix(0,c13,c21,j)]*normpot; 
	ptemp34[j]=phi[ix(0,c11,c22,j)]*normpot;     
	ptemp35[j]=phi[ix(0,c12,c22,j)]*normpot; 
	ptemp36[j]=phi[ix(0,c13,c22,j)]*normpot; 
	ptemp37[j]=phi[ix(0,c11,c23,j)]*normpot;     
	ptemp38[j]=phi[ix(0,c12,c23,j)]*normpot; 
	ptemp39[j]=phi[ix(0,c13,c23,j)]*normpot; 

      }
  fwrite(ptemp11, sizeof(double), size, probes11);
  fwrite(ptemp12, sizeof(double), size, probes12);
  fwrite(ptemp13, sizeof(double), size, probes13);  
  fwrite(ptemp14, sizeof(double), size, probes14);
  fwrite(ptemp15, sizeof(double), size, probes15);
  fwrite(ptemp16, sizeof(double), size, probes16);  
  fwrite(ptemp17, sizeof(double), size, probes17);
  fwrite(ptemp18, sizeof(double), size, probes18);
  fwrite(ptemp19, sizeof(double), size, probes19);

  fwrite(ptemp21, sizeof(double), size, probes21);
  fwrite(ptemp22, sizeof(double), size, probes22);
  fwrite(ptemp23, sizeof(double), size, probes23);  
  fwrite(ptemp24, sizeof(double), size, probes24);
  fwrite(ptemp25, sizeof(double), size, probes25);
  fwrite(ptemp26, sizeof(double), size, probes26);  
  fwrite(ptemp27, sizeof(double), size, probes27);
  fwrite(ptemp28, sizeof(double), size, probes28);
  fwrite(ptemp29, sizeof(double), size, probes29);

  fwrite(ptemp31, sizeof(double), size, probes31);
  fwrite(ptemp32, sizeof(double), size, probes32);
  fwrite(ptemp33, sizeof(double), size, probes33);  
  fwrite(ptemp34, sizeof(double), size, probes34);
  fwrite(ptemp35, sizeof(double), size, probes35);
  fwrite(ptemp36, sizeof(double), size, probes36);  
  fwrite(ptemp37, sizeof(double), size, probes37);
  fwrite(ptemp38, sizeof(double), size, probes38);
  fwrite(ptemp39, sizeof(double), size, probes39);
  
  //   free_dvecmem(temp1, 0, size);
  // free_dvecmem(temp2, 0, size); 
  // free_dvecmem(temp3, 0, size);
 
 }
  if(rank==0){
    fclose(probes11);
    fclose(probes12);
    fclose(probes13);
    fclose(probes14);
    fclose(probes15);
    fclose(probes16);
    fclose(probes17);
    fclose(probes18);
    fclose(probes19);

    fclose(probes21);
    fclose(probes22);
    fclose(probes23);
    fclose(probes24);
    fclose(probes25);
    fclose(probes26);
    fclose(probes27);
    fclose(probes28);
    fclose(probes29);

    fclose(probes31);
    fclose(probes32);
    fclose(probes33);
    fclose(probes34);
    fclose(probes35);
    fclose(probes36);
    fclose(probes37);
    fclose(probes38);
    fclose(probes39);

  }
}

/********probes for the potential************/
void pot_probes_init(void)
{
  int i1,i2,i3;
  double temp;
  int centr=(int)((ngy-1)/2);

  /*open files for potential probe arrays*/
  if(rank==0){
      probes11=my_file_open("./data/probes11.dat", "a");
      probes12=my_file_open("./data/probes12.dat", "a");	
      probes13=my_file_open("./data/probes13.dat", "a");
      probes14=my_file_open("./data/probes14.dat", "a");
      probes15=my_file_open("./data/probes15.dat", "a");	
      probes16=my_file_open("./data/probes16.dat", "a");
      probes17=my_file_open("./data/probes17.dat", "a");
      probes18=my_file_open("./data/probes18.dat", "a");	
      probes19=my_file_open("./data/probes19.dat", "a");

      probes21=my_file_open("./data/probes21.dat", "a");
      probes22=my_file_open("./data/probes22.dat", "a");	
      probes23=my_file_open("./data/probes23.dat", "a");
      probes24=my_file_open("./data/probes24.dat", "a");
      probes25=my_file_open("./data/probes25.dat", "a");	
      probes26=my_file_open("./data/probes26.dat", "a");
      probes27=my_file_open("./data/probes27.dat", "a");
      probes28=my_file_open("./data/probes28.dat", "a");	
      probes29=my_file_open("./data/probes29.dat", "a");

      probes31=my_file_open("./data/probes31.dat", "a");
      probes32=my_file_open("./data/probes32.dat", "a");	
      probes33=my_file_open("./data/probes33.dat", "a");
      probes34=my_file_open("./data/probes34.dat", "a");
      probes35=my_file_open("./data/probes35.dat", "a");	
      probes36=my_file_open("./data/probes36.dat", "a");
      probes37=my_file_open("./data/probes37.dat", "a");
      probes38=my_file_open("./data/probes38.dat", "a");	
      probes39=my_file_open("./data/probes39.dat", "a");

  }
  if(rank==0){
  i1=centr;
  i2=centr;
  i3=centr;
  int size=ngy;
  ptemp11=dvecmem(0,size);
  ptemp12=dvecmem(0,size);
  ptemp13=dvecmem(0,size);
  ptemp14=dvecmem(0,size);
  ptemp15=dvecmem(0,size);
  ptemp16=dvecmem(0,size);
  ptemp17=dvecmem(0,size);
  ptemp18=dvecmem(0,size);
  ptemp19=dvecmem(0,size);

  ptemp21=dvecmem(0,size);
  ptemp22=dvecmem(0,size);
  ptemp23=dvecmem(0,size);
  ptemp24=dvecmem(0,size);
  ptemp25=dvecmem(0,size);
  ptemp26=dvecmem(0,size);
  ptemp27=dvecmem(0,size);
  ptemp28=dvecmem(0,size);
  ptemp29=dvecmem(0,size);


  ptemp31=dvecmem(0,size);
  ptemp32=dvecmem(0,size);
  ptemp33=dvecmem(0,size);
  ptemp34=dvecmem(0,size);
  ptemp35=dvecmem(0,size);
  ptemp36=dvecmem(0,size);
  ptemp37=dvecmem(0,size);
  ptemp38=dvecmem(0,size);
  ptemp39=dvecmem(0,size);


  timerprobes=0;
  //print position of the slice
  temp=i1*dx*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);

  temp=i2*dx*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);
  temp=i3*dx*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);
  //print grid spacing, dx, dy, dz
  temp=dx*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);


  temp=dy*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);


  temp=dz*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);


  //print system dimensions, Lx, Ly, Lz
  temp=Lx*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);


  temp=Ly*normx;
  fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);


  temp=Lz*normx;
   fwrite(&temp, sizeof(double), 1, probes11);
  fwrite(&temp, sizeof(double), 1, probes12);
  fwrite(&temp, sizeof(double), 1, probes13);
  fwrite(&temp, sizeof(double), 1, probes14);
  fwrite(&temp, sizeof(double), 1, probes15);
  fwrite(&temp, sizeof(double), 1, probes16);
  fwrite(&temp, sizeof(double), 1, probes17);
  fwrite(&temp, sizeof(double), 1, probes18);
  fwrite(&temp, sizeof(double), 1, probes19);

  fwrite(&temp, sizeof(double), 1, probes21);
  fwrite(&temp, sizeof(double), 1, probes22);
  fwrite(&temp, sizeof(double), 1, probes23);
  fwrite(&temp, sizeof(double), 1, probes24);
  fwrite(&temp, sizeof(double), 1, probes25);
  fwrite(&temp, sizeof(double), 1, probes26);
  fwrite(&temp, sizeof(double), 1, probes27);
  fwrite(&temp, sizeof(double), 1, probes28);
  fwrite(&temp, sizeof(double), 1, probes29);

  fwrite(&temp, sizeof(double), 1, probes31);
  fwrite(&temp, sizeof(double), 1, probes32);
  fwrite(&temp, sizeof(double), 1, probes33);
  fwrite(&temp, sizeof(double), 1, probes34);
  fwrite(&temp, sizeof(double), 1, probes35);
  fwrite(&temp, sizeof(double), 1, probes36);
  fwrite(&temp, sizeof(double), 1, probes37);
  fwrite(&temp, sizeof(double), 1, probes38);
  fwrite(&temp, sizeof(double), 1, probes39);
  }
  if(rank==0){
    fclose(probes11);
    fclose(probes12);
    fclose(probes13);
    fclose(probes14);
    fclose(probes15);
    fclose(probes16);
    fclose(probes17);
    fclose(probes18);
    fclose(probes19);

    fclose(probes21);
    fclose(probes22);
    fclose(probes23);
    fclose(probes24);
    fclose(probes25);
    fclose(probes26);
    fclose(probes27);
    fclose(probes28);
    fclose(probes29);

    fclose(probes31);
    fclose(probes32);
    fclose(probes33);
    fclose(probes34);
    fclose(probes35);
    fclose(probes36);
    fclose(probes37);
    fclose(probes38);
    fclose(probes39);
  }

}


void print_current(int tid)
{
 	long int temp;
	int dno, sp;
if(rank==0)
{
	curr=my_file_open("./data/current.txt", "a");
	fprintf(curr, "%E\t", tid*dt*normtime);
//	printf("%E\n", tid*dt*normtime);
}

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	for(dno=0; dno<noofdusts; dno++)
		for(sp=0; sp<S; sp++)
		{
			MPI_Reduce(&current[dno][sp],&temp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
			if(rank==0)			
			  {
			    fprintf(curr, "%d\t", temp);
			    printf("reduced current %d %d %d\n", dno, sp, temp);
			  }			
		}
#else	
      		for(dno=0; dno<noofdusts; dno++)
		  for(sp=0; sp<S; sp++)
		    {
		      if(rank==0)
			{
			  fprintf(curr, "%d\t", current[dno][sp]);
			  printf("current %d %d %d\n", dno, sp, current[dno][sp]);
			}
		    }				
#endif	
	
if(rank==0)
{
//	printf("closing\n");
	fprintf(curr, "\n");
	fclose(curr);	
}
}


void printpotdistribution(int t)
{
	int i,j;
	double ptdmax, distrspacing;
	for(j=0; j<POTPOTS; j++) potdistrarray[j]=0;
	if(rank==0){
		potdistrmin=0;
		potdistrmax=0;
		ptdmax=distrspacing=0;
		potdistr=my_file_open("./data/potdistribution.dat", "a");
		for(i=0; i<(ngx*ngy*ngz); i++)
		{
				if(phi[i] <= potdistrmin) potdistrmin=phi[i];
				if(phi[i] > potdistrmax) potdistrmax=phi[i];
		}
		//create symmetric distribution
	//	printf("MAXMIN %E %E \n", potdistrmin*normpot, potdistrmax*normpot); getchar();
		if(potdistrmax > -potdistrmin)
			ptdmax=potdistrmax;
		else	
			ptdmax=-potdistrmin;
		if(ptdmax< 0) ptdmax*=-1;
		//find spacing
		distrspacing=2*ptdmax/POTPOTS;

		fprintf(potdistr, "%E\t%E\t%E\t%E\t%E\n", t*dt*normtime, potdistrmin*normpot, potdistrmax*normpot, -ptdmax*normpot, ptdmax*normpot);
		//create list
		for(i=0; i<(ngx*ngy*ngz); i++)	
			for(j=0; j<POTPOTS; j++)
				if((phi[i] >= (-ptdmax+distrspacing*j)) && (phi[i] < (-ptdmax+distrspacing*(j+1))))
					potdistrarray[j]++;
		//write list	
		int checkup=0;
		for(j=0; j<POTPOTS; j++)
		{
			fprintf(potdistr, "%d\t", potdistrarray[j]);
			checkup+=potdistrarray[j];
		//	printf("potdistr %d\t", potdistrarray[j]);

		}
		fprintf(potdistr, "\n");
		
		//printf("ptdmax %E, distrspacing %E checkup %d\n", ptdmax*normpot, distrspacing*normpot, checkup); getchar();
		fclose(potdistr);
	}
}
