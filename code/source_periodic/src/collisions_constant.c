/*DiP3D*/
/*Collisions.*/
/*Author: Wojciech Jacek Miloch*/
/*University of Oslo, Norway*/
/*2009*/
/*Last revised 08.03.10*/
/*Edited by Miroslav Horky 22.08.2013*/

#include <math.h>
#include "const.h"

//in const.h
/*
double nullcoll[S];
double nullcollfreq[S];
double massneutrals;
double vthneutr;
double ndesnity;
double * sigma;
*/

void collisions_init(void)
{
  double vmax, ndensitymax, potmax;
  double check, v, step, kener; 
  int i, ii, ii_max;
  double tempneutr;
  double fs, checksigma;
  potmax=Vpr_begin>Vpr_end ? Vpr_begin : Vpr_end;

  colltypes[0]=1;
  colltypes[1]=2;
#ifdef BEAM
  colltypes[2]=2;
  colltypes[3]=1;
#endif
  sigma=(double *)malloc(2*sizeof(double));

  ndensity=2E1;
    //22; //in m-3
  //kT=mV2, T=400, 1eV=11605.505 K
  tempneutr=0/11605.505;
  //in finite temp it was 30/11605.505
  massneutrals=mass[1]/ratio; //real mass
  vthneutr=sqrt(tempneutr*Q/massneutrals); //note factor 2 is missing
  vthneutr/=normvel; //normalized velocity
  //  printf("vthneutr %E", vthneutr*normvel);
  //  getchar();

  //homogeneous neutral density
  ndensitymax=ndensity;
  for(i=0; i<S; i++)
    {
      nullcollrest[i]=0.0;
      checksigma=0;
      fs=0;
      //frequency: velocity*sigma(Ek)*ndensity(local)
      //vmax in m/s
      vmax=sqrt(vdriftx[i]*vdriftx[i]+4*vthx[i]*vthx[i]+fabs(qm[i]*2*Vpr_begin))*normvel+2*vthneutr*normvel; //two accounts for moving neutrals
      nullcollfreq[i]=0;
      v=0; 
      ii_max=100;   
   step=1.0*vmax/ii_max;

   //printf("vmax %E\n", vmax);
   //getchar();
     for(ii=0; ii<ii_max; ii++)
	{
	  v+=step;
	  kener=0.5*mass[i]*v*v/(ratio*Q);
	  check=v*findsigma(i,kener)*ndensitymax;
	  printf("v %E enrer %E sigma %E\n", v, kener, findsigma(i,kener));
	  if(check > nullcollfreq[i])
	    nullcollfreq[i]=check;
	  checksigma=findsigma(i,kener);
	  if(checksigma > fs)
	    fs=checksigma;
	  printf("%E\n", fs);
	}
     if(i==0 || i==4) nullcollfreq[i]=0.0;
     if(i==1 || i==2) nullcollfreq[i]=0.0;
		
     nullcoll[i]=1-exp(-dt*normtime*nullcollfreq[i]);
     printf("nullcollisions; specie %d: fraction of collided particles %E freq %E\n normtime %E coll pr dt %E\n", i, nullcoll[i], nullcollfreq[i], normtime, nullcoll[i]*npart[i]);
     printf("MEAN FREE PATH FOR (old algorithm - not reliable) %d: %E m\n", i, 1/(fs*ndensity));
      getchar();
    }
}


void collisions(void)
{
  int i, j, ii, noofcoll, segment, type, colltype;
  double rr, rr2;
  double kener, pvel, collfreq, before, after, collfreqi0;
  double vvxn,vvyn,vvzn;
  double r2, lx,ly,lz;
   int  kii=313104;
  for(i=0; i<S; i++)
    {

   
	   {
	     printf("inside coll 3) 313104 x,y,z %E %E %E\n",spec[0].part[kii].vx, spec[0].part[kii].vy, spec[0].part[kii].vz);
	   }
	   noofcoll=(int)(nullcoll[i]*npart[i]);
	   nullcollrest[i]+=(nullcoll[i]*npart[i])-((int)(nullcoll[i]*npart[i]));
	   if(nullcollrest[i]>=1)
	     {
	       noofcoll++;
	       nullcollrest[i]-=1.0;
	     }
	   
       printf("collisions specie %d noofcoll %d null coll %E npart %d\n", i, noofcoll, nullcoll[i]*npart[i], npart[i]);
      // getchar();
      if(noofcoll!=0)
	segment=npart[i]/noofcoll;
      //choose particles for collisions
      for(ii=0; ii<noofcoll; ii++)
	{
	  //choosing the particle
	  rr=drand48();
	  j=rr*segment+ii*segment;
	  //choosing another random numer for finding the collision type.
	  rr2=primeroot();
	  printf("j %d segment %d", j, segment);
	  //calculating collision frequencies
	  //changing the frame of the reference for ions
	  if(i==1 || i==2)
	    {
	      //finding velocitites of neutrals, BOX MULLER
	      do
		{
		  lx=primeroot()*2-1;
		  ly=primeroot()*2-1;	     
		  r2=lx*lx+ly*ly;
		}
	      while(r2>=1);	      
	      vvxn=vthneutr*lx*sqrt(-2*log(r2)/r2);
	      vvyn=vthneutr*ly*sqrt(-2*log(r2)/r2);
	      
	      /*Find the Maxwellian velocity normalized in z direction*/	  
	      do
		{
		  lx=primeroot()*2-1;
		  lz=primeroot()*2-1;	     
		  r2=lx*lx+lz*lz;
		}
	      while(r2>=1);	      
	      vvzn=vthneutr*lz*sqrt(-2*log(r2)/r2);
	      
	      spec[i].part[j].vx-=vvxn;
	      spec[i].part[j].vy-=vvyn;
	      spec[i].part[j].vz-=vvzn;
	      kener=spec[i].part[j].vx*spec[i].part[j].vx+spec[i].part[j].vy*spec[i].part[j].vy+spec[i].part[j].vz*spec[i].part[j].vz;	      
	    }
	  else
	    {
	   	
	      kener=spec[i].part[j].vx*spec[i].part[j].vx+spec[i].part[j].vy*spec[i].part[j].vy+spec[i].part[j].vz*spec[i].part[j].vz;
	      // kener=spec[i].part[j].kenergy;
 	      // kener/=2; //averaging over two time steps
	    }
	  //velocity in m/s
	  pvel=sqrt(kener)*normvel;
	  //kenergy in eV.
	  kener*=0.5*mass[i]*normvel*normvel/(Q*ratio);	  //mass[i]=realmass*ratio
	  //find sigmas 
	  //	  printf("kener %E vel %E\n", kener, pvel);
	  //findsigma(i, kener);      //commented out MH
	  before=0;
	  after=0;
	  colltype=-1;
	  //find collision type

	if(i==0)
		{
		colltype=0;
		}
	if(i==1)
		{
		//if(rr2 <= collfreqi0/nullcollfreq[1])
		//colltype=0;
		//else
		colltype=1;
		}

	/*	for(type=0; type<colltypes[i]; type++)
	    {	
	      //collfreq=pvel*sigma[type]*ndensity;      //commented out MH
		  after+=collfreq[type];
	      if((rr2 > before/nullcollfreq[i]) && rr2 <= after/nullcollfreq[i])
		colltype=type;
	      before=after;
	    }
	*/
		
//	if(rr2 < nullcollfreq[i])	
    //only collide ions
		//colltype=1;  //commented out MH
	 if(i==0) 
	  if(colltype==0)
		{
		collide(i, colltype, j, kener, pvel);
		//	printf("collision %d %d", i, colltype);
		}
	 if(i==1){
	  if(colltype==0)
		{
		collide(i, colltype, j, kener, pvel);
		//	printf("collision %d %d", i, colltype);
		}
	  if(colltype==1) // there is a collision	
	  {
	    collide(i, colltype, j, kener, pvel);
	    //	  printf("collision %d %d", i, colltype);
	  }}
	  //going back to the proper reference frame for ions
	  if(i==1 || i==2)
	    {
	      spec[i].part[j].vx+=vvxn;
	      spec[i].part[j].vy+=vvyn;
	      spec[i].part[j].vz+=vvzn;	      
	    }
	}
    }
    
	   {
	//     printf("after coll 3) 313104 x,y,z %E %E %E\n",spec[0].part[kii].vx, spec[0].part[kii].vy, spec[0].part[kii].vz);
	   }
}

double findsigma(int i, double kinener)
{
  double sigmamax=0;
  double coeff=3.536*1E20;
double coeffel=1E20;  
//renormalize energy in EV - kinener is renormalized now!
  //0.25 is due to 0.5 from averaging over two timesteps, 0.5 from 0.5mv2 
  //  kinener=0.25*mass[i]*normvel*normvel;
  // kinener/=Q;

  //electrons
 if((i==0) || (i==3))
   {
     if(kinener < 0.35)
       sigma[0]=(3.71*exp(-(kinener/0.07338))+0.26554)/coeffel;
     else
       {
//printf("above 0.35 el \n");
	sigma[0]=(2.29192-0.35989*kinener+0.42984*kinener*kinener-0.02716*kinener*kinener*kinener-0.000439607*kinener*kinener*kinener*kinener)/coeffel;
} 
     sigmamax=sigma[0];
    }
  //ions or beam
  if(i==1 || i==2)
    {
      if(kinener < 0.35)
       {
	 //elastic scattering
	 sigma[0]=(93.72278+66.15386*exp(-kinener/6.057))/(coeff);
	 //inelastic scattering
	 sigma[1]=(119.42+51.14658*exp(-kinener/5.0849))/(coeff);
	 sigmamax=sigma[0]+sigma[1];
       }
      else
	{  
	  //elastic scattering 
	  sigma[0]=(69+38.5*exp(-kinener/62.92)+74*(-kinener/3.51))/coeff;
	  //inelastic scattering
	  sigma[1]=(70.85+57.68*exp(-kinener/105))/coeff;

	  //printf("kinetic energy too large for ion collisions\n");
	  //	  getchar();
	}
    }
  return sigmamax;
}


void collide(int i, int colltype, int j, double kener, double pvel)
{

  double coschi, sinchi, newkener, r, phi, sintheta, sss, scs;
  double pvelold=pvel;
 double vvx,vvy,vvz;
  //electrons
 if((i==0) || (i==3))
   {
      //electron-neutral elastic collision, argon
      if(colltype==0)
	{
	  r=primeroot();
	  coschi=(2+kener-2*pow((1+kener),r))/(kener);
	  if(coschi > 1) {coschi=(coschi-((int)coschi)); printf("collide 1 %E\n", (2+kener-2*pow((1+kener),r))/(kener)); getchar();}
	  sinchi=sqrt(1-coschi*coschi);
	  phi=2*M_PI*primeroot();
	  newkener=kener*(1-(2*(mass[i]/ratio)*(1-coschi))/massneutrals);
	  
	  if(newkener <= 0) {newkener=-newkener; printf("collide 2 %E, check dE %E\n", -newkener, (2*(mass[i]/ratio)*(1-coschi))/massneutrals); getchar();}
	  if(j==313104)
	    printf("new kinetic part %d energy %E and old %E cehck %E pvel %E\n", j, newkener, kener, (2*(mass[i]/ratio)*(1-coschi))/massneutrals, pvel);
	  vvx=spec[i].part[j].vx*normvel/pvel;
	  vvy=spec[i].part[j].vy*normvel/pvel;
	  vvz=spec[i].part[j].vz*normvel/pvel;
	  pvel=sqrt(2*newkener*ratio*Q/mass[i]);
       	  if(j==313104)
	  printf("vvx, vvy, vvz: %E %E %E\n", vvx, vvy, vvz);
	  //1 - cos2, where cos=(vx,vy,vz)*(1,0,0)
	  if(vvx*vvx >= 1.0) {printf("collide 3 %E\n", vvx); vvx=0.9; getchar();}
	  sintheta=sqrt(1-vvx*vvx);

	  sss=sinchi*sin(phi)/sintheta;
	  scs=sinchi*cos(phi)/sintheta;
	   if(j==313104)
	  printf("parameters sintheta, sss, scs %E %E %E check %E\n\n\n\n", sintheta, coschi, sinchi,sqrt(1-coschi*coschi) );
 
	  spec[i].part[j].vx=(vvx*coschi+(vvy*vvy+vvz*vvz)*scs)*pvel/normvel;
	  spec[i].part[j].vy=(vvy*coschi+vvz*sss-vvx*vvy*scs)*pvel/normvel;
	  spec[i].part[j].vz=(vvz*coschi-vvy*sss-vvx*vvz*scs)*pvel/normvel;
	   if(j==313104)
	  printf("vx %E, vy %E, vz %E, check %E\n",  spec[i].part[j].vx,  spec[i].part[j].vy,  spec[i].part[j].vz , pvel/normvel/sqrt( spec[i].part[j].vx*spec[i].part[j].vx+spec[i].part[j].vy*spec[i].part[j].vy+spec[i].part[j].vz* spec[i].part[j].vz));
	  //  printf("parameter check %E\n\n\n",(vvx*coschi+(vvy*vvy+vvz*vvz)*scs)*(vvx*coschi+(vvy*vvy+vvz*vvz)*scs)+(vvy*coschi+vvz*sss-vvx*vvy*scs)*(vvy*coschi+vvz*sss-vvx*vvy*scs)+(vvz*coschi-vvy*sss-vvx*vvz*scs)*(vvz*coschi-vvy*sss-vvx*vvz*scs));
	 
	}
    } 
  
  
  if((i==1) || (i==2)) 
    {
      //ion-neutral collsions, argon
      if(colltype==0) //elastic
	{
	   if(mass[i]/ratio!=massneutrals) printf("wrong mass[1]!=massneutrals in scatternig\n");
	  
	  coschi=sqrt(1-primeroot());
	  if(coschi > 1) {printf("collide 5 %E\n", coschi); coschi=(coschi-((int)coschi)); getchar();}

	  sinchi=sqrt(1-coschi*coschi);
	  phi=2*M_PI*primeroot();
	  newkener=kener*coschi*coschi;
	  
	  //account for changing the reference frame, neutral velocity
	  
	  vvx=spec[i].part[j].vx*normvel/pvel;
	  vvy=spec[i].part[j].vy*normvel/pvel;
	  vvz=spec[i].part[j].vz*normvel/pvel;
	  pvel=sqrt(2*newkener*ratio*Q/mass[i]);
	  //1 - cos2, where cos=(vx,vy,vz)*(1,0,0)
	  //  minus=vvx*vvx;
	  // if(minus>1) minus=1;
	  if(vvx*vvx >= 1.0) {printf("collide 6 %E\n", vvx); vvx=0.9; getchar();}
	  sintheta=sqrt(1-vvx*vvx);
	  sss=sinchi*sin(phi)/sintheta;
	  scs=sinchi*cos(phi)/sintheta;
	  
	  spec[i].part[j].vx=(vvx*coschi+(vvy*vvy+vvz*vvz)*scs)*pvel/normvel;
	  spec[i].part[j].vy=(vvy*coschi+vvz*sss-vvx*vvy*scs)*pvel/normvel;
	  spec[i].part[j].vz=(vvz*coschi-vvy*sss-vvx*vvz*scs)*pvel/normvel;
	}
      if(colltype==1) //charge exchange
	{
	  //neutral velocity assigned afterwards
	  spec[i].part[j].vx=0;
	  spec[i].part[j].vy=0;
	  spec[i].part[j].vz=0;
	}
      
    }
 
  //printf("\nAfter collision specie %d colltype %d old ener %E vel %E new vel %E pvel, %E \n", i, colltype, kener, pvelold, pvel, normvel*sqrt(spec[i].part[j].vx*spec[i].part[j].vx+spec[i].part[j].vy*spec[i].part[j].vy+spec[i].part[j].vz*spec[i].part[j].vz));
  //if(i==1) 
  // getchar();
}
