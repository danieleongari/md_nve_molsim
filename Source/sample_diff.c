#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"

// MAXT  = maximum timesteps for the correlation time
// MAXT0 = maximum number of time origins 
// FREQT0 = frequency with which a new time origin is selected

// time = Current time in units of Deltat
// t0time[TMAXT0] = Time of the stored time origin in units of Deltat
// t0Counter = Number of time origins stored so far
// t0index = The index of the current t0 in the array of stored time origins.
// CorrelTime = Time difference between the current time and the time origin
// SampleCounter[MAXT] = Number of samples taken at a given CorrelTime.

#define MAXT 1000
#define MAXT0 100
#define FREQT0 100

void SampleDiff(int Switch)
{
  int t0index,i,t,tmax,CorrelTime;
  double CumIntegration;
  static int time,t0time[MAXT0],t0Counter,SampleCounter[MAXT];
  static double Vacf[MAXT],R2[MAXT];
  static double Vxt0[MAXPART][MAXT0],Vyt0[MAXPART][MAXT0],Vzt0[MAXPART][MAXT0];
  static double Rx0[MAXPART][MAXT0],Ry0[MAXPART][MAXT0],Rz0[MAXPART][MAXT0];
  FILE *FilePtrMsd,*FilePtrVacf;

  switch(Switch)
  {
    // initialize everything
    case INITIALIZE:
      time=-1;
      t0Counter=-1;
      for(i=0;i<MAXT;i++) //MAXT-1 maximum time of time step (-1 because I decided to include 0)
      {
        SampleCounter[i]=0;  //
        R2[i]=0.0;           //zeros vector
        Vacf[i]=0.0;         //
      }
      break;

    case SAMPLE:
      time++;
      if((time%FREQT0)==0) //new time origin (t0) each FREQT0 time steps
      {
        // start modification: see algorithm 8 (page 91) of Frenkel/Smit
	t0Counter++;
	t0index=(t0Counter%MAXT0); //It's a way to store a maximum of MAXT0 time origins (t0), and if the index exceed MAXT0 become 1 again
	t0time[t0index]=time;
	for(i=0;i<NumberOfParticles;i++)  //t0index=(0 ... MAXT0) >>>> info about t0time, R_0[i], V_t0[i] 
	{
	  Rx0[i][t0index]=PositionsNONPDB[i].x;   //Store Positions @ t0 ---- I don't want wrapped ones: NONPDB
 	  Ry0[i][t0index]=PositionsNONPDB[i].y;	  
	  Rz0[i][t0index]=PositionsNONPDB[i].z;
	  Vxt0[i][t0index]=Velocities[i].x;       //Store Velocities @ t0
	  Vyt0[i][t0index]=Velocities[i].y;
	  Vzt0[i][t0index]=Velocities[i].z;
	}
        // end modification
      }//--------------------------end of collecting Pos,Vel @ t0

      // start modification: loop over all time origins that have been stored

      if (t0Counter<MAXT0)     //***
	tmax=t0Counter;        //***
      else                     //***
	tmax=MAXT0;            //***
 
      for(t=0;t<tmax;t++)      //*** tmax=min(t0Counter,MAXT0)
      { 
        CorrelTime=time-t0time[t]; //Time since t0
	if(CorrelTime<MAXT)
	{
	  SampleCounter[CorrelTime]++; //Counter to make the average later
	  for(i=0;i<NumberOfParticles;i++)
	  {
 	    Vacf[CorrelTime]+=( Velocities[i].x*Vxt0[i][t]+   // vector v@(t0+CorrelTime) * vector v@(t0)         
                                Velocities[i].y*Vyt0[i][t]+
                                Velocities[i].z*Vzt0[i][t]);           
	      R2[CorrelTime]+=SQR(PositionsNONPDB[i].x-Rx0[i][t])+
			      SQR(PositionsNONPDB[i].y-Ry0[i][t])+
			      SQR(PositionsNONPDB[i].z-Rz0[i][t]);
	  }
        }
      }  
   // end modification
      break;

    case WRITE_RESULTS:
      // write everything to disk
      CumIntegration=0.0;
      FilePtrMsd=fopen("msd.dat","w");
      FilePtrVacf=fopen("vacf.dat","w");
      for(i=0;i<MAXT;i++)
      {
        if(SampleCounter[i]>0)
        {
          Vacf[i]/=(double)(NumberOfParticles*SampleCounter[i]);  //To obtain an average from the big sum
          R2[i]/=(double)(NumberOfParticles*SampleCounter[i]);    //
        }
        else //if no samples of that CorrelTime t
        {
          Vacf[i]=0.0;
          R2[i]=0.0;
        }
        CumIntegration+=Deltat*Vacf[i]/3.0;
        fprintf(FilePtrVacf,"%lf %lf %lf\n",i*Deltat,Vacf[i],CumIntegration);

        if(i==0)
          fprintf(FilePtrMsd,"%lf %lf %lf\n",i*Deltat,R2[i],0.0);
        else
          fprintf(FilePtrMsd,"%lf %lf %lf\n",i*Deltat,R2[i],R2[i]/(6.0*i*Deltat));
      }
      printf("****************************************************\n");
      printf("Diffusion from velocity autocorrelation: D= %lf \n",CumIntegration);
      printf("Diffusion from mean square displacement: D= %lf \n",R2[MAXT-1]/(6.0*(MAXT-1)*Deltat));
      printf("****************************************************\n");
      fclose(FilePtrVacf);
      fclose(FilePtrMsd);
  }
}
