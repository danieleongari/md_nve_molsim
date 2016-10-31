#include <stdlib.h> 
#include <stdio.h> 
#include <math.h> 
#include "system.h"

void Mdloop(void)
{ 
  int i;
  double Av[6],Tempz,dUtot,Utot0;
  VECTOR Momentum;

  FILE *FilePtr;

  FilePtr=fopen("movie.pdb","w");
 
  for(i=0;i<6;i++)
    Av[i] = 0.0;

  dUtot = 0.0;
  Utot0 = 0.0;

  // initialize radial distribution function
  SampleRDF(INITIALIZE);

  // Initialize Diffusion Coefficient
  SampleDiff(INITIALIZE);

  // loop over all cycles

  for(i=0;i<NumberOfSteps;i++)
  {
    // Calculate The Force
    Force();
 
    // integrate the equations of motion
    Integrate(i,&Momentum);
 
    if(i==1)
    {
      printf("\n");
      printf("Momentum X-Dir.        : %lf\n",Momentum.x);
      printf("Momentum Y-Dir.        : %lf\n",Momentum.y);
      printf("Momentum Z-Dir.        : %lf\n",Momentum.z);
      printf("\n");
      printf("\n");
      printf("Step: %d UTotal: %lf UKinetic: %lf UPotential: %lf Temperature: %lf Pressure: %lf\n",
        i-1,UTotal,UKinetic,UPotential,Temperature,Pressure);
      printf("\n");
    }
 
    UTotal=UKinetic+UPotential;

    // Start Modification
    // Choose a starting point to calculate the energy drift
    // Use Utot0


    // End Modification

    Tempz=2.0*UKinetic/(3.0*NumberOfParticles-3.0);

    if(i>=NumberOfInitializationSteps&&((i%50)==0))
      printf("Step: %d UTotal: %lf UKinetic: %lf UPotential: %lf Temperature: %lf Pressure: %lf\n",
      i,UTotal,UKinetic,UPotential,Tempz,Pressure);
  
    if(i%10==0) WritePdb(FilePtr);

    if(i==NumberOfSteps-1)
    {
      printf("\n");
      printf("Momentum X-Dir.        : %lf\n",Momentum.x);
      printf("Momentum Y-Dir.        : %lf\n",Momentum.y);
      printf("Momentum Z-Dir.        : %lf\n",Momentum.z);
      printf("\n");
    }
 
    // calculate averages
 
    if(i>NumberOfInitializationSteps)
    {
      Av[0]+=Tempz;
      Av[1]+=Pressure;
      Av[2]+=UKinetic;
      Av[3]+=UPotential;
      Av[4]+=UTotal;
      Av[5]+=1.0;


    //  Start Modification
    //  update the energy drift
    //  (Use dUtot.)


    //  End Modification

    }

    // sample radial distribution function

    if((i%100)==0) SampleRDF(SAMPLE);
    SampleDiff(SAMPLE);
  }
 
  // print averages to screen
 
  if(Av[5]>0.5)
  {
    Av[5]=1.0/Av[5];
    for(i=0;i<5;i++)
      Av[i]*=Av[5];
 
    printf("\n");
    printf("Average Temperature   : %lf\n",Av[0]);
    printf("Average Pressure      : %lf\n",Av[1]);
    printf("Average UKinetic      : %lf\n",Av[2]);
    printf("Average UPotential    : %lf\n",Av[3]);
    printf("Average UTotal        : %lf\n",Av[4]);
  // print the average energy drift
    printf("Average DeltaU        : %lf\n",dUtot*Av[5]);

  }

  SampleRDF(WRITE_RESULTS);
  SampleDiff(WRITE_RESULTS);

  fclose(FilePtr);
}
