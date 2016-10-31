#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "system.h"

// read in system information
void Readdat(void)
{ 
  FILE *FilePtr;

  FilePtr=fopen("input","r");
  fscanf(FilePtr,"%lf %d %d %lf %lf %d\n",   
         &Box,&NumberOfParticles,&NumberOfSteps,&Temperature,
         &Deltat,&NumberOfInitializationSteps);
  fclose(FilePtr);

  if(NumberOfParticles>MAXPART)
  {
    printf("Maximum number of particles is : %d\n",MAXPART);
    exit(1);
  }
 
  // Calculate Some Parameters
  CutOff=0.49999*Box;
  Ecut=4.0*(pow(SQR(CutOff),-6.0)-pow(SQR(CutOff),-3.0));
 
  // print information to the screen
 
  printf("Molecular Dynamics Program\n");
  printf("\n");
  printf("Number of particles   : %d\n",NumberOfParticles);
  printf("Boxlength             : %f\n",Box);
  printf("Density               : %f\n",NumberOfParticles/(Box*Box*Box));
  printf("Temperature           : %f\n",Temperature);
  printf("Cut-Off radius        : %f\n",CutOff);
  printf("Cut-Off energy        : %f\n",Ecut);
  printf("Number of uteps       : %d\n",NumberOfSteps);
  printf("Number of init steps  : %d\n",NumberOfInitializationSteps);
  printf("Timestep              : %f\n",Deltat);
}
