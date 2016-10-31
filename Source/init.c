#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "system.h"
#include "ran_uniform.h"
 
// generates initial positions/velocities

void Init(void)
{
  int i,j,k,Number,Nplace;
 
  double Uold,tmp,Place,Size,scale;
  VECTOR OldForces[MAXPART];
  VECTOR Momentum;
 
  // generate velocities from a Gaussian; set impulse to zero

  Momentum.x=0.0;
  Momentum.y=0.0;
  Momentum.z=0.0;
  Uold=0.0;
  UKinetic=0.0;
 
  for(i=0;i<NumberOfParticles;i++)
  {
    Velocities[i].x=RandomGaussianNumber();
    Velocities[i].y=RandomGaussianNumber();
    Velocities[i].z=RandomGaussianNumber();
    Momentum.x+=Velocities[i].x;
    Momentum.y+=Velocities[i].y;
    Momentum.z+=Velocities[i].z;
  }
 
  Momentum.x/=NumberOfParticles;
  Momentum.y/=NumberOfParticles;
  Momentum.z/=NumberOfParticles;

  // calculate the kinetic energy UKinetic
 
  for(i=0;i<NumberOfParticles;i++)
  {
    Velocities[i].x-=Momentum.x;
    Velocities[i].y-=Momentum.y;
    Velocities[i].z-=Momentum.z;
 
    UKinetic+=(SQR(Velocities[i].x)+SQR(Velocities[i].y)+SQR(Velocities[i].z));
  }
 
  // scale all velocities to the correct temperature
 
  scale=sqrt(Temperature*(3.0*NumberOfParticles-3.0)/UKinetic);
  for(i=0;i<NumberOfParticles;i++)
  {
    Velocities[i].x*=scale;
    Velocities[i].y*=scale;
    Velocities[i].z*=scale;
  }
 
  // calculate initial positions on a lattice
 
  Number=(int)(pow(NumberOfParticles,1.0/3.0)+1.5);
  Size=Box/(Number+2.0);
 
  Place=0.2*Size;
 
  Nplace=0;
  for(i=0;i<Number;i++)
    for(j=0;j<Number;j++)
       for(k=0;k<Number;k++)
       {
         if(Nplace<NumberOfParticles)
         {
           Positions[Nplace].x=(i+0.01*(RandomNumber()-0.5))*Size;
           Positions[Nplace].y=(j+0.01*(RandomNumber()-0.5))*Size;
           Positions[Nplace].z=(k+0.01*(RandomNumber()-0.5))*Size;
         }
         Nplace++;
       }
 
  // calculate better positions using a steepest decent algorithm
 
  for(j=0;j<50;j++)
  {
    if(j==0)
    {
      Force();
      Uold=UPotential;
      printf("Initial Energy        : %lf\n",Uold);
    }
 
    // calculate maximum downhill gradient
    tmp=0.0;
    for(i=0;i<NumberOfParticles;i++)
    {
      OldPositions[i].x=Positions[i].x;
      OldPositions[i].y=Positions[i].y;
      OldPositions[i].z=Positions[i].z;
 
      OldForces[i].x=Forces[i].x;
      OldForces[i].y=Forces[i].y;
      OldForces[i].z=Forces[i].z;
 
      Momentum.x=fabs(Forces[i].x);
      Momentum.y=fabs(Forces[i].y);
      Momentum.z=fabs(Forces[i].z);
 
      if(Momentum.x>tmp) tmp=Momentum.x;
      if(Momentum.y>tmp) tmp=Momentum.y;
      if(Momentum.z>tmp) tmp=Momentum.z;
    }
 
    tmp=Place/tmp;

    // calculate improved positions
    for(i=0;i<NumberOfParticles;i++)
    {
      Positions[i].x+=tmp*Forces[i].x;
      Positions[i].y+=tmp*Forces[i].y;
      Positions[i].z+=tmp*Forces[i].z;
 
      // place particles back in the box
      if(Positions[i].x>=Box)
        Positions[i].x-=Box;
      else if(Positions[i].x<0.0)
        Positions[i].x+=Box;
 
      if(Positions[i].y>=Box)
        Positions[i].y-=Box;
      else if(Positions[i].y<0.0)
        Positions[i].y+=Box;
 
      if(Positions[i].z>=Box)
        Positions[i].z-=Box;
      else if(Positions[i].z<0.0)
        Positions[i].z+=Box;
    }
 
    // calculate new potential energy
    Force();
 
    // check if the new positions are acceptable
    if(UPotential<Uold)
    {
      Uold=UPotential;
      Place*=1.2;

      if(Place>0.5*Box) Place=0.5*Box;
    } 
    else
    {
      for(i=0;i<NumberOfParticles;i++)
      {
        Forces[i].x=OldForces[i].x;
        Forces[i].y=OldForces[i].y;
        Forces[i].z=OldForces[i].z;
        Positions[i].x=OldPositions[i].x;
        Positions[i].y=OldPositions[i].y;
        Positions[i].z=OldPositions[i].z;
      }
      Place*=0.1;
    }
  }
 
  // calculate previous position using the generated velocity
  for(i=0;i<NumberOfParticles;i++)
  {
    OldPositions[i].x=Positions[i].x-Deltat*Velocities[i].x;
    OldPositions[i].y=Positions[i].y-Deltat*Velocities[i].y;
    OldPositions[i].z=Positions[i].z-Deltat*Velocities[i].z;
    PositionsNONPDB[i].x=Positions[i].x;
    PositionsNONPDB[i].y=Positions[i].y;
    PositionsNONPDB[i].z=Positions[i].z;
  }
  printf("Final Energy          : %lf\n",Uold); 
} 
