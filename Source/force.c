#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "system.h"
 
// Calculate The Forces And Potential Energy
void Force(void)
{
  int i,j;
  double Ff,r2i,r6i;
  VECTOR dr;
 
  // set forces, potential energy and pressure to zero
  for(i=0;i<NumberOfParticles;i++)
  {
    Forces[i].x=0.0;
    Forces[i].y=0.0;
    Forces[i].z=0.0;
  }
 
  UPotential=0.0;
  Pressure=0.0;
 
  // loop over all particle pairs
  for(i=0;i<NumberOfParticles-1;i++)
  {
    for(j=i+1;j<NumberOfParticles;j++)
    {
      dr.x=Positions[i].x-Positions[j].x;
      dr.y=Positions[i].y-Positions[j].y;
      dr.z=Positions[i].z-Positions[j].z;

      // apply boundary conditions 
      dr.x-=Box*rint(dr.x/Box);
      dr.y-=Box*rint(dr.y/Box);
      dr.z-=Box*rint(dr.z/Box);

      r2i=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);

      // check if the distance is within the cutoff radius
      if(r2i<SQR(CutOff))
      {
        r2i=1.0/r2i;
        r6i=CUBE(r2i);
 
        UPotential+=4.0*r6i*(r6i-1.0)-Ecut;
        Ff=48.0*r2i*(r6i-0.5);
        Pressure+=Ff;
        Ff=Ff*r2i;
 
        Forces[i].x+=Ff*dr.x;
        Forces[i].y+=Ff*dr.y;
        Forces[i].z+=Ff*dr.z;
 
        Forces[j].x-=Ff*dr.x;
        Forces[j].y-=Ff*dr.y;
        Forces[j].z-=Ff*dr.z;
      }
    }
  }
 
  // scale the pressure
  Pressure/=3.0*CUBE(Box);
} 
