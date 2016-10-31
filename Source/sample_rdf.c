#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

#define Maxx 500

// samples the radial distribution function
void SampleRDF(int Ichoise)
{
  int i,j;
  double r2;
  static double Ggt,Gg[Maxx],Delta;
  VECTOR dr;
  FILE *FilePtr;

  switch(Ichoise)
  {
    case INITIALIZE:
      for(i=0;i<Maxx;i++)
        Gg[i]=0.0;

      Ggt=0.0;
      Delta=Box/(2.0*Maxx);
      break;
    case SAMPLE:
      Ggt+=1.0;
      // loop over all pairs
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

          r2=sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));

          // calculate in which bin this interaction is in
          if(r2<0.5*Box)
            Gg[(int)(r2/Delta)]+=2.0;
        }
      }
      break;
    case WRITE_RESULTS:
      // Write Results To Disk
      FilePtr=fopen("rdf.dat","w");
      for(i=0;i<Maxx-1;i++)
      {
        r2=(4.0*M_PI/3.0)*(NumberOfParticles/CUBE(Box))*CUBE(Delta)*(CUBE(i+1)-CUBE(i));
        fprintf(FilePtr,"%f %f\n",(i+0.5)*Delta,Gg[i]/(Ggt*NumberOfParticles*r2));
      }
      fclose(FilePtr);
      break;
  }
}
