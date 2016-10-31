#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..
void WritePdb(FILE *FilePtr)
{
  int i;
  static int Countmodel=0,Countatom=0;

  Countmodel++;
  fprintf(FilePtr,"%s %9d\n","MODEL",Countmodel);

  for(i=0;i<NumberOfParticles;i++)
  {
    Countatom++;
    fprintf(FilePtr,"%s%7d%s%12d%8.3lf%8.3lf%8.3lf\n",
      "ATOM",Countatom,"  H",Countatom,Positions[i].x*2.0,Positions[i].y*2.0,Positions[i].z*2.0);
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}
