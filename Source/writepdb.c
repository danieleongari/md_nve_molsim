#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "system.h"

// make a movie file of the simulation box
// use vmd to view it..

//Daniele: modified 11/10/16 because the printed PDB it was in a wrong format (added Box, adjusted xyz columns)

void WritePdb(FILE *FilePtr)
{
  int i;
  static int Countmodel=0,Countatom=0;

  Countmodel++;
  Countatom=0;

  fprintf(FilePtr,"%s%9.3lf%9.3lf%9.3lf%s\n",
	  "CRYST1",Box,Box,Box,"  90.00  90.00  90.00 P 1           1");

  for(i=0;i<NumberOfParticles;i++)
  {
    Countatom++;
    fprintf(FilePtr,"%s%7d%s%12d%s%8.3lf%8.3lf%8.3lf\n",
      "ATOM",Countatom,"  H",Countatom," xx ",Positions[i].x*1.0,Positions[i].y*1.0,Positions[i].z*1.0);
    Countatom=0;
  }
  fprintf(FilePtr,"%s\n","ENDMDL");
}
