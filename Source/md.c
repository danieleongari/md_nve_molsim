#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "ran_uniform.h"
#include "system.h"


int main(void)
{ 
  // Molecular Dynamics Program Of Argon,     C
  // Written By Thijs Vlugt On 6-10-1998      C

  printf("Starting\n");
 
  // initialize the random number generator with the system time
  InitializeRandomNumberGenerator(time(0l));

  // Read Data From Disk
  Readdat();
 
  // Generate Initial Coordinates/Velocities
  Init();
 
  // Finally, Perform An Md Simulation 
  Mdloop();

  return 0; 
}
