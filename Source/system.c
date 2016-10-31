#include "system.h"

VECTOR Positions[MAXPART];        // positions
VECTOR OldPositions[MAXPART];    // old positions
VECTOR Velocities[MAXPART];        // velocities
VECTOR Forces[MAXPART];           // forces
VECTOR PositionsNONPDB[MAXPART]; // positions that are not put back in the box

double Deltat;         // Timestep

double Box;           // Boxlengths
double CutOff;        // Cut-Off Radius
double Ecut;          // Cut-Off Energy

double UKinetic;      // Kinetic Energy
double UPotential;    // Potential Energy
double UTotal;        // Total Energy

double Temperature;   // Temperature
double Pressure;      // Pressure

int NumberOfSteps;
int NumberOfInitializationSteps;              
int NumberOfParticles;  

