#include <stdlib.h>
#include <stdio.h>

#define MAXPART 200

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define MIN(x,y) ((x)<(y)?(x):(y))

enum{INITIALIZE,SAMPLE,WRITE_RESULTS};

typedef struct
{
  double x;
  double y;
  double z;
} VECTOR;

extern VECTOR Forces[MAXPART];           // forces
extern VECTOR Positions[MAXPART];        // positions
extern VECTOR OldPositions[MAXPART];    // old positions
extern VECTOR Velocities[MAXPART];        // velocities
extern VECTOR PositionsNONPDB[MAXPART]; // positions that are not put back in the box

extern double Deltat;         // Timestep

extern double Box;           // Boxlengths
extern double CutOff;        // Cut-Off Radius
extern double Ecut;          // Cut-Off Energy

extern double UKinetic;      // Kinetic Energy
extern double UPotential;    // Potential Energy
extern double UTotal;        // Total Energy

extern double Temperature;   // Temperature
extern double Pressure;      // Pressure

extern int NumberOfSteps;
extern int NumberOfInitializationSteps;
extern int NumberOfParticles;

void WritePdb(FILE *FilePtr);
void SampleRDF(int Ichoise);
void SampleDiff(int Switch);
void Force(void);
void Integrate(double Step,VECTOR *Momentum);
void Readdat(void);
void Init(void);
void Mdloop(void);
