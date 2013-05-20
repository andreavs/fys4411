#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include "zignor.h"
#include "zigrandom.h"

// runs monte carlo integration with a method spesified on a wavefunction
// spesified by the user

class VMCSolver{
public:
    VMCSolver(WaveFunction *psi, int cycles, bool importance, bool incMolecule=false, bool saveR12=false);
    double runMonteCarloIntegration(double *e, double *es, bool saveResults=false, bool savePositions=false);
    double runDerivativeCalculation(double *dpda, double *dpdb, double *E=0);
    double run2DerivativeCalculation(double *dpda, double *dpdb, double *dpdaa, double *dpdbb, double *dpdab);

    double sign(double x){return ((x>0) - (x<0));}



private:
    bool saveR12;

    WaveFunction *waveFunction;

    bool useImportance;
    int nCycles;
    int nParticles;
    int nDimensions;

    bool stepAccepted;

    vec3 FOld;
    vec3 FNew;

    mat rOld;
    mat rNew;

    long idum;
    Random *rnd;

    bool molecule;

    double stepLength;

    void equilibriateSystem();
    void findOptimalStepLength();
    void setTrialPositions();
    void findNextPos(int i, double *waveFunctionOld);
    vec3 R;

};

#endif // VMCSOLVER_H
