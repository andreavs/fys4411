#ifndef VMCSOLVER_H
#define VMCSOLVER_H

// runs monte carlo integration with a method spesified on a wavefunction
// spesified by the user

class VMCSolver{
public:
    VMCSolver(WaveFunction *psi, int cycles, bool importance);
    double runMonteCarloIntegration(double *e, double *es);
    double sign(double x){return ((x>0) - (x<0));}


private:
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

    double stepLength;

    void equilibriateSystem();
    void findOptimalStepLength();
    void setTrialPositions();
    void findNextPos(int i, double *waveFunctionOld);


};

#endif // VMCSOLVER_H
