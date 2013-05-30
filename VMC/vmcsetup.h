#ifndef VMCSETUP_H
#define VMCSETUP_H

//sets up a variational monte carlo from config file, and runs dat shit yo.
//should be able to choose method (brute force vs gradient) here.

#include<armadillo>

class VMCSetup
{
public:
    VMCSetup(int nC, std::string atomType);
    vec runBruteForceSimulation(double firstAlpha, double lastAlpha, int nAlphas, double firstBeta, double lastBeta, int nBetas, bool incJas, bool incSelf, bool incPreComp);
    double runSingleSimulation(double alpha, double beta, bool incJas, bool incSelf, bool incPreComp, bool storeResults, bool storePositions, bool saveR12=false);
    void runConjGradSimulation(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp);
    vec runSteepestDescent(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp);
    int sign(double x){return ((x>0) - (x<0));}
    vec runBroydensMethod(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp);
    vec runNewtonsMethod(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp);
    int nCycles;
    Atom* atom;
    WaveFunction* wave;
    bool molecule;
    int nprocs;
    int my_rank;
};

#endif // VMCSETUP_H
