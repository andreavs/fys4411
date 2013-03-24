#ifndef VMCSETUP_H
#define VMCSETUP_H

//sets up a variational monte carlo from config file, and runs dat shit yo.
//should be able to choose method (brute force vs gradient) here.

class VMCSetup
{
public:
    VMCSetup(int nC, std::string atomType);
    void runBruteForceSimulation();
    void runSingleSimulation();
    void runConjGradSimulation();

private:
    int nCycles;
    Atom* atom;
    WaveFunction* wave;
};

#endif // VMCSETUP_H
