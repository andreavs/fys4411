#include "vmc.h"


using namespace std;
using namespace arma;

VMCSolver::VMCSolver(WaveFunction *psi, int cycles, bool importance)
{
    useImportance = importance;
    nCycles = cycles;
    waveFunction = psi;
    nParticles = waveFunction->getNParticles();
    nDimensions = waveFunction->getNDimensions();

    stepLength = 1.0;
    idum = -1;
    rnd = new Random(idum);

}

double VMCSolver::runMonteCarloIntegration(double *e, double *es){
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    setTrialPositions();

    double waveFunctionOld = waveFunction->waveFunction(rOld);
//    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // find step length and equilibriate

    equilibriateSystem();
    rNew = rOld;

    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {

        if(cycle % 100000 == 0){cout << cycle << endl;}

        // Store the current value of the wave function
        waveFunctionOld = waveFunction->waveFunction(rOld);
        waveFunction->slater->updateAll(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            findNextPos(i, &waveFunctionOld);
            // update energies
            deltaE = waveFunction->localEnergy(rNew);

            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    *e = energy;
    *es = energySquared;
}

void VMCSolver::equilibriateSystem(){
    int equiCycles;
    setTrialPositions();

    double waveFunctionOld = waveFunction->waveFunction(rOld);
    if(useImportance){
        equiCycles = 1e5;
    }
    else{
        findOptimalStepLength();
        equiCycles = 1e6;
    }
    for(int i=0; i<equiCycles; i++){
        for(int j=0; j<nParticles;j++){
            findNextPos(j, &waveFunctionOld);
        }
    }
}

void VMCSolver::findOptimalStepLength(){
    // Uses the secant method with 10 steps and the 1.0 max step rule to find an step
    // giving ~0.5 acceptance ratio.
    // First trial step lengths are set to 1.0 and 1.4, so that's cool and stuffs.
    int trialIterations = 10;
    double acceptRatio1 = 1.0;
    double acceptRatio2 = 1.4;
    int testCycles = 10000;
    double trialStepLength1 = 1.0;
    double trialStepLength2 = 1.1;
    double idealRatio = 0.5;
    double waveFunctionOld;
    double stepLengthDiff;
    double acceptRateDiff;
    double stepRule = 1.0;
    double changeInStepLength;

    // try first step length
    stepLength = trialStepLength1;
    // initial trial positions
    setTrialPositions();
    waveFunctionOld = waveFunction->waveFunction(rOld);
    for(int i=0;i<testCycles;i++){
        for(int j=0;j<nParticles;j++){
            findNextPos(j, &waveFunctionOld);
            if(stepAccepted){
                acceptRatio1++;
            }
        }
    }
    acceptRatio1 /= nParticles*testCycles;

    for(int iterator=0; iterator<trialIterations;iterator++){
        // try other step length
        stepLength = trialStepLength2;
        // initial trial positions
        setTrialPositions();
        waveFunctionOld = waveFunction->waveFunction(rOld);
        for(int i=0;i<testCycles;i++){
            for(int j=0;j<nParticles;j++){
                findNextPos(j, &waveFunctionOld);
                if(stepAccepted){
                    acceptRatio2++;
                }
            }
        }
        acceptRatio2 /= nParticles*testCycles;
        stepLengthDiff = trialStepLength2 - trialStepLength1;
        acceptRateDiff = acceptRatio2 - acceptRatio1;
        changeInStepLength = (idealRatio - acceptRatio2)*(stepLengthDiff)/(acceptRateDiff);
        changeInStepLength = min(abs(changeInStepLength),stepRule)*sign(changeInStepLength);
        trialStepLength1 = trialStepLength2;
        acceptRatio1 = acceptRatio2;
        trialStepLength2 = trialStepLength2 + changeInStepLength;
        acceptRatio2 = 0;
    }
}

void VMCSolver::setTrialPositions(){
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = stepLength * (rnd->nextDouble() - 0.5);
        }
    }
    rNew = rOld;
    waveFunction->slater->updateAll(rNew);
}

void VMCSolver::findNextPos(int i, double *waveFunctionOld){
    if(useImportance){
        stepAccepted = false;

        //importance step
        double dt = 0.02;
        FOld = waveFunction->drift(rOld,i);
        double D = 0.5;
        double stddev = sqrt(2*D*dt);
        for(int j=0; j<nDimensions;j++){
            rNew(i,j) = rOld(i,j) + stddev*randn() + D*FOld(j)*dt;
        }
        waveFunction->slater->updateAll(rNew);
        double waveFunctionNew = waveFunction->waveFunction(rNew);

        FNew = waveFunction->drift(rNew, i);
        // Check for step acceptance (if yes, update position, if no, reset position)
        double greensfunc = 0;

        for(int j=0;j<nDimensions;j++){
            greensfunc += 0.5*(FOld(j) + FNew(j))*(D*dt*0.5*(FOld(j) - FNew(j)) - rNew(i,j) + rOld(i,j));
        }

        greensfunc = exp(greensfunc);
        double testvar = greensfunc*(waveFunctionNew*waveFunctionNew) / ((*waveFunctionOld)*(*waveFunctionOld));

        if(rnd->nextDouble() <= testvar) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = rNew(i,j);
                *waveFunctionOld = waveFunctionNew;
            }
            stepAccepted = true;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j);

            }
            waveFunction->slater->updateAll(rNew);
        }

    }
    // brute force algo:
    else{
        stepAccepted = false;

        //random uniform step
        for(int j=0; j<nDimensions;j++){
            rNew(i,j) = rOld(i,j) + stepLength*(rnd->nextDouble()-0.5);
        }
        waveFunction->slater->updateAll(rNew);


        // Recalculate the value of the wave function
        double waveFunctionNew = waveFunction->waveFunction(rNew);
        // Check for step acceptance (if yes, update position, if no, reset position)
        if(rnd->nextDouble() <= (waveFunctionNew*waveFunctionNew) / ((*waveFunctionOld)*(*waveFunctionOld))) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = rNew(i,j);
                *waveFunctionOld = waveFunctionNew;
            }
            stepAccepted = true;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j);
            }
            waveFunction->slater->updateAll(rNew);
        }
    }
}
