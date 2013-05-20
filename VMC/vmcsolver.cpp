#include "vmc.h"


using namespace std;
using namespace arma;
#include "zignor.h"

VMCSolver::VMCSolver(WaveFunction *psi, int cycles, bool importance, bool incMolecule, bool saveR12)
{
    this->saveR12 = saveR12;
    useImportance = importance;
    nCycles = cycles;
    waveFunction = psi;
    nParticles = waveFunction->getNParticles();
    nDimensions = waveFunction->getNDimensions();

    stepLength = 1.0;
    idum = -1;
    rnd = new Random(idum);
    int seed = 10;
    RanNormalSetSeedZigVec(&seed, 100);
    molecule = incMolecule;

    R << 1.4 << 0 << 0;

}

double VMCSolver::runMonteCarloIntegration(double *e, double *es, bool saveResults, bool savePositions){
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
    waveFunction->slater->allNewPosUpdate(rNew);
    // loop over Monte Carlo cycles
    vec energies = zeros(nCycles);

    cube positions;
    cube xyzHistogram;
    mat xyHistogram;
    vec xHistogram;
    mat radii;
    double width;
    double r12 = 0;
    double binWidth;
    int nPoints;
    int numprocs, my_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    if(my_rank != 0){savePositions = false;}
    if(savePositions){
        nPoints = 200;
        width = 5;
        binWidth = (2.0*width)/nPoints;
        xyHistogram = zeros(nPoints,nPoints);
        xHistogram = zeros(nPoints);
        xyzHistogram = zeros(nPoints, nPoints, nPoints);
        positions = zeros(nParticles, nDimensions, nCycles);
        radii = zeros(nParticles,nCycles);
    }

    double acceptedCounter = 0;
    for(int cycle = 0; cycle < nCycles; cycle++) {

        if(cycle % 100000 == 0 && my_rank ==0){cout << "Now at cycle " << cycle << " out of " << nCycles << endl;}

        // Store the current value of the wave function
        waveFunctionOld = waveFunction->waveFunction(rOld);
        //


        // New position to test
        for(int i = 0; i < nParticles; i++) {

            findNextPos(i, &waveFunctionOld);
            if(stepAccepted) acceptedCounter++;
            // update energies
            deltaE = waveFunction->localEnergy(rNew);
            energies(cycle) = deltaE;
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;

        }
        if(savePositions){
            for(int i=0; i<nParticles; i++){
                radii(i,cycle) = sqrt(rNew(i,0)*rNew(i,0) + rNew(i,1)*rNew(i,1) + rNew(i,2)*rNew(i,2));
            }
            int xindex;
            int yindex;
            int zindex;
            for(int i=0; i<nParticles; i++){
                for(int j=0; j<nDimensions;j++){
                    positions(i,j,cycle) = rNew(i,j);
                }
                xindex = (rNew(i,0)+width)/binWidth;
                yindex = (rNew(i,1)+width)/binWidth;
                zindex = (rNew(i,2)+width)/binWidth;
                if (xindex < nPoints && xindex >= 0){
                    xHistogram(xindex) += 1;
                    if (yindex < nPoints && yindex >= 0){
                        xyHistogram(xindex, yindex) += 1;
                        if (zindex < nPoints && zindex >= 0){
                            xyzHistogram(xindex, yindex, zindex) += 1;
                        }
                    }
                }
            }
        }
        if(saveR12){
            for(int i=0; i<nParticles; i++){
                for(int j=i+1; j<nParticles; j++){
                    double rOne = 0;
                    for(int k=0; k<nDimensions; k++){
                        rOne += (rNew(i,k) - rNew(j,k))*(rNew(i,k) - rNew(j,k));

                    }
                    rOne = sqrt(rOne);
                    r12 += rOne;
                }
            }
        }
    }

    if(saveR12){
        r12 = r12/(nCycles*nParticles*(nParticles-1)/2.);
        if(my_rank == 0){
            cout << "avg distance between particles: " << r12 << endl;
        }
    }
    if(savePositions){
        string path1 = "results/" + waveFunction->myAtom->getAtomName() + "positions.mat";
        string path2 = "results/" + waveFunction->myAtom->getAtomName() + "radii.mat";
        string path3 = "results/" + waveFunction->myAtom->getAtomName() + "x.mat";
        string path4 = "results/" + waveFunction->myAtom->getAtomName() + "xy.mat";
        string path5 = "results/" + waveFunction->myAtom->getAtomName() + "xyz.mat";
        positions.save(path1);
        radii.save(path2);
        xHistogram.save(path3);
        xyHistogram.save(path4);
        xyzHistogram.save(path5);

    }
    double energy = energySum/(nCycles * nParticles);
    double energySquared = energySquaredSum/(nCycles * nParticles);
    //cout << "Energy: " << energy << " Energy (squared sum): " << energySquared << endl;
    *e = energy;
    *es = energySquared;
    if(saveResults){
        int numprocs, my_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        string dataPath = "results";
        energies.save(dataPath + "/data" + to_string(my_rank) + ".mat");
        //cout << "lol" << endl;
    }
    double tot_e;
    double tot_es;
    MPI_Allreduce(&energy, &tot_e, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&energySquared, &tot_es, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(my_rank ==0){
        tot_e = tot_e/numprocs;
        tot_es = tot_es/numprocs;
        cout << "Energy averaged over all nodes: " << tot_e << " Energy squared: " << tot_es << endl;
    }
    acceptedCounter = ((double) acceptedCounter)/ (nCycles*nParticles);
    cout<< "acceptance ratio: " << acceptedCounter << endl;
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

double VMCSolver::runDerivativeCalculation(double *dpda, double *dpdb, double *E){
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    setTrialPositions();

    double waveFunctionOld = waveFunction->waveFunction(rOld);
//    double waveFunctionNew = 0;

    double alphaSum = 0;
    double betaSum = 0;
    double ddbetaSum = 0;
    double alphaEnergySum = 0;
    double betaEnergySum = 0;
    double energySum = 0;

    double dAlpha = 0;
    double dBeta = 0;
    double dEnergy = 0;
    double ddAlpha = 0;
    double ddBeta = 0;




    // find step length and equilibriate

    equilibriateSystem();
    rNew = rOld;
    waveFunction->slater->allNewPosUpdate(rNew);
    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction->waveFunction(rOld);
        //


        // New position to test
        for(int i = 0; i < nParticles; i++) {
            findNextPos(i, &waveFunctionOld);
            // update energies
            dAlpha = waveFunction->dPhidAlpha(rNew);
            dBeta = waveFunction->dPhidBeta(rNew);
            dEnergy = waveFunction->localEnergy(rNew);

            alphaSum += dAlpha;
            betaSum += dBeta;
            energySum += dEnergy;
            alphaEnergySum += dAlpha*dEnergy;
            betaEnergySum += dBeta*dEnergy;
            //energySquaredSum += deltaE*deltaE;
        }
    }
    double alpha = alphaSum/(nCycles * nParticles);
    double beta = betaSum/(nCycles * nParticles);
    double energy = energySum/(nCycles * nParticles);
    double ae = alphaEnergySum/(nCycles * nParticles);
    double be = betaEnergySum/(nCycles * nParticles);
    cout << "Energy: " << energy << endl;
    *dpda = 2*(ae - energy*alpha);
    *dpdb = 2*(be - energy*beta);
    *E = energy;
}

double VMCSolver::run2DerivativeCalculation(double *dpda, double *dpdb, double *dpdaa, double *dpdbb, double *dpdab){
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    setTrialPositions();

    double waveFunctionOld = waveFunction->waveFunction(rOld);
//    double waveFunctionNew = 0;

    double alphaSum = 0;
    double betaSum = 0;
    double alphaEnergySum = 0;
    double betaEnergySum = 0;
    double aaSum = 0;
    double bbSum = 0;
    double abSum = 0;
    double abEnergySum = 0;
    double aaEnergySum = 0;
    double bbEnergySum = 0;
    double energySum = 0;
    double a2EnergySum = 0;
    double a2Sum = 0;
    double b2EnergySum = 0;
    double b2Sum = 0;

    double dAlpha = 0;
    double dBeta = 0;
    double ddAlpha = 0;
    double ddBeta = 0;
    double dAlphadBeta = 0;
    double dEnergy = 0;





    // find step length and equilibriate

    equilibriateSystem();
    rNew = rOld;
    waveFunction->slater->allNewPosUpdate(rNew);
    // loop over Monte Carlo cycles
    for(int cycle = 0; cycle < nCycles; cycle++) {
        // Store the current value of the wave function
        waveFunctionOld = waveFunction->waveFunction(rOld);
        //


        // New position to test
        for(int i = 0; i < nParticles; i++) {
            findNextPos(i, &waveFunctionOld);
            // update energies
            dAlpha = waveFunction->dPhidAlpha(rNew);
            dBeta = waveFunction->dPhidBeta(rNew);
            ddAlpha = waveFunction->d2PhidAlpha2(rNew);
            ddBeta = waveFunction->d2PhidBeta2(rNew);
            dEnergy = waveFunction->localEnergy(rNew);

            alphaSum += dAlpha;
            betaSum += dBeta;
            energySum += dEnergy;
            aaSum += ddAlpha;
            bbSum += ddBeta;
            abSum += dAlpha*dBeta;
            aaEnergySum += ddAlpha*dEnergy;
            bbEnergySum += ddBeta*dEnergy;
            abEnergySum += dAlpha*dBeta*dEnergy;
            alphaEnergySum += dAlpha*dEnergy;
            betaEnergySum += dBeta*dEnergy;
            a2EnergySum += dAlpha*dAlpha*dEnergy;
            a2Sum += dAlpha*dAlpha;
            b2EnergySum += dBeta*dBeta*dEnergy;
            b2Sum += dBeta*dBeta;
            //energySquaredSum += deltaE*deltaE;
        }
    }
    double alpha = alphaSum/(nCycles * nParticles);
    double beta = betaSum/(nCycles * nParticles);
    double energy = energySum/(nCycles * nParticles);
    double ae = alphaEnergySum/(nCycles * nParticles);
    double be = betaEnergySum/(nCycles * nParticles);
    double bbe = bbEnergySum/(nCycles * nParticles);
    double aae = aaEnergySum/(nCycles * nParticles);
    double abe = abEnergySum/(nCycles * nParticles);
    double ab = abSum/(nCycles * nParticles);
    double aa = aaSum/(nCycles * nParticles);
    double bb = bbSum/(nCycles * nParticles);
    double a2e = a2EnergySum/(nCycles * nParticles);
    double a2 = a2Sum/(nCycles * nParticles);
    double b2 = b2Sum/(nCycles * nParticles);
    double b2e = b2EnergySum/(nCycles * nParticles);
    cout << "Energy: " << energy << endl;
    *dpda = 2*(ae - energy*alpha);
    *dpdb = 2*(be - energy*beta);
    *dpdaa = 2*(aae + 3*a2e - aa*energy-3*a2*energy - 2*alpha*(*dpda));
    *dpdbb = 2*(bbe + 3*b2e - bb*energy-3*b2*energy - 2*beta*(*dpdb));
    *dpdab = 2*(4*abe - 4*ab*energy - alpha*(*dpdb) - beta*(*dpda));
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
    waveFunction->slater->allNewPosUpdate(rNew);
}

void VMCSolver::findNextPos(int i, double *waveFunctionOld){
    if(useImportance){

        stepAccepted = false;

        //importance step
        double dt = 0.02;
        FOld = waveFunction->drift(rOld,i);
        double D = 0.5;
        double stddev = sqrt(2*D*dt);
        int rand;
        for(int j=0; j<nDimensions;j++){
            rNew(i,j) = rOld(i,j) + stddev*randn() + D*FOld(j)*dt;
        }
        waveFunction->slater->updateAll(rNew, i);
        double waveFunctionNew = waveFunction->waveFunction(rNew);

        FNew = waveFunction->drift(rNew, i);
        // Check for step acceptance (if yes, update position, if no, reset position)
        double greensfunc = 0;

        for(int j=0;j<nDimensions;j++){
            greensfunc += 0.5*(FOld(j) + FNew(j))*(D*dt*0.5*(FOld(j) - FNew(j)) - rNew(i,j) + rOld(i,j));
        }

        greensfunc = exp(greensfunc);
        double testvar = greensfunc*(waveFunctionNew*waveFunctionNew) / ((*waveFunctionOld)*(*waveFunctionOld));

        if(randu() <= testvar) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = rNew(i,j);
                *waveFunctionOld = waveFunctionNew;
            }
            stepAccepted = true;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j);

            }
            waveFunction->slater->updateAll(rNew,i);
        }

    }
    // brute force algo:
    else{
        stepAccepted = false;

        //random uniform step
        for(int j=0; j<nDimensions;j++){
            rNew(i,j) = rOld(i,j) + stepLength*(randu()-0.5);
        }
        waveFunction->slater->updateAll(rNew,i);


        // Recalculate the value of the wave function
        double waveFunctionNew = waveFunction->waveFunction(rNew);
        // Check for step acceptance (if yes, update position, if no, reset position)
        if(randu() <= (waveFunctionNew*waveFunctionNew) / ((*waveFunctionOld)*(*waveFunctionOld))) {
            for(int j = 0; j < nDimensions; j++) {
                rOld(i,j) = rNew(i,j);
                *waveFunctionOld = waveFunctionNew;
            }
            stepAccepted = true;
        } else {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j);
            }
            waveFunction->slater->updateAll(rNew,i);
        }
    }
}
