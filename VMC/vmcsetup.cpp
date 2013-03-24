
#include "vmc.h"

using namespace std;
using namespace arma;

VMCSetup::VMCSetup(int nC, string atomType){
    nCycles = nC;
    Atom *a = new Atom("Be");
    atom = new Atom(atomType);
    wave = new WaveFunction(a, true, true, true, 3.96, 0.105);



}


void VMCSetup::runBruteForceSimulation(){
    double e;
    double es;
    #pragma omp parallel for
    for(int i=0;i<1;i++){
        for(int j=0;j<1;j++){
            VMCSolver *solver = new VMCSolver(wave, 1e6, true);
            solver->runMonteCarloIntegration(&e,&es);
            //energyMatrix(i,j) = e;
            //energySquaredMatrix(i,j) = es;
        }
    }
}

void VMCSetup::runSingleSimulation(){

}

void VMCSetup::runConjGradSimulation(){

}
