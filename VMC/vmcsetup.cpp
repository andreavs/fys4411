
#include "vmc.h"



using namespace std;
using namespace arma;

VMCSetup::VMCSetup(int nC, string atomType){
    nCycles = nC;
    //Atom *a = new Atom("Ne");
    atom = new Atom(atomType);
    molecule = atom->isMolecule;
    //wave = new WaveFunction(a, false, false, true, 10, 0.105);



}



vec VMCSetup::runBruteForceSimulation(double firstAlpha, double lastAlpha, int nAlphas, double firstBeta, double lastBeta, int nBetas, bool incJas, bool incSelf, bool incPreComp){
    double e;
    double es;
    mat energyMatrix = zeros<mat>(nAlphas, nBetas);
    mat energySquaredMatrix = zeros<mat>(nAlphas, nBetas);
    vec alpha = linspace(firstAlpha, lastAlpha, nAlphas);
    vec beta = linspace(firstBeta, lastBeta, nBetas);
    alpha.print("alpha");
    beta.print("beta");
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;
    WaveFunction *wave;
    VMCSolver *solver;
#pragma omp parallel for num_threads(3)
    for(int i=0;i<nAlphas;i++){
        for(int j=0;j<nBetas;j++){
            wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha[i], beta[j], molecule);
            solver = new VMCSolver(wave, nCycles, true, molecule);
            solver->runMonteCarloIntegration(&e,&es);
            energyMatrix(i,j) = e;
            energySquaredMatrix(i,j) = es;
        }
    }
    energyMatrix.print("Energy:");
    energySquaredMatrix.print("Energy squared: ");


    uword  row;
    uword col;
    double min_val = energyMatrix.min(row, col);
    cout << "the minimal energy is " << min_val << " with indices " << row << " " << col << endl;
    int i = int(row);
    int j = int(col);
    double bestAlpha = alpha(i);
    double bestBeta = beta(j);
    cout << "alpha: " << bestAlpha << " beta: " << bestBeta << endl;
    vec ret;
    ret << bestAlpha << bestBeta;
    return ret;
}

void VMCSetup::runSingleSimulation(double alpha, double beta, bool incJas, bool incSelf, bool incPreComp, bool storeResults, bool storePositions, bool saveR12){
    int numprocs, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand(my_rank);
    double e;
    double es;
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;
    WaveFunction *wave;

    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta, molecule);
    VMCSolver *solver = new VMCSolver(wave, nCycles, true, molecule, saveR12);
    if(storeResults)solver->runMonteCarloIntegration(&e,&es,storeResults, storePositions);
    else solver->runMonteCarloIntegration(&e,&es);


    cout << "Thread nr. " << my_rank << " has Energy: " << e << " Energy Squared; " << es << endl;
}

void VMCSetup::runConjGradSimulation(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp){
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;

    double alpha = firstAlpha;
    double beta = firstBeta;
    int counter = 0;
    vec r = zeros(2);
    r(0) = alpha; r(1) = beta;

    mat hessian = zeros(2,2);
    colvec gradient = zeros(2);
    WaveFunction *wave;
    VMCSolver *solver;
    double dpda;
    double dpdb;
    double error = tolerance + 1;
    //VMCSolver *solver = new VMCSolver(wave, nCycles, true);
    cout << counter << " " << maxIterations << endl;

    double c;
    while(counter < maxIterations && error > tolerance){
        wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta, molecule);
        solver = new VMCSolver(wave, nTestCycles, true, molecule);
        solver->runDerivativeCalculation(&dpda, &dpdb);
        gradient << dpda << dpdb;

        counter++;
    }
}

vec VMCSetup::runBroydensMethod(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp){
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;

    double alpha = firstAlpha;
    double beta = firstBeta;
    int counter = 0;
    colvec2 gradient = zeros(2);
    colvec2 aplusgradient = zeros(2);
    colvec2 aminusgradient = zeros(2);
    colvec2 bplusgradient = zeros(2);
    colvec2 bminusgradient = zeros(2);
    //colvec2 oldGradient = zeros(2);
    WaveFunction *wave;
    VMCSolver *solver;
    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta, molecule);
    double dpda;
    double dpdb;
    double error = tolerance + 1;

    solver = new VMCSolver(wave, 10*nTestCycles, true, molecule);
    solver->runDerivativeCalculation(&dpda, &dpdb);
    gradient << dpda << dpdb;
    vec gradientOld = gradient;
    double ah = 0.01;
    double bh = 0.01;

    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha+ah, beta);
    solver = new VMCSolver(wave, 10*nTestCycles, true);
    solver->runDerivativeCalculation(&dpda, &dpdb);
    aplusgradient << dpda << dpdb;

    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha-ah, beta);
    solver = new VMCSolver(wave, 10*nTestCycles, true);
    solver->runDerivativeCalculation(&dpda, &dpdb);
    aminusgradient << dpda << dpdb;

    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta+bh);
    solver = new VMCSolver(wave, 10*nTestCycles, true);
    solver->runDerivativeCalculation(&dpda, &dpdb);
    bplusgradient << dpda << dpdb;

    wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta-bh);
    solver = new VMCSolver(wave, 10*nTestCycles, true);
    solver->runDerivativeCalculation(&dpda, &dpdb);
    bminusgradient << dpda << dpdb;

    //VMCSolver *solver = new VMCSolver(wave, nCycles, true);
    cout << counter << " " << maxIterations << endl;

    mat B = zeros(2,2);
    B(0,0) = (aplusgradient(0) - aminusgradient(0))/(2*ah);
    B(1,1) = (bplusgradient(1) - bminusgradient(1))/(2*bh);
    B(1,0) = (aplusgradient(1) - aminusgradient(1))/(2*bh);
    B(0,1) = B(1,0); //(bplusgradient(0) - bminusgradient(0))/(2*ah);
    B.print("B");
    colvec2 p;
    colvec2 s;
    colvec2 y;
    //colvec2 gradientOld = ones(2);
    double step = 0.5;

    while(counter < maxIterations && error > tolerance){

        gradient.print("g");
        //c = min(2.0,1./((counter+1)*norm(gradient,2)));///(counter+1);

        p = -inv(B)*gradient;
        s = p*step;

        alpha = fabs(alpha + step*p(0));
        beta =fabs(beta + step*p(1));
        wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta);
        solver = new VMCSolver(wave, nTestCycles, true);
        solver->runDerivativeCalculation(&dpda, &dpdb);
        gradient << dpda << dpdb;
        y = gradient - gradientOld;
        y.print("y");
         ((1./dot(s,s))*(y - B*s)*s.t()).print("yoyo");
        B = B + (1./dot(s,s))*(y - B*s)*s.t();

        B.print("B");





        cout << alpha << " " << beta << " " << endl;


        if(dot(gradient,gradientOld) < 0){
            step /= 2.;
        }
        else{
            step *= 1.1;
        }
        gradientOld = gradient;
    }
    gradient << alpha << beta;
    return gradient;
}



vec VMCSetup::runSteepestDescent(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp){
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;

    double alpha = firstAlpha;
    double beta = firstBeta;
    int counter = 0;
    vec alphavals = zeros(maxIterations);
    vec betavec = zeros(maxIterations);
    vec energies = zeros(maxIterations);
    colvec gradient = zeros(2);
    colvec oldGradient = zeros(2);
    WaveFunction *wave;
    VMCSolver *solver;
    double dpda;
    double dpdb;
    double e;
    double error = tolerance + 1;

    int s1;
    int s2;
    //VMCSolver *solver = new VMCSolver(wave, nCycles, true);
    cout << counter << " " << maxIterations << endl;

    double alphaStep = 0.5;
    double betaStep = 0.1;

    while(counter < maxIterations && error > tolerance){
        alphavals(counter) = alpha;
        betavec(counter) = beta;
        wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta, molecule);
        solver = new VMCSolver(wave, nTestCycles, true, molecule);
        solver->runDerivativeCalculation(&dpda, &dpdb, &e);
        gradient << dpda << dpdb;
        gradient.print("g");
        energies(counter) = e;
        //c = min(2.0,1./((counter+1)*norm(gradient,2)));///(counter+1);
        alpha = 10.0;//max(alpha - alphaStep*sign(gradient(0)),0.01);
        beta = beta - betaStep*sign(gradient(1));
        if(beta<0){
            beta = 0.01;
            betaStep = 0.005;
        }

        s1 = sign(gradient(0));
        s2 = sign(oldGradient(0));
        if(s1 != s2){
            alphaStep /= 2.;
            cout << "alphastep: " << alphaStep << endl;
        }
        else{
            alphaStep *= 1.1;
            cout << "alphastep: " << alphaStep << endl;
        }
        cout << "balle" << sign(gradient(1)) << " " << sign(oldGradient(1)) << endl;
        s1 = sign(gradient(1));
        s2 = sign(oldGradient(1));
        if(s1 != s2){
           betaStep /= 2.;
            cout << "betastep: " << betaStep << endl;
        }
        else{
            betaStep *= 1.1;
            cout << "betastep: " << betaStep << endl;
        }
        cout << alpha << " " << beta << " " << endl;


        counter++;
        oldGradient = gradient;
    }
    gradient << alpha << beta;
    alphavals.print("alpha:");
    betavec.print("beta");
    energies.print("energies");
    return gradient;
}


vec VMCSetup::runNewtonsMethod(double firstAlpha, double firstBeta, int maxIterations, double tolerance, int nTestCycles, bool incJas, bool incSelf, bool incPreComp){
    bool includeJastrow = incJas;
    bool includeSelf = incSelf;
    bool includePreComputed = incPreComp;

    double alpha = firstAlpha;
    double beta = firstBeta;
    int counter = 0;
    colvec r = zeros(2);
    r(0) = alpha; r(1) = beta;

    mat hessian = zeros(2,2);
    colvec gradient = zeros(2);
    WaveFunction *wave;
    VMCSolver *solver;
    double dpda;
    double dpdb;
    double dpdaa;
    double dpdbb;
    double dpdab;
    double error = tolerance + 1;
    //VMCSolver *solver = new VMCSolver(wave, nCycles, true);
    cout << counter << " " << maxIterations << endl;

    double c;
    while(counter < maxIterations && error > tolerance){
        wave = new WaveFunction(atom, includeJastrow, includeSelf, includePreComputed, alpha, beta);
        solver = new VMCSolver(wave, nTestCycles, true);
        solver->run2DerivativeCalculation(&dpda, &dpdb, &dpdaa, &dpdbb, &dpdab);
        gradient << dpda << dpdb;
        hessian(0,0) = dpdaa; hessian(1,1) = dpdbb;
        hessian(1,0) = dpdab; hessian(0,1) = dpdab;
        r = r - inv(hessian)*gradient;
        hessian.print("hess");
        r.print("ab");
        counter++;
    }
}
