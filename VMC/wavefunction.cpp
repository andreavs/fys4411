#include "vmc.h"
#include "random.h"

using namespace std;
using namespace arma;

WaveFunction::WaveFunction(Atom* atom, bool incJas, bool incSelf, bool incPreComp, double a, double b)
{
    includeJastrow = incJas;
    includeSelfColoumb = incSelf;
    includePreComputedLocalEnergy = incPreComp;
    myAtom = atom;
    nParticles = myAtom->getNParticles();
    nDimensions = myAtom->getNDimensions();

    h = 1e-5;
    h2 = 1/(h*h);
    charge = nParticles;

    idum = -1;
    rnd = new Random(idum);


    alpha = a;
    beta = b;

    mat dummyr = randn<mat>(nParticles, nDimensions);
    slater = new SlaterDeterminant(nParticles, nDimensions, alpha, dummyr);
    jastrow = new Jastrow(nParticles, nDimensions, b);
}


double WaveFunction::waveFunction(const mat &r){
    double arg = slater->waveFunction(r);
    double jastrowArgument = 1;
    double a;
    if(includeJastrow){
        jastrowArgument = jastrow->waveFunction(r);
    }
    return arg*jastrowArgument;
}

double WaveFunction::localEnergy(const mat &r){
    if(includePreComputedLocalEnergy){
        double kineticTerms = 0;
        kineticTerms += -0.5*slater->laplace(r);
        double crossSum = 0;
        if(includeJastrow){
            kineticTerms -= 0.5*jastrow->laplace(r);
            for(int i=0; i<nParticles;i++){
                crossSum += dot(slater->gradient(r,i),jastrow->gradient(r,i));
                //cout << crossSum << endl;
            }
            kineticTerms -= crossSum;
        }

        double r1;
        double r12;
        double energy = kineticTerms;
        for(int i=0; i<nParticles; i++){
            r1 = sqrt(r(i,0)*r(i,0) + r(i,1)*r(i,1) + r(i,2)*r(i,2));
            energy -= (charge)*(1/r1);
            if(includeSelfColoumb){
                for(int j=i+1; j<nParticles; j++){
                    r12 = sqrt((r(i,0) - r(j,0))*(r(i,0) - r(j,0)) + (r(i,1) - r(j,1))*(r(i,1) - r(j,1)) + (r(i,2) - r(j,2))*(r(i,2) - r(j,2)));
                    energy += 1/r12;
                }
            }
        }
        //cout << energy << endl;
        return energy;


//    //precomputed local energy for helium atom, to speed up the code
//    double r1, r2, r12, r1dotr2;
//    double El1 =  - alpha*alpha;
//    double El2 = 0;
//    double jastrow = 0;
//    double self = 0;
//    double a;
//    for(int i=0; i<nParticles; i++){
//        r1 = sqrt(r(i,0)*r(i,0) + r(i,1)*r(i,1) + r(i,2)*r(i,2));
//        El1 += (alpha - 2)*(1/r1);
//        for(int j=i+1; j<nParticles; j++){
//            r2 = sqrt(r(j,0)*r(j,0) + r(j,1)*r(j,1) + r(j,2)*r(j,2));
//            r12 = sqrt((r(i,0) - r(j,0))*(r(i,0) - r(j,0)) + (r(i,1) - r(j,1))*(r(i,1) - r(j,1)) + (r(i,2) - r(j,2))*(r(i,2) - r(j,2)));
//            r1dotr2 = r(i,0)*r(j,0) + r(i,1)*r(j,1) + r(i,2)*r(j,2);
//            if(includeSelfColoumb){
//                self += 1/r12;
//            }
//            if(includeJastrow){
//                if((j+i)%2 == 1){a = 0.5;}
//                else{a = 0.25;}
//                jastrow += a/((beta*r12 + 1)*(beta*r12 + 1))*(alpha*(r1+r2)*(1 -
//                    r1dotr2/(r1*r2))/r12 - a/((beta*r12 + 1)*(beta*r12 + 1)) - 2/(r12*(beta*r12 + 1)));
//            }
//        }
//    }
//    El2 += El1 + self + jastrow;
//    return El2;
    }
    else{
        mat rPlus = zeros<mat>(nParticles, nDimensions);
        mat rMinus = zeros<mat>(nParticles, nDimensions);

        rPlus = rMinus = r;

        double waveFunctionMinus = 0;
        double waveFunctionPlus = 0;

        double waveFunctionCurrent = waveFunction(r);

        // Kinetic energy

        double kineticEnergy = 0;
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rPlus(i,j) += h;
                rMinus(i,j) -= h;
                waveFunctionMinus = waveFunction(rMinus);
                waveFunctionPlus = waveFunction(rPlus);
                kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
                rPlus(i,j) = r(i,j);
                rMinus(i,j) = r(i,j);
            }
        }
        kineticEnergy = 0.5 * h2 * kineticEnergy / waveFunctionCurrent;
//        kineticEnergy += 0.5*jastrow->laplace(r);
//        double crossSum = 0;
//        for(int i=0; i<nParticles;i++){
//            crossSum += dot(slater->gradient(r,i),jastrow->gradient(r,i));
//        }
//        kineticEnergy-=crossSum;
        // Potential energy
        double potentialEnergy = 0;
        double rSingleParticle = 0;
        for(int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r(i,j)*r(i,j);
            }
            potentialEnergy -= charge / sqrt(rSingleParticle);
        }


        // Contribution from electron-electron potential
        if(includeSelfColoumb){
        double r12 = 0;
            for(int i = 0; i < nParticles; i++) {
                for(int j = i + 1; j < nParticles; j++) {
                    r12 = 0;
                    for(int k = 0; k < nDimensions; k++) {
                        r12 += (r(i,k) - r(j,k)) * (r(i,k) - r(j,k));
                    }
                    potentialEnergy += 1 / sqrt(r12);
                }
            }
        }
        //cout << kineticEnergy << " " << potentialEnergy << endl;
        //cout << potentialEnergy << endl;
        return potentialEnergy + kineticEnergy;
    }
}

vec WaveFunction::drift(const mat &r, int i){
    double ri = sqrt(r(i,0)*r(i,0) + r(i,1)*r(i,1) + r(i,2)*r(i,2));
    double factor1 = -alpha/ri;
    double factor2 = 0;
    double a;
    if(includeJastrow){
        double r12n;
        for(int i=0; i< nParticles; i++){
            for(int j=i+1;j<nParticles;j++){
                if((i+j) % 2 == 1){a = 0.5;}
                else{a = 0.25;}
                r12n = norm(r.row(i) - r.row(j), 2);
                factor2 += a/(r12n*(beta*r12n+1)*(beta*r12n+1));
            }
        }
    }
    vec F = zeros(nDimensions);
    for(int j=0; j<nDimensions; j++){
        F(j) = 2*(factor1*r(i,j) + factor2*(r(i,j) + r((i+1)%2,j)));
    }
    return F;
}

void WaveFunction::oneBodyDensity(){
    //TODO LOLLOLOLOLOLOLOLOLOLOLOLOL
    //TODO sjekke jastroer
    int N = 101;
    int nCycles = 1e7;
    vec r1 = linspace(0,3,N);
    vec oneBody = zeros(N);
    double r12, exp1, exp2;
    double theta2 = 0;
    double r2 = 0;
    for(int j=0;j<N;j++){
        for(int i=0; i<nCycles; i++){
            r2 = 5*rnd->nextDouble();
            theta2 = 3.14159*rnd->nextDouble();
            r12 = sqrt(r1(j)*r1(j) -2*r1(j)*r2*cos(theta2) + r2*r2);
            exp1 = -2*alpha*(r1(j) + r2);
            exp2 = r12/(1+beta*r12);
            oneBody(j) += r1(j)*r1(j)*sin(theta2)*r2*r2*exp(exp1 + exp2);
        }
    }
    oneBody  = oneBody/sum(oneBody);

    for(int i=0;i<N;i++){
        cout << oneBody(i) << endl;
    }
}
