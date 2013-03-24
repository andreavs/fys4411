
#include "vmc.h"
using namespace std;
using namespace arma;

SlaterDeterminant::SlaterDeterminant(int nP, int nD, double a, const mat &r)
{
    nParticles = nP;
    nDimensions = nD;
    slater = zeros(nParticles, nDimensions);
    alpha = a;


    slater = zeros(nParticles, nParticles);
    slaterUp = zeros(nParticles/2, nParticles/2);
    slaterDown = zeros(nParticles/2, nParticles/2);
    slaterUpInverse = zeros(nParticles/2, nParticles/2);
    slaterDownInverse = zeros(nParticles/2, nParticles/2);


    // makin waves yo. yall aint never gun hold me down no mo
    orbitals.push_back(new Orbital(1,0,0,nD,nP,a));
    orbitals.push_back(new Orbital(1,0,0,nD,nP,a));
    orbitals.push_back(new Orbital(2,0,0,nD,nP,a));
    orbitals.push_back(new Orbital(2,0,0,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,-1,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,-1,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,0,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,0,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,1,nD,nP,a));
    orbitals.push_back(new Orbital(2,1,1,nD,nP,a));
    updateAll(r);

}

double SlaterDeterminant::getDeterminant(const mat &r){

    for(int i=0; i<nParticles/2; i++){
        for(int j=0; j<nParticles/2; j++){
            slaterUp(i,j) = orbitals[2*i]->waveFunction(r.row(j));
            slaterDown(i,j) = orbitals[2*i+1]->waveFunction(r.row(nParticles/2 + j));
        }
    }

    double ret = det(slaterUp)*det(slaterDown);
    return ret;
}

void SlaterDeterminant::updateAll(const mat &r){
    updateMatrix(r);
    updateInverse(r);
    updateDeterminant(r);
    updateLaplace(r);
}

void SlaterDeterminant::updateMatrix(const mat &r){
    for(int i=0; i<nParticles/2; i++){
        for(int j=0; j<nParticles/2; j++){
            slaterUp(i,j) = orbitals[2*i]->waveFunction(r.row(j));
            slaterDown(i,j) = orbitals[2*i+1]->waveFunction(r.row(nParticles/2 + j));
        }
    }
}

void SlaterDeterminant::updateInverse(const mat &){
    slaterUpInverse = inv(slaterUp);
    slaterDownInverse = inv(slaterDown);
}

void SlaterDeterminant::updateDeterminant(const mat &r){
    slaterUpDeterminant = det(slaterUp);
    slaterDownDeterminant = det(slaterDown);
}

void SlaterDeterminant::updateLaplace(const mat &r){
    updateInverse(r);
    slaterDownLaplace = 0;
    slaterUpLaplace = 0;


    for(int i=0; i<nParticles/2; i++){
        for(int j=0; j<nParticles/2;j++){
            slaterUpLaplace += orbitals[2*j]->laplace(r.row(i))*slaterUpInverse(i,j);
            slaterDownLaplace += orbitals[2*j+1]->laplace(r.row(nParticles/2 + i))*slaterDownInverse(i,j);

        }
    }

//    mat rPlus = zeros<mat>(nParticles, nDimensions);
//    mat rMinus = zeros<mat>(nParticles, nDimensions);

//    rPlus = rMinus = r;

//    double waveFunctionMinus = 0;
//    double waveFunctionPlus = 0;

//    double waveFunctionCurrent = waveFunction(r);
//    // Kinetic energy
//    double h = 1e-5;
//    double h2 = 1/(h*h);
//    double kineticEnergy = 0;
//    for(int i = 0; i < nParticles; i++) {
//        for(int j = 0; j < nDimensions; j++) {
//            rPlus(i,j) += h;
//            rMinus(i,j) -= h;
//            waveFunctionMinus = waveFunction(rMinus);
//            waveFunctionPlus = waveFunction(rPlus);
//            kineticEnergy -= (waveFunctionMinus + waveFunctionPlus - 2 * waveFunctionCurrent);
//            rPlus(i,j) = r(i,j);
//            rMinus(i,j) = r(i,j);
//        }
//    }
//    kineticEnergy =  -h2 * kineticEnergy / waveFunctionCurrent;
//    slaterDownLaplace = kineticEnergy;


}


double SlaterDeterminant::waveFunction(const mat &r){
    updateMatrix(r);
    updateInverse(r);
    updateDeterminant(r);
    return slaterUpDeterminant*slaterDownDeterminant;
}

vec SlaterDeterminant::gradient(const mat &r, int i){
    //returns grad(D)/D
    updateAll(r);
    slaterDownGradient = zeros(nDimensions);
    slaterUpGradient = zeros(nDimensions);

//    double h = 1e-6;
//    vec diff = slaterDownGradient;
//    mat rPlus = zeros<mat>(nParticles, nDimensions);
//    mat rMinus = zeros<mat>(nParticles, nDimensions);
//    rPlus = rMinus = r;
//    double waveFunctionPlus;
//    double waveFunctionMinus;
//    for(int j=0; j<nDimensions;j++){
//        rPlus(i,j) += h;
//        rMinus(i,j) -= h;
//        waveFunctionPlus = waveFunction(rPlus);
//        waveFunctionMinus = waveFunction(rMinus);
//        diff(j) = waveFunctionPlus - waveFunctionMinus;
//        rPlus(i,j) -= h;
//        rMinus(i,j) += h;
//    }
//    diff = diff/(2*h*waveFunction(r));
//    slaterUpGradient = diff;



    int dummyi;
    for(int j=0; j<nParticles/2;j++){
        if(i<nParticles/2){
            slaterUpGradient += orbitals[2*j]->gradient(r.row(i))*slaterUpInverse(i,j);
        }
        else{
            dummyi = i-nParticles/2;
            slaterDownGradient += orbitals[2*j+1]->gradient(r.row(i))*slaterDownInverse(dummyi,j);
        }
    }




    return slaterDownGradient + slaterUpGradient;

}

double SlaterDeterminant::laplace(const mat &r){
    updateAll(r);
    double ret = 0;
//    for(int i=0; i<nParticles;i++){
//        ret += dot(gradient(slaterUp, i), gradient(slaterDown,i));
//    }
//    ret = 2*ret;
//    r.print("r");
//    slaterUp.print("up");
//    slaterUpInverse.print("upinv");
//    slaterDownInverse.print("downinv");
//    slaterDown.print("down");

//    cout << slaterUpLaplace << " " << slaterDownLaplace << " " << slaterUpDeterminant << endl;
    ret = slaterUpLaplace + slaterDownLaplace;
    return ret;
}



