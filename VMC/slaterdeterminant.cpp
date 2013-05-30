
#include "vmc.h"
using namespace std;
using namespace arma;

SlaterDeterminant::SlaterDeterminant(Atom *atom, int nP, int nD, double a, const mat &r, bool preComp, bool molecule)
{
    nParticles = nP;
    nDimensions = nD;
    slater = zeros(nParticles, nDimensions);
    alpha = a;
    includePrecomputed = preComp;
    //cout << "lol: " << includePrecomputed << endl;


    slater = zeros(nParticles, nParticles);
    slaterUp = zeros(nParticles/2, nParticles/2);
    slaterDown = zeros(nParticles/2, nParticles/2);
    slaterUpInverse = zeros(nParticles/2, nParticles/2);
    slaterDownInverse = zeros(nParticles/2, nParticles/2);

    oldSlaterUp = zeros(nParticles/2, nParticles/2);
    oldSlaterDown = zeros(nParticles/2, nParticles/2);


    dPhidAlphaUp = zeros(nParticles/2, nParticles/2);
    dPhidAlphaDown = zeros(nParticles/2, nParticles/2);

    d2PhidAlpha2Up = zeros(nParticles/2, nParticles/2);
    d2PhidAlpha2Down = zeros(nParticles/2, nParticles/2);

    // makin waves yo. yall aint never gun hold me down no mo
    if(!molecule){
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
    }
    else{
        cout << "lolololo" << endl;
        orbitals.push_back(new MoleculeOrbital(1,0,0,nD,nP,a, atom,0));
        orbitals.push_back(new MoleculeOrbital(1,0,0,nD,nP,a, atom,0));
        orbitals.push_back(new MoleculeOrbital(1,0,0,nD,nP,a, atom,1));
        orbitals.push_back(new MoleculeOrbital(1,0,0,nD,nP,a, atom,1));
        orbitals.push_back(new MoleculeOrbital(2,0,0,nD,nP,a, atom,1));
        orbitals.push_back(new MoleculeOrbital(2,0,0,nD,nP,a, atom,1));
        orbitals.push_back(new MoleculeOrbital(2,0,0,nD,nP,a, atom,0));
        orbitals.push_back(new MoleculeOrbital(2,0,0,nD,nP,a, atom,0));
        orbitals.push_back(new MoleculeOrbital(2,1,-1,nD,nP,a, atom,1));
        orbitals.push_back(new MoleculeOrbital(2,1,-1,nD,nP,a, atom,1));
    }


    updateMatrix(r);
    oldSlaterUp = slaterUp;
    oldSlaterDown = slaterDown;
    ratioUp = 0;
    ratioDown = 0;
    //slaterUp.print("up");
    //slaterDown.print("down");
    slaterUpInverse = inv(slaterUp);
    slaterDownInverse = inv(slaterDown);
    slaterUpDeterminant = det(slaterUp);
    slaterDownDeterminant = det(slaterDown);
    updateLaplace(r);
    //updateAll(r);

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

void SlaterDeterminant::allNewPosUpdate(const mat &r){
    updateMatrix(r);
    slaterUpInverse = inv(slaterUp);
    slaterDownInverse = inv(slaterDown);
    slaterUpDeterminant = det(slaterUp);
    slaterDownDeterminant = det(slaterDown);
    updateLaplace(r);

}

void SlaterDeterminant::updateAll(const mat &r, int i){
    updateMatrix(r);
    //slaterUp.print("swags");
    //cout << "yolo" << endl;
    updateDeterminant(r,i);
    updateInverse(r,i);
    updateLaplace(r);
}

void SlaterDeterminant::updateMatrix(const mat &r){
    oldSlaterDown = slaterDown;
    oldSlaterUp = slaterUp;
    for(int i=0; i<nParticles/2; i++){
        for(int j=0; j<nParticles/2; j++){
            slaterUp(i,j) = orbitals[2*i]->waveFunction(r.row(j));
            //r.row(j).print("swag");
            slaterDown(i,j) = orbitals[2*i+1]->waveFunction(r.row(nParticles/2 + j));
        }
    }

}

void SlaterDeterminant::updateInverse(const mat &, int i){
    mat newSlaterUpInverse = zeros<mat>(nParticles/2,nParticles/2);
    mat newSlaterDownInverse = zeros<mat>(nParticles/2,nParticles/2);
    double sum;
    //cout << ratioUp << endl;
    if(i<nParticles/2){
        for(int j=0; j<nParticles/2; j++){
            sum = 0;
            if(j != i){
                for(int l=0; l<nParticles/2; l++){
                    sum += slaterUp.at(l,i)*slaterUpInverse.at(j,l);
                }

                for(int k=0; k<nParticles/2; k++){
                    newSlaterUpInverse.at(j,k) = slaterUpInverse.at(j,k) - slaterUpInverse.at(i,k)*sum/ratioUp;
                }
            }
            else{
                for(int l=0; l<nParticles/2; l++){
                    sum += oldSlaterUp.at(l,i)*slaterUpInverse.at(j,l);
                }

                for(int k=0; k<nParticles/2; k++){
                    newSlaterUpInverse.at(j,k) = slaterUpInverse.at(i,k)*sum/ratioUp;
                }
            }
        }
        slaterUpInverse = newSlaterUpInverse;
    }
    else{
        i = i-nParticles/2;
        for(int j=0; j<nParticles/2; j++){
            sum = 0;
            if(j != i){
                for(int l=0; l<nParticles/2; l++){
                    sum += slaterDown.at(l,i)*slaterDownInverse.at(j,l);
                }

                for(int k=0; k<nParticles/2; k++){
                    newSlaterDownInverse.at(j,k) = slaterDownInverse.at(j,k) - slaterDownInverse.at(i,k)*sum/ratioDown;
                }
            }
            else{
                for(int l=0; l<nParticles/2; l++){
                    sum += oldSlaterDown.at(l,i)*slaterDownInverse.at(j,l);
                }

                for(int k=0; k<nParticles/2; k++){
                    newSlaterDownInverse.at(j,k) = slaterDownInverse.at(i,k)*sum/ratioDown;
                }
            }
        }
        slaterDownInverse = newSlaterDownInverse;
    }

    //cout << "inv" << endl;
    //slaterUp.print();
    //slaterUpInverse = inv(slaterUp);
    //cout << "inv" << endl;
    //slaterDownInverse = inv(slaterDown);
}

void SlaterDeterminant::updateDeterminant(const mat &r, int i){
    ratioUp = 0;
    ratioDown = 0;
    if(i<nParticles/2){
        for(int j=0; j<nParticles/2; j++){
            ratioUp += slaterUp.at(j,i)*slaterUpInverse.at(i,j);
        }
        slaterUpDeterminant = ratioUp*slaterUpDeterminant;
    }
    else{
        i = i - nParticles/2;
        for(int j=0; j<nParticles/2; j++){
            ratioDown += slaterDown.at(j,i)*slaterDownInverse.at(i,j);
        }
        slaterDownDeterminant = ratioDown*slaterDownDeterminant;
    }

//    slaterUpDeterminant = det(slaterUp);
//    slaterDownDeterminant = det(slaterDown);
}

void SlaterDeterminant::updateLaplace(const mat &r){
    //updateInverse(r);
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
    //rupdateMatrix(r);
    //updateInverse(r);
    //updateDeterminant(r);
    double ret;
    if(includePrecomputed){
        ret =  slaterUpDeterminant*slaterDownDeterminant;
    }
    else{
        ret = getDeterminant(r);
        //cout << "lol" << endl;

    }
    return ret;
}

vec SlaterDeterminant::gradient(const mat &r, int i){
    //returns grad(D)/D
    //updateAll(r,i);
    slaterDownGradient = zeros(nDimensions);
    slaterUpGradient = zeros(nDimensions);

    if(includePrecomputed){
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
    }
    else{
        double h = 1e-6;
        vec diff = slaterDownGradient;
        mat rPlus = zeros<mat>(nParticles, nDimensions);
        mat rMinus = zeros<mat>(nParticles, nDimensions);
        rPlus = rMinus = r;
        double waveFunctionPlus;
        double waveFunctionMinus;
        for(int j=0; j<nDimensions;j++){
            rPlus(i,j) += h;
            rMinus(i,j) -= h;
            waveFunctionPlus = getDeterminant(rPlus);
            waveFunctionMinus = getDeterminant(rMinus);
            diff(j) = waveFunctionPlus - waveFunctionMinus;
            rPlus(i,j) -= h;
            rMinus(i,j) += h;
        }
        diff = diff/(2*h*getDeterminant(r));
        slaterUpGradient = diff;
    }

    return slaterDownGradient + slaterUpGradient;

}

double SlaterDeterminant::laplace(const mat &r){
    //updateAll(r);
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


double SlaterDeterminant::dPsidAlpha(const mat &r){
    for(int i=0; i<nParticles/2;i++){
        for(int j=0; j<nParticles/2; j++){
            dPhidAlphaUp(i,j) = orbitals[2*i]->dPhidAlpha(r.row(j));
            //r.row(j).print("swag");
            dPhidAlphaDown(i,j) = orbitals[2*i+1]->dPhidAlpha(r.row(nParticles/2 + j));
        }
    }

    double retUp = 0;
    double retDown = 0;

    for(int i=0; i<nParticles/2;i++){
        for(int j=0; j<nParticles/2; j++){
            retUp += slaterUpInverse(j,i)*dPhidAlphaUp(i,j);
            retDown += slaterDownInverse(j,i)*dPhidAlphaDown(i,j);
        }
    }
    double ret = retUp + retDown;
    return ret;
}

double SlaterDeterminant::d2PsidAlpha2(const mat &r){
    for(int i=0; i<nParticles/2;i++){
        for(int j=0; j<nParticles/2; j++){
            d2PhidAlpha2Up(i,j) = orbitals[2*i]->d2PhidAlpha2(r.row(j));
            //r.row(j).print("swag");
            d2PhidAlpha2Down(i,j) = orbitals[2*i+1]->d2PhidAlpha2(r.row(nParticles/2 + j));
        }
    }

    double retUp = 0;
    double retDown = 0;

    for(int i=0; i<nParticles/2;i++){
        for(int j=0; j<nParticles/2; j++){
            retUp += slaterUpInverse(j,i)*d2PhidAlpha2Up(i,j);
            retDown += slaterDownInverse(j,i)*d2PhidAlpha2Down(i,j);
        }
    }
    double ret = retUp + retDown;
    return ret;
}
