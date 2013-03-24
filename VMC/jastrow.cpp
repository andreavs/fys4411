#include "vmc.h"

Jastrow::Jastrow(int nP, int nD, double b)
{
    nParticles = nP;
    nDimensions = nD;
    beta = b;
}


double Jastrow::waveFunction(const mat &r){
    double sum = 0;
    double a;
    double r12;
    double r12n;
    for(int i=0; i<nParticles;i++){
        for(int j=i+1; j<nParticles;j++){
            if((j+i)%2 == 1){a = 0.5;}
            else{a = 0.25;}
            //r12n = norm(r.row(i) - r.row(j), 2);
            r12 = sqrt((r(i,0) - r(j,0))*(r(i,0) - r(j,0)) + (r(i,1) - r(j,1))*(r(i,1) - r(j,1)) + (r(i,2) - r(j,2))*(r(i,2) - r(j,2)));
            sum += a*r12/(1+beta*r12);
        }
    }

    return exp(sum);

}

vec Jastrow::gradient(const mat &r, int row){
    vec diff = zeros(3);
    double a;
    vec r12;
    double r12n;


    for(int j=0; j<nParticles;j++){
        if(j != row){
            if((j+row)%2 == 1){a = 0.5;}
            else{a = 0.25;}
            r12 = r.row(row) - r.row(j);
            r12n = sqrt((r(row,0) - r(j,0))*(r(row,0) - r(j,0)) + (r(row,1) - r(j,1))*(r(row,1) - r(j,1)) + (r(row,2) - r(j,2))*(r(row,2) - r(j,2)));
            diff = diff + a*r12/(r12n*(beta*r12n + 1)*(beta*r12n + 1));
        }
    }

//    double h = 1e-8;
//    mat rPlus = zeros<mat>(nParticles, nDimensions);
//    mat rMinus = zeros<mat>(nParticles, nDimensions);
//    rPlus = rMinus = r;
//    double waveFunctionPlus;
//    double waveFunctionMinus;
//    for(int i=0; i<nDimensions;i++){
//        rPlus(row,i) += h;
//        rMinus(row,i) -= h;
//        waveFunctionPlus = waveFunction(rPlus);
//        waveFunctionMinus = waveFunction(rMinus);
//        diff(i) = waveFunctionPlus - waveFunctionMinus;
//        rPlus(row,i) -= h;
//        rMinus(row,i) += h;
//    }
//    diff = diff/(2*h*waveFunction(r));

    return diff;

}

double Jastrow::laplace(const mat &r){
    double aik;
    vec rik;
    double rikn;
    double rijn;
    double partsum = 0;
    double sum = 0;


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
//    sum = kineticEnergy;

    for(int k=0; k<nParticles;k++){
        sum += dot(gradient(r,k),gradient(r,k));
        for(int i=0;i<nParticles;i++){
            if(i != k){
                if((i+k)%2 == 1){aik = 0.5;}
                else{aik = 0.25;}
                rik = r.row(i) - r.row(k);
                rikn = sqrt(rik(0)*rik(0) + rik(1)*rik(1) + rik(2)*rik(2));
                sum += (nDimensions-1)/rikn*aik/((1+beta*rikn)*(1+beta*rikn));
                sum -= 2*aik*beta/(pow(1+beta*rikn,3));
            }
        }
    }

//    double aij;
//    vec rij;
//    for(int i=0; i<nParticles;i++){
//        for(int j=0; j<nParticles;j++){
//            if(j != i){
//                if((j+i)%2 == 1){aij = 0.5;}
//                else{aij = 0.25;}
//                rij = r.row(i) - r.row(j);
//                rijn = sqrt((r(i,0) - r(j,0))*(r(i,0) - r(j,0)) + (r(i,1) - r(j,1))*(r(i,1) - r(j,1)) + (r(i,2) - r(j,2))*(r(i,2) - r(j,2)));

//                for(int k = 0; k<nParticles;k++){
//                    if(k != j && k != i){
//                        if((j+k)%2 == 1){aik = 0.5;}
//                        else{aik = 0.25;}
//                        rik = r.row(i) - r.row(k);
//                        rikn = sqrt(rik(0)*rik(0) + rik(1)*rik(1) + rik(2)*rik(2));
//                        partsum += aik*dot(rik, rij)/(rikn*rijn*(beta*rikn +1)*(beta*rikn +1));
//                    }
//                }
//                partsum = 2*partsum;

//                sum += aij/((beta*rijn + 1)*(beta*rijn + 1))*(partsum + aij/(((beta*rijn + 1)*(beta*rijn + 1))) + 2/(rijn*(beta*rijn + 1)));
//                partsum = 0.0;
//            }
//        }
//    }

    return sum;
}
