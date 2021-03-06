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
    vec3 diff = zeros(3);
    double a;
    vec3 r12 = zeros(3);
    double r12n;


    for(int j=0; j<nParticles;j++){
        if(j != row){
            if((j+row)%2 == 1){a = 0.5;}
            else{a = 0.25;}
            //r12 = r.row(row) - r.row(j);
            r12(0) = r.at(row,0) - r.at(j,0); r12.at(1) = r.at(row,1) - r.at(j,1); r12.at(2) = r.at(row,2) - r.at(j,2);
            //r12n = sqrt((r(row,0) - r(j,0))*(r(row,0) - r(j,0)) + (r(row,1) - r(j,1))*(r(row,1) - r(j,1)) + (r(row,2) - r(j,2))*(r(row,2) - r(j,2)));
            r12n = sqrt(r12.at(0)*r12.at(0) + r12.at(1)*r12.at(1) + r12.at(2)*r12.at(2));
            diff = diff + a*r12/(r12n*(beta*r12n + 1)*(beta*r12n + 1));
        }
    }

    // numerical (old) implementation:
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
    vec3 rik = zeros(3);
    double rikn;
    double rijn;
    double partsum = 0;
    double sum = 0;


    // numerical (old) implementation:
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
                rik.at(0) = r.at(i,0) - r.at(k,0); rik.at(1) = r.at(i,1) - r.at(k,1); rik.at(2) = r.at(i,2) - r.at(k,2);
                rikn = sqrt(rik.at(0)*rik.at(0) + rik.at(1)*rik.at(1) + rik.at(2)*rik.at(2));
                sum += (nDimensions-1)/rikn*aik/((1+beta*rikn)*(1+beta*rikn));
                sum -= 2*aik*beta/(pow(1+beta*rikn,3));
            }
        }
    }



    return sum;
}

double Jastrow::dPsidBeta(const mat &r){
    double aij;
    double rij;
    double ret = 0;
    for(int i=0; i<nParticles; i++){
        for(int j=i+1; j<nParticles; j++){
            if((j+i)%2 == 1){aij = 0.5;}
            else{aij = 0.25;}
            rij = (r(i,0)-r(j,0))*(r(i,0)-r(j,0)) + (r(i,1)-r(j,1))*(r(i,1)-r(j,1)) + (r(i,2)-r(j,2))*(r(i,2)-r(j,2));
            ret -= aij*rij*rij/((1+beta*rij)*(1+beta*rij));
        }
    }
    return ret;
}

double Jastrow::d2PsidBeta2(const mat &r){
    double aij;
    double rij;
    double ret = 0;
    for(int i=0; i<nParticles; i++){
        for(int j=i+1; j<nParticles; j++){
            if((j+i)%2 == 1){aij = 0.5;}
            else{aij = 0.25;}
            rij = (r(i,0)-r(j,0))*(r(i,0)-r(j,0)) + (r(i,1)-r(j,1))*(r(i,1)-r(j,1)) + (r(i,2)-r(j,2))*(r(i,2)-r(j,2));
            ret -= (aij*rij*rij/((1+beta*rij)*(1+beta*rij)))*(aij*rij*rij/((1+beta*rij)*(1+beta*rij))) - 3*rij*(aij*rij*rij/((1+beta*rij)*(1+beta*rij)*(1+beta*rij))) ;
        }
    }
    return ret;
}
