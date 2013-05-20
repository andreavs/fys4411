#include "vmc.h"
#include "random.h"

using namespace std;
using namespace arma;

WaveFunction::WaveFunction(Atom* atom, bool incJas, bool incSelf, bool incPreComp, double a, double b, bool molecule)
{
    includeJastrow = incJas;
    includeSelfColoumb = incSelf;
    includePreComputedLocalEnergy = incPreComp;
    myAtom = atom;
    nParticles = myAtom->getNParticles();
    nDimensions = myAtom->getNDimensions();

    h = 1e-6;
    h2 = 1/(h*h);
    charge = nParticles;

    idum = -1;
    rnd = new Random(idum);


    alpha = a;
    beta = b;

    mat dummyr = randn<mat>(nParticles, nDimensions);
    //cout << includePreComputedLocalEnergy << endl;
    slater = new SlaterDeterminant(atom, nParticles, nDimensions, alpha, dummyr, includePreComputedLocalEnergy, molecule);
    jastrow = new Jastrow(nParticles, nDimensions, b);
    this->molecule = molecule;
    R = atom->R;
}


double WaveFunction::waveFunction(const mat &r){
    double arg = slater->waveFunction(r);
    double jastrowArgument = 1;
    if(includeJastrow){
        jastrowArgument = jastrow->waveFunction(r);
    }
    return arg*jastrowArgument;
}

double WaveFunction::localEnergy(const mat &r){
    double kineticEnergy = 0;
    double potentialEnergy = 0;
    double energy = 0;
    if(includePreComputedLocalEnergy){

        kineticEnergy += -0.5*slater->laplace(r);
        double crossSum = 0;
        if(includeJastrow){
            kineticEnergy += -0.5*jastrow->laplace(r);
            for(int i=0; i<nParticles;i++){
                crossSum += dot(slater->gradient(r,i),jastrow->gradient(r,i));
                //cout << crossSum << endl;
            }
            kineticEnergy -= crossSum;
        }
    } else {
        mat rPlus = zeros<mat>(nParticles, nDimensions);
        mat rMinus = zeros<mat>(nParticles, nDimensions);

        rPlus = rMinus = r;

        double waveFunctionMinus = 0;
        double waveFunctionPlus = 0;

        double waveFunctionCurrent = waveFunction(r);

        // Kinetic energy

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
    }







//        double energy = kineticTerms;

//        for(int i=0; i<nParticles; i++){
//            r1 = sqrt(r(i,0)*r(i,0) + r(i,1)*r(i,1) + r(i,2)*r(i,2));
//            energy -= (charge)*(1/r1);
//            if(includeSelfColoumb){
//                for(int j=i+1; j<nParticles; j++){
//                    r12 = sqrt((r(i,0) - r(j,0))*(r(i,0) - r(j,0)) + (r(i,1) - r(j,1))*(r(i,1) - r(j,1)) + (r(i,2) - r(j,2))*(r(i,2) - r(j,2)));
//                    energy += 1/r12;
//                }
//            }
//        }
//        //cout << energy << endl;
//        return energy;


    if(!molecule){
        double rSingleParticle = 0;
        for(int i = 0; i < nParticles; i++) {
            rSingleParticle = 0;
            for(int j = 0; j < nDimensions; j++) {
                rSingleParticle += r(i,j)*r(i,j);
            }
            potentialEnergy -= charge / sqrt(rSingleParticle);
        }
    } else{
        double ch = myAtom->charge;
        vec3 r1p1; vec r1p2; //double r1p1n = 0; double r1p2n = 0;
        for(int i=0; i<nParticles; i++){
            r1p1 = r.row(i) + R/2;
            r1p2 = r.row(i) - R/2;
            double r1p1n = 0; double r1p2n = 0;
            for(int j=0; j<nDimensions; j++){
                r1p1n += r1p1(j)*r1p1(j);
                r1p2n += r1p2(j)*r1p2(j);
            }
            r1p1n = sqrt(r1p1n);
            r1p2n = sqrt(r1p2n);
            potentialEnergy -= ch/r1p1n;
            potentialEnergy -= ch/r1p2n;
        }
        double Rn = 0;
        for(int i=0; i<nDimensions; i++){
            Rn += R(i)*R(i);
        }
        Rn = sqrt(Rn);
        potentialEnergy += ch*ch/Rn;
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

vec WaveFunction::drift(const mat &r, int i){
//    double ri = sqrt(r(i,0)*r(i,0) + r(i,1)*r(i,1) + r(i,2)*r(i,2));
//    double factor1 = -alpha/ri;
//    double factor2 = 0;
//    double a;
//    if(includeJastrow){
//        double r12n;
//        for(int i=0; i< nParticles; i++){
//            for(int j=i+1;j<nParticles;j++){
//                if((i+j) % 2 == 1){a = 0.5;}
//                else{a = 0.25;}
//                r12n = norm(r.row(i) - r.row(j), 2);
//                factor2 += a/(r12n*(beta*r12n+1)*(beta*r12n+1));
//            }
//        }
//    }
    vec F = zeros(nDimensions);
//    for(int j=0; j<nDimensions; j++){
//        F(j) = 2*(factor1*r(i,j) + factor2*(r(i,j) + r((i+1)%2,j)));
//    }
    F += slater->gradient(r,i);
    if(includeJastrow){
        F += jastrow->gradient(r,i);
    }
    return F;
}

double WaveFunction::dPhidAlpha(const mat &r){
    return slater->dPsidAlpha(r);
}

double WaveFunction::dPhidBeta(const mat &r){
    if(includeJastrow){
        return jastrow->dPsidBeta(r);
    }
    else{
        return 0;
    }
}

double WaveFunction::d2PhidAlpha2(const mat &r){
    return slater->d2PsidAlpha2(r);
}

double WaveFunction::d2PhidBeta2(const mat &r){
    if(includeJastrow){
        return jastrow->d2PsidBeta2(r);
    }
    else{
        return 0;
    }
}


void WaveFunction::oneBodyDensity1D(string fn, double width, int nPoints, int nSamples, int dir){
    //1d 1b dens. dir = 0 for x axis etc.
    vec axis = linspace(-width, width, nPoints);
    vec values = zeros(nPoints);
    mat r = randu(nParticles, nDimensions);
    double factor = pow(2*width, nParticles*nDimensions-1)/nSamples;
    double sum = 0;
    double a;
    for(int i=0; i<nPoints; i++){
        sum = 0;
        for(int j=0; j<nSamples; j++){
            r = -width + 2*width*randu(nParticles, nDimensions);
            r(0,dir) = axis(i);
            slater->allNewPosUpdate(r);
            a = waveFunction(r);
            sum += a*a;
        }
        sum *= factor;
        values(i) = sum;

    }
    string s;
    stringstream out;
    out << dir;
    s = out.str();
    string dataPath = "results/oneBody/" + myAtom->getAtomName() + "res1d" + s +".mat";
    values.save(dataPath);

}

void WaveFunction::oneBodyDensity2D(string fn, double width, int nPoints, int nSamples, int dir){
    //2d 1b dens. dir = 0 for if yz plane etc.
    int dir1;
    int dir2;
    if(dir == 0){
        dir1 = 1;
        dir2 = 2;
    }
    else if(dir==1){
        dir1 = 0;
        dir2 = 2;
    }
    else{
        dir1 = 0;
        dir2 = 1;
    }

    vec axis1 = linspace(-width, width, nPoints);
    vec axis2 = linspace(-width, width, nPoints);
    mat values = zeros(nPoints, nPoints);
    mat r = randu(nParticles, nDimensions);
    double factor = pow(2*width, nParticles*nDimensions-2)/nSamples;
    double sum = 0;
    double a;
    for(int i=0; i<nPoints; i++){
        for(int k=0; k<nPoints; k++){
            sum = 0;
            for(int j=0; j<nSamples; j++){

                r = -width + 2*width*randu(nParticles, nDimensions);
                r(0,dir1) = axis1(i);
                r(0,dir2) = axis2(k);
                slater->allNewPosUpdate(r);
                a = waveFunction(r);
                sum += a*a;
            }
             sum *= factor;
             values(i,k) = sum;
        }
    }
    string s;
    stringstream out;
    out << dir;
    s = out.str();
    string dataPath = "results/oneBody/" + myAtom->getAtomName() + "res2d" + s +".mat";
    values.save(dataPath);

}

void WaveFunction::oneBodyDensity3D(string fn, double width, int nPoints, int nSamples){
    //2d 1b dens. dir = 0 for if yz plane etc.
    vec axis1 = linspace(-width, width, nPoints);
    vec axis2 = linspace(-width, width, nPoints);
    vec axis3 = linspace(-width, width, nPoints);
    cube values = zeros(nPoints, nPoints, nPoints);
    mat r = randu(nParticles, nDimensions);
    double factor = pow(2*width, nParticles*nDimensions-2)/nSamples;
    double sum = 0;
    double a;
    for(int i=0; i<nPoints; i++){
        cout << i << endl;
        for(int k=0; k<nPoints; k++){
            for(int h=0; h<nPoints; h++){

            sum = 0;
                for(int j=0; j<nSamples; j++){

                    r = -width + 2*width*randu(nParticles, nDimensions);
                    r(0,0) = axis1(i);
                    r(0,1) = axis2(k);
                    r(0,2) = axis3(h);
                    slater->allNewPosUpdate(r);
                    a = waveFunction(r);
                    sum += a*a;
                }
                sum *= factor;
                values(i,k,h) = sum;
            }
        }
    }
    string dataPath = "results/oneBody/" + myAtom->getAtomName() + "res3d.mat";
    values.save(dataPath);
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

    oneBody.print("onebody:");
}
