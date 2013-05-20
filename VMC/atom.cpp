#include "vmc.h"

using namespace std;
using namespace arma;

Atom::Atom()
{
    string defaultType = "He";
    constructor(defaultType);
}

Atom::Atom(string atomType)
{
    constructor(atomType);
}

void Atom::constructor(string atomType){
    atomName.assign(atomType);
    if(atomType == "He"){
        setParametersHelium();
    }
    else if(atomType == "Be"){
        setParametersBerylium();
    }
    else if(atomType == "Ne"){
        setParametersNeon();
    }
    else if(atomType == "H2"){
        setParametersH2();
    }
    else if(atomType == "Be2"){
        setParametersBe2();
    }
    else{
        cout << "Invalid atom; accepted ones are He, Be and Ne! Exiting.." << endl;
        exit(0);
    }
}

void Atom::setParametersHelium(){
    nParticles = 2;
    nDimensions = 3;
    isMolecule = false;
    charge = 2;
    optimalAlpha = 1.74;
    optimalBeta = 0.34;
}

void Atom::setParametersBerylium(){
    nParticles = 4;
    nDimensions = 3;
    isMolecule = false;
    charge = 4;
    optimalAlpha = 3.88;
    optimalBeta = 0.12;
}

void Atom::setParametersNeon(){
    nParticles = 10;
    nDimensions = 3;
    isMolecule = false;
    charge = 10;
    optimalAlpha = 10.33;
    optimalBeta = 0.073;
}

void Atom::setParametersH2(){
    nParticles = 2;
    nDimensions = 3;
    isMolecule = true;
    charge = 1;
    optimalAlpha = 1.29;
    optimalBeta = 0.39;
    R << 1.4 << 0 << 0;
}

void Atom::setParametersBe2(){
    nParticles = 8;
    nDimensions = 3;
    isMolecule = true;
    charge = 4;
    optimalAlpha = 3.91;
    optimalBeta = 0.25;
    // E = -19.41(2);
    // 4.63 // 2.46
    R << 4.63 << 0 << 0;
}
