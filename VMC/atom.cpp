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
    else{
        cout << "Invalid atom; accepted ones are He, Be and Ne! Exiting.." << endl;
        exit(0);
    }
}

void Atom::setParametersHelium(){
    nParticles = 2;
    nDimensions = 3;
}

void Atom::setParametersBerylium(){
    nParticles = 4;
    nDimensions = 3;
}

void Atom::setParametersNeon(){
    nParticles = 10;
    nDimensions = 3;
}
