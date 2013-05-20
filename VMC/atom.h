#ifndef ATOM_H
#define ATOM_H


// Atom contains information (essentially number of electrons) about
// the physical system.

#include <armadillo>

class Atom
{
public:
    Atom(std::string atomType);
    Atom();
    std::string getAtomName(){return atomName;}
    int getNParticles(){return nParticles;}
    int getNDimensions(){return nDimensions;}
    bool isMolecule;
    arma::vec R;
    double charge;
    double optimalAlpha;
    double optimalBeta;


private:
    void setParametersHelium();
    void setParametersBerylium();
    void setParametersNeon();
    void setParametersH2();
    void setParametersBe2();
    std::string atomName;
    void constructor(std::string atomType);
    int nParticles;
    int nDimensions;

    //optimal params

};

#endif // ATOM_H
