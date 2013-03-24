#ifndef ATOM_H
#define ATOM_H


// Atom contains information (essentially number of electrons) about
// the physical system.



class Atom
{
public:
    Atom(std::string atomType);
    Atom();
    std::string getAtomName(){return atomName;}
    int getNParticles(){return nParticles;}
    int getNDimensions(){return nDimensions;}


private:
    void setParametersHelium();
    void setParametersBerylium();
    void setParametersNeon();
    std::string atomName;
    void constructor(std::string atomType);
    int nParticles;
    int nDimensions;

    //optimal params
    double optimalAlpha;
    double optimalBeta;
};

#endif // ATOM_H
