#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
class SlaterDeterminant;

// The wavefunction class takes in an atom and creates different trial
// wave functions, some different properties such as gradient etc..
class Jastrow;
using namespace arma;

class WaveFunction
{
public:
    WaveFunction(Atom *atom, bool incJas, bool incSelf, bool incPreComp, double a, double b);
    double waveFunction(const mat &r);
    double localEnergy(const mat &r);
    double preComputedLocalEnergy(const mat &r);
    vec drift(const mat &r, int i);
    void oneBodyDensity();
    double getNParticles(){return nParticles;}
    double getNDimensions(){return nDimensions;}

    void setBeta(double newBeta){beta = newBeta;}
    double getBeta(){return beta;}

    void setAlpha(double newAlpha){alpha = newAlpha;}
    double getAlpha(){return alpha;}

    SlaterDeterminant *slater;
    Jastrow *jastrow;


private:
    bool includeJastrow;
    bool includeSelfColoumb;
    bool includePreComputedLocalEnergy;
    Atom *myAtom;

    //numerical differentiation things
    double h;
    double h2;

    //variational parameters
    double alpha;
    double beta;

    //physical parameters to be obtained from Atom
    int charge;
    int nParticles;
    int nDimensions;

    long idum;
    Random *rnd;


};

#endif // WAVEFUNCTION_H
