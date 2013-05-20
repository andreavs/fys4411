#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H
class Orbital;
class Atom;
#include <armadillo>
using namespace arma;
using namespace std;

class SlaterDeterminant
{
public:
    SlaterDeterminant(Atom *atom, int nP, int nD, double a, const mat &r, bool preComp=true, bool molecule=false);
    double getDeterminant(const mat &r);
    double waveFunction(const mat &r);
    vec gradient(const mat &r, int i);
    double laplace(const mat &r);
    void updateAll(const mat &r, int i);
    void updateDeterminant(const mat &r, int i);
    void updateLaplace(const mat &r);
    double dPsidAlpha(const mat &r);
    double d2PsidAlpha2(const mat &r);
    void allNewPosUpdate(const mat &r);


private:
    bool includePrecomputed;
    int nParticles;
    int nDimensions;
    mat slater;
    mat slaterUp;
    mat slaterDown;
    mat slaterUpInverse;
    mat slaterDownInverse;
    vec slaterDownGradient;
    vec slaterUpGradient;
    double slaterUpDeterminant;
    double slaterDownDeterminant;
    double slaterUpLaplace;
    double slaterDownLaplace;
    double alpha;
    double phi1s(const vec &r);
    double phi2s(const vec &r);
    double phi2p(const vec &r, int i);
    double phi2px(const vec &r);
    double phi2py(const vec &r);
    double phi2pz(const vec &r);
    void updateMatrix(const mat &r);
    void updateInverse(const mat &r, int i);

    double ratioUp;
    double ratioDown;

    mat oldSlaterUp;
    mat oldSlaterDown;


    vector<Orbital*> orbitals;//[nParticles];

    mat dPhidAlphaUp;
    mat dPhidAlphaDown;
    mat d2PhidAlpha2Up;
    mat d2PhidAlpha2Down;
//    typedef double(SlaterDeterminant::*wavefunc) (const vec &r);
//    wavefunc functions[10];
};

#endif // SLATERDETERMINANT_H
