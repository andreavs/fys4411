#ifndef SLATERDETERMINANT_H
#define SLATERDETERMINANT_H
class Orbital;

#include <armadillo>
using namespace arma;
using namespace std;

class SlaterDeterminant
{
public:
    SlaterDeterminant(int nP, int nD, double a, const mat &r);
    double getDeterminant(const mat &r);
    double waveFunction(const mat &r);
    vec gradient(const mat &r, int i);
    double laplace(const mat &r);
    void updateAll(const mat &r);
    void updateDeterminant(const mat &r);
    void updateLaplace(const mat &r);


private:
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
    void updateInverse(const mat &r);


    vector<Orbital*> orbitals;//[nParticles];


//    typedef double(SlaterDeterminant::*wavefunc) (const vec &r);
//    wavefunc functions[10];
};

#endif // SLATERDETERMINANT_H
