#ifndef ORBITAL_H
#define ORBITAL_H

class Orbital
{
public:
    Orbital(int n, int l, int m, int nD, int nP, double a);
    virtual double waveFunction(const vec &r);
    virtual vec gradient(const vec &r);
    virtual double laplace(const vec &r);
    virtual double dPhidAlpha(const vec &r);
    virtual double d2PhidAlpha2(const vec &r);


protected:
    double alpha;
    int nDimensions;
    int nParticles;
    int n;
    int l;
    int m;


};



#endif // ORBITAL_H
