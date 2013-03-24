#ifndef ORBITAL_H
#define ORBITAL_H

class Orbital
{
public:
    Orbital(int n, int l, int m, int nD, int nP, double a);
    double waveFunction(const vec &r);
    vec gradient(const vec &r);
    double laplace(const vec &r);


private:
    double alpha;
    int nDimensions;
    int nParticles;
    int n;
    int l;
    int m;


};

#endif // ORBITAL_H
