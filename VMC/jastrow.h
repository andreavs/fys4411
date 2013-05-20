#ifndef JASTROW_H
#define JASTROW_H

class Jastrow
{
public:
    Jastrow(int nP, int nD, double beta);
    double waveFunction(const mat &r);
    vec gradient(const mat &r, int row);
    double laplace(const mat &r);
    double dPsidBeta(const mat &r);
    double d2PsidBeta2(const mat &r);

private:
    int nParticles;
    int nDimensions;
    double beta;
};

#endif // JASTROW_H
