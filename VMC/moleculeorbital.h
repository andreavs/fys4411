#ifndef MOLECULEORBITAL_H
#define MOLECULEORBITAL_H

class Orbital;
class Atom;
class MoleculeOrbital : public Orbital
{
public:
    MoleculeOrbital(int n, int l, int m, int nD, int nP, double a, Atom *atom, int p) : Orbital(n, l, m, nD, nP, a){
    R = atom->R;
    this->p = p;
    }
    double waveFunction(const vec &r);
    vec gradient(const vec &r);
    double laplace(const vec &r);
    double dPhidAlpha(const vec &r);
    double d2PhidAlpha2(const vec &r);
    int p;



private:
    vec3 R;


};


#endif // MOLECULEORBITAL_H
