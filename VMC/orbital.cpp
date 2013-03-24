#include "vmc.h"

Orbital::Orbital(int n, int l, int m, int nD, int nP, double a)
{
    this->n = n;
    this->l = l;
    this->m = m;
    alpha = a;
    nDimensions = nD;
    nParticles = nP;

}


double Orbital::waveFunction(const vec &r){
    double r1 = 0;
    double ret;

    for(int i=0; i<nDimensions;i++){
        r1 += r(i)*r(i);
    }
    r1 = sqrt(r1);
    if(n==1){
        ret =  exp(-alpha*r1);
    }
    else if(n==2){
        if(l==0){
            return (1-alpha*r1/2)*exp(-alpha*r1/2);
        }
        else if(l==1){
            double ret = r(m+1)*alpha*exp(-alpha*r1/2);
            ret = ret;
        }
    }

    else{
        cout << "Error! tried to access unknown wavefunc" << endl;
    }
    return ret;
}

vec Orbital::gradient(const vec &r){
    double r1 = 0;
    for(int i=0; i<nDimensions;i++){
        r1 += r(i)*r(i);
    }
    r1 = sqrt(r1);
    vec F = r;
    if(n==1){
        F = -F*alpha*exp(-alpha*r1)/r1;
        return F;//waveFunction(r);
    }
    if(n==2){
        if(l==0){
            F = F*alpha*(alpha*r1-4)*exp(-alpha*r1/2)/(4*r1);
            return F;//waveFunction(r);
        }
        if(l==1){
            int ind = m+1;
            F = F*r(ind);
            F(ind) = alpha*r(ind)*r(ind)-2*r1;
            F = F*alpha*exp(-alpha*r1/2)/(2*r1);
            return F;//waveFunction(r);
        }
    }

    else{
        cout << "Error! tried to access unknown wavefunc" << endl;
    }
}

double Orbital::laplace(const vec &r){
    double r1 = 0;
    for(int i=0; i<nDimensions;i++){
        r1 += r(i)*r(i);
    }
    r1 = sqrt(r1);
    if(n==1){
        return alpha*(alpha*r1 - 2)*exp(-alpha*r1)/r1;//waveFunction(r);
    }
    if(n==2){
        if(l==0){
            return -alpha*(alpha*alpha*r1*r1 - 10*r1*alpha + 16)*exp(-alpha*r1/2)/(8*r1);//waveFunction(r);
        }
        if(l==1){
            return r(m+1)*alpha*alpha*(alpha*r1-8)*exp(-alpha*r1/2)/(4*r1);//waveFunction(r);
        }
    }

    else{
        cout << "Error! tried to access unknown wavefunc" << endl;
    }
}
