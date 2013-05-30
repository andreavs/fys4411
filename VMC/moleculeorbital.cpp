#include "vmc.h"

//MoleculeOrbital::MoleculeOrbital(int n, int l, int m, int nD, int nP, double a)
//{
//    this->n = n;
//    this->l = l;
//    this->m = m;
//    alpha = a;
//    nDimensions = nD;
//    nParticles = nP;
//    R << 1.4 << 0 << 0;

//}


double MoleculeOrbital::waveFunction(const vec &r){
    //cout << "lol" << endl;
    double r1p1 = 0;
    double r1p2 = 0;
    double ret = 0;
    double ret1 = 0;
    double ret2 = 0;


    for(int i=0; i<nDimensions;i++){
        r1p1 += (r(i)+R(i)/2)*(r(i)+R(i)/2);
        r1p2 += (r(i)-R(i)/2)*(r(i)-R(i)/2);
    }
    r1p1 = sqrt(r1p1);
    r1p2 = sqrt(r1p2);
    if(n==1){
         ret1 =  exp(-alpha*r1p1);
         ret2 =  exp(-alpha*r1p2);

    } else if(n==2){
        if(l==0){
            ret1 = (1-alpha*r1p1/2)*exp(-alpha*r1p1/2);
            ret2 =  (1-alpha*r1p2/2)*exp(-alpha*r1p2/2);
        }else if(l==1){
                ret1 = (r(m+1)+R(m+1)/2)*alpha*exp(-alpha*r1p1/2);
                ret2 = (r(m+1)-R(m+1)/2)*alpha*exp(-alpha*r1p2/2);
        }
    }

    else{
        cout << "Error! tried to access unknown wavefunction " << n << endl;
    }
    if(p==1){
        ret = ret1+ret2;
    }
    else if(p==0){
        ret = ret1-ret2;
    }
     //R.print("R:");
    return ret;
}

vec MoleculeOrbital::gradient(const vec &r){
    double r1p1 = 0;
    double r1p2 = 0;
    vec ret1 = r + R/2;
    vec ret2 = r - R/2;

    for(int i=0; i<nDimensions;i++){
        r1p1 += ret1(i)*ret1(i);
        r1p2 += ret2(i)*ret2(i);
    }
    r1p1 = sqrt(r1p1);
    r1p2 = sqrt(r1p2);
    if(n==1){
        if(p==1){
            ret1 =  -alpha*(ret1*exp(-alpha*r1p1)/r1p1 + ret2*exp(-alpha*r1p2)/r1p2);
        }
        else if(p==0){
            ret1 =  -alpha*(ret1*exp(-alpha*r1p1)/r1p1 - ret2*exp(-alpha*r1p2)/r1p2);
        }
    }
    else if(n==2){
        if(l==0){
            if(p==1){
                ret1 = ret1*alpha*(alpha*r1p1-4)*exp(-alpha*r1p1/2)/(4*r1p1) + ret2*alpha*(alpha*r1p2-4)*exp(-alpha*r1p2/2)/(4*r1p2);
                return ret1;//waveFunction(r);
            } else if(p==0){
                ret1 = ret1*alpha*(alpha*r1p1-4)*exp(-alpha*r1p1/2)/(4*r1p1) - ret2*alpha*(alpha*r1p2-4)*exp(-alpha*r1p2/2)/(4*r1p2);
                return ret1;//waveFunction(r);
            }
        }
        else if(l==1){

            int ind = m+1;
            double reff1 = r(ind)+R(ind)/2;
            double reff2 = r(ind)-R(ind)/2;
            ret1 = -alpha*alpha*ret1*reff1;
            ret1(ind) = -alpha*(alpha*reff1*reff1-2*r1p1);
            ret1 = ret1*exp(-alpha*r1p1/2)/(2*r1p1);

            ret2 = -alpha*alpha*ret2*reff2;
            ret2(ind) = -alpha*(alpha*reff2*reff2-2*r1p2);
            ret2 = ret2*exp(-alpha*r1p2/2)/(2*r1p2);
            if(p==1){
                return ret1 + ret2;//waveFunction(r);
            }
            else if(p==0){
                return ret1 - ret2;
            }
        }
    }
    else{
        cout << "Error! tried to access unknown wavefunct" << endl;
    }

    //cout << "lol " << endl;

//    ret1.zeros();
//    double h = 1e-8;
//    vec rPlus = r;
//    vec rMinus = r;
//    for(int i=0; i<nDimensions; i++){
//        rPlus(i) += h;
//        rMinus(i) -= h;
//        ret1(i) = waveFunction(rPlus) - waveFunction(rMinus);
//        rPlus(i) -= h;
//        rMinus(i) += h;
//    }
//    ret1 = ret1/(2*h);

    return ret1;
}

double MoleculeOrbital::laplace(const vec &r){
    double r1p1 = 0;
    double r1p2 = 0;
    vec ret1 = r + R/2;
    vec ret2 = r - R/2;
    double ret = 0;
    for(int i=0; i<nDimensions;i++){
        r1p1 += ret1(i)*ret1(i);
        r1p2 += ret2(i)*ret2(i);
    }
    r1p1 = sqrt(r1p1);
    r1p2 = sqrt(r1p2);
    double retd1 = 0;
    double retd2 = 0;
    if(n==1){
        retd1 += alpha*(alpha*r1p1 - 2)*exp(-alpha*r1p1)/r1p1;//waveFunction(r);
        retd2 += alpha*(alpha*r1p2 - 2)*exp(-alpha*r1p2)/r1p2;
    }
    else if(n==2){
        if(l==0){
            retd1 += -alpha*(alpha*alpha*r1p1*r1p1 - 10*r1p1*alpha + 16)*exp(-alpha*r1p1/2)/(8*r1p1);//waveFunction(r);
            retd2 += -alpha*(alpha*alpha*r1p2*r1p2 - 10*r1p2*alpha + 16)*exp(-alpha*r1p2/2)/(8*r1p2);//waveFunction(r);

        }
        else if(l==1){
            retd1 += (r(m+1)+R(m+1)/2)*alpha*alpha*(alpha*r1p1-8)*exp(-alpha*r1p1/2)/(4*r1p1);//waveFunction(r);
            retd2 += (r(m+1)-R(m+1)/2)*alpha*alpha*(alpha*r1p2-8)*exp(-alpha*r1p2/2)/(4*r1p2);
        }
    }
    else{
        cout << "Error! tried to access unknown wavefuncti" << endl;
    }

//    ret = 0;
//    vec rPlus = r;
//    vec rMinus = r;
//    double h = 1e-8;
//    for(int i=0; i<nDimensions; i++){
//        rPlus(i) = r(i)+h;
//        rMinus(i) = r(i)-h;
//        ret += waveFunction(rPlus) + waveFunction(rMinus) - 2*waveFunction(r);
//        rPlus(i) = r(i);
//        rMinus(i) = r(i);

//    }
//    ret = ret/(h*h);


    if(p==1) ret = retd1 + retd2;
    else if(p==0) ret = retd1 - retd2;
    return ret;
}

double MoleculeOrbital::dPhidAlpha(const vec &r){
    //cout << "lol" << endl;
    double r1p1 = 0;
    double r1p2 = 0;
    double ret;

    for(int i=0; i<nDimensions;i++){
        r1p1 += (r(i)+R(i)/2)*(r(i)+R(i)/2);
        r1p2 += (r(i)-R(i)/2)*(r(i)-R(i)/2);
    }
    double ret1 = 0;
    double ret2 = 0;
    r1p1 = sqrt(r1p1);
    r1p2 = sqrt(r1p2);
    if(n==1){

        ret1 =  -r1p1*exp(-alpha*r1p1);
        ret2 = - r1p2*exp(-alpha*r1p2);
    }
    else if(n==2){
        if(l==0){
            ret1 = (-r1p1 + alpha*r1p1*r1p1/4)*exp(-alpha*r1p1/2);
            ret2 = (-r1p2 + alpha*r1p2*r1p2/4)*exp(-alpha*r1p2/2);
        }
        else if(l==1){
            ret1 = (r(m+1)+R(m+1)/2)*(1-alpha*r1p1/2)*exp(-alpha*r1p1/2);
            ret2 = (r(m+1)-R(m+1)/2)*(1-alpha*r1p2/2)*exp(-alpha*r1p2/2);
        }
    }

    else{
        cout << "Error! tried to access unknown wavefunction " << n << endl;
    }


    if(p ==1){
        ret = ret1 + ret2;
    }
    else if(p==0){
        ret = ret1-ret2;
    }
    return ret;

    cout << "Error! tried to access unknown wavefunc" << endl;

}

double MoleculeOrbital::d2PhidAlpha2(const vec &r){
    //cout << "lol" << endl;
    double r1p1 = 0;
    double r1p2 = 0;
    double ret;

    for(int i=0; i<nDimensions;i++){
        r1p1 += (r(i)+R(i)/2)*(r(i)+R(i)/2);
        r1p2 += (r(i)-R(i)/2)*(r(i)-R(i)/2);
    }
    r1p1 = sqrt(r1p1);
    r1p2 = sqrt(r1p2);
    if(n==1){
        ret =  r1p1*r1p1*exp(-alpha*r1p1) - r1p2*r1p2*exp(-alpha*r1p2);
    }
    else{
        cout << "Error! tried to access unknown wavefunction " << n << endl;
    }
    return ret;

    cout << "Error! tried to access unknown wavefunc" << endl;

    cout << "Error! tried to access unknown wavefunc" << endl;

}
