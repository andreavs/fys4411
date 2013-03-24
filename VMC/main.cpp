#include <iostream>
#include "vmc.h"
using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    string atomType = "He";
    VMCSetup *mySim = new VMCSetup(1e6,atomType);
    mySim->runBruteForceSimulation();

    //SlaterDeterminant *slat = new SlaterDeterminant(3,3,1.0);

    /*
    string lol = "He";
    Atom *myAtom = new Atom(lol);
    lol = myAtom->getAtomName();
    cout << lol << endl;
    */
    return 0;
}

