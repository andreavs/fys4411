#include <iostream>
#include "vmc.h"
using namespace std;

int main(int nargs, char* args[])
{

    MPI_Init(&nargs, &args);
    int numprocs, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    cout << "Hello World!" << endl;

    string atomType = "H2";

    VMCSetup *mySim = new VMCSetup(1e6,atomType);
    bool incJas = true;
    bool incSelf = true;
    bool incPreComp = true;
    double alpha;
    double beta;
    VMCSetup *newSim;
    Atom *myAtom;
    WaveFunction *wave;

//    cout << "He" << endl;
//    atomType = "He"; alpha = 1.84; beta = 0.34;
//    newSim = new VMCSetup(1e7,atomType);
//    newSim->runSingleSimulation(alpha, beta, true, incSelf, incPreComp, true, true, true);
//    if(my_rank == 0){
//        newSim->runSteepestDescent(2.0,0.5,100,1,1e6,true,true,true);
//    }


//    cout << "Be" << endl;
//    atomType = "Be"; alpha = 3.34; beta = 0.12;
//    newSim = new VMCSetup(1e7,atomType);
//    newSim->runSingleSimulation(alpha, beta, inJas, incSelf, incPreComp, true, true, true);
//    if(my_rank == 0){
//        newSim->runSteepestDescent(2.0,0.5,100,1,1e5,true,true,true);
//    }


//    cout << "H2" << endl;
//    atomType = "H2"; alpha = 1.29; beta = 0.39;
//    newSim = new VMCSetup(1e6,atomType);
//    if(my_rank==0){
//        newSim->runSteepestDescent(2.0,0.5,50,1,1e5,true,true,true);
//        newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true, true, true);
//    }



//    cout << "Be2" << endl;
//    atomType = "Be2"; alpha = 3.65; beta = 0.55;
//    newSim = new VMCSetup(1e7,atomType);

//    newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true,true,true);
//    myAtom = new Atom(atomType);
//    wave = new WaveFunction(myAtom, true, true, true, alpha, beta);
//    wave->oneBodyDensity2D("test", 2, 50, 1e4);
//    wave->oneBodyDensity3D("test", 2, 50, 1e4);




    cout << "Ne" << endl;
    atomType = "Ne"; alpha = 10.0; beta = 0.104;
    newSim = new VMCSetup(1e6,atomType);
//        if(my_rank == 0){
//            newSim->runSteepestDescent(10.0,0.4,100,1,1e5,true,true,true);
//        }
    newSim->runSingleSimulation(alpha, beta, true, incSelf, incPreComp, true, true, true);
//    myAtom = new Atom(atomType);
//    wave = new WaveFunction(myAtom, true, true, true, alpha, beta);
//    wave->oneBodyDensity2D("test", 2, 50, 1e4);
//    wave->oneBodyDensity3D("test", 2, 50, 1e4);




    //newSim->runSteepestDescent(3,0.4,100,1,1e5,true, true, true);



//    Atom *atom = new Atom(atomType);
//    WaveFunction *wave = new WaveFunction(atom, true, true, true, alpha, beta);
//    wave->oneBodyDensity1D("test", 2, 200, 1e6,0);
//    wave->oneBodyDensity1D("test", 2, 200, 1e6,1);




    //newSim->runSteepestDescent(3,0.4,100,1,1e5,true, true, true);

//    atomType = "H2"; alpha = 1,29; beta = 0.39;
//    VMCSetup *newSim = new VMCSetup(1e5,atomType);
//    newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true,true);




    if(my_rank == 0){
        string dataPath = "results";
        Blocking *block = new Blocking(dataPath);
        block->doBlocking();
        string pythonPath = "python " + dataPath
                + "/plotBlocking.py "
                + dataPath + "/blocking.mat";
        system(pythonPath.c_str());
}


//    srand(1);
//    vec vec1 = randu(2);
//    vec1.print("1");

    MPI_Finalize();

    return 0;
}

