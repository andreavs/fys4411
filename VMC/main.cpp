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
//    newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true, true, true);
//    if(my_rank == 0){
//        newSim->runSteepestDescent(2.0,0.5,100,1,1e5,true,true,true);
//    }


    cout << "H2" << endl;
    atomType = "H2"; alpha = 1.29; beta = 0.39;
    int N = 30;
    newSim = new VMCSetup(1e7,atomType);
    newSim->atom->R(0) = 1.34;
    newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true, true, true);
//    if(my_rank == 0){
//        newSim->runSteepestDescent(2.0,0.5,100,1,1e5,true,true,true);
//    }


    //example of how to find optimal parameters for a set of R values:
//    vec Rvals = linspace(0.1,2.5,N);
//    double e;
//    vec avec = zeros(N);
//    vec bvec = zeros(N);
//    vec evals = zeros(N);
//    vec params = zeros(2);
//    if(my_rank==0){
//        for(int i=0; i<N;i++){
//            newSim->atom->R(0) = Rvals(i);
//            params = newSim->runSteepestDescent(2.0,0.5,50,1,1e5,true,true,true);
//            alpha = params(0);
//            avec(i) = alpha;
//            beta = params(1);
//            bvec(i) = beta;
//            evals(i) = newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true, true, true);
//        }

//        Rvals.print("R;");
//        evals.print("E:");
//        avec.print("alpha:");
//        bvec.print("beta:");
//    }



//    cout << "Be2" << endl;
//    atomType = "Be2"; alpha = 3.65; beta = 0.55;
//    newSim = new VMCSetup(1e7,atomType);

//    newSim->runSingleSimulation(alpha, beta, incJas, incSelf, incPreComp, true,true,true);
//    myAtom = new Atom(atomType);
//    wave = new WaveFunction(myAtom, true, true, true, alpha, beta);
//    wave->oneBodyDensity2D("test", 2, 50, 1e4);
//    wave->oneBodyDensity3D("test", 2, 50, 1e4);


    //run blocking on the results:


    if(my_rank == 0){
        string dataPath = "results";
        Blocking *block = new Blocking(dataPath);
        block->doBlocking();
        string pythonPath = "python " + dataPath
                + "/plotBlocking.py "
                + dataPath + "/blocking.mat";
        system(pythonPath.c_str());
}


    MPI_Finalize();

    return 0;
}

