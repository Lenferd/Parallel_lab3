#include <iostream>
#include <cmath>
#include "mpi.h"

typedef int ProcRank;
static const ProcRank ROOT_RANK = 0;
namespace Task{
    static double targetFunction(double x) {
        // x^4 + 3*x^3 + 6*x^3 + 8
//        return x*x*x*x + 6*x*x*x + 8;
        return x*x;
    }
    static double restriction1(double x) {
        return 2*x;
    }
    static double restriction2(double x) {
        return x * sin(x);
    }
}

struct TaskParametrs {
    double rangeStart;
    double rangeEnd;
    double accuracy;
};


void runParallel(int argc, char** argv) {
    TaskParametrs taskPar;
    taskPar.rangeStart = -10;
    taskPar.rangeEnd = 10;
    taskPar.accuracy = 0.005;


}

int main(int argc, char**argv) {
    MPI_Init(&argc, &argv);

    ProcRank rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time = MPI_Wtime();

    runParallel(argc, argv);

    if (rank == ROOT_RANK) {
        time = MPI_Wtime() - time;
        std::cout << "Runtime: " << time << std::endl;
    }

    MPI_Finalize();

    return 0;
}