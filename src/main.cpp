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

class TaskResult{
public:
    TaskResult() {}
    TaskResult(double pointX, double result) : pointX(pointX), result(result) {}

    void set(double pointX, double result) {
        this->pointX = pointX;
        this->result = result;
    }

    double getPointX() const {
        return pointX;
    }

    double getResult() const {
        return result;
    }

private:
    double pointX;
    double result;
};

class TaskParametrs {
public:
    TaskParametrs() {}

    TaskParametrs(double rangeStart, double rangeEnd, double accuracy, double stepSize) : rangeStart(rangeStart),
                                                                                          rangeEnd(rangeEnd),
                                                                                          accuracy(accuracy),
                                                                                          stepSize(stepSize) {}


    void setTask(double rangeStart, double rangeEnd, double accuracy){
        this->rangeStart = rangeStart;
        this->rangeEnd = rangeEnd;
        this->accuracy = accuracy;
    }

    void splitParametrs(double procRangeSize) {
        ProcRank rank = -1;
        int processCount = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &processCount);

        this->rangeStart += procRangeSize * rank;
        this->rangeEnd -= procRangeSize * (processCount - rank - 1);
    }

    void setAccuracy(double accuracy) {
        TaskParametrs::accuracy = accuracy;
    }

    double rangeStart;
    double rangeEnd;
    double accuracy;
    double stepSize;

};

class Solve {
public:
    static double splitTask(TaskParametrs taskParam, int &procCount) {
        double size = taskParam.rangeEnd - taskParam.rangeStart;
        double taskSize =  size / procCount;
        return taskSize;
    }
};

class Calculating {
public:
    static TaskResult runCalculating(TaskParametrs task){
        TaskResult minimum(task.rangeStart, Task::targetFunction(task.rangeStart));

        double tempNum;
        for (double x = task.rangeStart; x < task.rangeEnd; x+=task.stepSize) {
            tempNum = Task::targetFunction(x);
            if (tempNum < minimum.getResult()) {
                minimum.set(x, tempNum);
            }
        }
        return minimum;
    }
};



void runParallel(int argc, char** argv) {
    TaskParametrs taskParametrs(-10,10,0.005, 0.000001);

    int processCount = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    double taskSize;
    taskSize = Solve::splitTask(taskParametrs, processCount);

    MPI_Bcast(&taskSize, 1, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&taskParametrs, sizeof(taskParametrs), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);

    taskParametrs.splitParametrs(taskSize);

    // Start calculating for ROOT_PROCESS
    TaskResult localMin = Calculating::runCalculating(taskParametrs);

    // Get all result
//    double* resultsX = new double[processCount];
    TaskResult results[processCount];
    MPI_Gather(&localMin, sizeof(TaskResult), MPI_BYTE, results, sizeof(TaskResult), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);

//    double pX = localMin.getPointX();
//    MPI_Gather(&localMin, sizeof(TaskResult), MPI_BYTE, results, sizeof(TaskResult), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);
//    MPI_Gather(&pX, 1, MPI_DOUBLE, resultsX, 1, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);
//    std::cout << "point x: " << resultsX[0] << std::endl;

    //std::cout << "result: " << results[0]->getResult() << std::endl;
    for (int i = 0; i < processCount; ++i) {
        printf("x: %.15lf\tresult: %.15lf\n", results[i].getPointX(), results[i].getResult());
//        std::cout << "x: " << results[i].getPointX() << "\t\tresult: " << results[i].getResult() << std::endl;
    }

}

void runSubjects() {
    TaskParametrs taskParametrs;
    double taskSize;
    MPI_Bcast(&taskSize, 1, MPI_DOUBLE, ROOT_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&taskParametrs, sizeof(taskParametrs), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);
    taskParametrs.splitParametrs(taskSize);

    // Start calculating for other process
    TaskResult localMin = Calculating::runCalculating(taskParametrs);

    TaskResult results[0];
    // Send calculating data
    MPI_Gather(&localMin, sizeof(TaskResult), MPI_BYTE, results, sizeof(TaskResult), MPI_BYTE, ROOT_RANK, MPI_COMM_WORLD);

}

int main(int argc, char**argv) {
    MPI_Init(&argc, &argv);

    ProcRank rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time = MPI_Wtime();

    if (rank == ROOT_RANK) {
        runParallel(argc, argv);
    } else {
        runSubjects();
    }

    if (rank == ROOT_RANK) {
        time = MPI_Wtime() - time;
        std::cout << "Runtime: " << time << std::endl;
    }

    MPI_Finalize();

    return 0;
}