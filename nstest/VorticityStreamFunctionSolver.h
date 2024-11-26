#pragma once
#include "BasicSolver.h"
#include <vector>
#include <memory>

class VorticityStreamFunctionSolver : public BasicSolver {
public:
    VorticityStreamFunctionSolver(std::shared_ptr<ILogger> log);
    void solve(double time, double tau, double u_max) override;

protected:
    void init();
    void initializeBoundaryConditions(double u_max);
    void updateVorticityBoundaryConditions();
    void solveHelmholtzEquation(double tau);
    void updateStreamFunction();
    void calculateVelocities();

    double Q, epsilon;
    std::vector<std::vector<double>> stream_function, vorticity;
    std::vector<std::vector<double>> stream_function_new, vorticity_new;
    std::vector<std::vector<double>> obstruction;
};