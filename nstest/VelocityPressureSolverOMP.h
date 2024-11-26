#pragma once
#include "BasicSolver.h"
#include <vector>
#include <memory>

class VelocityPressureSolverOMP : public BasicSolver {
public:
    VelocityPressureSolverOMP(std::shared_ptr<ILogger> log);
    void solve(double time, double tau, double u_max) override;
private:
    void init();
    void initializeBoundaryConditions(double u_max);
    void updateBoundaryConditions();
    void updatePressure(double tau);
    void updateVelocityU(double tau);
    void updateVelocityV(double tau);

    double c;
    std::vector<std::vector<double>> u_new, v_new;
    std::vector<std::vector<double>> pressure, div_velocity;
};