#include "VelocityPressureSolverOMP.h"
#include <cmath>
#include <omp.h>

VelocityPressureSolverOMP::VelocityPressureSolverOMP(std::shared_ptr<ILogger> logger)
    : BasicSolver(std::move(logger)) {}

void VelocityPressureSolverOMP::init() {
    logger->log("start init", LogLevel::INFO);
    u.resize(nx, std::vector<double>(ny, 0.0));
    v.resize(nx, std::vector<double>(ny, 0.0));
    u_new.resize(nx, std::vector<double>(ny, 0.0));
    v_new.resize(nx, std::vector<double>(ny, 0.0));
    pressure.resize(nx, std::vector<double>(ny, 0.0));
    div_velocity.resize(nx, std::vector<double>(ny, 0.0));
    logger->log("init success", LogLevel::INFO);
}

void VelocityPressureSolverOMP::solve(double time, double tau, double u_max) {
    double t{};
    initializeBoundaryConditions(u_max);
    logger->log("Start solving", LogLevel::INFO);
    do {
        updatePressure(tau);
        updateBoundaryConditions();
        updateVelocityU(tau);
        updateVelocityV(tau);
        t += tau;
    } while (t <= time);
    logger->log("End solving", LogLevel::INFO);
}

void VelocityPressureSolverOMP::initializeBoundaryConditions(double u_max) {
    init();
#pragma omp parallel for
    for (int j = 1; j < ny - 1; j++) {
        u[0][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
        u[nx - 1][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
    }

    c = 1, density = 1, kinematic_viscosity = 1;
    logger->log("init boundary conditions success", LogLevel::INFO);
}

void VelocityPressureSolverOMP::updatePressure(double tau) {
#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            div_velocity[i][j] = ((u[i + 1][j + 1] + u[i + 1][j - 1]) - (u[i - 1][j - 1] + u[i - 1][j + 1])
                + (v[i - 1][j + 1] + v[i + 1][j + 1]) - (v[i - 1][j - 1] + v[i + 1][j - 1])) / (4 * step);
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            pressure[i][j] -= c * tau * div_velocity[i][j];
        }
    }
}

void VelocityPressureSolverOMP::updateBoundaryConditions() {
#pragma omp parallel for
    for (int i = 1; i < nx - 1; ++i) {
        pressure[i][0] = pressure[i][1];
        pressure[i][ny - 1] = pressure[i][ny - 2];
    }

#pragma omp parallel for
    for (int j = 1; j < ny - 1; ++j) {
        pressure[0][j] = 2 * pressure[1][j] - pressure[2][j];
        pressure[nx - 1][j] = 2 * pressure[nx - 2][j] - pressure[nx - 3][j];
    }
}

void VelocityPressureSolverOMP::updateVelocityU(double tau) {
#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            u_new[i][j] = u[i][j] + tau * (
                -((u[i][j] + std::abs(u[i][j])) / 2 * (u[i][j] - u[i - 1][j]) / step
                    + (u[i][j] - std::abs(u[i][j])) / 2 * (u[i + 1][j] - u[i][j]) / step)
                - ((v[i][j] + std::abs(v[i][j])) / 2 * (u[i][j] - u[i][j - 1]) / step
                    + (v[i][j] - std::abs(v[i][j])) / 2 * (u[i][j + 1] - u[i][j]) / step)
                - (pressure[i + 1][j + 1] + pressure[i + 1][j - 1] - pressure[i - 1][j + 1] - pressure[i - 1][j - 1]) / (4 * step * density)
                + kinematic_viscosity * (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j]) / (step * step)
                );
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i)
        for (int j = 1; j < ny - 1; ++j)
            u[i][j] = u_new[i][j];
}

void VelocityPressureSolverOMP::updateVelocityV(double tau) {
#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            v_new[i][j] = v[i][j] + tau * (
                -((u[i][j] + std::abs(u[i][j])) / 2 * (v[i][j] - v[i - 1][j]) / step
                    + (u[i][j] - std::abs(u[i][j])) / 2 * (v[i + 1][j] - v[i][j]) / step)
                - ((v[i][j] + std::abs(v[i][j])) / 2 * (v[i][j] - v[i][j - 1]) / step
                    + (v[i][j] - std::abs(v[i][j])) / 2 * (v[i][j + 1] - v[i][j]) / step)
                - (pressure[i + 1][j + 1] + pressure[i - 1][j + 1] - pressure[i + 1][j - 1] - pressure[i - 1][j - 1]) / (4 * step * density)
                + kinematic_viscosity * (v[i + 1][j] + v[i - 1][j] + v[i][j + 1] + v[i][j - 1] - 4 * v[i][j]) / (step * step)
                );
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i)
        for (int j = 1; j < ny - 1; ++j)
            v[i][j] = v_new[i][j];
}
