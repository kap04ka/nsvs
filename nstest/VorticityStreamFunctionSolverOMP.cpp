#include "VorticityStreamFunctionSolverOMP.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <omp.h>

VorticityStreamFunctionSolverOMP::VorticityStreamFunctionSolverOMP(std::shared_ptr<ILogger> logger)
    : BasicSolver(std::move(logger)) {}

void VorticityStreamFunctionSolverOMP::init() {
    logger->log("start init", LogLevel::INFO);
    u.resize(nx, std::vector<double>(ny, 0.0));
    v.resize(nx, std::vector<double>(ny, 0.0));
    stream_function.resize(nx, std::vector<double>(ny, 0.0));
    stream_function_new.resize(nx, std::vector<double>(ny, 0.0));
    vorticity.resize(nx, std::vector<double>(ny, 0.0));
    vorticity_new.resize(nx, std::vector<double>(ny, 0.0));
    obstruction.resize(nx, std::vector<double>(ny, 0.0));

#pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        obstruction[i][0] = 1;
        obstruction[i][ny - 1] = 1;
    }
    logger->log("init success", LogLevel::INFO);
}

void VorticityStreamFunctionSolverOMP::solve(double time, double tau, double u_max) {
    double t{};
    initializeBoundaryConditions(u_max);
    logger->log("Start solving", LogLevel::INFO);
    do {
        updateVorticityBoundaryConditions();
        solveHelmholtzEquation(tau);
        updateStreamFunction();
        calculateVelocities();
        t += tau;
    } while (t <= time);
    logger->log("End solving", LogLevel::INFO);
}

void VorticityStreamFunctionSolverOMP::initializeBoundaryConditions(double u_max) {
    init();

    // Параллелизация граничных значений для скорости
#pragma omp parallel for
    for (int j = 1; j < ny - 1; j++) {
        u[0][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
        u[nx - 1][j] = 4 * u_max * j * (ny - 1 - j) / ((ny - 1) * (ny - 1));
    }

#pragma omp parallel for
    for (int i = 0; i < nx; i++) {
        stream_function[i][0] = 0.0; // Нижняя граница
    }

    for (int j = 1; j < ny - 1; j++) {
        stream_function[0][j] = stream_function[0][j - 1] + u[0][j] * step;
        stream_function[nx - 1][j] = stream_function[nx - 1][j - 1] + u[nx - 1][j] * step;
    }

    stream_function[0][ny - 1] = stream_function[0][ny - 2];
    for (int i = 1; i < nx; i++) {
        stream_function[i][ny - 1] = stream_function[i - 1][ny - 1];
    }

    Q = stream_function[0][ny - 1] - stream_function[0][0];

    epsilon = 0.05;
    logger->log("init boundary conditions success", LogLevel::INFO);
    std::cout << Q << std::endl;
}

void VorticityStreamFunctionSolverOMP::updateVorticityBoundaryConditions() {
#pragma omp parallel for
    for (int i = 0; i < nx; ++i) {
        vorticity[i][ny - 1] = vorticity_new[i][ny - 1] = -(stream_function[i][ny - 1] - stream_function[i][ny - 2]) / (step * step);
        vorticity[i][0] = vorticity_new[i][0] = -(stream_function[i][1] - stream_function[i][0]) / (step * step);
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            if (obstruction[i][j] == 1) {
                if (obstruction[i][j - 1] == 0)
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j] - stream_function[i][j - 1]) / (step * step);
                if (obstruction[i][j + 1] == 0)
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j + 1] - stream_function[i][j]) / (step * step);
                if (obstruction[i - 1][j] == 0)
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i][j] - stream_function[i - 1][j]) / (step * step);
                if (obstruction[i + 1][j] == 0)
                    vorticity[i][j] = vorticity_new[i][j] = -(stream_function[i + 1][j] - stream_function[i][j]) / (step * step);
            }
        }
    }
}

void VorticityStreamFunctionSolverOMP::solveHelmholtzEquation(double tau) {
#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            if (obstruction[i][j] == 0) {
                vorticity_new[i][j] = vorticity[i][j] + tau * (
                    -((u[i][j] + fabs(u[i][j])) / 2.0) * (vorticity[i][j] - vorticity[i - 1][j]) / step
                    - ((u[i][j] - fabs(u[i][j])) / 2.0) * (vorticity[i + 1][j] - vorticity[i][j]) / step
                    - ((v[i][j] + fabs(v[i][j])) / 2.0) * (vorticity[i][j] - vorticity[i][j - 1]) / step
                    - ((v[i][j] - fabs(v[i][j])) / 2.0) * (vorticity[i][j + 1] - vorticity[i][j]) / step
                    + kinematic_viscosity * ((vorticity[i + 1][j] + vorticity[i - 1][j] + vorticity[i][j + 1] + vorticity[i][j - 1] - 4 * vorticity[i][j]) / (step * step))
                    );
            }
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i)
        for (int j = 1; j < ny - 1; ++j)
            vorticity[i][j] = vorticity_new[i][j];
}

void VorticityStreamFunctionSolverOMP::updateStreamFunction() {
    double eps;
    do {
        eps = -1.0;
#pragma omp parallel for reduction(max:eps) collapse(2)
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                if (obstruction[i][j] == 0) {
                    stream_function_new[i][j] = (step * step * vorticity[i][j] + stream_function[i + 1][j] + stream_function[i - 1][j] + stream_function[i][j + 1] + stream_function[i][j - 1]) / 4.0;
                    eps = std::max(eps, fabs(stream_function_new[i][j] - stream_function[i][j]));
                }
            }
        }

#pragma omp parallel for collapse(2)
        for (int i = 1; i < nx - 1; ++i)
            for (int j = 1; j < ny - 1; ++j)
                if (obstruction[i][j] == 0)
                    stream_function[i][j] = stream_function_new[i][j];
    } while (eps >= epsilon);
}

void VorticityStreamFunctionSolverOMP::calculateVelocities() {
#pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            if (obstruction[i][j] == 0) {
                u[i][j] = (stream_function[i + 1][j + 1] + stream_function[i - 1][j + 1] - stream_function[i + 1][j - 1] - stream_function[i - 1][j - 1]) / (4 * step);
                v[i][j] = -(stream_function[i + 1][j + 1] - stream_function[i - 1][j + 1] + stream_function[i + 1][j - 1] - stream_function[i - 1][j - 1]) / (4 * step);
            }
        }
    }
}