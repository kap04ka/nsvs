#include "BasicSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>

BasicSolver::BasicSolver(std::shared_ptr<ILogger> logger) : logger(logger) {}

void BasicSolver::setParams(double density, double kinematic_viscosity) {
    this->density = density;
    this->kinematic_viscosity = kinematic_viscosity;
}

void BasicSolver::setDimensions(double width, double height, double step) {
    this->width = width;
    this->height = height;
    this->step = step;
    this->nx = static_cast<int> (this->width / this->step) + 1;
    this->ny = static_cast<int> (this->height / this->step) + 1;
    logger->log("setDimensions passed", LogLevel::INFO);
}

void BasicSolver::printResults(const std::string& pathname_u, const std::string& pathname_v) {
    logger->log("Out u", LogLevel::WARNING);
    std::ofstream out_u;
    out_u.open(pathname_u);
    if (!out_u.is_open())  logger->log("Error open file", LogLevel::ERROR);
    for (int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            out_u << std::fixed << std::setprecision(6) << u[i][j] << "\t\t";
        }
        out_u << std::endl;
    }
    out_u.close();
    std::cout << std::endl;

    logger->log("Out v", LogLevel::WARNING);
    std::ofstream out_v;
    out_v.open(pathname_v);
    if (!out_v.is_open())  logger->log("Error open file", LogLevel::ERROR);
    for (int j = ny - 1; j >= 0; --j) {
        for (int i = 0; i < nx; ++i) {
            out_v << std::fixed << std::setprecision(6) << v[i][j] << "\t\t";
        }
        out_v << std::endl;
    }
    out_v.close();
    std::cout << std::endl;
}