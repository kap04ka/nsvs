#include "SolverFactory.h"
#include <stdexcept>

std::unique_ptr<BasicSolver> SolverFactory::createSolver(const std::string& method, std::shared_ptr<ILogger> logger) {
    if (method == "vorticity-stream") {
        return std::make_unique<VorticityStreamFunctionSolver>(logger);
    }
    else if (method == "velocity-pressure") {
        return std::make_unique<VelocityPressureSolver>(logger);
    }
    else if (method == "vorticity-stream-omp") {
        return std::make_unique<VorticityStreamFunctionSolverOMP>(logger);
    }
    else if (method == "velocity-pressure-omp") {
        return std::make_unique<VelocityPressureSolverOMP>(logger);
    }
    else if (method == "q" || method == "quit") {
        return nullptr;
    }
    else {
        throw std::invalid_argument("Invalid method.");
    }
}