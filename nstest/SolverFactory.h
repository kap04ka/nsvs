#pragma once
#include "BasicSolver.h"
#include "VelocityPressureSolver.h"
#include "VelocityPressureSolverOMP.h"
#include "VorticityStreamFunctionSolver.h"
#include "VorticityStreamFunctionSolverOMP.h"
#include <memory>
#include <string>

class SolverFactory {
public:
    static std::unique_ptr<BasicSolver> createSolver(const std::string& method, std::shared_ptr<ILogger> logger);
};