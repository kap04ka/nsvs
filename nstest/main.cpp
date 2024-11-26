#include "ConsoleLogger.h"
#include "SolverFactory.h"
#include <clocale>

int main() {
    setlocale(LC_ALL, "ru");

    std::shared_ptr<ILogger> logger = std::make_shared<ConsoleLogger>();

    double
        density = 1000.0,
        kinematic_viscosity = 0.009;

    double
        Lx = 15.0,
        Ly = 10.0,
        step = 1.0;

    double
        tmax = 1.5,
        tau = 0.0001,
        u_max = 2.0;

    while (true) {
        std::cout << "Methods:\n" << "vorticity-stream\n" << "velocity-pressure\n"
            << "vorticity-stream-omp\n" << "velocity-pressure-omp\n" << "q\\quit - exit\n";
        std::string option;
        std::cin >> option;

        auto solver = SolverFactory::createSolver(option, logger);
        if (solver == nullptr) return 1;

        std::string pathname_u = option + "_result_u.txt";
        std::string pathname_v = option + "_result_v.txt";

        solver->setParams(density, kinematic_viscosity);
        solver->setDimensions(Lx, Ly, step);
        solver->solve(tmax, tau, u_max);
        solver->printResults(pathname_u, pathname_v);
    }
    return 0;
}