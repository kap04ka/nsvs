#pragma once
#ifndef BASICSOLVER_HPP
#define BASICSOLVER_HPP

#include <memory>
#include <vector>
#include <string>
#include "ILogger.h"

class BasicSolver {
public:
    BasicSolver(std::shared_ptr<ILogger> logger);
    virtual void setDimensions(double width, double height, double step);
    virtual void setParams(double density, double kinematic_viscosity);
    virtual ~BasicSolver() = default;
    virtual void solve(double time, double tau, double u_max) = 0;
    virtual void printResults(const std::string& pathname_u, const std::string& pathname_v);
protected:
    std::shared_ptr<ILogger> logger;
    double width;  // Длина РО
    double height; // Высота РО
    double step; // Шаг по РО
    int nx; // Кол-во точек РО по длинне
    int ny; // Кол-во точек РО по высоте
    double density; // Плотность
    double kinematic_viscosity; // Кинематическая вязкость
    std::vector<std::vector<double>> turbulence_viscosity; // Турбулентная вязкость
    std::vector<std::vector<double>> u; // Горизонтальная компонента скорости
    std::vector<std::vector<double>> v; // Вертикальная компонента скорости
};

#endif // BASICSOLVER_HPP