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
    double width;  // ����� ��
    double height; // ������ ��
    double step; // ��� �� ��
    int nx; // ���-�� ����� �� �� ������
    int ny; // ���-�� ����� �� �� ������
    double density; // ���������
    double kinematic_viscosity; // �������������� ��������
    std::vector<std::vector<double>> turbulence_viscosity; // ������������ ��������
    std::vector<std::vector<double>> u; // �������������� ���������� ��������
    std::vector<std::vector<double>> v; // ������������ ���������� ��������
};

#endif // BASICSOLVER_HPP