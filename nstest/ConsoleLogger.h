#pragma once
#ifndef CONSOLELOGGER_HPP
#define CONSOLELOGGER_HPP

#include "ILogger.h"
#include <iostream>
#include <string>

class ConsoleLogger : public ILogger {
private:
    LogLevel currentLogLevel;

    // Можно добавить константы для цветов
    static constexpr const char* GREEN = "\033[32m";
    static constexpr const char* YELLOW = "\033[33m";
    static constexpr const char* RED = "\033[31m";
    static constexpr const char* RESET = "\033[0m";

public:
    ConsoleLogger();
    void log(const std::string& message, LogLevel level) override;
    void setLogLevel(LogLevel level) override;
};

#endif // CONSOLELOGGER_HPP