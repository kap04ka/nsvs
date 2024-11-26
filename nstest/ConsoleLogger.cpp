#include "ConsoleLogger.h"
#include <ctime>
#include <iomanip>

ConsoleLogger::ConsoleLogger() : currentLogLevel(LogLevel::INFO) {}

void ConsoleLogger::log(const std::string& message, LogLevel level) {
    if (level >= currentLogLevel) {
        // Получение текущего времени
        auto now = std::time(nullptr);
        struct tm tm;

        // Используем безопасную функцию localtime_s
        if (localtime_s(&tm, &now) == 0) {
            std::cout << std::put_time(&tm, "[%Y-%m-%d %H:%M:%S] ") << std::endl;
        }
        else {
            std::cerr << "Ошибка: невозможно получить локальное время." << std::endl;
        }

        switch (level) {
        case LogLevel::INFO:
            std::cout << "\033[32m[INFO]\033[0m "; // Зеленый цвет
            break;
        case LogLevel::WARNING:
            std::cout << "\033[33m[WARNING]\033[0m "; // Желтый цвет
            break;
        case LogLevel::ERROR:
            std::cout << "\033[31m[ERROR]\033[0m "; // Красный цвет
            break;
        }

        std::cout << message << std::endl;
    }
}

void ConsoleLogger::setLogLevel(LogLevel level) {
    currentLogLevel = level;
}