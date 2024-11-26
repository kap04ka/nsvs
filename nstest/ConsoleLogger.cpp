#include "ConsoleLogger.h"
#include <ctime>
#include <iomanip>

ConsoleLogger::ConsoleLogger() : currentLogLevel(LogLevel::INFO) {}

void ConsoleLogger::log(const std::string& message, LogLevel level) {
    if (level >= currentLogLevel) {
        // ��������� �������� �������
        auto now = std::time(nullptr);
        struct tm tm;

        // ���������� ���������� ������� localtime_s
        if (localtime_s(&tm, &now) == 0) {
            std::cout << std::put_time(&tm, "[%Y-%m-%d %H:%M:%S] ") << std::endl;
        }
        else {
            std::cerr << "������: ���������� �������� ��������� �����." << std::endl;
        }

        switch (level) {
        case LogLevel::INFO:
            std::cout << "\033[32m[INFO]\033[0m "; // ������� ����
            break;
        case LogLevel::WARNING:
            std::cout << "\033[33m[WARNING]\033[0m "; // ������ ����
            break;
        case LogLevel::ERROR:
            std::cout << "\033[31m[ERROR]\033[0m "; // ������� ����
            break;
        }

        std::cout << message << std::endl;
    }
}

void ConsoleLogger::setLogLevel(LogLevel level) {
    currentLogLevel = level;
}