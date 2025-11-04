QT += core widgets charts

CONFIG += c++11

SOURCES += \
    main.cpp \
    MainWindow.cpp \
    SimulationThread.cpp \
    simulation/fish_simulation.cpp

HEADERS += \
    MainWindow.hpp \
    SimulationThread.hpp \
    simulation/fish_simulation.hpp

# Настройки компиляции
QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2

# Определение для conditional compilation
!isEmpty(QT.charts) {
    DEFINES += USE_QT_CHARTS
}
