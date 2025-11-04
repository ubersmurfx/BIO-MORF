#!/bin/bash

echo "🔨 Компиляция Fish Simulation (Qt5)..."

# Проверяем наличие qmake
if ! command -v qmake &> /dev/null; then
    echo "❌ qmake не найден. Установите Qt5:"
    echo "   sudo apt install qt5-default qtbase5-dev"
    exit 1
fi

# Проверяем версию Qt
QT_VERSION=$(qmake -query QT_VERSION)
echo "📦 Используется Qt $QT_VERSION"

# Проверяем наличие Charts
if pkg-config --exists Qt5Charts; then
    echo "✅ Qt5 Charts доступен"
    QT_MODULES="core widgets charts"
else
    echo "⚠️  Qt5 Charts не доступен, компилируем без графиков"
    QT_MODULES="core widgets"
fi

# Создаем .pro файл
echo "📄 Создаю FishSimulation.pro..."
cat > FishSimulation.pro << EOF
QT += $QT_MODULES

CONFIG += c++11

SOURCES += \\
    main.cpp \\
    MainWindow.cpp \\
    SimulationThread.cpp \\
    simulation/fish_simulation.cpp

HEADERS += \\
    MainWindow.hpp \\
    SimulationThread.hpp \\
    simulation/fish_simulation.hpp

# Настройки компиляции
QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2

# Определение для conditional compilation
!isEmpty(QT.charts) {
    DEFINES += USE_QT_CHARTS
}
EOF

# Компиляция с qmake
echo "📦 Компиляция с qmake..."
qmake FishSimulation.pro
make

if [ $? -eq 0 ]; then
    echo "✅ Сборка завершена успешно!"
    echo "🚀 Запуск программы: ./FishSimulation"
else
    echo "❌ Ошибка сборки!"
    exit 1
fi