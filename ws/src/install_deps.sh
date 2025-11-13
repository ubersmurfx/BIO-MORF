#!/bin/bash

echo "📦 Установка зависимостей для Fish Simulation..."

# Обновление пакетов
sudo apt update

# Установка Qt5 и компилятора
sudo apt install -y qt5-default qtbase5-dev qtcharts5-dev
sudo apt install -y cmake build-essential g++
sudo apt-get install gnuplot -y
# Установить QT Charts
sudo apt install libqt5charts5-dev
sudo apt-get install gnuplot

echo "✅ Зависимости установлены!"
echo "🔨 Теперь можно запустить: ./build_sim.sh"