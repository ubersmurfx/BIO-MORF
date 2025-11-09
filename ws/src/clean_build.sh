#!/bin/bash

echo "🧹 Удаление сборки"


# Удаляем файлы
rm -f fish_simulation 2>/dev/null && echo "✅ Удален: fish_simulation"
rm -f trajectory_data.txt 2>/dev/null && echo "✅ Удален: trajectory_data.txt"
rm -f trajectory_plots.png 2>/dev/null && echo "✅ Удален: trajectory_plots.png"
rm -f 3d_trajectory.png 2>/dev/null && echo "✅ Удален: 3d_trajectory.png"
rm -f plot_trajectory.gnu 2>/dev/null && echo "✅ Удален: plot_trajectory.gnu"
rm -f plot_3d.gnu 2>/dev/null && echo "✅ Удален: plot_3d.gnu"
rm -f fish_simulation.o 2>/dev/null && echo "✅ Удален: fish_simulation.o"
rm -f FishSimulation 2>/dev/null && echo "✅ Удален: FishSimulation"
rm -f FishSimulation 2>/dev/null && echo "✅ Удален: FishSimulation"
rm -f FishSimulation.pro 2>/dev/null && echo "✅ Удален: FishSimulation.pro"
rm -f main.o 2>/dev/null && echo "✅ Удален: main.o"
rm -f -r .vscode 2>/dev/null && echo "✅ Удален: .vscode"
rm -f SimulationThread.o 2>/dev/null && echo "✅ Удален: SimulationThread.o"
rm -f moc_SimulationThread.o 2>/dev/null && echo "✅ Удален: moc_SimulationThread.o"
rm -f moc_MainWindow.o 2>/dev/null && echo "✅ Удален: moc_MainWindow.o"
rm -f MainWindow.o 2>/dev/null && echo "✅ Удален: MainWindow.o"


echo "🎯 Удаление завершено!"