#include "SimulationThread.hpp"
#include <QElapsedTimer>
#include <sstream>
#include <iomanip>
#include <atomic>
#include <iostream>

SimulationThread::SimulationThread(double x, double y, double vx, double vy,
                                 double theta, double omega, double duration, double dt,
                                 QObject *parent)
    : QThread(parent), x_(x), y_(y), vx_(vx), vy_(vy), theta_(theta), omega_(omega),
      duration_(duration), dt_(dt), stopRequested_(false)
{
}

void SimulationThread::stopSimulation()
{
    stopRequested_ = true;
}

void SimulationThread::run()
{
    try {
        FishTailSimulation simulation;
        FishTailSimulation::State initial_state(x_, y_, vx_, vy_, theta_, omega_, 0.0);
        
        QElapsedTimer timer;
        timer.start();
        
        auto results_ = simulation.simulate(duration_, dt_, initial_state);
        simulation.run_gnuplot();
        
        if (stopRequested_) {
            emit resultsReady("Симуляция остановлена пользователем");
            return;
        }
        
        // Prepare results text
        std::stringstream ss;
        ss << "Симуляция завершена!\n";
        ss << "Время выполнения: " << timer.elapsed() << " мс\n";
        ss << "Количество шагов: " << results_.size() << "\n\n";
        
        if (!results_.empty()) {
            const auto& last = results_.back();
            ss << "Итоговое состояние:\n";
            ss << "  Позиция: X = " << std::fixed << std::setprecision(4) << last.x 
               << ", Y = " << last.y << "\n";
            ss << "  Скорость: Vx = " << last.vx << ", Vy = " << last.vy << "\n";
            ss << "  Угол: θ = " << last.theta << " рад\n";
            ss << "  Пройденное расстояние: " << std::sqrt(last.x*last.x + last.y*last.y) << "\n";
        }
        
        // Сохраняем данные из основного потока
        QMetaObject::invokeMethod(this, [&simulation, &results_]() {
            simulation.save_data_to_file(results_, "trajectory_data.txt");
        }, Qt::BlockingQueuedConnection);
        
        emit resultsReady(QString::fromStdString(ss.str()));
        
    } catch (const std::exception& e) {
        emit simulationError(QString("Ошибка симуляции: %1").arg(e.what()));
    } catch (...) {
        emit simulationError("Неизвестная ошибка во время симуляции");
    }
}