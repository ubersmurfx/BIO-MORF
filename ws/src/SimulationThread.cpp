#include "SimulationThread.hpp"
#include <QMetaObject>

SimulationThread::SimulationThread(const SimulationParameters& params, QObject* parent)
    : QThread(parent), m_params(params), m_stopRequested(false)
{
}

void SimulationThread::run()
{
    m_stopRequested = false;
    
    try {
        FishTailSimulation simulation(m_params);
        FishTailSimulation::State initialState(
            m_params.initX, m_params.initY, 
            m_params.initVx, m_params.initVy,
            m_params.initTheta, m_params.initOmega
        );
        
        auto results = simulation.simulate(m_params.duration, m_params.dt, initialState);
        simulation.save_data_to_file(results, "trajectory_data.txt");
        simulation.run_gnuplot(); 
        
        
        int totalSteps = static_cast<int>(m_params.duration / m_params.dt);
        for (int step = 0; step < totalSteps && !m_stopRequested; ++step) {
            int progress = static_cast<int>((step * 100) / totalSteps);
            emit progressUpdated(progress);
            
            msleep(1);
        }
        
        if (!m_stopRequested && !results.empty()) {
            QString resultsText = processSimulationResults(results);
            emit resultsReady(resultsText);
            emit simulationFinished();
        }

    } catch (const std::exception& e) {
        emit simulationError(QString("Simulation error: %1").arg(e.what()));
    }
}

QString SimulationThread::processSimulationResults(const std::vector<FishTailSimulation::State>& results)
{
    if (results.empty()) {
        return "Симуляция не дала результатов";
    }
    
    const auto& finalState = results.back();
    
    // Анализ результатов
    double maxDepth = 0;
    double maxSpeed = 0;
    double totalDistance = 0;
    double minX = results[0].x, maxX = results[0].x;
    double minY = results[0].y, maxY = results[0].y;
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& state = results[i];
        
        double depth = std::abs(state.y);
        double speed = std::sqrt(state.vx * state.vx + state.vy * state.vy);
        
        if (depth > maxDepth) maxDepth = depth;
        if (speed > maxSpeed) maxSpeed = speed;
        if (state.x < minX) minX = state.x;
        if (state.x > maxX) maxX = state.x;
        if (state.y < minY) minY = state.y;
        if (state.y > maxY) maxY = state.y;
        
        if (i > 0) {
            const auto& prev = results[i-1];
            double distance = std::sqrt(std::pow(state.x - prev.x, 2) + std::pow(state.y - prev.y, 2));
            totalDistance += distance;
        }
    }
    
    double finalSpeed = std::sqrt(finalState.vx * finalState.vx + finalState.vy * finalState.vy);
    double finalAngleDeg = finalState.theta * 180.0 / M_PI;
    
    QString resultsText = QString(
        "=== РЕЗУЛЬТАТЫ СИМУЛЯЦИИ АНПА ===\n\n"
        "ОБЩАЯ ИНФОРМАЦИЯ:\n"
        "  Длительность симуляции: %1 с\n"
        "  Шаг времени: %2 с\n"
        "  Количество точек данных: %3\n"
        "  Общее время моделирования: %4 с\n\n"
        "ФИНАЛЬНОЕ СОСТОЯНИЕ:\n"
        "  Позиция: X=%5 м, Y=%6 м\n"
        "  Скорость: Vx=%7 м/с, Vy=%8 м/с\n"
        "  Угол: %9 рад (%10°)\n"
        "  Угловая скорость: %11 рад/с\n"
        "  Время: %12 с\n\n"
        "СТАТИСТИКА ДВИЖЕНИЯ:\n"
        "  Максимальная глубина: %13 м\n"
        "  Максимальная скорость: %14 м/с\n"
        "  Пройденное расстояние: %15 м\n"
        "  Финальная скорость: %16 м/с\n"
        "  Диапазон по X: %17 - %18 м\n"
        "  Диапазон по Y: %19 - %20 м\n\n"
        "ПАРАМЕТРЫ СИМУЛЯЦИИ:\n"
        "  Масса: %21 кг\n"
        "  Тяга: %22 Н\n"
        "  Коэффициент подъемной силы: %23\n"
        "  Демпфирование: Kx=%24, Ky=%25, Kθ=%26"
    ).arg(m_params.duration, 0, 'f', 1)
     .arg(m_params.dt, 0, 'f', 3)
     .arg(results.size())
     .arg(results.back().time, 0, 'f', 1)
     .arg(finalState.x, 0, 'f', 3)
     .arg(finalState.y, 0, 'f', 3)
     .arg(finalState.vx, 0, 'f', 3)
     .arg(finalState.vy, 0, 'f', 3)
     .arg(finalState.theta, 0, 'f', 3)
     .arg(finalAngleDeg, 0, 'f', 1)
     .arg(finalState.omega, 0, 'f', 3)
     .arg(finalState.time, 0, 'f', 1)
     .arg(maxDepth, 0, 'f', 3)
     .arg(maxSpeed, 0, 'f', 3)
     .arg(totalDistance, 0, 'f', 3)
     .arg(finalSpeed, 0, 'f', 3)
     .arg(minX, 0, 'f', 3)
     .arg(maxX, 0, 'f', 3)
     .arg(minY, 0, 'f', 3)
     .arg(maxY, 0, 'f', 3)
     .arg(m_params.mass, 0, 'f', 2)
     .arg(m_params.thrust, 0, 'f', 2)
     .arg(m_params.lambda, 0, 'f', 2)
     .arg(m_params.kx, 0, 'f', 1)
     .arg(m_params.ky, 0, 'f', 1)
     .arg(m_params.ktheta, 0, 'f', 1);
    
    return resultsText;
}

void SimulationThread::stopSimulation()
{
    m_stopRequested = true;
}