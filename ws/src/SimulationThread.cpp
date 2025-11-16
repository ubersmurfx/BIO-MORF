#include "SimulationThread.hpp"
#include <QMetaObject>
#include <cmath>
#include <iostream>

SimulationThread::SimulationThread(const SimulationParameters& params, QObject* parent)
    : QThread(parent), m_params(params), m_stopRequested(false)
{
}

void SimulationThread::run()
{
    m_stopRequested = false;
    
    try {
        FishTailSimulation simulation(m_params);
        auto analysis = simulation.analyze_equilibrium_scenarios();
        std::cout << "Статическое равновесие: " << analysis.y_static << std::endl;
        std::cout << "Динамическое равновесие: " << analysis.y_dynamic << std::endl;

        // Создаем начальное состояние согласно новой модели из статьи
        // z = -y (так как в статье z направлено вверх, а y в вашей системе вниз)
        // zdot = -vy
        FishTailSimulation::State initialState(
            0.0,                            // x (будет вычисляться автоматически)
            -m_params.initY,                // z = -y (переводим в систему координат статьи)
            -m_params.initVy,               // zdot = -vy
            m_params.initTheta,             // theta (такое же)
            m_params.initOmega,             // thetadot = omega
            0.0                             // time
        );
        
        std::cout << "Starting simulation with parameters:" << std::endl;
        simulation.print_parameters();
        
        // Анализ устойчивости перед симуляцией
        simulation.is_system_stable();
        
        auto results = simulation.simulate(m_params.duration, m_params.dt, initialState);
        simulation.save_data_to_file(results, "trajectory_data.txt");
        
        // Создаем графики
        QMetaObject::invokeMethod(this, [&simulation]() {
            simulation.create_gnuplot_script();
        }, Qt::BlockingQueuedConnection);
        
        // Эмуляция прогресса (можно адаптировать под реальный прогресс)
        int totalSteps = static_cast<int>(m_params.duration / m_params.dt);
        for (int step = 0; step < totalSteps && !m_stopRequested; ++step) {
            int progress = static_cast<int>((step * 100) / totalSteps);
            emit progressUpdated(progress);
            
            // Небольшая задержка для плавного обновления прогресса
            msleep(10);
        }
        
        if (!m_stopRequested && !results.empty()) {
            QString resultsText = processSimulationResults(results);
            emit resultsReady(resultsText);
            emit simulationFinished();
            
            // Запускаем GNUplot после завершения симуляции
            simulation.run_gnuplot();
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
    
    // Анализ результатов в системе координат статьи (z вверх)
    double maxDepth = 0;
    double maxSpeed = 0;
    double totalDistance = 0;
    double minX = results[0].x, maxX = results[0].x;
    double minZ = results[0].z, maxZ = results[0].z;
    
    // Для совместимости с UI вычисляем также в старой системе координат (y вниз)
    double minY = -results[0].z, maxY = -results[0].z;
    
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& state = results[i];
        
        // В системе статьи: z - высота (положительная вверх)
        // Глубина = -z (если z отрицательно - под водой)
        double depth = std::max(0.0, -state.z);
        double speed = std::sqrt(state.vx * state.vx + state.zdot * state.zdot);
        
        if (depth > maxDepth) maxDepth = depth;
        if (speed > maxSpeed) maxSpeed = speed;
        if (state.x < minX) minX = state.x;
        if (state.x > maxX) maxX = state.x;
        if (state.z < minZ) minZ = state.z;
        if (state.z > maxZ) maxZ = state.z;
        
        // Для старой системы координат
        double y = -state.z;
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
        
        if (i > 0) {
            const auto& prev = results[i-1];
            double distance = std::sqrt(std::pow(state.x - prev.x, 2) + std::pow(state.z - prev.z, 2));
            totalDistance += distance;
        }
    }
    
    double finalSpeed = std::sqrt(finalState.vx * finalState.vx + finalState.zdot * finalState.zdot);
    double finalAngleDeg = finalState.theta * 180.0 / M_PI;
    
    // Анализ устойчивости по конечному состоянию
    bool isStable = std::abs(finalState.zdot) < 0.01 && 
                   std::abs(finalState.thetadot) < 0.01 &&
                   std::abs(finalState.theta) < 0.1;
    
    QString stabilityStatus = isStable ? "УСТОЙЧИВО" : "НЕУСТОЙЧИВО";
    
    QString resultsText = QString(
        "=== РЕЗУЛЬТАТЫ СИМУЛЯЦИИ (Модель из статьи) ===\n\n"
        "ОБЩАЯ ИНФОРМАЦИЯ:\n"
        "  Длительность симуляции: %1 с\n"
        "  Шаг времени: %2 с\n"
        "  Количество точек данных: %3\n"
        "  Статус устойчивости: %4\n\n"
        "ФИНАЛЬНОЕ СОСТОЯНИЕ (система статьи):\n"
        "  Позиция: X=%5 м, Z=%6 м\n"
        "  Скорость: Vx=%7 м/с, Ż=%8 м/с\n"
        "  Угол тангажа: θ=%9 рад (%10°)\n"
        "  Угловая скорость: θ̇=%11 рад/с\n"
        "  Время: %12 с\n\n"
        "ФИНАЛЬНОЕ СОСТОЯНИЕ (старая система):\n"
        "  Позиция: X=%5 м, Y=%13 м\n"
        "  Скорость: Vx=%7 м/с, Vy=%14 м/с\n\n"
        "СТАТИСТИКА ДВИЖЕНИЯ:\n"
        "  Максимальная глубина: %15 м\n"
        "  Максимальная скорость: %16 м/с\n"
        "  Пройденное расстояние: %17 м\n"
        "  Финальная скорость: %18 м/с\n"
        "  Диапазон по X: %19 - %20 м\n"
        "  Диапазон по Z: %21 - %22 м\n"
        "  Диапазон по Y: %23 - %24 м\n\n"
        "ПАРАМЕТРЫ СИМУЛЯЦИИ:\n"
        "  Масса: %25 кг\n"
        "  Момент инерции: %26 кг·м²\n"
        "  Скорость U: %27 м/с\n"
        "  Коэффициент подъемной силы K: %28\n"
        "  Смещение плавучести rA: %29 м\n"
        "  Смещение подъемной силы rL: %30 м"
    ).arg(m_params.duration, 0, 'f', 1)
     .arg(m_params.dt, 0, 'f', 3)
     .arg(results.size())
     .arg(stabilityStatus)
     .arg(finalState.x, 0, 'f', 3)
     .arg(finalState.z, 0, 'f', 3)
     .arg(finalState.vx, 0, 'f', 3)
     .arg(finalState.zdot, 0, 'f', 3)
     .arg(finalState.theta, 0, 'f', 3)
     .arg(finalAngleDeg, 0, 'f', 1)
     .arg(finalState.thetadot, 0, 'f', 3)
     .arg(finalState.time, 0, 'f', 1)
     .arg(-finalState.z, 0, 'f', 3)  // Y = -Z
     .arg(-finalState.zdot, 0, 'f', 3) // Vy = -Ż
     .arg(maxDepth, 0, 'f', 3)
     .arg(maxSpeed, 0, 'f', 3)
     .arg(totalDistance, 0, 'f', 3)
     .arg(finalSpeed, 0, 'f', 3)
     .arg(minX, 0, 'f', 3)
     .arg(maxX, 0, 'f', 3)
     .arg(minZ, 0, 'f', 3)
     .arg(maxZ, 0, 'f', 3)
     .arg(minY, 0, 'f', 3)
     .arg(maxY, 0, 'f', 3)
     .arg(m_params.mass, 0, 'f', 2)
     .arg(m_params.inertia, 0, 'f', 2)
     .arg(m_params.u0, 0, 'f', 2)
     .arg(m_params.lambda, 0, 'f', 2)
     .arg(m_params.ra, 0, 'f', 2)
     .arg(m_params.rl, 0, 'f', 2);
    
    return resultsText;
}

void SimulationThread::stopSimulation()
{
    m_stopRequested = true;
}