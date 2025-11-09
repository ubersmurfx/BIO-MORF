#ifndef SIMULATIONTHREAD_HPP
#define SIMULATIONTHREAD_HPP

#include <QThread>
#include "simulation/fish_simulation.hpp"
#include "simulation_parameters.hpp"

class SimulationThread : public QThread
{
    Q_OBJECT

public:
    SimulationThread(const SimulationParameters& params, QObject* parent = nullptr);
    void stopSimulation();

signals:
    void simulationFinished();
    void resultsReady(const QString &results);
    void progressUpdated(int value);
    void simulationError(const QString &error);
    
protected:
    void run() override;

private:
    SimulationParameters m_params;
    bool m_stopRequested;
    QString processSimulationResults(const std::vector<FishTailSimulation::State>& results);
};

#endif // SIMULATIONTHREAD_HPP