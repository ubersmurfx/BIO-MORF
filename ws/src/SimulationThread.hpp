#ifndef SIMULATIONTHREAD_HPP
#define SIMULATIONTHREAD_HPP

#include <QThread>
#include <QString>
#include "simulation/fish_simulation.hpp"
#include "simulation_parameters.hpp"
#include <iostream>

class SimulationThread : public QThread
{
    Q_OBJECT

public:
    explicit SimulationThread(const SimulationParameters& params, QObject* parent = nullptr);
    void stopSimulation();

protected:
    void run() override;

private:
    SimulationParameters m_params;
    bool m_stopRequested;
    QString processSimulationResults(const std::vector<FishTailSimulation::State>& results);

signals:
    void progressUpdated(int progress);
    void resultsReady(const QString& results);
    void simulationFinished();
    void simulationError(const QString& error);
};

#endif // SIMULATIONTHREAD_HPP