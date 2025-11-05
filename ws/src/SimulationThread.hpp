#ifndef SIMULATIONTHREAD_HPP
#define SIMULATIONTHREAD_HPP

#include <QThread>
#include <QString>
#include <vector>
#include "simulation/fish_simulation.hpp"

class SimulationThread : public QThread
{
    Q_OBJECT

public:
    SimulationThread(double x, double y, double vx, double vy,
                    double theta, double omega, double duration, double dt,
                    QObject *parent = nullptr);
    
    void stopSimulation();
    const std::vector<FishTailSimulation::State>& getResults() const { return results_; }
signals:
    void progressUpdated(int percent);
    void resultsReady(const QString &results);
    void simulationError(const QString &error);

protected:
    void run() override;

private:
    double x_, y_, vx_, vy_, theta_, omega_, duration_, dt_;
    bool stopRequested_;
    std::vector<FishTailSimulation::State> results_; 
};

#endif // SIMULATIONTHREAD_HPP