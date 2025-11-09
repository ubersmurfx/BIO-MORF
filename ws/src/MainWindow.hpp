#ifndef MAINWINDOW_HPP
#define MAINWINDOW_HPP

#include <QMainWindow>
#include <QPushButton>
#include <QLineEdit>
#include <QLabel>
#include <QTextEdit>
#include <QGroupBox>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QProgressBar>
#include <QTabWidget>
#include <QTableWidget>

#include <vector>
#include <QStatusBar> 
#include "simulation_parameters.hpp"


class SimulationThread;
namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void startSimulation();
    void stopSimulation();
    void simulationFinished();
    void handleSimulationError(const QString &error);
    void showResults(const QString &results);
    void updateProgress(int value);
    void refreshPlot();
    void updatePlotImage();
    void saveParameters();
    void loadParameters();
    void resetParametersToDefault();
    void loadDataAndUpdateTable(const QString& filename);


private:
    void setupUI();
    QGroupBox* createInitialConditionsGroup();
    QGroupBox* createRobotParametersGroup();
    QGroupBox* createSimulationParametersGroup();
    QWidget* createParameterButtons();
    QWidget* createControlButtons();
    QWidget* createResultsTab();
    void loadParametersFromUI();
    void saveParametersToUI();
    void updateResultsTable();
    void updateResultsText(const QString& results);
    void fillResultsTable(const QVector<double>& times,
                         const QVector<double>& xPositions,
                         const QVector<double>& yPositions,
                         const QVector<double>& angles,
                         const QVector<double>& vxVelocities,
                         const QVector<double>& vyVelocities);

private:
    Ui::MainWindow *ui;
    SimulationParameters currentParams;  // Используем глобальную структуру со всеми параметрами симуляции
    SimulationThread *simulationThread;

    // UI elements
    QLineEdit *xEdit, *yEdit, *vxEdit, *vyEdit, *thetaEdit, *omegaEdit;
    QLineEdit *massEdit, *inertiaEdit, *lambdaEdit, *thrustEdit, *fYEdit, *gEdit;
    QLineEdit *raEdit, *rmgEdit, *rlEdit;
    QLineEdit *kxEdit, *kyEdit, *kthetaEdit;
    QLineEdit *yeqEdit, *u0Edit, *v0Edit;
    QLineEdit *rhoEdit, *p0Edit;
    QLineEdit *durationEdit, *dtEdit;
    
    QPushButton *saveParamsButton, *loadParamsButton, *resetParamsButton;
    QPushButton *startButton, *stopButton, *refreshPlotButton;
    QProgressBar *progressBar;
    
    QTextEdit *resultsText;
    QTableWidget *resultsTable;
    QLabel *plotLabel;
};

#endif // MAINWINDOW_HPP