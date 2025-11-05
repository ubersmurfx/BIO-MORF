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

// Только forward declaration
class SimulationThread;

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
    void updateProgress(int value);
    void showResults(const QString &results);
    void saveParameters();
    void loadParameters();

private:
    void setupUI();
    QGroupBox* createInitialConditionsGroup();
    QGroupBox* createSimulationParametersGroup();
    QWidget* createControlButtons();
    QWidget* createResultsTab();
    QWidget* createParameterButtons();
    
    // УБИРАЕМ методы с параметрами State, оставляем только простые методы
    void updateResultsTable();
    void updateResultsText(const QString& results);

    // Input fields
    QLineEdit *xEdit, *yEdit, *vxEdit, *vyEdit, *thetaEdit, *omegaEdit;
    QLineEdit *durationEdit, *dtEdit;
    
    // Controls
    QPushButton *startButton, *stopButton;
    QPushButton *saveParamsButton, *loadParamsButton;
    QProgressBar *progressBar;
    QTextEdit *resultsText;
    QTableWidget *resultsTable;
    
    SimulationThread *simulationThread;
};

#endif // MAINWINDOW_HPP