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
#include <QChart>
#include <QChartView>
#include <QLineSeries>
#include <QValueAxis>

QT_CHARTS_USE_NAMESPACE

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

private:
  void setupUI();
  void setupCharts();
  QGroupBox* createInitialConditionsGroup();
  QGroupBox* createSimulationParametersGroup();
  QWidget* createControlButtons();
  QWidget* createResultsTab();

  // Input fields
  QLineEdit *xEdit, *yEdit, *vxEdit, *vyEdit, *thetaEdit, *omegaEdit;
  QLineEdit *durationEdit, *dtEdit;
    
  // Controls
  QPushButton *startButton, *stopButton;
  QProgressBar *progressBar;
  QTextEdit *resultsText;
  QTableWidget *resultsTable;
    
  // Charts
  QChart *trajectoryChart, *positionChart, *velocityChart, *angleChart;
  QChartView *trajectoryView, *positionView, *velocityView, *angleView;
  
  SimulationThread *simulationThread;
};

#endif // MAINWINDOW_HPP