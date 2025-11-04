#include "MainWindow.hpp"
#include "SimulationThread.hpp"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QMessageBox>
#include <QHeaderView>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), simulationThread(nullptr)
{
    setupUI();
    setWindowTitle("Симуляция движения рыбы");
    setMinimumSize(800, 600);
}

MainWindow::~MainWindow()
{
    if (simulationThread && simulationThread->isRunning()) {
        simulationThread->stopSimulation();
        simulationThread->wait();
    }
}

void MainWindow::setupUI()
{
    QWidget *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    
    QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);
    
    // Tabs
    QTabWidget *tabWidget = new QTabWidget(this);
    
    // Setup tab
    QWidget *setupTab = new QWidget();
    QVBoxLayout *setupLayout = new QVBoxLayout(setupTab);
    
    setupLayout->addWidget(createInitialConditionsGroup());
    setupLayout->addWidget(createSimulationParametersGroup());
    setupLayout->addWidget(createControlButtons());
    setupLayout->addStretch();
    
    // Results tab
    QWidget *resultsTabWidget = createResultsTab();
    tabWidget->addTab(setupTab, "Настройка");
    tabWidget->addTab(resultsTabWidget, "Результаты");
    
    mainLayout->addWidget(tabWidget);
}

QGroupBox* MainWindow::createInitialConditionsGroup()
{
    QGroupBox *group = new QGroupBox("Начальные условия");
    QGridLayout *layout = new QGridLayout;
    
    // Create input fields with default values
    xEdit = new QLineEdit("0.0");
    yEdit = new QLineEdit("0.1");
    vxEdit = new QLineEdit("1.0");
    vyEdit = new QLineEdit("0.0");
    thetaEdit = new QLineEdit("0.0");
    omegaEdit = new QLineEdit("0.0");
    
    layout->addWidget(new QLabel("X:"), 0, 0);
    layout->addWidget(xEdit, 0, 1);
    layout->addWidget(new QLabel("Y:"), 1, 0);
    layout->addWidget(yEdit, 1, 1);
    layout->addWidget(new QLabel("Vx:"), 2, 0);
    layout->addWidget(vxEdit, 2, 1);
    layout->addWidget(new QLabel("Vy:"), 3, 0);
    layout->addWidget(vyEdit, 3, 1);
    layout->addWidget(new QLabel("θ:"), 4, 0);
    layout->addWidget(thetaEdit, 4, 1);
    layout->addWidget(new QLabel("ω:"), 5, 0);
    layout->addWidget(omegaEdit, 5, 1);
    
    group->setLayout(layout);
    return group;
}

QGroupBox* MainWindow::createSimulationParametersGroup()
{
    QGroupBox *group = new QGroupBox("Параметры симуляции");
    QGridLayout *layout = new QGridLayout;
    
    durationEdit = new QLineEdit("20.0");
    dtEdit = new QLineEdit("0.1");
    
    layout->addWidget(new QLabel("Длительность (с):"), 0, 0);
    layout->addWidget(durationEdit, 0, 1);
    layout->addWidget(new QLabel("Шаг по времени (с):"), 1, 0);
    layout->addWidget(dtEdit, 1, 1);
    
    group->setLayout(layout);
    return group;
}

QWidget* MainWindow::createControlButtons()
{
    QWidget *widget = new QWidget();
    QHBoxLayout *layout = new QHBoxLayout(widget);
    
    startButton = new QPushButton("Запуск симуляции");
    stopButton = new QPushButton("Остановить");
    stopButton->setEnabled(false);
    
    progressBar = new QProgressBar();
    progressBar->setVisible(false);
    
    layout->addWidget(startButton);
    layout->addWidget(stopButton);
    layout->addWidget(progressBar);
    layout->addStretch();
    
    connect(startButton, &QPushButton::clicked, this, &MainWindow::startSimulation);
    connect(stopButton, &QPushButton::clicked, this, &MainWindow::stopSimulation);
    
    return widget;
}

QWidget* MainWindow::createResultsTab()
{
    QWidget *widget = new QWidget();
    QVBoxLayout *layout = new QVBoxLayout(widget);
    
    resultsText = new QTextEdit();
    resultsText->setReadOnly(true);
    
    resultsTable = new QTableWidget();
    resultsTable->setColumnCount(6);
    resultsTable->setHorizontalHeaderLabels({"Время", "X", "Y", "θ", "Vx", "Vy"});
    
    layout->addWidget(new QLabel("Результаты:"));
    layout->addWidget(resultsText);
    layout->addWidget(new QLabel("Данные:"));
    layout->addWidget(resultsTable);
    
    return widget;
}

void MainWindow::startSimulation()
{
    // Get parameters from UI
    double x = xEdit->text().toDouble();
    double y = yEdit->text().toDouble();
    double vx = vxEdit->text().toDouble();
    double vy = vyEdit->text().toDouble();
    double theta = thetaEdit->text().toDouble();
    double omega = omegaEdit->text().toDouble();
    double duration = durationEdit->text().toDouble();
    double dt = dtEdit->text().toDouble();
    
    // Validate inputs
    if (duration <= 0 || dt <= 0) {
        QMessageBox::warning(this, "Ошибка", "Длительность и шаг по времени должны быть положительными");
        return;
    }
    
    // Create and start simulation thread
    simulationThread = new SimulationThread(x, y, vx, vy, theta, omega, duration, dt, this);
    
    connect(simulationThread, &SimulationThread::finished, this, &MainWindow::simulationFinished);
    connect(simulationThread, &SimulationThread::resultsReady, this, &MainWindow::showResults);
    connect(simulationThread, &SimulationThread::simulationError, this, &MainWindow::showResults);
    
    startButton->setEnabled(false);
    stopButton->setEnabled(true);
    progressBar->setVisible(true);
    progressBar->setValue(0);
    
    simulationThread->start();
}

void MainWindow::stopSimulation()
{
    if (simulationThread && simulationThread->isRunning()) {
        simulationThread->stopSimulation();
        simulationThread->wait();
        simulationFinished();
    }
}

void MainWindow::simulationFinished()
{
    startButton->setEnabled(true);
    stopButton->setEnabled(false);
    progressBar->setVisible(false);
    
    if (simulationThread) {
        simulationThread->deleteLater();
        simulationThread = nullptr;
    }
}

void MainWindow::updateProgress(int value)
{
    progressBar->setValue(value);
}

void MainWindow::showResults(const QString &results)
{
    resultsText->setPlainText(results);
}