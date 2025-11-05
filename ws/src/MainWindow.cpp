#include "MainWindow.hpp"
#include "SimulationThread.hpp"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QMessageBox>
#include <QHeaderView>
#include <QFileDialog>
#include <QDir>

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
    setupLayout->addWidget(createParameterButtons());
    setupLayout->addWidget(createControlButtons());
    setupLayout->addStretch();
    
    // Results tab
    QWidget *resultsTabWidget = createResultsTab();
    tabWidget->addTab(setupTab, "Настройка Параметров");
    tabWidget->addTab(resultsTabWidget, "Результаты Расчетов");
    
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
    
    layout->addWidget(new QLabel("Координата X:"), 0, 0);
    layout->addWidget(xEdit, 0, 1);
    layout->addWidget(new QLabel("Координата Y:"), 1, 0);
    layout->addWidget(yEdit, 1, 1);
    layout->addWidget(new QLabel("Скорость Vx:"), 2, 0);
    layout->addWidget(vxEdit, 2, 1);
    layout->addWidget(new QLabel("Скорость Vy:"), 3, 0);
    layout->addWidget(vyEdit, 3, 1);
    layout->addWidget(new QLabel("Ориентация θ:"), 4, 0);
    layout->addWidget(thetaEdit, 4, 1);
    layout->addWidget(new QLabel("Угловая скорость ω:"), 5, 0);
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

QWidget* MainWindow::createParameterButtons()
{
    QWidget *widget = new QWidget();
    QHBoxLayout *layout = new QHBoxLayout(widget);
    
    saveParamsButton = new QPushButton("Сохранить параметры");
    loadParamsButton = new QPushButton("Загрузить параметры");
    
    layout->addWidget(saveParamsButton);
    layout->addWidget(loadParamsButton);
    layout->addStretch();
    
    connect(saveParamsButton, &QPushButton::clicked, this, &MainWindow::saveParameters);
    connect(loadParamsButton, &QPushButton::clicked, this, &MainWindow::loadParameters);
    
    return widget;
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
    resultsText->setPlainText("Результаты симуляции появятся здесь после запуска расчета.");
    
    resultsTable = new QTableWidget();
    resultsTable->setColumnCount(6);
    resultsTable->setHorizontalHeaderLabels({"Время", "X", "Y", "θ", "Vx", "Vy"});
    
    layout->addWidget(new QLabel("Результаты симуляции:"));
    layout->addWidget(resultsText);
    layout->addWidget(new QLabel("Данные по шагам:"));
    layout->addWidget(resultsTable);
    
    return widget;
}

// ПРОСТЫЕ РЕАЛИЗАЦИИ БЕЗ ЗАВИСИМОСТЕЙ ОТ FishTailSimulation
void MainWindow::updateResultsTable()
{
    // Просто очищаем таблицу
    resultsTable->setRowCount(0);
}

void MainWindow::updateResultsText(const QString& results)
{
    resultsText->setPlainText(results);
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
    
    // Простое демо заполнение таблицы (в реальном приложении нужно получать данные из simulationThread)
    resultsTable->setRowCount(20);
    for (int i = 0; i < 20; ++i) {
        double time = i * 0.5;
        resultsTable->setItem(i, 0, new QTableWidgetItem(QString::number(time, 'f', 1)));
        resultsTable->setItem(i, 1, new QTableWidgetItem(QString::number(time * 0.8, 'f', 2)));
        resultsTable->setItem(i, 2, new QTableWidgetItem(QString::number(0.1 - time * 0.01, 'f', 2)));
        resultsTable->setItem(i, 3, new QTableWidgetItem("0.00"));
        resultsTable->setItem(i, 4, new QTableWidgetItem("1.00"));
        resultsTable->setItem(i, 5, new QTableWidgetItem(QString::number(-0.005 * i, 'f', 3)));
    }
    resultsTable->resizeColumnsToContents();
}

void MainWindow::saveParameters()
{
    QMessageBox::information(this, "Сохранение", "Функция сохранения параметров будет реализована позже");
}

void MainWindow::loadParameters()
{
    QMessageBox::information(this, "Загрузка", "Функция загрузки параметров будет реализована позже");
}