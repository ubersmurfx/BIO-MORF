#include "MainWindow.hpp"
#include "SimulationThread.hpp"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QMessageBox>
#include <QHeaderView>
#include <QFileDialog>
#include <QDir>
#include <QPixmap>
#include <QFile>
#include <QTimer>
#include <QMessageBox>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), simulationThread(nullptr)
{
    setupUI();
    setWindowTitle("Симуляция движения рыбы");
    setMinimumSize(1200, 700);

    qRegisterMetaType<SimulationParameters>();
    
    // Таймер для автообновления графика каждые 2 секунды
    QTimer *refreshTimer = new QTimer(this);
    connect(refreshTimer, &QTimer::timeout, this, &MainWindow::updatePlotImage);
    refreshTimer->start(2000); // 2000 мс = 2 секунды
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
    QHBoxLayout *mainLayout = new QHBoxLayout(centralWidget);
    
    // Левая панель - настройки и управление
    QWidget *leftPanel = new QWidget();
    QVBoxLayout *leftLayout = new QVBoxLayout(leftPanel);
    leftPanel->setMaximumWidth(400);
    
    leftLayout->addWidget(createInitialConditionsGroup());
    leftLayout->addWidget(createRobotParametersGroup()); 
    leftLayout->addWidget(createSimulationParametersGroup());
    leftLayout->addWidget(createParameterButtons());
    leftLayout->addWidget(createControlButtons());
    leftLayout->addStretch();
    
    // Правая панель - результаты и график
    QWidget *rightPanel = createResultsTab();
    
    mainLayout->addWidget(leftPanel);
    mainLayout->addWidget(rightPanel);
}

QGroupBox* MainWindow::createInitialConditionsGroup()
{
    QGroupBox *group = new QGroupBox("Начальные условия");
    QGridLayout *layout = new QGridLayout;
    
    // Create input fields with default values
    xEdit = new QLineEdit("0.0");
    yEdit = new QLineEdit("0.4");
    vxEdit = new QLineEdit("0.2");
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
    
    durationEdit = new QLineEdit("200.0");
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
    resetParamsButton = new QPushButton("Сбросить параметры"); 

    layout->addWidget(saveParamsButton);
    layout->addWidget(loadParamsButton);
    layout->addWidget(resetParamsButton);
    layout->addStretch();
    
    connect(saveParamsButton, &QPushButton::clicked, this, &MainWindow::saveParameters);
    connect(loadParamsButton, &QPushButton::clicked, this, &MainWindow::loadParameters);
    connect(resetParamsButton, &QPushButton::clicked, this, &MainWindow::resetParametersToDefault);

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
    
    // Вкладки для результатов
    QTabWidget *resultsTabs = new QTabWidget();
    
    // Вкладка с текстовыми результатами
    QWidget *textTab = new QWidget();
    QVBoxLayout *textLayout = new QVBoxLayout(textTab);
    
    resultsText = new QTextEdit();
    resultsText->setReadOnly(true);
    resultsText->setPlainText("Результаты симуляции появятся здесь после запуска расчета.");
    
    resultsTable = new QTableWidget();
    resultsTable->setColumnCount(6);
    resultsTable->setHorizontalHeaderLabels({"Время", "X", "Y", "θ", "Vx", "Vy"});
    
    textLayout->addWidget(new QLabel("Текстовые результаты:"));
    textLayout->addWidget(resultsText);
    textLayout->addWidget(new QLabel("Данные по шагам:"));
    textLayout->addWidget(resultsTable);
    
    // Вкладка с графиком
    QWidget *plotTab = new QWidget();
    QVBoxLayout *plotLayout = new QVBoxLayout(plotTab);
    
    // Панель управления графиком
    QHBoxLayout *plotControlLayout = new QHBoxLayout();
    refreshPlotButton = new QPushButton("Обновить график");
    QLabel *plotStatusLabel = new QLabel("График будет загружен после симуляции");
    
    plotControlLayout->addWidget(refreshPlotButton);
    plotControlLayout->addWidget(plotStatusLabel);
    plotControlLayout->addStretch();
    
    // Область для отображения графика
    plotLabel = new QLabel();
    plotLabel->setAlignment(Qt::AlignCenter);
    plotLabel->setMinimumSize(600, 800);
    plotLabel->setStyleSheet("QLabel { background-color: white; border: 1px solid gray; }");
    plotLabel->setText("Запустите симуляцию для генерации графика");
    
    plotLayout->addLayout(plotControlLayout);
    plotLayout->addWidget(plotLabel);
    
    connect(refreshPlotButton, &QPushButton::clicked, this, &MainWindow::refreshPlot);
    
    resultsTabs->addTab(textTab, "Текстовые результаты");
    resultsTabs->addTab(plotTab, "График траектории");
    
    layout->addWidget(resultsTabs);
    
    return widget;
}

void MainWindow::refreshPlot()
{
    updatePlotImage();
}

void MainWindow::updatePlotImage()
{
    QString imagePath = "trajectory_plots.png";
    
    if (QFile::exists(imagePath)) {
        QPixmap pixmap(imagePath);
        if (!pixmap.isNull()) {
            // Масштабируем изображение чтобы оно вписывалось в область
            plotLabel->setPixmap(pixmap.scaled(plotLabel->width() - 20, 
                                             plotLabel->height() - 20, 
                                             Qt::KeepAspectRatio, 
                                             Qt::SmoothTransformation));
        } else {
            plotLabel->setText("Ошибка загрузки изображения");
        }
    } else {
        plotLabel->setText("Файл trajectory_plots.png не найден\nЗапустите симуляцию для генерации графика");
    }
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
    loadParametersFromUI();
    
    // Валидация
    if (currentParams.duration <= 0 || currentParams.dt <= 0) {
        QMessageBox::warning(this, "Ошибка", "Длительность и шаг по времени должны быть положительными");
        return;
    }

    if (currentParams.mass <= 0) {
        QMessageBox::warning(this, "Ошибка", "Масса должна быть положительной");
        return;
    }
    
    FishTailSimulation::State initialState = {
        currentParams.initX,      // Make sure these variables are defined
        currentParams.initY,      // and have appropriate values
        currentParams.initVx,
        currentParams.initVy,
        currentParams.initTheta,
        currentParams.initOmega
    };

    // Create and start simulation thread
    simulationThread = new SimulationThread(currentParams, this);
    
    connect(simulationThread, &SimulationThread::finished, this, &MainWindow::simulationFinished);
    connect(simulationThread, &SimulationThread::resultsReady, this, &MainWindow::showResults);
    connect(simulationThread, &SimulationThread::simulationError, this, &MainWindow::handleSimulationError);
    connect(simulationThread, &SimulationThread::progressUpdated, this, &MainWindow::updateProgress);

    startButton->setEnabled(false);
    stopButton->setEnabled(true);
    progressBar->setVisible(true);
    progressBar->setValue(0);
    
    simulationThread->start();
}

void MainWindow::handleSimulationError(const QString &error)
{
    QMessageBox::critical(this, "Ошибка симуляции", error);
    simulationFinished();
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
    
    updatePlotImage();
    
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

QGroupBox* MainWindow::createRobotParametersGroup()
{
    QGroupBox *group = new QGroupBox("Параметры АНПА");
    QGridLayout *layout = new QGridLayout;
    
    // Создаем поля ввода
    massEdit = new QLineEdit("3.0");
    inertiaEdit = new QLineEdit("0.5");
    lambdaEdit = new QLineEdit("1.5");
    thrustEdit = new QLineEdit("0.1");
    fYEdit = new QLineEdit("0.2");
    gEdit = new QLineEdit("9.81");
    raEdit = new QLineEdit("0.06");
    rmgEdit = new QLineEdit("-0.05");
    rlEdit = new QLineEdit("-0.1");
    kxEdit = new QLineEdit("2.0");
    kyEdit = new QLineEdit("1.0");
    kthetaEdit = new QLineEdit("5.0");
    yeqEdit = new QLineEdit("0.344292");
    u0Edit = new QLineEdit("1.0");
    v0Edit = new QLineEdit("0.00305");
    rhoEdit = new QLineEdit("1000.0");
    p0Edit = new QLineEdit("101325.0");
    
    int row = 0;
    
    // Основные параметры
    layout->addWidget(new QLabel("Масса (кг):"), row, 0);
    layout->addWidget(massEdit, row++, 1);
    layout->addWidget(new QLabel("Момент инерции:"), row, 0);
    layout->addWidget(inertiaEdit, row++, 1);
    layout->addWidget(new QLabel("Коэф. подъемной силы:"), row, 0);
    layout->addWidget(lambdaEdit, row++, 1);
    
    // Параметры тяги
    layout->addWidget(new QLabel("Номинальная тяга (Н):"), row, 0);
    layout->addWidget(thrustEdit, row++, 1);
    layout->addWidget(new QLabel("Коэф. Fy:"), row, 0);
    layout->addWidget(fYEdit, row++, 1);
    layout->addWidget(new QLabel("Коэф. G:"), row, 0);
    layout->addWidget(gEdit, row++, 1);
    
    // Плечи сил
    layout->addWidget(new QLabel("Плечо Архимеда (м):"), row, 0);
    layout->addWidget(raEdit, row++, 1);
    layout->addWidget(new QLabel("Плечо тяжести (м):"), row, 0);
    layout->addWidget(rmgEdit, row++, 1);
    layout->addWidget(new QLabel("Плечо подъемной силы (м):"), row, 0);
    layout->addWidget(rlEdit, row++, 1);
    
    // Демпфирование
    layout->addWidget(new QLabel("Kx (гориз.):"), row, 0);
    layout->addWidget(kxEdit, row++, 1);
    layout->addWidget(new QLabel("Ky (вертик.):"), row, 0);
    layout->addWidget(kyEdit, row++, 1);
    layout->addWidget(new QLabel("Kθ (угловое):"), row, 0);
    layout->addWidget(kthetaEdit, row++, 1);
    
    // Равновесные параметры
    layout->addWidget(new QLabel("Равновесная глубина (м):"), row, 0);
    layout->addWidget(yeqEdit, row++, 1);
    layout->addWidget(new QLabel("Равновесная скорость (м/с):"), row, 0);
    layout->addWidget(u0Edit, row++, 1);
    layout->addWidget(new QLabel("Начальный объем (м³):"), row, 0);
    layout->addWidget(v0Edit, row++, 1);
    
    // Параметры среды
    layout->addWidget(new QLabel("Плотность воды (кг/м³):"), row, 0);
    layout->addWidget(rhoEdit, row++, 1);
    layout->addWidget(new QLabel("Атмосферное давление (Па):"), row, 0);
    layout->addWidget(p0Edit, row++, 1);
    
    group->setLayout(layout);
    return group;
}

void MainWindow::loadParametersFromUI()
{
    // Основные параметры
    currentParams.mass = massEdit->text().toDouble();
    currentParams.inertia = inertiaEdit->text().toDouble();
    currentParams.lambda = lambdaEdit->text().toDouble();
    currentParams.thrust = thrustEdit->text().toDouble();
    currentParams.fY = fYEdit->text().toDouble();
    currentParams.g = gEdit->text().toDouble();
    
    // Плечи сил
    currentParams.ra = raEdit->text().toDouble();
    currentParams.rmg = rmgEdit->text().toDouble();
    currentParams.rl = rlEdit->text().toDouble();
    
    // Демпфирование
    currentParams.kx = kxEdit->text().toDouble();
    currentParams.ky = kyEdit->text().toDouble();
    currentParams.ktheta = kthetaEdit->text().toDouble();
    
    // Равновесные параметры
    currentParams.yeq = yeqEdit->text().toDouble();
    currentParams.u0 = u0Edit->text().toDouble();
    currentParams.v0 = v0Edit->text().toDouble();
    
    // Параметры среды
    currentParams.rho = rhoEdit->text().toDouble();
    currentParams.p0 = p0Edit->text().toDouble();
    
    // Начальные условия
    currentParams.initX = xEdit->text().toDouble();
    currentParams.initY = yEdit->text().toDouble();
    currentParams.initVx = vxEdit->text().toDouble();
    currentParams.initVy = vyEdit->text().toDouble();
    currentParams.initTheta = thetaEdit->text().toDouble();
    currentParams.initOmega = omegaEdit->text().toDouble();
    
    // Параметры симуляции
    currentParams.duration = durationEdit->text().toDouble();
    currentParams.dt = dtEdit->text().toDouble();
}

void MainWindow::saveParametersToUI()
{
    // Основные параметры
    massEdit->setText(QString::number(currentParams.mass));
    inertiaEdit->setText(QString::number(currentParams.inertia));
    lambdaEdit->setText(QString::number(currentParams.lambda));
    thrustEdit->setText(QString::number(currentParams.thrust));
    fYEdit->setText(QString::number(currentParams.fY));
    gEdit->setText(QString::number(currentParams.g));
    
    // Плечи сил
    raEdit->setText(QString::number(currentParams.ra));
    rmgEdit->setText(QString::number(currentParams.rmg));
    rlEdit->setText(QString::number(currentParams.rl));
    
    // Демпфирование
    kxEdit->setText(QString::number(currentParams.kx));
    kyEdit->setText(QString::number(currentParams.ky));
    kthetaEdit->setText(QString::number(currentParams.ktheta));
    
    // Равновесные параметры
    yeqEdit->setText(QString::number(currentParams.yeq));
    u0Edit->setText(QString::number(currentParams.u0));
    v0Edit->setText(QString::number(currentParams.v0));
    
    // Параметры среды
    rhoEdit->setText(QString::number(currentParams.rho));
    p0Edit->setText(QString::number(currentParams.p0));
    
    // Начальные условия
    xEdit->setText(QString::number(currentParams.initX));
    yEdit->setText(QString::number(currentParams.initY));
    vxEdit->setText(QString::number(currentParams.initVx));
    vyEdit->setText(QString::number(currentParams.initVy));
    thetaEdit->setText(QString::number(currentParams.initTheta));
    omegaEdit->setText(QString::number(currentParams.initOmega));
    
    // Параметры симуляции
    durationEdit->setText(QString::number(currentParams.duration));
    dtEdit->setText(QString::number(currentParams.dt));
}

void MainWindow::resetParametersToDefault()
{
    currentParams = SimulationParameters(); // Сбрасываем к значениям по умолчанию
    saveParametersToUI();
    QMessageBox::information(this, "Сброс", "Параметры сброшены к значениям по умолчанию");
}

void MainWindow::fillResultsTable(const QVector<double>& times,
                                 const QVector<double>& xPositions,
                                 const QVector<double>& yPositions,
                                 const QVector<double>& angles,
                                 const QVector<double>& vxVelocities,
                                 const QVector<double>& vyVelocities)
{
    int rowCount = std::min(1000, times.size());
    resultsTable->setRowCount(rowCount);
    
    for (int i = 0; i < rowCount; ++i) {
        resultsTable->setItem(i, 0, new QTableWidgetItem(QString::number(times[i], 'f', 2)));
        resultsTable->setItem(i, 1, new QTableWidgetItem(QString::number(xPositions[i], 'f', 4)));
        resultsTable->setItem(i, 2, new QTableWidgetItem(QString::number(yPositions[i], 'f', 4)));
        resultsTable->setItem(i, 3, new QTableWidgetItem(QString::number(angles[i], 'f', 4)));
        resultsTable->setItem(i, 4, new QTableWidgetItem(QString::number(vxVelocities[i], 'f', 4)));
        resultsTable->setItem(i, 5, new QTableWidgetItem(QString::number(vyVelocities[i], 'f', 4)));
    }
    
    resultsTable->resizeColumnsToContents();
}

void MainWindow::loadDataAndUpdateTable(const QString& filename)
{
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::warning(this, "Ошибка", "Не удалось открыть файл: " + filename);
        return;
    }

    QTextStream in(&file);
    QVector<double> times, xPositions, yPositions, angles, vxVelocities, vyVelocities;
    
    // Пропускаем заголовок
    QString header = in.readLine();
    
    // Читаем данные
    while (!in.atEnd()) {
        QString line = in.readLine();
        QStringList fields = line.split(" ", Qt::SkipEmptyParts);
        
        if (fields.size() >= 6) {
            times.append(fields[0].toDouble());
            xPositions.append(fields[1].toDouble());
            yPositions.append(fields[2].toDouble());
            angles.append(fields[3].toDouble());
            vxVelocities.append(fields[4].toDouble());
            vyVelocities.append(fields[5].toDouble());
        }
    }
    file.close();

    // Заполняем таблицу
    fillResultsTable(times, xPositions, yPositions, angles, vxVelocities, vyVelocities);
    
    // Обновляем статус
    statusBar()->showMessage(QString("Загружено %1 строк данных").arg(times.size()));
}
