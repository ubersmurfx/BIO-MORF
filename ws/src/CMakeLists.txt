cmake_minimum_required(VERSION 3.5)
project(FishSimulationQt)

# Поиск Qt5
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Charts)

# Проверяем наличие Charts
find_package(Qt5Charts QUIET)
if(Qt5Charts_FOUND)
    set(USE_CHARTS TRUE)
    message(STATUS "Qt5 Charts found")
else()
    set(USE_CHARTS FALSE)
    message(STATUS "Qt5 Charts not found - disabling charts")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Автоматически подключаем moc, uic, rcc
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Исходные файлы
set(SOURCES
    main.cpp
    MainWindow.cpp
    SimulationThread.cpp
    simulation/fish_simulation.cpp
)

# Заголовочные файлы
set(HEADERS
    MainWindow.hpp
    SimulationThread.hpp
    simulation/fish_simulation.hpp
)

# Создаем исполняемый файл
add_executable(FishSimulationQt ${SOURCES} ${HEADERS})

# Подключаем библиотеки Qt5
target_link_libraries(FishSimulationQt Qt5::Core Qt5::Widgets)

# Если Charts доступны, подключаем их
if(USE_CHARTS)
    target_link_libraries(FishSimulationQt Qt5::Charts)
    target_compile_definitions(FishSimulationQt PRIVATE USE_QT_CHARTS)
endif()

# Настройки компиляции
target_compile_options(FishSimulationQt PRIVATE -Wall -Wextra -O2)