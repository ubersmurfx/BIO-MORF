#ifndef FISH_SIMULATION_HPP
#define FISH_SIMULATION_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include "simulation_parameters.hpp"


class FishTailSimulation {
public:


    struct State {
        double x, y, vx, vy, theta, omega;
        double time;

        State(double x=0, double y=0, double vx=0, double vy=0, 
              double theta=0, double omega=0, double t=0);
    };

    // Конструкторы
    FishTailSimulation();
    FishTailSimulation(const SimulationParameters& params);
    
    // Методы симуляции
    bool is_dynamically_stable(const State& state, double tolerance = 0.01) const;
    std::vector<double> equations_of_motion(double t, const std::vector<double>& state) const;
    std::vector<double> runge_kutta_4(double t, const std::vector<double>& state, double dt) const;
    std::vector<State> simulate(double duration, double dt, const State& initial_state);

    // Методы для работы с параметрами
    void setParameters(const SimulationParameters& new_params);
    
    // Методы ввода-вывода
    void save_data_to_file(const std::vector<State>& results, const std::string& filename);
    void create_gnuplot_script();
    void create_3d_gnuplot_script();
    void run_gnuplot();
    void print_simulation_results(const std::vector<State>& results);
    
    // Вспомогательные методы
    void print_parameters() const;
    void check_initial_moment();

private:
    SimulationParameters params;
    mutable double y_curr;
    mutable double adaptive_v0;
    double y_neutral;

    // Приватные методы
    double archimedes_force(double y) const;
    double adaptive_archimedes_force(double y, double vy) const;
    double thrust_force(double y, double theta) const;
    double lift_force(double vx) const;
    double find_static_equilibrium_depth() const;
    double calculate_dynamic_equilibrium(double vx) const;
    double calculate_critical_depth() const;
};

#endif // FISH_SIMULATION_HPP