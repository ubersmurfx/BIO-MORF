#pragma once
#include "simulation_parameters.hpp"
#include <vector>
#include <string>
#include <stdexcept>

class FishTailSimulation {
public:
    struct State {
        double x = 0;           // Horizontal position (m)
        double z = 0;           // Vertical position (m) - positive up (как в статье)
        double zdot = 0;        // Vertical velocity (m/s)
        double theta = 0;       // Pitch angle (rad) - positive nose up (как в статье)
        double thetadot = 0;    // Angular velocity (rad/s)
        double time = 0;        // Time (s)
        
        // Для совместимости с существующим кодом
        double y = 0;           // y = -z (для привычной системы координат)
        double vx = 0;          // Horizontal velocity
        double vy = 0;          // Vertical velocity (vy = -zdot)
        double omega = 0;       // Angular velocity (omega = thetadot)
        
        State() = default;
        State(double x, double z, double zdot, double theta, double thetadot, double t = 0)
            : x(x), z(z), zdot(zdot), theta(theta), thetadot(thetadot), time(t) {
            y = -z;
            vy = -zdot;
            omega = thetadot;
        }
    };
    
    FishTailSimulation();
    FishTailSimulation(const SimulationParameters& params);
    
    void setParameters(const SimulationParameters& new_params);
    void print_parameters() const;
    
    // Основные функции симуляции
    std::vector<State> simulate(double duration, double dt, const State& initial_state);
    bool is_system_stable() const;
    
    // Вспомогательные функции
    void save_data_to_file(const std::vector<State>& results, const std::string& filename = "trajectory_data.txt");
    void create_gnuplot_script() const;
    void run_gnuplot();
    void print_simulation_results(const std::vector<State>& results);
    
    // Анализ устойчивости
    void analyze_stability_window(double U_min = 0.1, double U_max = 5.0, int steps = 20);


    struct EquilibriumAnalysis {
        double y_static;      // Глубина статического равновесия
        double y_dynamic;     // Глубина динамического равновесия при текущей U
        bool can_stabilize_from_above; // Можно стабилизироваться при старте сверху
        bool can_stabilize_from_below; // Можно стабилизироваться при старте снизу
        double min_depth_for_stabilization; // Минимальная глубина для достижения динамического равновесия
    };
    double find_dynamic_equilibrium_depth() const;

private:
    SimulationParameters params;
    
    // Основные уравнения движения из статьи (уравнение 4)
    std::vector<double> equations_of_motion(double t, const std::vector<double>& state) const;
    std::vector<double> runge_kutta_4(double t, const std::vector<double>& state, double dt) const;
    
    // Расчет коэффициента жесткости плавучести k_b из статьи
    double calculate_buoyancy_stiffness() const;
    
    double calculate_U1() const;
    double calculate_U2() const;
    void calculate_optimal_geometry() const;
    double calculate_B_rL() const;
    bool check_delta3_condition(double U, double k_b, double K) const;
    void analyze_stability_map() const;
    double calculate_neutral_buoyancy_depth() const;



};