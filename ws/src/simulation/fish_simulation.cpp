#include "simulation/fish_simulation.hpp"


FishTailSimulation::State::State(double x, double y, double vx, double vy, 
                                double theta, double omega, double t) 
    : x(x), y(y), vx(vx), vy(vy), theta(theta), omega(omega), time(t) {}


FishTailSimulation::FishTailSimulation() : params(), adaptive_v0(params.v0) {
    params.validate();
    y_neutral = find_static_equilibrium_depth();
    check_initial_moment();
    print_parameters();
}

// Конструктор с параметрами
FishTailSimulation::FishTailSimulation(const SimulationParameters& params) 
    : params(params), adaptive_v0(params.v0) {
    params.validate();
    y_neutral = find_static_equilibrium_depth();
    check_initial_moment();
    print_parameters();
}

void FishTailSimulation::print_parameters() const {
    std::cout << "=== Параметры симуляции АНПА ===" << std::endl;
    std::cout << "Масса: " << params.mass << " кг" << std::endl;
    std::cout << "Момент инерции: " << params.inertia << " кг·м²" << std::endl;
    std::cout << "Тяга: " << params.thrust << " Н" << std::endl;
    std::cout << "Коэффициент Fy: " << params.fY << std::endl;
    std::cout << "Коэффициент G: " << params.g << std::endl;
    std::cout << "Плечо Архимеда: " << params.ra << " м" << std::endl;
    std::cout << "Плечо тяжести: " << params.rmg << " м" << std::endl;
    std::cout << "Плечо подъемной силы: " << params.rl << " м" << std::endl;
    std::cout << "Демпфирование: Kx=" << params.kx << ", Ky=" << params.ky << ", Kθ=" << params.ktheta << std::endl;
    std::cout << "Равновесная глубина: " << params.yeq << " м" << std::endl;
    std::cout << "Начальная скорость: " << params.u0 << " м/с" << std::endl;
    std::cout << "Начальный объем: " << params.v0 << " м³" << std::endl;
    std::cout << "Плотность воды: " << params.rho << " кг/м³" << std::endl;
    std::cout << "Давление: " << params.p0 << " Па" << std::endl;
    std::cout << "Длительность: " << params.duration << " с" << std::endl;
    std::cout << "Шаг времени: " << params.dt << " с" << std::endl;
    std::cout << "===============================" << std::endl;
}

void FishTailSimulation::setParameters(const SimulationParameters& new_params) {
    new_params.validate();
    params = new_params;
    adaptive_v0 = params.v0;
    y_neutral = find_static_equilibrium_depth();
}

double FishTailSimulation::archimedes_force(double y) const {
    double volume_at_depth = params.v0 * (params.p0 / (params.p0 + params.rho * params.g * y));
    return params.rho * params.g * volume_at_depth; 
}

bool FishTailSimulation::is_dynamically_stable(const State& state, double tolerance) const {
    return (std::abs(state.vy) < tolerance && 
            std::abs(state.omega) < tolerance &&
            std::abs(state.theta) < 0.1);
}

double FishTailSimulation::thrust_force(double y, double theta) const {
    double constrained_y = std::max(0.0, y);
    double delta_y = constrained_y - params.yeq;
    double result_thrust_force = params.thrust - params.fY * delta_y + params.g * theta;
    return result_thrust_force;
}

double FishTailSimulation::adaptive_archimedes_force(double y, double vy) const {
    double correction = 0.0;
    
    if (y < params.yeq && vy < -0.01) {
        correction = -0.00001;
    } else if (y > params.yeq && vy > 0.01) {
        correction = 0.00001;
    }
    
    adaptive_v0 += correction;
    adaptive_v0 = std::max(0.0025, std::min(0.0035, adaptive_v0));
    
    // std::cout << "ADAPTIVE_DEBUG: v0=" << adaptive_v0 << " correction=" << correction 
    //           << " y=" << y << " vy=" << vy << std::endl;
    
    return params.rho * params.g * adaptive_v0 * (params.p0 / (params.p0 + params.rho * params.g * y));
}

std::vector<double> FishTailSimulation::equations_of_motion(double t, const std::vector<double>& state) const {
    double x = state[0], y = state[1];
    double vx = state[2], vy = state[3];
    double theta = state[4], omega = state[5];
    
    double constrained_y = std::max(0.0, y);
    
    // вычисление сил
    double Ft = thrust_force(constrained_y, theta);
    double Fa = adaptive_archimedes_force(constrained_y, vy);
    double Fl = lift_force(vx);

    // ограничение сил
    Ft = std::max(-100.0, std::min(100.0, Ft));
    Fa = std::max(0.0, std::min(1000.0, Fa));
    Fl = std::max(0.0, std::min(1000.0, Fl));

    // Уравнения движения
    double ax = ((Ft * std::cos(theta)) - params.kx * vx * std::abs(vx)) / params.mass;
    double ay = (-Fa + params.mass * params.g - Fl - Ft * std::sin(theta) - 
                params.ky * vy * std::abs(vy)) / params.mass;

    // Вращательное движение
    double moment_archimedes = params.ra * Fa * std::cos(theta);
    double moment_gravity = params.rmg * params.mass * params.g * std::cos(theta);
    double moment_lift = params.rl * Fl;
    
    double alpha = (moment_archimedes + moment_gravity + moment_lift - 
                   params.ktheta * omega * std::abs(omega)) / params.inertia;
    
    y_curr = y;

    return {vx, vy, ax, ay, omega, alpha};
}

std::vector<double> FishTailSimulation::runge_kutta_4(double t, const std::vector<double>& state, double dt) const {
    auto k1 = equations_of_motion(t, state);
    
    std::vector<double> state2(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        state2[i] = state[i] + 0.5 * dt * k1[i];
    }
    auto k2 = equations_of_motion(t + 0.5 * dt, state2);
    
    std::vector<double> state3(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        state3[i] = state[i] + 0.5 * dt * k2[i];
    }
    auto k3 = equations_of_motion(t + 0.5 * dt, state3);
    
    std::vector<double> state4(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        state4[i] = state[i] + dt * k3[i];
    }
    auto k4 = equations_of_motion(t + dt, state4);
    
    std::vector<double> new_state(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        new_state[i] = state[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
    
    return new_state;
}

std::vector<FishTailSimulation::State> FishTailSimulation::simulate(double duration, double dt, const State& initial_state) {
    std::vector<State> results;
    
    std::vector<double> state = {
        initial_state.x, initial_state.y,
        initial_state.vx, initial_state.vy,
        initial_state.theta, initial_state.omega
    };
    
    int steps = static_cast<int>(duration / dt);
    results.reserve(steps + 1);
    results.push_back(initial_state);
    
    for (int i = 0; i < steps; ++i) {
        double t = i * dt;
        
        bool valid = true;
        for (double val : state) {
            if (std::isnan(val) || std::isinf(val)) {
                std::cout << "Numerical instability at step " << i << std::endl;
                valid = false;
                break;
            }
        }
        if (!valid) break;
        
        state = runge_kutta_4(t, state, dt);
        
        // Ограничения для стабильности
        state[1] = std::max(0.0, state[1]);
        results.emplace_back(state[0], state[1], state[2], state[3], state[4], state[5], t);
    }
    return results;
}

double FishTailSimulation::lift_force(double vx) const {
    return params.lambda * vx * vx;
}

// найти глубину нейтральной плавучести
double FishTailSimulation::find_static_equilibrium_depth() const {
    double y_left = 0.0;
    double y_right = 10.0;
    
    std::cout << "Поиск горизонта..." << std::endl;
    
    for (int i = 0; i < 50; ++i) {
        double y_mid = (y_left + y_right) / 2.0;
        double fa = archimedes_force(y_mid);
        double error = fa - params.mass * params.g;
        
        if (std::abs(error) < 1e-8) {
            return y_mid;
        }
        
        if (error > 0) {
            y_left = y_mid;
        } else {
            y_right = y_mid;
        }
        
        if (y_right - y_left < 1e-12) {
            std::cout << "Достигнута точность: " << y_mid << " м" << std::endl;
            return y_mid;
        }
    }
    
    double result = (y_left + y_right) / 2.0;
    std::cout << "Финальный результат: " << result << " м" << std::endl;
    return result;
}

double FishTailSimulation::calculate_dynamic_equilibrium(double vx) const {
    double target_force = params.mass * params.g - lift_force(vx);
    
    double y_left = -10.0;
    double y_right = 10.0;
    
    for (int i = 0; i < 50; ++i) {
        double y_mid = (y_left + y_right) / 2.0;
        double fa = archimedes_force(y_mid);
        
        if (std::abs(fa - target_force) < 1e-6) {
            return y_mid;
        }
        
        if (fa < target_force) {
            y_left = y_mid;
        } else {
            y_right = y_mid;
        }
    }
    
    return (y_left + y_right) / 2.0;
}

double FishTailSimulation::calculate_critical_depth() const {
    double v_max = 5.0;
    double FL_max = params.lambda * v_max * v_max;
    
    double y_crit = (params.mass * params.g * params.p0 / 
                    (params.rho * params.g * params.v0 * (params.mass * params.g - FL_max)) - 
                    params.p0) / (params.rho * params.g);
    
    std::cout << "Теоретическая критическая глубина: " << y_crit << " м" << std::endl;
    return y_crit;
}

void FishTailSimulation::save_data_to_file(const std::vector<State>& results, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    file << "Time X Y Theta Vx Vy\n";
    for (const auto& state : results) {
        file << std::fixed << std::setprecision(6)
             << state.time << " " << state.x << " " << state.y << " "
             << state.theta << " " << state.vx << " " << state.vy << "\n";
    }
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

void FishTailSimulation::create_gnuplot_script() {
  std::ofstream script("plot_trajectory.gnu");
  script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n";
  script << "set output 'trajectory_plots.png'\n\n";
  
  script << "set multiplot layout 2,2 title 'Fish Motion Analysis'\n\n";
  
  // График 1: Траектория X-Y
  script << "set title 'Trajectory X-Y'\n";
  script << "set xlabel 'X position'\n";
  script << "set ylabel 'Y position'\n";
  script << "plot 'trajectory_data.txt' using 2:3 with lines lw 2 title 'Trajectory', \\\n";
  script << "     'trajectory_data.txt' every 50 using 2:3 with points pt 7 ps 0.5 title 'Points'\n\n";
  
  // График 2: Положение от времени
  script << "set title 'Position vs Time'\n";
  script << "set xlabel 'Time'\n";
  script << "set ylabel 'Position'\n";
  script << "plot 'trajectory_data.txt' using 1:2 with lines lw 2 title 'X position', \\\n";
  script << "     'trajectory_data.txt' using 1:3 with lines lw 2 title 'Y position'\n\n";
  
  // График 3: Угол от времени
  script << "set title 'Angle vs Time'\n";
  script << "set xlabel 'Time'\n";
  script << "set ylabel 'Angle (rad)'\n";
  script << "plot 'trajectory_data.txt' using 1:4 with lines lw 2 title 'Theta'\n\n";
  
  // График 4: Скорости от времени
  script << "set title 'Velocity vs Time'\n";
  script << "set xlabel 'Time'\n";
  script << "set ylabel 'Velocity'\n";
  script << "plot 'trajectory_data.txt' using 1:5 with lines lw 2 title 'Vx', \\\n";
  script << "     'trajectory_data.txt' using 1:6 with lines lw 2 title 'Vy'\n\n";
  
  script << "unset multiplot\n";
  script << "set output\n";
  script.close();
  
  std::cout << "GNUplot script created: plot_trajectory.gnu" << std::endl;
}

void FishTailSimulation::create_3d_gnuplot_script() {
    std::ofstream script("plot_3d.gnu");
    script << "set terminal pngcairo size 800,600 enhanced font 'Arial,12'\n";
    script << "set output '3d_trajectory.png'\n\n";
    
    script << "set title '3D Trajectory: X-Y-Time'\n";
    script << "set xlabel 'X position'\n";
    script << "set ylabel 'Y position'\n";
    script << "set zlabel 'Time'\n";
    script << "splot 'trajectory_data.txt' using 2:3:1 with lines lw 2 title 'Trajectory'\n";
    script.close();
    
    std::cout << "3D GNUplot script created: plot_3d.gnu" << std::endl;
}

void FishTailSimulation::run_gnuplot() {
  std::cout << "Generating plots with GNUplot..." << std::endl;
  
  // Создаем скрипты
  create_gnuplot_script();
  create_3d_gnuplot_script();
  
  // Запускаем GNUplot
  int result1 = std::system("gnuplot plot_trajectory.gnu");
  int result2 = std::system("gnuplot plot_3d.gnu");
  
  if (result1 == 0 && result2 == 0) {
      std::cout << "Plots generated successfully!" << std::endl;
      std::cout << "Check files: trajectory_plots.png and 3d_trajectory.png" << std::endl;
  } else {
      std::cout << "Error running GNUplot. Make sure GNUplot is installed." << std::endl;
      std::cout << "You can install it with: sudo apt-get install gnuplot" << std::endl;
  }
}

void FishTailSimulation::print_simulation_results(const std::vector<State>& results) {
    std::cout << "Time\tX\tY\tTheta\tVx\tVy\n";
    std::cout << "--------------------------------------------------------\n";
    
    for (size_t i = 0; i < std::min(results.size(), size_t(20)); i++) {
        const auto& state = results[i];
        std::cout << std::fixed << std::setprecision(4);
        std::cout << state.time << "\t" << state.x << "\t" << state.y << "\t" 
                  << state.theta << "\t" << state.vx << "\t" << state.vy << "\n";
    }
    
    if (results.size() > 20) {
        std::cout << "... (showing first 20 of " << results.size() << " steps)\n";
    }
}


void FishTailSimulation::check_initial_moment() {
    double Fa = archimedes_force(y_curr);
    double M_arch = params.ra * Fa;
    double M_grav = params.rmg * params.mass * params.g;
    double total_M = M_arch + M_grav;
    
    std::cout << "Момент Архимеда: " << M_arch << " Н·м" << std::endl;
    std::cout << "Момент тяжести: " << M_grav << " Н·м" << std::endl; 
    std::cout << "Суммарный момент: " << total_M << " Н·м" << std::endl;
    std::cout << "Ожидаемое поведение: " 
         << (total_M < 0 ? "НОС ВНИЗ" : "НОС ВВЕРХ") << std::endl;
}