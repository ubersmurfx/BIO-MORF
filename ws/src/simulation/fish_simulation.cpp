#include "fish_simulation.hpp"

FishTailSimulation::State::State(double x, double y, double vx, double vy, 
                                double theta, double omega, double t) 
    : x(x), y(y), vx(vx), vy(vy), theta(theta), omega(omega), time(t) {}


FishTailSimulation::FishTailSimulation() {
    V0 = m / rho;
    
    // Находим глубину статического равновесия
    y_neutral = find_static_equilibrium_depth();
    
    std::cout << "=== ПАРАМЕТРЫ ДЛЯ ДИНАМИЧЕСКОЙ СТАБИЛИЗАЦИИ ===" << std::endl;
    std::cout << "Глубина статического равновесия: " << y_neutral << " м" << std::endl;
    std::cout << "Плечо архимедовой силы (ra): " << ra << " м" << std::endl;
    std::cout << "Плечо силы тяжести (rmg): " << rmg << " м" << std::endl;
    std::cout << "Плечо подъемной силы (rl): " << rl << " м" << std::endl;
    std::cout << "Коэффициент подъемной силы (Lambda): " << Lambda << std::endl;
    std::cout << "=============================================" << std::endl;
}

double FishTailSimulation::archimedes_force(double y) const {
    // При y=0 должна быть нейтральная плавучесть: Fa = m*g
    // Учитываем изменение объема с глубиной по закону Бойля-Мариотта
    double volume_at_depth = V0 * (p0 / (p0 + rho * g * y)); // учет изменения объема
    return rho * g * volume_at_depth; 
}

bool FishTailSimulation::is_dynamically_stable(const FishTailSimulation::State& state, double tolerance) const {
    return (std::abs(state.vy) < tolerance && 
            std::abs(state.omega) < tolerance &&
            std::abs(state.theta) < 0.1);
}

double FishTailSimulation::thrust_force(double y, double theta) const {
    double delta_y = y - yeq;
    double result_thrust_force = Ft0 + Fy * delta_y + G * theta;
    // std::cout << result_thrust_force << std::endl;
    return result_thrust_force;
}

std::vector<double> FishTailSimulation::equations_of_motion(double t, const std::vector<double>& state) const {
    double x = state[0], y = state[1];
    double vx = state[2], vy = state[3];
    double theta = state[4], omega = state[5];
    
    double Ft = thrust_force(y, theta);
    double Fa = archimedes_force(y);
    double Fl = lift_force(vx);

    // ОГРАНИЧИВАЕМ СИЛЫ
    Ft = std::max(-100.0, std::min(100.0, Ft));
    Fa = std::max(0.0, std::min(1000.0, Fa));
    Fl = std::max(0.0, std::min(1000.0, Fl));

    double ay = (Fa - m*g + Fl - Ft * std::sin(theta) - ky * vy * std::abs(vy)) / m;
    
    // Уравнение для X (без изменений)
    double ax = (Ft * std::cos(theta)) / m - kx * vx * std::abs(vx);
    
    std::cout << "t=" << t << " | Fa=" << Fa << " | mg=" << m*g << " | Fl=" << Fl 
                  << " | Ft=" << Ft << " | ay=" << ay << std::endl;

    
    // Вращательное движение
    double moment_archimedes = ra * Fa * std::cos(theta);
    double moment_gravity = rmg * m * g * std::cos(theta);
    double moment_lift = rl * Fl;
    
    double alpha = (moment_archimedes + moment_gravity + moment_lift - 
                   ktheta * omega * std::abs(omega)) / I;
    
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
        state[4] = std::max(-M_PI/4, std::min(M_PI/4, state[4]));
        state[5] = std::max(-5.0, std::min(5.0, state[5]));
        
        results.emplace_back(state[0], state[1], state[2], state[3], state[4], state[5], t);
    }
    return results;
}


double FishTailSimulation::lift_force(double vx) const {
    return Lambda * vx * vx;
}

double FishTailSimulation::find_static_equilibrium_depth() const {
    // Решаем уравнение: archimedes_force(y) = m*g
    // Используем метод бисекции для нахождения корня
    double y_left = -10.0;  // выше поверхности
    double y_right = 10.0;  // глубже
    
    double fa = archimedes_force(y_left) - m*g;
    double fb = archimedes_force(y_right) - m*g;
    
    if (fa * fb >= 0) {
        std::cout << "Предупреждение: не удалось найти статическое равновесие в диапазоне" << std::endl;
        return 0.0;
    }
    
    double y_mid;
    for (int i = 0; i < 50; ++i) {
        y_mid = (y_left + y_right) / 2.0;
        double fm = archimedes_force(y_mid) - m*g;
        
        if (std::abs(fm) < 1e-6) {
            break;
        }
        
        if (fa * fm < 0) {
            y_right = y_mid;
            fb = fm;
        } else {
            y_left = y_mid;
            fa = fm;
        }
    }
    
    return y_mid;
}

double FishTailSimulation::calculate_dynamic_equilibrium(double vx) const {
    // Для динамического равновесия: m*g = Fa(y) + Fl(vx)
    // Решаем уравнение: archimedes_force(y) = m*g - Lambda*vx*vx
    
    double target_force = m*g - lift_force(vx);
    
    // Используем метод бисекции
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

void FishTailSimulation::save_data_to_file(const std::vector<State>& results, 
                      const std::string& filename) {
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
