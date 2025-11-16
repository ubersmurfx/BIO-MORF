#include "simulation/fish_simulation.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

// Конструктор по умолчанию
FishTailSimulation::FishTailSimulation() : params() {
    params.validate();
}

// Конструктор с параметрами
FishTailSimulation::FishTailSimulation(const SimulationParameters& params) 
    : params(params) {
    params.validate();
    is_system_stable();
}

void FishTailSimulation::print_parameters() const {
    double k_b = calculate_buoyancy_stiffness();
    double U = params.u0;  // Используем u0 как скорость U из статьи
    
    std::cout << "=== Параметры симуляции (Модель из статьи) ===" << std::endl;
    std::cout << "Масса (m): " << params.mass << " кг" << std::endl;
    std::cout << "Момент инерции (J): " << params.inertia << " кг·м²" << std::endl;
    std::cout << "Скорость U: " << U << " м/с" << std::endl;
    std::cout << "Коэффициент подъемной силы K (lambda): " << params.lambda << " кг/м" << std::endl;
    std::cout << "Жесткость плавучести k_b: " << k_b << " Н/м" << std::endl;
    std::cout << "Структурная жесткость k_θ: " << params.ktheta << " Н·м/рад" << std::endl;
    std::cout << "Смещение плавучести rA: " << params.ra << " м" << std::endl;
    std::cout << "Смещение подъемной силы rL: " << params.rl << " м" << std::endl;
    std::cout << "Объем пузыря V0: " << params.v0 << " м³" << std::endl;
    std::cout << "Давление p0: " << params.p0 << " Па" << std::endl;
    std::cout << "===============================" << std::endl;
    std::cout << "Нейтральная глубина для заданных параметров" << calculate_neutral_buoyancy_depth() << std::endl;
}

void FishTailSimulation::setParameters(const SimulationParameters& new_params) {
    new_params.validate();
    params = new_params;
}

// Расчет жесткости плавучести k_b по формулe
double FishTailSimulation::calculate_buoyancy_stiffness() const {
    // k_b = (ρ_w² g² V_0) / (κ p_0)
    // Для упрощения принимаем κ = 1.4 (адиабатический процесс)
    double kappa = 1.4;
    double k_b = (params.rho * params.rho * params.g * params.g * params.v0) / (kappa * params.p0);
    return k_b;
}

// Основные уравнения движения из статьи (уравнение 4)
std::vector<double> FishTailSimulation::equations_of_motion(double t, const std::vector<double>& state) const {
    double z = state[0];
    double zdot = state[1];
    double theta = state[2];
    double thetadot = state[3];
    
    // Автоматически вычисляем статическое равновесие
    double y_static = (params.rho * params.v0 * params.p0 / params.mass - params.p0) / (params.rho * params.g);
    double z_equilibrium = -y_static; // z = -y
    double z_relative = z - z_equilibrium;
    
    // Остальной код без изменений...
    double U = params.u0;
    double k_b = calculate_buoyancy_stiffness();
    double K = params.lambda;
    
    double z_double_dot = (-K * U * zdot + 
                          k_b * z_relative + 
                          (k_b * params.ra + K * U * U) * theta) / params.mass;
    
    double net_pitch_stiffness = params.ktheta - params.ra * params.ra * k_b - params.rl * K * U * U;
    double theta_double_dot = (-params.rl * K * U * zdot + 
                              params.ra * k_b * z_relative + 
                              - net_pitch_stiffness * theta) / params.inertia;
    
    return {zdot, z_double_dot, thetadot, theta_double_dot};
}

// Метод Рунге-Кутты 4-го порядка
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

// Комплексный анализ устойчивости с учетом формул (5) из статьи
bool FishTailSimulation::is_system_stable() const {
    double U = params.u0;
    double k_b = calculate_buoyancy_stiffness();
    double K = params.lambda;
    
    // 1. Условие a₄ > 0 (U > U₁)
    double U1 = calculate_U1();
    bool condition_a4 = (params.rl > params.ra) ? (U > U1) : false;
    
    // 2. Условие Δ₃ > 0 ((KU² + kb·rA)·B(rL) < 0)
    bool condition_delta3 = check_delta3_condition(U, k_b, K);
    
    // 3. Геометрическое условие из Fig. 3
    bool geometric_condition = (params.ra < params.rl) && (params.rl < 0);
    
    bool fully_stable = condition_a4 && condition_delta3 && geometric_condition;
    
    std::cout << "ПОЛНЫЙ АНАЛИЗ УСТОЙЧИВОСТИ (Fig. 3):" << std::endl;
    std::cout << "  U = " << U << " м/с, rA = " << params.ra << ", rL = " << params.rl << std::endl;
    std::cout << "  a₄ > 0 (U > U₁): " << U << " > " << U1 << " = " << (condition_a4 ? "ДА" : "НЕТ") << std::endl;
    std::cout << "  Δ₃ > 0: " << (condition_delta3 ? "ДА" : "НЕТ") << std::endl;
    std::cout << "  Геометрия rA < rL < 0: " << (geometric_condition ? "ДА" : "НЕТ") << std::endl;
    std::cout << "  СТАТУС: " << (fully_stable ? "УСТОЙЧИВО" : "НЕУСТОЙЧИВО") << std::endl;
    
    return fully_stable;
}

// Проверка условия Δ₃ > 0: (KU² + kb·rA)·B(rL) < 0
bool FishTailSimulation::check_delta3_condition(double U, double k_b, double K) const {
    double term1 = K * U * U + k_b * params.ra;
    double B_rL = calculate_B_rL(); // Геометрический фактор B(rL)
    
    bool condition = (term1 * B_rL) < 0;
    
    std::cout << "  Δ₃ анализ: (K·U² + kb·rA) = " << term1;
    std::cout << ", B(rL) = " << B_rL;
    std::cout << ", произведение = " << (term1 * B_rL) << std::endl;
    
    return condition;
}

double FishTailSimulation::calculate_B_rL() const {
    // B(rL) - геометрический фактор, обращается в 0 на границе устойчивости
    // Упрощенная модель: B(rL) ≈ rL - r_critical
    double r_critical = -1.5; // Примерное значение из Fig. 3
    return params.rl - r_critical;
}

// Расчет оптимальной геометрии по формуле (6)
void FishTailSimulation::calculate_optimal_geometry() const {
    double k_b = calculate_buoyancy_stiffness();
    double kappa = (k_b * params.inertia) / (params.mass * params.ktheta);
    
    // Оптимальное соотношение плеч
    double alpha_optimal = std::sqrt(kappa / (1 + kappa));
    
    std::cout << "\nОПТИМАЛЬНАЯ ГЕОМЕТРИЯ (формула 6):" << std::endl;
    std::cout << "  κ = (kb·J)/(m·kθ) = " << kappa << std::endl;
    std::cout << "  α* = |rL|/|rA| = √[κ/(1+κ)] = " << alpha_optimal << std::endl;
    
    // Рекомендуемые значения для текущих параметров
    double recommended_rA = -2.0; // Пример
    double recommended_rL = recommended_rA * alpha_optimal;
    
    std::cout << "  Рекомендуемые значения: rA = " << recommended_rA;
    std::cout << ", rL = " << recommended_rL << std::endl;
}

// Анализ устойчивости для разных геометрий (как в Fig. 3)
void FishTailSimulation::analyze_stability_map() const {
    std::cout << "\nАНАЛИЗ КАРТЫ УСТОЙЧИВОСТИ (Fig. 3):" << std::endl;
    std::cout << "rA\\rL\t";
    
    // Заголовок
    std::vector<double> rL_values = {-3.0, -2.5, -2.0, -1.5, -1.0};
    for (double rL : rL_values) {
        std::cout << rL << "\t";
    }
    std::cout << std::endl;
    
    // Таблица устойчивости
    std::vector<double> rA_values = {-3.0, -2.5, -2.0, -1.5, -1.0};
    for (double rA : rA_values) {
        std::cout << rA << "\t";
        for (double rL : rL_values) {
            // Временно меняем параметры для анализа
            SimulationParameters temp_params = params;
            temp_params.ra = rA;
            temp_params.rl = rL;
            
            FishTailSimulation temp_sim(temp_params);
            bool stable = temp_sim.is_system_stable();
            
            std::cout << (stable ? "У" : "Н") << "\t";
        }
        std::cout << std::endl;
    }
}

// Улучшенный анализ окна устойчивости
void FishTailSimulation::analyze_stability_window(double U_min, double U_max, int steps) {
    std::cout << "\n=== АНАЛИЗ ОКНА УСТОЙЧИВОСТИ ===" << std::endl;
    std::cout << "Скорость U (м/с)\tЖесткость\tU>U1\tU<U2\tГеом.\tСтатус" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    
    double k_b = calculate_buoyancy_stiffness(); // ИСПРАВЛЕНО: buoyancy
    double K = params.lambda;
    double U1 = calculate_U1();
    double U2 = calculate_U2();
    bool geometric_ok = (params.ra < params.rl) && (params.rl < 0); // ИСПРАВЛЕНО: ra, rl
    
    for (int i = 0; i <= steps; ++i) {
        double U = U_min + (U_max - U_min) * i / steps;
        
        // Вычисляем все критерии
        double net_stiffness = params.ktheta - params.ra * params.ra * k_b - params.rl * K * U * U; // ИСПРАВЛЕНО: ra, rl
        bool condition1 = (params.rl > params.ra) ? (U > U1) : true; // ИСПРАВЛЕНО: rl
        bool condition2 = (params.rl < 0) ? (U < U2) : true; // ИСПРАВЛЕНО: rl
        bool stable = (net_stiffness > 0) && condition1 && condition2 && geometric_ok;
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << U << "\t\t" << std::setprecision(3) << net_stiffness << "\t";
        std::cout << (condition1 ? "ДА" : "НЕТ") << "\t";
        std::cout << (condition2 ? "ДА" : "НЕТ") << "\t";
        std::cout << (geometric_ok ? "ДА" : "НЕТ") << "\t";
        std::cout << (stable ? "УСТОЙЧИВО" : "НЕУСТОЙЧИВО") << std::endl;
    }
    
    // Вывод информации о порогах
    std::cout << "\nПороги устойчивости:" << std::endl;
    std::cout << "  U1 = " << U1 << " м/с" << std::endl;
    std::cout << "  U2 = " << U2 << " м/с" << std::endl;
    std::cout << "  Окно устойчивости: " << U1 << " < U < " << U2 << " м/с" << std::endl;
    
    if (!geometric_ok) {
        std::cout << "\nВНИМАНИЕ: Геометрия не соответствует устойчивой конфигурации!" << std::endl;
        std::cout << "Требуется: rA < rL < 0" << std::endl;
        std::cout << "Фактически: rA = " << params.ra << ", rL = " << params.rl << std::endl;
    }
}

// Основная функция симуляции
std::vector<FishTailSimulation::State> FishTailSimulation::simulate(double duration, double dt, const State& initial_state) {
    std::vector<State> results;
    


    // Состояние: [z, ż, θ, θ̇] - как в статье
    std::vector<double> state = {
        initial_state.z,
        initial_state.zdot, 
        initial_state.theta,
        initial_state.thetadot
    };
    
    int steps = static_cast<int>(duration / dt);
    results.reserve(steps + 1);
    
    // Добавляем начальное состояние
    State current_state = initial_state;
    current_state.time = 0;
    current_state.x = params.initX;
    results.push_back(current_state);
    
    std::cout << "Начало симуляции..." << std::endl;
    std::cout << "Начальные условия: z=" << initial_state.z << " theta=" << initial_state.theta << std::endl;
    std::cout << "Скорость U=" << params.u0 << " м/с" << std::endl;
    
    for (int i = 0; i < steps; ++i) {
        double t = (i + 1) * dt;
        
        // Проверка на численную неустойчивость
        bool valid = true;
        for (double val : state) {
            if (std::isnan(val) || std::isinf(val) || std::abs(val) > 1e6) {
                std::cout << "Численная неустойчивость на шаге " << i << std::endl;
                valid = false;
                break;
            }
        }
        if (!valid) break;
        
        // Интегрирование
        state = runge_kutta_4(t, state, dt);
        
        // Создаем новое состояние
        State new_state;
        new_state.z = state[0];
        new_state.zdot = state[1];
        new_state.theta = state[2];
        new_state.thetadot = state[3];
        new_state.time = t;
        
        // Вычисляем x положение (просто U*t для постоянной скорости)
        new_state.x = params.initX + params.u0 * t;
        new_state.y = -state[0]; // y = -z для привычной системы координат
        new_state.vx = params.u0;
        new_state.vy = -state[1];
        new_state.omega = state[3];
        
        results.push_back(new_state);
        
        // Периодический вывод прогресса
        if (i % static_cast<int>(10.0 / dt) == 0) { // Каждые 10 секунд
            std::cout << "t=" << t << "с: z=" << state[0] << "м, theta=" << state[2] << "рад" << std::endl;
        }
    }
    
    std::cout << "Симуляция завершена. Получено " << results.size() << " шагов." << std::endl;
    return results;
}

// Сохранение данных в файл
void FishTailSimulation::save_data_to_file(const std::vector<State>& results, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
        return;
    }
    
    file << "Time X Z Theta Zdot Thetadot Stability\n";
    for (const auto& state : results) {
        file << std::fixed << std::setprecision(6)
             << state.time << " " << state.x << " " << state.z << " "
             << state.theta << " " << state.zdot << " " << state.thetadot << "\n";
    }
    file.close();
    std::cout << "Данные сохранены в " << filename << std::endl;
}

// Создание скрипта для GNUplot
void FishTailSimulation::create_gnuplot_script() const {
    std::ofstream script("plot_results.gnu");
    script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'\n";
    script << "set output 'stability_analysis.png'\n\n";
    
    script << "set multiplot layout 2,2 title 'Stability Analysis (U=" << params.u0 << " m/s)'\n\n";
    
    // График 1: Глубина и угол от времени
    script << "set title 'Depth and Pitch vs Time'\n";
    script << "set xlabel 'Time (s)'\n";
    script << "set ylabel 'Depth z (m)'\n";
    script << "set y2label 'Pitch θ (rad)'\n";
    script << "set y2tics\n";
    script << "plot 'trajectory_data.txt' using 1:3 with lines lw 2 title 'Depth (z)', \\\n";
    script << "     'trajectory_data.txt' using 1:4 with lines lw 2 axes x1y2 title 'Pitch (θ)'\n\n";
    
    // График 2: Фазовая плоскость глубины
    script << "set title 'Phase Plane: Depth'\n";
    script << "set xlabel 'Depth z (m)'\n";
    script << "set ylabel 'Vertical velocity ż (m/s)'\n";
    script << "unset y2tics\n";
    script << "plot 'trajectory_data.txt' using 3:5 with lines lw 2 title 'Phase trajectory'\n\n";
    
    // График 3: Фазовая плоскость угла
    script << "set title 'Phase Plane: Pitch Angle'\n";
    script << "set xlabel 'Pitch θ (rad)'\n";
    script << "set ylabel 'Angular velocity θ̇ (rad/s)'\n";
    script << "plot 'trajectory_data.txt' using 4:6 with lines lw 2 title 'Phase trajectory'\n\n";
    
    // График 4: Траектория в пространстве
    script << "set title 'Spatial Trajectory'\n";
    script << "set xlabel 'X position (m)'\n";
    script << "set ylabel 'Depth z (m)'\n";
    script << "plot 'trajectory_data.txt' using 2:3 with lines lw 2 title 'Trajectory'\n\n";
    
    script << "unset multiplot\n";
    script << "set output\n";
    script.close();
    
    std::cout << "GNUplot script created: plot_results.gnu" << std::endl;
}

// Запуск GNUplot
void FishTailSimulation::run_gnuplot() {
    std::cout << "Generating plots with GNUplot..." << std::endl;
    
    create_gnuplot_script();
    
    int result = std::system("gnuplot plot_results.gnu");
    
    if (result == 0) {
        std::cout << "Plots generated successfully! Check stability_analysis.png" << std::endl;
    } else {
        std::cout << "Error running GNUplot." << std::endl;
    }
}

// Вывод результатов
void FishTailSimulation::print_simulation_results(const std::vector<State>& results) {
    if (results.empty()) return;
    
    std::cout << "\n=== РЕЗУЛЬТАТЫ СИМУЛЯЦИИ ===" << std::endl;
    std::cout << "Time\tZ\tTheta\tZdot\tThetadot\n";
    std::cout << "----------------------------------------\n";
    
    // Показываем начало, конец и несколько промежуточных точек
    size_t show_points = 8;
    size_t step = std::max(size_t(1), results.size() / show_points);
    
    for (size_t i = 0; i < results.size(); i += step) {
        const auto& state = results[i];
        std::cout << std::fixed << std::setprecision(3);
        std::cout << state.time << "\t" << state.z << "\t" << state.theta << "\t" 
                  << state.zdot << "\t" << state.thetadot << "\n";
        
        if (i >= show_points * step) break;
    }
    
    // Показываем конечное состояние
    const auto& final = results.back();
    std::cout << "...\nФинальное состояние:\n";
    std::cout << "z = " << final.z << " м, theta = " << final.theta << " рад\n";
    std::cout << "Система " << (is_system_stable() ? "УСТОЙЧИВА" : "НЕУСТОЙЧИВА") << " при U=" << params.u0 << " м/с" << std::endl;
}


double FishTailSimulation::calculate_U1() const {
    // U1 = sqrt(kθ / [K(rL - rA)]) при условии rL > rA
    if (params.rl <= params.ra) {
        return std::numeric_limits<double>::infinity();
    }
    double U1_squared = params.ktheta / (params.lambda * (params.rl - params.ra));
    return (U1_squared > 0) ? std::sqrt(U1_squared) : 0.0;
}

double FishTailSimulation::calculate_U2() const {
    // U2 = sqrt( k_b(J + m rL rA) / (-m K rL) ) при условии rL < 0
    if (params.rl >= 0) {
        return std::numeric_limits<double>::infinity();
    }
    double k_b = calculate_buoyancy_stiffness();
    double numerator = k_b * (params.inertia + params.mass * params.rl * params.ra);
    double denominator = -params.mass * params.lambda * params.rl;
    
    if (denominator <= 0 || numerator < 0) {
        return std::numeric_limits<double>::infinity();
    }
    
    double U2_squared = numerator / denominator;
    return (U2_squared > 0) ? std::sqrt(U2_squared) : 0.0;
}

double FishTailSimulation::calculate_neutral_buoyancy_depth() const {
    // y_neutral = (ρ · V₀ · p₀ / m - p₀) / (ρ · g)
    double numerator = params.rho * params.v0 * params.p0 / params.mass - params.p0;
    double denominator = params.rho * params.g;
    
    double y_neutral = numerator / denominator;
    
    return y_neutral;
}


FishTailSimulation::EquilibriumAnalysis FishTailSimulation::analyze_equilibrium_scenarios() const {
    EquilibriumAnalysis analysis;
    
    // 1. Глубина статического равновесия
    analysis.y_static = (params.rho * params.v0 * params.p0 / params.mass - params.p0) / (params.rho * params.g);
    
    // 2. Глубина динамического равновесия (решаем численно)
    analysis.y_dynamic = find_dynamic_equilibrium_depth();
    
    // 3. Анализ сценариев
    analysis.can_stabilize_from_above = analysis.y_dynamic > analysis.y_static;
    analysis.can_stabilize_from_below = false; // По гипотезе - невозможно
    
    // 4. Минимальная глубина для стабилизации
    analysis.min_depth_for_stabilization = analysis.y_static;
    
    return analysis;
}

double FishTailSimulation::find_dynamic_equilibrium_depth() const {
    // Правильное уравнение: m·g = ρ·g·V₀·p₀/(p₀ + ρ·g·y) + K·U²·α
    // Нужно учитывать угол атаки α, который зависит от глубины!
    
    double U = params.u0;
    double K = params.lambda;
    
    // Решаем численно: находим y, где сумма сил = 0
    double y_left = -10.0, y_right = 10.0;
    
    for (int i = 0; i < 50; ++i) {
        double y_mid = (y_left + y_right) / 2.0;
        
        // Выталкивающая сила на глубине y_mid
        double Fa = params.rho * params.g * params.v0 * (params.p0 / (params.p0 + params.rho * params.g * y_mid));

        double k_b = calculate_buoyancy_stiffness();
        double target_depth = 0.4; // Желаемая глубина стабилизации
        double theta_equilibrium = -k_b * (y_mid - target_depth) / (k_b * params.ra + K * U * U);
        
        double Fl = K * U * U * theta_equilibrium;
        
        // Суммарная сила
        double total_force = Fa + Fl - params.mass * params.g;
        
        if (std::abs(total_force) < 1e-6) {
            return y_mid;
        }
        
        if (total_force < 0) {
            y_left = y_mid; // Слишком много веса - увеличиваем глубину
        } else {
            y_right = y_mid; // Слишком много плавучести - уменьшаем глубину
        }
    }
    
    return (y_left + y_right) / 2.0;
}