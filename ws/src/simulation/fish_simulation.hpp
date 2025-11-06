#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstdlib>

class FishTailSimulation {
private:
    // Физические параметры
    double rho = 1000.0;    // плотность воды в кг на м3
    double g = 9.80665;        // ускорение свободного падения в м/с2
    double p0 = 101325.0;   // атмосферное давление в Па
    mutable double y_curr;
    // Параметры АНПА
    double Lambda = 1.5; // Коэффицент подъемной силы. Большая подъемная сила (хорошие плавники)
    double m = 3.0; // масса
    double I = 0.5; // момент инерции
    
    // Параметры тяги
    double Ft0 = 0.1; // номинальная тяга
    double Fy = 0.2; // коэффициент усиления по вертикали
    double G = 0.0; // коэффициент усиления по углу
    
    // статическая архимедова сила
    double Fa = 0.0;

    // аэродинамическая подъемная сила, H
    double Fl = 0.0;

    // плечи сил балансированный "нос вниз"
    double ra = 0.05; // плечо архимедовой силы
    double rmg = -0.06; // плечо силы тяжести
    double rl = -0.25; // плечо подъемной силы

    // коэффициенты
    double ky = 1.0; // коэффициенты квадратичного демпфирования
    double ktheta = 5;
    double kx = 2.0;    // крупный аппарат

    // эквивалентная глубина
    double yeq = 0.344292;
    double u0 = 1.0;   // равновесная скорость

    double V0;

    double y_neutral; // глубина статического равновесия или нулевой плавучести

    double archimedes_force(double y) const;
    double thrust_force(double y, double theta) const;
    double calculate_dynamic_equilibrium(double vx) const;
    double lift_force(double vx) const;
    double find_static_equilibrium_depth() const;
    double adaptive_archimedes_force(double y, double vy) const;
    mutable double adaptive_V0 = 0.00305;
    
public:
  struct State {
      double x, y, vx, vy, theta, omega;
      double time;
      State(double x=0, double y=0, double vx=0, double vy=0, double theta=0, double omega=0, double t=0);
  };
  bool is_dynamically_stable(const FishTailSimulation::State& state, double tolerance = 0.01) const;

  FishTailSimulation();
  std::vector<double> equations_of_motion(double t, const std::vector<double>& state) const;
  std::vector<double> runge_kutta_4(double t, const std::vector<double>& state, double dt) const;
  std::vector<State> simulate(double duration, double dt, const State& initial_state);

  void save_data_to_file(const std::vector<FishTailSimulation::State>& results, 
                      const std::string& filename);
                      
  void create_gnuplot_script();
  void create_3d_gnuplot_script();
  void run_gnuplot();
  void print_simulation_results(const std::vector<FishTailSimulation::State>& results);
  double calculate_critical_depth();
  void check_initial_moment();
};