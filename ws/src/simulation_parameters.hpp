#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include <QMetaType>

struct SimulationParameters {
    double mass;
    double inertia;
    double lambda;
    double thrust;
    double fY;
    double g;
    double ra;
    double rmg;
    double rl;
    double kx;
    double ky;
    double ktheta;
    double yeq;
    double u0;
    double v0;
    double rho;
    double p0;
    double initX;
    double initY;
    double initVx;
    double initVy;
    double initTheta;
    double initOmega;
    double duration;
    double dt;

    SimulationParameters() :
        mass(3.0), inertia(0.5), lambda(1.5), thrust(0.1), fY(0.2), g(9.81),
        ra(0.05), rmg(-0.06), rl(-0.25),
        kx(2.0), ky(1.0), ktheta(5.0),
        yeq(0.344292), u0(1.0), v0(0.00305),
        rho(1000.0), p0(101325.0),
        initX(0.0), initY(0.4), initVx(0.2), initVy(0.0), initTheta(0.0), initOmega(0.0),
        duration(200.0), dt(0.1)
    {}

    void validate() const {
      if (mass <= 0) throw std::invalid_argument("Масса должна быть положительной");
      if (inertia <= 0) throw std::invalid_argument("Момент инерции должен быть положительным");
      if (duration <= 0) throw std::invalid_argument("Длительность должна быть положительной");
      if (dt <= 0) throw std::invalid_argument("Шаг времени должен быть положительным");
      if (dt >= duration) throw std::invalid_argument("Шаг времени должен быть меньше длительности");
      if (rho <= 0) throw std::invalid_argument("Плотность должна быть положительной");
    }
};

Q_DECLARE_METATYPE(SimulationParameters)

#endif // SIMULATION_PARAMETERS_HPP