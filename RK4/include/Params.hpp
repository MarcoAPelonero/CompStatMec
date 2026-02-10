#ifndef PARAMS_HPP
#define PARAMS_HPP

#include "Vec3.hpp"

enum class FieldDir { Parallel, Perpendicular };

struct Params {
    // Time (must resolve spin if omega is realistic)
    double dt    = 1e-6;   // 10 microseconds
    double t_end = 0.20;   // 0.2 s flight (already long at 350 m/s)

    // Bullet-like mass/inertia (rough cylinder: m~8g, r~4.5mm, L~15mm)
    double m      = 0.008;     // kg
    double I_par  = 8.1e-8;    // kg m^2
    double I_perp = 1.9e-7;    // kg m^2

    // Gravity
    Vec3 g = Vec3(0.0, -9.81, 0.0);

    // Aerodynamic drag (quadratic; anisotropic)
    // Interpreted as: F = -k_par |v_par| v_par - k_perp_eff |v_perp| v_perp
    // k ~ 0.5*rho*Cd*A  (units kg/m)
    double k_par      = 1.15e-5;   // Cd_parallel ~0.30, A~6.36e-5 m^2
    double k_perp0    = 4.6e-5;    // Cd_perp ~1.2
    double k_perp_boost = 20.0;    // strong penalty when off-axis
    double v_c        = 25.0;      // m/s (≈ 4° at 350 m/s)
    double n_boost    = 3.0;       // sharpness

    // Optional Magnus-like force
    // Start at 0; if you want visible spirals, try 1e-6..1e-5.
    double k_magnus = 1.0e-6;

    // Aerodynamic torque (optional)
    double k_align = 0.0;
    double k_spin  = 0.0;

    // Magnetic model: Parallel (braking) or Perpendicular (deflection)
    FieldDir field_dir = FieldDir::Perpendicular;
    double field_strength = 20.0;           // Tesla (strong lab magnet)
    int impulse_steps = 2000;            // covers full 0.2 s flight

    // Magnetic force / torque coefficients (scaled by field_strength)
    // Parallel mode: F = -F0 * field_strength * v  (eddy-current braking, units: kg/(T·s))
    //   Time constant: tau = m / (F0 * field_strength)
    //   Default: tau = 0.008 / (0.002*5) = 0.8 s  →  ~22% speed loss in 0.2 s
    // Perpendicular mode: F = F0 * (v × B)  (Lorentz deflection, units: C)
    //   Radius: r = m*v / (F0*B) = 280 m  →  ~14° deflection in 0.2 s
    double F0   = 2.0e-3;     // Coupling strength — tune for desired effect
    double tau0 = 5.0e-3;     // Torque coupling (N·m per T per magnetic moment)
    double mu   = 1.0;        // Effective magnetic moment (dimensionless scale)

    // Adaptive timestep control (for improved stability)
    bool adaptive_dt = false;      // Enable adaptive sub-stepping
    double adaptive_tol = 1e-8;    // Error tolerance for adaptive mode

    // Initial conditions
    Vec3 x0 = Vec3(0, 0, 0);
    Vec3 v0 = Vec3(350.0, 0.0, 0.0);  // m/s pistol-like muzzle velocity
    Vec3 e0 = Vec3(1.0, 0.0, 0.0);
    double omega_spin0 = 9000.0;      // rad/s (≈ 86k rpm)
};

#endif // PARAMS_HPP
