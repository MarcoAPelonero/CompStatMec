#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

#include "Vec3.hpp"
#include "Params.hpp"
#include "State.hpp"
#include "Physics.hpp"
#include <algorithm>

inline Deriv compute_deriv(const Params& P, const State& S, double t, const ImpulseSchedule& sched) {
    Deriv D;

    Vec3 F, tau;
    forces_and_torques(P, S, t, sched, F, tau);

    Vec3 e = clamp_unit(S.e);
    Vec3 omega = omega_from_L(S.L, e, P.I_par, P.I_perp);

    D.dx = S.v;
    D.dv = F / P.m;
    D.de = cross(omega, e);  // kinematics of body axis in space
    D.dL = tau;

    return D;
}

inline State add_scaled(const State& S, const Deriv& K, double h) {
    State out = S;
    out.x += K.dx * h;
    out.v += K.dv * h;
    out.e += K.de * h;
    out.L += K.dL * h;
    return out;
}

inline void renormalize(State& S) {
    S.e = clamp_unit(S.e);
}

// Single RK4 step with given step size
inline State rk4_single_step(const Params& P, const State& S, double t, double h, const ImpulseSchedule& sched) {
    Deriv k1 = compute_deriv(P, S, t, sched);
    Deriv k2 = compute_deriv(P, add_scaled(S, k1, 0.5*h), t + 0.5*h, sched);
    Deriv k3 = compute_deriv(P, add_scaled(S, k2, 0.5*h), t + 0.5*h, sched);
    Deriv k4 = compute_deriv(P, add_scaled(S, k3, 1.0*h), t + 1.0*h, sched);

    State out = S;
    out.x += (k1.dx + k2.dx*2.0 + k3.dx*2.0 + k4.dx) * (h/6.0);
    out.v += (k1.dv + k2.dv*2.0 + k3.dv*2.0 + k4.dv) * (h/6.0);
    out.e += (k1.de + k2.de*2.0 + k3.de*2.0 + k4.de) * (h/6.0);
    out.L += (k1.dL + k2.dL*2.0 + k3.dL*2.0 + k4.dL) * (h/6.0);

    renormalize(out);
    return out;
}

// Estimate local truncation error using step-doubling
inline double estimate_error(const State& S_full, const State& S_half2) {
    // Compare positions (most sensitive to error typically)
    Vec3 diff_x = S_full.x - S_half2.x;
    Vec3 diff_v = S_full.v - S_half2.v;
    
    double err_x = norm(diff_x);
    double err_v = norm(diff_v);
    
    // Weighted error (position error scaled by typical length, velocity by typical speed)
    return std::max(err_x, err_v * 1e-3);
}

// Adaptive RK4 step with substeps for stability
// Returns the new state after advancing by P.dt, using substeps if needed
inline State rk4_step(const Params& P, const State& S, double t, const ImpulseSchedule& sched) {
    if (!P.adaptive_dt) {
        // Fixed timestep (original behavior)
        return rk4_single_step(P, S, t, P.dt, sched);
    }
    
    // Adaptive: subdivide the step if error is too large
    double h = P.dt;
    double t_target = t + P.dt;
    State S_current = S;
    double t_current = t;
    
    const double tol = P.adaptive_tol;
    const double min_h = P.dt * 1e-6;  // Don't go smaller than this
    const int max_substeps = 10000;
    int substep_count = 0;
    
    while (t_current < t_target - 1e-15 && substep_count < max_substeps) {
        // Don't overshoot
        if (t_current + h > t_target) {
            h = t_target - t_current;
        }
        
        // Full step
        State S_full = rk4_single_step(P, S_current, t_current, h, sched);
        
        // Two half steps
        State S_half = rk4_single_step(P, S_current, t_current, 0.5*h, sched);
        State S_half2 = rk4_single_step(P, S_half, t_current + 0.5*h, 0.5*h, sched);
        
        double err = estimate_error(S_full, S_half2);
        
        if (err < tol || h <= min_h) {
            // Accept the step (use the more accurate half-step result)
            S_current = S_half2;
            t_current += h;
            
            // Try to increase step size for next iteration
            if (err < tol * 0.1 && h < P.dt) {
                h = std::min(h * 2.0, P.dt);
            }
        } else {
            // Reject and try smaller step
            h = std::max(h * 0.5, min_h);
        }
        ++substep_count;
    }
    
    return S_current;
}

// Check if state appears unstable (NaN or extreme values)
inline bool is_unstable(const State& S, double max_pos = 1e6, double max_vel = 1e6) {
    if (std::isnan(S.x.x) || std::isnan(S.x.y) || std::isnan(S.x.z)) return true;
    if (std::isnan(S.v.x) || std::isnan(S.v.y) || std::isnan(S.v.z)) return true;
    if (norm(S.x) > max_pos || norm(S.v) > max_vel) return true;
    return false;
}

#endif // INTEGRATOR_HPP
