#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "Vec3.hpp"
#include "Params.hpp"
#include "State.hpp"
#include <cmath>

// Convert angular momentum to angular velocity for axisymmetric body
inline Vec3 omega_from_L(const Vec3& L, const Vec3& e_unit, double I_par, double I_perp) {
    double L_par = dot(L, e_unit);
    Vec3 L_par_vec = e_unit * L_par;
    Vec3 L_perp_vec = L - L_par_vec;
    Vec3 omega = e_unit * (L_par / I_par) + (L_perp_vec / I_perp);
    return omega;
}

struct ImpulseSchedule {
    double t_start = 0.0;
    double t_end = 0.0;
};

inline ImpulseSchedule make_impulse_schedule(const Params& P) {
    // Impulse begins after one full spin about axis at omega_spin0
    ImpulseSchedule S;
    double T = (P.omega_spin0 > 1e-12) ? (2.0 * M_PI / P.omega_spin0) : 0.0;
    S.t_start = T;
    S.t_end   = T + P.impulse_steps * P.dt;
    return S;
}

inline bool impulse_active(double t, const ImpulseSchedule& S) {
    return (t >= S.t_start && t < S.t_end);
}

// Auto-compute a fixed perpendicular B-field direction from initial velocity.
// Picks a direction perpendicular to v0 that gives lateral (not vertical) deflection.
inline Vec3 compute_perp_field_dir(const Params& P) {
    Vec3 v0hat = normalize(P.v0);
    // Try "up" direction (0,1,0) — gives horizontal/lateral Lorentz deflection
    Vec3 candidate(0, 1, 0);
    if (std::abs(dot(v0hat, candidate)) > 0.99) {
        candidate = Vec3(0, 0, 1);  // fallback if v0 is nearly vertical
    }
    // Remove component parallel to v0 to get a clean perpendicular
    Vec3 perp = candidate - v0hat * dot(candidate, v0hat);
    return normalize(perp);
}

inline void forces_and_torques(
    const Params& P,
    const State& S,
    double t,
    const ImpulseSchedule& sched,
    Vec3& F_total,
    Vec3& tau_total
) {
    Vec3 e = clamp_unit(S.e);

    // Gravity
    Vec3 Fg = P.g * P.m;

    // Aerodynamics: anisotropic quadratic drag
    Vec3 v_par, v_perp;
    split_parallel_perp(S.v, e, v_par, v_perp);

    double vpar_mag = norm(v_par);
    double vperp_mag = norm(v_perp);

    // Boosted transverse coefficient
    double boost = 0.0;
    if (P.v_c > 1e-12) {
        boost = P.k_perp_boost * std::pow(vperp_mag / P.v_c, P.n_boost);
    }
    double k_perp_eff = P.k_perp0 * (1.0 + boost);

    // Quadratic drag in both parallel and perpendicular directions
    Vec3 F_drag = v_par * (-P.k_par * vpar_mag) + v_perp * (-k_perp_eff * vperp_mag);

    // Optional Magnus-like term
    Vec3 omega = omega_from_L(S.L, e, P.I_par, P.I_perp);
    double sA = sin_alpha(S.v, e);
    Vec3 F_magnus = cross(omega, S.v) * (P.k_magnus * sA);

    // Aerodynamic torque (optional)
    Vec3 tau_air(0,0,0);
    if (P.k_align != 0.0 || P.k_spin != 0.0) {
        Vec3 vhat = normalize(S.v);
        if (norm(vhat) < 1e-12) vhat = e;
        tau_air += cross(e, vhat) * (-P.k_align * sA);
        tau_air += omega * (-P.k_spin);
    }

    // Magnetic interaction (active during impulse window)
    Vec3 F_mag(0,0,0);
    Vec3 tau_mag(0,0,0);
    if (impulse_active(t, sched)) {
        Vec3 bhat;  // unit direction of effective B field

        if (P.field_dir == FieldDir::Parallel) {
            // ── Ferromagnetic / eddy-current braking ──
            // F = -gamma * v  (velocity-proportional retarding force)
            // gamma = F0 * field_strength  (units: kg/s)
            // This gives exponential speed decay: v(t) = v0 * exp(-gamma*t/m)
            // so gravity progressively dominates the trajectory.
            double gamma = P.F0 * P.field_strength;
            F_mag = S.v * (-gamma);
            bhat = normalize(P.v0);  // fixed field direction along initial travel
        } else {
            // ── Lorentz-like lateral deflection ──
            // F = F0 * (v × B)   (perpendicular force, changes direction not speed)
            // B field is FIXED in space, perpendicular to initial velocity.
            bhat = compute_perp_field_dir(P);
            Vec3 B_vec = bhat * P.field_strength;
            F_mag = cross(S.v, B_vec) * P.F0;
        }

        // Magnetic alignment torque (both modes)
        tau_mag = cross(e, bhat) * (P.tau0 * P.field_strength * P.mu);
    }

    F_total = Fg + F_drag + F_magnus + F_mag;
    tau_total = tau_air + tau_mag;
}

#endif // PHYSICS_HPP
