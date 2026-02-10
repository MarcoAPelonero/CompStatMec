// projectile_aero_magnetic_rk4.cpp
// Minimal low-DOF aerodynamic cylinder + spin + gravity + magnetic impulse, integrated with RK4.
// State: position x, velocity v, body axis unit vector e, angular momentum L.
// Aerodynamics: anisotropic quadratic drag + optional Magnus force.
// Magnetic impulse: applied AFTER one full spin period; direction is PARALLEL or PERPENDICULAR to velocity.

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "include/Vec3.hpp"
#include "include/Params.hpp"
#include "include/State.hpp"
#include "include/Physics.hpp"
#include "include/Integrator.hpp"

FieldDir parse_field_dir(const std::string& s) {
    if (s == "parallel" || s == "PARALLEL") return FieldDir::Parallel;
    if (s == "perpendicular" || s == "PERPENDICULAR" || s == "perp") return FieldDir::Perpendicular;
    std::cerr << "Unknown field_dir '" << s << "'. Use: parallel | perpendicular\n";
    std::exit(1);
}

Vec3 parse_vec3(const std::string& s) {
    // Parse "x,y,z" format
    Vec3 v;
    size_t pos1 = s.find(',');
    size_t pos2 = s.find(',', pos1 + 1);
    if (pos1 == std::string::npos || pos2 == std::string::npos) {
        std::cerr << "Invalid Vec3 format. Use: x,y,z\n";
        std::exit(1);
    }
    v.x = std::stod(s.substr(0, pos1));
    v.y = std::stod(s.substr(pos1 + 1, pos2 - pos1 - 1));
    v.z = std::stod(s.substr(pos2 + 1));
    return v;
}

int main(int argc, char** argv) {
    Params P;

    // CLI (minimal):
    //   --dir parallel|perpendicular
    //   --strength <double>
    //   --impulse_steps <int>
    // Optional:
    //   --dt <double>
    //   --t_end <double>
    //   --save_every <int>        Number of steps between saves (default: 100)
    // Output:
    //   --out <file.csv>
    std::string out_path = "output.csv";
    int save_every = 100;  // Save every N steps (reduces output size drastically)

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (i + 1 >= argc) { std::cerr << "Missing value after " << name << "\n"; std::exit(1); }
            return std::string(argv[++i]);
        };

        if (a == "--dir") {
            P.field_dir = parse_field_dir(need("--dir"));
        } else if (a == "--strength") {
            P.field_strength = std::stod(need("--strength"));
        } else if (a == "--impulse_steps") {
            P.impulse_steps = std::stoi(need("--impulse_steps"));
        } else if (a == "--dt") {
            P.dt = std::stod(need("--dt"));
        } else if (a == "--t_end") {
            P.t_end = std::stod(need("--t_end"));
        } else if (a == "--out") {
            out_path = need("--out");
        } else if (a == "--save_every") {
            save_every = std::stoi(need("--save_every"));
            if (save_every < 1) save_every = 1;
        } else if (a == "--adaptive") {
            P.adaptive_dt = true;
        } else if (a == "--tol") {
            P.adaptive_tol = std::stod(need("--tol"));
        } else {
            std::cerr << "Unknown arg: " << a << "\n";
            return 1;
        }
    }

    // Initial angular momentum corresponding to spin around initial axis e0:
    State S;
    S.x = P.x0;
    S.v = P.v0;
    S.e = clamp_unit(P.e0);
    S.L = S.e * (P.I_par * P.omega_spin0);

    ImpulseSchedule sched = make_impulse_schedule(P);

    std::ofstream f(out_path);
    if (!f) {
        std::cerr << "Cannot open output file: " << out_path << "\n";
        return 1;
    }

    f << std::setprecision(10);
    f << "t,x,y,z,vhat_x,vhat_y,vhat_z,omega_x,omega_y,omega_z,ex,ey,ez,impulse\n";

    const int steps = (P.t_end > 0.0 && P.dt > 0.0) ? int(std::ceil(P.t_end / P.dt)) : 0;

    double t = 0.0;
    int saved_count = 0;
    for (int k = 0; k <= steps; ++k) {
        // Save data at specified intervals (or first/last step)
        bool should_save = (k % save_every == 0) || (k == steps);
        
        if (should_save) {
            int imp = impulse_active(t, sched) ? 1 : 0;
            
            // Compute velocity unit vector (speed versor)
            Vec3 vhat = normalize(S.v);
            
            // Compute angular velocity from angular momentum
            Vec3 e_unit = clamp_unit(S.e);
            Vec3 omega = omega_from_L(S.L, e_unit, P.I_par, P.I_perp);

            f << t << ","
              << S.x.x << "," << S.x.y << "," << S.x.z << ","
              << vhat.x << "," << vhat.y << "," << vhat.z << ","
              << omega.x << "," << omega.y << "," << omega.z << ","
              << S.e.x << "," << S.e.y << "," << S.e.z << ","
              << imp << "\n";
            ++saved_count;
        }

        if (k == steps) break;
        S = rk4_step(P, S, t, sched);
        t += P.dt;
        
        // Check for instability
        if (is_unstable(S)) {
            std::cerr << "WARNING: Simulation became unstable at t = " << t << "\n";
            std::cerr << "Consider using --adaptive flag or smaller --dt\n";
            break;
        }
    }
    
    std::cerr << "Saved " << saved_count << " data points (save_every=" << save_every << ")\n";

    std::cerr << "Wrote " << out_path << "\n";
    std::cerr << "Impulse starts at t = " << sched.t_start
              << " and ends at t = " << sched.t_end << "\n";
    return 0;
}

/*
Build:
  g++ -O3 -std=c++17 main.cpp -o sim

Run examples:
  ./sim --dir parallel --strength 5.0 --impulse_steps 200000 --dt 1e-6 --t_end 0.2 --out parallel.csv
  ./sim --dir perpendicular --strength 5.0 --impulse_steps 200000 --dt 1e-6 --t_end 0.2 --out perp.csv
  ./sim --adaptive --dt 1e-5 --t_end 0.1 --save_every 10 --out adaptive.csv

CLI Options:
  --dir parallel|perpendicular   Parallel = ferromagnetic braking; Perpendicular = Lorentz deflection
  --strength <double>            Field strength in Tesla (default: 5.0)
  --impulse_steps <int>          Duration of field interaction in timesteps
  --dt <double>                  Integration timestep
  --t_end <double>               Simulation end time
  --save_every <int>             Save data every N steps (default: 100)
  --adaptive                     Enable adaptive sub-stepping for improved stability
  --tol <double>                 Error tolerance for adaptive mode (default: 1e-8)
  --out <file.csv>               Output file path

Notes:
- The magnetic field activates after one full spin period T = 2*pi/omega_spin0.
- PARALLEL mode: eddy-current braking F = -gamma*v. Exponential slowdown, gravity curves trajectory.
- PERPENDICULAR mode: Lorentz deflection F = F0*(v x B). Lateral deviation, speed unchanged.
- B-field direction is automatic: parallel to v0 (parallel mode) or perpendicular to v0 (perp mode).
- Drag is anisotropic and sharply penalizes transverse velocity (v_perp) via k_perp_eff.
- Use save_every to reduce output file size and improve analysis performance.
- Use --adaptive if simulation becomes unstable (NaN or extreme values).
*/
