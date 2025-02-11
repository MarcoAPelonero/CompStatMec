#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

struct Particle {
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

double lj_potential(double r2) {
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;
    return 4.0 * (inv_r12 - inv_r6);
}

double lj_force(double r2) {
    double r = sqrt(r2);
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;
    return 24.0 * (2.0 * inv_r12 - inv_r6) / r;
}

void computeForces(vector<Particle>& particles) {
    int n = (int)particles.size();
    for (int i = 0; i < n; i++) {
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        particles[i].fz = 0.0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;
            double r2 = dx * dx + dy * dy + dz * dz;

            if (r2 < 1e-12) continue;
            double r = sqrt(r2);
            double fmag = lj_force(r2);
            double fx = fmag * dx/r;
            double fy = fmag * dy/r;
            double fz = fmag * dz/r;

            particles[i].fx += fx;
            particles[i].fy += fy;
            particles[i].fz += fz;

            particles[j].fx -= fx;
            particles[j].fy -= fy;
            particles[j].fz -= fz;
        }
    }
}

double computePotentialEnergy(const vector<Particle>& particles) {
    double U = 0.0;
    int n = (int)particles.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;
            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < 1e-12) continue;
            U += lj_potential(r2);
        }
    }
    return U;
}

double computeKineticEnergy(const vector<Particle>& particles) {
    double K = 0.0;
    for (auto &p : particles) {
        K += 0.5 * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    }
    return K;
}

double computeForceNormSquared(const vector<Particle>& particles) {
    double sum = 0.0;
    for (auto &p : particles) {
        sum += (p.fx * p.fx + p.fy * p.fy + p.fz * p.fz);
    }
    return sum;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <timestep>" << endl;
        return 1;
    }
    int N = 50;
    double dt = atof(argv[1]);
    int steps = atoi(argv[2]);
    vector<Particle> particles(N);
    int nside = (int)ceil(cbrt(N));
    double spacing = 1.5;
    int idx = 0;
    for (int i = 0; i < nside && idx < N; i++) {
        for (int j = 0; j < nside && idx < N; j++) {
            for (int k = 0; k < nside && idx < N; k++) {
                particles[idx].x = i * spacing;
                particles[idx].y = j * spacing;
                particles[idx].z = k * spacing;
                particles[idx].vx = 0.0;
                particles[idx].vy = 0.0;
                particles[idx].vz = 0.0;
                idx++;
            }
        }
    }

    computeForces(particles);

    ofstream file("energy.dat");
    file << "# time   H_standard   H_shadow\n";

    for (int step = 0; step <= steps; step++) {
        double U = computePotentialEnergy(particles);
        double K = computeKineticEnergy(particles);
        double H = U + K;

        double pdotF = 0.0;
        for (auto &p : particles) {
            pdotF += (p.vx * p.fx + p.vy * p.fy + p.vz * p.fz);
        }
        double H_shadow = H + 0.5 * dt * pdotF;

        file << (step * dt) << " " << H << " " << H_shadow << "\n";

        for (auto &p : particles) {
            p.vx += p.fx * dt;
            p.vy += p.fy * dt;
            p.vz += p.fz * dt;

            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
        }

        computeForces(particles);
    }

    file.close();
    return 0;
}