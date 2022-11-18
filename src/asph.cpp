#include "ParticleRenderer.hpp"
#include "ScopedTimer.hpp"
#include "Grid3.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/string_cast.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

template <int N> constexpr float pow(float x) { return x*pow<N-1>(x); }
template <> constexpr float pow<1>(float x) { return x; }
constexpr float sq(float x) { return x*x; }

struct Particle {
    float time;
    float deltaTime;
    int nextParticle;
    glm::vec3 position;
    glm::vec3 lastPosition;
    glm::vec3 velocity;
    glm::vec3 lastVelocity;
    glm::vec3 accel;
    float density;
    float lastDensity;
    float pressure;
    float lastPressure;


    Particle(const Particle& copy) = default;
    Particle(glm::vec3 position)
        : time(0.0f)
        , deltaTime(0.0f)
        , position(position)
        , lastPosition(position)
        , velocity(0.0f,0.0f,0.0f)
        , lastVelocity(velocity)
        , accel(0.0f,0.0f,0.0f)
        , density(rho0)
        , lastDensity(rho0)
        , pressure(0)
        , lastPressure(0)
        {}
    Particle()
        : Particle(glm::vec3(0.0f,0.0f,0.0f))
        {}

    static constexpr float mass = 0.524f;
    static constexpr float radius = 0.05f;
    static constexpr float cs = 1400.0f; // m/s^2
    static constexpr float rho0 = 1000; // water density
    static constexpr int gamma = 7;
};

Particle backtrace(const Particle& particle, float t)
{
    Particle res(particle);
    res.position = glm::mix(particle.lastPosition, particle.position, t);
    res.velocity = glm::mix(particle.lastVelocity, particle.velocity, t);
    res.density = glm::mix(particle.lastDensity, particle.density, t);
    res.pressure = glm::mix(particle.lastPressure, particle.pressure, t);
    return res;
}

std::vector<Particle> particles;
Grid3<int> cells(200,200,200,-1);
std::vector<int> temporalList;

float W(glm::vec3 disp) {
    float r = length(disp);
    float x = 1.0f - r / Particle::radius;
    return 315.0f / (64.0f * M_PI * pow<3>(Particle::radius)) * pow<3>(x);
}
glm::vec3 dW(glm::vec3 disp) {
    float r = length(disp);
    if (r > 0.0f) disp /= r;
    float x = 1.0f - r / Particle::radius;
    return disp * float(-45.0f / (M_PI * pow<4>(Particle::radius)) * sq(x));
}

void updateCells(void) {
    std::fill(cells.begin(), cells.end(), -1);
    for (int i = 0; i < particles.size(); i++) {
        glm::ivec3 index = particles[i].position / Particle::radius;
        particles[i].nextParticle = cells(index);
        cells(index) = i;
    }
}

void forNeighbors(const Particle& particle, std::function<void(const Particle &)> fn) {
    glm::ivec3 center = particle.position / Particle::radius;
    for (int dk = -1; dk <= 1; dk++)
    for (int dj = -1; dj <= 1; dj++)
    for (int di = -1; di <= 1; di++)
        for (int i = cells(center+glm::ivec3(di,dj,dk)); i != -1; i = particles[i].nextParticle)
            fn(particles[i]);
}

void step(int k)
{
    /*auto lessThan = [](int a, int b) {
        return particles[a].time < particles[b].time;
    };
    std::sort(temporalList.begin(), temporalList.end(), lessThan);*/

    // Determine possible time step dt
    float dt = std::min(
        0.01f*std::sqrt(Particle::radius/(0.001f+length(particles[k].accel))),
        0.05f*Particle::radius/(0.001f+length(particles[k].velocity))
    );

    // Reconstruct neighbor attributes
    bool escape = false;
    std::vector<Particle> neighbors;
    neighbors.push_back(particles[k]);
    forNeighbors(particles[k], [&](const Particle& neighbor){
        if (&particles[k] == &neighbor)
            return;

        float t = 0.0f;
        if (neighbor.deltaTime != 0)
            t = (neighbors.front().time - neighbor.time + neighbor.deltaTime)
                / (neighbor.deltaTime);

        auto backtrack = backtrace(neighbor, t);
        if (distance2(backtrack.position, neighbors.front().position) < sq(Particle::radius)) {
            if (neighbors.front().time > neighbor.time) {
                escape = true;
                return;
            }
            neighbors.push_back(backtrack);
        }
    });
    if (escape)
        return;
    
    Particle& particle = neighbors.front();
    particle.deltaTime = dt;
    particle.lastPosition = particle.position;
    particle.lastVelocity = particle.velocity;
    particle.lastDensity = particle.density;
    particle.lastPressure = particle.pressure;

    // Compute density
    particle.density = 0.0f;
    for (Particle& neighbor : neighbors) {
        auto xij = particle.position - neighbor.position;
        particle.density += Particle::mass*W(xij);
    }

    // Compute F* (Fvisc + Fext)
    constexpr float nu = 0.015f;
    particle.accel = glm::vec3(0.0f, -9.81f, 0.0f);
    for (const Particle& neighbor : neighbors) {
        auto vij = particle.velocity - neighbor.velocity;
        auto xij = particle.position - neighbor.position;
        if (length2(xij) > 0.0) {
            particle.accel += 2.0f * nu * Particle::mass / neighbor.density
                * vij * dot(xij, dW(xij))
                / (dot(xij, xij) + 0.01f * sq(Particle::radius));
        }
    }

    // Compute velocity using forces
    particle.velocity += dt*particle.accel;
    
    // Compute new density
    particle.density = 0.0f;
    for (const Particle& neighbor : neighbors) {
        auto xij = particle.position - neighbor.position;
        auto vij = particle.velocity - neighbor.velocity;
        particle.density += Particle::mass*W(xij);
        particle.density += dt*dot(dW(xij), vij);
    }

    // Compute pressure and pressure forces
    constexpr float kappa = 0.5f;
    auto accelP = glm::vec3{0.0f,0.0f,0.0f};
    particle.pressure = kappa*std::max(particle.density - Particle::rho0, 0.0f);
    for (const Particle& neighbor : neighbors) {
        auto xij = particle.position - neighbor.position;
        if (length2(xij) > 0.0) {
            accelP -= dW(xij) * Particle::mass
                * (particle.pressure / sq(particle.density)
                 + neighbor.pressure / sq(neighbor.density));
        }
    }
    particle.accel += accelP;

    // Integrate particle over time using dt
    particle.velocity += dt*accelP;
    particle.position += dt*particle.velocity + sq(dt)/2.0f*particle.accel;
    particle.time += dt;

    glm::vec3 r = glm::vec3(0.4f,0.3f,0.4f);
    glm::vec3 a = 5.0f-r, b = 5.0f+r;
    for (int d = 0; d < 3; d++) {
        if (particle.position[d] < a[d]) {
            particle.position[d] = a[d]+std::min(b[d]-a[d], a[d]-particle.position[d]);
            particle.velocity[d] *= -0.2f;
        }
        if (particle.position[d] > b[d]) {
            particle.position[d] = b[d]-std::min(b[d]-a[d], particle.position[d]-b[d]);
            particle.velocity[d] *= -0.2f;
        }
    }
    particles[k] = particle;
}

void backtrackAll(std::vector<glm::vec3>& result, float time)
{
    result.clear();
    for (int i = 0; i < particles.size(); i++) {
        float t = 0.0f;
        if (particles[i].deltaTime != 0)
            t = (time - particles[i].time + particles[i].deltaTime)
                / (particles[i].deltaTime);
        result.push_back(mix(particles[i].lastPosition, particles[i].position, t));
    }
}

void packSphere(const glm::vec3& center, float radius) {
    int r = (2.0f*radius) / Particle::radius;
    for (int z = -r; z <= r; z++)
    for (int y = -r; y <= r; y++)
    for (int x = -r; x <= r; x++) {
        auto pos = Particle::radius*(
            glm::vec3(0.0f,0.5f,0.5f)*float(x)
          + glm::vec3(0.5f,0.0f,0.5f)*float(y)
          + glm::vec3(0.5f,0.5f,0.0f)*float(z));
        if (length2(pos) < sq(radius)) {
            temporalList.push_back(particles.size());
            particles.emplace_back(center+pos);
        }
    }
}


int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <duration> <path>" << std::endl;
        return EXIT_FAILURE;
    }
    const float duration = std::stof(argv[1]);
    const std::string path = argv[2];

    ParticleRenderer renderer(512,512);
    packSphere({5.0f,5.0f,5.0f}, 0.25f);

    float time = 0.0f;
    int frameCount = 0;
    std::string filename = path+"XXXX.png";
    std::vector<glm::vec3> positions;
    while (time <= duration) {
        time += 1.0f/60.0f;
        float minTime;
        {
            ScopedTimer timer(std::to_string(frameCount));
            do {
                updateCells();
                #pragma omp parallel for
                for (int i = 0; i < particles.size(); i++)
                    step(i);
                
                minTime = particles.front().time;
                for (int i = 1; i < particles.size(); i++)
                    if (minTime > particles[i].time)
                        minTime = particles[i].time;
            } while (minTime < time);
        }
        backtrackAll(positions, time);
        renderer.render(positions);

        int rem = frameCount;
        for (auto p = filename.rbegin()+4; p != filename.rbegin()+8; ++p) {
            *p = '0' + (rem % 10);
            rem /= 10;
        }
        std::cerr << "writing " << filename << " at " << time << "s" << std::endl;

        stbi_flip_vertically_on_write(1);
        stbi_write_png(filename.c_str(), renderer.width(), renderer.height(),
            3, renderer.frontBuffer().data(), 0);
        
        frameCount++;
    }
    std::cerr << frameCount << " frames total." << std::endl;
    return EXIT_SUCCESS;
}
