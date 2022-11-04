#include "ParticleRenderer.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>
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
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 accel;
    float density;
    float pressure;


    Particle(const Particle& copy) = default;
    Particle(glm::vec3 position)
        : time(0.0f)
        , position(position)
        , velocity(0.0f,0.0f,0.0f)
        , accel(0.0f,0.0f,0.0f)
        , density(rho0)
        , pressure(0)
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

Particle mix(const Particle& a, const Particle& b, float t)
{
    Particle res;
    res.position = glm::mix(a.position, b.position, t);
    res.velocity = glm::mix(a.velocity, b.velocity, t);
    res.accel = glm::mix(a.accel, b.accel, t);
    res.density = glm::mix(a.density, b.density, t);
    res.pressure = glm::mix(a.pressure, b.pressure, t);
    return res;
}

std::vector<Particle> particles;
std::vector<Particle> lastParticles;
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


float step(void)
{
    auto lessThan = [](int a, int b) {
        return particles[a].time < particles[b].time;
    };
    std::sort(temporalList.begin(), temporalList.end(), lessThan);

    std::vector<Particle> neighbors;
    neighbors.push_back(particles[temporalList.front()]);

    // Determine possible time step dt
    float dt = std::min(
        0.01f*std::sqrt(Particle::radius/(0.001f+length(neighbors.front().accel))),
        0.05f*Particle::radius/(0.001f+length(neighbors.front().velocity))
    );;

    // Reconstruct neighbor attributes
    for (int i = 0; i < particles.size(); i++) {
        if (i == temporalList.front())
            continue;
        
        float t = 0.0f;
        if (particles[i].time != lastParticles[i].time)
            t = (neighbors.front().time - lastParticles[i].time)
                / (particles[i].time - lastParticles[i].time);

        auto backtrack = mix(lastParticles[i], particles[i], t);
        if (distance2(backtrack.position, neighbors.front().position) < sq(Particle::radius))
            neighbors.push_back(backtrack);
    }
    Particle& particle = neighbors.front();

    // Compute density
    particle.density = 0.0f;
    for (Particle& neighbor : neighbors) {
        auto xij = particle.position - neighbor.position;
        particle.density += Particle::mass*W(xij);
    }

    // Compute F* (Fvisc + Fext)
    constexpr float nu = 0.05f;
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
    constexpr float k = 0.1f;
    auto accelP = glm::vec3{0.0f,0.0f,0.0f};
    particle.pressure = k*std::max(particle.density - Particle::rho0, 0.0f);
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

    glm::vec3 r = glm::vec3(0.5f,1.0f,0.5f);
    glm::vec3 a = 5.0f-r, b = 5.0f+r;
    for (int d = 0; d < 3; d++) {
        if (particle.position[d] < a[d]) {
            particle.position[d] = a[d]+std::min(b[d]-a[d], a[d]-particle.position[d]);
            particle.velocity[d] = 0.0f;
        }
        if (particle.position[d] > b[d]) {
            particle.position[d] = b[d]-std::min(b[d]-a[d], particle.position[d]-b[d]);
            particle.velocity[d] = 0.0f;
        }
    }
    lastParticles[temporalList.front()] = particles[temporalList.front()];
    particles[temporalList.front()] = particle;

    return particle.time;
}

void backtrackAll(std::vector<glm::vec3>& result, float time)
{
    result.clear();
    for (int i = 0; i < particles.size(); i++) {
        float t = 0.0f;
        if (particles[i].time != lastParticles[i].time)
            t = (time - lastParticles[i].time)
                / (particles[i].time - lastParticles[i].time);
        result.push_back(mix(lastParticles[i].position, particles[i].position, t));
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
            lastParticles.push_back(particles.back());
        }
    }
}


int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <numFrames> <outPath>" << std::endl;
        return EXIT_FAILURE;
    }
    const int numFrames = std::stoi(argv[1]);
    const std::string filePrefix = argv[2];

    ParticleRenderer renderer(512,512);
    packSphere({5.0f,5.0f,5.0f}, 0.25f);

    float time = 0.0f;
    std::vector<glm::vec3> frame;
    for (int k = 0; k < numFrames; k++) {
        time += 1.0f/60.0f;
        while (step() < time)
            ;
        std::cerr << ".";
        backtrackAll(frame, time);
        renderer.render(frame);

        stbi_flip_vertically_on_write(1);
        stbi_write_png((filePrefix+std::to_string(k)+".png").c_str(),
            renderer.width(), renderer.height(),
            3, renderer.frontBuffer().data(), 0);
    }
    return EXIT_SUCCESS;
}
