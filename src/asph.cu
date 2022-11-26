#include "ParticleRenderer.hpp"
#include "ScopedTimer.hpp"
#include "Grid3.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <functional>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#define GLM_FORCE_CUDA
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#ifndef M_PI
#define M_PI 3.14159
#endif

template <int N> constexpr float pow(float x) { return x*pow<N-1>(x); }
template <> constexpr float pow<1>(float x) { return x; }
constexpr float sq(float x) { return x*x; }

struct Particle {
    int nextParticle;
    glm::vec3 position;
    glm::vec3 lastPosition;
    glm::vec3 velocity;
    glm::vec3 lastVelocity;
    glm::vec3 accel;
    float time;
    float density;
    float lastDensity;
    float pressure;
    float lastPressure;


    __host__ __device__
    Particle(const Particle& copy) = default;

    __host__ __device__
    Particle(glm::vec3 position)
        : nextParticle(-1)
        , position(position)
        , lastPosition(position)
        , velocity(0.0f,0.0f,0.0f)
        , lastVelocity(velocity)
        , accel(0.0f,0.0f,0.0f)
        , time(0.0f)
        , density(rho0)
        , lastDensity(rho0)
        , pressure(0)
        , lastPressure(0)
        {}
    __host__ __device__
    Particle()
        : Particle(glm::vec3(0.0f,0.0f,0.0f))
        {}

    static constexpr float mass = 0.524f;
    static constexpr float radius = 0.05f;
    static constexpr float cs = 1400.0f; // m/s^2
    static constexpr float rho0 = 1000; // water density
    static constexpr int gamma = 7;
};

__device__
Particle backtrace(const Particle& particle, float t)
{
    Particle res(particle);
    res.position = glm::mix(particle.lastPosition, particle.position, t);
    res.velocity = glm::mix(particle.lastVelocity, particle.velocity, t);
    res.density = glm::mix(particle.lastDensity, particle.density, t);
    res.pressure = glm::mix(particle.lastPressure, particle.pressure, t);
    return res;
}

__device__
float W(glm::vec3 disp) {
    float r = length(disp);
    if (r > Particle::radius)
        return 0.0f;
    float x = 1.0f - r / Particle::radius;
    return 315.0f / (64.0f * M_PI * pow<3>(Particle::radius)) * pow<3>(x);
}

__device__
glm::vec3 dW(glm::vec3 disp) {
    float r = length(disp);
    if (r > Particle::radius)
        return glm::vec3(0.0f,0.0f,0.0f);
    if (r > 0.0f) disp /= r;
    float x = 1.0f - r / Particle::radius;
    return disp * float(-45.0f / (M_PI * pow<4>(Particle::radius)) * sq(x));
}

struct Cell {
    int firstParticle;
    float time;
    float deltaTime;
    int updateCounter;

    Cell()
        : firstParticle(-1)
        , time(0.0f)
        , deltaTime(std::numeric_limits<float>::max()) {}
};

std::vector<Particle> particles;
std::vector<Cell> cells(200*200*200);

void updateCells() {
    for (int i = 0; i < cells.size(); i++) {
        cells[i].firstParticle = -1;
        cells[i].deltaTime = 0.01f;//std::numeric_limits<float>::max();
        cells[i].time = std::numeric_limits<float>::max();
        cells[i].updateCounter = 0;
    }
    
    for (int i = 0; i < particles.size(); i++) {
        glm::ivec3 index = particles[i].position / Particle::radius;
        index = glm::min(glm::max(index, glm::ivec3(0,0,0)), glm::ivec3(200,200,200)-1);
        int k = index.z*200*200+index.y*200+index.x;
        particles[i].nextParticle = cells[k].firstParticle;
        cells[k].firstParticle = i;

        float dt = min(
            0.02f*sqrt(Particle::radius/(0.001f+length(particles[i].accel))),
            0.05f*Particle::radius/(0.001f+length(particles[i].velocity))
        );
        cells[k].deltaTime = std::min(dt, cells[k].deltaTime);
        cells[k].time = std::min(particles[i].time, cells[k].time);
    }
}

__device__
Particle stepParticle(Particle& particle_, Particle* neighbors, int numNeighbors, float dt) {
    Particle particle = particle_;
    particle.lastPosition = particle.position;
    particle.lastVelocity = particle.velocity;
    particle.lastDensity = particle.density;
    particle.lastPressure = particle.pressure;
    
    // Compute density
    particle.density = 0.0f;
    for (int i = 0; i < numNeighbors; i++) {
        auto xij = particle.position - neighbors[i].position;
        particle.density += Particle::mass*W(xij);
    }

    // Compute F* (Fvisc + Fext)
    constexpr float nu = 0.015f;
    particle.accel = glm::vec3(0.0f, -9.81f, 0.0f);
    for (int i = 0; i < numNeighbors; i++) {
        auto vij = particle.velocity - neighbors[i].velocity;
        auto xij = particle.position - neighbors[i].position;
        if (length2(xij) > 0.0) {
            particle.accel += 2.0f * nu * Particle::mass / neighbors[i].density
                * vij * dot(xij, dW(xij))
                / (dot(xij, xij) + 0.01f * sq(Particle::radius));
        }
    }

    // Compute velocity using forces
    particle.velocity += dt*particle.accel;
    
    // Compute new density
    particle.density = 0.0f;
    for (int i = 0; i < numNeighbors; i++) {
        auto xij = particle.position - neighbors[i].position;
        auto vij = particle.velocity - neighbors[i].velocity;
        particle.density += Particle::mass*W(xij);
        particle.density += dt*dot(dW(xij), vij);
    }

    // Compute pressure and pressure forces
    constexpr float kappa = 0.5f;
    auto accelP = glm::vec3{0.0f,0.0f,0.0f};
    particle.pressure = kappa*std::max(particle.density - Particle::rho0, 0.0f);
    for (int i = 0; i < numNeighbors; i++) {
        auto xij = particle.position - neighbors[i].position;
        if (length2(xij) > 0.0) {
            accelP -= dW(xij) * Particle::mass
                * (particle.pressure / sq(particle.density)
                + neighbors[i].pressure / sq(neighbors[i].density));
        }
    }
    particle.accel += accelP;

    // Integrate particle over time using dt
    particle.velocity += dt*accelP;
    particle.position += dt*particle.velocity + sq(dt)/2.0f*particle.accel;

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
    particle.time += dt;
    return particle;
}

__global__
void stepCell(Cell* cells, Cell* nextCells, Particle* particles, Particle* nextParticles, float timeLimit) {
    constexpr static glm::ivec3 deltas[] = {
        glm::ivec3(-1,-1,-1),
        glm::ivec3(-1,-1, 0),
        glm::ivec3(-1,-1,+1),
        glm::ivec3(-1, 0,-1),
        glm::ivec3(-1, 0, 0),
        glm::ivec3(-1, 0,+1),
        glm::ivec3(-1,+1,-1),
        glm::ivec3(-1,+1, 0),
        glm::ivec3(-1,+1,+1),
        glm::ivec3( 0,-1,-1),
        glm::ivec3( 0,-1, 0),
        glm::ivec3( 0,-1,+1),
        glm::ivec3( 0, 0,-1),
        glm::ivec3( 0, 0,+1),
        glm::ivec3( 0,+1,-1),
        glm::ivec3( 0,+1, 0),
        glm::ivec3( 0,+1,+1),
        glm::ivec3(+1,-1,-1),
        glm::ivec3(+1,-1, 0),
        glm::ivec3(+1,-1,+1),
        glm::ivec3(+1, 0,-1),
        glm::ivec3(+1, 0, 0),
        glm::ivec3(+1, 0,+1),
        glm::ivec3(+1,+1,-1),
        glm::ivec3(+1,+1, 0),
        glm::ivec3(+1,+1,+1),
    };

    glm::ivec3 center(90+blockIdx.x, 90+blockIdx.y, 90+blockIdx.z);
    int centerIndex = center.x+center.y*200+center.z*200*200;

    nextCells[centerIndex] = cells[centerIndex];
    nextCells[centerIndex].updateCounter = 0;
    
    Cell& centerCell = cells[centerIndex];
    if (centerCell.firstParticle == -1 || centerCell.time > timeLimit)
        return;
    
    __shared__ Particle neighbors[160];
    __shared__ int numNeighbors;
    __shared__ float neighborTime;
    __shared__ bool doStep;

    if (threadIdx.x == 0) {
        doStep = true;
        numNeighbors = 0;

        for (int i = centerCell.firstParticle; i != -1; i = particles[i].nextParticle) {
            neighbors[numNeighbors++] = particles[i];
            assert(numNeighbors != 160);
        }

        neighborTime = std::numeric_limits<float>::max();
        for (glm::ivec3 delta : deltas) {
            const glm::ivec3 index = center+delta;
            const Cell& cell = cells[index.z*200*200+index.y*200+index.x];
            if (cell.firstParticle == -1)
                continue;
            if (cell.time < centerCell.time) {
                doStep = false;
            }
            neighborTime = min(neighborTime, cell.time);
            for (int i = cell.firstParticle; i != -1; i = particles[i].nextParticle) {
                float t = (centerCell.time - cell.time + cell.deltaTime) / cell.deltaTime;
                neighbors[numNeighbors++] = backtrace(particles[i], t);
                assert(numNeighbors != 160);
            }
        }
    }
Sync:
    __syncthreads();
   

    for (int k = 0, i = centerCell.firstParticle; i != -1; i = particles[i].nextParticle, k++) {
        if (k == threadIdx.x) {
            if (doStep)
                nextParticles[i] = stepParticle(neighbors[k], neighbors, numNeighbors, centerCell.deltaTime);
            else
                nextParticles[i] = particles[i];
            break;
        }
    }

    if (threadIdx.x == 0 && doStep) {
        nextCells[centerIndex].updateCounter = 1;
        nextCells[centerIndex].time = centerCell.time+centerCell.deltaTime;
    }
}

void backtrackAll(std::vector<glm::vec3>& result, float time)
{
    result.clear();
    for (int i = 0; i < particles.size(); i++) {
        glm::ivec3 index = particles[i].position / Particle::radius;
        if (index.x < 0)    index.x = 0;
        if (index.x >= 200) index.x = 199;
        if (index.y < 0)    index.y = 0;
        if (index.y >= 200) index.y = 199;
        if (index.z < 0)    index.z = 0;
        if (index.z >= 200) index.z = 199;
        int k = index.x+index.y*200+index.z*200*200;
        float t = (time - cells[k].time + cells[k].deltaTime) / (cells[k].deltaTime);
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
            particles.emplace_back(center+pos);
        }
    }
}

struct isUpdated {
    __device__ __host__
    bool operator()(Cell cell) {
        return cell.updateCounter;
    }
};


__device__
bool operator <(Particle a, Particle b) {
    glm::ivec3 indexA = a.position/Particle::radius;
    glm::ivec3 indexB = b.position/Particle::radius;
    return (indexA.z-indexB.z)*200*200
            +(indexA.y-indexB.y)*200
            +(indexA.x-indexB.x) < 0;
}

/*__global__
void update(Cell* cells, Cell* nextCells, Particle* particles, Particle* nextParticles, int numParticles, float time) {
    for (int k = 0; k < 50; k++) {
        __syncthreads();
        stepCell<<<dim3(20,20,20),1>>>(cells, nextCells, particles, nextParticles, time);
        cudaDeviceSynchronize();
        __syncthreads();
        memcpy(cells, nextCells, 200*200*200*sizeof(Cell));
        memcpy(particles, nextParticles, numParticles*sizeof(Particle));
    }
}*/

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

    int numParticles = particles.size();
    Particle *particles_dev;
    Particle *nextParticles_dev;
    Cell *cells_dev;
    Cell *nextCells_dev;
    cudaMalloc(&particles_dev, numParticles*sizeof(Particle));
    cudaMalloc(&nextParticles_dev, numParticles*sizeof(Particle));
    cudaMalloc(&cells_dev, 200*200*200*sizeof(Cell));
    cudaMalloc(&nextCells_dev, 200*200*200*sizeof(Cell));
    cudaMemcpy(particles_dev, particles.data(), numParticles*sizeof(Particle), cudaMemcpyHostToDevice);
    cudaMemcpy(nextParticles_dev, particles_dev, numParticles*sizeof(Particle), cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();

    float time = 0.0f;
    int frameCount = 0;
    std::string filename = path+"XXXX.png";
    std::vector<glm::vec3> positions;
    while (time <= duration) {
        time += 1.0f/60.0f;
        
        ScopedTimer timer(std::to_string(frameCount));
        updateCells();
        cudaMemcpyAsync(particles_dev, particles.data(), numParticles*sizeof(Particle), cudaMemcpyHostToDevice);
        cudaMemcpyAsync(nextParticles_dev, particles_dev, numParticles*sizeof(Particle), cudaMemcpyDeviceToDevice);
        cudaMemcpyAsync(cells_dev, cells.data(), 200*200*200*sizeof(Cell), cudaMemcpyHostToDevice);
        cudaMemcpyAsync(nextCells_dev, cells_dev, 200*200*200*sizeof(Cell), cudaMemcpyDeviceToDevice);
        //step<<<1,1>>>(res_dev, cells_dev, cellsNext_dev, particles_dev, particlesNext_dev, numParticles, time);
        
        for (;;) {
            for (int k = 0; k < 50; k++) {
                stepCell<<<dim3(20,20,20),160>>>(cells_dev, nextCells_dev, particles_dev, nextParticles_dev, time);
                //cudaMemcpyAsync(cells_dev, nextCells_dev, 200*200*200*sizeof(Cell), cudaMemcpyDeviceToDevice);
                std::swap(cells_dev, nextCells_dev);
                //cudaMemcpyAsync(particles_dev, nextParticles_dev, numParticles*sizeof(Particle), cudaMemcpyDeviceToDevice);
                std::swap(particles_dev, nextParticles_dev);
            }
            //update<<<1,1>>>(cells_dev, nextCells_dev, particles_dev, nextParticles_dev, numParticles, time);
            cudaDeviceSynchronize();

            thrust::device_ptr<Cell> p = thrust::device_pointer_cast(cells_dev);
            if (thrust::find_if(thrust::device, p, p+200*200*200, isUpdated()) == p+200*200*200)
                break;
        }
        /*{
            thrust::device_ptr<Particle> p = thrust::device_pointer_cast(particles_dev);
            thrust::sort(thrust::device, p, p+numParticles);
        }*/

        cudaMemcpyAsync(particles.data(), particles_dev, numParticles*sizeof(Particle), cudaMemcpyDeviceToHost);
        cudaMemcpyAsync(cells.data(), cells_dev, 200*200*200*sizeof(Cell), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();


        backtrackAll(positions, time);
        renderer.render(positions);

        int rem = frameCount;
        for (auto p = filename.rbegin()+4; p != filename.rbegin()+8; ++p) {
            *p = '0' + (rem % 10);
            rem /= 10;
        }
        //std::cerr << "writing " << filename << " at " << time << "s" << std::endl;

        stbi_flip_vertically_on_write(1);
        stbi_write_png(filename.c_str(), renderer.width(), renderer.height(),
            3, renderer.frontBuffer().data(), 0);
        
        frameCount++;
    }
    //std::cerr << frameCount << " frames total." << std::endl;
    return EXIT_SUCCESS;
}
