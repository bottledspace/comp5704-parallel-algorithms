#pragma once

#include "Field.hpp"
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#if SDF_USE_JULIA
#include <julia.h>
JULIA_DEFINE_FAST_TLS()
#endif

struct SDF {
    virtual float operator()(const glm::vec3&) const = 0;
};

#if SDF_USE_JULIA
struct JuliaSDF : public SDF {
    static bool init();
    JuliaSDF(const char*);
    float operator()(const glm::vec3&) const;
private:
    jl_value_t* m_func;
};

bool JuliaSDF::init() {
    jl_init_with_image("/lib/x86_64-linux-gnu/julia/", "sys.so");
    //atexit(jl_atexit_hook);
    return true;
}

JuliaSDF::JuliaSDF(const char* src)
: m_func(jl_eval_string(src)) {
    jl_set_global(jl_main_module, jl_symbol("func"), m_func);
}

float JuliaSDF::operator()(const glm::vec3& r) const {
    jl_value_t** args;
    JL_GC_PUSHARGS(args, 3);
    args[0] = jl_box_float32(r.x);
    args[1] = jl_box_float32(r.y);
    args[2] = jl_box_float32(r.z);
    args[3] = jl_call3(m_func, args[0], args[1], args[2]);
    float res = jl_unbox_float64(args[3]);
    JL_GC_POP();
    return res;
}
#endif // SDF_USE_JULIA


struct ParticleSDF : public SDF {
    struct Particle {
        glm::vec3 loc;
        int next;
    };

    ParticleSDF(const char*, float=0.01f, float=0.1f);
    float operator()(const glm::vec3&) const;
private:
    float m_particleSize;
    float m_smoothing;
    int m_gridSize;
    float m_scale;
    Grid<int,3> m_cells;
    std::vector<Particle> m_particles;
};

ParticleSDF::ParticleSDF(const char* filename, float particleSize, float smoothing)
: m_particleSize(particleSize)
, m_smoothing(smoothing) {
    std::ifstream is(filename);
    if (!is)
        return;
    std::string line;
    while (std::getline(is, line)) {
        std::stringstream ls(line);
        std::string cmd;
        if (ls >> cmd && cmd == "v") {
            float x,y,z;
            ls >> x >> y >> z;
            m_particles.push_back({glm::vec3(x,y,z),-1});
        }
    }
    m_gridSize = ceil(10.0f/(100.0f));
    m_cells.resize(m_gridSize,m_gridSize,m_gridSize);
    std::fill_n(m_cells.begin(), m_cells.size(), -1);
    for (int i = 0; i < m_particles.size(); i++) {
        glm::ivec3 pos = m_particles[i].loc/((10.0f));
        if (pos.x < 0 || pos.x >= m_gridSize
         || pos.y < 0 || pos.y >= m_gridSize
         || pos.z < 0 || pos.z >= m_gridSize)
            continue;
        m_particles[i].next = m_cells(pos);
        m_cells(pos) = i;
    }
}

// Modified from https://iquilezles.org/articles/distfunctions/
float opSmoothUnion(float d1, float d2, float k) {
    float h = glm::clamp(0.5 + 0.5*(d2-d1)/k, 0.0, 1.0);
    return glm::mix(d2, d1, h) - k*h*(1.0-h);
}

float ParticleSDF::operator()(const glm::vec3& r) const {
    /*glm::ivec3 cell = r/m_smoothing;
    float dist = distance(m_particles[m_cells(cell)].loc, r)-m_particleSize;
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++)
    for (int i = m_cells(cell+glm::ivec3(dx,dy,dz)); i != -1; i = m_particles[i].next)
        dist = opSmoothUnion(dist, distance(m_particles[i].loc, r)-m_particleSize, m_smoothing);*/
    
    glm::ivec3 center = r/((10.0f));
    int radius = std::ceil((10.0f)/m_smoothing);
    float dist = 100.0f;
    for (int dx = -radius; dx <= radius; dx++)
    for (int dy = -radius; dy <= radius; dy++)
    for (int dz = -radius; dz <= radius; dz++) {
        glm::ivec3 pos = center+glm::ivec3(dx,dy,dz);
        if (pos.x < 0 || pos.x >= m_gridSize
         || pos.y < 0 || pos.y >= m_gridSize
         || pos.z < 0 || pos.z >= m_gridSize)
            continue;
        for (int i = m_cells(pos); i != -1; i = m_particles[i].next)
            dist = opSmoothUnion(dist, distance(m_particles[i].loc, r)-m_particleSize, m_smoothing);
    }
    return dist;
}
