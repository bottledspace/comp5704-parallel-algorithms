#pragma once
#include <vector>
#include <algorithm>
#include <glm/glm.hpp>

template <typename T>
class Grid3 {
public:
    Grid3(const Grid3 &copy)
    : m_width(copy.m_width),
      m_height(copy.m_height),
      m_depth(copy.m_depth),
      m_grid(copy.m_grid) {
    }
    Grid3(Grid3 &&move)
    : m_width(std::exchange(move.m_width, 0)),
      m_height(std::exchange(move.m_height, 0)),
      m_depth(std::exchange(move.m_depth, 0)),
      m_grid(std::move(move.m_grid)) {
    }
    Grid3(int width, int height, int depth, const T& val = T(0.0f))
    : m_width(width),
      m_height(height),
      m_depth(depth),
      m_grid((m_width+2)*(m_height+2)*(m_depth+2), val) {
    }

    auto begin() { return m_grid.begin(); }
    auto end() { return m_grid.end(); }
    int width() const { return m_width; }
    int height() const { return m_height; }
    int depth() const { return m_depth; }
    template <typename Param>
    const T operator()(Param p) const = delete;
    const T operator()(const glm::vec3& pt) const {
        int   x = (int)pt.x, y = (int)pt.y, z = (int)pt.z;
        float u = pt.x - x,  v = pt.y - y,  w = pt.z - z;
        return        u*v*w        * (*this)(x+1,y+1,z+1)
             + (1.0f-u)*v*w        * (*this)(x,y+1,z+1)
             +        u*(1.0f-v)*w * (*this)(x+1,y,z+1)
             + (1.0f-u)*(1.0f-v)*w * (*this)(x,y,z+1)
             +        u*v*(1.0f-w)        * (*this)(x+1,y+1,z)
             + (1.0f-u)*v*(1.0f-w)        * (*this)(x,y+1,z)
             +        u*(1.0f-v)*(1.0f-w) * (*this)(x+1,y,z)
             + (1.0f-u)*(1.0f-v)*(1.0f-w) * (*this)(x,y,z);
    }
    const T &operator()(const glm::ivec3& pt) const
        { return (*this)(pt.x,pt.y,pt.z); }
    T &operator()(const glm::ivec3& pt)
        { return (*this)(pt.x,pt.y,pt.z); }
    const T &operator()(int x, int y, int z) const
        { return m_grid[y*(m_width+2) + x + (m_depth+2)*(m_width+2)*z]; }
    T &operator()(int x, int y, int z)
        { return m_grid[y*(m_width+2) + x + (m_depth+2)*(m_width+2)*z]; }
    const T *data() const { return m_grid.data(); }
    T *data() { return m_grid.data(); }
    
    Grid3 &operator =(const Grid3 &copy) {
        m_width = copy.m_width;
        m_height = copy.m_height;
        m_depth = copy.m_depth;
        m_grid = copy.m_grid;
        return *this;
    }
    
    friend void swap(Grid3<T>& a, Grid3<T>& b) {
        std::swap(a.m_width, b.m_width);
        std::swap(a.m_height, b.m_height);
        std::swap(a.m_depth, b.m_depth);
        std::swap(a.m_grid, b.m_grid);
    }
private:
    int m_width, m_height, m_depth;
    std::vector<T> m_grid;
};


void bounds(Grid3<glm::vec3>& out)
{
    const int w = out.width();
    const int h = out.height();
    const int d = out.depth();

    for (int z = 1; z <= d; z++)
    for (int y = 1; y <= h; y++) {
        out(  0,y,z) = glm::vec3(0.0f,0.2f,0.2f)*out(1,y,z);
        out(w+1,y,z) = glm::vec3(0.0f,0.2f,0.2f)*out(w,y,z);
    }
    for (int z = 1; z <= d; z++)
    for (int x = 1; x <= w; x++) {
        out(x,  0,z) = glm::vec3(0.2f,0.0f,0.2f)*out(x,1,z);
        out(x,h+1,z) = glm::vec3(0.2f,0.0f,0.2f)*out(x,h,z);
    }
    for (int y = 1; y <= h; y++)
    for (int x = 1; x <= w; x++) {
        out(x,y,0)   = glm::vec3(0.2f,0.2f,0.0f)*out(x,y,1);
        out(x,y,d+1) = glm::vec3(0.2f,0.2f,0.0f)*out(x,y,d);
    }
    out(0,0,0)     = 0.5f*(out(1,0,0)+out(0,1,0));
    out(w+1,0,0)   = 0.5f*(out(w,0,0)+out(w+1,1,0));
    out(w+1,h+1,0) = 0.5f*(out(w,h+1,0)+out(w+1,h,0));
    out(0,h+1,0)   = 0.5f*(out(1,h+1,0)+out(0,h,0));

    out(0,0,d+1)     = 0.5f*(out(1,0,d)+out(0,1,d));
    out(w+1,0,d+1)   = 0.5f*(out(w,0,d)+out(w+1,1,d));
    out(w+1,h+1,d+1) = 0.5f*(out(w,h+1,d)+out(w+1,h,d));
    out(0,h+1,d+1)   = 0.5f*(out(1,h+1,d)+out(0,h,d));
}

void bounds(Grid3<float>& out)
{
    const int w = out.width();
    const int h = out.height();
    const int d = out.depth();

    for (int z = 1; z <= d; z++)
    for (int y = 1; y <= h; y++) {
        out(  0,y,z) = out(1,y,z);
        out(w+1,y,z) = out(w,y,z);
    }
    for (int z = 1; z <= d; z++)
    for (int x = 1; x <= w; x++) {
        out(x,  0,z) = out(x,1,z);
        out(x,h+1,z) = out(x,h,z);
    }
    for (int y = 1; y <= h; y++)
    for (int x = 1; x <= w; x++) {
        out(x,y,0)   = out(x,y,1);
        out(x,y,d+1) = out(x,y,d);
    }
    out(0,0,0)     = 0.5f*(out(1,0,0)+out(0,1,0));
    out(w+1,0,0)   = 0.5f*(out(w,0,0)+out(w+1,1,0));
    out(w+1,h+1,0) = 0.5f*(out(w,h+1,0)+out(w+1,h,0));
    out(0,h+1,0)   = 0.5f*(out(1,h+1,0)+out(0,h,0));

    out(0,0,d+1)     = 0.5f*(out(1,0,d)+out(0,1,d));
    out(w+1,0,d+1)   = 0.5f*(out(w,0,d)+out(w+1,1,d));
    out(w+1,h+1,d+1) = 0.5f*(out(w,h+1,d)+out(w+1,h,d));
    out(0,h+1,d+1)   = 0.5f*(out(1,h+1,d)+out(0,h,d));
}