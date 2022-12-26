#pragma once
#include <iostream>
#include <vector>
#include <glm/glm.hpp>

template <typename T, int N>
class Grid {};

template <typename T>
class Grid<T, 2> {
public:
    Grid() = default;
    Grid(const Grid &copy) = default;
    Grid(Grid &&move) = default;
    Grid(int width, int height)
        : m_width(width)
        , m_height(height)
        , m_grid(m_width*m_height)
        {}
    Grid(int width, int height, const T& val)
        : m_width(width)
        , m_height(height)
        , m_grid(m_width*m_height, val)
        {}
    Grid &operator =(const Grid &copy) = default;

    void resize(int width, int height) {
        m_width = width;
        m_height = height;
        m_grid.resize(m_width*m_height);
    }

    auto begin() { return m_grid.begin(); }
    auto end() { return m_grid.end(); }
    int size() const { return m_grid.size(); }
    int width() const { return m_width; }
    int height() const { return m_height; }
    auto dim() const { return glm::ivec3(m_width,m_height,m_depth); }
    
    // Stop implicit cast from messing up overloads for ivec and vec
    template <typename U> U operator ()(U) = delete;
    template <typename U> U operator ()(U) const = delete;

    const T& operator()(const glm::ivec2& pt) const
        { return (*this)(pt.x, pt.y); }
    const T& operator()(int x, int y) const
        { return m_grid[(x) + (y)*m_width]; }
    T& operator()(const glm::ivec2& pt)
        { return (*this)(pt.x, pt.y); }
    T& operator()(int x, int y)
        { return m_grid[(x) + (y)*m_width]; }
    
    const T* data() const { return m_grid.data(); }
    T* data() { return m_grid.data(); }
private:
    int m_width, m_height, m_depth;
    std::vector<T> m_grid;
};

template <typename T>
class Grid<T, 3> {
public:
    Grid() = default;
    Grid(const Grid &copy) = default;
    Grid(Grid &&move) = default;
    Grid(int width, int height, int depth)
        : m_width(width)
        , m_height(height)
        , m_depth(depth)
        , m_grid(m_width*m_height*m_depth)
        {}
    Grid(int width, int height, int depth, const T& val)
        : m_width(width)
        , m_height(height)
        , m_depth(depth)
        , m_grid(m_width*m_height*m_depth, val)
        {}
    Grid &operator =(const Grid &copy) = default;

    void resize(int width, int height, int depth) {
        m_width = width;
        m_height = height;
        m_depth = depth;
        m_grid.resize(m_width*m_height*m_depth);
    }

    auto begin() { return m_grid.begin(); }
    auto end() { return m_grid.end(); }
    int size() const { return m_grid.size(); }
    int width() const { return m_width; }
    int height() const { return m_height; }
    int depth() const { return m_depth; }
    auto dim() const { return glm::ivec3(m_width,m_height,m_depth); }
    
    // Stop implicit cast from messing up overloads for ivec and vec
    template <typename U> U operator ()(U) = delete;
    template <typename U> U operator ()(U) const = delete;

    const T& operator()(const glm::ivec3& pt) const
        { return (*this)(pt.x, pt.y, pt.z); }
    const T& operator()(int x, int y, int z) const
        { return m_grid[(x) + ((y) + (z)*m_height)*m_width]; }
    T& operator()(const glm::ivec3& pt)
        { return (*this)(pt.x, pt.y, pt.z); }
    T& operator()(int x, int y, int z)
        { return m_grid[(x) + ((y) + (z)*m_height)*m_width]; }
    
    const T* data() const { return m_grid.data(); }
    T* data() { return m_grid.data(); }
private:
    int m_width, m_height, m_depth;
    std::vector<T> m_grid;
};

template <typename T, int N>
class Field {};

template <typename T>
class Field<T,2> : public Grid<T,2> {
public:
    Field() = default;
    Field(int width, int height)
        : Grid<T,2>(width, height) {}
    Field(int width, int height, const T& val)
        : Grid<T,2>(width, height, val) {}
    
    using Grid<T,2>::operator();

    const T operator()(const glm::vec3& pt) const {
        int   x = (int)pt.x, y = (int)pt.y;
        float u = pt.x - x,  v = pt.y - y;
        return        u*v        * (*this)(x+1,y+1)
             + (1.0f-u)*v        * (*this)(x,y+1)
             +        u*(1.0f-v) * (*this)(x+1,y)
             + (1.0f-u)*(1.0f-v) * (*this)(x,y);
    }
};

template <typename T>
class Field<T,3> : public Grid<T,3> {
public:
    Field() = default;
    Field(int width, int height, int depth)
        : Grid<T,3>(width, height, depth) {}
    Field(int width, int height, int depth, const T& val)
        : Grid<T,3>(width, height, depth, val) {}
    
    using Grid<T,3>::operator();

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
};

template <typename T, int N>
void deserialize(std::istream& is, Field<T,N>& field)
{
    // Read in SDF.
    int width, height, depth;
    is.read((char *)&width, sizeof(int));
    is.read((char *)&height, sizeof(int));
    is.read((char *)&depth, sizeof(int));
    field.resize(width,height,depth);
    is.read((char *)field.data(), sizeof(T)*field.size());
}