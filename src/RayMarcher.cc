#include "SDF.hpp"
#include "Field.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <sstream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cstring>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/norm.hpp>


constexpr float near = 0.1f;
constexpr float far = 5.0f;

template <typename T> auto sq(const T& a) { return a*a; }

glm::vec3 grad(const SDF& f, glm::vec3 position) {
    float centerValue = f(position);
    float sign = std::copysign(1.0, centerValue);
    glm::vec3 result;
    for (int d = 0; d < 3; d++) {
        glm::vec3 delta(0,0,0); delta[d] = 0.01f;
        float valA = f(position-delta);
        float valB = f(position+delta);
        if (sign*valA < sign*valB)
            result[d] = centerValue - valA;
        else
            result[d] = valB - centerValue;
    }
    return result / 0.01f;
}

enum ShadingMode {
    SM_NORMAL,
    SM_DEPTH
};

void raymarch(Field<glm::u8vec3,2>& image, const SDF& sdf,
              const glm::ivec2 &dim, const glm::mat4& mv, ShadingMode sm = SM_NORMAL)
{
    constexpr float thresh = 1e-3f;

    glm::mat4 proj = glm::perspective(45.0f, float(dim.x)/dim.y, near, far);

    glm::mat4 mvp_inv = inverse(proj*mv);
    for (int y = 0; y < dim.y; y++)
    for (int x = 0; x < dim.x; x++) {
        glm::vec2 ndc = glm::vec2(2.0f*float(x)/dim.x - 1.0f, 1.0f - 2.0f*float(y)/dim.y);
        glm::vec4 orig = mvp_inv * glm::vec4(ndc,-1.0f,1.0f);
        glm::vec4 dest = mvp_inv * glm::vec4(ndc,1.0f,1.0f);
        orig /= orig.w; dest /= dest.w; // Perspective divide
        glm::vec4 dir = normalize(dest-orig);

        float dist = near;
        glm::vec3 pos;
        for (;;) {
            pos = orig+dir*dist;
            float ddist = sdf(glm::vec3(pos));
            if (ddist < thresh) {
                glm::vec4 viewPos = proj*mv*glm::vec4(pos,1.0f);
                switch (sm) {
                case SM_NORMAL: image(x,y) = (normalize(glm::vec3(mv*glm::vec4(grad(sdf, pos),0.0f)))+1.0f)*128.0f; break;
                case SM_DEPTH:  image(x,y) = glm::u8vec3(1,1,1) * uint8_t(255 * viewPos.z); break;
                }
                break;
            }
            else if (ddist > 0.25f)
                ddist = 0.25f;

            dist += ddist;
            if (dist >= far) {
                image(x,y) = glm::u8vec3(255,255,255);
                break;
            }
        }
        
    }
}

glm::mat4 parse_transformation(const std::string& src)
{
	glm::mat4 result(1.0f);

    std::istringstream ss(src);
    std::string line;
    while (std::getline(ss, line, ';')) {
        std::istringstream ls(line);
        char command = ls.get();
        if (command == 't') {
            float x,y,z;
            ls >> x >> y >> z;
            result = result * glm::translate(glm::mat4(1.0f), glm::vec3(x,y,z));
        } else if (command == 'r') {
            float theta, x,y,z;
            ls >> theta >> x >> y >> z;
            result = result * glm::rotate(glm::mat4(1.0f), theta, glm::vec3(x,y,z));
        } else if (command == 's') {
            float x,y,z;
            ls >> x >> y >> z;
            result = result * glm::scale(glm::mat4(1.0f), glm::vec3(x,y,z));
        }
    }
    return result;
}

template <typename T, int Dim>
auto stov(const std::string& str)
{
    glm::vec<Dim,T> res;
    char sep;
    std::stringstream ss(str);
    for (int i = 0; i < Dim; i++)
        ss >> res[i] >> sep;
    return res;
}

int main(int argc, char** argv)
{
    glm::mat4 mv = parse_transformation(argv[1]);

    std::unique_ptr<SDF> sdf;
    if (std::strcmp(argv[2], "sdf") == 0) {
        assert(0);
    } else if (std::strcmp(argv[2], "particle") == 0)
        sdf = std::make_unique<ParticleSDF>(argv[3]);

    glm::ivec2 dim = stov<int,2>(argv[4]);

    ShadingMode sm = SM_NORMAL;
    if (std::strcmp(argv[5], "normal") == 0)
        sm = SM_NORMAL;
    else if (std::strcmp(argv[5], "depth") == 0)
        sm = SM_DEPTH;
    else {
        std::cerr << "Unknown mode." << std::endl;
        return EXIT_FAILURE;
    }
    
    Field<glm::u8vec3,2> image(dim.x, dim.y);
    raymarch(image, *sdf, dim, mv, sm);

    stbi_write_png(argv[6], dim.x, dim.y, 3, (void *)image.data(), dim.x*3);
}