#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/component_wise.hpp>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

enum ShadingMode {
    SM_NORMAL,
    SM_DEPTH
};

class Rasterizer {
public:
	static bool setup(int argc, char** argv);

	void rasterize(const glm::mat4& mvp, const std::vector<glm::vec3>& tris, ShadingMode sm = SM_NORMAL) const;
    void screenshot(const char *filename) const;

    static GLuint compileShader(GLenum type, const char* src);
    static GLuint createShaderProgram(const char* vertSrc, const char* fragSrc = nullptr);

	Rasterizer(int width, int height, const char *title = nullptr);

private:
	int m_width, m_height;
	int m_prog;
	GLuint m_vao, m_vbo;
	GLFWwindow* m_window;

	static const char* vertSrc;
	static const char* fragSrc;
};


// Shader based on https://math.hws.edu/graphicsbook/c7/s2.html
const char* Rasterizer::vertSrc = R"(#version 330
uniform mat4 mv;
uniform mat4 proj;
layout(location=0) in vec4 locIn;
layout(location=1) in vec3 normIn;
out vec3 color;
out float depth;
void main() {
	color = (normalize(vec3(mv*vec4(normIn,0.0)))+1.0)*0.5;
    gl_Position = proj*mv*locIn;
	depth = (gl_Position.z);
}
)";
const char* Rasterizer::fragSrc = R"(#version 330
uniform int mode;
in float depth;
in vec3 color;
out vec3 colorOut;
void main() {
	if (mode == 0)
		colorOut = color;
	else
		colorOut = vec3(depth,depth,depth);
}
)";

bool Rasterizer::setup(int argc, char** argv)
{
	return glfwInit();
}

Rasterizer::Rasterizer(int width, int height, const char *title)
	: m_width(width), m_height(height), m_prog(0), m_vao(0), m_vbo(0), m_window(nullptr)
{
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_FALSE);
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);

    m_window = glfwCreateWindow(m_width, m_height, "MeshRenderer", NULL, NULL);
    if (!m_window) {
        std::cerr << "Failed to create window." << std::endl;
        return;
    }
    glfwGetFramebufferSize(m_window, &m_width, &m_height);
    glfwMakeContextCurrent(m_window);

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW." << std::endl;
        return;
    }

	m_prog = createShaderProgram(vertSrc, fragSrc);

    // Create a Vertex Array Object to hold our mesh data.
    glGenVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);
    // The VAO will point to a single Vertex Buffer Object which will interleave
    // vertex information.
    glGenBuffers(1, &m_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    // First: Spatial coordinates.
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3)*2, (void *)0);
    // Second: Normal coordinates.
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3)*2, (void *)sizeof(glm::vec3));
}

GLuint Rasterizer::compileShader(GLenum type, const char* src)
{
	// Compile shader
	GLuint sh = glCreateShader(type);
	glShaderSource(sh, 1, &src, NULL);
	glCompileShader(sh);

	// Handle errors and warnings
	GLint loglen, stat;
	glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &loglen);
	if (loglen > 0) {
		std::vector<char> buffer(loglen + 1);
		glGetShaderInfoLog(sh, loglen, NULL, buffer.data());
		std::cerr << buffer.data();
	}
	glGetShaderiv(sh, GL_COMPILE_STATUS, &stat);
	if (!stat) {
		std::cerr << "Failed to compile shader program." << std::endl;
		glDeleteShader(sh);
		return 0;
	}
	return sh;
}

GLuint Rasterizer::createShaderProgram(const char* vertSrc, const char* fragSrc)
{
	GLuint prog, vert, frag;
	prog = glCreateProgram();
	vert = compileShader(GL_VERTEX_SHADER, vertSrc);
	if (fragSrc)
		frag = compileShader(GL_FRAGMENT_SHADER, fragSrc);
	glAttachShader(prog, vert);
	if (fragSrc)
		glAttachShader(prog, frag);
	glLinkProgram(prog);
	glDeleteShader(vert);
	if (fragSrc)
		glDeleteShader(frag);

	// Handle errors and warnings
	GLint loglen, stat;
	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &loglen);
	if (loglen > 0) {
		std::vector<char> buffer(loglen + 1);
		glGetProgramInfoLog(prog, loglen, NULL, buffer.data());
		std::cerr << buffer.data();
	}
	glGetProgramiv(prog, GL_LINK_STATUS, &stat);
	if (!stat) {
		std::cerr << "Failed to link shader program." << std::endl;
		glDeleteProgram(prog);
		return 0;
	}
	return prog;
}

void Rasterizer::screenshot(const char* filename) const
{
    std::vector<GLubyte> pixels(m_width*m_height*3);
    glReadPixels(0,0, m_width,m_height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    stbi_flip_vertically_on_write(1);
    stbi_write_png(filename, m_width, m_height, 3, pixels.data(), 0);
}

void Rasterizer::rasterize(const glm::mat4& mv, const std::vector<glm::vec3>& tris, ShadingMode sm) const
{
	constexpr int num_attrs = 2;

	// Update the buffer
    glBufferData(GL_ARRAY_BUFFER, tris.size()*sizeof(glm::vec3),
        tris.data(), GL_STATIC_DRAW);
    
    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(m_prog);
    GLint loc;
	loc = glGetUniformLocation(m_prog, "mode");
	glUniform1i(loc, sm);
    loc = glGetUniformLocation(m_prog, "mv");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mv));
    loc = glGetUniformLocation(m_prog, "proj");
    auto proj = glm::perspective(45.0f, 1.0f, 0.1f, 5.0f);
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(proj));
    glBindVertexArray(m_vao);
    glDrawArrays(GL_TRIANGLES, 0, tris.size() / num_attrs);
    glFinish();
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

void parse_obj(std::istream& is, std::vector<glm::vec3>& tris)
{
    std::vector<glm::vec3> verts;
	std::string line;
	while (std::getline(is, line)) {
		std::stringstream ls(line);
		std::string cmd;

		if (!(ls >> cmd))
			continue; // Empty line
		else if (cmd == "v") {
			float x, y, z;
			ls >> x >> y >> z;
			verts.emplace_back(x,y,z);
		} else if (cmd == "f") {
			int i;
			for (int n = 0; n < 3; n++) {
				ls >> i;
				tris.push_back(verts[i-1]);
			}
			if (ls >> i) {
				int n = tris.size();
				tris.insert(tris.end(), {tris[n-3],tris[n-1],verts[i-1]});
			}
		}
	}
}

void add_normals(std::vector<glm::vec3>& tris)
{
	std::vector<glm::vec3> tris2;
	tris2.reserve(tris.size()*2);
	for (int i = 0; i < tris.size(); i += 3) {
		glm::vec3 normal = normalize(cross(tris[i+1]-tris[i], tris[i+2]-tris[i]));
		tris2.insert(tris2.end(), {tris[i], normal, tris[i+1], normal, tris[i+2], normal});
	}
	std::swap(tris2,tris);
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

void normalize(std::vector<glm::vec3>& verts)
{
    glm::vec3 lb = glm::vec3(1.0f,1.0f,1.0f)*std::numeric_limits<float>::max();
    glm::vec3 ub = glm::vec3(1.0f,1.0f,1.0f)*std::numeric_limits<float>::min();
    for (auto& vert : verts) {
        lb = glm::min(lb, vert);
        ub = glm::max(ub, vert);
    }
    float scale = glm::compMax(ub - lb);
    glm::vec3 center = (ub + lb)/2.0f;
    for (auto& vert : verts)
        vert = 2.0f*(vert - center)/scale;
}

int main(int argc, char** argv)
{
	assert(Rasterizer::setup(argc, argv));

	glm::ivec2 dim = stov<int,2>(argv[1]);
	Rasterizer viewer(dim.x, dim.y);

	std::vector<glm::vec3> tris;
	parse_obj(std::cin, tris);
	add_normals(tris);

	glm::mat4 mvp = parse_transformation(argv[2]);
	ShadingMode sm;
    if (std::strcmp(argv[3], "normal") == 0)
        sm = SM_NORMAL;
    else if (std::strcmp(argv[3], "depth") == 0)
        sm = SM_DEPTH;
    else {
        std::cerr << "Unknown mode." << std::endl;
        return EXIT_FAILURE;
    }
	viewer.rasterize(mvp, tris, sm);
	viewer.screenshot(argv[4]);
}