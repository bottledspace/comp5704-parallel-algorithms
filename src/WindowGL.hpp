#pragma once
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include "stb_image_write.h"

class WindowGL {
public:
    WindowGL(const char *title, int width, int height, bool visible = true);

    void screenshot(const char *filename) const;
	bool shouldClose() const { return glfwWindowShouldClose(window); }

    static GLuint compileShader(GLenum type, const char* src);
    static GLuint createShaderProgram(const char* vertSrc, const char* fragSrc = nullptr);
protected:
    int m_width;
    int m_height;
    GLFWwindow* window;
};

WindowGL::WindowGL(const char *title, int width, int height, bool visible)
{
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW." << std::endl;
        return;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GL_TRUE);
    glfwWindowHint(GLFW_VISIBLE, visible);

    window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create window." << std::endl;
        return;
    }
    glfwGetFramebufferSize(window, &m_width, &m_height);
    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW." << std::endl;
        return;
    }
}

GLuint WindowGL::compileShader(GLenum type, const char* src)
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

GLuint WindowGL::createShaderProgram(const char* vertSrc, const char* fragSrc)
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

void WindowGL::screenshot(const char* filename) const
{
    std::vector<GLubyte> pixels(m_width*m_height*3);
    glReadPixels(0,0, m_width,m_height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    stbi_flip_vertically_on_write(1);
    stbi_write_png(filename, m_width, m_height, 3, pixels.data(), 0);
}