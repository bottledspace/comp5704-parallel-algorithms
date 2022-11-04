#include "WindowGL.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

const char* vertSrc = R"(#version 330
uniform mat4 proj;
uniform mat4 mv;
layout(location=0) in vec4 loc;
void main() {
	gl_Position = proj*mv*loc;
    gl_PointSize = 3.0;
}
)";
const char* fragSrc = R"(#version 330
out vec3 color;
void main() {
	color = vec3(1.0,1.0,1.0);
}
)";

template <typename Particle>
class ParticleRenderer : public WindowGL {
public:
    ParticleRenderer(int width, int height)
    : WindowGL("ParticleRenderer", width, height) {
        m_prog = createShaderProgram(vertSrc, fragSrc);
        
    	// Create a Vertex Array Object to hold our mesh data.
    	glGenVertexArrays(1, &m_particlesVAO);
    	glBindVertexArray(m_particlesVAO);
    	// The VAO will point to a single Vertex Buffer Object which will interleave
    	// vertex information.
    	glGenBuffers(1, &m_interleavedVBO);
    	glBindBuffer(GL_ARRAY_BUFFER, m_interleavedVBO);
    	// First: Spatial coordinates.
    	glEnableVertexAttribArray(0);
    	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Particle),
            (const GLvoid *)offsetof(Particle, position));
    }
    
    void render(const std::vector<Particle>& particles) const {
        glm::mat4 proj = glm::perspective(45.0f, 1.0f, 1.0f, 15.0f);
        glm::mat4 mv = glm::translate(glm::mat4(1.0f), glm::vec3(-5.0f,-5.0f,-7.0f))
			* glm::scale(glm::mat4(1.0f), glm::vec3(1.0f,1.0f,1.0f))
			* glm::rotate(glm::mat4(1.0f), -12.5f, glm::vec3(1.0f,0.0f,0.0f));
		    //* glm::rotate(glm::mat4(1.0f), -90.0f, glm::vec3(1.0f,0.0f,0.0f));

        glEnable(GL_PROGRAM_POINT_SIZE);
        
    	// Update the buffer
    	glBindBuffer(GL_ARRAY_BUFFER, m_interleavedVBO);
    	glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(Particle),
            particles.data(), GL_STATIC_DRAW);
            
		glUseProgram(m_prog);
        GLint loc;
		loc = glGetUniformLocation(m_prog, "proj");
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(proj));
		loc = glGetUniformLocation(m_prog, "mv");
		glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(mv));
		glClear(GL_COLOR_BUFFER_BIT);
		glBindVertexArray(m_particlesVAO);
		glDrawArrays(GL_POINTS, 0, particles.size());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
private:
    GLuint m_prog;
    GLuint m_particlesVAO;
    GLuint m_interleavedVBO;
};
