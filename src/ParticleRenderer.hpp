#include <vector>
#include <GL/osmesa.h>
#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class ParticleRenderer {
public:
	ParticleRenderer(int width, int height)
	: m_width(width)
	, m_height(height)
	, m_frontbuf(width*height*3)
	, m_context(OSMesaCreateContext(GL_RGB, NULL)) {
	}
	ParticleRenderer() {
		OSMesaDestroyContext(m_context);
	}
    void render(const std::vector<glm::vec3>& particles) {
		OSMesaMakeCurrent(m_context, m_frontbuf.data(), GL_UNSIGNED_BYTE,
			m_width, m_height);
        
		glMatrixMode(GL_PROJECTION);
		const glm::mat4 proj = glm::perspective(45.0f, 1.0f, 1.0f, 15.0f);
		glLoadMatrixf(value_ptr(proj));
		glMatrixMode(GL_MODELVIEW);

		glClear(GL_COLOR_BUFFER_BIT);
        const glm::mat4 mv =
			  glm::translate(glm::mat4(1.0f), glm::vec3(-5.0f,-5.0f,-7.0f))
			* glm::scale(glm::mat4(1.0f), glm::vec3(1.0f,1.0f,1.0f))
			* glm::rotate(glm::mat4(1.0f), -12.5f, glm::vec3(1.0f,0.0f,0.0f));
        glLoadMatrixf(value_ptr(mv));
		glPointSize(5);
		glColor3f(1.0f,1.0f,1.0f);
		glBegin(GL_POINTS);
		for (auto& particle : particles)
			glVertex3fv(value_ptr(particle));
		glEnd();
    	glFinish();
    }
	const std::vector<uint8_t>& frontBuffer() const { return m_frontbuf; }
	int width() const { return m_width; }
	int height() const { return m_height; }
private:
	int m_width, m_height;
	std::vector<uint8_t> m_frontbuf;
	OSMesaContext m_context;
};
