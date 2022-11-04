#include <vector>
#include <GL/osmesa.h>
#include <GL/gl.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class ParticleRenderer {
public:
	ParticleRenderer(int width, int height)
	: width(width)
	, height(height)
	, frontbuf(width*height*4)
	, context(OSMesaCreateContext(GL_RGBA, NULL)) {
	}
	ParticleRenderer() {
		OSMesaDestroyContext(context);
	}
    void render(const std::vector<glm::vec3>& particles) {
		OSMesaMakeCurrent(context, frontbuf.data(), GL_UNSIGNED_BYTE, width, height);
        
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
		glPointSize(3);
		glColor3f(1.0f,1.0f,1.0f);
		glBegin(GL_POINTS);
		for (auto& particle : particles)
			glVertex3fv(value_ptr(particle));
		glEnd();
    	glFinish();
    }
	const std::vector<uint8_t>& frontBuffer() const { return frontbuf; }
private:
	int width, height;
	std::vector<uint8_t> frontbuf;
	OSMesaContext context;
};
