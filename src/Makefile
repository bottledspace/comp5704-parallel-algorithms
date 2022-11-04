.PHONY : run
DURATION := 1
CFLAGS = -O3
LIBS = -lGLU -lOSMesa -lm -lavutil -lavformat \
    -lavcodec -lswscale -lavutil -lswresample
asph : asph.cpp ParticleRenderer.hpp
	g++ -o asph $(CFLAGS) asph.cpp $(LIBS)
frames :
	mkdir frames/
run : asph | frames
	./asph $(DURATION) frames/frame_ 
	ffmpeg -y -framerate 60 -pattern_type glob -i 'frames/frame_*.png' \
	-c:v libx264 -pix_fmt yuv420p out.mp4