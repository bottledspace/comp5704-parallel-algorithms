.PHONY : run
NUM_FRAMES := 10
CFLAGS = -O3
LIBS = -lGLU -lOSMesa -lm -lavutil -lavformat \
    -lavcodec -lswscale -lavutil -lswresample
asph : asph.cpp ParticleRenderer.hpp
	g++ -o asph $(CFLAGS) asph.cpp $(LIBS)
frames :
	mkdir frames/
run : asph | frames
	./asph $(NUM_FRAMES) frames/frame_ 
	ffmpeg -framerate 60 -pattern_type glob -i 'frames/frame_*.png' \
	-c:v libx264 -pix_fmt yuv420p out.mp4