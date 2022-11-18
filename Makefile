.PHONY : run
DURATION := 1
CFLAGS = -O3 -fopenmp
LIBS = -lOSMesa

bin :
	mkdir bin
bin/asph : src/asph.cpp src/ParticleRenderer.hpp | bin
	g++ -o bin/asph $(CFLAGS) src/asph.cpp $(LIBS)
out :
	mkdir out
run : bin/asph | out
	bin/asph $(DURATION) out/frame_ 
	ffmpeg -y -framerate 60 -pattern_type glob -i 'out/frame_*.png' \
	-c:v libx264 -pix_fmt yuv420p out.mp4
