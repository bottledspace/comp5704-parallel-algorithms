.PHONY : run
DURATION := 1
CXX := g++
CXXFLAGS := -O3 -fopenmp --std=c++17
LIBS := -lOSMesa

bin :
	mkdir bin
bin/asph : src/asph.cpp src/ParticleRenderer.hpp | bin
	$(CXX) -o bin/asph $(CXXFLAGS) src/asph.cpp $(LIBS)
bin/raymarcher : src/RayMarcher.cc src/SDF.hpp src/Field.hpp | bin
	$(CXX) $(CXXFLAGS) src/RayMarcher.cc $(LDFLAGS) -o bin/raymarcher
bin/rasterizer : src/Rasterizer.cc src/SDF.hpp src/Field.hpp | bin
	$(CXX) $(CXXFLAGS) src/RayMarcher.cc $(LDFLAGS) -o bin/raymarcher
bin/raymarcher : src/RayMarcher.cc src/SDF.hpp src/Field.hpp | bin
	$(CXX) $(CXXFLAGS) src/RayMarcher.cc $(LDFLAGS) -o bin/raymarcher
out :
	mkdir out
run : bin/asph | out
	bin/asph $(DURATION) out/frame_ 
	ffmpeg -y -framerate 60 -pattern_type glob -i 'out/frame_*.png' \
	-c:v libx264 -pix_fmt yuv420p out.mp4
