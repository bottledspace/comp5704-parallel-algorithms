CFLAGS = -O3
LIBS = -lGLU -lOSMesa -lm -lavutil -lavformat \
    -lavcodec -lswscale -lavutil -lswresample
asph : asph.cpp VideoStream.o
	g++ -o asph $(CFLAGS) asph.cpp VideoStream.o $(LIBS)
VideoStream.o : VideoStream.cpp VideoStream.hpp
	g++ -c $(CFLAGS) VideoStream.cpp $(LIBS)