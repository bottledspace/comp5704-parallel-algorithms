CFLAGS = -O3 -static-libgcc -static-libstdc++ -idirafter /usr/include/
LIBS = -lglew32 -lglfw3 -lgdi32 -lopengl32 -lavutil -lavformat \
    -lavcodec -lswscale -lavutil -lws2_32 -lswresample \
	-lsecur32 -lbcrypt -lole32 -lstrmiids
asph2.exe : asph2.cpp VideoStream.o
	x86_64-w64-mingw32-g++ -o asph.exe $(CFLAGS) asph2.cpp VideoStream.o $(LIBS)
VideoStream.o : VideoStream.cpp VideoStream.hpp
	x86_64-w64-mingw32-g++ -c $(CFLAGS) VideoStream.cpp $(LIBS)