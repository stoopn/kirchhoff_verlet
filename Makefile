CXX=clang++
CPPFLAGS2= -std=c++11 -I/Users/norbert/Code/LIBCONFIG/include -I/Users/norbert/Code/mili/
LDFLAGS2=   -L/Users/norbert/Code/LIBCONFIG/lib
LLFLAGS2 = -pthread -lconfig++ -lm -lstdc++

all : wire_sim   

wire_sim : matrix.cpp wire.h wire.cpp wire_sim.cpp quaternion.cpp wire_plot.cpp
	$(CXX) $(CPPFLAGS) $(CPPFLAGS2) $(CXXFLAGS) -O3 -o wire_sim matrix.cpp wire.cpp wire_sim.cpp quaternion.cpp wire_plot.cpp $(LDFLAGS) $(LDFLAGS2) $(LIBS)  -lz $(LLFLAGS2) 


