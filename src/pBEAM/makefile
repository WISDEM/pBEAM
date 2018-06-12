CXX = clang++
CFLAGS = -O3 -g -fPIC -std=c++11
INCL = -I ../include
OBJ = Poly.o myMath.o BeamFEA.o Beam.o

all: $(OBJ)

%.o : %.cpp
	$(CXX) $(CFLAGS) $(INCL) -c $<

clean:
	rm -f *.o
