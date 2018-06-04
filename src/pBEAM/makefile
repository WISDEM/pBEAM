CXX = clang++
CFLAGS = -O2 -g -fPIC -std=c++11
INCL = -I./ -I/opt/local/include
OBJ = Poly.o myMath.o BeamFEA.o Beam.o
#LIB = -framework Accelerate

all: $(OBJ)

%.o : %.cpp
	$(CXX) $(CFLAGS) $(INCL) -c $<

clean:
	rm -f *.o
