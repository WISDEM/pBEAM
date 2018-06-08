CXX = clang++
CFLAGS = -O2 -g -fwrapv -fPIC -std=c++11
INCL = -I../pBEAM -I../
LIB = -L/opt/local/lib -L/opt/local/lib/lapack -llapack -lblas -lpython2.7
CORE_OBJ = ../pBEAM/Poly.o ../pBEAM/myMath.o ../pBEAM/BeamFEA.o ../pBEAM/Beam.o
OBJ = runTests.o testBeam.o testBeamFEA.o testMath.o testPoly.o


all: $(OBJ)
	$(CXX) $(LIB) $(CORE_OBJ) $(OBJ) -o tests

%.o : %.cpp
	$(CXX) $(CFLAGS) $(INCL) -c $<

clean:
	rm -f *.o tests
