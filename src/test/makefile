CXX = clang++
CFLAGS = -O3 -g -fPIC -std=c++11
INCL = -I../pBEAM -I../
CORE_OBJ = ../pBEAM/Poly.o ../pBEAM/myMath.o ../pBEAM/BeamFEA.o ../pBEAM/Beam.o
OBJ = runTests.o testBeam.o testBeamFEA.o testMath.o testPoly.o


all: $(OBJ)
	$(CXX) $(CORE_OBJ) $(OBJ) -o runTests

%.o : %.cpp
	$(CXX) $(CFLAGS) $(INCL) -c $<

clean:
	rm -f *.o runTests
