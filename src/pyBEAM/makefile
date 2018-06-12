CXX = clang++
CFLAGS = -O2 -g -fPIC -std=c++11
#INCDIR = /opt/local/include
#LIBDIR = /opt/local/lib
#INCL = -I$(INCDIR) -I$(PY_FRAMEWORK)/include/python2.7 -I../pBEAM
INCL = -I../pBEAM -I../include `python-config --includes`
OBJ = ../pBEAM/Poly.o ../pBEAM/myMath.o ../pBEAM/BeamFEA.o ../pBEAM/Beam.o
#LIB = -L $(LIBDIR) -lboost_python-mt -lboost_numpy-mt -lpython2.7

pybeam: pyBEAM.o
#	$(CXX) -bundle -undefined dynamic_lookup -Wl,-headerpad_max_install_names $(LIB) pyBeam.o $(OBJ) -o _pBEAM.so
	$(CXX) $(INCL) $(CFLAGS) -dynamiclib -undefined dynamic_lookup -shared pyBeam.o $(OBJ) -o _pBEAM.so

.cpp.o:
	$(CXX) $(INCL) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o _pBEAM.so

# From setuptools
#C compiler: /usr/bin/clang++ -DNDEBUG -g -fwrapv -O3 -Wall -g -O2
#compile options: '-Isrc/pBEAM -I/opt/local/include -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include -I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7 -c'
#/usr/bin/clang++ -bundle -undefined dynamic_lookup -L/opt/local/lib -Wl,-headerpad_max_install_names -L/opt/local/lib/db48 -g -O2 build/temp.macosx-10.12-x86_64-2.7/src/pBEAM/Poly.o build/temp.macosx-10.12-x86_64-2.7/src/pBEAM/myMath.o build/temp.macosx-10.12-x86_64-2.7/src/pBEAM/BeamFEA.o build/temp.macosx-10.12-x86_64-2.7/src/pBEAM/Beam.o build/temp.macosx-10.12-x86_64-2.7/src/pyBEAM/pyBEAM.o -L/opt/local/lib -L/opt/local/lib/lapack/ -lboost_python-mt -lboost_numpy-mt -llapack -lblas -o build/lib.macosx-10.12-x86_64-2.7/_pBEAM.so
