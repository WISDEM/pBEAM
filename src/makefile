
all: py unit

py:
	cd pBEAM; make
	cd pyBEAM; make
	mv pyBEAM/_pBEAM.so .

unit:
	cd pBEAM; make
	cd test; make
	mv test/runTests .
	./runTests

clean:
	cd pBEAM; make clean
	cd pyBEAM; make clean
	cd test; make clean
	rm -f _pBEAM.so runTests
