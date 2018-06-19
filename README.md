A finite element method for beam-like structures.

Author: S. Andrew Ning
Revised June 2018 by Garrett E. Barter

## Prerequisites

C++ compiler (with c++11 support), NumPy

Note that any modern compiler (gcc, clang, etc) on Linux or MacOS will support the c++11 standard.  On Windows, recommend to use mingw through Anaconda or gcc on Cygwin.


## Installation

You should already have successfully installed Python (for a [single user only](http://bugs.python.org/issue5459#msg101098)), NumPy, and setuptools.

Install pBEAM with the following command.

    $ python setup.py install


To check if installation was successful run Python from the command line

    $ python

and import the module.  If no errors are issued, then the installation was successful.

    >>> import _pBEAM


## Unit Tests

pBEAM has a large range of unit tests, but they are only accessible through C++.  They are intended to test the integrity of the underying code for development purposes, rather than the python interface.  However, if you want to run the tests then change directory to `src` and run


    $ make CXX=g++

where the name of your C++ compiler should be inserted in the place of g++.  The script will build the test executable and run all tests.  The phrase "No errors detected" signifies that all the tests passed.  You can remove the remove the test executable and all object files by running

    $ make clean

For software issues please use <https://github.com/WISDEM/pBEAM/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).


## Detailed Documentation

Access the online version at <http://wisdem.github.io/pBEAM/>

