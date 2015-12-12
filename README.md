# CSE380_Project
## General structure
The C++ project is contained in the FEM directory.  The Matlab directory has a completely uncommented 2d FEM
prototype that was created to serve as a sort of road map.
Inside of the FEM directory there is a src/ directory for source code (cpp and hpp files) written for the
project.  There is a build directory where .o files are place.  The includes directory has the external header
files for the catch and Eigen libraries.  The bin/ directory is where the main and test binaries are output.
Also in this directory is the input file required to run main.

##Building
The program has been build and run on Stampede (in fact the grvy features don't compile on my local machine
yet).  Binaries for main and test are already built in the bin directory and should run with the masa and grvy
modules loaded.

To build, first make sure the following modules are loaded:
* gcc/4.7.1
* masa
* grvy

Navigate to the FEM folder which should contain two files, Makefile and Makefile\_test.

* Run make from the FEM folder to build main in the bin/ folder.  
* Run make -f Makefile\_test to build the catch tests into the bin/ folder
