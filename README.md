# CSE380_Project
## General structure
The C++ project is contained in the FEM directory.  The Matlab directory has a completely uncommented 2d FEM
prototype that was created to serve as a sort of road map.
Inside of the FEM directory there is a src/ directory for source code (cpp and hpp files) written for the
project.  There is a build directory where .o files are place.  The includes directory has the external header
files for the catch and Eigen libraries.  The bin/ directory is where the main and test binaries are output.
Also in this directory is the input file required to run main.

##Building
The program has been succesfully built on stampede with the gcc/4.7.1 module, along with the MASA and GRVY
modules.  Run make from the FEM folder to build main in the bin/ folder.  Run make -f Makefile\_Test to build
the catch tests into the bin/ folder
