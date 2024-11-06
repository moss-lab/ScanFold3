# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### Usage on Pronto:  
Before you run anything else, run the install_conda_environment.sh script. This will create a ScanFold3 environment for you, after this you don't need to run it again. If you want to recompile the library and install the conda environment at the same time, run setup.sh. This also only has to be run once.  
compile_and_run.sh will load necessary modules, compile fold.cpp as a library, and run the binary with an example python script   
run.sh will run the compiled library in an example python script      
compile.sh will compile fold.cpp but will not run it  
### Compiling manually on Pronto  
1. load necessary modules:  
   -module load miniconda3  
   if you haven't installed the environment yet, install it with:
   -cd env
   -conda env create --name ScanFold3 --file=ScanFold3.yml
   activate the environment:
   -conda activate ScanFold3
3. cd to the build directory  
4. enter these commands:  
   -cmake ..  
   -make  
5. compiled binary will be in the /lib/fold directory under root  
### Compiling manually on Linux (without conda)  
1. ensure boost is installed  
  -sudo apt-get install libboost-all-dev    
2. ensure python and G++ are installed  
  -sudo apt-get install G++  
  -https://docs.python.org/3/using/unix.html  
4. ensure cmake is installed  
  -https://cliutils.gitlab.io/modern-cmake/chapters/intro/installing.html  
5. cd into the build directory  
6. enter these commands:  
   -cmake ..  
   -make  
7. compiled binary will be in the /lib/fold directory    
### Requirements (to compile)  
Boost 1.74+  
CMake 3.15+  
Python 3.13+  
GCC version 8 or later  
