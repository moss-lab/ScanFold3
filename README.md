# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### Usage on Pronto:  
compile_and_run.sh will load necessary modules, compile fold.cpp, and run the binary on an example(has no email set currently)  
run.sh will run the compiled program on an example under test/coronaframeshift (make sure it's been compiled)    
compile.sh will compile fold.cpp but will not run it  
output for both run.sh and compile_and_run.sh will be in the /out/ folder  
fold takes the location of the .tsv file produced by scanfold-scan that contains all window data as an argument  
### Compiling manually on Pronto  
1. load necessary modules:  
   -module load boost  
   -module load gcc  
   -module load python  
   -module load cmake  
2. cd to the build directory  
3. enter these commands:  
   -cmake ..  
   -make  
4. compiled binary will be in the /bin directory  
### Compiling manually on Linux  
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
7. compiled binary will be in the /bin directory    
### Requirements (to compile)  
Boost 1.74+  
CMake 3.15+  
Python 3.9+  
GCC version 8 or later  
