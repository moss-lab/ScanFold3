# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### Usage on Pronto:  
If you want to re-compile the Python library, run the install_conda_environment.sh script. This will create a ScanFold3 environment for you which contains everything you need to compile ScanFold3, after this you don't need to run it again. compile.sh will compile fold.cpp but will not run it. If you want to recompile the library and haven't created the conda environment, run setup.sh.   
The ScanFold3 environment is necessary to compile the library, running ScanFold3 should be done inside the ScanFold2 environment. This is just due to how the environment is set up on Pronto, the release environment will merge the two.  
### Usage in Python:  
installation doesn't yet install the library, add ScanFold3 root to path then import with import lib.fold.fold  
### Compiling manually on Pronto  
1. load necessary modules:  
   -module load miniconda3  
   if you haven't installed the environment yet, install it with:
   -cd env
   -conda env create --name ScanFold3 --file=ScanFold3.yml
   activate the environment:
   -conda activate ScanFold3
3. cd to the build directory (create it if it doesn't already exist)  
4. enter these commands:  
   -cmake ..  
   -make  
5. compiled binary will be in the /lib/fold directory under root  
### Compiling manually on Linux (without conda)  
1. create the ScanFold3 environment from env/ScanFold3.yml, activate (see instructions for compiling on Pronto)      
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
Boost 1.74+ (will be removed in future versions)  
CMake 3.15+  
Python 3.13+  
PyBind11 2.13.6+  
GCC version 8 or later  
