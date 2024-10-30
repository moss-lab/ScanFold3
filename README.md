# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### Usage on Pronto:  
run.sh will load necessary modules, compile fold.cpp, and run the binary (has no email set currently)  
fold takes the location of the .tsv file produced by scanfold-scan that contains all window data as an argument, run.sh contains an example using the coronavirus frameshift element in the test folder
### Compiling manually  
1. ensure boost is installed  
on linux, sudo apt-get install libboost-all-dev    
on pronto, module load boost  
2. Compile on linux with g++ 11.4.0+  
on linux, sudo apt-get install g++  
on pronto, module load gcc  
3. Run this command:  
g++ fold.cpp -o fold  
### Requirements (to compile)  
Boost 1.74+  
CMake 3.22  
### Future Requirements  
Pybind11- will need to be added to the ScanFold environment  
