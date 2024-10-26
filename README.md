# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### usage on Pronto:  
run.sh will load necessary modules, compile fold.cpp, and run the binary (has no email set currently)
fold takes the location of a .tsv file produced by scanfold-scan that contains all window data as an argument, run.sh contains an example
### compiling manually  
1. ensure boost is installed  
on linux, sudo apt-get install libboost-all-dev    
on pronto, module load boost  
2. Compile on linux with g++ 11.4.0+  
on linux, sudo apt-get install g++  
on pronto, module load gcc  
3. Run this command:  
g++ fold.cpp -o fold  
