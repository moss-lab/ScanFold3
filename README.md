# ScanFold3
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

### new requirements:  
boost  
&emsp;on linux, sudo apt-get install libboost-all-dev    
&emsp;on pronto, module load boost  
### usage on Pronto:  
fold.cpp can't be compiled on Pronto afaik due to the version of gcc it uses, use the pre-compiled fold binary  
run.sh will run the binary and load necessary modules (has no email set currently)
### compiling
Compile on linux with g++ 11.4.0+
> g++ fold.cpp -o fold
