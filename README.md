# ScanFold3   
3.0 version of ScanFold, code from old version will be brought over as new code is finished  

## Requirements    
### To run   
Debian-based Operating System   
	Recommend Ubuntu 22+   
Anaconda 29+   
At least 7G of disc space   

Note: ScanFold3 does NOT currently have Windows support, though it is planned for the full release. If necessary, for help getting ScanFold3 to work on Windows contact us over email or Github (see page 1 for contact information). We cannot assist in getting ScanFold3 to work on MacOS.    
### To compile (optional)   
CMake 3.15+  
PyBind11 2.13+   
GCC Version 8 or later   

## Installation   
ScanFold3 uses Conda to install dependencies. The Conda environment can be created using the following steps:   
1. Download ScanFold3 from github (https://github.com/moss-lab/ScanFold3) to the desired folder   
2. Open a terminal and navigate to the env directory in ScanFold3 and install the conda environment from ScanFold3.yml using the following command:   
	```
	conda install -f ScanFold3.yml
	```   
3. When prompted, enter "y"   
	In order to run ScanFold3, the Conda environment must be active. This can be done using the command:   
	```
	conda activate ScanFold3
	```   
This only needs to be done once per session. In order to run ScanFold3 in a new terminal you must activate the environment in it.   
### Compiling ScanFold3      
ScanFold3 does not need to be compiled if you are using a compatible Debian-based operating system. If you are using a non-Debian OS or an old Debian OS, you may need to re-compile ScanFold3's libraries in order for it to work. If you encounter errors, particularly those caused by C++ name mangling (ImportError: undefined symbol: ...), re-compiling may fix them.     
In order to compile ScanFold3, first ensure you have installed the requirements listed in Requirements: To compile. Then follow these steps:   
1. In the ScanFold3 root directory, which contains the file CMakeLists.txt, create and navigate to a directory named "build":   
	```
	mkdir build    
	cd build    
	```   
2. In this directory, run the following commands:   
	```
	cmake ..   
	make
	```   
3. Navigate back to the ScanFold3 directory. (Optional) delete the build directory:   
	```
	cd ..   
	rm -r build    
	```
 
### Troubleshooting Installation   
If you are encountering errors running ScanFold3, try the following:   
1. Ensure the Conda environment is active   
2. Check Github for any updates   
3. Double check your input files to ensure they are the proper format, with no empty lines or invisible characters. If the file was created in Windows, you may need to convert the file to Unix.    
4. Re-compile ScanFold3's libraries (see Installation: Compiling ScanFold3)   
5. Delete the ScanFold3 environment and re-create it (this may be necessary if you accidentally update the environment)   
6. Check to ensure that the Conda environment is name correctly using the command:   
	```
	conda env list
	```   
    If it is not named ScanFold3, use the name listed here instead when activating the environment.   
   
## Usage   
### Running ScanFold   
1. Ensure the conda environment created during installation is active. If using default installation, use this command:   
	```
	conda activate ScanFold3
	```   
2. Run ScanFold using a fasta file input with the following command, with the path to where your fasta file and the ScanFold directory:   
	```
	python /path/to/ScanFold.py your_sequence.fasta
	```   
	   
ScanFold-Scan and ScanFold-Fold can optionally be run alone. ScanFold-Scan will produce only a .tsv file containing window data as an output. ScanFold-Fold must take a .tsv file in that format as an input.    
To run ScanFold-Scan alone:   
	```
	python /path/to/ScanFoldScan.py your_sequence.fasta
	```   
 
To run ScanFold-Fold alone:   
	```
	python /path/to/ScanFoldFold.py your_sequence.fasta --tsv scan_output.tsv
	```   
 
In the case that you have multiple sequences in one fasta, one .tsv must be given for each sequence. Ensure headers in the fasta match with the name of the .tsv file given; everything before the first "." in the .tsv file name must be the same as the header for its respective sequence. Headers that contain a "." will cause an error.    
	Example (two sequence):   
	```
	python /path/to/ScanFoldFold.py your_sequence.fasta --tsv header_1.tsv header_2.tsv
	```   
   
### Flags (optional)    
#### I/O    

--folder_name  
&nbsp;&nbsp;&nbsp;&nbsp;Name of output folder (defaults to date/time)    
--extract    
&nbsp;&nbsp;&nbsp;&nbsp;Extract structures from minus 1 or minus 2 dbn file (2 or 1); Default = 1    
--tsv    
&nbsp;&nbsp;&nbsp;&nbsp;Input tsv name from ScanFold-Scan (Only when directly running ScanFoldFold.py)	    
--id    
&nbsp;&nbsp;&nbsp;&nbsp;Name or ID of sequence being analyzed (default "UserInput")    
--es_path    
&nbsp;&nbsp;&nbsp;&nbsp;Name of extracted structures file (default "extracted_structures") 
--igv_path    
&nbsp;&nbsp;&nbsp;&nbsp;Name of IGV file (default "igv_files")    
--inforna_path    
&nbsp;&nbsp;&nbsp;&nbsp;Name of inforna file (default "inforna_structures")    
#### Scan Stage    
-s, --step    
&nbsp;&nbsp;&nbsp;&nbsp;Step size; default = 1    
-w, --window    
&nbsp;&nbsp;&nbsp;&nbsp;Window size; default = 120    
#### Fold Stage    
-f, --filter    
&nbsp;&nbsp;&nbsp;&nbsp;Z-score value for filtering output, default = -1    
-c, --competition    
&nbsp;&nbsp;&nbsp;&nbsp;Competition determine if each base must pair to one and only one base, or if bases may be reported as pairing to multiple other bases (1 for disallow competition, 0 for allow; 1 by default)    
--global_refold    
&nbsp;&nbsp;&nbsp;&nbsp;Global refold option. Refold full sequence using Zavg <-1 and <-2 base pairs    


