#  q-state Potts model    
### course project 'Computational Physics I', taught by Prof. Gattringer, University of Graz, Austria, Sep. 2009 
### Andreas Windisch, andreas.windisch@yahoo.com 

## 1. Files contained in this repository

- README.md - This file.
- Makefile - Makefile to produce the executables
- potts.cpp - This program generates the lattice configurations; output file: potts.dat
- analyzer.cpp - This program analyzes tyhe data produced by potts. Input: potts.dat; Output: finaldata.dat
- Potts_Modell_(in_German).pdf - Detailed description of the project (in German).

## 2. Short description of the project

This is a naive implementation of a q-sdtate Potts model that I implemented as an undergrad student at University of Graz.
The special case of q=2 reduces the system to the Ising model. Feel free to use to code as you please. If you have any questions regarding this problem, feel free to contact me using andreas.windisch (at) yahoo.com.

## 3. Building the program

On a linux machine, just open a console, cd into the directory where the files are strored, and run:
```bash
$ make
```
This should produce two executables, 'potts' and 'analyzer'.
Then, run the first program to produce configurations and measurements on those configurations:
```bash
$ ./potts
```
A file 'potts.dat' should have been created. Now you can run the analyzer on this data:
```bash
$ ./analyzer
```
You should get a summary of magnetizations, errors, etc, for every coupling value your used in the simulation. The data is also written to a file called 'finaldata.dat'.




