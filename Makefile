#=========================================
#=      Makefile for Potts model         =
#=---------------------------------------= 
#=         by Andreas Windisch           =
#=              May 2018                 =
#=========================================

CPP       = g++
CPPFLAGS  = -Wall -O3 

default:	all 

potts: potts.cpp
	$(CPP) $(CPPFLAGS) $^ -o $@ 

analyzer: analyzer.cpp
	$(CPP) $(CPPFLAGS) $^ -o $@

all: potts analyzer 
   

clean: 
	rm potts
	rm analyzer
