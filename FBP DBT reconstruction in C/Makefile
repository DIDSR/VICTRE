#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#								DIGITAL BREAST TOMOSYNTHESIS RECONSTRUCTION (FBP)
#
# BASED ON THE C RECONSTRUCTION CODE DEVELOPED BY LEESER, MUKHERJEE AND BROCK (BMC RESEARCH NOTES 7.1, p. 582, 2014).
# WHICH IN TURN IMPLEMENTS THE FDK RECONSTRUCTION ALGORITHM DEVELOPED BY FESSLER (http://web.eecs.umich.edu/~fessler/).
#
# AUTHOR: AUNNASHA SENGUPTA
#
#
# GITHUB LINK: 	https://github.com/DIDSR/VICTRE
#
#*********************************************************************************************************************************
#*********************************************************************************************************************************

##################################################################################################################################
#
#							DISCLAIMER
#
# This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. 
# Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. 
# Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, 
# and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, 
# its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. 
# Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. 
# Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, 
# and any modified versions bear some notice that they have been modified. 
#
##################################################################################################################################


################################################################################
#
# This program compiles the C codes for DBT reconstruction.
# Usage: make
 
################################################################################

################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= ./FBP_DBTrecon


################################################################################
# Rules and targets

# Compilers
CXX        := g++
LINK       := g++ -fPIC
CFLAGS     := -g -O0 -w

default: conebeam
	$(CXX) -o $(EXECUTABLE) FBP_DBTrecon.o -lm -lfftw3f

conebeam: FBP_DBTrecon.c
	$(CXX) $(CFLAGS) -o FBP_DBTrecon.o -c FBP_DBTrecon.c

clean: 
	rm *.o $(EXECUTABLE);

