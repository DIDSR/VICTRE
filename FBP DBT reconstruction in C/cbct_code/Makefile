################################################################################
#
# This program compiles the C codes for CT reconstruction.
# Usage: make
 
# Copyright (C) 2012 James Brock
 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= ./coneBeam
# C/C++ source files (compiled with gcc / c++)
#CCFILES		:= oclDeviceQuery.cpp


################################################################################
# Rules and targets

# Compilers
CXX        := g++
LINK       := g++ -fPIC
CFLAGS     := -g -O0 -w

default: conebeam
	$(CXX) -o $(EXECUTABLE) coneBeam.o -lm -lfftw3f

conebeam: coneBeamDBT.c
	$(CXX) $(CFLAGS) -o coneBeam.o -c coneBeamDBT.c

clean: 
	rm *.o $(EXECUTABLE);

