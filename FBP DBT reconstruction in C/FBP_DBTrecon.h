/*
 This program is a helper program for conebeam CT using C code.
 
 Copyright (C) 2012 James Brock
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>

#define MAXSTRLEN    1000
#define M_PI         3.141592653589793
// add to pi: 23846
#define FP_TH        1e-4

typedef float DAT_TYPE;

typedef enum proc_lvl_t {PREWEIGHT=0, PREFILTER, BACKPROJ};

proc_lvl_t runlvl =PREWEIGHT;

