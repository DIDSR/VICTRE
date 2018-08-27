#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#								DIGITAL BREAST TOMOSYNTHESIS RECONSTRUCTION (FBP)
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


## INPUT ARGUMENTS (GIVEN AS COMMAND LINE ARGUMENTS SEPARATED BY A BLANK SPACE): 
# Number of projection angles
# Number of detector elements in the direction of the x-ray tube
# Number of detector elements in the direction perpendicular to the x-ray tube movement
# Detector pixel size (in cm)
# Source to detector distance (in cm)
# Source to the rotation center distance (in cm) 
# Detector center offset along the x-ray tube movement direction relative to the x-ray source center (in pixels)
# Complete orbit of projection angles
# Number of voxels in the phantom X direction
# Number of voxels in the phantom Y direction
# Number of voxels in the phantom Z direction
# Phantom voxel size (in cm)
# Reconstruction pixel size (in cm)
# Reconstruction slice thickness in Z dimension (in cm)
# Offset of the volume center to the rotation center in the X direction (in pixels)
# First projection angle in the input dataset
# Difference between the first and second projection angles
# Random seed from the phantom file name


#!/bin/bash

thisSeed="GIVE THE SEED NUMBER FOR THE PHANTOM TO RUN - TAKEN FROM PHANTOM FILE NAME"

temp3=`sed -n 's/^ *DimSize *= *//p' "PATH TO THE FOLLOWING VICTRE PHANTOM FILE/pc_"$thisSeed"_crop.mhd"`
var5=$(echo $temp3 | cut -d ' ' -f1)
var6=$(echo $temp3 | cut -d ' ' -f2)
var7=$(echo $temp3 | cut -d ' ' -f3)

runCmd="time PATH TO THE FOLLOWING FILE/FBP_DBTrecon 25 3000 1500 0.0085 65.0 60.0 0.000 50.00 $var6 $var5 $var7 0.0050 0.0085 0.1 0 -25.0 2.083333333333333333 $thisSeed"
$runCmd
