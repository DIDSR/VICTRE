#*********************************************************************************************************************************
#*********************************************************************************************************************************
#			VICTRE (VIRTUAL CLINICAL TRIAL FOR REGULATORY EVALUATION)
#
#					ROI EXTRACTION FOR FFDM - SIGNAL ABSENT
#
# AUTHOR: 	DIKSHA SHARMA
#			DIKSHA.SHARMA@FDA.HHS.GOV
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
# Random seed taken from phantom file name
# Focal spot of the imaging detector - location X (mm)
# Focal spot of the imaging detector - location Y (mm)
# Focal spot of the imaging detector - location Z (mm)
# Number of Z planes cropped in projection FFDM images
# Voxel size for phantom (mm)
# Detector pixel size (in mm)
# Origin of the cropped phantom X (mm)
# Origin of the cropped phantom Y (mm)
# Origin of the cropped phantom Z (mm)
# Minimum voxel number in X (lowest number on phantom coordinate axes)
# Minimum voxel number in Y
# Minimum voxel number in Z 
# Total number of voxels in X (phantom)
# Total number of voxels in Y
# Total number of voxels in Z
# Lesion length (in voxels) for microcalcification cluster (considering lesion volume to be a cube)
# Lesion length (in voxels) for spiculated mass (considering lesion volume to be a cube)
# ROI length (in pixels) for the  microcalcification cluster
# ROI length (in pixels) for the spiculated mass
# Total number of ROIs to be extracted from this phantom
# Detector length in X (in pixels)
# Detector length in Y (in pixels)
# Z location for detector (in mm) - based on the air gap between bottom of the phantom and the detector plane
# Min limit on X dimension within the breast to make sure that lesions are not close to the chest wall



#!/bin/bash

thisSeed="GIVE THE SEED NUMBER FOR THE PHANTOM TO RUN - TAKEN FROM PHANTOM FILE NAME"

temp1=`sed -n 's/^ *ElementSpacing *= *//p' "PATH TO THE FOLLOWING FILE/pc_"$thisSeed"_crop.mhd"`
var1=$(echo $temp1 | cut -d ' ' -f1)

temp2=`sed -n 's/^ *Offset *= *//p' "PATH TO THE FOLLOWING FILE/pc_"$thisSeed"_crop.mhd"`
var2=$(echo $temp2 | cut -d ' ' -f1)
var3=$(echo $temp2 | cut -d ' ' -f2)
var4=$(echo $temp2 | cut -d ' ' -f3)

temp3=`sed -n 's/^ *DimSize *= *//p' "PATH TO THE FOLLOWING FILE/pc_"$thisSeed"_crop.mhd"`
var5=$(echo $temp3 | cut -d ' ' -f1)
var6=$(echo $temp3 | cut -d ' ' -f2)
var7=$(echo $temp3 | cut -d ' ' -f3)

runCmd="time python PATH TO THIS CODE/roiExtraction_mammo_SA.py $thisSeed 0 60.25 630 0 $var1 0.085 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 65 109 12 1500 3000 -20 115"
$runCmd
