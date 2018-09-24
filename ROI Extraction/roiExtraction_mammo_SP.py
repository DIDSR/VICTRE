#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#						ROI EXTRACTION (FFDM - SIGNAL PRESENT)
#
# AUTHOR: 	DIKSHA SHARMA
#		DIKSHA.SHARMA@FDA.HHS.GOV
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


# import numpy scientific computing package
import numpy as np
import argparse
import os

# Parsing input arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("RndSeed",type=int)			# phantom random seed - from the breast phantom filename
parser.add_argument("fs_x",type=float)			# focal spot locations X/Y/Z (in mm) (from x-ray detector coordinate system)
parser.add_argument("fs_y",type=float)
parser.add_argument("fs_z",type=float)
parser.add_argument("cr_z",type=int)			# number of Z planes to be cropped in MC-GPU output images
parser.add_argument("voxelsize",type=float)		# phantom voxel size (in mm)
parser.add_argument("pixelsize",type=float)		# detector pixel size (in mm)
parser.add_argument("orgX",type=float)          # origin X/Y/Z for cropped phantom (in mm)
parser.add_argument("orgY",type=float)
parser.add_argument("orgZ",type=float)
parser.add_argument("minvoxX",type=int)			# minimum voxel number in X,Y,Z dimensions (lowest number on cropped phantom coordinate axes)
parser.add_argument("minvoxY",type=int)
parser.add_argument("minvoxZ",type=int)
parser.add_argument("totvoxX",type=int)			# total number of voxels in X/Y/Z dimensions for cropped phantom
parser.add_argument("totvoxY",type=int)
parser.add_argument("totvoxZ",type=int)
parser.add_argument("roisizeCClus",type=int)	# ROI length in pixels for the  microcalcification cluster
parser.add_argument("roisizeSpic",type=int)		# ROI length in pixels for the spiculated mass
parser.add_argument("lendetX_pix",type=int)		# detector length in X/Y (in pixels)
parser.add_argument("lendetY_pix",type=int)		
parser.add_argument("det_z_loc",type=float)		# Z location for detector (in mm) - based on the air gap between bottom of phantom and the detector plane
args = parser.parse_args()

# change the current working directory
os.chdir('PATH TO DIRECTORY CONTAINING THIS CODE')

# read mammography projection image as float 32 bit (little endian order), MC-GPU output raw file contains 5 images but only the first one is needed
dt=np.dtype('<f4');	# little endian format for 32 bit float (real)

myimg=np.zeros(([args.lendetX_pix,args.lendetY_pix]),dtype=dt)  # flatfield corrected mammo image

fid = open("PATH TO MC-GPU OUTPUT RAW IMAGE/mcgpu_image_pcl_"+str(args.RndSeed)+"_crop.raw.gz-fatty_0000.raw","r")
tmpimgORG = np.fromfile(fid, dtype=dt, count=args.lendetX_pix*args.lendetY_pix)
myimgORG = tmpimgORG.reshape(args.lendetX_pix,args.lendetY_pix); 	# reshape for PROJECTION IMAGE X AND Y

## Read flatfield image
fidFF = open("PATH TO MC-GPU FLATFIELD RAW IMAGE/mcgpu_image-fatty-flatfield_0000.raw","r")
tmpimgFF = np.fromfile(fidFF, dtype=dt, count=args.lendetX_pix*args.lendetY_pix)
myimgFF = tmpimgFF.reshape(args.lendetX_pix,args.lendetY_pix)

## Generate flatfield corrected mammo image for breast phantom
for myxx in range(0,args.lendetX_pix):
	for myyy in range(0,args.lendetY_pix):
		myimg[myxx,myyy] = (float)(myimgFF[myxx,myyy]/myimgORG[myxx,myyy])

# read lesion locations
lesarr = np.loadtxt("PATH TO LESION LOCATIONS/pcl_"+str(args.RndSeed)+"_crop.loc", delimiter=' ', dtype = (int), unpack=True);	# voxel locations (in int)
lesloc = np.transpose(lesarr);	# because lesarr reads in col,row format. so transpose will give each row as (x,y,z,lesion type)
numloc = len(lesloc)	# number of rows

## file for writing ROI centers
floc = open("PATH TO THE FOLDER WHERE EXTRACTED ROI LOCATIONS ARE SAVED/roi_SP_"+str(args.RndSeed)+".loc",'w')

#### calculation pixels x1,y1 in projection space corresponding to the lesions

# Another way of determining these - using equation for a line
x1 = np.zeros(([1,numloc]),dtype=int);
y1 = np.zeros(([1,numloc]),dtype=int);

### ****************************************************************************************
## CONVERT EVERYTHING FROM THE PHANTOM SPACE TO THE DETECTOR SPACE
### ****************************************************************************************  
## AND IN THE END CONVERT X1,Y1 BACK TO THE PHANTOM SPACE
## The origin 0,0,0 for projection space (detector plane) is shifted by cr_z voxels. 
## this when converted to phantom space yields (dx, dy, dz) mm as the corresponding point.

for ii in range(0,numloc):
	
	lesZ = lesloc[ii,0];
	lesZ = lesZ - args.cr_z;	# REDUCE CR_Z NUMBER OF PLANES FROM lesZ

	lesX = lesloc[ii,1]      # In lesion insertion phantom is reshaped as [Z Y X] and lesions outputted as [Z X Y].  Thus X would be at [ii,1], and Y [ii,2]
	lesY = lesloc[ii,2]


	## LESION LOCATIONS FROM .LOC FILE HAS LOCATIONS IN VOXELS 		
	## CONVERT LESION VOXELS TO MM (CROPPING - PHANTOM COOR SYSTEM)
	## (lesion_vox - origin_vox)*voxelsize + origin_mm 
    	lesX_mm = (lesX-0)*args.voxelsize + args.orgX
	lesY_mm = (lesY-0)*args.voxelsize + args.orgY
	lesZ_mm = (lesZ-0)*args.voxelsize + args.orgZ


	## Focal spot (mm) and projection location (mm) both are in MCGPU coor system.
	## Translate lesion (mm) from phantom to MCGPU coor system
	## Phantom origin [0,0,0] vox or (X1,Y1,Z1) mm lies at lower left back corner of phantom.
	## Detector (mcgpu) origin lies in the same location with (0,0) mm. Remember this is not the top left corner (x=0 at chest wall) of detector.
    	## Detector origin (mm) coincides with the voxel origin (0,0,0 vox or X1,Y1,Z1 mm) of cropped phantom. 
	
	## so a point (a,b,c) mm in phantom is translated as (a-X1, b-Y1, c-Z1) mm in detector space.
	## The x,y,z axis directions are same for both phantom and detector, thus no flipping of xyz coordinates between the two systems.
	
	crop_phan_lenY_mm = args.totvoxY * args.voxelsize # cropped phantom length in Y dimension (mm)
    	det_lenY_mm = args.lendetY_pix * args.pixelsize # detector length in Y dimension (mm)

    	# translating from phantom to detector (mm)
	lesX_det_mm = lesX_mm - args.orgX
	lesY_det_mm = lesY_mm - args.orgY
	lesZ_det_mm = lesZ_mm - args.orgZ  # this indicates lesZ in detector space in mm
   

	## WITH Projection
    	alpha = (args.det_z_loc - args.fs_z)/(lesZ_det_mm - args.fs_z)   ## since detector lies at z=0. line equation is f+alpha(l-f) where alpha is (z_det-fz)/(lz-fz)
	projX_det_mm = args.fs_x + (alpha)*(lesX_det_mm - args.fs_x)   # mm
	projY_det_mm = args.fs_y + (alpha)*(lesY_det_mm - args.fs_y)   # mm

	## Calculate detector origin in pixels
	## Top left corner of detector is (0,0) pixels but this is not the origin
    	## Detector origin (0,0) mm is at x=0 pixels, y = {[(det_len_y_mm - crop_phan_len_y_mm)/2]/pixelsize} pixels
	det_orgX_pix = 0
	det_orgY_pix = ((det_lenY_mm - crop_phan_lenY_mm)*0.5)/args.pixelsize

	## Calculate distance in pixels between detector origin and projection
	dist_pix_proj_orgX = int((projX_det_mm - 0)/args.pixelsize)
	dist_pix_proj_orgY = int((projY_det_mm - 0)/args.pixelsize)

	
	## Calculate projection coordinates in pixels = detector origin pixels + distance in pixels between detector origin and projection
	proj_pix_X = dist_pix_proj_orgX + det_orgX_pix
	proj_pix_Y_tmp = dist_pix_proj_orgY + det_orgY_pix

	proj_pix_Y = args.lendetY_pix - proj_pix_Y_tmp # we figured out by looking at the voxels and pixels that Y dimension was flipped, so final Y_pix = len_det_Y - proj_Y_pix_tmp
	
	x1[0,ii] = proj_pix_X
	y1[0,ii] = proj_pix_Y
	
	print x1[0,ii], y1[0,ii]
	floc.write('%d %d %d\n' % (x1[0,ii],y1[0,ii],ii)) # write projection (or roi centers) to file x, y, index (0-7) indicating calc (0-3) or mass (4-7)



## extract roi
roipixel = np.zeros(([1,2]),dtype=int);		# number of pixels in the roi based on the lesion type {calc cluster (0) and spiculated mass (1)}

roipixel[0,0] = args.roisizeCClus;
roipixel[0,1] = args.roisizeSpic;


for ii in range(0,numloc):
	mypix = int((roipixel[0,lesloc[ii,3]] - 1)/2)	
        # roipixel[0,lesloc[ii,3]] is the number of pixels to be extracted for the roi; mypix is the number of pixels in each direction with x1,y1 at center (excluding x1,y1). for 5x5, mypix = 2.
	myroi = np.zeros(([roipixel[0,lesloc[ii,3]],roipixel[0,lesloc[ii,3]]]),dtype=dt)	# extracted roi image with the size based on the lesion type extracted

	for myx in range(x1[0,ii]-mypix,x1[0,ii]+mypix+1):	# we want x1+mypix to be included
		for myy in range(y1[0,ii]-mypix,y1[0,ii]+mypix+1):			
			myroi[abs(x1[0,ii]-mypix-myx),abs(y1[0,ii]-mypix-myy)] = myimg[myx,myy];

		

	# save as unzipped file
	if (lesloc[ii,3] == 0):		# cluster calcs
		myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR MICROCALCIFICATION CLUSTER/roiMAMMO_SP_"+str(ii)+"_"+str(args.RndSeed)+".raw")
	elif (lesloc[ii,3] == 1):	# spiculated mass
		myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR SPICULATED MASS/roiMAMMO_SP_"+str(ii)+"_"+str(args.RndSeed)+".raw")

