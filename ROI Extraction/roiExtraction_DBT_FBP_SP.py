#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#						ROI EXTRACTION (DBT - SIGNAL PRESENT)
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
parser.add_argument("RndSeed",type=int)					# phantom random seed - from the breast phantom filename
parser.add_argument("fs_x",type=float)					# focal spot locations X/Y/Z (in mm) (from x-ray detector coordinate system)
parser.add_argument("fs_y",type=float)
parser.add_argument("fs_z",type=float)
parser.add_argument("cr_z",type=int)					# number of Z planes to be cropped in MC-GPU output images
parser.add_argument("dim_x",type=int)					# number of pixels in X/Y/Z dimensions in the reconstructed raw image
parser.add_argument("dim_y",type=int)
parser.add_argument("dim_z",type=int)					# total number of slices
parser.add_argument("pixelsize",type=float)				# detector pixel size (in mm)
parser.add_argument("sl_thick",type=float)				# reconstruction slice thickness (in mm)
parser.add_argument("orgX",type=float)          		# phantom origin in X/Y/Z dimensions (in mm)
parser.add_argument("orgY",type=float)
parser.add_argument("orgZ",type=float)
parser.add_argument("minvoxX",type=int)         		# minimum voxel number in X,Y,Z dimensions (lowest number on cropped phantom coordinate axes)
parser.add_argument("minvoxY",type=int)
parser.add_argument("minvoxZ",type=int)
parser.add_argument("roisizeCClus",type=int)    		# ROI length in pixels for the  microcalcification cluster
parser.add_argument("roisizeSpic",type=int)     		# ROI length in pixels for the spiculated mass
parser.add_argument("slice3dCClus",type=int)			# number of slices needed in microcalification cluster VOI
parser.add_argument("slice3dSpic",type=int)				# number of slices needed in spiculated mass VOI
parser.add_argument("CClus_SliceOffsetRem",type=int)	# number of slices to deduct from the center slice to obtain lower bound of VOI for microcalfication cluster
parser.add_argument("CClus_SliceOffsetAdd",type=int)	# number of slices to add to the center slice to obtain the upper bound of VOI for microcalfication cluster
parser.add_argument("Spic_SliceOffsetRem",type=int)		# number of slices to deduct from the center slice to obtain lower bound of VOI for spiculated mass
parser.add_argument("Spic_SliceOffsetAdd",type=int)		# number of slices to add to the center slice to obtain lower bound of VOI for spiculated mass
parser.add_argument("voxelsize",type=float)       		# phantom voxel size (in mm)
args = parser.parse_args()			

# change the current working directory
os.chdir('PATH TO DIRECTORY CONTAINING THIS CODE')

# output file for saving roi locations
floc = open("PATH TO THE FOLDER WHERE EXTRACTED ROI LOCATIONS ARE SAVED/roi_SP_"+str(args.RndSeed)+".loc",'w')

# read DBT image set as float 64 bit real (little endian order), raw file contains args.dim_z images
dt=np.dtype('<f8');	# little endian format for 32 bit float (real)
fid = open("PATH TO THE RECONSTRUCTED VOLUME/DBT_fatty_pcl_"+str(args.RndSeed)+"_recon.raw","r")	# file open

myimg = np.zeros(([args.dim_x,args.dim_y,args.dim_z]),dtype=dt);
for jj in range(0,(args.dim_z-1)):
	myoffset = jj*args.dim_x*args.dim_y*8;	# offset for reading 64 bit image
	fid.seek(myoffset,0)	# offset the file position from where to read
	tmpimg = np.fromfile(fid, dtype=dt, count=args.dim_x*args.dim_y)
	myimg[:,:,jj] = tmpimg.reshape(args.dim_x,args.dim_y,order='F');

# read lesion locations
lesarr = np.loadtxt("PATH TO LESION LOCATIONS/pcl_"+str(args.RndSeed)+"_crop.loc", delimiter=' ', dtype = (int), unpack=True);	# these are in voxel locations (in int)
lesloc = np.transpose(lesarr);	# because lesarr reads in col,row format. so transpose will give each row as (x,y,z,lesion type)
numloc = len(lesloc)	# number of rows

myvoxel = args.pixelsize	# [in mm] voxel pitch = pixel pitch = 0.12 mm

x1 = np.zeros(([1,numloc]),dtype=int);
y1 = np.zeros(([1,numloc]),dtype=int);
z1 = np.zeros(([1,numloc]),dtype=int);

for ii in range(0,numloc):

	
	lesX = lesloc[ii,1];	# In lesion insertion phantom is reshaped as [Z Y X] and lesions outputted as [Z X Y].  Thus X would be at [ii,1], and Y [ii,2] 
	lesY = lesloc[ii,2];

	lesZ = lesloc[ii,0];
	lesZ = lesZ - args.cr_z;	# REDUCE CR_Z NUMBER OF PLANES FROM lesZ

	

    ## Since we know the lesion location in phantom voxels, so calculate distance (in voxels) from the origin (vox_loc - 0) [given origin vox in phantom = 0,0,0]
    ## convert distance from vox to mm = (vox_loc-0)* voxelsize
    ## convert this distance to DBT pixels by dividing with appropriate pixelsize (x=y=DBT pixelsize; z=slice thickness).
    ## This would be DBT pixel location assuming the origin point for DBT is also 0,0,0.

	x1[0,ii] = (int)((lesX * args.voxelsize)/(args.pixelsize)) # in recon voxels
	y1[0,ii] = (int)((lesY * args.voxelsize)/(args.pixelsize))
	z1[0,ii] = (int)((lesZ * args.voxelsize)/(args.sl_thick))


    ## MC-GPU coordinate system flips the Y axis, so find mirror (against center of y axis) of the lesion location in Y
    
	centerpix = (int)(args.dim_x/2)
	mydiffpix = (y1[0,ii] - centerpix)
	if (mydiffpix < 0): # y1[0,ii] lies to the left of center of phantom
		y1[0,ii] = centerpix + abs(mydiffpix)
	elif(mydiffpix >= 0): # lies to the right
		y1[0,ii] = centerpix - abs(mydiffpix)


    ## interchange x and y
	temp2 = x1[0,ii]
	x1[0,ii] = y1[0,ii]
	y1[0,ii] = temp2
	

	print x1[0,ii], y1[0,ii], z1[0,ii]
	floc.write('%d %d %d %d\n' % (x1[0,ii],y1[0,ii],z1[0,ii],ii))

## extract roi
roipixel = np.zeros(([1,2]),dtype=int);		# number of pixels in the roi based on the lesion type {calc cluster (0), spiculated mass (1)}

roipixel[0,0] = args.roisizeCClus;
roipixel[0,1] = args.roisizeSpic;

# 2D ROI extraction
for ii in range(0,numloc):

	mypix = int((roipixel[0,lesloc[ii,3]] - 1)/2)	# roipixel[0,lesloc[ii,3]] is the number of pixels to be extracted for the roi; mypix - number of pixels in each direction with x1,y1 at ctr (excluding x1,y1).
	myroi = np.zeros(([roipixel[0,lesloc[ii,3]],roipixel[0,lesloc[ii,3]]]),dtype=dt)	# extracted roi image with the size based on the lesion type extracted

	for myx in range(x1[0,ii]-mypix,x1[0,ii]+mypix+1):	# we want x1+mypix to be included
		for myy in range(y1[0,ii]-mypix,y1[0,ii]+mypix+1):
			myroi[(myx+mypix-x1[0,ii]),(myy+mypix-y1[0,ii])] = myimg[myx,myy,z1[0,ii]]; 

	# SAVE ROIs except for spiculated mass
	if (lesloc[ii,3] == 1):
		myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED 2D ROIs FOR SPICULATED MASS/roiDBT_SP_"+str(ii)+"_"+str(args.RndSeed)+"_2D.raw")


roipix3d = np.zeros(([1,2]),dtype=int);	# number of slices for 3d volume for each lesion type

roipix3d[0,0] = args.slice3dCClus;	# calc cluster
roipix3d[0,1] = args.slice3dSpic;	# spiculated mass


# 3D VOI extraction
for ii in range(0,numloc):
	mypix = int((roipixel[0,lesloc[ii,3]] - 1)/2)	# roipixel[0,lesloc[ii,3]] is the number of pixels to be extracted for the roi

	# find center slice for 3d volume
	if((roipix3d[0,lesloc[ii,3]] % 2) == 0):	# even	
		my3Dpix = roipix3d[0,lesloc[ii,3]]/2
	else:						# odd
		my3dpix = (roipix3d[0,lesloc[ii,3]] - 1)/2

	my2Droi = np.zeros(([roipixel[0,lesloc[ii,3]],roipixel[0,lesloc[ii,3]]]),dtype=dt)				# 2d roi for cluster calc - obtained from avg 3d voi
	my3Droi = np.zeros(([roipix3d[0,lesloc[ii,3]],roipixel[0,lesloc[ii,3]],roipixel[0,lesloc[ii,3]]]),dtype=dt)	# 3D volume - extracted voi with the size based on the lesion type extracted



	if (lesloc[ii,3] == 0):		# calc cluster - slice range [center-5, center+2] both inclusive
		for myx in range(x1[0,ii]-mypix,x1[0,ii]+mypix+1):	# we want x1+mypix to be included
			for myy in range(y1[0,ii]-mypix,y1[0,ii]+mypix+1):
				zcounter = 0
				for myz in range(z1[0,ii]-args.CClus_SliceOffsetRem,z1[0,ii]+args.CClus_SliceOffsetAdd+1):        
                                        my3Droi[zcounter,(myx+mypix-x1[0,ii]),(myy+mypix-y1[0,ii])] = myimg[myx,myy,myz]
                                        zcounter = zcounter + 1
		
		my2Droi = np.mean(my3Droi,axis=0)	# take mean over Z axis
		mip2Droi = np.max(my3Droi,axis=0)       # maximum intensity projection along z axis

		# SAVE 2D roi's - mean + MIP
		my2Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED 2D MEAN ROIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SP_"+str(ii)+"_"+str(args.RndSeed)+"_mean_2D.raw")
		mip2Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED 2D MIP ROIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SP_"+str(ii)+"_"+str(args.RndSeed)+"_mip_2D.raw")

	elif(lesloc[ii,3] == 1):	# spiculated masses, slice range (center +- 4) both inclusive
		for myx in range(x1[0,ii]-mypix,x1[0,ii]+mypix+1):	# we want x1+mypix to be included
			for myy in range(y1[0,ii]-mypix,y1[0,ii]+mypix+1):
				zcounter = 0
				for myz in range(z1[0,ii]-args.Spic_SliceOffsetRem,z1[0,ii]+args.Spic_SliceOffsetAdd+1):       				
                                        my3Droi[zcounter,(myx+mypix-x1[0,ii]),(myy+mypix-y1[0,ii])] = myimg[myx,myy,myz]
                                        zcounter = zcounter + 1


	# SAVE VOIs
	if (lesloc[ii,3] == 0):
		my3Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED VOIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SP_"+str(ii)+"_"+str(args.RndSeed)+"_3D.raw")
	elif (lesloc[ii,3] == 1):
		my3Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED VOIs FOR SPICULATED MASS/roiDBT_SP_"+str(ii)+"_"+str(args.RndSeed)+"_3D.raw")
