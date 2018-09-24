#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#						ROI EXTRACTION (FFDM - SIGNAL ABSENT)
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
import os.path
import gzip


# Change working directory
os.chdir('PATH TO DIRECTORY CONTAINING THIS CODE')


# Parsing input arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("RndSeed",type=int)				# phantom random seed - from the breast phantom filename
parser.add_argument("fs_x",type=float)  			# focal spot locations X/Y/Z (in mm) (from x-ray detector coordinate system)
parser.add_argument("fs_y",type=float)
parser.add_argument("fs_z",type=float)
parser.add_argument("cr_z",type=int)				# number of Z planes to be cropped in MC-GPU output images
parser.add_argument("voxelsize",type=float)			# phantom voxel size (in mm)
parser.add_argument("pixelsize",type=float)			# detector pixel size (in mm)
parser.add_argument("orgX",type=float)				# origin X/Y/Z for cropped phantom (in mm)
parser.add_argument("orgY",type=float)
parser.add_argument("orgZ",type=float)
parser.add_argument("minvoxX",type=int)				# minimum voxel number in X,Y,Z dimensions (lowest number on cropped phantom coordinate axes)
parser.add_argument("minvoxY",type=int)
parser.add_argument("minvoxZ",type=int)
parser.add_argument("totvoxX",type=int)				# total number of voxels in X/Y/Z dimensions for cropped phantom
parser.add_argument("totvoxY",type=int)
parser.add_argument("totvoxZ",type=int)
parser.add_argument("cclus_voxvol",type=int)		# length (in voxels) for microcalcification cluster (considering lesion volume to be a cube)
parser.add_argument("spic_voxvol",type=int)			# length (in voxels) for spiculated mass (considering lesion volume to be a cube)
parser.add_argument("roisizeCClus",type=int)   	 	# ROI length in pixels for the  microcalcification cluster
parser.add_argument("roisizeSpic",type=int)    		# ROI length in pixels for the spiculated mass
parser.add_argument("totnumROI",type=int)       	# total number of ROIs per phantom to be extracted (divided equally among diff lesion types)
parser.add_argument("lendetX_pix",type=int)     	# length of the detector in X/Y (in pixels)
parser.add_argument("lendetY_pix",type=int)     
parser.add_argument("det_z_loc",type=float)    		# Z location for detector (in mm) - based on the air gap between bottom of phantom and the detector plane
parser.add_argument("limit_chestX_vox",type=float)  # min X bound in phantom (in voxel number) to avoid overlap with the chest wall
args = parser.parse_args()

# Use system calls to unzip phantom to local temp for faster access (useful if phantom stored on a remote system)
fl1 = "PATH TO THE PHANTOM/pc_"
fl2 = str(args.RndSeed)
fl3 = "_crop.raw.gz"
fl4 = "_crop.raw"
fl5 = "pc_"

flname = fl1+fl2+fl3
flname2 = fl5+fl2+fl4

runcmd1 = "gunzip -c "+flname+" > /tmp/"+flname2 		# gunzip -c file.raw.gz > /tmp/file.raw
os.system(runcmd1)


# read binary file (compressed phantom) as unsigned int 8 bits and reshape in a matrix of size 1051x1657x713
phantom = np.fromfile("/tmp/pc_"+str(args.RndSeed)+"_crop.raw", dtype=np.uint8).reshape(args.totvoxZ,args.totvoxY,args.totvoxX);

px,py,pz = phantom.shape	# obtain phantom size

# read mammography projection image as florat 32 bit (little endian order), raw file contains 5 images (each being 850x850 pixels), we only need the first one
dt=np.dtype('<f4');     # little endian format for 32 bit float (real)

myimg=np.zeros(([args.lendetX_pix,args.lendetY_pix]),dtype=dt)  # define final flatfield corrected image

fid = open("PATH TO MC-GPU OUTPUT RAW IMAGE/mcgpu_image_pc_"+str(args.RndSeed)+"_crop.raw.gz-fatty_0000.raw","r")
tmpimgORG = np.fromfile(fid, dtype=dt, count=args.lendetX_pix*args.lendetY_pix)
myimgORG = tmpimgORG.reshape(args.lendetX_pix,args.lendetY_pix);      # reshape for PROJECTION IMAGE X AND Y

## Read flatfield image
fidFF = open("PATH TO MC-GPU FLATFIELD RAW IMAGE/mcgpu_results/mcgpu_image-fatty-flatfield_0000.raw","r")
tmpimgFF = np.fromfile(fidFF, dtype=dt, count=args.lendetX_pix*args.lendetY_pix)
myimgFF = tmpimgFF.reshape(args.lendetX_pix,args.lendetY_pix)

## Generate flatfield corrected mammo image for breast phantom
for myxx in range(0,args.lendetX_pix):
	for myyy in range(0,args.lendetY_pix):
		                myimg[myxx,myyy] = (float)(myimgFF[myxx,myyy]/myimgORG[myxx,myyy])
				

# create a mask for projection image for checking ROI overlap
maskimg = np.zeros((args.lendetX_pix,args.lendetY_pix),dtype=dt);

# store final ROI centers (projection)
x1 = np.zeros(([1,args.totnumROI]),dtype=int);
y1 = np.zeros(([1,args.totnumROI]),dtype=int);

# (phantom locations used to extract)				
phanX = np.zeros(([1,args.totnumROI]),dtype=int);
phanY = np.zeros(([1,args.totnumROI]),dtype=int);
phanZ = np.zeros(([1,args.totnumROI]),dtype=int);				


# read locations file (potential locations for extracting appropriate ROIs)
x, y, z = np.loadtxt("PATH TO FILE CONTAINING POTENTIAL ROI LOCATIONS/pc_"+str(args.RndSeed)+".loc", delimiter=',', usecols=(0, 1, 2), unpack=True);	# these are positions in mm, which needs to be changed to voxel number

numloc = len(x) # number of possible lesion locations

# create lesion locations array
lesloc = np.ndarray((numloc,3),dtype = float) 

# create lesion voxels array
lesvox = np.ndarray((numloc,4),dtype = int)  # array to store which voxels has been used or tried so far - x,y,z,flag - flag 0(available)/1(used or tried already)

lesvox[:,0:4] = 0


# copy x,y,z in respective elements of matrix {col 0-2 are float}, the fourth column is the flag {integer} (0 for now)
lesloc[:,0] = x # select all rows from column 0
lesloc[:,1] = y
lesloc[:,2] = z

# convert lesion location (in mm) to lesion voxel space (ASSUMING HARD CODED VALUES [taken from pc_*.VTI file] OF THE 0TH VOXEL FOR PRE-PILOT STUDY Origin="-9.44 -66.236709213 -33.089583664" Piece Extent="-20 1030 -276 1380 -60 652" Voxel "0.12 mm")
# this is based on 1051x1657x713, but python reads with x and z interchanged, so we interchange the location in x and z as well
# formula = Voxel number = {[curr pos - 0,0,0 pos (in mm)]/voxel size} - min voxel

myvoxel = args.voxelsize
loczeroX = args.orgX
loczeroY = args.orgY
loczeroZ = args.orgZ
voxminX = args.minvoxX
voxminY = args.minvoxY
voxminZ = args.minvoxZ


for myii in range(0,numloc):
	lesvox[myii,0] = np.floor(((lesloc[myii,0] - loczeroX)/myvoxel) - voxminX);			# x
	lesvox[myii,1] = np.floor(((lesloc[myii,1] - loczeroY)/myvoxel) - voxminY);               	# y
	lesvox[myii,2] = np.floor(((lesloc[myii,2] - loczeroZ)/myvoxel) - voxminZ);               	# z

## Calculate number of voxels in lesion model from number of pixels in ROI and pixelsize, voxelsize

numVox_CClus = int((args.roisizeCClus * args.pixelsize)/args.voxelsize);
numVox_Spic  = int((args.roisizeSpic * args.pixelsize)/args.voxelsize);


myoddflag = np.zeros(([1,2]),dtype=int);		# raise this if lesion is odd sized - each element correspond to each lesion type - calcCluster(0), spiculated (1)

if np.mod(numVox_CClus,2) == 0:      # even
        myoddflag[0,0] = 0;
else:                		    # odd
	myoddflag[0,0] = 1;


if np.mod(numVox_Spic,2) == 0:
	myoddflag[0,1] = 0;
else:
	myoddflag[0,1] = 1;


##### if loop to insert lesions
insertedLesions = 0	# variable to indicate number of inserted lesions so far
checkFail = 0	# flag to indicate if lesion passed the checks or not (0 pass or not tried yet, 1 fail)
myctr = 0 	# a counter to prevent while loop going in infinity

numLesToInsert = args.totnumROI  	# number of lesions to be inserted in total - in this case equal to total number of ROIs

if np.mod(args.totnumROI,2) == 0: # if total number of ROIs is even 
	myhalfnum = args.totnumROI/2
else: 
	myhalfnum = (args.totnumROI-1)/2
		
myLesionType = np.zeros(([args.totnumROI]),dtype=int);  

# initialize myLesionType - divide the total number of ROIs between how many to be extracted for sizes roisizeCClus and roisizeSpic
for ii in range (0,myhalfnum):
	myLesionType[ii] = 0	# indicating ROI size 'roisizeCClus'
	
for ii in range (myhalfnum,args.totnumROI):
	myLesionType[ii] = 1	# indicating ROI size 'roisizeSpic'
	
ctr1 = 0;	# counter for myLesionType

eliminated_locs = 0 # counter to check how many locations eliminated
maxattempt = 100;

# Initialize random seed generator (using /dev/urandom), goal to have reproducible runs if needed
# for reproducible runs, set (hardcode) urand_seed[0] which will then initialize np.random.seed sequence to generate same random sequence every time
finrnd=open("/dev/urandom","r")
urand_seed=np.fromfile(finrnd,dtype=int,count=1)
np.random.seed(abs(urand_seed[0]))

print urand_seed[0]

# variables to check percentage of overlap between same roi/lesion type
tot_vox_overlap = 0;
percent_overlap = 0.0;

while (insertedLesions < numLesToInsert) and (eliminated_locs < (numloc-1)):

	# pick a location for lesion insertion randomly
	myrnd = np.random.random_integers(numloc-1)
	currX = lesvox[myrnd,2];	# Z dim (713)
	currY = lesvox[myrnd,1];	# Y dim (1657)
	currZ = lesvox[myrnd,0];	# X dim (1051), interchanging positions of pos of x, y and z due to reading phantom in this way (lesion locations are based on 1051,1657,713, while we are reading phantom as 713,1657,1051

	# reset checkFail
	checkFail = 0;

	## insert lesion
	# choose type of lesion and accordingly change variables
        if (myLesionType[ctr1]==0) and (myoddflag[0,0] == 0):
                hfles = int(numVox_CClus*0.5);
                les = numVox_CClus;
	elif (myLesionType[ctr1]==1) and (myoddflag[0,1] == 0):
		hfles = int(numVox_Spic*0.5);
		les = numVox_Spic;
        elif (myLesionType[ctr1]==0) and (myoddflag[0,0] == 1):		# ODD SIZED LESIONS - add extra plane of zeros in all three dimensions
                hfles = int((numVox_CClus+1)*0.5);
                les = numVox_CClus+1;
        elif (myLesionType[ctr1]==1) and (myoddflag[0,1] == 1):
                hfles = int((numVox_Spic+1)*0.5);
                les = numVox_Spic+1;


	# is this lesion voxel available for use (flag == 0)?
	if lesvox[myrnd,3] == 0:

	        # location number available for use
	        print("Eliminated %d/%d \n") % (eliminated_locs,numloc)						


                ###### CHECK IF LESION IS TOO CLOSE TO THE CHEST WALL - HARDCODED X USING AMIDE VIEW OF RAW PHANTOM (diff (mm) from min x to edge of upper lip of paddle is converted to voxels)
                ## (assuming all breats in one category would have chest wall around same locations, thus same hardcoded X could be used)
                check_chestwall_X = args.limit_chestX_vox;   # in voxels

		if (currZ < check_chestwall_X):       # REMEMBER X AND Z HAVE BEEN INTERCHANGED, SO WE ARE ACTUALLY CHECKING IN X DIMENSION USING currZ
			lesvox[myrnd,3] = 1;          # lesion discarded due to being close to chest wall
			checkFail = 1;          # lesion failed the checks
			myctr = myctr + 1;
			print("check failed - close to chestwall \n")

		####### make changes for odd/even voxel lesion - FOR ODD SIZED LESIONS - CONSIDERING HFLES TO BE SIZE+1/2 -> SO THE BELOW SECTION WILL ADD ONE EXTRA VOXEL

		# check whether the lesion volume lies within the phantom boundaries	
		if (lesvox[myrnd,3] != 1) and ((int(currX-hfles) >= 0) and (int(currX+hfles-1) < px) and (int(currZ-hfles) >= 0) and (int(currZ+hfles-1) < pz) and (int(currY-hfles) >= 0) and (int(currY+hfles-1) < py)):

			## Perform required checks to determine if lesion can be inserted at current location

			# Create a temporary phantom 3d array - equal to size of lesion - contains phantom voxels corresponding to where lesion to be inserted
			tmp_phan=np.zeros((les,les,les),dtype=np.uint8) # assign 3D array (with lesion dimensions) of zeros with type uint8

			# Copy particular volume of phantom (where lesion to be placed) from phantom to tmp_phan
			for myx in range(int(currX-hfles),int(currX+hfles)):      # give the range [low,high) - high value is not inclusive, in x,y,z (3d volume) where lesion to be inserted
	                	for myy in range(int(currY-hfles),int(currY+hfles)):
	       	                	for myz in range(int(currZ-hfles),int(currZ+hfles)):
						tmp_phan[abs(currX-hfles-myx),abs(currY-hfles-myy),abs(currZ-hfles-myz)] = phantom[myx,myy,myz];	            


			# Preliminary check - check the eight corners of the lesion cube (in tmp_phan) to see if they lie in any of the foll (air,skin,nipple,muscle) - if yes, right away reject it because then the rest of the points will also not comply
			# if this is passed, then go for detailed check below
			corners=[0,les-1]
			for myxC in corners:
				for myyC in corners:
					for myzC in corners:
						if (tmp_phan[myxC,myyC,myzC] == 0) or (tmp_phan[myxC,myyC,myzC] == 2) or (tmp_phan[myxC,myyC,myzC] == 33) or (tmp_phan[myxC,myyC,myzC] == 40): # reject this lesion location
                                                        lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable for future due to not passing the checks
							checkFail = 1;          # lesion failed the checks
							myctr = myctr + 1;      # counter to prevent while loop going in infinity
							print("check failed - preliminary check %d,%d,%d \n") % (myxC,myyC,myzC)
                                                        break;


                        # Detailed check - Check if phantom where lesion to be inserted contains air(0), skin(2), nipple (33), or muscle (40) - if yes, do not insert lesion, choose another random location
                        # !!!! DEBUG !!!! CHECK ONLY SURFACES OF THE CUBE

                        if lesvox[myrnd,3] != 1:
                                # check for x = 0, les-1
                                for myxT in (0,les-1): # REMOVED 'RANGE()'
                                        for myyT in range(0,les):
                                                for myzT in range(0,les):
                                                        if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40): # re$
                                                                lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable for future due to not passing the checks
                                                                break;

                                # check for y = 0, les-1
                                for myyT in (0,les-1): # REMOVED 'RANGE()'
                                        for myxT in range(0,les):
                                                for myzT in range(0,les):
                                                        if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40):
                                                                lesvox[myrnd,3] = 1;
                                                                break;


                                # check for z = 0, les-1
                                for myzT in (0,les-1): # REMOVED 'RANGE()'
                                        for myxT in range(0,les):
                                                for myyT in range(0,les):
                                                        if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40):
                                                                lesvox[myrnd,3] = 1;
                                                                break;

				if lesvox[myrnd,3] == 1:
					checkFail = 1;          # lesion failed the checks
					myctr = myctr + 1;      # counter to prevent while loop going in infinity
					print("check failed - detailed check \n")



			# insert lesion at current location
			if checkFail == 0:	# Potential ROI available
				
				## Calculate the projection of these lesions
				##### MADE MODIFICATION OCT 19, 2016 - WHEN COMBINING LESION INSERTION AND ROI EXTRACTION CODES, THE X, Z SIMENSIONS WERE INTERCHANGED BETWEEN THE TWO CODES, SO
                ##### MODIFIED TO REPRESENT THIS BELOW
				
				lesZ = currX;	# these lesion voxel numbers start from 0, not negative voxels (that is how the lesion insertion code outputs them which is fine)
				lesY = currY;
				lesX = currZ;
				lesZ = lesZ - args.cr_z;        #!!!! DETECTOR MOVED CR_Z VOXELS UP, SO REDUCE THAT MANY Z PLANES (note that lesX contains the z voxel number for lesion)

			    ## LESION LOCATIONS FROM .LOC FILE HAS VOXEL NUMBERS FROM ORIGINAL PHANTOM SO CALCULATE LOCATIONS IN MM USING ORIGINAL OFFSET
			    ## ONCE IN MM, THEN IT DOES NOT MATTER WHETHER THE PHANTOM WAS CROPPED OR NOT BECAUSE THE PHANTOM COORDINATE SYSTEM REMAINS SAME FOR BOTH

			    ## CONVERT LESION VOXELS TO MM (CROPPING - PHANTOM COOR SYSTEM)
			    ## (lesion_vox - origin_vox)*voxelsize + origin_mm
			    lesX_mm = (lesX-0)*args.voxelsize + args.orgX
			    lesY_mm = (lesY-0)*args.voxelsize + args.orgY
			    lesZ_mm = (lesZ-0)*args.voxelsize + args.orgZ

			    ## Focal spot (mm) and projection location (mm) both are in MCGPU coor system.
			    ## Translate lesion (mm) from phantom to MCGPU coor system
			    ## Phantom origin [0,0,0] vox or (X1,Y1,Z1) mm lies at lower left back corner of phantom.
			    ## Detector (mcgpu) origin lies in the same location with (0,0) mm. Remember this is not the top left corner (x=0 at chest wall) of detector.
			    ## Detector origin (mm) coincides with the voxel origin (0,0,0 vox or X1,Y1,Z1 mm) of cropped phantom. so a point (a,b,c) mm in phantom is translated as (a-X1, b-Y1, c-Z1) mm in detector space.
				## The x,y,z axis directions are same for both phantom and detector, thus no flipping of xyz coordinates between the two systems.
				crop_phan_lenY_mm = args.totvoxY * args.voxelsize # cropped phantom length in Y dimension (mm)
				det_lenY_mm = args.lendetY_pix * args.pixelsize # detector length in Y dimension (mm)                                         

				# translating from phantom to detector (mm)
				lesX_det_mm = lesX_mm - args.orgX
				lesY_det_mm = lesY_mm - args.orgY
			        lesZ_det_mm = lesZ_mm - args.orgZ  # this indicates lesZ in detector space in mm

				# calculate alpha
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

			        ## Calculate pixel number in dim Y for origin
				## Calculate projection coordinates in pixels = detector origin pixels + distance in pixels between detector origin and projection
			        proj_pix_X = dist_pix_proj_orgX + det_orgX_pix
			        proj_pix_Y_tmp = dist_pix_proj_orgY + det_orgY_pix
				proj_pix_Y = args.lendetY_pix - proj_pix_Y_tmp # we figured out by looking at the voxels and pixels that Y dimension was flipped, so final Y_pix = len_det_Y - proj_Y_pix_tmp

				
				tmpX = int(proj_pix_X);		# projection in X
				tmpY = int(proj_pix_Y);		# in Y


				# calculate number of pixels in half ROI
				if (myLesionType[ctr1] == 0) and (np.mod(args.roisizeCClus,2) == 0):       # Calc clus, even number of pixels in ROI
				        mypix = int(args.roisizeCClus*0.5);
				elif (myLesionType[ctr1] == 0) and (np.mod(args.roisizeCClus,2) != 0):
				        mypix =int((args.roisizeCClus-1)*0.5);

				if (myLesionType[ctr1] == 1) and (np.mod(args.roisizeSpic,2) == 0): # Spiculated mass
				        mypix = int(args.roisizeSpic*0.5);
				elif (myLesionType[ctr1] == 1) and (np.mod(args.roisizeSpic,2) != 0):
				        mypix = int((args.roisizeSpic-1)*0.5);


				# determine if ROI is available (non-overlapping)	
				if (ctr1 == 0) or (ctr1 == 6):         # first lesion of each type - simply extract this ROI
				        x1[0,ctr1] = tmpX;
				        y1[0,ctr1] = tmpY;

                                        phanX[0,ctr1] = lesX;
                                        phanY[0,ctr1] = lesY;
                                        phanZ[0,ctr1] = lesZ + args.cr_z    # not including removing cr_z voxels
                                        
					# mark this region as '1's/'2's in the mask based on ROI type (calc clus - 1, spic - 2)
                                        if (myLesionType[ctr1] == 0):
                                                for myx in range (tmpX-mypix, tmpX+mypix+1):      # we want x1+mypix to be included
                                                        for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                maskimg[myx,myy] = 1;             # CALC CLUSTER
                                        elif (myLesionType[ctr1] == 1):
                                                for myx in range (tmpX-mypix, tmpX+mypix+1):                                     
                                                        for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                maskimg[myx,myy] = 2;             # SPICULATED MASS


	                                insertedLesions = insertedLesions + 1;  # mark lesion as inserted
	                                lesvox[myrnd,3] = 1                     # mark as unavailable
					myctr = 0;                              # reset myctr

					print("\n*************** check PASSED xy %d,%d ************** Inserted so far %d ***********\n\n") % (x1[0,ctr1],y1[0,ctr1], insertedLesions)

					ctr1 = ctr1 + 1;	# move to extract next ROI (next lesion type)

					# Extract ROI and save as raw image
					if (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) == 0):		# extracted roi image with the size based on the lesion type extracted
						myroi = np.zeros( ((args.roisizeCClus+1),(args.roisizeCClus+1)),dtype=dt ) # if ROI is even pixel, then add 1 pixel. we want tmpX,tmpY to lie in center
					elif (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) != 0):
						myroi = np.zeros( (args.roisizeCClus,args.roisizeCClus),dtype=dt )
						
					elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) == 0):
                                                myroi = np.zeros( ((args.roisizeSpic+1),(args.roisizeSpic+1)),dtype=dt )
					elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) != 0):
						myroi = np.zeros( (args.roisizeSpic,args.roisizeSpic),dtype=dt )

					#print myroi.shape	

					for myx in range(tmpX-mypix,tmpX+mypix+1):      # we want x1+mypix to be included
						for myy in range(tmpY-mypix,tmpY+mypix+1):
							myroi[abs(tmpX-mypix-myx),abs(tmpY-mypix-myy)] = myimg[myx,myy];

					# save as unzipped file
                                        if (myLesionType[ctr1-1] == 0):
                                                myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR MICROCALCIFICATION CLUSTER/roiMAMMO_SA_"+str(ctr1-1)+"_"+str(args.RndSeed)+".raw")
                                        elif (myLesionType[ctr1-1] == 1):
                                                myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR SPICULATED MASS/roiMAMMO_SA_"+str(ctr1-1)+"_"+str(args.RndSeed)+".raw")

					
									
					
				elif (ctr1 > 0):
				        # check whether current ROI projection would overlap any existing ROIs of same lesion/roi type (up to 10% overlap is allowed with same roi type)
                                        tot_vox_overlap = 0;# reset
                                        
                                        if (myLesionType[ctr1] == 0):
                                                for myx in range (tmpX-mypix, tmpX+mypix+1):	# check - for current ROI region, HOW MANY 1'S OR 2'S ARE THERE BASED ON ROI TYPE
                                                        for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                if(maskimg[myx,myy] == 1):      # CALC CLUSTER
                                                                        tot_vox_overlap = tot_vox_overlap + 1;

						
                                        elif (myLesionType[ctr1] == 1):
                                                for myx in range (tmpX-mypix, tmpX+mypix+1):    # SPIC MASS                          
                                                        for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                if(maskimg[myx,myy] == 2):
                                                                        tot_vox_overlap = tot_vox_overlap + 1;


                                        if (myLesionType[ctr1] == 0):
                                                percent_overlap = ( (float)(tot_vox_overlap) / (args.roisizeCClus * args.roisizeCClus) )                        
                                        elif (myLesionType[ctr1] == 1):
                                                percent_overlap = ( (float)(tot_vox_overlap) / (args.roisizeSpic * args.roisizeSpic) )
                                                
					

                                        if (percent_overlap > 0.1):      # ROI overlap greater than 10% between same ROI type
						lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable for future due to not passing the checks
						checkFail = 1;          # lesion failed the checks
						myctr = myctr + 1;      # counter to prevent while loop going in infinity - break while loop if consecutive 50 locations have been found inappropriate
						print("check failed - ROI mask overlap \n")
						


					if (checkFail == 0):		# ROI AVAILABLE
						x1[0,ctr1] = tmpX;      # projection voxels
	                                        y1[0,ctr1] = tmpY;


                                                phanX[0,ctr1] = lesX;   # phantom voxel location used
                                                phanY[0,ctr1] = lesY;
                                                phanZ[0,ctr1] = lesZ + args.cr_z    # not including removing cr_z voxels 

                                                # mark this region as '1's/'2's in the mask based on ROI type (calc clus - 1, spic - 2)                                                                    
                                                if (myLesionType[ctr1] == 0):
                                                        for myx in range (tmpX-mypix, tmpX+mypix+1):      # we want x1+mypix to be included                                                               
                                                                for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                        maskimg[myx,myy] = 1;             # CALC CLUSTER                                                                                  
                                                elif (myLesionType[ctr1] == 1):
                                                        for myx in range (tmpX-mypix, tmpX+mypix+1):
                                                                for myy in range(tmpY-mypix, tmpY+mypix+1):
                                                                        maskimg[myx,myy] = 2;             # SPICULATED MASS 


						insertedLesions = insertedLesions + 1;  # mark lesion as inserted
						lesvox[myrnd,3] = 1                     # mark as unavailable
						myctr = 0;                              # reset myctr

						print("\n*************** check PASSED xy %d,%d ************** Inserted so far %d ***********\n\n") % (x1[0,ctr1],y1[0,ctr1], insertedLesions)
						
						ctr1 = ctr1 + 1;        # move to extract next ROI (next lesion type) 
								
						# Extract ROI and save as raw image based on lesion type
						if (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) == 0):
	                                                myroi = np.zeros( ((args.roisizeCClus+1),(args.roisizeCClus+1)),dtype=dt ) # ROI is even pixel, add 1 pixel. we want tmpX,tmpY (to lie in center)
					    	elif (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) != 0):
						        myroi = np.zeros( (args.roisizeCClus,args.roisizeCClus),dtype=dt )

						elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) == 0):
						         myroi = np.zeros( ((args.roisizeSpic+1),(args.roisizeSpic+1)),dtype=dt )
						elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) != 0):
						         myroi = np.zeros( (args.roisizeSpic,args.roisizeSpic),dtype=dt )
																												
							

						for myx in range(tmpX-mypix,tmpX+mypix+1):      # we want x1+mypix to be included
							for myy in range(tmpY-mypix,tmpY+mypix+1):
								myroi[abs(tmpX-mypix-myx),abs(tmpY-mypix-myy)] = myimg[myx,myy];

						# save as unzipped file

                                                if (myLesionType[ctr1-1] == 0):
                                                        myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR MICROCALCIFICATION CLUSTER/roiMAMMO_SA_"+str(ctr1-1)+"_"+str(args.RndSeed)+".raw")
                                                elif (myLesionType[ctr1-1] == 1):        
                                                        myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR SPICULATED MASS/roiMAMMO_SA_"+str(ctr1-1)+"_"+str(args.RndSeed)+".raw")

		elif (lesvox[myrnd,3] != 1):	# OUTSIDE PHANTOM BOUNDARY - MARK UNAVAILABLE
			lesvox[myrnd,3] = 1;          # lesion discarded due to being close to chest wall
			checkFail = 1;          # lesion failed the checks
			print("check failed - OUTSIDE PHANTOM BOUNDARY \n")


		#####  MARK LOCATION UNAVAILABLE - EITHER INSERTED OR CHECK FAILED
		if (checkFail == 1) or (lesvox[myrnd,3] == 1):		# if either lesion was inserted or any of the checks failed - mark location unavailable
			eliminated_locs = eliminated_locs + 1;



## print message if ran out of potential locations
if(eliminated_locs >= (numloc-1)):
	print("\n\n\t WHILE LOOP BROKE DUE TO RUNNING OUT OF POTENTIAL LESION LOCATIONS!!!!...............................\n\n")

					
# Save locations - phantom locations (vox) [Z dim does not include removing cr_z voxels], projection locations (pix)
floc = open("PATH TO THE FOLDER WHERE EXTRACTED ROI LOCATIONS ARE SAVED/roi_SA_"+str(args.RndSeed)+".loc",'w')
for kl in range(0,numLesToInsert):
        floc.write('%d %d %d %d %d\n' % (phanX[0,kl], phanY[0,kl], phanZ[0,kl], x1[0,kl], y1[0,kl]))

# Delete pc_*.raw file under /tmp
runcmd3 = "/tmp/"+flname2
#print runcmd3
os.remove(runcmd3)
