#*********************************************************************************************************************************
#*********************************************************************************************************************************
#
#						VICTRE (VIRTUAL IMAGING CLINICAL TRIAL FOR REGULATORY EVALUATION)
#						ROI EXTRACTION (DBT - SIGNAL ABSENT)
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

# change the current working directory
os.chdir('PATH TO DIRECTORY CONTAINING THIS CODE')

# Parsing input arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("RndSeed",type=int)				# phantom random seed - from the breast phantom filename
parser.add_argument("fs_x",type=float)	        	# focal spot locations X/Y/Z (in mm) (from x-ray detector coordinate system)
parser.add_argument("fs_y",type=float)
parser.add_argument("fs_z",type=float)
parser.add_argument("cr_z",type=int)	        	# number of Z planes cropped in projecttion FFDM images - recon is constructed from FFDM images
parser.add_argument("dim_x",type=int)	        	# number of pixels in X/Y/Z dimensions in the reconstructed raw image
parser.add_argument("dim_y",type=int)
parser.add_argument("dim_z",type=int)	        	# total number of slices
parser.add_argument("recon_voxelsize",type=float)	# reconstructed volume voxel size (in mm)
parser.add_argument("sl_thick",type=float)			# slice thickness of reconstructed volume (in mm)
parser.add_argument("orgX",type=float)          	# origin X/Y/Z for cropped phantom (mm)
parser.add_argument("orgY",type=float)
parser.add_argument("orgZ",type=float)
parser.add_argument("minvoxX",type=int)        		# minimum voxel number in X,Y,Z dimensions (lowest number on cropped phantom coordinate axes)
parser.add_argument("minvoxY",type=int)
parser.add_argument("minvoxZ",type=int)
parser.add_argument("roisizeCClus",type=int)    	# ROI length in pixels for the  microcalcification cluster
parser.add_argument("roisizeSpic",type=int)     	# ROI length in pixels for the spiculated mass
parser.add_argument("slice3dCClus",type=int)		# number of slices needed in microcalification cluster VOI
parser.add_argument("slice3dSpic",type=int)			# number of slices needed in spiculated mass VOI
parser.add_argument("CClus_SliceOffsetRem",type=int)# number of slices to deduct from the center slice to obtain lower bound of VOI for microcalfication cluster
parser.add_argument("CClus_SliceOffsetAdd",type=int)# number of slices to add to the center slice to obtain the upper bound of VOI for microcalfication cluster
parser.add_argument("Spic_SliceOffsetRem",type=int)	# number of slices to deduct from the center slice to obtain lower bound of VOI for spiculated mass
parser.add_argument("Spic_SliceOffsetAdd",type=int) # number of slices to add to the center slice to obtain lower bound of VOI for spiculated mass
parser.add_argument("totnumROI",type=int)       	# total number of ROIs per phantom to be extracted (divided equally among diff lesion types)
parser.add_argument("totvoxX",type=int)         	# total number of phantom voxels in X/Y/Z
parser.add_argument("totvoxY",type=int)
parser.add_argument("totvoxZ",type=int)
parser.add_argument("voxelsize",type=float)			# phantom voxel size (in mm)
parser.add_argument("limit_chestX_vox",type=int)	# min X bound in phantom (in voxel number) to avoid overlap with the chest wall
args = parser.parse_args()			


# Use system calls to unzip phantom to local temp for faster access (useful if phantom stored on a remote system)
fl1 = "PATH TO THE PHANTOM/pc_"
fl2 = str(args.RndSeed)
fl3 = "_crop.raw.gz"
fl4 = "_crop.raw"
fl5 = "pc_"

flname = fl1+fl2+fl3
flname2 = fl5+fl2+fl4

runcmd1 = "gunzip -c "+flname+" > /tmp/"+flname2 # gunzip -c file.raw.gz > /tmp/file.raw
os.system(runcmd1)

# output file for saving roi locations
floc = open("PATH TO THE FOLDER WHERE EXTRACTED ROI LOCATIONS ARE SAVED/roi_SA_"+str(args.RndSeed)+".loc",'w')

# read binary file (compressed phantom) as unsigned int 8 bits and reshape in a matrix of size 1051x1657x713
phantom = np.fromfile("/tmp/pc_"+str(args.RndSeed)+"_crop.raw", dtype=np.uint8).reshape(args.totvoxZ,args.totvoxY,args.totvoxX);
px,py,pz = phantom.shape	# obtain phantom size

# read DBT image set as float 64 bit real (little endian order), raw file contains args.dim_z images
dt=np.dtype('<f8');	# little endian format for 32 bit float (real)
fid = open("PATH TO MC-GPU OUTPUT RAW IMAGE/DBT_fatty_pc_"+str(args.RndSeed)+"_recon.raw","r")	# file open

myimg = np.zeros(([args.dim_x,args.dim_y,args.dim_z]),dtype=dt);
for jj in range(0,(args.dim_z-1)):
	myoffset = jj*args.dim_x*args.dim_y*8; # MULTIPLY BY 8 FOR 64 BUT REAL IMAGE	
	fid.seek(myoffset,0)	# offset the file position from where to read
	tmpimg = np.fromfile(fid, dtype=dt, count=args.dim_x*args.dim_y)
	myimg[:,:,jj] = tmpimg.reshape(args.dim_x,args.dim_y,order='F');

# read locations file (potential locations for extracting appropriate ROIs)
x, y, z = np.loadtxt("PATH TO FILE CONTAINING POTENTIAL ROI LOCATIONS/pc_"+str(args.RndSeed)+".loc", delimiter=',', usecols=(0, 1, 2), unpack=True); # (in mm)
numloc = len(x)

# create lesion locations array
lesloc = np.ndarray((args.totnumROI,4),dtype = int)# randomly chosen potential lesion (ROI) locations in voxels + lesion type
temploc = np.ndarray((numloc,3),dtype = float) # all locations from potential lesion locations file

# create lesion voxels array
lesvox = np.ndarray((numloc,4),dtype = int)  # array to store which voxels has been used or tried so far - x,y,z,flag - flag 0(available)/1(used or tried already)
lesvox[:,0:4] = 0

# copy x,y,z in respective elements of matrix {col 0-2 are float}, the fourth column is the flag {integer} (0 for now)
temploc[:,0] = x # select all rows from column 0
temploc[:,1] = y
temploc[:,2] = z

loczeroX = args.orgX
loczeroY = args.orgY
loczeroZ = args.orgZ
voxminX = args.minvoxX
voxminY = args.minvoxY
voxminZ = args.minvoxZ

for myii in range(0,numloc):
        lesvox[myii,0] = np.floor(((temploc[myii,0] - loczeroX)/args.voxelsize) - voxminX);                     # x (in vox - phantom space)
        lesvox[myii,1] = np.floor(((temploc[myii,1] - loczeroY)/args.voxelsize) - voxminY);                     # y
        lesvox[myii,2] = np.floor(((temploc[myii,2] - loczeroZ)/args.voxelsize) - voxminZ);                     # z


## Calculate number of voxels in lesion model (phantom) from number of voxels in ROI and recon_voxelsize, phantom voxelsize

numVox_CClus = int((args.roisizeCClus * args.recon_voxelsize)/args.voxelsize);
numVox_Spic  = int((args.roisizeSpic * args.recon_voxelsize)/args.voxelsize);


x1 = np.zeros(([1,args.totnumROI]),dtype=int);
y1 = np.zeros(([1,args.totnumROI]),dtype=int);
z1 = np.zeros(([1,args.totnumROI]),dtype=int);

## Lesion model

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
maskflag = 0;	# flag used for checking overlap in maskimg, 1 if ROI available (overlap < 10%), 0 if ROI not available due to > 10% overlap


if np.mod(args.totnumROI,2) == 0: # if total number of ROIs is even 
	myhalfnum = args.totnumROI/2
else: 
	myhalfnum = (args.totnumROI-1)/2
		
myLesionType = np.zeros(([args.totnumROI]),dtype=int);  

# initialize myLesionType - divide the total number of ROIs between how many to be extracted for sizes roisizeCClus and roisizeSpic
for ii in range (0,myhalfnum):
	myLesionType[ii] = 0 	# indicating ROI size 'roisizeCClus'
	
for ii in range (myhalfnum,args.totnumROI):
	myLesionType[ii] = 1 	# indicating ROI size 'roisizeSpic'
	
ctr1 = 0;	# counter for myLesionType
eliminated_locs = 0     # counter to check how many locations eliminated

# array to store voxel number and type for the lesions inserted
InsLesions = np.ndarray((args.totnumROI,4),dtype = int) 
InsLesions[:,0:4] = 0;

maxattempt = 100;

# variables to check percentage of overlap between same roi/lesion type
tot_vox_overlap = 0;
percent_overlap = 0.0;

# create a mask for projection image for checking ROI overlap
maskimg = np.zeros(([args.dim_x,args.dim_y,args.dim_z]),dtype=dt);



## MAIN LOOP STARTS
while (insertedLesions < args.totnumROI) and (eliminated_locs < (numloc-1)):

	# pick a location for lesion insertion randomly
	myrnd = np.random.random_integers(numloc-1)
	lesloc[insertedLesions,0] = lesvox[myrnd,0] 
	lesloc[insertedLesions,1] = lesvox[myrnd,1]
	lesloc[insertedLesions,2] = lesvox[myrnd,2]
	lesloc[insertedLesions,3] = myLesionType[insertedLesions];


	#### START CHECKING IF CURRENT LOCATION PASSES THE CHECKS FOR LESION INSERTION
	currX = lesvox[myrnd,2];	 
	currY = lesvox[myrnd,1];	 
	currZ = lesvox[myrnd,0];	# interchaging positions of pos of x, y and z due to reading phantom in this way

	# reset checkFail
	checkFail = 0;
	selectROI = 0;

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

		###### CHECK IF LESION IS TOO CLOSE TO THE CHEST WALL - HARDCODED X USING AMIDE VIEW OF RAW PHANTOM (diff (mm) from min x to edge of upper lip of paddle is converted to voxels)
		## (assuming all breats in one category would have chest wall around same locations, thus same hardcoded X could be used)
                check_chestwall_X = args.limit_chestX_vox;   # in voxels

		if (currZ < check_chestwall_X):       # REMEMBER X AND Z HAVE BEEN INTERCHANGED, SO WE ARE ACTUALLY CHECKING IN X DIMENSION USING currZ
                        lesvox[myrnd,3] = 1;          # lesion discarded due to being too close to chest wall
                        eliminated_locs = eliminated_locs + 1;
									

		####### make changes for odd/even voxel lesion - FOR ODD SIZED LESIONS - CONSIDERING HFLES TO BE SIZE+1/2 -> SO THE BELOW SECTION WILL ADD ONE EXTRA VOXEL

		## length in X dimension correspond to 30 slices, so it is not same as Y,Z length. convert slice range to voxels (phantom)
		if (myLesionType[ctr1]==0):
			hfles_Xmin = int((args.CClus_SliceOffsetRem*args.sl_thick)/(args.voxelsize))
			hfles_Xmax = int((args.CClus_SliceOffsetAdd*args.sl_thick)/(args.voxelsize))
		elif (myLesionType[ctr1]==1):
			hfles_Xmin = int((args.Spic_SliceOffsetRem*args.sl_thick)/(args.voxelsize)) # (voxel number)
			hfles_Xmax = int((args.Spic_SliceOffsetAdd*args.sl_thick)/(args.voxelsize))


		# check whether the lesion volume lies within the phantom boundaries	
		if (int(currX-hfles_Xmin) >= 0) and (int(currX+hfles_Xmax-1) < px) and (int(currZ-hfles) >= 0) and (int(currZ+hfles-1) < pz) and (int(currY-hfles) >= 0) and (int(currY+hfles-1) < py):

			## Perform required checks to determine if lesion can be inserted at current location

			# Create a temporary phantom 3d array - equal to size of lesion - contains phantom voxels corresponding to where lesion to be inserted
			tmp_phan=np.zeros((((currX+hfles_Xmax+1)-(currX-hfles_Xmin)),les,les),dtype=np.uint8) # assign 3D array (with lesion dimensions) of zeros with type uint8

			# Copy particular volume of phantom (where lesion to be placed) from phantom to tmp_phan
			zcounter = 0
			for myx in range(int(currX-hfles_Xmin),int(currX+hfles_Xmax+1)):      # give the range [low,high) - high value is not inclusive, in x,y,z (3d volume) where lesion to be inserted
	                	for myy in range(int(currY-hfles),int(currY+hfles)):
	       	                	for myz in range(int(currZ-hfles),int(currZ+hfles)):
						tmp_phan[zcounter,abs(currY-hfles-myy),abs(currZ-hfles-myz)] = phantom[myx,myy,myz];	            
				zcounter = zcounter + 1;	# counter for X dim


			# Preliminary check - check the eight corners of the lesion cube (in tmp_phan) to see if they lie in any of the foll (air,skin,nipple,muscle) - if yes, right away reject it because then the rest of the points will also not comply
			# if this is passed, then go for detailed check below
			corners=[0,les-1]
			cornersX = [0,(hfles_Xmax-hfles_Xmin)]
			for myxC in cornersX:
				for myyC in corners:
					for myzC in corners:
						if (tmp_phan[myxC,myyC,myzC] == 0) or (tmp_phan[myxC,myyC,myzC] == 2) or (tmp_phan[myxC,myyC,myzC] == 33) or (tmp_phan[myxC,myyC,myzC] == 40): # reject this lesion location
                                                        lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable for future due to not passing the checks
                                                        break;



                        # Detailed check: Check if phantom where lesion to be inserted contains air(0), skin(2), nipple (33), or muscle (40) - if yes, do not insert lesion, choose another random location
                        # !!!! DEBUG !!!! CHECK ONLY SURFACES OF THE CUBE

                        if lesvox[myrnd,3] != 1:
				# check for x = 0, les-1
				for myxT in (0,(hfles_Xmax - hfles_Xmin + 1)): # REMOVED 'RANGE()'
					for myyT in range(0,les):
						for myzT in range(0,les):
							if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40): # reject this lesion location
                                                                lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable for future due to not passing the checks
                                                                break;

				# check for y = 0, les-1
				for myyT in (0,les-1): # REMOVED 'RANGE()'
					for myxT in range(0,(hfles_Xmax - hfles_Xmin + 1)):
						for myzT in range(0,les):
							if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40):
							        lesvox[myrnd,3] = 1;
							        break;


				# check for z = 0, les-1
				for myzT in (0,les-1): # REMOVED 'RANGE()'
					for myxT in range(0,(hfles_Xmax - hfles_Xmin + 1)):
						for myyT in range(0,les):
							if(tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40):
								lesvox[myrnd,3] = 1;
								break;
																																																																			

                        if lesvox[myrnd,3] == 1:                             
                                checkFail = 1;          # lesion failed the checks
				eliminated_locs = eliminated_locs + 1;
                                print("check failed - Air/Skin/Nipple/Muscle %d,%d,%d") % (myxC,myyC,myzC)


			# insert lesion at current location
			if checkFail == 0:	# Potential ROI available

				### Phantom is reshaped as ZYX
				lesZ = lesloc[insertedLesions,2];    # Z
				lesY = lesloc[insertedLesions,1];    # Y
				lesX = lesloc[insertedLesions,0];    # X

				lesZ = lesZ - args.cr_z;	#!!!! DETECTOR MOVED CR_Z VOXELS UP, SO REDUCE THAT MANY Z PLANES (note	that lesX contains the z voxel number for lesion)



			        ## Since we know the lesion location in phantom voxels, so calculate distance (in voxels) from the origin (vox_loc - 0) [given origin vox in phantom = 0,0,0]
				## convert distance from vox to mm = (vox_loc-0)* voxelsize
				## convert this distance to DBT pixels by dividing with appropriate pixelsize (x=y=DBT pixelsize; z=slice thickness).
			        ## This would be DBT pixel location assuming the origin point for DBT is also 0,0,0.

				lesvoxX_dbt = (int)((lesX * args.voxelsize)/(args.recon_voxelsize)) # in recon voxels
				lesvoxY_dbt = (int)((lesY * args.voxelsize)/(args.recon_voxelsize))
				lesvoxZ_dbt = (int)((lesZ * args.voxelsize)/(args.sl_thick))


				## since MC-GPU flips Y axis, we have to find mirror (against center of y axis) of the lesion location in Y
				## For us Y is the longer of the two dimensions (x and y).  Rongping inputs it as 1136x483x38, this dim_x (longer) is our Y
				centerpix = (int)(args.dim_x/2)
				mydiffpix = (lesvoxY_dbt - centerpix)
				if (mydiffpix < 0): # y1[0,ii] lies to the left of center of phantom
					lesvoxY_dbt = centerpix + abs(mydiffpix)
				elif(mydiffpix >= 0): # lies to the right
				        lesvoxY_dbt = centerpix - abs(mydiffpix)


			        ## interchange x and y
 		                temp2 = lesvoxX_dbt
				lesvoxX_dbt = lesvoxY_dbt
				lesvoxY_dbt = temp2


																	
				# MASKING - mark this region as '1's/'2's in the mask based on ROI type (calc clus - 1, spic - 2), and check for overlap
				roipixel = np.zeros(([1,2]),dtype=int);		# number of pixels in the roi based on the lesion type {calc cluster (0), spiculated mass (1)}
				roipixel[0,0] = args.roisizeCClus;
				roipixel[0,1] = args.roisizeSpic;

				myctrpix = int((roipixel[0,lesloc[insertedLesions,3]] - 1)/2)

                                if (myLesionType[ctr1] == 0) and (insertedLesions == 0):	## FIRST CC LESION INSERTED
                                	for myx in range (int(lesvoxX_dbt-myctrpix), int(lesvoxX_dbt+myctrpix+1)):     
                                        	for myy in range(int(lesvoxY_dbt-myctrpix), int(lesvoxY_dbt+myctrpix+1)):
							for myz in range(int(lesvoxZ_dbt-args.CClus_SliceOffsetRem), int(lesvoxZ_dbt+args.CClus_SliceOffsetAdd+1)):
								maskimg[myx,myy,myz] = 1;     # CC
								maskflag = 1;		# make ROI available

                                elif (myLesionType[ctr1] == 1) and (insertedLesions == int(args.totnumROI/2)):	## FIRST SPIC LESION INSERTED
                                        for myx in range (int(lesvoxX_dbt-myctrpix), int(lesvoxX_dbt+myctrpix+1)):                                     
                                                for myy in range(int(lesvoxY_dbt-myctrpix), int(lesvoxY_dbt+myctrpix+1)):
							for myz in range(int(lesvoxZ_dbt-args.Spic_SliceOffsetRem), int(lesvoxZ_dbt+args.Spic_SliceOffsetAdd+1)):
								maskimg[myx,myy,myz] = 2;     # SPIC
								maskflag = 1;

                                elif (myLesionType[ctr1] == 0) and (insertedLesions > 0):		## CC LESION already exist, calculate % overlap
					tot_vox_overlap = 0 # reset
					percent_overlap = 0
					
                                	for myx in range (int(lesvoxX_dbt-myctrpix), int(lesvoxX_dbt+myctrpix+1)): 	    
                                        	for myy in range(int(lesvoxY_dbt-myctrpix), int(lesvoxY_dbt+myctrpix+1)):
							for myz in range(int(lesvoxZ_dbt-args.CClus_SliceOffsetRem), int(lesvoxZ_dbt+args.CClus_SliceOffsetAdd+1)):
								if(maskimg[myx,myy,myz] == 1):      
									tot_vox_overlap = tot_vox_overlap + 1;

					percent_overlap = tot_vox_overlap/(args.roisizeCClus * args.roisizeCClus * args.roisizeCClus)

					if (percent_overlap > 0.1):
						maskflag = 0;	# not available
					else:
						maskflag = 1;	# available

                                elif (myLesionType[ctr1] == 1) and (insertedLesions > int(args.totnumROI/2)):		## SPIC LESION already exist, calculate % overlap
					tot_vox_overlap = 0 # reset
					percent_overlap = 0
					
                                        for myx in range (int(lesvoxX_dbt-myctrpix), int(lesvoxX_dbt+myctrpix+1)):                                     
                                                for myy in range(int(lesvoxY_dbt-myctrpix), int(lesvoxY_dbt+myctrpix+1)):
							for myz in range(int(lesvoxZ_dbt-args.Spic_SliceOffsetRem), int(lesvoxZ_dbt+args.Spic_SliceOffsetAdd+1)):
								if(maskimg[myx,myy,myz] == 2):
									tot_vox_overlap = tot_vox_overlap + 1;
					
					percent_overlap = tot_vox_overlap/(args.roisizeSpic * args.roisizeSpic * args.roisizeSpic)

					if (percent_overlap > 0.1):
						maskflag = 0;	# not available
					else:
						maskflag = 1;	# available





				## CHECK WHETHER (maskflag == 1) and (ROI VOLUME IS WITHIN THE BOUNDS OF DBT RECON)
				if (myLesionType[ctr1] == 0):
					midXY = int((args.roisizeCClus-1)/2)	# half the number of voxels in XY dim
					midZ = int((args.slice3dCClus-1)/2)	# half the number of slices
				elif (myLesionType[ctr1] == 1):
					midXY = int((args.roisizeSpic-1)/2)
					midZ = int((args.slice3dSpic-1)/2)

				if (int(lesvoxX_dbt-midXY) >= 0) and (int(lesvoxX_dbt+midXY) < args.dim_x) and (int(lesvoxY_dbt-midXY) >= 0) and (int(lesvoxY_dbt+midXY) < args.dim_y) and (int(lesvoxZ_dbt-midZ) >= 0) and (int(lesvoxZ_dbt+midZ) < args.dim_z) and (maskflag == 1):
					selectROI = 1;	# select this ROI
					ctr1 = ctr1 + 1;
				else:
					selectROI = 0;  # ROI extraction failed at this location


				
				if (selectROI == 1):				
					# DO NOT interchange x and y - x correspond to args.dim_x dim, and y to args.dim_y dim
					x1[0,insertedLesions] = lesvoxX_dbt   # 1300 dim
					y1[0,insertedLesions] = lesvoxY_dbt   # 1000 dim (recon volume in voxels)
					z1[0,insertedLesions] = lesvoxZ_dbt   # 30 slices

					print x1[0,insertedLesions], y1[0,insertedLesions], z1[0,insertedLesions]
					floc.write('%d %d %d %d\n' % (lesvoxX_dbt, lesvoxY_dbt, lesvoxZ_dbt, insertedLesions))

					## extract roi


					# 2D ROI extraction
					if (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) == 0):		# extracted roi image with the size based on the lesion type
						myroi = np.zeros( ((args.roisizeCClus+1),(args.roisizeCClus+1)),dtype=dt ) # if ROI is even pixel, then add 1 pixel. we want tmpX,tmpY to lie in center
					elif (myLesionType[ctr1-1] == 0) and (np.mod(args.roisizeCClus,2) != 0):
						myroi = np.zeros( (args.roisizeCClus,args.roisizeCClus),dtype=dt )
						
					elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) == 0):
                                                myroi = np.zeros( ((args.roisizeSpic+1),(args.roisizeSpic+1)),dtype=dt )
					elif (myLesionType[ctr1-1] == 1) and (np.mod(args.roisizeSpic,2) != 0):
						myroi = np.zeros( (args.roisizeSpic,args.roisizeSpic),dtype=dt )
					
					mypix = int((roipixel[0,lesloc[insertedLesions,3]] - 1)/2)	# roipixel[0,lesloc[ii,3]] is the number of pixels to be extracted for the roi; mypix - number of pixels in each direction with x1,y1 at ctr (excluding x1,y1).
					myroi = np.zeros(([roipixel[0,lesloc[insertedLesions,3]],roipixel[0,lesloc[insertedLesions,3]]]),dtype=dt)	# extracted roi image with the size based on the lesion type extracted

					for myx in range(x1[0,insertedLesions]-mypix,x1[0,insertedLesions]+mypix+1):	# we want x1+mypix to be included
						for myy in range(y1[0,insertedLesions]-mypix,y1[0,insertedLesions]+mypix+1):
							myroi[(myx+mypix-x1[0,insertedLesions]),(myy+mypix-y1[0,insertedLesions])] = myimg[myx,myy,z1[0,insertedLesions]]; 

					# SAVE ROIs except for spiculated mass
					if (lesloc[insertedLesions,3] == 1):
						myroi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED ROIs FOR SPICULATED MASS/roiDBT_SA_"+str(insertedLesions)+"_"+str(args.RndSeed)+"_2D.raw")

					roipix3d = np.zeros(([1,2]),dtype=int);	# number of slices for 3d volume for each lesion type
					roipix3d[0,0] = args.slice3dCClus;	# calc cluster
					roipix3d[0,1] = args.slice3dSpic;	# spiculated mass


					# 3D VOI extraction
					mypix = int((roipixel[0,lesloc[insertedLesions,3]] - 1)/2)	# roipixel[0,lesloc[ii,3]] is the number of pixels to be extracted for the roi

					# find center slice for 3d volume
					if((roipix3d[0,lesloc[insertedLesions,3]] % 2) == 0):	# even	
						my3Dpix = roipix3d[0,lesloc[insertedLesions,3]]/2
					else:						# odd
						my3dpix = (roipix3d[0,lesloc[insertedLesions,3]] - 1)/2

					my2Droi = np.zeros(([roipixel[0,lesloc[insertedLesions,3]],roipixel[0,lesloc[insertedLesions,3]]]),dtype=dt)				# 2d roi for cluster calc - obtained from avg 3d voi
					my3Droi = np.zeros(([roipix3d[0,lesloc[insertedLesions,3]],roipixel[0,lesloc[insertedLesions,3]],roipixel[0,lesloc[insertedLesions,3]]]),dtype=dt)	# 3D volume - extracted voi with the size based on the lesion type extracted
		
					if (lesloc[insertedLesions,3] == 0):		# calc cluster - slice range [center-5, center+2] both inclusive
						for myx in range(x1[0,insertedLesions]-mypix,x1[0,insertedLesions]+mypix+1):	# we want x1+mypix to be included
							for myy in range(y1[0,insertedLesions]-mypix,y1[0,insertedLesions]+mypix+1):
								zcounter = 0
								for myz in range(z1[0,insertedLesions]-args.CClus_SliceOffsetRem,z1[0,insertedLesions]+args.CClus_SliceOffsetAdd+1):        
				                                        my3Droi[zcounter,(myx+mypix-x1[0,insertedLesions]),(myy+mypix-y1[0,insertedLesions])] = myimg[myx,myy,myz]
			                                        	zcounter = zcounter + 1
							
						my2Droi = np.mean(my3Droi,axis=0)	# take mean over Z axis
						MIP2Droi = np.max(my3Droi,axis=0)       # maximum intensity projection along z axis

						# SAVE 2D roi
						my2Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED 2D MEAN ROIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SA_"+str(insertedLesions)+"_"+str(args.RndSeed)+"_mean_2D.raw")
						MIP2Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED 2D MIP ROIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SA_"+str(insertedLesions)+"_"+str(args.RndSeed)+"_mip_2D.raw")

					elif(lesloc[insertedLesions,3] == 1):	# spiculated masses, slice range (center +- 4) both inclusive
						for myx in range(x1[0,insertedLesions]-mypix,x1[0,insertedLesions]+mypix+1):	# we want x1+mypix to be included
							for myy in range(y1[0,insertedLesions]-mypix,y1[0,insertedLesions]+mypix+1):
								zcounter = 0
								for myz in range(z1[0,insertedLesions]-args.Spic_SliceOffsetRem,z1[0,insertedLesions]+args.Spic_SliceOffsetAdd+1):       
						      	                my3Droi[zcounter,(myx+mypix-x1[0,insertedLesions]),(myy+mypix-y1[0,insertedLesions])] = myimg[myx,myy,myz]
					                                zcounter = zcounter + 1


			


					# SAVE VOIs
					if (lesloc[insertedLesions,3] == 0):
						my3Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED VOIs FOR MICROCALCIFICATION CLUSTER/roiDBT_SA_"+str(insertedLesions)+"_"+str(args.RndSeed)+"_3D.raw")
					elif (lesloc[insertedLesions,3] == 1):
						my3Droi.astype(dt).tofile("PATH TO FOLDER CONTAINING EXTRACTED VOIs FOR SPICULATED MASS/roiDBT_SA_"+str(insertedLesions)+"_"+str(args.RndSeed)+"_3D.raw")

					insertedLesions = insertedLesions + 1;
                                        eliminated_locs = eliminated_locs + 1;

                                elif (selectROI == 0): # mask overlap failed
                                        eliminated_locs = eliminated_locs + 1;
                                        print "potential ROI failed mask overlap check...\n"

                else: # not within phantom boundaries 
                        eliminated_locs = eliminated_locs + 1
                        print "Lesion out of phantom bounds...\n"

        else:
                eliminated_locs = eliminated_locs + 1
                print "Lesion voxel unavailable for use...\n"
                                        
# Delete pc_*.raw file under /tmp
runcmd3 = "/tmp/"+flname2
os.remove(runcmd3)
