#################################################################

				LESION INSERTION
				FATTY PHANTOMS


AUTHOR: 	DIKSHA SHARMA
		DIKSHA.SHARMA@FDA.HHS.GOV

#################################################################

#################################################################

					DISCLAIMER

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

#################################################################
# import numpy scientific computing package
import numpy as np
import argparse
import os
import os.path
import gzip

# Change working directory
os.chdir('/raidb/VICTRE/pivotal/phantoms/withLesion/fatty/')


# Parsing input arguments using argparse
parser = argparse.ArgumentParser()
parser.add_argument("RndSeed",type=int)
parser.add_argument("fs_x",type=float)  	# focal spot locations (in mm)
parser.add_argument("fs_y",type=float)
parser.add_argument("fs_z",type=float)
parser.add_argument("voxelsize",type=float)	# voxel size (mm)
parser.add_argument("orgX",type=float)		# origin for phantom (mm)
parser.add_argument("orgY",type=float)
parser.add_argument("orgZ",type=float)
parser.add_argument("minvoxX",type=int)		# minimum voxel number
parser.add_argument("minvoxY",type=int)
parser.add_argument("minvoxZ",type=int)
parser.add_argument("totvoxX",type=int)		# total number of voxels
parser.add_argument("totvoxY",type=int)
parser.add_argument("totvoxZ",type=int)
parser.add_argument("cclus_voxlen",type=int)	# lesion length in voxels for calc cluster
parser.add_argument("spic_voxlen",type=int)	# spiculated mass
parser.add_argument("limit_chestX_vox",type=float)	# min X for lesions to make sure they are not close to chest wall
args = parser.parse_args()

# Use system calls to unzip phantom		
fl1 = "pc_"
fl2 = str(args.RndSeed)
fl3 = "_crop.raw.gz"
fl4 = "_crop.raw"
fl5 = "pcl_"

flname = fl1+fl2+fl3
flname2 = fl1+fl2+fl4
flname3 = fl5+fl2+fl4

print flname
print flname2
print flname3

runcmd1 = "gunzip -c "+flname+" > /tmp/"+flname2 # gunzip -c file.raw.gz > /tmp/file.raw
print runcmd1
os.system(runcmd1)


# read binary file (compressed phantom) as unsigned int 8 bits and reshape in a matrix of size 1051x1657x713
phantom = np.fromfile("/tmp/pc_"+str(args.RndSeed)+"_crop.raw", dtype=np.uint8).reshape(args.totvoxZ,args.totvoxY,args.totvoxX);

# check the size of matrix 
px,py,pz = phantom.shape

# read lesion locations file (potential locations for inserting lesions)
x, y, z = np.loadtxt("pc_"+str(args.RndSeed)+".loc", delimiter=',', usecols=(0, 1, 2), unpack=True);	# these are positions in mm, which needs to be changed to voxel number

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


# convert lesion location (in mm) to lesion voxel space 
# Voxel number = {[curr pos - 0,0,0 pos (in mm)]/voxel size} - min voxel

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


## Lesion model

lesionCalcClus=np.fromfile("/raidb/VICTRE/pivotal/phantoms/lesionModels/heteroCalc2_121_100.raw",dtype=np.uint8).reshape(args.cclus_voxlen,args.cclus_voxlen,args.cclus_voxlen) 	# microcalcification cluster
lesionSpic=np.fromfile("/raidb/VICTRE/pivotal/phantoms/lesionModels/mass_-308854003_cropped_166.raw",dtype=np.uint8).reshape(args.spic_voxlen,args.spic_voxlen,args.spic_voxlen) 	# spiculated lesion

myoddflag = np.zeros(([1,2]),dtype=int);		# raise this if lesion is odd sized - each element correspond to each lesion type - calcCluster(0), spiculated (1)

kk,ll,mm = lesionCalcClus.shape
lesCalcClus = kk                    # lesion size - assuming all dims are equal length
if np.mod(lesCalcClus,2) == 0:      # even
        hflesCalcClus=lesCalcClus/2
else:                		    # odd
        hflesCalcClus=(lesCalcClus+1)/2
	lesCalcClus = lesCalcClus+1;
	myoddflag[0,0] = 1;

kk1,ll1,mm1 = lesionSpic.shape
lesSpic = kk1  # lesion size - assuming all dims are equal length
if np.mod(lesSpic,2) == 0:
	hflesSpic=lesSpic/2
else:
	hflesSpic=(lesSpic+1)/2
	lesSpic = lesSpic+1;
	myoddflag[0,1] = 1;


# open file with locations for appending
floc = open("pcl_"+str(args.RndSeed)+"_crop.loc",'w')


## minimum distance (voxels) for checking projection overlap - if less than this distance then discard location
mindist_CC = args.cclus_voxlen+5   
mindist_CM = ( (int)((args.cclus_voxlen+args.spic_voxlen)/2) ) + 5 
mindist_MM = args.spic_voxlen   # between mass and mass


##### if loop to insert 8 lesions
insertedLesions = 0	# variable to indicate number of inserted lesions so far
checkFail = 0	# flag to indicate if lesion passed the checks or not (0 pass or not tried yet, 1 fail)
myctr = 0 	# a counter to prevent while loop going in infinity

numLesToInsert = 8  	# number of lesions to be inserted in total 

myLesionType = [0,0,0,0,1,1,1,1] # 8 lesions - 4 calcs ans 4 masses

ctr1 = 0;	# counter for myLesionType

eliminated_locs = 0     # counter to check how many locations eliminated

# array to store voxel number and type for the lesions inserted
InsLesions = np.ndarray((numLesToInsert,4),dtype = int) 
InsLesions[:,0:4] = 0;

maxattempt = 100;


# Initialize random seed generator (using /dev/urandom), goal to have reproducible runs if needed
# for reproducible runs, set (hardcode) urand_seed[0] which will then initialize np.random.seed sequence to generate same random sequence every time
finrnd=open("/dev/urandom","r")
urand_seed=np.fromfile(finrnd,dtype=int,count=1)
np.random.seed(abs(urand_seed[0]))

print urand_seed[0]

### check while lesions are yet to be inserted, or did not run out of potential locations in the list or myctr < maxattempt 
while (insertedLesions < numLesToInsert) and (eliminated_locs < (numloc-1)):
	
	# pick a location for lesion insertion randomly
	myrnd = np.random.random_integers(numloc-1)
	currX = lesvox[myrnd,2];
	currY = lesvox[myrnd,1];
	currZ = lesvox[myrnd,0];	# interchaging positions of pos of x, y and z due to reading phantom in this way (lesion locations are based on 1051,1657,713, while we are reading phantom as 713,1657,1051


	# reset checkFail
	checkFail = 0;

	## insert lesion
	# choose type of lesion and accordingly change variables
        if (myLesionType[ctr1]==0) and (myoddflag[0,0] == 0):
                hfles = hflesCalcClus;
                les = lesCalcClus;
                lesion = lesionCalcClus;
	elif (myLesionType[ctr1]==1) and (myoddflag[0,1] == 0):
		hfles = hflesSpic;
		les = lesSpic;
		lesion = lesionSpic;
        elif (myLesionType[ctr1]==0) and (myoddflag[0,0] == 1):		# ODD SIZED LESIONS - add extra plane of zeros in all three dimensions
                hfles = hflesCalcClus;
                les = lesCalcClus;
                newtmp = np.zeros((les,les,les), dtype=int);
                newtmp[0:les-1,0:les-1,0:les-1] = lesionCalcClus;
                lesion = newtmp;
        elif (myLesionType[ctr1]==1) and (myoddflag[0,1] == 1):
                hfles = hflesSpic;
                les = lesSpic;
                newtmp = np.zeros((les,les,les), dtype=int);
                newtmp[0:les-1,0:les-1,0:les-1] = lesionSpic;
                lesion = newtmp;
	
		
		
	# is this lesion voxel available for use (flag == 0)?
	if lesvox[myrnd,3] == 0:

		# location number available for use
		print("Eliminated %d/%d") % (eliminated_locs, numloc)

		###### CHECK IF LESION IS TOO CLOSE TO THE CHEST WALL - HARDCODED X USING AMIDE VIEW OF RAW PHANTOM (diff (mm) from min x to edge of upper lip of paddle is converted to voxels)
		## (assuming all breats in one category would have chest wall around same locations, thus same hardcoded X could be used)
		check_chestwall_X = args.limit_chestX_vox;   # in voxels

		if (currZ < check_chestwall_X):       # REMEMBER X AND Z HAVE BEEN INTERCHANGED, SO WE ARE ACTUALLY CHECKING IN X DIMENSION USING currZ
			lesvox[myrnd,3] = 1;          # lesion discarded due to being too close to chest wall
			checkFail = 1;
			myctr = myctr + 1;
			print("check failed - close to chestwall \n")

		####### make changes for odd/even voxel lesion - FOR ODD SIZED LESIONS - CONSIDERING HFLES TO BE SIZE+1/2 -> SO THE BELOW SECTION WILL ADD ONE EXTRA VOXEL

		#### IS LESION WITHIN PHANTOM BOUNDARIES? IF YES, PERFORM FURTHER CHECKS ELSE GO TO ELIF LOOP.	
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
			# CHECK ONLY SURFACES OF THE CUBE
			
			if lesvox[myrnd,3] != 1:
				# check for x = 0, les-1
				for myxT in (0,les-1): # REMOVED 'RANGE()'
					for myyT in range(0,les):
						for myzT in range(0,les):
							if (tmp_phan[myxT,myyT,myzT] == 0) or (tmp_phan[myxT,myyT,myzT] == 2) or (tmp_phan[myxT,myyT,myzT] == 33) or (tmp_phan[myxT,myyT,myzT] == 40): # reject this lesion location
								lesvox[myrnd,3] = 1;  	# mark this lesion to be unavailable for future due to not passing the checks
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
					print("check failed - detailed check - contains Air/Skin/Nipple/Muscle \n")



			# Check for overlapping lesions in Coronal plane (currY,currZ plane (1657,1051))
			# the overlapping distance varies based on the lesion type (ucalc cluster (args.cclus_voxlen), spiculated (args.spic_voxlen))
			if lesvox[myrnd,3] != 1:			# lesion insertion has not failed yet
				# calculate distance in yz between all inserted lesions and current yz, atleast > lesion type voxvolume, then insert, else look for another location
				if ctr1 > 0:				# atleast one lesion has been inserted
					dist = np.sqrt( np.power((InsLesions[ctr1-1,1] - currY),2) + np.power((InsLesions[ctr1-1,2] - currZ),2) );	# in voxels
					if (myLesionType[ctr1] == 0) and (dist < (2*args.cclus_voxlen)):                  	# ucalc cluster - if less than les volume then they are considered to be overlapping
                                                checkFail = 1;          # Lesion insertion failed!!
                                               	lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable
                                                myctr = myctr + 1;	# counter to prevent while loop going in infinity - break while loop if consecutive 50 locations have been found inappropriate to insert lesion
                                               	print("check failed - Coronal Overlapping %d,%d,%d - Calc Cluster") % (myxT,myyT,myzT)

                                        elif (myLesionType[ctr1] == 1) and (dist < (2*args.spic_voxlen)):                    	# spiculated - if less than les volume then they are considered to be overlapping
                                                checkFail = 1;          # Lesion insertion failed!!
                                                lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable
                                                myctr = myctr + 1;	# counter to prevent while loop going in infinity - break while loop if consecutive 50 locations have been found inappropriate to insert lesion
                                                print("check failed - Coronal Overlapping %d,%d,%d - Spiculated") % (myxT,myyT,myzT)

                              
                        ### PROJECTION OVERLAP CHECK
                        ## equation of line passing through F and L is f+alpha(l-f). So a point P on this line will be P = f+alpha(l-f).
                        ## We want shortest distance from potential lesion location 'a' to 'P'; (distance)^2 = (ax-px)^2 + (ay-py)^2 + (az-pz)^2.
                        ## But alpha is a variable here, so we solve for diffrential of alpha (it gives min/max values..in this case min is applicable since max could be infinity for alpha).
                        ## so we solve for d/d(alpha) = 0.  Diffrentiation of (ax-px)^2 or (ax-fx-alpha(lx-fx))^2 = 2.(ax-fx-alpha(lx-fx))*(fx-lx)...
                        ## Thus, entire equation would be 2.(ax-fx-alpha(lx-fx))*(fx-lx) + ... = 0.  we solve this to obtain 'alpha', which is then substituted to obtain P
                        ## and then we obtain distance between 'a' and 'P'.

                         
                        if lesvox[myrnd,3] != 1:
				if ctr1 > 0:                            # atleast one lesion has been inserted
					for myd in range(0,ctr1):       # check for all inseted lesions for far
						ax = currZ # potential lesion location
						ay = currY
						az = currX

						# focal spot from Andreu's system (mm) to phantom space (voxels)
						# Origin for both coincide at (0,0,0) vox = Andreu's (0,0) mm
						# Thus, conversion is as simple as divinding mm by voxelsize (mm) and type-casting to int
						fx = (int)(args.fs_x/args.voxelsize)
						fy = (int)(args.fs_y/args.voxelsize)
						fz = (int)(args.fs_z/args.voxelsize)


						lx = InsLesions[myd,2] # existing lesion
						ly = InsLesions[myd,1]
						lz = InsLesions[myd,0]

						tmpnume = (float)( (ax-fx)*(lx-fx) + (ay-fy)*(ly-fy) + (az-fz)*(lz-fz) )
						tmpdeno = (float)( (lx-fx)*(lx-fx) + (ly-fy)*(ly-fy) + (lz-fz)*(lz-fz) )
						alpha = (float)(tmpnume/tmpdeno)

						Px = (float)(fx+alpha*(lx-fx)) # Point P on line passing through F and L; distance between 'a' and 'P' is shortest distance.
						Py = (float)(fy+alpha*(ly-fy))
						Pz = (float)(fz+alpha*(lz-fz))

						mydist = (float)( np.sqrt( (ax-Px)*(ax-Px) + (ay-Py)*(ay-Py) + (az-Pz)*(az-Pz) ) )


						# based on which lesion types are being compared, set the min distance for checking overlap - HARDCODED FOR [C,C,C,C,M,M,M,M] order
						if(ctr1 < 4) and (myd < 3):
							mindistance = mindist_CC;
						elif(ctr1 >= 4) and (myd <= 3):
							mindistance = mindist_CM;
						elif(ctr1 >= 4) and (myd > 3):
							mindistance = mindist_MM;


						if(mydist < (mindistance)):    # lesion too close in projection plane, discard lesion location
							checkFail = 1;          # lesion insertion failed
							lesvox[myrnd,3] = 1;    # mark this lesion to be unavailable
							myctr = myctr + 1;      # counter to prevent while loop going in infinity - break while loop if consecutive 50 locations found inappropriate
							print("check failed - PROJECTION Overlapping  - loc xyz %d,%d,%d, prev loc xyz %d,%d,%d") % (ax,ay,az,lx,ly,lz)


			  
			# insert lesion at current location
			if (lesvox[myrnd,3] != 1) and (checkFail == 0):
				for myxP in range(int(currX-hfles),int(currX+hfles)):	# give the range [low,high) - high value is not inclusive, in x,y,z (3d volume) where lesion to be inserted
					for myyP in range(int(currY-hfles),int(currY+hfles)):
						for myzP in range(int(currZ-hfles),int(currZ+hfles)):
							if lesion[abs(currX-hfles-myxP),abs(currY-hfles-myyP),abs(currZ-hfles-myzP)] == 1:	# only replace where lesion=1, else keep original phantom data
								if (myLesionType[ctr1] == 0):	# ucalc cluster
									phantom[myxP,myyP,myzP] = 250;		# arbitrary value of lesion (constant - no lesion texture)
								elif (myLesionType[ctr1] == 1):	# spiculated mass
									phantom[myxP,myyP,myzP] = 200;		# different xray properties of calcs and masses, therefore different brightness

				print("\n*************** check PASSED zxy %d,%d,%d ******************* Inserted so far %d ****************\n\n") % (currX,currZ,currY,insertedLesions + 1)


				insertedLesions = insertedLesions + 1;  # mark lesion as inserted
				lesvox[myrnd,3] = 1 			# mark as unavailable
				myctr = 0; 				# reset myctr

				InsLesions[ctr1,0] = currX;		# copy voxel location and lesion type to InsLesions - used to check for overlapping lesions
				InsLesions[ctr1,1] = currY;
				InsLesions[ctr1,2] = currZ;
				InsLesions[ctr1,3] = myLesionType[ctr1];

			        floc.write('%d %d %d %d\n' % (currX,currZ,currY,myLesionType[ctr1]))  # append voxels (centers) where lesion inserted - output  z,x,y (713, 1051, 1657)
				ctr1 = ctr1 + 1;			# move to next lesion type

		elif (lesvox[myrnd,3] != 1):	# OUTSIDE PHANTOM BOUNDARY AND NOT TOO CLOSE TO CHESTWALL - MARK UNAVAILABLE
			lesvox[myrnd,3] = 1;    # lesion discarded due to being close to chest wall
			checkFail = 1;          # lesion failed the checks
			print("check failed - OUTSIDE PHANTOM BOUNDARY \n")



		#####  MARK LOCATION UNAVAILABLE - EITHER INSERTED OR CHECK FAILED
		if (checkFail == 1) or (lesvox[myrnd,3] == 1):		# if either lesion was inserted or any of the checks failed - mark location unavailable
			eliminated_locs = eliminated_locs + 1;
			


## print out if ran out of potential locations
if(eliminated_locs >= (numloc-1)):
	print("\n\n\t WHILE LOOP BROKE DUE TO RUNNING OUT OF POTENTIAL LESION LOCATIONS!!!!...............................\n\n")
   
## output pcl_*.raw.gz file
phantom.astype('uint8').tofile(flname3) # save as .raw

runcmd2 = "gzip "+flname3
print runcmd2
os.system(runcmd2) # gzip to raw.gz

# Delete pc_*.raw file under /tmp
runcmd3 = "/tmp/"+flname2
print runcmd3
os.remove(runcmd3)
