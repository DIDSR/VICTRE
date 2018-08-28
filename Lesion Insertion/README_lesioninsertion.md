# VICTRE - Lesion Insertion
*Author: Diksha Sharma (US Food and Drug Administration)*

*Questions should be directed to: diksha.sharma@fda.hhs.gov.*

This directory contains the lesion insertion python code and build script. This part of the pipeline comes after the phantom is compressed and cropped. Lesions are inserted in compressed breast phantoms to create cancer cases. The lesion insertion locations are randomly chosen from the list of possible locations given by the breast phantom generation code. The selected location is then passed through checks to ensure that the lesion is within the phantom boundaries, non-overlapping with tissues like air/muscle/nipple/skin, and non-overlapping with already inserted lesions.


Disclaimer
----------

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 


Folder contents
---------------
- lesionInsertion.py - lesion insertion python code
- README_lesioninsertion.md - this file.
- /build/lesionInsertion_script.sh - bash script for executing the lesion insertion code.
- /VICTRE_LesionModels/heteroCalc2_121_100.raw - sample microcalification cluster used for VICTRE.
- /VICTRE_LesionModels/mass_-308854003_cropped_166.raw - spiculated mass used for VICTRE.


Pre-requisites
--------------
- Python - this code has been tested on Python 2.7.5.
- Output files from phantom generation, compression and cropping - pc_SEED_crop.raw.gz, pc_SEED_crop.mhd (This SEED is unique to each phantom and is used as the initial random seed for phantom generation. We use this SEED to indicate which phantom to process).
- Lesion models (as raw files) - current version inserts two types of lesion models (microcalcification cluster and spiculated mass).
- The following input parameters are required to be passed as command line arguments - phantom seed, focal spot of the imaging system, voxel size of the phantom, phantom origin, lowest bounds for phantom coordinate axes, total number of voxels in the phantom in x/y/z, lesion length in voxels for the different types of lesions.


Output
------
- pcl_SEED_crop.raw.gz- gzip raw file containing lesion inserted phantoms.
- pcl_SEED_crop.loc - text file containing the inserted lesion locations (in voxels) for that phantom.


Caveat
------
The current version tries to insert maximum 8 lesions (two types of lesions (microcalcification cluster and spiculated mass); 4 of each type). To insert a different number of lesions, make the following edits to the code - look for line 'numLesToInsert = 8' (replace 8 with th total number of lesions you want);  'myLesionType = [0,0,0,0,1,1,1,1]' (flags 0 and 1 indicate the two types of lesions and 4 of each type; this order of 0's and 1's dictate what type of lesion is inserted when.  Edit this array to make sure the total number of 0's and 1's equal the number in numLesToInsert.)

To insert only one type of lesion, give same lesion file name under the section '## Lesion model'. To insert more than two types of lesions, the code will need further modifications.
