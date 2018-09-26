# VICTRE - ROI Extraction
*Author: Diksha Sharma (US Food and Drug Administration)*

*Questions should be directed to: diksha.sharma@fda.hhs.gov*

This code extracts the regions of interest (ROIs) for both, lesion-present and lesion-absent mammography projection images (or volumes of interest (VOIs) from the DBT volumes). For extracting lesion-absent ROIs, we applied rigorous checks including whether the subimages are within the reconstructed volume boundaries and non-overlapping, to find the appropriate locations.

Disclaimer
----------
This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

Folder contents
---------------
- Source codes - roiExtraction_mammo_SA.py (mammo signal absent), roiExtraction_mammo_SP.py (mammo signal present), roiExtraction_DBT_FBP_SA.py (DBT signal absent), roiExtraction_DBT_FBP_SP.py (DBT signal present).
- build/*.sh - bash scripts for execution
- README_roiextraction.md - this file

Pre-requisites
--------------
- Python - this code has been tested on Python 2.7.5.
- Needs (some of) the following inputs depending on the code to run: phantom raw file pc/pcl_SEED_crop.raw.gz (This SEED is unique to each phantom and is used as the initial random seed for phantom generation. We use this SEED to indicate which phantom to process), MC-GPU projection output (as .raw), flatfield projection output (as .raw).
- Output: ROIs (as .raw), ROI locations file (as .loc).  VOIs are output from roiExtraction_DBT_FBP_SA/SP.py (roiDBT_SA/SP_*_3D.raw).

Execution
---------
- The build/ folder contains scripts for execution - each script requires the phantom seed number (taken from its filename; this seed number is unique per phantom and thus determines which phantom to run) to be provided; user must modify other input arguments for their runs accordingly.
- The user needs to insert the path to input and output files (such as path to phantom raw image, potential lesion locations file, MC-GPU projection image as raw, flatfield image as raw, folder where extracted ROIs and their locations will be saved) in the source codes and build scripts.

Visualization of ROI/VOIs
-------------------------
The outputs are .raw files which can be visualized in any supporting tool. Image format - Mammo ROIs are in 32-bit real, little-endian byte order; DBT ROI/VOIs are in 64-bit real, little-endian byte order format.

One of the popular tools which can be used is ImageJ. 

ImageJ Example: 
	
		Visualizing Mammo ROIs containing microcalcification (ran with default parameters):

		Image type = 32-bit real
		Width = 65
		Height = 65
		Offset to first image = 0
		Number of images = 1
		Gap between images = 0
		Little-endian byte order - CHECKED
		
		
		Visualizing DBT VOI containing spiculated mass (ran with default parameters):
		
		Image type = 64-bit real
		Width = 109
		Height = 109
		Offset to first image = 0
		Number of images = 9
		Gap between images = 0
		Little-endian byte order - CHECKED
