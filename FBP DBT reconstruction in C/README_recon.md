
# VICTRE - Digital Breast Tomosynthesis (FBP) in C
*Authors: Aunnasha Sengupta, Rongping Zeng, Diksha Sharma, Aldo Badano (US Food and Drug Administration)*

*Questions should be directed to: diksha.sharma@fda.hhs.gov*

This freely available code performs 3D Filtered Back-Projection (FBP) reconstruction for Digital Breast Tomosynthesis (DBT) using single threaded C. 
It is based on the C reconstruction code developed by Leeser, Mukherjee and Brock (BMC RESEARCH NOTES 7.1, p. 582, 2014), 
which in turn implements the FDK reconstruction algorithm developed by Fessler (http://web.eecs.umich.edu/~fessler/).

Disclaimer
----------
This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

Folder contents
---------------

- Source codes - FBP_DBTrecon.c/.h
- Makefile - to compile the codes
- reconFBP_script.sh - bash script for execution
- cbct_code - Cone beam CBCT reconstruction code developed by Brock.  Modified by Aunnasha for use with VICTRE.
- Readme - this file

Pre-requisites
--------------

- FFTW3 package needs to be installed (on Ubuntu, use the command 'sudo apt-get install libfftw3-dev')
- Needs the following inputs: flatfield projection data, DBT projections for the phantom, and phantom dimensions as voxels (reconFBP_script.sh script needs pc_seed_crop.mhd file for obtaining this information for the phantom).
- Flatfield and DBT projections data format: RAW file, 32-bit real, little-endian byte order. All angle projections concatenated in a single raw file (one file each for flatfield and DBT projections).

Input arguments to the executable (given from command line separated by a blank space)
--------------------------------------------------------------------------------------
 
- Number of projection angles
- Number of detector elements in the direction of the x-ray tube
- Number of detector elements in the direction perpendicular to the x-ray tube movement
- Detector pixel size (in cm)
- Source to detector distance (in cm)
- Source to the rotation center distance (in cm) 
- Detector center offset along the x-ray tube movement direction relative to the x-ray source center (in pixels)
- Complete orbit of projection angles
- Number of voxels in the phantom X direction
- Number of voxels in the phantom Y direction
- Number of voxels in the phantom Z direction
- Phantom voxel size (in cm)
- Reconstruction pixel size (in cm)
- Reconstruction slice thickness in Z dimension (in cm)
- Offset of the volume center to the rotation center in the X direction (in pixels)
- First projection angle in the input dataset
- Difference between the first and second projection angles
- Random seed from the phantom file name (the code reconstructs the projections for this particular phantom)- this is a unique number per phantom and is used as the initial random seed for its generation. This is saved in the phantom file name and is used to tell the ROI extraction code which phantom to process.

Executing the code
-------------------

The Makefile in each subfolder compiles the FBP_DBTrecon.c and creates executable FBP_DBTrecon. 

$ make

$ ./reconFBP_script.sh 

This code outputs the reconstructed DBT image as a raw file (64 bit real, Little-endian byte format). 
We have a sample output placed under sample-input-output-files/output/DBT_pc_36008688_recon-example.raw, and named it so that the code does not overwrite it when run with default parameters.

Visualizing the output
----------------------

There are several ways to view the .raw reconstructed slices. The open source software, ImageJ, is the easiest way. 
To use ImageJ for visualization, it is important to know the dimesions of the output file, which is printed out during the execution of the code as nx, ny and nz. 

ImageJ input parameters: 
Image Type: 64-bit real
Width: nx
Height: ny
Number of images: nz
Check the box: "Little-endian byte order"

Use the slider to view the different reconstructed layers.

Changes from Leeser et. al's codes
----------------------------------

- Interpolation of DBT projections for DBT reconstruction
- Addition of hanning filter
- Corrected existing hanning filter function
- Corrected indexing problem in back-projection function
- corrected usage of fft package
- corrected formula in back-projection algorithm



