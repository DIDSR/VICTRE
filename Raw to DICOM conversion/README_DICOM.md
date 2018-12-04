# VICTRE - Raw to DICOM conversion

*Author: Eshan Dahal (eshan.dahal@fda.hhs.gov)*

The matlab code presented here reads the raw files from the VICTRE data set and creates the DICOM files for each raw file with the appropriate DICOM tags and UIDs. 


Disclaimer
----------

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

The matlab code presented here reads the raw files from the VICTRE data set and creates the DICOM files for each raw file with the appropriate DICOM tags and UIDs. 

Folder contents
---------------
- DICOM_DM.m - converts simulated digitial mamography raw images (32-bit) to DICOM format.
- DICOM_DBTpro.m - converts DBT projection raw images (32-bit) from VICTRE to DICOM format.
- DICOM_Recon.m - converts reconstructed raw images (64-bit) from VICTRE to DICOM format
- README_DICOM.md - this file.

Pre-requisites
--------------
 - Matlab
 - Input files in proper raw format.

Input
-----
1) Folder path to the raw files 
2) Folder path to write the dicom files 
3) Defining file naming convention or pattern.  
3) Metadata information from the users

Output
------
- 16-bit DICOM files for each raw file. 

Reference
---------
Badano A, Graff CG, Badal A, Sharma D, Zeng R, Samuelson F, Glick S, Myers K, "Evaluation of Digital Breast Tomosynthesis as Replacement of Full-Field Digital Mammography Using an In Silico Imaging Trial", JAMA Netw Open. 2018;1(7):e185474. *doi:10.1001/jamanetworkopen.2018.5474*

