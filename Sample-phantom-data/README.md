*Questions should be directed to: diksha.sharma@fda.hhs.gov.*

VICTRE sample phantom data
---------------------------

We have made available raw data for one phantom from each of the four breast density categories (dense, heterogeneously dense, scattered density, and fatty).  Due to github repository size limitations we are currently exploring alternative mechanisms for sharing additional phantoms. 

## Folder contents

Each folder contains multiple files. In the following filename descriptions *dddd* represents the density, *nnnnnnnn* represents the seed number:
- *dddd_p_nnnnnnnn.raw.gz* (Uncompressed)\
    The raw phantom volume (tissue labels) stored as 8-bit unsigned integers in a gzip archive. There is no file header. Additional information available at https://github.com/DIDSR/breastPhantom.  A lookup table for the tissue labels is available at https://breastphantom.readthedocs.io/en/latest/output.html
- *dddd_p_nnnnnnnn.mhd*\
    A metaImage header file containing information about the format of the file dddd_p_nnnnnnnn.raw.gz. Parsing this text header file will allow you to read and manipulate raw.gz phantom files. See https://itk.org/Wiki/MetaIO/Documentation for more information on the MetaImage format.
- *dddd_pc_nnnnnnnn_crop.raw.gz* (Compressed and cropped)\
    The raw compressed and cropped phantom volume (tissue labels) stored as 8-bit unsigned integers in a gzip archive.  More information available at https://github.com/DIDSR/breastCompress.
- *dddd_pc_nnnnnnnn_crop.mhd*\
    A metaImage header file containing information about the phantom data stored in the file dddd_pc_nnnnnnnn_crop.raw.gz
- *dddd_pcl_nnnnnnnn_crop.raw.gz* (Lesion inserted)\
    The raw lesion inserted, compressed and cropped phantom volume stored as 8-bit unsigned integers in a gzip archive.  More information available at https://github.com/DIDSR/VICTRE/tree/master/Lesion%20Insertion.
- *mcgpu_image_pcl_nnnnnnnn_crop.raw.gz-dddd_0000.raw* (X-ray image)\
    The raw MCGPU projection image stored as 32-bit real, little-endian byte order, with dimensions 3000x1500 pixels.  This contains two raw images concatenated - first is the primary+scatter (Compton+Rayleigh+multiple), second image is primary only.  More information available at https://github.com/DIDSR/VICTRE_MCGPU.
    
    
## Phantom dimensions
Phantom dimensions for various densities (this information is available in the .mhd files under *DimSize*). Lesion inserted phantoms have same dimensions as compressed and cropped.

**Dense** \
1010 1791 1434 (uncompressed) \
810  1920 745 (compressed+cropped)

**Heterogeneously dense** \
1495 1791 1794 (uncompressed) \
1280 1950 940 (compressed+cropped)

**Scattered density** \
1918 2227 2168 (uncompressed) \
1740 2415 1140 (compressed+cropped)

**Fatty** \
2440 2589 2198 (uncompressed) \
2250 2760 1240 (compressed+cropped)
   
## Unzipping and merging .gz files 
Following commands can be used on Linux,

**Unzipping:** 
```
$ gunzip file.gz
```

**Merging:** 
Due to the Github maximum file size limit, some of the phantom raw.gz files were split and need to be merged before unzipping. The split files end with file.gz.?. 
```
$ cat file.gz.? > file.gz
```
