VICTRE sample phantom data
---------------------------

We have made available raw data for one phantom from each of the four breast density categories (dense, heterogeneously dense, scattered desnity, and fatty). 

Each folder contains multiple files. In the following filename descriptions *dddd* represents the density, *nnnnnnnn* represents the seed number:
- *dddd_p_nnnnnnnn.raw.gz* (Uncompressed)\
    The raw phantom volume stored as 8-bit unsigned integers in a gzip archive. There is no file header. More informationa available at https://github.com/DIDSR/breastPhantom.
- *dddd_p_nnnnnnnn.mhd*\
    A metaImage header file containing information about the phantom data stored in the file p_nnnnnnnn.raw.gz. Parsing this header file will allow you to read and manipulate raw.gz phantom files. See https://itk.org/Wiki/MetaIO/Documentation for more information on the MetaImage format.
- *dddd_pc_nnnnnnnn_crop.raw.gz* (Compressed and cropped)\
    The raw compressed and cropped phantom volume stored as 8-bit unsigned integers in a gzip archive.  More information available at https://github.com/DIDSR/breastCompress.
- *dddd_pc_nnnnnnnn_crop.mhd*\
    A metaImage header file containing information about the phantom data stored in the file pc_nnnnnnnn_crop.raw.gz
- *dddd_pcl_nnnnnnnn_crop.raw.gz* (Lesion inserted)\
    The raw lesion inserted, compressed and cropped phantom volume stored as 8-bit unsigned integers in a gzip archive.  More information available at https://github.com/DIDSR/VICTRE/tree/master/Lesion%20Insertion.
- *mcgpu_image_pcl_nnnnnnnn_crop.raw.gz-dddd_0000.raw* (X-ray image)\
    The raw MCGPU projection image stored as 32-bit real, little-endian byte order.  More information available at https://github.com/DIDSR/VICTRE_MCGPU.
    
    
    Link to breast phantom/compress
    phantom dimensions for various densities
    MCGPU image format - 5 images (link to mcgpu)
    merging split files
    Viewing the data
