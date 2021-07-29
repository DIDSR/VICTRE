This folder contains phantom lesion locations in voxels and corresponding DBT locations.  The digital mammography ROI locations are available at https://github.com/DIDSR/VICTRE_DM_ROIs. Gzip files can be extracted using WinZip on Windows or from the command line in Linux with $tar -xvzf gzip-filename.

Phantom lesion locations files are named as pcl_PhantomSeed_crop.loc. It contains four columns (X Y Z lesionFlag), where X Y Z are voxel location in a compressed phantom and lesionFlag indicates type of lesion (0 - microcalcification, 1 - spiculated mass).

DBT location files are named as roi_SP/SA_PhantomSeed.loc. The DBT*.tar.gz contains both signal absent and present locations. The location files have four columns (X Y Z lesionNumber), where X Y Z are location in the DBT volume and lesion number ranges from 0-7 for signal present (0-3 calcification, 4-7 mass) and 0-11 for signal absent.

