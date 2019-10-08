
https://github.com/DIDSR/VICTRE/tree/master/ROI%20Extraction/build.

Focal spot parameters for the four breast densities:

**DENSE**\ 
$FX_x $FX_y $FX_z -> 0 48.25 630

**HETERO**\
$FX_x $FX_y $FX_z -> 0 49 630

**SCATTERED**\
$FX_x $FX_y $FX_z -> 0 60.25 630

**FATTY**\
$FX_x $FX_y $FX_z -> 0 69 630

**Extract the following information from the pc_$phantomSeed_crop.mhd**\
var1 -> ElementSpacing\
var2, var3, var4 -> Offset\
var5, var6, var7 -> DimSize

**Mammography Signal Present**\
time python roiExtraction_mammo_SP.py $phantom_seed $FX_x $FX_y $FX_z 0 $var1 0.085 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 65 109 1500 3000 -20

**Mammography Signal Absent**\
time python roiExtraction_mammo_SA.py $phantom_seed $FX_x $FX_y $FX_z 0 $var1 0.085 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 65 109 12 1500 3000 -20 115

**DBT Signal Present**\
time python roiExtraction_DBT_FBP_SP.py $thisSeed 0 60.25 630 0 1421 1024 57 0.085 1 $var2 $var3 $var4 0 0 0 65 109 5 9 2 2 4 4 $var1

**DBT Signal Absent**\
