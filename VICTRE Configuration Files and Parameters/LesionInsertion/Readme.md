Parameters used for inserting lesions in the four breast density phantoms in the VICTRE trial. For further details on the various parameters, visit https://github.com/DIDSR/VICTRE/tree/master/Lesion%20Insertion/build. 


**Extract the following information from the pc_$phantomSeed_crop.mhd**\
var1 -> ElementSpacing\
var2, var3, var4 -> Offset\
var5, var6, var7 -> DimSize

**DENSE**\
time python lesionInsertion.py $phantom_seed 0 48.25 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115

**HETERO**\
time python lesionInsertion.py $phantom_seed 0 49 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115

**SCATTERED**\
time python lesionInsertion.py $phantom_seed 0 60.25 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115

**FATTY**\
time python lesionInsertion.py $phantom_seed 0 69 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115

