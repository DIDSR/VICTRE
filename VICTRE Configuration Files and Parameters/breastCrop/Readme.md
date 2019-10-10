Parameters used for the four breast densities in the VICTRE trial. For further details on running this code, visit https://breastcrop.readthedocs.io/en/latest/#.

**DENSE**\
```
./phantomCrop -d $dense_phantom_directory -s $phantom_seed -g 1.0 -x 810 -y 1920 -z 745
```
**HETERO**\
```
./phantomCrop -d $hetero_phantom_directory -s $phantom_seed -g 1.0 -x 1280 -y 1950 -z 940
```
**SCATTERED**\
```
./phantomCrop -d $scattered_phantom_directory -s $phantom_seed -g 1.0 -x 1740 -y 2415 -z 1140
```
**FATTY**\
```
./phantomCrop -d $fatty_phantom_directory -s $phantom_seed -g 1.0 -x 2250 -y 2760 -z 1240
```
