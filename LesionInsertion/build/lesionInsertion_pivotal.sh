#################################################################

				LESION INSERTION SCRIPT

AUTHOR: 	DIKSHA SHARMA
		DIKSHA.SHARMA@FDA.HHS.GOV

#################################################################

#################################################################

					DISCLAIMER

This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other
parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified. 

#################################################################

## THIS SCRIPT INSERTS LESIONS IN THE FOUR CATEGORIES OF PHANTOMS (DENSE, HETERO, SCATTERED, FATTY) - DISTRIBUTING THE NUMBER OF PHANTOMS ON THE AVAILABLE CPU NODES.


#!/bin/bash

# start index from arg
numNodes="24" # number of CPU nodes to run on


# get list of phantom seeds in fatty category
seeds=($(cat /raidb/VICTRE/pivotal/phantoms/withLesion/fatty/seeds.txt))


for((i=$1; i<${#seeds[@]}; i+=$numNodes)) do
	thisSeed=${seeds[$i]}

	temp1=`sed -n 's/^ *ElementSpacing *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/fatty/pc_"$thisSeed"_crop.mhd"`
	var1=$(echo $temp1 | cut -d ' ' -f1)

	temp2=`sed -n 's/^ *Offset *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/fatty/pc_"$thisSeed"_crop.mhd"`
	var2=$(echo $temp2 | cut -d ' ' -f1)
	var3=$(echo $temp2 | cut -d ' ' -f2)
	var4=$(echo $temp2 | cut -d ' ' -f3)

	temp3=`sed -n 's/^ *DimSize *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/fatty/pc_"$thisSeed"_crop.mhd"`
	var5=$(echo $temp3 | cut -d ' ' -f1)
	var6=$(echo $temp3 | cut -d ' ' -f2)
	var7=$(echo $temp3 | cut -d ' ' -f3)

	lesinsertCmd="time python /raidb/VICTRE/pivotal/phantoms/code/LesionInsertion/lesionInsertion_pivotal_fatty_gzip_systemcall.py $thisSeed 0 69 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115"
	$lesinsertCmd

done


# get list of seeds in scattered category
seeds=($(cat /raidb/VICTRE/pivotal/phantoms/withLesion/scattered/seeds.txt))


for((i=$1; i<${#seeds[@]}; i+=$numNodes)) do
        thisSeed=${seeds[$i]}

        temp1=`sed -n 's/^ *ElementSpacing *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/scattered/pc_"$thisSeed"_crop.mhd"`
        var1=$(echo $temp1 | cut -d ' ' -f1)

        temp2=`sed -n 's/^ *Offset *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/scattered/pc_"$thisSeed"_crop.mhd"`
        var2=$(echo $temp2 | cut -d ' ' -f1)
        var3=$(echo $temp2 | cut -d ' ' -f2)
        var4=$(echo $temp2 | cut -d ' ' -f3)

        temp3=`sed -n 's/^ *DimSize *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/scattered/pc_"$thisSeed"_crop.mhd"`
        var5=$(echo $temp3 | cut -d ' ' -f1)
        var6=$(echo $temp3 | cut -d ' ' -f2)
        var7=$(echo $temp3 | cut -d ' ' -f3)

        lesinsertCmd="time python /raidb/VICTRE/pivotal/phantoms/code/LesionInsertion/lesionInsertion_pivotal_scattered_gzip_systemcall.py $thisSeed 0 60.25 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115"
        $lesinsertCmd
done



# get list of seeds in hetero category
seeds=($(cat /raidb/VICTRE/pivotal/phantoms/withLesion/hetero/seeds.txt))


for((i=$1; i<${#seeds[@]}; i+=$numNodes)) do
        thisSeed=${seeds[$i]}

        temp1=`sed -n 's/^ *ElementSpacing *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/hetero/pc_"$thisSeed"_crop.mhd"`
        var1=$(echo $temp1 | cut -d ' ' -f1)

        temp2=`sed -n 's/^ *Offset *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/hetero/pc_"$thisSeed"_crop.mhd"`
        var2=$(echo $temp2 | cut -d ' ' -f1)
        var3=$(echo $temp2 | cut -d ' ' -f2)
        var4=$(echo $temp2 | cut -d ' ' -f3)

        temp3=`sed -n 's/^ *DimSize *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/hetero/pc_"$thisSeed"_crop.mhd"`
        var5=$(echo $temp3 | cut -d ' ' -f1)
        var6=$(echo $temp3 | cut -d ' ' -f2)
        var7=$(echo $temp3 | cut -d ' ' -f3)

        lesinsertCmd="time python /raidb/VICTRE/pivotal/phantoms/code/LesionInsertion/lesionInsertion_pivotal_hetero_gzip_systemcall.py $thisSeed 0 49 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115"
        $lesinsertCmd

done




# get list of seeds in dense category
#seeds=($(cat /raidb/VICTRE/pivotal/phantoms/withLesion/dense/seeds.txt))


#for((i=$1; i<${#seeds[@]}; i+=$numNodes)) do
#        thisSeed=${seeds[$i]}

#        temp1=`sed -n 's/^ *ElementSpacing *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/dense/pc_"$thisSeed"_crop.mhd"`
#        var1=$(echo $temp1 | cut -d ' ' -f1)

#        temp2=`sed -n 's/^ *Offset *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/dense/pc_"$thisSeed"_crop.mhd"`
#        var2=$(echo $temp2 | cut -d ' ' -f1)
#        var3=$(echo $temp2 | cut -d ' ' -f2)
#        var4=$(echo $temp2 | cut -d ' ' -f3)

#        temp3=`sed -n 's/^ *DimSize *= *//p' "/raidb/VICTRE/pivotal/phantoms/withLesion/dense/pc_"$thisSeed"_crop.mhd"`
#        var5=$(echo $temp3 | cut -d ' ' -f1)
#        var6=$(echo $temp3 | cut -d ' ' -f2)
#        var7=$(echo $temp3 | cut -d ' ' -f3)

#        lesinsertCmd="time python /raidb/VICTRE/pivotal/phantoms/code/LesionInsertion/lesionInsertion_pivotal_dense.py $thisSeed 0 48.25 630 $var1 $var2 $var3 $var4 0 0 0 $var5 $var6 $var7 100 166 115"
#        $lesinsertCmd
#done
