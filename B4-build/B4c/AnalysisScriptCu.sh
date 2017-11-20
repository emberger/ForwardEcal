#!/bin/bash


Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor

FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1.75mm_absfirst_25x25mmGran/*.root

mkdir $Path/rmsx_over_E_Pb1.75mm_absfirst_25x25mmGran_analysis/

AnaPath=$Path/rmsx_over_E_Pb1.75mm_absfirst_25x25mmGran_analysis/


counter=0
for f in $FILES
do
    #echo "$f"
    FilePath=$f
    filename="${FilePath##*/}"

    foldername=${filename%_*}

    # mkdir $AnaPath

    ./Analysis 1000 $FilePath $AnaPath $foldername


done

# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_20x20mmGran/*.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_gapfirst_20x20mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_gapfirst_20x20mmGran_analysis/
#
#
# counter=0
# for f in $FILES
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#     #mkdir $AnaPath
#
#     ./Analysis 1000 $FilePath $AnaPath $foldername
#
#
# done
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_40x40mmGran/*.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_gapfirst_40x40mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_gapfirst_40x40mmGran_analysis/
#
#
# counter=0
# for f in $FILES
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#    # mkdir $AnaPath
#
#     ./Analysis 1000 $FilePath $AnaPath $foldername
#
#
# done
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_50x50mmGran/*.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_gapfirst_50x50mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_gapfirst_50x50mmGran_analysis/
#
#
# counter=0
# for f in $FILES
# do
#     #echo "$f"
#     FilePath=$f
#     filename="${FilePath##*/}"
#
#     foldername=${filename%_*}
#
#     # mkdir $AnaPath
#
#     ./Analysis 1000 $FilePath $AnaPath $foldername
#
#
# done
