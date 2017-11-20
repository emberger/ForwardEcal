#!/bin/bash


Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor

FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_absfirst_20x20mmGran_100layers/*.root

mkdir $Path/rmsx_over_E_Cu2mm_absfirst_20x20mmGran_100layers_analysis/

AnaPath=$Path/rmsx_over_E_Cu2mm_absfirst_20x20mmGran_100layers_analysis/


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

Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/

FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_absfirst_20x20mmGran_100layers/*.root

mkdir $Path/rmsx_over_E_Pb1mm_absfirst_20x20mmGran_100layers_analysis/

AnaPath=$Path/rmsx_over_E_Pb1mm_absfirst_20x20mmGran_100layers_analysis/


counter=0
for f in $FILES
do
    #echo "$f"
    FilePath=$f
    filename="${FilePath##*/}"

    foldername=${filename%_*}

    #mkdir $AnaPath

    ./Analysis 1000 $FilePath $AnaPath $foldername


done
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_absfirst_40x40mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_absfirst_40x40mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_absfirst_40x40mmGran_analysis/
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_absfirst_50x50mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_absfirst_50x50mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_absfirst_50x50mmGran_analysis/
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
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_gapfirst_10x10mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_gapfirst_10x10mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_gapfirst_10x10mmGran_analysis/
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
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_gapfirst_20x20mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_gapfirst_20x20mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_gapfirst_20x20mmGran_analysis/
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_gapfirst_40x40mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_gapfirst_40x40mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_gapfirst_40x40mmGran_analysis/
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Pb1mm_gapfirst_50x50mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Pb1mm_gapfirst_50x50mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Pb1mm_gapfirst_50x50mmGran_analysis/
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
#
#
#
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_absfirst_10x10mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_absfirst_10x10mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_absfirst_10x10mmGran_analysis/
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
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_absfirst_20x20mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_absfirst_20x20mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_absfirst_20x20mmGran_analysis/
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_absfirst_40x40mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_absfirst_40x40mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_absfirst_40x40mmGran_analysis/
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_absfirst_50x50mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_absfirst_50x50mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_absfirst_50x50mmGran_analysis/
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
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_10x10mmGran/gamma1500MeV_0.4:0.4_.root
#
# mkdir $Path/rmsx_over_E_Cu2mm_gapfirst_10x10mmGran_analysis/
#
# AnaPath=$Path/rmsx_over_E_Cu2mm_gapfirst_10x10mmGran_analysis/
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
#
# Path=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/
#
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_20x20mmGran/gamma1500MeV_0.4:0.4_.root
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_40x40mmGran/gamma1500MeV_0.4:0.4_.root
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
# FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Cu2mm_gapfirst_50x50mmGran/gamma1500MeV_0.4:0.4_.root
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
