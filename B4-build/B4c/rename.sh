#!/bin/bash

FILES=/home/iwsatlas1/emberger/Geant4/Data/ForwardCalor/rmsx_over_E_Copper2mm_absofirst_10x10mmGran_analysis/*.root

for f in $FILES
do

mv "$f" "${f/.root/_.root}"

done
