#!/bin/bash
#Specify needed variables
varType=double
object=Muon_
varList=(mujet_mindr mujet_pt muptSUmujetpt mujet_btagdisc)
varLast=nchi2_gt
varCount=p
#Kinematics
#ptError
#pTErrorSUpT
#
#dB bestTrack_pT pTerrorOVERbestTrackpT tunePBestTrack_pt
#ID
#isTrackerMuon isMediumMuon POGisGood
#Isolation
#iso_rho iso_nue_rel iso_ch_rel
#Track related variables
#chi2LocalPosition trkKink segmentCompatibility validFraction pixelLayersWithMeasurement qualityhighPurity
#dxy_it dz_it ip3d_val ip3d_err ip3d_sig nchi2_gt
#Other prop
#tunePBestTrackType
#Print info
echo " "
#Declare variables
echo -e "  vector<$varType> \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "$object${varList[$pos]}, \c"
  else
   echo "$object${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  $object${varList[$pos]}.push_back();"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo " AddBranch(&$object${varList[$pos]}               ,\"$object${varList[$pos]}\");"
  let pos=pos+1
done
echo " "
#Set clear 
pos=0
for count in ${varList[@]}; 
do
  echo " $object${varList[$pos]}.clear();"
  let pos=pos+1
done
echo " "
