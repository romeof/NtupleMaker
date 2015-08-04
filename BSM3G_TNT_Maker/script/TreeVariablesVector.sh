#!/bin/bash
#Specify needed variables
varType=double
varList=(Muon_tunePBestTrackType)
#Kinematics
#Muon_ptError
#Muon_dB Muon_bestTrack_pT Muon_pTerrorOVERbestTrackpT Muon_tunePBestTrack_pt
#ID
#Muon_isTrackerMuon Muon_isMediumMuon Muon_POGisGood
#Isolation
#Track related variables
#Muon_chi2LocalPosition Muon_trkKink Muon_segmentCompatibility Muon_validFraction Muon_pixelLayersWithMeasurement Muon_qualityhighPurity
#Other prop
#tunePBestTrackType
varLast=Muon_tunePBestTrackType
varCount=p
#Print info
echo " "
#Declare variables
echo -e "  vector<$varType> \c"
pos=0
for count in ${varList[@]}; 
do
  if [ "${varList[$pos]}" != "$varLast" ] 
  then
   echo -e "${varList[$pos]}, \c"
  else
   echo "${varList[$pos]};"
  fi
  let pos=pos+1
done
echo " "
echo " "
#Initialise
pos=0
for count in ${varList[@]}; 
do
  echo "  ${varList[$pos]}.push_back();"
  let pos=pos+1
done
echo " "
#Set branches
pos=0
for count in ${varList[@]}; 
do
  echo " AddBranch(&${varList[$pos]}               ,\"${varList[$pos]}\");"
  let pos=pos+1
done
echo " "
#Set clear 
pos=0
for count in ${varList[@]}; 
do
  echo " ${varList[$pos]}.clear();"
  let pos=pos+1
done
echo " "
