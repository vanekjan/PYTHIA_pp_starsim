#!/bin/bash
#submitdir=`/bin/pwd -P`
submitdir=$PWD
echo $submitdir

prodId=`date +%F_%H-%M`

nevents=10000

energy=200 #collisions energy

random=$RANDOM

config="pp:minbiasLambda"
#config="pp:dmeson"
#config="pp:dmesontune"
#config="pp:dmesontune_new"

#for output file name in xml file
tune="minbiasLambda"
#tune="dmeson"
#tune="dmesontune"
#tune="dmesontune_new"

mkdir -p ./SubmitInfo/


echo "output directories"
echo ./production/${prodId}
echo ./jobs/${prodId}

mkdir -p ./production/${prodId}
mkdir -p ./jobs/${prodId}

star-submit-template -template submitStarsimPYTHIA_8.xml -entities submitdir=$submitdir,nevents=$nevents,random=$random,productionId=${prodId},config=$config,tune=$tune,energy=$energy

