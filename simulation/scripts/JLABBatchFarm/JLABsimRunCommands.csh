#!/bin/tsch

setenv BUILD_DIR /home/tireman/simulation/e12-17-004/halla-rpgen/build/simulation
setenv pType $2
source $BUILD_DIR/scripts/JLABBatchFarm/JLABsetupRun.csh

setenv JOBNUMBER $1

cp -R $BUILD_DIR/gdml .
cp -R $BUILD_DIR/macros .
cp -R $BUILD_DIR/../../simulation/include .
cp -R $BUILD_DIR/../../npollib/include/*.hh .
cp $BUILD_DIR/../npollib/libNpolClasses.so .
cp $BUILD_DIR/../../simulation/include/*.hh .

source /site/12gev_phys/softenv.csh 2.2

echo "Starting up Job Number $1."	

$BUILD_DIR/Npolapp $BUILD_DIR/macros/$pType\.mac  
source $BUILD_DIR/scripts/JLABBatchFarm/JLABAnalysisRunCommands.csh $1 $2

