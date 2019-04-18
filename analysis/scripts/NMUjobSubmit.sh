#!/bin/bash

source /home/tireman/simulation/jlab/nmu-npol/env_setup/NMUnpolVariables.sh

for ((i=$1; i<=$2; i++))
do
    export JOBNUMBER=$i
    #$BUILD_DIR/../analysis/NpolAnalysis 1>$NPOLWORKDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.out 2>$NPOLWORKDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.err
	
	$BUILD_DIR/../analysis2/NpolProcessEvents 1>$NPOLDIR/dumpFiles/${NPOLBASENAME}Analysis2_$i.out 2>$NPOLDIR/dumpFiles/${NPOLBASENAME}Analysis2_$i.err 

	#root -l -q $BUILD_DIR/../../root_macros/ProcessElectrons.cxx
	
done    

