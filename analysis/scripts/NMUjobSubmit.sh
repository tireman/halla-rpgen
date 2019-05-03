#!/bin/bash

source /home/tireman/simulation/jlab/nmu-rpgen/env_setup/NMUnpolVariables.sh

for ((i=$1; i<=$2; i++))
do
    export JOBNUMBER=$i

	#export HistoOutputDir=$NPOLWORKDIR/histos
	#$BUILD_DIR/../analysis/NpolAnalysis 1>$NPOLWORKDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.out 2>$NPOLWORKDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.err
	
	#export CHARGE_TYPE=All
	#export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	#$BUILD_DIR/../analysis2/NpolProcessEvents

	#export CHARGE_TYPE=Neutral
	#export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	#$BUILD_DIR/../analysis2/NpolProcessEvents

	#export CHARGE_TYPE=Charged
	#export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	#$BUILD_DIR/../analysis2/NpolProcessEvents

	$BUILD_DIR/../analysis3/NpolRates

	
	#root -l -q $BUILD_DIR/../../root_macros/ProcessElectrons.cxx
	
done    

