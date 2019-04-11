#!/bin/bash

if [ $# -eq 0 ]
then
	START=1
	FINISH=1
elif [ $# -eq 1 ]
then
	START=1
	FINISH=$1
else
	START=$1
	FINISH=$(($1 + $2 - 1))
fi

for ((i=$START; i<=$FINISH; i++))
do
    export JOBNUMBER=$i
	export PType=QENeutron
	source /home/tireman/simulation/jlab/nmu-rpgen/env_setup/NMUnpolVariables.sh
	
    if [ -f "$NPOLDIR/${NPOLBASENAME}_$i.out" ]
    then
		rm $NPOLDIR/${NPOLBASENAME}_$i.out
		rm $NPOLDIR/${NPOLBASENAME}_$i.err
    fi

	echo "Starting up Job Number $i."
    $BUILD_DIR/Npolapp $BUILD_DIR/macros/RunEventGenerator.mac 1>$NPOLDIR/dumpFiles/${NPOLBASENAME}_$i.out 2>$NPOLDIR/dumpFiles/${NPOLBASENAME}_$i.err 

	hadd -f -k $NPOLDIR/root/source$PType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl_$i.root $NPOLDIR/root/source$PType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl_$i\_*.root

	rm $NPOLDIR/root/source$PType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl_$i\_*.root

	$BUILD_DIR/../analysis/NpolAnalysis 1>$NPOLDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.out 2>$NPOLDIR/dumpFiles/${NPOLBASENAME}Analysis1_$i.err

	rm $NPOLDIR/root/source$PType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl_$i.root
	
done
