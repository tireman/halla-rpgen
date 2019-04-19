#!/bin/tsch

setenv BUILD_DIR /home/tireman/simulation/e12-17-004/halla-rpgen/build/simulation
setenv pType $2

source $BUILD_DIR/../../env_setup/JLABnpolSetup.csh

cp -R $BUILD_DIR/gdml .
cp -R $BUILD_DIR/macros .
cp -R $BUILD_DIR/../../simulation/include .
cp -R $BUILD_DIR/../../npollib/include/*.hh .
cp $BUILD_DIR/../npollib/libNpolClasses.so .
cp $BUILD_DIR/../../simulation/include/*.hh .

source /site/12gev_phys/softenv.csh 2.2

@ NUM1 = ( $1 - 1 ) * 25 + 1
@ NUM2 = $1 * 25

foreach j (`seq $NUM1 1 $NUM2`)
  
  setenv JOBNUMBER $j

  echo "    Processing Job Number $j"

  if ( -e $NPOLDIR/root/$NPOLBASENAME\_$1.root ) then
	  
	  # This runs the NPOL efficiency code
	  $BUILD_DIR/../analysis/NpolAnalysis
	  
	  # The next 3 run the NPOL Process events code for neutral/charged particles types
	  setenv CHARGE_TYPE All
	  setenv HistoOutputDir $NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	  $BUILD_DIR/../analysis2/NpolProcessEvents
	  
	  setenv CHARGE_TYPE Neutral
	  setenv HistoOutputDir $NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	  $BUILD_DIR/../analysis2/NpolProcessEvents
	  
	  setenv CHARGE_TYPE Charged
	  setenv HistoOutputDir $NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
	  $BUILD_DIR/../analysis2/NpolProcessEvents
	  
	  # Uncomment this line if you DO NOT want to keep the raw ROOT file after run #
	  #rm $NPOLDIR/root/$NPOLBASENAME\_$1.root  
	  
	  endif 
	 
	 echo "Done processing job $j"
	 
end
