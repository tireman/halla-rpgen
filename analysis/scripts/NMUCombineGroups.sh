#!/bin/bash

export HistoOutputDir=$NPOLWORKDIR/histos
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx

export CHARGE_TYPE=All
export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
export RatesOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles
export RatesInputDir=$RatesOutputDir
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
root -l -q $BUILD_DIR/../../root_macros/RPgenRates.cxx

export CHARGE_TYPE=Neutral
export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
export RatesOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles
export RatesInputDir=$RatesOutputDir
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
root -l -q $BUILD_DIR/../../root_macros/RPgenRates.cxx

export CHARGE_TYPE=Charged
export HistoOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles/histos
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
export RatesOutputDir=$NPOLWORKDIR/$CHARGE_TYPE\Particles
export RatesInputDir=$RatesOutputDir
root -l -q $BUILD_DIR/../../root_macros/NpolCombineHistos.cxx
root -l -q $BUILD_DIR/../../root_macros/RPgenRates.cxx

	
