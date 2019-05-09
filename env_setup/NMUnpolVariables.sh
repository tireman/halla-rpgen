#!/bin/sh

export Lead=0
export Energy=4.4
export Bfield=2
export pType=ElectronBeam
export NPOLBASENAME=source$pType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl
export NPOLEVENTSPERFILE=10000

export BUILD_DIR=/home/tireman/simulation/jlab/nmu-rpgen/build/simulation
export NPOLLIB_DIR=/home/tireman/simulation/jlab/nmu-rpgen/build/npollib

#export NPOLDIR=/home/tireman/data1/DCStest/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2
#export NPOLWORKDIR=/home/tireman/data1/DCStest/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2

export NPOLDIR=/home2/tireman-h2/data1/volatile/hallc/cgen/tireman/RP-GEN/bckgnd/run15 #/ChargedParticles
export NPOLWORKDIR=/home2/tireman-h2/data1/volatile/hallc/cgen/tireman/RP-GEN/bckgnd/run15 #/AllParticles

#export NPOLDIR=/home/jmcmullen/data1/RP-GEN/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2/Run1
#export NPOLWORKDIR=/home/jmcmullen/data1/RP-GEN/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2/Run1

export RawDataDir=$NPOLDIR/root
export OutputDir=$NPOLDIR/Output
export InputDir=$NPOLDIR/root
export WorkOutputDir=$NPOLWORKDIR
export WorkInputDir=$NPOLWORKDIR
export HistoOutputDir=$NPOLWORKDIR/histos
export HistoInputDir=$NPOLWORKDIR/histos
export PlotsOutputDir=$NPOLWORKDIR/Plots
export RatesInputDir=$NPOLWORKDIR/rates
export RatesOutputDir=$NPOLWORKDIR/rates
export DumpFileDIR=$NPOLDIR/dumpFiles
 
if [ ! -e $NPOLDIR ]
then
	mkdir -p $NPOLDIR
fi
if [ ! -e $NPOLWORKDIR ]
then
	mkdir -p $NPOLWORKDIR
fi
if [ ! -e $NPOLDIR/root ]
then
	mkdir -p $NPOLDIR/root
fi
if [ ! -e $NPOLDIR/Output ]
then
	mkdir -p $NPOLDIR/Output
fi
if [ ! -e $NPOLWORKDIR/Output ]
then
	mkdir -p $NPOLWORKDIR/Output
fi
if [ ! -e $NPOLDIR/Plots ]
then
	mkdir -p $NPOLDIR/Plots
fi
if [ ! -e $NPOLWORKDIR/Plots ]
then
	mkdir -p $NPOLWORKDIR/Plots
fi
if [ ! -e $NPOLDIR/histos ]
then
	mkdir -p $NPOLDIR/histos
fi

if [ ! -e $NPOLDIR/dumpFiles ]
then
	mkdir -p $NPOLDIR/dumpFiles
fi

if [ ! -e $NPOLWORKDIR/dumpFiles ]
then
	mkdir -p $NPOLWORKDIR/dumpFiles
fi

if [ ! -e $NPOLWORKDIR/histos ]
then
	mkdir -p $NPOLWORKDIR/histos
fi

if [ ! -e $NPOLWORKDIR/rates ]
then
	mkdir -p $NPOLWORKDIR/rates
fi

if [ ! -e $NPOLWORKDIR/AllParticles/histos ]
then
	mkdir -p $NPOLWORKDIR/AllParticles/histos 
fi

if [ ! -e $NPOLWORKDIR/ChargedParticles/histos ]
then
	mkdir -p $NPOLWORKDIR/ChargedParticles/histos 
fi

if [ ! -e $NPOLWORKDIR/NeutralParticles/histos ]
then
	mkdir -p $NPOLWORKDIR/NeutralParticles/histos 
fi

if [ ! -e $NPOLWORKDIR/AllParticles/Output ]
then
	mkdir -p $NPOLWORKDIR/AllParticles/Output 
fi

if [ ! -e $NPOLWORKDIR/ChargedParticles/Output ]
then
	mkdir -p $NPOLWORKDIR/ChargedParticles/Output 
fi

if [ ! -e $NPOLWORKDIR/NeutralParticles/Output ]
then
	mkdir -p $NPOLWORKDIR/NeutralParticles/Output 
fi

if [ ! -e $NPOLWORKDIR/AllParticles/Plots ]
then
	mkdir -p $NPOLWORKDIR/AllParticles/Plots 
fi

if [ ! -e $NPOLWORKDIR/ChargedParticles/Plots ]
then
	mkdir -p $NPOLWORKDIR/ChargedParticles/Plots 
fi

if [ ! -e $NPOLWORKDIR/NeutralParticles/Plots ]
then
	mkdir -p $NPOLWORKDIR/NeutralParticles/Plots 
fi

if [ ! -e $NPOLWORKDIR/AllParticles/rates ]
then
	mkdir -p $NPOLWORKDIR/AllParticles/rates 
fi

if [ ! -e $NPOLWORKDIR/ChargedParticles/rates ]
then
	mkdir -p $NPOLWORKDIR/ChargedParticles/rates 
fi

if [ ! -e $NPOLWORKDIR/NeutralParticles/rates ]
then
	mkdir -p $NPOLWORKDIR/NeutralParticles/rates 
fi
