#!/bin/sh

export Lead=0
export Energy=4.4
export Bfield=2
export pType=ElectronBeam #QENeutron
export BUILD_DIR=/home/tireman/simulation/jlab/nmu-rpgen/build/simulation # nmu-npol/build/simulation
export NPOLLIB_DIR=/home/tireman/simulation/jlab/nmu-rpgen/build/simulation  #nmu-npol/build/npollib
export NPOLBASENAME=source$pType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl


#electronBeam_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl 
#source$PType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl

NPOLDIR=/home/tireman/data1/ElectronBeam/RP-GEN/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm
export NPOLWORKDIR=/home/tireman/data1/ElectronBeam/RP-GEN/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm

#export NPOLDIR=/home/tireman/data1/NeutronBeam/FixedBeam/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2
#export NPOLWORKDIR=/home/tireman/data1/NeutronBeam/FixedBeam/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2/Test_1

#export NPOLDIR=/home/tireman/data1/TargetTaggerSource/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2
#export NPOLWORKDIR=/home/tireman/data1/TargetTaggerSource/$pType\RawSim/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2/Test_5
#NPreal_Cut

#export NPOLDIR=/home/tireman/data1/PointSource/NeutronBeam/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2
#export NPOLWORKDIR=/home/tireman/data1/PointSource/NeutronBeam/$Energy\GeV/$Bfield\Bdl/Lead$Lead\cm/Location_2

export RawDataDir=$NPOLDIR/root
export OutputDir=$NPOLDIR/Output
export InputDir=$NPOLDIR/root
export WorkOutputDir=$NPOLWORKDIR
export WorkInputDir=$NPOLWORKDIR
export HistoOutputDir=$NPOLWORKDIR/histos
export HistoInputDir=$NPOLWORKDIR/histos
export PlotsOutputDir=$NPOLWORKDIR/Plots

if [ ! -e $NPOLDIR ]
then
	mkdir $NPOLDIR
fi
if [ ! -e $NPOLWORKDIR ]
then
	mkdir $NPOLWORKDIR
fi
if [ ! -e $NPOLDIR/Output ]
then
	mkdir $NPOLDIR/Output
fi
if [ ! -e $NPOLWORKDIR/Output ]
then
	mkdir $NPOLWORKDIR/Output
fi
if [ ! -e $NPOLDIR/Plots ]
then
	mkdir $NPOLDIR/Plots
fi
if [ ! -e $NPOLWORKDIR/Plots ]
then
	mkdir $NPOLWORKDIR/Plots
fi
if [ ! -e $NPOLDIR/histos ]
then
	mkdir $NPOLDIR/histos
fi

if [ ! -e $NPOLDIR/dumpFiles ]
then
	mkdir $NPOLDIR/dumpFiles
fi

if [ ! -e $NPOLWORKDIR/dumpFiles ]
then
	mkdir $NPOLWORKDIR/dumpFiles
fi

if [ ! -e $NPOLWORKDIR/histos ]
then
	mkdir $NPOLWORKDIR/histos
fi

