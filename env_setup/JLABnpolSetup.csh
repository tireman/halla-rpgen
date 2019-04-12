#!/bin/tsch

setenv Lead 0
setenv Energy 4.4
setenv Bfield 2

setenv NPOLLIB_DIR $BUILD_DIR/../npollib
setenv NPOLBASENAME source$pType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl
setenv NPOLDIR /volatile/hallc/cgen/tireman/RP-GEN/bckgnd/run4
setenv NPOLWORKDIR /volatile/hallc/cgen/tireman/RP-GEN/bckgnd/run4
setenv NPOLEVENTSPERFILE 100000

setenv RawDataDir $NPOLDIR/root
setenv OutputDir $NPOLWORKDIR/Output
setenv HistoInputDir $NPOLWORKDIR/histos
setenv HistoOutputDir $NPOLWORKDIR/histos
setenv WorkInputDir $NPOLWORKDIR
setenv WorkOutputDir $NPOLWORKDIR

if ( ! -e $NPOLDIR/root ) then
	mkdir $NPOLDIR/root
endif

if ( ! -e $NPOLDIR/dumpFiles ) then
	mkdir $NPOLDIR/dumpFiles
endif

if ( ! -e $NPOLWORKDIR/Plots ) then
	mkdir $NPOLWORKDIR/Plots
endif

if ( ! -e $NPOLWORKDIR/Output ) then
	mkdir $NPOLWORKDIR/Output
	endif
	
if ( ! -e $NPOLWORKDIR/histos ) then
	mkdir $NPOLWORKDIR/histos
endif
		
if ( ! -e $NPOLWORKDIR/AllParticles/histos ) then
	mkdir $NPOLWORKDIR/AllParticles/histos 
endif

if ( ! -e $NPOLWORKDIR/ChargedParticles/histos ) then
	mkdir $NPOLWORKDIR/ChargedParticles/histos 
endif

if ( ! -e $NPOLWORKDIR/NeutralParticles/histos ) then
	mkdir $NPOLWORKDIR/NeutralParticles/histos 
endif
