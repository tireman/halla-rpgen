#!/bin/tsch

setenv Lead 15
setenv Energy 4.4
setenv Bfield 4
#setenv pType Total

setenv NPOLLIB_DIR $BUILD_DIR/../npollib
setenv NPOLBASENAME source$pType\_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl 
#electronBeam_Lead$Lead\cm_$Energy\GeV_$Bfield\Bdl
setenv NPOLDIR /cache/hallc/cgen/simulation/tireman/TargetTaggerRuns/QENeutron/Lead$Lead\cm/$Energy\GeV/$Bfield\Bdl/Location_2
setenv NPOLWORKDIR /work/hallc/cgen/tireman/TargetTaggerRuns/QENeutron/Lead$Lead\cm/$Energy\GeV/$Bfield\Bdl/Location_2/NpolEffTesting/28-Jan-2018

setenv RawDataDir $NPOLDIR/root
setenv OutputDir $NPOLWORKDIR/Output
setenv HistoInputDir $NPOLWORKDIR/histos
setenv HistoOutputDir $NPOLWORKDIR/histos
setenv WorkInputDir $NPOLWORKDIR
setenv WorkOutputDir $NPOLWORKDIR

if ( ! -e $NPOLWORKDIR/Plots ) then
	mkdir $NPOLWORKDIR/Plots
endif

if ( ! -e $NPOLWORKDIR/Output ) then
	mkdir $NPOLWORKDIR/Output
endif

if ( ! -e $NPOLWORKDIR/histos ) then
	mkdir $NPOLWORKDIR/histos

endif
