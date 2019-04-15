#!/bin/csh

foreach i (`seq $1 1 $2`)

cat > jsubscan_$i << EOF1
PROJECT: cgen
TRACK : simulation
OS : centos7
JOBNAME : SCAN-RP-GEN_4.4GeV_2Bdl_0cm_$i
MAIL: tireman@jlab.org
TIME: 60
MEMORY: 1000 MB
COMMAND : source JLABscanRunCommands.csh $i $3
OTHER_FILES : /u/home/tireman/simulation/e12-17-004/halla-rpgen/analysis/scripts/JLABscanRunCommands.csh
EOF1
end

foreach j (`seq $1 1 $2`)

	jsub jsubscan_$j
	rm jsubscan_$j
end

