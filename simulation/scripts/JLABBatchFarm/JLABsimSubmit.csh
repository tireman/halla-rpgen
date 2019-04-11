#!/bin/csh

foreach i (`seq $1 1 $2`)

cat > jsubfile_SIM_$3_$i << EOF1
PROJECT: cgen
TRACK : simulation
OS : centos7
JOBNAME : RP-GEN_SIM_$3_$i
MAIL: tireman@jlab.org
TIME: 180
MEMORY: 1800 MB
COMMAND : source JLABsimRunCommands.csh $i $3
OTHER_FILES : /u/home/tireman/simulation/e12-17-004/halla-rpgen/build/simulation/scripts/JLABBatchFarm/JLABsimRunCommands.csh
EOF1
end

foreach j (`seq $1 1 $2`)

  jsub jsubfile_SIM_$3_$j
  rm jsubfile_SIM_$3_$j

end

