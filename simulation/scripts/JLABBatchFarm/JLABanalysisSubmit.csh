#!/bin/csh

foreach i (`seq $1 1 $2`)

cat > jsubfile_AY_$3_$i << EOF1
PROJECT: cgen
TRACK : simulation
OS : centos7
JOBNAME : TargetTagger_AY$3_$i
MAIL: tireman@jlab.org
TIME: 60
MEMORY: 1000 MB
COMMAND : source JLABAnalysisRunCommands.csh $i $3
OTHER_FILES : /u/home/tireman/simulation/e12-17-004/halla-rpgen/build/simulation/scripts/JLABBatchFarm/JLABAnalysisRunCommands.csh
EOF1
end

foreach j (`seq $1 1 $2`)

  jsub jsubfile_AY_$3_$j
  rm jsubfile_AY_$3_$j

end

