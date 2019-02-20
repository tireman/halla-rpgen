#!/bin/csh

foreach i (`seq $1 1 $2`)

cat > jsubfile_$i << EOF1
PROJECT: cgen
TRACK : simulation
OS : centos7
JOBNAME : Electron-11GeV-4Bdl_$i
MAIL: tireman@jlab.org
TIME: 900
MEMORY: 1800 MB
COMMAND : source JLABsimRunCommands.csh $i
OTHER_FILES : /u/home/tireman/simulation/e11_12_009/background/nmu-npol/build/simulation/scripts/JLABBatchFarm/JLABsimRunCommands.csh
EOF1
end

foreach j (`seq $1 1 $2`)

  jsub jsubfile_$j
  rm jsubfile_$j

end
