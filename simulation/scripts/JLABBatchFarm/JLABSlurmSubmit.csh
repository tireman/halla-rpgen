#!/bin/csh

#SBATCH -N1 -t3:00:00
#SBATCH -A cgen
#SBATCH -p simulation
OS : centos7
JOBNAME : RP-GEN_SIM_$3_$i
MAIL: tireman@jlab.org
TIME: 180
MEMORY: 1400 MB
COMMAND : source JLABsimRunCommands.csh $i $3
OTHER_FILES : /u/home/tireman/simulation/e12-17-004/halla-rpgen/build/simulation/scripts/JLABBatchFarm/JLABsimRunCommands.csh


