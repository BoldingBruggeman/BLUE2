#/usr/bin/bash

#ncpus="4 8 12 16 24 32 64"
ncpus="4 8 12 16"

for ncpu in $ncpus; do
   pygetm-subdiv optimize --legacy topo.nc $ncpu --pickle subdiv_$ncpu.pickle
done
