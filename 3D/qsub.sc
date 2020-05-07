#!/bin/bash
#COBALT -n 4
#COBALT -t 1:00:00
#COBALT -q debug-cache-quad --attrs mcdram=cache:numa=quad
#COBALT -A DynamicCS

PROC_PER_NODE=8 # currently I only use 1 rank per node since we don't have MPI support yet. We could run ensample job if needed 

aprun -n $((COBALT_JOBSIZE*PROC_PER_NODE)) -N ${PROC_PER_NODE} -d 4 -j 1 -e OMP_NUM_THREADS=8 ./sim_tv -d Data/1024_au_sto.h5 -m Measurement/1024_au_sto_measurement.h5 -o Output/1024_au_sto_recon.h5 -n 40
