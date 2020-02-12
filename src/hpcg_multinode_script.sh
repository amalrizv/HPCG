#!/bin/bash
mpirun -map-by node --hostfile myhosts -np 2 /home/arizvi/julia-1.3.1/bin/julia -L header.jl --np 2 --nx 16 --ny 16 --nz 16 --rt 200 --npx 0 --npy 0 --npz 0 --zl 16 --zu 0 --pz 0 
