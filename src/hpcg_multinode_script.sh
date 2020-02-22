#!/bin/bash
MPI=mpirun
MAPBY=host
HOSTFILE=/home/arizvi/HPCG/src/myhosts
JULIA_SCRIPT=/home/arizvi/julia-1.3.1/bin/julia
JULIA_DRIVER=/home/arizvi/HPCG/src/interim_driver.jl

for i in 14 28 56 112 224
do
	$MPI -map-by $MAPBY \
		-hostfile $HOSTFILE \
		-np $i \
		$JULIA_SCRIPT $JULIA_DRIVER \
		--np 14 \
		--nx 16 \
		--ny 16 \
		--nz 16 \
		--rt 200 \
		--npx 0 \
		--npy 0 \
		--npz 0 \
		--zl 0 \
		--zu 0 \
		--pz 0 
	if [[ "$?" -eq 0 ]]; then
		echo "OK: $i processes per node"
	else
		echo "FAIL: Could not run $i proc experiment"
		exit 2
	fi
done
