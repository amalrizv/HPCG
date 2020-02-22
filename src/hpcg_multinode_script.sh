#!/bin/bash
MPI=mpirun
MAPBY=node
HOSTFILE=/home/arizvi/HPCG/src/myhosts
JULIA_SCRIPT=/home/arizvi/julia-1.3.1/bin/julia
JULIA_DRIVER=/home/arizvi/HPCG/src/interim_driver.jl
BASE_NUM_HOSTS=$(wc -l $HOSTFILE | cut -f1 -d' ')

echo "Running HPCG Julia experiment for $BASE_NUM_HOSTS nodes"

for i in  1 2 4 8 16
do
	PROC_COUNT=$(( $i * $BASE_NUM_HOSTS ))
	$MPI -map-by $MAPBY -hostfile $HOSTFILE -np $PROC_COUNT $JULIA_SCRIPT $JULIA_DRIVER \
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
	if [[ "$?" -ne 0 ]]; then
		echo "FAIL: could not run $i MPI ranks per node ($PROC_COUNT total) experiment"
		exit 2
	else
		echo "OK: HPCG Julia $i MPI ranks per node ($PROC_COUNT total)"
	fi
done
