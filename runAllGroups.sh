#!/bin/bash
for mat in mat/*.matrix
do
	[ -e "$mat" ] || continue
	matName=$(echo $mat | cut -d'/' -f2 | cut -d'.' -f1)
#
#	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${matName}_X_${matName}_DATA.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm
#	
#	echo "likwid-perfctr -f -C 0 -g MEM_SP -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${matName}_X_${matName}_MEMSP.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm
#	
#	echo "likwid-perfctr -f -C 0 -g L2CACHE -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${matName}_X_${matName}_L2.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm
#	
#	echo "likwid-perfctr -f -C 0 -g L3CACHE -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${matName}_X_${matName}_L3.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm
#
	echo "likwid-perfctr -f -C 0 -g TOP_DOWN -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
	sbatch -n 1 -p intr.q -o ./output/${matName}_X_${matName}_TOPD.txt profileBase.slurm
	sed -i '$ d' profileBase.slurm
done

for massMat in mat/mass/*.matrix
do
	[ -e "$massMat" ] || continue
	massMatName=$(echo $massMat | cut -d'/' -f3 | cut -d'.' -f1)
	
	if [[ $massMatName = *DG* ]]; then
		laplaceMatName=${massMatName//"mass"/"ip-laplace"}
	else
		laplaceMatName=${massMatName//"mass"/"laplace"}
	fi
	
#	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${massMatName}_X_${laplaceMatName}_DATA.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm 
#
#	echo "likwid-perfctr -f -C 0 -g MEM_SP -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${massMatName}_X_${laplaceMatName}_MEMSP.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm 
#
#	echo "likwid-perfctr -f -C 0 -g L2CACHE -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
#	sbatch -n 1 -p intr.q -o ./output/${massMatName}_X_${laplaceMatName}_L2.txt profileBase.slurm
#	sed -i '$ d' profileBase.slurm 

	#echo "likwid-perfctr -f -C 0 -g L3CACHE -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
	#sbatch -n 1 -p intr.q -o ./output/${massMatName}_X_${laplaceMatName}_L3.txt profileBase.slurm
	#sed -i '$ d' profileBase.slurm 

	echo "likwid-perfctr -f -C 0 -g TOP_DOWN -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
	sbatch -n 1 -p intr.q -o ./output/${massMatName}_X_${laplaceMatName}_TOPD.txt profileBase.slurm
	sed -i '$ d' profileBase.slurm 
done
