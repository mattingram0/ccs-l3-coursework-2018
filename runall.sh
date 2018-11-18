#!/bin/bash
for mat in mat/*.matrix
do
	[ -e "$mat" ] || continue
	matName=$(echo $mat | cut -d'/' -f2 | cut -d'.' -f1)
	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileBase.slurm
	sbatch -n 1 -p test.q -o ./output/${matName}_X_${matName}.txt profileBase.slurm
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
	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileBase.slurm
	sbatch -n 1 -p test.q -o ./output/${massMatName}_X_${laplaceMatName}.txt profileBase.slurm
	sed -i '$ d' profileBase.slurm 
done
