#!/bin/bash
for mat in mat/*
do
	matName=$(echo $mat | cut -d'/' -f2 | cut -d'.' -f1)
	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${matName}_X_${matName}.matrix mat/${matName}.matrix mat/${matName}.matrix" >> profileall.slurm
	sbatch -n 1 -p -o output/${matName}_X_${matName}.txt test.q

for massMat in mat/mass/*
do
	massMatName=$(echo $massMat | cut -d'/' -f2 | cut -d'.' -f1)
	if grep -q massMatName <<<"D"; then
		laplaceMatName=${massMatName//"mass"/"ip-laplace"}
	else
		laplaceMatName=${massMatName//"mass"/"laplace"}
	fi
	echo "likwid-perfctr -f -C 0 -g DATA -m ./profile results/${massMatName}_X_${laplaceMatName}.matrix mat/${massMatName}.matrix mat/${laplaceMatName}.matrix" >> profileall.slurm
	sbatch -n 1 -p -o output/${massMatName}_X_${laplaceMatName}.txt test.q
