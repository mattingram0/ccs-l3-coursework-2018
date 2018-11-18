#!/bin/bash
mat1=$(echo $1 | cut -d'/' -f2 | cut -d'.' -f1)
mat2=$(echo $2 | cut -d'/' -f2 | cut -d'.' -f1)
echo "likwid-perfctr -f -C 0 -g $3 -m ./profile results/${mat1}_X_${mat2}.matrix mat/${mat1}.matrix mat/${mat2}.matrix" >> profileBase.slurm
sbatch -n 1 -p test.q -t $4 -o ./output/${mat1}_X_${mat2}.txt profileBase.slurm
sed -i '$ d' profileBase.slurm

#To Use: $1 - mat/(file.matrix), $2 - mat/(file.matrix), $3 - event
