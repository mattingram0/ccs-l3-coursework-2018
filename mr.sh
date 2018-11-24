#!/bin/bash
make sparsemm
./sparsemm test.matrix mat/DG1-mass-2D.matrix mat/DG1-mass-2D.matrix 
