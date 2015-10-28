#!/bin/bash

mpirun -np $1 ./run_knotted linear 3 neff sum_mnu planck lrg
