#!/bin/bash
mpicc -fopenmp -lm -O3 kmeans.c -o kmeans-lrt
