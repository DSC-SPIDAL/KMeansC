#!/bin/bash

wd=/N/u/sekanaya/sali/projects/flink/kmeans/data2
points=1000000
centers=500000
exp="$points"_"$centers"
p=$wd/$exp/points_LE.bin
c=$wd/$exp/centers_LE.bin
n=$points
d=2
k=$centers
m=100
t=0.00001
T=$1


:<<COMMENT
#original data
wd=$HOME/sali/git/github/esaliya/ccpp/KMeansC
c=$wd/data/1e6n_1000k/centers_LittleEndian.bin
p=$wd/data/1e6n_1000k/points_LittleEndian.bin
n=1000000
d=2
k=1000
m=1000
t=0.00001
T=$1
COMMENT

nodes=$4
hostfile=$3

tpp=$T
ppn=$2

cps=12
spn=2
cpn=$(($cps*$spn))

pe=$(($cpn/$ppn))

pat="$tpp"x"$ppn"x"$nodes"

explicitbind=$5
procbind=$6
verbose=$7

reportmpibindings=--report-bindings
#reportmpibindings=

btl="--mca btl ^openib"
if [ $procbind = "core" ]; then
    # with IB and bound to corresponding PEs
    mpirun $btl $reportmpibindings --map-by ppr:$ppn:node:PE=$pe --bind-to core -hostfile $hostfile -np $(($nodes*$ppn)) ./kmeans-fj -n$n -d$d -k$k -t$t -c$c -p$p -m$m -o"out.txt" -T$T -b$explicitbind $verbose 2>&1 | tee "$pat"_"$n"_"$k"_"$d"_"$m".txt
elif [ $procbind = "socket" ]; then
    # with IB but bound to socket
    mpirun $btl $reportmpibindings --map-by ppr:$ppn:node --bind-to socket -hostfile $hostfile -np $(($nodes*$ppn)) ./kmeans-fj -n$n -d$d -k$k -t$t -c$c -p$p -m$m -o"out.txt" -T$T -b$explicitbind $verbose 2>&1 | tee "$pat"_"$n"_"$k"_"$d"_"$m".txt
else
    # with IB but bound to none
    mpirun $btl $reportmpibindings --map-by ppr:$ppn:node --bind-to none -hostfile $hostfile -np $(($nodes*$ppn)) ./kmeans-fj -n$n -d$d -k$k -t$t -c$c -p$p -m$m -o"out.txt" -T$T -b$explicitbind $verbose 2>&1 | tee "$pat"_"$n"_"$k"_"$d"_"$m".txt
fi
