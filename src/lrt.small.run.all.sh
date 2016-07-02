#!/bin/bash

verbose=-v
#can be 1/0 to set thread pinning
explicitbind=$1
#can be core/socket/none
procbind=$2
nodes=1
name="$nodes"n
nodefile=nodes.$name.txt

#./lrt.small.run.generic.sh 1 2 $nodefile $nodes $explicitbind $procbind $verbose
./lrt.small.run.generic.sh 1 24 $nodefile $nodes $explicitbind $procbind $verbose
:<<COMMENT
./lrt.small.run.generic.sh 2 12 $nodefile $nodes $explicitbind $procbind $verbose
./lrt.small.run.generic.sh 3 8 $nodefile $nodes $explicitbind $procbind $verbose
COMMENT
#./lrt.small.run.generic.sh 4 6 $nodefile $nodes $explicitbind $procbind $verbose
:<<COMMENT
./lrt.small.run.generic.sh 6 4 $nodefile $nodes $explicitbind $procbind $verbose
./lrt.small.run.generic.sh 8 3 $nodefile $nodes $explicitbind $procbind $verbose
./lrt.small.run.generic.sh 12 2 $nodefile $nodes $explicitbind $procbind $verbose
COMMENT
#./lrt.small.run.generic.sh 24 1 $nodefile $nodes $explicitbind $procbind $verbose
