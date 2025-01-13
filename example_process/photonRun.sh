#!/bin/bash -l
#PBS -l nodes=1:ppn=32
#PBS -l walltime=24:00:00
#PBS -r n
#PBS -j oe
#PBS -q greenq

cd $PBS_O_WORKDIR 
module load gcc

# Record start time
start_time=$(date +%s)

./rundata 205 32

# Record end time
end_time=$(date +%s)

# Calculate and output the total runtime
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"

mv *.dat data
