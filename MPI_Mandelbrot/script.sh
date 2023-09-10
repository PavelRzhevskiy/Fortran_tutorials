#!/bin/bash

np_arr=(1 2 4 8 16)
exe=release
output="program_output.dat"
echo "Program output" > ${output}

for ((i=0; i<${#np_arr[@]}; i++))
do
    echo "np = " ${np_arr[$i]} >> ${output}
    mpirun --oversubscribe -np ${np_arr[$i]} ${exe} >> ${output}
    echo "" >> ${output}
    sleep 5
done
