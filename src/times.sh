#!/bin/bash

machine=0
ComputeCapab='sm_35'

CARGS="-I /usr/local/cuda/include -L /usr/local/cuda/lib64"

for n in `seq 1 10`; do
rm -fr times$n
mkdir times$n



cd bioinformatic

nvcc ${CARGS} -arch=$ComputeCapab -Xptxas -dlcm=ca -w -O2 SubSeqMax6.15.cu -o SubSeqMax
for i in 131072 262144 524288 1048576 2097152 4194304 8388608 16777216 33554432 67108864 134217728 268435456 536870912 1073741824 2147483648 4294967296; do
    nvprof --normalized-time-unit ms  ./SubSeqMax  $i $machine 2>> ../times$n/SubSeqMax.txt
done

cd ../matMul

nvcc ${CARGS} -arch=$ComputeCapab -Xptxas -dlcm=ca -w -O2 matMul_gpu_uncoalesced.cu -o matMul-Gm-Un
for i in  128 256 512 1024 2048 4096 8192 16384 32768 65536; do

    nvprof --normalized-time-unit ms  ./matMul-Gm-Un  $i $machine 2>> ../times$n/matMul-Gm-Un-SP.txt
done

nvcc ${CARGS} -arch=$ComputeCapab -Xptxas -dlcm=ca -w -O2 matMul_gpu.cu -o matMul-Gm
for i in  128 256 512 1024 2048 4096 8192 16384 32768 65536; do

    nvprof --normalized-time-unit ms -Xptxas -dlcm=ca ./matMul-Gm  $i $machine 2>> ../times$n/matMul-Gm-SP.txt
done


nvcc ${CARGS} -arch=$ComputeCapab -Xptxas -dlcm=ca -w -O2 matMul_gpu_sharedmem_uncoalesced.cu -o matMul-Sm-Un
for i in  128 256 512 1024 2048 4096 8192 16384 32768 65536; do

    nvprof --normalized-time-unit ms  ./matMul-Sm-Un  $i $machine 2>> ../times$n/matMul-Sm-Un-SP.txt
done


nvcc ${CARGS} -arch=$ComputeCapab -Xptxas -dlcm=ca -w -O2 matMul_gpu_sharedmem.cu -o matMul-Sm
for i in  128 256 512 1024 2048 4096 8192 16384 32768 65536; do

    nvprof --normalized-time-unit ms  ./matMul-Sm  $i $machine 2>> ../times$n/matMul-Sm-SP.txt
done


cd ../

mv times$n/ experiments/


done
