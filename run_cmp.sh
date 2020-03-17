#!/bin/bash
nvidia-smi --query-gpu=index,name,utilization.memory --format=csv
echo -n "Please select a GPU: "
read GPU_ID
export CUDA_VISIBLE_DEVICES=$GPU_ID
module load MATLAB/9.4
module load SuiteSparse/5.1.2-foss-2018b-METIS-5.1.0
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/datagrid/personal/lucivpav/gflags-2.2.2/build/lib:/datagrid/personal/lucivpav/InLocCIIRC_demo/functions/vlfeat/toolbox/mex/mexa64
cat startup.m inloc_demo.m | matlab -nodesktop 2>&1 | tee InLocCIIRC_demo.log