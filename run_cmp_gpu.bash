#!/bin/bash
bash nvidia-usage.sh
nvidia-smi --query-gpu=index,name,utilization.memory --format=csv
echo -n "Please select a GPU: "
read GPU_ID
export CUDA_VISIBLE_DEVICES=$GPU_ID
module load MATLAB/9.7
module load SuiteSparse/5.1.2-foss-2018b-METIS-5.1.0
module load LLVM/6.0.0-GCCcore-7.3.0 # avoid pyrender/OSMesa crash
module load CUDA/9.0.176-GCC-6.4.0-2.28
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/datagrid/personal/lucivpav/gflags-2.2.2/build/lib:"${PWD}"/functions/vlfeat/toolbox/mex/mexa64
export INLOC_HW="GPU"
echo "TODO: why does running this require ~28 GB RAM?"
cat inloc_demo.m | matlab -nodesktop -singleCompThread 2>&1 | tee /mnt/datagrid/personal/lucivpav/InLocCIIRC_dataset/logs/inloc_demo_gpu-"${INLOC_EXPERIMENT_NAME}".log