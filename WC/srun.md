# the following command request a p100 node in Wilson CLuster 
## (ignore the warning: bash: /nashome/w/wus/.bashrc: Permission denied)

[Go to official website for a full list of the nodes](https://computing.fnal.gov/wilsoncluster/hardware/)

srun --pty --nodes=1 --partition gpu_gce --gres=gpu:1 --nodelist wcgpu02 bash

**or**

[recommanded config]

srun --unbuffered --pty -A nova --partition=gpu_gce \
     --time=08:00:00 \
     --nodes=1 --ntasks-per-node=1 --gres=gpu:1 \
     --nodelist wcgpu06 /bin/bash

module load singularity

export SINGULARITY_CACHEDIR=/scratch/.singularity/cache


mkdir /scratch/work

singularity shell --userns --nv \
    --workdir=/scratch/work \
    --home=/work1/nova/wus/ \
    /work1/nova/singularity/scratch/singularity-ML-tf1.12-20191126/

# the mappped directory will store the change in container i.e. .bash_history, .bash_rc 

lstm_ee env:

Keras 2.2.4
TensorFlow 1.12
Cuda 9.0

Keras 2.2.5
TensorFlow 1.14
Cuda 10.0


