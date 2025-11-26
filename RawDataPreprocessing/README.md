## fast5 to pod5 conversion:
**pod5 installation steps:**
    - Create new environment with python version >= 3.8 
        - `mamba create -n pod5 python=3.8`
    - Activate env 
    - `python3 pip install --user pod5`

## Dorado Modified Basecalling:
    - model download: `dorado-0.7.3-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0`
        - moved `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` into `basecall_models` dir 


## GPU in ORCD cluster Information: 
    - `sinfo -o "%P %G %D %N"`
        - partitions and whether they have GPUs 
    - `scontrol show nodes`
        - show all nodes
        - have to manually filter out for GPU nodes
        - output file: gpu_node_info.txt
    - `scontrol show node [nodename]`
        - show info of node, e.g. node2804 (in partition mit_normal_gpu, has GPU)


ORCD GPU partitions/nodes:  
    - partition: mit_preemptable
        - gpu:4 2 node[2804,2906]
            - 4 GPUs each node
            - RealMemory=1031000
            - 2804=64 CPUs, 2906=128 CPUs
        - queue check: squeue -p mit_preemptable
    - partition: mit_normal_gpu
        - gpu:4 2 node[2804,2906]
        - same nodes as mit_preemptable partition ???
    - partition: sched_mit_hill
        - gpu:8 1 node235
            - RealMemory=128000
            - CPUTot=48
    - partition: sched_any
        - gpu:8 1 node235
            - same as sched_mit_hill gpu node

