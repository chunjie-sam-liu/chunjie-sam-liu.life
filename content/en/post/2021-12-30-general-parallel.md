---
title: "General Parallel on Cluster"
author: "Chun-Jie Liu"
date: "2021-12-30"
---

On a single machine, I can run a task in parallel with mulitple processes. I set up a simple script to run process in parallel with [`generalParallel`](https://github.com/chunjie-sam-liu/useful-scripts/blob/master/generalParallel). This script take advantage of fifo to pipe in pipe out the process. For example, if you have 1000 tasks, you can run them in parallel with 100 processes. After one process finish, the next process will start. The nuimber of running processes is always 100 until the last 100 processes.

Now, I'am working on a cluster with `Slurm` system. I create a [`generalParallelSlurm`](https://github.com/chunjie-sam-liu/useful-scripts/blob/master/generalParallelSlurm) script to submit jobs with in a similar way as single machine did. It's easy to use just like `generalParallel` script.

```
#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-12-30 12:44:36
# @DESCRIPTION:

# Number of input parameters


params=$#
script=${1}

# Set default threads as 20
ntasks=${2:-10}
nodes=${3:-2}
ntasks_per_node=${4:-5}


function usage {
  Description="Notice: The script aimed at running shell command in general parallel."
  Usage="Uage:  generalParallel shell_script.sh ntasks(default:20) nodes(default:5) ntasks_per_node(default:1)"
  ErrorNo="Error:  Number of arguments must be 1 or 2"
  ErrorScript="Error:  suffix of script file must be sh."
  if [ "$params" -lt 1 ]; then
    echo $ErrorNo
    echo $Description
    echo $Usage
    exit 1
  fi
  if [ "${script#*.}" != "sh" ]; then
    echo $ErrorScript
    echo $Description
    echo $Usage
    exit 1
  fi
  if [ $ntasks -gt 20 ]; then
    echo "Warnings: Your THREADS is $threads"
    echo "Warnings: Be careful! Your THREADS EXCEED 20!!!"
  fi
}

function run {
sbatch << EOF
#!/usr/bin/env bash
#SBATCH --job-name=${name}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks=${ntasks}
#SBATCH --ntasks-per-node=${ntasks_per_node}
#SBATCH --output=${HOME}/tmp/errout/${name}.%j.out
#SBATCH --error=${HOME}/tmp/errout/${name}.%j.err
#SBATCH --mail-user=chunjie.sam.liu@gmail.com

while read line;
do
  cmd="srun --exclusive --nodes 1 --ntasks 1 \${line} &"
  # echo \${cmd}
  eval \${cmd}
done < $script
wait
EOF
}

usage
name=`basename ${script%%.sh}`
echo "Notice - Your are running $script"
echo "Notice - Total number jobs `wc -l $script`"
run
```