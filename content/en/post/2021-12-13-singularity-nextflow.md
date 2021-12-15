---
title: "Portable Rmats Turbo: Part 2"
author: "Chun-Jie Liu"
date: "2021-12-13"
---

Nextflow and Singularity combined to run Rmats turbo on the HPC.

# Nextflow

The local executor is very useful for pipeline development and testing purposes, but for real world computational pipelines an HPC or cloud platform is often required.

In other words, Nextflow provides an abstraction between the pipeline’s functional logic and the underlying execution system. Thus it is possible to write a pipeline once and to seamlessly run it on your computer, a grid platform, or the cloud, without modifying it, by simply defining the target execution platform in the configuration file.

> Nextflow is based on the Dataflow programming model, which is a declarative programming model that enables the definition of a pipeline as a series of steps, each of which can be executed in parallel.

# Configure Singularity and Slurm for nextflow.config

Integrate Nextflow with rmats-turbo Singularity container, set up the Nextflow configuration.

```
params {
  alpha = 1
}

process {
  container = "/home/liuc9/sif/rmats-turbo_latest.sif"
  executor = "slurm"
}

executor {
  cpus = 1
  queue = "long"
  memory = "5G"

}

singularity {
  autoMount = true
  enabled = true
  runOptions = "--bind /scr1/users/liuc9/:/scr1/users/liuc9/,/mnt/isilon/xing_lab/liuc9/:/mnt/isilon/xing_lab/liuc9/"
}
```

# Submit jobs to the HPC with slurm

Test the process code chunk for user id in singularity container and hostname in slurm.

```
methods = ["prot", "dna", "rna"]
process bar {
    input:
    val x from methods

    output:
    stdout resultB

    """
    id
    hostname
    """
}
resultB.view { it.trim() }
```