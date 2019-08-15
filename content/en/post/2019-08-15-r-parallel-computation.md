---
title: "R parallel computation"
author: "Chun-Jie Liu"
date: "2019-08-15"
---

> For detailed parallelization, referring to [Beyond single core R](https://github.com/ljdursi/beyond-single-core-R)

For me, not majored in CS, alwasy hard to comprehend some terms of hardware and software.

## Hardware terms

- Node: A single motherboard with multiple sockets.
- Processor/Socket: the silicon containing multiple cores.
- Core: the unit of computation.
- Pseudo-cores: can appear to the OS as multiple cores.

![Sockets cores](./img/r-parallel/sockets-cores.png)

## Software terms

Interpreted languages generally you can only directly work with processes. You can call libraries that invoke threads (BLAS/LAPACK)

- Process: data and code in memory.
- One or more threads of execution within a process.
- Threads in the same process can see most of the same memory.
- Processes generally can not peer into another processes memory.

## The R parallel packages

- `parallel`: incorporates two pacakges below.
- `multicore`: **forking**, for using all processors on a single processor.
- `snow`: **spawning**, for using any group of processors, possibly across a cluster.
