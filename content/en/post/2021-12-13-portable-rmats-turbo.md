---
title: "Portable Rmats Turbo: Part 1"
author: "Chun-Jie Liu"
date: "2021-12-13"
---

> A portable Rmats turbo should be able to run on single machine, high-performance computing, and cloud computing.

The `portable rmats-turbo` was compiled into a container, it can be run on the local machine with `Docker`. However, on the HPC and cloud are not allowed to install `Docker`. I used `Singularity` to substitute `Docker`, and it can be run on the HPC and cloud. Next, `nextflow` is used to run the `portable rmats-turbo` on the cloud.

# Dockerfile for rmats-turbo

```
FROM debian:buster
LABEL maintainer="Chun-Jie Liu <chunjie.sam.liu@gmail.com>"

RUN apt-get update && apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common \
    cmake \
    curl \
    cython \
    g++ \
    gfortran \
    git \
    libblas-dev \
    libgsl-dev \
    liblapack-dev \
    make \
    python-dev \
    python-numpy \
    r-base \
    r-cran-nloptr \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    # Use a build dir to be removed after artifacts are extracted
    && mkdir /rmats_build \
    # && echo "140.82.114.4 github.com" >> /etc/hosts \
    && cd /rmats_build \
    && git clone https://github.com/Xinglab/rmats-turbo.git \
    && cd rmats-turbo \
    # The build will source setup_environment.sh which will source ~/.bashrc.
    # Skip that by truncating setup_environment.sh
    && echo '' > setup_environment.sh \
    && ./build_rmats \
    # Copy the build results
    && mkdir /rmats \
    && cd /rmats \
    && cp /rmats_build/rmats-turbo/rmats.py ./ \
    && cp /rmats_build/rmats-turbo/cp_with_prefix.py ./ \
    && cp /rmats_build/rmats-turbo/*.so ./ \
    && mkdir rMATS_C \
    && cp /rmats_build/rmats-turbo/rMATS_C/rMATSexe ./rMATS_C \
    && mkdir rMATS_P \
    && cp /rmats_build/rmats-turbo/rMATS_P/*.py ./rMATS_P \
    && mkdir rMATS_R \
    && cp /rmats_build/rmats-turbo/rMATS_R/*.R ./rMATS_R \
    # Remove build dir
    && rm -rf /rmats_build \
    && chmod 777 /rmats/rmats.py \
    && echo 'export PATH=/rmats/:${PATH}' >> /etc/bash.bashrc \
    && cp -r /rmats/* /usr/local/bin \
    # Build STAR
    && mkdir /star_build \
    && cd /star_build \
    && curl -L -O https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz \
    && tar -xvf 2.7.9a.tar.gz \
    && cd STAR-2.7.9a/source \
    && make STAR \
    && cp STAR /usr/local/bin

# Set defaults for running the image
WORKDIR /rmats
ENTRYPOINT ["python", "/rmats/rmats.py"]
CMD ["--help"]
```

# References GRCh38

1. Download the reference genome fasta files and gtf files from ensembl.

```
chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)

destdir=refdata/human-genome-ensembl-release-104

echo "Downloading human genome reference from Ensembl..."
for chromosome in ${chrs[@]}; do
    echo "Downloading chromosome ${chromosome}..."
    cmd="nohup curl -s -o ${destdir}/chr${chromosome}.fa.gz ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chromosome}.fa.gz 1>${destdir}/${chromosome}.nohup.out 2>&1 &"
    echo ${cmd}
    eval ${cmd}
done


chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
destdir=refdata/human-genome-ensembl-release-104

for chromosome in ${chrs[@]}; do
    echo "Processing chromosome ${chromosome}..."
    zcat ${destdir}/chr${chromosome}.fa.gz >> ${destdir}/Homo_sapiens.GRCh38.104.fa
done
```

2. Build STAR index

```
nohup docker run -v /workspace/liucj/refdata/star-genome-index-grch38:/refdata \
  chunjiesamliu/rmats-turbo:latest \
  STAR --runThreadN 80 \
  --runMode genomeGenerate \
  --genomeDir /refdata \
  --genomeFastaFiles /refdata/Homo_sapiens.GRCh38.104.fa \
  --sjdbGTFfile /refdata/Homo_sapiens.GRCh38.104.gtf \
  --sjdbOverhang 100 &
```