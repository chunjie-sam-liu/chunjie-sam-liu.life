---
title: "Use of jrocker in Singularity on Linux with Slurm"
author: "Chun-Jie Liu"
date: "2021-12-11"
---

The HPC was not allowed to use Docker, to facilitate the uniformity of running environment with R and rmats-turbo, I built the running environment with Dockerfile and pushed it to Docker Hub with GitHub Actions. Then pull the image from Docker Hub and run it with Singularity.

# Use Singularity with non-privileged user

## Install Go

```
# Install go
cd ~/tools
wget https://go.dev/dl/go1.17.5.linux-amd64.tar.gz
tar -xvf go1.17.5.linux-amd64.tar.gz
echo "export PATH=${HOME}/tools/go/bin/:\${PATH}" >> ~/.bashrc
source ~/.bashrc
```

## Install Singularity from CE source

```
cd ~/tools
wget https://github.com/sylabs/singularity/releases/download/v3.9.2/singularity-ce-3.9.2.tar.gz
tar -xvf singularity-ce-3.9.2.tar.gz
cd singularity-ce-3.9.2
./mconfig --without-suid --prefix=${HOME}/tools/singularity-ce-3.9.2 && \
    make -C ./builddir && \
    make -C ./builddir install
echo "export PATH=${HOME}/tools/singularity-ce-3.9.2/bin/:\${PATH}" >> ~/.bashrc
source ~/.bashrc
```

# Use jrocker in Singularity

```
singularity pull docker://chunjiesamliu/jrocker:latest
mkdir -p run var-lib-rstudio-server

printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf

PASSWORD="🙃" singularity exec --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf rstudio_4.0.4.sif rserver --auth-none=0  --auth-pam-helper-path=pam-helper --auth-timeout-minutes=0 --auth-stay-signed-in-days=30
```

Pointing your browser to http://hostname:8787, enter your local user ID on the system as the username, and the custom password specified in the PASSWORD environment variable.

# Use jrocker in Singularity with Slurm

Slurm job script

```
#!/usr/bin/env bash
#SBATCH --time=08:00:00
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8192
#SBATCH --output=/home/%u/tmp/errout/jrocker.job.%j
# customize --output path as appropriate (to a directory readable only by the user!)

# Create temporary directory to be populated with directories to bind-mount in the container
# where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(python -c 'import tempfile;from pathlib import Path; print(tempfile.mkdtemp(prefix="jrocker_", dir=str(Path.home()/"tmp/tmp")))')

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

# Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
# libraries used by R) from spawning more threads than the number of processors
# allocated to the job.
#
# Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
# personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<END
#!/usr/bin/env bash
export OMP_NUM_THREADS=${SLURM_JOB_CPUS_PER_NODE}
export R_LIBS_USER=${HOME}/R/jrocker
exec rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export SINGULARITY_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server,/scr1/users/liuc9:/scr1/users/liuc9,/mnt/isilon/xing_lab/liuc9:/mnt/isilon/xing_lab/liuc9"

# Do not suspend idle sessions.
# Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
# https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export SINGULARITYENV_RSTUDIO_SESSION_TIMEOUT=0

export SINGULARITYENV_USER=$(id -un)
export SINGULARITYENV_PASSWORD=$(openssl rand -base64 15)

# get unused socket per https://unix.stackexchange.com/a/132524
# tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:${HOSTNAME}:${PORT} ${SINGULARITYENV_USER}@hpc

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${SINGULARITYENV_USER}
   password: ${SINGULARITYENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f ${SLURM_JOB_ID}
END

singularity exec --cleanenv ${HOME}/sif/jrocker_latest.sif \
    rserver --www-port ${PORT} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh
printf 'rserver exited' 1>&2
```