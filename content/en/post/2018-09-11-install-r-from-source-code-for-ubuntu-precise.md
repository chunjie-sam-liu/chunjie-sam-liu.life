---
title: Install R from source code for ubuntu precise
author: Chun-Jie Liu
date: '2018-09-11'
---

The newest CRAN does not support ubuntu precise anymore.

But some packages need to upgrade the R to > 3.4.0.

So, you have to re-isntall the r from source code, and set repo as "cloud.r-project.org", then the install.package will work by install from source code.

```
options(repos = c(CRAN = "cloud.r-project.org"))
```

First, you should uninstall the r from the server.

```
apt-get --purge remove r-base
apt-get --purge remove r-base-core
apt-get --purge remove r-base-dev
apt-get autoremove
```
Then, download R source code from https://cran.r-project.org

```
$ tar -zcvr R-3.4.0.tar.gz
# cd R-3.4.0
# ./configure
```

It occurs errors.

```
apt-get install liblzma-dev
./configure --enable-R-shlib --with-blas --with-lapack
make
make check

cp R-3.4.- /usr/local/src
ln -s ../src/R-3.4.0/bin/R .
```

Third, install rstudio-server

```
rstudio-server verify-installation
rstudio-server start
```