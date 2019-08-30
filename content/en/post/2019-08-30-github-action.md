---
title: "GitHub Actions"
author: "Chun-Jie Liu"
date: "2019-08-30"
---

## GitHub Actions

GitHub actions is a automate build, test and deployment product. DockerHub provides automate building by connecting GitHub repo. However, the DockerHub limit the building storage and limit building time.

## Automate build and publish to DockerHub

Use [Jrocker](https://github.com/chunjie-sam-liu/jrocker) as an exmaple to automate build and publish DockerHub. The code chunk as follows.

```YAML
on: push
name: Build & Push containers
jobs:
  dockerRegistry:
    name: Docker Registry
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@master
    - name: Docker Registry
      uses: actions/docker/login@master
      env:
        DOCKER_PASSWORD: ${{ secrets.DOCKER_PASSWORD }}
        DOCKER_USERNAME: ${{ secrets.DOCKER_USERNAME }}
    - name: GitHub Action for Docker
      uses: actions/docker/cli@master
      with:
        args: build -t chunjiesamliu/jrocker:latest . && docker push chunjiesamliu/jrocker:latest
```
