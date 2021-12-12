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
name: jrocker

on:
  push:
    branches:
      - master

jobs:
  docker:
    runs-on: ubuntu-latest
    steps:
      -
        name: Set up QEMU
        uses: docker/setup-qemu-action@v1
      -
        name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v1
      -
        name: Login to DockerHub
        uses: docker/login-action@v1
        with:
          username: ${{ secrets.DOCKER_USERNAME }}
          password: ${{ secrets.DOCKER_PASSWORD }}
      -
        name: Build and push
        id: docker_build
        uses: docker/build-push-action@v2
        with:
          push: true
          tags: chunjiesamliu/jrocker:latest
```
