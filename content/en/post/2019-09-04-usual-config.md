---
title: "Configurations"
author: "Chun-Jie Liu"
date: "2019-09-04"
---

The configuration files were deposited into my [Gist](https://gist.github.com/chunjie-sam-liu), but the [Gist](https://gist.github.com/chunjie-sam-liu) access was blocked. In case of the not accessing, I copyed some configure files here.

## config in .ssh

```shell
Host *
  TCPKeepAlive yes
  ServerAliveInterval 120
  ServerAliveCountMax 5

Host guo-4
  HostName ip
  User username
  Port port

Host guo-10
  HostName ip
  User username
  Port port
  ProxyJump guo-4

Host guo-3
  HostName ip
  User username
  Port port
  ProxyJump guo-4

Host guo-1
  HostName ip
  User username
  Port port
  ProxyJump guo-4

Host mgt
  HostName ip
  User username
  Port port
```

## rstudio_bindings.json

```json
{
    "commentUncomment" : "Cmd+/",
    "fold" : "Cmd+[|Alt+Cmd+L|Cmd+F1",
    "foldAll" : "Alt+Cmd+[",
    "goToFileFunction" : "Ctrl+/",
    "unfold" : "Cmd+]|Shift+Alt+Cmd+L|Shift+Cmd+F1",
    "unfoldAll" : "Alt+Cmd+]"
}
```

## .gitconfig

```shell
[user]
	name = chunjie-sam-liu
	email = chunjie-sam-liu@foxmail.com
[color]
	diff = auto
	status = auto
	branch = auto
[alias]
	cl = clone
	ad = add
	co = checkout
	ss = status
	cm = commit -m
	br = branch -a -vv
	bd = branch -d
  bu = branch -u
	lg = log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit
	df = diff
	mg = merge --no-ff -m
	rm = rm --cached
	sv = status -v -v
[filter "lfs"]
	clean = git-lfs clean -- %f
	smudge = git-lfs smudge -- %f
	process = git-lfs filter-process
	required = true
```


## .Rproifle

```R
# Blogdown options --------------------------------------------------------
options(blogdown.author = "Chun-Jie Liu")
options(servr.daemon = FALSE)
options(blogdown.ext = ".Rmd")
options(blogdown.subdir = "post")
options(blogdown.yaml.empty = TRUE)

# General options ---------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(prompt = "Jrocker>", digits = 4, show.signif.stars = FALSE)
options(stringsAsFactors = FALSE)

# ggplot2 v3 options ------------------------------------------------------
options(
  ggplot2.continuous.color = "viridis",
  ggplot2.continuous.fill = "viridis"
)

# First -------------------------------------------------------------------

.First <- function(){
  library(magrittr)
}

# Last --------------------------------------------------------------------

.Last <- function(){}
```
