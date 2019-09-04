---
title: "git config"
author: "Chun-Jie Liu"
date: "2019-09-04"
---

## Git config shortcuts

The code chunk below is my daily used git config shortcuts.

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
	br = branch -a
	bd = branch -d
	lg = log --color --graph --pretty=format:'%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset' --abbrev-commit
	df = diff
	mg = merge --no-ff -m
	rm = rm -r --cached
	sv = status -v -v
[filter "lfs"]
	clean = git-lfs clean -- %f
	smudge = git-lfs smudge -- %f
	process = git-lfs filter-process
	required = true
```