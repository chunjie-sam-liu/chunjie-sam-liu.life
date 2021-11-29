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

## rstudio-prefs.json

```json
{
    "console_double_click_select": true,
    "save_workspace": "never",
    "load_workspace": false,
    "initial_working_directory": "~/github",
    "soft_wrap_r_files": true,
    "highlight_selected_line": true,
    "show_invisibles": true,
    "show_indent_guides": true,
    "blinking_cursor": false,
    "syntax_color_console": true,
    "scroll_past_end_of_document": true,
    "highlight_r_function_calls": true,
    "default_encoding": "UTF-8",
    "font_size_points": 14,
    "editor_theme": "Xcode",
    "posix_terminal_shell": "bash",
    "panes": {
        "quadrants": [
            "Source",
            "Console",
            "TabSet2",
            "TabSet1"
        ],
        "tabSet1": [
            "Files",
            "Connections",
            "Packages",
            "Help",
            "Tutorial",
            "Presentation"
        ],
        "tabSet2": [
            "Environment",
            "History",
            "Plots",
            "Build",
            "VCS",
            "Viewer"
        ],
        "hiddenTabSet": [],
        "console_left_on_top": false,
        "console_right_on_top": true,
        "additional_source_columns": 1
    },
    "show_last_dot_value": true,
    "always_save_history": false,
    "rainbow_parentheses": true,
    "tab_multiline_completion": true,
    "code_completion_delay": 100,
    "python_type": "virtualenv",
    "python_version": "3.8.10",
    "python_path": "/opt/virtualenvs/r-tensorflow/bin/python",
    "server_editor_font_enabled": true,
    "server_editor_font": "Menlo",
    "source_with_echo": true,
    "jobs_tab_visibility": "shown",
    "show_panel_focus_rectangle": true,
    "auto_append_newline": true,
    "strip_trailing_whitespace": true,
    "auto_save_on_blur": true,
    "show_help_tooltip_on_idle": true,
    "limit_visible_console": true,
    "auto_expand_error_tracebacks": true,
    "code_completion_characters": 2
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
[pager]
	branch = false
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
