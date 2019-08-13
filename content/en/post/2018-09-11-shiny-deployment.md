---
title: "Shiny deployment"
author: "Chun-Jie Liu"
date: "2018-09-11"
---


## Workspace mode

Before pointing workspace to the `shiny-server` srv directory, change the worksapce mode with `g+s`.

```
$chown -R liucj.shiny-apps GSCALite
$chmod -R g+s GSCALite
$ln -s GSCALite /srv/shiny-server/
```

## Apache2 configuration

The add `shiny-server` configuration to `Aapache2`. On Ubuntu system, the `Apache2` configuration file is `000-default.conf`.

The default http 80 port was used for default web directory `/var/wwww`. Change the proxy with provided port.

```
RewriteEngine on
RewriteCond %{HTTP:Upgrade} !=websocket
RewriteRule /web/(.*) http://localhost:port/$1 [P,L]
ProxyPass /web/ http://localhost:port/
ProxyPassReverse /web/ http://localhost:port/
```

## `shiny-server` configuration

Edit `shiny-server.conf` with root rights with following code.

Instruct Shiny Server to run applications as the user "shiny" run_as shiny; Define a server that listens on port port

```
server {
  listen port;

  # Define a location at the base URL
  location / {

    # Host the directory of Shiny Apps stored in this directory
    site_dir /srv/shiny-server;

    # Log all Shiny output to files in this directory
    log_dir /var/log/shiny-server;

    # When a user visits the base URL rather than a particular application,
    # an index of the applications available in this directory will be shown.
    directory_index on;
  }
}
```
