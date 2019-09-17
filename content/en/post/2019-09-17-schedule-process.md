---
title: "Schedule process on linux"
author: "Chun-Jie Liu"
date: "2019-09-17"
---

## Scheduling tasks

Automating system maintenance and administration of web server save my life for managing our web server. Scheduling certain tasks to run periodically or one time to backup data and start/stop server at certain time. `crontab` and `at` are to service processes to fulfill the requirement of a server manager.

## Run periodically `crontab`

![Time format](../../../img/misc-imgs/crontab-time.png)

The above picture is the `crontab` rule to schedule periodically tasks. Following code is an example to check `MySQL` status every five minutes.

```shell
# add in crontab
# crontab -e with root
*/5 * * * * /etc/mysql/listen.sh
```

```shell
#!/bin/bash
#apache
HTTP_CODE=`curl -o /dev/null -s -w "%{http_code}" "http://example.com"`
if [ $HTTP_CODE != 200 ]; then
    echo $HTTP_CODE>> /var/log/mysql_listen.log
    service apache2 start
else
    echo "lab server running">> /var/log/mysql_listen.log
fi

#mysql
pgrep mysqld &> /dev/null
if [ $? -gt 0 ]; then
    echo "`date` mysql is stop">> /var/log/mysql_listen.log
    service mysql start
else
    echo "`date` mysql running">> /var/log/mysql_listen.log
fi

#mongodb
netstat -anop | grep localhost:port
if [ $? -ne 1 ]; then
    echo "`date` mongodb running">> /var/log/mongodb_listen.log
else
    echo $(date +%T%n%F)" Restart mongodb Services " >> /var/log/mongodb_listen.log
    service mongod restart
fi
```

## Run one time task `at`

Schedule certain time to start/stop `apache` service.

```shell
# schedule one time process for stop apache2 at 06:00 2019-09-18
echo "/usr/sbin/service apache2 stop" | at -m 06:00 2019-09-18
# schedule one time process for start apache2 at 14:00 2019-09-20
echo "/usr/sbin/service apache2 start" | at -m 14:00 2019-09-20
```
