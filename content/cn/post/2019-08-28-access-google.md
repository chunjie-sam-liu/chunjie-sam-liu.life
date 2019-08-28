---
title: "科学上网"
author: "Chun-Jie Liu"
date: "2019-08-28"
---

科学上网很多个人博客都有分享，但是每隔一段时间会被删除，就算[GitHub shadowsocks](https://github.com/shadowsocks/shadowsocks)也需要做一个`rm branch`。应该遵守国家相关法律法规，科学上网。

但是，奈何`百度`有毒，阻碍有效搜索结果。并且，很多时候使用高效搜索都是因为需要`google scholar`高效搜索需要的文献。由于此Blog的受众应该是小众的科学工作者，在此分享一下在`VPS`上部署`ss`。

[代码](https://github.com/chunjie-sam-liu/useful-scripts/blob/master/ss-deployment.sh)是从多种不同的地方看到的，在此汇总一下，方便`VPS`的快速部署。在使用的过程中由于无法通过`ipv4`来访问`google scholar`，需要开启`ipv6`，同时将`ipv6`添加到`hosts`里面。

```shell
# get shadowsocks
wget https://raw.githubusercontent.com/chunjie-sam-liu/useful-scripts/master/ss.sh

# install shadowsocks
chmod a+x ss.sh
bash ss.sh
# specify password: password
# connecting port: 1070
# encryption: aes-256-cfb - 7

# open filewall
firewall-cmd --permanent --zone=public --add-port=1070/tcp
firewall-cmd --reload

# echo ipv4
echo "echo 3 > /proc/sys/net/ipv4/tcp_fastopen" >> /etc/rc.local
echo "net.ipv4.tcp_fastopen = 3" >> /etc/sysctl.conf

# restart shadowsocks
/etc/init.d/shadowsocks restart

# google scholar ipv6
echo "2404:6800:4008:c06::be scholar.google.com" >> /etc/hosts
echo "2404:6800:4008:c06::be scholar.google.com.hk" >> /etc/hosts
echo "2404:6800:4008:c06::be scholar.google.com.tw" >> /etc/hosts

# restart shadowsocks
/etc/init.d/shadowsocks restart

# run bbr
wget --no-check-certificate https://raw.githubusercontent.com/chunjie-sam-liu/useful-scripts/master/ss-bbr.sh
chmod +x ss-bbr.sh
bash ss-bbr.sh
```