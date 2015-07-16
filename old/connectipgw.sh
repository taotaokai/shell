#!/bin/bash

username=0006163488
passwd=notamatter1

wget --no-check-certificate -Y off -T 10 -t 3 -O /dev/null \
"http://162.105.67.5/ipgw/ipgw.ipgw?uid=$username&password=$passwd&timeout=1&range=2&operation=disconnectall"

sleep 1s

wget --no-check-certificate -Y off -T 10 -t 3 -O /dev/null \
"http://162.105.67.5/ipgw/ipgw.ipgw?uid=$username&password=$passwd&timeout=1&range=2&operation=connect"
