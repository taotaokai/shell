#!/bin/bash

# check if IRIS_fdsnws_dataselect.sh finishes OK. 

ls */IRIS_fdsnws_dataselect.log  | xargs tail -n 4 | grep saved | sed "s/.*IRIS\/\(.*\)\/mseed.*/\1/" > download.ok

grep -vf download.ok event.txt > download.failed
