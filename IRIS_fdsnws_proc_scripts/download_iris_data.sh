#!/bin/bash

#IRIS_fdsnws_dataselect.sh -s b/0/80/80/160 -t p,P/-500/1500 -o . event.txt

# re-download the failed events
#awk -F"|" '{printf "chmod u+w -R %s; rm -rf %s\n",$9,$9}' download.failed > rm_sh
#bash rm_sh
#IRIS_fdsnws_dataselect.sh -s b/0/80/80/160 -t p,P/-500/1500 -o . download.failed


# download additional events
IRIS_fdsnws_dataselect.sh -s b/0/80/80/160 -t p,P/-500/1500 -o . event_to_add.txt