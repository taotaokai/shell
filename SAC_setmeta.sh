#!/bin/bash

### Set sac header from channel.txt and event.txt 
# event: IRIS-fdsnws-event with format=text
# channel.txt: IRIS-fdsnws-station with level=channel,format=text

#====== read command line args
event_dir=${1:?[arg]need event_dir}

event_dir=$(readlink -f $event_dir)

event_id=${2:?[arg]need event_id}
sac_dir=$event_dir/${3:-dis}
channel_list=$event_dir/${4:-data/channel.txt}
event_list=$event_dir/${5:-data/event.txt}

# log file
log_file=$event_dir/SAC_setmeta.log

echo "# LOG created on $(date +%Y-%m-%dT%H:%M:%S)" > $log_file
echo "# PWD: $HOSTNAME:$(pwd)" >> $log_file
echo "# COMMAND: $0 $@" >> $log_file
echo "# Parameters: " >> $log_file
echo "#   event_dir=$event_dir " >> $log_file
echo "#   event_id=$event_id " >> $log_file
echo "#   channel_list=$channel_list " >> $log_file
echo "#   event_list=$event_list " >> $log_file
echo "#   sac_dir=$sac_dir " >> $log_file

#====== SAC
echo >> $log_file
echo "# Run SAC commands" >> $log_file
echo >> $log_file

tmp_log=$(mktemp)

#event=( $(grep $event_id $event_list | sed "s/|/ /g") )
event=$(grep $event_id $event_list)
evla=$(echo $event | awk -F"|" '{print $3}' )
evlo=$(echo $event | awk -F"|" '{print $4}' )
evdp=$(echo $event | awk -F"|" '{print $5}' )
mag=$(echo $event | awk -F"|" '{print $11}' )

origin_time=$(echo $event | awk -F"|" '{gsub(/T/," ",$2); print $2}' )
ogmt=$(date -ud "$origin_time" +"%Y %j %H %M %S %N" |\
  awk '{print $1,$2,$3,$4,$5,$6/1000000}')

# loop over echo channel in channel_list
grep -v ^# $channel_list |\
while IFS='|' read net sta loc chan stla stlo stel stdp cmpaz cmpdip dummy
do
    channel_id="${net}.${sta}.${loc}.${chan}"

    # get corresponding sac/sacpz files
    sac_file=(${sac_dir}/*${channel_id}*)
    if [ ! -f "${sac_file[0]}" ]; then
        echo "[WARN] SAC file not found for ${channel_id}, SKIP"
        continue
    fi
    if [ ${#sac_file[@]} -gt 1 ]; then
        echo "[WARN] Segmented SAC files found for ${channel_id}, SKIP"
        echo "${sac_file[@]}"
        continue
    fi

    # set sac header
    cmpinc=$(echo "90.0 + $cmpdip" | bc -l)
    sac<<EOF | tee $tmp_log
r ${sac_file[0]}
ch lcalda true
ch stla $stla stlo $stlo stel $stel stdp $stdp cmpaz $cmpaz cmpinc $cmpinc 
ch evla $evla evlo $evlo evdp $evdp mag $mag
ch o gmt $ogmt
lh
wh
q
EOF

    # check for SAC error
    if grep -Fiq ERROR $tmp_log; then
        echo "[WARN] SAC failed to set metadata for ${channel_id}, output REMOVED"
        cat $tmp_log
        continue
    fi
    
done >> $log_file 2>&1

rm $tmp_log


#====== print out error/WARN messages
echo >> $log_file
echo "# Finished on $(date +%Y-%m-%dT%H:%M:%S)" >> $log_file

echo "# ERROR/WARNING encoutered:"
echo
grep -i -e "error" -e "warn" $log_file