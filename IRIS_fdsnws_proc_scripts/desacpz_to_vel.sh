#!/bin/bash

### correct instrument response from SACPZ files using SAC

#====== read command line args
wkdir=${1:-.}

freqlim=${2:-"0.009 0.01 30.0 35.0"}

station_list=$wkdir/${3:-data/station.txt}
sac_dir=$wkdir/${4:-sac}
sacpz_file=$wkdir/${5:-data/sacpz.txt}

# output
out_dir=$wkdir/${6:-vel}

# log file
log_file=$wkdir/desacpz.log

echo "#LOG created on $(date +%Y-%m-%dT%H:%M:%S)" > $log_file
echo "#PWD: $HOSTNAME:$(pwd)" >> $log_file
echo "#COMMAND: $0 $@" >> $log_file
echo "#Parameter: wkdir= $wkdir " >> $log_file
echo "#Parameter: station_list= $station_list " >> $log_file
echo "#Parameter: freqlim= $freqlim " >> $log_file
echo "#Parameter: sac_dir= $sac_dir" >> $log_file
echo "#Parameter: sacpz_file= $sacpz_file" >> $log_file
echo "#Parameter: out_dir= $out_dir" >> $log_file


#====== Sanity check
if [ ! -d $out_dir ];then
    mkdir $out_dir
fi

if [ ! -f "${sacpz_file}" ]; then
    echo "[ERROR] SACPZ file not found: $sacpz_file"
    exit -1 
fi


#====== SAC
echo >> $log_file
echo "#Run SAC commands" >> $log_file
echo >> $log_file

tmp_log=$(mktemp)

# loop over echo channel in station_list
grep -v ^# $station_list |\
while IFS='|' read net sta loc chan dummy
do

    # get corresponding sac files
    sac_file=(${sac_dir}/*.${net}.${sta}.${loc}.${chan}.*SAC)
    if [ ! -f "${sac_file[0]}" ]; then
        echo "[WARNING] SKIP: SAC file not found for ${net}.${sta}.${loc}.${chan}"
        continue
    fi
    if [ ${#sac_file[@]} -gt 1 ]; then
        echo "[WARNING] SKIP: segmented SAC files found for ${net}.${sta}.${loc}.${chan}"
        echo "${sac_file[@]}"
        continue
    fi

    # sac
    out_file=$out_dir/${net}.${sta}.${loc}.${chan}
    sac<<EOF | tee $tmp_log
r ${sac_file[0]}
rtrend
trans from polezero s ${sacpz_file} to vel freq ${freqlim}
w $out_file
q
EOF
    # check for SAC error
    if grep -Fiq ERROR $tmp_log; then
        echo "[WARNING] SAC failed removing SACPZ for ${net}.${sta}.${loc}.${chan}"
        if [ -f $out_file ];then
            rm $out_file
        fi
        continue
    fi
    
done >> $log_file 2>&1

rm $tmp_log


#====== print out error/warning messages
echo >> $log_file
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)" >> $log_file

echo "# Error/Warning encoutered:"
echo
grep -i -e "error" -e "warn" $log_file