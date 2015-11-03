#!/bin/bash

### correct instrument response from SACPZ files using SAC

#====== read command line args
wkdir=${1:-.}

station_list=$wkdir/${2:-data/station.txt}
freqlim=${3:-"0.014 0.015 30.0 35.0"}

sac_dir=$wkdir/${4:-sac}
sacpz_dir=$wkdir/${5:-sacpz}

# output
disp_dir=$wkdir/${6:-disp}

# log file
log_file=$wkdir/desacpz.log

echo "#LOG created on $(date +%Y-%m-%dT%H:%M:%S)" > $log_file
echo "#PWD: $HOSTNAME:$(pwd)" >> $log_file
echo "#COMMAND: $0 $@" >> $log_file
echo "#Parameter: wkdir= $wkdir " >> $log_file
echo "#Parameter: station_list= $station_list " >> $log_file
echo "#Parameter: freqlim= $freqlim " >> $log_file
echo "#Parameter: sac_dir= $sac_dir" >> $log_file
echo "#Parameter: sacpz_dir= $sacpz_dir" >> $log_file
echo "#Parameter: disp_dir= $disp_dir" >> $log_file


#====== create directories if not exist
if [ ! -d $disp_dir ];then
    mkdir $disp_dir
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

    # get corresponding sac/sacpz files
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

    sacpz_file=(${sacpz_dir}/SAC_PZs_${net}_${sta}_${chan}_${loc}*)
    if [ ! -f "${sacpz_file[0]}" ]; then
        echo "[WARNING] SKIP: SACPZ file not found for ${net}.${sta}.${loc}.${chan}"
        continue
    fi
    if [ ${#sacpz_file[@]} -gt 1 ]; then
        echo "[WARNING] SKIP: more than one SACPZ file found for ${net}.${sta}.${loc}.${chan}"
        echo "${sacpz_file[@]}"
        continue
    fi

    # sac
    out_file=$disp_dir/${net}.${sta}.${loc}.${chan}
    sac<<EOF | tee $tmp_log
r ${sac_file[0]}
rtrend
trans from polezero s ${sacpz_file[0]} to none freq ${freqlim}
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