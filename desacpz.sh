#!/bin/bash

### correct instrument response from SACPZ files using SAC

#====== read command line args
wkdir=${1:?[arg] need event_dir}

outtype=${2:?[arg] need output type, one of dis,vel,acc}
if [ x${outtype} != xdis ] && \
   [ x${outtype} != xvel ] && \
   [ x${outtype} != xacc ]; then
    echo "[ERROR] outtype must be dis,vel,acc"
    exit -1
fi

freqlim=${3:-"0.009 0.01 30.0 35.0"}

channel_list=$wkdir/${4:-data/channel.txt}

sac_dir=$wkdir/${5:-sac}
sacpz_dir=$wkdir/${6:-sacpz}

# output
out_dir=$wkdir/${7:-${outtype}}

# log file
log_file=$wkdir/desacpz_to_${outtype}.log

echo "#LOG created on $(date +%Y-%m-%dT%H:%M:%S)" > $log_file
echo "#PWD: $HOSTNAME:$(pwd)" >> $log_file
echo "#COMMAND: $0 $@" >> $log_file
echo "#Parameter: wkdir= $wkdir " >> $log_file
echo "#Parameter: outtype= $outtype " >> $log_file
echo "#Parameter: freqlim= $freqlim " >> $log_file
echo "#Parameter: channel_list= $channel_list " >> $log_file
echo "#Parameter: sac_dir= $sac_dir" >> $log_file
echo "#Parameter: sacpz_dir= $sacpz_dir" >> $log_file
echo "#Parameter: out_dir= $out_dir" >> $log_file


#====== create directories if not exist
if [ ! -d $out_dir ];then
    mkdir $out_dir
fi

# set sac instrument type (transfer)
sac_inst=$outtype
if [ x$sac_inst == xdis ]; then
    sac_inst=none
fi

#====== SAC
echo >> $log_file
echo "#Run SAC commands" >> $log_file
echo >> $log_file

tmp_log=$(mktemp)

# loop over echo channel in channel_list
grep -v ^# $channel_list |\
while IFS='|' read net sta loc chan dummy
do

    channel_id="${net}.${sta}.${loc}.${chan}"

    # get corresponding sac/sacpz files
    sac_file=(${sac_dir}/*${channel_id}.*SAC)
    if [ ! -f "${sac_file[0]}" ]; then
        echo "[WARNING] SKIP: SAC file not found for ${channel_id}"
        continue
    fi
    if [ ${#sac_file[@]} -gt 1 ]; then
        echo "[WARNING] SKIP: segmented SAC files found for ${channel_id}"
        echo "${sac_file[@]}"
        continue
    fi

    #sacpz_file=(${sacpz_dir}/SAC_PZs_${net}_${sta}_${chan}_${loc}*)
    sacpz_file=(${sacpz_dir}/${channel_id}*)
    if [ ! -f "${sacpz_file[0]}" ]; then
        echo "[WARNING] SKIP: SACPZ file not found for ${channel_id}"
        continue
    fi
    if [ ${#sacpz_file[@]} -gt 1 ]; then
        echo "[WARNING] SKIP: more than one SACPZ file found for ${channel_id}"
        echo "${sacpz_file[@]}"
        continue
    fi

    # sac
    out_file=$out_dir/${channel_id}
    sac<<EOF | tee $tmp_log
r ${sac_file[0]}
rmean;rtrend
rmean;rtrend
rmean;rtrend
rmean;rtrend
trans from polezero s ${sacpz_file[0]} to ${sac_inst} freq ${freqlim}
w $out_file
q
EOF
    # check for SAC error
    if grep -Fiq ERROR $tmp_log; then
        echo "[WARNING] SAC failed removing SACPZ for ${channel_id}"
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