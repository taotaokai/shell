#!/bin/bash

# extract sac from mseed files with alternaltive response file (dataless seed).

#====== read command line args
wkdir=${1:-.}

mseed_dir=$wkdir/${2:-mseed}
dataless_seed=$mseed_dir/${3:-dataless.seed}

# output
sac_dir=$wkdir/${4:-sac}

# log file
log_file=$wkdir/extract_mseed.log

echo "#LOG created on $(date +%Y-%m-%dT%H:%M:%S)" > $log_file
echo "#PWD: $HOSTNAME:$(pwd)" >> $log_file
echo "#COMMAND: $0 $@" >> $log_file
echo "#Parameter: wkdir= $wkdir" >> $log_file
echo "#Parameter: mseed_dir= $mseed_dir" >> $log_file
echo "#Parameter: dataless_seed= $dataless_seed" >> $log_file
echo "#Parameter: sac_dir= $sac_dir" >> $log_file
echo "#Parameter: log_file= $log_file" >> $log_file


#====== create output directories if not exist 
if [ ! -d ${sac_dir} ]; then
  mkdir $sac_dir
fi


#====== rdseed
echo >> $log_file
echo "#running rdseed" >> $log_file
echo >> $log_file

find $mseed_dir -name "*.mseed" |\
while read mseed_file 
do

  rdseed -df $mseed_file -g $dataless_seed -q $sac_dir

done >> $log_file 2>&1


#====== print out error/warning messages
echo >> $log_file
echo "#Finished on $(date +%Y-%m-%dT%H:%M:%S)" >> $log_file

echo "# Error/Warning encoutered:"
echo
grep -i -e "error" -e "warn" $log_file