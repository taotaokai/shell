#!/bin/bash

# process all events
wkdir=$(pwd)
bin_dir=$wkdir/bin
gcmt_ndk=$wkdir/all.ndk

event_list=${1:?[arg] need fdsnws-event list}

for event_id in $(awk -F"|" '$1!~/#/{print $9}' $event_list)
do
    echo
    echo ====== processing $event_id
    echo
    event_dir=$wkdir/$event_id

    # if you have set a-w, but want to re-do the whole process
    chmod u+w -R $event_dir/*
    chmod a-w -R $event_dir/mseed
    chmod a-w -R $event_dir/sacpz

    echo ------ CMTSOLUTION
    #rm -rf $event_dir/data
    #mkdir -p $event_dir/data
    make_CMTSOLUTION.sh $event_id $gcmt_ndk $event_dir/data/CMTSOLUTION 

    #echo ------ make channel list 
    #$bin_dir/make_channel_list.py \
    #    $event_dir/data/CMTSOLUTION \
    #    $wkdir/metadata/channel.txt \
    #    $event_dir/data/channel.txt.1
    #sort $event_dir/data/channel.txt.1 > $event_dir/data/channel.txt
    #rm $event_dir/data/channel.txt.1
    
    echo ------ extract mseed
    rm -rf $event_dir/sac
    mkdir -p $event_dir/sac
    cd $event_dir/sac
    cat $event_dir/data/channel.txt |\
        awk -F"|" 'BEGIN{OFS="|"}; ($1!~/^#/ && $3==""){$3="--"}; {print $0}' \
        > metadata.txt 
    mseed2sac $event_dir/mseed/*.mseed -m metadata.txt > mseed2sac.log 2>&1

    echo ------ deconv sacpz 
    cd $event_dir
    #rm -rf sacpz vel dis
    #ln -sf $wkdir/metadata/sacpz
    $bin_dir/desacpz.sh $event_dir vel "0.005 0.008 5 8"
    $bin_dir/desacpz.sh $event_dir dis "0.005 0.008 5 8"

    echo ------ set sac header
    SAC_setmeta.sh $event_dir $event_id dis
    SAC_setmeta.sh $event_dir $event_id sac

    echo ------ reset permission
    chmod a-w -R $event_dir/*
done
