#!/bin/bash

# download IRIS data from a given event list
# dependecies:
# - IRIS Fetch scripts:
# - IRIS mseed2sac

#===== parameters

#-- events
event_file=${1:-event.txt}

#-- station/channel selection
staLatRange=15:55
staLonRange=90:145

#-- data time window referred to first arrival
twin_b=-300
twin_e=1500

#-- outputs
out_dir=${2:-events}

#===== get data for each event 

wkdir=$(pwd)
tmpdir=$wkdir

echo "########################"
echo "###### Parameters ######"
echo "########################"
echo "event_file=$event_file"
echo "staLatRange=$staLatRange "
echo "staLonRange=$staLonRange "
echo "twin_b=$twin_b"
echo "twin_e=$twin_e"
echo "out_dir=$out_dir"
echo "########################"

IFS=$'\n'
for event_line in $(cat $event_file)
do

  echo
  echo "#===== $event_line"
  echo

  #-- parse each line in the event list
  IFS='|' read -ra evt <<< "$event_line"
  evid=${evt[7]#*,} # use GCMT event ID
  evodate=${evt[1]/T/ } # replace T with a whitespace
  evla=${evt[2]//[[:space:]]/} # remove whitespaces
  evlo=${evt[3]//[[:space:]]/}
  evdp=${evt[4]//[[:space:]]/}

  evoepoch=$(date -u -d "$evodate UTC" +%s)
  evojday=$(date -u -d "$evodate UTC" +%Y,%j,%H:%M:%S.%N)

  #-- check potential stations/channels

  # data should exist 1 day before the earthquake
  starttime=$(date -u -d "$evodate UTC 1 days ago" +%Y-%m-%dT%H:%M:%S)

  tmp_station=$(mktemp --tmpdir=$tmpdir)

  FetchData.pl \
    -s $starttime \
    --lat $staLatRange --lon $staLonRange \
    -C 'BHZ' \
    -m $tmp_station

  # check if there are available stations
  nsta=$(awk 'END{print NR}' $tmp_station)
  nsta=$(( $nsta - 1 ))
  if [ $nsta -lt 1 ]; then
    echo "[WARN] No available stations found!"
    continue
  fi
  echo "# stations found: $nsta"

  #-- create directories
  evtdir=$wkdir/$out_dir/$evid

  if [ -d ${evtdir} ]; then
    rm -rf $evtdir
  fi

  mkdir -p $evtdir
  if [ -d $evtdir ]; then
    cd $evtdir
    mkdir mseed sac sacpz data
    cd $wkdir
  else
    echo "[ERROR] cannot create directory: $evtdir"
  fi

  #-- save event info
  echo $event_line > $evtdir/data/event.txt

  #-- generate data selection list
  cat /dev/null > $evtdir/data/dataselect.txt

  for station_line in $(grep -v "^#" $tmp_station)
  do

    IFS='|' read -ra sta <<< "$station_line"
    netwk=${sta[0]}
    stnm=${sta[1]}
    loc=${sta[2]}
    stla=${sta[4]}
    stlo=${sta[5]}

    # get traveltime
    ttp=$(aktimes $evla $evlo $evdp $stla $stlo $AK135_TABLE | \
          head -n1 | awk '{print $1}')

    if [ -z "$ttp" ]; then
      echo "[ERROR] wrong travaltime ttp=$ttp"
      exit
    fi

    # get time window
    t0=$(echo "$evoepoch + $twin_b + $ttp" | bc -l)
    t1=$(echo "$evoepoch + $twin_e + $ttp" | bc -l)
    starttime=$(date -u -d @$t0 +%Y-%m-%dT%H:%M:%S)
    endtime=$(date -u -d @$t1 +%Y-%m-%dT%H:%M:%S)

    # print data selection list 
    printf "%-6s %-6s %-4s  BH?  %s  %s\n" \
      $netwk $stnm $loc $starttime $endtime \
      >> $evtdir/data/dataselect.txt

  done

  rm $tmp_station

  #-- get data (waveform & sacpz)
  FetchData.pl \
    -l $evtdir/data/dataselect.txt \
    -o $evtdir/mseed/$evid.mseed \
    -m $evtdir/data/station.txt \
    -sd $evtdir/sacpz

  #-- convert mseed to sac
  cd $evtdir/sac
  mseed2sac \
    -m $evtdir/data/station.txt \
    -E "$evojday/$evla/$evlo/$evdp/${evt[9]}" \
    $evtdir/mseed/$evid.mseed

  cd $wkdir

done