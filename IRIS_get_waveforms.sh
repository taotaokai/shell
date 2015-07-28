#!/bin/bash

# download IRIS data from a given event list

# dependecies:
# - IRIS Fetch scripts:
# - IRIS mseed2sac
# - TAUP:taup_time

usage(){
  echo
  echo "Syntax:"
  echo "  $0 -s 15:55/90:145 -t p,P/-300:1500 -l- -o ./"
  echo
  echo "Parameters:"
  echo "  -s  geographic range of stations: lat0:lat1/lon0:lon1"
  echo "  -t  time window phase/begin/end. If multiple phases specified, \
use the smallest traveltime as the reference time."
  echo "  -l  list of events to download (default read from STDIN)."
  echo "  -o  output directory."
  echo 
  exit -1
}

#===== parameters

#defaults
event_list=-
geo_range=15:55/90:145
twin=p,P/-300/1500
out_dir=./

# parse options
while getopts s:t:l:o:h name
do
    case $name in
    s)  geo_range="$OPTARG";;
    t)  twin="$OPTARG";;
    l)  event_list="$OPTARG";;
    o)  out_dir="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

echo "##### The program is run as:"
echo "# $0 -s $geo_range -t $twin -l $event_list -o $out_dir"

staLatRange=${geo_range%/*}
staLonRange=${geo_range#*/}

if [ -z "$staLatRange" -o -z "$staLonRange" ]
then
  echo "[WARN] Wrong geographic range: lat=$staLatRange lon=$staLatRange "
  exit -1
fi

phase=$(echo ${twin%/*})
twin_b=$(echo ${twin#*/} | awk -F":" '{print $1}')
twin_e=$(echo ${twin#*/} | awk -F":" '{print $2}')

if [ -z "$phase" -o -z "$twin_b" -o -z "$twin_e" ]
then
  echo "[WARN] Wrong time window: phase=$phase b=$twin_b e=$twin_e"
  exit -1
fi

#===== get data for each event 

wkdir=$(pwd)

IFS=$'\n'
for event_line in $(cat $event_list)
do

  echo
  echo "===== EVENT: $event_line"
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
  echo "----- query for stations/channels"

  # data should exist 1 day before the earthquake
  starttime=$(date -u -d "$evodate UTC 1 days ago" +%Y-%m-%dT%H:%M:%S)

  tmp_station=$(mktemp)

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
    mkdir data mseed resp
  else
    echo "[ERROR] cannot create directory: $evtdir"
  fi

  #-- save event info
  echo $event_line > $evtdir/data/event.txt

  #-- generate data selection list

  echo "----- make data selection list"

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
    ttp=$(taup_time -evt $evla $evlo -h $evdp \
          -sta $stla $stlo --time -ph $phase | head -n1 | awk '{print $1}')

    if [ -z "$ttp" ]; then
      echo "[ERROR] taup_time failed to get travaltime ttp=$ttp"
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
  
  echo "----- query for waveform and sacpz files"

  FetchData.pl \
    -l $evtdir/data/dataselect.txt \
    -o $evtdir/mseed/$evid.mseed \
    -m $evtdir/data/station.txt \
    -rd $evtdir/resp

done