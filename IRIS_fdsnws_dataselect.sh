#!/bin/bash

wkdir=$(pwd)

# FDSN web services
#fdsnws_event="http://service.iris.edu/fdsnws/event/1/query"
fdsnws_station="http://service.iris.edu/fdsnws/station/1/query"
fdsnws_dataselect="http://service.iris.edu/fdsnws/dataselect/1/query"
irisws_sacpz="http://service.iris.edu/irisws/sacpz/1/query"

#======
usage(){
cat<<EOF
NAME

  IRIS_fdsnws_dataselect - a shell wrapper of IRIS web service to download waveform data.

SYNOPSIS

  IRIS_fdsnws_dataselect [-s <station_specs>] [-t <twin_specs>] [-o <out_dir>] <gcmt_catalog>

DESCRIPTION

  <gcmt_catalog>  a catalog list generated by fdsnws-event with catalog=GCMT

  Options:

  -s <station_specs> [b/20/60/90/150]
    box region: b/lat0/lat1/lon0/lon1
    radial region: r/lat0/lon0/radius0/radius1

  -t <twin_specs> [p,P/-500/1500]
    time window phase(s)/begin/end. If multiple phases specified,
    use the smallest traveltime as the reference time.

  -o <out_dir> [work directory]

NOTE
  1) TauP/bin/taup_time is called to calculate traveltime
  2) channel code is restricted to [FCHB]H?, i.e. only broadband high gain seismic data

EOF
}

#====== parameters

# default option args
station_specs=b/20/60/90/150
twin_specs=p,P/-500/1500
out_dir=$wkdir

# parse options
while getopts s:t:o:h name
do
  case $name in
  s)  station_specs="$OPTARG";;
  t)  twin_specs="$OPTARG";;
  o)  out_dir=$(readlink -f "$OPTARG");;
  [h,?])  usage; exit -1
  esac
done

# get all arguments
shift $((OPTIND-1))
event_list=${1:?[arg] need gcmt_catalog}
event_list=$(readlink -f $event_list)

# process arguments
geo_range=$(echo $station_specs | awk -F"/" '\
  $1=="b"{printf \
    "minlatitude=%s&maxlatitude=%s&minlongitude=%s&maxlongitude=%s",\
    $2,$3,$4,$5}; \
  $1=="r"{printf \
    "latitude=%s&longitude=%s&minradius=%s&maxradius=%s",$2,$3,$4,$5}')

phase_name=$(echo ${twin_specs} | cut -d / -f 1)
twin_begin=$(echo ${twin_specs} | cut -d / -f 2)
twin_end=$(echo ${twin_specs} | cut -d / -f 3)

if [ -z "$phase_name" -o -z "$twin_begin" -o -z "$twin_end" ]
then
  echo "[ERROR] time window: phase=$phase_name b=$twin_begin e=$twin_end"
  exit 1
fi

#====== get data for each event 

# LOG
echo "[LOG] wkdir=$wkdir"
echo "[LOG] $0 $event_list -s $station_specs -t $twin_specs -o $out_dir"

old_IFS="$IFS"
IFS=$'\n' 
for event_line in $(grep -v "^[ \t]*#" $event_list)
do
  echo
  echo "[LOG] ====== $event_line"
  echo

  #------ event info
  IFS='|' read -ra evt <<< "$event_line"
  evid=${evt[8]#*,} # use GCMT event ID
  evodate=${evt[1]/T/ } # replace T with a whitespace
  evla=${evt[2]//[[:space:]]/} # remove whitespaces
  evlo=${evt[3]//[[:space:]]/}
  evdp=${evt[4]//[[:space:]]/}

  evoepoch=$(date -ud "$evodate UTC" +%s)
  evojday=$(date -ud "$evodate UTC" +%Y,%j,%H:%M:%S.%N)

  #------ create event directories
  event_dir=$out_dir/$evid
  mkdir -p $event_dir
  cd $event_dir
  mkdir data mseed sacpz

  # save event info
  echo "$event_line" > $event_dir/data/event.txt

  #------ get available stations 
  echo
  echo "[LOG] ------ query seismic stations"
  echo
  station_list=$event_dir/data/station.txt
  # channel operating time range (start 1 day before and end 1 day after origin time)
  startbefore=$(date -ud "$evodate -1 days" +%Y-%m-%dT%H:%M:%S)
  endafter=$(date -ud "$evodate 1 days" +%Y-%m-%dT%H:%M:%S)
  str_twin="startbefore=${startbefore}&endafter=${endafter}"
  # query stations
  echo "[LOG] get station list"
  str_link="${fdsnws_station}?${str_twin}&${geo_range}&channel=FH?,CH?,HH?,BH?&level=station&format=text"
  echo "[LOG] wget $str_link -O $station_list"
  wget "$str_link" -O $station_list
  # check station number
  nsta=$(awk 'END{print NR}' $station_list)
  if [ $nsta -lt 1 ]; then
    echo "[WARN] No available stations! SKIP this event."
    continue
  fi
  echo "[LOG] stations found: $nsta"
  # query channels 
  channel_list=$event_dir/data/channel.txt
  echo "[LOG] get channel list"
  str_link="${fdsnws_station}?${str_twin}&${geo_range}&channel=FH?,CH?,HH?,BH?&level=channel&format=text"
  echo "[LOG] wget $str_link -O $channel_list"
  wget "$str_link" -O $channel_list

  #------ generate data selection list
  echo
  echo "[LOG] ------ make data selection list"
  echo
  dataselect_list=$event_dir/data/dataselect.txt
  cat /dev/null > $dataselect_list

  for station_line in $(grep -v "^#" $station_list)
  do
    # parse station info
    IFS='|' read -ra sta <<< "$station_line"
    netwk=${sta[0]}
    stnm=${sta[1]}
    stla=${sta[2]}
    stlo=${sta[3]}
    # get traveltime
    ttp=$(taup_time -evt $evla $evlo -h $evdp \
        -sta $stla $stlo --time -ph $phase_name | head -n1 | awk '{print $1}')
    if [ -z "$ttp" ]
    then
      echo "[ERROR] taup_time failed to get travaltime ttp=$ttp"
      exit -1
    fi
    # get time window
    t0=$(echo "$evoepoch + $twin_begin + $ttp" | bc -l)
    t1=$(echo "$evoepoch + $twin_end + $ttp" | bc -l)
    starttime=$(date -u -d @$t0 +%Y-%m-%dT%H:%M:%S)
    endtime=$(date -u -d @$t1 +%Y-%m-%dT%H:%M:%S)
    # print data selection list 
    printf "%-6s %-6s  *  BH?  %s  %s\n" \
      $netwk $stnm $starttime $endtime \
      >> $dataselect_list
  done

  #------ get waveform data
  echo 
  echo "[LOG] ------ get waveform data"
  echo
  mseed_file=$event_dir/mseed/${evid}.mseed
  # query waveforms 
  str_link="${fdsnws_dataselect}"
  echo "[LOG] wget --post-file=$dataselect_list $str_link -O $mseed_file"
  echo
  wget --post-file=$dataselect_list "$str_link" -O $mseed_file

  #------ get sacpz
  echo 
  echo "[LOG] ------ get sacpz"
  echo

  sacpz_dir=$event_dir/sacpz

  for line in $(grep -v "^#" $channel_list)
  do
    IFS='|' read -ra x <<< "$line"

    net="${x[0]}"
    sta="${x[1]}"
    loc="${x[2]}"
    cha="${x[3]}"

    IFS="$old_IFS" read -ra x <<< $(grep "${net}[ ]*${sta}" $dataselect_list)
    if [ "$?" -ne 0 ]
    then
      echo "[ERROR] ${net}.${sta} is not found in $dataselect_list! SKIP"
      continue
    fi
    starttime="${x[4]}"
    endtime="${x[5]}"

    echo
    echo "SACPZ ${net}.${sta}.${loc}.${cha} $starttime $endtime"
    echo
  
    url="${irisws_sacpz}?net=${net}&sta=${sta}&loc=${loc:---}&cha=${cha}&starttime=${starttime}&endtime=${endtime}"

    outfile="${sacpz_dir}/${net}.${sta}.${loc}.${cha}"
  
    wget "$url" -O $outfile
  
  done

done # for event_line in $(grep -v "^[ \t]*#" $event_list)

#END