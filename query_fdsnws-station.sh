#!/bin/bash

#get station information from http://service.iris.edu/fdsnws/station/1/query?

function usage(){
cat <<EOF

SYNOPSIS
  $0 -R[geographic region] -T[time range] -C[channel name] 
     -L[output level] -F[output format] -O[output file]

DESCRIPTION
  -R [b/122/126/43.5/46]
     * box region: [b<west>/<east>/<south>/<north>], Example: b/-122/-121.5/46.8/46.9
     * radial region: [r<lon0>/<lat0>/<minRadius>/<maxRadius>], Radius(degree), Example: r/-120/40/1.0/5.0
  -T [1900-01-01/2999-12-12]
  -C [BHZ] others: BH?
  -L [station] others: network, station, channel, response
  -F [text] others: xml 
  -O [station]

EOF
  exit 1
}

# default parameter values
R=b/122/126/43.5/46
T=1900-01-01/2999-12-12
C="BHZ"
L="station"
F="text"
O="station"

# parse options
while getopts R:T:C:L:F:O:h name
do
    case $name in
    R)  R="$OPTARG";;
    T)  T="$OPTARG";;
    C)  C="$OPTARG";;
    L)  L="$OPTARG";;
    F)  F="$OPTARG";;
    O)  O="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

# validate input arguments
echo "#This run: $0 -T$T -R$R -C$C -O$O" | tee $O

# form the query links 
#Query Usage
#
#/query? [channel-options] [geographic-constraints] [time-constraints] [misc-parameters] [nodata=404]
#
#channel-options   ::  [network=<network>] [station=<station>] [location=<location>] [channel=<channel>]
#geographic-constraints  ::  [boundaries-rect] OR [boundaries-circular]
#boundaries-rect   ::  [minlatitude=<degrees>] [maxlatitude=<degrees>] [minlongitude=<degrees>] [maxlongitude=<degrees>]
#boundaries-circular   ::  [latitude=<latitude>&longitude=<longitude>&maxradius=<degrees>] [minradius=<degrees>]
#time-constraints  ::  [starttime=<date>] [endtime=<date>] [startbefore=<date>] [startafter=<date>]
#        [endbefore=<date>] [endafter=<date>]
#misc-parameters   ::  [level=<station|network|channel|response>] [format=<xml|text>] [includeavailability=<false|true>] 
#        [updatedafter=<date>] [matchtimeseries=<false|true>] [includerestricted=<true|false>]
#                                [includecomments=<false|true>] 

strHead="http://service.iris.edu/fdsnws/station/1/query?"

strR=$(echo $R | awk -F"/" '$1=="b"{printf "minlongitude=%s&maxlongitude=%s&minlatitude=%s&maxlatitude=%s",$2,$3,$4,$5}; $1=="r"{printf "longitude=%s&latitude=%s&minradius=%s&maxradius=%s",$2,$3,$4,$5}')

strT=$(echo $T | awk -F"/" '{printf "starttime=%s&endtime=%s",$1,$2}')

strC="channel=$C"

strL="level=$L"

strF="format=$F"

strLink="$strHead&$strR&$strT&$strC&$strL&$strF"

echo "#link: $strLink" | tee $O

wget $strLink -O - >> $O

echo "Finished. "

# END
