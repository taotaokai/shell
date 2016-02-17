#!/bin/bash



usage(){
cat<<EOF
NAME 
    IRIS_fdsnws_station - a shell wrapper to IRIS web service to search seismic stations info.

SYNOPSIS
    IRIS_fdsnws_station -R [search region]  -T <start-time>/<end-time> -M <minMag>/<maxMag> -O <file>"

DESCRIPTION 

    -R Specify a box search region.  
        box region: [b<south>/<north>/<west>/<east>] e.g. b/-122/-121.5/46.8/46.9
        radial region: [r<lat0>/<lon0>/<minRadius>/<maxRadius>], Radius(degree), e.g. r/-120/40/1.0/5.0
    -T Specify the temporal range, e.g. 2012-11-29T00:00:00/2012-12-01T00:00:00
        2012-11-29 is equivalent to 2012-11-29T00:00:00
    -L Specify level of detail using network, station, channel.
    -O Specify the output file name, defaults to stdout

NOTE
  1) only [FCHB]H? channels are retrived (broadband high gain seismometer at decreasing sampling rate: F/C/H/B)

EXAMPLE
    IRIS_fdsnws_station -R b/20/60/90/150 -T 2009-01-01 -O -

EOF
}

# base url
fdsnws_station="http://service.iris.edu/fdsnws/station/1/query"

# default parameter values
R=b/-90/90/-180/180
T=2009-01-01
L=station
O=-

# parse options
while getopts R:T:L:O:h name
do
    case $name in
    R)	R="$OPTARG";;
    T)	T="$OPTARG";;
    L)	L="$OPTARG";;
    O)	O="$OPTARG";;
		[h,?])	usage; exit -1
    esac
done

# validate input arguments

echo "# Command: $0 -R$R -T$T -L$L -O$O"

# form the query links 
# keywords: [starttime,endtime,minlatitude,maxlatitude,minlongitude,maxlongitude,latitude,longitude,maxradius,minradius,mindepth,maxdepth,minmagnitude,maxmagnitude,magnitudetype]

strR=$(echo $R | awk -F"/" '\
    $1=="b"{printf \
        "minlatitude=%s&maxlatitude=%s&minlongitude=%s&maxlongitude=%s",\
        $2,$3,$4,$5}; \
    $1=="r"{printf \
        "latitude=%s&longitude=%s&minradius=%s&maxradius=%s",$2,$3,$4,$5}')

strT=$(echo $T | awk -F"/" '\
    NF==2{printf "starttime=%s&endtime=%s",$1,$2}; \
    NF==1{printf "starttime=%s",$1}')

strL="level=$L"

strLink="${fdsnws_station}?${strR}&${strT}&${strL}&channel=FH?,CH?,HH?,BH?&format=text"

echo "# Link: $strLink"

wget $strLink -O $O

echo "Finished. "

# END
