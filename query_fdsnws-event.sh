#!/bin/bash

#get earthquake catalog from http://service.iris.edu/fdsnws/event/1/query?
#USAGE: $0 -R lon0/lon1/lat0/lat1 -T start/end -M mag0/mag1

function usage(){
	echo ======================================
	echo "usage: $0 -R [search region]  -T <start-time>/<end-time> -M <minMag>/<maxMag> -O <file>"
	echo --------------------------------------
	echo "	-R Specify a box search region."
	echo "		box region: [b<west>/<east>/<south>/<north>], Example: b/-122/-121.5/46.8/46.9"
	echo "		radial region: [r<lon0>/<lat0>/<minRadius>/<maxRadius>], Radius(degree), Example: r/-120/40/1.0/5.0."
	echo "	-T Specify the temporal range. Example: 2012-11-29T00:00:00/2012-12-01T00:00:00, 2012-11-29 is equivalent to 2012-11-29T00:00:00"
	echo "	-M Specify the magnitude range. Example: -0.1/8.3"
	echo "	-O Specify the output file."
	echo ======================================
	echo "Default: $0 -R b/122/126/43.5/46 -T 2013-11-01/2013-12-01 -M -1.0/8.3 -O events"
	exit 1
}

# default parameter values
R=b/122/126/43.5/46
T=2013-11-01/2013-12-01
M=-1.0/8.3
O=events

# parse options
while getopts R:T:M:O:h name
do
    case $name in
    R)	R="$OPTARG";;
    T)	T="$OPTARG";;
    M)	M="$OPTARG";;
    O)	O="$OPTARG";;
		[h,?])	usage; exit -1
    esac
done

# validate input arguments

echo This run: $0 -T$T -M$M -R$R -O$O

# form the query links 
# keywords: [starttime,endtime,minlatitude,maxlatitude,minlongitude,maxlongitude,latitude,longitude,maxradius,minradius,mindepth,maxdepth,minmagnitude,maxmagnitude,magnitudetype]

strHead="http://service.iris.edu/fdsnws/event/1/query?"

strT=$(echo $T | awk -F"/" '{printf "starttime=%s&endtime=%s",$1,$2}')

strR=$(echo $R | awk -F"/" '$1=="b"{printf "minlongitude=%s&maxlongitude=%s&minlatitude=%s&maxlatitude=%s",$2,$3,$4,$5}; $1=="r"{printf "longitude=%s&latitude=%s&minradius=%s&maxradius=%s",$2,$3,$4,$5}')

strM=$(echo $M | awk -F"/" '{printf "minmagnitude=%s&maxmagnitude=%s",$1,$2}')

strLink="$strHead&$strT&$strM&$strR&format=text"

wget $strLink -O $O

# END
