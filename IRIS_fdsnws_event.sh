#!/bin/bash

usage(){
cat<<EOF
NAME 
    IRIS_fdsnws_event - a shell wrapper to IRIS web service to search event info.

SYNOPSIS
    IRIS_fdsnws_event -R [search region]  -T <start-time>/<end-time> -M <minMag>/<maxMag> -O <file>"

DESCRIPTION 

    -R Specify a box search region.  
        box region: [b/<south>/<north>/<west>/<east>] e.g. b/-122/-121.5/46.8/46.9
        radial region: [r/<lat0>/<lon0>/<minRadius>/<maxRadius>], Radius(degree), e.g. r/-120/40/1.0/5.0
    -T Specify the temporal range, e.g. 2012-11-29T00:00:00/2012-12-01T00:00:00
        2012-11-29 is equivalent to 2012-11-29T00:00:00
    -M Specify the magnitude range, e.g. -0.1/8.3
    -D Specify depth range, e.g. 100/1000
    -C Specify catalog, e.g. ANF GCMT TEST ISC UofW NEIC PDE
    -O Specify the output file name, defaults to stdout

EXAMPLE
    IRIS_fdsnws_event -R b/20/60/90/150 -T 2009-01-01 -M 5.5 -D 100 -C GCMT -O -

EOF
}

# base url
fdsnws_event="http://service.iris.edu/fdsnws/event/1/query"

# default parameter values
R=b/-90/90/-180/180
T=2009-01-01
M=5.5
depth=100
catalog=
O=-

# parse options
while getopts R:T:M:C:D:O:h name
do
    case $name in
    R)	R="$OPTARG";;
    T)	T="$OPTARG";;
    M)	M="$OPTARG";;
    D)	depth="$OPTARG";;
    C)	catalog="$OPTARG";;
    O)	O="$OPTARG";;
		[h,?])	usage; exit -1
    esac
done

# validate input arguments

#echo This run: $0 -T$T -M$M -R$R -O$O
echo "#This run: $0 -T$T -M$M -R$R -D$depth -C$catalog -O$O"

# form the query links 
# keywords: [starttime,endtime,minlatitude,maxlatitude,minlongitude,maxlongitude,latitude,longitude,maxradius,minradius,mindepth,maxdepth,minmagnitude,maxmagnitude,magnitudetype]

strT=$(echo $T | awk -F"/" '\
    NF==2{printf "starttime=%s&endtime=%s",$1,$2}; \
    NF==1{printf "starttime=%s",$1}')

strR=$(echo $R | awk -F"/" '\
    $1=="b"{printf \
        "minlatitude=%s&maxlatitude=%s&minlongitude=%s&maxlongitude=%s",\
        $2,$3,$4,$5}; \
    $1=="r"{printf \
        "latitude=%s&longitude=%s&minradius=%s&maxradius=%s",$2,$3,$4,$5}')

strM=$(echo $M | awk -F"/" '\
    NF==2{printf "minmagnitude=%s&maxmagnitude=%s",$1,$2}; \
    NF==1{printf "minmagnitude=%s",$1}')

strD=$(echo $depth | awk -F"/" '\
    NF==2{printf "mindepth=%s&maxdepth=%s",$1,$2}; \
    NF==1{printf "mindepth=%s",$1}')

if [ x"$catalog" != x ]
then
    strC="catalog=$catalog"
fi

strLink="${fdsnws_event}?${strT}&${strM}&${strR}&${strD}&${strC}&format=text"

echo "# $strLink"

wget $strLink -O $O

echo "Finished. "

# END
