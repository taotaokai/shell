#!/bin/bash

#get earthquake catalog from http://service.iris.edu/fdsnws/event/1/query?
#USAGE: $0 -R lon0/lon1/lat0/lat1 -T start/end -M mag0/mag1

function usage(){
    echo ======================================
    echo "usage: $0 -R [search region]  -T <start-time>/<end-time> -M <minMag>/<maxMag> -O <file>"
    echo --------------------------------------
    echo "  -R Specify a box search region."
    echo "      box region: [b<south>/<north>/<west>/<east>], Example: b/-122/-121.5/46.8/46.9"
    echo "      radial region: [r<lat0>/<lon0>/<minRadius>/<maxRadius>], Radius(degree), Example: r/-120/40/1.0/5.0."
    echo "  -T Specify the temporal range. Example: 2012-11-29T00:00:00/2012-12-01T00:00:00, 2012-11-29 is equivalent to 2012-11-29T00:00:00"
    echo "  -M Specify the magnitude range. Example: -0.1/8.3"
    echo "  -D Specify depth range. Example: 100/1000"
    echo "  -C Specify catalog, e.g. ANF GCMT TEST ISC UofW NEIC PDE"
    echo "  -O Specify the output file."
    echo ======================================
    echo "Example: $0 -R b/122/126/43.5/46 -T 2013-11-01/2013-12-01 -M 5.5/8.3 -D 100/1000 -C GCMT -O events"
    exit 1
}

# default parameter values
R=b/20/60/90/150
T=2009-01-01
M=5.5
depth=100
catalog=GCMT
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

strHead="http://service.iris.edu/fdsnws/event/1/query?"

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

strC="catalog=$catalog"

strLink="${strHead}&${strT}&${strM}&${strR}&${strD}&${strC}&format=text"

wget $strLink -O $O

echo "Finished. "

# END
