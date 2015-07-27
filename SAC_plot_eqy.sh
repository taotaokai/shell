#!/bin/bash

# plot traces within a specified time window (reduced time)
# sac header must be set correctly.
# dependencies: SAC, saclst, GMT-5

usage(){
  echo
  echo "Syntax:"
  echo "  $0 -b -50 -e 500 -s 10 -f \"bp co 0.01 0.2 n 4 p 2\" -d./ -t Profile -l list -o profile.ps"
  echo
  echo "Parameters:"
  echo "  -b/-e  begin and end time."
  echo "  -s  slowness for reduced time."
  echo "  -f  sac filter command."
  echo "  -d  sac data directory."
  echo "  -t  string of title."
  echo "  -l  list of sac files."
  echo "  -o  output PS file."
  echo
  exit -1
}

#====== read parameters

#defaults
b=-50
e=500
s=10
sac_filter="bp co 0.01 0.2 n 4 p 2"
sac_dir=./
title="Profile"
sac_list=list
ps=profile.ps

# parse options
while getopts b:e:s:f:d:t:l:o:h name
do
    case $name in
    b)  b="$OPTARG";;
    e)  e="$OPTARG";;
    s)  s="$OPTARG";;
    f)  sac_filter="$OPTARG";;
    d)  sac_dir="$OPTARG";;
    t)  title="$OPTARG";;
    l)  sac_list="$OPTARG";;
    o)  ps="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

if [ ! -f "$sac_list" ]
then
   echo "Error: $sac_list does not exist!"
   usage
fi

echo "#===== Parameters:"
echo "# $0 -b $b -e $e -s $s -f \"$sac_filter\" -d $sac_dir -t \"$title\" -l $sac_list -o $ps"

#====== cut sac, convert to txt format

wkdir=$(pwd)

# temporary directory
tmp_dir=$(mktemp -d)

# get station metadata 
data_list=$tmp_dir/data_list
cd $sac_dir
saclst gcarc o az knetwk kstnm f $(cat $sac_list) | sort -n -k2 > $data_list 
cd $wkdir
ntrace=$(cat $data_list | wc -l)

# cut each trace
n=0
cat $data_list |\
while read sacfn gcarc ot tmp
do

  n=$(( $n + 1 ))
  reduce_time=$(echo $s $gcarc $ot | awk '{print -1*($1*$2 + $3)}')

sac<<EOF
r dir $sac_dir $sacfn
ch allt $reduce_time
rtr;taper
$sac_filter
interpolate d 0.1
w dir $tmp_dir tmp.sac
cut $b $e
r dir $tmp_dir tmp.sac
w alpha dir $tmp_dir $n.txt
q
EOF

done

#====== GMT plot
ymin=0
ymax=$(( $ntrace + 1 ))

# make y-axis annotation 
awk '{printf "%d a %4.1f\n", NR,$2}' $data_list > $tmp_dir/y-annot.txt

# basemap
gmt psbasemap \
  -JX6i/9i \
  -BWSen+t"$title" \
  -By+l"dist (deg)" -Byc"$tmp_dir/y-annot.txt" \
  -Bxa+l"time - $s * dist (s)" \
  -R"$b/$e/$ymin/$ymax" \
  -P -K > $ps

#plot each trace
n=0
cat $data_list |\
while read sacfn gcarc ot az knetwk kstnm
do

  n=$(( $n + 1 ))

  dt=$(sed "1q;d" $tmp_dir/${n}.txt | awk '{print $1}')
  t0=$(sed "2q;d" $tmp_dir/${n}.txt | awk '{print $1}')
  maxamp=$(sed "1q;d" $tmp_dir/${n}.txt |\
    awk '{a=$2;b=$3;if(a<0)a=-a;if(b<0)b=-b;if(b<a)b=a; print b}')

  yscale=$(echo $maxamp | awk '{print 1/$1}')

  awk 'BEGIN{ORS=" "}; NR>30' $tmp_dir/${n}.txt | tr -s " " "\n" |\
    awk '{print t0+NR*dt, a*$1+y}' t0=$t0 dt=$dt a=$yscale y=$n |\
    gmt psxy -J -R -O -K -P >> $ps

   echo $e $n ${knetwk}.${kstnm} $az |\
    awk '{printf "%f %f %s (%03d) ",$1,$2,$3,$4}' |\
    gmt pstext -J -R -F+f8,+jLM -D0.3c/0 -N -P -O -K >> $ps

done

# clean up
rm -rf $tmp_dir

echo $ps
