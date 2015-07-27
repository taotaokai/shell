#!/bin/bash

# plot traces within a specified time window (reduced time)
# sac header must be set correctly.
# dependencies: SAC, saclst, GMT-5

usage(){
  echo
  echo "Syntax:"
  echo "  $0 -x -50:500 -n -50:500 -s 10 -f \"bp co 0.01 0.2 n 4 p 2\" -d./ -t Profile -l list -o profile.ps"
  echo
  echo "Parameters:"
  echo "  -x  plot time range(begin:end)."
  echo "  -n  time range to normalize amplitude(begin:end)."
  echo "  -s  slowness(s/deg) for reduced time."
  echo "  -a  align the traces at the specified SAC header, e.g. t1."
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
plot_range=-50:500
norm_range=$plot_range
slowness=10 # s/deg
align_sacHD=
SHIFT_ID=0
sac_filter="bp co 0.01 0.2 n 4 p 2"
sac_dir=./
title="Profile"
sac_list=list
ps=profile.ps

# parse options
while getopts x:n:s:a:f:d:t:l:o:h name
do
  case $name in
  x)  plot_range="$OPTARG";;
  n)  norm_range="$OPTARG";;
  s)  slowness="$OPTARG";;
  a)  align_sacHD="$OPTARG";;
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
  echo "[Error] $sac_list does not exist!"
  usage
fi

if [ ! -z "$align_sacHD" ]
then
  SHIFT_ID=1 # align time by sac header
else
  SHIFT_ID=0 # reduced time by linear moveout
fi

if [ $SHIFT_ID -eq 0 ]
then
  echo "#Command: $0 -x $plot_range -n $norm_range -s $s -f \"$sac_filter\" -d $sac_dir -t \"$title\" -l $sac_list -o $ps"
else
  echo "#Command: $0 -x $plot_range -n $norm_range -a $align_sacHD -f \"$sac_filter\" -d $sac_dir -t \"$title\" -l $sac_list -o $ps"
fi

plot_begin=${plot_range%:*}
plot_end=${plot_range#*:}

norm_begin=${norm_range%:*}
norm_end=${norm_range#*:}

#====== GMT plot

wkdir=$(pwd)

# temporary directory
tmp_dir=$(mktemp -d)

# get station metadata 
data_list=$tmp_dir/data_list

cd $sac_dir
if [ $SHIFT_ID -eq 0 ]
then
  saclst knetwk kstnm gcarc az o f $(cat $sac_list) | sort -k4 -n |\
    awk '{printf "%s %s %s %f %f %f\n",$1,$2,$3,$4,$5,-1*(v*$4 + $6)}' v=$slowness > $data_list 
else
  saclst knetwk kstnm gcarc az $align_sacHD f $(cat $sac_list) | sort -k4 -n |\
    awk '{printf "%s %s %s %f %f %f\n",$1,$2,$3,$4,$5,-1*$6}' > $data_list 
fi

cd $wkdir
ntrace=$(cat $data_list | wc -l)

ymin=0
ymax=$(( $ntrace + 1 ))

# make y-axis annotation 
awk '{printf "%d a %4.1f\n", NR,$4}' $data_list > $tmp_dir/y-annot.txt

# basemap
if [ $SHIFT_ID -eq 0 ]
then
  xlabel="time - $s * dist (s)"
else
  xlabel="time after $align_sacHD (s)"
fi

gmt psbasemap \
  -JX6i/9i \
  -BWSen+t"$title" \
  -By+l"dist (deg)" -Byc"$tmp_dir/y-annot.txt" \
  -Bxa+l"$xlabel" \
  -R"$plot_begin/$plot_end/$ymin/$ymax" \
  -P -K > $ps

#plot each trace
n=0
cat $data_list |\
while read sacfn knetwk kstnm gcarc az time_shift
do
  n=$(( $n + 1 ))

  # SAC cut plot range
sac<<EOF
r dir $sac_dir $sacfn
ch allt $time_shift
rtr;taper
$sac_filter
interpolate d 0.1
w dir $tmp_dir tmp.sac
cut ${plot_begin} ${plot_end}
r dir $tmp_dir tmp.sac
w alpha dir $tmp_dir sac.txt
q
EOF

  dt=$(sed "1q;d" ${tmp_dir}/sac.txt | awk '{print $1}')
  t0=$(sed "2q;d" ${tmp_dir}/sac.txt | awk '{print $1}')

  awk 'BEGIN{ORS=" "}; NR>30' $tmp_dir/sac.txt |\
    tr -s " " "\n" |\
    awk '{print t0+NR*dt, $1}' t0=$t0 dt=$dt > $tmp_dir/xy.txt

  yscale=$(awk '$1>=b&&$1<=e{if($2<0) $2=-1*$2; if($2>n) n=$2}; END{if(n==0) n=1; print 1/n}' b=$norm_begin e=$norm_end $tmp_dir/xy.txt)

  awk '{print $1, a*$2+y}' a=$yscale y=$n $tmp_dir/xy.txt |\
    gmt psxy -J -R -O -K -P >> $ps

  echo $plot_end $n ${knetwk}.${kstnm} $az |\
    awk '{printf "%f %f %s (%03d) ",$1,$2,$3,$4}' |\
    gmt pstext -J -R -F+f8,+jLM -D0.3c/0 -N -P -O -K >> $ps

done

# clean up
rm -rf $tmp_dir

echo $ps
