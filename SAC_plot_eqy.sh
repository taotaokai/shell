#!/bin/bash

# plot traces within a specified time window (reduced time)
# sac header must be set correctly.
# dependencies: SAC, saclst, GMT-5

usage(){
  echo
  echo "Syntax:"
  echo "  $0 -x o/-50:500 -n \"\" -s 0 -c -1 -f \"\" -l list -p \"\" -t \"\" -o profile.ps"
  echo
  echo "Parameters:"
  echo "  -x  time range(second) of plot (reference/begin:end). "
  echo "  -n  time range(second) to calculate the normalization amplitude (reference/begin:end). "
  echo "      If not specified the whole trace is normalized. "
  echo "  -s  slowness(s/deg) to calculate the reduce time. "
  echo "  -c  clip amplitude (no clip if < 0). "
  echo "  -f  SAC command to manipulate the data (e.g. bp co 0.02 0.2 n 4 p 2) "
  echo "  -l  list of SAC files to plot. "
  echo "  -p  prefix to the file names in file -l. "
  echo "  -t  title string of the plot. "
  echo "  -o  output PS file. "
  echo
  echo "Description:"
  exit -1
}

#====== read parameters

#defaults
plot_xlim=o/-50:500
norm_xlim=
slowness=0
clip=1
sac_macro=
sac_list=list
prefix=
plot_title="TITLE"
ps=profile.ps

# parse options
while getopts x:n:s:c:f:l:p:t:o:h name
do
  case $name in
  x)  plot_xlim="$OPTARG";;
  n)  norm_xlim="$OPTARG";;
  s)  slowness="$OPTARG";;
  c)  clip="$OPTARG";;
  f)  sac_macro="$OPTARG";;
  l)  sac_list="$OPTARG";;
  p)  prefix="$OPTARG";;
  t)  plot_title="$OPTARG";;
  o)  ps="$OPTARG";;
  [h,?])  usage; exit -1
  esac
done

if [ ! -f "$sac_list" ]
then
  echo "[Error] $sac_list does not exist!"
  usage
fi

echo "#Command: $0 -x $plot_xlim -n $norm_xlim -s $slowness -c $clip -f \"$sac_macro\" -l $sac_list -p \"$prefix\" -t \"$plot_title\" -o $ps"

# further parse arguments
plot_ref=${plot_xlim%/*}
tmpvar=${plot_xlim#*/}
plot_b=${tmpvar%:*}
plot_e=${tmpvar#*:}

if [ ! -z "$norm_xlim" ]
then
  NORM_WHOLE_TRACE=0
  norm_ref=${norm_xlim%/*}
  tmpvar=${norm_xlim#*/}
  norm_b=${tmpvar%:*}
  norm_e=${tmpvar#*:}
else
  NORM_WHOLE_TRACE=1
  norm_ref=o
  norm_b=0
  norm_e=0
fi

#====== GMT plot

# temporary directory
tmp_dir=$(mktemp -d)

# get station metadata 
data_list=$tmp_dir/data_list

saclst knetwk kstnm gcarc az $plot_ref $norm_ref \
  f $(awk '{printf "%s%s\n",v,$0}' v="$prefix" $sac_list) | sort -k4 -n |\
  awk '{printf "%s %s %s %f %f %f %f\n",$1,$2,$3,$4,$5,-1*(v*$4+$6),$7}' v=$slowness > $data_list 

ntrace=$(cat $data_list | wc -l)
ymin=-0.5
ymax=$(echo $ntrace 1.5 | awk '{print $1+$2}')

# make y-axis annotation 
awk '{printf "%d a %4.1f\n", NR,$4}' $data_list > $tmp_dir/y-annot.txt

# basemap
xlabel="time after $plot_ref - ${slowness}*dist (s)"

gmt psbasemap \
  -JX6i/9i \
  -BWSen+t"$plot_title" \
  -By+l"dist (deg)" -Byc"$tmp_dir/y-annot.txt" \
  -Bxa+l"$xlabel" \
  -R"$plot_b/$plot_e/$ymin/$ymax" \
  -P -K > $ps

#plot each trace
n=0
cat $data_list |\
while read sacfn knetwk kstnm gcarc az time_shift norm_t0
do
  n=$(( $n + 1 ))

  # SAC cut plot range
sac<<EOF
r $sacfn
ch allt $time_shift
rtr;taper
$sac_macro
interpolate d 0.1
w dir $tmp_dir tmp.sac
cut ${plot_b} ${plot_e}
r dir $tmp_dir tmp.sac
w alpha dir $tmp_dir sac.txt
q
EOF

  dt=$(sed "1q;d" ${tmp_dir}/sac.txt | awk '{print $1}')
  t0=$(sed "2q;d" ${tmp_dir}/sac.txt | awk '{print $1}')

  awk 'BEGIN{ORS=" "}; NR>30' $tmp_dir/sac.txt |\
    tr -s " " "\n" |\
    awk '{print t0+NR*dt, $1}' t0=$t0 dt=$dt > $tmp_dir/xy.txt

  if [ $NORM_WHOLE_TRACE -eq 0 ]
  then
    yscale=$(awk 'BEGIN{o=t0-v;b=b+o;e=e+o} \
      $1>=b&&$1<=e{if($2<0) $2=-1*$2; if($2>n) n=$2} \
      END{if(n==0) n=1; print 0.5/n}' \
      t0=$norm_t0 v=$time_shift b=$norm_b e=$norm_e $tmp_dir/xy.txt)
  else
    yscale=$(awk '{if($2<0) $2=-1*$2; if($2>n) n=$2} \
      END{if(n==0) n=1; print 0.5/n}' $tmp_dir/xy.txt)
  fi

  awk '{b=a*$2; if(c>0&&(b>c||b<-c)){print $1,"NAN"} else{print $1,b+y}}' \
    a=$yscale y=$n c=$clip $tmp_dir/xy.txt |\
    gmt psxy -J -R -O -K -P >> $ps

  echo $plot_e $n ${knetwk}.${kstnm} $az |\
    awk '{printf "%f %f %s (%03d) ",$1,$2,$3,$4}' |\
    gmt pstext -J -R -F+f8,+jLM -D0.3c/0 -N -P -O -K >> $ps

done

# clean up
rm -rf $tmp_dir

echo $ps
