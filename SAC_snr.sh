#!/bin/bash

# calculate SNR of SAC files
# SNR = A_signal / A_noise, where A_[] is the maximum amplitude within the
# corresponding time window

usage(){
  echo
  echo "Syntax:"
  echo "  $0 -n t0/-45/-15 -s t0/-10/20 -f \"bp co 0.01 0.2 n 4 p 2\" -l-"
  echo
  echo "Parameters:"
  echo "  -n  noise time window: tmarker/begin/end"
  echo "  -s  signal time window: tmarker/begin/end"
  echo "  -f  sac filter command."
  echo "  -l  sac file list. If not exist, read from STDIN"
  echo 
  exit -1
}

#====== read parameters

#defaults
tnoise=t0/-45/-15
tsignal=t0/-10/20
sac_filter="bp co 0.01 0.2 n 4 p 2"
sac_list=-

# parse options
while getopts n:s:f:l:h name
do
    case $name in
    n)  tnoise="$OPTARG";;
    s)  tsignal="$OPTARG";;
    f)  sac_filter="$OPTARG";;
    l)  sac_list="$OPTARG";;
    [h,?])  usage; exit -1
    esac
done

echo "#===== Parameters:"
echo "# $0 -n $tnoise -s $tsignal -f \"$sac_filter\" -l $sac_list"

tmark_noise=$(echo $tnoise | awk -F"/" '{print $1}')
noise_b=$(echo $tnoise | awk -F"/" '{print $2}')
noise_e=$(echo $tnoise | awk -F"/" '{print $3}')

tmark_signal=$(echo $tsignal | awk -F"/" '{print $1}')
signal_b=$(echo $tsignal | awk -F"/" '{print $2}')
signal_e=$(echo $tsignal | awk -F"/" '{print $3}')

#===== calculate SNR for each sac file

tmp_sac=$(mktemp)
tmp_txt=$(mktemp)

cat $sac_list |\
while read sacfn
do

sac<<EOF

r ${sacfn}
rtr;taper
${sac_filter}
abs
w $tmp_sac

cut $tmark_noise $noise_b $tmark $noise_e
r $tmp_sac
setbb amp0 &1,depmax

cut $tmark_signal $signal_b $tmark $signal_e
r $tmp_sac
sc echo %amp0 &1,depmax > $tmp_txt

q
EOF

  # calculate SNR
  cat $tmp_txt |\
    awk '{printf "%-16s | noise %s %5.1f %5.1f %9.2e | signal %s %5.1f %5.1f %9.2e | SNR %9.2e\n", \
          fn, tn,t0n,t1n,$1, ts,t0s,t1s,$2, $2/$1}' fn=$sacfn \
          tn=$tmark_noise t0n=$noise_b t1n=$noise_e \
          ts=$tmark_signal t0s=$signal_b t1s=$signal_e

done

#===== clean up
rm $tmp_sac $tmp_txt