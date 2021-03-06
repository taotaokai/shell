#!/bin/bash

usage(){
  cat<<EOF
Usage 
-----
$0 -L <SAC_List>  -D <SAC_DIR>

Default: $0 [-L-] -D.

Note
----
evdp(km) and gcarc(degree) should be correctly set in SAC header!

P: user0: rayp(s/km); user1: rayp(s/deg)
S: user2: rayp(s/km); user3: rayp(s/deg)

EOF

}

# Default parameters

saclst=-
sacdir=.

# parse options
while getopts L:D:P:h name
do
  case $name in
  L) saclst="$OPTARG";;
  D) sacdir="$OPTARG";;
  [h,?]) usage; exit -1
  esac
done

# loop each sac file
for sacf in $(cat $saclst)
do
  stninf=$(saclst gcarc evdp f $sacdir/$sacf)
  gcarc=$(echo $stninf | awk '{print $2}')
  evdp=$(echo $stninf | awk '{print $3}')
  
  # for P
  raypdeg_P=$(taup_time -h $evdp -deg $gcarc -ph P -rayp | sed -n "1p")
  raypkm_P=$(echo "scale=12;$raypdeg_P * 180 / 6371 / 4 / a(1)" | bc -l)

  raypdeg_S=$(taup_time -h $evdp -deg $gcarc -ph S -rayp | sed -n "1p")
  raypkm_S=$(echo "scale=12;$raypdeg_S * 180 / 6371 / 4 / a(1)" | bc -l)

  echo "===== $sacf gcarc $gcarc evdp $evdp"
  echo "----- P rayp(s/deg) $raypdeg_P rayp(s/km) $raypkm_P"
  echo "----- S rayp(s/deg) $raypdeg_S rayp(s/km) $raypkm_S"

sac<<EOF > /dev/null
r $sacdir/$sacf
ch user0 $raypkm_P user1 $raypdeg_P
ch user2 $raypkm_S user3 $raypdeg_S
wh
q
EOF

done
