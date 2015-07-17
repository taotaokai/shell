#!/bin/bash

#History:
#[20140614] STN_List (stla stlo) -> (stlo stla)
#[20140615] cat $STALST | sed '/^\s*$/d;/^#/d' |\ 
#[20140617] #merge sac: awk END{(if NR>=1)} in case no file exits
#[20140624] default -K.BHZ/.BHN/.BHE
#[20140624] output sac is named as ${gmt_btime}.${netwk}.${stnm}${cmpnm}
#[20140718] change order to "stlo stla"
#[20140719] 
#  the -K option now assume component name without preceding dot 
#  (e.g -K BHE/BHN/BHZ)
#[20140725] 
#  use bash build-in ${string:position:length} to cut the 
#  gmt time string instead of 
#  $(echo $gmt | cut -c1-4,5-6,7-8 --output-delimiter "-"), 
#  because sometimes the later won't add the output delimiter, 
#  which is strange and I cannot find out why.
#[20140726] used file descriptor to safely read lines from "while read -r"
#[20140726] change the time specification in CUT_LIST to "ref_time begin end"
#[20140726] remove the default behavior to read from stdin if -S not specified 
#[20140928] change the self document

selfdoc(){
cat <<EOF
NAME
  SACdb_wincut

DESCRIPTION
  extract data in a given time window from a SAC database

SYNOPSIS
  SACdb_wincut -C CUT.lst -S STA.lst -I SAC.index
               -K BHZ/BHN/BHE -D . -O .

PARAMETERS

  -C<CUT_List>  text file with each line consisting of 
    "ref_time begin end"
    ref_time:  yyyymmdd_HHMMSS(.SSS) and timezone UTC+00:00 is assumed
    begin/end are seconds relative to ref_time"

  -S<STA_List>  text file with each line consisting of 
    "netwk stnm stlo stla stel"

  -I<INDEX_File>  text file with each line consisting of
    "sacfn begin-time end-time"
    begin-time and end-time should be the epoch time (nanoseconds) 
    equivalent to the output of \$(date -d "date string" +%s%N)

  -K<CMP1>[/<CMP2>/<CMP3>] 
    slash separated strings to identify different components, 
    such as .Z/.N/.E, each component string will be used with 
    grep -F (treat as a fixed string!) to find out files of the same component 
    during the merging SAC process

  -D<DATA_DIR> SAC data directory 
    sac files should be stored as <DATA_DIR>/{stnm}/*.sac"

  -O<OUT_DIR> the output file path 
    {OUT_DIR}/{stnm}/{ref_time}.{netwk}.{stnm}.{CMP?}"

DEFAULT
  SACdb_wincut -C CUT.lst -S STA.lst -I SAC.index -K BHZ/BHN/BHE -D. -O.

EOF
  exit -1
}

# Default parameters

CUTLST=CUT.lst
STALST=STA.lst
CMPSTR="BHE/BHN/BHZ"
INDEX=SAC.index
OUTDIR=.
DATDIR=.

# parse options
while getopts C:S:K:I:D:O:h name
do
    case $name in
    C)  CUTLST="$OPTARG";;
    S)  STALST="$OPTARG";;
    K)  CMPSTR="$OPTARG";;
    I)  INDEX="$OPTARG";;
    D)  DATDIR="$OPTARG";;
    O)  OUTDIR="$OPTARG";;
    [h,?])  selfdoc;;
    esac
done

if [ ! -f "$CUTLST" ]
then
  echo "Error: $CUTLST does not exist!"
  selfdoc 
fi

if [ ! -f "$STALST" ]
then
  echo "Error: $STALST does not exist!"
  selfdoc
fi

if [ ! -f "$INDEX" ]
then
  echo "Warning: INDEX file [$INDEX] does not exist, use $DATDIR/SAC.index"
  INDEX=$DATDIR/SAC.index
fi

if [ ! -f "$INDEX" ]
then
  echo "Error: cannot find $DATDIR/SAC.index!"
  selfdoc
fi

echo "This run: $0 -C $CUTLST -S $STALST -K $CMPSTR -I $INDEX -D $DATDIR -O $OUTDIR"

# loop for each station
tmp_lst=$(mktemp -p. tmp.lst.XXX)
tmp_sacm=$(mktemp -p. tmp.sacm.XXX)
tmp_sta=$(mktemp -p. tmp.sta.XXX)
tmp_cut=$(mktemp -p. tmp.cut.XXX)

# sed: remove blank and commented lines ([20140615])
sed '/^\s*$/d;/^#/d' $STALST > $tmp_sta
sed '/^\s*$/d;/^#/d' $CUTLST > $tmp_cut

while read -r -u3 netwk stnm stlo stla stel tmp
do

  if [ ! -d $OUTDIR/$stnm ]
  then
      mkdir -p $OUTDIR/$stnm
  fi

  while read -r -u4 reftime begin end tmp
  do
  
    # begin/end cut time in epoch (nanoseconds)
    ymd=${reftime:0:4}-${reftime:4:2}-${reftime:6:2}
    hms=${reftime:9:2}:${reftime:11:2}:${reftime:13}
    otime=$(date -u -d "${ymd} ${hms}" +%s%N )

    btime=$(echo $otime $begin | awk '{print $1+10^9*$2}')
    etime=$(echo $otime $end | awk '{print $1+10^9*$2}')

    echo ===================
    echo Date: $(date)
    echo CUT: $reftime $begin $end
    echo STN: $netwk $stnm $stlo $stla $stel

    # find sac files that cover time range [btime,etime] in the INDEX file
    grep "$stnm" $INDEX | awk '$3>=b && $2<=e {print dir"/"$0}' \
      dir=$DATDIR e=$etime b=$btime | sort > $tmp_lst
    if [ ! -s $tmp_lst ]; then
      echo Warning: No file found!
      continue
    fi
    echo MERGE:
    cat $tmp_lst

    # generate sac macro
    cat /dev/null > $tmp_sacm # the content from the last loop must be cleaned

    #merge sac
    sacfn=$(mktemp -p. -u)
    for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
    do
      cat $tmp_lst | grep -F "$cmpnm" |\
      awk 'NR==1{printf "r %s\n",$1}
           NR>=2{printf "merge overlap average %s\n",$1}
           END{if(NR>=1) printf "ch kcmpnm %s;w %s.%s\n",cmpnm,sacfn,cmpnm}' \
           sacfn="$sacfn" cmpnm="$cmpnm" >> $tmp_sacm
    done

    #time cut and write event info
    tmp=$(echo "$otime/10^9" | bc -l)
    gmt_otime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S" "%N | awk '{printf "%s %s %s %s %s %3d",$1,$2,$3,$4,$5,$6/10^6}')
    tmp=$(echo "$btime/10^9" | bc -l)
    gmt_btime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S" "%N | awk '{printf "%s %s %s %s %s %3d",$1,$2,$3,$4,$5,$6/10^6}')
    tmp=$(echo "$etime/10^9" | bc -l)
    gmt_etime=$(date -u -d @$tmp +%Y" "%j" "%H" "%M" "%S" "%N | awk '{printf "%s %s %s %s %s %3d",$1,$2,$3,$4,$5,$6/10^6}')

    echo CUT: B $gmt_btime E $gmt_etime

    echo "r $sacfn*" >> $tmp_sacm
    echo "ch o  gmt $gmt_otime" >> $tmp_sacm
    echo "ch t1 gmt $gmt_btime" >> $tmp_sacm
    echo "ch t2 gmt $gmt_etime" >> $tmp_sacm
    echo "wh" >> $tmp_sacm
    echo "cut t1 0 t2 0" >> $tmp_sacm
    echo "r $sacfn*" >> $tmp_sacm
    echo "ch kevnm $reftime knetwk $netwk kstnm $stnm stla $stla stlo $stlo stel $stel" >> $tmp_sacm
    echo "w over" >> $tmp_sacm

    #set origin time to zero
    echo "cut off" >> $tmp_sacm
    for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
    do
      echo "r ${sacfn}.${cmpnm}" >> $tmp_sacm
      echo "eval to v -1 * &1,o" >> $tmp_sacm
      echo "ch allt %v" >> $tmp_sacm
      echo "wh" >> $tmp_sacm
    done

sac <<EOF
m $tmp_sacm
q
EOF

    # rename sac and move into the output directory
    for cmpnm in $(echo $CMPSTR | sed "s/\// /g")
    do
      #${reftime:0:15}: trim reftime to yyyymmdd_HHMMSS
      echo OUT: mv ${sacfn}.${cmpnm} $OUTDIR/$stnm/${reftime:0:15}.${netwk}.${stnm}.${cmpnm}
      mv ${sacfn}.${cmpnm} $OUTDIR/$stnm/${reftime:0:15}.${netwk}.${stnm}.${cmpnm}
    done

  done 4<$tmp_cut
done 3<$tmp_sta

# delete tmp files
rm  $tmp_lst $tmp_sacm $tmp_cut $tmp_sta

# END of script
