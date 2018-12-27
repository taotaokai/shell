#!/bin/bash

selfdoc(){
cat<<EOF
NAME
    make_sacpz_IRF - make impulse response function(IRF) of a SACPZ file

SYNOPSIS
    make_sacpz_IRF [-o <out_dir>]

DESCRIPTION
    check whether IRF is close to a low-pass filtered delta function

    read files from STDIN.
    -o  if used, output IRF files.
    -f  specify sac filter command before comparing, default "lp co 2 n 2 p 2"

NOTE
    Here the sacpz file is assumed to transform from ground displacement to analog signal.
    And a broadband velocity seismometer should have a flat velocity response
    within a wide frequency range, like 120s~100Hz (e.g. STS-2, Trillium-120,)
    That means, if the input is a delta function, the output should be very close to
    a broad band-pass filtered first derivative of delta in both time and waveform.  
    When I compare IRF with the first derivative of delta, both are low-passed below 2 Hz.

EOF
}

#====== get command line args
filter="lp co 2 n 2 p 2"
flag_out=0 # flag if output IRF files
while getopts "o:f:h" opt
do
  case $opt in
    o)
      out_dir="$OPTARG"
      flag_out=1
      ;;
    f)
      filter="$OPTARG"
      ;;
    [h,?])
      selfdoc
      exit 1
      ;;
  esac
done

cat <<EOF
#====================
#Parameters:
#   -o = "$out_dir"
#   -f = "$filter"
#====================
EOF

#====== create delta pulse(centered at 50s) and its derivative
pulse=$(mktemp)
pulse_d1=$(mktemp)
sac<<EOF
fg d 0.001 n 10001
w $pulse
fft
mulomega
ifft
w $pulse_d1
q
EOF

#====== create response 
tmp_list=$(mktemp)
tmp_sacm=$(mktemp)
while read sacpz_file
do

    echo ====== ${sacpz_file}

    sacpz_fn=${sacpz_file##.*/}

    # make list of all epochs
    grep -F "* START" $sacpz_file > /dev/null
    if [ $? -ne 0 ]
    then
        # not annotated sacpz, make a fake list
        echo "net|sta|loc|cha|2000-01-01T00:00:00|2000-01-01T00:00:00|" > $tmp_list
    else
        # annotated sacpz
        sed "/^[\s]*$/d" $sacpz_file |\
            grep -A6 "* NETWORK" |\
            sed "/^\*/s/$/|/;s/^\*[^:]*://" |\
            awk 'BEGIN{RS="--"}{printf "%s%s%s%s%s%s\n", $1,$2,$3,$4,$6,$7}' \
                > $tmp_list
    fi

    # make IRF for each epoch
    while IFS="|" read -u 10 net sta loc cha starttime endtime
    do
        channel_id=${net}.${sta}.${loc}.${cha}

        # get valid date: starttime + 1 second
        sacdate=$(date -ud "${starttime/T/ } 1 seconds" \
            +"nzyear %Y nzjday %j nzhour %H nzmin %M nzsec %S nzmsec 0")

        # make sac macro
        cat<<EOF > $tmp_sacm
r $pulse
ch knetwk $net kstnm $sta khole "$loc" kcmpnm $cha
ch $sacdate
trans from none to polezero s $sacpz_file
EOF

        if [ "$flag_out" -eq 1 ]
        then
            echo "w dir $out_dir ${sacpz_fn}.${channel_id}.${starttime}.${endtime}.sac" \
                >> $tmp_sacm
        fi

        cat<<EOF >> $tmp_sacm
r more $pulse_d1
${filter}
correlate master 2 normalized
setbb max &1,depmax
setbb maxtime ( gettime max )
sc echo HED: $sacpz_fn $channel_id ${starttime} ${endtime} cc_max %max shift %maxtime
*ppk
q
EOF

        # run sac macro
        sac<<EOF
m $tmp_sacm
q
EOF

    done 10<$tmp_list

done

# clean tmp files
rm $pulse $pulse_d1 $tmp_list $tmp_sacm

#END