#!/bin/bash

selfdoc(){
cat<<EOF
NAME
    make_sacpz_IRF - make impulse response function(IRF) of a SACPZ file

SYNOPSIS
    make_sacpz_IRF [-o <out_dir>]

DESCRIPTION
    read files from STDIN.
    -o  if used, output IRF files.
    check whether IRF is close to a low-pass filtered delta function
    a broadband velocity seismometer should have a flat velocity response
    within a wide frequency range, like 120s~100Hz (e.g. STS-2, Trillium-120,)
    Here, I only low-pass below 2 Hz.
EOF
}

#====== get command line args
flag_out=0 # flag if output IRF files
while getopts "o:h" opt
do
    case $opt in
        o)
            out_dir="$OPTARG"
            flag_out=1
            ;;
        [h,?])
            selfdoc
            exit 1
            ;;
    esac
done

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

    echo ${sacpz_file}

    sacpz_fn=${sacpz_file##.*/}

    # make list of all epochs
    grep -F "* START" $sacpz_file > /dev/null
    if [ $? -ne 0 ]
    then
        # not annotated sacpz, make a fake list
        echo "net|sta|loc|cha|2000-01-01T00:00:00|2000-01-01T00:00:00|" > $tmp_list
    else
        # annotated sacpz
        sed "/^[\s]*$/d" IU.INCN.00.BHZ |\
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
lp co 2 n 4 p 2
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