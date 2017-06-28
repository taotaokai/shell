#!/bin/bash

# make impulse response function(IRF) of SACPZ files created by rdseed, 
# the file naming style must follow: SAC_PZs_SN_XANT_BHZ_00_2007.182.00.00.00.0000_2223.187.14.13.20.99999

#====== get command line args
out_dir=${1:-.}

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
while read sacpz_file; do

    sacpz_fn=${sacpz_file##.*/}

    # parse net/sta/loc/cha
    # use '--' instead of empty location id, because bash array use spaces as delimiter
    x=( $(sed 's/__/_--_/g; s/_/ /g' <<< "$sacpz_fn") )
    net=${x[2]}
    sta=${x[3]}
    cha=${x[4]}
    loc=${x[5]/--/}

    seed_id=${net}.${sta}.${loc}.${cha}

    # get valid date: start time + 1 second
    d=( ${x[6]//./ } )
    yr=${d[0]}
    jd=${d[1]}
    hr=${d[2]}
    mi=${d[3]}
    se=${d[4]}
    sacdate=$(date -u -d "$yr-01-01 00:00:00 -1 days \
        $jd days $hr hours $mi minutes $se seconds +1 seconds" \
        +"nzyear %Y nzjday %j nzhour %H nzmin %M nzsec %S nzmsec 0")

sac<<EOF
r $pulse
ch knetwk $net kstnm $sta khole "$loc" kcmpnm $cha
ch $sacdate
trans from none to polezero s $sacpz_file
w dir $out_dir ${sacpz_fn}.FIR.sac
r more $pulse_d1
lp co 2 n 4 p 2
correlate master 2 normalized
setbb max &1,depmax
setbb maxtime ( gettime max )
sc echo HED: $seed_id $sacpz_fn cc_max %max shift %maxtime
*ppk
q
EOF

done

rm $pulse $pulse_d1