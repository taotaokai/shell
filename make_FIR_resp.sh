#!/bin/bash

# make response from the FIR stages in RESP files

# create delta pulse and its derivative
sac<<EOF
fg d 0.01 n 10001
w pulse
fft
mulomega
ifft
w pulse_d1
q
EOF

#  
tmp_file=$(mktemp)
while read resp_file; do

  # remove Poles&Zeros in the RESP file
    cat $resp_file > $tmp_file
    sed -i "/^B053F07/s/:.*/:   1.0/" $tmp_file
    sed -i "/^B053F08/s/:.*/:   1.0/" $tmp_file
    sed -i "/^B053F09/s/:.*/:   0/"   $tmp_file
    sed -i "/^B053F14/s/:.*/:   0/"   $tmp_file
    sed -i "/^B053F10-13/s/^/#/"      $tmp_file
    sed -i "/^B053F15-18/s/^/#/"      $tmp_file

sac<<EOF
r pulse
trans from none to evalresp fname $tmp_file
w ${resp_file}.FIR.sac
r more pulse_d1
correlate master 2
setbb maxtime ( gettime max )
sc echo FIR_SHIFT: $resp_file %maxtime
q
EOF

done

rm $tmp_file