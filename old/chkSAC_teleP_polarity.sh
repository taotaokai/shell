#!/bin/bash

# this script check if the polarity of the teleseismic P wave is consistent with the back-azimuth
# the ouput is the azimuthal deviation (in degree) of the horizontal polarity direction from the (back-azimuth-180)

# Note:
# the names of the sac data files of each teleseismic event should be in the format: [fsac].BHZ [fsac].BHE [fsac].BHN
# the P arrival time is assumed to be at time zero

# Usage:
# $ cat eventlist | bash ./check_polarity
# or  $ ls *.BHZ | sed 's/\.BHZ//' | bash ./check_polarity

while read fsac
do

sac<<EOF > /dev/null
r $fsac.BHZ $fsac.BHE $fsac.BHN
bp co 0.05 0.2 n 4 p 2
w z.sac e.sac n.sac
cut -10 20
r z.sac e.sac n.sac
correlate master 1
interpolate b 0 d 0.01
w over
cut -0.001 0.001
r z.sac e.sac n.sac
w over
q
EOF

baz=$(saclst baz f z.sac | awk '{print $2}')
Czz=$(saclst depmen f z.sac | awk '{print $2}')
Cez=$(saclst depmen f e.sac | awk '{print $2}')
Cnz=$(saclst depmen f n.sac | awk '{print $2}')

echo $fsac | awk -vbaz=$baz -vCzz=$Czz -vCez=$Cez -vCnz=$Cnz \
'{r2d=180/atan2(0,-1);pl=r2d*atan2(Cez/Czz,Cnz/Czz);az=(baz-180);if(pl<0)pl=pl+360;if(az<0)az=az+360;printf "%s polarity %6.1f radial %6.1f polarity-radial %6.1f\n",$1,pl,az,r2d*atan2(sin((pl-az)/r2d),cos((pl-az)/r2d))}'

done

rm z.sac e.sac n.sac
