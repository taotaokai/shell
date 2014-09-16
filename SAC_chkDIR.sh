#!/bin/bash

#Description
#	check SAC data under every directory {dirnm} listed in an index file:
# 	the paths of radial RFs are assumed to be {dirnm}/c_rfa/*.BHR
#   the direct P should be shifted to zero time
#		markers: bad(t5<0); needHPfilter(t6<0); reversedP(t7<0)
#	index file: dirnm ( | checked/skipped $(date) )
#		the index file must be supplied, and the first column is the directory name;
#		this program will append ( | checked/skipped $(date) ) to the end of each line after the corresponding directory is checked
#History
#	[20140721]tao: created

function usage(){
	echo "Description:"
	echo "	check (modify sac header via ppk) SAC data under every directory {dirnm} listed in an index file"
	echo "	the target phase (e.g. direct P) should be shifted to zero time"
	echo "	the program calls SAC and plot seismograms using ppk, the user can mark time markers on the seismogram for a certain purpos"
	echo "	example of time marker: bad(t5<0); needHPfilter(t6<0); reversedP(t7<0), you can used your only definition"
	echo
	echo "Usage: $0 -D<datDIR> -I<index> -S<subDIR> -R<sacfn> -M<sacMACRO>"
	echo 
	echo "	-D Specify the parent directory"
	echo "	-I Specify the index file, which is a list of <dirnm>"
	echo "	the index file must be supplied, and the first column is the directory name <dirnm>"
	echo "	this program will append ( | checked/skipped {date} ) to the end of each line after the corresponding directory is checked"
	echo "	-S Specify the subdirectory for SAC data, as will be {DIRNM}/{SUBDIR}/{CMPNM}"
	echo "	-R Specification of sac files to read, e.g. \"*.BH[ENZ]\" (quote the wildcard to avoid shell expansion)"
	echo "	the path of SAC data is <datDIR>/<dirnm>/<subDIR>/<sacfn>"
	echo "	-M Sepcify the SAC commands to execute after read in the sac data. Example: 'rmean;rtrend;taper;bp co 2 10 n 4 p2;'"
	echo 
	echo "Default: $0 -D. -I chk.index -S'' -R '*.BH[ENZ]' -M''" 
	exit -1
}

# Default parameters
INDEX=chk.index
DATDIR='.'
SUBDIR=''
SACNM='*.BH[ENZ]'
MACRO=''

# parse options
while getopts D:I:S:R:M:h name
do
    case $name in
    I)  INDEX="$OPTARG";;
    D)  DATDIR="$OPTARG";;
    S)  SUBDIR="$OPTARG";;
    R)  SACNM="$OPTARG";;
    M)  MACRO="$OPTARG";;
    [h,?])  usage;;
    esac
done

if [ ! -f "$INDEX" ]
then
	echo "Error: cannot find $INDEX!"
	usage
fi

echo "This run: $0 -D $DATDIR -I $INDEX -S \"$SUBDIR\" -R \"$SACNM\" -M \"$MACRO\""

#umask 002

wkdir=$(pwd)

DATDIR=$(readlink -f $DATDIR) #get the full path of the data directory 
INDEX=$(readlink -f $INDEX) #get the full path of the index file

# temporary list of events to check
tmp=$(mktemp)
awk '$0!~/^[\s\t]*#|checked|skipped/' $INDEX > $tmp # omit lines commented by # or marked as checked/skipped

echo !!!!!!!!!!!!!!!!!!!!
echo $(wc -l $tmp | awk '{printf "%d directories left",$1}')
echo !!!!!!!!!!!!!!!!!!!!

# loop each directory
IFS=$'\n'
for line in $(cat $tmp)
do

	dirnm=$(echo $line | awk '{print $1}')
	
	echo "...processing $line..."

	SACDIR=$DATDIR/$dirnm/$SUBDIR

sac<<END
cd $SACDIR
r ${SACNM}
fileid l ur t l KSTNM KCMPNM KEVNM GCARC MAG
qdp off
window x 0 1. y 0.2 1.
*xlim -100 100
${MACRO}
*ppk p 9 bell off
wh
q
END

	# parse user input: checked or skipped
	userinput=""
	
	echo 
	read -p "Please select([y]:done/[s]:skip/[r]:recheck/[q]:ignore&quit/Other:ignore&continue): " userinput
	
	datestr=$(date +%Y%m%d_%H%M%S)
	
	if [ "$userinput" == 'y' ]; then
		sed -i "/^$line/s/$/ | checked $datestr/" $INDEX
		#	record selection
# 	ls *.BHR | dumpSHD t5 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > bad.lst
# 	ls *.BHR | dumpSHD t6 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > needHPfilter.lst
# 	ls *.BHR | dumpSHD t7 | sed 's/\.BHR//' | awk '$3!=-12345 && $3<0 {print "'$evnm'",$1}' > reversedP.lst
	fi
	
	if [ "$userinput" == 's' ]; then
		sed -i "/^$line/s/$/ | skipped $datestr/" $INDEX
	fi

	if [ "$userinput" == 'r' ]; then
		sed -i "/^$line/s/$/ | recheck $datestr/" $INDEX
	fi
	
	if [ "$userinput" == 'q' ]; then
		break
	fi

done

# delete tempory file
if [ -f $tmp ]; then
	\rm $tmp
fi

# END of script
