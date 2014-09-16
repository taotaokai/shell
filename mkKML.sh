#!/bin/bash

# MKKML create .kml file for displaying in GoogleEarth
#--------------------------------------------------------------------------
# Usage: cat LIST | ./mkKML TITLE MARKER > FILE.kml
# LIST is a file containing "stlo stla stnm" on each row (like, 135 35 ST1).
# TITLE is the name showing in the GoogleEarth (if blank, default name used);
# MARKER is the symbol file name (if blank, default symbol used);
#--------------------------------------------------------------------------
# July 20, 2011: created
# Sep. 10, 2011: add Comments
# [2012-08-14] modified

title=$1
marker=$2

if [ -z "$title" ]
then
	title='title'
fi

if [ -z "$marker" ]
then
	marker='http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png'
fi

echo '<?xml version="1.0" encoding="UTF-8"?>'
echo '<kml xmlns="http://earth.google.com/kml/2.2">'
echo '<Document>'
echo "<name>$title</name>"

while read stlo stla stnm 
do
	echo '<Placemark>'
	echo "<Style><IconStyle><Icon><href>$marker</href></Icon></IconStyle></Style>"
	echo "<name>$stnm</name>"
	echo "<Point><coordinates>$stlo,$stla,0</coordinates></Point>"
	echo '</Placemark>'
done 

echo '</Document>'
echo '</kml>'
