#!/bin/bash

# MKKML create .kml file for displaying in GoogleEarth from given {name,lat,lon}
#--------------------------------------------------------------------------
# Usage: mkKML( kmlname, name, lat, lon, kmlfile )
#--------------------------------------------------------------------------
# Inputs:
#   kmlname: the folder name shown in the Places column of GoogleEarth;
#   name: the names shown for the markers in GoogleEarth;
#   lat: latitude of the markers;
#   lon: longitude of the markers;
#   kmlfile: output .kml file path;
#--------------------------------------------------------------------------
# July 20, 2011: created
# Sep. 10, 2011: add Comments

kmlname=$1

#--------------------------------------------------------------------------
# creat head
echo '<?xml version="1.0" encoding="UTF-8"?>'
echo '<kml xmlns="http://earth.google.com/kml/2.2">'
#--------------------------------------------------------------------------
# insert icon definition
# fprintf(fid,['<Style id="sn_ylw-pushpin">\n', ...
# 		'<IconStyle>\n', ...
# 			'<scale>1.1</scale>\n', ...
# 			'<Icon>\n', ...
# 				'<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>\n', ...
# 			'</Icon>\n', ...
# 			'<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>\n', ...
# 		'</IconStyle>\n', ...
# 		'<ListStyle>\n', ...
# 		'</ListStyle>\n', ...
# 	'</Style>\n', ...
#     '<StyleMap id="msn_ylw-pushpin">\n', ...
# 		'<Pair>\n', ...
# 			'<key>normal</key>\n', ...
# 			'<styleUrl>#sn_ylw-pushpin</styleUrl>\n', ...
# 		'</Pair>\n', ...
# 	'</StyleMap>\n']);
#--------------------------------------------------------------------------
# write placemark for each station
echo '<Document>'
echo "<name>$kmlname</name>"

while read line
do

	stnm=$( echo $line | awk '{print $1}' )
	stlo=$( echo $line | awk '{print $2}' )
	stla=$( echo $line | awk '{print $3}' )

    echo '<Placemark>'
    echo "<name>$stnm</name>"
# 	fprintf(fid,'<styleUrl>#sn_ylw-pushpin</styleUrl>\n');
    echo "<Point><coordinates>$stlo,$stla,0</coordinates></Point>"
    echo '</Placemark>'

done

#--------------------------------------------------------------------------
# write end of kml file
echo '</Document>'
echo '</kml>'
