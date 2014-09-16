#!/bin/bash

# mkKML_line create .kml file for displaying in GoogleEarth
#--------------------------------------------------------------------------
# Usage: mkKML < "x/longitude y/latitude z/height"
#--------------------------------------------------------------------------
# [2013-04-30] created 

#title=$1
#marker=$2

#if [ -z "$title" ]
#then
#	title='title'
#fi
#
#if [ -z "$marker" ]
#then
#	marker='http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png'
#fi

echo '<?xml version="1.0" encoding="UTF-8"?>'
echo '<kml xmlns="http://www.opengis.net/kml/2.2">'
echo '<Document>'
echo "  <name>Line</name>"
echo "  <description>Description</description>"
echo '  <Placemark>'
echo '      <name>LineString</name>'
echo '      <description>yellow lines</description>'
echo '      <Style>'
echo '      <LineStyle>'
echo '        <color>7f00ffff</color>'
echo '        <width>4</width>'
echo '      </LineStyle>'
echo '      </Style>'
echo '      <LineString>'
echo '<!--        <altitudeMode>relativeToGround</altitudeMode> -->'
echo '        <coordinates>'

while read lon lat height 
do
    echo "          $lon,$lat,$height"
done 

echo '        </coordinates>'
echo '      </LineString>'
echo '    </Placemark>'
echo '</Document>'
echo '</kml>'
