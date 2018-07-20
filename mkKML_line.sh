#!/bin/bash

selfdoc() {
cat<<EOF
NAME

SYNOPSIS
  mkKML_line.sh input_file group_name

DESCRIPTION

PARAMETERS"
  input_file: multi-segment lines of lon,lat seperated by >
  group_name: line group name

NOTES
EOF
}

#====== check input arguments
if [ "$#" -lt 1 ]
then
  selfdoc
  exit -1
fi

input_file=${1:?[arg]need intpu_file}
group_name=${2:-NAN}

#====== output KML file
cat<<EOF
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>${group_name}</name>
  <description></description>
  <Style id="yellowLine">
    <LineStyle>
      <color>7f00ffff</color>
      <width>4</width>
    </LineStyle>
  </Style>

  <Placemark> 
    <name>LineString</name> 
    <description>fault line</description> 
    <styleUrl>#yellowline</styleUrl> 
    <LineString> 
      <altitudeMode>relativeToGround</altitudeMode>
      <coordinates>
EOF

awk '$1==">"{printf "</coordinates> </LineString> </Placemark>\n<Placemark> <name>LineString</name> <description>fault line</description> <styleUrl>#yellowline</styleUrl> <LineString> <altitudeMode>relativeToGround</altitudeMode> <coordinates>\n"};
$1!=">"&&$1!~/#/{printf "%f,%f,0 \n", $1, $2}' $input_file

cat<<EOF
</coordinates> </LineString> </Placemark>
</Document>
</kml>
EOF
