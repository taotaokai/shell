#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Read stationXML file and print out text info at the channel level equivalent
to the FDSN-ws text format.

Network|Station|Location|Channel|Latitude|Longitude|Elevation|Depth|Azimuth|Dip|Senso
rDescription|Scale|ScaleFrequency|ScaleUnits|SampleRate|StartTime|EndTime
"""
import sys
from obspy import read_inventory

#====== parameters
#xml_filename = '201012080821291.BU.xml'
xml_filename = str(sys.argv[1])

#====== read xml file
inv = read_inventory(xml_filename, format='STATIONXML')

#====== loop each channel
for net in inv:
    net_code = net.code
    for sta in net:
        sta_code = sta.code
        for cha in sta:
            loc_code = cha.location_code
            cha_code = cha.code
            lat = cha.latitude
            lon = cha.longitude
            ele = cha.elevation
            dep = cha.depth
            az = cha.azimuth
            dip = cha.dip
            sensor_description = cha.sensor.description
            scale = cha.response.instrument_sensitivity.value
            scale_freq = cha.response.instrument_sensitivity.frequency
            scale_units = cha.response.instrument_sensitivity.input_units
            sample_rate = cha.sample_rate
            start_time = cha.start_date
            end_time = cha.end_date

            print '{:s}|{:s}|{:s}|{:s}|'\
                  '{:.4f}|{:.4f}|{:.1f}|{:.1f}|{:.1f}|{:.1f}|{:s}|'\
                  '{:E}|{:.1f}|{:s}|{:.1f}|'\
                  '{:s}|{:s}'.format(
                          net_code, sta_code, loc_code, cha_code, 
                          lat, lon, ele, dep, az, dip, sensor_description, 
                          scale, scale_freq, scale_units, sample_rate, 
                          start_time, end_time)

#END