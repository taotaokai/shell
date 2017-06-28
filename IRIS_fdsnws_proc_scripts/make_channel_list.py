#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""extract operating channels from a channel list  
    for a given event, and check channel infos

"""
import sys
import re
import numpy as np
from obspy import UTCDateTime

#====== parameters
cmt_file = str(sys.argv[1])
channel_file = str(sys.argv[2])
out_file = str(sys.argv[3])

#====== read CMTSOLUTION file
with open(cmt_file, 'r') as f:
    x = f.readline()

x = x.split()
event_time = UTCDateTime('%s-%s-%sT%s:%s:%s' % \
        (x[1], x[2], x[3], x[4], x[5], x[6]))

#====== parse channel list 
with open(channel_file, 'r') as f:
    lines = [x for x in f.readlines() if not(x.startswith('#'))]
lines = [x.replace('\n','').split('|') for x in lines]

channels = {}
for x in lines:
    date1 = [ int(a) for a in re.sub("\D", " ", x[15]).split() ]
    date2 = [ int(a) for a in re.sub("\D", " ", x[16]).split() ]
    t1 = UTCDateTime(date1[0], date1[1], date1[2]) \
            + 60.0*(60.0*date1[3] + date1[4]) + date1[5]
    t2 = UTCDateTime(date2[0], date2[1], date2[2]) \
            + 60.0*(60.0*date2[3] + date2[4]) + date2[5]

    # only get channels operating on event_time
    if t1 < event_time and event_time < t2:
        key = (x[0], x[1], x[2])
        channel = {'code':        x[3],
                   'latitude':    float(x[4]),
                   'longitude':   float(x[5]),
                   'elevation':   float(x[6]),
                   'depth':       float(x[7]),
                   'azimuth':     float(x[8]),
                   'dip':         float(x[9]),
                   'all':         x }
        if key not in channels:
            channels[key] = []
        channels[key].append(channel)

#====== check channel info
fp_out = open(out_file, 'w')

for key in channels:

    cha = channels[key]

    lats = [ x['latitude'] for x in cha ]
    lons = [ x['longitude'] for x in cha ]
    eles = [ x['elevation'] for x in cha ]
    deps = [ x['depth'] for x in cha ]
    Z_comp = [ (x['code'], x['azimuth'], x['dip'])
            for x in cha if x['code'][2] == 'Z']
    H_comp = [ (x['code'], x['azimuth'], x['dip']) \
            for x in cha if x['code'][2] != 'Z']

    flag_bad = False

    # 3 components 
    if len(cha) != 3: 
        print("# [WARNING] Not 3 components in ", key)
        flag_bad = True
    # same locations
    if lats.count(lats[0])!=len(lats) or lons.count(lons[0])!=len(lons) or \
            eles.count(eles[0])!=len(eles) or deps.count(deps[0])!=len(deps): 
        print("# [WARNING] Not the same coordinates in ", key)
        flag_bad = True 
    # Z-comp polarity 
    if len(Z_comp) != 1 or Z_comp[0][2] != -90.0: 
        print("# [WARNING] Problematic Z comp in ", key)
        flag_bad = True
    # H-comps polarity and orthogonality
    if len(H_comp) != 2 or abs(H_comp[0][2]) != 0.0 or \
            abs(H_comp[1][2]) != 0.0 or \
            abs(np.cos(np.deg2rad(H_comp[0][1] - H_comp[1][1]))) > 0.1: 
        print("# [WARNING] Problematic H comps in ", key)
        flag_bad = True

    # print out 3 componets
    for comp in cha:
        line = '|'.join(comp['all'])
        if flag_bad:
            line = "#" + line
        fp_out.write("%s\n" % line)

fp_out.close()

#END
