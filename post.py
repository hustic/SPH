#!/usr/bin/env python

import glob
from lxml import etree
import itertools
import pandas as pd
import numpy as np
import sys


files = sorted(glob.glob(sys.argv[1]+'*'+'.vtp'))
print(files)
iters = []
for f in files:
    doc = etree.parse(f)

    pointdata = doc.find('.//PointData')
    points = doc.find('.//Points')

    position = etree.tostring(points, encoding='unicode').split("\n")[2].strip().split(" ")


    pos = []
    for t in zip(*[iter(position)]*3):
        pos.append([float(t[0]), float(t[1]), float(t[2])])
    
    datas = etree.tostring(pointdata, encoding='unicode').strip().split("\n")

    pressure = datas[2].strip().split(" ")
    velocity = datas[5].strip().split(" ")
    types = datas[8].strip().split(" ")
    density = datas[11].strip().split(" ")

    vel = []
    for t in zip(*[iter(velocity)]*3):
        vel.append([float(t[0]), float(t[1]), float(t[2])])
    vel = np.array(vel)
    press = []
    for t in pressure:
        press.append(float(t))
    ty = []
    for t in types:
        ty.append(int(t)) 
    den = []
    for t in density:
        den.append(float(t))

#print(position)
    df = pd.DataFrame(pos, columns=['x', 'y', 'z'])
    df['vel_x'], df['vel_y'], df['vel_z'] = [vel[:, 0], vel[:, 1], vel[:, 2]]
    df['pressure'] = press
    df['pressure'] = df['pressure'].astype(float)
    df['types'] = types
    df['types'] = df['types'].astype(int)
    df['density'] = density
    if datas[14][0]!='<':
        spacingDensity = datas[14].strip().split(" ")
        sd = []
        for t in spacingDensity:
            sd.append(float(t))
        df['spacingDensity'] = sd
        iters.append(df)
  
hdf = pd.HDFStore(sys.argv[1]+'.h5', mode='w')
for i in range(len(iters)):
    #iters[i].to_hdf(sys.argv[1]+'.h5',key='s'+str(i), mode='w')
    hdf.put('s'+"{:05d}".format(i), iters[i])

hdf.close()
    
    
