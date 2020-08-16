#!/usr/bin/python

import argparse
from math import radians,sin,cos
import sys
import pylab as P
import numpy as np

def plot_line(ax,data,fx,fy,color='black',label=''):
  xx=[float(d[fx]) for d in data]
  yy=[float(d[fy]) for d in data]
  ax.plot(xx,yy,color=color,label=label)


def plot_times(ax, times):
  for t in times:
    ax.axvspan(t-0.1,t+0.1,facecolor='g',alpha=0.5)

def plot_checks(ax,data,fx,fy,checkGapFirst,checkGapMult,color='red'):
  if not checkGapFirst or not checkGapMult:
    return
  T = data[-1]['time']
  cx,cy=[],[]
  tX=[]
  gap = checkGapFirst
  t = checkGapFirst
  while(t < 0.5*T):
    tX.append(t)
    t = t + gap
    gap = gap * checkGapMult
  gap = checkGapFirst
  t = checkGapFirst
  while(t < 0.5*T):
    tX.append(T-t)
    t = t + gap
    gap = gap * checkGapMult

  for t in tX:
    done = False
    for i,d in enumerate(data):
      if data[i]['time'] >= t and (not done):
        cx.append(data[i][fx])
        cy.append(data[i][fy])
        done = True
  ax.plot(cx,cy,color=color,marker='o',markersize=10,linestyle='')


def plot(datas,labels,xmin,xmax,ymin,ymax,zmin,zmax,
         vxmin,vxmax,vymin,vymax,vzmin,vzmax,tmax=30,amult=1,askip=3,checkGapFirst=None,checkGapMult=None,minDescentAngle=None,
         amin=None,amax=None,
         filenames=[], thrust_times=[]):

  colors=['red','blue','green','black','pink','grey','purple','salmon']
  alim = [0,40]

  # altitude against time, and throttle
  fig = P.figure(1)

  P.subplot2grid((3,5),(0,0), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("x")
  ax.set_xlim([0,tmax])
  ax.set_ylim([xmin,xmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','x',color=col)
    plot_checks(ax,data,'time','x',checkGapFirst,checkGapMult)
  ax.grid()

  P.subplot2grid((3,5),(0,1), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("vx")
  ax.set_xlim([0,tmax])
  ax.set_ylim([vxmin,vxmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','vx',color=col)
    plot_checks(ax,data,'time','vx',checkGapFirst,checkGapMult)
  ax.grid()

  P.subplot2grid((3,5),(1,0), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("y")
  ax.set_xlim([0,tmax])
  ax.set_ylim([ymin,ymax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','y',color=col)
    plot_checks(ax,data,'time','y',checkGapFirst,checkGapMult)
  ax.grid()

  P.subplot2grid((3,5),(1,1), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("vy")
  ax.set_xlim([0,tmax])
  ax.set_ylim([vymin,vymax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','vy',color=col)
    plot_checks(ax,data,'time','vy',checkGapFirst,checkGapMult)
  ax.grid()

  P.subplot2grid((3,5),(2,0), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("z")
  ax.set_xlim([0,tmax])
  ax.set_ylim([zmin,zmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','z',color=col)
    plot_checks(ax,data,'time','z',checkGapFirst,checkGapMult)
  ax.grid()
  ax.legend()
  P.subplot2grid((3,5),(2,1), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("vz")
  ax.set_xlim([0,tmax])
  ax.set_ylim([vzmin,vzmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','vz',color=col)
    plot_checks(ax,data,'time','vz',checkGapFirst,checkGapMult)
  ax.grid()

  # Throttle
  P.subplot2grid((3,5),(0,2), colspan=2, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("mag(accel)")
  ax.set_xlim([0,tmax])
  ax.set_ylim(alim)

  # plot desired magnitude of acceleration
  tdata = []
  for col,data in zip(colors,datas):
    for d in data:
      T=np.array([d['ax'],d['ay'],d['az']])
      d['mag_accel'] = np.linalg.norm(T)
    plot_line(ax,data,'time','mag_accel',color=col)
    plot_checks(ax,data,'time','mag_accel',checkGapFirst,checkGapMult,color=col)
    if amin:
      ax.plot([0,data[-1]['time']],[amin,amin],color='blue',linestyle='--')
    if amax:
      ax.plot([0,data[-1]['time']],[amax,amax],color='blue',linestyle='--')
  plot_times(ax, thrust_times)
  ax.grid()

  # Attitude error
  P.subplot2grid((3,5),(0,4), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("attitude error /degrees")
  ax.set_xlim([0,tmax])
  ax.set_ylim([0,90])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','att_err',color=col)
  ax.grid()

  # XY
  P.subplot2grid((3,5),(1,2), colspan=3, rowspan=2)
  ax = P.gca()
  ax.set_xlabel("x")
  ax.set_ylabel("y")
  ax.set_xlim([xmin,xmax])
  ax.set_ylim([ymin,ymax])

  # plot side view of X,Y
  for di,data in enumerate(datas):
    xx,yy=[],[]
    for i,d in enumerate(data):
      xx = []
      yy = []
      xx.append(d['x'])
      yy.append(d['y'])
      xx.append(d['x']+d['ax']*amult)
      yy.append(d['y']+d['ay']*amult)
      ax.plot(xx,yy,color=colors[di],alpha=0.5)
    plot_line(ax,data,'x','y',color=colors[di],label=filenames[di])
    # Show checkpoints
    plot_checks(ax,data,'x','y',color=colors[di],checkGapFirst=checkGapFirst,checkGapMult=checkGapMult)

  # Draw min descent angle
  if minDescentAngle:
    fx = datas[0][-1]['x']
    fy = datas[0][-1]['y']
    fx = 0
    dy = 0
    s = sin(radians(minDescentAngle))
    c = cos(radians(minDescentAngle))
    d = (xmax-xmin) + (ymax-ymin)
    xx = [-d*c + fx,fx,d*c + fx]
    yy = [d*s + fy,fy,d*s + fy]
    ax.plot(xx,yy,color='blue',linestyle='--')

  ax.legend()
  ax.grid()
  P.show()

def extract_items(line):
  d = {}
  for kv in line.split(" "):
    if '=' in kv:
      k,v = kv.split('=',1)
      d[k] = v
  return d

def read_data(fname, d):
  """Reads column data file, values space or tab separated. First line in column names.
     Comments lines with hash can contain key=value pairs which will be returned in d"""
  fields=None
  dat=[]
  for line in file(fname):
    line=line.strip("\n\r")
    if line.startswith("#"):
      dd = extract_items(line[1:])
      print(dd)
      d.update(dd)
      continue
    if not fields:
      fields = line.split(None)
    else:
      try:
        data = [float(x) for x in line.split(" ")]
        if len(data)==len(fields):
          dat.append( dict(zip(fields,data)) )
      except:
        pass
  return dat

parser = argparse.ArgumentParser(description='Plot vessel data logs (or solutions) with X,Y,Z,VX,VY,VZ and Throttle in multiple plots')
parser.add_argument('filename', nargs='+',
                    help='Filename of TAB-separated data file, first line contains column names. Should contain time,x,y,z,vx,vy,vz,ax,ay,ax')
parser.add_argument('--xmin', type=float, help='Minimum x position', default=None)
parser.add_argument('--xmax', type=float, help='Maximum x position', default=None)
parser.add_argument('--ymin', type=float, help='Minimum y position', default=None)
parser.add_argument('--ymax', type=float, help='Maximum y position', default=None)
parser.add_argument('--zmin', type=float, help='Minimum z position', default=None)
parser.add_argument('--zmax', type=float, help='Maximum z position', default=None)
parser.add_argument('--vxmin', type=float, help='Minimum vx position', default=None)
parser.add_argument('--vxmax', type=float, help='Maximum vx position', default=None)
parser.add_argument('--vymin', type=float, help='Minimum vy position', default=None)
parser.add_argument('--vymax', type=float, help='Maximum vy position', default=None)
parser.add_argument('--vzmin', type=float, help='Minimum vz position', default=None)
parser.add_argument('--vzmax', type=float, help='Maximum vz position', default=None)
parser.add_argument('--tmax', type=float, help='Maximum time', default=None)
parser.add_argument('--amin', type=float, help='Show minimum acceleration constraint', default=None)
parser.add_argument('--amax', type=float, help='Show maximum acceleration constraint', default=None)
parser.add_argument('--checkGapFirst', type=float, help='Time of first check from start and end', default=None)
parser.add_argument('--checkGapMult', type=float, help='Multiplier to set time to next check (use >1)', default=None)
parser.add_argument('--minDescentAngle', type=float, help='Show cone for minimum descent angle', default=None)
parser.add_argument('--amult', type=float, help='Multiplier for scale up thrust acceleration lines', default=1)
parser.add_argument('--square', action='store_true', help='Make XY plot square (roughly as depends on window size)', default=False)

args = parser.parse_args()

datas=[]
info={}
for filename in args.filename:
  datas.append(read_data(filename,info))
thrust_times = []
try:
  args.minDescentAngle = float(info['minDescentAngle'])
except:
  pass
try:
  thrust_times = [float(t) for t in info['thrust_times'].split(',')]
except:
  pass

alldata = []
for data in datas:
  alldata = alldata + data

# Find min and max
if not args.xmin:
  args.xmin = min([d['x'] for d in alldata])
if not args.xmax:
  args.xmax = max([d['x'] for d in alldata])
if not args.ymin:
  args.ymin = min([d['y'] for d in alldata])
if not args.ymax:
  args.ymax = max([d['y'] for d in alldata])
if not args.zmin:
  args.zmin = min([d['z'] for d in alldata])
if not args.zmax:
  args.zmax = max([d['z'] for d in alldata])

if not args.vxmin:
  args.vxmin = min([d['vx'] for d in alldata])
if not args.vxmax:
  args.vxmax = max([d['vx'] for d in alldata])
if not args.vymin:
  args.vymin = min([d['vy'] for d in alldata])
if not args.vymax:
  args.vymax = max([d['vy'] for d in alldata])
if not args.vzmin:
  args.vzmin = min([d['vz'] for d in alldata])
if not args.vzmax:
  args.vzmax = max([d['vz'] for d in alldata])

if not args.tmax:
  args.tmax = max([d['time'] for d in alldata])

if args.square:
  width = args.xmax - args.xmin
  height = args.ymax - args.ymin
  ratio = width/height
  if ratio < 1:
    args.xmin = args.xmin - (height-width)*0.5
    args.xmax = args.xmax + (height-width)*0.5
  if ratio > 1:
    args.ymin = args.ymin - (width-height)*0.5
    args.ymax = args.ymax + (width-height)*0.5

plot(datas,args.filename,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,zmin=args.zmin,zmax=args.zmax,
     vxmin=args.vxmin,vxmax=args.vxmax,vymin=args.vymin,vymax=args.vymax,vzmin=args.vzmin,vzmax=args.vzmax,tmax=args.tmax,
     amult=args.amult,checkGapFirst=args.checkGapFirst,checkGapMult=args.checkGapMult,filenames=args.filename,minDescentAngle=args.minDescentAngle,amin=args.amin,amax=args.amax,thrust_times=thrust_times)
