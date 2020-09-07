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


def plot_times(ax, times, color):
  for t in times:
    ax.axvspan(t-0.1,t+0.1,facecolor=color,alpha=0.5)

def plot_checks(ax,data,fx,fy,times,color='black'):
  cx = []
  cy = []
  for t in times:
    done = False
    for i,d in enumerate(data):
      if data[i]['time'] >= t and (not done):
        cx.append(data[i][fx])
        cy.append(data[i][fy])
        done = True
  ax.plot(cx,cy,color=color,marker='o',markersize=4,linestyle='')

def plot_targets(ax,data,color='black'):
  tx = [d[0] for d in data]
  ty = [d[1] for d in data]
  ax.plot(tx,ty,color=color,marker='s',markersize=10,linestyle='')

class Vector3Time:
  def fromStr(self, s=''):
    v = s.split(":")[0]
    v = v.strip("[]")
    try:
      self.t = float(s.split(":")[1])
    except:
      self.t = -1
    self.x,self.y,self.z = [float(a) for a in v.split(",")]
    return self

  def __init__(self, x=0, y=0, z=0, t=-1):
    self.x = x
    self.y = y
    self.z = z
    self.t = t

def plot(labels,xmin,xmax,ymin,ymax,zmin,zmax,
         vxmin,vxmax,vymin,vymax,vzmin,vzmax,accelmax,tmax=30,amult=1,askip=3,
         filenames=[],showchecks=False):

  # Set up figures
  # altitude against time, and throttle
  fig = P.figure(1)

  colors=['red','blue','green','black','pink','grey','purple','salmon']

  P.subplot2grid((3,5),(0,0), colspan=1, rowspan=1)
  ax1 = P.gca()
  ax1.set_xlabel("time")
  ax1.set_ylabel("x")
  ax1.set_xlim([0,tmax])
  ax1.set_ylim([xmin,xmax])
  ax1.grid()

  P.subplot2grid((3,5),(0,1), colspan=1, rowspan=1)
  ax2 = P.gca()
  ax2.set_xlabel("time")
  ax2.set_ylabel("vx")
  ax2.set_xlim([0,tmax])
  ax2.set_ylim([vxmin,vxmax])
  ax2.grid()

  P.subplot2grid((3,5),(1,0), colspan=1, rowspan=1)
  ax3 = P.gca()
  ax3.set_xlabel("time")
  ax3.set_ylabel("y")
  ax3.set_xlim([0,tmax])
  ax3.set_ylim([ymin,ymax])
  ax3.grid()

  P.subplot2grid((3,5),(1,1), colspan=1, rowspan=1)
  ax4 = P.gca()
  ax4.set_xlabel("time")
  ax4.set_ylabel("vy")
  ax4.set_xlim([0,tmax])
  ax4.set_ylim([vymin,vymax])
  ax4.grid()

  P.subplot2grid((3,5),(2,0), colspan=1, rowspan=1)
  ax5 = P.gca()
  ax5.set_xlabel("time")
  ax5.set_ylabel("z")
  ax5.set_xlim([0,tmax])
  ax5.set_ylim([zmin,zmax])
  ax5.grid()

  P.subplot2grid((3,5),(2,1), colspan=1, rowspan=1)
  ax6 = P.gca()
  ax6.set_xlabel("time")
  ax6.set_ylabel("vz")
  ax6.set_xlim([0,tmax])
  ax6.set_ylim([vzmin,vzmax])
  ax6.grid()

  # Throttle
  P.subplot2grid((3,5),(0,2), colspan=2, rowspan=1)
  ax7 = P.gca()
  ax7.set_xlabel("time")
  ax7.set_ylabel("mag(accel)")
  ax7.set_xlim([0,tmax])
  ax7.set_ylim([0,accelmax])
  ax7.grid()

  # Attitude error
  P.subplot2grid((3,5),(0,4), colspan=1, rowspan=1)
  ax8 = P.gca()
  ax8.set_xlabel("time")
  ax8.set_ylabel("attitude error /degrees")
  ax8.set_xlim([0,tmax])
  ax8.set_ylim([0,90])
  ax8.grid()

  # XY
  P.subplot2grid((3,5),(1,2), colspan=3, rowspan=2)
  ax9 = P.gca()
  ax9.set_xlabel("x")
  ax9.set_ylabel("y")
  ax9.set_xlim([xmin,xmax])
  ax9.set_ylim([ymin,ymax])
  ax9.grid()


  for di,filename in enumerate(filenames):
    col = colors[di]
    data = read_data(filename,info)
    thrust_times = []
    check_times = []
    targets = []
    if 'minDescentAngle' in info:
      minDescentAngle = float(info['minDescentAngle'])
    else:
      minDescentAngle = None
    amin = 0
    if 'amin' in info:
      amin = float(info['amin'])
    amax = 30
    if 'amax' in info:
      amax = float(info['amax'])
    if 'thrust_times' in info:
      thrust_times = [float(t) for t in info['thrust_times'].split(",")]
    if 'check_times' in info:
      check_times = [float(t) for t in info['check_times'].split(",")]
    if not showchecks:
      check_times = []
    if 'target' in info:
      for s in info['target']:
        t=Vector3Time()
        targets.append(t.fromStr(s))
    if 'rf' in info:
      t=Vector3Time()
      targets.append(t.fromStr(info['rf']))

    plot_line(ax1,data,'time','x',color=col)
    plot_checks(ax1,data,'time','x',check_times,color=col)
    plot_targets(ax1,[(t.t,t.x) for t in targets])

    plot_line(ax2,data,'time','vx',color=col)
    plot_checks(ax2,data,'time','vx',check_times,color=col)

    plot_line(ax3,data,'time','y',color=col)
    plot_checks(ax3,data,'time','y',check_times,color=col)
    plot_targets(ax3,[(t.t,t.y) for t in targets])

    plot_line(ax4,data,'time','vy',color=col)
    plot_checks(ax4,data,'time','vy',check_times,color=col)

    plot_line(ax5,data,'time','z',color=col)
    plot_checks(ax5,data,'time','z',check_times,color=col)
    plot_targets(ax5,[(t.t,t.z) for t in targets])

    plot_line(ax6,data,'time','vz',color=col)
    plot_checks(ax6,data,'time','vz',check_times,color=col)

    # plot desired magnitude of acceleration
    tdata = []
    for d in data:
      T=np.array([d['ax'],d['ay'],d['az']])
      d['mag_accel'] = np.linalg.norm(T)
    plot_line(ax7,data,'time','mag_accel',color=col)
    plot_checks(ax7,data,'time','mag_accel',check_times,color=col)
    if amin:
      ax7.plot([0,data[-1]['time']],[amin,amin],color='blue',linestyle='--')
    if amax:
      ax7.plot([0,data[-1]['time']],[amax,amax],color='blue',linestyle='--')
    plot_times(ax7, thrust_times, color=col)

    plot_line(ax8,data,'time','att_err',color=col)


    # plot side view of X,Y
    xx,yy=[],[]
    for i,d in enumerate(data):
      xx = []
      yy = []
      xx.append(d['x'])
      yy.append(d['y'])
      xx.append(d['x']+d['ax']*amult)
      yy.append(d['y']+d['ay']*amult)
      ax9.plot(xx,yy,color=colors[di],alpha=0.5)
    plot_line(ax9,data,'x','y',color=colors[di],label=filenames[di])
    # Show checkpoints
    plot_checks(ax9,data,'x','y',check_times,color=colors[di])
    plot_targets(ax9,[(t.x,t.y) for t in targets])

  # Draw min descent angle
    if minDescentAngle is not None:
      fx = datas[0][-1]['x']
      fy = datas[0][-1]['y']
      fx = 0
      dy = 0
      s = sin(radians(minDescentAngle))
      c = cos(radians(minDescentAngle))
      d = (xmax-xmin) + (ymax-ymin)
      xx = [-d*c + fx,fx,d*c + fx]
      yy = [d*s + fy,fy,d*s + fy]
      ax9.plot(xx,yy,color=colors[di],linestyle='--')

  ax9.legend()
  P.show()

def extract_items(line, lists=[]):
  d = {}
  for kv in line.split(" "):
    if '=' in kv:
      k,v = kv.split('=',1)
      if k not in lists:
        d[k] = v
      if k in lists:
        try:
          d[k].append(v)
        except:
          d[k] = [v]
  return d

def read_data(fname, d):
  """Reads column data file, values space or tab separated. First line in column names.
     Comments lines with hash can contain key=value pairs which will be returned in d"""
  fields=None
  dat=[]
  for line in file(fname):
    line=line.strip("\n\r")
    if line.startswith("#"):
      dd = extract_items(line[1:], lists=['target'])
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
parser.add_argument('--accelmax', type=float, help='Maximum acceleration', default=None)
parser.add_argument('--tmax', type=float, help='Maximum time', default=None)
parser.add_argument('--amult', type=float, help='Multiplier for scale up thrust acceleration lines', default=1)
parser.add_argument('--square', action='store_true', help='Make XY plot square (roughly as depends on window size)', default=False)
parser.add_argument('--showchecks', action='store_true', help='Show time checks for max vel. and min descent angle', default=False)

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
if not args.accelmax:
  args.accelmax = 0
  for d in alldata:
    T=np.array([d['ax'],d['ay'],d['az']])
    args.accelmax = max(np.linalg.norm(T),args.accelmax)

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

plot(args.filename,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,zmin=args.zmin,zmax=args.zmax,
     vxmin=args.vxmin,vxmax=args.vxmax,vymin=args.vymin,vymax=args.vymax,vzmin=args.vzmin,vzmax=args.vzmax,tmax=args.tmax,
     amult=args.amult,filenames=args.filename,showchecks=args.showchecks,accelmax=args.accelmax)
