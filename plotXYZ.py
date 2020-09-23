#!/usr/bin/python

import argparse
from math import radians,sin,cos
import sys
import pylab as P
import numpy as np

def plot_line(ax,data,fx,fy,color='black',label='',linestyle='-'):
  xx=[float(d[fx]) for d in data]
  yy=[float(d[fy]) for d in data]
  ax.plot(xx,yy,color=color,label=label,linestyle=linestyle)


def plot_times(ax, times, color):
  for t in times:
    ax.axvspan(t-0.1,t+0.1,facecolor=color,alpha=0.5)

def plot_markers(ax,data,fx,fy,times,color='black',marker='o',markersize=4,alpha=1):
  cx = []
  cy = []
  for t in times:
    done = False
    for i,d in enumerate(data):
      if data[i]['time'] >= t and (not done):
        cx.append(data[i][fx])
        cy.append(data[i][fy])
        done = True
  ax.plot(cx,cy,color=color,marker=marker,markersize=markersize,linestyle='',alpha=alpha)

def plot_targets(ax,data,color='black'):
  tx = [d[0] for d in data]
  ty = [d[1] for d in data]
  ax.plot(tx,ty,color=color,marker='s',markersize=10,linestyle='')

def plot_scatter(ax,data,fx,fy,color='black'):
  tx = [d[fx] for d in data]
  ty = [d[fy] for d in data]
  ax.plot(tx,ty,color=color,marker='o',markersize=5,linestyle='')


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
         filenames=[],showchecks=False,savepng=None,
         marktime=None,timevsfuel=False):

  # Set up figures
  # altitude against time, and throttle
  fig = P.figure(1,figsize=(15,10))

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
  if timevsfuel:
    P.subplot2grid((3,5),(0,2), colspan=1, rowspan=1)
  else:
    P.subplot2grid((3,5),(0,2), colspan=2, rowspan=1)
  ax7 = P.gca()
  ax7.set_xlabel("time")
  ax7.set_ylabel("mag(accel)")
  ax7.set_xlim([0,tmax])
  ax7.set_ylim([0,accelmax*1.1])
  ax7.grid()

  # T vs fuel
  if timevsfuel:
    P.subplot2grid((3,5),(0,3), colspan=1, rowspan=1)
    ax8 = P.gca()
    ax8.set_xlabel("time")
    ax8.set_ylabel("fuel")
    ax8.grid()
  else:
    ax8 = None

  # Attitude error
  P.subplot2grid((3,5),(0,4), colspan=1, rowspan=1)
  ax9 = P.gca()
  ax9.set_xlabel("time")
  ax9.set_ylabel("attitude error /degrees")
  ax9.set_xlim([0,tmax])
  ax9.set_ylim([0,90])
  ax9.grid()

  # XY
  P.subplot2grid((3,5),(1,2), colspan=3, rowspan=2)
  ax10 = P.gca()
  ax10.set_xlabel("x")
  ax10.set_ylabel("y")
  ax10.set_xlim([xmin,xmax])
  ax10.set_ylim([ymin,ymax])
  ax10.grid()


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
    amax = 0
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
    rf = None
    if 'rf' in info:
      t=Vector3Time()
      rf = t.fromStr(info['rf'])
      targets.append(rf)
    solutions = []
    sln_Tmax = 0
    sln_fuelmin = 1000
    sln_fuelmax = 0
    if ('solution_time' in info) and ('solution_fuel' in info):
      for T,fuel in zip(info['solution_time'].split(","), info['solution_fuel'].split(",")):
        T = float(T)
        fuel = float(fuel)
        solutions.append({'T':T, 'fuel':fuel})
        sln_fuelmin = min(fuel,sln_fuelmin)
        sln_fuelmax = max(fuel,sln_fuelmax)
        sln_Tmax = max(T,sln_Tmax)

    plot_line(ax1,data,'time','x',color=col)
    plot_markers(ax1,data,'time','x',check_times,color=col)
    plot_targets(ax1,[(t.t,t.x) for t in targets])
    plot_markers(ax1,data,'time','x',[marktime],color=colors[di],markersize=10,alpha=0.5)

    plot_line(ax2,data,'time','vx',color=col)
    plot_markers(ax2,data,'time','vx',check_times,color=col)
    plot_markers(ax2,data,'time','vx',[marktime],color=colors[di],markersize=10,alpha=0.5)

    plot_line(ax3,data,'time','y',color=col)
    plot_markers(ax3,data,'time','y',check_times,color=col)
    plot_targets(ax3,[(t.t,t.y) for t in targets])
    plot_markers(ax3,data,'time','y',[marktime],color=colors[di],markersize=10,alpha=0.5)

    plot_line(ax4,data,'time','vy',color=col)
    plot_markers(ax4,data,'time','vy',check_times,color=col)
    plot_markers(ax4,data,'time','vy',[marktime],color=colors[di],markersize=10,alpha=0.5)

    plot_line(ax5,data,'time','z',color=col)
    plot_markers(ax5,data,'time','z',check_times,color=col)
    plot_targets(ax5,[(t.t,t.z) for t in targets])
    plot_markers(ax5,data,'time','z',[marktime],color=colors[di],markersize=10,alpha=0.5)

    plot_line(ax6,data,'time','vz',color=col)
    plot_markers(ax6,data,'time','vz',check_times,color=col)
    plot_markers(ax6,data,'time','vz',[marktime],color=colors[di],markersize=10,alpha=0.5)

    # plot desired magnitude of acceleration
    tdata = []
    for d in data:
      T=np.array([d['ax'],d['ay'],d['az']])
      d['mag_accel'] = np.linalg.norm(T)
    plot_line(ax7,data,'time','mag_accel',color=col)
    plot_markers(ax7,data,'time','mag_accel',check_times,color=col)
    plot_markers(ax7,data,'time','mag_accel',[marktime],color=colors[di],markersize=10,alpha=0.5)
    if 'amin' in data[0]: # continuos amin values
      plot_line(ax7,data,'time','amin',color=col,linestyle='--')
    elif amin:
      ax7.plot([0,data[-1]['time']],[amin,amin],color=col,linestyle='--')
    if 'amax' in data[0]: # continuos amax values
      plot_line(ax7,data,'time','amax',color=col,linestyle='--')
    elif amax:
      ax7.plot([0,data[-1]['time']],[amax,amax],color=col,linestyle='--')
    plot_times(ax7, thrust_times, color=col)

    if ax8:
      ax8.set_xlim([0,sln_Tmax])
      ax8.set_ylim([sln_fuelmin,sln_fuelmax])
      plot_scatter(ax8,solutions,'T','fuel',color=col)

    plot_line(ax9,data,'time','att_err',color=col)
    plot_markers(ax9,data,'time','att_err',[marktime],color=colors[di],markersize=10,alpha=0.5)


    # plot side view of X,Y
    xx,yy=[],[]
    for i,d in enumerate(data):
      xx = []
      yy = []
      xx.append(d['x'])
      yy.append(d['y'])
      xx.append(d['x']+d['ax']*amult)
      yy.append(d['y']+d['ay']*amult)
      ax10.plot(xx,yy,color=colors[di],alpha=0.5)
    plot_line(ax10,data,'x','y',color=colors[di],label=filenames[di])
    if marktime:
      plot_markers(ax10,data,'x','y',[marktime],color=colors[di],markersize=10,alpha=0.5)
    # Show checkpoints
    plot_markers(ax10,data,'x','y',check_times,color=colors[di])
    plot_targets(ax10,[(t.x,t.y) for t in targets])

  # Draw min descent angle
    if minDescentAngle is not None:
      if rf:
        fx = rf.x
        fy = rf.y
      else:
        fx = datas[0][-1]['x']
        fy = datas[0][-1]['y']
      fx = 0
      dy = 0
      s = sin(radians(minDescentAngle))
      c = cos(radians(minDescentAngle))
      d = (xmax-xmin) + (ymax-ymin)
      xx = [-d*c + fx,fx,d*c + fx]
      yy = [d*s + fy,fy,d*s + fy]
      ax10.plot(xx,yy,color=colors[di],linestyle='--')

  ax10.legend()

  fig.tight_layout(pad=0.5)
  if savepng:
    P.savefig(savepng)
  else:
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
parser.add_argument('--marktime', type=float, help='Put a marker a this time position', default=None)
parser.add_argument('--accelmax', type=float, help='Maximum acceleration', default=None)
parser.add_argument('--tmax', type=float, help='Maximum time', default=None)
parser.add_argument('--amult', type=float, help='Multiplier for scale up thrust acceleration lines', default=1)
parser.add_argument('--square', action='store_true', help='Make XY plot square (roughly as depends on window size)', default=False)
parser.add_argument('--showchecks', action='store_true', help='Show time checks for max vel. and min descent angle', default=False)
parser.add_argument('--timevsfuel', action='store_true', help='Include plot with time vs fuel', default=False)
parser.add_argument('--savepng', help='PNG filename to save plot to', default=None)

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
    args.zmin = args.zmin - (height-width)*0.5
    args.zmax = args.zmax + (height-width)*0.5
  if ratio > 1:
    args.ymin = args.ymin - (width-height)*0.5
    args.ymax = args.ymax + (width-height)*0.5

plot(args.filename,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,zmin=args.zmin,zmax=args.zmax,
     vxmin=args.vxmin,vxmax=args.vxmax,vymin=args.vymin,vymax=args.vymax,vzmin=args.vzmin,vzmax=args.vzmax,tmax=args.tmax,
     amult=args.amult,filenames=args.filename,showchecks=args.showchecks,accelmax=args.accelmax,savepng=args.savepng,
     marktime=args.marktime, timevsfuel=args.timevsfuel)
