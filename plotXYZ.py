#!/usr/bin/python

import argparse
import sys
import pylab as P
import numpy as np

def plot_line(ax,data,fx,fy,color='black',label=''):
  xx=[float(d[fx]) for d in data]
  yy=[float(d[fy]) for d in data]
  ax.plot(xx,yy,color=color,label=label)

def plot(datas,labels,xmin,xmax,ymin,ymax,zmin,zmax,
         vxmin,vxmax,vymin,vymax,vzmin,vzmax,tmax=30,amult=1,askip=3,numchecks=None,filenames=[]):

  colors=['red','blue','green','black','pink','grey','purple','salmon']
  amax = 40

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
  ax.grid()

  P.subplot2grid((3,5),(0,1), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("vx")
  ax.set_xlim([0,tmax])
  ax.set_ylim([vxmin,vxmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','vx',color=col)
  ax.grid()

  P.subplot2grid((3,5),(1,0), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("y")
  ax.set_xlim([0,tmax])
  ax.set_ylim([ymin,ymax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','y',color=col)
  ax.grid()

  P.subplot2grid((3,5),(1,1), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("vy")
  ax.set_xlim([0,tmax])
  ax.set_ylim([vymin,vymax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','vy',color=col)
  ax.grid()

  P.subplot2grid((3,5),(2,0), colspan=1, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("z")
  ax.set_xlim([0,tmax])
  ax.set_ylim([zmin,zmax])
  for col,data in zip(colors,datas):
    plot_line(ax,data,'time','z',color=col)
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
  ax.grid()

  # Throttle
  P.subplot2grid((3,5),(0,2), colspan=3, rowspan=1)
  ax = P.gca()
  ax.set_xlabel("time")
  ax.set_ylabel("mag(accel)")
  ax.set_xlim([0,tmax])
  ax.set_ylim([0,amax])

  # plot desired magnitude of acceleration
  for col,data in zip(colors,datas):
    tt=[]
    throttle=[]
    for d in data:
      T=np.array([d['ax'],d['ay'],d['az']])
      tt.append(d['time'])
      throttle.append( np.linalg.norm(T) )
    ax.plot(tt,throttle,color=col)
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
    throttle=[]
    for i,d in enumerate(data):
      xx = []
      yy = []
      xx.append(d['x'])
      yy.append(d['y'])
      xx.append(d['x']+d['ax']*amult)
      yy.append(d['y']+d['ay']*amult)
      ax.plot(xx,yy,color=col,alpha=0.5)
    plot_line(ax,data,'x','y',color=colors[di],label=filenames[di])
    # Show checkpoints
    if numchecks:
      cx,cy=[],[]
      gap = len(data)/(numchecks+1)
      for i in range(gap,len(data)-gap,gap):
        cx.append(data[i]['x'])
        cy.append(data[i]['y'])
      ax.plot(cx,cy,color=colors[di],marker='o',markersize=10,linestyle='')

  ax.legend()
  ax.grid()
  P.show()


def read_data(fname):
  fields=None
  dat=[]
  for line in file(fname):
    line=line.strip("\n\r")
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
parser.add_argument('--numchecks', type=int, help='How many checkpoints to show (spaced evenly accept at t={0,T}', default=None)
parser.add_argument('--amult', type=float, help='Multiplier for scale up thrust acceleration lines', default=1)

args = parser.parse_args()

datas=[]
for filename in args.filename:
  datas.append(read_data(filename))

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

plot(datas,args.filename,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,zmin=args.zmin,zmax=args.zmax,
     vxmin=args.vxmin,vxmax=args.vxmax,vymin=args.vymin,vymax=args.vymax,vzmin=args.vzmin,vzmax=args.vzmax,tmax=args.tmax,
     amult=args.amult,numchecks=args.numchecks,filenames=args.filename)
