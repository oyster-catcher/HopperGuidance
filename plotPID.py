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

def plot(data,xmin,xmax,ymin,ymax,zmin,zmax,
         vxmin,vxmax,vymin,vymax,vzmin,vzmax,
         axmin,axmax,aymin,aymax,azmin,azmax,
         tmin=0,tmax=30,
         savepng=None,timevsfuel=False):

  # Set up figures
  # altitude against time, and throttle
  fig = P.figure(1,figsize=(15,10))

  colors=['red','blue','green','black','pink','grey','purple','salmon']

  P.subplot2grid((3,3),(0,0), colspan=1, rowspan=1)
  ax1 = P.gca()
  ax1.set_xlabel("time")
  ax1.set_ylabel("x_err")
  ax1.set_xlim([tmin,tmax])
  ax1.set_ylim([xmin,xmax])
  ax1.grid()

  P.subplot2grid((3,3),(0,1), colspan=1, rowspan=1)
  ax2 = P.gca()
  ax2.set_xlabel("time")
  ax2.set_ylabel("y_err")
  ax2.set_xlim([tmin,tmax])
  ax2.set_ylim([ymin,ymax])
  ax2.grid()

  P.subplot2grid((3,3),(0,2), colspan=1, rowspan=1)
  ax3 = P.gca()
  ax3.set_xlabel("time")
  ax3.set_ylabel("z_err")
  ax3.set_xlim([tmin,tmax])
  ax3.set_ylim([zmin,zmax])
  ax3.grid()

  P.subplot2grid((3,3),(1,0), colspan=1, rowspan=1)
  ax4 = P.gca()
  ax4.set_xlabel("time")
  ax4.set_ylabel("vx_err")
  ax4.set_xlim([tmin,tmax])
  ax4.set_ylim([vxmin,vxmax])
  ax4.grid()

  P.subplot2grid((3,3),(1,1), colspan=1, rowspan=1)
  ax5 = P.gca()
  ax5.set_xlabel("time")
  ax5.set_ylabel("vy_err")
  ax5.set_xlim([tmin,tmax])
  ax5.set_ylim([vymin,vymax])
  ax5.grid()

  P.subplot2grid((3,3),(1,2), colspan=1, rowspan=1)
  ax6 = P.gca()
  ax6.set_xlabel("time")
  ax6.set_ylabel("vz_err")
  ax6.set_xlim([tmin,tmax])
  ax6.set_ylim([vzmin,vzmax])
  ax6.grid()

  P.subplot2grid((3,3),(2,0), colspan=1, rowspan=1)
  ax7 = P.gca()
  ax7.set_xlabel("time")
  ax7.set_ylabel("ax")
  ax7.set_xlim([tmin,tmax])
  ax7.set_ylim([axmin,axmax])
  ax7.grid()

  P.subplot2grid((3,3),(2,1), colspan=1, rowspan=1)
  ax8 = P.gca()
  ax8.set_xlabel("time")
  ax8.set_ylabel("ay")
  ax8.set_xlim([tmin,tmax])
  ax8.set_ylim([aymin,aymax])
  ax8.grid()

  P.subplot2grid((3,3),(2,2), colspan=1, rowspan=1)
  ax9 = P.gca()
  ax9.set_xlabel("time")
  ax9.set_ylabel("az")
  ax9.set_xlim([tmin,tmax])
  ax9.set_ylim([azmin,azmax])
  ax9.grid()

  col='red'
  plot_line(ax1,data,'time','x_err',color=col)

  col='red'
  plot_line(ax1,data,'time','x_err',color=col)
  plot_line(ax2,data,'time','y_err',color=col)
  plot_line(ax3,data,'time','z_err',color=col)

  plot_line(ax4,data,'time','vx_err',color=col)
  plot_line(ax5,data,'time','vy_err',color=col)
  plot_line(ax6,data,'time','vz_err',color=col)

  plot_line(ax7,data,'time','ax',color=col)
  plot_line(ax8,data,'time','ay',color=col)
  plot_line(ax9,data,'time','az',color=col)
  ax1.legend()

  fig.tight_layout(pad=0.5)
  if savepng:
    P.savefig(savepng)
  else:
    P.show()


def read_data(fname, d):
  """Reads column data file, values space or tab separated. First line in column names.
     Comments lines with hash can contain key=value pairs which will be returned in d"""
  fields=None
  dat=[]
  for line in open(fname):
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
parser.add_argument('filename',
                    help='Filename of TAB-separated data file, first line contains column names. Should contain time,x,y,z,vx,vy,vz,ax,ay,ax')
parser.add_argument('--xmin', type=float, help='Minimum x position', default=None)
parser.add_argument('--xmax', type=float, help='Maximum x position', default=None)
parser.add_argument('--ymin', type=float, help='Minimum y position', default=None)
parser.add_argument('--ymax', type=float, help='Maximum y position', default=None)
parser.add_argument('--zmin', type=float, help='Minimum z position', default=None)
parser.add_argument('--zmax', type=float, help='Maximum z position', default=None)
parser.add_argument('--tmin', type=float, help='Minimum time', default=None)
parser.add_argument('--tmax', type=float, help='Maximum time', default=None)
parser.add_argument('--vxmin', type=float, help='Minimum vx position', default=None)
parser.add_argument('--vxmax', type=float, help='Maximum vx position', default=None)
parser.add_argument('--vymin', type=float, help='Minimum vy position', default=None)
parser.add_argument('--vymax', type=float, help='Maximum vy position', default=None)
parser.add_argument('--vzmin', type=float, help='Minimum vz position', default=None)
parser.add_argument('--vzmax', type=float, help='Maximum vz position', default=None)
parser.add_argument('--axmin', type=float, help='Minimum ax position', default=None)
parser.add_argument('--axmax', type=float, help='Maximum ax position', default=None)
parser.add_argument('--aymin', type=float, help='Minimum ay position', default=None)
parser.add_argument('--aymax', type=float, help='Maximum ay position', default=None)
parser.add_argument('--azmin', type=float, help='Minimum az position', default=None)
parser.add_argument('--azmax', type=float, help='Maximum az position', default=None)
parser.add_argument('--savepng', help='PNG filename to save plot to', default=None)

args = parser.parse_args()

datas=[]
info={}
data = read_data(args.filename,info)

alldata = data

# Find min and max
if not args.xmin:
  args.xmin = min([d['x_err'] for d in alldata])
if not args.xmax:
  args.xmax = max([d['x_err'] for d in alldata])
if not args.ymin:
  args.ymin = min([d['y_err'] for d in alldata])
if not args.ymax:
  args.ymax = max([d['y_err'] for d in alldata])
if not args.zmin:
  args.zmin = min([d['z_err'] for d in alldata])
if not args.zmax:
  args.zmax = max([d['z_err'] for d in alldata])

if not args.vxmin:
  args.vxmin = min([d['vx_err'] for d in alldata])
if not args.vxmax:
  args.vxmax = max([d['vx_err'] for d in alldata])
if not args.vymin:
  args.vymin = min([d['vy_err'] for d in alldata])
if not args.vymax:
  args.vymax = max([d['vy_err'] for d in alldata])
if not args.vzmin:
  args.vzmin = min([d['vz_err'] for d in alldata])
if not args.vzmax:
  args.vzmax = max([d['vz_err'] for d in alldata])

if not args.axmin:
  args.axmin = min([d['ax'] for d in alldata])
if not args.axmax:
  args.axmax = max([d['ax'] for d in alldata])
if not args.aymin:
  args.aymin = min([d['ay'] for d in alldata])
if not args.aymax:
  args.aymax = max([d['ay'] for d in alldata])
if not args.azmin:
  args.azmin = min([d['az'] for d in alldata])
if not args.azmax:
  args.azmax = max([d['az'] for d in alldata])

if not args.tmin:
  args.tmin = min([d['time'] for d in alldata])
if not args.tmax:
  args.tmax = max([d['time'] for d in alldata])

plot(data,xmin=args.xmin,xmax=args.xmax,ymin=args.ymin,ymax=args.ymax,zmin=args.zmin,zmax=args.zmax,
     vxmin=args.vxmin,vxmax=args.vxmax,vymin=args.vymin,vymax=args.vymax,vzmin=args.vzmin,vzmax=args.vzmax,
     axmin=args.axmin,axmax=args.axmax,aymin=args.aymin,aymax=args.aymax,azmin=args.azmin,azmax=args.azmax,
     tmin=args.tmin,tmax=args.tmax,
     savepng=args.savepng)
