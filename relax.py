#!/usr/bin/python
#
# Call Solve.exe but relaxing constraints to try and find solution
#
# tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-100,828.23,0] v0=[0,-300,0] g=9.8137606716112 Tmin=0 Tmax=54.4644813537598 amin=37.6553230544797 amax=76.7035334738888 vmax=511 minDescentAngle=75 maxThrustAngle=20 maxLandingThrustAngle=5 extraTime=0  rf=[0,0,0]:-1 vf=[0.00,-1.40,0.00]:-1

import sys
import os

xfield="maxThrustAngle"
xrange=range(5,21,1)
yfield="maxLandingThrustAngle"
yrange=range(0,10,1)

for x in xrange:
  args=[]
  for kv in sys.argv[1:]:
    k,v = kv.split("=",1)
    if k==xfield:
      kv = "%s=%f" % (xfield,x)
    args.append(kv)
  args=" ".join(args)
  print("# %s" % args)
  ret = os.system("mono ./Solve.exe %s > tmp.dat" % args)
  print(ret)
