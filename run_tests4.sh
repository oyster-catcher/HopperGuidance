#!/bin/bash

# Why is v final still 20 m/s !!

# Try and range of start times
mono SolveTest.exe tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=0 r0=[0,6000,0] v0=[0,-200,0] g=9.8 Tmin=0 Tmax=82 amin=6 amax=14 vmax=500 minDescentAngle=75 maxThrustAngle=13 maxLandingThrustAngle=5 extraTime=0.8425344 rf=[0,0,0] vf=[0.00,0.00,0.00] TstartMax=0 TstartMin=0 > s.dat
tail s.dat

