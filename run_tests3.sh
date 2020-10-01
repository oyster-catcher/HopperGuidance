#!/bin/bash -e
# Tests for Falcon 9 landing

args="--square --showchecks"

# Simple slow landing on single central engine. TWR=1.8 (very rough)
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-30,500,0] v0=[0,-50,0] g=9.8 Tmin=-1 Tmax=-1 amin=4 amax=16 vmax=250 minDescentAngle=80 maxThrustAngle=10 maxLandingThrustAngle=5 rf=[0,0,0] vf=[0,-1,0] > s1.dat
./plotXYZ.py $args s1.dat

# Fast landing on 3 engine. TWR=4.5 (very rough). Min TWR>1 so impossible to come in slowly as can't hover
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-30.00,2500.00,0.00] v0=[0.00,-400.00,0.00] g=9.8 Tmin=0 Tmax=67.2460403442383 amin=7 amax=44 vmax=750 minDescentAngle=80 maxThrustAngle=10 maxLandingThrustAngle=5 rf=[0.00,0.00,0.00] vf=[0.00,-1.00,0.00] extraTime=0.5 > s2.dat
./plotXYZ.py $args s2.dat

# Same as above but more off-target
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-200.00,2500.00,0.00] v0=[0.00,-400.00,0.00] g=9.8 Tmin=0 Tmax=67.2460403442383 amin=11 amax=44 vmax=750 minDescentAngle=80 maxThrustAngle=10 maxLandingThrustAngle=5 rf=[0.00,0.00,0.00] vf=[0.00,-1.00,0.00] extraTime=0.5 > s3.dat
./plotXYZ.py $args s3.dat

# Same as above but reduce maxLandingAngle
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-200.00,2500.00,0.00] v0=[0.00,-400.00,0.00] g=9.8 Tmin=0 Tmax=67.2460403442383 amin=11 amax=44 vmax=750 minDescentAngle=80 maxThrustAngle=10 maxLandingThrustAngle=1 rf=[0.00,0.00,0.00] vf=[0.00,-1.00,0.00] extraTime=0.5 > s4.dat
./plotXYZ.py $args s4.dat

# Now use 5 engines for VERY high deceleration but min TWR very high which means this solution is only possible
# when altitude and vertical velocity are in a very narrow range since engines must be fired for the whole descent
for h in 4000 3400 2800 2400 2000 1900 1850 1800;
 do mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-100.00,$h,30] v0=[0.00,-400.00,0.00] g=9.8 Tmin=0 Tmax=67.2460403442383 amin=19 amax=66 vmax=750 minDescentAngle=80 maxThrustAngle=10 maxLandingThrustAngle=5 rf=[0.00,0.00,0.00] vf=[0.00,-1.00,0.00] extraTime=0.5 > s5.h=$h.dat
  ./plotXYZ.py $args s5.h=$h.dat
done

# Failure from ~/KSP_RO/KSP.log. amin=37! with 5 engines
mono ./SolveTest.exe tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-72.31,2028.23,88.66] v0=[5.08,-154.30,2.49] g=9.8137606716112 Tmin=0 Tmax=54.4644813537598 amin=37.6553230544797 amax=76.7035334738888 vmax=511 minDescentAngle=75 maxThrustAngle=7.19999997317791 maxLandingThrustAngle=5 extraTime=0.5  rf=[0.08,20.12,0.03] vf=[0.00,-1.40,0.00]
