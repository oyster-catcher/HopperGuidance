#!/bin/bash


# Test minimum descent angle
scales="--xmax 300 --ymax 200 --vymin -30 --vymax 30 --numchecks 5"
for ang in 1 10 20 30;
  do mono ./Solve.exe r0=[0,160,0] v0=[100,-10,0] Tmin=1 Tmax=300 N=5 g=9.8 amin=0 amax=20 minDescentAngle=$ang > test1.minDescentAngle=$ang.dat
done
./plotXYZ.py $scales test1.minDescentAngle={1,10,20,30}.dat

# Test maximum velocity
scales="--xmax 1200 --ymax 2000 --vymin -150 --vymax 10 --numchecks 5"
for vmax in 40 60 80 100 150;
  do mono ./Solve.exe r0=[300,2000,0] v0=[0,0,0] Tmin=1 Tmax=300 N=5 g=9.8 amin=0 amax=30 vmax=$vmax > test2.vmax=$vmax.dat
done
./plotXYZ.py $scales test2.vmax={40,60,80,100,150}.dat

# Test limit of horizontal acceleration at landing
scales="--xmax 30 --ymax 20 --numchecks 5 --amult 0.5"
for maxLandingThrustAngle in 0 10 20 40;
  do mono ./Solve.exe r0=[100,100,0] v0=[0,-10,0] Tmin=1 Tmax=300 N=5 g=9.8 amax=15 maxLandingThrustAngle=$maxLandingThrustAngle > test3.maxLandingThrustAngle=$maxLandingThrustAngle.dat
done
./plotXYZ.py $scales test3.maxLandingThrustAngle={0,10,20,40}.dat

# Test max thrust angle
scales="--xmax 1200 --ymax 900 --numchecks 5"
for ang in 10 20 60 120;
  do mono ./Solve.exe r0=[400,400,0] v0=[40,10,0] Tmin=1 Tmax=300 N=5 g=9.8 amax=15 maxThrustAngle=$ang maxLandingThrustAngle=30 > test4.maxThrustAngle=$ang.dat
done
./plotXYZ.py $scales test4.*.dat
