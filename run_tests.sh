#!/bin/bash


# Test minimum descent angle
for ang in 1 10 20 30;
  do mono ./SolveTest.exe r0=[200,160,0] v0=[0,-10,0] Tmin=1 Tmax=30 g=9.8 amin=0 amax=20 minDescentAngle=$ang rf=[0,0,0] vf=[0,0,0] maxThrustsBetweenTargets=4 > test1.minDescentAngle=$ang.dat
done
./plotXYZ.py --showchecks --square test1.minDescentAngle={1,10,20,30}.dat
exit

# Test maximum velocity
for vmax in 40 60 80 100 150;
  do mono ./SolveTest.exe r0=[300,2000,0] v0=[0,0,0] Tmin=1 Tmax=300 g=9.8 amin=0 amax=30 vmax=$vmax rf=[0,0,0] vf=[0,0,0] > test2.vmax=$vmax.dat
done
./plotXYZ.py --square test2.vmax={40,60,80,100,150}.dat

# Test limit of horizontal acceleration at landing
for maxLandingThrustAngle in 0 10 20 40;
  do mono ./SolveTest.exe r0=[100,100,0] v0=[0,-10,0] Tmin=1 Tmax=300 g=9.8 amax=15 maxLandingThrustAngle=$maxLandingThrustAngle rf=[0,0,0] vf=[0,0,0] > test3.maxLandingThrustAngle=$maxLandingThrustAngle.dat
done
./plotXYZ.py --square test3.maxLandingThrustAngle={0,10,20,40}.dat

# Test max thrust angle
for ang in 10 20 60 120;
  do mono ./SolveTest.exe r0=[400,400,0] v0=[40,10,0] Tmin=1 Tmax=300 g=9.8 amax=15 maxThrustAngle=$ang maxLandingThrustAngle=30 rf=[0,0,0] vf=[0,0,0] > test4.maxThrustAngle=$ang.dat
done
./plotXYZ.py --square test4.*.dat

# Test min thrust
for amin in 0 2 4 6 8;
  do mono ./SolveTest.exe r0=[400,400,0] v0=[40,10,0] Tmin=1 Tmax=300 g=9.8 amin=$amin amax=15 maxLandingThrustAngle=30 rf=[0,0,0] vf=[0,0,0] > test5.amin=$amin.dat
done
./plotXYZ.py --square test5.*.dat

# Test multiple targets
for amin in 0 2 4 6 8;
  do mono ./SolveTest.exe r0=[400,400,0] v0=[40,10,0] Tmin=1 Tmax=300 g=9.8 amin=$amin amax=15 maxLandingThrustAngle=30 rf=[0,0,0] vf=[0,0,0] > test5.amin=$amin.dat
done
./plotXYZ.py --square test5.*.dat
