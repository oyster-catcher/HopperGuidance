#!/bin/bash


# Test minimum descent angle
scales="--xmax 300 --ymax 200 --vymin -30 --vymax 30 --numchecks 5"
for ang in 1 10 20 30;
  do mono ./Solve.exe r0=[0,160,0] v0=[100,-10,0] att0=[0,1,0] Tmin=1 Tmax=300 N=5 g=9.8 amin=0 amax=20 minDescentAngle=$ang > test1.$ang.dat
done
#./plotXYZ.py $scales test1.{1,10,20,30}.dat

# Test maximum velocity
scales="--xmax 1200 --ymax 2000 --vymin -150 --vymax 10 --numchecks 5"
for vmax in 40 60 80 100 150;
  do mono ./Solve.exe r0=[300,2000,0] v0=[0,0,0] att0=[0,1,0] Tmin=1 Tmax=300 N=5 g=9.8 amin=0 amax=30 vmax=$vmax > test2.$vmax.dat
done
#./plotXYZ.py $scales test2.{40,60,80,100,150}.dat

# Test limit of horizontal acceleration at landing
scales="--xmax 30 --ymax 20 --numchecks 5"
#TODO finalsideamax >10 doesn't work - seems to do final 0.5*(a+b) step and fail
for finalsideamax in 5 2 1 0;
  do mono ./Solve.exe r0=[100,100,0] v0=[0,-10,0] att0=[0,1,0] Tmin=1 Tmax=300 N=5 g=9.8 amax=15 finalHorizontalAMax=$finalsideamax > test3.$finalsideamax.dat
done
./plotXYZ.py $scales test3.{10,5,2,1,0}.dat
