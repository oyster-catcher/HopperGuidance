#!/bin/bash
# Test VesselSim
#  This is a 'fuller' simulation of the vessel trying to follow the trajectory
#  is uses the PID controllers and the Controller class which will also implement
#  the special final descent model set when distance to final target < finalDescentDistance

# Seem to have a problem steering to trajectory when well off it

# This is OK
#mono ./VesselSim.exe r0offset=[100,0,0] r0=[-500,1000,0] rf=[0,0,0] vf=[0,-2,0] maxT=250 finalDescentDistance=50 touchdownSpeed=2 kp1=0.2 kp2=0.6 > s1.dat

# This is not
# Problem is because target position is at r0 when stationary
# vessel tries to get to that position but gets stuck with a balance in position, velocity and
# acceleration targets that balance out to leave vessel hovering
#mono ./VesselSim.exe r0offset=[200,0,0] r0=[-500,1000,0] rf=[0,0,0] vf=[0,-2,0] maxT=250 finalDescentDistance=50 touchdownSpeed=2 kp1=0.2 kp2=0.6 > s2.dat

# This is OK because falling
#mono ./VesselSim.exe r0offset=[200,0,0] v0offset=[0,-20,0] r0=[-500,1000,0] rf=[0,0,0] vf=[0,-2,0] maxT=250 finalDescentDistance=50 touchdownSpeed=2 kp1=0.2 kp2=0.6 > s3.dat

#./plotXYZ.py --square solution.dat s1.dat s2.dat s3.dat

# Try and reproduce oscillation
mono ./VesselSim.exe tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[17.34,1101.97,49.52] v0=[1.27,20.48,3.96] g=9.81659694510968 Tmin=0 Tmax=56.3245032757322 amin=3.31267338625936 amax=16.1175997992514 vmax=150 minDescentAngle=20 maxThrustAngle=11.9999999552965 maxLandingThrustAngle=5 extraTime=0.6756973 rf=[0.49,23.22,0.04] vf=[0.00,-2.50,0.00] > vessel.dat

./plotXYZ.py --square solution.dat vessel.dat
