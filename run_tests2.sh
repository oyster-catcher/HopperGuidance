#!/bin/bash -e
# More difficult cases with multi-target trajectories and seen failing within KSP

args="--showchecks"

# Large craft seems to want to make a long route
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[30,30,0] v0=[0.00,0.00,0.00] g=9.80735028176503 Tmin=1 Tmax=260.169830322266 amin=0 amax=18 vmax=51 minDescentAngle=10 maxThrustAngle=45 maxLandingThrustAngle=2.00000002235174 target=[20,20,0] rf=[0.00,0,0.00] vf=[0.00,-1.40,0.00] > s2a.dat
./plotXYZ.py $args s2a.dat

# First part
mono ./SolveTest.exe tol=1 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[0.00,1.00,0.00] v0=[0.00,0.00,0.00] g=9.80735028176503 Tmin=1 Tmax=260.169830322266 amin=0 amax=18 vmax=51 minDescentAngle=10 maxThrustAngle=15 maxLandingThrustAngle=2 rf=[100.00,100.00,0.00] > s2b.dat
./plotXYZ.py $args s2b.dat

# Failure to compute trajectory when on ground
mono ./SolveTest.exe tol=0.1 minDurationPerThrust=4 maxThrustsBetweenTargets=4 r0=[0,1,0] v0=[0.00,0.00,0.00] g=9.8 Tmin=1 Tmax=40 amin=0 amax=18 vmax=51 minDescentAngle=10 maxThrustAngle=15 maxLandingThrustAngle=2 target=[100,100,0] rf=[0.00,0,0.00] vf=[0.00,-1.40,0.00] > s2c.dat
./plotXYZ.py $args s2c.dat

# Route to get to height and travel to 3 targets at same height then drop
# Unfortunately seems to get almost to landing then go high in order to reach the right maxLandingThrustAngle constraint
mono ./SolveTest.exe tol=0.2 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[6.41,11.82,32.82] v0=[0.00,0.00,0.00] g=9.80740002422221 Tmin=1 Tmax=289.309143066406 amin=0 amax=10.1867346139428 vmax=51 minDescentAngle=10 maxThrustAngle=19.9999999254942 maxLandingThrustAngle=4.99999998137355 target=[92.96,80.15,57.46] target=[128.31,79.73,159.57] target=[-34.09,80.08,195.12]  rf=[0.00,7.83,0.00] vf=[0.00,-1.40,0.00] > s3d.dat
./plotXYZ.py $args s3d.dat

# Failed SN5 trajectory
# looks like problem was previous uncleared targets
mono ./SolveTest.exe minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[125.19,11.83,36.43] v0=[0.00,0.00,0.00] g=9.80740369934396 Tmin=1 Tmax=632.105651855469 amin=0 amax=10.1867175012821 vmax=51 minDescentAngle=20 maxThrustAngle=19.9999999254942 maxLandingThrustAngle=4.99999998137355 rf=[0.00,7.80,0.00] vf=[0.00,-1.40,0.00] > s3e.dat
./plotXYZ.py $args s3e.dat

# Tmax was too high, even setting it to 40 was too high even though final solution is 15
mono ./SolveTest.exe minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[125.19,110.83,36.43] v0=[0.00,0.00,0.00] g=9.80740369934396 Tmin=1 Tmax=35 amin=0 amax=10.1867175012821 vmax=51 minDescentAngle=20 maxThrustAngle=19.9999999254942 maxLandingThrustAngle=20 rf=[0.00,7.80,0.00] vf=[0.00,-1.40,0.00] > s3f.dat
./plotXYZ.py $args s3f.dat

# Failure with Test Craft Heavy
# Low trajectory with excessiive extra long loop
mono ./SolveTest.exe tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[67.64,13.33,19.77] v0=[0.00,0.00,0.00] g=9.80735191908658 Tmin=295.046055189007 Tmax=295.984191143588 amin=0 amax=10.1700578514072 vmax=51 minDescentAngle=5.12181663513184 maxThrustAngle=19.9999999254942 maxLandingThrustAngle=4.99999998137355 target=[821.83,136.52,779.85] target=[-556.90,144.91,1406.68]  rf=[0.00,9.31,0.00] vf=[0.00,-1.40,0.00] > s3g.dat
./plotXYZ.py $args s3g.dat

# Failure with test craft on mun to close landing site
mono ./SolveTest.exe tol=0.5 minDurationPerThrust=4 maxThrustsBetweenTargets=3 r0=[-0.65,23.31,-7.04] v0=[-0.16,16.49,-5.28] g=1.53671856318373 Tmin=10.5800775593291 Tmax=17.1189250946045 amin=0 amax=42.7053677878649 vmax=150 minDescentAngle=20 maxThrustAngle=63.9999997615814 maxLandingThrustAngle=8.00000011920929 extraTime=0.5  rf=[0.00,4.97,0.00]:-1 vf=[0.00,-1.40,0.00]:-1
