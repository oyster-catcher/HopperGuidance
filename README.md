HopperGuidance v0.2.3
=====================

A KSP mod which enables the calculation of a fuel efficient trajectory in-flight through a number of targets to make a soft landing at a designated location respecting the capabilities of the craft. The craft will be steered automonously to follow this trajectory. You can choose many parameters of the trajectory to make it safe or more exciting. This has been a culmination of years of effort on and off, particularly to give to the ability to find the optimum path to steer quickly through a sequence of targets.

Features
- Pick targets to pass through ending on a ground target
- Shows the trajectory to landing and required thrust vectors
- Lets you tune the parameters to the craft

HopperGuidance is designed handle trajectories close to the landing site but you can see how far you can push it. It may lose accuracy beyond 10km. It can't handle descent from orbit well as it treats the planet surface as flat and doesn't take into account atmospheric drag, though can correct for it while flying. It may waste quite a bit of fuel trying to keep to the pre-calculated trajectory too. You can give it a try though, its all good fun.

The algorithm used is aimed to be an implementation of the G-FOLD algorithm for fuel optimal diverts
reputably used by Lars Blackmore for landing the SpaceX Falcon-9 rocket. However my algorithm is many ways a simplified version of it, though I've changed the constraints quite a bit and have extended it to compute trajectories though intermediate targets which involved some compromises, but it works pretty well.

Please see the videos from earlier versions here https://www.youtube.com/watch?v=Ekywp-P6EBQ&t=301
and https://www.youtube.com/watch?v=ouIQUfHsZac&t=101s
For the latest version, which is much improved, see https://www.youtube.com/watch?v=VJ_4bBZ_jT4

References

-  "G-FOLD: A Real-Time Implementable Fuel Optimal Large Divert Guidance Algorithm for Planetary Pinpoint Landing" https://www.researchgate.net/publication/258676350_G-FOLD_A_Real-Time_Implementable_Fuel_Optimal_Large_Divert_Guidance_Algorithm_for_Planetary_Pinpoint_Landing.

- "Convex Programming Approach to Powered Descent Guidance for Mars Landing" https://arc.aiaa.org/doi/abs/10.2514/1.27553

- "Minimum-Landing-Error Powered-Descent Guidance for Mars Landing Using Convex Optimization" http://www.larsblackmore.com/BlackmoreEtAlJGCD10.pdf

Pre-requistes
=============

Works with KSP versions 1.8.1 to 1.10 and possibly beyond.
If you want to run the solver outside of KSP you need python and matplotlib to show the trajectories.

Installing
==========

Copy GameData on top of the existing GameData folder in your KSP installation.

Using HopperGuidance
====================

The mod includes shiny cuboid HopperGuidanceCore part in the Command & Control category and two test craft, Test Craft and Test Craft Heavy. Right click to the part to show its UI.

Setting the target position. Click "Pick target" to select a ground position with the mouse pointer.
You can changes its height via the part dialog box which always sets the height of the last target.
A yellow target, inspired a little by the SpaceX drone ship markings will be shown. If the height is greater than zero then the centre cross is extended on a stick above the ground.

Once in flight click "Enable guidance" to automatically calculate a trajectory to safely land at the target point. A trajectory will be shown in green with the calculated thrust direction and magnitude in red. The craft will be steered to attempt to follow this trajectory. You may get the message "Failed: try again" in which case the craft is likely to be either: too far away, beyond max velocity or below the minimum descent angle. If you tweak these numbers the trajectory will be automatically recalculated. However, if the solution has failed because is the position of the craft, may be its too low than you will need to manually press click "Enable guidance" again.

You can cancel the autopilot with "Cancel guidance" and take control again. You might want to do this to recalculate the trajectory if the craft has diverged too much from it.

Control parameters. The algorithm is in two parts.

1. Calculating a trajectory
2. Following the trajectory

Parameters for (1) are used to compute the trajectory and can be used to control how 'tame' the trajectory is or how exciting is. If you want your craft to come slowly and gently down from above like Blue Origin you can choose the slightly more exciting Space X Falcon 9 landing, or more extreme doomed to failure trajectories coming in fast and very low.

- Min. descent angle - imagine an upside down cone with the apex at the landing point. The cone describes the safe area of descent. A high angle to the ground means the craft must descent steeply, only getting near the ground at the landing point.
- Max velocity - Keep velocity below this limit
- Max thrust angle - the maximum allow angle the craft will tilt. If 90 degrees then the craft can tip from vertical to the horizontal position but no more. If beyond >90 this value is ignored completely.

Click "Pick Target" then click on the main view to select the target, and then sets if height on the part UI if you wish. You will need to increase the target size to see far away targets. Note that the target is horizontal so it will get hidden on slopes.

Parameters for (2) describe what to do when the craft is off the trajectory. The nearest point on the trajectory is marked by a blue line from the craft to the trajectory. This nearest point takes into account and position and velocity with a little more weight for velocity. This makes the craft behave more smoothly rather than blindly aiming for the nearest point in position. So the craft tries to match the position and velocity of the nearest point. A correct to the thrust vector is calculated minimise the discrepancy. Six PID controllers are used to achieve this. 3 for the X,Y and Z positions and 3 for X, Y and Z velocities.

- Idle attitude angle - If the craft is pointing more this angle away from the required direction of thrust then idle the engine at 1% thrust. This prevent the craft thrusting in the wrong direction and waits of it to point correctly
- Correction factor - Sets the velocity to aim for to close the position error. If set of 1 then when 1m aware aim for a velocity of 1m/s towards the target. 0.15 to 0.25 are good values. Use down to 0.15 for large vessels that take longer to change their attitude. Lower this if the craft oscillates around the trajectory.

Finally there are a few toggle switches

- Keep engine ignited - Use in Realism Overhaul where engines have limited ignition and take time so start up
- Show track - Display track, align and steer vectors
- Logging - If enabled the files vessel.dat, target.dat and solution.dat are logged to the same directory as KSP.log

Here's a diagram to explain some of the parameters

![](docs/Diagram.png)

General tips:

- Have sufficient RCS or reaction wheel for your craft so it can change its attitude quickly enough (look at whether the steer vector is often misaligned with the crafts attitude)
- Enable RCS
- It can take longer to compute highly constrained trajectories for heavy craft or not find a solution at all. Retry, try something easier or reduce the number of targets.
- If your craft comes in too low to land either reduce the thrust angle to land more vertically or increase min descent angle
- If your craft is way off the trajectory its best to cancel and re-enable guidance again to compute a new trajectory

Failure to find solution:

Usually caused by one of the following

- The craft is going faster than max velocity
- The craft is below the minimum descent angle (to final target)
- You have run out of fuel
- Its impossible not to hit the ground
- The targets are too close and near impossible to fly between
- There are too many targets. About 10 is the maximum
- You have run out of ignitions (Realism Overhaul)
- Other. Sadly sometimes the solver just fails to find a solution. Simplify things or just try again

Fun things to try
=================

Position targets for a slalom course around some buildings
Find the bridges at KSC and try to fly under them
Try and reproduce the SN5 Starship hop in Realism Overhaul
Try and reproduce a Falcon 9 landing in Realism Overhaul
Try setting targets round and round in a circle to make a loop

Dealing with Small Craft
========================

Small craft are exciting, like a sports car.
You can try and push of max velocity, max thrust angle to high values. The craft will then be able to move quickly and accelerate quickly leads to exciting but possibly risky landings.

Dealing with Large Craft
========================

Large craft are hard to maneuver since they have high mass and will take time to change its attitude. This leads to it pointing is the wrong direction and either thrusting in that wrong direction or waiting to change its attitude. Both these mean it cannot following the track. Another problem is that atmospheric drag is not included when solving for the best trajectory. So here's some tips for dealing with large craft

- Keep max thrust angle low, says <15 degrees as this will keep the craft more upright and mean the craft is never pointing in the wrong direction for long
- Lower correction factor down to as low as 0.15. This will mean if the craft is off trajectory it will only move slowly towards the correct trajectory. If this value is too high it will overshoot. Observe whether the red steer attitude is far from the actual crafts attitude
- You can compensate for the lack of ability to steer the craft by adding more RCS thrusters or overly powerful reaction wheels
- Don't position targets too close to each other, especially if large velocity changes are required

Dealing with Realism Overhaul:

- Choose an engine that is throttleable
- Switch on keepIgnited otherwise you may use up all your ignitions very quickly
- Note that the engine has a minimum thrust and if this results in craft rising then you have too much thrust to hover. So choose your engine carefully
- Reaction wheels don't work too well and RCS might be limited. Thrust vectoring of the main engine will be needed to steer the craft, but that can work well if the engine is always ignited and has enough thrust

Hints for a Falcon 9 landing:

HopperGuidance will only use the currently activated engines, as it can't read the min/max thrust otherwise. With engines that can't throttle to 0% this mean they will be providing thrust which is a nuisance. To emulate a Falcon 9 landing as in the latest video I shutdown the engines while falling and enable them at minimum thrust at 2500m with a brief press of SHIFT. This only slows the rocket at about 2m/s/s and then click "Enable Guidance". Also note that with 3 merlin engines the minimum thrust is greater than gravity so a Falcon 9 can't hover or descend slowly with three engines. If you find you are hovering above the ground you will need to cancel guidance and cut the engine otherwise you will rise again. In reality SpaceX only use a single engine for the final part of the descent, shutting down the other two, and I plan to implement this in a future version.

Debugging
=========

If you have a problem look in the KSP directory for KSP.log. Log lines are proceeded by [HopperGuidance]
You will be able to find all the parameters when either a solution was requested and found or it failed. You can rerun by taking the parameters from the log line and using them directly on the command line of Solve.exe outside of KSP. See below.

Testing
=======

There is a way to run the algorithm outside of KSP via Solve.exe.
I'm running on the Mac so I run
```
mono ./Solve.exe r0=[-340,100,20] v0=[10,10,10] target=[-200,120,0] target=[150,50,0] rf=[0,0,0] vf=[0,0,0] g=9.8 maxThrustsBetweenTargets=2 > solution.dat; ./plotXYZ.py solution.dat
```
for example and this produces a data file with a trajectory.
This can be plotted (via matplotlib) by running

```
./plotXYZ.py --square solution.dat

sage: plotXYZ.py [-h] [--xmin XMIN] [--xmax XMAX] [--ymin YMIN] [--ymax YMAX]
                  [--zmin ZMIN] [--zmax ZMAX] [--vxmin VXMIN] [--vxmax VXMAX]
                  [--vymin VYMIN] [--vymax VYMAX] [--vzmin VZMIN]
                  [--vzmax VZMAX] [--tmax TMAX] [--amult AMULT] [--square]
                  [--showchecks]
                  filename [filename ...]

Plot vessel data logs (or solutions) with X,Y,Z,VX,VY,VZ and Throttle in
multiple plots

positional arguments:
  filename       Filename of TAB-separated data file, first line contains
                 column names. Should contain time,x,y,z,vx,vy,vz,ax,ay,ax

optional arguments:
  -h, --help     show this help message and exit
  --xmin XMIN    Minimum x position
  --xmax XMAX    Maximum x position
  --ymin YMIN    Minimum y position
  --ymax YMAX    Maximum y position
  --zmin ZMIN    Minimum z position
  --zmax ZMAX    Maximum z position
  --vxmin VXMIN  Minimum vx position
  --vxmax VXMAX  Maximum vx position
  --vymin VYMIN  Minimum vy position
  --vymax VYMAX  Maximum vy position
  --vzmin VZMIN  Minimum vz position
  --vzmax VZMAX  Maximum vz position
  --tmax TMAX    Maximum time
  --amult AMULT  Multiplier for scale up thrust acceleration lines
  --square       Make XY plot square (roughly as depends on window size)
  --showchecks   Show time checks for max vel. and min descent angle
```

![](docs/plotXYZ.png)

Note: The maximum acceleration is shown as the dotted line on the time vs mag(accel) plot and this is exceeded temporarily. The craft will not be able to produce this much acceleration and the craft with diverge from the trajectory. This is because the constraint is actually applied independently to each X, Y, Z axis so the total magnitude can be larger than the axis limit. This is a simplification to help speed up solving time but there is scope for improvement.

You can also run all the tests by running
```
./run_tests.sh
and
./run_test2.sh # multi-target tests
```
This is all designed to run on a Mac where you have a bash shell. It should be failing easy to convert though.

Logging
=======

If you switch on logging in the part UI window then the file vessel.dat and solution.dat with by logged to the KSP directory. solution.dat is what would be generated by Solve.exe but running inside the Game. vessel.dat is a log of the actual vessel parameters while trying to follow the trajectory.

Heres an example:

![](docs/logging.png)

You can see the vessel followed trajectory reasonably well though couldn't make the turn half way through. The vessel is Test Craft Heavy and the problem is that its attitude could not match that required quickly enough so it didn't slowdown faster enough. Look at the attitude error plot.

Compilation
===========

Running on a Mac I used Mono to compile C#.
I couldn't get the Mac Visual Studio to use the correct version of Mono for KSP.
But I did manage to get it to work with the Makefile.
Just type
```
make
```
But I include the Visual Studio project as it may be easy to get running on a Windows machine or on a Mac with some configuration hacks. Any help on this much appreciated.

ALGLIB
======

I make use of the free version of ALGLIB for convex optimisation. I use the free version provided until the GPL. Here's a description of ALGLIB from alglib.net

ALGLIB is a cross-platform numerical analysis and data processing library. It supports several programming languages (C++, C#, Delphi) and several operating systems (Windows and POSIX, including Linux). ALGLIB features include:

- Data analysis (classification/regression, statistics)
- Optimization and nonlinear solvers
- Interpolation and linear/nonlinear least-squares fitting
- Linear algebra (direct algorithms, EVD/SVD), direct and iterative linear solvers
- Fast Fourier Transform and many other algorithms
- ALGLIB Project offers you several editions of ALGLIB:

Details of the Algorithm
========================

My implementation is not quite the same as the 'proper' G-FOLD algorithm. I made some simplifications due to the convex solver code I had available and in some cases due to my lack of understanding.
The full solution makes account of the consumption of fuel and the improved acceleration when mass drops due to this. Mine does not. If you run out of fuel then sorry.

The proper G-FOLD uses a Second Order Cone Programming (SOCP) solver. This is fancy stuff, but basically it means that various constraints such a thrust directions and minimum descent angle are all cones. You can see that they are since they are a maximum angle from the vertical in all directions.
I couldn't find a SOCP solver but I could find a Quadratic Programming (QP) solver in the shape of alglib. See alglib.net. This means I use 4 planes rather than a cone so its an approximation that you probably won't notice.

To follow the trajectory I use 6 PID controllers. There are a pair of controller for each X, Y and Z. The controllers are very simple using only the proportional element, kP. The 1st controller say for the X axis sets a target (or setpoint) for the velocity, VX given the error in the X position. To start with we find the closest position on the trajectory by considering both position and velocity of the craft and trajectory and finding the nearest point. kP is given by "Correction factor". To this we add the actual desired velocity from the nearest point on the trajectory. So this says aim to move towards the X position at this speed while also moving along the trajectory. Then the second PID for X set the acceleration in the X direction given the mismatch between the velocity desired given the target velocity and the current velocity. kP in this controller is set from "Acceleration gain". This is done for the three dimensions and gives a desired correction thrust vector. This correction is added to the thrust vector already computed from the solution. Finally the thrust vector in limited by the max thrust angle from vertical which will prevent the craft leaning too much. The direction of this thrust vector is shown by the red bar drawn from the centre of the craft.

The recent changes in v0.2.3 was to enable multiple targets. Each target has a position whereas the final target also has a velocity constraint of a small vertical descent. The G-FOLD algoritm has to test out solutions over a range of durations with a golden search. But with multiple targets there are multiples times so the search would quickly become high dimensional and take for other. To get around this the solution is done iteratively, computing the solution to the first target with minimum fuel and then pinning down the time of the first target. Then the next target is adding, trying out a range of durations to travel to it, and continuing. The trajectory from the start point to target 1 will generally needs changing once target 2 is known so the whole trajectory is computed again as each target is added, but there is a problem if target 2 requires a radical change of direction it might mean the craft can't get to target 1 so quickly, so a slight extra time allowance is the made to get to each target to allow for some slack, and usually more fuel can be burned to get to the target more quickly. So the multi-target solution was not always be optimal but generally it works well, if sometimes taking a few seconds to compute.

Adrian Skilling

adrian.skilling@gmail.com
