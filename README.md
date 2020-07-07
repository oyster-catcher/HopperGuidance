HopperGuidance
==============

A KSP mod which enables can calculate an optimal trajectory in-flight which attempts to minimise fuel usage
to make a precision landing. The craft will then be steered automonously to follow this trajectory.

Features
- Choose a landing latitude, longitude and altitude
- Shows the trajectory to landing and require thrust vectors
- Lets you tune the parameters to the craft

You can choose a particular latitude and longitude and altitude to land at, although currently this is
designed to be used close to the landing site. It certainly can't currently handle descent from orbit
as it suffers from Flat Earthism and treats the planet surface as flat so it can't handle orbits or
atmospheric drag. It has however great for getting to your landing site within say 1km away. I plan
to extend this. Currently it works best for small craft which can steer quickly, I'm working on this too.

The algorithm used is aimed to be an implementation of the G-FOLD algorithm for fuel optimal diverts
reputably used by Lars Blackmore for landing the SpaceX Falcon-9 rocket. However this algorithm
currently has the following limitations.

- It used the Quadratic Solver in alglib rather than a Second Order Cone Programming (SOCP) solution.
  which means!! its not quite fuel optimal since the squared magnitude of rocket thrust vectors is
  minimized for each trajectory of a given duration rather than the sum of magnitudes (SOCP lets you
  use magnitudes)
- It doesn't to work out whether you have enough fuel to reach the landing point so you might run out.

References
-  "G-FOLD: A Real-Time Implementable Fuel Optimal Large Divert Guidance Algorithm for Planetary Pinpoint Landing" https://www.researchgate.net/publication/258676350_G-FOLD_A_Real-Time_Implementable_Fuel_Optimal_Large_Divert_Guidance_Algorithm_for_Planetary_Pinpoint_Landing.

- "Convex Programming Approach to Powered Descent Guidance for Mars Landing" https://arc.aiaa.org/doi/abs/10.2514/1.27553

- "Minimum-Landing-Error Powered-Descent Guidance for Mars Landing Using Convex Optimization" http://www.larsblackmore.com/BlackmoreEtAlJGCD10.pdf

Prerequistes
============

Tested with KSP version 1.8.1. May work with later versions.

Installing
==========

Copy GameData to your KSP installation

Using HopperGuidance
====================

Creates a new part, HopperGuidanceCore. Right-click to get its UI up.

Setting the target position. Click "Set Target Position" to set and landing target to where you are now.
A flickery red cross will be shown at that position.

Once in flight click "Enable autopilot" to automatically calculate a trajectory to safely land at the target point. Its not recommended that you do this more than 5km from the target currently. A trajectory will be shown in green with the calculated thrust direction and magnitude in red. The craft will be steered to attempt to follow this trajectory. You may get the message "Failed: try again" in which case the craft is likely to be either: too far away, beyond max velocity or below the minimum descent angle. Try to resolve these and try again.

You can cancel the autopilot with "Disable autopilot" and take control again. You might want to do this to recalculate the trajectory if the craft has diverged too much from it.

Control parameters. The algorithm is in two parts.

1. Calculating a trajectory
2. Following the trajectory

Parameters for (1) are used to compute the trajectory and can be used to control how 'tame' the trajectory is or how exciting is. If you want your craft to come slowly and gently down from above like Blue Origin you can choose the slightly more exciting Space X Falcon 9 landing, or more extreme doomed to failure trajectories coming in fast and very low.

- Target latitude
- Target longitude
- Target altitude
- Final max side acceleration - set this low or zero and the craft must descent vertically or increase to allow descent from the side
- Min. descent angle - imagine an upside down cone with the apex at the landing point. The cone describes to safe area of descent. A high angle to the ground means the craft must descent steeply, only getting near the ground at the landing point.
- Max thrust - how much of the available thrust to use. Set this slow and the descent will be slow and careful. You can use more than 100% since the calculation is done when calculating the trajectory you might have lost weight due to fuel usage and you will have spare acceleration by landing.
- Max velocity - Keep velocity below this limit

Parameters for (2) describe what to do when the craft is off the trajectory. The nearest point on the trajectory is marked by a blue line from the craft to the trajectory. This nearest point takes into account and position and velocity with a little more weight for velocity. This makes the craft behave more smoothly rather than blindy aiming for the nearest point in position. So the craft tries to match the position and velocity of the nearest point. If calculates a correct to the thrust vector to try and minimise the discrepancy. Six PID controllers are used to achieve this. 3 for the X,Y and Z positions and 3 for X, Y and Z velocities.

- Idle attitude angle - If the craft is pointing more this angle away from the required direction of thrust then idle the engine at 0.1% thrust. This prevent the craft thrusting in the wrong direction and waits of it to point correctly
- Err. Position gain - Sets the velocity to aim to close the position error. If set of 1 then when 1m aware aim for a velocity of 1m/s towards the target. 0.2 to 0.4 are good values.
- Velocity gain - If set of 1 then when the target velocity is out by 1m/s than accelerate at 1m/s2
- Err. max thrust - How much thrust to use to minimise the error


Issues:





Testing
=======
