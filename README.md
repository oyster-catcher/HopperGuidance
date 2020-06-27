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
