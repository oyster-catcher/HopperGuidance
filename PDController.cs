// Powered Descent Controller
//   - calculates a trajectory to follow using GFold algorithm and attempts to follow it

using System;
using System.Collections.Generic;

using UnityEngine;
using UnityEngine.UI;

namespace HopperGuidance
{
  public class PDController
  {
    // Parameters
    public PID3d _pid3d = new PID3d();
    public float maxThrustAngle = 45;
    public float idleAngle = 90;
    public float amin = 0;
    public float amax = 30;
    public float touchdownSpeed = 1;
    public float powerDescentAlt = 5000; // Powered descent below this altitude
    public float aeroDescentAlt = 70000; // Aerodynamic descent below this altitude
    public float finalDescentAlt = 0; // Ignore trajectory under this distance to final landing spot
    public double attitudeTimeConstant = 0; // Instanteous - speed of attitude change
    public double throttleTimeConstant = 0; // Instanteous - speed of throttle change

    // Live values
    public bool ignited = false;
    public Vector3d currentThrustDir = Vector3d.zero;
    public double currentThrustMag = 0;
    float last_t = 0;
    bool finalDescent = false;

    public PDController(float kp1=1, float ki1=0, float kd1=0, float kp2=1, float ki2=0, float kd2=0,
                        float kp1y=1, float ki1y=0, float kd1y=0,
                        float kp2y=1, float ki2y=0, float kd2y=0,
                        float a_amin=0, float a_amax=30, float a_maxThrustAngle=45, string logFilename="")
    {
      _pid3d = new PID3d(kp1,ki1,kd1,kp2,ki2,kd2,kp1y,ki1y,kd1y,kp2y,ki2y,kd2y,1000,amax,logFilename);
      amin = a_amin;
      amax = a_amax;
      maxThrustAngle = a_maxThrustAngle;
      ignited = false;
      finalDescent = false;
    }

    public void ResetPIDs(float kp1, float ki1, float kd1, float kp2, float ki2, float kd2,
                          float kp1y, float ki1y, float kd1y,
                          float kp2y, float ki2y, float kd2y)
    {
      _pid3d.kp1 = kp1;
      _pid3d.ki1 = ki1;
      _pid3d.kd1 = kd1;
      _pid3d.kp1y = kp1y;
      _pid3d.ki1y = ki1y;
      _pid3d.kd1y = kd1y;
      _pid3d.kp2 = kp2;
      _pid3d.ki2 = ki2;
      _pid3d.kd2 = kd2;
      _pid3d.kp2y = kp2y;
      _pid3d.ki2y = ki2y;
      _pid3d.kd2y = kd2y;
      _pid3d.Reset();
    }

    public double MinHeightAtMinThrust(double y, double vy,double g)
    {
      double minHeight = 0;
      double lower_amin = amin*0.8;
      if (lower_amin < g)
        return -float.MaxValue;
      double tHover = -vy/lower_amin; // time to come to hover
      minHeight = y + vy*tHover + 0.5*lower_amin*tHover*tHover - 0.5*g*tHover*tHover;
      return minHeight;
    }

    public bool Set(string k, string v) {
      if (k=="kp1")
        _pid3d.kp1 = (float)Convert.ToDouble(v);
      else if (k=="ki1")
        _pid3d.ki1 = (float)Convert.ToDouble(v);
      else if (k=="kd1")
        _pid3d.kd1 = (float)Convert.ToDouble(v);
      else if (k=="kp2")
        _pid3d.kp2 = (float)Convert.ToDouble(v);
      else if (k=="ki2")
        _pid3d.ki2 = (float)Convert.ToDouble(v);
      else if (k=="kd2")
        _pid3d.kd2 = (float)Convert.ToDouble(v);
      //TODO: Add y PIDs params
      else if (k=="kp1y")
        _pid3d.kp1y = (float)Convert.ToDouble(v);
      else if (k=="ki1y")
        _pid3d.ki1y = (float)Convert.ToDouble(v);
      else if (k=="kd1y")
        _pid3d.kd1y = (float)Convert.ToDouble(v);
      else if (k=="kp2y")
        _pid3d.kp2y = (float)Convert.ToDouble(v);
      else if (k=="ki2y")
        _pid3d.kp2y = (float)Convert.ToDouble(v);
      else if (k=="kd2y")
        _pid3d.kp2y = (float)Convert.ToDouble(v);
      else if (k=="amin")
        amin = (float)Convert.ToDouble(v);
      else if (k=="amax")
        amax = (float)Convert.ToDouble(v);
      else if (k=="maxThrustAngle")
      {
        maxThrustAngle = (float)Convert.ToDouble(v);
      }
      else if (k=="finalDescentAlt")
        finalDescentAlt = (float)Convert.ToDouble(v);
      else if (k=="touchdownSpeed")
        touchdownSpeed = (float)Convert.ToDouble(v);
      else if (k=="throttleTimeConstant")
        throttleTimeConstant = (float)Convert.ToDouble(v);
      else if (k=="attitudeTimeConstant")
        attitudeTimeConstant = (float)Convert.ToDouble(v);
      else if (k=="controllerLogFilename")
        _pid3d.SetLogFilename(v);
      else
        return false;
      _pid3d.Reset();
      return true;
    }

    // returns the thrust cone limited vector
    // accepted original thrust vector as F, and returns unlimited F after adjustment
    Vector3d GetThrustVector(Vector3d tr, Vector3d tv, Vector3d dr, Vector3d dv,
                             float t,
                             Vector3d F, out Vector3d unlimF)
    {
      Vector3d F2 = _pid3d.Update(tr,tv,dr,dv,t);
      unlimF = F + F2;
      // Reduce sideways components
      return ConeUtils.ClosestThrustInsideCone(maxThrustAngle,(float)amin,(float)amax,unlimF);
    }

    // F should be initial acceleration before any corrections
    // For forced landing this is the opposite vector to gravity (to hover)
    // For an autopilot path this is the acceleration vector from the computed solution
    public void AutopilotStepToTarget(
            Vector3d tatt, Vector3d tr, Vector3d tv,
            Vector3d dr, Vector3d dv, Vector3d da,
            Vector3d g, float t,
            out double throttle, out Vector3d thrustV, out float att_err, out bool shutdownEnginesNow)
    {
      float deltaTime = t - last_t;
      last_t = t;
      shutdownEnginesNow = false;
      // In Realism Overhaul where engines have limited throttle keep throttle >= 0.01 once
      // engines have been ignited and engines have non 0% minThrust (suggests Realism Overhaul)
      float minThrottle = (ignited && (amin>0)) ? 0.01f : 0;
      throttle = minThrottle;

      Vector3d da2 = Vector3d.zero; // actual thrust vector
      att_err = 0;
      if (ignited)
      {
        Vector3d unlimda;
        da2 = GetThrustVector(tr,tv,dr,dv,t,da, out unlimda);
        throttle = Mathf.Clamp((float)(da2.magnitude-amin)/(amax-amin+0.01f),minThrottle,1);
        // Decide to shutdown engines for final touch down? (within 3 secs)
        // Criteria should be if amin maintained with current engines we could not reach next target
        // height
        double minHeight  = MinHeightAtMinThrust(tr.y,tv.y,g.magnitude);
        // Criteria for shutting down engines
        // - we could not reach ground at minimum thrust (would ascend)
        // - falling less than 20m/s (otherwise can decide to shutdown engines when still high and travelling fast)
        bool cant_reach_ground = (minHeight > 0) && (tv.y > -10);

        // Other criteria?
        // - amin > required acceleration for rest of trajectory?
        if (cant_reach_ground)
          shutdownEnginesNow = true;
        // Shutoff throttle if pointing in wrong direction
        float ddot = (float)Vector3d.Dot(Vector3d.Normalize(tatt),Vector3d.Normalize(da2));
        att_err = Mathf.Acos(ddot)*180/Mathf.PI;
        if ((att_err >= idleAngle) && (da2.magnitude>0.01))
          throttle = 0.01f; // some throttle to steer? (if no RCS and main thruster gimbals)
      }
      else
      {
        // Engines not ignited
        // Only ignite engines if solution has non-zero thrust
        if (da.magnitude > 0.1f)
        {
          throttle = 0.01f; // will get corrected on next iteration
          ignited = true;
        }
      }

      // If no clear steer direction point upwards
      if (da2.magnitude < 0.01f)
        da2 = new Vector3d(0,0.01f,0);

      // Update current thrust with smoothing
      // decompose into throttle setting and direction
      double alphaThrottle = 1, alphaAttitude = 1;
      if (throttleTimeConstant > 0)
        alphaThrottle = 1 - Math.Exp(-deltaTime/throttleTimeConstant);
      if (attitudeTimeConstant > 0)
        alphaAttitude = 1 - Math.Exp(-deltaTime/attitudeTimeConstant);
      currentThrustMag = currentThrustMag*(1-alphaThrottle) + da2.magnitude*alphaThrottle;
      currentThrustDir = currentThrustDir*(1-alphaAttitude) + Vector3d.Normalize(da2)*alphaAttitude;
      thrustV = currentThrustMag * currentThrustDir;
    }


    public void GetControlOutputs(Trajectory traj,
                    Vector3d r,
                    Vector3d v,
                    Vector3d att,
                    Vector3d g,
                    float t,
                    out Vector3d dr,
                    out Vector3d dv,
                    out Vector3d da,
                    out double throttle, out Vector3d thrustV, out float att_err,
                    out bool shutdownEnginesNow)
    {
      double desired_t; // closest time in trajectory (desired)
      shutdownEnginesNow = false;
      traj.FindClosest(r, v, out dr, out dv, out da, out desired_t, 0.5f, 0.5f);
      double distToLast = (traj.r[traj.Length()-1] - r).magnitude;
      //System.Console.Error.WriteLine("r="+r+" distToLast="+distToLast);
      Vector3d touchdownv = new Vector3d(0,-touchdownSpeed,0);
      if (distToLast < finalDescentAlt) // If near ground just descent vertically
      {
        // TODO: Assumes touching down in ground downwards
        double h = r.y - traj.r[traj.Length()-1].y;
        // p=0 when h==distToLast
        // p=1 when h<=0.5*distToLast
        float p = Mathf.Clamp((float)(h/finalDescentAlt),0,1);
        // Smoothed from trajectory to touchdown from 0.5*finalDescentAlt to finalDescentAlt
        p = Mathf.Clamp(HGUtils.LinearMap(p,0.5f,1,0,1),0,1);
        // Interpolate between computed trajectory and touchdown speed for smooth transition
        r = dr*p + r*(1-p);
        dv = dv*p + touchdownv*(1-p);
        da = da*p + (-g)*(1-p);
        finalDescent = true;
      }
      AutopilotStepToTarget(att, r, v, dr, dv, da, g, t, out throttle, out thrustV, out att_err, out shutdownEnginesNow);
    }

    public bool IsFinalDescent()
    {
      return finalDescent;
    }
  }
}
