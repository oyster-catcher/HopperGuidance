using System;
using System.Collections.Generic;

using UnityEngine;
using UnityEngine.UI;

namespace HopperGuidance
{
  public class Controller
  {
    // Parameters
    public PID3d _pid3d = new PID3d();
    public float minThrust = 1;
    public float maxThrustAngle = 25;
    public float idleAngle = 90;
    public float amin = 0;
    public float amax = 30;
    public float touchdownSpeed = 1;
    public float finalDescentHeight = 0; // Ignore trajectory under this height and descent to ground
    public double throttleTimeConstant = 0; // Instanteous

    // Live values
    public bool ignited = false;
    public Vector3d currentThrustV = Vector3d.zero;


    public Controller(float kp1=1, float ki1=0, float kd1=0, float kp2=1, float ki2=0, float kd2=0, float amax=1000, float ymult=2)
    {
      _pid3d.kp1 = kp1;
      _pid3d.ki1 = ki1;
      _pid3d.kd1 = kd1;
      _pid3d.kp2 = kp2;
      _pid3d.ki2 = ki2;
      _pid3d.kd2 = kd2;
      _pid3d.amax = amax;
      _pid3d.ymult = ymult;
      _pid3d.Reset();
    }

    double MinHeightAtMinThrust(double y, double vy,double amin,double g)
    {
      double minHeight = 0;
      if (amin < g)
        return -float.MaxValue;
      double tHover = -vy/amin; // time to come to hover
      minHeight = y + vy*tHover + 0.5*amin*tHover*tHover - 0.5*g*tHover*tHover;
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
      else if (k=="ymult")
        _pid3d.ymult = (float)Convert.ToDouble(v);
      else if (k=="finalDescentHeight")
        finalDescentHeight = (float)Convert.ToDouble(v);
      else if (k=="touchdownSpeed")
        touchdownSpeed = (float)Convert.ToDouble(v);
      else if (k=="throttleTimeConstant")
        throttleTimeConstant = (float)Convert.ToDouble(v);
      else
        return false;
      _pid3d.Reset();
      return true;
    }

    // returns the thrust cone limited vector
    // accepted original thrust vector as F, and returns unlimited F after adjustment
    Vector3d GetThrustVector(Vector3d tr, Vector3d tv, Vector3d dr, Vector3d dv,
                             float maxThrustAngle,
                             float deltaTime,
                             Vector3d F, out Vector3d unlimF)
    {
      //System.Console.WriteLine("r="+tr+" deltaTime="+deltaTime);
      Vector3d F2 = _pid3d.Update(tr,tv,dr,dv,deltaTime);
      unlimF = F + F2;
      // Reduce sideways components
      return ConeUtils.ClosestThrustInsideCone((float)maxThrustAngle,(float)amin,(float)amax,unlimF);
    }

    // F should be initial acceleration before any corrections
    // For forced landing this is the opposite vector to gravity (to hover)
    // For an autopilot path this is the acceleration vector from the computed solution
    public void AutopilotStepToTarget(
            Vector3d tatt, Vector3d tr, Vector3d tv,
            Vector3d dr, Vector3d dv, Vector3d da,
            Vector3d g, float deltaTime,
            out double throttle, out Vector3d thrustV)
    {
      //ComputeMinMaxThrust(out _minThrust,out _maxThrust);
      //float amax = (float)(maxThrust/mass);
      //float amin = (float)(minThrust/mass);
      // In Realism Overhaul where engines have limited throttle keep throttle >= 0.01 once
      // engines have been ignited and engines have non 0% minThrust (suggests Realism Overhaul)
      float minThrottle = (ignited && (minThrust>0)) ? 0.01f : 0;
      throttle = minThrottle;

      Vector3d da2 = Vector3d.zero; // actual thrust vector
      float att_err = 0;
      if (ignited)
      {
        Vector3d unlimda;
        da2 = GetThrustVector(tr,tv,dr,dv,maxThrustAngle,deltaTime,da, out unlimda);
        throttle = Mathf.Clamp((float)(da2.magnitude-amin)/(amax-amin+0.01f),minThrottle,1);
        // Decide to shutdown engines for final touch down? (within 3 secs)
        // Criteria should be if amin maintained with current engines we could not reach next target
        // height
        double minHeight  = MinHeightAtMinThrust(tr.y,tv.y,amin,g.magnitude);
        // Criteria for shutting down engines
        // - we could not reach ground at minimum thrust (would ascend)
        // - engines not already shutdown
        // - needed acceleration less than the minimum
        // - falling less than 20m/s (otherwise can decide to shutdown engines when still high and travelling fast)
        bool cant_reach_ground = (minHeight > 0) && (tv.y > -20);
        //Log("minHeight="+minHeight+" tv.y="+tv.y+" unlimda.y="+unlimda.y);
        // Other criteria?
        // - amin > required acceleration for rest of trajectory?
        if (cant_reach_ground)
        {
          // Check if thrust beyond this point is higher again?
          // TODO: Add this back in
          //int numShutdown = ShutdownOuterEngines(unlimda.y*vessel.totalMass); // set thrust in range for desired accel
          //if (numShutdown > 0)
          //Log("Shutdown "+numShutdown+" engines for controlled touchdown", true);
        }
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

      // TODO: Open log files
      //if ((_vesselWriter == null) && (_logging))
      //{
      //  _vesselWriter = new System.IO.StreamWriter(_vesselLogFilename);
      //  _vesselWriter.WriteLine("time x y z vx vy vz ax ay az att_err amin amax");
      //  _tgtWriter = new System.IO.StreamWriter(_tgtLogFilename);
      //  _tgtWriter.WriteLine("time x y z vx vy vz ax ay az att_err amin amax");
      //}

      //double tRel = Time.time - _startTime;
      //if ((_logging) && (tRel >= last_t+log_interval))
      //{
      //  LogData(_tgtWriter, tRel, dr, dv, da, 0, solver.amin, solver.amax);
      //  double thrustProp = 0;
      //  if ((_maxThrust > 0) && (da2.magnitude > 0))
      //    thrustProp = (amax*GetCurrentThrust())/(da2.magnitude*_maxThrust);
      //  LogData(_vesselWriter, tRel, tr, tv, da2*thrustProp, att_err, amin, amax);
      //  last_t = tRel;
      //}

      //DrawAlign(tr,dr,_transform,aligncol);
      // If no clear steer direction point upwards
      if (da2.magnitude < 0.01f)
        da2 = new Vector3d(0,0.01f,0);
      //DrawSteer(tr, tr+3*vessel.vesselSize.x*Vector3d.Normalize(da2), _transform, thrustcol);
      //da2 = _transform.TransformVector(da2); // transform back to world co-ordinates
      //vessel.Autopilot.SAS.lockedMode = false;
      //vessel.Autopilot.SAS.SetTargetOrientation(da2,false);
      //Log("da2="+da2+" throttle="+throttle+" minThrottle="+minThrottle+" amax="+amax+" ignited="+ignited);
      //
      // Update current thrust with smoothing
      double alpha = 1;
      if (throttleTimeConstant > 0)
        alpha = 1 - Math.Exp(-deltaTime/throttleTimeConstant);
      currentThrustV = currentThrustV*(1-alpha) + da2*alpha;
      thrustV = currentThrustV;
    }


    public void GetControlOutputs(Trajectory traj,
                    Vector3d r,
                    Vector3d v,
                    Vector3d g,
                    float deltaTime,
                    out Vector3d dr,
                    out Vector3d dv,
                    out Vector3d da,
                    out double throttle, out Vector3d thrustV)
    {
        Vector3d att = new Vector3d(0,1,0);
        double desired_t; // closest time in trajectory (desired)
        if (r.y < finalDescentHeight) // If near ground just descent vertically
        {
          System.Console.Error.WriteLine("FINAL DESCENT!!!");
          desired_t = 0;
          dr = r;
          dv = v;
          da = -g - new Vector3d(v.x*0.4,v.y*0.4,v.z*0.4);
        }
        else
          traj.FindClosest(r, v, out dr, out dv, out da, out desired_t, 0.5f, 0.5f);
        AutopilotStepToTarget(att, r, v, dr, dv, da, g, deltaTime, out throttle, out thrustV);
    }
  }
}
