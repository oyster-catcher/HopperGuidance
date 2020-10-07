using KSPAssets;
using UnityEngine;
using System;
using System.Collections.Generic;
using HopperGuidance;
using Trajectories;

namespace HopperGuidance
{
  public class Trajectory
  {
    public double dt = 0;
    public Vector3d [] r = null;
    public Vector3d [] v = null;
    public Vector3d [] a = null;
    public double T = 0;

    // Take a set of thrust vectors at positions 0, dt, 2*dt, etc...
    // and return a higher fidelity track of positions and thrust directions: r, rrr
    // at intervals of simdt
    public void Simulate(double a_T, ThrustVectorTime [] thrusts, Vector3d r0, Vector3d v0, Vector3d g, double a_dt)
    {
      if (thrusts == null)
        return;
      T = a_T;
      List<float> thrust_times = new List<float>();
      for(int i=0; i<thrusts.Length; i++)
        thrust_times.Add(thrusts[i].t);

      // Simulate
      dt = a_dt;
      Vector3d cr = r0;
      Vector3d cv = v0;
      int start = 0;
      if (r != null)
        start = r.Length;
      int j = start;
      int N = thrusts.Length;
      int M = (int)(T/dt+1);
      double t=0;
      if (r == null)
      {
        r = new Vector3d[M];
        v = new Vector3d[M];
        a = new Vector3d[M];
      }
      else
      {
        Array.Resize(ref r, start+M);
        Array.Resize(ref v, start+M);
        Array.Resize(ref a, start+M);
      }
      while(j < start+M)
      {
        Vector3d ca = Vector3d.zero;
        cr = r0;
        cv = v0;
        double [] wr;
        double [] wv;
        // TODO: This is slow and expensive in RVWeightsToTime but accurate and
        // avoid accumulation of error when computing forwards only
        Solve.RVWeightsToTime(t, dt, thrust_times, out wr, out wv);
        double [] w = Solve.BasisWeights(t, thrust_times);
        for(int i = 0 ; i < N ; i++)
        {
          ca += w[i]  * thrusts[i].v;
          cr += wr[i] * thrusts[i].v;
          cv += wv[i] * thrusts[i].v;
        }
        cr = cr + v0*t + 0.5*t*t*g;
        cv = cv + g*t;
        a[j] = ca;
        r[j] = cr;
        v[j] = cv;
        j++;
        t += dt;
      }
    }

    // Simulates an aerodynamic descent trajectory at current attitude with no thrust
    // applied. The reference from in the vessels main body transform
    public void SimulateAero(Vessel vessel, float tgtAlt)
    {
      Vector3d wr = vessel.GetWorldPos3D();
      CelestialBody b = vessel.mainBody;
      Vector3d g = FlightGlobals.getGeeForceAtPosition(wr);
      Vector3d wv = vessel.GetObtVelocity();
      T = 0;
      dt = 1;
      int newsize=50;
      int i = 0;

      Array.Resize<Vector3d>(ref r, newsize);
      Array.Resize<Vector3d>(ref v, newsize);
      Array.Resize<Vector3d>(ref a, newsize); 
      double altitude = (wr - b.transform.position).magnitude - b.Radius;
      while(altitude > tgtAlt)
      {
        if (i > newsize-1)
        {
          newsize = newsize + 50;
          Array.Resize<Vector3d>(ref r, newsize);
          Array.Resize<Vector3d>(ref v, newsize);
          Array.Resize<Vector3d>(ref a, newsize); 
          g = FlightGlobals.getGeeForceAtPosition(wr);
        }
        r[i] = b.transform.InverseTransformPoint(wr);
        v[i] = b.transform.InverseTransformVector(wv);
        //Debug.Log("wr="+wr+" wv="+wv+" alt="+altitude);
        Vector3d f = Trajectories.StockAeroUtil.SimAeroForce(vessel, wv, wr);
        altitude = (wr - b.position).magnitude - b.Radius; // NOTE: latitude of last point
        wr = wr + wv*dt;
        wv = wv - (f/vessel.totalMass)*dt;
        wv = wv + g*dt;
        T = T + dt;
        i++;
      }
      Array.Resize<Vector3d>(ref r, i);
      Array.Resize<Vector3d>(ref v, i);
      Array.Resize<Vector3d>(ref a, i);
      // Reset last point by interpolating to altitude=0
      double alt1 = r[i-2].magnitude - b.Radius;
      double alt2 = r[i-1].magnitude - b.Radius;
      double p = (alt1 - tgtAlt)/(alt1 - alt2);
      r[i-1] = r[i-2] + (0.5*v[i-2] + 0.5*v[i-1]) * p;
      double altf = r[i-1].magnitude - b.Radius;
      Debug.Log("alt1="+alt1+" alt2="+alt2+" altf="+altf+" p="+p);
      Debug.Log("SimulateAero() finished. Trajectory size="+i+" alt="+altitude);
    }

    public void CompensateForBodyRotation(CelestialBody body)
    {
      double t=0;
      for(int i=0; i<r.Length; i++)
      {
        float ang = (float)((-t) * body.angularVelocity.magnitude / Math.PI * 180.0);
        Quaternion bodyRotation = Quaternion.AngleAxis(ang, body.angularVelocity.normalized);
        r[i] = bodyRotation * r[i];
        t = t + dt;
      }
    }

    // Extend the trajectory at the same velocity for a period of time
    // compute the thrust vector to cancel gravity
    public void Extend(Vector3d a_v, Vector3d a_g, double a_extendTime)
    {
      int oldsize = r.Length;
      int newsize = r.Length + (int)(a_extendTime/dt);
      Array.Resize<Vector3d>(ref r, newsize);
      Array.Resize<Vector3d>(ref v, newsize);
      Array.Resize<Vector3d>(ref a, newsize); 
      Vector3d cr = r[oldsize-1];
      System.Console.Error.WriteLine("g="+a_g);
      for( int i=oldsize; i<newsize; i++ )
      {
        r[i] = cr;
        v[i] = a_v;
        a[i] = -a_g;
        cr = cr + a_v*dt;
      }
    }

    // A,B describes a line from A to B and P is the point
    // return the closest position and the parameter to from 0 to 1 (A to B)
    public static Vector3d FindClosestPointOnLine(Vector3d A, Vector3d B, Vector3d P, out double t)
    {
      Vector3d a_to_p = P - A;
      Vector3d a_to_b = B - A;
      double atb2 = a_to_b.x*a_to_b.x + a_to_b.y*a_to_b.y + a_to_b.z*a_to_b.z;
      double atp_dot_atb = Vector3d.Dot(a_to_b,a_to_p);
      t = atp_dot_atb / atb2; // The normalized "distance" from a to your closest point
      if (t<0) {t=0;}
      if (t>1) {t=1;}
      return (A + a_to_b*t);
    }

    // rA,rB describes a line from rA to rB and R is the position
    // vA,vB describes a line from vA to vB and V is the velocity
    // return the closest position and the parameter to from 0 to 1 (A to B)
    public static double FindClosestPointOnLineRV(Vector3d rA,Vector3d rB, Vector3d vA, Vector3d vB,
                                                  Vector3d r,Vector3d v,
                                                  out Vector3d c_r, out Vector3d c_v, out double t,
                                                  double rWeight=1, double vWeight=1)
    {
      //Vector3d a_to_p = P - A;
      //Vector3d a_to_b = B - A;
      Vector3d srA = rA * rWeight;
      Vector3d srB = rB * rWeight;
      Vector3d svA = vA * vWeight;
      Vector3d svB = vB * vWeight;
      Vector3d sr = r * rWeight;
      Vector3d sv = v * vWeight;
      double [] a_to_p = new double [] {sr.x-srA.x,  sr.y-srA.y,  sr.z-srA.z,  sv.x-svA.x, sv.y-svA.y,  sv.z-svA.z};
      double [] a_to_b = new double [] {srB.x-srA.x, srB.y-srA.y, srB.z-srA.z, svB.x-svA.x, svB.y-svA.y, svB.z-svA.z};
      double atb2 = a_to_b[0]*a_to_b[0] + a_to_b[1]*a_to_b[1] + a_to_b[2]*a_to_b[2] +
                    a_to_b[3]*a_to_b[3] + a_to_b[4]*a_to_b[4] + a_to_b[5]*a_to_b[5];
      //double atp_dot_atb = Vector3d.Dot(a_to_b,a_to_p);
      double atp_dot_atb = a_to_b[0]*a_to_p[0] + a_to_b[1]*a_to_p[1] + a_to_b[2]*a_to_p[2] +
                           a_to_b[3]*a_to_p[3] + a_to_b[4]*a_to_p[4] + a_to_b[5]*a_to_p[5];
      t = atp_dot_atb / atb2; // The normalized "distance" from a to your closest point
      if (t<0) {t=0;}
      if (t>1) {t=1;}
      c_r = rA + (rB - rA) * t;
      c_v = vA + (vB - vA) * t;
      return rWeight*(c_r-r).magnitude + vWeight*(c_v-v).magnitude;
    }

    // Correct to go exactly through the targets
    public void CorrectForTargets(List<SolveTarget> targets)
    {
      int start=0, j=0;
      Vector3d start_r_err = Vector3d.zero;
      Vector3d start_v_err = Vector3d.zero;
      for( int i = 0 ; i < targets.Count ; i++ )
      {
        int end = (int)((targets[i].t/T) * Length() + 0.5f);
        end = Math.Min(Length()-1,end); // inclusive
        Vector3d end_r_err = r[end] - targets[i].r;
        Vector3d end_v_err = Vector3d.zero;
        if (targets[i].vaxes != 0)
        {
          end_v_err = v[end] - targets[i].v;
        }
        for ( j=start; j<=end ; j++)
        {
          double p = (j-start)/(double)(end-start);
          r[j] = r[j]  - start_r_err * (1-p) - end_r_err * p;
          v[j] = v[j]  - start_v_err * (1-p) - end_v_err * p;
          //System.Console.Error.WriteLine("Correcting "+j+" by "+(start_r_err * (1-p))+" and "+(end_r_err * p)+" to "+r[j]+" p="+p);
        }
        start = j;
        start_r_err = end_r_err;
        start_v_err = end_v_err;
      }
    }

    public void Init(double a_dt, Vector3d [] ar, Vector3d [] av, Vector3d [] a_a)
    {
      dt = a_dt;
      r = ar;
      v = av;
      a = a_a;
    }

    public int Length()
    {
      if (r != null)
        return r.Length;
      else
        return 0;
    }

    // Weight distance on the position and velocity vector from 1 to 2
    public double WeightedDist(Vector3d r1,Vector3d v1,Vector3d r2,Vector3d v2,double rWeight,double vWeight)
    {
        double dr = (r2-r1).magnitude;
        double dv = (v2-v1).magnitude;
        return dr*rWeight + dv*vWeight;
    }

    public double FindClosest(Vector3d a_r, Vector3d a_v,
                            out Vector3d closest_r, out Vector3d closest_v, out Vector3d closest_a, out double closest_t,
                            double rWeight = 1.0, double vWeight = 1.0)
    {
      // Find closest point in speed and position and returns the index
      int ci = -1;
      double t;
      double cdist = 99999999;
      for(int i=0;i<r.Length; i++)
      {
        double dr = (a_r-r[i]).magnitude;
        double dv = (a_v-v[i]).magnitude;
        double d = dr*rWeight + dv*vWeight;
        if (d < cdist)
        {
          ci = i;
          cdist = d;
        }
      }
      if (ci!=-1)
      {
        //Debug.Log("Found ci="+ci+" out of "+r.Length);
        closest_r = r[ci];
        closest_v = v[ci];
        closest_a = a[ci];
        closest_t = ci*dt;
        // Now check line-segment either side
        if (ci-1 >= 0)
        {
          double d1 = FindClosestPointOnLineRV(r[ci-1],r[ci],v[ci-1],v[ci],a_r,a_v,
                                               out Vector3d c_r,out Vector3d c_v,out t,rWeight,vWeight);
          if (d1 < cdist)
          {
            closest_r = c_r;
            closest_v = c_v;
            closest_a = a[ci-1] + (a[ci] - a[ci-1])*t;
            cdist = d1;
          }
        }
        if (ci+1 < r.Length)
        {
          double d2 = FindClosestPointOnLineRV(r[ci],r[ci+1],v[ci],v[ci+1],a_r,a_v,
                                               out Vector3d c_r,out Vector3d c_v,out t,rWeight,vWeight);
          if (d2 < cdist)
          {
            closest_r = c_r;
            closest_v = c_v;
            closest_a = a[ci] + (a[ci+1] - a[ci])*t;
            cdist = d2;
          }
        }
        return (closest_r - a_r).magnitude;
      }
      else
      {
        closest_r = Vector3d.zero;
        closest_v = Vector3d.zero;
        closest_a = Vector3d.zero;
        closest_t = 0;
        return (closest_r - a_r).magnitude;
      }
    }

    public void Write(string filename = null, List<string> comments = null)
    {
      System.IO.StreamWriter f;
      if (filename != null)
        f = new System.IO.StreamWriter(filename);
      else
        f = new System.IO.StreamWriter(System.Console.OpenStandardOutput());
      double t = 0;
      Vector3d tr = Vector3d.zero;
      Vector3d tv,ta;
      if (comments != null)
      {
        foreach( string comment in comments)
          f.WriteLine("# "+comment);
      }
      f.WriteLine("time x y z vx vy vz ax ay az att_err");
      for(int i = 0 ; i < Length() ; i++)
      {
        tr = r[i];
        tv = v[i];
        ta = a[i];
        f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1} 0",t,tr.x,tr.y,tr.z,tv.x,tv.y,tv.z,ta.x,ta.y,ta.z));
        t += dt;
      }
      f.Close();
    }

/*
    static int Main(string[] argv)
    {
      Trajectory traj = new Trajectory();

      if (argv.Length != 8) {
        System.Console.Error.WriteLine("usage: Trajectory.exe <weight_r> <weight_v> <x> <y> <z> <vx> <vy> <vz>");
        return(1);
      }
      double wr = Convert.ToDouble(argv[0]);
      double wv = Convert.ToDouble(argv[1]);
      double x = Convert.ToDouble(argv[2]);
      double y = Convert.ToDouble(argv[3]);
      double z = Convert.ToDouble(argv[4]);
      double vx = Convert.ToDouble(argv[5]);
      double vy = Convert.ToDouble(argv[6]);
      double vz = Convert.ToDouble(argv[7]);

      Vector3d [] trackr = {new Vector3d(0,2,0),new Vector3d(0,1,0), new Vector3d(0,0,0), new Vector3d(0,-1,0)};
      Vector3d [] trackv = {new Vector3d(0,-5,0),new Vector3d(0,-4,0), new Vector3d(0,-3,0), new Vector3d(0,-3,0)};
      Vector3d [] tracka = {new Vector3d(0,13,0),new Vector3d(0,13,0), new Vector3d(0,9,0), new Vector3d(0,9,0)};

      traj.Init(0.1,trackr,trackv,tracka);

      Vector3d r = new Vector3d(x,y,z);
      Vector3d v = new Vector3d(vx,vy,vz);
      Vector3d close_r,close_v,close_a;
      double close_t;
      System.Console.WriteLine("Closest:"+traj.FindClosest(r,v,out close_r,out close_v,out close_a,out close_t,wr,wv));
      System.Console.WriteLine("close_r="+close_r+" close_v="+close_v+" close_a="+close_a);

      return 0;
    }
*/
  }
}
