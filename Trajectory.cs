using KSPAssets;
using UnityEngine;
using System;
using System.Collections.Generic;
using HopperGuidance;

namespace HopperGuidance
{
  public class Trajectory
  {
    public double dt = 0;
    public Vector3d [] r = null;
    public Vector3d [] v = null;
    public Vector3d [] a = null;

    // Take a set of thrust vectors at positions 0, dt, 2*dt, etc...
    // and return a higher fidelity track of positions and thrust directions: r, rrr
    // at intervals of simdt
    public void Simulate(double T, ThrustVectorTime [] thrusts, Vector3d r0, Vector3d v0, Vector3d g, double a_dt, double extendTime)
    {
      if (thrusts == null)
        return;
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
      int extendM = (int)((T+extendTime)/dt+1);
      double t=0;
      if (r == null)
      {
        r = new Vector3d[extendM];
        v = new Vector3d[extendM];
        a = new Vector3d[extendM];
      }
      else
      {
        Array.Resize(ref r, start+extendM);
        Array.Resize(ref v, start+extendM);
        Array.Resize(ref a, start+extendM);
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
      // Continue with same position and velocity and acceleration
      // to equal gravity
      while(j < start+extendM)
      {
        r[j] = cr;
        v[j] = cv;
        a[j] = -g;
        cr += cv*dt;
        t += dt;
        j++;
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

    // Correct final position in trajectory to cover up numeric errors

    // Correct final position in trajectory to cover up numeric errors
    public void CorrectFinal(Vector3d rf, Vector3d vf, bool correct_r, bool correct_v)
    {
      Vector3d r_err = r[Length()-1] - rf;
      Vector3d v_err = v[Length()-1] - vf;
      for( int i = 0 ; i < Length() ; i++ )
      {
        if (correct_r)
          r[i] = r[i] - r_err*((double)i/Length());
        if (correct_v)
          v[i] = v[i] - v_err*((double)i/Length());
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
      double rdist = 99999999;
      for(int i=0;i<r.Length; i++)
      {
        double dr = (a_r-r[i]).magnitude;
        double dv = (a_v-v[i]).magnitude;
        double d = dr*rWeight + dv*vWeight;
        if (d < cdist)
        {
          ci = i;
          cdist = d;
          rdist = dr;
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
        return rdist;
      }
      else
      {
        closest_r = Vector3d.zero;
        closest_v = Vector3d.zero;
        closest_a = Vector3d.zero;
        closest_t = 0;
        return rdist;
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

//    static int Main(string[] argv)
//    {
//      Trajectory traj = new Trajectory();
//
//      Vector3d [] trackr = {new Vector3d(0,0,0),new Vector3d(3,10,4), new Vector3d(8,12,8)};
//
//      traj.Init(0.1,trackr,trackr,trackr);
//
//      Vector3d r0 = new Vector3d(8,8,8);
//      Vector3d closer,closev,close_a;
//      System.Console.WriteLine("Closest:"+traj.FindClosest(r0,r0,out closer,out closev,out close_a));
//      System.Console.WriteLine("closer="+closer);
//
//      return 0;
//    }
  }
}
