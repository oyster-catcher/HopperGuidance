using KSPAssets;
using UnityEngine;

namespace HopperGuidance
{
  public class Trajectory
  {
    double dt = 0;
    public Vector3d [] r;
    public Vector3d [] v;
    public Vector3d [] a;

    // Take a set of thrust vectors at positions 0, dt, 2*dt, etc...
    // and return a higher fidelity track of positions and thrust directions: r, rrr
    // at intervals of simdt
    public void Simulate(double T, Vector3d[] thrusts, Vector3d r0, Vector3d v0, Vector3d g, double a_dt, double extendTime)
    {
      // Simulate
      dt = a_dt;
      Vector3d cr = r0;
      Vector3d cv = v0;
      int j = 0;
      int N = thrusts.Length;
      int M = (int)(T/dt+1);
      int extendM = (int)((T+extendTime)/dt+1);
      double t=0;
      r = new Vector3d[extendM];
      v = new Vector3d[extendM];
      a = new Vector3d[extendM];
      while(j < M)
      {
        r[j] = cr;
        v[j] = cv;
        Vector3d ca = Vector3d.zero;
        double [] w = Solve.BasisWeights(t,T,N);
        for( int i = 0 ; i < N ; i++ )
        {
          ca = ca + w[i] * thrusts[i];
        }
        a[j] = ca;
        cr += cv*dt + 0.5*ca*dt*dt + 0.5*g*dt*dt;
        cv += ca*dt + g*dt;
        t += dt;
        j++;
      }
      // Continue with same position and velocity and acceleration
      // to equal gravity
      while(j < extendM)
      {
        r[j] = cr;
        v[j] = cv;
        a[j] = -g*dt;
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
      double atp_dot_atb = Vector3d.Dot(a_to_p,a_to_p);
      t = atp_dot_atb / atb2; // The normalized "distance" from a to your closest point
      if (t<0) {t=0;}
      if (t>1) {t=1;}
      return (A + a_to_b*t);
    }

    // Correct final position in trajectory to cover up numeric errors
    public void CorrectFinal(Vector3d rf, Vector3d vf)
    {
      Vector3d r_err = r[Length()-1] - rf;
      Vector3d v_err = v[Length()-1] - vf;
      for( int i = 0 ; i < Length() ; i++ )
      {
        r[i] = r[i] - r_err*((double)i/Length());
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
      return r.Length;
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
                            double rWeight = 0.3, double vWeight = 0.7)
    {
      // Find closest point in speed and position and returns the index
      int ci = -1;
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
      Debug.Log("Found "+ci+" in "+r.Length);
      if (ci!=-1)
      {
        closest_r = r[ci];
        closest_v = v[ci];
        closest_a = a[ci];
        closest_t = ci*dt;
        // Now check line-segment either side
        if (ci-1 >= 0)
        {
          Vector3d r1 = FindClosestPointOnLine(r[ci-1],r[ci],a_r,out double t);
          double dr = (a_r-r1).magnitude;
          if (dr < cdist)
          {
            closest_r = r1;
            closest_v = v[ci-1]*(1-t) + v[ci]*t;
            closest_a = a[ci-1]*(1-t) + a[ci]*t;
            cdist = dr;
          }
        }

        // Now check line-segment either side
        if (ci+1 < r.Length)
        {
          Vector3d r2 = FindClosestPointOnLine(r[ci],r[ci+1],a_r,out double t);
          double dr = (a_r-r2).magnitude;
          if (dr < cdist)
          {
            closest_r = r2;
            closest_v = v[ci]*(1-t) + v[ci+1]*t;
            closest_a = a[ci]*(1-t) + a[ci+1]*t;
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

    public void WriteLog(string filename, Transform transform)
    {
      System.IO.StreamWriter f = new System.IO.StreamWriter(filename);
      double t = 0;
      Vector3d tr = Vector3.zero;
      Vector3d tv,ta;
      f.WriteLine("time x y z vx vy vz ax ay az");
      for(int i = 0 ; i < Length() ; i++)
      {
        tr = transform.InverseTransformPoint(r[i]);
        tv = transform.InverseTransformVector(v[i]);
        ta = transform.InverseTransformVector(a[i]);
        f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1}",t,tr.x,tr.y,tr.z,tv.x,tv.y,tv.z,ta.x,ta.y,ta.z));
        t += dt;
      }
      f.Close();
      //Debug.Log("Final target: "+tr);
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
