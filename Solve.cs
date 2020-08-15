#define MINDESCENTANGLE
#define MAXTHRUSTANGLE
#define MAXVELOCITY
#define MINTHRUST
//#define DUMP
#define DEBUG

using System;
using System.IO;
using System.Collections.Generic;
using KSPAssets;
using UnityEngine;

namespace HopperGuidance
{
  public enum SolveTargetType
  {
    Position,
    Velocity,
    Both
  }
  public class SolveTarget
  {
    public Vector3d r;
    public Vector3d v;
    public float t;
    public SolveTargetType type;
  }

  public class Solve
  {
    // Parameters to control solution
    public double Tmin = 0;
    public double Tmax = 300;
    public double minDurationPerThrust = 4; // Insert extra thrust vector between targets
    public int maxThrustsBetweenTargets = 2;
    public double checkGapFirst = 1; // duration between checks (none at T=0)
    public double checkGapMult = 1.5; // increase time from ends by this proportion
    public double g = 9.8;
    public double amin = 0;
    public double amax = 30;
    public double vmax = 1000;
    public double minDescentAngle = 10;
    public double tol = 0.5;
    public double maxThrustAngle = 180;
    public double maxLandingThrustAngle = 20;
    public int fidelity = 20;
    public double timePenalty = 0; // If zero minimize fuel, as 1 then 1 extra second cost 1 fuel
    public float extraTime = 0.5f; // Add this on to minimum fuel solution (helps larger craft?)
    public Vector3d apex; // Apex position of cone of minimum descent

    // Last stored inputs to GFold() - in natural space for solution
    // with Y as the up direction
    public Vector3d r0;
    public Vector3d v0;
    public List<SolveTarget> targets;

    // Outputs of solution (stored and returned)
    public int N; // choosen N given duration of flight
    public double T;
    public double dt;
    public Vector3d [] thrusts;
    int retval;

#if (DUMP)
    public void WriteMatrix(string name,double [,] a,int rows,int cols)
    {
      string tab="";
      System.Console.Write("{0}=",name);
      for(int i=0; i<rows; i++) {
        System.Console.Write("{0}[",tab);
        tab="  ";
        for(int j=0; j<cols; j++) {
          System.Console.Write("{0:F2} ",a[i,j]);
        }
        System.Console.WriteLine("]");
      }
    }

    void WriteVector(string name,double [] a,int size)
    {
      System.Console.Write("{0}=[",name);
      for(int i=0; i<size; i++) {
        System.Console.Write("{0:F2} ",a[i]);
      }
      System.Console.WriteLine("]");
    }

    void WriteVector(string name,int [] a,int size)
    {
      System.Console.Write("{0}=[",name);
      for(int i=0; i<size; i++) {
        System.Console.Write("{0:F2} ",a[i]);
      }
      System.Console.WriteLine("]");
    }
#endif
    public string Vec2Str(Vector3d v)
    {
      return string.Format("[{0:F2},{1:F2},{2:F2}]",v.x,v.y,v.z);
    }

    public string DumpString()
    {
      string msg = ((retval>=1)&&(retval<=5))?"SUCCEED":"FAIL";
      string stargets = "";
      foreach( SolveTarget tgt in targets )
      {
        stargets = stargets + "target={";
        if(( tgt.type == SolveTargetType.Position ) || (tgt.type == SolveTargetType.Both))
        {
          stargets = stargets + String.Format("r="+Vec2Str(tgt.r)) + ":";
        }
        if(( tgt.type == SolveTargetType.Velocity ) || (tgt.type == SolveTargetType.Both))
        {
          stargets = stargets + String.Format("v="+Vec2Str(tgt.v)) + ":";
        }
        stargets = stargets + String.Format("t={0:F2}",tgt.t) + "} ";
      }
      // TODO - Missing constraints?
      return string.Format("HopperGuidance: "+msg+" tol="+tol+" minDurationPerThrust="+minDurationPerThrust+" maxThrustsBetweenTargets="+maxThrustsBetweenTargets+" N="+N+" r0="+Vec2Str(r0)+" v0="+Vec2Str(v0)+" g="+g+" Tmin="+Tmin+" Tmax="+Tmax+" amin="+amin+" amax="+amax+" vmax="+vmax+" minDescentAngle="+minDescentAngle+" maxThrustAngle="+maxThrustAngle+" maxLandingThrustAngle="+maxLandingThrustAngle+" apex="+Vec2Str(apex)+" {0}",stargets);
    }
 
    public static double [] BasisWeights(double t, List<float> thrust_t)
    {
      // Returns vector of weights[N] for each acceleration vector at a given time, t
      // for time t in range 0...T, with N vectors
      // T divided into N-1 parts, 0, T/(N-1), 2T/(N-1), (N-1)T/(N-1) == T
      int N = thrust_t.Count;
      double [] w = new double[N];
      for( int j = 0; j < N-1; j++ )
      {
        // Find which two thrust vectors are closest
        if(( t >= thrust_t[j]) && (t <= thrust_t[j+1]))
        {
          double d = (t-thrust_t[j])/(thrust_t[j+1]-thrust_t[j]); // range 0 to 1
          double b = Math.Cos(d*0.5*Math.PI);
          w[j] = b*b;
          w[j+1] = 1-b*b;
          break;
        }
      }
      // Check if t beyond last, assume numerical error and use last thrust
      int k = thrust_t.Count-1;
      if (t > thrust_t[k])
        w[k] = 1;
        
      return w;
    }

    // Calculate weights on position and velocity from thrust vectors up to time tX
    static void RVWeightsToTime(double tX, double dt, List<float> thrusts_t, out double[] wr, out double[] wv)
    {
      int N = thrusts_t.Count;
      double T = thrusts_t[N-1];
      wr = new double[N];
      wv = new double[N];
      for(double t = 0; t < tX; t += dt)
      {
        double [] w = BasisWeights(t,thrusts_t); // Vector for all N weights at time, t
        //System.Console.Error.WriteLine("t="+t);
        //for(int j=0; j<w.Length; j++)
        //  System.Console.Error.WriteLine("  w["+j+"]="+w[j]+" at "+thrusts_t[j]);
        double tr = tX - t - dt; // time remaining
        for(int i = 0 ; i < N ; i++)
        {
          // extra. vel after accel to time, tr + movement during acceleration over dt
          wr[i] += tr*w[i]*dt + 0.5*w[i]*dt*dt;
          // additional velocity over time, dt
          wv[i] += w[i]*dt;
        }
      }
    }

    // Are the constraints satisfied?
    void CheckSolution(double [,] c, int [] ct, double [] x)
    {
      int M = x.Length;
      int L = ct.Length; // Number of constraints
   
      double [] lhs = new double[L];
      for(int i=0; i<L ; i++)
      { 
        for(int j=0; j<M ; j++)
        { 
          lhs[i] += c[i,j] * x[j];
        }
        string cmp = "=";
        if (ct[i] > 0)
          cmp=">";
        if (ct[i] < 0)
          cmp="<";
        {
          System.Console.Error.WriteLine("[{0:F0}] {1:F2} {2} {3:F2}", i, lhs[i], cmp, c[i,M]);
        }
      }
    }

    // It is expected that g is position and acts in the direction Y downwards
    // returns fuel used
    public double GFold(Vector3d r0, Vector3d v0,
                        List<SolveTarget> a_targets, double T,
                        out List<float> o_thrusts_t, // Times of thrusts returned
                        out Vector3d[] a_thrusts, out int a_retval)
    {
      // TODO: Ignoring Nmin and Nmax
      a_retval = 0;
      o_thrusts_t = new List<float>();
      // Put a thrust vector at every target
      o_thrusts_t.Add(0);
      float last_t = 0;
      foreach( SolveTarget tgt in a_targets )
      {
        float tgt_t = tgt.t;
        if (tgt_t < 0)
          tgt_t = (float)T;
        if (tgt_t - last_t > minDurationPerThrust)
        {
          System.Console.Error.WriteLine("ADDING "+last_t+" to "+tgt_t+" minDur="+minDurationPerThrust);
          int num = (int)((tgt_t - last_t)/minDurationPerThrust); // rounds down
          num = Math.Min(num, maxThrustsBetweenTargets);
          float dt = (tgt_t - last_t)/(num+1);
          // if num=1
          // t=last_t, 0.5*(tgt.t+last_t), tgt.t
          for(int i=1 ; i<num ; i++)
          {
            System.Console.Error.WriteLine("i="+i);
            o_thrusts_t.Add(last_t + dt*(1+i));
          }
        }
        o_thrusts_t.Add(tgt_t);
        last_t = tgt_t;
      }
      foreach( float t in o_thrusts_t )
      {
#if (DEBUG)
        System.Console.Error.WriteLine("thrust_t="+t);
#endif
      }
      N = o_thrusts_t.Count;
      a_thrusts = new Vector3d[N];

      List<double> mdCheckTimes = new List<double>();
      double tX;
      if (minDescentAngle >= 0)
      {
        // from start
        for( tX=checkGapFirst; tX<0.5*T; tX=tX*checkGapMult )
          mdCheckTimes.Add(tX);
        // from end
        for( tX=T-checkGapFirst; tX>0.5*T; tX=T-(T-tX)*checkGapMult )
          mdCheckTimes.Add(tX);
      }

      alglib.minqpstate state;
      alglib.minqpreport rep;

      // Constraints
      //
      // Minimise fuel
      // f(x) = SUM |Ti|   // all thrust magnitudes
      // Exact constraints
      //   Final position = rf
      //   Final velocity = vf
      //   Start position = r0
      //   Start velocity = v0
      //   Gravity = g (acts downloads in Y direction)
      //   Initial mass = m
      // Minimise
      //   SUM(thrust)

      dt = T/(N*fidelity);

      // Coefficients of function to minimize of squared terms (none)
      // Note that a is made of N (X,Y,Z) acceleration vectors
      // followed by N thrust magnitudes
      // the acceleration vectors must be < thrust magnitude is every X,Y,Z component
      // so the constraint is to a box current. We could cut the corners off the box
      double[,] a = new double[N*4,N*4];

      // Coefficent of minimise function of linear terms
      // this is the sum of all the thrust magnitudes
      double[] b = new double[N*4];
      for(int i=0;i<N;i++) // thrust magnitudes
        b[N*3+i]=1;

      // Accelation for (ax0,ay0,ax1,ay1...)
      double[] x; // dimensionality of N*4. Each thrust vector as x,y,z and magnitude

      alglib.minqpcreate(N*4, out state);
      alglib.minqpsetquadraticterm(state, a);
      alglib.minqpsetlinearterm(state, b);

      double[] bndl = new double[N*4];
      double[] bndu = new double[N*4];
      for(int i=0;i<N;i++)
      {
        bndl[i*3]   = -amax;
        bndl[i*3+1] = -amax;
        bndl[i*3+2] = -amax;
        bndu[i*3]   = amax;
        bndu[i*3+1] = amax;
        bndu[i*3+2] = amax;
        // thrust magnitudes
        bndl[N*3+i] = amin;
        bndu[N*3+i] = amax;
      }

      alglib.minqpsetbc(state, bndl, bndu);

      int constraints = 6*N; // thrust magnitude limits
      foreach( SolveTarget tgt in a_targets )
      {
        if (tgt.type == SolveTargetType.Position)
          constraints += 3;
        if (tgt.type == SolveTargetType.Velocity)
          constraints += 3;
        if (tgt.type == SolveTargetType.Both)
          constraints += 6;
      }

      constraints++; // For hit the ground in future constraint

      // for minDescentAngle, N points for 4 planes to make square 'cones'
#if (MINDESCENTANGLE)
      if (minDescentAngle >= 0)
        constraints += 4*(mdCheckTimes.Count);
#endif

#if (MAXVELOCITY)
      constraints += 6*N;
#endif

#if (MAXTHRUSTANGLE)
      // Can't handle >90 since the space of possible thrust directions becomes
      // non-convex :-(
      if (maxThrustAngle<90)
        constraints += 4*(N-1);
      if (maxLandingThrustAngle<90)
        constraints += 4;
#endif

#if (MINTHRUST)
      constraints += N;
#endif
      int k=0;
#if (DEBUG)
      System.Console.Error.WriteLine("Constraints="+constraints+" N="+N+" T="+T);
#endif

      int rhs = N*4;
      double [,] c = new double[constraints,rhs+1]; // zeroed?
      int [] ct = new int[constraints]; // type of constraint, =, > or <  (default to 0 -> =)

      // Constrain thrust vectors to be below thrust magnitudes
      for( int i = 0 ; i < N ; i++ )
      {
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector X axis less than thrust magnitude
        k++;
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector X axis greater than -thrust magnitude
        k++;

        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector Y axis less than thrust magnitude
        k++;
        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector Y axis greater than -thrust magnitude
        k++;

        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector Z axis less than thrust magnitude
        k++;
        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector Z axis greater than -thrust magnitude
        k++;
      }

      foreach( SolveTarget tgt in a_targets )
      {
        tX = tgt.t;
        if (tX < 0)
          tX = T;
        RVWeightsToTime(tX,dt,o_thrusts_t,out double[] wr,out double[] wv);
        // Position contraint
        if(( tgt.type == SolveTargetType.Position ) || (tgt.type == SolveTargetType.Both) )
        {
          for(int i = 0 ; i < N ; i++)
          {
            c[k+0,i*3+0] = wr[i]; // X
            c[k+1,i*3+1] = wr[i]; // Y
            c[k+2,i*3+2] = wr[i]; // Z
          }
          c[k+0,rhs] = tgt.r.x - (r0.x + v0.x*tX);
          c[k+1,rhs] = tgt.r.y - (r0.y + v0.y*tX - 0.5*tX*tX*g);
          c[k+2,rhs] = tgt.r.z - (r0.z + v0.z*tX);
          // Must equal
          ct[k+0] = 0;
          ct[k+1] = 0;
          ct[k+2] = 0;
          k += 3;
        }
        // Velocity contraint
        if(( tgt.type == SolveTargetType.Velocity ) || (tgt.type == SolveTargetType.Both) )
        {
          for(int i = 0 ; i < N ; i++)
          {
            c[k+0,i*3+0] = wv[i]; // X
            c[k+1,i*3+1] = wv[i]; // Y
            c[k+2,i*3+2] = wv[i]; // Z
          }
          c[k+0,rhs] = tgt.v.x - v0.x;
          c[k+1,rhs] = tgt.v.y - (v0.y - tX*g);
          c[k+2,rhs] = tgt.v.z - v0.z;
          // Must equal
          ct[k+0] = 0;
          ct[k+1] = 0;
          ct[k+2] = 0;
          k += 3;
        }
      }

      // Find time to hit ground with just gravity
      // 0 = a*x*x + b*x + c;
      // height = 0.5*amax*t*t + v*t;

      // distance = (init vel) * (time) * 0.5

      // Use last Y position target
      // TODO: Calculate final position
      Vector3d tgt_r = a_targets[a_targets.Count-1].r;
      float height = (float)tgt_r.y;
      // Crude height above descent angle
      //float d = (float)Math.Sqrt(tgt_r.x*tgt_r.x + tgt_r.y*tgt_r.y);
      float maxv = Mathf.Sqrt((float)amax*1.7f*height); // should really be 2.0
      RVWeightsToTime(T,dt,o_thrusts_t,out double[] wr2,out double[] wv2);
      for(int i = 0 ; i < N ; i++)
        c[k,i*3+1] = wv2[i]; // Y
      c[k,rhs] = -maxv - (v0.y - T*g);
      ct[k] = 1; // LHS > RHS
      k++;

#if (MINDESCENTANGLE)
      // Constrain N intermediate positions to be within minimumDescentAngle
      foreach( float mdt in mdCheckTimes )
      {
#if (DEBUG)
        System.Console.Error.WriteLine("minDescent angle t="+mdt+" k="+k+"/"+constraints);
#endif
        // No check at t=T
        //double tX = T*((float)(j+1)/(numchecks+1));
        // Get whole weight vector up to time t
        RVWeightsToTime(mdt,dt,o_thrusts_t,out double[] wr,out double[] wv);

        // Calculate Normal for plane to be above (like an upside down pyramid)
        double vx = Math.Sin(minDescentAngle*Math.PI/180.0);
        double vy = Math.Cos(minDescentAngle*Math.PI/180.0);
        double [] V1 = new double [] {vx,vy,0}; // Normal vector of plane to be above
        double [] V2 = new double [] {-vx,vy,0}; // Normal vector of plane to be above
        double [] V3 = new double [] {0,vy,vx}; // Normal vector of plane to be above
        double [] V4 = new double [] {0,vy,-vx}; // Normal vector of plane to be above
        for(int i = 0 ; i < N ; i++)
        {
          // proportions of thrusts[i] for XYZ for position
          // 45 degrees when X<0
          c[k+0,i*3+0] = V1[0] * wr[i]; // X
          c[k+0,i*3+1] = V1[1] * wr[i]; // Y
          c[k+0,i*3+2] = V1[2] * wr[i]; // Z
          // proportions of thrusts[i] for XYZ for position
          // 45 degrees when X<0
          c[k+1,i*3+0] = V2[0] * wr[i]; // X
          c[k+1,i*3+1] = V2[1] * wr[i]; // Y
          c[k+1,i*3+2] = V2[2] * wr[i]; // Z
          // proportions of thrusts[i] for XYZ for position
          // 45 degrees when X<0
          c[k+2,i*3+0] = V3[0] * wr[i]; // X
          c[k+2,i*3+1] = V3[1] * wr[i]; // Y
          c[k+2,i*3+2] = V3[2] * wr[i]; // Z
          // proportions of thrusts[i] for XYZ for position
          // 45 degrees when X<0
          c[k+3,i*3+0] = V4[0] * wr[i]; // X
          c[k+3,i*3+1] = V4[1] * wr[i]; // Y
          c[k+3,i*3+2] = V4[2] * wr[i]; // Z
        }
        // LHS factors to tX[i] and tY[i]
        // Final X + Y + Z
        double rfx = apex.x;
        double rfy = apex.y;
        double rfz = apex.z;
        c[k+0,rhs] = V1[0] * (rfx - (r0[0] + v0[0]*mdt))
                   + V1[1] * (rfy - (r0[1] + v0[1]*mdt - 0.5*mdt*mdt*g))
                   + V1[2] * (rfz - (r0[2] + v0[2]*mdt)); // RHS
        ct[k+0] = 1; // LHS > RHS
        c[k+1,rhs] = V2[0] * (rfx - (r0[0] + v0[0]*mdt))
                   + V2[1] * (rfy - (r0[1] + v0[1]*mdt - 0.5*mdt*mdt*g))
                   + V2[2] * (rfz - (r0[2] + v0[2]*mdt)); // RHS
        ct[k+1] = 1; // LHS > RHS
        c[k+2,rhs] = V3[0] * (rfx - (r0[0] + v0[0]*mdt))
                   + V3[1] * (rfy - (r0[1] + v0[1]*mdt - 0.5*mdt*mdt*g))
                   + V3[2] * (rfz - (r0[2] + v0[2]*mdt)); // RHS
        ct[k+2] = 1; // LHS > RHS
        c[k+3,rhs] = V4[0] * (rfx - (r0[0] + v0[0]*mdt))
                   + V4[1] * (rfy - (r0[1] + v0[1]*mdt - 0.5*mdt*mdt*g))
                   + V4[2] * (rfz - (r0[2] + v0[2]*mdt)); // RHS
        ct[k+3] = 1; // LHS > RHS
        k += 4;
      }
#endif

#if (MAXVELOCITY)
      // Constrain N intermediate positions to be within vmax
      for( int j=0; j<N; j++ )
      {
        // No check at t=0 and t=T
        tX = T*((double)(j+1))/(N+1);
        // Get whole weight vector up to time t
        RVWeightsToTime(tX,dt,o_thrusts_t,out double[] wr,out double[] wv);

        for(int i = 0 ; i < N ; i++)
        {
          c[k+0,i*3+0] = wv[i]; // vx increase by tX
          c[k+1,i*3+0] = wv[i]; // vx increase by tX
          c[k+2,i*3+1] = wv[i]; // vy increase by tX
          c[k+3,i*3+1] = wv[i]; // vy increase by tX
          c[k+4,i*3+2] = wv[i]; // vz increase by tX
          c[k+5,i*3+2] = wv[i]; // vy increase by tX
        }
        c[k+0,rhs] = - v0[0] - vmax;
        ct[k+0] = 1; // incV@tx + v0 - g*tX > -vmax
        c[k+1,rhs] = - v0[0] + vmax;
        ct[k+1] = -1; // incV@tx + v0 - g*tX < -vmax
        c[k+2,rhs] = - v0[1] + g*tX - vmax;
        ct[k+2] = 1; // incV@tx + v0 - g*tX > -vmax
        c[k+3,rhs] = - v0[1] + g*tX + vmax;
        ct[k+3] = -1; // incV@tx + v0 -g*tX < +vmax
        c[k+4,rhs] = - v0[2] - vmax;
        ct[k+4] = 1; // incV@tx + v0 -g*tX > -vmax
        c[k+5,rhs] = - v0[2] + vmax;
        ct[k+5] = -1; // incV@tx + v0 -g*tX < +vmax
        k+=6;
      }
#endif

#if (MAXTHRUSTANGLE)
      // Constrain thrust directions to be within angle of vertical
      for( int i=0; i<N; i++ )
      {
        // Calculate Normal for plane to be above (like an upside down pyramid)
        double ang = maxThrustAngle;
        if (i==N-1)
          ang = maxLandingThrustAngle;
        if (ang < 90)
        {
          double vx = Math.Cos(ang*Math.PI/180.0);
          double vy = Math.Sin(ang*Math.PI/180.0);
          double [] V1 = new double [] {vx,vy,0}; // Normal vector of plane to be above
          double [] V2 = new double [] {-vx,vy,0}; // Normal vector of plane to be above
          double [] V3 = new double [] {0,vy,vx}; // Normal vector of plane to be above
          double [] V4 = new double [] {0,vy,-vx}; // Normal vector of plane to be above
          // proportions of thrusts[i] for XYZ for position
          c[k+0,i*3+0] = V1[0]; // X
          c[k+0,i*3+1] = V1[1]; // Y
          c[k+0,i*3+2] = V1[2]; // Z
          // proportions of thrusts[i] for XYZ for position
          c[k+1,i*3+0] = V2[0]; // X
          c[k+1,i*3+1] = V2[1]; // Y
          c[k+1,i*3+2] = V2[2]; // Z
          // proportions of thrusts[i] for XYZ for position
          c[k+2,i*3+0] = V3[0]; // X
          c[k+2,i*3+1] = V3[1]; // Y
          c[k+2,i*3+2] = V3[2]; // Z
          // proportions of thrusts[i] for XYZ for position
          c[k+3,i*3+0] = V4[0]; // X
          c[k+3,i*3+1] = V4[1]; // Y
          c[k+3,i*3+2] = V4[2]; // Z
          ct[k+0] = 1; // LHS > RHS
          ct[k+1] = 1; // LHS > RHS
          ct[k+2] = 1; // LHS > RHS
          ct[k+3] = 1; // LHS > RHS
          k += 4;
        }
      }
#endif

#if (MINTHRUST)
      // Constrain thrust to be at least amin in any direction
      for( int i=0; i<N; i++ )
      {
          c[k,i*3+1] = 1.0; // thrust weight of Y component
          c[k,rhs] = amin; // may mean thrust is TOO upright
          ct[k] = 1; // LHS > RHS
          k++;
      }
#endif

      if (k!=constraints)
      {
        System.Console.Error.WriteLine("Initialised contraints not equal to allocated: "+k+"!="+constraints);
        return 0;
      }

      // zeroes for equality constraints
      alglib.minqpsetlc(state, c, ct);

      double[] s = new double[N*4];
      for(int i=0;i<N*4;i++)
      {
        s[i] = 1;
      }
      alglib.minqpsetscale(state, s);

#if (DUMP)
      WriteMatrix("a",a,N*3,N*3);
      WriteVector("b",b,N*3);
      WriteVector("bndl",bndl,N*3);
      WriteVector("bndu",bndu,N*3);
      WriteMatrix("c",c,k,N*3+1);
      WriteVector("ct",ct,k);
#endif

      alglib.minqpsetalgodenseipm(state, 0.001);

      alglib.minqpoptimize(state);
      alglib.minqpresults(state, out x, out rep);
      double fuel=0;
      for(int i=0;i<N;i++)
      {
        double tx = x[i*3+0];
        double ty = x[i*3+1];
        double tz = x[i*3+2];
        fuel = fuel + Math.Sqrt(tx*tx+ty*ty+tz*tz)*(T/N);
      }
      //CheckSolution(c,ct,x);

      a_retval = rep.terminationtype;
      if ((a_retval>=1) && (a_retval<=5))
      {
        for(int i=0;i<N;i++)
        {
          a_thrusts[i].x = x[i*3+0];
          a_thrusts[i].y = x[i*3+1];
          a_thrusts[i].z = x[i*3+2];
        }
        System.Console.Error.WriteLine("PASS: T={0:F4} FUEL={1:F2} retval={2:F0}", T, fuel, a_retval);
        retval = a_retval;
      }
      else
      {
        System.Console.Error.WriteLine("FAIL: T={0:F4} FUEL=inf retval={1}", T, a_retval);
        retval = a_retval;
        fuel = 9e+20;
      }
      return fuel + T*timePenalty;
  }

   static double[] convToDouble3(Vector3d v)
   {
     return new double[] { v.x, v.y, v.z };
   }

   // Note that time of the final position constraint a_tr_t and a_tv_t
   // will be adjusted to the end of the time period tested, T
   public double GoldenSearchGFold(Vector3d a_r0, Vector3d a_v0,
                                   List<SolveTarget> a_targets,
                                   out List<float> o_thrusts_t,
                                   out Vector3d [] o_thrusts, out double o_fuel, out int o_retval)
   {
      // Store best solution values
      r0 = a_r0;
      v0 = a_v0;
      targets = a_targets;
      retval = -1;

      // golden section search
      // to find the minimum of f on [a,b]
      // f: a strictly unimodal function on [a,b]
      //
      // example:
      // >>> f = lambda x: (x-2)**2
      // >>> x = golden_search(f, 1, 5)
      // >>> x
      // 2.000009644875678
      double gr = (Math.Sqrt(5) + 1) / 2;
      double fc = 0;
      double fd = 0;
      double a = Tmin;
      double b = Tmax;
      double c,d;
      double last_c,last_d;
      int retvalc=-1,retvald=-1;

      o_thrusts = new Vector3d[N];
      o_thrusts_t = new List<float>();
      Vector3d [] thrustsc = new Vector3d[N];
      Vector3d [] thrustsd = new Vector3d[N];

      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      last_c = c;
      last_d = d;
      while (Math.Abs(c - d) > tol)
      {
        last_c = c;
        last_d = d;
        // TODO - Adjust times on a_tr_t and a_tv_t for be fc and fd
        fc = GFold(a_r0,a_v0,a_targets,c,out o_thrusts_t,out thrustsc,out retvalc);
        fd = GFold(a_r0,a_v0,a_targets,d,out o_thrusts_t,out thrustsd,out retvald);
        if (fc < fd)
            {b = d;}
        else
            {a = c;}

        // We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
      }

      double bestT = 0.5*(a+b);
      // TODO - Adjust times on a_tr_t and a_tv_t for be fc and fd
      o_fuel = GFold(a_r0,a_v0,a_targets,bestT,out o_thrusts_t,out o_thrusts,out o_retval);
      retval = o_retval;
      r0 = a_r0;
      v0 = a_v0;
      if ((o_retval<1) || (o_retval>5))
      {
        System.Console.Error.WriteLine("Fallback "+last_c+" "+last_d+" "+fc+","+fd+","+retvalc+","+retvald);
        // Fall back to c or d positions for 0.5*(a+b) failed
        if ((fc < fd) && ((retvalc>=1) && (retvalc<=5)))
        {
          o_fuel = fc;
          o_thrusts = thrustsc;
          o_retval = retvalc;
          retval = retvalc;
          targets[targets.Count-1].t = (float)last_c + extraTime;
          return last_c;
        }
        if ((fd < fc) && ((retvald>=1) && (retvald<=5)))
        {
          o_fuel = fd;
          o_thrusts = thrustsd;
          o_retval = retvald;
          retval = retvald;
          targets[targets.Count-1].t = (float)last_d + extraTime;;
          return last_d;
        }
        System.Console.Error.WriteLine("FAILED AT T={0}",0.5*(a+b));
        return bestT;
      }
      targets[targets.Count-1].t = (float)bestT + extraTime;
      return bestT;
    }


  }
  public class MainProg
  {
    // All positions and velocities supplied in local space relative to landing target at (0,0,0)
    // but the trajectory, world_traj, is also computed using the transform, given in transform
    static public double MultiPartSolve(ref Solve solver,
                                        ref Trajectory local_traj,
                                        Vector3d local_r, Vector3d local_v,
                                        List<SolveTarget> a_targets,
                                        float g,
                                        out double o_fuel, out int retval)
    {
      //double dt = 0;
      double T = 0, bestT = 0;
      o_fuel = 0;
      retval = -1;

      // Compute trajectory to landing spot
      double fuel;
      Vector3d [] local_thrusts = null;

      // Solve for best time adding each intermediate constraint at a time
      bool done = false;
      List<float> thrusts_t = new List<float>();
      while(!done)
      {
        List<SolveTarget> partial_targets = new List<SolveTarget>();

        int last = 0;
        // Create list of targets up to final undefined target
        foreach( SolveTarget tgt in a_targets )
        {
          partial_targets.Add(tgt);
          if (tgt.t < 0)
            break; // Make this last
          last++;
        }

        // Currently uses intermediate positions, ir[]
        bestT = solver.GoldenSearchGFold(local_r, local_v, partial_targets, out thrusts_t, out local_thrusts, out fuel, out retval);
        System.Console.Error.WriteLine(solver.DumpString());

        if ((retval>=1) && (retval<=5))
        {
          T += solver.T;
          // Fill in T for position target
          a_targets[last].t = (float)bestT;
          // All targets added?
          done = (a_targets.Count == partial_targets.Count);
          o_fuel += fuel;
        }
        else
        {
          o_fuel = 0;
          solver.targets = a_targets; // fill in times
          local_traj = new Trajectory();
          local_traj.Simulate(bestT, local_thrusts, thrusts_t, local_r, local_v, new Vector3d(0,-g,0), solver.dt, 0);
          return 0; // failed part way
        }
      }
      local_traj = new Trajectory();
      local_traj.Simulate(bestT, local_thrusts, thrusts_t, local_r, local_v, new Vector3d(0,-g,0), solver.dt, 0);
      solver.targets = a_targets; // fill in times
      return T;
    }

    static int RunTest(string[] args)
    {
      Solve solver = new Solve();
      Vector3d r0 = Vector3d.zero;
      Vector3d v0 = Vector3d.zero;
      Vector3d rf = Vector3d.zero;
      Vector3d vf = Vector3d.zero;
      bool rfset = false;
      bool vfset = false;
      // Intermediate positions and velocities
      List<SolveTarget> targets = new List<SolveTarget>();
      Vector3d c;
      double d;
      for(int i=0;i<args.Length;i++)
      {
        char [] delim={'='};
        string k = args[i].Split(delim,2)[0];
        string v = args[i].Split(delim,2)[1];
        if (v.StartsWith("[")) {
          v = v.Replace("[","").Replace("]","");
          double x = Convert.ToDouble(v.Split(',')[0]);
          double y = Convert.ToDouble(v.Split(',')[1]);
          double z = Convert.ToDouble(v.Split(',')[2]);
          c = new Vector3d(x,y,z);
          if (k=="r0")
            r0 = c;
          else if (k=="v0")
            v0 = c;
          else if (k=="rf")
          {
            rf = c;
            rfset = true;
          }
          else if (k=="vf")
          {
           vf = c;
           vfset = true;
          }
          else if (k=="apex")
            solver.apex = c;
          else
          {
            System.Console.Error.WriteLine("No such parameter name: "+k);
            return 1;
          }
        } else if (v.StartsWith("T")||v.StartsWith("F"))
        {
          //bool flag = v.StartsWith("T")?true:false;
          System.Console.Error.WriteLine("No such parameter name: "+k);
          return 1;
        } else if (v.StartsWith("{") && (v.EndsWith("}"))) // compound
        {
          v = v.Replace("{","").Replace("}","");
          string [] parts = v.Split(':');
          SolveTarget tgt = new SolveTarget();
          bool rset = false, vset = false;
          for(int j=0; j<parts.Length; j++)
          {
            tgt.t = -1;
            k = parts[j].Split('=')[0];
            v = parts[j].Split('=')[1];
            c = Vector3d.zero;
            d = -1;
            if (v.StartsWith("[") && v.EndsWith("]"))
            {
              v = v.Replace("[","").Replace("]","");
              double x = Convert.ToDouble(v.Split(',')[0]);
              double y = Convert.ToDouble(v.Split(',')[1]);
              double z = Convert.ToDouble(v.Split(',')[2]);
              c = new Vector3d(x,y,z);
            } else {
              d = Convert.ToDouble(v);
            }
            if( k=="r" )
            {
              tgt.r = c;
              tgt.type = (vset)?SolveTargetType.Both:SolveTargetType.Position;
              rset = true;
            }
            else if( k=="v" )
            {
              tgt.v = c;
              tgt.type = (rset)?SolveTargetType.Both:SolveTargetType.Velocity;
            }
            else if( k=="t" )
              tgt.t = (float)d;
            else
              System.Console.Error.WriteLine("No such element {0} in target. Only r, v or t",k);
          }
          targets.Add(tgt);
        } else {
          d = Convert.ToDouble(v);
          if (k=="N")
            solver.N = (int)d; // ignored - output in DumpString() for info
          else if (k=="maxThrustsBetweenTargets")
            solver.maxThrustsBetweenTargets = (int)d;
          else if (k=="minDurationPerThrust")
            solver.minDurationPerThrust = d;
          else if (k=="amin")
            solver.amin = d;
          else if (k=="amax")
            solver.amax = d;
          else if (k=="g")
            solver.g = d;
          else if (k=="minDescentAngle")
            solver.minDescentAngle = d;
          else if (k=="tol")
            solver.tol = d;
          else if (k=="vmax")
            solver.vmax = d;
          else if (k=="Tmin")
            solver.Tmin = d;
          else if (k=="Tmax")
            solver.Tmax = d;
          else if (k=="maxLandingThrustAngle")
            solver.maxLandingThrustAngle = d;
          else if (k=="maxThrustAngle")
            solver.maxThrustAngle = d;
          else if (k=="timePenalty")
            solver.timePenalty = d;
          else if (k=="extraTime")
            solver.extraTime = (float)d;
          else if (k=="checkGapFirst")
            solver.checkGapFirst = d;
          else if (k=="checkGapMult")
            solver.checkGapMult = d;
          else
          {
            System.Console.Error.WriteLine("No such parameter name: "+k);
            return(1);
          }
        }
      }

      // Add final target if set
      if (rfset || vfset)
      {
        SolveTarget final = new SolveTarget();
        if (rfset && vfset)
          final.type = SolveTargetType.Both;
        if (rfset && !vfset)
          final.type = SolveTargetType.Position;
        if (!rfset && vfset)
          final.type = SolveTargetType.Velocity;
        final.r = rf;
        final.v = vf;
        final.t = -1; // unset
        targets.Add(final);
      }

      double fuel;
      int retval;
      Trajectory traj = new Trajectory();
      MultiPartSolve(ref solver, ref traj, r0, v0, targets, (float)solver.g, out fuel, out retval);
      List<string> comments = new List<string>();
      comments.Add(solver.DumpString());
      traj.Write(null, comments);
      System.Console.Error.WriteLine(solver.DumpString());
      return(0);
    }
    static int Main(string[] args)
    {
      if (args.Length == 0) {
        System.Console.Error.WriteLine("usage: Solve.exe k=v ... - set of key=value pairs from r0,v0,rf,vf,Tmin,Tmax,amin,amax,g,minDescentAngle,tol,vmax,ir,iv which also have defaults (ir,iv are intermediates which be specified multiple times)");
        return(1);
      } else {
        return RunTest(args);
      }
    }
  }  
}
