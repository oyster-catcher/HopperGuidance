#define MINDESCENTANGLE
#define MAXTHRUSTANGLE
#define MAXVELOCITY
#define MINTHRUST
//#define DUMP
//#define DEBUG
//#define UNITYDEBUG

using System;
using System.IO;
using System.Collections.Generic;
using KSPAssets;
using UnityEngine;

namespace HopperGuidance
{
  public class Solution
  {
    public Solution(double a_T=0, double a_fuel=0)
    {
      T = a_T;
      fuel = a_fuel;
    }
    public double T;
    public double fuel;
  }

  public class SolveResult
  {
    public double T;
    public double Tstart;
    public double dt;
    public double fuel;
    public List<double> checktimes;
    public ThrustVectorTime [] thrusts;
    public int retval; // return value of Convex optimiser for last call to GFold()
    public Solve inputs;
    public List<Solution> solutions;

    public SolveResult()
    {
      T=0;
      Tstart=0;
      dt=0;
      fuel=float.MaxValue;
      thrusts=null;
      retval=0;
      solutions=new List<Solution>();
    }

    public bool isSolved()
    {
      return (retval>=1) && (retval<=5);
    }

    public string DumpString()
    {
      string msg = isSolved()?"SUCCEED":"FAIL";
      if (thrusts != null)
        msg = "N="+thrusts.Length;
      else
        msg = "N=0";
      return string.Format(msg+" fuel="+fuel+" Tstart="+Tstart+" T="+T+" retval="+retval);
    }
  }

  public class Solve
  {
    public const double toDegrees = 180/Math.PI;
    public const double toRadians = Math.PI/180;

    // Parameters to control solution
    public double Tmin = -1;
    public double Tmax = -1;
    public double TstartMin = 0;
    public double TstartMax = 0;
    public int TstartIntervals = 5;
    public double minDurationPerThrust = 4; // Insert extra thrust vector between targets
    public int maxThrustsBetweenTargets = 1;
    public double checkGapStep = 1; // seconds between each height check
    public int maxChecks = 50; // Maximum number of height checks
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
    public Vector3d apex = Vector3d.zero; // min descent relative to this point
    public bool full = false;
    public float extendTime = 0;

    // Stored input values
    Vector3d r0;
    Vector3d v0;
    List<SolveTarget> targets;
 
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
    public bool Set(string key, string value)
    {
      if (key=="maxThrustsBetweenTargets")
        maxThrustsBetweenTargets = Convert.ToInt32(value);
      else if (key=="minDurationPerThrust")
        minDurationPerThrust = Convert.ToDouble(value);
      else if (key=="minDescentAngle")
         minDescentAngle = Convert.ToDouble(value);
      else if (key=="maxThrustsBetweenTargets")
        maxThrustsBetweenTargets = Convert.ToInt32(value);
      else if (key=="amin")
        amin = Convert.ToDouble(value);
      else if (key=="amax")
        amax = Convert.ToDouble(value);
      else if (key=="vmax")
        vmax = Convert.ToDouble(value);
      else if (key=="tol")
         tol = Convert.ToDouble(value);
      else if (key=="Tmin")
         Tmin = Convert.ToDouble(value);
      else if (key=="Tmax")
         Tmax = Convert.ToDouble(value);
      else if (key=="TstartMin")
         TstartMin = Convert.ToDouble(value);
      else if (key=="TstartMax")
         TstartMax = Convert.ToDouble(value);
      else if (key=="maxThrustAngle")
         maxThrustAngle = Convert.ToDouble(value);
      else if (key=="maxLandingThrustAngle")
         maxLandingThrustAngle = Convert.ToDouble(value);
      else if (key=="timePenalty")
         timePenalty = Convert.ToDouble(value);
      else if (key=="g")
         g = Convert.ToDouble(value);
      else if (key=="extraTime")
         extraTime = (float)Convert.ToDouble(value);
      else if (key=="maxChecks")
         maxChecks = Convert.ToInt32(value);
      else if (key=="checkGapStep")
         checkGapStep = Convert.ToDouble(value);
      else if (key=="extendTime")
         extendTime = (float)Convert.ToDouble(value);
      else
        return false;
      return true;
    }

    // axes it bitfield of X, Y, Z flags
    public static string Vec2Str(Vector3d v, int axes=7)
    {
      string s = "";
      if ((axes & SolveTarget.X) != 0)
        s = s + string.Format("{0:F2}",v.x);
      else
        s = s + "*";
      s = s + ",";
      if ((axes & SolveTarget.Y) != 0)
        s = s + string.Format("{0:F2}",v.y);
      else
        s = s + "*";
      s = s + ",";
      if ((axes & SolveTarget.Z) != 0)
        s = s + string.Format("{0:F2}",v.z);
      else
        s = s + "*";
      return "[" + s + "]";
    }

    public string DumpString()
    {
      string stargets = "";
      Vector3d rf = Vector3d.zero;
      Vector3d vf = Vector3d.zero;
      int rfaxes = 0;
      int vfaxes = 0;
      float ft = -1;
      for(int i=0; i<targets.Count; i++)
      {
        SolveTarget tgt = targets[i];
        if (i == targets.Count-1)
        {
          rf = tgt.r;
          rfaxes = tgt.raxes;
          vf = tgt.v;
          vfaxes = tgt.vaxes;
          ft = tgt.t;
        }
        else
        {
          stargets = stargets + String.Format("target="+Vec2Str(tgt.r,tgt.raxes)) + ":" + tgt.t + " ";
        }
      }
      if (rfaxes != 0)
        stargets = stargets + " rf="+Vec2Str(rf,rfaxes)+":"+ft;
      if (vfaxes != 0)
        stargets = stargets + " vf="+Vec2Str(vf,vfaxes)+":"+ft;
      // TODO - Missing constraints?
      return string.Format("tol="+tol+" minDurationPerThrust="+minDurationPerThrust+" maxThrustsBetweenTargets="+maxThrustsBetweenTargets+" r0="+Vec2Str(r0)+" v0="+Vec2Str(v0)+" g="+g+" Tmin="+Tmin+" Tmax="+Tmax+" amin="+amin+" amax="+amax+" vmax="+vmax+" minDescentAngle="+minDescentAngle+" maxThrustAngle="+maxThrustAngle+" maxLandingThrustAngle="+maxLandingThrustAngle+" extraTime="+extraTime+" TstartMin="+TstartMin+" TstartMax="+TstartMax+" {0}",stargets);
    }
 
    public static double [] BasisWeights(double t, List<float> times)
    {
      // Returns vector of weights[N] for each acceleration vector at a given time, t
      // for time t in range 0...T, with N vectors
      // T divided into N-1 parts, 0, T/(N-1), 2T/(N-1), (N-1)T/(N-1) == T
      int N = times.Count;
      double [] w = new double[N];
      for( int j = 0; j < N-1; j++ )
      {
        // Find which two thrust vectors are closest
        if(( t >= times[j]) && (t <= times[j+1]))
        {
          double d = (t-times[j])/(times[j+1]-times[j]); // range 0 to 1
          double b = Math.Cos(d*0.5*Math.PI);
          w[j] = b*b;
          w[j+1] = 1-b*b;
          break;
        }
      }
      // Check if t beyond last, assume numerical error and use last thrust
      if (t > times[N-1])
        w[N-1] = 1;
        
      return w;
    }

    // Calculate weights on position and velocity from thrust vectors up to time tX
    public static void RVWeightsToTime(double tX, double dt, List<float> times, out double[] wr, out double[] wv)
    {
      // TODO
      int N = times.Count;
      wr = new double[N];
      wv = new double[N];
      for(double t = 0; t < tX; t += dt)
      {
        double [] w = BasisWeights(t,times); // Vector for all N weights at time, t
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

    // Assumes starting from stationary, acceleration by amin(negative) or amax(positive) for half way
    // and then deacceleration. amin is -v.e and amax is +v.e
    static float AxisTime(float d, float amin, float amax, float vmax)
    {
      float t1,t2;

      amin = Mathf.Abs(amin);
      d = Mathf.Abs(d);
      t1 = Mathf.Sqrt(d/(0.5f*amax+0.5f*(amax*amax)/amin));
      t2 = t1*(amax/amin);

      float t1max = vmax/amax;
      float t = t1 + t2;
      if( t1 > t1max )  // can accerate to max. vel
      {
        // remaining distance at vmax
        float rd = d - 0.5f*vmax*(vmax/amax) - 0.5f*vmax*(vmax/amin);
        t = vmax/amax + vmax/amin + rd/vmax;
      }

      System.Console.Error.WriteLine("t1="+t1+" t2="+t2+" d="+d+" amin="+amin+" amax="+amax+" vmax="+vmax+" t1max="+t1max);
      return t;
    }

    public static float EstimateTimeBetweenTargets(Vector3d r0, Vector3d v0, List<SolveTarget> tgts, float amax, float g, float vmax, float maxThrustAngle)
    {
      if (g > amax)
        return -1; // can't even hover
      float t=0;
      // Estimate time to go from stationary at one target to stationary at next, to provide
      // an upper estimate on the solution time
      // Find angle at which vertical acceleration equals gravity
      float ang = (float)(Math.Acos(g/amax) * toDegrees);
      float maxt_sideamax = (float)(Math.Sin(ang * toRadians) * amax);
      System.Console.Error.WriteLine("Estimate: angle (g=vert.accel) = "+ang+" sideamax="+maxt_sideamax+" amax="+amax);
      // Use maxAngle at a lower throttle where vertical accel = g
      if (ang > maxThrustAngle)
      {
        maxt_sideamax = (float)(g*Math.Tan(maxThrustAngle * toRadians));
        System.Console.Error.WriteLine("Estimate: angle>maxThrustAngle ang="+maxThrustAngle+" sideamax="+maxt_sideamax);
      }

      float hor_min = -maxt_sideamax;
      float hor_max = +maxt_sideamax; 
      float ver_min = -g;
      float ver_max = amax-g;

      // Compute position with zero velocity
      if (v0.magnitude > 0.1f)
      {
        Vector3d ca = -(amax-g) * v0/v0.magnitude;
        r0 = r0 + v0*t + 0.5*ca*t*t;
        float stop_hor_t = 0;
        float stop_ver_t = 0;
        if (v0.y > 0)
        {
          // move up to stop
          stop_ver_t = (float)(v0.y/g);
          r0.y = r0.y + v0.y*stop_ver_t - 0.5f*g*stop_ver_t*stop_ver_t;
        }
        else
        {
          // move down to stop
          stop_ver_t = (float)(-v0.y/ver_max);
          r0.y = r0.y + v0.y*stop_ver_t + 0.5f*ver_max*stop_ver_t*stop_ver_t;
        }
        // move sideways to stop
        v0.y = 0;
        if (v0.magnitude > 0.1f)
        {
          stop_hor_t = Mathf.Sqrt((float)(v0.x*v0.x + v0.z*v0.z))/maxt_sideamax;
          r0 = r0 + v0*stop_hor_t - 0.5f*maxt_sideamax*stop_hor_t*stop_hor_t*(v0/v0.magnitude); 
        }
        t = stop_hor_t + stop_ver_t;
        System.Console.Error.WriteLine("Position with zero velocity "+r0+" hor_stop_t="+stop_hor_t+" ver_top_t="+stop_ver_t);
      }
      // r0 represents stationary position after velocity cancelled out

      foreach( SolveTarget tgt in tgts )
      {
        float dx = (float)(tgt.r.x - r0.x);
        float dy = (float)(tgt.r.y - r0.y);
        float dz = (float)(tgt.r.z - r0.z);
        float d = Mathf.Sqrt(dx*dx+dz*dz);
        // Compute time to move in each orthogonal axes: horizontally or vertically
        float hor_t = AxisTime(d, hor_min*0.5f, hor_max*0.5f, vmax);
        float ver_t = AxisTime(dy, ver_min, ver_max, vmax);
        t = t + hor_t + ver_t;
      }
      return t;
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
    // returns fuel used or float.MaxFloat is no solution
    public SolveResult GFold(Vector3d a_r0, Vector3d a_v0,
                             List<SolveTarget> a_targets, double T, double Tstart)
    {
      // Save inputs
      r0 = a_r0;
      v0 = a_v0;
      targets = a_targets;
      SolveResult result = new SolveResult();

      result.inputs = this;
      // Apex is for the final target for min descent angle
      // If not in min descent region use ground and min descent angle=0 to keep above ground
      Vector3d ground = r0;
      foreach(SolveTarget tgt in a_targets)
      {
        if (tgt.r.y < ground.y)
          ground = tgt.r;
      }
      if( apex.y < ground.y)
        ground = apex;
      // Hack to ensure starting position is off ground by moving starting
      // position up by 1m
      if (ground.y > r0.y - 1)
        r0.y = r0.y + 1;

      // Quick check if Tstart is large than time to hit ground
      double tHit = -(r0.y-ground.y)/(a_v0.y-0.5*g*g);
      if (tHit < Tstart)
        return result;

      // TODO: Ignoring Nmin and Nmax
      List<float> thrust_times = new List<float>();
      // Put a thrust vector at every target
      thrust_times.Add((float)Tstart);
      float last_t = (float)Tstart;
      float last_target_t = 0;
      foreach( SolveTarget tgt in a_targets )
      {
        float tgt_t = tgt.t;
        if ((tgt.raxes != 0) && (tgt_t > 0))
          last_target_t = tgt.t;
        if (tgt_t < 0)
          tgt_t = (float)T;
        if (tgt_t - last_t > minDurationPerThrust)
        {
#if (DEBUG)
          System.Console.Error.WriteLine("ADDING "+last_t+" to "+tgt_t+" minDur="+minDurationPerThrust);
#endif
           int num = (int)((tgt_t - last_t)/minDurationPerThrust); // rounds down
           num = Math.Min(num, maxThrustsBetweenTargets);
           float gap = (tgt_t - last_t)/(num+1);
           for(int i=0 ; i<num ; i++)
           {
             float t = last_t + gap*(1+i);
             if (t > Tstart)
               thrust_times.Add(t);
           }
        }
        if (tgt_t > Tstart)
          thrust_times.Add(tgt_t);
        last_t = tgt_t;
      }
      // Final target isn't a landing point so this is a multi-part solution, don't restrict
      // this part to the minDescentAngle
      if ((a_targets[a_targets.Count-1].raxes == 0) || (a_targets[a_targets.Count-1].vaxes == 0))
        last_target_t = (float)T;

      // Final target as position and velocity constrained (i.e. land target)
      if( minDescentAngle > 0)
      {
        if ((a_targets[a_targets.Count-1].raxes != 0) && (a_targets[a_targets.Count-1].vaxes != 0))
        {
          Vector3 rl = r0;
          Vector3 rf = a_targets[a_targets.Count-1].r;
          if( a_targets.Count > 1 )
            rl = a_targets[a_targets.Count-2].r; // target before final
          // Find descent angle from this point to final target
          float fx = rf.x - rl.x;
          float fz = rf.z - rl.z;
          float adj = Mathf.Sqrt(fx*fx+fz*fz);
          float hyp = (rf-rl).magnitude;
          float ang = (float)(Mathf.Acos(adj/hyp)*180/Math.PI);
          if (ang < minDescentAngle)
          {
            System.Console.Error.WriteLine("Minimum descent possible is "+ang);
            minDescentAngle = ang;
          }
        }
      }

      int N = thrust_times.Count;
      result.thrusts = new ThrustVectorTime[N];
      for(int i=0; i<thrust_times.Count; i++)
      {
#if (DEBUG)
        System.Console.Error.WriteLine("thrust_t="+thrust_times[i]);
#endif
#if (UNITYDEBUG)
        Debug.Log("thrust_t="+thrust_times[i]);
#endif
        result.thrusts[i] = new ThrustVectorTime();
        result.thrusts[i].v = Vector3d.zero;
        result.thrusts[i].t = thrust_times[i];
      }

      List<double> checktimes = new List<double>();
      double tX;
      if (minDescentAngle >= 0)
      {
        // check time (for below ground and min descent angle)
        // every checkGapStep seconds, but limit number of checks
        // to checkMax
        double step = checkGapStep;
        if ((T/step) > maxChecks)
          step = (double)T/maxChecks;
        for( tX=step; tX<=T-step; tX=tX+step )
          checktimes.Add(tX);
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

      double dt = T/(N*fidelity);

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
        if ((tgt.raxes & SolveTarget.X) != 0)
          constraints++;
        if ((tgt.raxes & SolveTarget.Y) != 0)
          constraints++;
        if ((tgt.raxes & SolveTarget.Z) != 0)
          constraints++;
        if ((tgt.vaxes & SolveTarget.X) != 0)
          constraints++;
        if ((tgt.vaxes & SolveTarget.Y) != 0)
          constraints++;
        if ((tgt.vaxes & SolveTarget.Z) != 0)
          constraints++;
      }

      if (a_targets[a_targets.Count-1].vaxes == 0)
        constraints++; // For hit the ground in future constraint

      // for minDescentAngle, N points for 4 planes to make square 'cones'
#if (MINDESCENTANGLE)
      foreach( float mdt in checktimes )
      {
        if( mdt < last_target_t )
          constraints++; // just check Y above ground
        else
          constraints += 4;
      }
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
      double[] s = new double[N*4]; // scale

      // set default scaling
      for(int i=0;i<N*4;i++)
        s[i] = 1;

      // Constrain thrust vectors to be below thrust magnitudes
      for( int i = 0 ; i < N ; i++ )
      {
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector X axes less than thrust magnitude
        k++;
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector X axes greater than -thrust magnitude
        k++;

        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector Y axes less than thrust magnitude
        k++;
        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector Y axes greater than -thrust magnitude
        k++;

        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Means thrust vector Z axes less than thrust magnitude
        k++;
        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Means thrust vector Z axes greater than -thrust magnitude
        k++;
      }

      foreach( SolveTarget tgt in a_targets )
      {
        tX = tgt.t;
        if (tX < 0)
          tX = T;
#if (DEBUG)
        System.Console.Error.WriteLine("tX="+tX+" tgt.r="+tgt.r+" tgt.v="+tgt.v+" T="+T);
#endif
        RVWeightsToTime(tX,dt,thrust_times,out double[] wr,out double[] wv);
        // Position contraint
        if( tgt.raxes !=0 )
        {
          int k_start = k;
          for(int i = 0 ; i < N ; i++)
          {
            int kk = k;
            if ((tgt.raxes & SolveTarget.X) != 0)
              c[kk++,i*3+0] = wr[i]; // X
            if ((tgt.raxes & SolveTarget.Y) != 0)
              c[kk++,i*3+1] = wr[i]; // Y
            if ((tgt.raxes & SolveTarget.Z) != 0)
              c[kk++,i*3+2] = wr[i]; // Z
          }
          k = k_start;
          // Must equal
          if ((tgt.raxes & SolveTarget.X) != 0)
          {
            c[k,rhs] = tgt.r.x - (r0.x + v0.x*tX);
            ct[k++] = 0;
          }
          if ((tgt.raxes & SolveTarget.Y) != 0)
          {
            c[k,rhs] = tgt.r.y - (r0.y + v0.y*tX - 0.5*tX*tX*g);
            ct[k++] = 0;
          }
          if ((tgt.raxes & SolveTarget.Z) != 0)
          {
            c[k,rhs] = tgt.r.z - (r0.z + v0.z*tX);
            ct[k++] = 0;
          }
        }
        // Velocity contraint
        if( tgt.vaxes != 0)
        {
          int k_start = k;
          for(int i = 0 ; i < N ; i++)
          {
            int kk = k;
            if ((tgt.vaxes & SolveTarget.X) != 0)
              c[kk++,i*3+0] = wv[i]; // X
            if ((tgt.vaxes & SolveTarget.Y) != 0)
              c[kk++,i*3+1] = wv[i]; // Y
            if ((tgt.vaxes & SolveTarget.Z) != 0)
              c[kk++,i*3+2] = wv[i]; // Z
          }
          k = k_start;
          // Must equal
          if ((tgt.vaxes & SolveTarget.X) != 0)
          {
            c[k,rhs] = tgt.v.x - v0.x;
            ct[k++] = 0;
          }
          if ((tgt.vaxes & SolveTarget.Y) != 0)
          {
            c[k,rhs] = tgt.v.y - (v0.y - tX*g);
            ct[k++] = 0;
          }
          if ((tgt.vaxes & SolveTarget.Z) != 0)
          {
            c[k,rhs] = tgt.v.z - v0.z;
            ct[k++] = 0;
          }
        }
      }

      if (a_targets[a_targets.Count-1].vaxes == 0)
      {
        // Use lowest Y position as ground - not true but this only to avoid
        // making partial trajectories don't allow time to not hit ground
        // TODO: Calculate final position
        Vector3d tgt_r = a_targets[a_targets.Count-1].r;
        float height = (float)tgt_r.y - (float)ground.y;
        float maxv = Mathf.Sqrt((float)(amax-g)*height*2);
        RVWeightsToTime(T,dt,thrust_times,out double[] wr2,out double[] wv2);
        for(int i = 0 ; i < N ; i++)
          c[k,i*3+1] = wv2[i]; // Y
        c[k,rhs] = -maxv - (v0.y - T*g);
        ct[k] = 1; // LHS > RHS
        k++;
      }

#if (MINDESCENTANGLE)
      // Constrain N intermediate positions to be within minimumDescentAngle
      foreach( float mdt in checktimes )
      {
        // No check at t=T
        //double tX = T*((float)(j+1)/(numchecks+1));
        // Get whole weight vector up to time t
        RVWeightsToTime(mdt,dt,thrust_times,out double[] wr,out double[] wv);

        // Calculate Normal for plane to be above (like an upside down pyramid)
        // ensure within 45 degress of landing in last 3 seconds
        // TODO: And beyond previous target
        if( mdt < last_target_t )
        {
          for(int i = 0 ; i < N ; i++)
            c[k,i*3+1] = wr[i]; // Y
          c[k,rhs] = ground.y - (r0[1] + v0[1]*mdt - 0.5*mdt*mdt*g);
          ct[k] = 1; // LHS > RHS
          k++;
        }
        else
        {
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
          // scale down to avoid too much weight
          for(int ii = 0 ; ii < N*4+1 ; ii++) 
          {
            c[k+0,ii] = c[k+0,ii] * 0.01f;
            c[k+1,ii] = c[k+1,ii] * 0.01f;
            c[k+2,ii] = c[k+2,ii] * 0.01f;
            c[k+3,ii] = c[k+3,ii] * 0.01f;
          }
          k += 4;
        }
      }
#endif

#if (MAXVELOCITY)
      // Constrain N intermediate positions to be within vmax
      for( int j=0 ; j<N ; j++ )
      {
        // No check at t=0 and t=T
        tX = T*((double)(j+1))/(N+1);
        // Get whole weight vector up to time t
#if (DEBUG)
        System.Console.Error.WriteLine("tX="+tX+" T="+T+" N="+N);
#endif
        RVWeightsToTime(tX,dt,thrust_times,out double[] wr,out double[] wv);

        for(int i = 0 ; i < N ; i++)
        {
          c[k+0,i*3+0] = wv[i]; // vx increase by tX
          c[k+1,i*3+0] = wv[i]; // vx increase by tX
          c[k+2,i*3+1] = wv[i]; // vy increase by tX
          c[k+3,i*3+1] = wv[i]; // vy increase by tX
          c[k+4,i*3+2] = wv[i]; // vz increase by tX
          c[k+5,i*3+2] = wv[i]; // vz increase by tX
        }
        c[k+0,rhs] = - v0[0] - vmax;
        ct[k+0] = 1; // incV@tx + v0 > -vmax
        c[k+1,rhs] = - v0[0] + vmax;
        ct[k+1] = -1; // incV@tx + v0 < -vmax
        c[k+2,rhs] = - v0[1] + g*tX - vmax;
        ct[k+2] = 1; // incV@tx + v0 - g*tX > -vmax
        c[k+3,rhs] = - v0[1] + g*tX + vmax;
        ct[k+3] = -1; // incV@tx + v0 -g*tX < +vmax
        c[k+4,rhs] = - v0[2] - vmax;
        ct[k+4] = 1; // incV@tx + v0 > -vmax
        c[k+5,rhs] = - v0[2] + vmax;
        ct[k+5] = -1; // incV@tx + v0 < +vmax
        // scale down to avoid too much weight
        for(int i = 0 ; i < N*4+1 ; i++) 
        {
          c[k+0,i] = c[k+0,i] * 0.01f;
          c[k+1,i] = c[k+1,i] * 0.01f;
          c[k+2,i] = c[k+2,i] * 0.01f;
          c[k+3,i] = c[k+3,i] * 0.01f;
          c[k+4,i] = c[k+4,i] * 0.01f;
          c[k+5,i] = c[k+5,i] * 0.01f;
        }
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
        return result;
      }

      // zeroes for equality constraints
      alglib.minqpsetlc(state, c, ct);

      alglib.minqpsetscale(state, s);

#if (DUMP)
      WriteMatrix("a",a,N*3,N*3);
      WriteVector("b",b,N*3);
      WriteVector("bndl",bndl,N*3);
      WriteVector("bndu",bndu,N*3);
      WriteMatrix("c",c,k,N*4+1);
      WriteVector("ct",ct,k);
#endif

      alglib.minqpsetalgodenseipm(state, 0.00001);
      //alglib.minqpsetscaleautodiag(state);
      alglib.minqpoptimize(state);
      alglib.minqpresults(state, out x, out rep);
      result.dt = dt;
      result.fuel = 0;
      result.T = T;
      for(int i=0;i<N;i++)
      {
        double tx = x[i*3+0];
        double ty = x[i*3+1];
        double tz = x[i*3+2];
        result.fuel += Math.Sqrt(tx*tx+ty*ty+tz*tz)*(T/N);
      }
      //CheckSolution(c,ct,x);
      result.T = T;
      result.retval = rep.terminationtype;
      result.checktimes = checktimes;
      if (result.isSolved()) // uses retval
      {
        for(int i=0;i<N;i++)
        {
          result.thrusts[i].v.x = x[i*3+0];
          result.thrusts[i].v.y = x[i*3+1];
          result.thrusts[i].v.z = x[i*3+2];
#if (DEBUG)
          System.Console.Error.WriteLine("T["+i+"]="+(Vector3)result.thrusts[i].v);
#endif
        }
        System.Console.Error.WriteLine("PASS: T={0:F2} Tstart={1:F2} FUEL={2:F2} retval={3:F0}", T, Tstart, result.fuel, result.retval);
      }
      else
      {
        System.Console.Error.WriteLine("FAIL: T={0:F2} Tstart={1:F2} FUEL=inf retval={2}", T, Tstart, result.retval);
        result.fuel = 9e+20;
      }
      result.Tstart = Tstart;
      return result;
  }

   static double[] convToDouble3(Vector3d v)
   {
     return new double[] { v.x, v.y, v.z };
   }

   public SolveResult FullSearchGFold(Vector3d a_r0, Vector3d a_v0,
                                     List<SolveTarget> a_targets, double Tstart)
   {
     SolveResult result = new SolveResult();
     double t = Tmin;

     while(t<Tmax)
     {
       SolveResult res = GFold(a_r0,a_v0,a_targets,t,Tstart);
       System.Console.Error.WriteLine("best_T="+res.T+" best_fuel="+res.fuel+" retval="+res.retval);
       if (res.isSolved())
       {
         if (res.fuel < result.fuel)
           // Save values
           result = res;
         else
           return result;
       }
       else
         // If failure to solve and already have solution stop
         if (result.isSolved())
           return result;
       t = t + tol;
     }
     return result;
   }

   // Note that time of the final position constraint a_tr_t and a_tv_t
   // will be adjusted to the end of the time period tested, T
   public SolveResult GoldenSearchGFold(Vector3d a_r0, Vector3d a_v0,
                                   List<SolveTarget> a_targets, double Tstart)
   {
      SolveResult result = null;
      List<Solution> solutions = new List<Solution>();
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
      double a = Tmin;
      double b = Tmax;
      double c,d;
      double last_c,last_d;

      SolveResult resc = new SolveResult();
      SolveResult resd = new SolveResult();
      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      last_c = c;
      last_d = d;
      while (Math.Abs(c - d) > tol)
      {
        last_c = c;
        last_d = d;
        // TODO - Adjust times on a_tr_t and a_tv_t for be fc and fd
        resc = GFold(a_r0,a_v0,a_targets,c,Tstart);
        resd = GFold(a_r0,a_v0,a_targets,d,Tstart);
        // Store duration and fuel
        if (resc.isSolved())
          solutions.Add(new Solution(c,resc.fuel));
        if (resd.isSolved())
          solutions.Add(new Solution(d,resd.fuel));
        if ((!resc.isSolved()) && (!resd.isSolved()))
        {
          // Where to go - nowhere is good!
          double step = (b-a)*0.1f;
          c = a;
          //System.Console.Error.WriteLine("STEPPING a="+a+" b="+b+" step="+step);
          while(c <= b+0.01f)
          {
            resc = GFold(a_r0,a_v0,a_targets,c,Tstart);
            //System.Console.Error.WriteLine("STEPPING: T="+c+" fuel="+resc.fuel+" retval="+resc.retval);
            if (resc.isSolved())
            {
              if (resc.fuel < resd.fuel) // best so far
              {
                d = c; // best time
                resd = resc; // keep best
              }
            }
            c = c + step;
          }
          if (!resd.isSolved()) // no solutions found
            return resd;
          // New search interval
          a = Math.Max(d - step,Tmin);
          b = Math.Min(d + step,Tmax);
        }
        else
        {
          if (resc.fuel < resd.fuel)
            b = d;
          else
            a = c;
        }

        // We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
      }

      double bestT = 0.5*(a+b);
      // TODO - Adjust times on a_tr_t and a_tv_t for be fc and fd
      result = GFold(a_r0,a_v0,a_targets,bestT,Tstart);
      
      if (!result.isSolved())
      {
        System.Console.Error.WriteLine("Fallback "+last_c+" "+last_d+" "+resc.fuel+","+resd.fuel+","+resc.retval+","+resd.retval);
        // Fall back to c or d positions for 0.5*(a+b) failed
        if ((resc.fuel < resd.fuel) && (resc.isSolved()))
        {
          a_targets[a_targets.Count-1].t = (float)resc.T + extraTime;
          resc.solutions = solutions;
          return resc;
        }
        if ((resd.fuel < resc.fuel) && (resd.isSolved()))
        {
          a_targets[a_targets.Count-1].t = (float)resd.T + extraTime;
          resd.solutions = solutions;
          return resd;
        }
        System.Console.Error.WriteLine("FAILED AT T={0}",0.5*(a+b));
        return result;
      }
      else
      {
        solutions.Add(new Solution(bestT, result.fuel));
        a_targets[a_targets.Count-1].t = (float)result.T + extraTime;
        result.solutions = solutions;
        return result;
      }
    }
  }

  public class MainProg
  {
    // Loop over start times. No thrust below start times
    static public SolveResult MultiPartSolve(ref Solve solver,
                                            ref Trajectory local_traj,
                                            Vector3d local_r, Vector3d local_v,
                                            ref List<SolveTarget> a_targets,
                                            float g,
                                            float extendTime, bool correct=false)
    {
      double tStart = solver.TstartMin;
      double Tmin = solver.Tmin;
      double Tmax = solver.Tmax;
      SolveResult result = new SolveResult();
      while(tStart <= solver.TstartMax)
      {
        List<SolveTarget> tTargets = new List<SolveTarget>(a_targets);
        solver.Tmin = Tmin;
        solver.Tmax = Tmax;
        result = MultiPartSolve2(ref solver,ref local_traj,local_r,local_v,ref tTargets,g,extendTime,tStart,correct);
        if (result.isSolved())
        {
          a_targets = tTargets;
          return result;
        }
        tStart = tStart + solver.minDurationPerThrust;
      }
      return result;
    }

    // All positions and velocities supplied in local space relative to landing target at (0,0,0)
    // but the trajectory, world_traj, is also computed using the transform, given in transform
    static public SolveResult MultiPartSolve2(ref Solve solver,
                                              ref Trajectory local_traj,
                                              Vector3d local_r, Vector3d local_v,
                                              ref List<SolveTarget> a_targets,
                                              float g,
                                              float extendTime, double Tstart,
                                              bool correct=false)
    {
      SolveResult result = new SolveResult();
#if (UNITYDEBBUG)
      Debug.Log("local_r="+(Vector3)local_r);
      Debug.Log("local_v="+(Vector3)local_v);
#endif
      solver.apex = a_targets[a_targets.Count-1].r;
      solver.apex -= a_targets[a_targets.Count-1].v * extendTime * 0.5f;
      double Tmin = solver.Tmin;
      double Tmax = solver.Tmax;

      // Solve for best time adding each intermediate constraint at a time
      bool done = false;
      double tol = solver.tol;
      // Current last position in trajectory
      Vector3d lr = local_r;
      Vector3d lv = local_v;
      while(!done)
      {
        List<SolveTarget> partial_targets = new List<SolveTarget>();

        int last = 0;
        // Create list of targets up to final undefined target
        List<SolveTarget> next_targets = new List<SolveTarget>();
        foreach( SolveTarget tgt in a_targets )
        {
#if (UNITYDEBUG)
          Debug.Log("target="+(Vector3)tgt.r);
#endif
          partial_targets.Add(tgt);
          if (tgt.t < 0)
          {
            next_targets.Add(tgt);
            break; // Make this last
          }
          last++;
        }

        // Move final target off ground a little as extend time was extend trajectory vertically downloads
        if (a_targets.Count == partial_targets.Count)
          partial_targets[partial_targets.Count-1].r -= partial_targets[partial_targets.Count-1].v * extendTime * 0.5f;
 
        solver.Tmin = (Tmin>=0) ? Tmin : result.T;
        if (Tmax < 0)
        { 
          float estT = Solve.EstimateTimeBetweenTargets(lr, lv, next_targets, (float)solver.amax, (float)solver.g, (float)solver.vmax, (float)solver.maxThrustAngle);
          solver.Tmax = result.T + estT*1.2f;
          System.Console.Error.WriteLine("Estimated Tmax="+solver.Tmax);
        }
        else
          solver.Tmax = Tmax;

        // Currently uses intermediate positions, ir[]
        // use rough tolerance for intermediate targets and reduce to actual set value for
        // last target when complete final trajectory is computed
        solver.tol = (a_targets.Count == partial_targets.Count)?tol:tol*4;
        if (solver.full)
          result = solver.FullSearchGFold(local_r, local_v, partial_targets, Tstart);
        else
          result = solver.GoldenSearchGFold(local_r, local_v, partial_targets, Tstart);
        System.Console.Error.WriteLine(solver.DumpString()+" "+result.DumpString());

        if (result.isSolved())
        {
          solver.Tmin = result.T;
          // All targets added?
          done = (a_targets.Count == partial_targets.Count);
        }
        else
        {
          local_traj = new Trajectory();
          local_traj.Simulate(result.T, result.thrusts, local_r, local_v, new Vector3d(0,-g,0), result.dt);
          return result; // failed part way through
        }
        lr = partial_targets[partial_targets.Count-1].r; // assume got to last target
      }
      local_traj = new Trajectory();
      local_traj.Simulate(result.T, result.thrusts, local_r, local_v, new Vector3d(0,-g,0), result.dt);
      Vector3d vf = a_targets[a_targets.Count-1].v;
      if (correct)
      {
        local_traj.CorrectForTargets(a_targets);
      }
      local_traj.Extend(vf, new Vector3d(0,-g,0), extendTime);
      return result;
    }


    static public bool RunTest(Dictionary<string,string> args, ref Solve solver, out SolveResult result, out Trajectory traj)
    {
      bool correct = false;
      result = new SolveResult();
      traj = new Trajectory();
      solver.Tmax = -1; // means estimate if not set
      Vector3d r0 = Vector3d.zero;
      Vector3d v0 = Vector3d.zero;
      Vector3d rf = Vector3d.zero;
      Vector3d vf = Vector3d.zero;
      int rfaxes = 0;
      int vfaxes = 0;
      // Intermediate positions and velocities
      List<SolveTarget> targets = new List<SolveTarget>();
      float t = -1;

      foreach(var arg in args)
      {
        if (arg.Key=="r0")
          r0 = HGUtils.ToVector3d(arg.Value);
        if (arg.Key=="v0")
          v0 = HGUtils.ToVector3d(arg.Value);
        if (arg.Key=="rf") {
          rf = HGUtils.ToVector3d(arg.Value);
          rfaxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
        }
        if (arg.Key=="vf") {
          vf = HGUtils.ToVector3d(arg.Value);
          vfaxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
        }
        // TODO: Handle : to specify time
        if (arg.Key=="target")
        {
          SolveTarget tgt = new SolveTarget();
          tgt.r = HGUtils.ToVector3d(arg.Value);
          tgt.raxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
          tgt.vaxes = 0;
          tgt.t = t;
          targets.Add(tgt);
        }
      }

      // Add final target if set
      if ((rfaxes != 0)||(vfaxes != 0))
      {
        SolveTarget final = new SolveTarget();
        final.r = rf;
        final.raxes = rfaxes;
        final.v = vf;
        final.vaxes = vfaxes;
        final.t = -1; // unset
        targets.Add(final);
      }
      else
      {
        if (targets.Count == 0)
        {
          System.Console.Error.WriteLine("Must set targets with target= or set rf= and/or vf=");
          return false;
        }
      }

      // targets gets times filled in if set to -1
      result = MainProg.MultiPartSolve(ref solver, ref traj, r0, v0, ref targets, (float)solver.g, solver.extendTime, correct);
      return result.isSolved();
    }
  }  
}
