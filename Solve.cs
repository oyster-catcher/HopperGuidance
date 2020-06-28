using System;
using UnityEngine;
using KSPAssets;

namespace HopperGuidance
{
  public class Solve
  {
    public static double [] BasisWeights(double t, double T, int N)
    {
      // Returns vector of weights[N] for each acceleration vector at a given time, t
      // for time t in range 0...T, with N vectors
      // T divided into N-1 parts, 0, T/(N-1), 2T/(N-1), (N-1)T/(N-1) == T
      double [] w = new double[N];
      double dt = T/(N-1);   // distance between basis points in time
      for(int j=0; j<N; j++)
      {
        double d = j - (t/dt);
        double b = 0; // when outside local range
        if ((d>-1) && (d<1))
        {
          b = Math.Cos(d*0.5*Math.PI);
        }
        w[j] = b*b;
      }
      return w;
    }

    // Calculate weights on position and velocity from thrust vectors up to time tX
    public static void RVWeightsToTime(double tX, double T,int N, double dt, out double[] wr, out double[] wv)
    {
      wr = new double[N];
      wv = new double[N];
      for(double t = 0; t < tX; t += dt)
      {
        double tr = tX - t - dt; // remaining time after applying acceleration
        double [] w = BasisWeights(t,T,N); // Vector for all N weights at time, t
        for(int i = 0 ; i < N ; i++)
        {
          wr[i] += tr*w[i]*dt + 0.5*w[i]*dt*dt;
          wv[i] += w[i]*dt;
        }
      }
    }

    public static void OrthogonalVectors(Vector3d v, out Vector3d a, out Vector3d b)
    {
      a = new Vector3d(-v.y,v.x,0);
      if (Vector3d.Dot(v,a)<0.001) // if looks parallel try another
      {
        a = new Vector3d(0,v.z,-v.y);
      }
      b = Vector3d.Cross(v,a);
      System.Console.WriteLine("v="+v+" a="+a+" b="+b);
      System.Console.WriteLine("v.a="+Vector3d.Dot(v,a)+" v.b="+Vector3d.Dot(v,b));
    }

    // It is expected that g is position and acts in the direction Y downwards
    public static double gfold(double[] r0, double[] v0, double[] att0, double[] rf, double[] vf, double T, int N, double g, double amax, double minDescentAngle, double minFinalDescentAngle, double amin, out double[,] thrusts, out int retval)
    {
      int fidelity = 20;
      thrusts = null;

      alglib.minqpstate state;
      alglib.minqpreport rep;

      // Constraints
      //
      // Minimise fuel
      // f(x) = SUM |Ti|   // all thrusts
      // Exact constraints
      //   Final position = rf
      //   Final velocity = vf
      //   Start position = r0
      //   Start velocity = v0
      //   Gravity = g (acts downloads in Y direction)
      //   Initial mass = m
      // Ideally minimise
      //   SUM(|thrust|)
      // but due to quadratic constraints
      // Minimise
      //   SUM(thrust^2)
      // Special constraints
      //   Initial thrust constraints to be parallel to initial attitude

      double dt = (T/(double)N)/fidelity;

      // Coefficients of function to minimize of squared terms (a bad approximation)
      // we really want to minimise magnitude of the accelerations
      double[,] a = new double[N*3,N*3];
      for(int i=0;i<N*3;i++)
      {
        a[i,i]=1;
      }
      // Coefficent of minimise function of linear terms, zeroes
      double[] b = new double[N*3];

      // Accelation for (ax0,ay0,ax1,ay1...)
      double[] x; // dimensionality of N*2

      alglib.minqpcreate(N*3, out state);
      // Minimise thrust^2
      alglib.minqpsetquadraticterm(state, a);
      alglib.minqpsetlinearterm(state, b);

      double[] bndl = new double[N*3];
      double[] bndu = new double[N*3];
      for(int i=0;i<N;i++)
      {
        bndl[i*3]   = -amax;
        bndl[i*3+1] = -amax;
        bndl[i*3+2] = -amax;
        bndu[i*3]   = amax;
        bndu[i*3+1] = amax;
        bndu[i*3+2] = amax;
      }

      //TODO: Need?
      //alglib.minqpsetstartingpoint(state, x0);
      alglib.minqpsetbc(state, bndl, bndu);

      // 1 constraint on final X position
      // 1 constraint on final Y position
      // 1 constraint on final X vel.
      // 1 constraint on final Y vel.

      int constraints = 3+3; // rf, vf

      // Constrain initial attitude
      //constraints += 3;

      // for minDescentAngle, N points for 4 planes to make square 'cones'
      constraints += 4*N;
      int k=0;

      // for amin constraint to keep thrust downwards
      constraints += N;

      double [,] c = new double[constraints,N*3+1]; // zeroed?
      int [] ct = new int[constraints]; // type of constrain, =, > or <  (default to 0 -> =)

      RVWeightsToTime(T,T,N,dt,out double [] wr,out double [] wv);
      for(int i = 0 ; i < N ; i++)
      {
        // proportions of thrusts[i] for XYZ for position
        c[k+0,i*3+0]  = wr[i];
        c[k+1,i*3+1]  = wr[i];
        c[k+2,i*3+2]  = wr[i];
        // final r function sums to this
        c[k+0,N*3] = rf[0] - (r0[0] + v0[0]*T); // X constant
        c[k+1,N*3] = rf[1] - (r0[1] + v0[1]*T - 0.5*T*T*g); // Y constant
        c[k+2,N*3] = rf[2] - (r0[2] + v0[2]*T); // Z constant
        // proportions of thrusts[i] for XYZ for velocity
        c[k+3,i*3+0]  = wv[i];
        c[k+4,i*3+1]  = wv[i];
        c[k+5,i*3+2]  = wv[i];
        // final v function sums to this
        c[k+3,N*3] = vf[0] - v0[0];
        c[k+4,N*3] = vf[1] - (v0[1] - T*g);
        c[k+5,N*3] = vf[2] - v0[2];
      }
      k += 6;

      // Constrain initial thrust to be the same as att0
      // which will be small by will point craft in the correct direction
#if (false)
      c[k,0] = 1.0; // thrust[0].x
      c[k,N*3] = att0[0]; // =att0.x
      ct[k] = 0;
      c[k+1,1] = 1.0; // thrust[1].y
      c[k+1,N*3] = att0[1]; // =att0.x
      ct[k+1] = 0;
      c[k+2,2] = 1.0; // thrust[2].z
      c[k+2,N*3] = att0[2]; // =att0.x
      ct[k+2] = 0;
      k+=3;
#endif

      // Constrain N intermediate positions to be within minimumDescentAngle
      for( int j=0; j<N; j++ )
      {
        double tX = T*((double)(j+1)/(N+2));
        // Get whole weight vector up to time t
        RVWeightsToTime(tX,T,N,dt,out wr,out wv);
       
        double vx = Math.Sin(minDescentAngle*Math.PI/180.0); 
        double vy = Math.Cos(minDescentAngle*Math.PI/180.0); 
        //vx = 1.0;
        //vy = 1.0;
        double [] V1 = new double [] {vx,vy,0}; // Normal vector of plane to be above
        double [] V2 = new double [] {-vx,vy,0}; // Normal vector of plane to be above
        double [] V3 = new double [] {0,vy,vx}; // Normal vector of plane to be above
        double [] V4 = new double [] {0,vy,-vx}; // Normal vector of plane to be above
        for(int i = 0 ; i < N ; i++)
        {
          if (i == N-1) // Use constraint for final descent angle
          {
            vx = Math.Sin(minFinalDescentAngle*Math.PI/180.0); 
            vy = Math.Cos(minFinalDescentAngle*Math.PI/180.0); 
            V1 = new double [] {vx,vy,0}; // Normal vector of plane to be above
            V2 = new double [] {-vx,vy,0}; // Normal vector of plane to be above
            V3 = new double [] {0,vy,vx}; // Normal vector of plane to be above
            V4 = new double [] {0,vy,-vx}; // Normal vector of plane to be above
          }
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
        c[k+0,N*3] = V1[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V1[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V1[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+0] = 1; // LHS > RHS
        c[k+1,N*3] = V2[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V2[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V2[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+1] = 1; // LHS > RHS
        c[k+2,N*3] = V3[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V3[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V3[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+2] = 1; // LHS > RHS
        c[k+3,N*3] = V4[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V4[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V4[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+3] = 1; // LHS > RHS
        k += 4;
      }

      // constraint thrust to be downwards > amin
#if (true)
      for( int i=0; i<N; i++)
      {
        c[k,i*3+1] = 1.0; // weight on T[i].y
        c[k,N*3] = amin; // RHS
        ct[k] = 1; // LHS > RHS
        k += 1;
      }
#endif

      // zeroes for equality constraints
      alglib.minqpsetlc(state, c, ct);

      double[] s = new double[N*3];
      for(int i=0;i<N*3;i++)
      {
        s[i] = 1;
      }
      alglib.minqpsetscale(state, s);

      alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
      alglib.minqpoptimize(state);
      alglib.minqpresults(state, out x, out rep);
      double fuel=0;
      for(int i=0;i<N;i++)
      {
        double tx = x[i*3+0];
        double ty = x[i*3+1];
        double tz = x[i*3+2];
        fuel = fuel + Math.Sqrt(tx*tx+ty*ty+tz*tz)*dt;
      }

      retval = rep.terminationtype;
      if (retval>0)
      {
        thrusts = new double[N,3];
        for(int i=0;i<N;i++)
        {
          thrusts[i,0] = x[i*3+0];
          thrusts[i,1] = x[i*3+1];
          thrusts[i,2] = x[i*3+2];
        }
        System.Console.WriteLine("T={0:F2} FUEL={1:F2}", T, fuel);
      }
      else
      {
        fuel = 9e+20;
      }

      return fuel;
  }

   static double[] convToDouble3(Vector3d v)
   {
       return new double[] { v.x, v.y, v.z };
   }

   // If retval>0 than success. If retval<0 then various kind of failure. See https://www.alglib.net/translator/man/manual.csharp.html#minqpreportclass
   public static double gfold(Vector3d r0, Vector3d v0, Vector3d att0, Vector3d rf, Vector3d vf, double T, int N, double g, double amax, double minDescentAngle, double minFinalDescentAngle, double amin, out int retval, out Vector3d [] thrusts)
   {
      double [] _r0 = convToDouble3(r0);
      double [] _v0 = convToDouble3(v0);
      double [] _att0 = convToDouble3(att0);
      double [] _rf = convToDouble3(rf);
      double [] _vf = convToDouble3(vf);
      double [,] _thrusts = new double[N,3];
      double fuel = gfold(_r0, _v0, _att0, _rf, _vf, T, N, g, amax, minDescentAngle, minFinalDescentAngle, amin, out _thrusts, out retval);
      thrusts = null;
      if (retval > 0)
      {
        thrusts = new Vector3d[N];
        for(int i=0; i<N; i++)
        {
          thrusts[i].x = _thrusts[i,0];
          thrusts[i].y = _thrusts[i,1];
          thrusts[i].z = _thrusts[i,2];
        }
      }
      return fuel;
   }

#if (UNITY)
  public static double golden_search_gfold(Vector3d r0, Vector3d v0, Vector3d att0, Vector3d rf, Vector3d vf, double Tmin, double Tmax, int N, double g, Transform transform, double amax, double minDescentAngle, double minFinalDescentAngle, double amin, double tol, out int retval, out Vector3d[] thrusts)
  {
    if (transform != null)
    {
      r0 = transform.InverseTransformPoint(r0);
      v0 = transform.InverseTransformVector(v0);
      att0 = transform.InverseTransformVector(att0);
      rf = transform.InverseTransformPoint(rf);
      vf = transform.InverseTransformVector(vf);
    }
    double fuel = golden_search_gfold(r0, v0, att0, rf, vf, Tmin, Tmax, N, g, amax, minDescentAngle, minFinalDescentAngle, amin, tol, out retval, out thrusts);
    if (transform)
    {
      for(int i=0; i<N; i++)
      {
        thrusts[i] = transform.TransformVector(thrusts[i]);
      }
    }
    return fuel;
  }
#endif

  public static double golden_search_gfold(Vector3d r0, Vector3d v0, Vector3d att0, Vector3d rf, Vector3d vf, double Tmin, double Tmax, int N, double g, double amax, double minDescentAngle, double minFinalDescentAngle, double amin, double tol, out int retval, out Vector3d[] thrusts)
    {

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

      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      while (Math.Abs(c - d) > tol)
      {
        fc=gfold(r0,v0,att0,rf,vf,c,N,g,amax,minDescentAngle,minFinalDescentAngle,amin,out retval,out thrusts);
        fd=gfold(r0,v0,att0,rf,vf,d,N,g,amax,minDescentAngle,minFinalDescentAngle,amin,out retval,out thrusts);
        if (fc < fd)
            {b = d;}
        else
            {a = c;}

        // We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
      }
      gfold(r0,v0,att0,rf,vf,(b+a)/2,N,g,amax,minDescentAngle,minFinalDescentAngle,amin,out retval,out thrusts);
      return (b + a) / 2;
    }

    static int Main(string[] args)
    {
      // Initial position
      Vector3d r0 = new Vector3d(500,600,0);
      // Initial velocity
      Vector3d v0 = new Vector3d(-20,10,0);
      // Initial attitude
      Vector3d att0 = new Vector3d(0,1,0);
      // Final position
      Vector3d rf = new Vector3d(0,0,0);
      // Final velocity
      Vector3d vf = new Vector3d(0,0,0);
      // gravity
      Vector3d g = new Vector3d(0,-9.8,0);
      // maximum acceleration
      double amax = 30;
      // Number of discrete steps
      int N = 5;
      Vector3d [] thrusts;
      double Tmin = 2;
      double Tmax = 30;
      double tol = 0.5; // find best duration down to this tolerance
      double minDescentAngle = 40;
      double minFinalDescentAngle = 40;
      double amin = 1.0;
      int retval;
      double bestT = golden_search_gfold(r0,v0,att0,rf,vf,Tmin,Tmax,N,g.magnitude,amax,minDescentAngle,minFinalDescentAngle,amin,tol,out retval,out thrusts);
      if (retval > 0)
      {
        double dt = 0.5;
        System.Console.WriteLine("bestT="+bestT);
        Trajectory traj = new Trajectory();
        traj.Simulate(bestT,thrusts,r0,v0,g,dt);
        traj.CorrectFinal(rf, vf);

        System.IO.StreamWriter _slnWriter = new System.IO.StreamWriter("solution.dat");
        double t = 0;
        _slnWriter.WriteLine("time x y z vx vy vz ax ay az");
        for( int i = 0 ; i < traj.Length() ; i++ )
        {
          _slnWriter.WriteLine(t+" "+traj.r[i].x+" "+traj.r[i].y+" "+traj.r[i].z+" "+traj.v[i].x+" "+traj.v[i].y+" "+traj.v[i].z+" "+traj.a[i].x+" "+traj.a[i].y+" "+traj.a[i].z);
          t += dt;
        }
        _slnWriter.Close();
        System.Console.WriteLine("Written solution.dat. View with plotXYZ.py solution.dat");

        traj.FindClosest(new Vector3d(0,500,0), new Vector3d(0,-80,0), out Vector3d close_r, out Vector3d close_v, out Vector3d close_a, out double close_t, 1, 1);
        OrthogonalVectors(new Vector3d(10,-10,5),out Vector3d a,out Vector3d b);
        OrthogonalVectors(new Vector3d(0,0,10),out a,out b);
      }
      else
      {
        System.Console.WriteLine("Failed to find a solution: error code "+retval);
      }

      return(0);
    }
  }
}
