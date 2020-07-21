#define MINDESCENTANGLE
#define MAXTHRUSTANGLE
#define FINALTHRUSTDIR
#define MAXVELOCITY
#define MINTHRUST
//#define DUMP

using System;
using System.IO;
using KSPAssets;


namespace HopperGuidance
{
  public class Solve
  {
    // Parameters to control solution
    public double Tmin = 0;
    public double Tmax = 30;
    public int Nmin = 2;
    public int Nmax = 10;
    public double minDurationPerThrust = 4;
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

    // Last stored inputs to GFold() - in natural space for solution
    // with Y as the up direction
    public Vector3d r0;
    public Vector3d v0;
    public Vector3d rf;
    public Vector3d vf;

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
      return string.Format("HopperGuidance: "+msg+" Nmin="+Nmin+ " Nmax="+Nmax+ " minDurationPerThrust="+minDurationPerThrust+" N="+N+" r0="+Vec2Str(r0)+" v0="+Vec2Str(v0)+" rf="+Vec2Str(rf)+" vf="+Vec2Str(vf)+" g="+g+" Tmin="+Tmin+" Tmax="+Tmax+" amax="+amax+" vmax="+vmax+" minDescentAngle="+minDescentAngle+" maxThrustAngle="+maxThrustAngle+" maxLandingThrustAngle="+maxLandingThrustAngle);
    }
 
    public static double [] BasisWeights(double t, double a_T, int N)
    {
      // Returns vector of weights[N] for each acceleration vector at a given time, t
      // for time t in range 0...T, with N vectors
      // T divided into N-1 parts, 0, T/(N-1), 2T/(N-1), (N-1)T/(N-1) == T
      double [] w = new double[N];
      double dt = a_T/(N-1);   // distance between basis points in time
      for(int j=0; j<N; j++)
      {
        double d = j - (t/dt);
        double b = 0; // when outside local range
        if ((d>-1) && (d<1))
          b = Math.Cos(d*0.5*Math.PI);
        w[j] = b*b;
      }
      return w;
    }

    // Calculate weights on position and velocity from thrust vectors up to time tX
    static void RVWeightsToTime(double tX, double T,int N, double dt, out double[] wr, out double[] wv)
    {
      wr = new double[N];
      wv = new double[N];
      for(double t = 0; t < tX; t += dt)
      {
        double tr = tX - t - dt; // remaining time after applying acceleration for dt
        double [] w = BasisWeights(t,T,N); // Vector for all N weights at time, t
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
    public double GFold(double[] r0, double[] v0, double[] rf, double[] vf, double T,
                        out double[,] a_thrusts, out int a_retval)
    {
      double fidelity = 10; // this many steps inbetween thrust positions
      a_thrusts = null;

      N = (int)(T/minDurationPerThrust);
      if (N < Nmin)
        N = Nmin;
      if (N > Nmax)
        N = Nmax;
      int numchecks = (int)(0.5*T); // number of checks for descent angle

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
        bndl[N*3+i] = 0;
        bndu[N*3+i] = amax;
      }

      alglib.minqpsetbc(state, bndl, bndu);

      // 1 constraint on final X position
      // 1 constraint on final Y position
      // 1 constraint on final X vel.
      // 1 constraint on final Y vel.
      // 1 constraint for each Ty[i]>0

      int constraints = 3 + 3 + 6*N;

      // for minDescentAngle, N points for 4 planes to make square 'cones'
#if (MINDESCENTANGLE)
      constraints += 4*numchecks;
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

      int rhs = N*4;
      double [,] c = new double[constraints,rhs+1]; // zeroed?
      int [] ct = new int[constraints]; // type of constraint, =, > or <  (default to 0 -> =)

      // Constrain thrust vectors to be below thrust magnitudes
      for( int i = 0 ; i < N ; i++ )
      {
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Mean thrust vector X axis less than thrust magnitude
        k++;
        c[k,i*3+0] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Mean thrust vector X axis greater than -thrust magnitude
        k++;

        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Mean thrust vector Y axis less than thrust magnitude
        k++;
        c[k,i*3+1] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Mean thrust vector Y axis greater than -thrust magnitude
        k++;

        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = -1.0;
        ct[k] = -1; // LHS < 0. Mean thrust vector Z axis less than thrust magnitude
        k++;
        c[k,i*3+2] = 1.0;
        c[k,N*3+i] = 1.0;
        ct[k] = 1; // LHS > 0. Mean thrust vector Z axis greater than -thrust magnitude
        k++;
      }

      for(double t = 0; t < T; t += dt)
      {
        double tr = T - t - dt; // remaining time after applying acceleration
        double [] w = BasisWeights(t,T,N); // Vector for all N weights at time, t
        for(int i = 0 ; i < N ; i++)
        {
          // constrain final r=[0,0,0]
          c[k+0,i*3+0] += tr*w[i]*dt + 0.5*w[i]*dt*dt;
          c[k+1,i*3+1] += tr*w[i]*dt + 0.5*w[i]*dt*dt;
          c[k+2,i*3+2] += tr*w[i]*dt + 0.5*w[i]*dt*dt;
          // constrain final v=[0,0,0]
          c[k+3,i*3+0] += w[i]*dt;
          c[k+4,i*3+1] += w[i]*dt;
          c[k+5,i*3+2] += w[i]*dt;
        }
      }
      for( int i = 0 ; i < N ; i++ )
      {
        // final r function sums to this
        c[k+0,rhs] = rf[0] - (r0[0] + v0[0]*T); // X constant
        c[k+1,rhs] = rf[1] - (r0[1] + v0[1]*T - 0.5*T*T*g); // Y constant
        c[k+2,rhs] = rf[2] - (r0[2] + v0[2]*T); // Z constant
        // final v function sums to this
        c[k+3,rhs] = vf[0] - v0[0];
        c[k+4,rhs] = vf[1] - (v0[1] - T*g);
        c[k+5,rhs] = vf[2] - v0[2];
      }
      k+=6;

#if (MINDESCENTANGLE)
      // Constrain N intermediate positions to be within minimumDescentAngle
      for( int j=0; j<numchecks; j++ )
      {
        // No check at t=0 and t=T
        double tX = (T*(j+1))/(numchecks+1);
        // Get whole weight vector up to time t
        RVWeightsToTime(tX,T,N,dt,out double[] wr,out double[] wv);

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
        c[k+0,rhs] = V1[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V1[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V1[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+0] = 1; // LHS > RHS
        c[k+1,rhs] = V2[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V2[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V2[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+1] = 1; // LHS > RHS
        c[k+2,rhs] = V3[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V3[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V3[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+2] = 1; // LHS > RHS
        c[k+3,rhs] = V4[0] * (rf[0] - (r0[0] + v0[0]*tX))
                   + V4[1] * (rf[1] - (r0[1] + v0[1]*tX - 0.5*tX*tX*g))
                   + V4[2] * (rf[2] - (r0[2] + v0[2]*tX)); // RHS
        ct[k+3] = 1; // LHS > RHS
        k += 4;
      }
#endif

#if (MAXVELOCITY)
      // Constrain N intermediate positions to be within minimumDescentAngle
      for( int j=0; j<N; j++ )
      {
        // No check at t=0 and t=T
        double tX = T*((double)(j+1))/(N+1);
        // Get whole weight vector up to time t
        RVWeightsToTime(tX,T,N,dt,out double[] wr,out double[] wv);

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
          c[k,N*3+i] = 1.0; // thrust weight
          c[k,rhs] = amin;
          ct[k] = 1; // LHS > RHS
          k++;
      }
#endif

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
      a_thrusts = new double[N,3];
      if ((a_retval>=1) && (a_retval<=5))
      {
        for(int i=0;i<N;i++)
        {
          a_thrusts[i,0] = x[i*3+0];
          a_thrusts[i,1] = x[i*3+1];
          a_thrusts[i,2] = x[i*3+2];
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

   // If retval>0 than success. If retval<0 then various kind of failure. See https://www.alglib.net/translator/man/manual.csharp.html#minqpreportclass
   public double GFold(Vector3d a_r0, Vector3d a_v0, Vector3d a_rf, Vector3d a_vf, double a_T,
                       out Vector3d [] a_thrusts, out int a_retval)
   {
      T = a_T;
      double [] _r0 = convToDouble3(a_r0);
      double [] _v0 = convToDouble3(a_v0);
      double [] _rf = convToDouble3(a_rf);
      double [] _vf = convToDouble3(a_vf);
      double [,] _thrusts = new double[N,3];
      double fuel = GFold(_r0, _v0, _rf, _vf, T, out _thrusts, out a_retval);
      a_thrusts = null;
      if (retval > 0)
      {
        a_thrusts = new Vector3d[N];
        for(int i=0; i<N; i++)
        {
          a_thrusts[i].x = _thrusts[i,0];
          a_thrusts[i].y = _thrusts[i,1];
          a_thrusts[i].z = _thrusts[i,2];
        }
      }
      return fuel;
   }


    public double GoldenSearchGFold(Vector3d a_r0, Vector3d a_v0, Vector3d a_rf, Vector3d a_vf,
                                    out Vector3d [] o_thrusts, out double o_fuel, out int o_retval)
    {
      // Store best solution values
      r0 = Vector3d.zero;
      v0 = Vector3d.zero;
      rf = Vector3d.zero;
      vf = Vector3d.zero;
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
      int retvalc=-1,retvald=-1;

      o_thrusts = new Vector3d[N];

      c = b - (b - a) / gr;
      d = a + (b - a) / gr;
      while (Math.Abs(c - d) > tol)
      {
        fc = GFold(a_r0,a_v0,a_rf,a_vf,c,out o_thrusts,out retvalc);
        fd = GFold(a_r0,a_v0,a_rf,a_vf,d,out o_thrusts,out retvald);
        if (fc < fd)
            {b = d;}
        else
            {a = c;}

        // We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
      }

      double bestT = 0.5*(a+b);
      o_fuel = GFold(a_r0,a_v0,a_rf,a_vf,bestT,out o_thrusts,out o_retval);
      retval = o_retval;
      if ((o_retval<1) || (o_retval>5))
      {
        System.Console.Error.WriteLine("FAILED AT T={0}",0.5*(a+b));
      }
      else
      {
        r0 = a_r0;
        v0 = a_v0;
        rf = a_rf;
        vf = a_vf;
      }
      
      return bestT;
    }

    static int RunTest(string[] args)
    {
      Solve solver = new Solve();
      Vector3d r0 = Vector3d.zero;
      Vector3d v0 = Vector3d.zero;
      Vector3d rf = Vector3d.zero;
      Vector3d vf = Vector3d.zero;
      Vector3d c;
      double d;
      for(int i=0;i<args.Length;i++)
      {
        string k = args[i].Split('=')[0];
        string v = args[i].Split('=')[1];
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
            rf = c;
          else if (k=="vf")
            vf = c;
        } else {
          d = Convert.ToDouble(v);
          System.Console.Error.WriteLine(k+"="+d);
          if (k=="Nmin")
            solver.Nmin = (int)d;
          else if (k=="Nmax")
            solver.Nmax = (int)d;
          else if (k=="N")
          {
            solver.Nmin = (int)d;
            solver.Nmax = (int)d;
          }
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
          else
          {
            System.Console.Error.WriteLine("No such parameter: {0}",k);
            return(1);
          }
        }
      }
      // Now run Solver
      double fuel;
      int retval;
      Vector3d [] thrusts;
      double bestT = solver.GoldenSearchGFold(r0,v0,rf,vf,out thrusts,out fuel,out retval);

      System.Console.Error.WriteLine(solver.DumpString());
      if ((retval>=1) && (retval<=5))
      {
        double final_r_err = 0;
        Trajectory traj = new Trajectory();
        Vector3d vg = new Vector3d(0,-solver.g,0);
        traj.Simulate(bestT, thrusts, r0, v0, vg, solver.dt, 0);
        traj.CorrectFinal(Vector3d.zero,Vector3d.zero);
        traj.Write(null);
        final_r_err = (traj.r[traj.r.Length-1] - rf).magnitude;
      }
      return(0);
    }

    static int Main(string[] args)
    {
      if (args.Length == 0) {
        System.Console.Error.WriteLine("usage: Solve.exe k=v ... - set of key=value pairs from r0,v0,rf,vf,Tmin,Tmax,amin,amax,g,minDescentAngle,tol,vmax which also have defaults");
        return(1);
      } else {
        return RunTest(args);
      }
    }
  }
}
