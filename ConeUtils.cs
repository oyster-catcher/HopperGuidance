using System;
using UnityEngine;

namespace HopperGuidance
{
  public class ConeUtils
  {
    public static Vector2 FindClosestPointOnLine(Vector2 P, Vector2 A, Vector2 B, out float dist)
    {
      Vector2 a_to_p = P - A;
      Vector2 a_to_b = B - A;
      float atb2 = a_to_b.x*a_to_b.x + a_to_b.y*a_to_b.y;
      float atp_dot_atb = Vector2.Dot(a_to_b,a_to_p);
      float t = atp_dot_atb / atb2; // The normalized "distance" from a to your closest point
      if (t<0) {t=0;}
      if (t>1) {t=1;}
      Vector2 Z = A + a_to_b*t;
      dist = (P-Z).magnitude;
      return Z;
    }


    public static bool RightOfLine(Vector2 P, Vector2 L1, Vector2 L2)
    {
      // Normal is line vector rotated by 90 degrees to right
      Vector2 N = new Vector2(L2.y - L1.y, - (L2.x - L1.x));

      return Vector2.Dot(P-L1,N) > 0;
    }

    // shape should be described clockwise
    public static bool InsideConvexPolygon(Vector2 P, Vector2 [] shape)
    {
      bool inside = true;
      //Console.WriteLine("P="+P);
      for(int i=0; i<shape.Length-1; i++)
      {
        //Console.WriteLine(shape[i]+" to "+shape[i+1]+" "+RightOfLine(P,shape[i],shape[i+1]));
        inside = inside & RightOfLine(P,shape[i],shape[i+1]);
      }
      //Console.WriteLine(shape[shape.Length-1]+" to "+shape[0]+" "+RightOfLine(P,shape[shape.Length-1],shape[0]));
      return inside & RightOfLine(P,shape[shape.Length-1],shape[0]);
    }

    public static Vector2 FindClosestPointOnConvexPolygon(Vector2 P, Vector2 [] shape)
    {
      if (InsideConvexPolygon(P,shape))
      {
        //Console.WriteLine("inside");
        return P;
      }
      //Console.WriteLine("outside");
      // Find closest point on boundary 
      float dmin = 1.0e+10f;
      Vector2 pmin = Vector2.zero;
      Vector2 p = Vector2.zero;
      float d;
      for(int i=0; i<shape.Length-1; i++)
      {
        p = FindClosestPointOnLine(P, shape[i], shape[i+1], out d);
        if (d < dmin)
          { dmin = d; pmin = p; }
      }
      p = FindClosestPointOnLine(P, shape[shape.Length-1], shape[0], out d);
      if (d < dmin)
        { dmin = d; pmin = p; }
      return pmin;
    }


    public static Vector2 FindClosestPointOnLines(Vector2 P, Vector2 [] shape)
    {
      // Find closest point on boundary 
      float dmin = 1.0e+10f;
      Vector2 pmin = Vector2.zero;
      Vector2 p = Vector2.zero;
      float d;
      for(int i=0; i<shape.Length-1; i++)
      {
        p = FindClosestPointOnLine(P, shape[i], shape[i+1], out d);
        if (d < dmin)
          { dmin = d; pmin = p; }
      }
      return pmin;
    }


    public static Vector3d ClosestThrustInsideCone(float maxAngle, float amin, float amax, Vector3 F)
    {
      if (maxAngle>90)
        return F;

      Vector3 N = new Vector3(F.x,0,F.z); // plane of vector F
      N = Vector3.Normalize(N);

      //Console.WriteLine("N="+N);

      maxAngle = maxAngle*Mathf.PI/180;
      // Define boundaries of allow thrusts in a truncated approximate cone (without min thrust cone)
      //    ..---- P3 ----..
      //   P2              P4
      //    \             /
      //     \           /
      //      \         /
      //           P0 
      // and add some extra points!

      // Need separate thrust cones to keep them convex

      // Thrust cone thrust 0 to amax
      Vector2 [] conemax = {new Vector2(0,0),
                      new Vector2(Mathf.Sin(-maxAngle)*amax,Mathf.Cos(-maxAngle)*amax),
                      new Vector2(Mathf.Sin(-maxAngle*0.5f)*amax,Mathf.Cos(-maxAngle*0.5f)*amax),
                      new Vector2(0,amax),
                      new Vector2(Mathf.Sin(maxAngle*0.5f)*amax,Mathf.Cos(maxAngle*0.5f)*amax),
                      new Vector2(Mathf.Sin(maxAngle)*amax,Mathf.Cos(maxAngle)*amax)};

      // Put on 3D plane
      Vector2 F2 = new Vector2(Vector3.Dot(N,F),F.y);
      Console.WriteLine("F2="+F2);

      // Find closest ( can be inside)
      Vector2 pmin = FindClosestPointOnConvexPolygon(F2, conemax);
      Console.WriteLine("pmin="+pmin);

      // Now consider thrust cone 0 to amin. This is better check than polygon
      // cos we will be on edge which causes problems
      bool inside_min = (pmin.magnitude < amin);
      Console.WriteLine("inside_min="+inside_min);

      // Find closest point on boundary with at amin
      if (inside_min)
      {
        // The central point is lower than it should be. This is a hack to that when
        // the desired thrust in strongly downwards we don't use a nearest point on the
        // edge of the min thrust cone which is lower in Y but highly positive in X
        //Vector2 [] bound = {new Vector2(Mathf.Sin(-maxAngle)*amin,Mathf.Cos(-maxAngle)*amin),
        //                    new Vector2(0,amin*0.7f),
        //                    new Vector2(Mathf.Sin(maxAngle*0.5f)*amin,Mathf.Cos(maxAngle*0.5f)*amin)};
        //pmin = FindClosestPointOnLines(F2, bound);
        // Actually it seems to work better just to use a vector is the same direction
        // but with greater magnitude
        Console.WriteLine("mag="+pmin.magnitude);
        if (pmin.magnitude > 0)
          pmin = amin*(pmin/pmin.magnitude);
        else
        {
          pmin.x = 0;
          pmin.y = amin;
        }
      }

      //Console.WriteLine("pmin="+pmin);
      // Transform back into 3D
      F = pmin.x*N + (new Vector3(0,pmin.y,0));

      return F;
    }

    public static int Main(string[] args)
    {
      if (args.Length != 6) {
        Console.Error.WriteLine("usage: ConeUtils.exe <maxAngle> <amin> <amax> <Fx> <Fy> <Fz>");
        return 1;
      }
      float maxAngle = (float)Convert.ToDouble(args[0]);
      float amin = (float)Convert.ToDouble(args[1]);
      float amax = (float)Convert.ToDouble(args[2]);
      float x = (float)Convert.ToDouble(args[3]);
      float y = (float)Convert.ToDouble(args[4]);
      float z = (float)Convert.ToDouble(args[5]);

      Vector3 F = new Vector3(x,y,z);
      Vector3 Fp = ClosestThrustInsideCone(maxAngle,amin,amax,F);
      System.Console.WriteLine("F="+F+" F^="+Fp);

      return 0;
    } 
  }
}
