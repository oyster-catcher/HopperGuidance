using System;
using System.IO;
using System.Collections.Generic;

namespace HopperGuidance
{
  public class VesselSim
  {
    static void Simulate(Controller controller, ref Trajectory traj, Vector3d r, Vector3d v, Vector3d g, double dt, double maxT, string tgtLogFilename="target.dat")
    {
      double t = 0;
      double throttle;
      Vector3d thrustV = Vector3d.zero;
      Vector3d dr,dv,da; // desired values on trajectory

      var tgtWriter = new System.IO.StreamWriter(tgtLogFilename);
      tgtWriter.WriteLine("time x y z vx vy vz ax ay az att_err");

      Vector3d up = new Vector3d(0,1,0);
      System.Console.WriteLine("time x y z vx vy vz ax ay az att_err");
      while(( t < maxT ) && (r.y > 0))
      {
        System.Console.WriteLine("{0:F2} {1:F2} {2:F2} {3:F2} {4:F2} {5:F2} {6:F2} {7:F2} {8:F2} {9:F2} 0",t,r.x,r.y,r.z,v.x,v.y,v.z,thrustV.x, thrustV.y, thrustV.z);
        // Equations of motion
        r = r + v*dt;
        v = v + g*dt;

        // Control loop
        float att_err;
        controller.GetControlOutputs(traj, r, v, g, (float)dt, out dr, out dv, out da, out throttle, out thrustV,out att_err);
        tgtWriter.WriteLine("{0:F2} {1:F2} {2:F2} {3:F2} {4:F2} {5:F2} {6:F2} {7:F2} {8:F2} {9:F2} {10:F2}",t,dr.x,dr.y,dr.z,dv.x,dv.y,dv.z,da.x,da.y,da.z,att_err);
        //System.Console.WriteLine("throttle="+throttle+" thrustV="+thrustV);
        //
        // Make engine thrust have effect
        v = v + thrustV*dt;

        t = t + dt;
      }
      tgtWriter.Close();
    }

       static bool SetTargets(string key, string value, ref Vector3d r0, ref Vector3d v0, ref Vector3d rf, ref Vector3d vf, ref int rfaxes, ref int vfaxes, ref List<SolveTarget> targets, ref bool correct) {
      float t = -1;

      if (key=="--correct")
        correct = true;
      else if (key=="r0")
        r0 = HGUtils.ToVector3d(value);
      else if (key=="v0")
        v0 = HGUtils.ToVector3d(value);
      else if (key=="rf") {
        rf = HGUtils.ToVector3d(value);
        rfaxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
      }
      else if (key=="vf") {
        vf = HGUtils.ToVector3d(value);
        vfaxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
      }
      // TODO: Handle : to specify time
      else if (key=="target")
      {
        SolveTarget tgt = new SolveTarget();
        tgt.r = HGUtils.ToVector3d(value);
        tgt.raxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
        tgt.vaxes = 0;
        tgt.t = t;
        targets.Add(tgt);
      }
      else
        return false;
      return true;
    }

    static void SetFinalTarget(Vector3d rf, Vector3d vf, int rfaxes, int vfaxes,ref List<SolveTarget> targets) {
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
    }

    static bool Set(string k, string v, ref double dt, ref double maxT)
    {
      if (k=="dt") {
        dt = Convert.ToDouble(v);
      } else if (k=="maxT") {
        maxT = Convert.ToDouble(v);
      } else
        return false;
      return true;
    }

    static int Main(string[] mainargs)
    {
      if (mainargs.Length == 0)
      {
        System.Console.Error.WriteLine("usage: VesselSim.exe  - computes a solution given initial conditions like SolveText.exe but also simulates a vessel like running inside KSP");
        System.Console.Error.WriteLine("vessel options:  drag=<float>");
        System.Console.Error.WriteLine("                 start-up-time=<float>");
        System.Console.Error.WriteLine("                 dt=<float>");
        return 1;
      }
      var solver = new Solve();
      var result = new SolveResult();
      var traj = new Trajectory();
      var controller = new Controller();
      bool correct = false;
      Vector3d r0 = Vector3d.zero;
      Vector3d v0 = Vector3d.zero;
      Vector3d rf = Vector3d.zero;
      Vector3d vf = Vector3d.zero;
      int rfaxes = 0;
      int vfaxes = 0;
      double dt = 0.1;
      double maxT = 0;
      List<SolveTarget> targets = new List<SolveTarget>();

      foreach(var arg in HGUtils.ToKeyValuePairs(mainargs))
      {
        bool used = solver.Set(arg.Item1, arg.Item2);
        used = used || controller.Set(arg.Item1, arg.Item2);
        used = used || SetTargets(arg.Item1, arg.Item2, ref r0, ref v0, ref rf, ref vf, ref rfaxes, ref vfaxes, ref targets, ref correct);
        used = used || Set(arg.Item1, arg.Item2, ref dt, ref maxT);
        if (!used) {
          throw new System.ArgumentException("Unexpected argument: "+arg.Item1+"="+arg.Item2);
        }
      }
      SetFinalTarget(rf, vf, rfaxes, vfaxes, ref targets);

      if (targets.Count == 0) {
        System.Console.Error.WriteLine("Error! You must specify targets with either target=[x,y,z], rf=[x,y,z] or vf=[x,y,z]");
        return 1;
      }

      result = MainProg.MultiPartSolve(ref solver, ref traj, r0, v0, ref targets, (float)solver.g, solver.extendTime, correct);
      if (result.isSolved()) {
        List<string> comments = new List<string>();
        System.Console.Error.WriteLine("Writing solution to: solution.dat");
        traj.Write("solution.dat", comments);
        if (traj.T > maxT)
          maxT = traj.T;
        Simulate(controller, ref traj, traj.r[0], traj.v[0], new Vector3d(0,-solver.g,0), dt, maxT);
      }
      return 0;
    }
  }
}
