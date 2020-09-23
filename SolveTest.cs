using System;
using System.IO;
using System.Collections.Generic;

namespace HopperGuidance
{
  public class SolveTest
  {
    static void ListArgs(List<Tuple<string,string>> args) {
      foreach(var arg in args) {
        if (arg.Item2=="")
          System.Console.Error.WriteLine(arg.Item1);
        else
          System.Console.Error.WriteLine(arg.Item1+"="+arg.Item2);
      }
    }

    static bool Set(string key, string value, ref Vector3d r0, ref Vector3d v0, ref Vector3d rf, ref Vector3d vf, ref int rfaxes, ref int vfaxes, ref List<SolveTarget> targets, ref bool correct) {
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

    static void SetFinal(Vector3d rf, Vector3d vf, int rfaxes, int vfaxes,ref List<SolveTarget> targets) {
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

    // MAIN PROGRAM
    static int Main(string[] mainargs)
    {
      if (mainargs.Length == 0) {
        System.Console.Error.WriteLine("usage: Solve.exe [--correct] k=v ... - set of key=value pairs from r0,v0,rf,vf,Tmin,Tmax,amin,amax,g,minDescentAngle,tol,vmax,ir,iv which also have defaults (ir,iv are intermediates which be specified multiple times)");
        return(1);
      } else {
        Solve solver = new Solve();
        SolveResult result = new SolveResult();;
        Trajectory traj = new Trajectory();
        List<SolveTarget> targets = new List<SolveTarget>();
        Vector3d r0 = Vector3d.zero;
        Vector3d v0 = Vector3d.zero;
        Vector3d rf = Vector3d.zero;
        Vector3d vf = Vector3d.zero;
        int rfaxes = 0;
        int vfaxes = 0;
        bool correct = false;

        // Argument processing - arguments use Set() on different objects they might apply to
        foreach(var arg in HGUtils.ToKeyValuePairs(mainargs)) {
          bool used = Set(arg.Item1, arg.Item2, ref r0, ref v0, ref rf, ref vf, ref rfaxes, ref vfaxes, ref targets, ref correct);
          used = used || solver.Set(arg.Item1, arg.Item2);
          if (!used) {
            System.Console.Error.WriteLine("Unexpected argument: "+arg.Item1+"="+arg.Item2);
            return 1;
          }
        }
        SetFinal(rf, vf, rfaxes, vfaxes, ref targets);

        if (targets.Count == 0) {
          System.Console.Error.WriteLine("Error! You must specify targets with either target=[x,y,z], rf=[x,y,z] or vf=[x,y,z]");
          return 1;
        }
        result = MainProg.MultiPartSolve(ref solver, ref traj, r0, v0, ref targets, (float)solver.g, solver.extendTime, correct);

        if (result.isSolved()) {
          List<string> comments = new List<string>();
          comments.Add(solver.DumpString());
          comments.Add(result.DumpString());
          // Thrusts
          List<float> thrust_times = new List<float>();
          for( int i=0; i<result.thrusts.Length; i++)
            thrust_times.Add(result.thrusts[i].t);
          comments.Add("thrust_times="+String.Join(",",thrust_times));
          if( result.checktimes != null )
            comments.Add("check_times="+String.Join(",",result.checktimes));
          List<double> sln_Ts = new List<double>();
          List<double> sln_fuels = new List<double>();
          foreach(var sln in result.solutions) {
            sln_Ts.Add(sln.T);
            sln_fuels.Add(sln.fuel);
          }
          comments.Add("solution_time="+String.Join(",",sln_Ts));
          comments.Add("solution_fuel="+String.Join(",",sln_fuels));
      
          traj.Write(null, comments);
          System.Console.Error.WriteLine(result.DumpString());
          return(0);
        }
        return(1);
      }
    }
  }  
}
