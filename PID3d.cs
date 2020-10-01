//using UnityEngine;

namespace HopperGuidance
{
  public class PID3d
  {
    public float kp1;
    public float ki1;
    public float kd1;
    public float kp2;
    public float ki2;
    public float kd2;
    public float vmax;
    public float amax;
    public float ymult;
    PIDclamp pid1x,pid1y,pid1z;
    PIDclamp pid2x,pid2y,pid2z;
    float last_t;
    System.IO.StreamWriter f = null;

    public PID3d(float a_kp1=1, float a_ki1=0, float a_kd1=0,
                 float a_kp2=1, float a_ki2=0, float a_kd2=0,
                 float a_vmax=1000, float a_amax=100, float a_ymult=1.0f,
                 string logFilename="")
    {
      // Probably set I1=0, D1=0, I2=0, D2=0
      // P1 determines proportion of position error to velocity 1, to close 1m as 1m/s
      // P2 determines maximum spare acceleration available, 2+ times higher than P2 to avoid little oscillation
      // I1>0 would help overcome something like a constant wind or gravity (but can compensate exactly for gravity)
      // ymult scales up P1 for make the PID focus on hitting the height target
      // (this is easier since we have gravity to pull us down if overshooting upward)
      kp1 = a_kp1;
      ki1 = a_ki1;
      kd1 = a_kd1;
      kp2 = a_kp2;
      ki2 = a_ki2;
      kd2 = a_kd2;
      vmax = a_vmax;
      amax = a_amax;
      ymult = a_ymult;
      last_t = 0;
      SetLogFilename(logFilename);
      Reset();
    }

    public void SetLogFilename(string logFilename)
    {
      if (f != null)
        f.Close();
      if (logFilename != "")
      {
        f = new System.IO.StreamWriter(logFilename);
        f.WriteLine("time x_err y_err z_err vx_err vy_err vz_err ax ay az");
      }
    }

    public void Reset()
    {
      pid1x = new PIDclamp("pid1x",kp1,ki1,kd1,vmax);
      pid1y = new PIDclamp("pid1y",kp1*ymult,ki1,kd1,vmax);
      pid1z = new PIDclamp("pid1z",kp1,ki1,kd1,vmax);
      // Scale up amax so its effectively ignored
      pid2x = new PIDclamp("pid2x",kp2,ki2,kd2,amax);
      pid2y = new PIDclamp("pid2y",kp2*ymult,ki2,kd2,amax);
      pid2z = new PIDclamp("pid2z",kp2,ki2,kd2,amax);
    }

    public Vector3d Update(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv, float t)
    {
      Vector3d tgt_v1;
      float dt = t - last_t;
      last_t = t;

      // Errors
      Vector3d err_r = dr - r;
      Vector3d err_v = dv - v;

      // Compute modified velocity targets given position and velocity error
      tgt_v1 = new Vector3d(pid1x.Update(err_r.x, dt),
                            pid1y.Update(err_r.y, dt),
                            pid1z.Update(err_r.z, dt));

      tgt_v1 = tgt_v1 + err_v;

      // Execute PID 2's to give acceleration
      Vector3d a = new Vector3d(pid2x.Update(tgt_v1.x, dt),
                                pid2y.Update(tgt_v1.y, dt),
                                pid2z.Update(tgt_v1.z, dt));

      if (f != null)
        f.WriteLine("{0:F3} {1:F2} {2:F2} {3:F2} {4:F2} {5:F2} {6:F2} {7:F2} {8:F2} {9:F2}",t,err_r.x,err_r.y,err_r.z,err_v.x,err_v.y,err_v.z,a.x,a.y,a.z);

      return a;
    }
  }
}
