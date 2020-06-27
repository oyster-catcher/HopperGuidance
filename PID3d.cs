namespace HopperGuidance
{
  public class PID3d
  {
    PIDclamp pid1x,pid1y,pid1z;
    PIDclamp pid2x,pid2y,pid2z;

    public void Init(float kp1, float ki1, float kd1, float kp2, float ki2, float kd2, float vmax, float amax)
    {
      // Probably set I1=0, D1=0, I2=0, D2=0
      // P1 determines proportion of position error to velocity 1, to close 1m as 1m/s
      // P2 determines maximum spare acceleration available, 2+ times higher than P2 to avoid little oscillation
      // I1>0 would help overcome something like a constant wind or gravity (but can compensate exactly for gravity)
      pid1x = new PIDclamp("pid1x",kp1,ki1,kd1,vmax);
      pid1y = new PIDclamp("pid1y",kp1,ki1,kd1,vmax);
      pid1z = new PIDclamp("pid1y",kp1,ki1,kd1,vmax);
      pid2x = new PIDclamp("pid2x",kp2,ki2,kd2,amax);
      pid2y = new PIDclamp("pid2y",kp2,ki2,kd2,amax);
      pid2z = new PIDclamp("pid2y",kp2,ki2,kd2,amax);
      //Debug.Log("kp1="+kp1+" ki1="+ki1+" kd1="+kd1+" kp2="+kp2+" ki2="+ki2+" kd2="+kd2);
    }

    public void Reset()
    {
      pid1x.Reset();
      pid1y.Reset();
      pid1z.Reset();
    }

    public Vector3d Update(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv, float dt)
    {
      Vector3d tgt_v1;

      // Errors
      Vector3d err_r = dr - r;
      Vector3d err_v = dv - v;

      // Compute modified velocity targets given position and velocity error
      tgt_v1 = new Vector3d(pid1x.Update(err_r.x + err_v.x, dt),
                            pid1y.Update(err_r.y + err_v.y, dt),
                            pid1z.Update(err_r.z + err_v.z, dt));


      // Execute PID 2's to give acceleration
      Vector3d a = new Vector3d(pid2x.Update(tgt_v1.x, dt),
                                pid2y.Update(tgt_v1.y, dt),
                                pid2z.Update(tgt_v1.z, dt));

      //Debug.Log("tgt_v1="+tg1_v1+ "a="+a);
      return a;
    }
  }
}