using System;
using System.IO;
using System.Collections.Generic;

using UnityEngine;
using System.Collections;
using System.Numerics;

public class Vessel
{
  // Current attitude
  System.Numerics.Quaternion attitude = System.Numerics.Quaternion.Identity;
  // Current angular rotation
  System.Numerics.Quaternion rot_rate = System.Numerics.Quaternion.Identity;

  // Steer attitude to this vector
  public class Steer(Vector3d v)
  {


  }

  public class Update(float dt)
  {
    attitude = attitude * System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, rot_rate,dt);
  }
}

public class Spin
{
  public static int Main(string[] args)
  {
    float kD = 0.25f; // Turn speed is 25% of difference from attitude to steer direction each second
    //var pid = new PIDclamp("pid",1,0,0,100);
    System.Numerics.Quaternion attitude = System.Numerics.Quaternion.CreateFromYawPitchRoll(0,0,0);
    //System.Numerics.Quaternion steer = System.Numerics.Quaternion.CreateFromYawPitchRoll(90,0,0);
    System.Numerics.Quaternion steer = System.Numerics.Quaternion.FromVectorToVector(new Vector3(0,1,0),new Vector3(2,2,0));
    System.Numerics.Quaternion rot_rate = System.Numerics.Quaternion.Identity;
    //Vector3 v = Vector3.zero;
    System.Numerics.Quaternion attitude_error;
    System.Numerics.Quaternion desired_rot_rate;
    System.Numerics.Quaternion error_rot_rate;
    System.Numerics.Quaternion applied_rot_rate_change;
    float torque = 0.8f; // sort off, amount of rotation / sec as proportion of error
    float dt = 0.1f;

    System.Console.WriteLine("time pitch yaw roll desired_pitch_rate desired_yaw_rate desired_roll_rate pitch_rate yaw_rate roll_rate");
    float t = 0;
    //System.Console.WriteLine("attitude="+attitude.ToEulerAngles());
    while(t<=50)
    {
      // Calculate desired rot speed to turn to steer direction
      attitude_error = System.Numerics.Quaternion.Inverse(attitude) * steer;
      desired_rot_rate = System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, attitude_error, kD);

      // Calculate difference in current rotation to desired rotation
      error_rot_rate = System.Numerics.Quaternion.Inverse(rot_rate) * desired_rot_rate;

      // TODO: Limit by max rotational force
      applied_rot_rate_change = System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, error_rot_rate, torque);

      // Modify rotation speed by applied rotation
      rot_rate = rot_rate * System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, applied_rot_rate_change,dt);


      System.Console.WriteLine("{0:F1} {1:F1} {2:F1} {3:F2} {4:F2} {5:F2} {6:F2} {7:F2} {8:F2} {9:F2}",t,attitude.ToEulerAngles().x,attitude.ToEulerAngles().y,attitude.ToEulerAngles().z,desired_rot_rate.ToEulerAngles().x,desired_rot_rate.ToEulerAngles().y,desired_rot_rate.ToEulerAngles().z,rot_rate.ToEulerAngles().x,rot_rate.ToEulerAngles().y,rot_rate.ToEulerAngles().z);

      // Change attitude by rotation rate
      attitude = attitude * System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, rot_rate,dt);
      t = t + dt;
    }

    return 0;
  }
}
