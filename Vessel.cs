using System;
using System.IO;
using System.Collections.Generic;

using UnityEngine;
using System.Collections;
using System.Numerics;

public class RotBody
{
  // Current attitude
  System.Numerics.Quaternion attitude = System.Numerics.Quaternion.Identity;
  // Current angular rotation
  System.Numerics.Quaternion rot_speed = System.Numerics.Quaternion.Identity;
  // Steer to Quaternion
  System.Numerics.Quaternion steer = System.Numerics.Quaternion.Identity;
  float last_t = 0;
  float kD = 0.5f; // speed to turn towards steer vector
  // Vessel parameters
  Vector3d rot_inertia = new Vector3d(2,2,2); // rotational inertia of pitch, yaw, roll
  float torque_gain = 5; // sensitivity of thrusters to attitude error
  float max_torque = 10; // max turning power of thrusters

  public static double Clamp(double val, double min, double max)
  {
    if (val<min)
      return min;
    if (val>max)
      return max;
    return val;
  }

  // Steer attitude to this vector
  public void Steer(Vector3d v)
  {
    steer = System.Numerics.Quaternion.FromVectorToVector(Vector3.up, v);
  }

  public void Update(float t, bool logging=false)
  {
    System.Numerics.Quaternion torque;
    Vector3 rot_accel_euler;

    float dt = t - last_t;
    last_t = t;

    // Steer to steer vector
    System.Numerics.Quaternion attitude_error;
    System.Numerics.Quaternion desired_rot_speed;
    System.Numerics.Quaternion error_rot_speed;
    System.Numerics.Quaternion rot_accel;

    // Calculate desired rot speed to turn to steer direction
    attitude_error = System.Numerics.Quaternion.Inverse(attitude) * steer;
    desired_rot_speed = System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, attitude_error, kD);

    // Calculate difference in current rotation to desired rotation
    error_rot_speed = System.Numerics.Quaternion.Inverse(rot_speed) * desired_rot_speed;

    // TODO: Limit by max rotational force
    torque = System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, error_rot_speed, torque_gain);
    Vector3d torque_euler = torque.ToEulerAngles();
    torque_euler.x = Clamp(torque_euler.x,-max_torque,max_torque);
    torque_euler.y = Clamp(torque_euler.y,-max_torque,max_torque);
    torque_euler.z = Clamp(torque_euler.z,-max_torque,max_torque);
    rot_accel_euler = new Vector3d(torque_euler.x/rot_inertia.x,torque_euler.y/rot_inertia.y,torque_euler.z/rot_inertia.z);
    rot_accel = System.Numerics.Quaternion.CreateFromYawPitchRoll((float)rot_accel_euler.x,(float)rot_accel_euler.z,(float)rot_accel_euler.y);

    // Modify rotation speed by applied rotation
    rot_speed = rot_speed * System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, rot_accel,dt);

    attitude = attitude * System.Numerics.Quaternion.Lerp(System.Numerics.Quaternion.Identity, rot_speed,dt);

    // Logging
    if (logging)
      System.Console.WriteLine("{0:F1} {1:F1} {2:F1} {3:F2} {4:F2} {5:F2} {6:F2} {7:F2} {8:F2} {9:F2}",t,attitude.ToEulerAngles().x,attitude.ToEulerAngles().y,attitude.ToEulerAngles().z,desired_rot_speed.ToEulerAngles().x,desired_rot_speed.ToEulerAngles().y,desired_rot_speed.ToEulerAngles().z,rot_speed.ToEulerAngles().x,rot_speed.ToEulerAngles().y,rot_speed.ToEulerAngles().z);
  }
}

public class Spin
{
  public static int Main(string[] args)
  {
    RotBody vessel = new RotBody();

    Vector3d steer = new Vector3(0.1f,1,0.5f);
    float t=0;
    float dt=0.2f;

    System.Console.WriteLine("time pitch yaw roll desired_pitch_rate desired_yaw_rate desired_roll_rate pitch_rate yaw_rate roll_rate");
    while(t<25)
    {
      // Change attitude by rotation rate
      vessel.Steer(steer);
      vessel.Update(t, true);
      t = t + dt;
    }

    return 0;
  }
}
