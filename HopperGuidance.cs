﻿#define LIMIT_ATTITUDE

using System;
using System.Collections.Generic;

using UnityEngine;
using UnityEngine.UI;

using KSPAssets;


namespace HopperGuidance
{
    enum AutoMode
    {
      Off,
      LandAtTarget,
      LandHere,
      Failed
    }

    public class HopperGuidance : PartModule
    {
        // Constants
        static Color trackcol = new Color(0,1,0,0.3f); // transparent green
        static Color targetcol = new Color(1,1,0,0.5f); // solid yellow
        static Color thrustcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        static Color idlecol = new Color(1,0.2f,0.2f,0.9f); // right red (idling as off attitude target)
        static Color aligncol = new Color(0,0.1f,1.0f,0.3f); // blue

        GameObject _tgt_obj = null; // new GameObject("Target");
        GameObject _track_obj = null; // new GameObject("Track");
        GameObject _align_obj = null; // new GameObject("Track");
        GameObject _steer_obj = null;
        LineRenderer _align_line = null;
        LineRenderer _steer_line = null;
        LinkedList<GameObject> thrusts = new LinkedList<GameObject>();
        PID3d _pid3d = new PID3d();
        bool checkingLanded = false; // only check once in flight to avoid failure to start when already on ground
        AutoMode autoMode = AutoMode.Off;
        float errMargin = 0.1f; // margin of error in solution to allow for headroom in amax (use double for maxThrustAngle)
        float predictTime = 0.2f;
        float finalDescentHeight = 0; // Descent straight down below this height
        double lowestY = 0; // Y position of bottom of craft relative to centre
        float _minThrust, _maxThrust;
        Solve solver; // Stores solution inputs, output and trajectory
        Trajectory _traj; // trajectory of solution in world space?
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _transform = null;
        double last_t = -1; // last time data was logged
        double log_interval = 0.05f; // Interval between logging
        System.IO.StreamWriter _tgtWriter = null; // Actual vessel
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        double extendTime = 1; // duration to extend trajectory to slowly descent to touch down and below at touchDownSpeed
        //double touchDownSpeed = 1.4f;
        double touchDownSpeed = 0;
        bool setShowTrack = true;
        float lastTgtLatitude, lastTgtLongitude, lastTgtAltitude;
        bool pickingPositionTarget = false;
        string _vesselLogFilename = "vessel.dat";
        string _tgtLogFilename = "target.dat";
        string _solutionLogFilename = "solution.dat";

        [UI_FloatRange(minValue = 5, maxValue = 500, stepIncrement = 5)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target size", guiFormat = "F0", isPersistant = false)]
        float tgtSize = 10;
        float setTgtSize;

        [UI_FloatRange(minValue = -90.0f, maxValue = 90.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Latitude", guiFormat = "F7", isPersistant = true, guiUnits="°")]
        float tgtLatitude = -0.0968071692165f; // H-Pad
        float setTgtLatitude;

        [UI_FloatRange(minValue = -180.0f, maxValue = 180.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Longitude", guiFormat = "F7", isPersistant = true, guiUnits = "°")]
        float tgtLongitude = -74.6172808614f; // H-Pad
        float setTgtLongitude;

        [UI_FloatRange(minValue = 0.1f, maxValue = 10000.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Altitude", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float tgtAltitude = 176; // H-Pad
        float setTgtAltitude;

        [UI_FloatRange(minValue = -1, maxValue = 90.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min descent angle", guiFormat = "F0", isPersistant = true, guiUnits = "°")]
        float minDescentAngle = 20.0f;
        float setMinDescentAngle;

        [UI_FloatRange(minValue = 1, maxValue = 1500, stepIncrement = 10f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max velocity", guiFormat = "F0", isPersistant = true, guiUnits = "m/s")]
        float maxV = 150f; // Max. vel to add to get towards target - not too large that vessel can't turn around
        float setMaxV;

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust angle", guiFormat = "F1", isPersistant = true, guiUnits="°")]
        float maxThrustAngle = 45f; // Max. thrust angle from vertical
        float setMaxThrustAngle;

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max landing thrust angle", guiFormat = "F1", isPersistant = true, guiUnits="°")]
        float maxLandingThrustAngle = 5f; // Max. final thrust angle from vertical
        float setMaxLandingThrustAngle;

        [UI_FloatRange(minValue = 0f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Idle attitude angle", guiFormat = "F0", isPersistant = true, guiUnits = "°")]
        float idleAngle = 90.0f;

        [UI_FloatRange(minValue = 0.01f, maxValue = 1f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Correction factor", guiFormat = "F2", isPersistant = true)]
        float corrFactor = 0.2f; // If 1 then at 1m error aim to close a 1m/s
        float setCorrFactor;
        // Set this down to 0.05 for large craft and up to 0.4 for very agile small craft

        [UI_FloatRange(minValue = 0, maxValue = 1, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Acceleration gain", guiFormat = "F1", isPersistant = true)]
        float accGain = 0.8f;
        float setAccGain;

        float yMult = 2; // Extra weight for Y to try and minimise height error over other errors

        // Can't get any improvement by raising these coeffs from zero
        float ki1 = 0;
        float ki2 = 0;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Keep engine ignited", isPersistant = false)]
        bool _keepIgnited = false;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Show track", isPersistant = false)]
        bool showTrack = true;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Logging", isPersistant = false)]
        bool _logging = false;

        // Quad should be described a,b,c,d in anti-clockwise order when looking at it
        public void AddQuad(Vector3[] vertices, int vi, int[] triangles, int ti,
                            Vector3d a, Vector3d b, Vector3d c, Vector3d d)
        {
          vertices[vi+0] = a;
          vertices[vi+1] = b;
          vertices[vi+2] = c;
          vertices[vi+3] = d;
          triangles[ti++] = vi;
          triangles[ti++] = vi+2;
          triangles[ti++] = vi+1;
          triangles[ti++] = vi;
          triangles[ti++] = vi+3;
          triangles[ti++] = vi+2;
        }

        public void DrawTarget(Vector3d pos, Transform transform, Color color, double size)
        {
          double[] r = new double[]{size*0.5,size*0.55,size*0.95,size};
          if (_tgt_obj != null)
          {
            Destroy(_tgt_obj);
          }
          pos = transform.TransformPoint(pos); // convert to World Pos

          Vector3d vx = new Vector3d(1,0,0);
          Vector3d vy = new Vector3d(0,1,0);
          Vector3d vz = new Vector3d(0,0,1);

          vx = transform.TransformVector(vx);
          vy = transform.TransformVector(vy);
          vz = transform.TransformVector(vz);

          _tgt_obj = new GameObject("Target");
          MeshFilter meshf = _tgt_obj.AddComponent<MeshFilter>();
          MeshRenderer meshr = _tgt_obj.AddComponent<MeshRenderer>();
          meshr.transform.parent = transform;
          meshr.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          meshr.material.color = color;

          Mesh mesh = new Mesh();
          Vector3[] vertices = new Vector3[36*4+4+4+4];
          int[] triangles = new int[(36*2*2-8+2+2+2)*3]; // take away gaps
          int i,j;
          int v=0,t=0;
          for(j=0;j<4;j++) // four concentric rings
          {
            for(i=0;i<36;i++)
            {
              float a = -(i*10)*Mathf.PI/180.0f;
              vertices[v++] = pos + vx*Mathf.Sin(a)*r[j] + vz*Mathf.Cos(a)*r[j];
            }
          }
          for(j=0;j<2;j++)
          {
            int start = j*72;
            for(i=0;i<36;i++)
            {
              if ((j==1) || (i%9!=0)) // make 4 gaps in inner ring
              {
                triangles[t++] = start+i;
                triangles[t++] = start+(i+1)%36;
                triangles[t++] = start+36+i%36;

                triangles[t++] = start+(i+1)%36;
                triangles[t++] = start+36+(i+1)%36;
                triangles[t++] = start+36+i%36;
              }
            }
          }
          // Add cross across centre
          Vector3 cx = vx*size*0.03;
          Vector3 cz = vz*size*0.03;
          float cs=8;
          AddQuad(vertices,v,triangles,t,
                  pos-cx*cs-cz,pos+cx*cs-cz,pos+cx*cs+cz,pos-cx*cs+cz);
          v+=4; t+=6;
          // One side
          AddQuad(vertices,v,triangles,t,
                  pos-cx+cz,pos+cx+cz,pos+cx+cz*cs,pos-cx+cz*cs);
          v+=4; t+=6;
          // Other size
          AddQuad(vertices,v,triangles,t,
                  pos-cx-cz*cs,pos+cx-cz*cs,pos+cx-cz,pos-cx-cz);
          v+=4; t+=6;

          mesh.vertices = vertices;
          mesh.triangles = triangles;
          mesh.RecalculateNormals();
          meshf.mesh = mesh;
        }

        public void DrawTrack(Trajectory traj, Transform transform, Color color, bool pretransform=true, float amult=1)
        {
            if (_track_obj != null)
            {
              Destroy(_track_obj); // delete old track
              // delete old thrusts
              LinkedListNode<GameObject> node = thrusts.First;
              while(node != null)
              {
                Destroy(node.Value);
                node = node.Next;
              }
              thrusts.Clear();
            }
            if (!showTrack)
              return;

            _track_obj = new GameObject("Track");
            LineRenderer line = _track_obj.AddComponent<LineRenderer>();
            line.transform.parent = transform;
            line.useWorldSpace = false;
            line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            line.material.color = color;
            line.startWidth = 0.4f;
            line.endWidth = 0.4f;
            line.positionCount = traj.Length();
            int j = 0;
            for (int i = 0; i < traj.Length(); i++)
            {
              if (pretransform)
                line.SetPosition(j++,transform.TransformPoint(traj.r[i]));
              else
                line.SetPosition(j++,traj.r[i]);
              // Draw accelerations
              GameObject obj = new GameObject("Accel");
              // TODO - Work out why only one LineRender per GameObject - it seems wrong!
              thrusts.AddLast(obj);
              LineRenderer line2 = obj.AddComponent<LineRenderer>();
              line2.transform.parent = transform;
              line2.useWorldSpace = false;
              line2.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
              line2.material.color = thrustcol;
              line2.startWidth = 0.4f;
              line2.endWidth = 0.4f;
              line2.positionCount = 2;
              if (pretransform)
              {
                line2.SetPosition(0,transform.TransformPoint(traj.r[i]));
                line2.SetPosition(1,transform.TransformPoint(traj.r[i] + traj.a[i]*amult));
              }
              else
              {
                line2.SetPosition(0,traj.r[i]);
                line2.SetPosition(1,traj.r[i] + traj.a[i]*amult);
              }
            }
        }

        public void DrawAlign(Vector3 r_from,Vector3 r_to, Transform transform, Color color)
        {
            if (!showTrack)
            {
              if (_align_obj != null)
              {
                Destroy(_align_obj);
                _align_obj = null;
                _align_line = null;
              }
              return;
            }
            if (_align_line == null)
            {
              _align_obj = new GameObject("Align");
              _align_line= _align_obj.AddComponent<LineRenderer>();
            }
            _align_line.transform.parent = transform;
            _align_line.useWorldSpace = true;
            _align_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _align_line.material.color = color;
            _align_line.startWidth = 0.3f;
            _align_line.endWidth = 0.3f;
            _align_line.positionCount = 2;
            _align_line.SetPosition(0,transform.TransformPoint(r_from));
            _align_line.SetPosition(1,transform.TransformPoint(r_to));
        }

        public void DrawSteer(Vector3 r_from,Vector3 r_to, Transform transform, Color color)
        {
            if (!showTrack)
            {
              if (_steer_obj != null)
              {
                Destroy(_steer_obj);
                _steer_obj = null;
                _steer_line = null;
              }
              return;
            }

            if (_steer_line == null)
            {
              _steer_obj = new GameObject("Steer");
              _steer_line= _steer_obj.AddComponent<LineRenderer>();
            }
            _steer_line.transform.parent = transform;
            _steer_line.useWorldSpace = true;
            _steer_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _steer_line.material.color = color;
            _steer_line.startWidth = 0.3f;
            _steer_line.endWidth = 0.3f;
            _steer_line.positionCount = 2;
            _steer_line.SetPosition(0,transform.TransformPoint(r_from));
            _steer_line.SetPosition(1,transform.TransformPoint(r_to));
        }

        public void OnDestroy()
        {
          DisableLand();
        }

        public void SetUpTransform(float lat, float lon, float alt)
        {
          // Set up transform so Y is up and (0,0,0) is target position
          CelestialBody body = vessel.mainBody;
          Vector3d origin = body.GetWorldSurfacePosition(lat, lon, alt);
          Vector3d vEast = body.GetWorldSurfacePosition(lat, lon-0.1, alt) - origin;
          Vector3d vUp = body.GetWorldSurfacePosition(lat, lon, alt+1) - origin;
          // Convert to body co-ordinates
          origin = body.transform.InverseTransformPoint(origin);
          vEast = body.transform.InverseTransformVector(vEast);
          vUp = body.transform.InverseTransformVector(vUp);

          GameObject go = new GameObject();
          // Need to rotation that converts (0,1,0) to vUp in the body transform
          Quaternion quat = Quaternion.FromToRotation(new Vector3(0,1,0),vUp);

          _transform = go.transform;
          _transform.SetPositionAndRotation(origin,quat);
          _transform.SetParent(body.transform,false);
        }

        public void LogStop()
        {
          if (_vesselWriter != null)
            _vesselWriter.Close();
          _vesselWriter = null;
          if (_tgtWriter != null)
            _tgtWriter.Close();
          _tgtWriter = null;
          last_t = -1;
        }

        void LogData(System.IO.StreamWriter f, double t, Vector3 r, Vector3 v, Vector3 a, float att_err)
        {
          f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1} {10:F1}",t,r.x,r.y,r.z,v.x,v.y,v.z,a.x,a.y,a.z,att_err));
        }

        // Find Y offset to lowest part from origin of the vessel
        double FindLowestPointOnVessel()
        {
          Vector3 CoM, up;

          CoM = vessel.localCoM;
          Vector3 bottom = Vector3.zero; // Offset from CoM
          up = FlightGlobals.getUpAxis(CoM); //Gets up axis
          Vector3 pos = vessel.GetWorldPos3D();
          Vector3 distant = pos - 1000*up; // distant below craft
          double miny = 0;
          foreach (Part p in vessel.parts)
          {
            if (p.collider != null) //Makes sure the part actually has a collider to touch ground
            {
              Vector3 pbottom = p.collider.ClosestPointOnBounds(distant); //Gets the bottom point
              double y = Vector3.Dot(up,pbottom-pos); // relative to centre of vessel
              if (y < miny)
              {
                bottom = pbottom;
                miny = y;
              }
            }
          }
          return miny;
        }

        public void ComputeMinMaxThrust(out float minThrust, out float maxThrust)
        {
          minThrust = 0;
          maxThrust = 0;
          foreach (Part part in vessel.GetActiveParts())
          {
              Debug.Log("Part: "+part);
              part.isEngine(out List<ModuleEngines> engines);
              foreach (ModuleEngines engine in engines)
              {
                  Debug.Log("Engine: "+engine);
                  //engine.Activate(); // must be active to get thrusts or else realIsp=0
                  for(float throttle=0; throttle<=1; throttle+=0.1f)
                    Debug.Log("isp="+engine.realIsp+" throttle="+throttle+" Thrust="+engine.GetEngineThrust(engine.realIsp,throttle));
                  // I think this will get the correct thrust given throttle in atmosphere (or wherever)
                  minThrust += engine.GetEngineThrust(engine.realIsp, 0);
                  maxThrust += engine.GetEngineThrust(engine.realIsp, 1);
              }
          }
        }

        public void ShutdownAllEngines()
        {
          foreach (Part part in vessel.GetActiveParts())
          {
            part.isEngine(out List<ModuleEngines> engines);
            foreach (ModuleEngines engine in engines)
              engine.Shutdown();
          }
        }

        // All positions and velocities supplied in local space relative to landing target at (0,0,0)
        // but the trajectory, world_traj, is also computed using the transform, given in transform
        public double ComputeTrajectory(ref Trajectory local_traj, Transform transform,
                                        ref Trajectory world_traj,
                                        Vector3d local_r, Vector3d local_v,
                                        List<Vector3d> local_int_r, List<Vector3d> local_int_v,
                                        Vector3d local_rf, Vector3d local_vf, float g,
                                        out double o_fuel, out int retval)
        {
          double dt = 0;
          double T = 0;
          Vector3d world_r = transform.TransformPoint(local_r);
          Vector3d world_v = transform.TransformVector(local_v);
          Vector3d world_g = transform.TransformVector(new Vector3d(0,-g,0));
          o_fuel = 0;
          retval = -1;

          local_int_r.Add(local_rf);
          local_int_v.Add(local_vf);
          solver.minDescentAngle = -1;
          for(int i=0 ; i < local_int_r.Count ; i++ )
          {
            // Compute trajectory to landing spot
            double fuel;
            Vector3d [] local_thrusts;
            // Currently uses intermediate positions, ir[], but ignores iv[
            double bestT = solver.GoldenSearchGFold(local_r, local_v, local_int_r[i], true, local_int_v[i], (i!=0), out local_thrusts, out fuel, out retval);
            Debug.Log(solver.DumpString());
            solver.minDescentAngle = minDescentAngle;
            if ((retval>=1) && (retval<=5))
            {
               T += solver.T;
               if (i==0)
                 dt = solver.dt; // set from first segment. TODO - might not be best
              // Simulate local trajectory given thrusts
              local_traj.Simulate(bestT, local_thrusts, local_r, local_v, new Vector3d(0,-world_g.magnitude,0), dt, extendTime);
              local_r = local_traj.r[local_traj.r.Length-1];
              local_v = local_traj.v[local_traj.v.Length-1];
              // Simulate world trajectory given transformed thrusts
              Vector3d [] world_thrusts = new Vector3d[local_thrusts.Length];
              for(int j=0; j<local_thrusts.Length; j++)
                world_thrusts[j] = transform.TransformVector(local_thrusts[j]);
              world_traj.Simulate(bestT, world_thrusts, world_r, world_v, world_g, dt, extendTime);
              world_r = world_traj.r[world_traj.r.Length-1];
              world_v = world_traj.v[world_traj.v.Length-1];
              o_fuel += fuel;
            }
            else
              return 0;
          }
          return T;
        }

        public void EnableLandAtTarget()
        {
          checkingLanded = false; // stop trajectory being cancelled while on ground
          lowestY = FindLowestPointOnVessel();
          //Vector3d[] thrusts;
          Vector3d r0 = vessel.GetWorldPos3D();
          Vector3d v0 = vessel.GetSrfVelocity();
          Vector3d g = FlightGlobals.getGeeForceAtPosition(r0);
          Vector3d rf = new Vector3d(0,-lowestY -0.5f*touchDownSpeed*extendTime,0);
          Vector3d vf = new Vector3d(0,-touchDownSpeed,0);
          ComputeMinMaxThrust(out _minThrust,out _maxThrust); // This might be including RCS (i.e. non main Throttle)
          _startTime = Time.time;
          double amin = _minThrust/vessel.totalMass;
          double amax = _maxThrust/vessel.totalMass;

          solver = new Solve();
          solver.Tmin = 1;
          solver.Tmax = 300; // 5 mins
          solver.tol = 0.1;
          solver.vmax = maxV;
          solver.amin = amin*(1+errMargin);
          solver.amax = amax*(1-errMargin);
          solver.Nmin = 2;
          solver.Nmax = 6;
          solver.minDurationPerThrust = 2;
          solver.g = g.magnitude;
          solver.minDescentAngle = minDescentAngle;
          solver.maxThrustAngle = maxThrustAngle*(1-2*errMargin);
          solver.maxLandingThrustAngle = maxLandingThrustAngle*(1-2*errMargin);

          int retval;
          // Predict into future since solution makes 0.1-0.3 secs to compute
          Vector3 att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          r0 = r0 + v0*predictTime + 0.5*g*predictTime*predictTime + 0.5f*(float)amin*att*predictTime*predictTime;
          v0 = v0 + g*predictTime + (float)amin*att*predictTime;

          // Shut-off throttle
          FlightCtrlState ctrl = new FlightCtrlState();
          vessel.GetControlState(ctrl);
          ctrl.mainThrottle = (_keepIgnited)?0.01f:0;

          // Compute trajectory to landing spot
          double fuel;
          List<Vector3d> ir = new List<Vector3d>();
          Vector3d tr0 = _transform.InverseTransformPoint(r0);
          Vector3d tv0 = _transform.InverseTransformVector(v0);
          List<Vector3d> iv = new List<Vector3d>();
          _traj = new Trajectory();
          Trajectory traj2 = new Trajectory();
          ComputeTrajectory(ref _traj, _transform, ref traj2, tr0, tv0, ir, iv, rf, vf, (float)g.magnitude, out fuel, out retval);
          Debug.Log(solver.DumpString());
          if ((retval>=1) && (retval<=5))
          {
            // Draw track computed in world space
            DrawTrack(traj2, _transform, trackcol, false);
            // Enable autopilot
            _pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult);
            // TODO - Testing out using in solution co-ordinates
            DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            vessel.OnFlyByWire += new FlightInputCallback(Fly);
            // Write solution
            if (_logging) {_traj.Write(_solutionLogFilename);}
            autoMode = AutoMode.LandAtTarget;
            Events["ToggleLandAtTarget"].guiName = "Cancel land at target";
          }
          else
          {
            DisableLand();
            Events["ToggleLandAtTarget"].guiName = "Failed! - Cancel land at target";
            string msg = "HopperGuidance: Failure to find solution because ";
            Vector3d r = tr0 - rf;
            double cos_descentAng = Math.Sqrt(r.x*r.x + r.z*r.z) / r.magnitude;
            // Do some checks
            if (v0.magnitude > maxV)
              msg = msg + " velocity over "+maxV+" m/s";
            else if (cos_descentAng < Math.Cos(Math.PI*minDescentAngle))
              msg = msg + "below min. descent angle "+minDescentAngle+"°";
            else if (amax < g.magnitude)
              msg = msg + "engine has insufficient thrust, no engines active or no fuel";
            else
              msg = msg + "impossible to reach target within constraints";
            Debug.Log(msg);
            ScreenMessages.PostScreenMessage(msg, 3.0f, ScreenMessageStyle.UPPER_CENTER);
            autoMode = AutoMode.Failed;
          }
        }

        public void DisableLand()
        {
          autoMode = AutoMode.Off;
          if (_track_obj != null) {Destroy(_track_obj); _track_obj=null;}
          if (_tgt_obj   != null) {Destroy(_tgt_obj);     _tgt_obj=null;}
          if (_align_obj != null) {Destroy(_align_obj); _align_obj=null;}
          if (_steer_obj != null) {Destroy(_steer_obj); _steer_obj=null;}
          LinkedListNode<GameObject> node = thrusts.First;
          while(node != null)
          {
            Destroy(node.Value);
            node = node.Next;
          }
          thrusts.Clear();
          LogStop();
          vessel.OnFlyByWire -= new FlightInputCallback(Fly);
          Events["ToggleLandHere"].guiName = "Land here";
          Events["ToggleLandAtTarget"].guiName = "Land at target";
          vessel.Autopilot.Disable();
        }

        ~HopperGuidance()
        {
          DisableLand();
        }

        public override void OnInactive()
        {
          Debug.Log("OnInactive()");
          base.OnInactive();
          DisableLand();
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Pick target", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void PickTarget()
        {
          pickingPositionTarget = true;
          string message = "Click to select a target";
          ScreenMessages.PostScreenMessage(message, 3.0f, ScreenMessageStyle.UPPER_CENTER);
          lastTgtLongitude = tgtLongitude;
          lastTgtLatitude = tgtLatitude;
          lastTgtAltitude = tgtAltitude;
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Land at target", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void ToggleLandAtTarget()
        {
          if (autoMode == AutoMode.Failed)
            DisableLand();
          else
          {
            if (autoMode != AutoMode.LandAtTarget)
              EnableLandAtTarget();
            else
              DisableLand();
          }
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Land here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void ToggleLandHere()
        {
          if (autoMode != AutoMode.LandHere)
          {
            DisableLand();
            lowestY = FindLowestPointOnVessel();
            ComputeMinMaxThrust(out _minThrust,out _maxThrust); // This might be including RCS (i.e. non main Throttle)
            float amax = (float)(_maxThrust/vessel.totalMass); // F = m*a. We want, unit of throttle for each 1m/s/s
            _pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult);
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            vessel.OnFlyByWire += new FlightInputCallback(Fly);
            Events["ToggleLandHere"].guiName = "Cancel land here";
            SetUpTransform((float)vessel.latitude,(float)vessel.longitude,(float)vessel.mainBody.TerrainAltitude(vessel.latitude,vessel.longitude));
            _startTime = Time.time;
            autoMode = AutoMode.LandHere;
          }
          else
            DisableLand();
        }

        Vector3d GetLandHereDesiredVel(float height, Vector3d tv, float g, float amax)
        {
          // Find time to hit ground with just gravity
          // 0 = a*x*x + b*x + c;
          // height = 0.5*amax*t*t + v*t; 
          float a = 0.5f*g;
          float b = (float)tv.magnitude;
          float c = -height;
          float x1 = (-b + Mathf.Sqrt(b*b-4*a*c))/(2*a);
          float x2 = (-b - Mathf.Sqrt(b*b-4*a*c))/(2*a);
          float t = x1;
          if (x2 > 0)
            t = x2;
          float maxv = t*(amax-g)*0.75f;
          return new Vector3d(0,-maxv,0);
        }

        // F should be initial acceleration before any corrections
        // For forced landing this is the opposite vector to gravity (to hover)
        // For an autopilot path this is the acceleration vector from the computed solution
        public void AutopilotStepToTarget(FlightCtrlState state, Vector3d tr, Vector3d tv, Vector3d dr, Vector3d dv, Vector3d da, float g)
        {
          float throttle = 0;

          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          Vector3d tatt = _transform.InverseTransformVector(att);
          ComputeMinMaxThrust(out _minThrust,out _maxThrust);
          float amax = (float)(_maxThrust/vessel.totalMass);
          float amin = (float)(_minThrust/vessel.totalMass);

          Vector3d da2 = GetThrustVector(tr,tv,dr,dv,amin,amax,maxThrustAngle,da, out Vector3d unlimda);
          if (_keepIgnited)
            throttle = Mathf.Clamp((float)(da2.magnitude-amin)/(amax-amin),0.01f,1);
          else
            throttle = Mathf.Clamp((float)(da2.magnitude-amin)/(amax-amin),0,1);

          // Shutoff throttle if pointing in wrong direction
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(tatt),Vector3d.Normalize(da2));
          float att_err = Mathf.Acos(ddot)*180/Mathf.PI;

          // Open log files
          if ((_vesselWriter == null) && (_logging))
          {
            _vesselWriter = new System.IO.StreamWriter(_vesselLogFilename);
            _vesselWriter.WriteLine("time x y z vx vy vz ax ay az att_err");
            _tgtWriter = new System.IO.StreamWriter(_tgtLogFilename);
            _tgtWriter.WriteLine("time x y z vx vy vz ax ay az att_err");
          }

          double t = Time.time - _startTime;
          if ((_logging)&&(t >= last_t+log_interval))
          {
            LogData(_tgtWriter, t, dr, dv, da, 0);
            LogData(_vesselWriter, t, tr, tv, da2, att_err);
            last_t = t;
          }

          if ((att_err >= idleAngle) && (da2.magnitude>0.01))
          {
            throttle = 0.01f; // some throttle to steer? (if no RCS and main thruster gimbals)
            // Draw steer vector
            DrawSteer(tr, tr+3*vessel.vesselSize.x*Vector3d.Normalize(unlimda), _transform, idlecol);
          }
          else
          {
            // Draw steer vector
            DrawSteer(tr, tr+3*vessel.vesselSize.x*Vector3d.Normalize(unlimda), _transform, thrustcol);
          }
          DrawAlign(tr,dr,_transform,aligncol);
          // If no clear steer direction point upwards
          if (da2.magnitude < 0.1f)
            da2 = new Vector3d(0,1,0);
          da2 = _transform.TransformVector(da2); // transform back to world co-ordinates
          vessel.Autopilot.SAS.lockedMode = false;
          vessel.Autopilot.SAS.SetTargetOrientation(da2,false);
          state.mainThrottle = throttle;
        }

        // returns the thrust cone limited vector
        // accepted original thrust vector as F, and returns unlimited F after adjustment
        public Vector3d GetThrustVector(Vector3d tr, Vector3d tv, Vector3d dr, Vector3d dv, float amin, float amax, float maxThrustAngle, Vector3d F, out Vector3d unlimF)
        {
          Vector3 F2 = _pid3d.Update(tr,tv,dr,dv,Time.deltaTime);
          unlimF = F + F2;
          // Reduce sideways components
          return ConeUtils.ClosestThrustInsideCone((float)maxThrustAngle,(float)amin,(float)amax,unlimF);
        }

        public void DesiredPosVelAccForLandHere(Vector3d tr, Vector3d tv, out Vector3d dr, out Vector3d dv, out Vector3d da, float g)
        {
          // Find height via ray
          RaycastHit hit = new RaycastHit();
          Vector3d up = FlightGlobals.getUpAxis(vessel.localCoM); //Gets up axis
          Vector3d origin = vessel.GetWorldPos3D() + up*(float)lowestY;
          bool isHit = Physics.Raycast(origin, -up, out hit, Mathf.Infinity, 1 << 15);
          if (!isHit)
          {
            dr = tr;
            dv = tv;
            da = Vector3.zero;
            return;
          }

          ComputeMinMaxThrust(out _minThrust,out _maxThrust);
          float amax = (float)(_maxThrust/vessel.totalMass);

          // Find height. Vertical distance from bottom of craft to intersection
          float rheight = Vector3.Dot(up,origin - hit.point);
          Vector3 tg = new Vector3(0,-g,0);
          dv = GetLandHereDesiredVel(rheight, tv, g, amax);
          da = -tg; // compensate for gravity
          Debug.Log("HopperGuidance: h="+rheight+" v="+(Vector3d)tv+" dv="+(Vector3d)dv);

          //if (dv.magnitude > tv.magnitude) // desired/max velocity is downwards
          //  dv = tv;
          dr = tr; // Don't change position
        } 
        
        public void DesiredPosVelAccForLandAtTarget(Vector3d tr, Vector3d tv, out Vector3d dr, out Vector3d dv, out Vector3d da, float g)
        {
          double desired_t; // closest time in trajectory (desired)
          // TODO - Transform with _transform
          _traj.FindClosest(tr, tv, out dr, out dv, out da, out desired_t, 0.5f, 0.5f);
        }

        public void Fly(FlightCtrlState state)
        {
          if ((vessel == null) || (vessel.checkLanded() && checkingLanded) )
          {
            ScreenMessages.PostScreenMessage("Landed!", 3.0f, ScreenMessageStyle.UPPER_CENTER);
            DisableLand();
            // Shut-off throttle
            FlightCtrlState ctrl = new FlightCtrlState();
            vessel.GetControlState(ctrl);
            ctrl.mainThrottle = 0;
            // Necessary on Realism Overhaul to shutdown engine as at throttle=0 the engine may still have
            // a lot of thrust
            if (_keepIgnited)
              ShutdownAllEngines();
            autoMode = AutoMode.Off;
            return;
          }
          // Only start checking if landed when taken off
          if (!vessel.checkLanded())
            checkingLanded = true;

          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3 tr = _transform.InverseTransformPoint(r);
          Vector3 tv = _transform.InverseTransformVector(v);
          float g = (float)FlightGlobals.getGeeForceAtPosition(r).magnitude;

          Vector3d dr,dv,da;
          // Switch to land here when close to ground
          if ((autoMode == AutoMode.LandHere) && (tr.y < finalDescentHeight - lowestY))
            DesiredPosVelAccForLandHere(tr,tv,out dr,out dv,out da,g);
          else
            DesiredPosVelAccForLandAtTarget(tr,tv,out dr,out dv,out da,g);

          // Uses transformed positions and vectors
          // sets throttle and desired attitude based on targets
          // F is the inital force (acceleration) vector
          AutopilotStepToTarget(state, tr, tv, dr, dv, da, g);
        }
        
        public override void OnUpdate()
        {
          base.OnUpdate();
          // Check for changes to slider. This can trigger one of many actions
          // - If autopilot enabled, recompute trajectory
          // - Reset PID controller
          // - Redraw target
          // - Draw/delete track
          bool recomputeTrajectory = false;
          bool redrawTarget = false;
          bool resetPID = false;
          bool resetHop = false;
          if ((tgtLatitude != setTgtLatitude) || (tgtLongitude != setTgtLongitude)) 
          {
            redrawTarget = true;
            recomputeTrajectory = true;
            resetHop = true;
            setTgtLatitude = tgtLatitude;
            setTgtLongitude = tgtLongitude;
          }
          if (tgtSize != setTgtSize)
          {
            setTgtSize = tgtSize;
            redrawTarget = true;
          }
          if ((tgtAltitude != setTgtAltitude) || (tgtSize != setTgtSize))
          {
            setTgtAltitude = tgtAltitude;
            redrawTarget = true;
            recomputeTrajectory = true;
            resetHop = true;
          }
          if ((corrFactor != setCorrFactor) || (accGain != setAccGain))
          {
            setCorrFactor = corrFactor;
            setAccGain = accGain;
            resetPID = true;
          }
          if ((minDescentAngle != setMinDescentAngle)||(maxThrustAngle != setMaxThrustAngle)||(maxLandingThrustAngle != setMaxLandingThrustAngle))
          {
            setMinDescentAngle = minDescentAngle;
            setMaxThrustAngle = maxThrustAngle;
            setMaxLandingThrustAngle = maxLandingThrustAngle;
            recomputeTrajectory = true;
          }
          if (maxV != setMaxV)
          {
            setMaxV = maxV;
            resetPID = true;
            recomputeTrajectory = true;
          }  
          if (showTrack != setShowTrack)
          {
            DrawTrack(_traj, _transform, trackcol);
            setShowTrack = showTrack;
          }
          if (pickingPositionTarget)
          {
            if (Input.GetKeyDown(KeyCode.Escape))
            {
              // Previous position
              tgtLatitude = lastTgtLatitude;
              tgtLongitude = lastTgtLongitude;
              tgtAltitude = lastTgtAltitude;
              redrawTarget = true;
              pickingPositionTarget = false;
            }
            RaycastHit hit;
            if (GuiUtils.GetMouseHit(vessel.mainBody,out hit,part))
            {
              // Picked
              double lat, lon, alt;
              vessel.mainBody.GetLatLonAlt(hit.point, out lat, out lon, out alt);
              tgtLatitude = (float)lat;
              tgtLongitude = (float)lon;
              tgtAltitude = (float)(alt+0.5);
              redrawTarget = true;
              // If clicked stop picking
              if (Input.GetMouseButtonDown(0))
              {
                recomputeTrajectory = true;
                pickingPositionTarget = false;
              }
            }
          }
          // Activate the required updates
          if (redrawTarget)
          {
            setTgtLatitude = tgtLatitude;
            setTgtLongitude = tgtLongitude;
            setTgtAltitude = tgtAltitude;
            setTgtSize = tgtSize;
            SetUpTransform(tgtLatitude, tgtLongitude, tgtAltitude);
            DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
          }
          if ((recomputeTrajectory)&&((autoMode == AutoMode.LandAtTarget)||(autoMode == AutoMode.Failed)))
            EnableLandAtTarget();
          if (resetPID)
          {
            ComputeMinMaxThrust(out _minThrust,out _maxThrust); // This might be including RCS (i.e. non main Throttle)
            double amax = _maxThrust/vessel.totalMass;
            _pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult);
          }
          if ((!_logging) && (_vesselWriter != null))
          {
            LogStop();
          }
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Set target here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void SetTargetHere()
        {
            // Find vessel co-ordinates
            tgtLatitude = (float)vessel.latitude;
            tgtLongitude = (float)vessel.longitude;
            // Note: compensate for height of vessel by getting bottom Y of vessel
            lowestY = FindLowestPointOnVessel();
            Vector3 up = FlightGlobals.getUpAxis(vessel.GetWorldPos3D());
            tgtAltitude = (float)(FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D() + up*(float)lowestY));
            SetUpTransform(tgtLatitude, tgtLongitude, tgtAltitude);
            DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
            setTgtLatitude = tgtLatitude;
            setTgtLongitude = tgtLongitude;
            setTgtAltitude = tgtAltitude;
        }
    }
}
