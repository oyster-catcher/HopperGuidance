using System.Collections.Generic;

using UnityEngine;
using UnityEngine.UI;

using KSPAssets;

//#define FLYPID

namespace HopperGuidance
{
    public class HopperGuidance : PartModule
    {
        GameObject _tgt_obj = null; // new GameObject("Target");
        LineRenderer _tgt_line = null;
        GameObject _track_obj = null; // new GameObject("Track");
        LineRenderer _track_line = null;
        GameObject _align_obj = null; // new GameObject("Track");
        LineRenderer _align_line = null;
        PID3d _pid3d = new PID3d();
        bool _enabled = false;
        float _maxThrust;
        Color trackcol = new Color(0,1,0,0.5f); // transparent green
        Color targetcol = new Color(1,0.2f,0.2f,0.5f); // transparent red
        Color aligncol = new Color(0,0.1f,1.0f,1.0f); // blue
        Trajectory _traj; // trajectory of solutio
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _logTransform = null;
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        Vector3d _targetPos;
        float _accelFactor = 1.21f; // max. accel = maxThrust/(totalMass*_accelFactor)
        bool _logging = true;
        

        //[UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Latitude", guiFormat = "F7", isPersistant = true)]
        double tgtLatitude = -0.0972078;

        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Longitude", guiFormat = "F7", isPersistant = true)]
        double tgtLongitude = -74.5576822;

        [UI_FloatRange(minValue = 0.1f, maxValue = 1000.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Altitude", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float tgtAltitude = 74.7f;

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min. descent angle", guiFormat = "F0", isPersistant = true, guiUnits = "m")]
        float minDescentAngle = 30.0f;

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min. final descent angle", guiFormat = "F0", isPersistant = true, guiUnits = "m")]
        float minFinalDescentAngle = 30.0f;

        [UI_FloatRange(minValue = 0f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Idle attitude angle", guiFormat = "F0", isPersistant = true)]
        float idleAngle = 45.0f;

        [UI_FloatRange(minValue = 0, maxValue = 100f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust %", isPersistant = true, guiUnits = "%")]
        float maxPercentThrust = 50f;

        [UI_FloatRange(minValue = 0, maxValue = 100f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min thrust (downwards) %", isPersistant = true, guiUnits = "%")]
        float minPercentThrust = 1f;

        [UI_FloatRange(minValue = 0.01f, maxValue = 1f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Position gain", guiFormat = "F2", isPersistant = true)]
        float kP1 = 0.3f; // If 1 then at 1m error aim to close a 1m/s

        [UI_FloatRange(minValue = 1f, maxValue = 100f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Max velocity", guiFormat = "F1", isPersistant = true)]
        float maxErrV = 50f; // Max. vel to add to get towards target

        [UI_FloatRange(minValue = 0.0f, maxValue = 2f, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Velocity gain", guiFormat = "F1", isPersistant = true)]
        float kP2 = 1.0f; // If 1 then at 1m/s error in velocity acceleration at an extra 1m/s/s

        [UI_FloatRange(minValue = 0.0f, maxValue = 100.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Max thrust", guiFormat = "F1", isPersistant = true)]
        float maxErrPercentThrust = 20.0f; // Max acceleration to add to get towards target


        public void DrawTarget(Vector3d pos, Transform transform, Color color, double size=10)
        {
          if (_tgt_line == null)
          {
            _tgt_obj = new GameObject("Target");
            _tgt_line= _tgt_obj.AddComponent<LineRenderer>();
          }
          _tgt_line.transform.parent = transform;
          _tgt_line.useWorldSpace = false;
          _tgt_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          _tgt_line.material.color = color;
          _tgt_line.startWidth = 0.4f;
          _tgt_line.endWidth = 0.4f;
          _tgt_line.positionCount = 8;
          _tgt_line.SetPosition(0,pos-transform.up);
          _tgt_line.SetPosition(1,pos+transform.up);
          _tgt_line.SetPosition(2,pos);
          _tgt_line.SetPosition(3,pos-transform.forward);
          _tgt_line.SetPosition(4,pos+transform.forward);
          _tgt_line.SetPosition(5,pos);
          _tgt_line.SetPosition(6,pos-transform.right);
          _tgt_line.SetPosition(7,pos+transform.right);
        }

        public void DrawTrack(Trajectory traj, Transform transform, Color color, float amult=1)
        {
            Debug.Log("DrawTrack N="+traj.Length());
            if (_track_line == null)
            {
              Debug.Log("Making gameObject");
              _track_obj = new GameObject("Track");
              _track_line= _track_obj.AddComponent<LineRenderer>();
            }
            _track_line.transform.parent = transform;
            _track_line.useWorldSpace = false;
            _track_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _track_line.material.color = color;
            _track_line.startWidth = 0.4f;
            _track_line.endWidth = 0.4f;
            _track_line.positionCount = traj.Length();
            int j = 0;
            for (int i = 0; i < traj.Length(); i++)
            {
              _track_line.SetPosition(j++,_traj.r[i]);
              //_track_line.SetPosition(j++,_traj.r[i] + _traj.a[i]*amult);
              //_track_line.SetPosition(j++,_traj.r[i]);
            }
        }

        public void DrawAlign(Vector3 r_from,Vector3 r_to, Transform transform, Color color)
        {
            if (_align_line == null)
            {
              Debug.Log("Making gameObject");
              _align_obj = new GameObject("Track");
              _align_line= _align_obj.AddComponent<LineRenderer>();
            }
            _align_line.transform.parent = transform;
            _align_line.useWorldSpace = false;
            _align_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _align_line.material.color = color;
            _align_line.startWidth = 0.4f;
            _align_line.endWidth = 0.4f;
            _align_line.positionCount = 2;
            _align_line.SetPosition(0,r_from);
            _align_line.SetPosition(1,r_to);
        }


        public override void OnStart(StartState state)
        {
            base.OnStart(state);
            Debug.Log("OnStart");
        }

        public void LogSetUpTransform()
        {
          // Set up transform so Y is up and (0,0,0) is target position
          GameObject go = new GameObject();
          CelestialBody body = vessel.mainBody;
          Vector3d origin = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
          Vector3d vEast = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude-0.1, tgtAltitude) - origin;
          Vector3d vUp = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude+1) - origin;
          Quaternion quat = Quaternion.LookRotation(vEast,vUp);
          _logTransform = go.transform;
          _logTransform.SetPositionAndRotation(origin,quat);
        }

        public void LogStart(string filename)
        {
          if (_vesselWriter == null)
          {
            _vesselWriter = new System.IO.StreamWriter(filename);
          }
          _vesselWriter.WriteLine("time x y z vx vy vz ax ay az");
        }

        public void LogStop()
        {
          if (_vesselWriter != null)
          {
            _vesselWriter.Close();
          }
          _vesselWriter = null;
        }

        void LogData(System.IO.StreamWriter f, double t, Vector3 r, Vector3 v, Vector3 a)
        {
          if (f != null)
          {
            r = _logTransform.InverseTransformPoint(r);
            v = _logTransform.InverseTransformVector(v);
            a = _logTransform.InverseTransformVector(a);
            f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1}",t,r.x,r.y,r.z,v.x,v.y,v.z,a.x,a.y,a.z));
          }
        }

        public float ComputeMaxThrust()
        {
          float maxThrust = 0;
          foreach (Part part in vessel.GetActiveParts())
          {
              part.isEngine(out List<ModuleEngines> engines);
              foreach (ModuleEngines engine in engines)
              {
                  maxThrust += engine.GetMaxThrust();
              }
          }
          return maxThrust;
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Enable autopilot", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void OnToggle()
        {
            if (!_enabled)
            {

                CelestialBody body = vessel.mainBody;
                Vector3d targetPos = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);

                // Compute trajectory to landing spot
                double Tmin = 1;
                double Tmax = 60;
                double tol = 0.5;
                int N = 5;
                Vector3d[] thrusts;
                Vector3d r0 = vessel.GetWorldPos3D();
                Vector3d v0 = vessel.GetSrfVelocity();
                Vector3d att0 = vessel.transform.up; // pointing direction
                Vector3d g = FlightGlobals.getGeeForceAtPosition(r0);
                _maxThrust = ComputeMaxThrust();
                double amax = (_maxThrust/vessel.totalMass)/_accelFactor;
                _startTime = Time.time;
                int retval;
                // Use _logTransform so that Y is vertical direction, and the gravity which acts downwards in the Y direction
                LogSetUpTransform(); // initialises _logTransform
                double bestT = Solve.golden_search_gfold(r0, v0, att0, targetPos, Vector3d.zero, Tmin, Tmax, N, g.magnitude, _logTransform, 0.01*maxPercentThrust*amax, tol, minDescentAngle, minFinalDescentAngle, 0.01, out retval,out thrusts);
                if (retval > 0)
                {
                  for(int i=0; i<N; i++)
                  {
                    Debug.Log("a["+i+"] = "+thrusts[i]);
                  }
                  Debug.Log("Found Solution: bestT=" + bestT + " retval="+retval);
                  double dt = 0.05;

                  _traj = new Trajectory();
                  _traj.Simulate(bestT, thrusts, r0, v0, g, dt, 5.0);
                  _traj.CorrectFinal(targetPos, Vector3.zero);
                  Debug.Log("Traj length="+_traj.Length());

                  // Enable autopilot
                  _pid3d.Init(kP1,0,0,kP2,0,0,maxErrV,(float)(0.01f*maxErrPercentThrust*amax));
                  _targetPos = vessel.mainBody.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
                  DrawTarget(_targetPos, body.transform, targetcol);
                  DrawTrack(_traj, body.transform, trackcol);
                  vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
                  Debug.Log("Enabled");
                  vessel.OnFlyByWire += new FlightInputCallback(Fly);
                  if (_logging) {LogStart("vessel.dat");}
                  // Write solution
                  if (_logging) {_traj.WriteLog("solution.dat", _logTransform);}
                  _enabled = true;
                  Events["OnToggle"].guiName = "Disable autopilot";
                }
                else
                {
                  Events["OnToggle"].guiName = "Enable fail: Try again";
                  Debug.Log("Failure to find solution: error code "+retval);
                }
            }
            else
            {
                Debug.Log("Disabled");
                if (_logging) {LogStop();}
                vessel.Autopilot.Disable();
                vessel.OnFlyByWire -= new FlightInputCallback(Fly);
                _enabled = false;
                Events["OnToggle"].guiName = "Enable autopilot";
            } 
        }

        // Use 3-D pairs of PIDs to calculate adjustment to force vector and optionally transfrom
        // space where Y is vertical
        public Vector3d GetFAdjust(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv, float dt, Transform transform=null)
        {
            if (transform)
            {
              Vector3d r2 = transform.InverseTransformPoint(r);
              Vector3d v2 = transform.InverseTransformVector(v);
              Vector3d dr2 = transform.InverseTransformPoint(dr);
              Vector3d dv2 = transform.InverseTransformVector(dv);
              Vector3d F2 = _pid3d.Update(r2,
                                          v2,
                                          dr2,
                                          dv2,
                                          Time.deltaTime);
              Debug.Log("Err. r="+(Vector3)(dr2-r2)+" Err. v="+(Vector3)(dv2-v2)+ " F2="+(Vector3)F2);
              return transform.TransformVector(F2);
            }
            else
            {
              return _pid3d.Update(r,v,dr,dv,Time.deltaTime);
            }
        }

        public void Fly(FlightCtrlState state)
        {
          // Get nearest point on trajectory in position and velocity
          // use both so we can handle when the trajectory crosses itself,
          // and it works a bit better anyway if velocity and position are not in step
          if (vessel == null)
          {
            if (_logging) {LogStop();}
            vessel.OnFlyByWire -= new FlightInputCallback(Fly);
            _enabled = false;
            Events["OnToggle"].guiName = "Enable autopilot";
          }

          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3d dr,dv,da;
          double desired_t; // closest time in trajectory (desired)
          bool found = _traj.FindClosest(r, v, out dr, out dv, out da, out desired_t);
          DrawAlign(r,dr,vessel.mainBody.transform,aligncol);
          float throttle=0;
          Vector3 F = vessel.mainBody.transform.up; // default to up
          
          if (found)
          {
            if (da.magnitude > 0.01)
            {
              F = da;
            }

            // Adjust thrust direction when off trajectory
            Vector3d F2 = GetFAdjust(r,v,dr,dv,Time.deltaTime,_logTransform);
            //F = F + (dr-r)*_rGain + (dv-v)*_vGain;
            F = F + F2;
            Debug.Log("F="+(Vector3d)F+"  F(trans)="+(Vector3d)_logTransform.InverseTransformPoint(F));
            float amax = (float)(_maxThrust/vessel.totalMass);  // F = m.a. We want, unit of throttle for each 1m/s/s
            Debug.Log("amax="+amax+ "maxThrust="+_maxThrust+" mass="+vessel.totalMass+ "mag(F)="+F.magnitude);
            //Vector3 g = FlightGlobals.getGeeForceAtPosition(r); // since F2 excludes g
            throttle = _accelFactor*F.magnitude/(amax+0.001f); // protect against divide by zero // 2* is HACK!!!
            //throttle = (float)(F.magnitude/(amax+0.001f));
            Debug.Log("t="+(Time.time-_startTime)+ "close_t="+desired_t+" vel="+v.magnitude+" desired_vel="+dv.magnitude+" spec.accel="+vessel.specificAcceleration);
          }
          else
          {
            Debug.Log("Not found");
          }

          // Shutoff throttle if pointing in wrong direction
          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(att),Vector3d.Normalize(F));
          Debug.Log("ddot="+ddot+" at 45="+Mathf.Cos(45*Mathf.PI/180.0f));
          if (ddot < Mathf.Cos((float)idleAngle*(Mathf.PI/180.0f)))
          {
            throttle = 0.01f;
          }

          // Set throttle and direction
          if (F.magnitude > 0.01)
          {
            // Only steer for significant direction
            vessel.Autopilot.SAS.SetTargetOrientation(F,false);
          }
          vessel.Autopilot.SAS.lockedMode = false;
          state.mainThrottle = throttle;
          Debug.Log("throttle="+throttle);

          // Logging
          if (_logging) {LogData(_vesselWriter,Time.time-_startTime,vessel.GetWorldPos3D(),vessel.GetSrfVelocity(),F);}
        }

//        public override void OnUpdate()
//        {
//          base.OnUpdate();
//          _targetPos = vessel.mainBody.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
//          DrawTarget(_targetPos,vessel.mainBody.transform,targetcol);
//        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Set Target Here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void SetTargetHere()
        {
            // Find vessel co-ordinates
            tgtLatitude = vessel.latitude;
            tgtLongitude = vessel.longitude;
            tgtAltitude = (float)FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D());
            CelestialBody body = vessel.mainBody;
            Vector3d targetPos = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
            Debug.Log("Error in target pos="+(targetPos-vessel.GetWorldPos3D()).magnitude);
            DrawTarget(targetPos,body.transform,targetcol);
        }
    }
}
