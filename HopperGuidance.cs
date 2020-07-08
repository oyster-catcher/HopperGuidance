using System.Collections.Generic;

using UnityEngine;
using UnityEngine.UI;

using KSPAssets;

namespace HopperGuidance
{
    public class HopperGuidance : PartModule
    {
        GameObject _tgt_obj = null; // new GameObject("Target");
        GameObject _track_obj = null; // new GameObject("Track");
        GameObject _align_obj = null; // new GameObject("Track");
        LineRenderer _align_line = null;
        LinkedList<GameObject> thrusts = new LinkedList<GameObject>();
        PID3d _pid3d = new PID3d();
        bool _enabled = false;
        double predictTime = 0.2;
        float _maxThrust;
        Color trackcol = new Color(0,1,0,0.3f); // transparent green
        Color targetcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        Color thrustcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        Color aligncol = new Color(0,0.1f,1.0f,0.3f); // blue
        Trajectory _traj; // trajectory of solution
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _logTransform = null;
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        Vector3d _targetPos;
        float _accelFactor = 1.0f; // max. accel = maxThrust/(totalMass*_accelFactor)
        bool _logging = true;
        double disengageDistance = 10000;
        double extendTime = 5; // extend trajectory to slowly descent to touch down

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Latitude", guiFormat = "F7", isPersistant = false)]
        float tgtLatitude = -0.0972078f;

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Longitude", guiFormat = "F7", isPersistant = false)]
        float tgtLongitude = -74.5576822f;

        [UI_FloatRange(minValue = 0.1f, maxValue = 1000.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Altitude", guiFormat = "F1", isPersistant = false, guiUnits = "m")]
        float tgtAltitude = 74.7f;

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min descent angle", guiFormat = "F0", isPersistant = false, guiUnits = "m")]
        float minDescentAngle = 20.0f;

        [UI_FloatRange(minValue = 1f, maxValue = 100f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max velocity", guiFormat = "F1", isPersistant = false)]
        float maxErrV = 50f; // Max. vel to add to get towards target

        [UI_FloatRange(minValue = 5f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust angle", guiFormat = "F1", isPersistant = false)]
        float maxThrustAngle = 45f; // Max. thrust angle from vertical

        [UI_FloatRange(minValue = 5f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max final thrust angle", guiFormat = "F1", isPersistant = false)]
        float maxFinalThrustAngle = 20f; // Max. final thrust angle from vertical

        [UI_FloatRange(minValue = 0, maxValue = 200f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust %", isPersistant = false, guiUnits = "%")]
        float maxPercentThrust = 100f;

        [UI_FloatRange(minValue = 0f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Idle attitude angle", guiFormat = "F0", isPersistant = false)]
        float idleAngle = 90.0f;

        [UI_FloatRange(minValue = 0.01f, maxValue = 1f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Position gain", guiFormat = "F2", isPersistant = false)]
        float kP1 = 0.3f; // If 1 then at 1m error aim to close a 1m/s


        [UI_FloatRange(minValue = 0.0f, maxValue = 2f, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Velocity gain", guiFormat = "F1", isPersistant = false)]
        float kP2 = 1.0f; // If 1 then at 1m/s error in velocity acceleration at an extra 1m/s/s

        [UI_FloatRange(minValue = 0.0f, maxValue = 100.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Max thrust", guiFormat = "F1", isPersistant = false)]
        float maxErrPercentThrust = 20.0f; // Max acceleration to add to get towards target


        public void DrawTarget(Vector3d pos, Transform transform, Color color, double size=10)
        {
          if (_tgt_obj != null)
          {
            Destroy(_tgt_obj);
          }
          _tgt_obj = new GameObject("Target");
          LineRenderer line= _tgt_obj.AddComponent<LineRenderer>();
          line.transform.parent = transform;
          line.useWorldSpace = false;
          line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          line.material.color = color;
          line.startWidth = 0.4f;
          line.endWidth = 0.4f;
          line.positionCount = 8;
          line.SetPosition(0,pos-transform.up);
          line.SetPosition(1,pos+transform.up);
          line.SetPosition(2,pos);
          line.SetPosition(3,pos-transform.forward);
          line.SetPosition(4,pos+transform.forward);
          line.SetPosition(5,pos);
          line.SetPosition(6,pos-transform.right);
          line.SetPosition(7,pos+transform.right);
        }

        public void DrawTrack(Trajectory traj, Transform transform, Color color, float amult=1)
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
              line.SetPosition(j++,_traj.r[i]);
              // Draw accelerations
              GameObject obj = new GameObject("Accel");
              thrusts.AddLast(obj);
              LineRenderer line2 = obj.AddComponent<LineRenderer>();
              line2.transform.parent = transform;
              line2.useWorldSpace = false;
              line2.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
              line2.material.color = thrustcol;
              line2.startWidth = 0.4f;
              line2.endWidth = 0.4f;
              line2.positionCount = 2;
              line2.SetPosition(0,_traj.r[i]);
              line2.SetPosition(1,_traj.r[i] + _traj.a[i]*amult);
            }
        }

        public void DrawAlign(Vector3 r_from,Vector3 r_to, Transform transform, Color color)
        {
            if (_align_line == null)
            {
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
            Disable();
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

        public void Disable()
        {
          _enabled = false;
          if (_track_obj != null) {Destroy(_track_obj); _track_obj=null;}
          if (_align_obj != null) {Destroy(_align_obj); _align_obj=null;}
          LinkedListNode<GameObject> node = thrusts.First;
          while(node != null)
          {
            Destroy(node.Value);
            node = node.Next;
          }
          thrusts.Clear();

          LogStop();
          vessel.Autopilot.Disable();
          vessel.OnFlyByWire -= new FlightInputCallback(Fly);
          Events["OnToggle"].guiName = "Enable autopilot";
        }

        ~HopperGuidance()
        {
          Disable();
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Enable autopilot", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void OnToggle()
        {
            if (!_enabled)
            {

                CelestialBody body = vessel.mainBody;
                Vector3d targetPos = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);

                Vector3d[] thrusts;
                Vector3d r0 = vessel.GetWorldPos3D();
                Vector3d v0 = vessel.GetSrfVelocity();
                Vector3d g = FlightGlobals.getGeeForceAtPosition(r0);
                Vector3d vf = new Vector3d(0,0,0);
                _maxThrust = ComputeMaxThrust();
                _startTime = Time.time;
                double amax = (_maxThrust/vessel.totalMass)/_accelFactor;

                Solve solver = new Solve();
                solver.Tmin = 1;
                solver.Tmax = 300; // 5 mins
                solver.tol = 0.5;
                solver.vmax = 50;
                solver.amax = amax*maxPercentThrust*0.01;
                solver.N = 5;
                solver.g = g.magnitude;
                solver.minDescentAngle = minDescentAngle;
                solver.maxThrustAngle = maxThrustAngle;
                solver.maxFinalThrustAngle = maxFinalThrustAngle;

                int retval;
                // Use _logTransform so that Y is vertical direction, and the gravity which acts downwards in the Y direction
                LogSetUpTransform(); // initialises _logTransform
                // Predict into future since solution makes 0.2-2 secs to compute
                r0 = r0 + v0*predictTime + 0.5*g*predictTime*predictTime;
                v0 = v0 + g*predictTime;

                // Shut-off throttle
                FlightCtrlState ctrl = new FlightCtrlState();
                vessel.GetControlState(ctrl);
                ctrl.mainThrottle = 0;

                // Compute trajectory to landing spot
                double fuel;
                double bestT = solver.GoldenSearchGFold(r0, v0, targetPos, vf, _logTransform, out thrusts, out fuel, out retval);
                Debug.Log(solver.DumpString());
                if ((retval>=1) && (retval<=5))
                {
                  double dt = solver.dt;

                  _traj = new Trajectory();
                  _traj.Simulate(bestT, thrusts, r0, v0, g, dt, extendTime);
                  bool corrected = _traj.CorrectFinal(targetPos, Vector3.zero);
                  if (!corrected)
                  {
                    Debug.Log("HopperGuidance: Target not corrected!");
                  }

                  // Enable autopilot
                  _pid3d.Init(kP1,0,0,kP2,0,0,maxErrV,(float)(0.01f*maxErrPercentThrust*amax),2.0f);
                  _targetPos = vessel.mainBody.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
                  DrawTarget(_targetPos, body.transform, targetcol);
                  DrawTrack(_traj, body.transform, trackcol);
                  vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
                  vessel.OnFlyByWire += new FlightInputCallback(Fly);
                  if (_logging) {LogStart("vessel.dat");}
                  // Write solution
                  if (_logging) {_traj.Write("solution.dat",_logTransform);}
                  _enabled = true;
                  Events["OnToggle"].guiName = "Disable autopilot";
                }
                else
                {
                  Disable();
                  Events["OnToggle"].guiName = "Enable fail: Try again";
                  Debug.Log("Failure to find solution: error code "+retval);
                }
            }
            else
            {
                Disable();
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
          if ((vessel == null) || (vessel.state != Vessel.State.ACTIVE) ||
              (vessel.situation == Vessel.Situations.LANDED) || (vessel.situation == Vessel.Situations.SPLASHED))
          {
            Disable();
            return;
          }

          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3d dr,dv,da;
          double desired_t; // closest time in trajectory (desired)
          double dist = _traj.FindClosest(r, v, out dr, out dv, out da, out desired_t);
          if (dist > disengageDistance)
          {
            Debug.Log("Disengaging autopilot as trajectory too distance. "+dist+" > "+disengageDistance);
            Disable();
            return;
          }
          DrawAlign(r,dr,vessel.mainBody.transform,aligncol);
          float throttle=0;
          Vector3 F = vessel.mainBody.transform.up; // default to up
          
          if (da.magnitude > 0.01)
          {
            F = da;
          }

          // Adjust thrust direction when off trajectory
          Vector3d F2 = GetFAdjust(r,v,dr,dv,Time.deltaTime,_logTransform);
          F = F + F2;
          float amax = (float)(_maxThrust/vessel.totalMass);  // F = m.a. We want, unit of throttle for each 1m/s/s
          throttle = _accelFactor*F.magnitude/(amax+0.001f); // protect against divide by zero

          // Shutoff throttle if pointing in wrong direction
          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(att),Vector3d.Normalize(F));
          if ((ddot < Mathf.Cos((float)idleAngle*(Mathf.PI/180.0f))) && (throttle>0.01))
          {
            throttle = 0.001f; // some throttle to steer (if no RCS and main thruster gimbals)
          }

          // Set throttle and direction
          if (F.magnitude > 0.01)
          {
            // Only steer for significant direction
            vessel.Autopilot.SAS.SetTargetOrientation(F,false);
          }
          vessel.Autopilot.SAS.lockedMode = false;
          state.mainThrottle = throttle;

          // Logging
          if (_logging) {LogData(_vesselWriter,Time.time-_startTime,vessel.GetWorldPos3D(),vessel.GetSrfVelocity(),F);}
        }

        public override void OnUpdate()
        {
          base.OnUpdate();
          //_targetPos = vessel.mainBody.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
          //DrawTarget(_targetPos,vessel.mainBody.transform,targetcol);
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Set Target Here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void SetTargetHere()
        {
            // Find vessel co-ordinates
            tgtLatitude = (float)vessel.latitude;
            tgtLongitude = (float)vessel.longitude;
            tgtAltitude = (float)FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D());
            CelestialBody body = vessel.mainBody;
            Vector3d targetPos = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
            DrawTarget(targetPos,body.transform,targetcol);
        }
    }
}
