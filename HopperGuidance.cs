#define LIMIT_ATTITUDE

using System;
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
        double predictTime = 0.5f;
        float _maxThrust;
        Color trackcol = new Color(0,1,0,0.3f); // transparent green
        Color targetcol = new Color(1,1,0,0.5f); // solid yellow
        Color thrustcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        Color aligncol = new Color(0,0.1f,1.0f,0.3f); // blue
        Solve solver; // Stores solution inputs, output and trajectory
        Trajectory _traj; // trajectory of solution in world space?
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _logTransform = null;
        double last_t = -1; // last time data was logged
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        bool _logging = false;
        double extendTime = 2.0f; // extend trajectory to slowly descent to touch down
        double setTgtLatitude, setTgtLongitude, setTgtAltitude, setTgtSize;
        bool setShowTrack = true;
        float lastTgtLatitude, lastTgtLongitude, lastTgtAltitude;
        bool pickingPositionTarget = false;
        double lowestY; // lowest Y part of vessel (recalculated at SetTargetHere(), PickTarget() and Enable())

        [UI_FloatRange(minValue = 5, maxValue = 500, stepIncrement = 5)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target size", guiFormat = "F0", isPersistant = false)]
        float tgtSize = 10;

        [UI_FloatRange(minValue = -90.0f, maxValue = 90.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Latitude", guiFormat = "F7", isPersistant = true)]
        float tgtLatitude = -0.0968071692165f; // H-Pad

        [UI_FloatRange(minValue = -180.0f, maxValue = 180.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Longitude", guiFormat = "F7", isPersistant = true)]
        float tgtLongitude = -74.6172808614f; // H-Pad

        [UI_FloatRange(minValue = 0.1f, maxValue = 10000.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Altitude", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float tgtAltitude = 176f; // H-Pad

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min descent angle", guiFormat = "F0", isPersistant = true, guiUnits = "m")]
        float minDescentAngle = 20.0f;

        [UI_FloatRange(minValue = 1, maxValue = 1500, stepIncrement = 10f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max velocity", guiFormat = "F0", isPersistant = true)]
        float maxV = 150f; // Max. vel to add to get towards target - not too large that vessel can't turn around

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust angle", guiFormat = "F1", isPersistant = true)]
        float maxThrustAngle = 45f; // Max. thrust angle from vertical

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max landing thrust angle", guiFormat = "F1", isPersistant = true)]
        float maxLandingThrustAngle = 5f; // Max. final thrust angle from vertical

        [UI_FloatRange(minValue = 0, maxValue = 100, stepIncrement = 1)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min thrust %", guiFormat = "F0", isPersistant = true)]
        float minPercentThrust = 10;

        [UI_FloatRange(minValue = 0, maxValue = 300f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust %", isPersistant = true, guiUnits = "%")]
        float maxPercentThrust = 100f;

        [UI_FloatRange(minValue = 0.0f, maxValue = 10.0f, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Time penalty", guiFormat = "F1", isPersistant = false)]
        float timePenalty = 0.0f; // Fuel penalty to every extra second

        [UI_FloatRange(minValue = 0f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Idle attitude angle", guiFormat = "F0", isPersistant = true)]
        float idleAngle = 90.0f;


        [UI_FloatRange(minValue = 0.01f, maxValue = 1f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Position gain", guiFormat = "F2", isPersistant = true)]
        float kP1 = 0.3f; // If 1 then at 1m error aim to close a 1m/s

        [UI_FloatRange(minValue = 0.0f, maxValue = 2f, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Velocity gain", guiFormat = "F1", isPersistant = true)]
        float kP2 = 1.0f; // If 1 then at 1m/s error in velocity acceleration at an extra 1m/s/s

        [UI_FloatRange(minValue = 0.0f, maxValue = 90.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Err: Extra thrust angle", guiFormat = "F1", isPersistant = true)]
        float errExtraThrustAngle = 10.0f; // Additional thrust angle from vertical allowed to correct for error

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Show track", isPersistant = false)]
        bool showTrack = true;

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
          //meshr.shadowCastingMode = renderer.shadowCastingMode.Off;

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
              line.SetPosition(j++,transform.TransformPoint(traj.r[i]));
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
              line2.SetPosition(0,transform.TransformPoint(traj.r[i]));
              line2.SetPosition(1,transform.TransformPoint(traj.r[i] + traj.a[i]*amult));
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
            _align_line.startWidth = 0.3f;
            _align_line.endWidth = 0.3f;
            _align_line.positionCount = 2;
            _align_line.SetPosition(0,transform.TransformPoint(r_from));
            _align_line.SetPosition(1,transform.TransformPoint(r_to));
        }

        public void OnDestroy()
        {
            Debug.Log("HopperGuidance: OnDestroy()");
            Disable();
        }

        public void LogSetUpTransform()
        {
          // Keep track so we can discover if sliders are moved
          setTgtLatitude = tgtLatitude;
          setTgtLongitude = tgtLongitude;
          setTgtAltitude = tgtAltitude;

          // Set up transform so Y is up and (0,0,0) is target position
          CelestialBody body = vessel.mainBody;
          Vector3d origin = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude);
          Vector3d vEast = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude-0.1, tgtAltitude) - origin;
          Vector3d vUp = body.GetWorldSurfacePosition(tgtLatitude, tgtLongitude, tgtAltitude+1) - origin;
          // Convert to body co-ordinates
          origin = body.transform.InverseTransformPoint(origin);
          vEast = body.transform.InverseTransformVector(vEast);
          vUp = body.transform.InverseTransformVector(vUp);

          GameObject go = new GameObject();
          // Need to rotation that converts (0,1,0) to vUp in the body transform
          Quaternion quat = Quaternion.FromToRotation(new Vector3(0,1,0),vUp);

          _logTransform = go.transform;
          _logTransform.SetPositionAndRotation(origin,quat);
          _logTransform.SetParent(body.transform,false);
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
          last_t = -1;
        }

        void LogData(System.IO.StreamWriter f, double t, Vector3 r, Vector3 v, Vector3 a)
        {
          if (f != null)
          {
            f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1}",t,r.x,r.y,r.z,v.x,v.y,v.z,a.x,a.y,a.z));
            last_t = t;
          }
        }

        // Find Y offset to lowest part from origin of the vessel
        double FindAltLowestPointOnVessel(out double miny)
        {
          Vector3 CoM, up;

          CoM = vessel.localCoM;
          Vector3 bottom = Vector3.zero; // Offset from CoM
          up = FlightGlobals.getUpAxis(CoM); //Gets up axis
          miny = 0;
          double alt = FlightGlobals.getAltitudeAtPos(CoM);
          foreach (Part p in vessel.parts)
          {
            if (p.collider != null) //Makes sure the part actually has a collider to touch ground
            {
              Vector3 pbottom = p.collider.ClosestPointOnBounds(vessel.mainBody.position); //Gets the bottom point
              double y = Vector3.Dot(up,pbottom);
              if (FlightGlobals.getAltitudeAtPos(bottom) < alt)
                alt = FlightGlobals.getAltitudeAtPos(bottom);
              if (y < miny)
              {
                bottom = pbottom;
                miny = y;
              }
            }
          }
          miny = miny - Vector3.Dot(up,CoM); // Add on offset of CoM
          return alt;
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

        public void Enable()
        {
          lowestY = -vessel.vesselSize.y*0.4; // approximate if CoG at centre of vessel
          //Debug.Log("lowestY="+lowestY);
          Vector3d[] thrusts;
          Vector3d r0 = vessel.GetWorldPos3D();
          Vector3d v0 = vessel.GetSrfVelocity();
          Vector3d g = FlightGlobals.getGeeForceAtPosition(r0);
          Vector3d rf = new Vector3d(0,-lowestY,0); // a little above the surface
          Vector3d vf = new Vector3d(0,-0.1,0);
          _maxThrust = ComputeMaxThrust();
          _startTime = Time.time;
          double amax = _maxThrust/vessel.totalMass;

          solver = new Solve();
          solver.Tmin = 1;
          solver.Tmax = 300; // 5 mins
          solver.tol = 0.1;
          solver.vmax = maxV;
          solver.amax = amax*maxPercentThrust*0.01;
          solver.amin = amax*minPercentThrust*0.01;
          solver.Nmin = 2;
          solver.Nmax = 6;
          solver.minDurationPerThrust = 4;
          solver.g = g.magnitude;
          solver.minDescentAngle = minDescentAngle;
          solver.maxThrustAngle = maxThrustAngle;
          solver.maxLandingThrustAngle = maxLandingThrustAngle;
          solver.timePenalty = timePenalty;

          int retval;
          // Predict into future since solution makes 0.1-0.3 secs to compute
          r0 = r0 + v0*predictTime + 0.5*g*predictTime*predictTime;
          v0 = v0 + g*predictTime;

          // Shut-off throttle
          FlightCtrlState ctrl = new FlightCtrlState();
          vessel.GetControlState(ctrl);
          ctrl.mainThrottle = 0.01f*minPercentThrust;
          vessel.Autopilot.SAS.SetTargetOrientation(-g,false);

          // Compute trajectory to landing spot
          double fuel;
          Vector3d tr0 = _logTransform.InverseTransformPoint(r0);
          Vector3d tv0 = _logTransform.InverseTransformVector(v0);
          double bestT = solver.GoldenSearchGFold(tr0, tv0, rf, vf, out thrusts, out fuel, out retval);
          Debug.Log(solver.DumpString());
          if ((retval>=1) && (retval<=5))
          {
             double dt = solver.dt;
            _traj = new Trajectory(); // Transformed into world space
            // Use g vector from solution calculation
            _traj.Simulate(bestT, thrusts, tr0, tv0, new Vector3d(0,-g.magnitude,0), dt, extendTime);
            double fdist = (_traj.r[_traj.r.Length-1] - rf).magnitude;
            Debug.Log("HopperGuidance: Final pos error = "+fdist);
            _traj.CorrectFinal(rf,vf);

            // Enable autopilot
            _pid3d.Init(kP1,0,0,kP2,0,0,maxV,(float)(0.01f*maxPercentThrust*amax),2.0f);
            // TODO - Testing out using in solution co-ordinates
            DrawTrack(_traj, _logTransform, trackcol);
            DrawTarget(Vector3d.zero,_logTransform,targetcol,tgtSize);
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            vessel.OnFlyByWire += new FlightInputCallback(Fly);
            if (_logging) {LogStart("vessel.dat");}
            // Write solution
            if (_logging) {_traj.Write("solution.dat");}
            _enabled = true;
            Events["OnToggle"].guiName = "Disable autopilot";
          }
          else
          {
            Disable();
            Events["OnToggle"].guiName = "Enable fail: Try again";
            Debug.Log("HopperGuidance: Failure to find solution: error code "+retval);
            string message = "Failure to find solution within constraints - try relaxing them";
            ScreenMessages.PostScreenMessage(message, 3.0f, ScreenMessageStyle.UPPER_CENTER);
          }
        }

        public void Disable()
        {
          _enabled = false;
          if (_track_obj != null) {Destroy(_track_obj); _track_obj=null;}
          if (_tgt_obj != null) {Destroy(_tgt_obj); _tgt_obj=null;}
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

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Enable autopilot", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void OnToggle()
        {
          if (!_enabled)
            Enable();
          else
            Disable();
        }

        // Use 3-D pairs of PIDs to calculate adjustment to force vector and optionally transfrom
        // space where Y is vertical
        public Vector3d GetFAdjust(Vector3d r, Vector3d v, Vector3d dr, Vector3d dv, float dt, Transform transform=null)
        {
            if (transform)
            {
              Vector3d r2 = transform.TransformPoint(r);
              Vector3d v2 = transform.TransformVector(v);
              Vector3d dr2 = transform.TransformPoint(dr);
              Vector3d dv2 = transform.TransformVector(dv);
              Vector3d F2 = _pid3d.Update(r2,
                                          v2,
                                          dr2,
                                          dv2,
                                          Time.deltaTime);
              return transform.InverseTransformVector(F2);
            }
            else
            {
              return _pid3d.Update(r,v,dr,dv,Time.deltaTime);
            }
        }

        public void Fly(FlightCtrlState state)
        {
          //LogSetUpTransform();
          // Get nearest point on trajectory in position and velocity
          // use both so we can handle when the trajectory crosses itself,
          // and it works a bit better anyway if velocity and position are not in step
          if ((vessel == null) || vessel.checkLanded())
          {
            Disable();
            return;
          }

          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3d dr,dv,da;
          double desired_t; // closest time in trajectory (desired)
          // TODO - Transform with _logTransform
          Vector3d tr = _logTransform.InverseTransformPoint(r);
          Vector3d tv = _logTransform.InverseTransformVector(v);
          _traj.FindClosest(tr, tv, out dr, out dv, out da, out desired_t);
          DrawAlign(tr,dr,_logTransform,aligncol);
          float throttle=0;
          
          // Adjust thrust direction when off trajectory
          Vector3d F2 = GetFAdjust(tr,tv,dr,dv,Time.deltaTime,_logTransform);
          Vector3d F = da + F2;

#if (LIMIT_ATTITUDE)
          // Restrict angle from vertical to thrust angle + allowed error
          double maxAngle = maxThrustAngle + errExtraThrustAngle;

          // Calculate normal for plane at maxAngle (make it 2D)
          if (true)
          {
            // TODO: Working upside down?
            // Note tan at 90 degrees goes infinite
            double R = Math.Abs(F.y)*Math.Tan(maxAngle*Math.PI/180); // cone radius at this point
            double x = Math.Sqrt(F.x*F.x + F.z*F.z);
            if (F.y > 0)
            {
              if ((x > R)&&(maxAngle < 90)) // upright and thrust cone upwards
              {
                // Need to reduce fx (horizontal component of thrust until within cone)
                // 0 = nx*fhor2 + ny*F.y;
                // fhor2 = p*Math.Sqrt(F.x*F.x+F.z*F.z)
                // fhor2 = Math.Sqrt(p1*F.x*F.x+p1*F.z*F.z)
                // fhor2*fhor2 = p1*F.x*F.x + p1*F.z*F.z
                // fhor2*Fhor2 = p1*(F.x*F.x + F.z*F.z)
                // p1 = (fhor2*fhor2)/(F.x*F.x + F.z*F.z)
                // need to scale back F.x and F.y to be on boundary
                F.x = (float)((R/x)*F.x);
                F.z = (float)((R/x)*F.z);
                //Debug.Log("Outside of thrust cone - R="+R+" F.x="+F.x+" F.z="+F.z+" scaling to "+F);
              }
              // If F.y > 0 and maxAngle > 90 then by definite we are inside upwards hemisphere
            } else {
              // If F.y < 0 and maxAngle < 90 we are definite outside thrust cone
              if (maxAngle < 90)
              {
                // Project to nearest point inside thrust cone
                //Debug.Log("Zero thrust since in wrong hemisphere");
                F.x = 0;
                F.y = 0;
                F.z = 0;
              }
              if ((x < R)&&(maxAngle > 90)) // Now we need to be outside of thrust cone
              {
                //Debug.Log("Inside of upside down thrust cone -> zero thrust");
                F.x = 0;
                F.y = 0;
                F.z = 0;
              }
            }
          }

#endif
          // Logging
          double t = Time.time - _startTime;
          if ((_logging)&&(t+solver.dt >= last_t))
          {
            LogData(_vesselWriter,t,tr,tv,F); // Updates last_t
          }

          F = _logTransform.TransformVector(F);

          float amax = (float)(_maxThrust/vessel.totalMass); // F = m*a. We want, unit of throttle for each 1m/s/s
          throttle = (float)F.magnitude/(amax+0.001f); // protect against divide by zero

          // Shutoff throttle if pointing in wrong direction
          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(att),Vector3d.Normalize(F));
          if ((ddot < Mathf.Cos((float)idleAngle*(Mathf.PI/180.0f))) && (F.magnitude>0.01))
          {
            throttle = 0.01f*minPercentThrust; // some throttle to steer? (if no RCS and main thruster gimbals)
          }

          // Set throttle and direction
          if (throttle > 0)
          {
            // Only steer for significant direction
            vessel.Autopilot.SAS.SetTargetOrientation(F,false);
          }
          vessel.Autopilot.SAS.lockedMode = false;
          state.mainThrottle = Mathf.Max(throttle,0.01f*minPercentThrust);

        }

        public override void OnUpdate()
        {
          base.OnUpdate();
          if ((tgtLatitude != setTgtLatitude) || (tgtLongitude != setTgtLongitude) || (tgtAltitude != setTgtAltitude) || (tgtSize != setTgtSize))
          {
            LogSetUpTransform();
            DrawTarget(Vector3d.zero,_logTransform,targetcol,tgtSize);
            setTgtSize = tgtSize;
          }
          if (showTrack != setShowTrack)
          {
            DrawTrack(_traj, _logTransform, trackcol);
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
              LogSetUpTransform();
              DrawTarget(Vector3d.zero,_logTransform,targetcol,tgtSize);
              pickingPositionTarget = false;
              return;
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
              LogSetUpTransform();
              DrawTarget(Vector3d.zero,_logTransform,targetcol,tgtSize);
              // TODO - Recompute trajectory is autopilot on

              // If clicked stop picking
              if (Input.GetMouseButtonDown(0))
              {
                pickingPositionTarget = false;
                if (_enabled)
                  Enable(); // recalculate trajectory
                pickingPositionTarget = false;
              }
            }
          }
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Set Target Here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void SetTargetHere()
        {
            // Fire ray down to ground
            //Ray ray = FlightCamera.fetch.mainCamera.ScreenPointToRay(Input.mousePosition);
            //Vector3d g = FlightGlobals.getGeeForceAtPosition(vessel.GetWorldPos3D());
            // Hit ground within 50 units
            //RaycastHit hit;
            //bool isHit = Physics.Raycast(vessel.GetWorldPos3D(), g, out hit, 50, 1 << 15);

            // Find vessel co-ordinates
            tgtLatitude = (float)vessel.latitude;
            tgtLongitude = (float)vessel.longitude;
            // Note: compensate of height of vessel to get height at bottom of vessel
            tgtAltitude = (float)(FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D() - vessel.transform.up*vessel.vesselSize.y*0.5f ));
            LogSetUpTransform();
            DrawTarget(Vector3d.zero,_logTransform,targetcol,tgtSize);
        }
    }
}
