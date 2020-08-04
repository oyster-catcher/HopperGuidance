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
        GameObject _steer_obj = null;
        LineRenderer _align_line = null;
        LineRenderer _steer_line = null;
        LinkedList<GameObject> thrusts = new LinkedList<GameObject>();
        PID3d _pid3d = new PID3d();
        bool _enabled = false;
        float errMargin = 0.1f; // margin of error in solution to allow for headroom in amax (use double for maxThrustAngle)
        float predictTime = 0.2f;
        double lowestY = 0; // Y position of bottom of craft relative to centre
        float _minThrust, _maxThrust;
        Color trackcol = new Color(0,1,0,0.3f); // transparent green
        Color targetcol = new Color(1,1,0,0.5f); // solid yellow
        Color thrustcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        Color idlecol = new Color(1,0.2f,0.2f,0.9f); // right red (idling as off attitude target)
        Color aligncol = new Color(0,0.1f,1.0f,0.3f); // blue
        Solve solver; // Stores solution inputs, output and trajectory
        Trajectory _traj; // trajectory of solution in world space?
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _transform = null;
        double last_t = -1; // last time data was logged
        double log_interval = 0.05f; // Interval between logging
        System.IO.StreamWriter _tgtWriter = null; // Actual vessel
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        double extendTime = 1; // duration to extend trajectory to slowly descent to touch down and below at touchDownSpeed
        double touchDownSpeed = 1.4f;
        double setTgtLatitude, setTgtLongitude, setTgtAltitude, setTgtSize;
        bool setShowTrack = true;
        float lastTgtLatitude, lastTgtLongitude, lastTgtAltitude;
        bool pickingPositionTarget = false;
        string _vesselLogFilename = "vessel.dat";
        string _tgtLogFilename = "target.dat";
        string _solutionLogFilename = "solution.dat";



        [UI_FloatRange(minValue = 5, maxValue = 500, stepIncrement = 5)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target size", guiFormat = "F0", isPersistant = false)]
        float tgtSize = 10;

        [UI_FloatRange(minValue = -90.0f, maxValue = 90.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Latitude", guiFormat = "F7", isPersistant = true, guiUnits="°")]
        float tgtLatitude = -0.0968071692165f; // H-Pad

        [UI_FloatRange(minValue = -180.0f, maxValue = 180.0f, stepIncrement = 0.0001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Longitude", guiFormat = "F7", isPersistant = true, guiUnits = "°")]
        float tgtLongitude = -74.6172808614f; // H-Pad

        [UI_FloatRange(minValue = 0.1f, maxValue = 10000.0f, stepIncrement = 0.001f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target Altitude", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float tgtAltitude = 176f; // H-Pad

        [UI_FloatRange(minValue = 0.1f, maxValue = 90.0f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min descent angle", guiFormat = "F0", isPersistant = true, guiUnits = "°")]
        float minDescentAngle = 20.0f;

        [UI_FloatRange(minValue = 1, maxValue = 1500, stepIncrement = 10f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max velocity", guiFormat = "F0", isPersistant = true, guiUnits = "m/s")]
        float maxV = 150f; // Max. vel to add to get towards target - not too large that vessel can't turn around

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust angle", guiFormat = "F1", isPersistant = true, guiUnits="°")]
        float maxThrustAngle = 45f; // Max. thrust angle from vertical

        [UI_FloatRange(minValue = 0f, maxValue = 180f, stepIncrement = 1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max landing thrust angle", guiFormat = "F1", isPersistant = true, guiUnits="°")]
        float maxLandingThrustAngle = 5f; // Max. final thrust angle from vertical

        //[UI_FloatRange(minValue = 0, maxValue = 100, stepIncrement = 0.1f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "Min thrust %", guiFormat = "F1", isPersistant = true, guiUnits="%")]
        //float minPercentThrust = 0.01f; // raise for Realism Overhaul to engine doesn't cut out and can steer

        //[UI_FloatRange(minValue = 0, maxValue = 300f, stepIncrement = 5f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "Max thrust %", isPersistant = true, guiUnits = "%")]
        //float maxPercentThrust = 90f;

        [UI_FloatRange(minValue = 0.0f, maxValue = 5.0f, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Time penalty", guiFormat = "F1", isPersistant = false)]
        float timePenalty = 0.0f; // Fuel penalty to every extra second

        [UI_FloatRange(minValue = 0f, maxValue = 90f, stepIncrement = 5f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Idle attitude angle", guiFormat = "F0", isPersistant = true, guiUnits = "°")]
        float idleAngle = 90.0f;


        [UI_FloatRange(minValue = 0.01f, maxValue = 1f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "correction factor", guiFormat = "F2", isPersistant = true)]
        float corrFactor = 0.2f; // If 1 then at 1m error aim to close a 1m/s
        // Set this down to 0.05 for large craft and up to 0.4 for very agile small craft

        [UI_FloatRange(minValue = 0, maxValue = 10, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "accel gain", guiFormat = "F1", isPersistant = true)]
        float accGain = 0.8f;

        //[UI_FloatRange(minValue = 1, maxValue = 5, stepIncrement = 0.1f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "Acceleration mutiplier", guiFormat = "F1", isPersistant = true)]
        //float vMult = 1;
        float yMult = 2; // Extra weight for Y to try and minimise height error over other errors

        float ki1 = 0;
        float ki2 = 0;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Keep engine ignited", isPersistant = false)]
        bool _keepIgnited = false;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Show track", isPersistant = false)]
        bool showTrack = true;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Enable logging", isPersistant = false)]
        bool _logging = false;

        bool _forceLanding = false; // make a soft-landing wherever you are

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
            Debug.Log("HopperGuidance: OnDestroy()");
            Disable();
        }

        public void SetUpTransform(float lat,float lon,float alt)
        {
          // Keep track so we can discover if sliders are moved
          setTgtLatitude = lat;
          setTgtLongitude = lon;
          setTgtAltitude = alt;

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
              part.isEngine(out List<ModuleEngines> engines);
              foreach (ModuleEngines engine in engines)
              {
                  engine.Activate(); // must be active to get thrusts or else realIsp=0
                  //for(float throttle=0; throttle<=1; throttle+=0.1f)
                  //  Debug.Log("isp="+engine.realIsp+" throttle="+throttle+" Thrust="+engine.GetEngineThrust(engine.realIsp,throttle));
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

        public void Enable()
        {
          lowestY = FindLowestPointOnVessel();
          Vector3d[] thrusts;
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
          solver.timePenalty = timePenalty;

          Debug.Log("amin="+amin+" amax="+amax);

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
          Vector3d tr0 = _transform.InverseTransformPoint(r0);
          Vector3d tv0 = _transform.InverseTransformVector(v0);
          double t1 = Time.time;
          double bestT = solver.GoldenSearchGFold(tr0, tv0, rf, vf, out thrusts, out fuel, out retval);
          double t2 = Time.time;
          Debug.Log("GoldenSearchGFold time = "+(t2-t1));
          Debug.Log(solver.DumpString());
          if ((retval>=1) && (retval<=5))
          {
             double dt = solver.dt;
            _traj = new Trajectory(); // Transformed into world space
            // Use g vector from solution calculation
            _traj.Simulate(bestT, thrusts, tr0, tv0, new Vector3d(0,-g.magnitude,0), dt, extendTime);

            // Simulate in world space - to get better precision
            Trajectory traj2 = new Trajectory();
            Vector3d [] thrusts2 = new Vector3d[thrusts.Length];
            for(int i=0;i<thrusts.Length;i++)
              thrusts2[i] = _transform.TransformVector(thrusts[i]);
            traj2.Simulate(bestT, thrusts2, r0, v0, g, dt, extendTime);
            traj2.CorrectFinal(_transform.TransformPoint(rf),_transform.TransformVector(vf));
            DrawTrack(traj2, _transform, trackcol, false);

            double fdist = (_traj.r[_traj.r.Length-1] - rf).magnitude;
            Debug.Log("HopperGuidance: Final pos error = "+fdist);
            _traj.CorrectFinal(rf,vf);

            // Enable autopilot
            //_pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult,vMult);
            _pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult);
            // TODO - Testing out using in solution co-ordinates
            DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            vessel.OnFlyByWire += new FlightInputCallback(Fly);
            // Write solution
            if (_logging) {_traj.Write(_solutionLogFilename);}
            _enabled = true;
            Events["OnToggle"].guiName = "Disable autopilot";
          }
          else
          {
            Disable();
            Events["OnToggle"].guiName = "Enable fail: Try again";
            string msg = "HopperGuidance: Failure to find solution because ";
            // Do some checks
            if (tv0.magnitude > maxV)
              msg = msg + " velocity over "+maxV+" m/s";
            else
            {
              Vector3d r = tr0 - rf;
              double cos_descentAng = Math.Sqrt(r.x*r.x + r.z*r.z) / r.magnitude;
              if (cos_descentAng < Math.Cos(Math.PI*minDescentAngle))
                msg = msg + "below min. descent angle "+minDescentAngle+"°";
              else
                msg = msg + "impossible to reach target within constraints";
            }
            Debug.Log(msg);
            ScreenMessages.PostScreenMessage(msg, 3.0f, ScreenMessageStyle.UPPER_CENTER);
          }
        }

        public void Disable()
        {
          _enabled = false;
          _forceLanding = false;
          if (_track_obj != null) {Destroy(_track_obj); _track_obj=null;}
          if (_tgt_obj != null) {Destroy(_tgt_obj); _tgt_obj=null;}
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


        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Force landing", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void ToggleForceLanding()
        {
          if (!_forceLanding)
          {
            Disable();
            lowestY = FindLowestPointOnVessel();
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            ComputeMinMaxThrust(out _minThrust,out _maxThrust); // This might be including RCS (i.e. non main Throttle)
            float amax = (float)(_maxThrust/vessel.totalMass); // F = m*a. We want, unit of throttle for each 1m/s/s
            //_pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult,vMult);
            _pid3d.Init(corrFactor,ki1,0,accGain,ki2,0,maxV,(float)amax,yMult);
            vessel.OnFlyByWire += new FlightInputCallback(FlyForceLanding);
            Events["ToggleForceLanding"].guiName = "Disable forced landing";
            SetUpTransform((float)vessel.latitude,(float)vessel.longitude,0);
            _forceLanding = true;
          }
          else
          {
            vessel.Autopilot.Disable();
            vessel.OnFlyByWire -= new FlightInputCallback(FlyForceLanding);
            Events["ToggleForceLanding"].guiName = "Force landing";
            _forceLanding = false;
          }
        }

        public void FlyForceLanding(FlightCtrlState state)
        {
          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3 tv = _transform.InverseTransformVector(v);
          Vector3 tr = _transform.InverseTransformPoint(r);
          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          Vector3d tatt = _transform.InverseTransformVector(att);

          Vector3 up = FlightGlobals.getUpAxis(vessel.GetWorldPos3D());
          float h1 = (float)FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D());
          float h2 = (float)(FlightGlobals.getAltitudeAtPos(vessel.GetWorldPos3D() + up*(float)lowestY));
          float height = vessel.GetHeightFromTerrain() - (h1-h2);
          
          // Find height via ray
          RaycastHit hit = new RaycastHit();
          Vector3d origin = vessel.GetWorldPos3D() + up*(float)lowestY;
          bool isHit = Physics.Raycast(origin, -up, out hit, Mathf.Infinity, 1 << 15);
          if (!isHit)
            return;
          // Find height. Vertical distance from bottom of craft to intersection
          float rheight = Vector3.Dot(up,hit.point - origin);
 
          Debug.Log("h(terrain)="+height+" h(surface)="+vessel.GetHeightFromSurface()+" rheight="+rheight+" rdist="+hit.distance);

          if ((vessel == null) || vessel.checkLanded() || (height<0))
          {
            Disable();
            // Shut-off throttle
            FlightCtrlState ctrl = new FlightCtrlState();
            vessel.GetControlState(ctrl);
            ctrl.mainThrottle = 0;
            // Necessary on Realism Overhaul to shutdown engine as at throttle=0 the engine may still have
            // a lot of thrust
            if (_keepIgnited)
              ShutdownAllEngines();
            vessel.Autopilot.Disable();
            vessel.OnFlyByWire -= new FlightInputCallback(FlyForceLanding);
            Events["ToggleForceLanding"].guiName = "Force landing";
            _forceLanding = false;
            return;
          }

          ComputeMinMaxThrust(out _minThrust,out _maxThrust); // This might be including RCS (i.e. non main Throttle)
          if (_minThrust > 0)
            return;
          float amax = (float)(_maxThrust/vessel.totalMass); // F = m*a. We want, unit of throttle for each 1m/s/s
          float amin = (float)(_minThrust/vessel.totalMass);


          // Time to hit surface
          Vector3 g = FlightGlobals.getGeeForceAtPosition(r);

          // Find time to hit ground with just gravity
          // 0 = a*x*x + b*x + c;
          // height = 0.5*amax*t*t + v*t; 
          float a = 0.5f*g.magnitude;
          float b = v.magnitude;
          float c = -height;
          float x1 = (-b + Mathf.Sqrt(b*b-4*a*c))/(2*a);
          float x2 = (-b - Mathf.Sqrt(b*b-4*a*c))/(2*a);
          float t = x1;
          if (x2 > 0)
            t = x2;
          // Final velocity after gravity > max deceleration over time period
          float maxv = t*(amax-g.magnitude);
          Vector3 dv = new Vector3(0,-maxv,0);
          float throttle = (_keepIgnited)?0.01f:0;
          Vector3 tg = new Vector3(0,-g.magnitude,0);
          Vector3 F = -tg; // compensate for gravity
          if (v.magnitude > maxv)
          {
            // Now compute again
            F = -tg;

            // Note: r is not transformed but the same r is used for actual and desired position
            // because we want zero error as we just correct for velocity as we don't worry
            // about where we are landing
            Vector3 F2 = _pid3d.Update(tr,tv,tr,dv,Time.deltaTime);
            F = F + F2;
            // Reduce sideways components
            F = ConeUtils.ClosestThrustInsideCone((float)maxThrustAngle,(float)amin,(float)amax,F);
            if (_keepIgnited)
              throttle = Mathf.Clamp((float)(F.magnitude-amin)/(amax-amin),0.01f,1);
            else
              throttle = Mathf.Clamp((float)(F.magnitude-amin)/(amax-amin),0,1);
            Debug.Log("HopperGuidance: h="+height+" tv="+(Vector3d)tv+" dv="+(Vector3d)dv+" t_hit="+t+" F="+(Vector3d)F+" F2="+(Vector3d)F2+" throttle="+throttle);
          }


          // Shutoff throttle if pointing in wrong direction
          F = _transform.TransformVector(F); // transform back to world co-ordinates
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(att),Vector3d.Normalize(F));
          if ((ddot < Mathf.Cos((float)idleAngle*(Mathf.PI/180.0f))) && (F.magnitude>0.01))
          {
            Debug.Log("HopperGuidance: Attitude wrong");
            throttle = 0.01f; // some throttle to steer? (if no RCS and main thruster gimbals)
            // Draw steer vector
            DrawSteer(tr, tr+3*vessel.vesselSize.y*Vector3d.Normalize(F), _transform, idlecol);
          }
          else
          {
            // Draw steer vector
            DrawSteer(tr, tr+3*vessel.vesselSize.y*Vector3d.Normalize(F), _transform, thrustcol);
          }
          vessel.Autopilot.SAS.lockedMode = false;
          vessel.Autopilot.SAS.SetTargetOrientation(F,false);
          state.mainThrottle = throttle;
        }

        public void Fly(FlightCtrlState state)
        {
          // Get nearest point on trajectory in position and velocity
          // use both so we can handle when the trajectory crosses itself,
          // and it works a bit better anyway if velocity and position are not in step
          if ((vessel == null) || vessel.checkLanded())
          {
            Disable();
            // Shut-off throttle
            FlightCtrlState ctrl = new FlightCtrlState();
            vessel.GetControlState(ctrl);
            ctrl.mainThrottle = 0;
            // Necessary on Realism Overhaul to shutdown engine as at throttle=0 the engine may still have
            // a lot of thrust
            if (_keepIgnited)
              ShutdownAllEngines();
            return;
          }

          // Open log files
          if ((_vesselWriter == null) && (_logging))
          {
            _vesselWriter = new System.IO.StreamWriter(_vesselLogFilename);
            _vesselWriter.WriteLine("time x y z vx vy vz ax ay az att_err");
            _tgtWriter = new System.IO.StreamWriter(_tgtLogFilename);
            _tgtWriter.WriteLine("time x y z vx vy vz ax ay az att_err");
          }

          // Logging
          double t = Time.time - _startTime;

          Vector3 r = vessel.GetWorldPos3D();
          Vector3 v = vessel.GetSrfVelocity();
          Vector3d dr,dv,da;
          double desired_t; // closest time in trajectory (desired)
          // TODO - Transform with _transform
          Vector3d tr = _transform.InverseTransformPoint(r);
          Vector3d tv = _transform.InverseTransformVector(v);
          _traj.FindClosest(tr, tv, out dr, out dv, out da, out desired_t, 0.5f, 0.5f);
          DrawAlign(tr,dr,_transform,aligncol);
          float throttle=0;
          
          // Adjust thrust direction when off trajectory
          Vector3d F2 = _pid3d.Update(tr,tv,dr,dv,Time.deltaTime);
          Vector3d F = da + F2;

          float amax = (float)(_maxThrust/vessel.totalMass); // F = m*a. We want, unit of throttle for each 1m/s/s
          float amin = (float)(_minThrust/vessel.totalMass);

          // Show before limiting - flash if limiting?
          Vector3 sF = F; // save unlimited F

          // Return the closest possible thrust vector to the desired one but that
          // is inside the thrust cone constrainted by maxThrustAngle from vertical
          // and the minimum acceleration amin and max acceleration amax
          // Note that currently of maxThrustAngle > 90 this function has no effect
          F = ConeUtils.ClosestThrustInsideCone((float)maxThrustAngle,(float)amin,(float)amax,F);

          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          Vector3d tatt = _transform.InverseTransformVector(att);
          float ddot = (float)Vector3d.Dot(Vector3d.Normalize(tatt),Vector3d.Normalize(F));
          float att_err = Mathf.Acos(ddot)*180/Mathf.PI;

          if (_keepIgnited)
            throttle = Mathf.Clamp((float)(F.magnitude-amin)/(amax-amin),0.01f,1);
          else
            throttle = Mathf.Clamp((float)(F.magnitude-amin)/(amax-amin),0,1);

          // Shutoff throttle if pointing in wrong direction
          Vector3d logF = F; // report actual thrust vector used in case in idle region
          Debug.Log("vesselSize="+vessel.vesselSize);
          if ((att_err >= idleAngle) && (F.magnitude>0.01))
          {
           throttle = 0.01f; // some throttle to steer? (if no RCS and main thruster gimbals)
           logF = Vector3d.Normalize(F)*amin;
            DrawSteer(tr, tr+3*vessel.vesselSize.y*Vector3d.Normalize(sF), _transform, idlecol);
          }
          else
          {
            DrawSteer(tr, tr+3*vessel.vesselSize.y*Vector3d.Normalize(sF), _transform, thrustcol);
          }

          if ((_logging)&&(t >= last_t+log_interval))
          {
            LogData(_tgtWriter, t, dr, dv, da, 0);
            LogData(_vesselWriter, t, tr, tv, logF, att_err);
            last_t = t;
          }

          F = _transform.TransformVector(F);
          if (throttle > 0)
          {
            // Only steer for significant direction
            vessel.Autopilot.SAS.SetTargetOrientation(F,false);
          }
          vessel.Autopilot.SAS.lockedMode = false;
          state.mainThrottle = throttle;
        }

        public override void OnUpdate()
        {
          base.OnUpdate();
          if ((tgtLatitude != setTgtLatitude) || (tgtLongitude != setTgtLongitude) || (tgtAltitude != setTgtAltitude) || (tgtSize != setTgtSize))
          {
            SetUpTransform(tgtLatitude, tgtLongitude, tgtAltitude);
            DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
            setTgtSize = tgtSize;
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
              SetUpTransform(tgtLatitude, tgtLongitude, tgtAltitude);
              DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
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
              SetUpTransform(tgtLatitude, tgtLongitude, tgtAltitude);
              DrawTarget(Vector3d.zero,_transform,targetcol,tgtSize);
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
          if ((!_logging) && (_vesselWriter != null))
          {
            LogStop();
          }
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Set Target Here", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
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
        }
    }
}
