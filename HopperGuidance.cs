using System;
using System.Collections.Generic;
using System.Reflection;

using UnityEngine;
using UnityEngine.UI;

using KSPAssets;

namespace HopperGuidance
{
    enum AutoMode
    {
      Off,
      LandAtTarget,
      Failed
    }

    public class Target
    {
      public Target(float a_lat,float a_lon,float a_alt,float a_height)
      {
        lat = a_lat;
        lon = a_lon;
        alt = a_alt;
        height = a_height;
      }
      public float lat,lon,alt,height;
    }

    public class HopperGuidance : PartModule
    {
        // Constants
        static Color trackcol = new Color(0,1,0,0.3f); // transparent green
        static Color targetcol = new Color(1,1,0,0.5f); // solid yellow
        static Color thrustcol = new Color(1,0.2f,0.2f,0.3f); // transparent red
        static Color aligncol = new Color(0,0.1f,1.0f,0.3f); // blue

        List<GameObject> _tgt_objs = new List<GameObject>(); // new GameObject("Target");
        GameObject _track_obj = null;
        GameObject _thrusts_obj = null;
        GameObject _align_obj = null;
        GameObject _steer_obj = null;
        LineRenderer _align_line = null; // so it can be updated
        LineRenderer _steer_line = null; // so it can be updated
        bool checkingLanded = false; // only check once in flight to avoid failure to start when already on ground
        AutoMode autoMode = AutoMode.Off;
        Controller controller;
        float errMargin = 0.1f; // margin of error in solution to allow for headroom in amax (use double for maxThrustAngle)
        double lowestY = 0; // Y position of bottom of craft relative to centre
        float _minThrust, _maxThrust;
        Solve solver; // Stores solution inputs, output and trajectory
        SolveResult _result;
        Trajectory _traj = null; // trajectory of solution in local space
        double _startTime = 0; // Start solution starts to normalize vessel times to start at 0
        Transform _transform = null;
        double last_t = -1; // last time data was logged
        double log_interval = 0.05f; // Interval between logging
        System.IO.StreamWriter _tgtWriter = null; // Actual vessel
        System.IO.StreamWriter _vesselWriter = null; // Actual vessel
        float extendTime = 0; // duration to extend trajectory to slowly descent to touch down and below at touchdownSpeed
        bool pickingPositionTarget = false;
        string _vesselLogFilename = "vessel.dat";
        string _tgtLogFilename = "target.dat";
        string _solutionLogFilename = "solution.dat";

        [UI_FloatRange(minValue = 0.1f, maxValue = 5, stepIncrement = 0.1f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Touchdown speed", guiFormat = "F0", isPersistant = false, guiUnits = "m/s")]
        float touchdownSpeed = 2.5f;

        [UI_FloatRange(minValue = 0, maxValue = 100, stepIncrement = 1)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Final descent distance", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float finalDescentDistance = 10;

        // Specials for Realism Overhaul
        List<ModuleEngines> allEngines = new List<ModuleEngines>();
        List<ModuleEngines> shutdownEngines = new List<ModuleEngines>();
        bool finalDescentLogged = false;

        [UI_FloatRange(minValue = 5, maxValue = 500, stepIncrement = 5)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target size", guiFormat = "F0", isPersistant = false)]
        float tgtSize = 10;
        float setTgtSize;

        List<Target> _tgts = new List<Target>();

        [UI_FloatRange(minValue = 0, maxValue = 1000, stepIncrement = 1)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Target height", guiFormat = "F1", isPersistant = true, guiUnits = "m")]
        float tgtHeight;
        float setTgtHeight;

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

        // This factor just affects X and Z PIDs (horizontal)
        [UI_FloatRange(minValue = 0.1f, maxValue = 0.5f, stepIncrement = 0.01f)]
        [KSPField(guiActive = true, guiActiveEditor = true, guiName = "Correction factor", guiFormat = "F2", isPersistant = true)]
        // Set this down to 0.05 for large craft and up to 0.4 for very agile small craft
        float corrFactor = 0.25f; // If 1 then at 1m error aim to close a 1m/s
        float setCorrFactor;

        // Fixed values for Y PIDs (vertical)
        // - these must be quite high otherwise we can get into hovering above landing spot
        //   if vessel has slow too much but trajectory says high thrust required
        float kp1y = 0.5f, ki1y = 0, kd1y = 0;
        float kp2y = 0.6f, ki2y = 0, kd2y = 0;

        //[UI_FloatRange(minValue = 0.01f, maxValue = 1.0f, stepIncrement = 0.01f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "Accel factor", guiFormat = "F2", isPersistant = true)]
        //float accelFactor = 0.6f; // Gain on thrust vector
        float corrFactorP1 = 0.15f, corrFactorP2 = 0.4f;
        float accelFactorP1 = 0.2f, accelFactorP2 = 0.6f;

        //[UI_FloatRange(minValue = 0.00f, maxValue = 2.0f, stepIncrement = 0.1f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "kd1", guiFormat = "F1", isPersistant = true)]
        float kd1 = 0;
        //float setKd1;

        //[UI_FloatRange(minValue = 0.00f, maxValue = 2.0f, stepIncrement = 0.1f)]
        //[KSPField(guiActive = true, guiActiveEditor = true, guiName = "kd2", guiFormat = "F1", isPersistant = true)]
        float kd2 = 0;
        //float setKd2;

        // Can't get any improvement by raising these coeffs from zero
        float ki1 = 0;
        float ki2 = 0;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Show track", isPersistant = false)]
        bool showTrack = true;
        bool setShowTrack = true;

        [UI_Toggle(disabledText = "Off", enabledText = "On")]
        [KSPField(guiActive = true, guiActiveEditor = false, guiName = "Logging", isPersistant = false)]
        bool _logging = false;

        void Log(string msg, bool onScreen=false)
        {
          Debug.Log("[HopperGuidance] "+msg);
          if (onScreen)
            ScreenMessages.PostScreenMessage(msg, 3.0f, ScreenMessageStyle.UPPER_CENTER);
        }

        // Quad should be described a,b,c,d in anti-clockwise order when looking at it
        void AddQuad(Vector3[] vertices, ref int vi, int[] triangles, ref int ti,
                      Vector3d a, Vector3d b, Vector3d c, Vector3d d,
                      bool double_sided = false)
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
          if (double_sided)
          {
            triangles[ti++] = vi;
            triangles[ti++] = vi+1;
            triangles[ti++] = vi+2;
            triangles[ti++] = vi;
            triangles[ti++] = vi+2;
            triangles[ti++] = vi+3;
          }
          vi += 4;
        }

        void AddLine(Vector3[] vertices, ref int vi, int[] triangles, ref int ti,
                            Vector3d a, Vector3d b, float width,
                            bool double_sided = false)
        {
          Vector3 Y = new Vector3(0,1,0);
          Vector3 C = Vector3.Cross(b-a,Y);
          Vector3 B = Vector3.Cross(b-a,C);
          if (B.magnitude < 0.01f)
          {
            Vector3 X = new Vector3(1,0,0);
            C = Vector3.Cross(b-a,X);
            B = Vector3.Cross(b-a,C);
          }
          B = B/B.magnitude * width;
          C = C/C.magnitude * width;
          AddQuad(vertices,ref vi,triangles,ref ti,a-C,a+C,b+C,b-C,double_sided);
          AddQuad(vertices,ref vi,triangles,ref ti,a-B,a+B,b+B,b-B,double_sided);
        }

        double MinHeightAtMinThrust(double y, double vy,double amin,double g)
        {
          double minHeight = 0;
          if (amin < g)
            return -float.MaxValue;
          double tHover = -vy/amin; // time to come to hover
          minHeight = y + vy*tHover + 0.5*amin*tHover*tHover - 0.5*g*tHover*tHover;
          return minHeight;
        }

        // pos is ground position, but draw up to height
        void DrawTarget(Vector3d pos, Transform a_transform, Color color, double size, float height)
        {
          double[] r = new double[]{size*0.5,size*0.55,size*0.95,size};
          Vector3d gpos = pos;
          Vector3d tpos = pos + new Vector3d(0,height,0);

          Vector3d vx = new Vector3d(1,0,0);
          Vector3d vz = new Vector3d(0,0,1);

          GameObject o = new GameObject();
          o.transform.SetParent(a_transform, false);
          _tgt_objs.Add(o);
          MeshFilter meshf = o.AddComponent<MeshFilter>();
          MeshRenderer meshr = o.AddComponent<MeshRenderer>();
          meshr.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          meshr.material.color = color;
          meshr.receiveShadows = false;

          Mesh mesh = new Mesh();
          Vector3 [] vertices = new Vector3[36*4+4+4+4+8];
          int [] triangles = new int[(36*2*2-8+2+2+2+4+4)*3]; // take away gaps
          int i,j;
          int v=0,t=0;
          for(j=0;j<4;j++) // four concentric rings
          {
            for(i=0;i<36;i++)
            {
              float a = -(i*10)*Mathf.PI/180.0f;
              vertices[v++] = gpos + vx*Mathf.Sin(a)*r[j] + vz*Mathf.Cos(a)*r[j];
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
          AddQuad(vertices,ref v,triangles,ref t,
                  tpos-cx*cs-cz,tpos+cx*cs-cz,tpos+cx*cs+cz,tpos-cx*cs+cz);
          // One side
          AddQuad(vertices,ref v,triangles,ref t,
                  tpos-cx+cz,tpos+cx+cz,tpos+cx+cz*cs,tpos-cx+cz*cs);
          // Other size
          AddQuad(vertices,ref v,triangles,ref t,
                  tpos-cx-cz*cs,tpos+cx-cz*cs,tpos+cx-cz,tpos-cx-cz);

          // Draw quads from cross at actual height to the rings on the ground
          cx = vx*size*0.01;
          cz = vz*size*0.01;
          AddQuad(vertices,ref v,triangles,ref t,
                  gpos-cx,gpos+cx,tpos+cx,tpos-cx,true);
          AddQuad(vertices,ref v,triangles,ref t,
                  gpos-cz,gpos+cz,tpos+cz,tpos-cz,true);

          mesh.vertices = vertices;
          mesh.triangles = triangles;
          meshf.mesh = mesh;
          mesh.RecalculateNormals();
        }

        void DrawTargets(List<Target> tgts, Transform a_transform, Color color, double size)
        {
          foreach (GameObject obj in _tgt_objs)
          {
            if (obj != null )
              Destroy(obj);
          }
          _tgt_objs.Clear();
          foreach (Target t in tgts)
          {
            Vector3d pos = vessel.mainBody.GetWorldSurfacePosition(t.lat, t.lon, t.alt);
            pos = a_transform.InverseTransformPoint(pos); // convert to local (for orientation)
            DrawTarget(pos,a_transform,color,size,t.height);
          }
        }

        void DrawTrack(Trajectory traj, Transform a_transform, float amult=1)
        {
          if (traj == null)
            return;
          if (_track_obj != null)
          {
            Destroy(_track_obj); // delete old track
            _track_obj = null;
          }
          if (_thrusts_obj != null)
          {
            Destroy(_thrusts_obj); // delete old track
            _thrusts_obj = null;
          }
          if (!showTrack)
            return;

          // Track
          _track_obj = new GameObject("Track");
          _track_obj.transform.SetParent(a_transform, false);
          MeshFilter meshf = _track_obj.AddComponent<MeshFilter>();
          MeshRenderer meshr = _track_obj.AddComponent<MeshRenderer>();
          meshr.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          meshr.material.color = trackcol;
          meshr.receiveShadows = false;
          meshr.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
          Mesh mesh = new Mesh();
          Vector3 [] vertices = new Vector3[(traj.Length()-1)*4*2];
          int [] triangles = new int[(traj.Length()-1)*4*2*3]; // number of vertices in tris
          float lineWidth = 0.2f;
          int v=0,t=0;
          for (int i = 0; i < traj.Length()-1; i++)
          {
            Vector3 p1 = traj.r[i];
            Vector3 p2 = traj.r[i+1];
            AddLine(vertices,ref v,triangles,ref t,p1,p2,lineWidth,true);
          }
          mesh.vertices = vertices;
          mesh.triangles = triangles;
          meshf.mesh = mesh;
          mesh.RecalculateNormals();

          // Thrust vectors
          _thrusts_obj = new GameObject("Thrusts");
          _thrusts_obj.transform.SetParent(a_transform, false);
          meshf = _thrusts_obj.AddComponent<MeshFilter>();
          meshr = _thrusts_obj.AddComponent<MeshRenderer>();
          mesh = new Mesh();
          meshr.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
          meshr.material.color = thrustcol;
          meshr.receiveShadows = false;
          meshr.shadowCastingMode = UnityEngine.Rendering.ShadowCastingMode.Off;
          vertices = new Vector3[traj.Length()*4*2];
          triangles = new int[traj.Length()*4*2*3]; // number of vertices in tris
          v=0;
          t=0;
          for (int i = 0; i < traj.Length(); i++)
          {
            Vector3 p1 = traj.r[i];
            Vector3 p2 = traj.r[i] + traj.a[i]*amult;
            AddLine(vertices,ref v,triangles,ref t,p1,p2,lineWidth,true);
          }
          mesh.vertices = vertices;
          mesh.triangles = triangles;
          meshf.mesh = mesh;
          mesh.RecalculateNormals();
        }

        void DrawAlign(Vector3 r_from,Vector3 r_to, Transform a_transform, Color color)
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
            _align_line.transform.parent = a_transform;
            _align_line.useWorldSpace = true;
            _align_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _align_line.material.color = color;
            _align_line.startWidth = 0.3f;
            _align_line.endWidth = 0.3f;
            _align_line.positionCount = 2;
            _align_line.SetPosition(0,a_transform.TransformPoint(r_from));
            _align_line.SetPosition(1,a_transform.TransformPoint(r_to));
        }

        void DrawSteer(Vector3 r_from,Vector3 r_to, Transform a_transform, Color color)
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
            _steer_line.transform.parent = a_transform;
            _steer_line.useWorldSpace = true;
            _steer_line.material = new Material(Shader.Find("KSP/Alpha/Unlit Transparent"));
            _steer_line.material.color = color;
            _steer_line.startWidth = 0.3f;
            _steer_line.endWidth = 0.3f;
            _steer_line.positionCount = 2;
            _steer_line.SetPosition(0,a_transform.TransformPoint(r_from));
            _steer_line.SetPosition(1,a_transform.TransformPoint(r_to));
        }

        public void OnDestroy()
        {
          DisableLand();
          // Remove targets as they are not removed on DisableLand()
          foreach (GameObject obj in _tgt_objs)
          {
            if (obj != null)
              Destroy(obj);
          }
          _tgt_objs.Clear();
        }

        Transform SetUpTransform(Target final)
        {
          // Set up transform so Y is up and (0,0,0) is target position
          CelestialBody body = vessel.mainBody;
          Vector3d origin = body.GetWorldSurfacePosition(final.lat, final.lon, final.alt);
          Vector3d vEast = body.GetWorldSurfacePosition(final.lat, final.lon-0.1, final.alt) - origin;
          Vector3d vUp = body.GetWorldSurfacePosition(final.lat, final.lon, final.alt+1) - origin;
          // Convert to body co-ordinates
          origin = body.transform.InverseTransformPoint(origin);
          vEast = body.transform.InverseTransformVector(vEast);
          vUp = body.transform.InverseTransformVector(vUp);

          GameObject go = new GameObject();
          // Need to rotation that converts (0,1,0) to vUp in the body transform
          Quaternion quat = Quaternion.FromToRotation(new Vector3(0,1,0),vUp);

          Transform o_transform = go.transform;
          o_transform.SetPositionAndRotation(origin,quat);
          o_transform.SetParent(body.transform,false);
          return o_transform;
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

        static void LogData(System.IO.StreamWriter f, double t, Vector3 r, Vector3 v, Vector3 a, float att_err, double amin, double amax)
        {
          f.WriteLine(string.Format("{0} {1:F5} {2:F5} {3:F5} {4:F5} {5:F5} {6:F5} {7:F1} {8:F1} {9:F1} {10:F1} {11:F2} {12:F2}",t,r.x,r.y,r.z,v.x,v.y,v.z,a.x,a.y,a.z,att_err,amin,amax));
        }

        static void LogSolution(SolveResult result, Trajectory traj, string filename)
        {
          List<string> comments = new List<string>();
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
          foreach(var sln in result.solutions)
          {
            sln_Ts.Add(sln.T);
            sln_fuels.Add(sln.fuel);
          }
          comments.Add("solution_time="+String.Join(",",sln_Ts));
          comments.Add("solution_fuel="+String.Join(",",sln_fuels));
          traj.Write(filename, comments);
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

        int ComputeMinMaxThrust(out float minThrust, out float maxThrust, bool log=false)
        {
          allEngines.Clear();
          int numEngines = 0;
          minThrust = 0;
          maxThrust = 0;
          foreach (Part part in vessel.parts)
          {
              part.isEngine(out List<ModuleEngines> engines);
              foreach (ModuleEngines engine in engines)
              {
                  Vector3 relpos = vessel.transform.InverseTransformPoint(part.transform.position);
                  float isp = (engine.realIsp>0)? engine.realIsp : 280; // guess!
                  float pressure = (float)FlightGlobals.getStaticPressure()*0.01f; // so 1.0 at Kerbin sea level? 
                  float atmMaxThrust = engine.MaxThrustOutputAtm(true, true, pressure, FlightGlobals.getExternalTemperature());
                  if (log)
                    Log("  engine="+engine+" relpos="+relpos+" isp="+isp+" MinThrust="+engine.GetEngineThrust(isp,0)+" MaxThrust="+atmMaxThrust+" operational="+engine.isOperational);
                  if (engine.isOperational)
                  {
                      minThrust += engine.GetEngineThrust(isp, 0); // can't get atmMinThrust (this ignore throttle limiting but thats ok)
                      maxThrust += atmMaxThrust; // this uses throttle limiting and should give vac thrust as pressure/temp specified too
                      allEngines.Add(engine);
                      numEngines++;
                  }
              }
          }
          return numEngines;
        }

        // Compute engine thrust if one set of symmetrical engines is shutdown
        // (primarily for a Falcon 9 landing to shutdown engines for slow touchdown)
        List<ModuleEngines> ShutdownOuterEngines(float desiredThrust, bool log=false)
        {
          List<ModuleEngines> shutdown = new List<ModuleEngines>();

          // Find engine parts and sort by closest to centre first
          List<Tuple<double,ModuleEngines>> allEngines = new List<Tuple<double,ModuleEngines>>();
          foreach (Part part in vessel.GetActiveParts())
          {
            Vector3 relpos = vessel.transform.InverseTransformPoint(part.transform.position);
            part.isEngine(out List<ModuleEngines> engines);
            double dist = Math.Sqrt(relpos.x*relpos.x+relpos.z*relpos.z);
            foreach (ModuleEngines engine in engines)
              allEngines.Add(new Tuple<double,ModuleEngines>(dist,engine));
          }
          allEngines.Sort();

          // Loop through engines starting a closest to axis
          // Accumulate minThrust, once minThrust exceeds desiredThrust shutdown this and all
          // further out engines
          float minThrust=0, maxThrust=0;
          double shutdownDist = float.MaxValue;
          foreach (var engDist in allEngines)
          {
            ModuleEngines engine = engDist.Item2;
            if (engine.isOperational)
            {
              minThrust += engine.GetEngineThrust(engine.realIsp, 0);
              maxThrust += engine.GetEngineThrust(engine.realIsp, 1);
              if (shutdownDist == float.MaxValue)
              {
                if ((minThrust < desiredThrust) && (desiredThrust < maxThrust)) // good amount of thrust
                  shutdownDist = engDist.Item1 + 0.1f;
                if (minThrust > desiredThrust)
                  shutdownDist = engDist.Item1 - 0.1f;
              }

              if (engDist.Item1 > shutdownDist)
              {
                if (log)
                  Log("ComputeShutdownMinMaxThrust(): minThrust="+minThrust+" desiredThrust="+desiredThrust+" SHUTDOWN");
                engine.Shutdown();
                shutdown.Add(engine);
              }
              else
                if (log)
                  Log("ComputeShutdownMinMaxThrust(): minThrust="+minThrust+" desiredThrust="+desiredThrust+" KEEP");
            }
          }
          Debug.Log(shutdown.Count+" engines shutdown");
          return shutdown;
        }

        float GetCurrentThrust()
        {
          float thrust = 0;
          foreach (ModuleEngines engine in allEngines)
            thrust += engine.GetCurrentThrust();
          return thrust;
        }


        float AxisTime(float dist, float amin, float amax)
        {
          float t;
          if (dist > 0)
            t = Mathf.Sqrt(dist/amax);
          else
            t = Mathf.Sqrt(dist/amin);
          return 2*t;
        }


        float EstimateTimeBetweenTargets(Vector3d r0, Vector3d v0, List<SolveTarget> tgts, float amax, float g, float vmax)
        {
          if (g > amax)
            return -1; // can't even hover
          float t=0;
          // Estimate time to go from stationary at one target to stationary at next, to provide
          // an upper estimate on the solution time
          float xmax = amax - g;
          float xmin = -(amax-g);
          float ymin = -g;
          float ymax = amax - g;
          float zmax = amax - g;
          float zmin = -(amax-g);

          // Compute position with zero velocity
          t = (float)v0.magnitude / (amax-g);
       
          Vector3d ca = -(amax-g) * v0/v0.magnitude; 
          r0 = r0 + v0*t + 0.5*ca*t*t;

          // r0 and v0 represent stationary after velocity cancelled out
          v0 = Vector3d.zero;

          foreach( SolveTarget tgt in tgts )
          {
            float dx = (float)(tgt.r.x - r0.x);
            float dy = (float)(tgt.r.y - r0.y);
            float dz = (float)(tgt.r.z - r0.z);
            // Compute time to move in each orthogonal axis
            float tx = AxisTime(dx, xmin, xmax);
            float ty = AxisTime(dy, ymin, ymax);
            float tz = AxisTime(dz, zmin, zmax);
            t = t + tx + ty + tz;
          }
          return t;
        }


        public void EnableLandAtTarget()
        {

          Type type = Type.GetType("Mono.Runtime");
          if (type != null)
          {
             MethodInfo displayName = type.GetMethod("GetDisplayName", BindingFlags.NonPublic | BindingFlags.Static);
             if (displayName != null)
               Log("Mono version: "+displayName.Invoke(null, null));
          }

          if (_tgts.Count == 0)
          {
            ScreenMessages.PostScreenMessage("No targets - adding one on the ground directly below", 3.0f, ScreenMessageStyle.UPPER_CENTER);
            double alt = vessel.mainBody.TerrainAltitude(vessel.latitude, vessel.longitude);
            Target tgt = new Target((float)vessel.latitude,(float)vessel.longitude,(float)alt,0);
            _tgts.Add(tgt);
            _transform = SetUpTransform(_tgts[_tgts.Count-1]);
            DrawTargets(_tgts,_transform,targetcol,tgtSize);
          }
          // Reset
          finalDescentLogged = false;
          shutdownEngines.Clear();
          checkingLanded = false; // stop trajectory being cancelled while on ground
          lowestY = FindLowestPointOnVessel();
          Vector3d r0 = vessel.GetWorldPos3D();
          Vector3d v0 = vessel.GetSrfVelocity();
          Vector3d g = FlightGlobals.getGeeForceAtPosition(r0);
          Vector3d vf = new Vector3d(0,-touchdownSpeed,0);

          ComputeMinMaxThrust(out _minThrust,out _maxThrust, true);
          if( _maxThrust == 0 )
          {
            Log("No engine thrust (activate engines)", true);
            autoMode = AutoMode.Off;
            return;
          }
          _startTime = Time.time;
          // Shutdown amin/amax will be recomputed in Fly()
          double amin = _minThrust/vessel.totalMass;
          double amax = _maxThrust/vessel.totalMass;
          float accelFactor  = HGUtils.LinearMap(corrFactor, corrFactorP1, corrFactorP2, accelFactorP1, accelFactorP2);
          Log("accelFactor="+accelFactor);
          controller = new Controller(corrFactor,ki1,kd1,accelFactor,ki2,kd2,kp1y,ki1y,kd1y,kp2y,ki2y,kd2y,(float)amin,(float)amax,maxThrustAngle);
          controller.touchdownSpeed = touchdownSpeed;
          controller.finalDescentDistance = finalDescentDistance;
          controller.attitudeTimeConstant = 0.2f;
          controller.throttleTimeConstant = 0.2f;

          if( amin > amax*0.95 )
          {
            Log("Engine doesn't appear to be throttleable. This makes precision guidance impossible", true);
            autoMode = AutoMode.Off;
            return;
          }

          solver = new Solve();
          solver.Tmin = -1;
          solver.tol = 0.5;
          solver.vmax = maxV;
          solver.amin = amin*(1+errMargin);
          solver.amax = amax*(1-errMargin);
          solver.minDurationPerThrust = 4;
          solver.maxThrustsBetweenTargets = 3;
          solver.g = g.magnitude;
          solver.minDescentAngle = minDescentAngle;
          solver.maxThrustAngle = maxThrustAngle*(1-2*errMargin);
          solver.maxLandingThrustAngle = Math.Max(5,0.1f*maxThrustAngle); // 10% of max thrust angle down to 10 degrees
          solver.TstartMin = 0;
          solver.TstartMax = 0;
          if (amin > 0)
            solver.TstartMax = 20; //start engines up to 20 secs later
          // hack for large craft to allow extra slowdown time at target to prepare for next target
          // where thrust is just over gravity give 5 seconds extra time
          // where thrust is double gravity than use 0.5 secs extra time
          solver.extraTime = (float)(2.5 - 2 * Math.Min(0.5*(amax/g.magnitude),1));

          // Compute trajectory to landing spot
          List<SolveTarget> targets = new List<SolveTarget>();
          Vector3d tr0 = _transform.InverseTransformPoint(r0);
          tr0.y += 0.1f; // move up slightly to ensure above ground plane
          Vector3d tv0 = _transform.InverseTransformVector(v0);

          // Create list of solve targets
          double d = 0;
          Vector3 cr = tr0;
          Vector3d wrf = Vector3d.zero;
          Vector3d wvf = Vector3d.zero;
          for(int i=0; i<_tgts.Count; i++)
          {
            SolveTarget tgt = new SolveTarget();
            Vector3d pos = vessel.mainBody.GetWorldSurfacePosition(_tgts[i].lat, _tgts[i].lon, _tgts[i].alt + _tgts[i].height);
            tgt.r = _transform.InverseTransformPoint(pos); // convert to local (for orientation)
            tgt.raxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
            if (i==_tgts.Count-1) // final target
            {
              tgt.r.y += - lowestY;
              wrf = _transform.TransformPoint(tgt.r);
              tgt.vaxes = SolveTarget.X | SolveTarget.Y | SolveTarget.Z;
              tgt.v = vf;
              wvf = _transform.TransformVector(tgt.v); // to correct final later
            }
            tgt.t = -1;
            d = d + (cr-tgt.r).magnitude;
            cr = tgt.r;
            targets.Add(tgt);
          }

          solver.Tmax = -1; // Forces estimation given initial position, velocity and targets

          solver.apex = targets[targets.Count-1].r;
          _traj = new Trajectory();

          _result = MainProg.MultiPartSolve(ref solver, ref _traj, tr0, tv0, ref targets, (float)g.magnitude, extendTime, true);
          Log(solver.DumpString()+" "+_result.DumpString());
          if (_result.isSolved()) // solved for complete path?
          {
            string msg = String.Format("Found solution T={0:F1} Tstart={1:F1} Fuel={2:F1}",_result.T,_result.Tstart,_result.fuel);
            Log(msg, true);
            // Enable autopilot
            // TODO: Use New functions on Controller()
            DrawTargets(_tgts,_transform,targetcol,tgtSize);
            vessel.Autopilot.Enable(VesselAutopilot.AutopilotMode.StabilityAssist);
            vessel.OnFlyByWire += new FlightInputCallback(Fly);
            // Write solution
            if (_logging)
              LogSolution(_result, _traj, _solutionLogFilename);
            autoMode = AutoMode.LandAtTarget;
            Events["ToggleGuidance"].guiName = "Cancel guidance";
          }
          else
          {
            DisableLand();
            Events["ToggleGuidance"].guiName = "Failed! - Cancel guidance";
            string msg = "Failure to find solution as ";
            double minHeight = MinHeightAtMinThrust(tr0.y, tv0.y, amin, g.magnitude);
            // Do some checks
            if (v0.magnitude > maxV)
              msg = msg + " velocity over "+maxV+" m/s";
            else if (amax < g.magnitude)
              msg = msg + "engine has insufficient thrust, no engines active or no fuel";
            else if (minHeight > 0)
              msg = msg + "can't throttle engine low enough to descend to ground";
            else
              msg = msg + "impossible to reach target within constraints";
            Log(msg, true);
            autoMode = AutoMode.Failed;
          }
          if (_result.isSolved()) // solved for complete path? - show partial?
            DrawTrack(_traj, _transform);
        }

        public void DisableLand()
        {
          autoMode = AutoMode.Off;
          if (_track_obj != null)   {Destroy(_track_obj); _track_obj=null;}
          if (_thrusts_obj != null) {Destroy(_thrusts_obj); _thrusts_obj=null;}
          if (_align_obj != null)   {Destroy(_align_obj); _align_obj=null;}
          if (_steer_obj != null)   {Destroy(_steer_obj); _steer_obj=null;}
          _traj = null;
          LogStop();
          if (vessel != null)
          {
            if (vessel.OnFlyByWire != null)
              vessel.OnFlyByWire -= new FlightInputCallback(Fly);
            if (vessel.Autopilot != null)
              vessel.Autopilot.Disable();
          }
          Events["ToggleGuidance"].guiName = "Enable guidance";
        }

        ~HopperGuidance()
        {
          DisableLand();
        }

        public override void OnInactive()
        {
          base.OnInactive();
          DisableLand();
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Clear targets", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        public void ClearTargets()
        {
          _tgts.Clear();
          DrawTargets(_tgts,_transform,targetcol,tgtSize);
          DisableLand();
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Pick target", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        void PickTarget()
        {
          pickingPositionTarget = true;
          string message = "Click to select a target";
          ScreenMessages.PostScreenMessage(message, 3.0f, ScreenMessageStyle.UPPER_CENTER);
        }

        [KSPEvent(guiActive = true, guiActiveEditor = true, guiName = "Enable guidance", active = true, guiActiveUnfocused = true, unfocusedRange = 1000)]
        void ToggleGuidance()
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
            autoMode = AutoMode.Off;
            // Re-activate shutdown engines to make it easier for the player
            foreach( ModuleEngines engine in shutdownEngines )
              engine.Activate();
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
          // Update controller vessel parameters
          ComputeMinMaxThrust(out _minThrust,out _maxThrust);
          // New thrust after shutdown should span a range that enables
          // a hover
          controller.amin = (float)(_minThrust/vessel.totalMass);
          controller.amax = (float)(_maxThrust/vessel.totalMass);
          controller.maxThrustAngle = maxThrustAngle;
          Vector3d dr, dv, da, tthrustV;
          float att_err;
          double throttle;
          double tRel = Time.time - _startTime;
          Vector3d att = new Vector3d(vessel.transform.up.x,vessel.transform.up.y,vessel.transform.up.z);
          Vector3d tatt = _transform.InverseTransformVector(att);
          bool shutdownEnginesNow;
          controller.GetControlOutputs(_traj, tr, tv, tatt, new Vector3d(0,-FlightGlobals.getGeeForceAtPosition(r).magnitude,0), (float)tRel, out dr, out dv, out da, out throttle, out tthrustV, out att_err, out shutdownEnginesNow);
          if ((controller.IsFinalDescent()) && (!finalDescentLogged))
          {
            Log("Final descent activated at "+(int)(tr.y+lowestY)+"m",true);
            finalDescentLogged = true;
          }

          Vector3d thrustV = _transform.TransformVector(tthrustV); // transform back to world co-ordinates
          vessel.Autopilot.SAS.lockedMode = false;
          vessel.Autopilot.SAS.SetTargetOrientation(thrustV ,false);
          state.mainThrottle = (float)throttle;

          if (shutdownEnginesNow)
            shutdownEngines = ShutdownOuterEngines(g*(float)vessel.totalMass, true);

          // Drawing
          DrawAlign(tr,dr,_transform,aligncol);
          DrawSteer(tr, tr+3*vessel.vesselSize.x*Vector3d.Normalize(tthrustV), _transform, thrustcol);

          // Open log files
          if ((_vesselWriter == null) && (_logging))
          {
            _vesselWriter = new System.IO.StreamWriter(_vesselLogFilename);
            _vesselWriter.WriteLine("time x y z vx vy vz ax ay az att_err amin amax");
            _tgtWriter = new System.IO.StreamWriter(_tgtLogFilename);
            _tgtWriter.WriteLine("time x y z vx vy vz ax ay az att_err amin amax");
          }

          // Logging
          if ((_logging) && (tRel >= last_t+log_interval))
          {
            LogData(_tgtWriter, tRel, dr, dv, da, 0, solver.amin, solver.amax);
            double thrustProp = 0;
            if ((_maxThrust > 0) && (thrustV.magnitude > 0))
              thrustProp = (controller.amax*GetCurrentThrust())/(tthrustV.magnitude*_maxThrust);
            LogData(_vesselWriter, tRel, tr, tv, tthrustV*thrustProp, att_err, controller.amin, controller.amax);
            last_t = tRel;
          }
        }

        public override void OnUpdate()
        {
          base.OnUpdate();

          List<Target> tgts = new List<Target>(_tgts);

          // Check for changes to slider. This can trigger one of many actions
          // - If autopilot enabled, recompute trajectory
          // - Reset PID controller
          // - Redraw target
          // - Draw/delete track
          bool recomputeTrajectory = false;
          bool redrawTargets = false;
          bool resetPID = false;
          if (tgtHeight != setTgtHeight)
          {
            redrawTargets = true;
            recomputeTrajectory = true;
            setTgtHeight = tgtHeight;
          }
          if (tgtSize != setTgtSize)
          {
            setTgtSize = tgtSize;
            redrawTargets = true;
          }
          if (corrFactor != setCorrFactor)
          {
            setCorrFactor = corrFactor;
            resetPID = true;
          }
          //if (kd1 != setKd1)
          //{
          //  setKd1 = kd1;
          //  resetPID = true;
          //}
          //if (kd2 != setKd2)
          //{
          //  setKd2 = kd2;
          //  resetPID = true;
          //}
          if ((minDescentAngle != setMinDescentAngle)||(maxThrustAngle != setMaxThrustAngle))
          {
            setMinDescentAngle = minDescentAngle;
            setMaxThrustAngle = maxThrustAngle;
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
            DrawTrack(_traj, _transform);
            setShowTrack = showTrack;
          }
          if (pickingPositionTarget)
          {
            if (Input.GetKeyDown(KeyCode.Escape))
            {
              // Previous position
              redrawTargets = true;
              pickingPositionTarget = false;
            }
            RaycastHit hit;
            if (GuiUtils.GetMouseHit(vessel.mainBody,out hit,part))
            {
              // Picked
              double lat, lon, alt;
              vessel.mainBody.GetLatLonAlt(hit.point, out lat, out lon, out alt);
              Target tgt = new Target((float)lat,(float)lon,(float)alt+0.2f,tgtHeight);
              tgts.Add(tgt); // Add temporarily to end of list
              redrawTargets = true;
              // If clicked stop picking
              if (Input.GetMouseButtonDown(0))
              {
                _tgts = new List<Target>(tgts); // Copy to final list of targets
                _transform = SetUpTransform(_tgts[_tgts.Count-1]);
                recomputeTrajectory = true;
                pickingPositionTarget = false;
              }
            }
          }
          // Activate the required updates
          if (redrawTargets)
          {
            setTgtSize = tgtSize;
            if (tgts.Count > 0)
            {
              // Reset height of last target
              tgts[tgts.Count-1].height = tgtHeight;
              _transform = SetUpTransform(tgts[tgts.Count-1]);
              DrawTargets(tgts,_transform,targetcol,tgtSize);
            }
          }
          if ((recomputeTrajectory)&&((autoMode == AutoMode.LandAtTarget)||(autoMode == AutoMode.Failed)))
            EnableLandAtTarget();
          if (resetPID)
          {
            float accelFactor  = HGUtils.LinearMap(corrFactor, corrFactorP1, corrFactorP2, accelFactorP1, accelFactorP2);
            Log("accelFactor="+accelFactor);
            controller.ResetPIDs(corrFactor,ki1,kd1,accelFactor,ki2,kd2,kp1y,ki1y,kd1y,kp2y,ki2y,kd2y);
          }
          if ((!_logging) && (_vesselWriter != null))
            LogStop();
        }
    }
}
