// Simple Vessel class using by VesselSim

namespace HopperGuidance
{
  public class Engine
  {
    Vector3d pos;
    double maxGimbal;
    double minThrust;
    double maxThrust;

    public Engine(Vector3d a_pos, double a_maxGimbal, double a_minThrust, double a_maxThrust)
    {
      pos = a_pos;
      maxGimbal = a_maxGimbal;
      minThrust = a_minThrust;
      maxThrust = a_maxThrust;
    }
  }

  public class Vessel
  {
    public Vector3d pos; // position
    public Vector3d vel; // velocity
    public Vector3d att; // attitude
    public double mass;
    public double height;
    public double diameter;
    public double fuel; // engines consume 1 unit of fuel at max thrust per second
    List<Engines> engines;

    public Vessel(double a_mass, double a_height, double a_diameter)
    {
      mass = a_mass;
      height = a_height;
      diameter = a_diameter;
    }

    public AddEngine(Vector3d pos, double maxGimbal, double minThrust, double maxThrust)
    {
      engines.Add(new Engine(pos, maxGimbal, minThrust, maxThrust));
    }
  }
}

