namespace HopperGuidance
{
    public class ThrustVectorTime
    {
        public Vector3d v;
        public float t;
    }

    public class SolveTarget
    {
        public const int X = 1;
        public const int Y = 2;
        public const int Z = 4;

        public Vector3d r;
        public int raxes; // combination of X, Y, Z
        public Vector3d v;
        public int vaxes; // combination of X, Y, Z
        public float t;
    }
}
