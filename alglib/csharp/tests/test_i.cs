
#pragma warning disable 219
#pragma warning disable 162
using System;
public class MemoryLeaksTest : System.Runtime.ConstrainedExecution.CriticalFinalizerObject
{
    public int dummy;
    public MemoryLeaksTest()
    {
        dummy = 0;
    }
    ~MemoryLeaksTest()
    {
        alglib.free_disposed_items();
        long cnt = alglib.alloc_counter();
        System.Console.WriteLine("Allocation counter checked... "+(cnt==0 ? "OK" : "FAILED"));
        if( cnt!=0 )
            System.Environment.ExitCode = 1;
    }
}
public class MainTest
{

    public static bool doc_test_bool(bool val, bool test_val)
    { return (val && test_val) || (!val && !test_val); }

    public static bool doc_test_int(int val, int test_val)
    { return val==test_val; }

    public static bool doc_test_real(double val, double test_val, double _threshold)
    {
        double s = _threshold>=0 ? 1.0 : Math.Abs(test_val); 
        double threshold = Math.Abs(_threshold);
        return Math.Abs(val-test_val)/s<=threshold; 
    }

    public static bool doc_test_complex(alglib.complex val, alglib.complex test_val, double _threshold)
    { 
        double s = _threshold>=0 ? 1.0 : alglib.math.abscomplex(test_val);
        double threshold = Math.Abs(_threshold);
        return alglib.math.abscomplex(val-test_val)/s<=threshold;
    }


    public static bool doc_test_bool_vector(bool[] val, bool[] test_val)
    {
        int i;
        if( alglib.ap.len(val)!=alglib.ap.len(test_val) )
            return false;
        for(i=0; i<alglib.ap.len(val); i++)
            if( val[i]!=test_val[i] )
                return false;
        return true;
    }

    public static bool doc_test_bool_matrix(bool[,] val, bool[,] test_val)
    {
        int i, j;
        if( alglib.ap.rows(val)!=alglib.ap.rows(test_val) )
            return false;
        if( alglib.ap.cols(val)!=alglib.ap.cols(test_val) )
            return false;
        for(i=0; i<alglib.ap.rows(val); i++)
            for(j=0; j<alglib.ap.cols(val); j++)
                if( val[i,j]!=test_val[i,j] )
                    return false;
        return true;
    }

    public static bool doc_test_int_vector(int[] val, int[] test_val)
    {
        int i;
        if( alglib.ap.len(val)!=alglib.ap.len(test_val) )
            return false;
        for(i=0; i<alglib.ap.len(val); i++)
            if( val[i]!=test_val[i] )
                return false;
        return true;
    }

    public static bool doc_test_int_matrix(int[,] val, int[,] test_val)
    {
        int i, j;
        if( alglib.ap.rows(val)!=alglib.ap.rows(test_val) )
            return false;
        if( alglib.ap.cols(val)!=alglib.ap.cols(test_val) )
            return false;
        for(i=0; i<alglib.ap.rows(val); i++)
            for(j=0; j<alglib.ap.cols(val); j++)
                if( val[i,j]!=test_val[i,j] )
                    return false;
        return true;
    }

    public static bool doc_test_real_vector(double[] val, double[] test_val, double _threshold)
    {
        int i;
        if( alglib.ap.len(val)!=alglib.ap.len(test_val) )
            return false;
        for(i=0; i<alglib.ap.len(val); i++)
        {
            double s = _threshold>=0 ? 1.0 : Math.Abs(test_val[i]); 
            double threshold = Math.Abs(_threshold);
            if( Math.Abs(val[i]-test_val[i])/s>threshold )
                return false;
        }
        return true;
    }

    public static bool doc_test_real_matrix(double[,] val, double[,] test_val, double _threshold)
    {
        int i, j;
        if( alglib.ap.rows(val)!=alglib.ap.rows(test_val) )
            return false;
        if( alglib.ap.cols(val)!=alglib.ap.cols(test_val) )
            return false;
        for(i=0; i<alglib.ap.rows(val); i++)
            for(j=0; j<alglib.ap.cols(val); j++)
            {
                double s = _threshold>=0 ? 1.0 : Math.Abs(test_val[i,j]);
                double threshold = Math.Abs(_threshold);
                if( Math.Abs(val[i,j]-test_val[i,j])/s>threshold )
                    return false;
            }
        return true;
    }

    public static bool doc_test_complex_vector(alglib.complex[] val, alglib.complex[] test_val, double _threshold)
    {
        int i;
        if( alglib.ap.len(val)!=alglib.ap.len(test_val) )
            return false;
        for(i=0; i<alglib.ap.len(val); i++)
        {
            double s = _threshold>=0 ? 1.0 : alglib.math.abscomplex(test_val[i]);
            double threshold = Math.Abs(_threshold);
            if( alglib.math.abscomplex(val[i]-test_val[i])/s>threshold )
                return false;
        }
        return true;
    }

    public static bool doc_test_complex_matrix(alglib.complex[,] val, alglib.complex[,] test_val, double _threshold)
    {
        int i, j;
        if( alglib.ap.rows(val)!=alglib.ap.rows(test_val) )
            return false;
        if( alglib.ap.cols(val)!=alglib.ap.cols(test_val) )
            return false;
        for(i=0; i<alglib.ap.rows(val); i++)
            for(j=0; j<alglib.ap.cols(val); j++)
            {
                double s = _threshold>=0 ? 1.0 : alglib.math.abscomplex(test_val[i,j]);
                double threshold = Math.Abs(_threshold);
                if( alglib.math.abscomplex(val[i,j]-test_val[i,j])/s>threshold )
                    return false;
            }
        return true;
    }

    public static void spoil_vector_by_adding_element<T>(ref T[] x) where T : new()
    {
        int i;
        T[] y = x;
        x = new T[y.Length+1];
        for(i=0; i<y.Length; i++)
            x[i] = y[i];
        x[y.Length] = new T();
    }

    public static void spoil_vector_by_deleting_element<T>(ref T[] x) where T : new()
    {
        int i;
        T[] y = x;
        x = new T[y.Length-1];
        for(i=0; i<y.Length-1; i++)
            x[i] = y[i];
    }

    public static void spoil_matrix_by_adding_row<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0)+1,y.GetLength(1)];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
        for(j=0; j<y.GetLength(1); j++)
            x[y.GetLength(0),j] = new T();
    }

    public static void spoil_matrix_by_deleting_row<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0)-1,y.GetLength(1)];
        for(i=0; i<y.GetLength(0)-1; i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
    }

    public static void spoil_matrix_by_adding_col<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0), y.GetLength(1)+1];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1); j++)
                x[i,j] = y[i,j];
        for(i=0; i<y.GetLength(0); i++)
            x[i,y.GetLength(1)] = new T();
    }

    public static void spoil_matrix_by_deleting_col<T>(ref T[,] x) where T : new()
    {
        int i, j;
        T[,] y = x;
        x = new T[y.GetLength(0), y.GetLength(1)-1];
        for(i=0; i<y.GetLength(0); i++)
            for(j=0; j<y.GetLength(1)-1; j++)
                x[i,j] = y[i,j];
    }

    public static void spoil_vector_by_value<T>(ref T[] x, T val)
    {
        if( x.Length!=0 )
            x[alglib.math.randominteger(x.Length)] = val;
    }
    public static void spoil_matrix_by_value<T>(ref T[,] x, T val)
    {
        if( x.GetLength(0)!=0 && x.GetLength(1)!=0 )
            x[alglib.math.randominteger(x.GetLength(0)),alglib.math.randominteger(x.GetLength(1))] = val;
    }

    public static void function1_func(double[] x, ref double func, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
    }
    public static void function1_grad(double[] x, ref double func, double[] grad, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        // and its derivatives df/d0 and df/dx1
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 400*System.Math.Pow(x[0]+3,3);
        grad[1] = 4*System.Math.Pow(x[1]-3,3);
    }
    public static void function1_hess(double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        // this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        // its derivatives df/d0 and df/dx1
        // and its Hessian.
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 400*System.Math.Pow(x[0]+3,3);
        grad[1] = 4*System.Math.Pow(x[1]-3,3);
        hess[0,0] = 1200*System.Math.Pow(x[0]+3,2);
        hess[0,1] = 0;
        hess[1,0] = 0;
        hess[1,1] = 12*System.Math.Pow(x[1]-3,2);
    }
    public static void  function1_fvec(double[] x, double[] fi, object obj)
    {
        //
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        //
        fi[0] = 10*System.Math.Pow(x[0]+3,2);
        fi[1] = System.Math.Pow(x[1]-3,2);
    }
    public static void  function1_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        // and Jacobian matrix J = [dfi/dxj]
        fi[0] = 10*System.Math.Pow(x[0]+3,2);
        fi[1] = System.Math.Pow(x[1]-3,2);
        jac[0,0] = 20*(x[0]+3);
        jac[0,1] = 0;
        jac[1,0] = 0;
        jac[1,1] = 2*(x[1]-3);
    }
    public static void function2_func(double[] x, ref double func, object obj)
    {
        //
        // this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        //
        func = System.Math.Pow(x[0]*x[0]+1,2) + System.Math.Pow(x[1]-1,2);
    }
    public static void function2_grad(double[] x, ref double func, double[] grad, object obj)
    {
        //
        // this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        // and its derivatives df/d0 and df/dx1
        //
        func = System.Math.Pow(x[0]*x[0]+1,2) + System.Math.Pow(x[1]-1,2);
        grad[0] = 4*(x[0]*x[0]+1)*x[0];
        grad[1] = 2*(x[1]-1);
    }
    public static void function2_hess(double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        //
        // this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        // its gradient and Hessian
        //
        func = System.Math.Pow(x[0]*x[0]+1,2) + System.Math.Pow(x[1]-1,2);
        grad[0] = 4*(x[0]*x[0]+1)*x[0];
        grad[1] = 2*(x[1]-1);
        hess[0,0] = 12*x[0]*x[0]+4;
        hess[0,1] = 0;
        hess[1,0] = 0;
        hess[1,1] = 2;
    }
    public static void  function2_fvec(double[] x, double[] fi, object obj)
    {
        //
        // this callback calculates
        // f0(x0,x1) = 100*(x0+3)^4,
        // f1(x0,x1) = (x1-3)^4
        //
        fi[0] = x[0]*x[0]+1;
        fi[1] = x[1]-1;
    }
    public static void  function2_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates
        // f0(x0,x1) = x0^2+1
        // f1(x0,x1) = x1-1
        // and Jacobian matrix J = [dfi/dxj]
        //
        fi[0] = x[0]*x[0]+1;
        fi[1] = x[1]-1;
        jac[0,0] = 2*x[0];
        jac[0,1] = 0;
        jac[1,0] = 0;
        jac[1,1] = 1;
    }
    public static void  nlcfunc1_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates
        //
        //     f0(x0,x1) = -x0+x1
        //     f1(x0,x1) = x0^2+x1^2-1
        //
        // and Jacobian matrix J = [dfi/dxj]
        //
        fi[0] = -x[0]+x[1];
        fi[1] = x[0]*x[0] + x[1]*x[1] - 1.0;
        jac[0,0] = -1.0;
        jac[0,1] = +1.0;
        jac[1,0] = 2*x[0];
        jac[1,1] = 2*x[1];
    }
    public static void  nlcfunc2_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates
        //
        //     f0(x0,x1,x2) = x0+x1
        //     f1(x0,x1,x2) = x2-exp(x0)
        //     f2(x0,x1,x2) = x0^2+x1^2-1
        //
        // and Jacobian matrix J = [dfi/dxj]
        //
        fi[0] = x[0]+x[1];
        fi[1] = x[2]-System.Math.Exp(x[0]);
        fi[2] = x[0]*x[0] + x[1]*x[1] - 1.0;
        jac[0,0] = 1.0;
        jac[0,1] = 1.0;
        jac[0,2] = 0.0;
        jac[1,0] = -System.Math.Exp(x[0]);
        jac[1,1] = 0.0;
        jac[1,2] = 1.0;
        jac[2,0] = 2*x[0];
        jac[2,1] = 2*x[1];
        jac[2,2] = 0.0;
    }
    public static void  nsfunc1_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates
        //
        //     f0(x0,x1) = 2*|x0|+|x1|
        //
        // and Jacobian matrix J = [ df0/dx0, df0/dx1 ]
        //
        fi[0] = 2*System.Math.Abs(x[0])+System.Math.Abs(x[1]);
        jac[0,0] = 2*System.Math.Sign(x[0]);
        jac[0,1] = System.Math.Sign(x[1]);
    }
    public static void  nsfunc1_fvec(double[] x, double[] fi, object obj)
    {
        //
        // this callback calculates
        //
        //     f0(x0,x1) = 2*|x0|+|x1|
        //
        fi[0] = 2*System.Math.Abs(x[0])+System.Math.Abs(x[1]);
    }
    public static void  nsfunc2_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates function vector
        //
        //     f0(x0,x1) = 2*|x0|+x1
        //     f1(x0,x1) = x0-1
        //     f2(x0,x1) = -x1-1
        //
        // and Jacobian matrix J
        //
        //         [ df0/dx0   df0/dx1 ]
        //     J = [ df1/dx0   df1/dx1 ]
        //         [ df2/dx0   df2/dx1 ]
        //
        fi[0] = 2*System.Math.Abs(x[0])+System.Math.Abs(x[1]);
        jac[0,0] = 2*System.Math.Sign(x[0]);
        jac[0,1] = System.Math.Sign(x[1]);
        fi[1] = x[0]-1;
        jac[1,0] = 1;
        jac[1,1] = 0;
        fi[2] = -x[1]-1;
        jac[2,0] = 0;
        jac[2,1] = -1;
    }
    public static void bad_func(double[] x, ref double func, object obj)
    {
        //
        // this callback calculates 'bad' function,
        // i.e. function with incorrectly calculated derivatives
        //
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
    }
    public static void bad_grad(double[] x, ref double func, double[] grad, object obj)
    {
        //
        // this callback calculates 'bad' function,
        // i.e. function with incorrectly calculated derivatives
        //
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 40*System.Math.Pow(x[0]+3,3);
        grad[1] = 40*System.Math.Pow(x[1]-3,3);
    }
    public static void bad_hess(double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        //
        // this callback calculates 'bad' function,
        // i.e. function with incorrectly calculated derivatives
        //
        func = 100*System.Math.Pow(x[0]+3,4) + System.Math.Pow(x[1]-3,4);
        grad[0] = 40*System.Math.Pow(x[0]+3,3);
        grad[1] = 40*System.Math.Pow(x[1]-3,3);
        hess[0,0] = 120*System.Math.Pow(x[0]+3,2);
        hess[0,1] = 1;
        hess[1,0] = 1;
        hess[1,1] = 120*System.Math.Pow(x[1]-3,2);
    }
    public static void  bad_fvec(double[] x, double[] fi, object obj)
    {
        //
        // this callback calculates 'bad' function,
        // i.e. function with incorrectly calculated derivatives
        //
        fi[0] = 10*System.Math.Pow(x[0]+3,2);
        fi[1] = System.Math.Pow(x[1]-3,2);
    }
    public static void  bad_jac(double[] x, double[] fi, double[,] jac, object obj)
    {
        //
        // this callback calculates 'bad' function,
        // i.e. function with incorrectly calculated derivatives
        //
        fi[0] = 10*System.Math.Pow(x[0]+3,2);
        fi[1] = System.Math.Pow(x[1]-3,2);
        jac[0,0] = 20*(x[0]+3);
        jac[0,1] = 0;
        jac[1,0] = 1;
        jac[1,1] = 20*(x[1]-3);
    }
    public static void function_cx_1_func(double[] c, double[] x, ref double func, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0))
        // where x is a position on X-axis and c is adjustable parameter
        func = System.Math.Exp(-c[0]*x[0]*x[0]);
    }
    public static void function_cx_1_grad(double[] c, double[] x, ref double func, double[] grad, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient G={df/dc[i]}
        // where x is a position on X-axis and c is adjustable parameter.
        // IMPORTANT: gradient is calculated with respect to C, not to X
        func = System.Math.Exp(-c[0]*System.Math.Pow(x[0],2));
        grad[0] = -System.Math.Pow(x[0],2)*func;
    }
    public static void function_cx_1_hess(double[] c, double[] x, ref double func, double[] grad, double[,] hess, object obj)
    {
        // this callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])}
        // where x is a position on X-axis and c is adjustable parameter.
        // IMPORTANT: gradient/Hessian are calculated with respect to C, not to X
        func = System.Math.Exp(-c[0]*System.Math.Pow(x[0],2));
        grad[0] = -System.Math.Pow(x[0],2)*func;
        hess[0,0] = System.Math.Pow(x[0],4)*func;
    }
    public static void ode_function_1_diff(double[] y, double x, double[] dy, object obj)
    {
        // this callback calculates f(y[],x)=-y[0]
        dy[0] = -y[0];
    }
    public static void int_function_1_func(double x, double xminusa, double bminusx, ref double y, object obj)
    {
        // this callback calculates f(x)=exp(x)
        y = Math.Exp(x);
    }
    public static void function_debt_func(double[] c, double[] x, ref double func, object obj)
    {
        //
        // this callback calculates f(c,x)=c[0]*(1+c[1]*(pow(x[0]-1999,c[2])-1))
        //
        func = c[0]*(1+c[1]*(System.Math.Pow(x[0]-1999,c[2])-1));
    }
    public static void s1_grad(double[] x, ref double func, double[] grad, object obj)
    {
        //
        // this callback calculates f(x) = (1+x)^(-0.2) + (1-x)^(-0.3) + 1000*x and its gradient.
        //
        // function is trimmed when we calculate it near the singular points or outside of the [-1,+1].
        // Note that we do NOT calculate gradient in this case.
        //
        if( (x[0]<=-0.999999999999) || (x[0]>=+0.999999999999) )
        {
            func = 1.0E+300;
            return;
        }
        func = System.Math.Pow(1+x[0],-0.2) + System.Math.Pow(1-x[0],-0.3) + 1000*x[0];
        grad[0] = -0.2*System.Math.Pow(1+x[0],-1.2) +0.3*System.Math.Pow(1-x[0],-1.3) + 1000;
    }

    public static void Main(string[] args)
    {
        bool _TotalResult = true;
        bool _TestResult;
        int _spoil_scenario;
        System.Console.WriteLine("C# interface tests. Please wait...");
        alglib.alloc_counter_activate();
        System.Console.WriteLine("Allocation counter activated...");
        try
        {
            //
            // TEST nneighbor_d_1
            //      Nearest neighbor search, KNN queries
            //
            System.Console.WriteLine("0/151");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{0,0},{0,1},{1,0},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    int nx = 2;
                    int ny = 0;
                    int normtype = 2;
                    alglib.kdtree kdt;
                    double[] x;
                    double[,] r = new double[0,0];
                    int k;
                    alglib.kdtreebuild(a, nx, ny, normtype, out kdt);
                    x = new double[]{-1,0};
                    k = alglib.kdtreequeryknn(kdt, x, 1);
                    _TestResult = _TestResult && doc_test_int(k, 1);
                    alglib.kdtreequeryresultsx(kdt, ref r);
                    _TestResult = _TestResult && doc_test_real_matrix(r, new double[,]{{0,0}}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nneighbor_t_2
            //      Subsequent queries; buffered functions must use previously allocated storage (if large enough), so buffer may contain some info from previous call
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{0,0},{0,1},{1,0},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    int nx = 2;
                    int ny = 0;
                    int normtype = 2;
                    alglib.kdtree kdt;
                    double[] x;
                    double[,] rx = new double[0,0];
                    int k;
                    alglib.kdtreebuild(a, nx, ny, normtype, out kdt);
                    x = new double[]{+2,0};
                    k = alglib.kdtreequeryknn(kdt, x, 2, true);
                    _TestResult = _TestResult && doc_test_int(k, 2);
                    alglib.kdtreequeryresultsx(kdt, ref rx);
                    _TestResult = _TestResult && doc_test_real_matrix(rx, new double[,]{{1,0},{1,1}}, 0.05);
                    x = new double[]{-2,0};
                    k = alglib.kdtreequeryknn(kdt, x, 1, true);
                    _TestResult = _TestResult && doc_test_int(k, 1);
                    alglib.kdtreequeryresultsx(kdt, ref rx);
                    _TestResult = _TestResult && doc_test_real_matrix(rx, new double[,]{{0,0},{1,1}}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nneighbor_d_2
            //      Serialization of KD-trees
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{0,0},{0,1},{1,0},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    int nx = 2;
                    int ny = 0;
                    int normtype = 2;
                    alglib.kdtree kdt0;
                    alglib.kdtree kdt1;
                    string s;
                    double[] x;
                    double[,] r0 = new double[0,0];
                    double[,] r1 = new double[0,0];

                    //
                    // Build tree and serialize it
                    //
                    alglib.kdtreebuild(a, nx, ny, normtype, out kdt0);
                    alglib.kdtreeserialize(kdt0, out s);
                    alglib.kdtreeunserialize(s, out kdt1);

                    //
                    // Compare results from KNN queries
                    //
                    x = new double[]{-1,0};
                    alglib.kdtreequeryknn(kdt0, x, 1);
                    alglib.kdtreequeryresultsx(kdt0, ref r0);
                    alglib.kdtreequeryknn(kdt1, x, 1);
                    alglib.kdtreequeryresultsx(kdt1, ref r1);
                    _TestResult = _TestResult && doc_test_real_matrix(r0, new double[,]{{0,0}}, 0.05);
                    _TestResult = _TestResult && doc_test_real_matrix(r1, new double[,]{{0,0}}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST odesolver_d1
            //      Solving y'=-y with ODE solver
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] x = new double[]{0,1,2,3};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double eps = 0.00001;
                    if( _spoil_scenario==7 )
                        eps = (double)System.Double.NaN;
                    if( _spoil_scenario==8 )
                        eps = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        eps = (double)System.Double.NegativeInfinity;
                    double h = 0;
                    if( _spoil_scenario==10 )
                        h = (double)System.Double.NaN;
                    if( _spoil_scenario==11 )
                        h = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==12 )
                        h = (double)System.Double.NegativeInfinity;
                    alglib.odesolverstate s;
                    int m;
                    double[] xtbl;
                    double[,] ytbl;
                    alglib.odesolverreport rep;
                    alglib.odesolverrkck(y, x, eps, h, out s);
                    alglib.odesolversolve(s, ode_function_1_diff, null);
                    alglib.odesolverresults(s, out m, out xtbl, out ytbl, out rep);
                    _TestResult = _TestResult && doc_test_int(m, 4);
                    _TestResult = _TestResult && doc_test_real_vector(xtbl, new double[]{0,1,2,3}, 0.005);
                    _TestResult = _TestResult && doc_test_real_matrix(ytbl, new double[,]{{1},{0.367},{0.135},{0.050}}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "odesolver_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST sparse_d_1
            //      Basic operations with sparse matrices
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<1; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates creation/initialization of the sparse matrix
                    // and matrix-vector multiplication.
                    //
                    // First, we have to create matrix and initialize it. Matrix is initially created
                    // in the Hash-Table format, which allows convenient initialization. We can modify
                    // Hash-Table matrix with sparseset() and sparseadd() functions.
                    //
                    // NOTE: Unlike CRS format, Hash-Table representation allows you to initialize
                    // elements in the arbitrary order. You may see that we initialize a[0][0] first,
                    // then move to the second row, and then move back to the first row.
                    //
                    alglib.sparsematrix s;
                    alglib.sparsecreate(2, 2, out s);
                    alglib.sparseset(s, 0, 0, 2.0);
                    alglib.sparseset(s, 1, 1, 1.0);
                    alglib.sparseset(s, 0, 1, 1.0);

                    alglib.sparseadd(s, 1, 1, 4.0);

                    //
                    // Now S is equal to
                    //   [ 2 1 ]
                    //   [   5 ]
                    // Lets check it by reading matrix contents with sparseget().
                    // You may see that with sparseget() you may read both non-zero
                    // and zero elements.
                    //
                    double v;
                    v = alglib.sparseget(s, 0, 0);
                    _TestResult = _TestResult && doc_test_real(v, 2.0000, 0.005);
                    v = alglib.sparseget(s, 0, 1);
                    _TestResult = _TestResult && doc_test_real(v, 1.0000, 0.005);
                    v = alglib.sparseget(s, 1, 0);
                    _TestResult = _TestResult && doc_test_real(v, 0.0000, 0.005);
                    v = alglib.sparseget(s, 1, 1);
                    _TestResult = _TestResult && doc_test_real(v, 5.0000, 0.005);

                    //
                    // After successful creation we can use our matrix for linear operations.
                    //
                    // However, there is one more thing we MUST do before using S in linear
                    // operations: we have to convert it from HashTable representation (used for
                    // initialization and dynamic operations) to CRS format with sparseconverttocrs()
                    // call. If you omit this call, ALGLIB will generate exception on the first
                    // attempt to use S in linear operations. 
                    //
                    alglib.sparseconverttocrs(s);

                    //
                    // Now S is in the CRS format and we are ready to do linear operations.
                    // Lets calculate A*x for some x.
                    //
                    double[] x = new double[]{1,-1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[0];
                    alglib.sparsemv(s, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{1.000,-5.000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "sparse_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST sparse_d_crs
            //      Advanced topic: creation in the CRS format.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<2; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates creation/initialization of the sparse matrix in the
                    // CRS format.
                    //
                    // Hash-Table format used by default is very convenient (it allows easy
                    // insertion of elements, automatic memory reallocation), but has
                    // significant memory and performance overhead. Insertion of one element 
                    // costs hundreds of CPU cycles, and memory consumption is several times
                    // higher than that of CRS.
                    //
                    // When you work with really large matrices and when you can tell in 
                    // advance how many elements EXACTLY you need, it can be beneficial to 
                    // create matrix in the CRS format from the very beginning.
                    //
                    // If you want to create matrix in the CRS format, you should:
                    // * use sparsecreatecrs() function
                    // * know row sizes in advance (number of non-zero entries in the each row)
                    // * initialize matrix with sparseset() - another function, sparseadd(), is not allowed
                    // * initialize elements from left to right, from top to bottom, each
                    //   element is initialized only once.
                    //
                    alglib.sparsematrix s;
                    int[] row_sizes = new int[]{2,2,2,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_deleting_element(ref row_sizes);
                    alglib.sparsecreatecrs(4, 4, row_sizes, out s);
                    alglib.sparseset(s, 0, 0, 2.0);
                    alglib.sparseset(s, 0, 1, 1.0);
                    alglib.sparseset(s, 1, 1, 4.0);
                    alglib.sparseset(s, 1, 2, 2.0);
                    alglib.sparseset(s, 2, 2, 3.0);
                    alglib.sparseset(s, 2, 3, 1.0);
                    alglib.sparseset(s, 3, 3, 9.0);

                    //
                    // Now S is equal to
                    //   [ 2 1     ]
                    //   [   4 2   ]
                    //   [     3 1 ]
                    //   [       9 ]
                    //
                    // We should point that we have initialized S elements from left to right,
                    // from top to bottom. CRS representation does NOT allow you to do so in
                    // the different order. Try to change order of the sparseset() calls above,
                    // and you will see that your program generates exception.
                    //
                    // We can check it by reading matrix contents with sparseget().
                    // However, you should remember that sparseget() is inefficient on
                    // CRS matrices (it may have to pass through all elements of the row 
                    // until it finds element you need).
                    //
                    double v;
                    v = alglib.sparseget(s, 0, 0);
                    _TestResult = _TestResult && doc_test_real(v, 2.0000, 0.005);
                    v = alglib.sparseget(s, 2, 3);
                    _TestResult = _TestResult && doc_test_real(v, 1.0000, 0.005);

                    // you may see that you can read zero elements (which are not stored) with sparseget()
                    v = alglib.sparseget(s, 3, 2);
                    _TestResult = _TestResult && doc_test_real(v, 0.0000, 0.005);

                    //
                    // After successful creation we can use our matrix for linear operations.
                    // Lets calculate A*x for some x.
                    //
                    double[] x = new double[]{1,-1,1,-1};
                    if( _spoil_scenario==1 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[0];
                    alglib.sparsemv(s, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{1.000,-2.000,2.000,-9}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "sparse_d_crs");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ablas_d_gemm
            //      Matrix multiplication (single-threaded)
            //
            _TestResult = true;
            try
            {
                double[,] a = new double[,]{{2,1},{1,3}};
                double[,] b = new double[,]{{2,1},{0,1}};
                double[,] c = new double[,]{{0,0},{0,0}};

                //
                // rmatrixgemm() function allows us to calculate matrix product C:=A*B or
                // to perform more general operation, C:=alpha*op1(A)*op2(B)+beta*C,
                // where A, B, C are rectangular matrices, op(X) can be X or X^T,
                // alpha and beta are scalars.
                //
                // This function:
                // * can apply transposition and/or multiplication by scalar to operands
                // * can use arbitrary part of matrices A/B (given by submatrix offset)
                // * can store result into arbitrary part of C
                // * for performance reasons requires C to be preallocated
                //
                // Parameters of this function are:
                // * M, N, K            -   sizes of op1(A) (which is MxK), op2(B) (which
                //                          is KxN) and C (which is MxN)
                // * Alpha              -   coefficient before A*B
                // * A, IA, JA          -   matrix A and offset of the submatrix
                // * OpTypeA            -   transformation type:
                //                          0 - no transformation
                //                          1 - transposition
                // * B, IB, JB          -   matrix B and offset of the submatrix
                // * OpTypeB            -   transformation type:
                //                          0 - no transformation
                //                          1 - transposition
                // * Beta               -   coefficient before C
                // * C, IC, JC          -   preallocated matrix C and offset of the submatrix
                //
                // Below we perform simple product C:=A*B (alpha=1, beta=0)
                //
                // IMPORTANT: this function works with preallocated C, which must be large
                //            enough to store multiplication result.
                //
                int m = 2;
                int n = 2;
                int k = 2;
                double alpha = 1.0;
                int ia = 0;
                int ja = 0;
                int optypea = 0;
                int ib = 0;
                int jb = 0;
                int optypeb = 0;
                double beta = 0.0;
                int ic = 0;
                int jc = 0;
                alglib.rmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{4,3},{2,4}}, 0.0001);

                //
                // Now we try to apply some simple transformation to operands: C:=A*B^T
                //
                optypeb = 1;
                alglib.rmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{5,1},{5,3}}, 0.0001);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ablas_d_gemm");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ablas_d_syrk
            //      Symmetric rank-K update (single-threaded)
            //
            _TestResult = true;
            try
            {
                //
                // rmatrixsyrk() function allows us to calculate symmetric rank-K update
                // C := beta*C + alpha*A'*A, where C is square N*N matrix, A is square K*N
                // matrix, alpha and beta are scalars. It is also possible to update by
                // adding A*A' instead of A'*A.
                //
                // Parameters of this function are:
                // * N, K       -   matrix size
                // * Alpha      -   coefficient before A
                // * A, IA, JA  -   matrix and submatrix offsets
                // * OpTypeA    -   multiplication type:
                //                  * 0 - A*A^T is calculated
                //                  * 2 - A^T*A is calculated
                // * Beta       -   coefficient before C
                // * C, IC, JC  -   preallocated input/output matrix and submatrix offsets
                // * IsUpper    -   whether upper or lower triangle of C is updated;
                //                  this function updates only one half of C, leaving
                //                  other half unchanged (not referenced at all).
                //
                // Below we will show how to calculate simple product C:=A'*A
                //
                // NOTE: beta=0 and we do not use previous value of C, but still it
                //       MUST be preallocated.
                //
                int n = 2;
                int k = 1;
                double alpha = 1.0;
                int ia = 0;
                int ja = 0;
                int optypea = 2;
                double beta = 0.0;
                int ic = 0;
                int jc = 0;
                bool isupper = true;
                double[,] a = new double[,]{{1,2}};

                // preallocate space to store result
                double[,] c = new double[,]{{0,0},{0,0}};

                // calculate product, store result into upper part of c
                alglib.rmatrixsyrk(n, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);

                // output result.
                // IMPORTANT: lower triangle of C was NOT updated!
                _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{1,2},{0,4}}, 0.0001);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ablas_d_syrk");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ablas_t_complex
            //      Basis test for complex matrix functions (correctness and presence of SMP support)
            //
            _TestResult = true;
            try
            {
                alglib.complex[,] a;
                alglib.complex[,] b;
                alglib.complex[,] c;

                // test cmatrixgemm()
                a = new alglib.complex[,]{{new alglib.complex(0,2),new alglib.complex(0,1)},{1,3}};
                b = new alglib.complex[,]{{2,1},{0,1}};
                c = new alglib.complex[,]{{0,0},{0,0}};
                alglib.cmatrixgemm(2, 2, 2, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, ref c, 0, 0);
                _TestResult = _TestResult && doc_test_complex_matrix(c, new alglib.complex[,]{{new alglib.complex(0,4),new alglib.complex(0,3)},{2,4}}, 0.0001);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ablas_t_complex");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_r1
            //      Real matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{1,-1},{1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.rmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_matrix(a, new double[,]{{0.5,0.5},{-0.5,0.5}}, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.r1, 0.5, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.rinf, 0.5, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_r1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_c1
            //      Complex matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] a = new alglib.complex[,]{{new alglib.complex(0,1),-1},{new alglib.complex(0,1),1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.cmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_complex_matrix(a, new alglib.complex[,]{{new alglib.complex(0,-0.5),new alglib.complex(0,-0.5)},{-0.5,0.5}}, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.r1, 0.5, 0.00005);
                    _TestResult = _TestResult && doc_test_real(rep.rinf, 0.5, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_c1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_spd1
            //      SPD matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] a = new double[,]{{2,1},{1,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.spdmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_matrix(a, new double[,]{{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_spd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_d_hpd1
            //      HPD matrix inverse
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] a = new alglib.complex[,]{{2,1},{1,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref a);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref a);
                    int info;
                    alglib.matinvreport rep;
                    alglib.hpdmatrixinverse(ref a, out info, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_complex_matrix(a, new alglib.complex[,]{{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_hpd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_t_r1
            //      Real matrix inverse: singular matrix
            //
            _TestResult = true;
            try
            {
                double[,] a = new double[,]{{1,-1},{-2,2}};
                int info;
                alglib.matinvreport rep;
                alglib.rmatrixinverse(ref a, out info, out rep);
                _TestResult = _TestResult && doc_test_int(info, -3);
                _TestResult = _TestResult && doc_test_real(rep.r1, 0.0, 0.00005);
                _TestResult = _TestResult && doc_test_real(rep.rinf, 0.0, 0.00005);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_r1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_t_c1
            //      Complex matrix inverse: singular matrix
            //
            _TestResult = true;
            try
            {
                alglib.complex[,] a = new alglib.complex[,]{{new alglib.complex(0,1),new alglib.complex(0,-1)},{-2,2}};
                int info;
                alglib.matinvreport rep;
                alglib.cmatrixinverse(ref a, out info, out rep);
                _TestResult = _TestResult && doc_test_int(info, -3);
                _TestResult = _TestResult && doc_test_real(rep.r1, 0.0, 0.00005);
                _TestResult = _TestResult && doc_test_real(rep.rinf, 0.0, 0.00005);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_c1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_e_spd1
            //      Attempt to use SPD function on nonsymmetrix matrix
            //
            _TestResult = true;
            try
            {
                double[,] a = new double[,]{{1,0},{1,1}};
                int info;
                alglib.matinvreport rep;
                alglib.spdmatrixinverse(ref a, out info, out rep);
                _TestResult = false;
            }
            catch(alglib.alglibexception)
            {}
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_spd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matinv_e_hpd1
            //      Attempt to use SPD function on nonsymmetrix matrix
            //
            _TestResult = true;
            try
            {
                alglib.complex[,] a = new alglib.complex[,]{{1,0},{1,1}};
                int info;
                alglib.matinvreport rep;
                alglib.hpdmatrixinverse(ref a, out info, out rep);
                _TestResult = false;
            }
            catch(alglib.alglibexception)
            {}
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_hpd1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlbfgs_d_1
            //      Nonlinear optimization by L-BFGS
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // using LBFGS method, with:
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minlbfgssetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsg = 0;
                    if( _spoil_scenario==6 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==9 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlbfgsstate state;
                    alglib.minlbfgscreate(1, x, out state);
                    alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlbfgssetscale(state, s);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target function (C0 continuity violation)
                    // * nonsmoothness of the target function (C1 continuity violation)
                    // * erroneous analytic gradient, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minlbfgsoptguardsmoothness(state);
                    alglib.minlbfgsoptguardgradient(state, 0.001);

                    //
                    // Optimize and examine results.
                    //
                    alglib.minlbfgsreport rep;
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minlbfgsoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlbfgs_d_2
            //      Nonlinear optimization with additional settings and restarts
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // using LBFGS method.
                    //
                    // Several advanced techniques are demonstrated:
                    // * upper limit on step size
                    // * restart from new point
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsg = 0;
                    if( _spoil_scenario==6 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==9 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double stpmax = 0.1;
                    if( _spoil_scenario==15 )
                        stpmax = (double)System.Double.NaN;
                    if( _spoil_scenario==16 )
                        stpmax = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        stpmax = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlbfgsstate state;
                    alglib.minlbfgsreport rep;

                    // create and tune optimizer
                    alglib.minlbfgscreate(1, x, out state);
                    alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlbfgssetstpmax(state, stpmax);
                    alglib.minlbfgssetscale(state, s);

                    // Set up OptGuard integrity checker which catches errors
                    // like nonsmooth targets or errors in the analytic gradient.
                    //
                    // OptGuard is essential at the early prototyping stages.
                    //
                    // NOTE: gradient verification needs 3*N additional function
                    //       evaluations; DO NOT USE IT IN THE PRODUCTION CODE
                    //       because it leads to unnecessary slowdown of your app.
                    alglib.minlbfgsoptguardsmoothness(state);
                    alglib.minlbfgsoptguardgradient(state, 0.001);

                    // first run
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // second run - algorithm is restarted
                    x = new double[]{10,10};
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.minlbfgsrestartfrom(state, x);
                    alglib.minlbfgsoptimize(state, function1_grad, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // check OptGuard integrity report. Why do we need it at all?
                    // Well, try breaking the gradient by adding 1.0 to some
                    // of its components - OptGuard should report it as error.
                    // And it may also catch unintended errors too :)
                    alglib.optguardreport ogrep;
                    alglib.minlbfgsoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlbfgs_numdiff
            //      Nonlinear optimization by L-BFGS with numerical differentiation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // using numerical differentiation to calculate gradient.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double epsg = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==6 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0;
                    if( _spoil_scenario==9 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double diffstep = 1.0e-6;
                    if( _spoil_scenario==12 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        diffstep = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlbfgsstate state;
                    alglib.minlbfgsreport rep;

                    alglib.minlbfgscreatef(1, x, diffstep, out state);
                    alglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits);
                    alglib.minlbfgsoptimize(state, function1_func, null, null);
                    alglib.minlbfgsresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_numdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST linlsqr_d_1
            //      Solution of sparse linear systems with CG
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example illustrates solution of sparse linear least squares problem
                    // with LSQR algorithm.
                    // 
                    // Suppose that we have least squares problem min|A*x-b| with sparse A
                    // represented by sparsematrix object
                    //         [ 1 1 ]
                    //         [ 1 1 ]
                    //     A = [ 2 1 ]
                    //         [ 1   ]
                    //         [   1 ]
                    // and right part b
                    //     [ 4 ]
                    //     [ 2 ]
                    // b = [ 4 ]
                    //     [ 1 ]
                    //     [ 2 ]
                    // and we want to solve this system in the least squares sense using
                    // LSQR algorithm. In order to do so, we have to create left part
                    // (sparsematrix object) and right part (dense array).
                    //
                    // Initially, sparse matrix is created in the Hash-Table format,
                    // which allows easy initialization, but do not allow matrix to be
                    // used in the linear solvers. So after construction you should convert
                    // sparse matrix to CRS format (one suited for linear operations).
                    //
                    alglib.sparsematrix a;
                    alglib.sparsecreate(5, 2, out a);
                    alglib.sparseset(a, 0, 0, 1.0);
                    alglib.sparseset(a, 0, 1, 1.0);
                    alglib.sparseset(a, 1, 0, 1.0);
                    alglib.sparseset(a, 1, 1, 1.0);
                    alglib.sparseset(a, 2, 0, 2.0);
                    alglib.sparseset(a, 2, 1, 1.0);
                    alglib.sparseset(a, 3, 0, 1.0);
                    alglib.sparseset(a, 4, 1, 1.0);

                    //
                    // Now our matrix is fully initialized, but we have to do one more
                    // step - convert it from Hash-Table format to CRS format (see
                    // documentation on sparse matrices for more information about these
                    // formats).
                    //
                    // If you omit this call, ALGLIB will generate exception on the first
                    // attempt to use A in linear operations. 
                    //
                    alglib.sparseconverttocrs(a);

                    //
                    // Initialization of the right part
                    //
                    double[] b = new double[]{4,2,4,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref b);

                    //
                    // Now we have to create linear solver object and to use it for the
                    // solution of the linear system.
                    //
                    alglib.linlsqrstate s;
                    alglib.linlsqrreport rep;
                    double[] x;
                    alglib.linlsqrcreate(5, 2, out s);
                    alglib.linlsqrsolvesparse(s, a, b);
                    alglib.linlsqrresults(s, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{1.000,2.000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "linlsqr_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minbleic_d_1
            //      Nonlinear optimization with bound constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // subject to box constraints
                    //
                    //     -1<=x<=+1, -1<=y<=+1
                    //
                    // using BLEIC optimizer with:
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minbleicsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties:
                    // * set box constraints
                    // * set variable scales
                    // * set stopping criteria
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double epsg = 0;
                    if( _spoil_scenario==11 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==12 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==14 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==15 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==16 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==17 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==18 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minbleicstate state;
                    alglib.minbleiccreate(x, out state);
                    alglib.minbleicsetbc(state, bndl, bndu);
                    alglib.minbleicsetscale(state, s);
                    alglib.minbleicsetcond(state, epsg, epsf, epsx, maxits);

                    //
                    // Then we activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target function (C0 continuity violation)
                    // * nonsmoothness of the target function (C1 continuity violation)
                    // * erroneous analytic gradient, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minbleicoptguardsmoothness(state);
                    alglib.minbleicoptguardgradient(state, 0.001);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.minbleicreport rep;
                    alglib.minbleicoptimize(state, function1_grad, null, null);
                    alglib.minbleicresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,1}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minbleicoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minbleic_d_2
            //      Nonlinear optimization with linear inequality constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // subject to inequality constraints
                    //
                    // * x>=2 (posed as general linear constraint),
                    // * x+y>=6
                    //
                    // using BLEIC optimizer with
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minbleicsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties:
                    // * set linear constraints
                    // * set variable scales
                    // * set stopping criteria
                    //
                    double[] x = new double[]{5,5};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref s);
                    double[,] c = new double[,]{{1,0,2},{1,1,6}};
                    if( _spoil_scenario==7 )
                        spoil_matrix_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_matrix_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==9 )
                        spoil_matrix_by_value(ref c, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==10 )
                        spoil_matrix_by_deleting_row(ref c);
                    if( _spoil_scenario==11 )
                        spoil_matrix_by_deleting_col(ref c);
                    int[] ct = new int[]{1,1};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref ct);
                    alglib.minbleicstate state;
                    double epsg = 0;
                    if( _spoil_scenario==13 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==16 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==17 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==18 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==20 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==21 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;

                    alglib.minbleiccreate(x, out state);
                    alglib.minbleicsetlc(state, c, ct);
                    alglib.minbleicsetscale(state, s);
                    alglib.minbleicsetcond(state, epsg, epsf, epsx, maxits);

                    //
                    // Then we activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target function (C0 continuity violation)
                    // * nonsmoothness of the target function (C1 continuity violation)
                    // * erroneous analytic gradient, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minbleicoptguardsmoothness(state);
                    alglib.minbleicoptguardgradient(state, 0.001);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.minbleicreport rep;
                    alglib.minbleicoptimize(state, function1_grad, null, null);
                    alglib.minbleicresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2,4}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minbleicoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minbleic_numdiff
            //      Nonlinear optimization with bound constraints and numerical differentiation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // subject to box constraints
                    //
                    //     -1<=x<=+1, -1<=y<=+1
                    //
                    // using BLEIC optimizer with:
                    // * numerical differentiation being used
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minbleicsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties:
                    // * set box constraints
                    // * set variable scales
                    // * set stopping criteria
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_deleting_element(ref bndu);
                    alglib.minbleicstate state;
                    double epsg = 0;
                    if( _spoil_scenario==11 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==12 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==14 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==15 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==16 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==17 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==18 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    double diffstep = 1.0e-6;
                    if( _spoil_scenario==20 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==21 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==22 )
                        diffstep = (double)System.Double.NegativeInfinity;

                    alglib.minbleiccreatef(x, diffstep, out state);
                    alglib.minbleicsetbc(state, bndl, bndu);
                    alglib.minbleicsetscale(state, s);
                    alglib.minbleicsetcond(state, epsg, epsf, epsx, maxits);

                    //
                    // Then we activate OptGuard integrity checking.
                    //
                    // Numerical differentiation always produces "correct" gradient
                    // (with some truncation error, but unbiased). Thus, we just have
                    // to check smoothness properties of the target: C0 and C1 continuity.
                    //
                    // Sometimes user accidentally tries to solve nonsmooth problems
                    // with smooth optimizer. OptGuard helps to detect such situations
                    // early, at the prototyping stage.
                    //
                    alglib.minbleicoptguardsmoothness(state);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.minbleicreport rep;
                    alglib.minbleicoptimize(state, function1_func, null, null);
                    alglib.minbleicresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,1}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // Want to challenge OptGuard? Try to make your problem
                    // nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    // re-run optimizer.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minbleicoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_numdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minqp_d_u1
            //      Unconstrained dense quadratic programming
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<17; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    //
                    // Exact solution is [x0,x1] = [3,2]
                    //
                    // We provide algorithm with starting point, although in this case
                    // (dense matrix, no constraints) it can work without such information.
                    //
                    // Several QP solvers are tried: QuickQP, BLEIC, DENSE-AUL.
                    //
                    // IMPORTANT: this solver minimizes  following  function:
                    //     f(x) = 0.5*x'*A*x + b'*x.
                    // Note that quadratic term has 0.5 before it. So if you want to minimize
                    // quadratic function, you should rewrite it in such way that quadratic term
                    // is multiplied by 0.5 too.
                    //
                    // For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    //     f(x) = 0.5*(2*x0^2+2*x1^2) + .... 
                    // and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    //
                    double[,] a = new double[,]{{2,0},{0,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref a);
                    double[] b = new double[]{-6,-4};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref b);
                    double[] x0 = new double[]{0,1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref x0);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] x;
                    alglib.minqpstate state;
                    alglib.minqpreport rep;

                    // create solver, set quadratic/linear terms
                    alglib.minqpcreate(2, out state);
                    alglib.minqpsetquadraticterm(state, a);
                    alglib.minqpsetlinearterm(state, b);
                    alglib.minqpsetstartingpoint(state, x0);

                    // Set scale of the parameters.
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    //       which automatically determines variable scales.
                    alglib.minqpsetscale(state, s);

                    //
                    // Solve problem with QuickQP solver.
                    //
                    // This solver is intended for medium and large-scale problems with box
                    // constraints (general linear constraints are not supported), but it can
                    // also be efficiently used on unconstrained problems.
                    //
                    // Default stopping criteria are used, Newton phase is active.
                    //
                    alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{3,2}, 0.005);

                    //
                    // Solve problem with BLEIC-based QP solver.
                    //
                    // This solver is intended for problems with moderate (up to 50) number
                    // of general linear constraints and unlimited number of box constraints.
                    // Of course, unconstrained problems can be solved too.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{3,2}, 0.005);

                    //
                    // Solve problem with DENSE-AUL solver.
                    //
                    // This solver is optimized for problems with up to several thousands of
                    // variables and large amount of general linear constraints. Problems with
                    // less than 50 general linear constraints can be efficiently solved with
                    // BLEIC, problems with box-only constraints can be solved with QuickQP.
                    // However, DENSE-AUL will work in any (including unconstrained) case.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{3,2}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_u1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minqp_d_bc1
            //      Bound constrained dense quadratic programming
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    // subject to bound constraints 0<=x0<=2.5, 0<=x1<=2.5
                    //
                    // Exact solution is [x0,x1] = [2.5,2]
                    //
                    // We provide algorithm with starting point. With such small problem good starting
                    // point is not really necessary, but with high-dimensional problem it can save us
                    // a lot of time.
                    //
                    // Several QP solvers are tried: QuickQP, BLEIC, DENSE-AUL.
                    //
                    // IMPORTANT: this solver minimizes  following  function:
                    //     f(x) = 0.5*x'*A*x + b'*x.
                    // Note that quadratic term has 0.5 before it. So if you want to minimize
                    // quadratic function, you should rewrite it in such way that quadratic term
                    // is multiplied by 0.5 too.
                    // For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    //     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    // and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    //
                    double[,] a = new double[,]{{2,0},{0,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref a);
                    double[] b = new double[]{-6,-4};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref b);
                    double[] x0 = new double[]{0,1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref x0);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{0.0,0.0};
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{2.5,2.5};
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double[] x;
                    alglib.minqpstate state;
                    alglib.minqpreport rep;

                    // create solver, set quadratic/linear terms
                    alglib.minqpcreate(2, out state);
                    alglib.minqpsetquadraticterm(state, a);
                    alglib.minqpsetlinearterm(state, b);
                    alglib.minqpsetstartingpoint(state, x0);
                    alglib.minqpsetbc(state, bndl, bndu);

                    // Set scale of the parameters.
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    //       which automatically determines variable scales.
                    alglib.minqpsetscale(state, s);

                    //
                    // Solve problem with QuickQP solver.
                    //
                    // This solver is intended for medium and large-scale problems with box
                    // constraints (general linear constraints are not supported).
                    //
                    // Default stopping criteria are used, Newton phase is active.
                    //
                    alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 4);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2.5,2}, 0.005);

                    //
                    // Solve problem with BLEIC-based QP solver.
                    //
                    // This solver is intended for problems with moderate (up to 50) number
                    // of general linear constraints and unlimited number of box constraints.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2.5,2}, 0.005);

                    //
                    // Solve problem with DENSE-AUL solver.
                    //
                    // This solver is optimized for problems with up to several thousands of
                    // variables and large amount of general linear constraints. Problems with
                    // less than 50 general linear constraints can be efficiently solved with
                    // BLEIC, problems with box-only constraints can be solved with QuickQP.
                    // However, DENSE-AUL will work in any (including unconstrained) case.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2.5,2}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_bc1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minqp_d_lc1
            //      Linearly constrained dense quadratic programming
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    // subject to linear constraint x0+x1<=2
                    //
                    // Exact solution is [x0,x1] = [1.5,0.5]
                    //
                    // IMPORTANT: this solver minimizes  following  function:
                    //     f(x) = 0.5*x'*A*x + b'*x.
                    // Note that quadratic term has 0.5 before it. So if you want to minimize
                    // quadratic function, you should rewrite it in such way that quadratic term
                    // is multiplied by 0.5 too.
                    // For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    //     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    // and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    //
                    double[,] a = new double[,]{{2,0},{0,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref a);
                    double[] b = new double[]{-6,-4};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref b);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref s);
                    double[,] c = new double[,]{{1.0,1.0,2.0}};
                    if( _spoil_scenario==13 )
                        spoil_matrix_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_matrix_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==15 )
                        spoil_matrix_by_value(ref c, (double)System.Double.NegativeInfinity);
                    int[] ct = new int[]{-1};
                    double[] x;
                    alglib.minqpstate state;
                    alglib.minqpreport rep;

                    // create solver, set quadratic/linear terms
                    alglib.minqpcreate(2, out state);
                    alglib.minqpsetquadraticterm(state, a);
                    alglib.minqpsetlinearterm(state, b);
                    alglib.minqpsetlc(state, c, ct);

                    // Set scale of the parameters.
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    //       which automatically determines variable scales.
                    alglib.minqpsetscale(state, s);

                    //
                    // Solve problem with BLEIC-based QP solver.
                    //
                    // This solver is intended for problems with moderate (up to 50) number
                    // of general linear constraints and unlimited number of box constraints.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{1.500,0.500}, 0.05);

                    //
                    // Solve problem with DENSE-AUL solver.
                    //
                    // This solver is optimized for problems with up to several thousands of
                    // variables and large amount of general linear constraints. Problems with
                    // less than 50 general linear constraints can be efficiently solved with
                    // BLEIC, problems with box-only constraints can be solved with QuickQP.
                    // However, DENSE-AUL will work in any (including unconstrained) case.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{1.500,0.500}, 0.05);

                    //
                    // Solve problem with QuickQP solver.
                    //
                    // This solver is intended for medium and large-scale problems with box
                    // constraints, and...
                    //
                    // ...Oops! It does not support general linear constraints, -5 returned as completion code!
                    //
                    alglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, -5);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_lc1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minqp_d_u2
            //      Unconstrained sparse quadratic programming
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1,
                    // with quadratic term given by sparse matrix structure.
                    //
                    // Exact solution is [x0,x1] = [3,2]
                    //
                    // We provide algorithm with starting point, although in this case
                    // (dense matrix, no constraints) it can work without such information.
                    //
                    // IMPORTANT: this solver minimizes  following  function:
                    //     f(x) = 0.5*x'*A*x + b'*x.
                    // Note that quadratic term has 0.5 before it. So if you want to minimize
                    // quadratic function, you should rewrite it in such way that quadratic term
                    // is multiplied by 0.5 too.
                    //
                    // For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    //     f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    // and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    //
                    alglib.sparsematrix a;
                    double[] b = new double[]{-6,-4};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref b);
                    double[] x0 = new double[]{0,1};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref x0);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] x;
                    alglib.minqpstate state;
                    alglib.minqpreport rep;

                    // initialize sparsematrix structure
                    alglib.sparsecreate(2, 2, 0, out a);
                    alglib.sparseset(a, 0, 0, 2.0);
                    alglib.sparseset(a, 1, 1, 2.0);

                    // create solver, set quadratic/linear terms
                    alglib.minqpcreate(2, out state);
                    alglib.minqpsetquadratictermsparse(state, a, true);
                    alglib.minqpsetlinearterm(state, b);
                    alglib.minqpsetstartingpoint(state, x0);

                    // Set scale of the parameters.
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    // NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    //       which automatically determines variable scales.
                    alglib.minqpsetscale(state, s);

                    //
                    // Solve problem with BLEIC-based QP solver.
                    //
                    // This solver is intended for problems with moderate (up to 50) number
                    // of general linear constraints and unlimited number of box constraints.
                    // It also supports sparse problems.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{3,2}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_u2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minqp_d_nonconvex
            //      Nonconvex quadratic programming
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of nonconvex function
                    //     F(x0,x1) = -(x0^2+x1^2)
                    // subject to constraints x0,x1 in [1.0,2.0]
                    // Exact solution is [x0,x1] = [2,2].
                    //
                    // Non-convex problems are harded to solve than convex ones, and they
                    // may have more than one local minimum. However, ALGLIB solves may deal
                    // with such problems (altough they do not guarantee convergence to
                    // global minimum).
                    //
                    // IMPORTANT: this solver minimizes  following  function:
                    //     f(x) = 0.5*x'*A*x + b'*x.
                    // Note that quadratic term has 0.5 before it. So if you want to minimize
                    // quadratic function, you should rewrite it in such way that quadratic term
                    // is multiplied by 0.5 too.
                    //
                    // For example, our function is f(x)=-(x0^2+x1^2), but we rewrite it as 
                    //     f(x) = 0.5*(-2*x0^2-2*x1^2)
                    // and pass diag(-2,-2) as quadratic term - NOT diag(-1,-1)!
                    //
                    double[,] a = new double[,]{{-2,0},{0,-2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref a);
                    double[] x0 = new double[]{1,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref x0);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{1.0,1.0};
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{2.0,2.0};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double[] x;
                    alglib.minqpstate state;
                    alglib.minqpreport rep;

                    // create solver, set quadratic/linear terms, constraints
                    alglib.minqpcreate(2, out state);
                    alglib.minqpsetquadraticterm(state, a);
                    alglib.minqpsetstartingpoint(state, x0);
                    alglib.minqpsetbc(state, bndl, bndu);

                    // Set scale of the parameters.
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    // NOTE: there also exists minqpsetscaleautodiag() function
                    //       which automatically determines variable scales; however,
                    //       it does NOT work for non-convex problems.
                    alglib.minqpsetscale(state, s);

                    //
                    // Solve problem with BLEIC-based QP solver.
                    //
                    // This solver is intended for problems with moderate (up to 50) number
                    // of general linear constraints and unlimited number of box constraints.
                    //
                    // It may solve non-convex problems as long as they are bounded from
                    // below under constraints.
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2,2}, 0.005);

                    //
                    // Solve problem with DENSE-AUL solver.
                    //
                    // This solver is optimized for problems with up to several thousands of
                    // variables and large amount of general linear constraints. Problems with
                    // less than 50 general linear constraints can be efficiently solved with
                    // BLEIC, problems with box-only constraints can be solved with QuickQP.
                    // However, DENSE-AUL will work in any (including unconstrained) case.
                    //
                    // Algorithm convergence is guaranteed only for convex case, but you may
                    // expect that it will work for non-convex problems too (because near the
                    // solution they are locally convex).
                    //
                    // Default stopping criteria are used.
                    //
                    alglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{2,2}, 0.005);

                    // Hmm... this problem is bounded from below (has solution) only under constraints.
                    // What it we remove them?
                    //
                    // You may see that BLEIC algorithm detects unboundedness of the problem, 
                    // -4 is returned as completion code. However, DENSE-AUL is unable to detect
                    // such situation and it will cycle forever (we do not test it here).
                    double[] nobndl = new double[]{-System.Double.PositiveInfinity,-System.Double.PositiveInfinity};
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref nobndl, (double)System.Double.NaN);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_deleting_element(ref nobndl);
                    double[] nobndu = new double[]{System.Double.PositiveInfinity,+System.Double.PositiveInfinity};
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref nobndu, (double)System.Double.NaN);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref nobndu);
                    alglib.minqpsetbc(state, nobndl, nobndu);
                    alglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
                    alglib.minqpoptimize(state);
                    alglib.minqpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, -4);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_nonconvex");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlp_basic
            //      Basic linear programming example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates how to minimize
                    //
                    //     F(x0,x1) = -0.1*x0 - x1
                    //
                    // subject to box constraints
                    //
                    //     -1 <= x0,x1 <= +1 
                    //
                    // and general linear constraints
                    //
                    //     x0 - x1 >= -1
                    //     x0 + x1 <=  1
                    //
                    // We use dual simplex solver provided by ALGLIB for this task. Box
                    // constraints are specified by means of constraint vectors bndl and
                    // bndu (we have bndl<=x<=bndu). General linear constraints are
                    // specified as AL<=A*x<=AU, with AL/AU being 2x1 vectors and A being
                    // 2x2 matrix.
                    //
                    // NOTE: some/all components of AL/AU can be +-INF, same applies to
                    //       bndl/bndu. You can also have AL[I]=AU[i] (as well as
                    //       BndL[i]=BndU[i]).
                    //
                    double[,] a = new double[,]{{1,-1},{1,+1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_deleting_row(ref a);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_deleting_col(ref a);
                    double[] al = new double[]{-1,-System.Double.PositiveInfinity};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref al, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref al);
                    double[] au = new double[]{System.Double.PositiveInfinity,+1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref au, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref au);
                    double[] c = new double[]{-0.1,-1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref c);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double[] x;
                    alglib.minlpstate state;
                    alglib.minlpreport rep;

                    alglib.minlpcreate(2, out state);

                    //
                    // Set cost vector, box constraints, general linear constraints.
                    //
                    // Box constraints can be set in one call to minlpsetbc() or minlpsetbcall()
                    // (latter sets same constraints for all variables and accepts two scalars
                    // instead of two vectors).
                    //
                    // General linear constraints can be specified in several ways:
                    // * minlpsetlc2dense() - accepts dense 2D array as input; sometimes this
                    //   approach is more convenient, although less memory-efficient.
                    // * minlpsetlc2() - accepts sparse matrix as input
                    // * minlpaddlc2dense() - appends one row to the current set of constraints;
                    //   row being appended is specified as dense vector
                    // * minlpaddlc2() - appends one row to the current set of constraints;
                    //   row being appended is specified as sparse set of elements
                    // Independently from specific function being used, LP solver uses sparse
                    // storage format for internal representation of constraints.
                    //
                    alglib.minlpsetcost(state, c);
                    alglib.minlpsetbc(state, bndl, bndu);
                    alglib.minlpsetlc2dense(state, a, al, au, 2);

                    //
                    // Set scale of the parameters.
                    //
                    // It is strongly recommended that you set scale of your variables.
                    // Knowing their scales is essential for evaluation of stopping criteria
                    // and for preconditioning of the algorithm steps.
                    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    //
                    alglib.minlpsetscale(state, s);

                    // Solve
                    alglib.minlpoptimize(state);
                    alglib.minlpresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{0,1}, 0.0005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlp_basic");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minnlc_d_inequality
            //      Nonlinearly constrained optimization (inequality constraints)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = -x0+x1
                    //
                    // subject to box constraints
                    //
                    //    x0>=0, x1>=0
                    //
                    // and nonlinear inequality constraint
                    //
                    //    x0^2 + x1^2 - 1 <= 0
                    //
                    double[] x0 = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    double[] bndl = new double[]{0,0};
                    double[] bndu = new double[]{System.Double.PositiveInfinity,+System.Double.PositiveInfinity};
                    alglib.minnlcstate state;

                    //
                    // Create optimizer object and tune its settings:
                    // * epsx=0.000001  stopping condition for inner iterations
                    // * s=[1,1]        all variables have unit scale; it is important to
                    //                  tell optimizer about scales of your variables - it
                    //                  greatly accelerates convergence and helps to perform
                    //                  some important integrity checks.
                    //
                    alglib.minnlccreate(2, x0, out state);
                    alglib.minnlcsetcond(state, epsx, maxits);
                    alglib.minnlcsetscale(state, s);

                    //
                    // Choose one of the nonlinear programming solvers supported by minnlc
                    // optimizer:
                    // * SQP - sequential quadratic programming NLP solver
                    // * AUL - augmented Lagrangian NLP solver
                    // * SLP - successive linear programming NLP solver
                    //
                    // Different solvers have different properties:
                    // * SQP needs less function evaluations than any other solver, but it
                    //   has much higher iteration cost than other solvers (a QP subproblem
                    //   has to be solved during each step)
                    // * AUL solver has cheaper iterations, but needs more target function
                    //   evaluations
                    // * SLP is the most robust solver provided by ALGLIB, but it performs
                    //   order of magnitude more iterations than SQP.
                    //
                    // In the code below we set solver to be AUL but then override it with SLP,
                    // and then with SQP, so the effective choice is to use SLP. We recommend
                    // you to use SQP at least for early prototyping stages, and then switch
                    // to AUL if possible.
                    //
                    double rho = 1000.0;
                    int outerits = 5;
                    alglib.minnlcsetalgoaul(state, rho, outerits);
                    alglib.minnlcsetalgoslp(state);
                    alglib.minnlcsetalgosqp(state);

                    //
                    // Set constraints:
                    //
                    // 1. boundary constraints are passed with minnlcsetbc() call
                    //
                    // 2. nonlinear constraints are more tricky - you can not "pack" general
                    //    nonlinear function into double precision array. That's why
                    //    minnlcsetnlc() does not accept constraints itself - only constraint
                    //    counts are passed: first parameter is number of equality constraints,
                    //    second one is number of inequality constraints.
                    //
                    //    As for constraining functions - these functions are passed as part
                    //    of problem Jacobian (see below).
                    //
                    // NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    //       linear and general nonlinear constraints. This example does not
                    //       show how to work with general linear constraints, but you can
                    //       easily find it in documentation on minnlcsetlc() function.
                    //
                    alglib.minnlcsetbc(state, bndl, bndu);
                    alglib.minnlcsetnlc(state, 0, 1);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target/constraints (C0 continuity violation)
                    // * nonsmoothness of the target/constraints (C1 continuity violation)
                    // * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minnlcoptguardsmoothness(state);
                    alglib.minnlcoptguardgradient(state, 0.001);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints.
                    //
                    // So, our vector function has form
                    //
                    //     {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
                    //
                    // with Jacobian
                    //
                    //         [  -1    +1  ]
                    //     J = [            ]
                    //         [ 2*x0  2*x1 ]
                    //
                    // with f0 being target function, f1 being constraining function. Number
                    // of equality/inequality constraints is specified by minnlcsetnlc(),
                    // with equality ones always being first, inequality ones being last.
                    //
                    alglib.minnlcreport rep;
                    double[] x1;
                    alglib.minnlcoptimize(state, nlcfunc1_jac, null, null);
                    alglib.minnlcresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{1.0000,0.0000}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minnlcoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_inequality");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minnlc_d_equality
            //      Nonlinearly constrained optimization (equality constraints)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = -x0+x1
                    //
                    // subject to nonlinear equality constraint
                    //
                    //    x0^2 + x1^2 - 1 = 0
                    //
                    double[] x0 = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnlcstate state;

                    //
                    // Create optimizer object and tune its settings:
                    // * epsx=0.000001  stopping condition for inner iterations
                    // * s=[1,1]        all variables have unit scale
                    //
                    alglib.minnlccreate(2, x0, out state);
                    alglib.minnlcsetcond(state, epsx, maxits);
                    alglib.minnlcsetscale(state, s);

                    //
                    // Choose one of the nonlinear programming solvers supported by minnlc
                    // optimizer:
                    // * SLP - successive linear programming NLP solver
                    // * AUL - augmented Lagrangian NLP solver
                    //
                    // Different solvers have different properties:
                    // * SLP is the most robust solver provided by ALGLIB: it can solve both
                    //   convex and nonconvex optimization problems, it respects box and
                    //   linear constraints (after you find feasible point it won't move away
                    //   from the feasible area) and tries to respect nonlinear constraints
                    //   as much as possible. It also usually needs less function evaluations
                    //   to converge than AUL.
                    //   However, it solves LP subproblems at each iterations which adds
                    //   significant overhead to its running time. Sometimes it can be as much
                    //   as 7x times slower than AUL.
                    // * AUL solver is less robust than SLP - it can violate box and linear
                    //   constraints at any moment, and it is intended for convex optimization
                    //   problems (although in many cases it can deal with nonconvex ones too).
                    //   Also, unlike SLP it needs some tuning (penalty factor and number of
                    //   outer iterations).
                    //   However, it is often much faster than the current version of SLP.
                    //
                    // In the code below we set solver to be AUL but then override it with SLP,
                    // so the effective choice is to use SLP. We recommend you to use SLP at
                    // least for early prototyping stages.
                    //
                    // You can comment out line with SLP if you want to solve your problem with
                    // AUL solver.
                    //
                    double rho = 1000.0;
                    int outerits = 5;
                    alglib.minnlcsetalgoaul(state, rho, outerits);
                    alglib.minnlcsetalgoslp(state);

                    //
                    // Set constraints:
                    //
                    // Nonlinear constraints are tricky - you can not "pack" general
                    // nonlinear function into double precision array. That's why
                    // minnlcsetnlc() does not accept constraints itself - only constraint
                    // counts are passed: first parameter is number of equality constraints,
                    // second one is number of inequality constraints.
                    //
                    // As for constraining functions - these functions are passed as part
                    // of problem Jacobian (see below).
                    //
                    // NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    //       linear and general nonlinear constraints. This example does not
                    //       show how to work with general linear constraints, but you can
                    //       easily find it in documentation on minnlcsetbc() and
                    //       minnlcsetlc() functions.
                    //
                    alglib.minnlcsetnlc(state, 1, 0);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target/constraints (C0 continuity violation)
                    // * nonsmoothness of the target/constraints (C1 continuity violation)
                    // * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minnlcoptguardsmoothness(state);
                    alglib.minnlcoptguardgradient(state, 0.001);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints.
                    //
                    // So, our vector function has form
                    //
                    //     {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
                    //
                    // with Jacobian
                    //
                    //         [  -1    +1  ]
                    //     J = [            ]
                    //         [ 2*x0  2*x1 ]
                    //
                    // with f0 being target function, f1 being constraining function. Number
                    // of equality/inequality constraints is specified by minnlcsetnlc(),
                    // with equality ones always being first, inequality ones being last.
                    //
                    alglib.minnlcreport rep;
                    double[] x1;
                    alglib.minnlcoptimize(state, nlcfunc1_jac, null, null);
                    alglib.minnlcresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{0.70710,-0.70710}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minnlcoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_equality");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minnlc_d_mixed
            //      Nonlinearly constrained optimization with mixed equality/inequality constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = x0+x1
                    //
                    // subject to nonlinear inequality constraint
                    //
                    //    x0^2 + x1^2 - 1 <= 0
                    //
                    // and nonlinear equality constraint
                    //
                    //    x2-exp(x0) = 0
                    //
                    double[] x0 = new double[]{0,0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnlcstate state;
                    alglib.minnlcreport rep;
                    double[] x1;

                    //
                    // Create optimizer object and tune its settings:
                    // * epsx=0.000001  stopping condition for inner iterations
                    // * s=[1,1]        all variables have unit scale
                    // * upper limit on step length is specified (to avoid probing locations where exp() is large)
                    //
                    alglib.minnlccreate(3, x0, out state);
                    alglib.minnlcsetcond(state, epsx, maxits);
                    alglib.minnlcsetscale(state, s);
                    alglib.minnlcsetstpmax(state, 10.0);

                    //
                    // Choose one of the nonlinear programming solvers supported by minnlc
                    // optimizer:
                    // * SLP - successive linear programming NLP solver
                    // * AUL - augmented Lagrangian NLP solver
                    //
                    // Different solvers have different properties:
                    // * SLP is the most robust solver provided by ALGLIB: it can solve both
                    //   convex and nonconvex optimization problems, it respects box and
                    //   linear constraints (after you find feasible point it won't move away
                    //   from the feasible area) and tries to respect nonlinear constraints
                    //   as much as possible. It also usually needs less function evaluations
                    //   to converge than AUL.
                    //   However, it solves LP subproblems at each iterations which adds
                    //   significant overhead to its running time. Sometimes it can be as much
                    //   as 7x times slower than AUL.
                    // * AUL solver is less robust than SLP - it can violate box and linear
                    //   constraints at any moment, and it is intended for convex optimization
                    //   problems (although in many cases it can deal with nonconvex ones too).
                    //   Also, unlike SLP it needs some tuning (penalty factor and number of
                    //   outer iterations).
                    //   However, it is often much faster than the current version of SLP.
                    //
                    // In the code below we set solver to be AUL but then override it with SLP,
                    // so the effective choice is to use SLP. We recommend you to use SLP at
                    // least for early prototyping stages.
                    //
                    // You can comment out line with SLP if you want to solve your problem with
                    // AUL solver.
                    //
                    double rho = 1000.0;
                    int outerits = 5;
                    alglib.minnlcsetalgoaul(state, rho, outerits);
                    alglib.minnlcsetalgoslp(state);

                    //
                    // Set constraints:
                    //
                    // Nonlinear constraints are tricky - you can not "pack" general
                    // nonlinear function into double precision array. That's why
                    // minnlcsetnlc() does not accept constraints itself - only constraint
                    // counts are passed: first parameter is number of equality constraints,
                    // second one is number of inequality constraints.
                    //
                    // As for constraining functions - these functions are passed as part
                    // of problem Jacobian (see below).
                    //
                    // NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    //       linear and general nonlinear constraints. This example does not
                    //       show how to work with boundary or general linear constraints, but you
                    //       can easily find it in documentation on minnlcsetbc() and
                    //       minnlcsetlc() functions.
                    //
                    alglib.minnlcsetnlc(state, 1, 1);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target/constraints (C0 continuity violation)
                    // * nonsmoothness of the target/constraints (C1 continuity violation)
                    // * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minnlcoptguardsmoothness(state);
                    alglib.minnlcoptguardgradient(state, 0.001);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints.
                    //
                    // So, our vector function has form
                    //
                    //     {f0,f1,f2} = { x0+x1 , x2-exp(x0) , x0^2+x1^2-1 }
                    //
                    // with Jacobian
                    //
                    //         [  +1      +1       0 ]
                    //     J = [-exp(x0)  0        1 ]
                    //         [ 2*x0    2*x1      0 ]
                    //
                    // with f0 being target function, f1 being equality constraint "f1=0",
                    // f2 being inequality constraint "f2<=0". Number of equality/inequality
                    // constraints is specified by minnlcsetnlc(), with equality ones always
                    // being first, inequality ones being last.
                    //
                    alglib.minnlcoptimize(state, nlcfunc2_jac, null, null);
                    alglib.minnlcresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{-0.70710,-0.70710,0.49306}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minnlcoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_mixed");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minbc_d_1
            //      Nonlinear optimization with box constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // subject to box constraints
                    //
                    //     -1<=x<=+1, -1<=y<=+1
                    //
                    // using MinBC optimizer with:
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minbcsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties:
                    // * set box constraints
                    // * set variable scales
                    // * set stopping criteria
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_deleting_element(ref bndu);
                    alglib.minbcstate state;
                    double epsg = 0;
                    if( _spoil_scenario==11 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==12 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==14 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==15 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==16 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==17 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==18 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minbccreate(x, out state);
                    alglib.minbcsetbc(state, bndl, bndu);
                    alglib.minbcsetscale(state, s);
                    alglib.minbcsetcond(state, epsg, epsf, epsx, maxits);

                    //
                    // Then we activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target function (C0 continuity violation)
                    // * nonsmoothness of the target function (C1 continuity violation)
                    // * erroneous analytic gradient, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.minbcoptguardsmoothness(state);
                    alglib.minbcoptguardgradient(state, 0.001);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.minbcreport rep;
                    alglib.minbcoptimize(state, function1_grad, null, null);
                    alglib.minbcresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,1}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minbcoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minbc_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minbc_numdiff
            //      Nonlinear optimization with bound constraints and numerical differentiation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // subject to box constraints
                    //
                    //    -1<=x<=+1, -1<=y<=+1
                    //
                    // using MinBC optimizer with:
                    // * numerical differentiation being used
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see minbcsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_deleting_element(ref s);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_deleting_element(ref bndu);
                    alglib.minbcstate state;
                    double epsg = 0;
                    if( _spoil_scenario==11 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==12 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==14 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==15 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==16 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.000001;
                    if( _spoil_scenario==17 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==18 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    double diffstep = 1.0e-6;
                    if( _spoil_scenario==20 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==21 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==22 )
                        diffstep = (double)System.Double.NegativeInfinity;

                    //
                    // Now we are ready to actually optimize something:
                    // * first we create optimizer
                    // * we add boundary constraints
                    // * we tune stopping conditions
                    // * and, finally, optimize and obtain results...
                    //
                    alglib.minbccreatef(x, diffstep, out state);
                    alglib.minbcsetbc(state, bndl, bndu);
                    alglib.minbcsetscale(state, s);
                    alglib.minbcsetcond(state, epsg, epsf, epsx, maxits);

                    //
                    // Then we activate OptGuard integrity checking.
                    //
                    // Numerical differentiation always produces "correct" gradient
                    // (with some truncation error, but unbiased). Thus, we just have
                    // to check smoothness properties of the target: C0 and C1 continuity.
                    //
                    // Sometimes user accidentally tries to solve nonsmooth problems
                    // with smooth optimizer. OptGuard helps to detect such situations
                    // early, at the prototyping stage.
                    //
                    alglib.minbcoptguardsmoothness(state);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.minbcreport rep;
                    alglib.minbcoptimize(state, function1_func, null, null);
                    alglib.minbcresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,1}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // Want to challenge OptGuard? Try to make your problem
                    // nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    // re-run optimizer.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minbcoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minbc_numdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minns_d_unconstrained
            //      Nonsmooth unconstrained optimization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = 2*|x0|+|x1|
                    //
                    // using nonsmooth nonlinear optimizer.
                    //
                    double[] x0 = new double[]{1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.00001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double radius = 0.1;
                    if( _spoil_scenario==9 )
                        radius = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        radius = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        radius = (double)System.Double.NegativeInfinity;
                    double rho = 0.0;
                    if( _spoil_scenario==12 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        rho = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnsstate state;
                    alglib.minnsreport rep;
                    double[] x1;

                    //
                    // Create optimizer object, choose AGS algorithm and tune its settings:
                    // * radius=0.1     good initial value; will be automatically decreased later.
                    // * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    //                  because we do not have such constraints
                    // * epsx=0.000001  stopping conditions
                    // * s=[1,1]        all variables have unit scale
                    //
                    alglib.minnscreate(2, x0, out state);
                    alglib.minnssetalgoags(state, radius, rho);
                    alglib.minnssetcond(state, epsx, maxits);
                    alglib.minnssetscale(state, s);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints
                    // (box/linear ones are passed separately by means of minnssetbc() and
                    // minnssetlc() calls).
                    //
                    // If you do not have nonlinear constraints (exactly our situation), then
                    // you will have one-component function vector and 1xN Jacobian matrix.
                    //
                    // So, our vector function has form
                    //
                    //     {f0} = { 2*|x0|+|x1| }
                    //
                    // with Jacobian
                    //
                    //         [                       ]
                    //     J = [ 2*sign(x0)   sign(x1) ]
                    //         [                       ]
                    //
                    // NOTE: nonsmooth optimizer requires considerably more function
                    //       evaluations than smooth solver - about 2N times more. Using
                    //       numerical differentiation introduces additional (multiplicative)
                    //       2N speedup.
                    //
                    //       It means that if smooth optimizer WITH user-supplied gradient
                    //       needs 100 function evaluations to solve 50-dimensional problem,
                    //       then AGS solver with user-supplied gradient will need about 10.000
                    //       function evaluations, and with numerical gradient about 1.000.000
                    //       function evaluations will be performed.
                    //
                    // NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    //       optimization problems. It has convergence guarantees, i.e. it will
                    //       converge to stationary point of the function after running for some
                    //       time.
                    //
                    //       However, it is important to remember that "stationary point" is not
                    //       equal to "solution". If your problem is convex, everything is OK.
                    //       But nonconvex optimization problems may have "flat spots" - large
                    //       areas where gradient is exactly zero, but function value is far away
                    //       from optimal. Such areas are stationary points too, and optimizer
                    //       may be trapped here.
                    //
                    //       "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    //       orders of magnitude worse properties - they may be quite large and
                    //       hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    //       problem, because it is impossible to automatically distinguish "flat
                    //       spot" from true solution.
                    //
                    //       This note is here to warn you that you should be very careful when
                    //       you solve nonsmooth optimization problems. Visual inspection of
                    //       results is essential.
                    //
                    alglib.minnsoptimize(state, nsfunc1_jac, null, null);
                    alglib.minnsresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{0.0000,0.0000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_unconstrained");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minns_d_diff
            //      Nonsmooth unconstrained optimization with numerical differentiation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = 2*|x0|+|x1|
                    //
                    // using nonsmooth nonlinear optimizer with numerical
                    // differentiation provided by ALGLIB.
                    //
                    // NOTE: nonsmooth optimizer requires considerably more function
                    //       evaluations than smooth solver - about 2N times more. Using
                    //       numerical differentiation introduces additional (multiplicative)
                    //       2N speedup.
                    //
                    //       It means that if smooth optimizer WITH user-supplied gradient
                    //       needs 100 function evaluations to solve 50-dimensional problem,
                    //       then AGS solver with user-supplied gradient will need about 10.000
                    //       function evaluations, and with numerical gradient about 1.000.000
                    //       function evaluations will be performed.
                    //
                    double[] x0 = new double[]{1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.00001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double diffstep = 0.000001;
                    if( _spoil_scenario==9 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        diffstep = (double)System.Double.NegativeInfinity;
                    double radius = 0.1;
                    if( _spoil_scenario==12 )
                        radius = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        radius = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        radius = (double)System.Double.NegativeInfinity;
                    double rho = 0.0;
                    if( _spoil_scenario==15 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==16 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        rho = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnsstate state;
                    alglib.minnsreport rep;
                    double[] x1;

                    //
                    // Create optimizer object, choose AGS algorithm and tune its settings:
                    // * radius=0.1     good initial value; will be automatically decreased later.
                    // * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    //                  because we do not have such constraints
                    // * epsx=0.000001  stopping conditions
                    // * s=[1,1]        all variables have unit scale
                    //
                    alglib.minnscreatef(2, x0, diffstep, out state);
                    alglib.minnssetalgoags(state, radius, rho);
                    alglib.minnssetcond(state, epsx, maxits);
                    alglib.minnssetscale(state, s);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function, with first component
                    // being target function, and next components being nonlinear equality
                    // and inequality constraints (box/linear ones are passed separately
                    // by means of minnssetbc() and minnssetlc() calls).
                    //
                    // If you do not have nonlinear constraints (exactly our situation), then
                    // you will have one-component function vector.
                    //
                    // So, our vector function has form
                    //
                    //     {f0} = { 2*|x0|+|x1| }
                    //
                    alglib.minnsoptimize(state, nsfunc1_fvec, null, null);
                    alglib.minnsresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{0.0000,0.0000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_diff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minns_d_bc
            //      Nonsmooth box constrained optimization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<17; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = 2*|x0|+|x1|
                    //
                    // subject to box constraints
                    //
                    //        1 <= x0 < +INF
                    //     -INF <= x1 < +INF
                    //
                    // using nonsmooth nonlinear optimizer.
                    //
                    double[] x0 = new double[]{1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double[] bndl = new double[]{1,-System.Double.PositiveInfinity};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    double[] bndu = new double[]{System.Double.PositiveInfinity,+System.Double.PositiveInfinity};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    double epsx = 0.00001;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==9 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double radius = 0.1;
                    if( _spoil_scenario==11 )
                        radius = (double)System.Double.NaN;
                    if( _spoil_scenario==12 )
                        radius = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        radius = (double)System.Double.NegativeInfinity;
                    double rho = 0.0;
                    if( _spoil_scenario==14 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==15 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==16 )
                        rho = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnsstate state;
                    alglib.minnsreport rep;
                    double[] x1;

                    //
                    // Create optimizer object, choose AGS algorithm and tune its settings:
                    // * radius=0.1     good initial value; will be automatically decreased later.
                    // * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    //                  because we do not have such constraints
                    // * epsx=0.000001  stopping conditions
                    // * s=[1,1]        all variables have unit scale
                    //
                    alglib.minnscreate(2, x0, out state);
                    alglib.minnssetalgoags(state, radius, rho);
                    alglib.minnssetcond(state, epsx, maxits);
                    alglib.minnssetscale(state, s);

                    //
                    // Set box constraints.
                    //
                    // General linear constraints are set in similar way (see comments on
                    // minnssetlc() function for more information).
                    //
                    // You may combine box, linear and nonlinear constraints in one optimization
                    // problem.
                    //
                    alglib.minnssetbc(state, bndl, bndu);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints
                    // (box/linear ones are passed separately by means of minnssetbc() and
                    // minnssetlc() calls).
                    //
                    // If you do not have nonlinear constraints (exactly our situation), then
                    // you will have one-component function vector and 1xN Jacobian matrix.
                    //
                    // So, our vector function has form
                    //
                    //     {f0} = { 2*|x0|+|x1| }
                    //
                    // with Jacobian
                    //
                    //         [                       ]
                    //     J = [ 2*sign(x0)   sign(x1) ]
                    //         [                       ]
                    //
                    // NOTE: nonsmooth optimizer requires considerably more function
                    //       evaluations than smooth solver - about 2N times more. Using
                    //       numerical differentiation introduces additional (multiplicative)
                    //       2N speedup.
                    //
                    //       It means that if smooth optimizer WITH user-supplied gradient
                    //       needs 100 function evaluations to solve 50-dimensional problem,
                    //       then AGS solver with user-supplied gradient will need about 10.000
                    //       function evaluations, and with numerical gradient about 1.000.000
                    //       function evaluations will be performed.
                    //
                    // NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    //       optimization problems. It has convergence guarantees, i.e. it will
                    //       converge to stationary point of the function after running for some
                    //       time.
                    //
                    //       However, it is important to remember that "stationary point" is not
                    //       equal to "solution". If your problem is convex, everything is OK.
                    //       But nonconvex optimization problems may have "flat spots" - large
                    //       areas where gradient is exactly zero, but function value is far away
                    //       from optimal. Such areas are stationary points too, and optimizer
                    //       may be trapped here.
                    //
                    //       "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    //       orders of magnitude worse properties - they may be quite large and
                    //       hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    //       problem, because it is impossible to automatically distinguish "flat
                    //       spot" from true solution.
                    //
                    //       This note is here to warn you that you should be very careful when
                    //       you solve nonsmooth optimization problems. Visual inspection of
                    //       results is essential.
                    //
                    //
                    alglib.minnsoptimize(state, nsfunc1_jac, null, null);
                    alglib.minnsresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{1.0000,0.0000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_bc");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minns_d_nlc
            //      Nonsmooth nonlinearly constrained optimization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x0,x1) = 2*|x0|+|x1|
                    //
                    // subject to combination of equality and inequality constraints
                    //
                    //      x0  =  1
                    //      x1 >= -1
                    //
                    // using nonsmooth nonlinear optimizer. Although these constraints
                    // are linear, we treat them as general nonlinear ones in order to
                    // demonstrate nonlinearly constrained optimization setup.
                    //
                    double[] x0 = new double[]{1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.00001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double radius = 0.1;
                    if( _spoil_scenario==9 )
                        radius = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        radius = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        radius = (double)System.Double.NegativeInfinity;
                    double rho = 50.0;
                    if( _spoil_scenario==12 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        rho = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minnsstate state;
                    alglib.minnsreport rep;
                    double[] x1;

                    //
                    // Create optimizer object, choose AGS algorithm and tune its settings:
                    // * radius=0.1     good initial value; will be automatically decreased later.
                    // * rho=50.0       penalty coefficient for nonlinear constraints. It is your
                    //                  responsibility to choose good one - large enough that it
                    //                  enforces constraints, but small enough in order to avoid
                    //                  extreme slowdown due to ill-conditioning.
                    // * epsx=0.000001  stopping conditions
                    // * s=[1,1]        all variables have unit scale
                    //
                    alglib.minnscreate(2, x0, out state);
                    alglib.minnssetalgoags(state, radius, rho);
                    alglib.minnssetcond(state, epsx, maxits);
                    alglib.minnssetscale(state, s);

                    //
                    // Set general nonlinear constraints.
                    //
                    // This part is more tricky than working with box/linear constraints - you
                    // can not "pack" general nonlinear function into double precision array.
                    // That's why minnssetnlc() does not accept constraints itself - only
                    // constraint COUNTS are passed: first parameter is number of equality
                    // constraints, second one is number of inequality constraints.
                    //
                    // As for constraining functions - these functions are passed as part
                    // of problem Jacobian (see below).
                    //
                    // NOTE: MinNS optimizer supports arbitrary combination of boundary, general
                    //       linear and general nonlinear constraints. This example does not
                    //       show how to work with general linear constraints, but you can
                    //       easily find it in documentation on minnlcsetlc() function.
                    //
                    alglib.minnssetnlc(state, 1, 1);

                    //
                    // Optimize and test results.
                    //
                    // Optimizer object accepts vector function and its Jacobian, with first
                    // component (Jacobian row) being target function, and next components
                    // (Jacobian rows) being nonlinear equality and inequality constraints
                    // (box/linear ones are passed separately by means of minnssetbc() and
                    // minnssetlc() calls).
                    //
                    // Nonlinear equality constraints have form Gi(x)=0, inequality ones
                    // have form Hi(x)<=0, so we may have to "normalize" constraints prior
                    // to passing them to optimizer (right side is zero, constraints are
                    // sorted, multiplied by -1 when needed).
                    //
                    // So, our vector function has form
                    //
                    //     {f0,f1,f2} = { 2*|x0|+|x1|,  x0-1, -x1-1 }
                    //
                    // with Jacobian
                    //
                    //         [ 2*sign(x0)   sign(x1) ]
                    //     J = [     1           0     ]
                    //         [     0          -1     ]
                    //
                    // which means that we have optimization problem
                    //
                    //     min{f0} subject to f1=0, f2<=0
                    //
                    // which is essentially same as
                    //
                    //     min { 2*|x0|+|x1| } subject to x0=1, x1>=-1
                    //
                    // NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    //       optimization problems. It has convergence guarantees, i.e. it will
                    //       converge to stationary point of the function after running for some
                    //       time.
                    //
                    //       However, it is important to remember that "stationary point" is not
                    //       equal to "solution". If your problem is convex, everything is OK.
                    //       But nonconvex optimization problems may have "flat spots" - large
                    //       areas where gradient is exactly zero, but function value is far away
                    //       from optimal. Such areas are stationary points too, and optimizer
                    //       may be trapped here.
                    //
                    //       "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    //       orders of magnitude worse properties - they may be quite large and
                    //       hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    //       problem, because it is impossible to automatically distinguish "flat
                    //       spot" from true solution.
                    //
                    //       This note is here to warn you that you should be very careful when
                    //       you solve nonsmooth optimization problems. Visual inspection of
                    //       results is essential.
                    //
                    alglib.minnsoptimize(state, nsfunc2_jac, null, null);
                    alglib.minnsresults(state, out x1, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x1, new double[]{1.0000,0.0000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_nlc");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mincg_d_1
            //      Nonlinear optimization by CG
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<15; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // using nonlinear conjugate gradient method with:
                    // * initial point x=[0,0]
                    // * unit scale being set for all variables (see mincgsetscale for more info)
                    // * stopping criteria set to "terminate after short enough step"
                    // * OptGuard integrity check being used to check problem statement
                    //   for some common errors like nonsmoothness or bad analytic gradient
                    //
                    // First, we create optimizer object and tune its properties
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsg = 0;
                    if( _spoil_scenario==6 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==9 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.mincgstate state;
                    alglib.mincgcreate(x, out state);
                    alglib.mincgsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.mincgsetscale(state, s);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to catch common coding and problem statement
                    // issues, like:
                    // * discontinuity of the target function (C0 continuity violation)
                    // * nonsmoothness of the target function (C1 continuity violation)
                    // * erroneous analytic gradient, i.e. one inconsistent with actual
                    //   change in the target/constraints
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    //
                    //            Other OptGuard checks add moderate overhead, but anyway
                    //            it is better to turn them off when they are not needed.
                    //
                    alglib.mincgoptguardsmoothness(state);
                    alglib.mincgoptguardgradient(state, 0.001);

                    //
                    // Optimize and evaluate results
                    //
                    alglib.mincgreport rep;
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    //       1.0 to some of its components.
                    //
                    alglib.optguardreport ogrep;
                    alglib.mincgoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mincg_d_2
            //      Nonlinear optimization with additional settings and restarts
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    // with nonlinear conjugate gradient method.
                    //
                    // Several advanced techniques are demonstrated:
                    // * upper limit on step size
                    // * restart from new point
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsg = 0;
                    if( _spoil_scenario==6 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==9 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double stpmax = 0.1;
                    if( _spoil_scenario==15 )
                        stpmax = (double)System.Double.NaN;
                    if( _spoil_scenario==16 )
                        stpmax = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        stpmax = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.mincgstate state;
                    alglib.mincgreport rep;

                    // create and tune optimizer
                    alglib.mincgcreate(x, out state);
                    alglib.mincgsetscale(state, s);
                    alglib.mincgsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.mincgsetstpmax(state, stpmax);

                    // Set up OptGuard integrity checker which catches errors
                    // like nonsmooth targets or errors in the analytic gradient.
                    //
                    // OptGuard is essential at the early prototyping stages.
                    //
                    // NOTE: gradient verification needs 3*N additional function
                    //       evaluations; DO NOT USE IT IN THE PRODUCTION CODE
                    //       because it leads to unnecessary slowdown of your app.
                    alglib.mincgoptguardsmoothness(state);
                    alglib.mincgoptguardgradient(state, 0.001);

                    // first run
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // second run - algorithm is restarted with mincgrestartfrom()
                    x = new double[]{10,10};
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.mincgrestartfrom(state, x);
                    alglib.mincgoptimize(state, function1_grad, null, null);
                    alglib.mincgresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    // check OptGuard integrity report. Why do we need it at all?
                    // Well, try breaking the gradient by adding 1.0 to some
                    // of its components - OptGuard should report it as error.
                    // And it may also catch unintended errors too :)
                    alglib.optguardreport ogrep;
                    alglib.mincgoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mincg_numdiff
            //      Nonlinear optimization by CG with numerical differentiation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<18; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of
                    //
                    //     f(x,y) = 100*(x+3)^4+(y-3)^4
                    //
                    // using numerical differentiation to calculate gradient.
                    //
                    // We also show how to use OptGuard integrity checker to catch common
                    // problem statement errors like accidentally specifying nonsmooth target
                    // function.
                    //
                    // First, we set up optimizer...
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsg = 0;
                    if( _spoil_scenario==6 )
                        epsg = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsg = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsg = (double)System.Double.NegativeInfinity;
                    double epsf = 0;
                    if( _spoil_scenario==9 )
                        epsf = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        epsf = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        epsf = (double)System.Double.NegativeInfinity;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double diffstep = 1.0e-6;
                    if( _spoil_scenario==15 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==16 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        diffstep = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.mincgstate state;
                    alglib.mincgcreatef(x, diffstep, out state);
                    alglib.mincgsetcond(state, epsg, epsf, epsx, maxits);
                    alglib.mincgsetscale(state, s);

                    //
                    // Then, we activate OptGuard integrity checking.
                    //
                    // Numerical differentiation always produces "correct" gradient
                    // (with some truncation error, but unbiased). Thus, we just have
                    // to check smoothness properties of the target: C0 and C1 continuity.
                    //
                    // Sometimes user accidentally tried to solve nonsmooth problems
                    // with smooth optimizer. OptGuard helps to detect such situations
                    // early, at the prototyping stage.
                    //
                    alglib.mincgoptguardsmoothness(state);

                    //
                    // Now we are ready to run the optimization
                    //
                    alglib.mincgreport rep;
                    alglib.mincgoptimize(state, function1_func, null, null);
                    alglib.mincgresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,3}, 0.005);

                    //
                    // ...and to check OptGuard integrity report.
                    //
                    // Want to challenge OptGuard? Try to make your problem
                    // nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    // re-run optimizer.
                    //
                    alglib.optguardreport ogrep;
                    alglib.mincgoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc0suspected, false);
                    _TestResult = _TestResult && doc_test_bool(ogrep.nonc1suspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mincg_numdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_v
            //      Nonlinear least squares optimization using function vector only
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    //
                    //     f0(x0,x1) = 10*(x0+3)^2
                    //     f1(x0,x1) = (x1-3)^2
                    //
                    // using "V" mode of the Levenberg-Marquardt optimizer.
                    //
                    // Optimization algorithm uses:
                    // * function vector f[] = {f1,f2}
                    //
                    // No other information (Jacobian, gradient, etc.) is needed.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    //
                    // Create optimizer, tell it to:
                    // * use numerical differentiation with step equal to 0.0001
                    // * use unit scale for all variables (s is a unit vector)
                    // * stop after short enough step (less than epsx)
                    //
                    alglib.minlmcreatev(2, x, 0.0001, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmsetscale(state, s);

                    //
                    // Optimize
                    //
                    alglib.minlmoptimize(state, function1_fvec, null, null);

                    //
                    // Test optimization results
                    //
                    // NOTE: because we use numerical differentiation, we do not
                    //       verify Jacobian correctness - it is always "correct".
                    //       However, if you switch to analytic gradient, consider
                    //       checking it with OptGuard (see other examples).
                    //
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_v");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_vj
            //      Nonlinear least squares optimization using function vector and Jacobian
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    //
                    //     f0(x0,x1) = 10*(x0+3)^2
                    //     f1(x0,x1) = (x1-3)^2
                    //
                    // using "VJ" mode of the Levenberg-Marquardt optimizer.
                    //
                    // Optimization algorithm uses:
                    // * function vector f[] = {f1,f2}
                    // * Jacobian matrix J = {dfi/dxj}.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==6 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;

                    //
                    // Create optimizer, tell it to:
                    // * use analytic gradient provided by user
                    // * use unit scale for all variables (s is a unit vector)
                    // * stop after short enough step (less than epsx)
                    //
                    alglib.minlmcreatevj(2, x, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmsetscale(state, s);

                    //
                    // Activate OptGuard integrity checking.
                    //
                    // OptGuard monitor helps to detect erroneous analytic Jacobian,
                    // i.e. one inconsistent with actual change in the target function.
                    //
                    // OptGuard is essential for early prototyping stages because such
                    // problems often result in premature termination of the optimizer
                    // which is really hard to distinguish from the correct termination.
                    //
                    // IMPORTANT: JACOBIAN VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    //            DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    //
                    alglib.minlmoptguardgradient(state, 0.001);

                    //
                    // Optimize
                    //
                    alglib.minlmoptimize(state, function1_fvec, function1_jac, null, null);

                    //
                    // Test optimization results
                    //
                    alglib.minlmreport rep;
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);

                    //
                    // Check that OptGuard did not report errors
                    //
                    // NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    //       1.0 to some of its components.
                    //
                    // NOTE: unfortunately, specifics of LM optimization do not allow us
                    //       to detect errors like nonsmoothness (like we do with other
                    //       optimizers). So, only Jacobian correctness is verified.
                    //
                    alglib.optguardreport ogrep;
                    alglib.minlmoptguardresults(state, out ogrep);
                    _TestResult = _TestResult && doc_test_bool(ogrep.badgradsuspected, false);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_vj");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_fgh
            //      Nonlinear Hessian-based optimization for general functions
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = 100*(x0+3)^4+(x1-3)^4
                    // using "FGH" mode of the Levenberg-Marquardt optimizer.
                    //
                    // F is treated like a monolitic function without internal structure,
                    // i.e. we do NOT represent it as a sum of squares.
                    //
                    // Optimization algorithm uses:
                    // * function value F(x0,x1)
                    // * gradient G={dF/dxi}
                    // * Hessian H={d2F/(dxi*dxj)}
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    alglib.minlmcreatefgh(x, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_grad, function1_hess, null, null);
                    alglib.minlmresults(state, out x, out rep);

                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_fgh");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_vb
            //      Bound constrained nonlinear least squares optimization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    //
                    //     f0(x0,x1) = 10*(x0+3)^2
                    //     f1(x0,x1) = (x1-3)^2
                    //
                    // with boundary constraints
                    //
                    //     -1 <= x0 <= +1
                    //     -1 <= x1 <= +1
                    //
                    // using "V" mode of the Levenberg-Marquardt optimizer.
                    //
                    // Optimization algorithm uses:
                    // * function vector f[] = {f1,f2}
                    //
                    // No other information (Jacobian, gradient, etc.) is needed.
                    //
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[] s = new double[]{1,1};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    double[] bndl = new double[]{-1,-1};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{+1,+1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==10 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==11 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==12 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;

                    //
                    // Create optimizer, tell it to:
                    // * use numerical differentiation with step equal to 1.0
                    // * use unit scale for all variables (s is a unit vector)
                    // * stop after short enough step (less than epsx)
                    // * set box constraints
                    //
                    alglib.minlmcreatev(2, x, 0.0001, out state);
                    alglib.minlmsetbc(state, bndl, bndu);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmsetscale(state, s);

                    //
                    // Optimize
                    //
                    alglib.minlmoptimize(state, function1_fvec, null, null);

                    //
                    // Test optimization results
                    //
                    // NOTE: because we use numerical differentiation, we do not
                    //       verify Jacobian correctness - it is always "correct".
                    //       However, if you switch to analytic gradient, consider
                    //       checking it with OptGuard (see other examples).
                    //
                    alglib.minlmreport rep;
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-1,+1}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_vb");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_d_restarts
            //      Efficient restarts of LM optimizer
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    //
                    //     f0(x0,x1) = 10*(x0+3)^2
                    //     f1(x0,x1) = (x1-3)^2
                    //
                    // using several starting points and efficient restarts.
                    //
                    double[] x;
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==0 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==1 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==2 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;

                    //
                    // create optimizer using minlmcreatev()
                    //
                    x = new double[]{10,10};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.minlmcreatev(2, x, 0.0001, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmoptimize(state, function1_fvec, null, null);
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);

                    //
                    // restart optimizer using minlmrestartfrom()
                    //
                    // we can use different starting point, different function,
                    // different stopping conditions, but problem size
                    // must remain unchanged.
                    //
                    x = new double[]{4,4};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.minlmrestartfrom(state, x);
                    alglib.minlmoptimize(state, function2_fvec, null, null);
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{0,1}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_restarts");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_t_1
            //      Nonlinear least squares optimization, FJ scheme (obsolete, but supported)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;
                    alglib.minlmcreatefj(2, x, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_jac, null, null);
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_t_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST minlm_t_2
            //      Nonlinear least squares optimization, FGJ scheme (obsolete, but supported)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double epsx = 0.0000000001;
                    if( _spoil_scenario==3 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    alglib.minlmstate state;
                    alglib.minlmreport rep;
                    alglib.minlmcreatefgj(2, x, out state);
                    alglib.minlmsetcond(state, epsx, maxits);
                    alglib.minlmoptimize(state, function1_func, function1_grad, function1_jac, null, null);
                    alglib.minlmresults(state, out x, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{-3,+3}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "minlm_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_d_base
            //      Basic functionality (moments, adev, median, percentile)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double mean;
                    double variance;
                    double skewness;
                    double kurtosis;
                    double adev;
                    double p;
                    double v;

                    //
                    // Here we demonstrate calculation of sample moments
                    // (mean, variance, skewness, kurtosis)
                    //
                    alglib.samplemoments(x, out mean, out variance, out skewness, out kurtosis);
                    _TestResult = _TestResult && doc_test_real(mean, 28.5, 0.01);
                    _TestResult = _TestResult && doc_test_real(variance, 801.1667, 0.01);
                    _TestResult = _TestResult && doc_test_real(skewness, 0.5751, 0.01);
                    _TestResult = _TestResult && doc_test_real(kurtosis, -1.2666, 0.01);

                    //
                    // Average deviation
                    //
                    alglib.sampleadev(x, out adev);
                    _TestResult = _TestResult && doc_test_real(adev, 23.2, 0.01);

                    //
                    // Median and percentile
                    //
                    alglib.samplemedian(x, out v);
                    _TestResult = _TestResult && doc_test_real(v, 20.5, 0.01);
                    p = 0.5;
                    if( _spoil_scenario==3 )
                        p = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        p = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        p = (double)System.Double.NegativeInfinity;
                    alglib.samplepercentile(x, p, out v);
                    _TestResult = _TestResult && doc_test_real(v, 20.5, 0.01);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_base");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_d_c2
            //      Correlation (covariance) between two random variables
            //
            System.Console.WriteLine("50/151");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    //
                    // We have two samples - x and y, and want to measure dependency between them
                    //
                    double[] x = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double v;

                    //
                    // Three dependency measures are calculated:
                    // * covariation
                    // * Pearson correlation
                    // * Spearman rank correlation
                    //
                    v = alglib.cov2(x, y);
                    _TestResult = _TestResult && doc_test_real(v, 82.5, 0.001);
                    v = alglib.pearsoncorr2(x, y);
                    _TestResult = _TestResult && doc_test_real(v, 0.9627, 0.001);
                    v = alglib.spearmancorr2(x, y);
                    _TestResult = _TestResult && doc_test_real(v, 1.000, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_c2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_d_cm
            //      Correlation (covariance) between components of random vector
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // X is a sample matrix:
                    // * I-th row corresponds to I-th observation
                    // * J-th column corresponds to J-th variable
                    //
                    double[,] x = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[,] c;

                    //
                    // Three dependency measures are calculated:
                    // * covariation
                    // * Pearson correlation
                    // * Spearman rank correlation
                    //
                    // Result is stored into C, with C[i,j] equal to correlation
                    // (covariance) between I-th and J-th variables of X.
                    //
                    alglib.covm(x, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{1.80,0.60,-1.40},{0.60,0.70,-0.80},{-1.40,-0.80,14.70}}, 0.01);
                    alglib.pearsoncorrm(x, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{1.000,0.535,-0.272},{0.535,1.000,-0.249},{-0.272,-0.249,1.000}}, 0.01);
                    alglib.spearmancorrm(x, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{1.000,0.556,-0.306},{0.556,1.000,-0.750},{-0.306,-0.750,1.000}}, 0.01);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_cm");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_d_cm2
            //      Correlation (covariance) between two random vectors
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // X and Y are sample matrices:
                    // * I-th row corresponds to I-th observation
                    // * J-th column corresponds to J-th variable
                    //
                    double[,] x = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    double[,] y = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double[,] c;

                    //
                    // Three dependency measures are calculated:
                    // * covariation
                    // * Pearson correlation
                    // * Spearman rank correlation
                    //
                    // Result is stored into C, with C[i,j] equal to correlation
                    // (covariance) between I-th variable of X and J-th variable of Y.
                    //
                    alglib.covm2(x, y, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{4.100,-3.250},{2.450,-1.500},{13.450,-5.750}}, 0.01);
                    alglib.pearsoncorrm2(x, y, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0.519,-0.699},{0.497,-0.518},{0.596,-0.433}}, 0.01);
                    alglib.spearmancorrm2(x, y, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0.541,-0.649},{0.216,-0.433},{0.433,-0.135}}, 0.01);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_cm2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_t_base
            //      Tests ability to detect errors in inputs
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<34; _spoil_scenario++)
            {
                try
                {
                    double mean;
                    double variance;
                    double skewness;
                    double kurtosis;
                    double adev;
                    double p;
                    double v;

                    //
                    // first, we test short form of functions
                    //
                    double[] x1 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x1, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x1, (double)System.Double.NegativeInfinity);
                    alglib.samplemoments(x1, out mean, out variance, out skewness, out kurtosis);
                    double[] x2 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref x2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NegativeInfinity);
                    alglib.sampleadev(x2, out adev);
                    double[] x3 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x3, (double)System.Double.NaN);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref x3, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref x3, (double)System.Double.NegativeInfinity);
                    alglib.samplemedian(x3, out v);
                    double[] x4 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref x4, (double)System.Double.NaN);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref x4, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref x4, (double)System.Double.NegativeInfinity);
                    p = 0.5;
                    if( _spoil_scenario==12 )
                        p = (double)System.Double.NaN;
                    if( _spoil_scenario==13 )
                        p = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==14 )
                        p = (double)System.Double.NegativeInfinity;
                    alglib.samplepercentile(x4, p, out v);

                    //
                    // and then we test full form
                    //
                    double[] x5 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref x5, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref x5, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref x5, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_deleting_element(ref x5);
                    alglib.samplemoments(x5, 10, out mean, out variance, out skewness, out kurtosis);
                    double[] x6 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref x6, (double)System.Double.NaN);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref x6, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref x6, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_deleting_element(ref x6);
                    alglib.sampleadev(x6, 10, out adev);
                    double[] x7 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==23 )
                        spoil_vector_by_value(ref x7, (double)System.Double.NaN);
                    if( _spoil_scenario==24 )
                        spoil_vector_by_value(ref x7, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==25 )
                        spoil_vector_by_value(ref x7, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==26 )
                        spoil_vector_by_deleting_element(ref x7);
                    alglib.samplemedian(x7, 10, out v);
                    double[] x8 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==27 )
                        spoil_vector_by_value(ref x8, (double)System.Double.NaN);
                    if( _spoil_scenario==28 )
                        spoil_vector_by_value(ref x8, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==29 )
                        spoil_vector_by_value(ref x8, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==30 )
                        spoil_vector_by_deleting_element(ref x8);
                    p = 0.5;
                    if( _spoil_scenario==31 )
                        p = (double)System.Double.NaN;
                    if( _spoil_scenario==32 )
                        p = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==33 )
                        p = (double)System.Double.NegativeInfinity;
                    alglib.samplepercentile(x8, 10, p, out v);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_t_base");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST basestat_t_covcorr
            //      Tests ability to detect errors in inputs
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<126; _spoil_scenario++)
            {
                try
                {
                    double v;
                    double[,] c;

                    //
                    // 2-sample short-form cov/corr are tested
                    //
                    double[] x1 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x1, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x1, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x1);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x1);
                    double[] y1 = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y1, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y1, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y1);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y1);
                    v = alglib.cov2(x1, y1);
                    double[] x2 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref x2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref x2);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref x2);
                    double[] y2 = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref y2, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref y2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref y2, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_adding_element(ref y2);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref y2);
                    v = alglib.pearsoncorr2(x2, y2);
                    double[] x3 = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref x3, (double)System.Double.NaN);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref x3, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_value(ref x3, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==23 )
                        spoil_vector_by_adding_element(ref x3);
                    if( _spoil_scenario==24 )
                        spoil_vector_by_deleting_element(ref x3);
                    double[] y3 = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==25 )
                        spoil_vector_by_value(ref y3, (double)System.Double.NaN);
                    if( _spoil_scenario==26 )
                        spoil_vector_by_value(ref y3, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==27 )
                        spoil_vector_by_value(ref y3, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==28 )
                        spoil_vector_by_adding_element(ref y3);
                    if( _spoil_scenario==29 )
                        spoil_vector_by_deleting_element(ref y3);
                    v = alglib.spearmancorr2(x3, y3);

                    //
                    // 2-sample full-form cov/corr are tested
                    //
                    double[] x1a = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==30 )
                        spoil_vector_by_value(ref x1a, (double)System.Double.NaN);
                    if( _spoil_scenario==31 )
                        spoil_vector_by_value(ref x1a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==32 )
                        spoil_vector_by_value(ref x1a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==33 )
                        spoil_vector_by_deleting_element(ref x1a);
                    double[] y1a = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==34 )
                        spoil_vector_by_value(ref y1a, (double)System.Double.NaN);
                    if( _spoil_scenario==35 )
                        spoil_vector_by_value(ref y1a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==36 )
                        spoil_vector_by_value(ref y1a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==37 )
                        spoil_vector_by_deleting_element(ref y1a);
                    v = alglib.cov2(x1a, y1a, 10);
                    double[] x2a = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==38 )
                        spoil_vector_by_value(ref x2a, (double)System.Double.NaN);
                    if( _spoil_scenario==39 )
                        spoil_vector_by_value(ref x2a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==40 )
                        spoil_vector_by_value(ref x2a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==41 )
                        spoil_vector_by_deleting_element(ref x2a);
                    double[] y2a = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==42 )
                        spoil_vector_by_value(ref y2a, (double)System.Double.NaN);
                    if( _spoil_scenario==43 )
                        spoil_vector_by_value(ref y2a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==44 )
                        spoil_vector_by_value(ref y2a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==45 )
                        spoil_vector_by_deleting_element(ref y2a);
                    v = alglib.pearsoncorr2(x2a, y2a, 10);
                    double[] x3a = new double[]{0,1,4,9,16,25,36,49,64,81};
                    if( _spoil_scenario==46 )
                        spoil_vector_by_value(ref x3a, (double)System.Double.NaN);
                    if( _spoil_scenario==47 )
                        spoil_vector_by_value(ref x3a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==48 )
                        spoil_vector_by_value(ref x3a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==49 )
                        spoil_vector_by_deleting_element(ref x3a);
                    double[] y3a = new double[]{0,1,2,3,4,5,6,7,8,9};
                    if( _spoil_scenario==50 )
                        spoil_vector_by_value(ref y3a, (double)System.Double.NaN);
                    if( _spoil_scenario==51 )
                        spoil_vector_by_value(ref y3a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==52 )
                        spoil_vector_by_value(ref y3a, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==53 )
                        spoil_vector_by_deleting_element(ref y3a);
                    v = alglib.spearmancorr2(x3a, y3a, 10);

                    //
                    // vector short-form cov/corr are tested.
                    //
                    double[,] x4 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==54 )
                        spoil_matrix_by_value(ref x4, (double)System.Double.NaN);
                    if( _spoil_scenario==55 )
                        spoil_matrix_by_value(ref x4, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==56 )
                        spoil_matrix_by_value(ref x4, (double)System.Double.NegativeInfinity);
                    alglib.covm(x4, out c);
                    double[,] x5 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==57 )
                        spoil_matrix_by_value(ref x5, (double)System.Double.NaN);
                    if( _spoil_scenario==58 )
                        spoil_matrix_by_value(ref x5, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==59 )
                        spoil_matrix_by_value(ref x5, (double)System.Double.NegativeInfinity);
                    alglib.pearsoncorrm(x5, out c);
                    double[,] x6 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==60 )
                        spoil_matrix_by_value(ref x6, (double)System.Double.NaN);
                    if( _spoil_scenario==61 )
                        spoil_matrix_by_value(ref x6, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==62 )
                        spoil_matrix_by_value(ref x6, (double)System.Double.NegativeInfinity);
                    alglib.spearmancorrm(x6, out c);

                    //
                    // vector full-form cov/corr are tested.
                    //
                    double[,] x7 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==63 )
                        spoil_matrix_by_value(ref x7, (double)System.Double.NaN);
                    if( _spoil_scenario==64 )
                        spoil_matrix_by_value(ref x7, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==65 )
                        spoil_matrix_by_value(ref x7, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==66 )
                        spoil_matrix_by_deleting_row(ref x7);
                    if( _spoil_scenario==67 )
                        spoil_matrix_by_deleting_col(ref x7);
                    alglib.covm(x7, 5, 3, out c);
                    double[,] x8 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==68 )
                        spoil_matrix_by_value(ref x8, (double)System.Double.NaN);
                    if( _spoil_scenario==69 )
                        spoil_matrix_by_value(ref x8, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==70 )
                        spoil_matrix_by_value(ref x8, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==71 )
                        spoil_matrix_by_deleting_row(ref x8);
                    if( _spoil_scenario==72 )
                        spoil_matrix_by_deleting_col(ref x8);
                    alglib.pearsoncorrm(x8, 5, 3, out c);
                    double[,] x9 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==73 )
                        spoil_matrix_by_value(ref x9, (double)System.Double.NaN);
                    if( _spoil_scenario==74 )
                        spoil_matrix_by_value(ref x9, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==75 )
                        spoil_matrix_by_value(ref x9, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==76 )
                        spoil_matrix_by_deleting_row(ref x9);
                    if( _spoil_scenario==77 )
                        spoil_matrix_by_deleting_col(ref x9);
                    alglib.spearmancorrm(x9, 5, 3, out c);

                    //
                    // cross-vector short-form cov/corr are tested.
                    //
                    double[,] x10 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==78 )
                        spoil_matrix_by_value(ref x10, (double)System.Double.NaN);
                    if( _spoil_scenario==79 )
                        spoil_matrix_by_value(ref x10, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==80 )
                        spoil_matrix_by_value(ref x10, (double)System.Double.NegativeInfinity);
                    double[,] y10 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==81 )
                        spoil_matrix_by_value(ref y10, (double)System.Double.NaN);
                    if( _spoil_scenario==82 )
                        spoil_matrix_by_value(ref y10, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==83 )
                        spoil_matrix_by_value(ref y10, (double)System.Double.NegativeInfinity);
                    alglib.covm2(x10, y10, out c);
                    double[,] x11 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==84 )
                        spoil_matrix_by_value(ref x11, (double)System.Double.NaN);
                    if( _spoil_scenario==85 )
                        spoil_matrix_by_value(ref x11, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==86 )
                        spoil_matrix_by_value(ref x11, (double)System.Double.NegativeInfinity);
                    double[,] y11 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==87 )
                        spoil_matrix_by_value(ref y11, (double)System.Double.NaN);
                    if( _spoil_scenario==88 )
                        spoil_matrix_by_value(ref y11, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==89 )
                        spoil_matrix_by_value(ref y11, (double)System.Double.NegativeInfinity);
                    alglib.pearsoncorrm2(x11, y11, out c);
                    double[,] x12 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==90 )
                        spoil_matrix_by_value(ref x12, (double)System.Double.NaN);
                    if( _spoil_scenario==91 )
                        spoil_matrix_by_value(ref x12, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==92 )
                        spoil_matrix_by_value(ref x12, (double)System.Double.NegativeInfinity);
                    double[,] y12 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==93 )
                        spoil_matrix_by_value(ref y12, (double)System.Double.NaN);
                    if( _spoil_scenario==94 )
                        spoil_matrix_by_value(ref y12, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==95 )
                        spoil_matrix_by_value(ref y12, (double)System.Double.NegativeInfinity);
                    alglib.spearmancorrm2(x12, y12, out c);

                    //
                    // cross-vector full-form cov/corr are tested.
                    //
                    double[,] x13 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==96 )
                        spoil_matrix_by_value(ref x13, (double)System.Double.NaN);
                    if( _spoil_scenario==97 )
                        spoil_matrix_by_value(ref x13, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==98 )
                        spoil_matrix_by_value(ref x13, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==99 )
                        spoil_matrix_by_deleting_row(ref x13);
                    if( _spoil_scenario==100 )
                        spoil_matrix_by_deleting_col(ref x13);
                    double[,] y13 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==101 )
                        spoil_matrix_by_value(ref y13, (double)System.Double.NaN);
                    if( _spoil_scenario==102 )
                        spoil_matrix_by_value(ref y13, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==103 )
                        spoil_matrix_by_value(ref y13, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==104 )
                        spoil_matrix_by_deleting_row(ref y13);
                    if( _spoil_scenario==105 )
                        spoil_matrix_by_deleting_col(ref y13);
                    alglib.covm2(x13, y13, 5, 3, 2, out c);
                    double[,] x14 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==106 )
                        spoil_matrix_by_value(ref x14, (double)System.Double.NaN);
                    if( _spoil_scenario==107 )
                        spoil_matrix_by_value(ref x14, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==108 )
                        spoil_matrix_by_value(ref x14, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==109 )
                        spoil_matrix_by_deleting_row(ref x14);
                    if( _spoil_scenario==110 )
                        spoil_matrix_by_deleting_col(ref x14);
                    double[,] y14 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==111 )
                        spoil_matrix_by_value(ref y14, (double)System.Double.NaN);
                    if( _spoil_scenario==112 )
                        spoil_matrix_by_value(ref y14, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==113 )
                        spoil_matrix_by_value(ref y14, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==114 )
                        spoil_matrix_by_deleting_row(ref y14);
                    if( _spoil_scenario==115 )
                        spoil_matrix_by_deleting_col(ref y14);
                    alglib.pearsoncorrm2(x14, y14, 5, 3, 2, out c);
                    double[,] x15 = new double[,]{{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}};
                    if( _spoil_scenario==116 )
                        spoil_matrix_by_value(ref x15, (double)System.Double.NaN);
                    if( _spoil_scenario==117 )
                        spoil_matrix_by_value(ref x15, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==118 )
                        spoil_matrix_by_value(ref x15, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==119 )
                        spoil_matrix_by_deleting_row(ref x15);
                    if( _spoil_scenario==120 )
                        spoil_matrix_by_deleting_col(ref x15);
                    double[,] y15 = new double[,]{{2,3},{2,1},{-1,6},{-9,9},{7,1}};
                    if( _spoil_scenario==121 )
                        spoil_matrix_by_value(ref y15, (double)System.Double.NaN);
                    if( _spoil_scenario==122 )
                        spoil_matrix_by_value(ref y15, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==123 )
                        spoil_matrix_by_value(ref y15, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==124 )
                        spoil_matrix_by_deleting_row(ref y15);
                    if( _spoil_scenario==125 )
                        spoil_matrix_by_deleting_col(ref y15);
                    alglib.spearmancorrm2(x15, y15, 5, 3, 2, out c);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "basestat_t_covcorr");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ssa_d_basic
            //      Simple SSA analysis demo
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate SSA trend/noise separation for some toy problem:
                    // small monotonically growing series X are analyzed with 3-tick window
                    // and "top-K" version of SSA, which selects K largest singular vectors
                    // for analysis, with K=1.
                    //
                    alglib.ssamodel s;
                    double[] x = new double[]{0,0.5,1,1,1.5,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);

                    //
                    // First, we create SSA model, set its properties and add dataset.
                    //
                    // We use window with width=3 and configure model to use direct SSA
                    // algorithm - one which runs exact O(N*W^2) analysis - to extract
                    // one top singular vector. Well, it is toy problem :)
                    //
                    // NOTE: SSA model may store and analyze more than one sequence
                    //       (say, different sequences may correspond to data collected
                    //       from different devices)
                    //
                    alglib.ssacreate(out s);
                    alglib.ssasetwindow(s, 3);
                    alglib.ssaaddsequence(s, x);
                    alglib.ssasetalgotopkdirect(s, 1);

                    //
                    // Now we begin analysis. Internally SSA model stores everything it needs:
                    // data, settings, solvers and so on. Right after first call to analysis-
                    // related function it will analyze dataset, build basis and perform analysis.
                    //
                    // Subsequent calls to analysis functions will reuse previously computed
                    // basis, unless you invalidate it by changing model settings (or dataset).
                    //
                    double[] trend;
                    double[] noise;
                    alglib.ssaanalyzesequence(s, x, out trend, out noise);
                    _TestResult = _TestResult && doc_test_real_vector(trend, new double[]{0.3815,0.5582,0.7810,1.0794,1.5041,2.0105}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_basic");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ssa_d_forecast
            //      Simple SSA forecasting demo
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate SSA forecasting on some toy problem with clearly
                    // visible linear trend and small amount of noise.
                    //
                    alglib.ssamodel s;
                    double[] x = new double[]{0.05,0.96,2.04,3.11,3.97,5.03,5.98,7.02,8.02};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);

                    //
                    // First, we create SSA model, set its properties and add dataset.
                    //
                    // We use window with width=3 and configure model to use direct SSA
                    // algorithm - one which runs exact O(N*W^2) analysis - to extract
                    // two top singular vectors. Well, it is toy problem :)
                    //
                    // NOTE: SSA model may store and analyze more than one sequence
                    //       (say, different sequences may correspond to data collected
                    //       from different devices)
                    //
                    alglib.ssacreate(out s);
                    alglib.ssasetwindow(s, 3);
                    alglib.ssaaddsequence(s, x);
                    alglib.ssasetalgotopkdirect(s, 2);

                    //
                    // Now we begin analysis. Internally SSA model stores everything it needs:
                    // data, settings, solvers and so on. Right after first call to analysis-
                    // related function it will analyze dataset, build basis and perform analysis.
                    //
                    // Subsequent calls to analysis functions will reuse previously computed
                    // basis, unless you invalidate it by changing model settings (or dataset).
                    //
                    // In this example we show how to use ssaforecastlast() function, which
                    // predicts changed in the last sequence of the dataset. If you want to
                    // perform prediction for some other sequence, use ssaforecastsequence().
                    //
                    double[] trend;
                    alglib.ssaforecastlast(s, 3, out trend);

                    //
                    // Well, we expected it to be [9,10,11]. There exists some difference,
                    // which can be explained by the artificial noise in the dataset.
                    //
                    _TestResult = _TestResult && doc_test_real_vector(trend, new double[]{9.0005,9.9322,10.8051}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_forecast");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST ssa_d_realtime
            //      Real-time SSA algorithm with fast incremental updates
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    //
                    // Suppose that you have a constant stream of incoming data, and you want
                    // to regularly perform singular spectral analysis of this stream.
                    //
                    // One full run of direct algorithm costs O(N*Width^2) operations, so
                    // the more points you have, the more it costs to rebuild basis from
                    // scratch.
                    // 
                    // Luckily we have incremental SSA algorithm which can perform quick
                    // updates of already computed basis in O(K*Width^2) ops, where K
                    // is a number of singular vectors extracted. Usually it is orders of
                    // magnitude faster than full update of the basis.
                    //
                    // In this example we start from some initial dataset x0. Then we
                    // start appending elements one by one to the end of the last sequence.
                    //
                    // NOTE: direct algorithm also supports incremental updates, but
                    //       with O(Width^3) cost. Typically K<<Width, so specialized
                    //       incremental algorithm is still faster.
                    //
                    alglib.ssamodel s1;
                    double[,] a1;
                    double[] sv1;
                    int w;
                    int k;
                    double[] x0 = new double[]{0.009,0.976,1.999,2.984,3.977,5.002};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x0, (double)System.Double.NegativeInfinity);
                    alglib.ssacreate(out s1);
                    alglib.ssasetwindow(s1, 3);
                    alglib.ssaaddsequence(s1, x0);

                    // set algorithm to the real-time version of top-K, K=2
                    alglib.ssasetalgotopkrealtime(s1, 2);

                    // one more interesting feature of the incremental algorithm is "power-up" cycle.
                    // even with incremental algorithm initial basis calculation costs O(N*Width^2) ops.
                    // if such startup cost is too high for your real-time app, then you may divide
                    // initial basis calculation across several model updates. It results in better
                    // latency at the price of somewhat lesser precision during first few updates.
                    alglib.ssasetpoweruplength(s1, 3);

                    // now, after we prepared everything, start to add incoming points one by one;
                    // in the real life, of course, we will perform some work between subsequent update
                    // (analyze something, predict, and so on).
                    //
                    // After each append we perform one iteration of the real-time solver. Usually
                    // one iteration is more than enough to update basis. If you have REALLY tight
                    // performance constraints, you may specify fractional amount of iterations,
                    // which means that iteration is performed with required probability.
                    double updateits = 1.0;
                    if( _spoil_scenario==3 )
                        updateits = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        updateits = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        updateits = (double)System.Double.NegativeInfinity;
                    alglib.ssaappendpointandupdate(s1, 5.951, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 7.074, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 7.925, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 8.992, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 9.942, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 11.051, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 11.965, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 13.047, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    alglib.ssaappendpointandupdate(s1, 13.970, updateits);
                    alglib.ssagetbasis(s1, out a1, out sv1, out w, out k);

                    // Ok, we have our basis in a1[] and singular values at sv1[].
                    // But is it good enough? Let's print it.
                    _TestResult = _TestResult && doc_test_real_matrix(a1, new double[,]{{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}}, 0.0005);

                    // Ok, two vectors with 3 components each.
                    // But how to understand that is it really good basis?
                    // Let's compare it with direct SSA algorithm on the entire sequence.
                    alglib.ssamodel s2;
                    double[,] a2;
                    double[] sv2;
                    double[] x2 = new double[]{0.009,0.976,1.999,2.984,3.977,5.002,5.951,7.074,7.925,8.992,9.942,11.051,11.965,13.047,13.970};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NaN);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref x2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref x2, (double)System.Double.NegativeInfinity);
                    alglib.ssacreate(out s2);
                    alglib.ssasetwindow(s2, 3);
                    alglib.ssaaddsequence(s2, x2);
                    alglib.ssasetalgotopkdirect(s2, 2);
                    alglib.ssagetbasis(s2, out a2, out sv2, out w, out k);

                    // it is exactly the same as one calculated with incremental approach!
                    _TestResult = _TestResult && doc_test_real_matrix(a2, new double[,]{{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}}, 0.0005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_realtime");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST linreg_d_basic
            //      Linear regression used to build the very basic model and unpack coefficients
            //
            _TestResult = true;
            try
            {
                //
                // In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
                //
                // We have:
                // * xy - matrix of basic function values (exp(0.5*x)) and expected values
                //
                double[,] xy = new double[,]{{0.606531,1.133719},{0.670320,1.306522},{0.740818,1.504604},{0.818731,1.554663},{0.904837,1.884638},{1.000000,2.072436},{1.105171,2.257285},{1.221403,2.534068},{1.349859,2.622017},{1.491825,2.897713},{1.648721,3.219371}};
                int info;
                int nvars;
                alglib.linearmodel model;
                alglib.lrreport rep;
                double[] c;

                alglib.lrbuildz(xy, 11, 1, out info, out model, out rep);
                _TestResult = _TestResult && doc_test_int(info, 1);
                alglib.lrunpack(model, out c, out nvars);
                _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.98650,0.00000}, 0.00005);
            }
            catch(alglib.alglibexception)
            { _TestResult = false; }
            catch
            { throw; }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "linreg_d_basic");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST filters_d_sma
            //      SMA(k) filter
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate SMA(k) filtering for time series.
                    //
                    double[] x = new double[]{5,6,7,8};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);

                    //
                    // Apply filter.
                    // We should get [5, 5.5, 6.5, 7.5] as result
                    //
                    alglib.filtersma(ref x, 2);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{5,5.5,6.5,7.5}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_sma");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST filters_d_ema
            //      EMA(alpha) filter
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate EMA(0.5) filtering for time series.
                    //
                    double[] x = new double[]{5,6,7,8};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);

                    //
                    // Apply filter.
                    // We should get [5, 5.5, 6.25, 7.125] as result
                    //
                    alglib.filterema(ref x, 0.5);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{5,5.5,6.25,7.125}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_ema");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST filters_d_lrma
            //      LRMA(k) filter
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate LRMA(3) filtering for time series.
                    //
                    double[] x = new double[]{7,8,8,9,12,12};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);

                    //
                    // Apply filter.
                    // We should get [7.0000, 8.0000, 8.1667, 8.8333, 11.6667, 12.5000] as result
                    //    
                    alglib.filterlrma(ref x, 3);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{7.0000,8.0000,8.1667,8.8333,11.6667,12.5000}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_lrma");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mcpd_simple1
            //      Simple unconstrained MCPD model (no entry/exit states)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple MCPD example
                    //
                    // We have a loan portfolio. Our loans can be in one of two states:
                    // * normal loans ("good" ones)
                    // * past due loans ("bad" ones)
                    //
                    // We assume that:
                    // * loans can transition from any state to any other state. In 
                    //   particular, past due loan can become "good" one at any moment 
                    //   with same (fixed) probability. Not realistic, but it is toy example :)
                    // * portfolio size does not change over time
                    //
                    // Thus, we have following model
                    //     state_new = P*state_old
                    // where
                    //         ( p00  p01 )
                    //     P = (          )
                    //         ( p10  p11 )
                    //
                    // We want to model transitions between these two states using MCPD
                    // approach (Markov Chains for Proportional/Population Data), i.e.
                    // to restore hidden transition matrix P using actual portfolio data.
                    // We have:
                    // * poportional data, i.e. proportion of loans in the normal and past 
                    //   due states (not portfolio size measured in some currency, although 
                    //   it is possible to work with population data too)
                    // * two tracks, i.e. two sequences which describe portfolio
                    //   evolution from two different starting states: [1,0] (all loans 
                    //   are "good") and [0.8,0.2] (only 80% of portfolio is in the "good"
                    //   state)
                    //
                    alglib.mcpdstate s;
                    alglib.mcpdreport rep;
                    double[,] p;
                    double[,] track0 = new double[,]{{1.00000,0.00000},{0.95000,0.05000},{0.92750,0.07250},{0.91738,0.08263},{0.91282,0.08718}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.NegativeInfinity);
                    double[,] track1 = new double[,]{{0.80000,0.20000},{0.86000,0.14000},{0.88700,0.11300},{0.89915,0.10085}};
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.NegativeInfinity);

                    alglib.mcpdcreate(2, out s);
                    alglib.mcpdaddtrack(s, track0);
                    alglib.mcpdaddtrack(s, track1);
                    alglib.mcpdsolve(s);
                    alglib.mcpdresults(s, out p, out rep);

                    //
                    // Hidden matrix P is equal to
                    //         ( 0.95  0.50 )
                    //         (            )
                    //         ( 0.05  0.50 )
                    // which means that "good" loans can become "bad" with 5% probability, 
                    // while "bad" loans will return to good state with 50% probability.
                    //
                    _TestResult = _TestResult && doc_test_real_matrix(p, new double[,]{{0.95,0.50},{0.05,0.50}}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mcpd_simple1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST mcpd_simple2
            //      Simple MCPD model (no entry/exit states) with equality constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // Simple MCPD example
                    //
                    // We have a loan portfolio. Our loans can be in one of three states:
                    // * normal loans
                    // * past due loans
                    // * charged off loans
                    //
                    // We assume that:
                    // * normal loan can stay normal or become past due (but not charged off)
                    // * past due loan can stay past due, become normal or charged off
                    // * charged off loan will stay charged off for the rest of eternity
                    // * portfolio size does not change over time
                    // Not realistic, but it is toy example :)
                    //
                    // Thus, we have following model
                    //     state_new = P*state_old
                    // where
                    //         ( p00  p01    )
                    //     P = ( p10  p11    )
                    //         (      p21  1 )
                    // i.e. four elements of P are known a priori.
                    //
                    // Although it is possible (given enough data) to In order to enforce 
                    // this property we set equality constraints on these elements.
                    //
                    // We want to model transitions between these two states using MCPD
                    // approach (Markov Chains for Proportional/Population Data), i.e.
                    // to restore hidden transition matrix P using actual portfolio data.
                    // We have:
                    // * poportional data, i.e. proportion of loans in the current and past 
                    //   due states (not portfolio size measured in some currency, although 
                    //   it is possible to work with population data too)
                    // * two tracks, i.e. two sequences which describe portfolio
                    //   evolution from two different starting states: [1,0,0] (all loans 
                    //   are "good") and [0.8,0.2,0.0] (only 80% of portfolio is in the "good"
                    //   state)
                    //
                    alglib.mcpdstate s;
                    alglib.mcpdreport rep;
                    double[,] p;
                    double[,] track0 = new double[,]{{1.000000,0.000000,0.000000},{0.950000,0.050000,0.000000},{0.927500,0.060000,0.012500},{0.911125,0.061375,0.027500},{0.896256,0.060900,0.042844}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref track0, (double)System.Double.NegativeInfinity);
                    double[,] track1 = new double[,]{{0.800000,0.200000,0.000000},{0.860000,0.090000,0.050000},{0.862000,0.065500,0.072500},{0.851650,0.059475,0.088875},{0.838805,0.057451,0.103744}};
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref track1, (double)System.Double.NegativeInfinity);

                    alglib.mcpdcreate(3, out s);
                    alglib.mcpdaddtrack(s, track0);
                    alglib.mcpdaddtrack(s, track1);
                    alglib.mcpdaddec(s, 0, 2, 0.0);
                    alglib.mcpdaddec(s, 1, 2, 0.0);
                    alglib.mcpdaddec(s, 2, 2, 1.0);
                    alglib.mcpdaddec(s, 2, 0, 0.0);
                    alglib.mcpdsolve(s);
                    alglib.mcpdresults(s, out p, out rep);

                    //
                    // Hidden matrix P is equal to
                    //         ( 0.95 0.50      )
                    //         ( 0.05 0.25      )
                    //         (      0.25 1.00 ) 
                    // which means that "good" loans can become past due with 5% probability, 
                    // while past due loans will become charged off with 25% probability or
                    // return back to normal state with 50% probability.
                    //
                    _TestResult = _TestResult && doc_test_real_matrix(p, new double[,]{{0.95,0.50,0.00},{0.05,0.25,0.00},{0.00,0.25,1.00}}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "mcpd_simple2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_regr
            //      Regression problem with one output (2=>1)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple example on neural network: network is trained to reproduce
                    // small 2x2 multiplication table.
                    //
                    // NOTE: we use network with excessive amount of neurons, which guarantees
                    //       almost exact reproduction of the training set. Generalization ability
                    //       of such network is rather low, but we are not concerned with such
                    //       questions in this basic demo.
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpreport rep;

                    //
                    // Training set:
                    // * one row corresponds to one record A*B=C in the multiplication table
                    // * first two columns store A and B, last column stores C
                    //
                    // [1 * 1 = 1]
                    // [1 * 2 = 2]
                    // [2 * 1 = 2]
                    // [2 * 2 = 4]
                    //
                    double[,] xy = new double[,]{{1,1,1},{1,2,2},{2,1,2},{2,2,4}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    // Network is created.
                    // Trainer object is created.
                    // Dataset is attached to trainer object.
                    //
                    alglib.mlpcreatetrainer(2, 1, out trn);
                    alglib.mlpcreate1(2, 5, 1, out network);
                    alglib.mlpsetdataset(trn, xy, 4);

                    //
                    // Network is trained with 5 restarts from random positions
                    //
                    alglib.mlptrainnetwork(trn, network, 5, out rep);

                    //
                    // 2*2=?
                    //
                    double[] x = new double[]{2,2};
                    double[] y = new double[]{0};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{4.000}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_regr");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_regr_n
            //      Regression problem with multiple outputs (2=>2)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Network with 2 inputs and 2 outputs is trained to reproduce vector function:
                    //     (x0,x1) => (x0+x1, x0*x1)
                    //
                    // Informally speaking, we want neural network to simultaneously calculate
                    // both sum of two numbers and their product.
                    //
                    // NOTE: we use network with excessive amount of neurons, which guarantees
                    //       almost exact reproduction of the training set. Generalization ability
                    //       of such network is rather low, but we are not concerned with such
                    //       questions in this basic demo.
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpreport rep;

                    //
                    // Training set. One row corresponds to one record [A,B,A+B,A*B].
                    //
                    // [ 1   1  1+1  1*1 ]
                    // [ 1   2  1+2  1*2 ]
                    // [ 2   1  2+1  2*1 ]
                    // [ 2   2  2+2  2*2 ]
                    //
                    double[,] xy = new double[,]{{1,1,2,1},{1,2,3,2},{2,1,3,2},{2,2,4,4}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    // Network is created.
                    // Trainer object is created.
                    // Dataset is attached to trainer object.
                    //
                    alglib.mlpcreatetrainer(2, 2, out trn);
                    alglib.mlpcreate1(2, 5, 2, out network);
                    alglib.mlpsetdataset(trn, xy, 4);

                    //
                    // Network is trained with 5 restarts from random positions
                    //
                    alglib.mlptrainnetwork(trn, network, 5, out rep);

                    //
                    // 2+1=?
                    // 2*1=?
                    //
                    double[] x = new double[]{2,1};
                    double[] y = new double[]{0,0};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{3.000,2.000}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_regr_n");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_cls2
            //      Binary classification problem
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Suppose that we want to classify numbers as positive (class 0) and negative
                    // (class 1). We have training set which includes several strictly positive
                    // or negative numbers - and zero.
                    //
                    // The problem is that we are not sure how to classify zero, so from time to
                    // time we mark it as positive or negative (with equal probability). Other
                    // numbers are marked in pure deterministic setting. How will neural network
                    // cope with such classification task?
                    //
                    // NOTE: we use network with excessive amount of neurons, which guarantees
                    //       almost exact reproduction of the training set. Generalization ability
                    //       of such network is rather low, but we are not concerned with such
                    //       questions in this basic demo.
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpreport rep;
                    double[] x = new double[]{0};
                    double[] y = new double[]{0,0};

                    //
                    // Training set. One row corresponds to one record [A => class(A)].
                    //
                    // Classes are denoted by numbers from 0 to 1, where 0 corresponds to positive
                    // numbers and 1 to negative numbers.
                    //
                    // [ +1  0]
                    // [ +2  0]
                    // [ -1  1]
                    // [ -2  1]
                    // [  0  0]   !! sometimes we classify 0 as positive, sometimes as negative
                    // [  0  1]   !!
                    //
                    double[,] xy = new double[,]{{+1,0},{+2,0},{-1,1},{-2,1},{0,0},{0,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    //
                    // When we solve classification problems, everything is slightly different from
                    // the regression ones:
                    //
                    // 1. Network is created. Because we solve classification problem, we use
                    //    mlpcreatec1() function instead of mlpcreate1(). This function creates
                    //    classifier network with SOFTMAX-normalized outputs. This network returns
                    //    vector of class membership probabilities which are normalized to be
                    //    non-negative and sum to 1.0
                    //
                    // 2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
                    //    create trainer object. Trainer object process dataset and neural network
                    //    slightly differently to account for specifics of the classification
                    //    problems.
                    //
                    // 3. Dataset is attached to trainer object. Note that dataset format is slightly
                    //    different from one used for regression.
                    //
                    alglib.mlpcreatetrainercls(1, 2, out trn);
                    alglib.mlpcreatec1(1, 5, 2, out network);
                    alglib.mlpsetdataset(trn, xy, 6);

                    //
                    // Network is trained with 5 restarts from random positions
                    //
                    alglib.mlptrainnetwork(trn, network, 5, out rep);

                    //
                    // Test our neural network on strictly positive and strictly negative numbers.
                    //
                    // IMPORTANT! Classifier network returns class membership probabilities instead
                    // of class indexes. Network returns two values (probabilities) instead of one
                    // (class index).
                    //
                    // Thus, for +1 we expect to get [P0,P1] = [1,0], where P0 is probability that
                    // number is positive (belongs to class 0), and P1 is probability that number
                    // is negative (belongs to class 1).
                    //
                    // For -1 we expect to get [P0,P1] = [0,1]
                    //
                    // Following properties are guaranteed by network architecture:
                    // * P0>=0, P1>=0   non-negativity
                    // * P0+P1=1        normalization
                    //
                    x = new double[]{1};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{1.000,0.000}, 0.05);
                    x = new double[]{-1};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,1.000}, 0.05);

                    //
                    // But what our network will return for 0, which is between classes 0 and 1?
                    //
                    // In our dataset it has two different marks assigned (class 0 AND class 1).
                    // So network will return something average between class 0 and class 1:
                    //     0 => [0.5, 0.5]
                    //
                    x = new double[]{0};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.500,0.500}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_cls2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_cls3
            //      Multiclass classification problem
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // Suppose that we want to classify numbers as positive (class 0) and negative
                    // (class 1). We also have one more class for zero (class 2).
                    //
                    // NOTE: we use network with excessive amount of neurons, which guarantees
                    //       almost exact reproduction of the training set. Generalization ability
                    //       of such network is rather low, but we are not concerned with such
                    //       questions in this basic demo.
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpreport rep;
                    double[] x = new double[]{0};
                    double[] y = new double[]{0,0,0};

                    //
                    // Training set. One row corresponds to one record [A => class(A)].
                    //
                    // Classes are denoted by numbers from 0 to 2, where 0 corresponds to positive
                    // numbers, 1 to negative numbers, 2 to zero
                    //
                    // [ +1  0]
                    // [ +2  0]
                    // [ -1  1]
                    // [ -2  1]
                    // [  0  2]
                    //
                    double[,] xy = new double[,]{{+1,0},{+2,0},{-1,1},{-2,1},{0,2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    //
                    // When we solve classification problems, everything is slightly different from
                    // the regression ones:
                    //
                    // 1. Network is created. Because we solve classification problem, we use
                    //    mlpcreatec1() function instead of mlpcreate1(). This function creates
                    //    classifier network with SOFTMAX-normalized outputs. This network returns
                    //    vector of class membership probabilities which are normalized to be
                    //    non-negative and sum to 1.0
                    //
                    // 2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
                    //    create trainer object. Trainer object process dataset and neural network
                    //    slightly differently to account for specifics of the classification
                    //    problems.
                    //
                    // 3. Dataset is attached to trainer object. Note that dataset format is slightly
                    //    different from one used for regression.
                    //
                    alglib.mlpcreatetrainercls(1, 3, out trn);
                    alglib.mlpcreatec1(1, 5, 3, out network);
                    alglib.mlpsetdataset(trn, xy, 5);

                    //
                    // Network is trained with 5 restarts from random positions
                    //
                    alglib.mlptrainnetwork(trn, network, 5, out rep);

                    //
                    // Test our neural network on strictly positive and strictly negative numbers.
                    //
                    // IMPORTANT! Classifier network returns class membership probabilities instead
                    // of class indexes. Network returns three values (probabilities) instead of one
                    // (class index).
                    //
                    // Thus, for +1 we expect to get [P0,P1,P2] = [1,0,0],
                    // for -1 we expect to get [P0,P1,P2] = [0,1,0],
                    // and for 0 we will get [P0,P1,P2] = [0,0,1].
                    //
                    // Following properties are guaranteed by network architecture:
                    // * P0>=0, P1>=0, P2>=0    non-negativity
                    // * P0+P1+P2=1             normalization
                    //
                    x = new double[]{1};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{1.000,0.000,0.000}, 0.05);
                    x = new double[]{-1};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,1.000,0.000}, 0.05);
                    x = new double[]{0};
                    alglib.mlpprocess(network, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,0.000,1.000}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_cls3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_trainerobject
            //      Advanced example on trainer object
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // Trainer object is used to train network. It stores dataset, training settings,
                    // and other information which is NOT part of neural network. You should use
                    // trainer object as follows:
                    // (1) you create trainer object and specify task type (classification/regression)
                    //     and number of inputs/outputs
                    // (2) you add dataset to the trainer object
                    // (3) you may change training settings (stopping criteria or weight decay)
                    // (4) finally, you may train one or more networks
                    //
                    // You may interleave stages 2...4 and repeat them many times. Trainer object
                    // remembers its internal state and can be used several times after its creation
                    // and initialization.
                    //
                    alglib.mlptrainer trn;

                    //
                    // Stage 1: object creation.
                    //
                    // We have to specify number of inputs and outputs. Trainer object can be used
                    // only for problems with same number of inputs/outputs as was specified during
                    // its creation.
                    //
                    // In case you want to train SOFTMAX-normalized network which solves classification
                    // problems,  you  must  use  another  function  to  create  trainer  object:
                    // mlpcreatetrainercls().
                    //
                    // Below we create trainer object which can be used to train regression networks
                    // with 2 inputs and 1 output.
                    //
                    alglib.mlpcreatetrainer(2, 1, out trn);

                    //
                    // Stage 2: specification of the training set
                    //
                    // By default trainer object stores empty dataset. So to solve your non-empty problem
                    // you have to set dataset by passing to trainer dense or sparse matrix.
                    //
                    // One row of the matrix corresponds to one record A*B=C in the multiplication table.
                    // First two columns store A and B, last column stores C
                    //
                    //     [1 * 1 = 1]   [ 1 1 1 ]
                    //     [1 * 2 = 2]   [ 1 2 2 ]
                    //     [2 * 1 = 2] = [ 2 1 2 ]
                    //     [2 * 2 = 4]   [ 2 2 4 ]
                    //
                    double[,] xy = new double[,]{{1,1,1},{1,2,2},{2,1,2},{2,2,4}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.mlpsetdataset(trn, xy, 4);

                    //
                    // Stage 3: modification of the training parameters.
                    //
                    // You may modify parameters like weights decay or stopping criteria:
                    // * we set moderate weight decay
                    // * we choose iterations limit as stopping condition (another condition - step size -
                    //   is zero, which means than this condition is not active)
                    //
                    double wstep = 0.000;
                    if( _spoil_scenario==3 )
                        wstep = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        wstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        wstep = (double)System.Double.NegativeInfinity;
                    int maxits = 100;
                    alglib.mlpsetdecay(trn, 0.01);
                    alglib.mlpsetcond(trn, wstep, maxits);

                    //
                    // Stage 4: training.
                    //
                    // We will train several networks with different architecture using same trainer object.
                    // We may change training parameters or even dataset, so different networks are trained
                    // differently. But in this simple example we will train all networks with same settings.
                    //
                    // We create and train three networks:
                    // * network 1 has 2x1 architecture     (2 inputs, no hidden neurons, 1 output)
                    // * network 2 has 2x5x1 architecture   (2 inputs, 5 hidden neurons, 1 output)
                    // * network 3 has 2x5x5x1 architecture (2 inputs, two hidden layers, 1 output)
                    //
                    // NOTE: these networks solve regression problems. For classification problems you
                    //       should use mlpcreatec0/c1/c2 to create neural networks which have SOFTMAX-
                    //       normalized outputs.
                    //
                    alglib.multilayerperceptron net1;
                    alglib.multilayerperceptron net2;
                    alglib.multilayerperceptron net3;
                    alglib.mlpreport rep;

                    alglib.mlpcreate0(2, 1, out net1);
                    alglib.mlpcreate1(2, 5, 1, out net2);
                    alglib.mlpcreate2(2, 5, 5, 1, out net3);

                    alglib.mlptrainnetwork(trn, net1, 5, out rep);
                    alglib.mlptrainnetwork(trn, net2, 5, out rep);
                    alglib.mlptrainnetwork(trn, net3, 5, out rep);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_trainerobject");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_crossvalidation
            //      Cross-validation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example shows how to perform cross-validation with ALGLIB
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpreport rep;

                    //
                    // Training set: f(x)=1/(x^2+1)
                    // One row corresponds to one record [x,f(x)]
                    //
                    double[,] xy = new double[,]{{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    // Trainer object is created.
                    // Dataset is attached to trainer object.
                    //
                    // NOTE: it is not good idea to perform cross-validation on sample
                    //       as small as ours (13 examples). It is done for demonstration
                    //       purposes only. Generalization error estimates won't be
                    //       precise enough for practical purposes.
                    //
                    alglib.mlpcreatetrainer(1, 1, out trn);
                    alglib.mlpsetdataset(trn, xy, 13);

                    //
                    // The key property of the cross-validation is that it estimates
                    // generalization properties of neural ARCHITECTURE. It does NOT
                    // estimates generalization error of some specific network which
                    // is passed to the k-fold CV routine.
                    //
                    // In our example we create 1x4x1 neural network and pass it to
                    // CV routine without training it. Original state of the network
                    // is not used for cross-validation - each round is restarted from
                    // random initial state. Only geometry of network matters.
                    //
                    // We perform 5 restarts from different random positions for each
                    // of the 10 cross-validation rounds.
                    //
                    alglib.mlpcreate1(1, 4, 1, out network);
                    alglib.mlpkfoldcv(trn, network, 5, 10, out rep);

                    //
                    // Cross-validation routine stores estimates of the generalization
                    // error to MLP report structure. You may examine its fields and
                    // see estimates of different errors (RMS, CE, Avg).
                    //
                    // Because cross-validation is non-deterministic, in our manual we
                    // can not say what values will be stored to rep after call to
                    // mlpkfoldcv(). Every CV round will return slightly different
                    // estimates.
                    //
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_crossvalidation");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_ensembles_es
            //      Early stopping ensembles
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example shows how to train early stopping ensebles.
                    //
                    alglib.mlptrainer trn;
                    alglib.mlpensemble ensemble;
                    alglib.mlpreport rep;

                    //
                    // Training set: f(x)=1/(x^2+1)
                    // One row corresponds to one record [x,f(x)]
                    //
                    double[,] xy = new double[,]{{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    //
                    // Trainer object is created.
                    // Dataset is attached to trainer object.
                    //
                    // NOTE: it is not good idea to use early stopping ensemble on sample
                    //       as small as ours (13 examples). It is done for demonstration
                    //       purposes only. Ensemble training algorithm won't find good
                    //       solution on such small sample.
                    //
                    alglib.mlpcreatetrainer(1, 1, out trn);
                    alglib.mlpsetdataset(trn, xy, 13);

                    //
                    // Ensemble is created and trained. Each of 50 network is trained
                    // with 5 restarts.
                    //
                    alglib.mlpecreate1(1, 4, 1, 50, out ensemble);
                    alglib.mlptrainensemblees(trn, ensemble, 5, out rep);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_ensembles_es");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST nn_parallel
            //      Parallel training
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example shows how to use parallel functionality of ALGLIB.
                    // We generate simple 1-dimensional regression problem and show how
                    // to use parallel training, parallel cross-validation, parallel
                    // training of neural ensembles.
                    //
                    // We assume that you already know how to use ALGLIB in serial mode
                    // and concentrate on its parallel capabilities.
                    //
                    // NOTE: it is not good idea to use parallel features on sample as small
                    //       as ours (13 examples). It is done only for demonstration purposes.
                    //
                    alglib.mlptrainer trn;
                    alglib.multilayerperceptron network;
                    alglib.mlpensemble ensemble;
                    alglib.mlpreport rep;
                    double[,] xy = new double[,]{{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.mlpcreatetrainer(1, 1, out trn);
                    alglib.mlpsetdataset(trn, xy, 13);
                    alglib.mlpcreate1(1, 4, 1, out network);
                    alglib.mlpecreate1(1, 4, 1, 50, out ensemble);

                    //
                    // Below we demonstrate how to perform:
                    // * parallel training of individual networks
                    // * parallel cross-validation
                    // * parallel training of neural ensembles
                    //
                    // In order to use multithreading, you have to:
                    // 1) Install SMP edition of ALGLIB.
                    // 2) This step is specific for C++ users: you should activate OS-specific
                    //    capabilities of ALGLIB by defining AE_OS=AE_POSIX (for *nix systems)
                    //    or AE_OS=AE_WINDOWS (for Windows systems).
                    //    C# users do not have to perform this step because C# programs are
                    //    portable across different systems without OS-specific tuning.
                    // 3) Tell ALGLIB that you want it to use multithreading by means of
                    //    setnworkers() call:
                    //          * alglib::setnworkers(0)  = use all cores
                    //          * alglib::setnworkers(-1) = leave one core unused
                    //          * alglib::setnworkers(-2) = leave two cores unused
                    //          * alglib::setnworkers(+2) = use 2 cores (even if you have more)
                    //    During runtime ALGLIB will automatically determine whether it is
                    //    feasible to start worker threads and split your task between cores.
                    //
                    alglib.setnworkers(+2);

                    //
                    // First, we perform parallel training of individual network with 5
                    // restarts from random positions. These 5 rounds of  training  are
                    // executed in parallel manner,  with  best  network  chosen  after
                    // training.
                    //
                    // ALGLIB can use additional way to speed up computations -  divide
                    // dataset   into   smaller   subsets   and   process these subsets
                    // simultaneously. It allows us  to  efficiently  parallelize  even
                    // single training round. This operation is performed automatically
                    // for large datasets, but our toy dataset is too small.
                    //
                    alglib.mlptrainnetwork(trn, network, 5, out rep);

                    //
                    // Then, we perform parallel 10-fold cross-validation, with 5 random
                    // restarts per each CV round. I.e., 5*10=50  networks  are trained
                    // in total. All these operations can be parallelized.
                    //
                    // NOTE: again, ALGLIB can parallelize  calculation   of   gradient
                    //       over entire dataset - but our dataset is too small.
                    //
                    alglib.mlpkfoldcv(trn, network, 5, 10, out rep);

                    //
                    // Finally, we train early stopping ensemble of 50 neural networks,
                    // each  of them is trained with 5 random restarts. I.e.,  5*50=250
                    // networks aretrained in total.
                    //
                    alglib.mlptrainensemblees(trn, ensemble, 5, out rep);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "nn_parallel");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST clst_ahc
            //      Simple hierarchical clusterization with Euclidean distance function
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple clusterization example
                    //
                    // We have a set of points in 2D space:
                    //     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    //
                    //  |
                    //  |     P3
                    //  |
                    //  | P1          
                    //  |             P4
                    //  | P0          P2
                    //  |-------------------------
                    //
                    // We want to perform Agglomerative Hierarchic Clusterization (AHC),
                    // using complete linkage (default algorithm) and Euclidean distance
                    // (default metric).
                    //
                    // In order to do that, we:
                    // * create clusterizer with clusterizercreate()
                    // * set points XY and metric (2=Euclidean) with clusterizersetpoints()
                    // * run AHC algorithm with clusterizerrunahc
                    //
                    // You may see that clusterization itself is a minor part of the example,
                    // most of which is dominated by comments :)
                    //
                    alglib.clusterizerstate s;
                    alglib.ahcreport rep;
                    double[,] xy = new double[,]{{1,1},{1,2},{4,1},{2,3},{4,1.5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.clusterizercreate(out s);
                    alglib.clusterizersetpoints(s, xy, 2);
                    alglib.clusterizerrunahc(s, out rep);

                    //
                    // Now we've built our clusterization tree. Rep.z contains information which
                    // is required to build dendrogram. I-th row of rep.z represents one merge
                    // operation, with first cluster to merge having index rep.z[I,0] and second
                    // one having index rep.z[I,1]. Merge result has index NPoints+I.
                    //
                    // Clusters with indexes less than NPoints are single-point initial clusters,
                    // while ones with indexes from NPoints to 2*NPoints-2 are multi-point
                    // clusters created during merges.
                    //
                    // In our example, Z=[[2,4], [0,1], [3,6], [5,7]]
                    //
                    // It means that:
                    // * first, we merge C2=(P2) and C4=(P4),    and create C5=(P2,P4)
                    // * then, we merge  C2=(P0) and C1=(P1),    and create C6=(P0,P1)
                    // * then, we merge  C3=(P3) and C6=(P0,P1), and create C7=(P0,P1,P3)
                    // * finally, we merge C5 and C7 and create C8=(P0,P1,P2,P3,P4)
                    //
                    // Thus, we have following dendrogram:
                    //  
                    //      ------8-----
                    //      |          |
                    //      |      ----7----
                    //      |      |       |
                    //   ---5---   |    ---6---
                    //   |     |   |    |     |
                    //   P2   P4   P3   P0   P1
                    //
                    _TestResult = _TestResult && doc_test_int_matrix(rep.z, new int[,]{{2,4},{0,1},{3,6},{5,7}});

                    //
                    // We've built dendrogram above by reordering our dataset.
                    //
                    // Without such reordering it would be impossible to build dendrogram without
                    // intersections. Luckily, ahcreport structure contains two additional fields
                    // which help to build dendrogram from your data:
                    // * rep.p, which contains permutation applied to dataset
                    // * rep.pm, which contains another representation of merges 
                    //
                    // In our example we have:
                    // * P=[3,4,0,2,1]
                    // * PZ=[[0,0,1,1,0,0],[3,3,4,4,0,0],[2,2,3,4,0,1],[0,1,2,4,1,2]]
                    //
                    // Permutation array P tells us that P0 should be moved to position 3,
                    // P1 moved to position 4, P2 moved to position 0 and so on:
                    //
                    //   (P0 P1 P2 P3 P4) => (P2 P4 P3 P0 P1)
                    //
                    // Merges array PZ tells us how to perform merges on the sorted dataset.
                    // One row of PZ corresponds to one merge operations, with first pair of
                    // elements denoting first of the clusters to merge (start index, end
                    // index) and next pair of elements denoting second of the clusters to
                    // merge. Clusters being merged are always adjacent, with first one on
                    // the left and second one on the right.
                    //
                    // For example, first row of PZ tells us that clusters [0,0] and [1,1] are
                    // merged (single-point clusters, with first one containing P2 and second
                    // one containing P4). Third row of PZ tells us that we merge one single-
                    // point cluster [2,2] with one two-point cluster [3,4].
                    //
                    // There are two more elements in each row of PZ. These are the helper
                    // elements, which denote HEIGHT (not size) of left and right subdendrograms.
                    // For example, according to PZ, first two merges are performed on clusterization
                    // trees of height 0, while next two merges are performed on 0-1 and 1-2
                    // pairs of trees correspondingly.
                    //
                    _TestResult = _TestResult && doc_test_int_vector(rep.p, new int[]{3,4,0,2,1});
                    _TestResult = _TestResult && doc_test_int_matrix(rep.pm, new int[,]{{0,0,1,1,0,0},{3,3,4,4,0,0},{2,2,3,4,0,1},{0,1,2,4,1,2}});
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "clst_ahc");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST clst_kmeans
            //      Simple k-means clusterization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple clusterization example
                    //
                    // We have a set of points in 2D space:
                    //     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    //
                    //  |
                    //  |     P3
                    //  |
                    //  | P1          
                    //  |             P4
                    //  | P0          P2
                    //  |-------------------------
                    //
                    // We want to perform k-means++ clustering with K=2.
                    //
                    // In order to do that, we:
                    // * create clusterizer with clusterizercreate()
                    // * set points XY and metric (must be Euclidean, distype=2) with clusterizersetpoints()
                    // * (optional) set number of restarts from random positions to 5
                    // * run k-means algorithm with clusterizerrunkmeans()
                    //
                    // You may see that clusterization itself is a minor part of the example,
                    // most of which is dominated by comments :)
                    //
                    alglib.clusterizerstate s;
                    alglib.kmeansreport rep;
                    double[,] xy = new double[,]{{1,1},{1,2},{4,1},{2,3},{4,1.5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.clusterizercreate(out s);
                    alglib.clusterizersetpoints(s, xy, 2);
                    alglib.clusterizersetkmeanslimits(s, 5, 0);
                    alglib.clusterizerrunkmeans(s, 2, out rep);

                    //
                    // We've performed clusterization, and it succeeded (completion code is +1).
                    //
                    // Now first center is stored in the first row of rep.c, second one is stored
                    // in the second row. rep.cidx can be used to determine which center is
                    // closest to some specific point of the dataset.
                    //
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);

                    // We called clusterizersetpoints() with disttype=2 because k-means++
                    // algorithm does NOT support metrics other than Euclidean. But what if we
                    // try to use some other metric?
                    //
                    // We change metric type by calling clusterizersetpoints() one more time,
                    // and try to run k-means algo again. It fails.
                    //
                    alglib.clusterizersetpoints(s, xy, 0);
                    alglib.clusterizerrunkmeans(s, 2, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, -5);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "clst_kmeans");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST clst_linkage
            //      Clusterization with different linkage types
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // We have a set of points in 1D space:
                    //     (P0,P1,P2,P3,P4) = (1, 3, 10, 16, 20)
                    //
                    // We want to perform Agglomerative Hierarchic Clusterization (AHC),
                    // using either complete or single linkage and Euclidean distance
                    // (default metric).
                    //
                    // First two steps merge P0/P1 and P3/P4 independently of the linkage type.
                    // However, third step depends on linkage type being used:
                    // * in case of complete linkage P2=10 is merged with [P0,P1]
                    // * in case of single linkage P2=10 is merged with [P3,P4]
                    //
                    alglib.clusterizerstate s;
                    alglib.ahcreport rep;
                    double[,] xy = new double[,]{{1},{3},{10},{16},{20}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    int[] cidx;
                    int[] cz;

                    alglib.clusterizercreate(out s);
                    alglib.clusterizersetpoints(s, xy, 2);

                    // use complete linkage, reduce set down to 2 clusters.
                    // print clusterization with clusterizergetkclusters(2).
                    // P2 must belong to [P0,P1]
                    alglib.clusterizersetahcalgo(s, 0);
                    alglib.clusterizerrunahc(s, out rep);
                    alglib.clusterizergetkclusters(rep, 2, out cidx, out cz);
                    _TestResult = _TestResult && doc_test_int_vector(cidx, new int[]{1,1,1,0,0});

                    // use single linkage, reduce set down to 2 clusters.
                    // print clusterization with clusterizergetkclusters(2).
                    // P2 must belong to [P2,P3]
                    alglib.clusterizersetahcalgo(s, 1);
                    alglib.clusterizerrunahc(s, out rep);
                    alglib.clusterizergetkclusters(rep, 2, out cidx, out cz);
                    _TestResult = _TestResult && doc_test_int_vector(cidx, new int[]{0,0,1,1,1});
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "clst_linkage");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST clst_distance
            //      Clusterization with different metric types
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // We have three points in 4D space:
                    //     (P0,P1,P2) = ((1, 2, 1, 2), (6, 7, 6, 7), (7, 6, 7, 6))
                    //
                    // We want to try clustering them with different distance functions.
                    // Distance function is chosen when we add dataset to the clusterizer.
                    // We can choose several distance types - Euclidean, city block, Chebyshev,
                    // several correlation measures or user-supplied distance matrix.
                    //
                    // Here we'll try three distances: Euclidean, Pearson correlation,
                    // user-supplied distance matrix. Different distance functions lead
                    // to different choices being made by algorithm during clustering.
                    //
                    alglib.clusterizerstate s;
                    alglib.ahcreport rep;
                    int disttype;
                    double[,] xy = new double[,]{{1,2,1,2},{6,7,6,7},{7,6,7,6}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.clusterizercreate(out s);

                    // With Euclidean distance function (disttype=2) two closest points
                    // are P1 and P2, thus:
                    // * first, we merge P1 and P2 to form C3=[P1,P2]
                    // * second, we merge P0 and C3 to form C4=[P0,P1,P2]
                    disttype = 2;
                    alglib.clusterizersetpoints(s, xy, disttype);
                    alglib.clusterizerrunahc(s, out rep);
                    _TestResult = _TestResult && doc_test_int_matrix(rep.z, new int[,]{{1,2},{0,3}});

                    // With Pearson correlation distance function (disttype=10) situation
                    // is different - distance between P0 and P1 is zero, thus:
                    // * first, we merge P0 and P1 to form C3=[P0,P1]
                    // * second, we merge P2 and C3 to form C4=[P0,P1,P2]
                    disttype = 10;
                    alglib.clusterizersetpoints(s, xy, disttype);
                    alglib.clusterizerrunahc(s, out rep);
                    _TestResult = _TestResult && doc_test_int_matrix(rep.z, new int[,]{{0,1},{2,3}});

                    // Finally, we try clustering with user-supplied distance matrix:
                    //     [ 0 3 1 ]
                    // P = [ 3 0 3 ], where P[i,j] = dist(Pi,Pj)
                    //     [ 1 3 0 ]
                    //
                    // * first, we merge P0 and P2 to form C3=[P0,P2]
                    // * second, we merge P1 and C3 to form C4=[P0,P1,P2]
                    double[,] d = new double[,]{{0,3,1},{3,0,3},{1,3,0}};
                    alglib.clusterizersetdistances(s, d, true);
                    alglib.clusterizerrunahc(s, out rep);
                    _TestResult = _TestResult && doc_test_int_matrix(rep.z, new int[,]{{0,2},{1,3}});
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "clst_distance");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST clst_kclusters
            //      Obtaining K top clusters from clusterization tree
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // We have a set of points in 2D space:
                    //     (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    //
                    //  |
                    //  |     P3
                    //  |
                    //  | P1          
                    //  |             P4
                    //  | P0          P2
                    //  |-------------------------
                    //
                    // We perform Agglomerative Hierarchic Clusterization (AHC) and we want
                    // to get top K clusters from clusterization tree for different K.
                    //
                    alglib.clusterizerstate s;
                    alglib.ahcreport rep;
                    double[,] xy = new double[,]{{1,1},{1,2},{4,1},{2,3},{4,1.5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    int[] cidx;
                    int[] cz;

                    alglib.clusterizercreate(out s);
                    alglib.clusterizersetpoints(s, xy, 2);
                    alglib.clusterizerrunahc(s, out rep);

                    // with K=5, every points is assigned to its own cluster:
                    // C0=P0, C1=P1 and so on...
                    alglib.clusterizergetkclusters(rep, 5, out cidx, out cz);
                    _TestResult = _TestResult && doc_test_int_vector(cidx, new int[]{0,1,2,3,4});

                    // with K=1 we have one large cluster C0=[P0,P1,P2,P3,P4,P5]
                    alglib.clusterizergetkclusters(rep, 1, out cidx, out cz);
                    _TestResult = _TestResult && doc_test_int_vector(cidx, new int[]{0,0,0,0,0});

                    // with K=3 we have three clusters C0=[P3], C1=[P2,P4], C2=[P0,P1]
                    alglib.clusterizergetkclusters(rep, 3, out cidx, out cz);
                    _TestResult = _TestResult && doc_test_int_vector(cidx, new int[]{2,2,1,0,1});
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "clst_kclusters");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST randomforest_cls
            //      Simple classification with random forests
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple classification example: classify points (x,y) in 2D space
                    // as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
                    // has to find out it).
                    //
                    // First, we have to create decision forest builder object, load dataset and
                    // specify training settings. Our dataset is specified as matrix, which has
                    // following format:
                    //
                    //     x0 y0 class0
                    //     x1 y1 class1
                    //     x2 y2 class2
                    //     ....
                    //
                    // Here xi and yi can be any values (and in fact you can have any number of
                    // independent variables), and classi MUST be integer number in [0,NClasses)
                    // range. In our example we denote points with x>=0 as class #0, and
                    // ones with negative xi as class #1.
                    //
                    // NOTE: if you want to solve regression problem, specify NClasses=1. In
                    //       this case last column of xy can be any numeric value.
                    //
                    // For the sake of simplicity, our example includes only 4-point dataset.
                    // However, random forests are able to cope with extremely large datasets
                    // having millions of examples.
                    //
                    alglib.decisionforestbuilder builder;
                    int nvars = 2;
                    int nclasses = 2;
                    int npoints = 4;
                    double[,] xy = new double[,]{{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.dfbuildercreate(out builder);
                    alglib.dfbuildersetdataset(builder, xy, npoints, nvars, nclasses);

                    // in our example we train decision forest using full sample - it allows us
                    // to get zero classification error. However, in practical applications smaller
                    // values are used: 50%, 25%, 5% or even less.
                    alglib.dfbuildersetsubsampleratio(builder, 1.0);

                    // we train random forest with just one tree; again, in real life situations
                    // you typically need from 50 to 500 trees.
                    int ntrees = 1;
                    alglib.decisionforest forest;
                    alglib.dfreport rep;
                    alglib.dfbuilderbuildrandomforest(builder, ntrees, out forest, out rep);

                    // with such settings (100% of the training set is used) you can expect
                    // zero classification error. Beautiful results, but remember - in real life
                    // you do not need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult && doc_test_real(rep.relclserror, 0.0000, 0.00005);

                    // now, let's perform some simple processing with dfprocess()
                    double[] x = new double[]{+1,0};
                    double[] y = new double[0];
                    alglib.dfprocess(forest, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{+1,0}, 0.0005);

                    // another option is to use dfprocess0() which returns just first component
                    // of the output vector y. ideal for regression problems and binary classifiers.
                    double y0;
                    y0 = alglib.dfprocess0(forest, x);
                    _TestResult = _TestResult && doc_test_real(y0, 1.000, 0.0005);

                    // finally, you can use dfclassify() which returns most probable class index (i.e. argmax y[i]).
                    int i;
                    i = alglib.dfclassify(forest, x);
                    _TestResult = _TestResult && doc_test_int(i, 0);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "randomforest_cls");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST randomforest_reg
            //      Simple classification with decision forest
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple regression example: model f(x,y)=x+y
                    //
                    // First, we have to create DF builder object, load dataset and specify
                    // training settings. Our dataset is specified as matrix, which has following
                    // format:
                    //
                    //     x0 y0 f0
                    //     x1 y1 f1
                    //     x2 y2 f2
                    //     ....
                    //
                    // Here xi and yi can be any values, and fi is a dependent function value.
                    //
                    // NOTE: you can also solve classification problems with DF models, see
                    //       another example for this unit.
                    //
                    alglib.decisionforestbuilder builder;
                    int nvars = 2;
                    int nclasses = 1;
                    int npoints = 4;
                    double[,] xy = new double[,]{{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.dfbuildercreate(out builder);
                    alglib.dfbuildersetdataset(builder, xy, npoints, nvars, nclasses);

                    // in our example we train decision forest using full sample - it allows us
                    // to get zero classification error. However, in practical applications smaller
                    // values are used: 50%, 25%, 5% or even less.
                    alglib.dfbuildersetsubsampleratio(builder, 1.0);

                    // we train random forest with just one tree; again, in real life situations
                    // you typically need from 50 to 500 trees.
                    int ntrees = 1;
                    alglib.decisionforest model;
                    alglib.dfreport rep;
                    alglib.dfbuilderbuildrandomforest(builder, ntrees, out model, out rep);

                    // with such settings (full sample is used) you can expect zero RMS error on the
                    // training set. Beautiful results, but remember - in real life you do not
                    // need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult && doc_test_real(rep.rmserror, 0.0000, 0.00005);

                    // now, let's perform some simple processing with dfprocess()
                    double[] x = new double[]{+1,+1};
                    double[] y = new double[0];
                    alglib.dfprocess(model, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{+2}, 0.0005);

                    // another option is to use dfprocess0() which returns just first component
                    // of the output vector y. ideal for regression problems and binary classifiers.
                    double y0;
                    y0 = alglib.dfprocess0(model, x);
                    _TestResult = _TestResult && doc_test_real(y0, 2.000, 0.0005);

                    // there also exist another convenience function, dfclassify(),
                    // but it does not work for regression problems - it always returns -1.
                    int i;
                    i = alglib.dfclassify(model, x);
                    _TestResult = _TestResult && doc_test_int(i, -1);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "randomforest_reg");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST knn_cls
            //      Simple classification with KNN model
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple classification example: classify points (x,y) in 2D space
                    // as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
                    // has to find out it).
                    //
                    // First, we have to create KNN builder object, load dataset and specify
                    // training settings. Our dataset is specified as matrix, which has following
                    // format:
                    //
                    //     x0 y0 class0
                    //     x1 y1 class1
                    //     x2 y2 class2
                    //     ....
                    //
                    // Here xi and yi can be any values (and in fact you can have any number of
                    // independent variables), and classi MUST be integer number in [0,NClasses)
                    // range. In our example we denote points with x>=0 as class #0, and
                    // ones with negative xi as class #1.
                    //
                    // NOTE: if you want to solve regression problem, specify dataset in similar
                    //       format, but with dependent variable(s) instead of class labels. You
                    //       can have dataset with multiple dependent variables, by the way!
                    //
                    // For the sake of simplicity, our example includes only 4-point dataset and
                    // really simple K=1 nearest neighbor search. Industrial problems typically
                    // need larger values of K.
                    //
                    alglib.knnbuilder builder;
                    int nvars = 2;
                    int nclasses = 2;
                    int npoints = 4;
                    double[,] xy = new double[,]{{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.knnbuildercreate(out builder);
                    alglib.knnbuildersetdatasetcls(builder, xy, npoints, nvars, nclasses);

                    // we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
                    int k = 1;
                    double eps = 0;
                    alglib.knnmodel model;
                    alglib.knnreport rep;
                    alglib.knnbuilderbuildknnmodel(builder, k, eps, out model, out rep);

                    // with such settings (k=1 is used) you can expect zero classification
                    // error on training set. Beautiful results, but remember - in real life
                    // you do not need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult && doc_test_real(rep.relclserror, 0.0000, 0.00005);

                    // now, let's perform some simple processing with knnprocess()
                    double[] x = new double[]{+1,0};
                    double[] y = new double[0];
                    alglib.knnprocess(model, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{+1,0}, 0.0005);

                    // another option is to use knnprocess0() which returns just first component
                    // of the output vector y. ideal for regression problems and binary classifiers.
                    double y0;
                    y0 = alglib.knnprocess0(model, x);
                    _TestResult = _TestResult && doc_test_real(y0, 1.000, 0.0005);

                    // finally, you can use knnclassify() which returns most probable class index (i.e. argmax y[i]).
                    int i;
                    i = alglib.knnclassify(model, x);
                    _TestResult = _TestResult && doc_test_int(i, 0);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "knn_cls");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST knn_reg
            //      Simple classification with KNN model
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // The very simple regression example: model f(x,y)=x+y
                    //
                    // First, we have to create KNN builder object, load dataset and specify
                    // training settings. Our dataset is specified as matrix, which has following
                    // format:
                    //
                    //     x0 y0 f0
                    //     x1 y1 f1
                    //     x2 y2 f2
                    //     ....
                    //
                    // Here xi and yi can be any values, and fi is a dependent function value.
                    // By the way, with KNN algorithm you can even model functions with multiple
                    // dependent variables!
                    //
                    // NOTE: you can also solve classification problems with KNN models, see
                    //       another example for this unit.
                    //
                    // For the sake of simplicity, our example includes only 4-point dataset and
                    // really simple K=1 nearest neighbor search. Industrial problems typically
                    // need larger values of K.
                    //
                    alglib.knnbuilder builder;
                    int nvars = 2;
                    int nout = 1;
                    int npoints = 4;
                    double[,] xy = new double[,]{{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);

                    alglib.knnbuildercreate(out builder);
                    alglib.knnbuildersetdatasetreg(builder, xy, npoints, nvars, nout);

                    // we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
                    int k = 1;
                    double eps = 0;
                    alglib.knnmodel model;
                    alglib.knnreport rep;
                    alglib.knnbuilderbuildknnmodel(builder, k, eps, out model, out rep);

                    // with such settings (k=1 is used) you can expect zero RMS error on the
                    // training set. Beautiful results, but remember - in real life you do not
                    // need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult && doc_test_real(rep.rmserror, 0.0000, 0.00005);

                    // now, let's perform some simple processing with knnprocess()
                    double[] x = new double[]{+1,+1};
                    double[] y = new double[0];
                    alglib.knnprocess(model, x, ref y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{+2}, 0.0005);

                    // another option is to use knnprocess0() which returns just first component
                    // of the output vector y. ideal for regression problems and binary classifiers.
                    double y0;
                    y0 = alglib.knnprocess0(model, x);
                    _TestResult = _TestResult && doc_test_real(y0, 2.000, 0.0005);

                    // there also exist another convenience function, knnclassify(),
                    // but it does not work for regression problems - it always returns -1.
                    int i;
                    i = alglib.knnclassify(model, x);
                    _TestResult = _TestResult && doc_test_int(i, -1);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "knn_reg");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST autogk_d1
            //      Integrating f=exp(x) by adaptive integrator
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates integration of f=exp(x) on [0,1]:
                    // * first, autogkstate is initialized
                    // * then we call integration function
                    // * and finally we obtain results with autogkresults() call
                    //
                    double a = 0;
                    if( _spoil_scenario==0 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==1 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==2 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = 1;
                    if( _spoil_scenario==3 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==4 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        b = (double)System.Double.NegativeInfinity;
                    alglib.autogkstate s;
                    double v;
                    alglib.autogkreport rep;

                    alglib.autogksmooth(a, b, out s);
                    alglib.autogkintegrate(s, int_function_1_func, null);
                    alglib.autogkresults(s, out v, out rep);

                    _TestResult = _TestResult && doc_test_real(v, 1.7182, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "autogk_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_d1
            //      Complex FFT: simple example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1i,1i,1i,1i] is converted to [4i, 0, 0, 0]
                    //
                    alglib.complex[] z = new alglib.complex[]{new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1)};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NegativeInfinity);
                    alglib.fftc1d(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(0,4),0,0,0}, 0.0001);

                    //
                    // now we convert [4i, 0, 0, 0] back to [1i,1i,1i,1i]
                    // with backward FFT
                    //
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1),new alglib.complex(0,1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_d2
            //      Complex FFT: advanced example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [0,1,0,1i] is converted to [1+1i, -1-1i, -1-1i, 1+1i]
                    //
                    alglib.complex[] z = new alglib.complex[]{0,1,0,new alglib.complex(0,1)};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NegativeInfinity);
                    alglib.fftc1d(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{new alglib.complex(1,+1),new alglib.complex(-1,-1),new alglib.complex(-1,-1),new alglib.complex(1,+1)}, 0.0001);

                    //
                    // now we convert result back with backward FFT
                    //
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{0,1,0,new alglib.complex(0,1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_real_d1
            //      Real FFT: simple example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1,1,1,1] is converted to [4, 0, 0, 0]
                    //
                    double[] x = new double[]{1,1,1,1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.complex[] f;
                    double[] x2;
                    alglib.fftr1d(x, out f);
                    _TestResult = _TestResult && doc_test_complex_vector(f, new alglib.complex[]{4,0,0,0}, 0.0001);

                    //
                    // now we convert [4, 0, 0, 0] back to [1,1,1,1]
                    // with backward FFT
                    //
                    alglib.fftr1dinv(f, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,1,1,1}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_real_d2
            //      Real FFT: advanced example
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // first we demonstrate forward FFT:
                    // [1,2,3,4] is converted to [10, -2+2i, -2, -2-2i]
                    //
                    // note that output array is self-adjoint:
                    // * f[0] = conj(f[0])
                    // * f[1] = conj(f[3])
                    // * f[2] = conj(f[2])
                    //
                    double[] x = new double[]{1,2,3,4};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.complex[] f;
                    double[] x2;
                    alglib.fftr1d(x, out f);
                    _TestResult = _TestResult && doc_test_complex_vector(f, new alglib.complex[]{10,new alglib.complex(-2,+2),-2,new alglib.complex(-2,-2)}, 0.0001);

                    //
                    // now we convert [10, -2+2i, -2, -2-2i] back to [1,2,3,4]
                    //
                    alglib.fftr1dinv(f, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,2,3,4}, 0.0001);

                    //
                    // remember that F is self-adjoint? It means that we can pass just half
                    // (slightly larger than half) of F to inverse real FFT and still get our result.
                    //
                    // I.e. instead [10, -2+2i, -2, -2-2i] we pass just [10, -2+2i, -2] and everything works!
                    //
                    // NOTE: in this case we should explicitly pass array length (which is 4) to ALGLIB;
                    // if not, it will automatically use array length to determine FFT size and
                    // will erroneously make half-length FFT.
                    //
                    f = new alglib.complex[]{10,new alglib.complex(-2,+2),-2};
                    alglib.fftr1dinv(f, 4, out x2);
                    _TestResult = _TestResult && doc_test_real_vector(x2, new double[]{1,2,3,4}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST fft_complex_e1
            //      error detection in backward FFT
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[] z = new alglib.complex[]{0,2,0,-2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref z, (alglib.complex)System.Double.NegativeInfinity);
                    alglib.fftc1dinv(ref z);
                    _TestResult = _TestResult && doc_test_complex_vector(z, new alglib.complex[]{0,new alglib.complex(0,1),0,new alglib.complex(0,-1)}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_e1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST idw_d_mstab
            //      Simple model built with IDW-MSTAB algorithm
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example illustrates basic concepts of the IDW models:
                    // creation and evaluation.
                    // 
                    // Suppose that we have set of 2-dimensional points with associated
                    // scalar function values, and we want to build an IDW model using
                    // our data.
                    // 
                    // NOTE: we can work with N-dimensional models and vector-valued functions too :)
                    // 
                    // Typical sequence of steps is given below:
                    // 1. we create IDW builder object
                    // 2. we attach our dataset to the IDW builder and tune algorithm settings
                    // 3. we generate IDW model
                    // 4. we use IDW model instance (evaluate, serialize, etc.)
                    //
                    double v;

                    //
                    // Step 1: IDW builder creation.
                    //
                    // We have to specify dimensionality of the space (2 or 3) and
                    // dimensionality of the function (scalar or vector).
                    //
                    // New builder object is empty - it has not dataset and uses
                    // default model construction settings
                    //
                    alglib.idwbuilder builder;
                    alglib.idwbuildercreate(2, 1, out builder);

                    //
                    // Step 2: dataset addition
                    //
                    // XY contains two points - x0=(-1,0) and x1=(+1,0) -
                    // and two function values f(x0)=2, f(x1)=3.
                    //
                    double[,] xy = new double[,]{{-1,0,2},{+1,0,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.idwbuildersetpoints(builder, xy);

                    //
                    // Step 3: choose IDW algorithm and generate model
                    //
                    // We use modified stabilized IDW algorithm with following parameters:
                    // * SRad - set to 5.0 (search radius must be large enough)
                    //
                    // IDW-MSTAB algorithm is a state-of-the-art implementation of IDW which
                    // is competitive with RBFs and bicubic splines. See comments on the
                    // idwbuildersetalgomstab() function for more information.
                    //
                    alglib.idwmodel model;
                    alglib.idwreport rep;
                    alglib.idwbuildersetalgomstab(builder, 5.0);
                    alglib.idwfit(builder, out model, out rep);

                    //
                    // Step 4: model was built, evaluate its value
                    //
                    v = alglib.idwcalc2(model, 1.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 3.000, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "idw_d_mstab");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST idw_d_serialize
            //      IDW model serialization/unserialization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example shows how to serialize and unserialize IDW model.
                    // 
                    // Suppose that we have set of 2-dimensional points with associated
                    // scalar function values, and we have built an IDW model using
                    // our data.
                    //
                    // This model can be serialized to string or stream. ALGLIB supports
                    // flexible (un)serialization, i.e. you can move serialized model
                    // representation between different machines (32-bit or 64-bit),
                    // different CPU architectures (x86/64, ARM) or even different
                    // programming languages supported by ALGLIB (C#, C++, ...).
                    //
                    // Our first step is to build model, evaluate it at point (1,0),
                    // and serialize it to string.
                    //
                    string s;
                    double v;
                    double[,] xy = new double[,]{{-1,0,2},{+1,0,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.idwbuilder builder;
                    alglib.idwmodel model;
                    alglib.idwmodel model2;
                    alglib.idwreport rep;
                    alglib.idwbuildercreate(2, 1, out builder);
                    alglib.idwbuildersetpoints(builder, xy);
                    alglib.idwbuildersetalgomstab(builder, 5.0);
                    alglib.idwfit(builder, out model, out rep);
                    v = alglib.idwcalc2(model, 1.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 3.000, 0.005);

                    //
                    // Serialization + unserialization to a different instance
                    // of the model class.
                    //
                    alglib.idwserialize(model, out s);
                    alglib.idwunserialize(s, out model2);

                    //
                    // Evaluate unserialized model at the same point
                    //
                    v = alglib.idwcalc2(model2, 1.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 3.000, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "idw_d_serialize");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline1d_d_linear
            //      Piecewise linear spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use piecewise linear spline to interpolate f(x)=x^2 sampled 
                    // at 5 equidistant nodes on [-1,+1].
                    //
                    double[] x = new double[]{-1.0,-0.5,0.0,+0.5,+1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{+1.0,0.25,0.0,0.25,+1.0};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = 0.25;
                    if( _spoil_scenario==10 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    alglib.spline1dinterpolant s;

                    // build spline
                    alglib.spline1dbuildlinear(x, y, out s);

                    // calculate S(0.25) - it is quite different from 0.25^2=0.0625
                    v = alglib.spline1dcalc(s, t);
                    _TestResult = _TestResult && doc_test_real(v, 0.125, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_linear");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline1d_d_cubic
            //      Cubic spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use cubic spline to interpolate f(x)=x^2 sampled 
                    // at 5 equidistant nodes on [-1,+1].
                    //
                    // First, we use default boundary conditions ("parabolically terminated
                    // spline") because cubic spline built with such boundary conditions 
                    // will exactly reproduce any quadratic f(x).
                    //
                    // Then we try to use natural boundary conditions
                    //     d2S(-1)/dx^2 = 0.0
                    //     d2S(+1)/dx^2 = 0.0
                    // and see that such spline interpolated f(x) with small error.
                    //
                    double[] x = new double[]{-1.0,-0.5,0.0,+0.5,+1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{+1.0,0.25,0.0,0.25,+1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = 0.25;
                    if( _spoil_scenario==8 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    alglib.spline1dinterpolant s;
                    int natural_bound_type = 2;
                    //
                    // Test exact boundary conditions: build S(x), calculare S(0.25)
                    // (almost same as original function)
                    //
                    alglib.spline1dbuildcubic(x, y, out s);
                    v = alglib.spline1dcalc(s, t);
                    _TestResult = _TestResult && doc_test_real(v, 0.0625, 0.00001);

                    //
                    // Test natural boundary conditions: build S(x), calculare S(0.25)
                    // (small interpolation error)
                    //
                    alglib.spline1dbuildcubic(x, y, 5, natural_bound_type, 0.0, natural_bound_type, 0.0, out s);
                    v = alglib.spline1dcalc(s, t);
                    _TestResult = _TestResult && doc_test_real(v, 0.0580, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_cubic");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline1d_d_monotone
            //      Monotone interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    //
                    // Spline built witn spline1dbuildcubic() can be non-monotone even when
                    // Y-values form monotone sequence. Say, for x=[0,1,2] and y=[0,1,1]
                    // cubic spline will monotonically grow until x=1.5 and then start
                    // decreasing.
                    //
                    // That's why ALGLIB provides special spline construction function
                    // which builds spline which preserves monotonicity of the original
                    // dataset.
                    //
                    // NOTE: in case original dataset is non-monotonic, ALGLIB splits it
                    // into monotone subsequences and builds piecewise monotonic spline.
                    //
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,1,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    alglib.spline1dinterpolant s;

                    // build spline
                    alglib.spline1dbuildmonotone(x, y, out s);

                    // calculate S at x = [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
                    // you may see that spline is really monotonic
                    double v;
                    v = alglib.spline1dcalc(s, -0.5);
                    _TestResult = _TestResult && doc_test_real(v, 0.0000, 0.00005);
                    v = alglib.spline1dcalc(s, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.0000, 0.00005);
                    v = alglib.spline1dcalc(s, +0.5);
                    _TestResult = _TestResult && doc_test_real(v, 0.5000, 0.00005);
                    v = alglib.spline1dcalc(s, 1.0);
                    _TestResult = _TestResult && doc_test_real(v, 1.0000, 0.00005);
                    v = alglib.spline1dcalc(s, 1.5);
                    _TestResult = _TestResult && doc_test_real(v, 1.0000, 0.00005);
                    v = alglib.spline1dcalc(s, 2.0);
                    _TestResult = _TestResult && doc_test_real(v, 1.0000, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_monotone");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline1d_d_griddiff
            //      Differentiation on the grid using cubic splines
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use cubic spline to do grid differentiation, i.e. having
                    // values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
                    // we calculate derivatives of cubic spline at nodes WITHOUT
                    // CONSTRUCTION OF SPLINE OBJECT.
                    //
                    // There are efficient functions spline1dgriddiffcubic() and
                    // spline1dgriddiff2cubic() for such calculations.
                    //
                    // We use default boundary conditions ("parabolically terminated
                    // spline") because cubic spline built with such boundary conditions 
                    // will exactly reproduce any quadratic f(x).
                    //
                    // Actually, we could use natural conditions, but we feel that 
                    // spline which exactly reproduces f() will show us more 
                    // understandable results.
                    //
                    double[] x = new double[]{-1.0,-0.5,0.0,+0.5,+1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{+1.0,0.25,0.0,0.25,+1.0};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] d1;
                    double[] d2;

                    //
                    // We calculate first derivatives: they must be equal to 2*x
                    //
                    alglib.spline1dgriddiffcubic(x, y, out d1);
                    _TestResult = _TestResult && doc_test_real_vector(d1, new double[]{-2.0,-1.0,0.0,+1.0,+2.0}, 0.0001);

                    //
                    // Now test griddiff2, which returns first AND second derivatives.
                    // First derivative is 2*x, second is equal to 2.0
                    //
                    alglib.spline1dgriddiff2cubic(x, y, out d1, out d2);
                    _TestResult = _TestResult && doc_test_real_vector(d1, new double[]{-2.0,-1.0,0.0,+1.0,+2.0}, 0.0001);
                    _TestResult = _TestResult && doc_test_real_vector(d2, new double[]{2.0,2.0,2.0,2.0,2.0}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_griddiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline1d_d_convdiff
            //      Resampling using cubic splines
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use cubic spline to do resampling, i.e. having
                    // values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
                    // we calculate values/derivatives of cubic spline on 
                    // another grid (equidistant with 9 nodes on [-1,+1])
                    // WITHOUT CONSTRUCTION OF SPLINE OBJECT.
                    //
                    // There are efficient functions spline1dconvcubic(),
                    // spline1dconvdiffcubic() and spline1dconvdiff2cubic() 
                    // for such calculations.
                    //
                    // We use default boundary conditions ("parabolically terminated
                    // spline") because cubic spline built with such boundary conditions 
                    // will exactly reproduce any quadratic f(x).
                    //
                    // Actually, we could use natural conditions, but we feel that 
                    // spline which exactly reproduces f() will show us more 
                    // understandable results.
                    //
                    double[] x_old = new double[]{-1.0,-0.5,0.0,+0.5,+1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x_old, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x_old, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x_old, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x_old);
                    double[] y_old = new double[]{+1.0,0.25,0.0,0.25,+1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y_old, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y_old, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y_old, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y_old);
                    double[] x_new = new double[]{-1.00,-0.75,-0.50,-0.25,0.00,+0.25,+0.50,+0.75,+1.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref x_new, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref x_new, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref x_new, (double)System.Double.NegativeInfinity);
                    double[] y_new;
                    double[] d1_new;
                    double[] d2_new;

                    //
                    // First, conversion without differentiation.
                    //
                    //
                    alglib.spline1dconvcubic(x_old, y_old, x_new, out y_new);
                    _TestResult = _TestResult && doc_test_real_vector(y_new, new double[]{1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001);

                    //
                    // Then, conversion with differentiation (first derivatives only)
                    //
                    //
                    alglib.spline1dconvdiffcubic(x_old, y_old, x_new, out y_new, out d1_new);
                    _TestResult = _TestResult && doc_test_real_vector(y_new, new double[]{1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001);
                    _TestResult = _TestResult && doc_test_real_vector(d1_new, new double[]{-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0}, 0.0001);

                    //
                    // Finally, conversion with first and second derivatives
                    //
                    //
                    alglib.spline1dconvdiff2cubic(x_old, y_old, x_new, out y_new, out d1_new, out d2_new);
                    _TestResult = _TestResult && doc_test_real_vector(y_new, new double[]{1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001);
                    _TestResult = _TestResult && doc_test_real_vector(d1_new, new double[]{-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0}, 0.0001);
                    _TestResult = _TestResult && doc_test_real_vector(d2_new, new double[]{2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0}, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_convdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST parametric_rdp
            //      Parametric Ramer-Douglas-Peucker approximation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use RDP algorithm to approximate parametric 2D curve given by
                    // locations in t=0,1,2,3 (see below), which form piecewise linear
                    // trajectory through D-dimensional space (2-dimensional in our example).
                    // 
                    //     |
                    //     |
                    //     -     *     *     X2................X3
                    //     |                .
                    //     |               .
                    //     -     *     *  .  *     *     *     *
                    //     |             .
                    //     |            .
                    //     -     *     X1    *     *     *     *
                    //     |      .....
                    //     |  ....
                    //     X0----|-----|-----|-----|-----|-----|---
                    //
                    int npoints = 4;
                    int ndimensions = 2;
                    double[,] x = new double[,]{{0,0},{2,1},{3,3},{6,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);

                    //
                    // Approximation of parametric curve is performed by another parametric curve
                    // with lesser amount of points. It allows to work with "compressed"
                    // representation, which needs smaller amount of memory. Say, in our example
                    // (we allow points with error smaller than 0.8) approximation will have
                    // just two sequential sections connecting X0 with X2, and X2 with X3.
                    // 
                    //     |
                    //     |
                    //     -     *     *     X2................X3
                    //     |               . 
                    //     |             .  
                    //     -     *     .     *     *     *     *
                    //     |         .    
                    //     |       .     
                    //     -     .     X1    *     *     *     *
                    //     |   .       
                    //     | .    
                    //     X0----|-----|-----|-----|-----|-----|---
                    //
                    //
                    double[,] y;
                    int[] idxy;
                    int nsections;
                    int limitcnt = 0;
                    double limiteps = 0.8;
                    if( _spoil_scenario==5 )
                        limiteps = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==6 )
                        limiteps = (double)System.Double.NegativeInfinity;
                    alglib.parametricrdpfixed(x, npoints, ndimensions, limitcnt, limiteps, out y, out idxy, out nsections);
                    _TestResult = _TestResult && doc_test_int(nsections, 2);
                    _TestResult = _TestResult && doc_test_int_vector(idxy, new int[]{0,2,3});
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "parametric_rdp");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline3d_trilinear
            //      Trilinear spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use trilinear spline to interpolate f(x,y,z)=x+xy+z sampled 
                    // at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
                    //
                    // We store x, y and z-values at local arrays with same names.
                    // Function values are stored in the array F as follows:
                    //     f[0]     (x,y,z) = (0,0,0)
                    //     f[1]     (x,y,z) = (1,0,0)
                    //     f[2]     (x,y,z) = (0,1,0)
                    //     f[3]     (x,y,z) = (1,1,0)
                    //     f[4]     (x,y,z) = (0,0,1)
                    //     f[5]     (x,y,z) = (1,0,1)
                    //     f[6]     (x,y,z) = (0,1,1)
                    //     f[7]     (x,y,z) = (1,1,1)
                    //
                    double[] x = new double[]{0.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] z = new double[]{0.0,1.0};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref z, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref z, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref z, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref z);
                    double[] f = new double[]{0,1,0,2,1,2,1,3};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_deleting_element(ref f);
                    double vx = 0.50;
                    if( _spoil_scenario==16 )
                        vx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        vx = (double)System.Double.NegativeInfinity;
                    double vy = 0.50;
                    if( _spoil_scenario==18 )
                        vy = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        vy = (double)System.Double.NegativeInfinity;
                    double vz = 0.50;
                    if( _spoil_scenario==20 )
                        vz = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==21 )
                        vz = (double)System.Double.NegativeInfinity;
                    double v;
                    alglib.spline3dinterpolant s;

                    // build spline
                    alglib.spline3dbuildtrilinearv(x, 2, y, 2, z, 2, f, 1, out s);

                    // calculate S(0.5,0.5,0.5)
                    v = alglib.spline3dcalc(s, vx, vy, vz);
                    _TestResult = _TestResult && doc_test_real(v, 1.2500, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline3d_trilinear");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline3d_vector
            //      Vector-valued trilinear spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<22; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use trilinear vector-valued spline to interpolate {f0,f1}={x+xy+z,x+xy+yz+z}
                    // sampled at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
                    //
                    // We store x, y and z-values at local arrays with same names.
                    // Function values are stored in the array F as follows:
                    //     f[0]     f0, (x,y,z) = (0,0,0)
                    //     f[1]     f1, (x,y,z) = (0,0,0)
                    //     f[2]     f0, (x,y,z) = (1,0,0)
                    //     f[3]     f1, (x,y,z) = (1,0,0)
                    //     f[4]     f0, (x,y,z) = (0,1,0)
                    //     f[5]     f1, (x,y,z) = (0,1,0)
                    //     f[6]     f0, (x,y,z) = (1,1,0)
                    //     f[7]     f1, (x,y,z) = (1,1,0)
                    //     f[8]     f0, (x,y,z) = (0,0,1)
                    //     f[9]     f1, (x,y,z) = (0,0,1)
                    //     f[10]    f0, (x,y,z) = (1,0,1)
                    //     f[11]    f1, (x,y,z) = (1,0,1)
                    //     f[12]    f0, (x,y,z) = (0,1,1)
                    //     f[13]    f1, (x,y,z) = (0,1,1)
                    //     f[14]    f0, (x,y,z) = (1,1,1)
                    //     f[15]    f1, (x,y,z) = (1,1,1)
                    //
                    double[] x = new double[]{0.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] z = new double[]{0.0,1.0};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref z, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref z, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref z, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref z);
                    double[] f = new double[]{0,0,1,1,0,0,2,2,1,1,2,2,1,2,3,4};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_deleting_element(ref f);
                    double vx = 0.50;
                    if( _spoil_scenario==16 )
                        vx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==17 )
                        vx = (double)System.Double.NegativeInfinity;
                    double vy = 0.50;
                    if( _spoil_scenario==18 )
                        vy = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        vy = (double)System.Double.NegativeInfinity;
                    double vz = 0.50;
                    if( _spoil_scenario==20 )
                        vz = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==21 )
                        vz = (double)System.Double.NegativeInfinity;
                    alglib.spline3dinterpolant s;

                    // build spline
                    alglib.spline3dbuildtrilinearv(x, 2, y, 2, z, 2, f, 2, out s);

                    // calculate S(0.5,0.5,0.5) - we have vector of values instead of single value
                    double[] v;
                    alglib.spline3dcalcv(s, vx, vy, vz, out v);
                    _TestResult = _TestResult && doc_test_real_vector(v, new double[]{1.2500,1.5000}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline3d_vector");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_calcdiff
            //      Interpolation and differentiation using barycentric representation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate polynomial interpolation and differentiation
                    // of y=x^2-x sampled at [0,1,2]. Barycentric representation of polynomial is used.
                    //
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==10 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    double dv;
                    double d2v;
                    alglib.barycentricinterpolant p;

                    // barycentric model is created
                    alglib.polynomialbuild(x, y, out p);

                    // barycentric interpolation is demonstrated
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);

                    // barycentric differentation is demonstrated
                    alglib.barycentricdiff1(p, t, out v, out dv);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dv, -3.0, 0.00005);

                    // second derivatives with barycentric representation
                    alglib.barycentricdiff1(p, t, out v, out dv);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dv, -3.0, 0.00005);
                    alglib.barycentricdiff2(p, t, out v, out dv, out d2v);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dv, -3.0, 0.00005);
                    _TestResult = _TestResult && doc_test_real(d2v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_calcdiff");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_conv
            //      Conversion between power basis and barycentric representation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    //
                    // Here we demonstrate conversion of y=x^2-x
                    // between power basis and barycentric representation.
                    //
                    double[] a = new double[]{0,-1,+1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref a, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref a, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref a, (double)System.Double.NegativeInfinity);
                    double t = 2;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double[] a2;
                    double v;
                    alglib.barycentricinterpolant p;

                    //
                    // a=[0,-1,+1] is decomposition of y=x^2-x in the power basis:
                    //
                    //     y = 0 - 1*x + 1*x^2
                    //
                    // We convert it to the barycentric form.
                    //
                    alglib.polynomialpow2bar(a, out p);

                    // now we have barycentric interpolation; we can use it for interpolation
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.005);

                    // we can also convert back from barycentric representation to power basis
                    alglib.polynomialbar2pow(p, out a2);
                    _TestResult = _TestResult && doc_test_real_vector(a2, new double[]{0,-1,+1}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_conv");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_d_spec
            //      Polynomial interpolation on special grids (equidistant, Chebyshev I/II)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    //
                    // Temporaries:
                    // * values of y=x^2-x sampled at three special grids:
                    //   * equdistant grid spanning [0,2],     x[i] = 2*i/(N-1), i=0..N-1
                    //   * Chebyshev-I grid spanning [-1,+1],  x[i] = 1 + Cos(PI*(2*i+1)/(2*n)), i=0..N-1
                    //   * Chebyshev-II grid spanning [-1,+1], x[i] = 1 + Cos(PI*i/(n-1)), i=0..N-1
                    // * barycentric interpolants for these three grids
                    // * vectors to store coefficients of quadratic representation
                    //
                    double[] y_eqdist = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y_eqdist, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y_eqdist, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y_eqdist, (double)System.Double.NegativeInfinity);
                    double[] y_cheb1 = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref y_cheb1, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y_cheb1, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y_cheb1, (double)System.Double.NegativeInfinity);
                    double[] y_cheb2 = new double[]{0,0,2};
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y_cheb2, (double)System.Double.NaN);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y_cheb2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref y_cheb2, (double)System.Double.NegativeInfinity);
                    alglib.barycentricinterpolant p_eqdist;
                    alglib.barycentricinterpolant p_cheb1;
                    alglib.barycentricinterpolant p_cheb2;
                    double[] a_eqdist;
                    double[] a_cheb1;
                    double[] a_cheb2;

                    //
                    // First, we demonstrate construction of barycentric interpolants on
                    // special grids. We unpack power representation to ensure that
                    // interpolant was built correctly.
                    //
                    // In all three cases we should get same quadratic function.
                    //
                    alglib.polynomialbuildeqdist(0.0, 2.0, y_eqdist, out p_eqdist);
                    alglib.polynomialbar2pow(p_eqdist, out a_eqdist);
                    _TestResult = _TestResult && doc_test_real_vector(a_eqdist, new double[]{0,-1,+1}, 0.00005);

                    alglib.polynomialbuildcheb1(-1, +1, y_cheb1, out p_cheb1);
                    alglib.polynomialbar2pow(p_cheb1, out a_cheb1);
                    _TestResult = _TestResult && doc_test_real_vector(a_cheb1, new double[]{0,-1,+1}, 0.00005);

                    alglib.polynomialbuildcheb2(-1, +1, y_cheb2, out p_cheb2);
                    alglib.polynomialbar2pow(p_cheb2, out a_cheb2);
                    _TestResult = _TestResult && doc_test_real_vector(a_cheb2, new double[]{0,-1,+1}, 0.00005);

                    //
                    // Now we demonstrate polynomial interpolation without construction 
                    // of the barycentricinterpolant structure.
                    //
                    // We calculate interpolant value at x=-2.
                    // In all three cases we should get same f=6
                    //
                    double t = -2;
                    if( _spoil_scenario==9 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalceqdist(0.0, 2.0, y_eqdist, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);

                    v = alglib.polynomialcalccheb1(-1, +1, y_cheb1, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);

                    v = alglib.polynomialcalccheb2(-1, +1, y_cheb2, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_spec");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_1
            //      Polynomial interpolation, full list of parameters.
            //
            System.Console.WriteLine("100/151");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0,1,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==8 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        t = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuild(x, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_2
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildeqdist(0.0, 2.0, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_3
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb1(-1.0, +1.0, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_4
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -2;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb2(a, b, y, 3, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_5
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalceqdist(0.0, 2.0, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_6
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -1;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb1(a, b, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_6");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_7
            //      Polynomial interpolation, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref y);
                    double t = -2;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==5 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==8 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        b = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb2(a, b, y, 3, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_7");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_8
            //      Polynomial interpolation: y=x^2-x, equidistant grid, barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildeqdist(0.0, 2.0, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_8");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_9
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind), barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb1(a, b, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_9");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_10
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind), barycentric form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -2;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.NegativeInfinity;
                    alglib.barycentricinterpolant p;
                    double v;
                    alglib.polynomialbuildcheb2(a, b, y, out p);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_10");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_11
            //      Polynomial interpolation: y=x^2-x, equidistant grid
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalceqdist(0.0, 2.0, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_11");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_12
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{-0.116025,0.000000,1.616025};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -1;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb1(a, b, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_12");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST polint_t_13
            //      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind)
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<11; _spoil_scenario++)
            {
                try
                {
                    double[] y = new double[]{0,0,2};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    double t = -2;
                    if( _spoil_scenario==3 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==4 )
                        t = (double)System.Double.NegativeInfinity;
                    double a = -1;
                    if( _spoil_scenario==5 )
                        a = (double)System.Double.NaN;
                    if( _spoil_scenario==6 )
                        a = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==7 )
                        a = (double)System.Double.NegativeInfinity;
                    double b = +1;
                    if( _spoil_scenario==8 )
                        b = (double)System.Double.NaN;
                    if( _spoil_scenario==9 )
                        b = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==10 )
                        b = (double)System.Double.NegativeInfinity;
                    double v;
                    v = alglib.polynomialcalccheb2(a, b, y, t);
                    _TestResult = _TestResult && doc_test_real(v, 6.0, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_13");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_nlf
            //      Nonlinear fitting using function value only
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<24; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate exponential fitting
                    // by f(x) = exp(-c*x^2)
                    // using function value only.
                    //
                    // Gradient is estimated using combination of numerical differences
                    // and secant updates. diffstep variable stores differentiation step 
                    // (we have to tell algorithm what step to use).
                    //
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    double diffstep = 0.0001;
                    if( _spoil_scenario==16 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==17 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==18 )
                        diffstep = (double)System.Double.NegativeInfinity;

                    //
                    // Fitting without weights
                    //
                    alglib.lsfitcreatef(x, y, c, diffstep, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);

                    //
                    // Fitting with weights
                    // (you can change weights and see how it changes result)
                    //
                    double[] w = new double[]{1,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==19 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==23 )
                        spoil_vector_by_deleting_element(ref w);
                    alglib.lsfitcreatewf(x, y, w, c, diffstep, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlf");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_nlfg
            //      Nonlinear fitting using gradient
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate exponential fitting
                    // by f(x) = exp(-c*x^2)
                    // using function value and gradient (with respect to c).
                    //
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;

                    //
                    // Fitting without weights
                    //
                    alglib.lsfitcreatefg(x, y, c, true, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);

                    //
                    // Fitting with weights
                    // (you can change weights and see how it changes result)
                    //
                    double[] w = new double[]{1,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref w);
                    alglib.lsfitcreatewfg(x, y, w, c, true, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfg");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_nlfgh
            //      Nonlinear fitting using gradient and Hessian
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<21; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate exponential fitting
                    // by f(x) = exp(-c*x^2)
                    // using function value, gradient and Hessian (with respect to c)
                    //
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, (double)System.Double.NegativeInfinity);
                    double epsx = 0.000001;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;

                    //
                    // Fitting without weights
                    //
                    alglib.lsfitcreatefgh(x, y, c, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);

                    //
                    // Fitting with weights
                    // (you can change weights and see how it changes result)
                    //
                    double[] w = new double[]{1,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref w);
                    alglib.lsfitcreatewfgh(x, y, w, c, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, function_cx_1_grad, function_cx_1_hess, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.5}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfgh");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_nlfb
            //      Bound contstrained nonlinear fitting using function value only
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate exponential fitting by
                    //     f(x) = exp(-c*x^2)
                    // subject to bound constraints
                    //     0.0 <= c <= 1.0
                    // using function value only.
                    //
                    // Gradient is estimated using combination of numerical differences
                    // and secant updates. diffstep variable stores differentiation step 
                    // (we have to tell algorithm what step to use).
                    //
                    // Unconstrained solution is c=1.5, but because of constraints we should
                    // get c=1.0 (at the boundary).
                    //
                    double[,] x = new double[,]{{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{0.3};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, (double)System.Double.NegativeInfinity);
                    double[] bndl = new double[]{0.0};
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{1.0};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double epsx = 0.000001;
                    if( _spoil_scenario==17 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==18 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==19 )
                        epsx = (double)System.Double.NegativeInfinity;
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    double diffstep = 0.0001;
                    if( _spoil_scenario==20 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==21 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==22 )
                        diffstep = (double)System.Double.NegativeInfinity;

                    alglib.lsfitcreatef(x, y, c, diffstep, out state);
                    alglib.lsfitsetbc(state, bndl, bndu);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitfit(state, function_cx_1_func, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.0}, 0.05);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfb");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_nlscale
            //      Nonlinear fitting with custom scaling and bound constraints
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<27; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate fitting by
                    //     f(x) = c[0]*(1+c[1]*((x-1999)^c[2]-1))
                    // subject to bound constraints
                    //     -INF  < c[0] < +INF
                    //      -10 <= c[1] <= +10
                    //      0.1 <= c[2] <= 2.0
                    // Data we want to fit are time series of Japan national debt
                    // collected from 2000 to 2008 measured in USD (dollars, not
                    // millions of dollars).
                    //
                    // Our variables are:
                    //     c[0] - debt value at initial moment (2000),
                    //     c[1] - direction coefficient (growth or decrease),
                    //     c[2] - curvature coefficient.
                    // You may see that our variables are badly scaled - first one 
                    // is order of 10^12, and next two are somewhere about 1 in 
                    // magnitude. Such problem is difficult to solve without some
                    // kind of scaling.
                    // That is exactly where lsfitsetscale() function can be used.
                    // We set scale of our variables to [1.0E12, 1, 1], which allows
                    // us to easily solve this problem.
                    //
                    // You can try commenting out lsfitsetscale() call - and you will 
                    // see that algorithm will fail to converge.
                    //
                    double[,] x = new double[,]{{2000},{2001},{2002},{2003},{2004},{2005},{2006},{2007},{2008}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref x);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref x);
                    double[] y = new double[]{4323239600000.0,4560913100000.0,5564091500000.0,6743189300000.0,7284064600000.0,7050129600000.0,7092221500000.0,8483907600000.0,8625804400000.0};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] c = new double[]{1.0e+13,1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref c, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref c, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref c, (double)System.Double.NegativeInfinity);
                    double epsx = 1.0e-5;
                    if( _spoil_scenario==13 )
                        epsx = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        epsx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        epsx = (double)System.Double.NegativeInfinity;
                    double[] bndl = new double[]{-System.Double.PositiveInfinity,-10,0.1};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref bndl, (double)System.Double.NaN);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_deleting_element(ref bndl);
                    double[] bndu = new double[]{System.Double.PositiveInfinity,+10,2.0};
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref bndu, (double)System.Double.NaN);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref bndu);
                    double[] s = new double[]{1.0e+12,1,1};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref s, (double)System.Double.NaN);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref s, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_value(ref s, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==23 )
                        spoil_vector_by_deleting_element(ref s);
                    int maxits = 0;
                    int info;
                    alglib.lsfitstate state;
                    alglib.lsfitreport rep;
                    double diffstep = 1.0e-5;
                    if( _spoil_scenario==24 )
                        diffstep = (double)System.Double.NaN;
                    if( _spoil_scenario==25 )
                        diffstep = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==26 )
                        diffstep = (double)System.Double.NegativeInfinity;

                    alglib.lsfitcreatef(x, y, c, diffstep, out state);
                    alglib.lsfitsetcond(state, epsx, maxits);
                    alglib.lsfitsetbc(state, bndl, bndu);
                    alglib.lsfitsetscale(state, s);
                    alglib.lsfitfit(state, function_debt_func, null, null);
                    alglib.lsfitresults(state, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 2);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{4.142560e+12,0.434240,0.565376}, -0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlscale");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_lin
            //      Unconstrained (general) linear least squares fitting with and without weights
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<13; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
                    //
                    // We have:
                    // * y - vector of experimental data
                    // * fmatrix -  matrix of basis functions calculated at sample points
                    //              Actually, we have only one basis function F0 = exp(0.5*x).
                    //
                    double[,] fmatrix = new double[,]{{0.606531},{0.670320},{0.740818},{0.818731},{0.904837},{1.000000},{1.105171},{1.221403},{1.349859},{1.491825},{1.648721}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.NegativeInfinity);
                    double[] y = new double[]{1.133719,1.306522,1.504604,1.554663,1.884638,2.072436,2.257285,2.534068,2.622017,2.897713,3.219371};
                    if( _spoil_scenario==3 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;

                    //
                    // Linear fitting without weights
                    //
                    alglib.lsfitlinear(y, fmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.98650}, 0.00005);

                    //
                    // Linear fitting with individual weights.
                    // Slightly different result is returned.
                    //
                    double[] w = new double[]{1.414213,1,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_deleting_element(ref w);
                    alglib.lsfitlinearw(y, w, fmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{1.983354}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_lin");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_linc
            //      Constrained (general) linear least squares fitting with and without weights
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate linear fitting by f(x|a,b) = a*x+b
                    // with simple constraint f(0)=0.
                    //
                    // We have:
                    // * y - vector of experimental data
                    // * fmatrix -  matrix of basis functions sampled at [0,1] with step 0.2:
                    //                  [ 1.0   0.0 ]
                    //                  [ 1.0   0.2 ]
                    //                  [ 1.0   0.4 ]
                    //                  [ 1.0   0.6 ]
                    //                  [ 1.0   0.8 ]
                    //                  [ 1.0   1.0 ]
                    //              first column contains value of first basis function (constant term)
                    //              second column contains second basis function (linear term)
                    // * cmatrix -  matrix of linear constraints:
                    //                  [ 1.0  0.0  0.0 ]
                    //              first two columns contain coefficients before basis functions,
                    //              last column contains desired value of their sum.
                    //              So [1,0,0] means "1*constant_term + 0*linear_term = 0" 
                    //
                    double[] y = new double[]{0.072436,0.246944,0.491263,0.522300,0.714064,0.921929};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref y);
                    double[,] fmatrix = new double[,]{{1,0.0},{1,0.2},{1,0.4},{1,0.6},{1,0.8},{1,1.0}};
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_matrix_by_value(ref fmatrix, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_matrix_by_adding_row(ref fmatrix);
                    if( _spoil_scenario==9 )
                        spoil_matrix_by_adding_col(ref fmatrix);
                    if( _spoil_scenario==10 )
                        spoil_matrix_by_deleting_row(ref fmatrix);
                    if( _spoil_scenario==11 )
                        spoil_matrix_by_deleting_col(ref fmatrix);
                    double[,] cmatrix = new double[,]{{1,0,0}};
                    if( _spoil_scenario==12 )
                        spoil_matrix_by_value(ref cmatrix, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_matrix_by_value(ref cmatrix, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_matrix_by_value(ref cmatrix, (double)System.Double.NegativeInfinity);
                    int info;
                    double[] c;
                    alglib.lsfitreport rep;

                    //
                    // Constrained fitting without weights
                    //
                    alglib.lsfitlinearc(y, fmatrix, cmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{0,0.932933}, 0.0005);

                    //
                    // Constrained fitting with individual weights
                    //
                    double[] w = new double[]{1,1.414213,1,1,1,1};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref w);
                    alglib.lsfitlinearwc(y, w, fmatrix, cmatrix, out info, out c, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && doc_test_real_vector(c, new double[]{0,0.938322}, 0.0005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_linc");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_pol
            //      Unconstrained polynomial fitting
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<20; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates polynomial fitting.
                    //
                    // Fitting is done by two (M=2) functions from polynomial basis:
                    //     f0 = 1
                    //     f1 = x
                    // Basically, it just a linear fit; more complex polynomials may be used
                    // (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
                    // us to demonstrate polynomialfit() function in action.
                    //
                    // We have:
                    // * x      set of abscissas
                    // * y      experimental data
                    //
                    // Additionally we demonstrate weighted fitting, where second point has
                    // more weight than other ones.
                    //
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==10 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==11 )
                        t = (double)System.Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;

                    //
                    // Fitting without individual weights
                    //
                    // NOTE: result is returned as barycentricinterpolant structure.
                    //       if you want to get representation in the power basis,
                    //       you can use barycentricbar2pow() function to convert
                    //       from barycentric to power representation (see docs for 
                    //       POLINT subpackage for more info).
                    //
                    alglib.polynomialfit(x, y, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.011, 0.002);

                    //
                    // Fitting with individual weights
                    //
                    // NOTE: slightly different result is returned
                    //
                    double[] w = new double[]{1,1.414213562,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[0];
                    if( _spoil_scenario==17 )
                        spoil_vector_by_adding_element(ref xc);
                    double[] yc = new double[0];
                    if( _spoil_scenario==18 )
                        spoil_vector_by_adding_element(ref yc);
                    int[] dc = new int[0];
                    if( _spoil_scenario==19 )
                        spoil_vector_by_adding_element(ref dc);
                    alglib.polynomialfitwc(x, y, w, xc, yc, dc, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.023, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_pol");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_polc
            //      Constrained polynomial fitting
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<29; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates polynomial fitting.
                    //
                    // Fitting is done by two (M=2) functions from polynomial basis:
                    //     f0 = 1
                    //     f1 = x
                    // with simple constraint on function value
                    //     f(0) = 0
                    // Basically, it just a linear fit; more complex polynomials may be used
                    // (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
                    // us to demonstrate polynomialfit() function in action.
                    //
                    // We have:
                    // * x      set of abscissas
                    // * y      experimental data
                    // * xc     points where constraints are placed
                    // * yc     constraints on derivatives
                    // * dc     derivative indices
                    //          (0 means function itself, 1 means first derivative)
                    //
                    double[] x = new double[]{1.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.9,1.1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1};
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_adding_element(ref w);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{0};
                    if( _spoil_scenario==15 )
                        spoil_vector_by_value(ref xc, (double)System.Double.NaN);
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref xc, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref xc, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_adding_element(ref xc);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref xc);
                    double[] yc = new double[]{0};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_value(ref yc, (double)System.Double.NaN);
                    if( _spoil_scenario==21 )
                        spoil_vector_by_value(ref yc, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==22 )
                        spoil_vector_by_value(ref yc, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==23 )
                        spoil_vector_by_adding_element(ref yc);
                    if( _spoil_scenario==24 )
                        spoil_vector_by_deleting_element(ref yc);
                    int[] dc = new int[]{0};
                    if( _spoil_scenario==25 )
                        spoil_vector_by_adding_element(ref dc);
                    if( _spoil_scenario==26 )
                        spoil_vector_by_deleting_element(ref dc);
                    double t = 2;
                    if( _spoil_scenario==27 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==28 )
                        t = (double)System.Double.NegativeInfinity;
                    int m = 2;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;

                    alglib.polynomialfitwc(x, y, w, xc, yc, dc, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.000, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_polc");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_d_spline
            //      Unconstrained fitting by penalized regression spline
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<19; _spoil_scenario++)
            {
                try
                {
                    //
                    // In this example we demonstrate penalized spline fitting of noisy data
                    //
                    // We have:
                    // * x - abscissas
                    // * y - vector of experimental data, straight line with small noise
                    //
                    double[] x = new double[]{0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_adding_element(ref x);
                    if( _spoil_scenario==4 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.10,0.00,0.30,0.40,0.30,0.40,0.62,0.68,0.75,0.95};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_adding_element(ref y);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_deleting_element(ref y);
                    int info;
                    double v;
                    alglib.spline1dinterpolant s;
                    alglib.spline1dfitreport rep;
                    double rho;

                    //
                    // Fit with VERY small amount of smoothing (rho = -5.0)
                    // and large number of basis functions (M=50).
                    //
                    // With such small regularization penalized spline almost fully reproduces function values
                    //
                    rho = -5.0;
                    if( _spoil_scenario==10 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==11 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==12 )
                        rho = (double)System.Double.NegativeInfinity;
                    alglib.spline1dfitpenalized(x, y, 50, rho, out info, out s, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    v = alglib.spline1dcalc(s, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.10, 0.01);

                    //
                    // Fit with VERY large amount of smoothing (rho = 10.0)
                    // and large number of basis functions (M=50).
                    //
                    // With such regularization our spline should become close to the straight line fit.
                    // We will compare its value in x=1.0 with results obtained from such fit.
                    //
                    rho = +10.0;
                    if( _spoil_scenario==13 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==14 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        rho = (double)System.Double.NegativeInfinity;
                    alglib.spline1dfitpenalized(x, y, 50, rho, out info, out s, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    v = alglib.spline1dcalc(s, 1.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.969, 0.001);

                    //
                    // In real life applications you may need some moderate degree of fitting,
                    // so we try to fit once more with rho=3.0.
                    //
                    rho = +3.0;
                    if( _spoil_scenario==16 )
                        rho = (double)System.Double.NaN;
                    if( _spoil_scenario==17 )
                        rho = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==18 )
                        rho = (double)System.Double.NegativeInfinity;
                    alglib.spline1dfitpenalized(x, y, 50, rho, out info, out s, out rep);
                    _TestResult = _TestResult && doc_test_int(info, 1);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_spline");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_t_polfit_1
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<10; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==8 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==9 )
                        t = (double)System.Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfit(x, y, 11, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.011, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_t_polfit_2
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<14; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1.414213562,1,1,1,1,1,1,1,1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[0];
                    double[] yc = new double[0];
                    int[] dc = new int[0];
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==12 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        t = (double)System.Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, 11, xc, yc, dc, 0, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.023, 0.002);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_t_polfit_3
            //      Polynomial fitting, full list of parameters.
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<23; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{1.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.9,1.1};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] w = new double[]{1,1};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref w, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref w, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref w, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref w);
                    double[] xc = new double[]{0};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref xc, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref xc, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref xc, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_deleting_element(ref xc);
                    double[] yc = new double[]{0};
                    if( _spoil_scenario==16 )
                        spoil_vector_by_value(ref yc, (double)System.Double.NaN);
                    if( _spoil_scenario==17 )
                        spoil_vector_by_value(ref yc, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==18 )
                        spoil_vector_by_value(ref yc, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==19 )
                        spoil_vector_by_deleting_element(ref yc);
                    int[] dc = new int[]{0};
                    if( _spoil_scenario==20 )
                        spoil_vector_by_deleting_element(ref dc);
                    int m = 2;
                    double t = 2;
                    if( _spoil_scenario==21 )
                        t = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==22 )
                        t = (double)System.Double.NegativeInfinity;
                    int info;
                    alglib.barycentricinterpolant p;
                    alglib.polynomialfitreport rep;
                    double v;
                    alglib.polynomialfitwc(x, y, w, 2, xc, yc, dc, 1, m, out info, out p, out rep);
                    v = alglib.barycentriccalc(p, t);
                    _TestResult = _TestResult && doc_test_real(v, 2.000, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_t_4pl
            //      4-parameter logistic fitting
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<8; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{1,2,3,4,5,6,7,8};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.06313223,0.44552624,0.61838364,0.71385108,0.77345838,0.81383140,0.84280033,0.86449822};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int n = 8;
                    double a;
                    double b;
                    double c;
                    double d;
                    alglib.lsfitreport rep;

                    //
                    // Test logisticfit4() on carefully designed data with a priori known answer.
                    //
                    alglib.logisticfit4(x, y, n, out a, out b, out c, out d, out rep);
                    _TestResult = _TestResult && doc_test_real(a, -1.000, 0.01);
                    _TestResult = _TestResult && doc_test_real(b, 1.200, 0.01);
                    _TestResult = _TestResult && doc_test_real(c, 0.900, 0.01);
                    _TestResult = _TestResult && doc_test_real(d, 1.000, 0.01);

                    //
                    // Evaluate model at point x=0.5
                    //
                    double v;
                    v = alglib.logisticcalc4(0.5, a, b, c, d);
                    _TestResult = _TestResult && doc_test_real(v, -0.33874308, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_4pl");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lsfit_t_5pl
            //      5-parameter logistic fitting
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<8; _spoil_scenario++)
            {
                try
                {
                    double[] x = new double[]{1,2,3,4,5,6,7,8};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.1949776139,0.5710060208,0.726002637,0.8060434158,0.8534547965,0.8842071579,0.9054773317,0.9209088299};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    int n = 8;
                    double a;
                    double b;
                    double c;
                    double d;
                    double g;
                    alglib.lsfitreport rep;

                    //
                    // Test logisticfit5() on carefully designed data with a priori known answer.
                    //
                    alglib.logisticfit5(x, y, n, out a, out b, out c, out d, out g, out rep);
                    _TestResult = _TestResult && doc_test_real(a, -1.000, 0.01);
                    _TestResult = _TestResult && doc_test_real(b, 1.200, 0.01);
                    _TestResult = _TestResult && doc_test_real(c, 0.900, 0.01);
                    _TestResult = _TestResult && doc_test_real(d, 1.000, 0.01);
                    _TestResult = _TestResult && doc_test_real(g, 1.200, 0.01);

                    //
                    // Evaluate model at point x=0.5
                    //
                    double v;
                    v = alglib.logisticcalc5(0.5, a, b, c, d, g);
                    _TestResult = _TestResult && doc_test_real(v, -0.2354656824, 0.001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_5pl");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_bilinear
            //      Bilinear spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled 
                    // at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
                    //
                    double[] x = new double[]{0.0,0.5,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] f = new double[]{0.00,0.25,1.00,2.00,2.25,3.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref f);
                    double vx = 0.25;
                    if( _spoil_scenario==12 )
                        vx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        vx = (double)System.Double.NegativeInfinity;
                    double vy = 0.50;
                    if( _spoil_scenario==14 )
                        vy = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        vy = (double)System.Double.NegativeInfinity;
                    double v;
                    alglib.spline2dinterpolant s;

                    // build spline
                    alglib.spline2dbuildbilinearv(x, 3, y, 2, f, 1, out s);

                    // calculate S(0.25,0.50)
                    v = alglib.spline2dcalc(s, vx, vy);
                    _TestResult = _TestResult && doc_test_real(v, 1.1250, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_bilinear");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_bicubic
            //      Bilinear spline interpolation
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled 
                    // at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
                    //
                    double[] x = new double[]{0.0,0.5,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] f = new double[]{0.00,0.25,1.00,2.00,2.25,3.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref f);
                    double vx = 0.25;
                    if( _spoil_scenario==12 )
                        vx = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==13 )
                        vx = (double)System.Double.NegativeInfinity;
                    double vy = 0.50;
                    if( _spoil_scenario==14 )
                        vy = (double)System.Double.PositiveInfinity;
                    if( _spoil_scenario==15 )
                        vy = (double)System.Double.NegativeInfinity;
                    double v;
                    double dx;
                    double dy;
                    double dxy;
                    alglib.spline2dinterpolant s;

                    // build spline
                    alglib.spline2dbuildbicubicv(x, 3, y, 2, f, 1, out s);

                    // calculate S(0.25,0.50)
                    v = alglib.spline2dcalc(s, vx, vy);
                    _TestResult = _TestResult && doc_test_real(v, 1.0625, 0.00005);

                    // calculate derivatives
                    alglib.spline2ddiff(s, vx, vy, out v, out dx, out dy, out dxy);
                    _TestResult = _TestResult && doc_test_real(v, 1.0625, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dx, 0.5000, 0.00005);
                    _TestResult = _TestResult && doc_test_real(dy, 2.0000, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_bicubic");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_fit_blocklls
            //      Fitting bicubic spline to irregular data
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    //
                    // We use bicubic spline to reproduce f(x,y)=1/(1+x^2+2*y^2) sampled
                    // at irregular points (x,y) from [-1,+1]*[-1,+1]
                    //
                    // We have 5 such points, located approximately at corners of the area
                    // and its center -  but not exactly at the grid. Thus, we have to FIT
                    // the spline, i.e. to solve least squares problem
                    //
                    double[,] xy = new double[,]{{-0.987,-0.902,0.359},{0.948,-0.992,0.347},{-1.000,1.000,0.333},{1.000,0.973,0.339},{0.017,0.180,0.968}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref xy);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref xy);

                    //
                    // First step is to create spline2dbuilder object and set its properties:
                    // * d=1 means that we create vector-valued spline with 1 component
                    // * we specify dataset xy
                    // * we rely on automatic selection of interpolation area
                    // * we tell builder that we want to use 5x5 grid for an underlying spline
                    // * we choose least squares solver named BlockLLS and configure it by
                    //   telling that we want to apply zero nonlinearity penalty.
                    //
                    // NOTE: you can specify non-zero lambdav if you want to make your spline
                    //       more "rigid", i.e. to penalize nonlinearity.
                    //
                    // NOTE: ALGLIB has two solvers which fit bicubic splines to irregular data,
                    //       one of them is BlockLLS and another one is FastDDM. Former is
                    //       intended for moderately sized grids (up to 512x512 nodes, although
                    //       it may take up to few minutes); it is the most easy to use and
                    //       control spline fitting function in the library. Latter, FastDDM,
                    //       is intended for efficient solution of large-scale problems
                    //       (up to 100.000.000 nodes). Both solvers can be parallelized, but
                    //       FastDDM is much more efficient. See comments for more information.
                    //
                    alglib.spline2dbuilder builder;
                    int d = 1;
                    double lambdav = 0.000;
                    alglib.spline2dbuildercreate(d, out builder);
                    alglib.spline2dbuildersetpoints(builder, xy, 5);
                    alglib.spline2dbuildersetgrid(builder, 5, 5);
                    alglib.spline2dbuildersetalgoblocklls(builder, lambdav);

                    //
                    // Now we are ready to fit and evaluate our results
                    //
                    alglib.spline2dinterpolant s;
                    alglib.spline2dfitreport rep;
                    alglib.spline2dfit(builder, out s, out rep);

                    // evaluate results - function value at the grid is reproduced exactly
                    double v;
                    v = alglib.spline2dcalc(s, -1, 1);
                    _TestResult = _TestResult && doc_test_real(v, 0.333000, 0.005);

                    // check maximum error - it must be nearly zero
                    _TestResult = _TestResult && doc_test_real(rep.maxerror, 0.000, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_fit_blocklls");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_unpack
            //      Unpacking bilinear spline
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // We build bilinear spline for f(x,y)=x+2*y+3*xy for (x,y) in [0,1].
                    // Then we demonstrate how to unpack it.
                    //
                    double[] x = new double[]{0.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] f = new double[]{0.00,1.00,2.00,6.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref f);
                    double[,] c;
                    int m;
                    int n;
                    int d;
                    alglib.spline2dinterpolant s;

                    // build spline
                    alglib.spline2dbuildbilinearv(x, 2, y, 2, f, 1, out s);

                    // unpack and test
                    alglib.spline2dunpackv(s, out m, out n, out d, out c);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0,1,0,1,0,2,0,0,1,3,0,0,0,0,0,0,0,0,0,0}}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_unpack");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_copytrans
            //      Copy and transform
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<16; _spoil_scenario++)
            {
                try
                {
                    //
                    // We build bilinear spline for f(x,y)=x+2*y for (x,y) in [0,1].
                    // Then we apply several transformations to this spline.
                    //
                    double[] x = new double[]{0.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] f = new double[]{0.00,1.00,2.00,3.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref f);
                    alglib.spline2dinterpolant s;
                    alglib.spline2dinterpolant snew;
                    double v;
                    alglib.spline2dbuildbilinearv(x, 2, y, 2, f, 1, out s);

                    // copy spline, apply transformation x:=2*xnew, y:=4*ynew
                    // evaluate at (xnew,ynew) = (0.25,0.25) - should be same as (x,y)=(0.5,1.0)
                    alglib.spline2dcopy(s, out snew);
                    alglib.spline2dlintransxy(snew, 2.0, 0.0, 4.0, 0.0);
                    v = alglib.spline2dcalc(snew, 0.25, 0.25);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.00005);

                    // copy spline, apply transformation SNew:=2*S+3
                    alglib.spline2dcopy(s, out snew);
                    alglib.spline2dlintransf(snew, 2.0, 3.0);
                    v = alglib.spline2dcalc(snew, 0.5, 1.0);
                    _TestResult = _TestResult && doc_test_real(v, 8.000, 0.00005);

                    //
                    // Same example, but for vector spline (f0,f1) = {x+2*y, 2*x+y}
                    //
                    double[] f2 = new double[]{0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00};
                    if( _spoil_scenario==12 )
                        spoil_vector_by_value(ref f2, (double)System.Double.NaN);
                    if( _spoil_scenario==13 )
                        spoil_vector_by_value(ref f2, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==14 )
                        spoil_vector_by_value(ref f2, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==15 )
                        spoil_vector_by_deleting_element(ref f2);
                    double[] vr;
                    alglib.spline2dbuildbilinearv(x, 2, y, 2, f2, 2, out s);

                    // copy spline, apply transformation x:=2*xnew, y:=4*ynew
                    alglib.spline2dcopy(s, out snew);
                    alglib.spline2dlintransxy(snew, 2.0, 0.0, 4.0, 0.0);
                    alglib.spline2dcalcv(snew, 0.25, 0.25, out vr);
                    _TestResult = _TestResult && doc_test_real_vector(vr, new double[]{2.500,2.000}, 0.00005);

                    // copy spline, apply transformation SNew:=2*S+3
                    alglib.spline2dcopy(s, out snew);
                    alglib.spline2dlintransf(snew, 2.0, 3.0);
                    alglib.spline2dcalcv(snew, 0.5, 1.0, out vr);
                    _TestResult = _TestResult && doc_test_real_vector(vr, new double[]{8.000,7.000}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_copytrans");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST spline2d_vector
            //      Copy and transform
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<12; _spoil_scenario++)
            {
                try
                {
                    //
                    // We build bilinear vector-valued spline (f0,f1) = {x+2*y, 2*x+y}
                    // Spline is built using function values at 2x2 grid: (x,y)=[0,1]*[0,1]
                    // Then we perform evaluation at (x,y)=(0.1,0.3)
                    //
                    double[] x = new double[]{0.0,1.0};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref x);
                    double[] y = new double[]{0.0,1.0};
                    if( _spoil_scenario==4 )
                        spoil_vector_by_value(ref y, (double)System.Double.NaN);
                    if( _spoil_scenario==5 )
                        spoil_vector_by_value(ref y, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==6 )
                        spoil_vector_by_value(ref y, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==7 )
                        spoil_vector_by_deleting_element(ref y);
                    double[] f = new double[]{0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00};
                    if( _spoil_scenario==8 )
                        spoil_vector_by_value(ref f, (double)System.Double.NaN);
                    if( _spoil_scenario==9 )
                        spoil_vector_by_value(ref f, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==10 )
                        spoil_vector_by_value(ref f, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==11 )
                        spoil_vector_by_deleting_element(ref f);
                    alglib.spline2dinterpolant s;
                    double[] vr;
                    alglib.spline2dbuildbilinearv(x, 2, y, 2, f, 2, out s);
                    alglib.spline2dcalcv(s, 0.1, 0.3, out vr);
                    _TestResult = _TestResult && doc_test_real_vector(vr, new double[]{0.700,0.500}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_vector");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST rbf_d_hrbf
            //      Simple model built with HRBF algorithm
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example illustrates basic concepts of the RBF models: creation, modification,
                    // evaluation.
                    // 
                    // Suppose that we have set of 2-dimensional points with associated
                    // scalar function values, and we want to build a RBF model using
                    // our data.
                    // 
                    // NOTE: we can work with 3D models too :)
                    // 
                    // Typical sequence of steps is given below:
                    // 1. we create RBF model object
                    // 2. we attach our dataset to the RBF model and tune algorithm settings
                    // 3. we rebuild RBF model using QNN algorithm on new data
                    // 4. we use RBF model (evaluate, serialize, etc.)
                    //
                    double v;

                    //
                    // Step 1: RBF model creation.
                    //
                    // We have to specify dimensionality of the space (2 or 3) and
                    // dimensionality of the function (scalar or vector).
                    //
                    // New model is empty - it can be evaluated,
                    // but we just get zero value at any point.
                    //
                    alglib.rbfmodel model;
                    alglib.rbfcreate(2, 1, out model);

                    v = alglib.rbfcalc2(model, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.000, 0.005);

                    //
                    // Step 2: we add dataset.
                    //
                    // XY contains two points - x0=(-1,0) and x1=(+1,0) -
                    // and two function values f(x0)=2, f(x1)=3.
                    //
                    // We added points, but model was not rebuild yet.
                    // If we call rbfcalc2(), we still will get 0.0 as result.
                    //
                    double[,] xy = new double[,]{{-1,0,2},{+1,0,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.rbfsetpoints(model, xy);

                    v = alglib.rbfcalc2(model, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.000, 0.005);

                    //
                    // Step 3: rebuild model
                    //
                    // After we've configured model, we should rebuild it -
                    // it will change coefficients stored internally in the
                    // rbfmodel structure.
                    //
                    // We use hierarchical RBF algorithm with following parameters:
                    // * RBase - set to 1.0
                    // * NLayers - three layers are used (although such simple problem
                    //   does not need more than 1 layer)
                    // * LambdaReg - is set to zero value, no smoothing is required
                    //
                    alglib.rbfreport rep;
                    alglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0);
                    alglib.rbfbuildmodel(model, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);

                    //
                    // Step 4: model was built
                    //
                    // After call of rbfbuildmodel(), rbfcalc2() will return
                    // value of the new model.
                    //
                    v = alglib.rbfcalc2(model, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_hrbf");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST rbf_d_vector
            //      Working with vector functions
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    //
                    // Suppose that we have set of 2-dimensional points with associated VECTOR
                    // function values, and we want to build a RBF model using our data.
                    // 
                    // Typical sequence of steps is given below:
                    // 1. we create RBF model object
                    // 2. we attach our dataset to the RBF model and tune algorithm settings
                    // 3. we rebuild RBF model using new data
                    // 4. we use RBF model (evaluate, serialize, etc.)
                    //
                    double[] x;
                    double[] y;

                    //
                    // Step 1: RBF model creation.
                    //
                    // We have to specify dimensionality of the space (equal to 2) and
                    // dimensionality of the function (2-dimensional vector function).
                    //
                    // New model is empty - it can be evaluated,
                    // but we just get zero value at any point.
                    //
                    alglib.rbfmodel model;
                    alglib.rbfcreate(2, 2, out model);

                    x = new double[]{+1,+1};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref x, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref x, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref x, (double)System.Double.NegativeInfinity);
                    alglib.rbfcalc(model, x, out y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,0.000}, 0.005);

                    //
                    // Step 2: we add dataset.
                    //
                    // XY arrays containt four points:
                    // * (x0,y0) = (+1,+1), f(x0,y0)=(0,-1)
                    // * (x1,y1) = (+1,-1), f(x1,y1)=(-1,0)
                    // * (x2,y2) = (-1,-1), f(x2,y2)=(0,+1)
                    // * (x3,y3) = (-1,+1), f(x3,y3)=(+1,0)
                    //
                    double[,] xy = new double[,]{{+1,+1,0,-1},{+1,-1,-1,0},{-1,-1,0,+1},{-1,+1,+1,0}};
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.rbfsetpoints(model, xy);

                    // We added points, but model was not rebuild yet.
                    // If we call rbfcalc(), we still will get 0.0 as result.
                    alglib.rbfcalc(model, x, out y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,0.000}, 0.005);

                    //
                    // Step 3: rebuild model
                    //
                    // We use hierarchical RBF algorithm with following parameters:
                    // * RBase - set to 1.0
                    // * NLayers - three layers are used (although such simple problem
                    //   does not need more than 1 layer)
                    // * LambdaReg - is set to zero value, no smoothing is required
                    //
                    // After we've configured model, we should rebuild it -
                    // it will change coefficients stored internally in the
                    // rbfmodel structure.
                    //
                    alglib.rbfreport rep;
                    alglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0);
                    alglib.rbfbuildmodel(model, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);

                    //
                    // Step 4: model was built
                    //
                    // After call of rbfbuildmodel(), rbfcalc() will return
                    // value of the new model.
                    //
                    alglib.rbfcalc(model, x, out y);
                    _TestResult = _TestResult && doc_test_real_vector(y, new double[]{0.000,-1.000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_vector");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST rbf_d_polterm
            //      RBF models - working with polynomial term
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example show how to work with polynomial term
                    // 
                    // Suppose that we have set of 2-dimensional points with associated
                    // scalar function values, and we want to build a RBF model using
                    // our data.
                    //
                    // We use hierarchical RBF algorithm with following parameters:
                    // * RBase - set to 1.0
                    // * NLayers - three layers are used (although such simple problem
                    //   does not need more than 1 layer)
                    // * LambdaReg - is set to zero value, no smoothing is required
                    //
                    double v;
                    alglib.rbfmodel model;
                    double[,] xy = new double[,]{{-1,0,2},{+1,0,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.rbfreport rep;

                    alglib.rbfcreate(2, 1, out model);
                    alglib.rbfsetpoints(model, xy);
                    alglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0);

                    //
                    // By default, RBF model uses linear term. It means that model
                    // looks like
                    //     f(x,y) = SUM(RBF[i]) + a*x + b*y + c
                    // where RBF[i] is I-th radial basis function and a*x+by+c is a
                    // linear term. Having linear terms in a model gives us:
                    // (1) improved extrapolation properties
                    // (2) linearity of the model when data can be perfectly fitted
                    //     by the linear function
                    // (3) linear asymptotic behavior
                    //
                    // Our simple dataset can be modelled by the linear function
                    //     f(x,y) = 0.5*x + 2.5
                    // and rbfbuildmodel() with default settings should preserve this
                    // linearity.
                    //
                    int nx;
                    int ny;
                    int nc;
                    int modelversion;
                    double[,] xwr;
                    double[,] c;
                    alglib.rbfbuildmodel(model, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);
                    alglib.rbfunpack(model, out nx, out ny, out xwr, out nc, out c, out modelversion);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0.500,0.000,2.500}}, 0.005);

                    // asymptotic behavior of our function is linear
                    v = alglib.rbfcalc2(model, 1000.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 502.50, 0.05);

                    //
                    // Instead of linear term we can use constant term. In this case
                    // we will get model which has form
                    //     f(x,y) = SUM(RBF[i]) + c
                    // where RBF[i] is I-th radial basis function and c is a constant,
                    // which is equal to the average function value on the dataset.
                    //
                    // Because we've already attached dataset to the model the only
                    // thing we have to do is to call rbfsetconstterm() and then
                    // rebuild model with rbfbuildmodel().
                    //
                    alglib.rbfsetconstterm(model);
                    alglib.rbfbuildmodel(model, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);
                    alglib.rbfunpack(model, out nx, out ny, out xwr, out nc, out c, out modelversion);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0.000,0.000,2.500}}, 0.005);

                    // asymptotic behavior of our function is constant
                    v = alglib.rbfcalc2(model, 1000.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.005);

                    //
                    // Finally, we can use zero term. Just plain RBF without polynomial
                    // part:
                    //     f(x,y) = SUM(RBF[i])
                    // where RBF[i] is I-th radial basis function.
                    //
                    alglib.rbfsetzeroterm(model);
                    alglib.rbfbuildmodel(model, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);
                    alglib.rbfunpack(model, out nx, out ny, out xwr, out nc, out c, out modelversion);
                    _TestResult = _TestResult && doc_test_real_matrix(c, new double[,]{{0.000,0.000,0.000}}, 0.005);

                    // asymptotic behavior of our function is just zero constant
                    v = alglib.rbfcalc2(model, 1000.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.000, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_polterm");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST rbf_d_serialize
            //      Serialization/unserialization
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<3; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example show how to serialize and unserialize RBF model
                    // 
                    // Suppose that we have set of 2-dimensional points with associated
                    // scalar function values, and we want to build a RBF model using
                    // our data. Then we want to serialize it to string and to unserialize
                    // from string, loading to another instance of RBF model.
                    //
                    // Here we assume that you already know how to create RBF models.
                    //
                    string s;
                    double v;
                    alglib.rbfmodel model0;
                    alglib.rbfmodel model1;
                    double[,] xy = new double[,]{{-1,0,2},{+1,0,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref xy, (double)System.Double.NegativeInfinity);
                    alglib.rbfreport rep;

                    // model initialization
                    alglib.rbfcreate(2, 1, out model0);
                    alglib.rbfsetpoints(model0, xy);
                    alglib.rbfsetalgohierarchical(model0, 1.0, 3, 0.0);
                    alglib.rbfbuildmodel(model0, out rep);
                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);

                    //
                    // Serialization - it looks easy,
                    // but you should carefully read next section.
                    //
                    alglib.rbfserialize(model0, out s);
                    alglib.rbfunserialize(s, out model1);

                    // both models return same value
                    v = alglib.rbfcalc2(model0, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.005);
                    v = alglib.rbfcalc2(model1, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.005);

                    //
                    // Previous section shows that model state is saved/restored during
                    // serialization. However, some properties are NOT serialized.
                    //
                    // Serialization saves/restores RBF model, but it does NOT saves/restores
                    // settings which were used to build current model. In particular, dataset
                    // which was used to build model, is not preserved.
                    //
                    // What does it mean in for us?
                    //
                    // Do you remember this sequence: rbfcreate-rbfsetpoints-rbfbuildmodel?
                    // First step creates model, second step adds dataset and tunes model
                    // settings, third step builds model using current dataset and model
                    // construction settings.
                    //
                    // If you call rbfbuildmodel() without calling rbfsetpoints() first, you
                    // will get empty (zero) RBF model. In our example, model0 contains
                    // dataset which was added by rbfsetpoints() call. However, model1 does
                    // NOT contain dataset - because dataset is NOT serialized.
                    //
                    // This, if we call rbfbuildmodel(model0,rep), we will get same model,
                    // which returns 2.5 at (x,y)=(0,0). However, after same call model1 will
                    // return zero - because it contains RBF model (coefficients), but does NOT
                    // contain dataset which was used to build this model.
                    //
                    // Basically, it means that:
                    // * serialization of the RBF model preserves anything related to the model
                    //   EVALUATION
                    // * but it does NOT creates perfect copy of the original object.
                    //
                    alglib.rbfbuildmodel(model0, out rep);
                    v = alglib.rbfcalc2(model0, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 2.500, 0.005);

                    alglib.rbfbuildmodel(model1, out rep);
                    v = alglib.rbfcalc2(model1, 0.0, 0.0);
                    _TestResult = _TestResult && doc_test_real(v, 0.000, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_serialize");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_1
            //      Determinant calculation, real matrix, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    double[,] b = new double[,]{{1,2},{2,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    double a;
                    a = alglib.rmatrixdet(b);
                    _TestResult = _TestResult && doc_test_real(a, -3, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_2
            //      Determinant calculation, real matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double[,] b = new double[,]{{5,4},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    double a;
                    a = alglib.rmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_real(a, 9, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_3
            //      Determinant calculation, complex matrix, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(1,+1),2},{2,new alglib.complex(1,-1)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    alglib.complex a;
                    a = alglib.cmatrixdet(b);
                    _TestResult = _TestResult && doc_test_complex(a, -2, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_4
            //      Determinant calculation, complex matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(0,5),4},{new alglib.complex(0,4),5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_complex(a, new alglib.complex(0,9), 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_d_5
            //      Determinant calculation, complex matrix with zero imaginary part, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<7; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{9,1},{2,1}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b);
                    _TestResult = _TestResult && doc_test_complex(a, 7, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_0
            //      Determinant calculation, real matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{3,4},{-4,3}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.rmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_real(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_0");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_1
            //      Determinant calculation, real matrix, LU, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{1,2},{2,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{1,1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_adding_element(ref p);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.rmatrixludet(b, p);
                    _TestResult = _TestResult && doc_test_real(a, -5, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_2
            //      Determinant calculation, real matrix, LU, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    double a;
                    double[,] b = new double[,]{{5,4},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{0,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.rmatrixludet(b, p, 2);
                    _TestResult = _TestResult && doc_test_real(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_2");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_3
            //      Determinant calculation, complex matrix, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<5; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{new alglib.complex(0,5),4},{-4,new alglib.complex(0,5)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    a = alglib.cmatrixdet(b, 2);
                    _TestResult = _TestResult && doc_test_complex(a, -9, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_3");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_4
            //      Determinant calculation, complex matrix, LU, short form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<9; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{1,2},{2,new alglib.complex(0,5)}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_adding_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_adding_col(ref b);
                    if( _spoil_scenario==5 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==6 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{1,1};
                    if( _spoil_scenario==7 )
                        spoil_vector_by_adding_element(ref p);
                    if( _spoil_scenario==8 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.cmatrixludet(b, p);
                    _TestResult = _TestResult && doc_test_complex(a, new alglib.complex(0,-5), 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_4");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST matdet_t_5
            //      Determinant calculation, complex matrix, LU, full form
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<6; _spoil_scenario++)
            {
                try
                {
                    alglib.complex a;
                    alglib.complex[,] b = new alglib.complex[,]{{5,new alglib.complex(0,4)},{4,5}};
                    if( _spoil_scenario==0 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_matrix_by_value(ref b, (alglib.complex)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_matrix_by_deleting_row(ref b);
                    if( _spoil_scenario==4 )
                        spoil_matrix_by_deleting_col(ref b);
                    int[] p = new int[]{0,1};
                    if( _spoil_scenario==5 )
                        spoil_vector_by_deleting_element(ref p);
                    a = alglib.cmatrixludet(b, p, 2);
                    _TestResult = _TestResult && doc_test_complex(a, 25, 0.0001);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_5");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST solvesks_d_1
            //      Solving positive definite sparse system using Skyline (SKS) solver
            //
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example demonstrates creation/initialization of the sparse matrix
                    // in the SKS (Skyline) storage format and solution using SKS-based direct
                    // solver.
                    //
                    // First, we have to create matrix and initialize it. Matrix is created
                    // in the SKS format, using fixed bandwidth initialization function.
                    // Several points should be noted:
                    //
                    // 1. SKS sparse storage format also allows variable bandwidth matrices;
                    //    we just do not want to overcomplicate this example.
                    //
                    // 2. SKS format requires you to specify matrix geometry prior to
                    //    initialization of its elements with sparseset(). If you specified
                    //    bandwidth=1, you can not change your mind afterwards and call
                    //    sparseset() for non-existent elements.
                    // 
                    // 3. Because SKS solver need just one triangle of SPD matrix, we can
                    //    omit initialization of the lower triangle of our matrix.
                    //
                    int n = 4;
                    int bandwidth = 1;
                    alglib.sparsematrix s;
                    alglib.sparsecreatesksband(n, n, bandwidth, out s);
                    alglib.sparseset(s, 0, 0, 2.0);
                    alglib.sparseset(s, 0, 1, 1.0);
                    alglib.sparseset(s, 1, 1, 3.0);
                    alglib.sparseset(s, 1, 2, 1.0);
                    alglib.sparseset(s, 2, 2, 3.0);
                    alglib.sparseset(s, 2, 3, 1.0);
                    alglib.sparseset(s, 3, 3, 2.0);

                    //
                    // Now we have symmetric positive definite 4x4 system width bandwidth=1:
                    //
                    //     [ 2 1     ]   [ x0]]   [  4 ]
                    //     [ 1 3 1   ]   [ x1 ]   [ 10 ]
                    //     [   1 3 1 ] * [ x2 ] = [ 15 ]
                    //     [     1 2 ]   [ x3 ]   [ 11 ]
                    //
                    // After successful creation we can call SKS solver.
                    //
                    double[] b = new double[]{4,10,15,11};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref b);
                    alglib.sparsesolverreport rep;
                    double[] x;
                    bool isuppertriangle = true;
                    alglib.sparsesolvesks(s, n, isuppertriangle, b, out rep, out x);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{1.0000,2.0000,3.0000,4.0000}, 0.00005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "solvesks_d_1");
            _TotalResult = _TotalResult && _TestResult;


            //
            // TEST lincg_d_1
            //      Solution of sparse linear systems with CG
            //
            System.Console.WriteLine("150/151");
            _TestResult = true;
            for(_spoil_scenario=-1; _spoil_scenario<4; _spoil_scenario++)
            {
                try
                {
                    //
                    // This example illustrates solution of sparse linear systems with
                    // conjugate gradient method.
                    // 
                    // Suppose that we have linear system A*x=b with sparse symmetric
                    // positive definite A (represented by sparsematrix object)
                    //         [ 5 1       ]
                    //         [ 1 7 2     ]
                    //     A = [   2 8 1   ]
                    //         [     1 4 1 ]
                    //         [       1 4 ]
                    // and right part b
                    //     [  7 ]
                    //     [ 17 ]
                    // b = [ 14 ]
                    //     [ 10 ]
                    //     [  6 ]
                    // and we want to solve this system using sparse linear CG. In order
                    // to do so, we have to create left part (sparsematrix object) and
                    // right part (dense array).
                    //
                    // Initially, sparse matrix is created in the Hash-Table format,
                    // which allows easy initialization, but do not allow matrix to be
                    // used in the linear solvers. So after construction you should convert
                    // sparse matrix to CRS format (one suited for linear operations).
                    //
                    // It is important to note that in our example we initialize full
                    // matrix A, both lower and upper triangles. However, it is symmetric
                    // and sparse solver needs just one half of the matrix. So you may
                    // save about half of the space by filling only one of the triangles.
                    //
                    alglib.sparsematrix a;
                    alglib.sparsecreate(5, 5, out a);
                    alglib.sparseset(a, 0, 0, 5.0);
                    alglib.sparseset(a, 0, 1, 1.0);
                    alglib.sparseset(a, 1, 0, 1.0);
                    alglib.sparseset(a, 1, 1, 7.0);
                    alglib.sparseset(a, 1, 2, 2.0);
                    alglib.sparseset(a, 2, 1, 2.0);
                    alglib.sparseset(a, 2, 2, 8.0);
                    alglib.sparseset(a, 2, 3, 1.0);
                    alglib.sparseset(a, 3, 2, 1.0);
                    alglib.sparseset(a, 3, 3, 4.0);
                    alglib.sparseset(a, 3, 4, 1.0);
                    alglib.sparseset(a, 4, 3, 1.0);
                    alglib.sparseset(a, 4, 4, 4.0);

                    //
                    // Now our matrix is fully initialized, but we have to do one more
                    // step - convert it from Hash-Table format to CRS format (see
                    // documentation on sparse matrices for more information about these
                    // formats).
                    //
                    // If you omit this call, ALGLIB will generate exception on the first
                    // attempt to use A in linear operations. 
                    //
                    alglib.sparseconverttocrs(a);

                    //
                    // Initialization of the right part
                    //
                    double[] b = new double[]{7,17,14,10,6};
                    if( _spoil_scenario==0 )
                        spoil_vector_by_value(ref b, (double)System.Double.NaN);
                    if( _spoil_scenario==1 )
                        spoil_vector_by_value(ref b, (double)System.Double.PositiveInfinity);
                    if( _spoil_scenario==2 )
                        spoil_vector_by_value(ref b, (double)System.Double.NegativeInfinity);
                    if( _spoil_scenario==3 )
                        spoil_vector_by_deleting_element(ref b);

                    //
                    // Now we have to create linear solver object and to use it for the
                    // solution of the linear system.
                    //
                    // NOTE: lincgsolvesparse() accepts additional parameter which tells
                    //       what triangle of the symmetric matrix should be used - upper
                    //       or lower. Because we've filled both parts of the matrix, we
                    //       can use any part - upper or lower.
                    //
                    alglib.lincgstate s;
                    alglib.lincgreport rep;
                    double[] x;
                    alglib.lincgcreate(5, out s);
                    alglib.lincgsolvesparse(s, a, true, b);
                    alglib.lincgresults(s, out x, out rep);

                    _TestResult = _TestResult && doc_test_int(rep.terminationtype, 1);
                    _TestResult = _TestResult && doc_test_real_vector(x, new double[]{1.000,2.000,1.000,2.000,1.000}, 0.005);
                    _TestResult = _TestResult && (_spoil_scenario==-1);
                }
                catch(alglib.alglibexception)
                { _TestResult = _TestResult && (_spoil_scenario!=-1); }
                catch
                { throw; }
            }
            if( !_TestResult)
                System.Console.WriteLine("{0,-32} FAILED", "lincg_d_1");
            _TotalResult = _TotalResult && _TestResult;


            System.Console.WriteLine("151/151");
        }
        catch
        {
            System.Console.WriteLine("Unhandled exception was raised!");
            System.Environment.ExitCode = 1;
            return;
        }
        
        // This object is descendant of CriticalFinalizerObject class,
        // which guarantees that it will be finalized AFTER all other
        // ALGLIB objects which hold pointers to unmanaged memory.
        //
        // Tests for memory leaks aredone within object's destructor.
        //
        // IMPORTANT: these tests report results to console, but do
        //            not set application exit code because it can
        //            not be reliably set from finalizer!
        MemoryLeaksTest _test_object = new MemoryLeaksTest();
        if( !_TotalResult )
            System.Environment.ExitCode = 1;
    }
}
