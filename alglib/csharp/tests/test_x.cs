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
public class XTest
{
    public static void Main(string[] args)
    {
        bool _TotalResult = true;
        bool _TestResult;
        System.Console.WriteLine("x-tests. Please wait...");
        alglib.alloc_counter_activate();
        System.Console.WriteLine("Allocation counter activated...");
        try
        {
            const int max1d = 70;
            const int max2d = 40;
            
            System.Console.WriteLine("Basic tests:");
            {
                // deallocateimmediately()
                alglib.minlbfgsstate s;
                double[] x = new double[100];
                long cnt0, cnt1;
                cnt0 = alglib.alloc_counter();
                alglib.minlbfgscreate(x.Length, 10, x, out s);
                alglib.deallocateimmediately(ref s);
                cnt1 = alglib.alloc_counter();
                _TestResult = cnt1<=cnt0;
                System.Console.WriteLine("* deallocateimmediately()    "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // boolean 1D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int n, i, cnt;
                _TestResult = true;
                for(n=0; n<=max1d; n++)
                {
                    bool[] arr0 = new bool[n];
                    bool[] arr1 = new bool[n];
                    bool[] arr2 = new bool[n];
                    bool[] arr3 = null;
                    cnt = 0;
                    for(i=0; i<n; i++)
                    {
                        arr0[i] = alglib.math.randomreal()>0.5;
                        arr1[i] = arr0[i];
                        arr2[i] = arr0[i];
                        if( arr0[i] )
                            cnt++;
                    }
                    _TestResult = _TestResult && (alglib.xdebugb1count(arr0)==cnt);
                    alglib.xdebugb1not(ref arr1);
                    if( alglib.ap.len(arr1)==n )
                    {
                        for(i=0; i<n; i++)
                            _TestResult = _TestResult && (arr1[i]==!arr0[i]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugb1appendcopy(ref arr2);
                    if( alglib.ap.len(arr2)==2*n )
                    {
                        for(i=0; i<2*n; i++)
                            _TestResult = _TestResult && (arr2[i]==arr0[i%n]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugb1outeven(n, out arr3);
                    if( alglib.ap.len(arr3)==n )
                    {
                        for(i=0; i<n; i++)
                            _TestResult = _TestResult && (arr3[i]==(i%2==0));
                    }
                    else
                        _TestResult = false;
                }
                System.Console.WriteLine("* boolean 1D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // integer 1D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int n, i, sum;
                _TestResult = true;
                for(n=0; n<=max1d; n++)
                {
                    int[] arr0 = new int[n];
                    int[] arr1 = new int[n];
                    int[] arr2 = new int[n];
                    int[] arr3 = null;
                    sum = 0;
                    for(i=0; i<n; i++)
                    {
                        arr0[i] = alglib.math.randominteger(10);
                        arr1[i] = arr0[i];
                        arr2[i] = arr0[i];
                        sum+=arr0[i];
                    }
                    _TestResult = _TestResult && (alglib.xdebugi1sum(arr0)==sum);
                    alglib.xdebugi1neg(ref arr1);
                    if( alglib.ap.len(arr1)==n )
                    {
                        for(i=0; i<n; i++)
                            _TestResult = _TestResult && (arr1[i]==-arr0[i]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugi1appendcopy(ref arr2);
                    if( alglib.ap.len(arr2)==2*n )
                    {
                        for(i=0; i<2*n; i++)
                            _TestResult = _TestResult && (arr2[i]==arr0[i%n]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugi1outeven(n,out arr3);
                    if( alglib.ap.len(arr3)==n )
                    {
                        for(i=0; i<n; i++)
                            if( i%2==0 )
                                _TestResult = _TestResult && (arr3[i]==i);
                            else
                                _TestResult = _TestResult && (arr3[i]==0);
                    }
                    else
                        _TestResult = false;
                }
                System.Console.WriteLine("* integer 1D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // real 1D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int n, i;
                double sum;
                _TestResult = true;
                for(n=0; n<=max1d; n++)
                {
                    double[] arr0 = new double[n];
                    double[] arr1 = new double[n];
                    double[] arr2 = new double[n];
                    double[] arr3 = null;
                    sum = 0;
                    for(i=0; i<n; i++)
                    {
                        arr0[i] = alglib.math.randomreal()-0.5;
                        arr1[i] = arr0[i];
                        arr2[i] = arr0[i];
                        sum+=arr0[i];
                    }
                    _TestResult = _TestResult && (Math.Abs(alglib.xdebugr1sum(arr0)-sum)<1.0E-10);
                    alglib.xdebugr1neg(ref arr1);
                    if( alglib.ap.len(arr1)==n )
                    {
                        for(i=0; i<n; i++)
                            _TestResult = _TestResult && (Math.Abs(arr1[i]+arr0[i])<1.0E-10);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugr1appendcopy(ref arr2);
                    if( alglib.ap.len(arr2)==2*n )
                    {
                        for(i=0; i<2*n; i++)
                            _TestResult = _TestResult && (arr2[i]==arr0[i%n]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugr1outeven(n,out arr3);
                    if( alglib.ap.len(arr3)==n )
                    {
                        for(i=0; i<n; i++)
                            if( i%2==0 )
                                _TestResult = _TestResult && (arr3[i]==i*0.25);
                            else
                                _TestResult = _TestResult && (arr3[i]==0);
                    }
                    else
                        _TestResult = false;
                }
                System.Console.WriteLine("* real 1D arrays             "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // complex 1D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int n, i;
                alglib.complex sum;
                _TestResult = true;
                for(n=0; n<=max1d; n++)
                {
                    alglib.complex[] arr0 = new alglib.complex[n];
                    alglib.complex[] arr1 = new alglib.complex[n];
                    alglib.complex[] arr2 = new alglib.complex[n];
                    alglib.complex[] arr3 = null;
                    sum = 0;
                    for(i=0; i<n; i++)
                    {
                        arr0[i].x = alglib.math.randomreal()-0.5;
                        arr0[i].y = alglib.math.randomreal()-0.5;
                        arr1[i] = arr0[i];
                        arr2[i] = arr0[i];
                        sum+=arr0[i];
                    }
                    _TestResult = _TestResult && (alglib.math.abscomplex(alglib.xdebugc1sum(arr0)-sum)<1.0E-10);
                    alglib.xdebugc1neg(ref arr1);
                    if( alglib.ap.len(arr1)==n )
                    {
                        for(i=0; i<n; i++)
                            _TestResult = _TestResult && (alglib.math.abscomplex(arr1[i]+arr0[i])<1.0E-10);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugc1appendcopy(ref arr2);
                    if( alglib.ap.len(arr2)==2*n )
                    {
                        for(i=0; i<2*n; i++)
                            _TestResult = _TestResult && (arr2[i]==arr0[i%n]);
                    }
                    else
                        _TestResult = false;
                    alglib.xdebugc1outeven(n,out arr3);
                    if( alglib.ap.len(arr3)==n )
                    {
                        for(i=0; i<n; i++)
                            if( i%2==0 )
                            {
                                _TestResult = _TestResult && (arr3[i].x==i*0.250);
                                _TestResult = _TestResult && (arr3[i].y==i*0.125);
                            }
                            else
                                _TestResult = _TestResult && (arr3[i]==0);
                    }
                    else
                        _TestResult = false;
                }
                System.Console.WriteLine("* complex 1D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // boolean 2D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int m, n, i, j, cnt;
                _TestResult = true;
                for(n=0; n<=max2d; n++)
                    for(m=0; m<=max2d; m++)
                    {
                        // skip situations when n*m==0, but n!=0 or m!=0
                        if( n*m==0 && (n!=0 || m!=0) )
                            continue;
                        
                        // proceed to testing
                        bool[,] arr0 = new bool[m,n];
                        bool[,] arr1 = new bool[m,n];
                        bool[,] arr2 = new bool[m,n];
                        bool[,] arr3 = null;
                        cnt = 0;
                        for(i=0; i<m; i++)
                            for(j=0; j<n; j++)
                            {
                                arr0[i,j] = alglib.math.randomreal()>0.5;
                                arr1[i,j] = arr0[i,j];
                                arr2[i,j] = arr0[i,j];
                                if( arr0[i,j] )
                                    cnt++;
                            }
                        _TestResult = _TestResult && (alglib.xdebugb2count(arr0)==cnt);
                        alglib.xdebugb2not(ref arr1);
                        if( alglib.ap.rows(arr1)==m && alglib.ap.cols(arr1)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr1[i,j]==!arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugb2transpose(ref arr2);
                        if( alglib.ap.rows(arr2)==n && alglib.ap.cols(arr2)==m )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr2[j,i]==arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugb2outsin(m, n, out arr3);
                        if( alglib.ap.rows(arr3)==m && alglib.ap.cols(arr3)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr3[i,j]==(Math.Sin(3*i+5*j)>0));
                        }
                        else
                            _TestResult = false;
                    }
                System.Console.WriteLine("* boolean 2D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // integer 2D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int m, n, i, j;
                int sum;
                _TestResult = true;
                for(n=0; n<=max2d; n++)
                    for(m=0; m<=max2d; m++)
                    {
                        // skip situations when n*m==0, but n!=0 or m!=0
                        if( n*m==0 && (n!=0 || m!=0) )
                            continue;
                        
                        // proceed to testing
                        int[,] arr0 = new int[m,n];
                        int[,] arr1 = new int[m,n];
                        int[,] arr2 = new int[m,n];
                        int[,] arr3 = null;
                        sum = 0;
                        for(i=0; i<m; i++)
                            for(j=0; j<n; j++)
                            {
                                arr0[i,j] = alglib.math.randominteger(10);
                                arr1[i,j] = arr0[i,j];
                                arr2[i,j] = arr0[i,j];
                                sum += arr0[i,j];
                            }
                        _TestResult = _TestResult && (alglib.xdebugi2sum(arr0)==sum);
                        alglib.xdebugi2neg(ref arr1);
                        if( alglib.ap.rows(arr1)==m && alglib.ap.cols(arr1)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr1[i,j]==-arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugi2transpose(ref arr2);
                        if( alglib.ap.rows(arr2)==n && alglib.ap.cols(arr2)==m )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr2[j,i]==arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugi2outsin(m, n, out arr3);
                        if( alglib.ap.rows(arr3)==m && alglib.ap.cols(arr3)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr3[i,j]==System.Math.Sign(Math.Sin(3*i+5*j)));
                        }
                        else
                            _TestResult = false;
                    }
                System.Console.WriteLine("* integer 2D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // real 2D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int m, n, i, j;
                double sum;
                _TestResult = true;
                for(n=0; n<=max2d; n++)
                    for(m=0; m<=max2d; m++)
                    {
                        // skip situations when n*m==0, but n!=0 or m!=0
                        if( n*m==0 && (n!=0 || m!=0) )
                            continue;
                        
                        // proceed to testing
                        double[,] arr0 = new double[m,n];
                        double[,] arr1 = new double[m,n];
                        double[,] arr2 = new double[m,n];
                        double[,] arr3 = null;
                        sum = 0;
                        for(i=0; i<m; i++)
                            for(j=0; j<n; j++)
                            {
                                arr0[i,j] = alglib.math.randomreal()-0.5;
                                arr1[i,j] = arr0[i,j];
                                arr2[i,j] = arr0[i,j];
                                sum += arr0[i,j];
                            }
                        _TestResult = _TestResult && (System.Math.Abs(alglib.xdebugr2sum(arr0)-sum)<1.0E-10);
                        alglib.xdebugr2neg(ref arr1);
                        if( alglib.ap.rows(arr1)==m && alglib.ap.cols(arr1)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr1[i,j]==-arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugr2transpose(ref arr2);
                        if( alglib.ap.rows(arr2)==n && alglib.ap.cols(arr2)==m )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr2[j,i]==arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugr2outsin(m, n, out arr3);
                        if( alglib.ap.rows(arr3)==m && alglib.ap.cols(arr3)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (System.Math.Abs(arr3[i,j]-Math.Sin(3*i+5*j))<1E-10);
                        }
                        else
                            _TestResult = false;
                    }
                System.Console.WriteLine("* real 2D arrays             "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // real 2D arrays (this test checks both interface and ref/out conventions used by ALGLIB)
                int m, n, i, j;
                alglib.complex sum;
                _TestResult = true;
                for(n=0; n<=max2d; n++)
                    for(m=0; m<=max2d; m++)
                    {
                        // skip situations when n*m==0, but n!=0 or m!=0
                        if( n*m==0 && (n!=0 || m!=0) )
                            continue;
                        
                        // proceed to testing
                        alglib.complex[,] arr0 = new alglib.complex[m,n];
                        alglib.complex[,] arr1 = new alglib.complex[m,n];
                        alglib.complex[,] arr2 = new alglib.complex[m,n];
                        alglib.complex[,] arr3 = null;
                        sum = 0;
                        for(i=0; i<m; i++)
                            for(j=0; j<n; j++)
                            {
                                arr0[i,j].x = alglib.math.randomreal()-0.5;
                                arr0[i,j].y = alglib.math.randomreal()-0.5;
                                arr1[i,j] = arr0[i,j];
                                arr2[i,j] = arr0[i,j];
                                sum += arr0[i,j];
                            }
                        _TestResult = _TestResult && (alglib.math.abscomplex(alglib.xdebugc2sum(arr0)-sum)<1.0E-10);
                        alglib.xdebugc2neg(ref arr1);
                        if( alglib.ap.rows(arr1)==m && alglib.ap.cols(arr1)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr1[i,j]==-arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugc2transpose(ref arr2);
                        if( alglib.ap.rows(arr2)==n && alglib.ap.cols(arr2)==m )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                    _TestResult = _TestResult && (arr2[j,i]==arr0[i,j]);
                        }
                        else
                            _TestResult = false;
                        alglib.xdebugc2outsincos(m, n, out arr3);
                        if( alglib.ap.rows(arr3)==m && alglib.ap.cols(arr3)==n )
                        {
                            for(i=0; i<m; i++)
                                for(j=0; j<n; j++)
                                {
                                    _TestResult = _TestResult && (System.Math.Abs(arr3[i,j].x-Math.Sin(3*i+5*j))<1E-10);
                                    _TestResult = _TestResult && (System.Math.Abs(arr3[i,j].y-Math.Cos(3*i+5*j))<1E-10);
                                }
                        }
                        else
                            _TestResult = false;
                    }
                System.Console.WriteLine("* complex 2D arrays          "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                // "biased product / sum" test
                int m, n, i, j;
                double sum;
                _TestResult = true;
                for(n=1; n<=max2d; n++)
                    for(m=1; m<=max2d; m++)
                    {
                        // proceed to testing
                        double[,] a = new double[m,n];
                        double[,] b = new double[m,n];
                        bool[,]   c = new bool[m,n];
                        sum = 0;
                        for(i=0; i<m; i++)
                            for(j=0; j<n; j++)
                            {
                                a[i,j] = alglib.math.randomreal()-0.5;
                                b[i,j] = alglib.math.randomreal()-0.5;
                                c[i,j] = alglib.math.randomreal()>0.5;
                                if( c[i,j] )
                                    sum += a[i,j]*(1+b[i,j]);
                            }
                        _TestResult = _TestResult && (Math.Abs(alglib.xdebugmaskedbiasedproductsum(m,n,a,b,c)-sum)<1.0E-10);
                    }
                System.Console.WriteLine("* multiple arrays            "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            {
                //
                // Test multithreading-related settings
                //
                // For this test we measure performance of large NxNxN GEMMs
                // with different threading settings.
                //
                System.Console.WriteLine("SMP settings vs GEMM speedup:");
                if( alglib.smp.cores_count>1 )
                {
                    bool passed = true;
                    int n = 800, mintime = 2000, t0;
                    ulong default_global_threading = alglib.ae_get_global_threading();
                    int default_nworkers = alglib.getnworkers();
                    double time_default          = 0,
                           time_glob_ser         = 0,
                           time_glob_smp         = 0,
                           time_glob_ser_loc_ser = 0,
                           time_glob_ser_loc_smp = 0,
                           time_glob_smp_loc_ser = 0,
                           time_glob_smp_loc_smp = 0,
                           time_glob_smp_nw1     = 0;
                    try
                    {
                        // allocate temporary matrices
                        double[,] a = new double[n,n], b = new double[n,n], c = new double[n,n];
                        int i, j;
                        for(i=0; i<n; i++)
                            for(j=0; j<n; j++)
                            {
                                a[i,j] = alglib.math.randomreal()-0.5;
                                b[i,j] = alglib.math.randomreal()-0.5;
                                c[i,j] = 0.0;
                            }
                        
                        // measure time; interleave measurements with different settings in order to
                        // reduce variance of results
                        while(time_default<mintime)
                        {
                            // default threading
                            t0 = System.Environment.TickCount;
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0);
                            time_default += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global serial
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.serial);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0);
                            time_glob_ser += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global parallel
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.parallel);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0);
                            time_glob_smp += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global serial, local serial
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.serial);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0,
                                alglib.serial);
                            time_glob_ser_loc_ser += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global serial, local parallel
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.serial);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0,
                                alglib.parallel);
                            time_glob_ser_loc_smp += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global parallel, local serial
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.parallel);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0,
                                alglib.serial);
                            time_glob_smp_loc_ser += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global parallel, local parallel
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.parallel);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0,
                                alglib.parallel);
                            time_glob_smp_loc_smp += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            
                            // global parallel, nworkers 1
                            t0 = System.Environment.TickCount;
                            alglib.setglobalthreading(alglib.parallel);
                            alglib.setnworkers(1);
                            alglib.rmatrixgemm(
                                n, n, n,
                                1.0,
                                a, 0, 0, 0,
                                b, 0, 0, 0,
                                0.0,
                                ref c, 0, 0);
                            time_glob_smp_nw1 += System.Environment.TickCount-t0;
                            alglib.ae_set_global_threading(default_global_threading); // restore
                            alglib.setnworkers(default_nworkers);
                        }
                    }
                    catch
                    { passed = false; }
                    System.Console.WriteLine("* default speedup         {0,5:F1}x", 1.0);
                    System.Console.WriteLine("* serial (global)         {0,5:F1}x", time_glob_ser/time_default);
                    System.Console.WriteLine("* serial (local)          {0,5:F1}x", time_glob_ser/time_glob_ser_loc_ser);
                    System.Console.WriteLine("* serial (nworkers=1)     {0,5:F1}x", time_glob_ser/time_glob_smp_nw1);
                    System.Console.WriteLine("* parallel (global)       {0,5:F1}x", time_glob_ser/time_glob_smp);
                    System.Console.WriteLine("* parallel (local) v1     {0,5:F1}x", time_glob_ser/time_glob_ser_loc_smp);
                    passed = passed && (time_glob_ser/time_default         >0.85) && (time_glob_ser/time_default         <1.15);
                    passed = passed && (time_glob_ser/time_glob_ser        >0.85) && (time_glob_ser/time_glob_ser        <1.15);
                    passed = passed && (time_glob_ser/time_glob_ser_loc_ser>0.85) && (time_glob_ser/time_glob_ser_loc_ser<1.15);
                    passed = passed && (time_glob_ser/time_glob_smp_loc_ser>0.85) && (time_glob_ser/time_glob_smp_loc_ser<1.15);
                    passed = passed && (time_glob_ser/time_glob_smp_nw1    >0.85) && (time_glob_ser/time_glob_smp_nw1    <1.15);
                    passed = passed && (time_glob_ser/time_glob_smp        >1.30);
                    passed = passed && (time_glob_ser/time_glob_ser_loc_smp>1.30);
                    passed = passed && (time_glob_ser/time_glob_smp_loc_smp>1.30);
                    System.Console.WriteLine("* test result                 "+(passed ? "OK" : "FAILED (soft failure)"));
                    //
                    // soft failure:
                    // // if( !passed )
                    // //   return 1;
                    //
                }
                else
                    System.Console.WriteLine("* test skipped (no SMP)       "+"??");
            }
            
            
            //////////////////////////////////
            // Advanced tests
            //////
            System.Console.WriteLine("Advanced tests:");

            //
            // Testing CSV functionality
            //
            {
                string csv_name = "alglib-tst-35252-ndg4sf.csv";
                _TestResult = true;
                try
                {
                    // CSV_DEFAULT must be zero
                    _TestResult = _TestResult && alglib.CSV_DEFAULT==0;
                    
                    // absent file - must fail
                    try
                    {
                        double[,] arr;
                        alglib.read_csv("nonexistent123foralgtestinglib", '\t', alglib.CSV_DEFAULT, out arr);
                        _TestResult = false;
                    }
                    catch
                    { }
                    
                    // non-rectangular file - must fail
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "a,b,c\r\n1,2");
                        alglib.read_csv(csv_name, ',', alglib.CSV_SKIP_HEADERS, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = false;
                    }
                    catch
                    { }
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "a,b,c\r\n1,2,3,4");
                        alglib.read_csv(csv_name, ',', alglib.CSV_SKIP_HEADERS, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = false;
                    }
                    catch
                    { }
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "1,2,3,4\n1,2,3\n1,2,3");
                        alglib.read_csv(csv_name, ',', alglib.CSV_DEFAULT, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = false;
                    }
                    catch
                    { }
                    
                    // empty file
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "");
                        alglib.read_csv(csv_name, '\t', alglib.CSV_DEFAULT, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = _TestResult && arr.GetLength(0)==0 && arr.GetLength(1)==0;
                    }
                    catch
                    { _TestResult = false; }
                    
                    // one row with header, tab separator
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "a\tb\tc\n");
                        alglib.read_csv(csv_name, '\t', alglib.CSV_SKIP_HEADERS, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = _TestResult && arr.GetLength(0)==0 && arr.GetLength(1)==0;
                    }
                    catch
                    { _TestResult = false; }
                    
                    // no header, comma-separated, full stop as decimal point
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "1.5,2,3.25\n4,5,6");
                        alglib.read_csv(csv_name, ',', alglib.CSV_DEFAULT, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = _TestResult && alglib.ap.format(arr,2)=="{{1.50,2.00,3.25},{4.00,5.00,6.00}}";
                    }
                    catch
                    { _TestResult = false; }
                    
                    // header, tab-separated, mixed use of comma and full stop as decimal points
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, "a\tb\tc\n1.5\t2\t3,25\n4\t5.25\t6,1\n");
                        alglib.read_csv(csv_name, '\t', alglib.CSV_SKIP_HEADERS, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = _TestResult && alglib.ap.format(arr,2)=="{{1.50,2.00,3.25},{4.00,5.25,6.10}}";
                    }
                    catch
                    { _TestResult = false; }
                    
                    // header, tab-separated, fixed/exponential, spaces, mixed use of comma and full stop as decimal points
                    try
                    {
                        double[,] arr;
                        System.IO.File.WriteAllText(csv_name, " a\t b \tc\n1,1\t 2.9\t -3.5  \n  1.1E1  \t 2.0E-1 \t-3E+1 \n+1  \t -2\t 3.    \n.1\t-.2\t+.3\n");
                        alglib.read_csv(csv_name, '\t', alglib.CSV_SKIP_HEADERS, out arr);
                        System.IO.File.Delete(csv_name);
                        _TestResult = _TestResult && alglib.ap.format(arr,2)=="{{1.10,2.90,-3.50},{11.00,0.20,-30.00},{1.00,-2.00,3.00},{0.10,-0.20,0.30}}";
                    }
                    catch
                    { _TestResult = false; }
                }
                catch
                {
                    _TestResult = false;
                }
                
                //
                // Report
                //
                System.Console.WriteLine("* CSV support                "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }

            //
            // Testing serialization functionality (using kd-trees as playground)
            //
            {
                
                _TestResult = true;
                try
                {
                    // prepare data
                    alglib.hqrndstate rs;
                    alglib.kdtree tree0;
                    double[,] xy, rxy0 = new double[0,0], rxy1 = new double[0,0];
                    double[]  qx;
                    const int npts = 50;
                    const int nx = 2;
                    const int ny = 1;
                    int cnt0, cnt1;
                    alglib.hqrndrandomize(out rs);
                    xy = new double[npts,nx+ny];
                    for(int i=0; i<npts; i++)
                        for(int j=0; j<nx+ny; j++)
                            xy[i,j] = alglib.hqrndnormal(rs);
                    alglib.kdtreebuild(xy, npts, nx, ny, 2, out tree0);
                    qx = new double[nx];
                    
                    try
                    {
                        // test string serialization/unserialization
                        alglib.kdtree tree1;
                        string s;
                        alglib.kdtreeserialize(tree0, out s);
                        alglib.kdtreeunserialize(s, out tree1);
                        for(int i=0; i<100; i++)
                        {
                            for(int j=0; j<nx; j++)
                                qx[j] = alglib.hqrndnormal(rs);
                            cnt0 = alglib.kdtreequeryknn(tree0, qx, 1, true);
                            cnt1 = alglib.kdtreequeryknn(tree1, qx, 1, true);
                            if( (cnt0!=1) || (cnt1!=1) )
                            {
                                _TestResult = false;
                                break;
                            }
                            alglib.kdtreequeryresultsxy(tree0, ref rxy0);
                            alglib.kdtreequeryresultsxy(tree1, ref rxy1);
                            for(int j=0; j<nx+ny; j++)
                                _TestResult = _TestResult && (rxy0[0,j]==rxy1[0,j]);
                        }
                    }
                    catch
                    { _TestResult = false; }
                    
                    try
                    {
                        // test stream serialization/unserialization
                        //
                        // NOTE: we add a few symbols at the beginning and after the end of the data
                        //       in order to test algorithm ability to work in the middle of the stream
                        alglib.kdtree tree1;
                        System.IO.MemoryStream s = new System.IO.MemoryStream();
                        s.WriteByte((byte)'b');
                        s.WriteByte((byte)' ');
                        s.WriteByte((byte)'e');
                        s.WriteByte((byte)'g');
                        alglib.kdtreeserialize(tree0, s);
                        s.WriteByte((byte)'@');
                        s.WriteByte((byte)' ');
                        s.WriteByte((byte)'n');
                        s.WriteByte((byte)'d');
                        s.Seek(0, System.IO.SeekOrigin.Begin);
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'b');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)' ');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'e');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'g');
                        alglib.kdtreeunserialize(s, out tree1);
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'@');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)' ');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'n');
                        _TestResult = _TestResult && (s.ReadByte()==(byte)'d');
                        for(int i=0; i<100; i++)
                        {
                            for(int j=0; j<nx; j++)
                                qx[j] = alglib.hqrndnormal(rs);
                            cnt0 = alglib.kdtreequeryknn(tree0, qx, 1, true);
                            cnt1 = alglib.kdtreequeryknn(tree1, qx, 1, true);
                            if( (cnt0!=1) || (cnt1!=1) )
                            {
                                _TestResult = false;
                                break;
                            }
                            alglib.kdtreequeryresultsxy(tree0, ref rxy0);
                            alglib.kdtreequeryresultsxy(tree1, ref rxy1);
                            for(int j=0; j<nx+ny; j++)
                                _TestResult = _TestResult && (rxy0[0,j]==rxy1[0,j]);
                        }
                    }
                    catch
                    { _TestResult = false; }
                    
                    try
                    {
                        // test string-to-stream serialization/unserialization
                        alglib.kdtree tree1;
                        string s0;
                        alglib.kdtreeserialize(tree0, out s0);
                        System.IO.MemoryStream s1 = new System.IO.MemoryStream(System.Text.Encoding.UTF8.GetBytes(s0));
                        alglib.kdtreeunserialize(s1, out tree1);
                        for(int i=0; i<100; i++)
                        {
                            for(int j=0; j<nx; j++)
                                qx[j] = alglib.hqrndnormal(rs);
                            cnt0 = alglib.kdtreequeryknn(tree0, qx, 1, true);
                            cnt1 = alglib.kdtreequeryknn(tree1, qx, 1, true);
                            if( (cnt0!=1) || (cnt1!=1) )
                            {
                                _TestResult = false;
                                break;
                            }
                            alglib.kdtreequeryresultsxy(tree0, ref rxy0);
                            alglib.kdtreequeryresultsxy(tree1, ref rxy1);
                            for(int j=0; j<nx+ny; j++)
                                _TestResult = _TestResult && (rxy0[0,j]==rxy1[0,j]);
                        }
                    }
                    catch
                    { _TestResult = false; }
                    
                    try
                    {
                        // test stream-to-string serialization/unserialization
                        alglib.kdtree tree1;
                        System.IO.MemoryStream s0 = new System.IO.MemoryStream();
                        alglib.kdtreeserialize(tree0, s0);
                        s0.Seek(0, System.IO.SeekOrigin.Begin);
                        string s1 = System.Text.Encoding.UTF8.GetString(s0.ToArray());
                        alglib.kdtreeunserialize(s1, out tree1);
                        for(int i=0; i<100; i++)
                        {
                            for(int j=0; j<nx; j++)
                                qx[j] = alglib.hqrndnormal(rs);
                            cnt0 = alglib.kdtreequeryknn(tree0, qx, 1, true);
                            cnt1 = alglib.kdtreequeryknn(tree1, qx, 1, true);
                            if( (cnt0!=1) || (cnt1!=1) )
                            {
                                _TestResult = false;
                                break;
                            }
                            alglib.kdtreequeryresultsxy(tree0, ref rxy0);
                            alglib.kdtreequeryresultsxy(tree1, ref rxy1);
                            for(int j=0; j<nx+ny; j++)
                                _TestResult = _TestResult && (rxy0[0,j]==rxy1[0,j]);
                        }
                    }
                    catch
                    { _TestResult = false; }

                }
                catch
                {
                    _TestResult = false;
                }
                
                //
                // Report
                //
                System.Console.WriteLine("* Serialization (kd-tree)    "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            
            //////////////////////////////////
            // Test issues from Mantis
            //////
            System.Console.WriteLine("Testing issies from Mantis:");
                
            
            //
            // Task #594 (http://bugs.alglib.net/view.php?id=594) - additional
            // test for correctness of copying of objects. When we copy ALGLIB
            // object, indenendent new copy is created.
            //
            {
                //
                // First, test copying of alglib.multilayerperceptron, which
                // is an "opaque object".
                //
                // Test copy constructors:
                // * copy object with make_copy()
                // * process vector with original network
                // * randomize original network
                // * process vector with copied networks and compare
                //
                alglib.multilayerperceptron net0, net1;
                double[] x  = new double[]{1,2};
                double[] y0 = new double[]{0,0};
                double[] y1 = new double[]{0,0};
                double[] y2 = new double[]{0,0};
                _TestResult = true;
                alglib.mlpcreate0(2, 2, out net0);
                alglib.mlpprocess(net0, x, ref y0);
                net1 = (alglib.multilayerperceptron)net0.make_copy();
                alglib.mlprandomize(net0);
                alglib.mlpprocess(net1, x, ref y1);
                _TestResult = _TestResult && (Math.Abs(y0[0]-y1[0])<1.0E-9) && (Math.Abs(y0[1]-y1[1])<1.0E-9);
                
                //
                // Then, test correctness of copying "records", i.e.
                // objects with publicly visible fields.
                //
                alglib.xdebugrecord1 r0, r1;
                alglib.xdebuginitrecord1(out r0);
                r1 = (alglib.xdebugrecord1)r0.make_copy();
                _TestResult = _TestResult && (r1.i==r0.i);
                _TestResult = _TestResult && (r1.c==r0.c);
                
                _TestResult = _TestResult && (r1.a.Length==2);
                _TestResult = _TestResult && (r0.a.Length==2);
                _TestResult = _TestResult && (r1.a!=r0.a);
                _TestResult = _TestResult && (r1.a[0]==r0.a[0]);
                _TestResult = _TestResult && (r1.a[1]==r0.a[1]);
                
                //
                // Test result
                //
                System.Console.WriteLine("* issue 594                  "+(_TestResult ? " OK" : " FAILED"));
                _TotalResult = _TotalResult && _TestResult;
            }
            
        }
        catch
        {
            System.Console.WriteLine("Unhandled exception was raised!");
            System.Environment.ExitCode = 1;
            return;
        }
        
            
        //////////////////////////////////
        // Backward compatibility tests
        //////
        System.Console.WriteLine("Backward compatibility tests:");

        //
        // Testing RBF storage format
        //
        {
            double eps = 0.0000000001;
            double[] ref_val = new double[]{
                -0.042560546916643,
                 0.942523544654062,
                 0.875197036560778,
                 0.0656948997826632,
                -0.743065973803404,
                -0.8903682039297,
                -0.26994815318748,
                 0.602248517290195,
                 0.980011992233124,
                 0.436594293214176
                };
            string _ss = @"50000000000 00000000000 20000000000 10000000000 A0000000000
30000000000 20000000000 00000000000 A0000000000 30000000000
00000000000 20000000000 A0000000000 60000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000m_3 00000000000 00000000000 00000000m_3 00000000000
00000000000 00000000004 00000000000 00000000000 00000000004
00000000000 00000000000 00000000804 00000000000 00000000000
00000000804 00000000000 00000000000 00000000G04 00000000000
00000000000 00000000G04 00000000000 00000000000 00000000O04
00000000000 00000000000 00000000O04 00000000000 00000000000
00000000S04 00000000000 00000000000 00000000S04 00000000000
00000000000 00000000W04 00000000000 00000000000 00000000W04
00000000000 00000000000 00000000Y04 00000000000 00000000000
00000000Y04 00000000000 00000000000 00000000K04 00000000000
00000000000 00000000K04 00000000000 00000000000 A0000000000
00000000000 10000000000 20000000000 30000000000 40000000000
60000000000 70000000000 80000000000 90000000000 50000000000
30000000000 00000000000 00000000000 00000000000 30000000000
00000000Y04 00000000000 00000000000 u1000000000 00000000000
00000000000 00000000000 60000000000 80000000000 00000000000
50000000000 00000000000 50000000000 50000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 K0000000000
00000000I04 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
00000000000 00000000000 00000000000 00000000000 00000000000
A0000000000 30000000000 00000000000 00000000000 00000000000
00000000m_3 00000000000 00000000000 00000000004 00000000000
00000000000 00000000804 00000000000 00000000000 00000000G04
00000000000 00000000000 00000000K04 00000000000 00000000000
00000000O04 00000000000 00000000000 00000000S04 00000000000
00000000000 00000000W04 00000000000 00000000000 00000000Y04
00000000000 00000000000 A0000000000 40000000000 00000000q04
-pAGQnQBI14 UqUWierJ91C esm8ag6G61C 00000000q04 4wcFMyCtu04
oPDvwHqst04 CExQXp8Ct04 00000000q04 litzPFhRb0C oKJvjcct314
5-fT-X8w614 00000000q04 3HSOsPVH11C vZWf4dgfv04 GbZg4MTJn04
00000000q04 iv7rMhuR71C hRtixp15r_3 EvCEDtLu-0C 00000000q04
41CXzA_q71C umRYLK2yp0C 1zzY3Zqd91C 00000000q04 JvxJzDeI21C
TVbyd7Ygz0C JLywRdR1n0C 00000000q04 KmFarhc4g0C 1ehrn2tUt0C
AECfwTIX814 00000000q04 Big__6hwt04 nSPzmAQrh_B 2H3o-KftH14
00000000q04 n1b9361vI14 mhJhviUE114 54a_qyBrH1C 00000000q04
10000000000 40000000000 StLCgor39-3 00000000000 00000000000
6qTG7Ae-1_3
";

            // test string unserialization without trailing dot symbol (end-of-stream marker); must work
            try
            {
                string s = _ss;
                alglib.rbfmodel model;
                alglib.rbfunserialize(s, out model);
                _TestResult = true;
                for(int i=0; i<ref_val.Length; i++)
                    _TestResult = _TestResult && (Math.Abs(alglib.rbfcalc2(model,i,0)-ref_val[i])<eps);
            }
            catch
            { _TestResult = false; }

            // test string unserialization with trailing dot symbol (end-of-stream marker); must work
            try
            {
                string s = _ss+".";
                alglib.rbfmodel model;
                alglib.rbfunserialize(s, out model);
                _TestResult = true;
                for(int i=0; i<ref_val.Length; i++)
                    _TestResult = _TestResult && (Math.Abs(alglib.rbfcalc2(model,i,0)-ref_val[i])<eps);
            }
            catch
            { _TestResult = false; }

            // test stream unserialization with trailing dot symbol (end-of-stream marker); must work
            try
            {
                System.IO.Stream s = new System.IO.MemoryStream(System.Text.Encoding.UTF8.GetBytes(_ss+"."));
                alglib.rbfmodel model;
                alglib.rbfunserialize(s, out model);
                _TestResult = true;
                for(int i=0; i<ref_val.Length; i++)
                    _TestResult = _TestResult && (Math.Abs(alglib.rbfcalc2(model,i,0)-ref_val[i])<eps);
            }
            catch
            { _TestResult = false; }

            // test stream unserialization with trailing dot symbol (end-of-stream marker); MUST FAIL
            try
            {
                System.IO.Stream s = new System.IO.MemoryStream(System.Text.Encoding.UTF8.GetBytes(_ss));
                alglib.rbfmodel model;
                alglib.rbfunserialize(s, out model);
                _TestResult = false; // must FAIL!
            }
            catch
            {  }
            
            //
            // Report
            //
            System.Console.WriteLine("* RBF bwd compatibility      "+(_TestResult ? " OK" : " FAILED"));
            _TotalResult = _TotalResult && _TestResult;
        }
        
        //
        // Performance tests
        //
        System.Console.WriteLine("Performance:");
        {
            {
                int[] _n = new int[]{ 16, 32, 64, 256, 1024, 0};
                int i, j, k, t, nidx;
                for(nidx=0; _n[nidx]!=0; nidx++)
                {
                    //
                    // Settings:
                    // * n - matrix size
                    // * nrepeat - number of repeated multiplications, always divisible by 4
                    //
                    int n = _n[nidx];
                    double desiredflops = n>64 ? 1.0E10 : 1.0E9;
                    int nrepeat = (int)(desiredflops/(2*System.Math.Pow(n,3.0)));
                    nrepeat = 4*(nrepeat/4+1);
                    
                    //
                    // Actual processing
                    //
                    double[,] a, b, c;
                    double perf0, perf1, perf2;
                    a = new double[n,n];
                    b = new double[n,n];
                    c = new double[n,n];
                    for(i=0; i<n; i++)
                        for(j=0; j<n; j++)
                        {
                            a[i,j] = alglib.math.randomreal()-0.5;
                            b[i,j] = alglib.math.randomreal()-0.5;
                            c[i,j] = 0.0;
                        }
                    
                    alglib.setnworkers(0);
                    t = System.Environment.TickCount;
                    for(k=0; k<nrepeat; k++)
                        alglib.rmatrixgemm(
                            n, n, n,
                            1.0,
                            a, 0, 0, k%2,
                            b, 0, 0, (k/2)%2,
                            0.0,
                            ref c, 0, 0);
                    t = System.Environment.TickCount-t;
                    perf0 = 1.0E-6*System.Math.Pow(n,3)*2.0*nrepeat/(0.001*t);
                    System.Console.WriteLine("* RGEMM-SEQ-{0,-4:D} (MFLOPS)  {1,5:F0}", n, perf0);
                    
                    alglib.setnworkers(0);
                    t = System.Environment.TickCount;
                    for(k=0; k<nrepeat; k++)
                        alglib.rmatrixgemm(
                            n, n, n,
                            1.0,
                            a, 0, 0, k%2,
                            b, 0, 0, (k/2)%2,
                            0.0,
                            ref c, 0, 0,
                            alglib.parallel);
                    t = System.Environment.TickCount-t;
                    perf2 = 1.0E-6*System.Math.Pow(n,3)*2.0*nrepeat/(0.001*t);
                    System.Console.WriteLine("* RGEMM-MTN-{0,-4:D}           {1,4:F1}x", n, perf2/perf0);
                }
            }
        }
        
        
        //
        // Test below creates instance of MemoryLeaksTest object.
        //
        // This object is descendant of CriticalFinalizerObject class,
        // which guarantees that it will be finalized AFTER all other
        // ALGLIB objects which hold pointers to unmanaged memory.
        //
        // Tests for memory leaks are done within object's destructor.
        //
        MemoryLeaksTest _test_object = new MemoryLeaksTest();
        if( !_TotalResult )
            System.Environment.ExitCode = 1;
    }
}
