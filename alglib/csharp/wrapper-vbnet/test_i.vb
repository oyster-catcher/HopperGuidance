
Module MainTest

    Public Function doc_test_bool(v As Boolean, t As Boolean) As Boolean
        Return (v AndAlso t) OrElse (Not v AndAlso Not t)
    End Function

    Public Function doc_test_int(v As Integer, t As Integer) As Boolean
        Return v = t
    End Function

    Public Function doc_test_real(val As Double, test_val As Double, _threshold As Double) As Boolean
        Dim s As Double
        If _threshold >= 0 Then
            s = 1.0
        Else
            s =Math.Abs(test_val)
        End If
        Dim threshold As Double = Math.Abs(_threshold)
        Return Math.Abs(val - test_val) / s <= threshold
    End Function

    Public Function doc_test_complex(val As alglib.complex, test_val As alglib.complex, _threshold As Double) As Boolean
        Dim s As Double
        If _threshold >= 0 Then
            s = 1.0
        Else
            s =alglib.math.abscomplex(test_val)
        End If
        Dim threshold As Double = Math.Abs(_threshold)
        Return alglib.math.abscomplex(val - test_val) / s <= threshold
    End Function

    Public Function doc_test_bool_vector(v As Boolean(), t As Boolean()) As Boolean
        Dim i As Integer
        If alglib.ap.len(v) <> alglib.ap.len(t) Then
            Return False
        End If
        For i = 0 To alglib.ap.len(v) - 1
            If v(i) <> t(i) Then
                Return False
            End If
        Next
        Return True
    End Function

    Public Function doc_test_bool_matrix(v As Boolean(,), t As Boolean(,)) As Boolean
        Dim i As Integer, j As Integer
        If alglib.ap.rows(v) <> alglib.ap.rows(t) Then
            Return False
        End If
        If alglib.ap.cols(v) <> alglib.ap.cols(t) Then
            Return False
        End If
        For i = 0 To alglib.ap.rows(v) - 1
            For j = 0 To alglib.ap.cols(v) - 1
                If v(i, j) <> t(i, j) Then
                    Return False
                End If
            Next
        Next
        Return True
    End Function

    Public Function doc_test_int_vector(v As Integer(), t As Integer()) As Boolean
        Dim i As Integer
        If alglib.ap.len(v) <> alglib.ap.len(t) Then
            Return False
        End If
        For i = 0 To alglib.ap.len(v) - 1
            If v(i) <> t(i) Then
                Return False
            End If
        Next
        Return True
    End Function

    Public Function doc_test_int_matrix(v As Integer(,), t As Integer(,)) As Boolean
        Dim i As Integer, j As Integer
        If alglib.ap.rows(v) <> alglib.ap.rows(t) Then
            Return False
        End If
        If alglib.ap.cols(v) <> alglib.ap.cols(t) Then
            Return False
        End If
        For i = 0 To alglib.ap.rows(v) - 1
            For j = 0 To alglib.ap.cols(v) - 1
                If v(i, j) <> t(i, j) Then
                    Return False
                End If
            Next
        Next
        Return True
    End Function

    Public Function doc_test_real_vector(val As Double(), test_val As Double(), _threshold As Double) As Boolean
        Dim i As Integer
        If alglib.ap.len(val) <> alglib.ap.len(test_val) Then
            Return False
        End If
        For i = 0 To alglib.ap.len(val) - 1
            Dim s As Double
            If _threshold >= 0 Then
                s = 1.0
            Else
                s =Math.Abs(test_val(i))
            End If
            Dim threshold As Double = Math.Abs(_threshold)
            If Math.Abs(val(i) - test_val(i)) / s > threshold Then
                Return False
            End If
        Next
        Return True
    End Function

    Public Function doc_test_real_matrix(val As Double(,), test_val As Double(,), _threshold As Double) As Boolean
        Dim i As Integer, j As Integer
        If alglib.ap.rows(val) <> alglib.ap.rows(test_val) Then
            Return False
        End If
        If alglib.ap.cols(val) <> alglib.ap.cols(test_val) Then
            Return False
        End If
        For i = 0 To alglib.ap.rows(val) - 1
            For j = 0 To alglib.ap.cols(val) - 1
                Dim s As Double
                If _threshold >= 0 Then
                    s = 1.0
                Else
                    s =Math.Abs(test_val(i,j))
                End If
                Dim threshold As Double = Math.Abs(_threshold)
                If Math.Abs(val(i, j) - test_val(i, j)) / s > threshold Then
                    Return False
                End If
            Next
        Next
        Return True
    End Function
    
    Public Function doc_test_complex_vector(val As alglib.complex(), test_val As alglib.complex(), _threshold As Double) As Boolean
        Dim i As Integer
        If alglib.ap.len(val) <> alglib.ap.len(test_val) Then
            Return False
        End If
        For i = 0 To alglib.ap.len(val) - 1
            Dim s As Double
            If _threshold >= 0 Then
                s = 1.0
            Else
                s =alglib.math.abscomplex(test_val(i))
            End If
            Dim threshold As Double = Math.Abs(_threshold)
            If alglib.math.abscomplex(val(i) - test_val(i)) / s > threshold Then
                Return False
            End If
        Next
        Return True
    End Function

    Public Function doc_test_complex_matrix(val As alglib.complex(,), test_val As alglib.complex(,), _threshold As Double) As Boolean
        Dim i As Integer, j As Integer
        If alglib.ap.rows(val) <> alglib.ap.rows(test_val) Then
            Return False
        End If
        If alglib.ap.cols(val) <> alglib.ap.cols(test_val) Then
            Return False
        End If
        For i = 0 To alglib.ap.rows(val) - 1
            For j = 0 To alglib.ap.cols(val) - 1
                Dim s As Double
                If _threshold >= 0 Then
                    s = 1.0
                Else
                    s =alglib.math.abscomplex(test_val(i,j))
                End If
                Dim threshold As Double = Math.Abs(_threshold)
                If alglib.math.abscomplex(val(i, j) - test_val(i, j)) / s > threshold Then
                    Return False
                End If
            Next
        Next
        Return True
    End Function

    Public Sub spoil_vector_by_adding_element(Of T As New)(ByRef x As T(), val As T)
        Dim i As Integer
        Dim y As T() = x
        x = New T(y.Length) {}
        For i = 0 To y.Length - 1
            x(i) = y(i)
        Next
        x(y.Length) = val
    End Sub

    Public Sub spoil_vector_by_deleting_element(Of T As New)(ByRef x As T())
        Dim i As Integer
        Dim y As T() = x
        x = New T(y.Length - 2) {}
        For i = 0 To y.Length - 2
            x(i) = y(i)
        Next
    End Sub

    Public Sub spoil_matrix_by_adding_row(Of T As New)(ByRef x As T(,), val As T)
        Dim i As Integer, j As Integer
        Dim y As T(,) = x
        x = New T(y.GetLength(0), y.GetLength(1) - 1) {}
        For i = 0 To y.GetLength(0) - 1
            For j = 0 To y.GetLength(1) - 1
                x(i, j) = y(i, j)
            Next
        Next
        For j = 0 To y.GetLength(1) - 1
            x(y.GetLength(0), j) = val
        Next
    End Sub

    Public Sub spoil_matrix_by_deleting_row(Of T As New)(ByRef x As T(,))
        Dim i As Integer, j As Integer
        Dim y As T(,) = x
        x = New T(y.GetLength(0) - 2, y.GetLength(1) - 1) {}
        For i = 0 To y.GetLength(0) - 2
            For j = 0 To y.GetLength(1) - 1
                x(i, j) = y(i, j)
            Next
        Next
    End Sub

    Public Sub spoil_matrix_by_adding_col(Of T As New)(ByRef x As T(,), val As T)
        Dim i As Integer, j As Integer
        Dim y As T(,) = x
        x = New T(y.GetLength(0) - 1, y.GetLength(1)) {}
        For i = 0 To y.GetLength(0) - 1
            For j = 0 To y.GetLength(1) - 1
                x(i, j) = y(i, j)
            Next
        Next
        For i = 0 To y.GetLength(0) - 1
            x(i, y.GetLength(1)) = val
        Next
    End Sub

    Public Sub spoil_matrix_by_deleting_col(Of T As New)(ByRef x As T(,))
        Dim i As Integer, j As Integer
        Dim y As T(,) = x
        x = New T(y.GetLength(0) - 1, y.GetLength(1) - 2) {}
        For i = 0 To y.GetLength(0) - 1
            For j = 0 To y.GetLength(1) - 2
                x(i, j) = y(i, j)
            Next
        Next
    End Sub

    Public Sub spoil_vector_by_value(Of T)(ByRef x As T(), val As T)
        If x.Length <> 0 Then
            x(alglib.math.randominteger(x.Length)) = val
        End If
    End Sub

    Public Sub spoil_matrix_by_value(Of T)(ByRef x As T(,), val As T)
        If x.GetLength(0) <> 0 AndAlso x.GetLength(1) <> 0 Then
            x(alglib.math.randominteger(x.GetLength(0)), alglib.math.randominteger(x.GetLength(1))) = val
        End If
    End Sub

    Public Sub function1_func(x As Double(), ByRef func As Double, obj As Object)
        '
        ' this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
    End Sub

    Public Sub function1_grad(x As Double(), ByRef func As Double, grad As Double(), obj As Object)
        '
        ' this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        ' and its derivatives df/d0 and df/dx1
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
        grad(0) = 400 * System.Math.Pow(x(0) + 3, 3)
        grad(1) = 4 * System.Math.Pow(x(1) - 3, 3)
    End Sub

    Public Sub function1_hess(x As Double(), ByRef func As Double, grad As Double(), hess As Double(,), obj As Object)
        '
        ' this callback calculates f(x0,x1) = 100*(x0+3)^4 + (x1-3)^4
        ' its derivatives df/d0 and df/dx1
        ' and its Hessian.
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
        grad(0) = 400 * System.Math.Pow(x(0) + 3, 3)
        grad(1) = 4 * System.Math.Pow(x(1) - 3, 3)
        hess(0, 0) = 1200 * System.Math.Pow(x(0) + 3, 2)
        hess(0, 1) = 0
        hess(1, 0) = 0
        hess(1, 1) = 12 * System.Math.Pow(x(1) - 3, 2)
    End Sub

    Public Sub function1_fvec(x As Double(), fi As Double(), obj As Object)
        '
        ' this callback calculates
        ' f0(x0,x1) = 100*(x0+3)^4,
        ' f1(x0,x1) = (x1-3)^4
        '
        fi(0) = 10 * System.Math.Pow(x(0) + 3, 2)
        fi(1) = System.Math.Pow(x(1) - 3, 2)
    End Sub

    Public Sub function1_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates
        ' f0(x0,x1) = 100*(x0+3)^4,
        ' f1(x0,x1) = (x1-3)^4
        ' and Jacobian matrix J = [dfi/dxj]
        '
        fi(0) = 10 * System.Math.Pow(x(0) + 3, 2)
        fi(1) = System.Math.Pow(x(1) - 3, 2)
        jac(0, 0) = 20 * (x(0) + 3)
        jac(0, 1) = 0
        jac(1, 0) = 0
        jac(1, 1) = 2 * (x(1) - 3)
    End Sub

    Public Sub function2_func(x As Double(), ByRef func As Double, obj As Object)
        '
        ' this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        '
        func = System.Math.Pow(x(0) * x(0) + 1, 2) + System.Math.Pow(x(1) - 1, 2)
    End Sub

    Public Sub function2_grad(x As Double(), ByRef func As Double, grad As Double(), obj As Object)
        '
        ' this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        ' and its derivatives df/d0 and df/dx1
        '
        func = System.Math.Pow(x(0) * x(0) + 1, 2) + System.Math.Pow(x(1) - 1, 2)
        grad(0) = 4 * (x(0) * x(0) + 1) * x(0)
        grad(1) = 2 * (x(1) - 1)
    End Sub

    Public Sub function2_hess(x As Double(), ByRef func As Double, grad As Double(), hess As Double(,), obj As Object)
        '
        ' this callback calculates f(x0,x1) = (x0^2+1)^2 + (x1-1)^2
        ' its gradient and Hessian
        '
        func = System.Math.Pow(x(0) * x(0) + 1, 2) + System.Math.Pow(x(1) - 1, 2)
        grad(0) = 4 * (x(0) * x(0) + 1) * x(0)
        grad(1) = 2 * (x(1) - 1)
        hess(0, 0) = 12 * x(0) * x(0) + 4
        hess(0, 1) = 0
        hess(1, 0) = 0
        hess(1, 1) = 2
    End Sub

    Public Sub function2_fvec(x As Double(), fi As Double(), obj As Object)
        '
        ' this callback calculates
        ' f0(x0,x1) = 100*(x0+3)^4,
        ' f1(x0,x1) = (x1-3)^4
        '
        fi(0) = x(0) * x(0) + 1
        fi(1) = x(1) - 1
    End Sub

    Public Sub function2_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates
        ' f0(x0,x1) = x0^2+1
        ' f1(x0,x1) = x1-1
        ' and Jacobian matrix J = [dfi/dxj]
        '
        fi(0) = x(0) * x(0) + 1
        fi(1) = x(1) - 1
        jac(0, 0) = 2 * x(0)
        jac(0, 1) = 0
        jac(1, 0) = 0
        jac(1, 1) = 1
    End Sub

    Public Sub nlcfunc1_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates
        '
        '     f0(x0,x1) = -x0+x1
        '     f1(x0,x1) = x0^2+x1^2-1
        '
        ' and Jacobian matrix J = [dfi/dxj]
        '
        fi(0) = -x(0)+x(1)
        fi(1) = x(0)*x(0) + x(1)*x(1) - 1.0
        jac(0,0) = -1.0
        jac(0,1) = +1.0
        jac(1,0) = 2*x(0)
        jac(1,1) = 2*x(1)
    End Sub

    Public Sub nlcfunc2_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates
        '
        '     f0(x0,x1,x2) = x0+x1
        '     f1(x0,x1,x2) = x2-exp(x0)
        '     f2(x0,x1,x2) = x0^2+x1^2-1
        '
        ' and Jacobian matrix J = [dfi/dxj]
        '
        fi(0) = x(0)+x(1)
        fi(1) = x(2)-System.Math.Exp(x(0))
        fi(2) = x(0)*x(0) + x(1)*x(1) - 1.0
        jac(0,0) = 1.0
        jac(0,1) = 1.0
        jac(0,2) = 0.0
        jac(1,0) = -System.Math.Exp(x(0))
        jac(1,1) = 0.0
        jac(1,2) = 1.0
        jac(2,0) = 2*x(0)
        jac(2,1) = 2*x(1)
        jac(2,2) = 0.0
    End Sub

    Public Sub nsfunc1_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates
        '
        '     f0(x0,x1) = 2*|x0|+|x1|
        '
        ' and Jacobian matrix J = [df0/dx0 df0/dx1]
        '
        fi(0) = 2*System.Math.Abs(x(0))+System.Math.Abs(x(1))
        jac(0,0) = 2*System.Math.Sign(x(0))
        jac(0,1) = System.Math.Sign(x(1))
    End Sub

    Public Sub nsfunc1_fvec(x As Double(), fi As Double(), obj As Object)
        '
        ' this callback calculates
        '
        '     f0(x0,x1) = 2*|x0|+|x1|
        '
        ' and Jacobian matrix J = [df0/dx0 df0/dx1]
        '
        fi(0) = 2*System.Math.Abs(x(0))+System.Math.Abs(x(1))
    End Sub

    Public Sub nsfunc2_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates function vector
        '
        '     f0(x0,x1) = 2*|x0|+x1
        '     f1(x0,x1) = x0-1
        '     f2(x0,x1) = -x1-1
        '
        ' and Jacobian matrix J
        '
        '         [ df0/dx0   df0/dx1 ]
        '     J = [ df1/dx0   df1/dx1 ]
        '         [ df2/dx0   df2/dx1 ]
        '
        fi(0) = 2*System.Math.Abs(x(0))+System.Math.Abs(x(1))
        jac(0,0) = 2*System.Math.Sign(x(0))
        jac(0,1) = System.Math.Sign(x(1))
        fi(1) = x(0)-1
        jac(1,0) = 1
        jac(1,1) = 0
        fi(2) = -x(1)-1
        jac(2,0) = 0
        jac(2,1) = -1
    End Sub

    Public Sub bad_func(x As Double(), ByRef func As Double, obj As Object)
        '
        ' this callback calculates 'bad' function,
        ' i.e. function with incorrectly calculated derivatives
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
    End Sub

    Public Sub bad_grad(x As Double(), ByRef func As Double, grad As Double(), obj As Object)
        '
        ' this callback calculates 'bad' function,
        ' i.e. function with incorrectly calculated derivatives
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
        grad(0) = 40 * System.Math.Pow(x(0) + 3, 3)
        grad(1) = 40 * System.Math.Pow(x(1) - 3, 3)
    End Sub

    Public Sub bad_hess(x As Double(), ByRef func As Double, grad As Double(), hess As Double(,), obj As Object)
        '
        ' this callback calculates 'bad' function,
        ' i.e. function with incorrectly calculated derivatives
        '
        func = 100 * System.Math.Pow(x(0) + 3, 4) + System.Math.Pow(x(1) - 3, 4)
        grad(0) = 40 * System.Math.Pow(x(0) + 3, 3)
        grad(1) = 40 * System.Math.Pow(x(1) - 3, 3)
        hess(0, 0) = 120 * System.Math.Pow(x(0) + 3, 2)
        hess(0, 1) = 1
        hess(1, 0) = 1
        hess(1, 1) = 120 * System.Math.Pow(x(1) - 3, 2)
    End Sub

    Public Sub bad_fvec(x As Double(), fi As Double(), obj As Object)
        '
        ' this callback calculates 'bad' function,
        ' i.e. function with incorrectly calculated derivatives
        '
        fi(0) = 10 * System.Math.Pow(x(0) + 3, 2)
        fi(1) = System.Math.Pow(x(1) - 3, 2)
    End Sub

    Public Sub bad_jac(x As Double(), fi As Double(), jac As Double(,), obj As Object)
        '
        ' this callback calculates 'bad' function,
        ' i.e. function with incorrectly calculated derivatives
        '
        fi(0) = 10 * System.Math.Pow(x(0) + 3, 2)
        fi(1) = System.Math.Pow(x(1) - 3, 2)
        jac(0, 0) = 20 * (x(0) + 3)
        jac(0, 1) = 0
        jac(1, 0) = 1
        jac(1, 1) = 20 * (x(1) - 3)
    End Sub

    Public Sub function_cx_1_func(c As Double(), x As Double(), ByRef func As Double, obj As Object)
        '
        ' this callback calculates f(c,x)=exp(-c0*sqr(x0))
        ' where x is a position on X-axis and c is adjustable parameter
        '
        func = System.Math.Exp(-c(0) * x(0) * x(0))
    End Sub

    Public Sub function_cx_1_grad(c As Double(), x As Double(), ByRef func As Double, grad As Double(), obj As Object)
        '
        ' this callback calculates f(c,x)=exp(-c0*sqr(x0)) and gradient G={df/dc[i]}
        ' where x is a position on X-axis and c is adjustable parameter.
        ' IMPORTANT: gradient is calculated with respect to C, not to X
        '
        func = System.Math.Exp(-c(0) * System.Math.Pow(x(0), 2))
        grad(0) = -System.Math.Pow(x(0), 2) * func
    End Sub

    Public Sub function_cx_1_hess(c As Double(), x As Double(), ByRef func As Double, grad As Double(), hess As Double(,), obj As Object)
        '
        ' this callback calculates f(c,x)=exp(-c0*sqr(x0)), gradient G={df/dc[i]} and Hessian H={d2f/(dc[i]*dc[j])}
        ' where x is a position on X-axis and c is adjustable parameter.
        ' IMPORTANT: gradient/Hessian are calculated with respect to C, not to X
        '
        func = System.Math.Exp(-c(0) * System.Math.Pow(x(0), 2))
        grad(0) = -System.Math.Pow(x(0), 2) * func
        hess(0, 0) = System.Math.Pow(x(0), 4) * func
    End Sub

    Public Sub ode_function_1_diff(y As Double(), x As Double, dy As Double(), obj As Object)
        '
        ' this callback calculates f(y[],x)=-y[0]
        '
        dy(0) = -y(0)
    End Sub

    Public Sub int_function_1_func(x As Double, xminusa As Double, bminusx As Double, ByRef y As Double, obj As Object)
        '
        ' this callback calculates f(x)=exp(x)
        '
        y = Math.Exp(x)
    End Sub

    Public Sub function_debt_func(c As Double(), x As Double(), ByRef func As Double, obj As Object)
        '
        ' this callback calculates f(c,x)=c[0]*(1+c[1]*(pow(x[0]-1999,c[2])-1))
        '
        func = c(0)*(1+c(1)*(System.Math.Pow(x(0)-1999,c(2))-1))
    End Sub

    Public Sub s1_grad(x As Double(), ByRef func As Double, grad As Double(), obj As Object)
        '
        ' this callback calculates f(x) = (1+x)^(-0.2) + (1-x)^(-0.3) + 1000*x and its gradient.
        '
        ' function is trimmed when we calculate it near the singular points or outside of the [-1,+1].
        ' Note that we do NOT calculate gradient in this case.
        '
        If x(0)<=-0.999999999999 Or x(0)>=+0.999999999999 Then
            func = 1.0E+300
            Return
        End If
        func = System.Math.Pow(1+x(0),-0.2) + System.Math.Pow(1-x(0),-0.3) + 1000*x(0)
        grad(0) = -0.2*System.Math.Pow(1+x(0),-1.2) +0.3*System.Math.Pow(1-x(0),-1.3) + 1000
    End Sub

    Sub Main()
        Dim _TotalResult As Boolean = True
        Dim _TestResult As Boolean
        Dim _spoil_scenario As Integer
		Dim v_spoil_bool As Boolean
		Dim v_spoil_int As Integer
		Dim v_spoil_real As Double
		Dim v_spoil_complex As alglib.Complex
        System.Console.WriteLine("VB.NET interface tests. Please wait...")
        Try
            '
            ' TEST nneighbor_d_1
            '      Nearest neighbor search, KNN queries
            '
            System.Console.WriteLine("0/151")
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    Dim a(,) As Double = New Double(,){{0,0},{0,1},{1,0},{1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    Dim nx As Integer = 2
                    Dim ny As Integer = 0
                    Dim normtype As Integer = 2
                    Dim kdt As kdtree = New XAlglib.kdtree() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){}
                    Dim r(,) As Double = New Double(,){{}}
                    Dim k As Integer
                    xalglib.kdtreebuild(a, nx, ny, normtype, kdt)
                    x = New Double(){-1,0}
                    k = xalglib.kdtreequeryknn(kdt, x, 1)
                    _TestResult = _TestResult And doc_test_int(k, 1)
                    xalglib.kdtreequeryresultsx(kdt, r)
                    _TestResult = _TestResult And doc_test_real_matrix(r, New Double(,){{0,0}}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nneighbor_t_2
            '      Subsequent queries; buffered functions must use previously allocated storage (if large enough), so buffer may contain some info from previous call
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    Dim a(,) As Double = New Double(,){{0,0},{0,1},{1,0},{1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    Dim nx As Integer = 2
                    Dim ny As Integer = 0
                    Dim normtype As Integer = 2
                    Dim kdt As kdtree = New XAlglib.kdtree() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){}
                    Dim rx(,) As Double = New Double(,){{}}
                    Dim k As Integer
                    xalglib.kdtreebuild(a, nx, ny, normtype, kdt)
                    x = New Double(){+2,0}
                    k = xalglib.kdtreequeryknn(kdt, x, 2, true)
                    _TestResult = _TestResult And doc_test_int(k, 2)
                    xalglib.kdtreequeryresultsx(kdt, rx)
                    _TestResult = _TestResult And doc_test_real_matrix(rx, New Double(,){{1,0},{1,1}}, 0.05)
                    x = New Double(){-2,0}
                    k = xalglib.kdtreequeryknn(kdt, x, 1, true)
                    _TestResult = _TestResult And doc_test_int(k, 1)
                    xalglib.kdtreequeryresultsx(kdt, rx)
                    _TestResult = _TestResult And doc_test_real_matrix(rx, New Double(,){{0,0},{1,1}}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_t_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nneighbor_d_2
            '      Serialization of KD-trees
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    Dim a(,) As Double = New Double(,){{0,0},{0,1},{1,0},{1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    Dim nx As Integer = 2
                    Dim ny As Integer = 0
                    Dim normtype As Integer = 2
                    Dim kdt0 As kdtree = New XAlglib.kdtree() ' initializer can be dropped, but compiler will issue warning
                    Dim kdt1 As kdtree = New XAlglib.kdtree() ' initializer can be dropped, but compiler will issue warning
                    Dim s As String = "" ' initializer can be dropped, but compiler will issue warning 
                    Dim x() As Double = New Double(){}
                    Dim r0(,) As Double = New Double(,){{}}
                    Dim r1(,) As Double = New Double(,){{}}

                    ' 
                    '  Build tree and serialize it
                    ' 
                    xalglib.kdtreebuild(a, nx, ny, normtype, kdt0)
                    xalglib.kdtreeserialize(kdt0, s)
                    xalglib.kdtreeunserialize(s, kdt1)

                    ' 
                    '  Compare results from KNN queries
                    ' 
                    x = New Double(){-1,0}
                    xalglib.kdtreequeryknn(kdt0, x, 1)
                    xalglib.kdtreequeryresultsx(kdt0, r0)
                    xalglib.kdtreequeryknn(kdt1, x, 1)
                    xalglib.kdtreequeryresultsx(kdt1, r1)
                    _TestResult = _TestResult And doc_test_real_matrix(r0, New Double(,){{0,0}}, 0.05)
                    _TestResult = _TestResult And doc_test_real_matrix(r1, New Double(,){{0,0}}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nneighbor_d_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST odesolver_d1
            '      Solving y'=-y with ODE solver
            '
            _TestResult = true
            For _spoil_scenario = -1 To 12
                Try
                    Dim y() As Double = New Double(){1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim x() As Double = New Double(){0,1,2,3}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim eps As Double = 0.00001
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        eps = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.PositiveInfinity
                        eps = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NegativeInfinity
                        eps = v_spoil_real
                    End If
                    Dim h As Double = 0
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        h = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        h = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        h = v_spoil_real
                    End If
                    Dim s As odesolverstate = New XAlglib.odesolverstate() ' initializer can be dropped, but compiler will issue warning
                    Dim m As Integer
                    Dim xtbl() As Double = New Double(){}
                    Dim ytbl(,) As Double = New Double(,){{}}
                    Dim rep As odesolverreport = New XAlglib.odesolverreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.odesolverrkck(y, x, eps, h, s)
                    xalglib.odesolversolve(s, AddressOf ode_function_1_diff, Nothing)
                    xalglib.odesolverresults(s, m, xtbl, ytbl, rep)
                    _TestResult = _TestResult And doc_test_int(m, 4)
                    _TestResult = _TestResult And doc_test_real_vector(xtbl, New Double(){0,1,2,3}, 0.005)
                    _TestResult = _TestResult And doc_test_real_matrix(ytbl, New Double(,){{1},{0.367},{0.135},{0.050}}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "odesolver_d1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST sparse_d_1
            '      Basic operations with sparse matrices
            '
            _TestResult = true
            For _spoil_scenario = -1 To 0
                Try
                    ' 
                    '  This example demonstrates creation/initialization of the sparse matrix
                    '  and matrix-vector multiplication.
                    ' 
                    '  First, we have to create matrix and initialize it. Matrix is initially created
                    '  in the Hash-Table format, which allows convenient initialization. We can modify
                    '  Hash-Table matrix with sparseset() and sparseadd() functions.
                    ' 
                    '  NOTE: Unlike CRS format, Hash-Table representation allows you to initialize
                    '  elements in the arbitrary order. You may see that we initialize a[0][0] first,
                    '  then move to the second row, and then move back to the first row.
                    ' 
                    Dim s As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    xalglib.sparsecreate(2, 2, s)
                    xalglib.sparseset(s, 0, 0, 2.0)
                    xalglib.sparseset(s, 1, 1, 1.0)
                    xalglib.sparseset(s, 0, 1, 1.0)

                    xalglib.sparseadd(s, 1, 1, 4.0)

                    ' 
                    '  Now S is equal to
                    '    [ 2 1 ]
                    '    [   5 ]
                    '  Lets check it by reading matrix contents with sparseget().
                    '  You may see that with sparseget() you may read both non-zero
                    '  and zero elements.
                    ' 
                    Dim v As Double
                    v = xalglib.sparseget(s, 0, 0)
                    _TestResult = _TestResult And doc_test_real(v, 2.0000, 0.005)
                    v = xalglib.sparseget(s, 0, 1)
                    _TestResult = _TestResult And doc_test_real(v, 1.0000, 0.005)
                    v = xalglib.sparseget(s, 1, 0)
                    _TestResult = _TestResult And doc_test_real(v, 0.0000, 0.005)
                    v = xalglib.sparseget(s, 1, 1)
                    _TestResult = _TestResult And doc_test_real(v, 5.0000, 0.005)

                    ' 
                    '  After successful creation we can use our matrix for linear operations.
                    ' 
                    '  However, there is one more thing we MUST do before using S in linear
                    '  operations: we have to convert it from HashTable representation (used for
                    '  initialization and dynamic operations) to CRS format with sparseconverttocrs()
                    '  call. If you omit this call, ALGLIB will generate exception on the first
                    '  attempt to use S in linear operations. 
                    ' 
                    xalglib.sparseconverttocrs(s)

                    ' 
                    '  Now S is in the CRS format and we are ready to do linear operations.
                    '  Lets calculate A*x for some x.
                    ' 
                    Dim x() As Double = New Double(){1,-1}
                    If _spoil_scenario=0 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){}
                    xalglib.sparsemv(s, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){1.000,-5.000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "sparse_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST sparse_d_crs
            '      Advanced topic: creation in the CRS format.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 1
                Try
                    ' 
                    '  This example demonstrates creation/initialization of the sparse matrix in the
                    '  CRS format.
                    ' 
                    '  Hash-Table format used by default is very convenient (it allows easy
                    '  insertion of elements, automatic memory reallocation), but has
                    '  significant memory and performance overhead. Insertion of one element 
                    '  costs hundreds of CPU cycles, and memory consumption is several times
                    '  higher than that of CRS.
                    ' 
                    '  When you work with really large matrices and when you can tell in 
                    '  advance how many elements EXACTLY you need, it can be beneficial to 
                    '  create matrix in the CRS format from the very beginning.
                    ' 
                    '  If you want to create matrix in the CRS format, you should:
                    '  * use sparsecreatecrs() function
                    '  * know row sizes in advance (number of non-zero entries in the each row)
                    '  * initialize matrix with sparseset() - another function, sparseadd(), is not allowed
                    '  * initialize elements from left to right, from top to bottom, each
                    '    element is initialized only once.
                    ' 
                    Dim s As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    Dim row_sizes() As Integer = New Integer(){2,2,2,1}
                    If _spoil_scenario=0 Then
                        spoil_vector_by_deleting_element(row_sizes)
                    End If
                    xalglib.sparsecreatecrs(4, 4, row_sizes, s)
                    xalglib.sparseset(s, 0, 0, 2.0)
                    xalglib.sparseset(s, 0, 1, 1.0)
                    xalglib.sparseset(s, 1, 1, 4.0)
                    xalglib.sparseset(s, 1, 2, 2.0)
                    xalglib.sparseset(s, 2, 2, 3.0)
                    xalglib.sparseset(s, 2, 3, 1.0)
                    xalglib.sparseset(s, 3, 3, 9.0)

                    ' 
                    '  Now S is equal to
                    '    [ 2 1     ]
                    '    [   4 2   ]
                    '    [     3 1 ]
                    '    [       9 ]
                    ' 
                    '  We should point that we have initialized S elements from left to right,
                    '  from top to bottom. CRS representation does NOT allow you to do so in
                    '  the different order. Try to change order of the sparseset() calls above,
                    '  and you will see that your program generates exception.
                    ' 
                    '  We can check it by reading matrix contents with sparseget().
                    '  However, you should remember that sparseget() is inefficient on
                    '  CRS matrices (it may have to pass through all elements of the row 
                    '  until it finds element you need).
                    ' 
                    Dim v As Double
                    v = xalglib.sparseget(s, 0, 0)
                    _TestResult = _TestResult And doc_test_real(v, 2.0000, 0.005)
                    v = xalglib.sparseget(s, 2, 3)
                    _TestResult = _TestResult And doc_test_real(v, 1.0000, 0.005)

                    '  you may see that you can read zero elements (which are not stored) with sparseget()
                    v = xalglib.sparseget(s, 3, 2)
                    _TestResult = _TestResult And doc_test_real(v, 0.0000, 0.005)

                    ' 
                    '  After successful creation we can use our matrix for linear operations.
                    '  Lets calculate A*x for some x.
                    ' 
                    Dim x() As Double = New Double(){1,-1,1,-1}
                    If _spoil_scenario=1 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){}
                    xalglib.sparsemv(s, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){1.000,-2.000,2.000,-9}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "sparse_d_crs")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ablas_d_gemm
            '      Matrix multiplication (single-threaded)
            '
            _TestResult = True
            try
                Dim a(,) As Double = New Double(,){{2,1},{1,3}}
                Dim b(,) As Double = New Double(,){{2,1},{0,1}}
                Dim c(,) As Double = New Double(,){{0,0},{0,0}}

                ' 
                '  rmatrixgemm() function allows us to calculate matrix product C:=A*B or
                '  to perform more general operation, C:=alpha*op1(A)*op2(B)+beta*C,
                '  where A, B, C are rectangular matrices, op(X) can be X or X^T,
                '  alpha and beta are scalars.
                ' 
                '  This function:
                '  * can apply transposition and/or multiplication by scalar to operands
                '  * can use arbitrary part of matrices A/B (given by submatrix offset)
                '  * can store result into arbitrary part of C
                '  * for performance reasons requires C to be preallocated
                ' 
                '  Parameters of this function are:
                '  * M, N, K            -   sizes of op1(A) (which is MxK), op2(B) (which
                '                           is KxN) and C (which is MxN)
                '  * Alpha              -   coefficient before A*B
                '  * A, IA, JA          -   matrix A and offset of the submatrix
                '  * OpTypeA            -   transformation type:
                '                           0 - no transformation
                '                           1 - transposition
                '  * B, IB, JB          -   matrix B and offset of the submatrix
                '  * OpTypeB            -   transformation type:
                '                           0 - no transformation
                '                           1 - transposition
                '  * Beta               -   coefficient before C
                '  * C, IC, JC          -   preallocated matrix C and offset of the submatrix
                ' 
                '  Below we perform simple product C:=A*B (alpha=1, beta=0)
                ' 
                '  IMPORTANT: this function works with preallocated C, which must be large
                '             enough to store multiplication result.
                ' 
                Dim m As Integer = 2
                Dim n As Integer = 2
                Dim k As Integer = 2
                Dim alpha As Double = 1.0
                Dim ia As Integer = 0
                Dim ja As Integer = 0
                Dim optypea As Integer = 0
                Dim ib As Integer = 0
                Dim jb As Integer = 0
                Dim optypeb As Integer = 0
                Dim beta As Double = 0.0
                Dim ic As Integer = 0
                Dim jc As Integer = 0
                xalglib.rmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc)
                _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{4,3},{2,4}}, 0.0001)

                ' 
                '  Now we try to apply some simple transformation to operands: C:=A*B^T
                ' 
                optypeb = 1
                xalglib.rmatrixgemm(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc)
                _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{5,1},{5,3}}, 0.0001)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ablas_d_gemm")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ablas_d_syrk
            '      Symmetric rank-K update (single-threaded)
            '
            _TestResult = True
            try
                ' 
                '  rmatrixsyrk() function allows us to calculate symmetric rank-K update
                '  C := beta*C + alpha*A'*A, where C is square N*N matrix, A is square K*N
                '  matrix, alpha and beta are scalars. It is also possible to update by
                '  adding A*A' instead of A'*A.
                ' 
                '  Parameters of this function are:
                '  * N, K       -   matrix size
                '  * Alpha      -   coefficient before A
                '  * A, IA, JA  -   matrix and submatrix offsets
                '  * OpTypeA    -   multiplication type:
                '                   * 0 - A*A^T is calculated
                '                   * 2 - A^T*A is calculated
                '  * Beta       -   coefficient before C
                '  * C, IC, JC  -   preallocated input/output matrix and submatrix offsets
                '  * IsUpper    -   whether upper or lower triangle of C is updated;
                '                   this function updates only one half of C, leaving
                '                   other half unchanged (not referenced at all).
                ' 
                '  Below we will show how to calculate simple product C:=A'*A
                ' 
                '  NOTE: beta=0 and we do not use previous value of C, but still it
                '        MUST be preallocated.
                ' 
                Dim n As Integer = 2
                Dim k As Integer = 1
                Dim alpha As Double = 1.0
                Dim ia As Integer = 0
                Dim ja As Integer = 0
                Dim optypea As Integer = 2
                Dim beta As Double = 0.0
                Dim ic As Integer = 0
                Dim jc As Integer = 0
                Dim isupper As Boolean = true
                Dim a(,) As Double = New Double(,){{1,2}}

                '  preallocate space to store result
                Dim c(,) As Double = New Double(,){{0,0},{0,0}}

                '  calculate product, store result into upper part of c
                xalglib.rmatrixsyrk(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper)

                '  output result.
                '  IMPORTANT: lower triangle of C was NOT updated!
                _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{1,2},{0,4}}, 0.0001)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ablas_d_syrk")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ablas_t_complex
            '      Basis test for complex matrix functions (correctness and presence of SMP support)
            '
            _TestResult = True
            try
                Dim a(,) As alglib.complex = New alglib.complex(,){{}}
                Dim b(,) As alglib.complex = New alglib.complex(,){{}}
                Dim c(,) As alglib.complex = New alglib.complex(,){{}}

                '  test cmatrixgemm()
                a = New alglib.complex(,){{New alglib.complex(0,2),New alglib.complex(0,1)},{1,3}}
                b = New alglib.complex(,){{2,1},{0,1}}
                c = New alglib.complex(,){{0,0},{0,0}}
                xalglib.cmatrixgemm(2, 2, 2, 1.0, a, 0, 0, 0, b, 0, 0, 0, 0.0, c, 0, 0)
                _TestResult = _TestResult And doc_test_complex_matrix(c, New alglib.complex(,){{New alglib.complex(0,4),New alglib.complex(0,3)},{2,4}}, 0.0001)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ablas_t_complex")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_d_r1
            '      Real matrix inverse
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim a(,) As Double = New Double(,){{1,-1},{1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_row(a, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_col(a, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim info As Integer
                    Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.rmatrixinverse(a, info, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_matrix(a, New Double(,){{0.5,0.5},{-0.5,0.5}}, 0.00005)
                    _TestResult = _TestResult And doc_test_real(rep.r1, 0.5, 0.00005)
                    _TestResult = _TestResult And doc_test_real(rep.rinf, 0.5, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_r1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_d_c1
            '      Complex matrix inverse
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim a(,) As alglib.complex = New alglib.complex(,){{New alglib.complex(0,1),-1},{New alglib.complex(0,1),1}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_row(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_col(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim info As Integer
                    Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.cmatrixinverse(a, info, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_complex_matrix(a, New alglib.complex(,){{New alglib.complex(0,-0.5),New alglib.complex(0,-0.5)},{-0.5,0.5}}, 0.00005)
                    _TestResult = _TestResult And doc_test_real(rep.r1, 0.5, 0.00005)
                    _TestResult = _TestResult And doc_test_real(rep.rinf, 0.5, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_c1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_d_spd1
            '      SPD matrix inverse
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim a(,) As Double = New Double(,){{2,1},{1,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_row(a, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_col(a, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim info As Integer
                    Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.spdmatrixinverse(a, info, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_matrix(a, New Double(,){{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_spd1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_d_hpd1
            '      HPD matrix inverse
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim a(,) As alglib.complex = New alglib.complex(,){{2,1},{1,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_row(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_col(a, v_spoil_complex)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim info As Integer
                    Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.hpdmatrixinverse(a, info, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_complex_matrix(a, New alglib.complex(,){{0.666666,-0.333333},{-0.333333,0.666666}}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_d_hpd1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_t_r1
            '      Real matrix inverse: singular matrix
            '
            _TestResult = True
            try
                Dim a(,) As Double = New Double(,){{1,-1},{-2,2}}
                Dim info As Integer
                Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                xalglib.rmatrixinverse(a, info, rep)
                _TestResult = _TestResult And doc_test_int(info, -3)
                _TestResult = _TestResult And doc_test_real(rep.r1, 0.0, 0.00005)
                _TestResult = _TestResult And doc_test_real(rep.rinf, 0.0, 0.00005)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_r1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_t_c1
            '      Complex matrix inverse: singular matrix
            '
            _TestResult = True
            try
                Dim a(,) As alglib.complex = New alglib.complex(,){{New alglib.complex(0,1),New alglib.complex(0,-1)},{-2,2}}
                Dim info As Integer
                Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                xalglib.cmatrixinverse(a, info, rep)
                _TestResult = _TestResult And doc_test_int(info, -3)
                _TestResult = _TestResult And doc_test_real(rep.r1, 0.0, 0.00005)
                _TestResult = _TestResult And doc_test_real(rep.rinf, 0.0, 0.00005)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_t_c1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_e_spd1
            '      Attempt to use SPD function on nonsymmetrix matrix
            '
            _TestResult = true
            Try
                Dim a(,) As Double = New Double(,){{1,0},{1,1}}
                Dim info As Integer
                Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                xalglib.spdmatrixinverse(a, info, rep)
                _TestResult = False
            Catch E As AlglibException
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_spd1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matinv_e_hpd1
            '      Attempt to use SPD function on nonsymmetrix matrix
            '
            _TestResult = true
            Try
                Dim a(,) As alglib.complex = New alglib.complex(,){{1,0},{1,1}}
                Dim info As Integer
                Dim rep As matinvreport = New XAlglib.matinvreport() ' initializer can be dropped, but compiler will issue warning
                xalglib.hpdmatrixinverse(a, info, rep)
                _TestResult = False
            Catch E As AlglibException
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matinv_e_hpd1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlbfgs_d_1
            '      Nonlinear optimization by L-BFGS
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  using LBFGS method, with:
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minlbfgssetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlbfgsstate = New XAlglib.minlbfgsstate() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlbfgscreate(1, x, state)
                    xalglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.minlbfgssetscale(state, s)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target function (C0 continuity violation)
                    '  * nonsmoothness of the target function (C1 continuity violation)
                    '  * erroneous analytic gradient, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minlbfgsoptguardsmoothness(state)
                    xalglib.minlbfgsoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and examine results.
                    ' 
                    Dim rep As minlbfgsreport = New XAlglib.minlbfgsreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlbfgsoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minlbfgsresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlbfgsoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlbfgs_d_2
            '      Nonlinear optimization with additional settings and restarts
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    '  using LBFGS method.
                    ' 
                    '  Several advanced techniques are demonstrated:
                    '  * upper limit on step size
                    '  * restart from new point
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim stpmax As Double = 0.1
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        stpmax = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        stpmax = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        stpmax = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlbfgsstate = New XAlglib.minlbfgsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlbfgsreport = New XAlglib.minlbfgsreport() ' initializer can be dropped, but compiler will issue warning

                    '  create and tune optimizer
                    xalglib.minlbfgscreate(1, x, state)
                    xalglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.minlbfgssetstpmax(state, stpmax)
                    xalglib.minlbfgssetscale(state, s)

                    '  Set up OptGuard integrity checker which catches errors
                    '  like nonsmooth targets or errors in the analytic gradient.
                    ' 
                    '  OptGuard is essential at the early prototyping stages.
                    ' 
                    '  NOTE: gradient verification needs 3*N additional function
                    '        evaluations; DO NOT USE IT IN THE PRODUCTION CODE
                    '        because it leads to unnecessary slowdown of your app.
                    xalglib.minlbfgsoptguardsmoothness(state)
                    xalglib.minlbfgsoptguardgradient(state, 0.001)

                    '  first run
                    xalglib.minlbfgsoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minlbfgsresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    '  second run - algorithm is restarted
                    x = New Double(){10,10}
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    xalglib.minlbfgsrestartfrom(state, x)
                    xalglib.minlbfgsoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minlbfgsresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    '  check OptGuard integrity report. Why do we need it at all?
                    '  Well, try breaking the gradient by adding 1.0 to some
                    '  of its components - OptGuard should report it as error.
                    '  And it may also catch unintended errors too :)
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlbfgsoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_d_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlbfgs_numdiff
            '      Nonlinear optimization by L-BFGS with numerical differentiation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    '  using numerical differentiation to calculate gradient.
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim epsg As Double = 0.0000000001
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim diffstep As Double = 1.0e-6
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlbfgsstate = New XAlglib.minlbfgsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlbfgsreport = New XAlglib.minlbfgsreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.minlbfgscreatef(1, x, diffstep, state)
                    xalglib.minlbfgssetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.minlbfgsoptimize(state, AddressOf function1_func, Nothing, Nothing)
                    xalglib.minlbfgsresults(state, x, rep)

                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlbfgs_numdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST linlsqr_d_1
            '      Solution of sparse linear systems with CG
            '
            _TestResult = true
            For _spoil_scenario = -1 To 3
                Try
                    ' 
                    '  This example illustrates solution of sparse linear least squares problem
                    '  with LSQR algorithm.
                    '  
                    '  Suppose that we have least squares problem min|A*x-b| with sparse A
                    '  represented by sparsematrix object
                    '          [ 1 1 ]
                    '          [ 1 1 ]
                    '      A = [ 2 1 ]
                    '          [ 1   ]
                    '          [   1 ]
                    '  and right part b
                    '      [ 4 ]
                    '      [ 2 ]
                    '  b = [ 4 ]
                    '      [ 1 ]
                    '      [ 2 ]
                    '  and we want to solve this system in the least squares sense using
                    '  LSQR algorithm. In order to do so, we have to create left part
                    '  (sparsematrix object) and right part (dense array).
                    ' 
                    '  Initially, sparse matrix is created in the Hash-Table format,
                    '  which allows easy initialization, but do not allow matrix to be
                    '  used in the linear solvers. So after construction you should convert
                    '  sparse matrix to CRS format (one suited for linear operations).
                    ' 
                    Dim a As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    xalglib.sparsecreate(5, 2, a)
                    xalglib.sparseset(a, 0, 0, 1.0)
                    xalglib.sparseset(a, 0, 1, 1.0)
                    xalglib.sparseset(a, 1, 0, 1.0)
                    xalglib.sparseset(a, 1, 1, 1.0)
                    xalglib.sparseset(a, 2, 0, 2.0)
                    xalglib.sparseset(a, 2, 1, 1.0)
                    xalglib.sparseset(a, 3, 0, 1.0)
                    xalglib.sparseset(a, 4, 1, 1.0)

                    ' 
                    '  Now our matrix is fully initialized, but we have to do one more
                    '  step - convert it from Hash-Table format to CRS format (see
                    '  documentation on sparse matrices for more information about these
                    '  formats).
                    ' 
                    '  If you omit this call, ALGLIB will generate exception on the first
                    '  attempt to use A in linear operations. 
                    ' 
                    xalglib.sparseconverttocrs(a)

                    ' 
                    '  Initialization of the right part
                    ' 
                    Dim b() As Double = New Double(){4,2,4,1,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(b)
                    End If

                    ' 
                    '  Now we have to create linear solver object and to use it for the
                    '  solution of the linear system.
                    ' 
                    Dim s As linlsqrstate = New XAlglib.linlsqrstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As linlsqrreport = New XAlglib.linlsqrreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){}
                    xalglib.linlsqrcreate(5, 2, s)
                    xalglib.linlsqrsolvesparse(s, a, b)
                    xalglib.linlsqrresults(s, x, rep)

                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){1.000,2.000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "linlsqr_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minbleic_d_1
            '      Nonlinear optimization with bound constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 19
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  subject to box constraints
                    ' 
                    '      -1<=x<=+1, -1<=y<=+1
                    ' 
                    '  using BLEIC optimizer with:
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minbleicsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties:
                    '  * set box constraints
                    '  * set variable scales
                    '  * set stopping criteria
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minbleicstate = New XAlglib.minbleicstate() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleiccreate(x, state)
                    xalglib.minbleicsetbc(state, bndl, bndu)
                    xalglib.minbleicsetscale(state, s)
                    xalglib.minbleicsetcond(state, epsg, epsf, epsx, maxits)

                    ' 
                    '  Then we activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target function (C0 continuity violation)
                    '  * nonsmoothness of the target function (C1 continuity violation)
                    '  * erroneous analytic gradient, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minbleicoptguardsmoothness(state)
                    xalglib.minbleicoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As minbleicreport = New XAlglib.minbleicreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minbleicresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-1,1}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minbleic_d_2
            '      Nonlinear optimization with linear inequality constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 21
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  subject to inequality constraints
                    ' 
                    '  * x>=2 (posed as general linear constraint),
                    '  * x+y>=6
                    ' 
                    '  using BLEIC optimizer with
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minbleicsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties:
                    '  * set linear constraints
                    '  * set variable scales
                    '  * set stopping criteria
                    ' 
                    Dim x() As Double = New Double(){5,5}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim c(,) As Double = New Double(,){{1,0,2},{1,1,6}}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_matrix_by_deleting_row(c)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_matrix_by_deleting_col(c)
                    End If
                    Dim ct() As Integer = New Integer(){1,1}
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(ct)
                    End If
                    Dim state As minbleicstate = New XAlglib.minbleicstate() ' initializer can be dropped, but compiler will issue warning
                    Dim epsg As Double = 0
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0

                    xalglib.minbleiccreate(x, state)
                    xalglib.minbleicsetlc(state, c, ct)
                    xalglib.minbleicsetscale(state, s)
                    xalglib.minbleicsetcond(state, epsg, epsf, epsx, maxits)

                    ' 
                    '  Then we activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target function (C0 continuity violation)
                    '  * nonsmoothness of the target function (C1 continuity violation)
                    '  * erroneous analytic gradient, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minbleicoptguardsmoothness(state)
                    xalglib.minbleicoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As minbleicreport = New XAlglib.minbleicreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minbleicresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2,4}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_d_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minbleic_numdiff
            '      Nonlinear optimization with bound constraints and numerical differentiation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 22
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  subject to box constraints
                    ' 
                    '      -1<=x<=+1, -1<=y<=+1
                    ' 
                    '  using BLEIC optimizer with:
                    '  * numerical differentiation being used
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minbleicsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties:
                    '  * set box constraints
                    '  * set variable scales
                    '  * set stopping criteria
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim state As minbleicstate = New XAlglib.minbleicstate() ' initializer can be dropped, but compiler will issue warning
                    Dim epsg As Double = 0
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim diffstep As Double = 1.0e-6
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If

                    xalglib.minbleiccreatef(x, diffstep, state)
                    xalglib.minbleicsetbc(state, bndl, bndu)
                    xalglib.minbleicsetscale(state, s)
                    xalglib.minbleicsetcond(state, epsg, epsf, epsx, maxits)

                    ' 
                    '  Then we activate OptGuard integrity checking.
                    ' 
                    '  Numerical differentiation always produces "correct" gradient
                    '  (with some truncation error, but unbiased). Thus, we just have
                    '  to check smoothness properties of the target: C0 and C1 continuity.
                    ' 
                    '  Sometimes user accidentally tries to solve nonsmooth problems
                    '  with smooth optimizer. OptGuard helps to detect such situations
                    '  early, at the prototyping stage.
                    ' 
                    xalglib.minbleicoptguardsmoothness(state)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As minbleicreport = New XAlglib.minbleicreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptimize(state, AddressOf function1_func, Nothing, Nothing)
                    xalglib.minbleicresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-1,1}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  Want to challenge OptGuard? Try to make your problem
                    '  nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    '  re-run optimizer.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbleicoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minbleic_numdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minqp_d_u1
            '      Unconstrained dense quadratic programming
            '
            _TestResult = true
            For _spoil_scenario = -1 To 16
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    ' 
                    '  Exact solution is [x0,x1] = [3,2]
                    ' 
                    '  We provide algorithm with starting point, although in this case
                    '  (dense matrix, no constraints) it can work without such information.
                    ' 
                    '  Several QP solvers are tried: QuickQP, BLEIC, DENSE-AUL.
                    ' 
                    '  IMPORTANT: this solver minimizes  following  function:
                    '      f(x) = 0.5*x'*A*x + b'*x.
                    '  Note that quadratic term has 0.5 before it. So if you want to minimize
                    '  quadratic function, you should rewrite it in such way that quadratic term
                    '  is multiplied by 0.5 too.
                    ' 
                    '  For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    '      f(x) = 0.5*(2*x0^2+2*x1^2) + .... 
                    '  and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    ' 
                    Dim a(,) As Double = New Double(,){{2,0},{0,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim b() As Double = New Double(){-6,-4}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(b)
                    End If
                    Dim x0() As Double = New Double(){0,1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(x0)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim x() As Double = New Double(){}
                    Dim state As minqpstate = New XAlglib.minqpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minqpreport = New XAlglib.minqpreport() ' initializer can be dropped, but compiler will issue warning

                    '  create solver, set quadratic/linear terms
                    xalglib.minqpcreate(2, state)
                    xalglib.minqpsetquadraticterm(state, a)
                    xalglib.minqpsetlinearterm(state, b)
                    xalglib.minqpsetstartingpoint(state, x0)

                    '  Set scale of the parameters.
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    '  NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    '        which automatically determines variable scales.
                    xalglib.minqpsetscale(state, s)

                    ' 
                    '  Solve problem with QuickQP solver.
                    ' 
                    '  This solver is intended for medium and large-scale problems with box
                    '  constraints (general linear constraints are not supported), but it can
                    '  also be efficiently used on unconstrained problems.
                    ' 
                    '  Default stopping criteria are used, Newton phase is active.
                    ' 
                    xalglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){3,2}, 0.005)

                    ' 
                    '  Solve problem with BLEIC-based QP solver.
                    ' 
                    '  This solver is intended for problems with moderate (up to 50) number
                    '  of general linear constraints and unlimited number of box constraints.
                    '  Of course, unconstrained problems can be solved too.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){3,2}, 0.005)

                    ' 
                    '  Solve problem with DENSE-AUL solver.
                    ' 
                    '  This solver is optimized for problems with up to several thousands of
                    '  variables and large amount of general linear constraints. Problems with
                    '  less than 50 general linear constraints can be efficiently solved with
                    '  BLEIC, problems with box-only constraints can be solved with QuickQP.
                    '  However, DENSE-AUL will work in any (including unconstrained) case.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){3,2}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_u1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minqp_d_bc1
            '      Bound constrained dense quadratic programming
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    '  subject to bound constraints 0<=x0<=2.5, 0<=x1<=2.5
                    ' 
                    '  Exact solution is [x0,x1] = [2.5,2]
                    ' 
                    '  We provide algorithm with starting point. With such small problem good starting
                    '  point is not really necessary, but with high-dimensional problem it can save us
                    '  a lot of time.
                    ' 
                    '  Several QP solvers are tried: QuickQP, BLEIC, DENSE-AUL.
                    ' 
                    '  IMPORTANT: this solver minimizes  following  function:
                    '      f(x) = 0.5*x'*A*x + b'*x.
                    '  Note that quadratic term has 0.5 before it. So if you want to minimize
                    '  quadratic function, you should rewrite it in such way that quadratic term
                    '  is multiplied by 0.5 too.
                    '  For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    '      f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    '  and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    ' 
                    Dim a(,) As Double = New Double(,){{2,0},{0,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim b() As Double = New Double(){-6,-4}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(b)
                    End If
                    Dim x0() As Double = New Double(){0,1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(x0)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){0.0,0.0}
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){2.5,2.5}
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim x() As Double = New Double(){}
                    Dim state As minqpstate = New XAlglib.minqpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minqpreport = New XAlglib.minqpreport() ' initializer can be dropped, but compiler will issue warning

                    '  create solver, set quadratic/linear terms
                    xalglib.minqpcreate(2, state)
                    xalglib.minqpsetquadraticterm(state, a)
                    xalglib.minqpsetlinearterm(state, b)
                    xalglib.minqpsetstartingpoint(state, x0)
                    xalglib.minqpsetbc(state, bndl, bndu)

                    '  Set scale of the parameters.
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    '  NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    '        which automatically determines variable scales.
                    xalglib.minqpsetscale(state, s)

                    ' 
                    '  Solve problem with QuickQP solver.
                    ' 
                    '  This solver is intended for medium and large-scale problems with box
                    '  constraints (general linear constraints are not supported).
                    ' 
                    '  Default stopping criteria are used, Newton phase is active.
                    ' 
                    xalglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 4)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2.5,2}, 0.005)

                    ' 
                    '  Solve problem with BLEIC-based QP solver.
                    ' 
                    '  This solver is intended for problems with moderate (up to 50) number
                    '  of general linear constraints and unlimited number of box constraints.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2.5,2}, 0.005)

                    ' 
                    '  Solve problem with DENSE-AUL solver.
                    ' 
                    '  This solver is optimized for problems with up to several thousands of
                    '  variables and large amount of general linear constraints. Problems with
                    '  less than 50 general linear constraints can be efficiently solved with
                    '  BLEIC, problems with box-only constraints can be solved with QuickQP.
                    '  However, DENSE-AUL will work in any (including unconstrained) case.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2.5,2}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_bc1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minqp_d_lc1
            '      Linearly constrained dense quadratic programming
            '
            _TestResult = true
            For _spoil_scenario = -1 To 15
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1
                    '  subject to linear constraint x0+x1<=2
                    ' 
                    '  Exact solution is [x0,x1] = [1.5,0.5]
                    ' 
                    '  IMPORTANT: this solver minimizes  following  function:
                    '      f(x) = 0.5*x'*A*x + b'*x.
                    '  Note that quadratic term has 0.5 before it. So if you want to minimize
                    '  quadratic function, you should rewrite it in such way that quadratic term
                    '  is multiplied by 0.5 too.
                    '  For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    '      f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    '  and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    ' 
                    Dim a(,) As Double = New Double(,){{2,0},{0,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim b() As Double = New Double(){-6,-4}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(b)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim c(,) As Double = New Double(,){{1.0,1.0,2.0}}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(c, v_spoil_real)
                    End If
                    Dim ct() As Integer = New Integer(){-1}
                    Dim x() As Double = New Double(){}
                    Dim state As minqpstate = New XAlglib.minqpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minqpreport = New XAlglib.minqpreport() ' initializer can be dropped, but compiler will issue warning

                    '  create solver, set quadratic/linear terms
                    xalglib.minqpcreate(2, state)
                    xalglib.minqpsetquadraticterm(state, a)
                    xalglib.minqpsetlinearterm(state, b)
                    xalglib.minqpsetlc(state, c, ct)

                    '  Set scale of the parameters.
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    '  NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    '        which automatically determines variable scales.
                    xalglib.minqpsetscale(state, s)

                    ' 
                    '  Solve problem with BLEIC-based QP solver.
                    ' 
                    '  This solver is intended for problems with moderate (up to 50) number
                    '  of general linear constraints and unlimited number of box constraints.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){1.500,0.500}, 0.05)

                    ' 
                    '  Solve problem with DENSE-AUL solver.
                    ' 
                    '  This solver is optimized for problems with up to several thousands of
                    '  variables and large amount of general linear constraints. Problems with
                    '  less than 50 general linear constraints can be efficiently solved with
                    '  BLEIC, problems with box-only constraints can be solved with QuickQP.
                    '  However, DENSE-AUL will work in any (including unconstrained) case.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){1.500,0.500}, 0.05)

                    ' 
                    '  Solve problem with QuickQP solver.
                    ' 
                    '  This solver is intended for medium and large-scale problems with box
                    '  constraints, and...
                    ' 
                    '  ...Oops! It does not support general linear constraints, -5 returned as completion code!
                    ' 
                    xalglib.minqpsetalgoquickqp(state, 0.0, 0.0, 0.0, 0, true)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, -5)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_lc1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minqp_d_u2
            '      Unconstrained sparse quadratic programming
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = x0^2 + x1^2 -6*x0 - 4*x1,
                    '  with quadratic term given by sparse matrix structure.
                    ' 
                    '  Exact solution is [x0,x1] = [3,2]
                    ' 
                    '  We provide algorithm with starting point, although in this case
                    '  (dense matrix, no constraints) it can work without such information.
                    ' 
                    '  IMPORTANT: this solver minimizes  following  function:
                    '      f(x) = 0.5*x'*A*x + b'*x.
                    '  Note that quadratic term has 0.5 before it. So if you want to minimize
                    '  quadratic function, you should rewrite it in such way that quadratic term
                    '  is multiplied by 0.5 too.
                    ' 
                    '  For example, our function is f(x)=x0^2+x1^2+..., but we rewrite it as 
                    '      f(x) = 0.5*(2*x0^2+2*x1^2) + ....
                    '  and pass diag(2,2) as quadratic term - NOT diag(1,1)!
                    ' 
                    Dim a As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    Dim b() As Double = New Double(){-6,-4}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(b)
                    End If
                    Dim x0() As Double = New Double(){0,1}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(x0)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim x() As Double = New Double(){}
                    Dim state As minqpstate = New XAlglib.minqpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minqpreport = New XAlglib.minqpreport() ' initializer can be dropped, but compiler will issue warning

                    '  initialize sparsematrix structure
                    xalglib.sparsecreate(2, 2, 0, a)
                    xalglib.sparseset(a, 0, 0, 2.0)
                    xalglib.sparseset(a, 1, 1, 2.0)

                    '  create solver, set quadratic/linear terms
                    xalglib.minqpcreate(2, state)
                    xalglib.minqpsetquadratictermsparse(state, a, true)
                    xalglib.minqpsetlinearterm(state, b)
                    xalglib.minqpsetstartingpoint(state, x0)

                    '  Set scale of the parameters.
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    '  NOTE: for convex problems you may try using minqpsetscaleautodiag()
                    '        which automatically determines variable scales.
                    xalglib.minqpsetscale(state, s)

                    ' 
                    '  Solve problem with BLEIC-based QP solver.
                    ' 
                    '  This solver is intended for problems with moderate (up to 50) number
                    '  of general linear constraints and unlimited number of box constraints.
                    '  It also supports sparse problems.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){3,2}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_u2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minqp_d_nonconvex
            '      Nonconvex quadratic programming
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  This example demonstrates minimization of nonconvex function
                    '      F(x0,x1) = -(x0^2+x1^2)
                    '  subject to constraints x0,x1 in [1.0,2.0]
                    '  Exact solution is [x0,x1] = [2,2].
                    ' 
                    '  Non-convex problems are harded to solve than convex ones, and they
                    '  may have more than one local minimum. However, ALGLIB solves may deal
                    '  with such problems (altough they do not guarantee convergence to
                    '  global minimum).
                    ' 
                    '  IMPORTANT: this solver minimizes  following  function:
                    '      f(x) = 0.5*x'*A*x + b'*x.
                    '  Note that quadratic term has 0.5 before it. So if you want to minimize
                    '  quadratic function, you should rewrite it in such way that quadratic term
                    '  is multiplied by 0.5 too.
                    ' 
                    '  For example, our function is f(x)=-(x0^2+x1^2), but we rewrite it as 
                    '      f(x) = 0.5*(-2*x0^2-2*x1^2)
                    '  and pass diag(-2,-2) as quadratic term - NOT diag(-1,-1)!
                    ' 
                    Dim a(,) As Double = New Double(,){{-2,0},{0,-2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim x0() As Double = New Double(){1,1}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(x0)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){1.0,1.0}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){2.0,2.0}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim x() As Double = New Double(){}
                    Dim state As minqpstate = New XAlglib.minqpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minqpreport = New XAlglib.minqpreport() ' initializer can be dropped, but compiler will issue warning

                    '  create solver, set quadratic/linear terms, constraints
                    xalglib.minqpcreate(2, state)
                    xalglib.minqpsetquadraticterm(state, a)
                    xalglib.minqpsetstartingpoint(state, x0)
                    xalglib.minqpsetbc(state, bndl, bndu)

                    '  Set scale of the parameters.
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    '  NOTE: there also exists minqpsetscaleautodiag() function
                    '        which automatically determines variable scales; however,
                    '        it does NOT work for non-convex problems.
                    xalglib.minqpsetscale(state, s)

                    ' 
                    '  Solve problem with BLEIC-based QP solver.
                    ' 
                    '  This solver is intended for problems with moderate (up to 50) number
                    '  of general linear constraints and unlimited number of box constraints.
                    ' 
                    '  It may solve non-convex problems as long as they are bounded from
                    '  below under constraints.
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2,2}, 0.005)

                    ' 
                    '  Solve problem with DENSE-AUL solver.
                    ' 
                    '  This solver is optimized for problems with up to several thousands of
                    '  variables and large amount of general linear constraints. Problems with
                    '  less than 50 general linear constraints can be efficiently solved with
                    '  BLEIC, problems with box-only constraints can be solved with QuickQP.
                    '  However, DENSE-AUL will work in any (including unconstrained) case.
                    ' 
                    '  Algorithm convergence is guaranteed only for convex case, but you may
                    '  expect that it will work for non-convex problems too (because near the
                    '  solution they are locally convex).
                    ' 
                    '  Default stopping criteria are used.
                    ' 
                    xalglib.minqpsetalgodenseaul(state, 1.0e-9, 1.0e+4, 5)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){2,2}, 0.005)

                    '  Hmm... this problem is bounded from below (has solution) only under constraints.
                    '  What it we remove them?
                    ' 
                    '  You may see that BLEIC algorithm detects unboundedness of the problem, 
                    '  -4 is returned as completion code. However, DENSE-AUL is unable to detect
                    '  such situation and it will cycle forever (we do not test it here).
                    Dim nobndl() As Double = New Double(){-Double.PositiveInfinity,-Double.PositiveInfinity}
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(nobndl, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        spoil_vector_by_deleting_element(nobndl)
                    End If
                    Dim nobndu() As Double = New Double(){Double.PositiveInfinity,+Double.PositiveInfinity}
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(nobndu, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        spoil_vector_by_deleting_element(nobndu)
                    End If
                    xalglib.minqpsetbc(state, nobndl, nobndu)
                    xalglib.minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0)
                    xalglib.minqpoptimize(state)
                    xalglib.minqpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, -4)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minqp_d_nonconvex")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlp_basic
            '      Basic linear programming example
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates how to minimize
                    ' 
                    '      F(x0,x1) = -0.1*x0 - x1
                    ' 
                    '  subject to box constraints
                    ' 
                    '      -1 <= x0,x1 <= +1 
                    ' 
                    '  and general linear constraints
                    ' 
                    '      x0 - x1 >= -1
                    '      x0 + x1 <=  1
                    ' 
                    '  We use dual simplex solver provided by ALGLIB for this task. Box
                    '  constraints are specified by means of constraint vectors bndl and
                    '  bndu (we have bndl<=x<=bndu). General linear constraints are
                    '  specified as AL<=A*x<=AU, with AL/AU being 2x1 vectors and A being
                    '  2x2 matrix.
                    ' 
                    '  NOTE: some/all components of AL/AU can be +-INF, same applies to
                    '        bndl/bndu. You can also have AL[I]=AU[i] (as well as
                    '        BndL[i]=BndU[i]).
                    ' 
                    Dim a(,) As Double = New Double(,){{1,-1},{1,+1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        spoil_matrix_by_deleting_row(a)
                    End If
                    If _spoil_scenario=2 Then
                        spoil_matrix_by_deleting_col(a)
                    End If
                    Dim al() As Double = New Double(){-1,-Double.PositiveInfinity}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(al, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(al)
                    End If
                    Dim au() As Double = New Double(){Double.PositiveInfinity,+1}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(au, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(au)
                    End If
                    Dim c() As Double = New Double(){-0.1,-1}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(c)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim x() As Double = New Double(){}
                    Dim state As minlpstate = New XAlglib.minlpstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlpreport = New XAlglib.minlpreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.minlpcreate(2, state)

                    ' 
                    '  Set cost vector, box constraints, general linear constraints.
                    ' 
                    '  Box constraints can be set in one call to minlpsetbc() or minlpsetbcall()
                    '  (latter sets same constraints for all variables and accepts two scalars
                    '  instead of two vectors).
                    ' 
                    '  General linear constraints can be specified in several ways:
                    '  * minlpsetlc2dense() - accepts dense 2D array as input; sometimes this
                    '    approach is more convenient, although less memory-efficient.
                    '  * minlpsetlc2() - accepts sparse matrix as input
                    '  * minlpaddlc2dense() - appends one row to the current set of constraints;
                    '    row being appended is specified as dense vector
                    '  * minlpaddlc2() - appends one row to the current set of constraints;
                    '    row being appended is specified as sparse set of elements
                    '  Independently from specific function being used, LP solver uses sparse
                    '  storage format for internal representation of constraints.
                    ' 
                    xalglib.minlpsetcost(state, c)
                    xalglib.minlpsetbc(state, bndl, bndu)
                    xalglib.minlpsetlc2dense(state, a, al, au, 2)

                    ' 
                    '  Set scale of the parameters.
                    ' 
                    '  It is strongly recommended that you set scale of your variables.
                    '  Knowing their scales is essential for evaluation of stopping criteria
                    '  and for preconditioning of the algorithm steps.
                    '  You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
                    ' 
                    xalglib.minlpsetscale(state, s)

                    '  Solve
                    xalglib.minlpoptimize(state)
                    xalglib.minlpresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){0,1}, 0.0005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlp_basic")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minnlc_d_inequality
            '      Nonlinearly constrained optimization (inequality constraints)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = -x0+x1
                    ' 
                    '  subject to box constraints
                    ' 
                    '     x0>=0, x1>=0
                    ' 
                    '  and nonlinear inequality constraint
                    ' 
                    '     x0^2 + x1^2 - 1 <= 0
                    ' 
                    Dim x0() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim bndl() As Double = New Double(){0,0}
                    Dim bndu() As Double = New Double(){Double.PositiveInfinity,+Double.PositiveInfinity}
                    Dim state As minnlcstate = New XAlglib.minnlcstate() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Create optimizer object and tune its settings:
                    '  * epsx=0.000001  stopping condition for inner iterations
                    '  * s=[1,1]        all variables have unit scale; it is important to
                    '                   tell optimizer about scales of your variables - it
                    '                   greatly accelerates convergence and helps to perform
                    '                   some important integrity checks.
                    ' 
                    xalglib.minnlccreate(2, x0, state)
                    xalglib.minnlcsetcond(state, epsx, maxits)
                    xalglib.minnlcsetscale(state, s)

                    ' 
                    '  Choose one of the nonlinear programming solvers supported by minnlc
                    '  optimizer:
                    '  * SQP - sequential quadratic programming NLP solver
                    '  * AUL - augmented Lagrangian NLP solver
                    '  * SLP - successive linear programming NLP solver
                    ' 
                    '  Different solvers have different properties:
                    '  * SQP needs less function evaluations than any other solver, but it
                    '    has much higher iteration cost than other solvers (a QP subproblem
                    '    has to be solved during each step)
                    '  * AUL solver has cheaper iterations, but needs more target function
                    '    evaluations
                    '  * SLP is the most robust solver provided by ALGLIB, but it performs
                    '    order of magnitude more iterations than SQP.
                    ' 
                    '  In the code below we set solver to be AUL but then override it with SLP,
                    '  and then with SQP, so the effective choice is to use SLP. We recommend
                    '  you to use SQP at least for early prototyping stages, and then switch
                    '  to AUL if possible.
                    ' 
                    Dim rho As Double = 1000.0
                    Dim outerits As Integer = 5
                    xalglib.minnlcsetalgoaul(state, rho, outerits)
                    xalglib.minnlcsetalgoslp(state)
                    xalglib.minnlcsetalgosqp(state)

                    ' 
                    '  Set constraints:
                    ' 
                    '  1. boundary constraints are passed with minnlcsetbc() call
                    ' 
                    '  2. nonlinear constraints are more tricky - you can not "pack" general
                    '     nonlinear function into double precision array. That's why
                    '     minnlcsetnlc() does not accept constraints itself - only constraint
                    '     counts are passed: first parameter is number of equality constraints,
                    '     second one is number of inequality constraints.
                    ' 
                    '     As for constraining functions - these functions are passed as part
                    '     of problem Jacobian (see below).
                    ' 
                    '  NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    '        linear and general nonlinear constraints. This example does not
                    '        show how to work with general linear constraints, but you can
                    '        easily find it in documentation on minnlcsetlc() function.
                    ' 
                    xalglib.minnlcsetbc(state, bndl, bndu)
                    xalglib.minnlcsetnlc(state, 0, 1)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target/constraints (C0 continuity violation)
                    '  * nonsmoothness of the target/constraints (C1 continuity violation)
                    '  * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minnlcoptguardsmoothness(state)
                    xalglib.minnlcoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [  -1    +1  ]
                    '      J = [            ]
                    '          [ 2*x0  2*x1 ]
                    ' 
                    '  with f0 being target function, f1 being constraining function. Number
                    '  of equality/inequality constraints is specified by minnlcsetnlc(),
                    '  with equality ones always being first, inequality ones being last.
                    ' 
                    Dim rep As minnlcreport = New XAlglib.minnlcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}
                    xalglib.minnlcoptimize(state, AddressOf nlcfunc1_jac, Nothing, Nothing)
                    xalglib.minnlcresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){1.0000,0.0000}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minnlcoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_inequality")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minnlc_d_equality
            '      Nonlinearly constrained optimization (equality constraints)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = -x0+x1
                    ' 
                    '  subject to nonlinear equality constraint
                    ' 
                    '     x0^2 + x1^2 - 1 = 0
                    ' 
                    Dim x0() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnlcstate = New XAlglib.minnlcstate() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Create optimizer object and tune its settings:
                    '  * epsx=0.000001  stopping condition for inner iterations
                    '  * s=[1,1]        all variables have unit scale
                    ' 
                    xalglib.minnlccreate(2, x0, state)
                    xalglib.minnlcsetcond(state, epsx, maxits)
                    xalglib.minnlcsetscale(state, s)

                    ' 
                    '  Choose one of the nonlinear programming solvers supported by minnlc
                    '  optimizer:
                    '  * SLP - successive linear programming NLP solver
                    '  * AUL - augmented Lagrangian NLP solver
                    ' 
                    '  Different solvers have different properties:
                    '  * SLP is the most robust solver provided by ALGLIB: it can solve both
                    '    convex and nonconvex optimization problems, it respects box and
                    '    linear constraints (after you find feasible point it won't move away
                    '    from the feasible area) and tries to respect nonlinear constraints
                    '    as much as possible. It also usually needs less function evaluations
                    '    to converge than AUL.
                    '    However, it solves LP subproblems at each iterations which adds
                    '    significant overhead to its running time. Sometimes it can be as much
                    '    as 7x times slower than AUL.
                    '  * AUL solver is less robust than SLP - it can violate box and linear
                    '    constraints at any moment, and it is intended for convex optimization
                    '    problems (although in many cases it can deal with nonconvex ones too).
                    '    Also, unlike SLP it needs some tuning (penalty factor and number of
                    '    outer iterations).
                    '    However, it is often much faster than the current version of SLP.
                    ' 
                    '  In the code below we set solver to be AUL but then override it with SLP,
                    '  so the effective choice is to use SLP. We recommend you to use SLP at
                    '  least for early prototyping stages.
                    ' 
                    '  You can comment out line with SLP if you want to solve your problem with
                    '  AUL solver.
                    ' 
                    Dim rho As Double = 1000.0
                    Dim outerits As Integer = 5
                    xalglib.minnlcsetalgoaul(state, rho, outerits)
                    xalglib.minnlcsetalgoslp(state)

                    ' 
                    '  Set constraints:
                    ' 
                    '  Nonlinear constraints are tricky - you can not "pack" general
                    '  nonlinear function into double precision array. That's why
                    '  minnlcsetnlc() does not accept constraints itself - only constraint
                    '  counts are passed: first parameter is number of equality constraints,
                    '  second one is number of inequality constraints.
                    ' 
                    '  As for constraining functions - these functions are passed as part
                    '  of problem Jacobian (see below).
                    ' 
                    '  NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    '        linear and general nonlinear constraints. This example does not
                    '        show how to work with general linear constraints, but you can
                    '        easily find it in documentation on minnlcsetbc() and
                    '        minnlcsetlc() functions.
                    ' 
                    xalglib.minnlcsetnlc(state, 1, 0)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target/constraints (C0 continuity violation)
                    '  * nonsmoothness of the target/constraints (C1 continuity violation)
                    '  * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minnlcoptguardsmoothness(state)
                    xalglib.minnlcoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0,f1} = { -x0+x1 , x0^2+x1^2-1 }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [  -1    +1  ]
                    '      J = [            ]
                    '          [ 2*x0  2*x1 ]
                    ' 
                    '  with f0 being target function, f1 being constraining function. Number
                    '  of equality/inequality constraints is specified by minnlcsetnlc(),
                    '  with equality ones always being first, inequality ones being last.
                    ' 
                    Dim rep As minnlcreport = New XAlglib.minnlcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}
                    xalglib.minnlcoptimize(state, AddressOf nlcfunc1_jac, Nothing, Nothing)
                    xalglib.minnlcresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){0.70710,-0.70710}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minnlcoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_equality")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minnlc_d_mixed
            '      Nonlinearly constrained optimization with mixed equality/inequality constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = x0+x1
                    ' 
                    '  subject to nonlinear inequality constraint
                    ' 
                    '     x0^2 + x1^2 - 1 <= 0
                    ' 
                    '  and nonlinear equality constraint
                    ' 
                    '     x2-exp(x0) = 0
                    ' 
                    Dim x0() As Double = New Double(){0,0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnlcstate = New XAlglib.minnlcstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minnlcreport = New XAlglib.minnlcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}

                    ' 
                    '  Create optimizer object and tune its settings:
                    '  * epsx=0.000001  stopping condition for inner iterations
                    '  * s=[1,1]        all variables have unit scale
                    '  * upper limit on step length is specified (to avoid probing locations where exp() is large)
                    ' 
                    xalglib.minnlccreate(3, x0, state)
                    xalglib.minnlcsetcond(state, epsx, maxits)
                    xalglib.minnlcsetscale(state, s)
                    xalglib.minnlcsetstpmax(state, 10.0)

                    ' 
                    '  Choose one of the nonlinear programming solvers supported by minnlc
                    '  optimizer:
                    '  * SLP - successive linear programming NLP solver
                    '  * AUL - augmented Lagrangian NLP solver
                    ' 
                    '  Different solvers have different properties:
                    '  * SLP is the most robust solver provided by ALGLIB: it can solve both
                    '    convex and nonconvex optimization problems, it respects box and
                    '    linear constraints (after you find feasible point it won't move away
                    '    from the feasible area) and tries to respect nonlinear constraints
                    '    as much as possible. It also usually needs less function evaluations
                    '    to converge than AUL.
                    '    However, it solves LP subproblems at each iterations which adds
                    '    significant overhead to its running time. Sometimes it can be as much
                    '    as 7x times slower than AUL.
                    '  * AUL solver is less robust than SLP - it can violate box and linear
                    '    constraints at any moment, and it is intended for convex optimization
                    '    problems (although in many cases it can deal with nonconvex ones too).
                    '    Also, unlike SLP it needs some tuning (penalty factor and number of
                    '    outer iterations).
                    '    However, it is often much faster than the current version of SLP.
                    ' 
                    '  In the code below we set solver to be AUL but then override it with SLP,
                    '  so the effective choice is to use SLP. We recommend you to use SLP at
                    '  least for early prototyping stages.
                    ' 
                    '  You can comment out line with SLP if you want to solve your problem with
                    '  AUL solver.
                    ' 
                    Dim rho As Double = 1000.0
                    Dim outerits As Integer = 5
                    xalglib.minnlcsetalgoaul(state, rho, outerits)
                    xalglib.minnlcsetalgoslp(state)

                    ' 
                    '  Set constraints:
                    ' 
                    '  Nonlinear constraints are tricky - you can not "pack" general
                    '  nonlinear function into double precision array. That's why
                    '  minnlcsetnlc() does not accept constraints itself - only constraint
                    '  counts are passed: first parameter is number of equality constraints,
                    '  second one is number of inequality constraints.
                    ' 
                    '  As for constraining functions - these functions are passed as part
                    '  of problem Jacobian (see below).
                    ' 
                    '  NOTE: MinNLC optimizer supports arbitrary combination of boundary, general
                    '        linear and general nonlinear constraints. This example does not
                    '        show how to work with boundary or general linear constraints, but you
                    '        can easily find it in documentation on minnlcsetbc() and
                    '        minnlcsetlc() functions.
                    ' 
                    xalglib.minnlcsetnlc(state, 1, 1)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target/constraints (C0 continuity violation)
                    '  * nonsmoothness of the target/constraints (C1 continuity violation)
                    '  * erroneous analytic Jacobian, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minnlcoptguardsmoothness(state)
                    xalglib.minnlcoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0,f1,f2} = { x0+x1 , x2-exp(x0) , x0^2+x1^2-1 }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [  +1      +1       0 ]
                    '      J = [-exp(x0)  0        1 ]
                    '          [ 2*x0    2*x1      0 ]
                    ' 
                    '  with f0 being target function, f1 being equality constraint "f1=0",
                    '  f2 being inequality constraint "f2<=0". Number of equality/inequality
                    '  constraints is specified by minnlcsetnlc(), with equality ones always
                    '  being first, inequality ones being last.
                    ' 
                    xalglib.minnlcoptimize(state, AddressOf nlcfunc2_jac, Nothing, Nothing)
                    xalglib.minnlcresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){-0.70710,-0.70710,0.49306}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minnlcoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minnlc_d_mixed")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minbc_d_1
            '      Nonlinear optimization with box constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 19
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  subject to box constraints
                    ' 
                    '      -1<=x<=+1, -1<=y<=+1
                    ' 
                    '  using MinBC optimizer with:
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minbcsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties:
                    '  * set box constraints
                    '  * set variable scales
                    '  * set stopping criteria
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim state As minbcstate = New XAlglib.minbcstate() ' initializer can be dropped, but compiler will issue warning
                    Dim epsg As Double = 0
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    xalglib.minbccreate(x, state)
                    xalglib.minbcsetbc(state, bndl, bndu)
                    xalglib.minbcsetscale(state, s)
                    xalglib.minbcsetcond(state, epsg, epsf, epsx, maxits)

                    ' 
                    '  Then we activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target function (C0 continuity violation)
                    '  * nonsmoothness of the target function (C1 continuity violation)
                    '  * erroneous analytic gradient, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.minbcoptguardsmoothness(state)
                    xalglib.minbcoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As minbcreport = New XAlglib.minbcreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbcoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.minbcresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-1,1}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbcoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minbc_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minbc_numdiff
            '      Nonlinear optimization with bound constraints and numerical differentiation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 22
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  subject to box constraints
                    ' 
                    '     -1<=x<=+1, -1<=y<=+1
                    ' 
                    '  using MinBC optimizer with:
                    '  * numerical differentiation being used
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see minbcsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim state As minbcstate = New XAlglib.minbcstate() ' initializer can be dropped, but compiler will issue warning
                    Dim epsg As Double = 0
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim diffstep As Double = 1.0e-6
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If

                    ' 
                    '  Now we are ready to actually optimize something:
                    '  * first we create optimizer
                    '  * we add boundary constraints
                    '  * we tune stopping conditions
                    '  * and, finally, optimize and obtain results...
                    ' 
                    xalglib.minbccreatef(x, diffstep, state)
                    xalglib.minbcsetbc(state, bndl, bndu)
                    xalglib.minbcsetscale(state, s)
                    xalglib.minbcsetcond(state, epsg, epsf, epsx, maxits)

                    ' 
                    '  Then we activate OptGuard integrity checking.
                    ' 
                    '  Numerical differentiation always produces "correct" gradient
                    '  (with some truncation error, but unbiased). Thus, we just have
                    '  to check smoothness properties of the target: C0 and C1 continuity.
                    ' 
                    '  Sometimes user accidentally tries to solve nonsmooth problems
                    '  with smooth optimizer. OptGuard helps to detect such situations
                    '  early, at the prototyping stage.
                    ' 
                    xalglib.minbcoptguardsmoothness(state)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As minbcreport = New XAlglib.minbcreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbcoptimize(state, AddressOf function1_func, Nothing, Nothing)
                    xalglib.minbcresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-1,1}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  Want to challenge OptGuard? Try to make your problem
                    '  nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    '  re-run optimizer.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minbcoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minbc_numdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minns_d_unconstrained
            '      Nonsmooth unconstrained optimization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = 2*|x0|+|x1|
                    ' 
                    '  using nonsmooth nonlinear optimizer.
                    ' 
                    Dim x0() As Double = New Double(){1,1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.00001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim radius As Double = 0.1
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        radius = v_spoil_real
                    End If
                    Dim rho As Double = 0.0
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnsstate = New XAlglib.minnsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minnsreport = New XAlglib.minnsreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}

                    ' 
                    '  Create optimizer object, choose AGS algorithm and tune its settings:
                    '  * radius=0.1     good initial value; will be automatically decreased later.
                    '  * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    '                   because we do not have such constraints
                    '  * epsx=0.000001  stopping conditions
                    '  * s=[1,1]        all variables have unit scale
                    ' 
                    xalglib.minnscreate(2, x0, state)
                    xalglib.minnssetalgoags(state, radius, rho)
                    xalglib.minnssetcond(state, epsx, maxits)
                    xalglib.minnssetscale(state, s)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints
                    '  (box/linear ones are passed separately by means of minnssetbc() and
                    '  minnssetlc() calls).
                    ' 
                    '  If you do not have nonlinear constraints (exactly our situation), then
                    '  you will have one-component function vector and 1xN Jacobian matrix.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0} = { 2*|x0|+|x1| }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [                       ]
                    '      J = [ 2*sign(x0)   sign(x1) ]
                    '          [                       ]
                    ' 
                    '  NOTE: nonsmooth optimizer requires considerably more function
                    '        evaluations than smooth solver - about 2N times more. Using
                    '        numerical differentiation introduces additional (multiplicative)
                    '        2N speedup.
                    ' 
                    '        It means that if smooth optimizer WITH user-supplied gradient
                    '        needs 100 function evaluations to solve 50-dimensional problem,
                    '        then AGS solver with user-supplied gradient will need about 10.000
                    '        function evaluations, and with numerical gradient about 1.000.000
                    '        function evaluations will be performed.
                    ' 
                    '  NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    '        optimization problems. It has convergence guarantees, i.e. it will
                    '        converge to stationary point of the function after running for some
                    '        time.
                    ' 
                    '        However, it is important to remember that "stationary point" is not
                    '        equal to "solution". If your problem is convex, everything is OK.
                    '        But nonconvex optimization problems may have "flat spots" - large
                    '        areas where gradient is exactly zero, but function value is far away
                    '        from optimal. Such areas are stationary points too, and optimizer
                    '        may be trapped here.
                    ' 
                    '        "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    '        orders of magnitude worse properties - they may be quite large and
                    '        hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    '        problem, because it is impossible to automatically distinguish "flat
                    '        spot" from true solution.
                    ' 
                    '        This note is here to warn you that you should be very careful when
                    '        you solve nonsmooth optimization problems. Visual inspection of
                    '        results is essential.
                    ' 
                    xalglib.minnsoptimize(state, AddressOf nsfunc1_jac, Nothing, Nothing)
                    xalglib.minnsresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){0.0000,0.0000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_unconstrained")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minns_d_diff
            '      Nonsmooth unconstrained optimization with numerical differentiation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 17
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = 2*|x0|+|x1|
                    ' 
                    '  using nonsmooth nonlinear optimizer with numerical
                    '  differentiation provided by ALGLIB.
                    ' 
                    '  NOTE: nonsmooth optimizer requires considerably more function
                    '        evaluations than smooth solver - about 2N times more. Using
                    '        numerical differentiation introduces additional (multiplicative)
                    '        2N speedup.
                    ' 
                    '        It means that if smooth optimizer WITH user-supplied gradient
                    '        needs 100 function evaluations to solve 50-dimensional problem,
                    '        then AGS solver with user-supplied gradient will need about 10.000
                    '        function evaluations, and with numerical gradient about 1.000.000
                    '        function evaluations will be performed.
                    ' 
                    Dim x0() As Double = New Double(){1,1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.00001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim diffstep As Double = 0.000001
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If
                    Dim radius As Double = 0.1
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        radius = v_spoil_real
                    End If
                    Dim rho As Double = 0.0
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnsstate = New XAlglib.minnsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minnsreport = New XAlglib.minnsreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}

                    ' 
                    '  Create optimizer object, choose AGS algorithm and tune its settings:
                    '  * radius=0.1     good initial value; will be automatically decreased later.
                    '  * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    '                   because we do not have such constraints
                    '  * epsx=0.000001  stopping conditions
                    '  * s=[1,1]        all variables have unit scale
                    ' 
                    xalglib.minnscreatef(2, x0, diffstep, state)
                    xalglib.minnssetalgoags(state, radius, rho)
                    xalglib.minnssetcond(state, epsx, maxits)
                    xalglib.minnssetscale(state, s)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function, with first component
                    '  being target function, and next components being nonlinear equality
                    '  and inequality constraints (box/linear ones are passed separately
                    '  by means of minnssetbc() and minnssetlc() calls).
                    ' 
                    '  If you do not have nonlinear constraints (exactly our situation), then
                    '  you will have one-component function vector.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0} = { 2*|x0|+|x1| }
                    ' 
                    xalglib.minnsoptimize(state, AddressOf nsfunc1_fvec, Nothing, Nothing)
                    xalglib.minnsresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){0.0000,0.0000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_diff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minns_d_bc
            '      Nonsmooth box constrained optimization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 16
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = 2*|x0|+|x1|
                    ' 
                    '  subject to box constraints
                    ' 
                    '         1 <= x0 < +INF
                    '      -INF <= x1 < +INF
                    ' 
                    '  using nonsmooth nonlinear optimizer.
                    ' 
                    Dim x0() As Double = New Double(){1,1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim bndl() As Double = New Double(){1,-Double.PositiveInfinity}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    Dim bndu() As Double = New Double(){Double.PositiveInfinity,+Double.PositiveInfinity}
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.00001
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim radius As Double = 0.1
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NaN
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        radius = v_spoil_real
                    End If
                    Dim rho As Double = 0.0
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnsstate = New XAlglib.minnsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minnsreport = New XAlglib.minnsreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}

                    ' 
                    '  Create optimizer object, choose AGS algorithm and tune its settings:
                    '  * radius=0.1     good initial value; will be automatically decreased later.
                    '  * rho=0.0        penalty coefficient for nonlinear constraints; can be zero
                    '                   because we do not have such constraints
                    '  * epsx=0.000001  stopping conditions
                    '  * s=[1,1]        all variables have unit scale
                    ' 
                    xalglib.minnscreate(2, x0, state)
                    xalglib.minnssetalgoags(state, radius, rho)
                    xalglib.minnssetcond(state, epsx, maxits)
                    xalglib.minnssetscale(state, s)

                    ' 
                    '  Set box constraints.
                    ' 
                    '  General linear constraints are set in similar way (see comments on
                    '  minnssetlc() function for more information).
                    ' 
                    '  You may combine box, linear and nonlinear constraints in one optimization
                    '  problem.
                    ' 
                    xalglib.minnssetbc(state, bndl, bndu)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints
                    '  (box/linear ones are passed separately by means of minnssetbc() and
                    '  minnssetlc() calls).
                    ' 
                    '  If you do not have nonlinear constraints (exactly our situation), then
                    '  you will have one-component function vector and 1xN Jacobian matrix.
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0} = { 2*|x0|+|x1| }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [                       ]
                    '      J = [ 2*sign(x0)   sign(x1) ]
                    '          [                       ]
                    ' 
                    '  NOTE: nonsmooth optimizer requires considerably more function
                    '        evaluations than smooth solver - about 2N times more. Using
                    '        numerical differentiation introduces additional (multiplicative)
                    '        2N speedup.
                    ' 
                    '        It means that if smooth optimizer WITH user-supplied gradient
                    '        needs 100 function evaluations to solve 50-dimensional problem,
                    '        then AGS solver with user-supplied gradient will need about 10.000
                    '        function evaluations, and with numerical gradient about 1.000.000
                    '        function evaluations will be performed.
                    ' 
                    '  NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    '        optimization problems. It has convergence guarantees, i.e. it will
                    '        converge to stationary point of the function after running for some
                    '        time.
                    ' 
                    '        However, it is important to remember that "stationary point" is not
                    '        equal to "solution". If your problem is convex, everything is OK.
                    '        But nonconvex optimization problems may have "flat spots" - large
                    '        areas where gradient is exactly zero, but function value is far away
                    '        from optimal. Such areas are stationary points too, and optimizer
                    '        may be trapped here.
                    ' 
                    '        "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    '        orders of magnitude worse properties - they may be quite large and
                    '        hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    '        problem, because it is impossible to automatically distinguish "flat
                    '        spot" from true solution.
                    ' 
                    '        This note is here to warn you that you should be very careful when
                    '        you solve nonsmooth optimization problems. Visual inspection of
                    '        results is essential.
                    ' 
                    ' 
                    xalglib.minnsoptimize(state, AddressOf nsfunc1_jac, Nothing, Nothing)
                    xalglib.minnsresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){1.0000,0.0000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_bc")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minns_d_nlc
            '      Nonsmooth nonlinearly constrained optimization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x0,x1) = 2*|x0|+|x1|
                    ' 
                    '  subject to combination of equality and inequality constraints
                    ' 
                    '       x0  =  1
                    '       x1 >= -1
                    ' 
                    '  using nonsmooth nonlinear optimizer. Although these constraints
                    '  are linear, we treat them as general nonlinear ones in order to
                    '  demonstrate nonlinearly constrained optimization setup.
                    ' 
                    Dim x0() As Double = New Double(){1,1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.00001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim radius As Double = 0.1
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        radius = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        radius = v_spoil_real
                    End If
                    Dim rho As Double = 50.0
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minnsstate = New XAlglib.minnsstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minnsreport = New XAlglib.minnsreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x1() As Double = New Double(){}

                    ' 
                    '  Create optimizer object, choose AGS algorithm and tune its settings:
                    '  * radius=0.1     good initial value; will be automatically decreased later.
                    '  * rho=50.0       penalty coefficient for nonlinear constraints. It is your
                    '                   responsibility to choose good one - large enough that it
                    '                   enforces constraints, but small enough in order to avoid
                    '                   extreme slowdown due to ill-conditioning.
                    '  * epsx=0.000001  stopping conditions
                    '  * s=[1,1]        all variables have unit scale
                    ' 
                    xalglib.minnscreate(2, x0, state)
                    xalglib.minnssetalgoags(state, radius, rho)
                    xalglib.minnssetcond(state, epsx, maxits)
                    xalglib.minnssetscale(state, s)

                    ' 
                    '  Set general nonlinear constraints.
                    ' 
                    '  This part is more tricky than working with box/linear constraints - you
                    '  can not "pack" general nonlinear function into double precision array.
                    '  That's why minnssetnlc() does not accept constraints itself - only
                    '  constraint COUNTS are passed: first parameter is number of equality
                    '  constraints, second one is number of inequality constraints.
                    ' 
                    '  As for constraining functions - these functions are passed as part
                    '  of problem Jacobian (see below).
                    ' 
                    '  NOTE: MinNS optimizer supports arbitrary combination of boundary, general
                    '        linear and general nonlinear constraints. This example does not
                    '        show how to work with general linear constraints, but you can
                    '        easily find it in documentation on minnlcsetlc() function.
                    ' 
                    xalglib.minnssetnlc(state, 1, 1)

                    ' 
                    '  Optimize and test results.
                    ' 
                    '  Optimizer object accepts vector function and its Jacobian, with first
                    '  component (Jacobian row) being target function, and next components
                    '  (Jacobian rows) being nonlinear equality and inequality constraints
                    '  (box/linear ones are passed separately by means of minnssetbc() and
                    '  minnssetlc() calls).
                    ' 
                    '  Nonlinear equality constraints have form Gi(x)=0, inequality ones
                    '  have form Hi(x)<=0, so we may have to "normalize" constraints prior
                    '  to passing them to optimizer (right side is zero, constraints are
                    '  sorted, multiplied by -1 when needed).
                    ' 
                    '  So, our vector function has form
                    ' 
                    '      {f0,f1,f2} = { 2*|x0|+|x1|,  x0-1, -x1-1 }
                    ' 
                    '  with Jacobian
                    ' 
                    '          [ 2*sign(x0)   sign(x1) ]
                    '      J = [     1           0     ]
                    '          [     0          -1     ]
                    ' 
                    '  which means that we have optimization problem
                    ' 
                    '      min{f0} subject to f1=0, f2<=0
                    ' 
                    '  which is essentially same as
                    ' 
                    '      min { 2*|x0|+|x1| } subject to x0=1, x1>=-1
                    ' 
                    '  NOTE: AGS solver used by us can handle nonsmooth and nonconvex
                    '        optimization problems. It has convergence guarantees, i.e. it will
                    '        converge to stationary point of the function after running for some
                    '        time.
                    ' 
                    '        However, it is important to remember that "stationary point" is not
                    '        equal to "solution". If your problem is convex, everything is OK.
                    '        But nonconvex optimization problems may have "flat spots" - large
                    '        areas where gradient is exactly zero, but function value is far away
                    '        from optimal. Such areas are stationary points too, and optimizer
                    '        may be trapped here.
                    ' 
                    '        "Flat spots" are nonsmooth equivalent of the saddle points, but with
                    '        orders of magnitude worse properties - they may be quite large and
                    '        hard to avoid. All nonsmooth optimizers are prone to this kind of the
                    '        problem, because it is impossible to automatically distinguish "flat
                    '        spot" from true solution.
                    ' 
                    '        This note is here to warn you that you should be very careful when
                    '        you solve nonsmooth optimization problems. Visual inspection of
                    '        results is essential.
                    ' 
                    xalglib.minnsoptimize(state, AddressOf nsfunc2_jac, Nothing, Nothing)
                    xalglib.minnsresults(state, x1, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x1, New Double(){1.0000,0.0000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minns_d_nlc")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST mincg_d_1
            '      Nonlinear optimization by CG
            '
            _TestResult = true
            For _spoil_scenario = -1 To 14
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  using nonlinear conjugate gradient method with:
                    '  * initial point x=[0,0]
                    '  * unit scale being set for all variables (see mincgsetscale for more info)
                    '  * stopping criteria set to "terminate after short enough step"
                    '  * OptGuard integrity check being used to check problem statement
                    '    for some common errors like nonsmoothness or bad analytic gradient
                    ' 
                    '  First, we create optimizer object and tune its properties
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As mincgstate = New XAlglib.mincgstate() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgcreate(x, state)
                    xalglib.mincgsetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.mincgsetscale(state, s)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to catch common coding and problem statement
                    '  issues, like:
                    '  * discontinuity of the target function (C0 continuity violation)
                    '  * nonsmoothness of the target function (C1 continuity violation)
                    '  * erroneous analytic gradient, i.e. one inconsistent with actual
                    '    change in the target/constraints
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: GRADIENT VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION. DO NOT USE IT IN PRODUCTION CODE!!!!!!!
                    ' 
                    '             Other OptGuard checks add moderate overhead, but anyway
                    '             it is better to turn them off when they are not needed.
                    ' 
                    xalglib.mincgoptguardsmoothness(state)
                    xalglib.mincgoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize and evaluate results
                    ' 
                    Dim rep As mincgreport = New XAlglib.mincgreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.mincgresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the gradient - say, add
                    '        1.0 to some of its components.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST mincg_d_2
            '      Nonlinear optimization with additional settings and restarts
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  This example demonstrates minimization of f(x,y) = 100*(x+3)^4+(y-3)^4
                    '  with nonlinear conjugate gradient method.
                    ' 
                    '  Several advanced techniques are demonstrated:
                    '  * upper limit on step size
                    '  * restart from new point
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim stpmax As Double = 0.1
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        stpmax = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        stpmax = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        stpmax = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As mincgstate = New XAlglib.mincgstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mincgreport = New XAlglib.mincgreport() ' initializer can be dropped, but compiler will issue warning

                    '  create and tune optimizer
                    xalglib.mincgcreate(x, state)
                    xalglib.mincgsetscale(state, s)
                    xalglib.mincgsetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.mincgsetstpmax(state, stpmax)

                    '  Set up OptGuard integrity checker which catches errors
                    '  like nonsmooth targets or errors in the analytic gradient.
                    ' 
                    '  OptGuard is essential at the early prototyping stages.
                    ' 
                    '  NOTE: gradient verification needs 3*N additional function
                    '        evaluations; DO NOT USE IT IN THE PRODUCTION CODE
                    '        because it leads to unnecessary slowdown of your app.
                    xalglib.mincgoptguardsmoothness(state)
                    xalglib.mincgoptguardgradient(state, 0.001)

                    '  first run
                    xalglib.mincgoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.mincgresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    '  second run - algorithm is restarted with mincgrestartfrom()
                    x = New Double(){10,10}
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    xalglib.mincgrestartfrom(state, x)
                    xalglib.mincgoptimize(state, AddressOf function1_grad, Nothing, Nothing)
                    xalglib.mincgresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    '  check OptGuard integrity report. Why do we need it at all?
                    '  Well, try breaking the gradient by adding 1.0 to some
                    '  of its components - OptGuard should report it as error.
                    '  And it may also catch unintended errors too :)
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "mincg_d_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST mincg_numdiff
            '      Nonlinear optimization by CG with numerical differentiation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 17
                Try
                    ' 
                    '  This example demonstrates minimization of
                    ' 
                    '      f(x,y) = 100*(x+3)^4+(y-3)^4
                    ' 
                    '  using numerical differentiation to calculate gradient.
                    ' 
                    '  We also show how to use OptGuard integrity checker to catch common
                    '  problem statement errors like accidentally specifying nonsmooth target
                    '  function.
                    ' 
                    '  First, we set up optimizer...
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsg As Double = 0
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsg = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsg = v_spoil_real
                    End If
                    Dim epsf As Double = 0
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsf = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsf = v_spoil_real
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim diffstep As Double = 1.0e-6
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As mincgstate = New XAlglib.mincgstate() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgcreatef(x, diffstep, state)
                    xalglib.mincgsetcond(state, epsg, epsf, epsx, maxits)
                    xalglib.mincgsetscale(state, s)

                    ' 
                    '  Then, we activate OptGuard integrity checking.
                    ' 
                    '  Numerical differentiation always produces "correct" gradient
                    '  (with some truncation error, but unbiased). Thus, we just have
                    '  to check smoothness properties of the target: C0 and C1 continuity.
                    ' 
                    '  Sometimes user accidentally tried to solve nonsmooth problems
                    '  with smooth optimizer. OptGuard helps to detect such situations
                    '  early, at the prototyping stage.
                    ' 
                    xalglib.mincgoptguardsmoothness(state)

                    ' 
                    '  Now we are ready to run the optimization
                    ' 
                    Dim rep As mincgreport = New XAlglib.mincgreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgoptimize(state, AddressOf function1_func, Nothing, Nothing)
                    xalglib.mincgresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,3}, 0.005)

                    ' 
                    '  ...and to check OptGuard integrity report.
                    ' 
                    '  Want to challenge OptGuard? Try to make your problem
                    '  nonsmooth by replacing 100*(x+3)^4 by 100*|x+3| and
                    '  re-run optimizer.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.mincgoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc0suspected, false)
                    _TestResult = _TestResult And doc_test_bool(ogrep.nonc1suspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "mincg_numdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_d_v
            '      Nonlinear least squares optimization using function vector only
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    ' 
                    '      f0(x0,x1) = 10*(x0+3)^2
                    '      f1(x0,x1) = (x1-3)^2
                    ' 
                    '  using "V" mode of the Levenberg-Marquardt optimizer.
                    ' 
                    '  Optimization algorithm uses:
                    '  * function vector f[] = {f1,f2}
                    ' 
                    '  No other information (Jacobian, gradient, etc.) is needed.
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Create optimizer, tell it to:
                    '  * use numerical differentiation with step equal to 0.0001
                    '  * use unit scale for all variables (s is a unit vector)
                    '  * stop after short enough step (less than epsx)
                    ' 
                    xalglib.minlmcreatev(2, x, 0.0001, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmsetscale(state, s)

                    ' 
                    '  Optimize
                    ' 
                    xalglib.minlmoptimize(state, AddressOf function1_fvec, Nothing, Nothing)

                    ' 
                    '  Test optimization results
                    ' 
                    '  NOTE: because we use numerical differentiation, we do not
                    '        verify Jacobian correctness - it is always "correct".
                    '        However, if you switch to analytic gradient, consider
                    '        checking it with OptGuard (see other examples).
                    ' 
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_v")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_d_vj
            '      Nonlinear least squares optimization using function vector and Jacobian
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    ' 
                    '      f0(x0,x1) = 10*(x0+3)^2
                    '      f1(x0,x1) = (x1-3)^2
                    ' 
                    '  using "VJ" mode of the Levenberg-Marquardt optimizer.
                    ' 
                    '  Optimization algorithm uses:
                    '  * function vector f[] = {f1,f2}
                    '  * Jacobian matrix J = {dfi/dxj}.
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Create optimizer, tell it to:
                    '  * use analytic gradient provided by user
                    '  * use unit scale for all variables (s is a unit vector)
                    '  * stop after short enough step (less than epsx)
                    ' 
                    xalglib.minlmcreatevj(2, x, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmsetscale(state, s)

                    ' 
                    '  Activate OptGuard integrity checking.
                    ' 
                    '  OptGuard monitor helps to detect erroneous analytic Jacobian,
                    '  i.e. one inconsistent with actual change in the target function.
                    ' 
                    '  OptGuard is essential for early prototyping stages because such
                    '  problems often result in premature termination of the optimizer
                    '  which is really hard to distinguish from the correct termination.
                    ' 
                    '  IMPORTANT: JACOBIAN VERIFICATION IS PERFORMED BY MEANS OF NUMERICAL
                    '             DIFFERENTIATION, THUS DO NOT USE IT IN PRODUCTION CODE!
                    ' 
                    xalglib.minlmoptguardgradient(state, 0.001)

                    ' 
                    '  Optimize
                    ' 
                    xalglib.minlmoptimize(state, AddressOf function1_fvec, AddressOf function1_jac, Nothing, Nothing)

                    ' 
                    '  Test optimization results
                    ' 
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)

                    ' 
                    '  Check that OptGuard did not report errors
                    ' 
                    '  NOTE: want to test OptGuard? Try breaking the Jacobian - say, add
                    '        1.0 to some of its components.
                    ' 
                    '  NOTE: unfortunately, specifics of LM optimization do not allow us
                    '        to detect errors like nonsmoothness (like we do with other
                    '        optimizers). So, only Jacobian correctness is verified.
                    ' 
                    Dim ogrep As optguardreport = New XAlglib.optguardreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlmoptguardresults(state, ogrep)
                    _TestResult = _TestResult And doc_test_bool(ogrep.badgradsuspected, false)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_vj")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_d_fgh
            '      Nonlinear Hessian-based optimization for general functions
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = 100*(x0+3)^4+(x1-3)^4
                    '  using "FGH" mode of the Levenberg-Marquardt optimizer.
                    ' 
                    '  F is treated like a monolitic function without internal structure,
                    '  i.e. we do NOT represent it as a sum of squares.
                    ' 
                    '  Optimization algorithm uses:
                    '  * function value F(x0,x1)
                    '  * gradient G={dF/dxi}
                    '  * Hessian H={d2F/(dxi*dxj)}
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.minlmcreatefgh(x, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmoptimize(state, AddressOf function1_func, AddressOf function1_grad, AddressOf function1_hess, Nothing, Nothing)
                    xalglib.minlmresults(state, x, rep)

                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_fgh")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_d_vb
            '      Bound constrained nonlinear least squares optimization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 12
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    ' 
                    '      f0(x0,x1) = 10*(x0+3)^2
                    '      f1(x0,x1) = (x1-3)^2
                    ' 
                    '  with boundary constraints
                    ' 
                    '      -1 <= x0 <= +1
                    '      -1 <= x1 <= +1
                    ' 
                    '  using "V" mode of the Levenberg-Marquardt optimizer.
                    ' 
                    '  Optimization algorithm uses:
                    '  * function vector f[] = {f1,f2}
                    ' 
                    '  No other information (Jacobian, gradient, etc.) is needed.
                    ' 
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim s() As Double = New Double(){1,1}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    Dim bndl() As Double = New Double(){-1,-1}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){+1,+1}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Create optimizer, tell it to:
                    '  * use numerical differentiation with step equal to 1.0
                    '  * use unit scale for all variables (s is a unit vector)
                    '  * stop after short enough step (less than epsx)
                    '  * set box constraints
                    ' 
                    xalglib.minlmcreatev(2, x, 0.0001, state)
                    xalglib.minlmsetbc(state, bndl, bndu)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmsetscale(state, s)

                    ' 
                    '  Optimize
                    ' 
                    xalglib.minlmoptimize(state, AddressOf function1_fvec, Nothing, Nothing)

                    ' 
                    '  Test optimization results
                    ' 
                    '  NOTE: because we use numerical differentiation, we do not
                    '        verify Jacobian correctness - it is always "correct".
                    '        However, if you switch to analytic gradient, consider
                    '        checking it with OptGuard (see other examples).
                    ' 
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-1,+1}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_vb")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_d_restarts
            '      Efficient restarts of LM optimizer
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  This example demonstrates minimization of F(x0,x1) = f0^2+f1^2, where 
                    ' 
                    '      f0(x0,x1) = 10*(x0+3)^2
                    '      f1(x0,x1) = (x1-3)^2
                    ' 
                    '  using several starting points and efficient restarts.
                    ' 
                    Dim x() As Double = New Double(){}
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  create optimizer using minlmcreatev()
                    ' 
                    x = New Double(){10,10}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    xalglib.minlmcreatev(2, x, 0.0001, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmoptimize(state, AddressOf function1_fvec, Nothing, Nothing)
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)

                    ' 
                    '  restart optimizer using minlmrestartfrom()
                    ' 
                    '  we can use different starting point, different function,
                    '  different stopping conditions, but problem size
                    '  must remain unchanged.
                    ' 
                    x = New Double(){4,4}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    xalglib.minlmrestartfrom(state, x)
                    xalglib.minlmoptimize(state, AddressOf function2_fvec, Nothing, Nothing)
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){0,1}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_d_restarts")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_t_1
            '      Nonlinear least squares optimization, FJ scheme (obsolete, but supported)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlmcreatefj(2, x, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmoptimize(state, AddressOf function1_func, AddressOf function1_jac, Nothing, Nothing)
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_t_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST minlm_t_2
            '      Nonlinear least squares optimization, FGJ scheme (obsolete, but supported)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim x() As Double = New Double(){0,0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.0000000001
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim state As minlmstate = New XAlglib.minlmstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As minlmreport = New XAlglib.minlmreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.minlmcreatefgj(2, x, state)
                    xalglib.minlmsetcond(state, epsx, maxits)
                    xalglib.minlmoptimize(state, AddressOf function1_func, AddressOf function1_grad, AddressOf function1_jac, Nothing, Nothing)
                    xalglib.minlmresults(state, x, rep)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){-3,+3}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "minlm_t_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_d_base
            '      Basic functionality (moments, adev, median, percentile)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim x() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim mean As Double
                    Dim variance As Double
                    Dim skewness As Double
                    Dim kurtosis As Double
                    Dim adev As Double
                    Dim p As Double
                    Dim v As Double

                    ' 
                    '  Here we demonstrate calculation of sample moments
                    '  (mean, variance, skewness, kurtosis)
                    ' 
                    xalglib.samplemoments(x, mean, variance, skewness, kurtosis)
                    _TestResult = _TestResult And doc_test_real(mean, 28.5, 0.01)
                    _TestResult = _TestResult And doc_test_real(variance, 801.1667, 0.01)
                    _TestResult = _TestResult And doc_test_real(skewness, 0.5751, 0.01)
                    _TestResult = _TestResult And doc_test_real(kurtosis, -1.2666, 0.01)

                    ' 
                    '  Average deviation
                    ' 
                    xalglib.sampleadev(x, adev)
                    _TestResult = _TestResult And doc_test_real(adev, 23.2, 0.01)

                    ' 
                    '  Median and percentile
                    ' 
                    xalglib.samplemedian(x, v)
                    _TestResult = _TestResult And doc_test_real(v, 20.5, 0.01)
                    p = 0.5
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        p = v_spoil_real
                    End If
                    xalglib.samplepercentile(x, p, v)
                    _TestResult = _TestResult And doc_test_real(v, 20.5, 0.01)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_base")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_d_c2
            '      Correlation (covariance) between two random variables
            '
            System.Console.WriteLine("50/151")
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    ' 
                    '  We have two samples - x and y, and want to measure dependency between them
                    ' 
                    Dim x() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim v As Double

                    ' 
                    '  Three dependency measures are calculated:
                    '  * covariation
                    '  * Pearson correlation
                    '  * Spearman rank correlation
                    ' 
                    v = xalglib.cov2(x, y)
                    _TestResult = _TestResult And doc_test_real(v, 82.5, 0.001)
                    v = xalglib.pearsoncorr2(x, y)
                    _TestResult = _TestResult And doc_test_real(v, 0.9627, 0.001)
                    v = xalglib.spearmancorr2(x, y)
                    _TestResult = _TestResult And doc_test_real(v, 1.000, 0.001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_c2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_d_cm
            '      Correlation (covariance) between components of random vector
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  X is a sample matrix:
                    '  * I-th row corresponds to I-th observation
                    '  * J-th column corresponds to J-th variable
                    ' 
                    Dim x(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    Dim c(,) As Double = New Double(,){{}}

                    ' 
                    '  Three dependency measures are calculated:
                    '  * covariation
                    '  * Pearson correlation
                    '  * Spearman rank correlation
                    ' 
                    '  Result is stored into C, with C[i,j] equal to correlation
                    '  (covariance) between I-th and J-th variables of X.
                    ' 
                    xalglib.covm(x, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{1.80,0.60,-1.40},{0.60,0.70,-0.80},{-1.40,-0.80,14.70}}, 0.01)
                    xalglib.pearsoncorrm(x, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{1.000,0.535,-0.272},{0.535,1.000,-0.249},{-0.272,-0.249,1.000}}, 0.01)
                    xalglib.spearmancorrm(x, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{1.000,0.556,-0.306},{0.556,1.000,-0.750},{-0.306,-0.750,1.000}}, 0.01)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_cm")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_d_cm2
            '      Correlation (covariance) between two random vectors
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  X and Y are sample matrices:
                    '  * I-th row corresponds to I-th observation
                    '  * J-th column corresponds to J-th variable
                    ' 
                    Dim x(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    Dim y(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y, v_spoil_real)
                    End If
                    Dim c(,) As Double = New Double(,){{}}

                    ' 
                    '  Three dependency measures are calculated:
                    '  * covariation
                    '  * Pearson correlation
                    '  * Spearman rank correlation
                    ' 
                    '  Result is stored into C, with C[i,j] equal to correlation
                    '  (covariance) between I-th variable of X and J-th variable of Y.
                    ' 
                    xalglib.covm2(x, y, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{4.100,-3.250},{2.450,-1.500},{13.450,-5.750}}, 0.01)
                    xalglib.pearsoncorrm2(x, y, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0.519,-0.699},{0.497,-0.518},{0.596,-0.433}}, 0.01)
                    xalglib.spearmancorrm2(x, y, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0.541,-0.649},{0.216,-0.433},{0.433,-0.135}}, 0.01)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_d_cm2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_t_base
            '      Tests ability to detect errors in inputs
            '
            _TestResult = true
            For _spoil_scenario = -1 To 33
                Try
                    Dim mean As Double
                    Dim variance As Double
                    Dim skewness As Double
                    Dim kurtosis As Double
                    Dim adev As Double
                    Dim p As Double
                    Dim v As Double

                    ' 
                    '  first, we test short form of functions
                    ' 
                    Dim x1() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    xalglib.samplemoments(x1, mean, variance, skewness, kurtosis)
                    Dim x2() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    xalglib.sampleadev(x2, adev)
                    Dim x3() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    xalglib.samplemedian(x3, v)
                    Dim x4() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x4, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x4, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x4, v_spoil_real)
                    End If
                    p = 0.5
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        p = v_spoil_real
                    End If
                    xalglib.samplepercentile(x4, p, v)

                    ' 
                    '  and then we test full form
                    ' 
                    Dim x5() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x5, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x5, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x5, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        spoil_vector_by_deleting_element(x5)
                    End If
                    xalglib.samplemoments(x5, 10, mean, variance, skewness, kurtosis)
                    Dim x6() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x6, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x6, v_spoil_real)
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x6, v_spoil_real)
                    End If
                    If _spoil_scenario=22 Then
                        spoil_vector_by_deleting_element(x6)
                    End If
                    xalglib.sampleadev(x6, 10, adev)
                    Dim x7() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=23 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=24 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=25 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=26 Then
                        spoil_vector_by_deleting_element(x7)
                    End If
                    xalglib.samplemedian(x7, 10, v)
                    Dim x8() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=27 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=28 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=29 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=30 Then
                        spoil_vector_by_deleting_element(x8)
                    End If
                    p = 0.5
                    If _spoil_scenario=31 Then
                        v_spoil_real = Double.NaN
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=32 Then
                        v_spoil_real = Double.PositiveInfinity
                        p = v_spoil_real
                    End If
                    If _spoil_scenario=33 Then
                        v_spoil_real = Double.NegativeInfinity
                        p = v_spoil_real
                    End If
                    xalglib.samplepercentile(x8, 10, p, v)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_t_base")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST basestat_t_covcorr
            '      Tests ability to detect errors in inputs
            '
            _TestResult = true
            For _spoil_scenario = -1 To 125
                Try
                    Dim v As Double
                    Dim c(,) As Double = New Double(,){{}}

                    ' 
                    '  2-sample short-form cov/corr are tested
                    ' 
                    Dim x1() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x1, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x1)
                    End If
                    Dim y1() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y1, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y1, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y1, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y1, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y1)
                    End If
                    v = xalglib.cov2(x1, y1)
                    Dim x2() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        spoil_vector_by_deleting_element(x2)
                    End If
                    Dim y2() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y2, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y2, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y2, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y2, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        spoil_vector_by_deleting_element(y2)
                    End If
                    v = xalglib.pearsoncorr2(x2, y2)
                    Dim x3() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=23 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x3, v_spoil_real)
                    End If
                    If _spoil_scenario=24 Then
                        spoil_vector_by_deleting_element(x3)
                    End If
                    Dim y3() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=25 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y3, v_spoil_real)
                    End If
                    If _spoil_scenario=26 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y3, v_spoil_real)
                    End If
                    If _spoil_scenario=27 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y3, v_spoil_real)
                    End If
                    If _spoil_scenario=28 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y3, v_spoil_real)
                    End If
                    If _spoil_scenario=29 Then
                        spoil_vector_by_deleting_element(y3)
                    End If
                    v = xalglib.spearmancorr2(x3, y3)

                    ' 
                    '  2-sample full-form cov/corr are tested
                    ' 
                    Dim x1a() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=30 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x1a, v_spoil_real)
                    End If
                    If _spoil_scenario=31 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x1a, v_spoil_real)
                    End If
                    If _spoil_scenario=32 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x1a, v_spoil_real)
                    End If
                    If _spoil_scenario=33 Then
                        spoil_vector_by_deleting_element(x1a)
                    End If
                    Dim y1a() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=34 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y1a, v_spoil_real)
                    End If
                    If _spoil_scenario=35 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y1a, v_spoil_real)
                    End If
                    If _spoil_scenario=36 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y1a, v_spoil_real)
                    End If
                    If _spoil_scenario=37 Then
                        spoil_vector_by_deleting_element(y1a)
                    End If
                    v = xalglib.cov2(x1a, y1a, 10)
                    Dim x2a() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=38 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x2a, v_spoil_real)
                    End If
                    If _spoil_scenario=39 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x2a, v_spoil_real)
                    End If
                    If _spoil_scenario=40 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x2a, v_spoil_real)
                    End If
                    If _spoil_scenario=41 Then
                        spoil_vector_by_deleting_element(x2a)
                    End If
                    Dim y2a() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=42 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y2a, v_spoil_real)
                    End If
                    If _spoil_scenario=43 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y2a, v_spoil_real)
                    End If
                    If _spoil_scenario=44 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y2a, v_spoil_real)
                    End If
                    If _spoil_scenario=45 Then
                        spoil_vector_by_deleting_element(y2a)
                    End If
                    v = xalglib.pearsoncorr2(x2a, y2a, 10)
                    Dim x3a() As Double = New Double(){0,1,4,9,16,25,36,49,64,81}
                    If _spoil_scenario=46 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x3a, v_spoil_real)
                    End If
                    If _spoil_scenario=47 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x3a, v_spoil_real)
                    End If
                    If _spoil_scenario=48 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x3a, v_spoil_real)
                    End If
                    If _spoil_scenario=49 Then
                        spoil_vector_by_deleting_element(x3a)
                    End If
                    Dim y3a() As Double = New Double(){0,1,2,3,4,5,6,7,8,9}
                    If _spoil_scenario=50 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y3a, v_spoil_real)
                    End If
                    If _spoil_scenario=51 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y3a, v_spoil_real)
                    End If
                    If _spoil_scenario=52 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y3a, v_spoil_real)
                    End If
                    If _spoil_scenario=53 Then
                        spoil_vector_by_deleting_element(y3a)
                    End If
                    v = xalglib.spearmancorr2(x3a, y3a, 10)

                    ' 
                    '  vector short-form cov/corr are tested.
                    ' 
                    Dim x4(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=54 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x4, v_spoil_real)
                    End If
                    If _spoil_scenario=55 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x4, v_spoil_real)
                    End If
                    If _spoil_scenario=56 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x4, v_spoil_real)
                    End If
                    xalglib.covm(x4, c)
                    Dim x5(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=57 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x5, v_spoil_real)
                    End If
                    If _spoil_scenario=58 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x5, v_spoil_real)
                    End If
                    If _spoil_scenario=59 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x5, v_spoil_real)
                    End If
                    xalglib.pearsoncorrm(x5, c)
                    Dim x6(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=60 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x6, v_spoil_real)
                    End If
                    If _spoil_scenario=61 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x6, v_spoil_real)
                    End If
                    If _spoil_scenario=62 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x6, v_spoil_real)
                    End If
                    xalglib.spearmancorrm(x6, c)

                    ' 
                    '  vector full-form cov/corr are tested.
                    ' 
                    Dim x7(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=63 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=64 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=65 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x7, v_spoil_real)
                    End If
                    If _spoil_scenario=66 Then
                        spoil_matrix_by_deleting_row(x7)
                    End If
                    If _spoil_scenario=67 Then
                        spoil_matrix_by_deleting_col(x7)
                    End If
                    xalglib.covm(x7, 5, 3, c)
                    Dim x8(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=68 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=69 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=70 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x8, v_spoil_real)
                    End If
                    If _spoil_scenario=71 Then
                        spoil_matrix_by_deleting_row(x8)
                    End If
                    If _spoil_scenario=72 Then
                        spoil_matrix_by_deleting_col(x8)
                    End If
                    xalglib.pearsoncorrm(x8, 5, 3, c)
                    Dim x9(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=73 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x9, v_spoil_real)
                    End If
                    If _spoil_scenario=74 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x9, v_spoil_real)
                    End If
                    If _spoil_scenario=75 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x9, v_spoil_real)
                    End If
                    If _spoil_scenario=76 Then
                        spoil_matrix_by_deleting_row(x9)
                    End If
                    If _spoil_scenario=77 Then
                        spoil_matrix_by_deleting_col(x9)
                    End If
                    xalglib.spearmancorrm(x9, 5, 3, c)

                    ' 
                    '  cross-vector short-form cov/corr are tested.
                    ' 
                    Dim x10(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=78 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x10, v_spoil_real)
                    End If
                    If _spoil_scenario=79 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x10, v_spoil_real)
                    End If
                    If _spoil_scenario=80 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x10, v_spoil_real)
                    End If
                    Dim y10(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=81 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y10, v_spoil_real)
                    End If
                    If _spoil_scenario=82 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y10, v_spoil_real)
                    End If
                    If _spoil_scenario=83 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y10, v_spoil_real)
                    End If
                    xalglib.covm2(x10, y10, c)
                    Dim x11(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=84 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x11, v_spoil_real)
                    End If
                    If _spoil_scenario=85 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x11, v_spoil_real)
                    End If
                    If _spoil_scenario=86 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x11, v_spoil_real)
                    End If
                    Dim y11(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=87 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y11, v_spoil_real)
                    End If
                    If _spoil_scenario=88 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y11, v_spoil_real)
                    End If
                    If _spoil_scenario=89 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y11, v_spoil_real)
                    End If
                    xalglib.pearsoncorrm2(x11, y11, c)
                    Dim x12(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=90 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x12, v_spoil_real)
                    End If
                    If _spoil_scenario=91 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x12, v_spoil_real)
                    End If
                    If _spoil_scenario=92 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x12, v_spoil_real)
                    End If
                    Dim y12(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=93 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y12, v_spoil_real)
                    End If
                    If _spoil_scenario=94 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y12, v_spoil_real)
                    End If
                    If _spoil_scenario=95 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y12, v_spoil_real)
                    End If
                    xalglib.spearmancorrm2(x12, y12, c)

                    ' 
                    '  cross-vector full-form cov/corr are tested.
                    ' 
                    Dim x13(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=96 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x13, v_spoil_real)
                    End If
                    If _spoil_scenario=97 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x13, v_spoil_real)
                    End If
                    If _spoil_scenario=98 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x13, v_spoil_real)
                    End If
                    If _spoil_scenario=99 Then
                        spoil_matrix_by_deleting_row(x13)
                    End If
                    If _spoil_scenario=100 Then
                        spoil_matrix_by_deleting_col(x13)
                    End If
                    Dim y13(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=101 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y13, v_spoil_real)
                    End If
                    If _spoil_scenario=102 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y13, v_spoil_real)
                    End If
                    If _spoil_scenario=103 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y13, v_spoil_real)
                    End If
                    If _spoil_scenario=104 Then
                        spoil_matrix_by_deleting_row(y13)
                    End If
                    If _spoil_scenario=105 Then
                        spoil_matrix_by_deleting_col(y13)
                    End If
                    xalglib.covm2(x13, y13, 5, 3, 2, c)
                    Dim x14(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=106 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x14, v_spoil_real)
                    End If
                    If _spoil_scenario=107 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x14, v_spoil_real)
                    End If
                    If _spoil_scenario=108 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x14, v_spoil_real)
                    End If
                    If _spoil_scenario=109 Then
                        spoil_matrix_by_deleting_row(x14)
                    End If
                    If _spoil_scenario=110 Then
                        spoil_matrix_by_deleting_col(x14)
                    End If
                    Dim y14(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=111 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y14, v_spoil_real)
                    End If
                    If _spoil_scenario=112 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y14, v_spoil_real)
                    End If
                    If _spoil_scenario=113 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y14, v_spoil_real)
                    End If
                    If _spoil_scenario=114 Then
                        spoil_matrix_by_deleting_row(y14)
                    End If
                    If _spoil_scenario=115 Then
                        spoil_matrix_by_deleting_col(y14)
                    End If
                    xalglib.pearsoncorrm2(x14, y14, 5, 3, 2, c)
                    Dim x15(,) As Double = New Double(,){{1,0,1},{1,1,0},{-1,1,0},{-2,-1,1},{-1,0,9}}
                    If _spoil_scenario=116 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x15, v_spoil_real)
                    End If
                    If _spoil_scenario=117 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x15, v_spoil_real)
                    End If
                    If _spoil_scenario=118 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x15, v_spoil_real)
                    End If
                    If _spoil_scenario=119 Then
                        spoil_matrix_by_deleting_row(x15)
                    End If
                    If _spoil_scenario=120 Then
                        spoil_matrix_by_deleting_col(x15)
                    End If
                    Dim y15(,) As Double = New Double(,){{2,3},{2,1},{-1,6},{-9,9},{7,1}}
                    If _spoil_scenario=121 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(y15, v_spoil_real)
                    End If
                    If _spoil_scenario=122 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(y15, v_spoil_real)
                    End If
                    If _spoil_scenario=123 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(y15, v_spoil_real)
                    End If
                    If _spoil_scenario=124 Then
                        spoil_matrix_by_deleting_row(y15)
                    End If
                    If _spoil_scenario=125 Then
                        spoil_matrix_by_deleting_col(y15)
                    End If
                    xalglib.spearmancorrm2(x15, y15, 5, 3, 2, c)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "basestat_t_covcorr")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ssa_d_basic
            '      Simple SSA analysis demo
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Here we demonstrate SSA trend/noise separation for some toy problem:
                    '  small monotonically growing series X are analyzed with 3-tick window
                    '  and "top-K" version of SSA, which selects K largest singular vectors
                    '  for analysis, with K=1.
                    ' 
                    Dim s As ssamodel = New XAlglib.ssamodel() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){0,0.5,1,1,1.5,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If

                    ' 
                    '  First, we create SSA model, set its properties and add dataset.
                    ' 
                    '  We use window with width=3 and configure model to use direct SSA
                    '  algorithm - one which runs exact O(N*W^2) analysis - to extract
                    '  one top singular vector. Well, it is toy problem :)
                    ' 
                    '  NOTE: SSA model may store and analyze more than one sequence
                    '        (say, different sequences may correspond to data collected
                    '        from different devices)
                    ' 
                    xalglib.ssacreate(s)
                    xalglib.ssasetwindow(s, 3)
                    xalglib.ssaaddsequence(s, x)
                    xalglib.ssasetalgotopkdirect(s, 1)

                    ' 
                    '  Now we begin analysis. Internally SSA model stores everything it needs:
                    '  data, settings, solvers and so on. Right after first call to analysis-
                    '  related function it will analyze dataset, build basis and perform analysis.
                    ' 
                    '  Subsequent calls to analysis functions will reuse previously computed
                    '  basis, unless you invalidate it by changing model settings (or dataset).
                    ' 
                    Dim trend() As Double = New Double(){}
                    Dim noise() As Double = New Double(){}
                    xalglib.ssaanalyzesequence(s, x, trend, noise)
                    _TestResult = _TestResult And doc_test_real_vector(trend, New Double(){0.3815,0.5582,0.7810,1.0794,1.5041,2.0105}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_basic")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ssa_d_forecast
            '      Simple SSA forecasting demo
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Here we demonstrate SSA forecasting on some toy problem with clearly
                    '  visible linear trend and small amount of noise.
                    ' 
                    Dim s As ssamodel = New XAlglib.ssamodel() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){0.05,0.96,2.04,3.11,3.97,5.03,5.98,7.02,8.02}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If

                    ' 
                    '  First, we create SSA model, set its properties and add dataset.
                    ' 
                    '  We use window with width=3 and configure model to use direct SSA
                    '  algorithm - one which runs exact O(N*W^2) analysis - to extract
                    '  two top singular vectors. Well, it is toy problem :)
                    ' 
                    '  NOTE: SSA model may store and analyze more than one sequence
                    '        (say, different sequences may correspond to data collected
                    '        from different devices)
                    ' 
                    xalglib.ssacreate(s)
                    xalglib.ssasetwindow(s, 3)
                    xalglib.ssaaddsequence(s, x)
                    xalglib.ssasetalgotopkdirect(s, 2)

                    ' 
                    '  Now we begin analysis. Internally SSA model stores everything it needs:
                    '  data, settings, solvers and so on. Right after first call to analysis-
                    '  related function it will analyze dataset, build basis and perform analysis.
                    ' 
                    '  Subsequent calls to analysis functions will reuse previously computed
                    '  basis, unless you invalidate it by changing model settings (or dataset).
                    ' 
                    '  In this example we show how to use ssaforecastlast() function, which
                    '  predicts changed in the last sequence of the dataset. If you want to
                    '  perform prediction for some other sequence, use ssaforecastsequence().
                    ' 
                    Dim trend() As Double = New Double(){}
                    xalglib.ssaforecastlast(s, 3, trend)

                    ' 
                    '  Well, we expected it to be [9,10,11]. There exists some difference,
                    '  which can be explained by the artificial noise in the dataset.
                    ' 
                    _TestResult = _TestResult And doc_test_real_vector(trend, New Double(){9.0005,9.9322,10.8051}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_forecast")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST ssa_d_realtime
            '      Real-time SSA algorithm with fast incremental updates
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    ' 
                    '  Suppose that you have a constant stream of incoming data, and you want
                    '  to regularly perform singular spectral analysis of this stream.
                    ' 
                    '  One full run of direct algorithm costs O(N*Width^2) operations, so
                    '  the more points you have, the more it costs to rebuild basis from
                    '  scratch.
                    '  
                    '  Luckily we have incremental SSA algorithm which can perform quick
                    '  updates of already computed basis in O(K*Width^2) ops, where K
                    '  is a number of singular vectors extracted. Usually it is orders of
                    '  magnitude faster than full update of the basis.
                    ' 
                    '  In this example we start from some initial dataset x0. Then we
                    '  start appending elements one by one to the end of the last sequence.
                    ' 
                    '  NOTE: direct algorithm also supports incremental updates, but
                    '        with O(Width^3) cost. Typically K<<Width, so specialized
                    '        incremental algorithm is still faster.
                    ' 
                    Dim s1 As ssamodel = New XAlglib.ssamodel() ' initializer can be dropped, but compiler will issue warning
                    Dim a1(,) As Double = New Double(,){{}}
                    Dim sv1() As Double = New Double(){}
                    Dim w As Integer
                    Dim k As Integer
                    Dim x0() As Double = New Double(){0.009,0.976,1.999,2.984,3.977,5.002}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x0, v_spoil_real)
                    End If
                    xalglib.ssacreate(s1)
                    xalglib.ssasetwindow(s1, 3)
                    xalglib.ssaaddsequence(s1, x0)

                    '  set algorithm to the real-time version of top-K, K=2
                    xalglib.ssasetalgotopkrealtime(s1, 2)

                    '  one more interesting feature of the incremental algorithm is "power-up" cycle.
                    '  even with incremental algorithm initial basis calculation costs O(N*Width^2) ops.
                    '  if such startup cost is too high for your real-time app, then you may divide
                    '  initial basis calculation across several model updates. It results in better
                    '  latency at the price of somewhat lesser precision during first few updates.
                    xalglib.ssasetpoweruplength(s1, 3)

                    '  now, after we prepared everything, start to add incoming points one by one;
                    '  in the real life, of course, we will perform some work between subsequent update
                    '  (analyze something, predict, and so on).
                    ' 
                    '  After each append we perform one iteration of the real-time solver. Usually
                    '  one iteration is more than enough to update basis. If you have REALLY tight
                    '  performance constraints, you may specify fractional amount of iterations,
                    '  which means that iteration is performed with required probability.
                    Dim updateits As Double = 1.0
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        updateits = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        updateits = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        updateits = v_spoil_real
                    End If
                    xalglib.ssaappendpointandupdate(s1, 5.951, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 7.074, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 7.925, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 8.992, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 9.942, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 11.051, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 11.965, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 13.047, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    xalglib.ssaappendpointandupdate(s1, 13.970, updateits)
                    xalglib.ssagetbasis(s1, a1, sv1, w, k)

                    '  Ok, we have our basis in a1[] and singular values at sv1[].
                    '  But is it good enough? Let's print it.
                    _TestResult = _TestResult And doc_test_real_matrix(a1, New Double(,){{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}}, 0.0005)

                    '  Ok, two vectors with 3 components each.
                    '  But how to understand that is it really good basis?
                    '  Let's compare it with direct SSA algorithm on the entire sequence.
                    Dim s2 As ssamodel = New XAlglib.ssamodel() ' initializer can be dropped, but compiler will issue warning
                    Dim a2(,) As Double = New Double(,){{}}
                    Dim sv2() As Double = New Double(){}
                    Dim x2() As Double = New Double(){0.009,0.976,1.999,2.984,3.977,5.002,5.951,7.074,7.925,8.992,9.942,11.051,11.965,13.047,13.970}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x2, v_spoil_real)
                    End If
                    xalglib.ssacreate(s2)
                    xalglib.ssasetwindow(s2, 3)
                    xalglib.ssaaddsequence(s2, x2)
                    xalglib.ssasetalgotopkdirect(s2, 2)
                    xalglib.ssagetbasis(s2, a2, sv2, w, k)

                    '  it is exactly the same as one calculated with incremental approach!
                    _TestResult = _TestResult And doc_test_real_matrix(a2, New Double(,){{0.510607,0.753611},{0.575201,0.058445},{0.639081,-0.654717}}, 0.0005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "ssa_d_realtime")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST linreg_d_basic
            '      Linear regression used to build the very basic model and unpack coefficients
            '
            _TestResult = True
            try
                ' 
                '  In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
                ' 
                '  We have:
                '  * xy - matrix of basic function values (exp(0.5*x)) and expected values
                ' 
                Dim xy(,) As Double = New Double(,){{0.606531,1.133719},{0.670320,1.306522},{0.740818,1.504604},{0.818731,1.554663},{0.904837,1.884638},{1.000000,2.072436},{1.105171,2.257285},{1.221403,2.534068},{1.349859,2.622017},{1.491825,2.897713},{1.648721,3.219371}}
                Dim info As Integer
                Dim nvars As Integer
                Dim model As linearmodel = New XAlglib.linearmodel() ' initializer can be dropped, but compiler will issue warning
                Dim rep As lrreport = New XAlglib.lrreport() ' initializer can be dropped, but compiler will issue warning
                Dim c() As Double = New Double(){}

                xalglib.lrbuildz(xy, 11, 1, info, model, rep)
                _TestResult = _TestResult And doc_test_int(info, 1)
                xalglib.lrunpack(model, c, nvars)
                _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.98650,0.00000}, 0.00005)
            Catch E As AlglibException
                _TestResult = False
            Catch E As Exception
                Throw
            End Try
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "linreg_d_basic")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST filters_d_sma
            '      SMA(k) filter
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Here we demonstrate SMA(k) filtering for time series.
                    ' 
                    Dim x() As Double = New Double(){5,6,7,8}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If

                    ' 
                    '  Apply filter.
                    '  We should get [5, 5.5, 6.5, 7.5] as result
                    ' 
                    xalglib.filtersma(x, 2)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){5,5.5,6.5,7.5}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_sma")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST filters_d_ema
            '      EMA(alpha) filter
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Here we demonstrate EMA(0.5) filtering for time series.
                    ' 
                    Dim x() As Double = New Double(){5,6,7,8}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If

                    ' 
                    '  Apply filter.
                    '  We should get [5, 5.5, 6.25, 7.125] as result
                    ' 
                    xalglib.filterema(x, 0.5)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){5,5.5,6.25,7.125}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_ema")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST filters_d_lrma
            '      LRMA(k) filter
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Here we demonstrate LRMA(3) filtering for time series.
                    ' 
                    Dim x() As Double = New Double(){7,8,8,9,12,12}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If

                    ' 
                    '  Apply filter.
                    '  We should get [7.0000, 8.0000, 8.1667, 8.8333, 11.6667, 12.5000] as result
                    '     
                    xalglib.filterlrma(x, 3)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){7.0000,8.0000,8.1667,8.8333,11.6667,12.5000}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "filters_d_lrma")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST mcpd_simple1
            '      Simple unconstrained MCPD model (no entry/exit states)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  The very simple MCPD example
                    ' 
                    '  We have a loan portfolio. Our loans can be in one of two states:
                    '  * normal loans ("good" ones)
                    '  * past due loans ("bad" ones)
                    ' 
                    '  We assume that:
                    '  * loans can transition from any state to any other state. In 
                    '    particular, past due loan can become "good" one at any moment 
                    '    with same (fixed) probability. Not realistic, but it is toy example :)
                    '  * portfolio size does not change over time
                    ' 
                    '  Thus, we have following model
                    '      state_new = P*state_old
                    '  where
                    '          ( p00  p01 )
                    '      P = (          )
                    '          ( p10  p11 )
                    ' 
                    '  We want to model transitions between these two states using MCPD
                    '  approach (Markov Chains for Proportional/Population Data), i.e.
                    '  to restore hidden transition matrix P using actual portfolio data.
                    '  We have:
                    '  * poportional data, i.e. proportion of loans in the normal and past 
                    '    due states (not portfolio size measured in some currency, although 
                    '    it is possible to work with population data too)
                    '  * two tracks, i.e. two sequences which describe portfolio
                    '    evolution from two different starting states: [1,0] (all loans 
                    '    are "good") and [0.8,0.2] (only 80% of portfolio is in the "good"
                    '    state)
                    ' 
                    Dim s As mcpdstate = New XAlglib.mcpdstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mcpdreport = New XAlglib.mcpdreport() ' initializer can be dropped, but compiler will issue warning
                    Dim p(,) As Double = New Double(,){{}}
                    Dim track0(,) As Double = New Double(,){{1.00000,0.00000},{0.95000,0.05000},{0.92750,0.07250},{0.91738,0.08263},{0.91282,0.08718}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    Dim track1(,) As Double = New Double(,){{0.80000,0.20000},{0.86000,0.14000},{0.88700,0.11300},{0.89915,0.10085}}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If

                    xalglib.mcpdcreate(2, s)
                    xalglib.mcpdaddtrack(s, track0)
                    xalglib.mcpdaddtrack(s, track1)
                    xalglib.mcpdsolve(s)
                    xalglib.mcpdresults(s, p, rep)

                    ' 
                    '  Hidden matrix P is equal to
                    '          ( 0.95  0.50 )
                    '          (            )
                    '          ( 0.05  0.50 )
                    '  which means that "good" loans can become "bad" with 5% probability, 
                    '  while "bad" loans will return to good state with 50% probability.
                    ' 
                    _TestResult = _TestResult And doc_test_real_matrix(p, New Double(,){{0.95,0.50},{0.05,0.50}}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "mcpd_simple1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST mcpd_simple2
            '      Simple MCPD model (no entry/exit states) with equality constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  Simple MCPD example
                    ' 
                    '  We have a loan portfolio. Our loans can be in one of three states:
                    '  * normal loans
                    '  * past due loans
                    '  * charged off loans
                    ' 
                    '  We assume that:
                    '  * normal loan can stay normal or become past due (but not charged off)
                    '  * past due loan can stay past due, become normal or charged off
                    '  * charged off loan will stay charged off for the rest of eternity
                    '  * portfolio size does not change over time
                    '  Not realistic, but it is toy example :)
                    ' 
                    '  Thus, we have following model
                    '      state_new = P*state_old
                    '  where
                    '          ( p00  p01    )
                    '      P = ( p10  p11    )
                    '          (      p21  1 )
                    '  i.e. four elements of P are known a priori.
                    ' 
                    '  Although it is possible (given enough data) to In order to enforce 
                    '  this property we set equality constraints on these elements.
                    ' 
                    '  We want to model transitions between these two states using MCPD
                    '  approach (Markov Chains for Proportional/Population Data), i.e.
                    '  to restore hidden transition matrix P using actual portfolio data.
                    '  We have:
                    '  * poportional data, i.e. proportion of loans in the current and past 
                    '    due states (not portfolio size measured in some currency, although 
                    '    it is possible to work with population data too)
                    '  * two tracks, i.e. two sequences which describe portfolio
                    '    evolution from two different starting states: [1,0,0] (all loans 
                    '    are "good") and [0.8,0.2,0.0] (only 80% of portfolio is in the "good"
                    '    state)
                    ' 
                    Dim s As mcpdstate = New XAlglib.mcpdstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mcpdreport = New XAlglib.mcpdreport() ' initializer can be dropped, but compiler will issue warning
                    Dim p(,) As Double = New Double(,){{}}
                    Dim track0(,) As Double = New Double(,){{1.000000,0.000000,0.000000},{0.950000,0.050000,0.000000},{0.927500,0.060000,0.012500},{0.911125,0.061375,0.027500},{0.896256,0.060900,0.042844}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(track0, v_spoil_real)
                    End If
                    Dim track1(,) As Double = New Double(,){{0.800000,0.200000,0.000000},{0.860000,0.090000,0.050000},{0.862000,0.065500,0.072500},{0.851650,0.059475,0.088875},{0.838805,0.057451,0.103744}}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(track1, v_spoil_real)
                    End If

                    xalglib.mcpdcreate(3, s)
                    xalglib.mcpdaddtrack(s, track0)
                    xalglib.mcpdaddtrack(s, track1)
                    xalglib.mcpdaddec(s, 0, 2, 0.0)
                    xalglib.mcpdaddec(s, 1, 2, 0.0)
                    xalglib.mcpdaddec(s, 2, 2, 1.0)
                    xalglib.mcpdaddec(s, 2, 0, 0.0)
                    xalglib.mcpdsolve(s)
                    xalglib.mcpdresults(s, p, rep)

                    ' 
                    '  Hidden matrix P is equal to
                    '          ( 0.95 0.50      )
                    '          ( 0.05 0.25      )
                    '          (      0.25 1.00 ) 
                    '  which means that "good" loans can become past due with 5% probability, 
                    '  while past due loans will become charged off with 25% probability or
                    '  return back to normal state with 50% probability.
                    ' 
                    _TestResult = _TestResult And doc_test_real_matrix(p, New Double(,){{0.95,0.50,0.00},{0.05,0.25,0.00},{0.00,0.25,1.00}}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "mcpd_simple2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_regr
            '      Regression problem with one output (2=>1)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple example on neural network: network is trained to reproduce
                    '  small 2x2 multiplication table.
                    ' 
                    '  NOTE: we use network with excessive amount of neurons, which guarantees
                    '        almost exact reproduction of the training set. Generalization ability
                    '        of such network is rather low, but we are not concerned with such
                    '        questions in this basic demo.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Training set:
                    '  * one row corresponds to one record A*B=C in the multiplication table
                    '  * first two columns store A and B, last column stores C
                    ' 
                    '  [1 * 1 = 1]
                    '  [1 * 2 = 2]
                    '  [2 * 1 = 2]
                    '  [2 * 2 = 4]
                    ' 
                    Dim xy(,) As Double = New Double(,){{1,1,1},{1,2,2},{2,1,2},{2,2,4}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    '  Network is created.
                    '  Trainer object is created.
                    '  Dataset is attached to trainer object.
                    ' 
                    xalglib.mlpcreatetrainer(2, 1, trn)
                    xalglib.mlpcreate1(2, 5, 1, network)
                    xalglib.mlpsetdataset(trn, xy, 4)

                    ' 
                    '  Network is trained with 5 restarts from random positions
                    ' 
                    xalglib.mlptrainnetwork(trn, network, 5, rep)

                    ' 
                    '  2*2=?
                    ' 
                    Dim x() As Double = New Double(){2,2}
                    Dim y() As Double = New Double(){0}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){4.000}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_regr")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_regr_n
            '      Regression problem with multiple outputs (2=>2)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Network with 2 inputs and 2 outputs is trained to reproduce vector function:
                    '      (x0,x1) => (x0+x1, x0*x1)
                    ' 
                    '  Informally speaking, we want neural network to simultaneously calculate
                    '  both sum of two numbers and their product.
                    ' 
                    '  NOTE: we use network with excessive amount of neurons, which guarantees
                    '        almost exact reproduction of the training set. Generalization ability
                    '        of such network is rather low, but we are not concerned with such
                    '        questions in this basic demo.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Training set. One row corresponds to one record [A,B,A+B,A*B].
                    ' 
                    '  [ 1   1  1+1  1*1 ]
                    '  [ 1   2  1+2  1*2 ]
                    '  [ 2   1  2+1  2*1 ]
                    '  [ 2   2  2+2  2*2 ]
                    ' 
                    Dim xy(,) As Double = New Double(,){{1,1,2,1},{1,2,3,2},{2,1,3,2},{2,2,4,4}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    '  Network is created.
                    '  Trainer object is created.
                    '  Dataset is attached to trainer object.
                    ' 
                    xalglib.mlpcreatetrainer(2, 2, trn)
                    xalglib.mlpcreate1(2, 5, 2, network)
                    xalglib.mlpsetdataset(trn, xy, 4)

                    ' 
                    '  Network is trained with 5 restarts from random positions
                    ' 
                    xalglib.mlptrainnetwork(trn, network, 5, rep)

                    ' 
                    '  2+1=?
                    '  2*1=?
                    ' 
                    Dim x() As Double = New Double(){2,1}
                    Dim y() As Double = New Double(){0,0}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){3.000,2.000}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_regr_n")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_cls2
            '      Binary classification problem
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Suppose that we want to classify numbers as positive (class 0) and negative
                    '  (class 1). We have training set which includes several strictly positive
                    '  or negative numbers - and zero.
                    ' 
                    '  The problem is that we are not sure how to classify zero, so from time to
                    '  time we mark it as positive or negative (with equal probability). Other
                    '  numbers are marked in pure deterministic setting. How will neural network
                    '  cope with such classification task?
                    ' 
                    '  NOTE: we use network with excessive amount of neurons, which guarantees
                    '        almost exact reproduction of the training set. Generalization ability
                    '        of such network is rather low, but we are not concerned with such
                    '        questions in this basic demo.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){0}
                    Dim y() As Double = New Double(){0,0}

                    ' 
                    '  Training set. One row corresponds to one record [A => class(A)].
                    ' 
                    '  Classes are denoted by numbers from 0 to 1, where 0 corresponds to positive
                    '  numbers and 1 to negative numbers.
                    ' 
                    '  [ +1  0]
                    '  [ +2  0]
                    '  [ -1  1]
                    '  [ -2  1]
                    '  [  0  0]   !! sometimes we classify 0 as positive, sometimes as negative
                    '  [  0  1]   !!
                    ' 
                    Dim xy(,) As Double = New Double(,){{+1,0},{+2,0},{-1,1},{-2,1},{0,0},{0,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    ' 
                    '  When we solve classification problems, everything is slightly different from
                    '  the regression ones:
                    ' 
                    '  1. Network is created. Because we solve classification problem, we use
                    '     mlpcreatec1() function instead of mlpcreate1(). This function creates
                    '     classifier network with SOFTMAX-normalized outputs. This network returns
                    '     vector of class membership probabilities which are normalized to be
                    '     non-negative and sum to 1.0
                    ' 
                    '  2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
                    '     create trainer object. Trainer object process dataset and neural network
                    '     slightly differently to account for specifics of the classification
                    '     problems.
                    ' 
                    '  3. Dataset is attached to trainer object. Note that dataset format is slightly
                    '     different from one used for regression.
                    ' 
                    xalglib.mlpcreatetrainercls(1, 2, trn)
                    xalglib.mlpcreatec1(1, 5, 2, network)
                    xalglib.mlpsetdataset(trn, xy, 6)

                    ' 
                    '  Network is trained with 5 restarts from random positions
                    ' 
                    xalglib.mlptrainnetwork(trn, network, 5, rep)

                    ' 
                    '  Test our neural network on strictly positive and strictly negative numbers.
                    ' 
                    '  IMPORTANT! Classifier network returns class membership probabilities instead
                    '  of class indexes. Network returns two values (probabilities) instead of one
                    '  (class index).
                    ' 
                    '  Thus, for +1 we expect to get [P0,P1] = [1,0], where P0 is probability that
                    '  number is positive (belongs to class 0), and P1 is probability that number
                    '  is negative (belongs to class 1).
                    ' 
                    '  For -1 we expect to get [P0,P1] = [0,1]
                    ' 
                    '  Following properties are guaranteed by network architecture:
                    '  * P0>=0, P1>=0   non-negativity
                    '  * P0+P1=1        normalization
                    ' 
                    x = New Double(){1}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){1.000,0.000}, 0.05)
                    x = New Double(){-1}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,1.000}, 0.05)

                    ' 
                    '  But what our network will return for 0, which is between classes 0 and 1?
                    ' 
                    '  In our dataset it has two different marks assigned (class 0 AND class 1).
                    '  So network will return something average between class 0 and class 1:
                    '      0 => [0.5, 0.5]
                    ' 
                    x = New Double(){0}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.500,0.500}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_cls2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_cls3
            '      Multiclass classification problem
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  Suppose that we want to classify numbers as positive (class 0) and negative
                    '  (class 1). We also have one more class for zero (class 2).
                    ' 
                    '  NOTE: we use network with excessive amount of neurons, which guarantees
                    '        almost exact reproduction of the training set. Generalization ability
                    '        of such network is rather low, but we are not concerned with such
                    '        questions in this basic demo.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){0}
                    Dim y() As Double = New Double(){0,0,0}

                    ' 
                    '  Training set. One row corresponds to one record [A => class(A)].
                    ' 
                    '  Classes are denoted by numbers from 0 to 2, where 0 corresponds to positive
                    '  numbers, 1 to negative numbers, 2 to zero
                    ' 
                    '  [ +1  0]
                    '  [ +2  0]
                    '  [ -1  1]
                    '  [ -2  1]
                    '  [  0  2]
                    ' 
                    Dim xy(,) As Double = New Double(,){{+1,0},{+2,0},{-1,1},{-2,1},{0,2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    ' 
                    '  When we solve classification problems, everything is slightly different from
                    '  the regression ones:
                    ' 
                    '  1. Network is created. Because we solve classification problem, we use
                    '     mlpcreatec1() function instead of mlpcreate1(). This function creates
                    '     classifier network with SOFTMAX-normalized outputs. This network returns
                    '     vector of class membership probabilities which are normalized to be
                    '     non-negative and sum to 1.0
                    ' 
                    '  2. We use mlpcreatetrainercls() function instead of mlpcreatetrainer() to
                    '     create trainer object. Trainer object process dataset and neural network
                    '     slightly differently to account for specifics of the classification
                    '     problems.
                    ' 
                    '  3. Dataset is attached to trainer object. Note that dataset format is slightly
                    '     different from one used for regression.
                    ' 
                    xalglib.mlpcreatetrainercls(1, 3, trn)
                    xalglib.mlpcreatec1(1, 5, 3, network)
                    xalglib.mlpsetdataset(trn, xy, 5)

                    ' 
                    '  Network is trained with 5 restarts from random positions
                    ' 
                    xalglib.mlptrainnetwork(trn, network, 5, rep)

                    ' 
                    '  Test our neural network on strictly positive and strictly negative numbers.
                    ' 
                    '  IMPORTANT! Classifier network returns class membership probabilities instead
                    '  of class indexes. Network returns three values (probabilities) instead of one
                    '  (class index).
                    ' 
                    '  Thus, for +1 we expect to get [P0,P1,P2] = [1,0,0],
                    '  for -1 we expect to get [P0,P1,P2] = [0,1,0],
                    '  and for 0 we will get [P0,P1,P2] = [0,0,1].
                    ' 
                    '  Following properties are guaranteed by network architecture:
                    '  * P0>=0, P1>=0, P2>=0    non-negativity
                    '  * P0+P1+P2=1             normalization
                    ' 
                    x = New Double(){1}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){1.000,0.000,0.000}, 0.05)
                    x = New Double(){-1}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,1.000,0.000}, 0.05)
                    x = New Double(){0}
                    xalglib.mlpprocess(network, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,0.000,1.000}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_cls3")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_trainerobject
            '      Advanced example on trainer object
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  Trainer object is used to train network. It stores dataset, training settings,
                    '  and other information which is NOT part of neural network. You should use
                    '  trainer object as follows:
                    '  (1) you create trainer object and specify task type (classification/regression)
                    '      and number of inputs/outputs
                    '  (2) you add dataset to the trainer object
                    '  (3) you may change training settings (stopping criteria or weight decay)
                    '  (4) finally, you may train one or more networks
                    ' 
                    '  You may interleave stages 2...4 and repeat them many times. Trainer object
                    '  remembers its internal state and can be used several times after its creation
                    '  and initialization.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Stage 1: object creation.
                    ' 
                    '  We have to specify number of inputs and outputs. Trainer object can be used
                    '  only for problems with same number of inputs/outputs as was specified during
                    '  its creation.
                    ' 
                    '  In case you want to train SOFTMAX-normalized network which solves classification
                    '  problems,  you  must  use  another  function  to  create  trainer  object:
                    '  mlpcreatetrainercls().
                    ' 
                    '  Below we create trainer object which can be used to train regression networks
                    '  with 2 inputs and 1 output.
                    ' 
                    xalglib.mlpcreatetrainer(2, 1, trn)

                    ' 
                    '  Stage 2: specification of the training set
                    ' 
                    '  By default trainer object stores empty dataset. So to solve your non-empty problem
                    '  you have to set dataset by passing to trainer dense or sparse matrix.
                    ' 
                    '  One row of the matrix corresponds to one record A*B=C in the multiplication table.
                    '  First two columns store A and B, last column stores C
                    ' 
                    '      [1 * 1 = 1]   [ 1 1 1 ]
                    '      [1 * 2 = 2]   [ 1 2 2 ]
                    '      [2 * 1 = 2] = [ 2 1 2 ]
                    '      [2 * 2 = 4]   [ 2 2 4 ]
                    ' 
                    Dim xy(,) As Double = New Double(,){{1,1,1},{1,2,2},{2,1,2},{2,2,4}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.mlpsetdataset(trn, xy, 4)

                    ' 
                    '  Stage 3: modification of the training parameters.
                    ' 
                    '  You may modify parameters like weights decay or stopping criteria:
                    '  * we set moderate weight decay
                    '  * we choose iterations limit as stopping condition (another condition - step size -
                    '    is zero, which means than this condition is not active)
                    ' 
                    Dim wstep As Double = 0.000
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        wstep = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        wstep = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        wstep = v_spoil_real
                    End If
                    Dim maxits As Integer = 100
                    xalglib.mlpsetdecay(trn, 0.01)
                    xalglib.mlpsetcond(trn, wstep, maxits)

                    ' 
                    '  Stage 4: training.
                    ' 
                    '  We will train several networks with different architecture using same trainer object.
                    '  We may change training parameters or even dataset, so different networks are trained
                    '  differently. But in this simple example we will train all networks with same settings.
                    ' 
                    '  We create and train three networks:
                    '  * network 1 has 2x1 architecture     (2 inputs, no hidden neurons, 1 output)
                    '  * network 2 has 2x5x1 architecture   (2 inputs, 5 hidden neurons, 1 output)
                    '  * network 3 has 2x5x5x1 architecture (2 inputs, two hidden layers, 1 output)
                    ' 
                    '  NOTE: these networks solve regression problems. For classification problems you
                    '        should use mlpcreatec0/c1/c2 to create neural networks which have SOFTMAX-
                    '        normalized outputs.
                    ' 
                    Dim net1 As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim net2 As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim net3 As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.mlpcreate0(2, 1, net1)
                    xalglib.mlpcreate1(2, 5, 1, net2)
                    xalglib.mlpcreate2(2, 5, 5, 1, net3)

                    xalglib.mlptrainnetwork(trn, net1, 5, rep)
                    xalglib.mlptrainnetwork(trn, net2, 5, rep)
                    xalglib.mlptrainnetwork(trn, net3, 5, rep)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_trainerobject")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_crossvalidation
            '      Cross-validation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example shows how to perform cross-validation with ALGLIB
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Training set: f(x)=1/(x^2+1)
                    '  One row corresponds to one record [x,f(x)]
                    ' 
                    Dim xy(,) As Double = New Double(,){{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    '  Trainer object is created.
                    '  Dataset is attached to trainer object.
                    ' 
                    '  NOTE: it is not good idea to perform cross-validation on sample
                    '        as small as ours (13 examples). It is done for demonstration
                    '        purposes only. Generalization error estimates won't be
                    '        precise enough for practical purposes.
                    ' 
                    xalglib.mlpcreatetrainer(1, 1, trn)
                    xalglib.mlpsetdataset(trn, xy, 13)

                    ' 
                    '  The key property of the cross-validation is that it estimates
                    '  generalization properties of neural ARCHITECTURE. It does NOT
                    '  estimates generalization error of some specific network which
                    '  is passed to the k-fold CV routine.
                    ' 
                    '  In our example we create 1x4x1 neural network and pass it to
                    '  CV routine without training it. Original state of the network
                    '  is not used for cross-validation - each round is restarted from
                    '  random initial state. Only geometry of network matters.
                    ' 
                    '  We perform 5 restarts from different random positions for each
                    '  of the 10 cross-validation rounds.
                    ' 
                    xalglib.mlpcreate1(1, 4, 1, network)
                    xalglib.mlpkfoldcv(trn, network, 5, 10, rep)

                    ' 
                    '  Cross-validation routine stores estimates of the generalization
                    '  error to MLP report structure. You may examine its fields and
                    '  see estimates of different errors (RMS, CE, Avg).
                    ' 
                    '  Because cross-validation is non-deterministic, in our manual we
                    '  can not say what values will be stored to rep after call to
                    '  mlpkfoldcv(). Every CV round will return slightly different
                    '  estimates.
                    '
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_crossvalidation")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_ensembles_es
            '      Early stopping ensembles
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example shows how to train early stopping ensebles.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim ensemble As mlpensemble = New XAlglib.mlpensemble() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Training set: f(x)=1/(x^2+1)
                    '  One row corresponds to one record [x,f(x)]
                    ' 
                    Dim xy(,) As Double = New Double(,){{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    ' 
                    '  Trainer object is created.
                    '  Dataset is attached to trainer object.
                    ' 
                    '  NOTE: it is not good idea to use early stopping ensemble on sample
                    '        as small as ours (13 examples). It is done for demonstration
                    '        purposes only. Ensemble training algorithm won't find good
                    '        solution on such small sample.
                    ' 
                    xalglib.mlpcreatetrainer(1, 1, trn)
                    xalglib.mlpsetdataset(trn, xy, 13)

                    ' 
                    '  Ensemble is created and trained. Each of 50 network is trained
                    '  with 5 restarts.
                    ' 
                    xalglib.mlpecreate1(1, 4, 1, 50, ensemble)
                    xalglib.mlptrainensemblees(trn, ensemble, 5, rep)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_ensembles_es")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST nn_parallel
            '      Parallel training
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example shows how to use parallel functionality of ALGLIB.
                    '  We generate simple 1-dimensional regression problem and show how
                    '  to use parallel training, parallel cross-validation, parallel
                    '  training of neural ensembles.
                    ' 
                    '  We assume that you already know how to use ALGLIB in serial mode
                    '  and concentrate on its parallel capabilities.
                    ' 
                    '  NOTE: it is not good idea to use parallel features on sample as small
                    '        as ours (13 examples). It is done only for demonstration purposes.
                    ' 
                    Dim trn As mlptrainer = New XAlglib.mlptrainer() ' initializer can be dropped, but compiler will issue warning
                    Dim network As multilayerperceptron = New XAlglib.multilayerperceptron() ' initializer can be dropped, but compiler will issue warning
                    Dim ensemble As mlpensemble = New XAlglib.mlpensemble() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As mlpreport = New XAlglib.mlpreport() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{-2.0,0.2},{-1.6,0.3},{-1.3,0.4},{-1,0.5},{-0.6,0.7},{-0.3,0.9},{0,1},{2.0,0.2},{1.6,0.3},{1.3,0.4},{1,0.5},{0.6,0.7},{0.3,0.9}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.mlpcreatetrainer(1, 1, trn)
                    xalglib.mlpsetdataset(trn, xy, 13)
                    xalglib.mlpcreate1(1, 4, 1, network)
                    xalglib.mlpecreate1(1, 4, 1, 50, ensemble)

                    ' 
                    '  Below we demonstrate how to perform:
                    '  * parallel training of individual networks
                    '  * parallel cross-validation
                    '  * parallel training of neural ensembles
                    ' 
                    '  In order to use multithreading, you have to:
                    '  1) Install SMP edition of ALGLIB.
                    '  2) This step is specific for C++ users: you should activate OS-specific
                    '     capabilities of ALGLIB by defining AE_OS=AE_POSIX (for *nix systems)
                    '     or AE_OS=AE_WINDOWS (for Windows systems).
                    '     C# users do not have to perform this step because C# programs are
                    '     portable across different systems without OS-specific tuning.
                    '  3) Tell ALGLIB that you want it to use multithreading by means of
                    '     setnworkers() call:
                    '           * alglib::setnworkers(0)  = use all cores
                    '           * alglib::setnworkers(-1) = leave one core unused
                    '           * alglib::setnworkers(-2) = leave two cores unused
                    '           * alglib::setnworkers(+2) = use 2 cores (even if you have more)
                    '     During runtime ALGLIB will automatically determine whether it is
                    '     feasible to start worker threads and split your task between cores.
                    ' 
                    xalglib.setnworkers(+2)

                    ' 
                    '  First, we perform parallel training of individual network with 5
                    '  restarts from random positions. These 5 rounds of  training  are
                    '  executed in parallel manner,  with  best  network  chosen  after
                    '  training.
                    ' 
                    '  ALGLIB can use additional way to speed up computations -  divide
                    '  dataset   into   smaller   subsets   and   process these subsets
                    '  simultaneously. It allows us  to  efficiently  parallelize  even
                    '  single training round. This operation is performed automatically
                    '  for large datasets, but our toy dataset is too small.
                    ' 
                    xalglib.mlptrainnetwork(trn, network, 5, rep)

                    ' 
                    '  Then, we perform parallel 10-fold cross-validation, with 5 random
                    '  restarts per each CV round. I.e., 5*10=50  networks  are trained
                    '  in total. All these operations can be parallelized.
                    ' 
                    '  NOTE: again, ALGLIB can parallelize  calculation   of   gradient
                    '        over entire dataset - but our dataset is too small.
                    ' 
                    xalglib.mlpkfoldcv(trn, network, 5, 10, rep)

                    ' 
                    '  Finally, we train early stopping ensemble of 50 neural networks,
                    '  each  of them is trained with 5 random restarts. I.e.,  5*50=250
                    '  networks aretrained in total.
                    ' 
                    xalglib.mlptrainensemblees(trn, ensemble, 5, rep)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "nn_parallel")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST clst_ahc
            '      Simple hierarchical clusterization with Euclidean distance function
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple clusterization example
                    ' 
                    '  We have a set of points in 2D space:
                    '      (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    ' 
                    '   |
                    '   |     P3
                    '   |
                    '   | P1          
                    '   |             P4
                    '   | P0          P2
                    '   |-------------------------
                    ' 
                    '  We want to perform Agglomerative Hierarchic Clusterization (AHC),
                    '  using complete linkage (default algorithm) and Euclidean distance
                    '  (default metric).
                    ' 
                    '  In order to do that, we:
                    '  * create clusterizer with clusterizercreate()
                    '  * set points XY and metric (2=Euclidean) with clusterizersetpoints()
                    '  * run AHC algorithm with clusterizerrunahc
                    ' 
                    '  You may see that clusterization itself is a minor part of the example,
                    '  most of which is dominated by comments :)
                    ' 
                    Dim s As clusterizerstate = New XAlglib.clusterizerstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As ahcreport = New XAlglib.ahcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{1,1},{1,2},{4,1},{2,3},{4,1.5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.clusterizercreate(s)
                    xalglib.clusterizersetpoints(s, xy, 2)
                    xalglib.clusterizerrunahc(s, rep)

                    ' 
                    '  Now we've built our clusterization tree. Rep.z contains information which
                    '  is required to build dendrogram. I-th row of rep.z represents one merge
                    '  operation, with first cluster to merge having index rep.z[I,0] and second
                    '  one having index rep.z[I,1]. Merge result has index NPoints+I.
                    ' 
                    '  Clusters with indexes less than NPoints are single-point initial clusters,
                    '  while ones with indexes from NPoints to 2*NPoints-2 are multi-point
                    '  clusters created during merges.
                    ' 
                    '  In our example, Z=[[2,4], [0,1], [3,6], [5,7]]
                    ' 
                    '  It means that:
                    '  * first, we merge C2=(P2) and C4=(P4),    and create C5=(P2,P4)
                    '  * then, we merge  C2=(P0) and C1=(P1),    and create C6=(P0,P1)
                    '  * then, we merge  C3=(P3) and C6=(P0,P1), and create C7=(P0,P1,P3)
                    '  * finally, we merge C5 and C7 and create C8=(P0,P1,P2,P3,P4)
                    ' 
                    '  Thus, we have following dendrogram:
                    '   
                    '       ------8-----
                    '       |          |
                    '       |      ----7----
                    '       |      |       |
                    '    ---5---   |    ---6---
                    '    |     |   |    |     |
                    '    P2   P4   P3   P0   P1
                    ' 
                    _TestResult = _TestResult And doc_test_int_matrix(rep.z, New Integer(,){{2,4},{0,1},{3,6},{5,7}})

                    ' 
                    '  We've built dendrogram above by reordering our dataset.
                    ' 
                    '  Without such reordering it would be impossible to build dendrogram without
                    '  intersections. Luckily, ahcreport structure contains two additional fields
                    '  which help to build dendrogram from your data:
                    '  * rep.p, which contains permutation applied to dataset
                    '  * rep.pm, which contains another representation of merges 
                    ' 
                    '  In our example we have:
                    '  * P=[3,4,0,2,1]
                    '  * PZ=[[0,0,1,1,0,0],[3,3,4,4,0,0],[2,2,3,4,0,1],[0,1,2,4,1,2]]
                    ' 
                    '  Permutation array P tells us that P0 should be moved to position 3,
                    '  P1 moved to position 4, P2 moved to position 0 and so on:
                    ' 
                    '    (P0 P1 P2 P3 P4) => (P2 P4 P3 P0 P1)
                    ' 
                    '  Merges array PZ tells us how to perform merges on the sorted dataset.
                    '  One row of PZ corresponds to one merge operations, with first pair of
                    '  elements denoting first of the clusters to merge (start index, end
                    '  index) and next pair of elements denoting second of the clusters to
                    '  merge. Clusters being merged are always adjacent, with first one on
                    '  the left and second one on the right.
                    ' 
                    '  For example, first row of PZ tells us that clusters [0,0] and [1,1] are
                    '  merged (single-point clusters, with first one containing P2 and second
                    '  one containing P4). Third row of PZ tells us that we merge one single-
                    '  point cluster [2,2] with one two-point cluster [3,4].
                    ' 
                    '  There are two more elements in each row of PZ. These are the helper
                    '  elements, which denote HEIGHT (not size) of left and right subdendrograms.
                    '  For example, according to PZ, first two merges are performed on clusterization
                    '  trees of height 0, while next two merges are performed on 0-1 and 1-2
                    '  pairs of trees correspondingly.
                    ' 
                    _TestResult = _TestResult And doc_test_int_vector(rep.p, New Integer(){3,4,0,2,1})
                    _TestResult = _TestResult And doc_test_int_matrix(rep.pm, New Integer(,){{0,0,1,1,0,0},{3,3,4,4,0,0},{2,2,3,4,0,1},{0,1,2,4,1,2}})
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "clst_ahc")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST clst_kmeans
            '      Simple k-means clusterization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple clusterization example
                    ' 
                    '  We have a set of points in 2D space:
                    '      (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    ' 
                    '   |
                    '   |     P3
                    '   |
                    '   | P1          
                    '   |             P4
                    '   | P0          P2
                    '   |-------------------------
                    ' 
                    '  We want to perform k-means++ clustering with K=2.
                    ' 
                    '  In order to do that, we:
                    '  * create clusterizer with clusterizercreate()
                    '  * set points XY and metric (must be Euclidean, distype=2) with clusterizersetpoints()
                    '  * (optional) set number of restarts from random positions to 5
                    '  * run k-means algorithm with clusterizerrunkmeans()
                    ' 
                    '  You may see that clusterization itself is a minor part of the example,
                    '  most of which is dominated by comments :)
                    ' 
                    Dim s As clusterizerstate = New XAlglib.clusterizerstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As kmeansreport = New XAlglib.kmeansreport() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{1,1},{1,2},{4,1},{2,3},{4,1.5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.clusterizercreate(s)
                    xalglib.clusterizersetpoints(s, xy, 2)
                    xalglib.clusterizersetkmeanslimits(s, 5, 0)
                    xalglib.clusterizerrunkmeans(s, 2, rep)

                    ' 
                    '  We've performed clusterization, and it succeeded (completion code is +1).
                    ' 
                    '  Now first center is stored in the first row of rep.c, second one is stored
                    '  in the second row. rep.cidx can be used to determine which center is
                    '  closest to some specific point of the dataset.
                    ' 
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)

                    '  We called clusterizersetpoints() with disttype=2 because k-means++
                    '  algorithm does NOT support metrics other than Euclidean. But what if we
                    '  try to use some other metric?
                    ' 
                    '  We change metric type by calling clusterizersetpoints() one more time,
                    '  and try to run k-means algo again. It fails.
                    ' 
                    xalglib.clusterizersetpoints(s, xy, 0)
                    xalglib.clusterizerrunkmeans(s, 2, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, -5)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "clst_kmeans")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST clst_linkage
            '      Clusterization with different linkage types
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  We have a set of points in 1D space:
                    '      (P0,P1,P2,P3,P4) = (1, 3, 10, 16, 20)
                    ' 
                    '  We want to perform Agglomerative Hierarchic Clusterization (AHC),
                    '  using either complete or single linkage and Euclidean distance
                    '  (default metric).
                    ' 
                    '  First two steps merge P0/P1 and P3/P4 independently of the linkage type.
                    '  However, third step depends on linkage type being used:
                    '  * in case of complete linkage P2=10 is merged with [P0,P1]
                    '  * in case of single linkage P2=10 is merged with [P3,P4]
                    ' 
                    Dim s As clusterizerstate = New XAlglib.clusterizerstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As ahcreport = New XAlglib.ahcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{1},{3},{10},{16},{20}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    Dim cidx() As Integer = New Integer(){}
                    Dim cz() As Integer = New Integer(){}

                    xalglib.clusterizercreate(s)
                    xalglib.clusterizersetpoints(s, xy, 2)

                    '  use complete linkage, reduce set down to 2 clusters.
                    '  print clusterization with clusterizergetkclusters(2).
                    '  P2 must belong to [P0,P1]
                    xalglib.clusterizersetahcalgo(s, 0)
                    xalglib.clusterizerrunahc(s, rep)
                    xalglib.clusterizergetkclusters(rep, 2, cidx, cz)
                    _TestResult = _TestResult And doc_test_int_vector(cidx, New Integer(){1,1,1,0,0})

                    '  use single linkage, reduce set down to 2 clusters.
                    '  print clusterization with clusterizergetkclusters(2).
                    '  P2 must belong to [P2,P3]
                    xalglib.clusterizersetahcalgo(s, 1)
                    xalglib.clusterizerrunahc(s, rep)
                    xalglib.clusterizergetkclusters(rep, 2, cidx, cz)
                    _TestResult = _TestResult And doc_test_int_vector(cidx, New Integer(){0,0,1,1,1})
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "clst_linkage")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST clst_distance
            '      Clusterization with different metric types
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  We have three points in 4D space:
                    '      (P0,P1,P2) = ((1, 2, 1, 2), (6, 7, 6, 7), (7, 6, 7, 6))
                    ' 
                    '  We want to try clustering them with different distance functions.
                    '  Distance function is chosen when we add dataset to the clusterizer.
                    '  We can choose several distance types - Euclidean, city block, Chebyshev,
                    '  several correlation measures or user-supplied distance matrix.
                    ' 
                    '  Here we'll try three distances: Euclidean, Pearson correlation,
                    '  user-supplied distance matrix. Different distance functions lead
                    '  to different choices being made by algorithm during clustering.
                    ' 
                    Dim s As clusterizerstate = New XAlglib.clusterizerstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As ahcreport = New XAlglib.ahcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim disttype As Integer
                    Dim xy(,) As Double = New Double(,){{1,2,1,2},{6,7,6,7},{7,6,7,6}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.clusterizercreate(s)

                    '  With Euclidean distance function (disttype=2) two closest points
                    '  are P1 and P2, thus:
                    '  * first, we merge P1 and P2 to form C3=[P1,P2]
                    '  * second, we merge P0 and C3 to form C4=[P0,P1,P2]
                    disttype = 2
                    xalglib.clusterizersetpoints(s, xy, disttype)
                    xalglib.clusterizerrunahc(s, rep)
                    _TestResult = _TestResult And doc_test_int_matrix(rep.z, New Integer(,){{1,2},{0,3}})

                    '  With Pearson correlation distance function (disttype=10) situation
                    '  is different - distance between P0 and P1 is zero, thus:
                    '  * first, we merge P0 and P1 to form C3=[P0,P1]
                    '  * second, we merge P2 and C3 to form C4=[P0,P1,P2]
                    disttype = 10
                    xalglib.clusterizersetpoints(s, xy, disttype)
                    xalglib.clusterizerrunahc(s, rep)
                    _TestResult = _TestResult And doc_test_int_matrix(rep.z, New Integer(,){{0,1},{2,3}})

                    '  Finally, we try clustering with user-supplied distance matrix:
                    '      [ 0 3 1 ]
                    '  P = [ 3 0 3 ], where P[i,j] = dist(Pi,Pj)
                    '      [ 1 3 0 ]
                    ' 
                    '  * first, we merge P0 and P2 to form C3=[P0,P2]
                    '  * second, we merge P1 and C3 to form C4=[P0,P1,P2]
                    Dim d(,) As Double = New Double(,){{0,3,1},{3,0,3},{1,3,0}}
                    xalglib.clusterizersetdistances(s, d, true)
                    xalglib.clusterizerrunahc(s, rep)
                    _TestResult = _TestResult And doc_test_int_matrix(rep.z, New Integer(,){{0,2},{1,3}})
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "clst_distance")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST clst_kclusters
            '      Obtaining K top clusters from clusterization tree
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  We have a set of points in 2D space:
                    '      (P0,P1,P2,P3,P4) = ((1,1),(1,2),(4,1),(2,3),(4,1.5))
                    ' 
                    '   |
                    '   |     P3
                    '   |
                    '   | P1          
                    '   |             P4
                    '   | P0          P2
                    '   |-------------------------
                    ' 
                    '  We perform Agglomerative Hierarchic Clusterization (AHC) and we want
                    '  to get top K clusters from clusterization tree for different K.
                    ' 
                    Dim s As clusterizerstate = New XAlglib.clusterizerstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As ahcreport = New XAlglib.ahcreport() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{1,1},{1,2},{4,1},{2,3},{4,1.5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    Dim cidx() As Integer = New Integer(){}
                    Dim cz() As Integer = New Integer(){}

                    xalglib.clusterizercreate(s)
                    xalglib.clusterizersetpoints(s, xy, 2)
                    xalglib.clusterizerrunahc(s, rep)

                    '  with K=5, every points is assigned to its own cluster:
                    '  C0=P0, C1=P1 and so on...
                    xalglib.clusterizergetkclusters(rep, 5, cidx, cz)
                    _TestResult = _TestResult And doc_test_int_vector(cidx, New Integer(){0,1,2,3,4})

                    '  with K=1 we have one large cluster C0=[P0,P1,P2,P3,P4,P5]
                    xalglib.clusterizergetkclusters(rep, 1, cidx, cz)
                    _TestResult = _TestResult And doc_test_int_vector(cidx, New Integer(){0,0,0,0,0})

                    '  with K=3 we have three clusters C0=[P3], C1=[P2,P4], C2=[P0,P1]
                    xalglib.clusterizergetkclusters(rep, 3, cidx, cz)
                    _TestResult = _TestResult And doc_test_int_vector(cidx, New Integer(){2,2,1,0,1})
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "clst_kclusters")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST randomforest_cls
            '      Simple classification with random forests
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple classification example: classify points (x,y) in 2D space
                    '  as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
                    '  has to find out it).
                    ' 
                    '  First, we have to create decision forest builder object, load dataset and
                    '  specify training settings. Our dataset is specified as matrix, which has
                    '  following format:
                    ' 
                    '      x0 y0 class0
                    '      x1 y1 class1
                    '      x2 y2 class2
                    '      ....
                    ' 
                    '  Here xi and yi can be any values (and in fact you can have any number of
                    '  independent variables), and classi MUST be integer number in [0,NClasses)
                    '  range. In our example we denote points with x>=0 as class #0, and
                    '  ones with negative xi as class #1.
                    ' 
                    '  NOTE: if you want to solve regression problem, specify NClasses=1. In
                    '        this case last column of xy can be any numeric value.
                    ' 
                    '  For the sake of simplicity, our example includes only 4-point dataset.
                    '  However, random forests are able to cope with extremely large datasets
                    '  having millions of examples.
                    ' 
                    Dim builder As decisionforestbuilder = New XAlglib.decisionforestbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim nvars As Integer = 2
                    Dim nclasses As Integer = 2
                    Dim npoints As Integer = 4
                    Dim xy(,) As Double = New Double(,){{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.dfbuildercreate(builder)
                    xalglib.dfbuildersetdataset(builder, xy, npoints, nvars, nclasses)

                    '  in our example we train decision forest using full sample - it allows us
                    '  to get zero classification error. However, in practical applications smaller
                    '  values are used: 50%, 25%, 5% or even less.
                    xalglib.dfbuildersetsubsampleratio(builder, 1.0)

                    '  we train random forest with just one tree; again, in real life situations
                    '  you typically need from 50 to 500 trees.
                    Dim ntrees As Integer = 1
                    Dim forest As decisionforest = New XAlglib.decisionforest() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As dfreport = New XAlglib.dfreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.dfbuilderbuildrandomforest(builder, ntrees, forest, rep)

                    '  with such settings (100% of the training set is used) you can expect
                    '  zero classification error. Beautiful results, but remember - in real life
                    '  you do not need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult And doc_test_real(rep.relclserror, 0.0000, 0.00005)

                    '  now, let's perform some simple processing with dfprocess()
                    Dim x() As Double = New Double(){+1,0}
                    Dim y() As Double = New Double(){}
                    xalglib.dfprocess(forest, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){+1,0}, 0.0005)

                    '  another option is to use dfprocess0() which returns just first component
                    '  of the output vector y. ideal for regression problems and binary classifiers.
                    Dim y0 As Double
                    y0 = xalglib.dfprocess0(forest, x)
                    _TestResult = _TestResult And doc_test_real(y0, 1.000, 0.0005)

                    '  finally, you can use dfclassify() which returns most probable class index (i.e. argmax y[i]).
                    Dim i As Integer
                    i = xalglib.dfclassify(forest, x)
                    _TestResult = _TestResult And doc_test_int(i, 0)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "randomforest_cls")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST randomforest_reg
            '      Simple classification with decision forest
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple regression example: model f(x,y)=x+y
                    ' 
                    '  First, we have to create DF builder object, load dataset and specify
                    '  training settings. Our dataset is specified as matrix, which has following
                    '  format:
                    ' 
                    '      x0 y0 f0
                    '      x1 y1 f1
                    '      x2 y2 f2
                    '      ....
                    ' 
                    '  Here xi and yi can be any values, and fi is a dependent function value.
                    ' 
                    '  NOTE: you can also solve classification problems with DF models, see
                    '        another example for this unit.
                    ' 
                    Dim builder As decisionforestbuilder = New XAlglib.decisionforestbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim nvars As Integer = 2
                    Dim nclasses As Integer = 1
                    Dim npoints As Integer = 4
                    Dim xy(,) As Double = New Double(,){{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.dfbuildercreate(builder)
                    xalglib.dfbuildersetdataset(builder, xy, npoints, nvars, nclasses)

                    '  in our example we train decision forest using full sample - it allows us
                    '  to get zero classification error. However, in practical applications smaller
                    '  values are used: 50%, 25%, 5% or even less.
                    xalglib.dfbuildersetsubsampleratio(builder, 1.0)

                    '  we train random forest with just one tree; again, in real life situations
                    '  you typically need from 50 to 500 trees.
                    Dim ntrees As Integer = 1
                    Dim model As decisionforest = New XAlglib.decisionforest() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As dfreport = New XAlglib.dfreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.dfbuilderbuildrandomforest(builder, ntrees, model, rep)

                    '  with such settings (full sample is used) you can expect zero RMS error on the
                    '  training set. Beautiful results, but remember - in real life you do not
                    '  need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult And doc_test_real(rep.rmserror, 0.0000, 0.00005)

                    '  now, let's perform some simple processing with dfprocess()
                    Dim x() As Double = New Double(){+1,+1}
                    Dim y() As Double = New Double(){}
                    xalglib.dfprocess(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){+2}, 0.0005)

                    '  another option is to use dfprocess0() which returns just first component
                    '  of the output vector y. ideal for regression problems and binary classifiers.
                    Dim y0 As Double
                    y0 = xalglib.dfprocess0(model, x)
                    _TestResult = _TestResult And doc_test_real(y0, 2.000, 0.0005)

                    '  there also exist another convenience function, dfclassify(),
                    '  but it does not work for regression problems - it always returns -1.
                    Dim i As Integer
                    i = xalglib.dfclassify(model, x)
                    _TestResult = _TestResult And doc_test_int(i, -1)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "randomforest_reg")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST knn_cls
            '      Simple classification with KNN model
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple classification example: classify points (x,y) in 2D space
                    '  as ones with x>=0 and ones with x<0 (y is ignored, but our classifier
                    '  has to find out it).
                    ' 
                    '  First, we have to create KNN builder object, load dataset and specify
                    '  training settings. Our dataset is specified as matrix, which has following
                    '  format:
                    ' 
                    '      x0 y0 class0
                    '      x1 y1 class1
                    '      x2 y2 class2
                    '      ....
                    ' 
                    '  Here xi and yi can be any values (and in fact you can have any number of
                    '  independent variables), and classi MUST be integer number in [0,NClasses)
                    '  range. In our example we denote points with x>=0 as class #0, and
                    '  ones with negative xi as class #1.
                    ' 
                    '  NOTE: if you want to solve regression problem, specify dataset in similar
                    '        format, but with dependent variable(s) instead of class labels. You
                    '        can have dataset with multiple dependent variables, by the way!
                    ' 
                    '  For the sake of simplicity, our example includes only 4-point dataset and
                    '  really simple K=1 nearest neighbor search. Industrial problems typically
                    '  need larger values of K.
                    ' 
                    Dim builder As knnbuilder = New XAlglib.knnbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim nvars As Integer = 2
                    Dim nclasses As Integer = 2
                    Dim npoints As Integer = 4
                    Dim xy(,) As Double = New Double(,){{1,1,0},{1,-1,0},{-1,1,1},{-1,-1,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.knnbuildercreate(builder)
                    xalglib.knnbuildersetdatasetcls(builder, xy, npoints, nvars, nclasses)

                    '  we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
                    Dim k As Integer = 1
                    Dim eps As Double = 0
                    Dim model As knnmodel = New XAlglib.knnmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As knnreport = New XAlglib.knnreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.knnbuilderbuildknnmodel(builder, k, eps, model, rep)

                    '  with such settings (k=1 is used) you can expect zero classification
                    '  error on training set. Beautiful results, but remember - in real life
                    '  you do not need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult And doc_test_real(rep.relclserror, 0.0000, 0.00005)

                    '  now, let's perform some simple processing with knnprocess()
                    Dim x() As Double = New Double(){+1,0}
                    Dim y() As Double = New Double(){}
                    xalglib.knnprocess(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){+1,0}, 0.0005)

                    '  another option is to use knnprocess0() which returns just first component
                    '  of the output vector y. ideal for regression problems and binary classifiers.
                    Dim y0 As Double
                    y0 = xalglib.knnprocess0(model, x)
                    _TestResult = _TestResult And doc_test_real(y0, 1.000, 0.0005)

                    '  finally, you can use knnclassify() which returns most probable class index (i.e. argmax y[i]).
                    Dim i As Integer
                    i = xalglib.knnclassify(model, x)
                    _TestResult = _TestResult And doc_test_int(i, 0)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "knn_cls")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST knn_reg
            '      Simple classification with KNN model
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  The very simple regression example: model f(x,y)=x+y
                    ' 
                    '  First, we have to create KNN builder object, load dataset and specify
                    '  training settings. Our dataset is specified as matrix, which has following
                    '  format:
                    ' 
                    '      x0 y0 f0
                    '      x1 y1 f1
                    '      x2 y2 f2
                    '      ....
                    ' 
                    '  Here xi and yi can be any values, and fi is a dependent function value.
                    '  By the way, with KNN algorithm you can even model functions with multiple
                    '  dependent variables!
                    ' 
                    '  NOTE: you can also solve classification problems with KNN models, see
                    '        another example for this unit.
                    ' 
                    '  For the sake of simplicity, our example includes only 4-point dataset and
                    '  really simple K=1 nearest neighbor search. Industrial problems typically
                    '  need larger values of K.
                    ' 
                    Dim builder As knnbuilder = New XAlglib.knnbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim nvars As Integer = 2
                    Dim nout As Integer = 1
                    Dim npoints As Integer = 4
                    Dim xy(,) As Double = New Double(,){{1,1,+2},{1,-1,0},{-1,1,0},{-1,-1,-2}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If

                    xalglib.knnbuildercreate(builder)
                    xalglib.knnbuildersetdatasetreg(builder, xy, npoints, nvars, nout)

                    '  we build KNN model with k=1 and eps=0 (exact k-nn search is performed)
                    Dim k As Integer = 1
                    Dim eps As Double = 0
                    Dim model As knnmodel = New XAlglib.knnmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As knnreport = New XAlglib.knnreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.knnbuilderbuildknnmodel(builder, k, eps, model, rep)

                    '  with such settings (k=1 is used) you can expect zero RMS error on the
                    '  training set. Beautiful results, but remember - in real life you do not
                    '  need zero TRAINING SET error, you need good generalization.

                    _TestResult = _TestResult And doc_test_real(rep.rmserror, 0.0000, 0.00005)

                    '  now, let's perform some simple processing with knnprocess()
                    Dim x() As Double = New Double(){+1,+1}
                    Dim y() As Double = New Double(){}
                    xalglib.knnprocess(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){+2}, 0.0005)

                    '  another option is to use knnprocess0() which returns just first component
                    '  of the output vector y. ideal for regression problems and binary classifiers.
                    Dim y0 As Double
                    y0 = xalglib.knnprocess0(model, x)
                    _TestResult = _TestResult And doc_test_real(y0, 2.000, 0.0005)

                    '  there also exist another convenience function, knnclassify(),
                    '  but it does not work for regression problems - it always returns -1.
                    Dim i As Integer
                    i = xalglib.knnclassify(model, x)
                    _TestResult = _TestResult And doc_test_int(i, -1)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "knn_reg")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST autogk_d1
            '      Integrating f=exp(x) by adaptive integrator
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  This example demonstrates integration of f=exp(x) on [0,1]:
                    '  * first, autogkstate is initialized
                    '  * then we call integration function
                    '  * and finally we obtain results with autogkresults() call
                    ' 
                    Dim a As Double = 0
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = 1
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim s As autogkstate = New XAlglib.autogkstate() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    Dim rep As autogkreport = New XAlglib.autogkreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.autogksmooth(a, b, s)
                    xalglib.autogkintegrate(s, AddressOf int_function_1_func, Nothing)
                    xalglib.autogkresults(s, v, rep)

                    _TestResult = _TestResult And doc_test_real(v, 1.7182, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "autogk_d1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST fft_complex_d1
            '      Complex FFT: simple example
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  first we demonstrate forward FFT:
                    '  [1i,1i,1i,1i] is converted to [4i, 0, 0, 0]
                    ' 
                    Dim z() As alglib.complex = New alglib.complex(){New alglib.complex(0,1),New alglib.complex(0,1),New alglib.complex(0,1),New alglib.complex(0,1)}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    xalglib.fftc1d(z)
                    _TestResult = _TestResult And doc_test_complex_vector(z, New alglib.complex(){New alglib.complex(0,4),0,0,0}, 0.0001)

                    ' 
                    '  now we convert [4i, 0, 0, 0] back to [1i,1i,1i,1i]
                    '  with backward FFT
                    ' 
                    xalglib.fftc1dinv(z)
                    _TestResult = _TestResult And doc_test_complex_vector(z, New alglib.complex(){New alglib.complex(0,1),New alglib.complex(0,1),New alglib.complex(0,1),New alglib.complex(0,1)}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST fft_complex_d2
            '      Complex FFT: advanced example
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  first we demonstrate forward FFT:
                    '  [0,1,0,1i] is converted to [1+1i, -1-1i, -1-1i, 1+1i]
                    ' 
                    Dim z() As alglib.complex = New alglib.complex(){0,1,0,New alglib.complex(0,1)}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    xalglib.fftc1d(z)
                    _TestResult = _TestResult And doc_test_complex_vector(z, New alglib.complex(){New alglib.complex(1,+1),New alglib.complex(-1,-1),New alglib.complex(-1,-1),New alglib.complex(1,+1)}, 0.0001)

                    ' 
                    '  now we convert result back with backward FFT
                    ' 
                    xalglib.fftc1dinv(z)
                    _TestResult = _TestResult And doc_test_complex_vector(z, New alglib.complex(){0,1,0,New alglib.complex(0,1)}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_d2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST fft_real_d1
            '      Real FFT: simple example
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  first we demonstrate forward FFT:
                    '  [1,1,1,1] is converted to [4, 0, 0, 0]
                    ' 
                    Dim x() As Double = New Double(){1,1,1,1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim f() As alglib.complex = New alglib.complex(){}
                    Dim x2() As Double = New Double(){}
                    xalglib.fftr1d(x, f)
                    _TestResult = _TestResult And doc_test_complex_vector(f, New alglib.complex(){4,0,0,0}, 0.0001)

                    ' 
                    '  now we convert [4, 0, 0, 0] back to [1,1,1,1]
                    '  with backward FFT
                    ' 
                    xalglib.fftr1dinv(f, x2)
                    _TestResult = _TestResult And doc_test_real_vector(x2, New Double(){1,1,1,1}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST fft_real_d2
            '      Real FFT: advanced example
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  first we demonstrate forward FFT:
                    '  [1,2,3,4] is converted to [10, -2+2i, -2, -2-2i]
                    ' 
                    '  note that output array is self-adjoint:
                    '  * f[0] = conj(f[0])
                    '  * f[1] = conj(f[3])
                    '  * f[2] = conj(f[2])
                    ' 
                    Dim x() As Double = New Double(){1,2,3,4}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    Dim f() As alglib.complex = New alglib.complex(){}
                    Dim x2() As Double = New Double(){}
                    xalglib.fftr1d(x, f)
                    _TestResult = _TestResult And doc_test_complex_vector(f, New alglib.complex(){10,New alglib.complex(-2,+2),-2,New alglib.complex(-2,-2)}, 0.0001)

                    ' 
                    '  now we convert [10, -2+2i, -2, -2-2i] back to [1,2,3,4]
                    ' 
                    xalglib.fftr1dinv(f, x2)
                    _TestResult = _TestResult And doc_test_real_vector(x2, New Double(){1,2,3,4}, 0.0001)

                    ' 
                    '  remember that F is self-adjoint? It means that we can pass just half
                    '  (slightly larger than half) of F to inverse real FFT and still get our result.
                    ' 
                    '  I.e. instead [10, -2+2i, -2, -2-2i] we pass just [10, -2+2i, -2] and everything works!
                    ' 
                    '  NOTE: in this case we should explicitly pass array length (which is 4) to ALGLIB;
                    '  if not, it will automatically use array length to determine FFT size and
                    '  will erroneously make half-length FFT.
                    ' 
                    f = New alglib.complex(){10,New alglib.complex(-2,+2),-2}
                    xalglib.fftr1dinv(f, 4, x2)
                    _TestResult = _TestResult And doc_test_real_vector(x2, New Double(){1,2,3,4}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "fft_real_d2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST fft_complex_e1
            '      error detection in backward FFT
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    Dim z() As alglib.complex = New alglib.complex(){0,2,0,-2}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_vector_by_value(z, v_spoil_complex)
                    End If
                    xalglib.fftc1dinv(z)
                    _TestResult = _TestResult And doc_test_complex_vector(z, New alglib.complex(){0,New alglib.complex(0,1),0,New alglib.complex(0,-1)}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "fft_complex_e1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST idw_d_mstab
            '      Simple model built with IDW-MSTAB algorithm
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example illustrates basic concepts of the IDW models:
                    '  creation and evaluation.
                    '  
                    '  Suppose that we have set of 2-dimensional points with associated
                    '  scalar function values, and we want to build an IDW model using
                    '  our data.
                    '  
                    '  NOTE: we can work with N-dimensional models and vector-valued functions too :)
                    '  
                    '  Typical sequence of steps is given below:
                    '  1. we create IDW builder object
                    '  2. we attach our dataset to the IDW builder and tune algorithm settings
                    '  3. we generate IDW model
                    '  4. we use IDW model instance (evaluate, serialize, etc.)
                    ' 
                    Dim v As Double

                    ' 
                    '  Step 1: IDW builder creation.
                    ' 
                    '  We have to specify dimensionality of the space (2 or 3) and
                    '  dimensionality of the function (scalar or vector).
                    ' 
                    '  New builder object is empty - it has not dataset and uses
                    '  default model construction settings
                    ' 
                    Dim builder As idwbuilder = New XAlglib.idwbuilder() ' initializer can be dropped, but compiler will issue warning
                    xalglib.idwbuildercreate(2, 1, builder)

                    ' 
                    '  Step 2: dataset addition
                    ' 
                    '  XY contains two points - x0=(-1,0) and x1=(+1,0) -
                    '  and two function values f(x0)=2, f(x1)=3.
                    ' 
                    Dim xy(,) As Double = New Double(,){{-1,0,2},{+1,0,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.idwbuildersetpoints(builder, xy)

                    ' 
                    '  Step 3: choose IDW algorithm and generate model
                    ' 
                    '  We use modified stabilized IDW algorithm with following parameters:
                    '  * SRad - set to 5.0 (search radius must be large enough)
                    ' 
                    '  IDW-MSTAB algorithm is a state-of-the-art implementation of IDW which
                    '  is competitive with RBFs and bicubic splines. See comments on the
                    '  idwbuildersetalgomstab() function for more information.
                    ' 
                    Dim model As idwmodel = New XAlglib.idwmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As idwreport = New XAlglib.idwreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.idwbuildersetalgomstab(builder, 5.0)
                    xalglib.idwfit(builder, model, rep)

                    ' 
                    '  Step 4: model was built, evaluate its value
                    ' 
                    v = xalglib.idwcalc2(model, 1.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 3.000, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "idw_d_mstab")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST idw_d_serialize
            '      IDW model serialization/unserialization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example shows how to serialize and unserialize IDW model.
                    '  
                    '  Suppose that we have set of 2-dimensional points with associated
                    '  scalar function values, and we have built an IDW model using
                    '  our data.
                    ' 
                    '  This model can be serialized to string or stream. ALGLIB supports
                    '  flexible (un)serialization, i.e. you can move serialized model
                    '  representation between different machines (32-bit or 64-bit),
                    '  different CPU architectures (x86/64, ARM) or even different
                    '  programming languages supported by ALGLIB (C#, C++, ...).
                    ' 
                    '  Our first step is to build model, evaluate it at point (1,0),
                    '  and serialize it to string.
                    ' 
                    Dim s As String = "" ' initializer can be dropped, but compiler will issue warning 
                    Dim v As Double
                    Dim xy(,) As Double = New Double(,){{-1,0,2},{+1,0,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    Dim builder As idwbuilder = New XAlglib.idwbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim model As idwmodel = New XAlglib.idwmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim model2 As idwmodel = New XAlglib.idwmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As idwreport = New XAlglib.idwreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.idwbuildercreate(2, 1, builder)
                    xalglib.idwbuildersetpoints(builder, xy)
                    xalglib.idwbuildersetalgomstab(builder, 5.0)
                    xalglib.idwfit(builder, model, rep)
                    v = xalglib.idwcalc2(model, 1.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 3.000, 0.005)

                    ' 
                    '  Serialization + unserialization to a different instance
                    '  of the model class.
                    ' 
                    xalglib.idwserialize(model, s)
                    xalglib.idwunserialize(s, model2)

                    ' 
                    '  Evaluate unserialized model at the same point
                    ' 
                    v = xalglib.idwcalc2(model2, 1.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 3.000, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "idw_d_serialize")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline1d_d_linear
            '      Piecewise linear spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    ' 
                    '  We use piecewise linear spline to interpolate f(x)=x^2 sampled 
                    '  at 5 equidistant nodes on [-1,+1].
                    ' 
                    Dim x() As Double = New Double(){-1.0,-0.5,0.0,+0.5,+1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){+1.0,0.25,0.0,0.25,+1.0}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = 0.25
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    Dim s As spline1dinterpolant = New XAlglib.spline1dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline1dbuildlinear(x, y, s)

                    '  calculate S(0.25) - it is quite different from 0.25^2=0.0625
                    v = xalglib.spline1dcalc(s, t)
                    _TestResult = _TestResult And doc_test_real(v, 0.125, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_linear")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline1d_d_cubic
            '      Cubic spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    ' 
                    '  We use cubic spline to interpolate f(x)=x^2 sampled 
                    '  at 5 equidistant nodes on [-1,+1].
                    ' 
                    '  First, we use default boundary conditions ("parabolically terminated
                    '  spline") because cubic spline built with such boundary conditions 
                    '  will exactly reproduce any quadratic f(x).
                    ' 
                    '  Then we try to use natural boundary conditions
                    '      d2S(-1)/dx^2 = 0.0
                    '      d2S(+1)/dx^2 = 0.0
                    '  and see that such spline interpolated f(x) with small error.
                    ' 
                    Dim x() As Double = New Double(){-1.0,-0.5,0.0,+0.5,+1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){+1.0,0.25,0.0,0.25,+1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = 0.25
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    Dim s As spline1dinterpolant = New XAlglib.spline1dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim natural_bound_type As Integer = 2
                    ' 
                    '  Test exact boundary conditions: build S(x), calculare S(0.25)
                    '  (almost same as original function)
                    ' 
                    xalglib.spline1dbuildcubic(x, y, s)
                    v = xalglib.spline1dcalc(s, t)
                    _TestResult = _TestResult And doc_test_real(v, 0.0625, 0.00001)

                    ' 
                    '  Test natural boundary conditions: build S(x), calculare S(0.25)
                    '  (small interpolation error)
                    ' 
                    xalglib.spline1dbuildcubic(x, y, 5, natural_bound_type, 0.0, natural_bound_type, 0.0, s)
                    v = xalglib.spline1dcalc(s, t)
                    _TestResult = _TestResult And doc_test_real(v, 0.0580, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_cubic")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline1d_d_monotone
            '      Monotone interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    ' 
                    '  Spline built witn spline1dbuildcubic() can be non-monotone even when
                    '  Y-values form monotone sequence. Say, for x=[0,1,2] and y=[0,1,1]
                    '  cubic spline will monotonically grow until x=1.5 and then start
                    '  decreasing.
                    ' 
                    '  That's why ALGLIB provides special spline construction function
                    '  which builds spline which preserves monotonicity of the original
                    '  dataset.
                    ' 
                    '  NOTE: in case original dataset is non-monotonic, ALGLIB splits it
                    '  into monotone subsequences and builds piecewise monotonic spline.
                    ' 
                    Dim x() As Double = New Double(){0,1,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0,1,1}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim s As spline1dinterpolant = New XAlglib.spline1dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline1dbuildmonotone(x, y, s)

                    '  calculate S at x = [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
                    '  you may see that spline is really monotonic
                    Dim v As Double
                    v = xalglib.spline1dcalc(s, -0.5)
                    _TestResult = _TestResult And doc_test_real(v, 0.0000, 0.00005)
                    v = xalglib.spline1dcalc(s, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.0000, 0.00005)
                    v = xalglib.spline1dcalc(s, +0.5)
                    _TestResult = _TestResult And doc_test_real(v, 0.5000, 0.00005)
                    v = xalglib.spline1dcalc(s, 1.0)
                    _TestResult = _TestResult And doc_test_real(v, 1.0000, 0.00005)
                    v = xalglib.spline1dcalc(s, 1.5)
                    _TestResult = _TestResult And doc_test_real(v, 1.0000, 0.00005)
                    v = xalglib.spline1dcalc(s, 2.0)
                    _TestResult = _TestResult And doc_test_real(v, 1.0000, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_monotone")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline1d_d_griddiff
            '      Differentiation on the grid using cubic splines
            '
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    ' 
                    '  We use cubic spline to do grid differentiation, i.e. having
                    '  values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
                    '  we calculate derivatives of cubic spline at nodes WITHOUT
                    '  CONSTRUCTION OF SPLINE OBJECT.
                    ' 
                    '  There are efficient functions spline1dgriddiffcubic() and
                    '  spline1dgriddiff2cubic() for such calculations.
                    ' 
                    '  We use default boundary conditions ("parabolically terminated
                    '  spline") because cubic spline built with such boundary conditions 
                    '  will exactly reproduce any quadratic f(x).
                    ' 
                    '  Actually, we could use natural conditions, but we feel that 
                    '  spline which exactly reproduces f() will show us more 
                    '  understandable results.
                    ' 
                    Dim x() As Double = New Double(){-1.0,-0.5,0.0,+0.5,+1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){+1.0,0.25,0.0,0.25,+1.0}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim d1() As Double = New Double(){}
                    Dim d2() As Double = New Double(){}

                    ' 
                    '  We calculate first derivatives: they must be equal to 2*x
                    ' 
                    xalglib.spline1dgriddiffcubic(x, y, d1)
                    _TestResult = _TestResult And doc_test_real_vector(d1, New Double(){-2.0,-1.0,0.0,+1.0,+2.0}, 0.0001)

                    ' 
                    '  Now test griddiff2, which returns first AND second derivatives.
                    '  First derivative is 2*x, second is equal to 2.0
                    ' 
                    xalglib.spline1dgriddiff2cubic(x, y, d1, d2)
                    _TestResult = _TestResult And doc_test_real_vector(d1, New Double(){-2.0,-1.0,0.0,+1.0,+2.0}, 0.0001)
                    _TestResult = _TestResult And doc_test_real_vector(d2, New Double(){2.0,2.0,2.0,2.0,2.0}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_griddiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline1d_d_convdiff
            '      Resampling using cubic splines
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    ' 
                    '  We use cubic spline to do resampling, i.e. having
                    '  values of f(x)=x^2 sampled at 5 equidistant nodes on [-1,+1]
                    '  we calculate values/derivatives of cubic spline on 
                    '  another grid (equidistant with 9 nodes on [-1,+1])
                    '  WITHOUT CONSTRUCTION OF SPLINE OBJECT.
                    ' 
                    '  There are efficient functions spline1dconvcubic(),
                    '  spline1dconvdiffcubic() and spline1dconvdiff2cubic() 
                    '  for such calculations.
                    ' 
                    '  We use default boundary conditions ("parabolically terminated
                    '  spline") because cubic spline built with such boundary conditions 
                    '  will exactly reproduce any quadratic f(x).
                    ' 
                    '  Actually, we could use natural conditions, but we feel that 
                    '  spline which exactly reproduces f() will show us more 
                    '  understandable results.
                    ' 
                    Dim x_old() As Double = New Double(){-1.0,-0.5,0.0,+0.5,+1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x_old, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x_old, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x_old, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x_old)
                    End If
                    Dim y_old() As Double = New Double(){+1.0,0.25,0.0,0.25,+1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y_old, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y_old, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y_old, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y_old)
                    End If
                    Dim x_new() As Double = New Double(){-1.00,-0.75,-0.50,-0.25,0.00,+0.25,+0.50,+0.75,+1.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x_new, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x_new, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x_new, v_spoil_real)
                    End If
                    Dim y_new() As Double = New Double(){}
                    Dim d1_new() As Double = New Double(){}
                    Dim d2_new() As Double = New Double(){}

                    ' 
                    '  First, conversion without differentiation.
                    ' 
                    ' 
                    xalglib.spline1dconvcubic(x_old, y_old, x_new, y_new)
                    _TestResult = _TestResult And doc_test_real_vector(y_new, New Double(){1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001)

                    ' 
                    '  Then, conversion with differentiation (first derivatives only)
                    ' 
                    ' 
                    xalglib.spline1dconvdiffcubic(x_old, y_old, x_new, y_new, d1_new)
                    _TestResult = _TestResult And doc_test_real_vector(y_new, New Double(){1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001)
                    _TestResult = _TestResult And doc_test_real_vector(d1_new, New Double(){-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0}, 0.0001)

                    ' 
                    '  Finally, conversion with first and second derivatives
                    ' 
                    ' 
                    xalglib.spline1dconvdiff2cubic(x_old, y_old, x_new, y_new, d1_new, d2_new)
                    _TestResult = _TestResult And doc_test_real_vector(y_new, New Double(){1.0000,0.5625,0.2500,0.0625,0.0000,0.0625,0.2500,0.5625,1.0000}, 0.0001)
                    _TestResult = _TestResult And doc_test_real_vector(d1_new, New Double(){-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0}, 0.0001)
                    _TestResult = _TestResult And doc_test_real_vector(d2_new, New Double(){2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0}, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline1d_d_convdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST parametric_rdp
            '      Parametric Ramer-Douglas-Peucker approximation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    ' 
                    '  We use RDP algorithm to approximate parametric 2D curve given by
                    '  locations in t=0,1,2,3 (see below), which form piecewise linear
                    '  trajectory through D-dimensional space (2-dimensional in our example).
                    '  
                    '      |
                    '      |
                    '      -     *     *     X2................X3
                    '      |                .
                    '      |               .
                    '      -     *     *  .  *     *     *     *
                    '      |             .
                    '      |            .
                    '      -     *     X1    *     *     *     *
                    '      |      .....
                    '      |  ....
                    '      X0----|-----|-----|-----|-----|-----|---
                    ' 
                    Dim npoints As Integer = 4
                    Dim ndimensions As Integer = 2
                    Dim x(,) As Double = New Double(,){{0,0},{2,1},{3,3},{6,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If

                    ' 
                    '  Approximation of parametric curve is performed by another parametric curve
                    '  with lesser amount of points. It allows to work with "compressed"
                    '  representation, which needs smaller amount of memory. Say, in our example
                    '  (we allow points with error smaller than 0.8) approximation will have
                    '  just two sequential sections connecting X0 with X2, and X2 with X3.
                    '  
                    '      |
                    '      |
                    '      -     *     *     X2................X3
                    '      |               . 
                    '      |             .  
                    '      -     *     .     *     *     *     *
                    '      |         .    
                    '      |       .     
                    '      -     .     X1    *     *     *     *
                    '      |   .       
                    '      | .    
                    '      X0----|-----|-----|-----|-----|-----|---
                    ' 
                    ' 
                    Dim y(,) As Double = New Double(,){{}}
                    Dim idxy() As Integer = New Integer(){}
                    Dim nsections As Integer
                    Dim limitcnt As Integer = 0
                    Dim limiteps As Double = 0.8
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        limiteps = v_spoil_real
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        limiteps = v_spoil_real
                    End If
                    xalglib.parametricrdpfixed(x, npoints, ndimensions, limitcnt, limiteps, y, idxy, nsections)
                    _TestResult = _TestResult And doc_test_int(nsections, 2)
                    _TestResult = _TestResult And doc_test_int_vector(idxy, New Integer(){0,2,3})
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "parametric_rdp")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline3d_trilinear
            '      Trilinear spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 21
                Try
                    ' 
                    '  We use trilinear spline to interpolate f(x,y,z)=x+xy+z sampled 
                    '  at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
                    ' 
                    '  We store x, y and z-values at local arrays with same names.
                    '  Function values are stored in the array F as follows:
                    '      f[0]     (x,y,z) = (0,0,0)
                    '      f[1]     (x,y,z) = (1,0,0)
                    '      f[2]     (x,y,z) = (0,1,0)
                    '      f[3]     (x,y,z) = (1,1,0)
                    '      f[4]     (x,y,z) = (0,0,1)
                    '      f[5]     (x,y,z) = (1,0,1)
                    '      f[6]     (x,y,z) = (0,1,1)
                    '      f[7]     (x,y,z) = (1,1,1)
                    ' 
                    Dim x() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim z() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(z)
                    End If
                    Dim f() As Double = New Double(){0,1,0,2,1,2,1,3}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim vx As Double = 0.50
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        vx = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        vx = v_spoil_real
                    End If
                    Dim vy As Double = 0.50
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        vy = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        vy = v_spoil_real
                    End If
                    Dim vz As Double = 0.50
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.PositiveInfinity
                        vz = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.NegativeInfinity
                        vz = v_spoil_real
                    End If
                    Dim v As Double
                    Dim s As spline3dinterpolant = New XAlglib.spline3dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline3dbuildtrilinearv(x, 2, y, 2, z, 2, f, 1, s)

                    '  calculate S(0.5,0.5,0.5)
                    v = xalglib.spline3dcalc(s, vx, vy, vz)
                    _TestResult = _TestResult And doc_test_real(v, 1.2500, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline3d_trilinear")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline3d_vector
            '      Vector-valued trilinear spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 21
                Try
                    ' 
                    '  We use trilinear vector-valued spline to interpolate {f0,f1}={x+xy+z,x+xy+yz+z}
                    '  sampled at (x,y,z) from [0.0, 1.0] X [0.0, 1.0] X [0.0, 1.0].
                    ' 
                    '  We store x, y and z-values at local arrays with same names.
                    '  Function values are stored in the array F as follows:
                    '      f[0]     f0, (x,y,z) = (0,0,0)
                    '      f[1]     f1, (x,y,z) = (0,0,0)
                    '      f[2]     f0, (x,y,z) = (1,0,0)
                    '      f[3]     f1, (x,y,z) = (1,0,0)
                    '      f[4]     f0, (x,y,z) = (0,1,0)
                    '      f[5]     f1, (x,y,z) = (0,1,0)
                    '      f[6]     f0, (x,y,z) = (1,1,0)
                    '      f[7]     f1, (x,y,z) = (1,1,0)
                    '      f[8]     f0, (x,y,z) = (0,0,1)
                    '      f[9]     f1, (x,y,z) = (0,0,1)
                    '      f[10]    f0, (x,y,z) = (1,0,1)
                    '      f[11]    f1, (x,y,z) = (1,0,1)
                    '      f[12]    f0, (x,y,z) = (0,1,1)
                    '      f[13]    f1, (x,y,z) = (0,1,1)
                    '      f[14]    f0, (x,y,z) = (1,1,1)
                    '      f[15]    f1, (x,y,z) = (1,1,1)
                    ' 
                    Dim x() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim z() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(z, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(z)
                    End If
                    Dim f() As Double = New Double(){0,0,1,1,0,0,2,2,1,1,2,2,1,2,3,4}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim vx As Double = 0.50
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        vx = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        vx = v_spoil_real
                    End If
                    Dim vy As Double = 0.50
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        vy = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        vy = v_spoil_real
                    End If
                    Dim vz As Double = 0.50
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.PositiveInfinity
                        vz = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.NegativeInfinity
                        vz = v_spoil_real
                    End If
                    Dim s As spline3dinterpolant = New XAlglib.spline3dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline3dbuildtrilinearv(x, 2, y, 2, z, 2, f, 2, s)

                    '  calculate S(0.5,0.5,0.5) - we have vector of values instead of single value
                    Dim v() As Double = New Double(){}
                    xalglib.spline3dcalcv(s, vx, vy, vz, v)
                    _TestResult = _TestResult And doc_test_real_vector(v, New Double(){1.2500,1.5000}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline3d_vector")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_d_calcdiff
            '      Interpolation and differentiation using barycentric representation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    ' 
                    '  Here we demonstrate polynomial interpolation and differentiation
                    '  of y=x^2-x sampled at [0,1,2]. Barycentric representation of polynomial is used.
                    ' 
                    Dim x() As Double = New Double(){0,1,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    Dim dv As Double
                    Dim d2v As Double
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  barycentric model is created
                    xalglib.polynomialbuild(x, y, p)

                    '  barycentric interpolation is demonstrated
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)

                    '  barycentric differentation is demonstrated
                    xalglib.barycentricdiff1(p, t, v, dv)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And doc_test_real(dv, -3.0, 0.00005)

                    '  second derivatives with barycentric representation
                    xalglib.barycentricdiff1(p, t, v, dv)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And doc_test_real(dv, -3.0, 0.00005)
                    xalglib.barycentricdiff2(p, t, v, dv, d2v)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And doc_test_real(dv, -3.0, 0.00005)
                    _TestResult = _TestResult And doc_test_real(d2v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_calcdiff")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_d_conv
            '      Conversion between power basis and barycentric representation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    ' 
                    '  Here we demonstrate conversion of y=x^2-x
                    '  between power basis and barycentric representation.
                    ' 
                    Dim a() As Double = New Double(){0,-1,+1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(a, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(a, v_spoil_real)
                    End If
                    Dim t As Double = 2
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a2() As Double = New Double(){}
                    Dim v As Double
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  a=[0,-1,+1] is decomposition of y=x^2-x in the power basis:
                    ' 
                    '      y = 0 - 1*x + 1*x^2
                    ' 
                    '  We convert it to the barycentric form.
                    ' 
                    xalglib.polynomialpow2bar(a, p)

                    '  now we have barycentric interpolation; we can use it for interpolation
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.005)

                    '  we can also convert back from barycentric representation to power basis
                    xalglib.polynomialbar2pow(p, a2)
                    _TestResult = _TestResult And doc_test_real_vector(a2, New Double(){0,-1,+1}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_conv")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_d_spec
            '      Polynomial interpolation on special grids (equidistant, Chebyshev I/II)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    ' 
                    '  Temporaries:
                    '  * values of y=x^2-x sampled at three special grids:
                    '    * equdistant grid spanning [0,2],     x[i] = 2*i/(N-1), i=0..N-1
                    '    * Chebyshev-I grid spanning [-1,+1],  x[i] = 1 + Cos(PI*(2*i+1)/(2*n)), i=0..N-1
                    '    * Chebyshev-II grid spanning [-1,+1], x[i] = 1 + Cos(PI*i/(n-1)), i=0..N-1
                    '  * barycentric interpolants for these three grids
                    '  * vectors to store coefficients of quadratic representation
                    ' 
                    Dim y_eqdist() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y_eqdist, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y_eqdist, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y_eqdist, v_spoil_real)
                    End If
                    Dim y_cheb1() As Double = New Double(){-0.116025,0.000000,1.616025}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y_cheb1, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y_cheb1, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y_cheb1, v_spoil_real)
                    End If
                    Dim y_cheb2() As Double = New Double(){0,0,2}
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y_cheb2, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y_cheb2, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y_cheb2, v_spoil_real)
                    End If
                    Dim p_eqdist As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim p_cheb1 As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim p_cheb2 As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim a_eqdist() As Double = New Double(){}
                    Dim a_cheb1() As Double = New Double(){}
                    Dim a_cheb2() As Double = New Double(){}

                    ' 
                    '  First, we demonstrate construction of barycentric interpolants on
                    '  special grids. We unpack power representation to ensure that
                    '  interpolant was built correctly.
                    ' 
                    '  In all three cases we should get same quadratic function.
                    ' 
                    xalglib.polynomialbuildeqdist(0.0, 2.0, y_eqdist, p_eqdist)
                    xalglib.polynomialbar2pow(p_eqdist, a_eqdist)
                    _TestResult = _TestResult And doc_test_real_vector(a_eqdist, New Double(){0,-1,+1}, 0.00005)

                    xalglib.polynomialbuildcheb1(-1, +1, y_cheb1, p_cheb1)
                    xalglib.polynomialbar2pow(p_cheb1, a_cheb1)
                    _TestResult = _TestResult And doc_test_real_vector(a_cheb1, New Double(){0,-1,+1}, 0.00005)

                    xalglib.polynomialbuildcheb2(-1, +1, y_cheb2, p_cheb2)
                    xalglib.polynomialbar2pow(p_cheb2, a_cheb2)
                    _TestResult = _TestResult And doc_test_real_vector(a_cheb2, New Double(){0,-1,+1}, 0.00005)

                    ' 
                    '  Now we demonstrate polynomial interpolation without construction 
                    '  of the barycentricinterpolant structure.
                    ' 
                    '  We calculate interpolant value at x=-2.
                    '  In all three cases we should get same f=6
                    ' 
                    Dim t As Double = -2
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalceqdist(0.0, 2.0, y_eqdist, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)

                    v = xalglib.polynomialcalccheb1(-1, +1, y_cheb1, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)

                    v = xalglib.polynomialcalccheb2(-1, +1, y_cheb2, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_d_spec")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_1
            '      Polynomial interpolation, full list of parameters.
            '
            System.Console.WriteLine("100/151")
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    Dim x() As Double = New Double(){0,1,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuild(x, y, 3, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_2
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildeqdist(0.0, 2.0, y, 3, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_3
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim y() As Double = New Double(){-0.116025,0.000000,1.616025}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildcheb1(-1.0, +1.0, y, 3, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_3")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_4
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -2
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildcheb2(a, b, y, 3, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_4")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_5
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalceqdist(0.0, 2.0, y, 3, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_5")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_6
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    Dim y() As Double = New Double(){-0.116025,0.000000,1.616025}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalccheb1(a, b, y, 3, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_6")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_7
            '      Polynomial interpolation, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim t As Double = -2
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalccheb2(a, b, y, 3, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_7")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_8
            '      Polynomial interpolation: y=x^2-x, equidistant grid, barycentric form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildeqdist(0.0, 2.0, y, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_8")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_9
            '      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind), barycentric form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    Dim y() As Double = New Double(){-0.116025,0.000000,1.616025}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildcheb1(a, b, y, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_9")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_10
            '      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind), barycentric form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -2
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialbuildcheb2(a, b, y, p)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_10")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_11
            '      Polynomial interpolation: y=x^2-x, equidistant grid
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalceqdist(0.0, 2.0, y, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_11")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_12
            '      Polynomial interpolation: y=x^2-x, Chebyshev grid (first kind)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    Dim y() As Double = New Double(){-0.116025,0.000000,1.616025}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -1
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalccheb1(a, b, y, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_12")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST polint_t_13
            '      Polynomial interpolation: y=x^2-x, Chebyshev grid (second kind)
            '
            _TestResult = true
            For _spoil_scenario = -1 To 10
                Try
                    Dim y() As Double = New Double(){0,0,2}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    Dim t As Double = -2
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim a As Double = -1
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        a = v_spoil_real
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        a = v_spoil_real
                    End If
                    Dim b As Double = +1
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        b = v_spoil_real
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        b = v_spoil_real
                    End If
                    Dim v As Double
                    v = xalglib.polynomialcalccheb2(a, b, y, t)
                    _TestResult = _TestResult And doc_test_real(v, 6.0, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "polint_t_13")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_nlf
            '      Nonlinear fitting using function value only
            '
            _TestResult = true
            For _spoil_scenario = -1 To 23
                Try
                    ' 
                    '  In this example we demonstrate exponential fitting
                    '  by f(x) = exp(-c*x^2)
                    '  using function value only.
                    ' 
                    '  Gradient is estimated using combination of numerical differences
                    '  and secant updates. diffstep variable stores differentiation step 
                    '  (we have to tell algorithm what step to use).
                    ' 
                    Dim x(,) As Double = New Double(,){{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If
                    Dim y() As Double = New Double(){0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim c() As Double = New Double(){0.3}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim info As Integer
                    Dim state As lsfitstate = New XAlglib.lsfitstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim diffstep As Double = 0.0001
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If

                    ' 
                    '  Fitting without weights
                    ' 
                    xalglib.lsfitcreatef(x, y, c, diffstep, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)

                    ' 
                    '  Fitting with weights
                    '  (you can change weights and see how it changes result)
                    ' 
                    Dim w() As Double = New Double(){1,1,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=23 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    xalglib.lsfitcreatewf(x, y, w, c, diffstep, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlf")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_nlfg
            '      Nonlinear fitting using gradient
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  In this example we demonstrate exponential fitting
                    '  by f(x) = exp(-c*x^2)
                    '  using function value and gradient (with respect to c).
                    ' 
                    Dim x(,) As Double = New Double(,){{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If
                    Dim y() As Double = New Double(){0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim c() As Double = New Double(){0.3}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim info As Integer
                    Dim state As lsfitstate = New XAlglib.lsfitstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Fitting without weights
                    ' 
                    xalglib.lsfitcreatefg(x, y, c, true, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, AddressOf function_cx_1_grad, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)

                    ' 
                    '  Fitting with weights
                    '  (you can change weights and see how it changes result)
                    ' 
                    Dim w() As Double = New Double(){1,1,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    xalglib.lsfitcreatewfg(x, y, w, c, true, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, AddressOf function_cx_1_grad, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfg")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_nlfgh
            '      Nonlinear fitting using gradient and Hessian
            '
            _TestResult = true
            For _spoil_scenario = -1 To 20
                Try
                    ' 
                    '  In this example we demonstrate exponential fitting
                    '  by f(x) = exp(-c*x^2)
                    '  using function value, gradient and Hessian (with respect to c)
                    ' 
                    Dim x(,) As Double = New Double(,){{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If
                    Dim y() As Double = New Double(){0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim c() As Double = New Double(){0.3}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim info As Integer
                    Dim state As lsfitstate = New XAlglib.lsfitstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Fitting without weights
                    ' 
                    xalglib.lsfitcreatefgh(x, y, c, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, AddressOf function_cx_1_grad, AddressOf function_cx_1_hess, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)

                    ' 
                    '  Fitting with weights
                    '  (you can change weights and see how it changes result)
                    ' 
                    Dim w() As Double = New Double(){1,1,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=20 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    xalglib.lsfitcreatewfgh(x, y, w, c, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, AddressOf function_cx_1_grad, AddressOf function_cx_1_hess, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.5}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfgh")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_nlfb
            '      Bound contstrained nonlinear fitting using function value only
            '
            _TestResult = true
            For _spoil_scenario = -1 To 22
                Try
                    ' 
                    '  In this example we demonstrate exponential fitting by
                    '      f(x) = exp(-c*x^2)
                    '  subject to bound constraints
                    '      0.0 <= c <= 1.0
                    '  using function value only.
                    ' 
                    '  Gradient is estimated using combination of numerical differences
                    '  and secant updates. diffstep variable stores differentiation step 
                    '  (we have to tell algorithm what step to use).
                    ' 
                    '  Unconstrained solution is c=1.5, but because of constraints we should
                    '  get c=1.0 (at the boundary).
                    ' 
                    Dim x(,) As Double = New Double(,){{-1},{-0.8},{-0.6},{-0.4},{-0.2},{0},{0.2},{0.4},{0.6},{0.8},{1.0}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If
                    Dim y() As Double = New Double(){0.223130,0.382893,0.582748,0.786628,0.941765,1.000000,0.941765,0.786628,0.582748,0.382893,0.223130}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim c() As Double = New Double(){0.3}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    Dim bndl() As Double = New Double(){0.0}
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){1.0}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim epsx As Double = 0.000001
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=19 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim maxits As Integer = 0
                    Dim info As Integer
                    Dim state As lsfitstate = New XAlglib.lsfitstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim diffstep As Double = 0.0001
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If

                    xalglib.lsfitcreatef(x, y, c, diffstep, state)
                    xalglib.lsfitsetbc(state, bndl, bndu)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitfit(state, AddressOf function_cx_1_func, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.0}, 0.05)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlfb")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_nlscale
            '      Nonlinear fitting with custom scaling and bound constraints
            '
            _TestResult = true
            For _spoil_scenario = -1 To 26
                Try
                    ' 
                    '  In this example we demonstrate fitting by
                    '      f(x) = c[0]*(1+c[1]*((x-1999)^c[2]-1))
                    '  subject to bound constraints
                    '      -INF  < c[0] < +INF
                    '       -10 <= c[1] <= +10
                    '       0.1 <= c[2] <= 2.0
                    '  Data we want to fit are time series of Japan national debt
                    '  collected from 2000 to 2008 measured in USD (dollars, not
                    '  millions of dollars).
                    ' 
                    '  Our variables are:
                    '      c[0] - debt value at initial moment (2000),
                    '      c[1] - direction coefficient (growth or decrease),
                    '      c[2] - curvature coefficient.
                    '  You may see that our variables are badly scaled - first one 
                    '  is order of 10^12, and next two are somewhere about 1 in 
                    '  magnitude. Such problem is difficult to solve without some
                    '  kind of scaling.
                    '  That is exactly where lsfitsetscale() function can be used.
                    '  We set scale of our variables to [1.0E12, 1, 1], which allows
                    '  us to easily solve this problem.
                    ' 
                    '  You can try commenting out lsfitsetscale() call - and you will 
                    '  see that algorithm will fail to converge.
                    ' 
                    Dim x(,) As Double = New Double(,){{2000},{2001},{2002},{2003},{2004},{2005},{2006},{2007},{2008}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(x)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(x)
                    End If
                    Dim y() As Double = New Double(){4323239600000.0,4560913100000.0,5564091500000.0,6743189300000.0,7284064600000.0,7050129600000.0,7092221500000.0,8483907600000.0,8625804400000.0}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim c() As Double = New Double(){1.0e+13,1,1}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(c, v_spoil_real)
                    End If
                    Dim epsx As Double = 1.0e-5
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        epsx = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        epsx = v_spoil_real
                    End If
                    Dim bndl() As Double = New Double(){-Double.PositiveInfinity,-10,0.1}
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndl, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        spoil_vector_by_deleting_element(bndl)
                    End If
                    Dim bndu() As Double = New Double(){Double.PositiveInfinity,+10,2.0}
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(bndu, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        spoil_vector_by_deleting_element(bndu)
                    End If
                    Dim s() As Double = New Double(){1.0e+12,1,1}
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(s, v_spoil_real)
                    End If
                    If _spoil_scenario=23 Then
                        spoil_vector_by_deleting_element(s)
                    End If
                    Dim maxits As Integer = 0
                    Dim info As Integer
                    Dim state As lsfitstate = New XAlglib.lsfitstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim diffstep As Double = 1.0e-5
                    If _spoil_scenario=24 Then
                        v_spoil_real = Double.NaN
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=25 Then
                        v_spoil_real = Double.PositiveInfinity
                        diffstep = v_spoil_real
                    End If
                    If _spoil_scenario=26 Then
                        v_spoil_real = Double.NegativeInfinity
                        diffstep = v_spoil_real
                    End If

                    xalglib.lsfitcreatef(x, y, c, diffstep, state)
                    xalglib.lsfitsetcond(state, epsx, maxits)
                    xalglib.lsfitsetbc(state, bndl, bndu)
                    xalglib.lsfitsetscale(state, s)
                    xalglib.lsfitfit(state, AddressOf function_debt_func, Nothing, Nothing)
                    xalglib.lsfitresults(state, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 2)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){4.142560e+12,0.434240,0.565376}, -0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_nlscale")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_lin
            '      Unconstrained (general) linear least squares fitting with and without weights
            '
            _TestResult = true
            For _spoil_scenario = -1 To 12
                Try
                    ' 
                    '  In this example we demonstrate linear fitting by f(x|a) = a*exp(0.5*x).
                    ' 
                    '  We have:
                    '  * y - vector of experimental data
                    '  * fmatrix -  matrix of basis functions calculated at sample points
                    '               Actually, we have only one basis function F0 = exp(0.5*x).
                    ' 
                    Dim fmatrix(,) As Double = New Double(,){{0.606531},{0.670320},{0.740818},{0.818731},{0.904837},{1.000000},{1.105171},{1.221403},{1.349859},{1.491825},{1.648721}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    Dim y() As Double = New Double(){1.133719,1.306522,1.504604,1.554663,1.884638,2.072436,2.257285,2.534068,2.622017,2.897713,3.219371}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim info As Integer
                    Dim c() As Double = New Double(){}
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Linear fitting without weights
                    ' 
                    xalglib.lsfitlinear(y, fmatrix, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.98650}, 0.00005)

                    ' 
                    '  Linear fitting with individual weights.
                    '  Slightly different result is returned.
                    ' 
                    Dim w() As Double = New Double(){1.414213,1,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    xalglib.lsfitlinearw(y, w, fmatrix, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){1.983354}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_lin")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_linc
            '      Constrained (general) linear least squares fitting with and without weights
            '
            _TestResult = true
            For _spoil_scenario = -1 To 19
                Try
                    ' 
                    '  In this example we demonstrate linear fitting by f(x|a,b) = a*x+b
                    '  with simple constraint f(0)=0.
                    ' 
                    '  We have:
                    '  * y - vector of experimental data
                    '  * fmatrix -  matrix of basis functions sampled at [0,1] with step 0.2:
                    '                   [ 1.0   0.0 ]
                    '                   [ 1.0   0.2 ]
                    '                   [ 1.0   0.4 ]
                    '                   [ 1.0   0.6 ]
                    '                   [ 1.0   0.8 ]
                    '                   [ 1.0   1.0 ]
                    '               first column contains value of first basis function (constant term)
                    '               second column contains second basis function (linear term)
                    '  * cmatrix -  matrix of linear constraints:
                    '                   [ 1.0  0.0  0.0 ]
                    '               first two columns contain coefficients before basis functions,
                    '               last column contains desired value of their sum.
                    '               So [1,0,0] means "1*constant_term + 0*linear_term = 0" 
                    ' 
                    Dim y() As Double = New Double(){0.072436,0.246944,0.491263,0.522300,0.714064,0.921929}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim fmatrix(,) As Double = New Double(,){{1,0.0},{1,0.2},{1,0.4},{1,0.6},{1,0.8},{1,1.0}}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_row(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_col(fmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        spoil_matrix_by_deleting_row(fmatrix)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_matrix_by_deleting_col(fmatrix)
                    End If
                    Dim cmatrix(,) As Double = New Double(,){{1,0,0}}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(cmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(cmatrix, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(cmatrix, v_spoil_real)
                    End If
                    Dim info As Integer
                    Dim c() As Double = New Double(){}
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Constrained fitting without weights
                    ' 
                    xalglib.lsfitlinearc(y, fmatrix, cmatrix, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){0,0.932933}, 0.0005)

                    ' 
                    '  Constrained fitting with individual weights
                    ' 
                    Dim w() As Double = New Double(){1,1.414213,1,1,1,1}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    xalglib.lsfitlinearwc(y, w, fmatrix, cmatrix, info, c, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And doc_test_real_vector(c, New Double(){0,0.938322}, 0.0005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_linc")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_pol
            '      Unconstrained polynomial fitting
            '
            _TestResult = true
            For _spoil_scenario = -1 To 19
                Try
                    ' 
                    '  This example demonstrates polynomial fitting.
                    ' 
                    '  Fitting is done by two (M=2) functions from polynomial basis:
                    '      f0 = 1
                    '      f1 = x
                    '  Basically, it just a linear fit; more complex polynomials may be used
                    '  (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
                    '  us to demonstrate polynomialfit() function in action.
                    ' 
                    '  We have:
                    '  * x      set of abscissas
                    '  * y      experimental data
                    ' 
                    '  Additionally we demonstrate weighted fitting, where second point has
                    '  more weight than other ones.
                    ' 
                    Dim x() As Double = New Double(){0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim m As Integer = 2
                    Dim t As Double = 2
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim info As Integer
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As polynomialfitreport = New XAlglib.polynomialfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double

                    ' 
                    '  Fitting without individual weights
                    ' 
                    '  NOTE: result is returned as barycentricinterpolant structure.
                    '        if you want to get representation in the power basis,
                    '        you can use barycentricbar2pow() function to convert
                    '        from barycentric to power representation (see docs for 
                    '        POLINT subpackage for more info).
                    ' 
                    xalglib.polynomialfit(x, y, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.011, 0.002)

                    ' 
                    '  Fitting with individual weights
                    ' 
                    '  NOTE: slightly different result is returned
                    ' 
                    Dim w() As Double = New Double(){1,1.414213562,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    Dim xc() As Double = New Double(){}
                    If _spoil_scenario=17 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(xc, v_spoil_real)
                    End If
                    Dim yc() As Double = New Double(){}
                    If _spoil_scenario=18 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(yc, v_spoil_real)
                    End If
                    Dim dc() As Integer = New Integer(){}
                    If _spoil_scenario=19 Then
                        v_spoil_int = 0
                        spoil_vector_by_adding_element(dc, v_spoil_int)
                    End If
                    xalglib.polynomialfitwc(x, y, w, xc, yc, dc, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.023, 0.002)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_pol")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_polc
            '      Constrained polynomial fitting
            '
            _TestResult = true
            For _spoil_scenario = -1 To 28
                Try
                    ' 
                    '  This example demonstrates polynomial fitting.
                    ' 
                    '  Fitting is done by two (M=2) functions from polynomial basis:
                    '      f0 = 1
                    '      f1 = x
                    '  with simple constraint on function value
                    '      f(0) = 0
                    '  Basically, it just a linear fit; more complex polynomials may be used
                    '  (e.g. parabolas with M=3, cubic with M=4), but even such simple fit allows
                    '  us to demonstrate polynomialfit() function in action.
                    ' 
                    '  We have:
                    '  * x      set of abscissas
                    '  * y      experimental data
                    '  * xc     points where constraints are placed
                    '  * yc     constraints on derivatives
                    '  * dc     derivative indices
                    '           (0 means function itself, 1 means first derivative)
                    ' 
                    Dim x() As Double = New Double(){1.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.9,1.1}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim w() As Double = New Double(){1,1}
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(w, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    Dim xc() As Double = New Double(){0}
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        spoil_vector_by_deleting_element(xc)
                    End If
                    Dim yc() As Double = New Double(){0}
                    If _spoil_scenario=20 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=23 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=24 Then
                        spoil_vector_by_deleting_element(yc)
                    End If
                    Dim dc() As Integer = New Integer(){0}
                    If _spoil_scenario=25 Then
                        v_spoil_int = 0
                        spoil_vector_by_adding_element(dc, v_spoil_int)
                    End If
                    If _spoil_scenario=26 Then
                        spoil_vector_by_deleting_element(dc)
                    End If
                    Dim t As Double = 2
                    If _spoil_scenario=27 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=28 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim m As Integer = 2
                    Dim info As Integer
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As polynomialfitreport = New XAlglib.polynomialfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double

                    xalglib.polynomialfitwc(x, y, w, xc, yc, dc, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.000, 0.001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_polc")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_d_spline
            '      Unconstrained fitting by penalized regression spline
            '
            _TestResult = true
            For _spoil_scenario = -1 To 18
                Try
                    ' 
                    '  In this example we demonstrate penalized spline fitting of noisy data
                    ' 
                    '  We have:
                    '  * x - abscissas
                    '  * y - vector of experimental data, straight line with small noise
                    ' 
                    Dim x() As Double = New Double(){0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(x, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.10,0.00,0.30,0.40,0.30,0.40,0.62,0.68,0.75,0.95}
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=8 Then
                        v_spoil_real = 0
                        spoil_vector_by_adding_element(y, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim info As Integer
                    Dim v As Double
                    Dim s As spline1dinterpolant = New XAlglib.spline1dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As spline1dfitreport = New XAlglib.spline1dfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim rho As Double

                    ' 
                    '  Fit with VERY small amount of smoothing (rho = -5.0)
                    '  and large number of basis functions (M=50).
                    ' 
                    '  With such small regularization penalized spline almost fully reproduces function values
                    ' 
                    rho = -5.0
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=11 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    xalglib.spline1dfitpenalized(x, y, 50, rho, info, s, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    v = xalglib.spline1dcalc(s, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.10, 0.01)

                    ' 
                    '  Fit with VERY large amount of smoothing (rho = 10.0)
                    '  and large number of basis functions (M=50).
                    ' 
                    '  With such regularization our spline should become close to the straight line fit.
                    '  We will compare its value in x=1.0 with results obtained from such fit.
                    ' 
                    rho = +10.0
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    xalglib.spline1dfitpenalized(x, y, 50, rho, info, s, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    v = xalglib.spline1dcalc(s, 1.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.969, 0.001)

                    ' 
                    '  In real life applications you may need some moderate degree of fitting,
                    '  so we try to fit once more with rho=3.0.
                    ' 
                    rho = +3.0
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        rho = v_spoil_real
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        rho = v_spoil_real
                    End If
                    xalglib.spline1dfitpenalized(x, y, 50, rho, info, s, rep)
                    _TestResult = _TestResult And doc_test_int(info, 1)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_d_spline")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_t_polfit_1
            '      Polynomial fitting, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 9
                Try
                    Dim x() As Double = New Double(){0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim m As Integer = 2
                    Dim t As Double = 2
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim info As Integer
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As polynomialfitreport = New XAlglib.polynomialfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialfit(x, y, 11, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.011, 0.002)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_t_polfit_2
            '      Polynomial fitting, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 13
                Try
                    Dim x() As Double = New Double(){0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.00,0.05,0.26,0.32,0.33,0.43,0.60,0.60,0.77,0.98,1.02}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim w() As Double = New Double(){1,1.414213562,1,1,1,1,1,1,1,1,1}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    Dim xc() As Double = New Double(){}
                    Dim yc() As Double = New Double(){}
                    Dim dc() As Integer = New Integer(){}
                    Dim m As Integer = 2
                    Dim t As Double = 2
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim info As Integer
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As polynomialfitreport = New XAlglib.polynomialfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialfitwc(x, y, w, 11, xc, yc, dc, 0, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.023, 0.002)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_t_polfit_3
            '      Polynomial fitting, full list of parameters.
            '
            _TestResult = true
            For _spoil_scenario = -1 To 22
                Try
                    Dim x() As Double = New Double(){1.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.9,1.1}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim w() As Double = New Double(){1,1}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(w, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(w)
                    End If
                    Dim xc() As Double = New Double(){0}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(xc, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        spoil_vector_by_deleting_element(xc)
                    End If
                    Dim yc() As Double = New Double(){0}
                    If _spoil_scenario=16 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=17 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=18 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(yc, v_spoil_real)
                    End If
                    If _spoil_scenario=19 Then
                        spoil_vector_by_deleting_element(yc)
                    End If
                    Dim dc() As Integer = New Integer(){0}
                    If _spoil_scenario=20 Then
                        spoil_vector_by_deleting_element(dc)
                    End If
                    Dim m As Integer = 2
                    Dim t As Double = 2
                    If _spoil_scenario=21 Then
                        v_spoil_real = Double.PositiveInfinity
                        t = v_spoil_real
                    End If
                    If _spoil_scenario=22 Then
                        v_spoil_real = Double.NegativeInfinity
                        t = v_spoil_real
                    End If
                    Dim info As Integer
                    Dim p As barycentricinterpolant = New XAlglib.barycentricinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As polynomialfitreport = New XAlglib.polynomialfitreport() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.polynomialfitwc(x, y, w, 2, xc, yc, dc, 1, m, info, p, rep)
                    v = xalglib.barycentriccalc(p, t)
                    _TestResult = _TestResult And doc_test_real(v, 2.000, 0.001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_polfit_3")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_t_4pl
            '      4-parameter logistic fitting
            '
            _TestResult = true
            For _spoil_scenario = -1 To 7
                Try
                    Dim x() As Double = New Double(){1,2,3,4,5,6,7,8}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.06313223,0.44552624,0.61838364,0.71385108,0.77345838,0.81383140,0.84280033,0.86449822}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim n As Integer = 8
                    Dim a As Double
                    Dim b As Double
                    Dim c As Double
                    Dim d As Double
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Test logisticfit4() on carefully designed data with a priori known answer.
                    ' 
                    xalglib.logisticfit4(x, y, n, a, b, c, d, rep)
                    _TestResult = _TestResult And doc_test_real(a, -1.000, 0.01)
                    _TestResult = _TestResult And doc_test_real(b, 1.200, 0.01)
                    _TestResult = _TestResult And doc_test_real(c, 0.900, 0.01)
                    _TestResult = _TestResult And doc_test_real(d, 1.000, 0.01)

                    ' 
                    '  Evaluate model at point x=0.5
                    ' 
                    Dim v As Double
                    v = xalglib.logisticcalc4(0.5, a, b, c, d)
                    _TestResult = _TestResult And doc_test_real(v, -0.33874308, 0.001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_4pl")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lsfit_t_5pl
            '      5-parameter logistic fitting
            '
            _TestResult = true
            For _spoil_scenario = -1 To 7
                Try
                    Dim x() As Double = New Double(){1,2,3,4,5,6,7,8}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.1949776139,0.5710060208,0.726002637,0.8060434158,0.8534547965,0.8842071579,0.9054773317,0.9209088299}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim n As Integer = 8
                    Dim a As Double
                    Dim b As Double
                    Dim c As Double
                    Dim d As Double
                    Dim g As Double
                    Dim rep As lsfitreport = New XAlglib.lsfitreport() ' initializer can be dropped, but compiler will issue warning

                    ' 
                    '  Test logisticfit5() on carefully designed data with a priori known answer.
                    ' 
                    xalglib.logisticfit5(x, y, n, a, b, c, d, g, rep)
                    _TestResult = _TestResult And doc_test_real(a, -1.000, 0.01)
                    _TestResult = _TestResult And doc_test_real(b, 1.200, 0.01)
                    _TestResult = _TestResult And doc_test_real(c, 0.900, 0.01)
                    _TestResult = _TestResult And doc_test_real(d, 1.000, 0.01)
                    _TestResult = _TestResult And doc_test_real(g, 1.200, 0.01)

                    ' 
                    '  Evaluate model at point x=0.5
                    ' 
                    Dim v As Double
                    v = xalglib.logisticcalc5(0.5, a, b, c, d, g)
                    _TestResult = _TestResult And doc_test_real(v, -0.2354656824, 0.001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lsfit_t_5pl")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_bilinear
            '      Bilinear spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 15
                Try
                    ' 
                    '  We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled 
                    '  at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
                    ' 
                    Dim x() As Double = New Double(){0.0,0.5,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim f() As Double = New Double(){0.00,0.25,1.00,2.00,2.25,3.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim vx As Double = 0.25
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        vx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        vx = v_spoil_real
                    End If
                    Dim vy As Double = 0.50
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        vy = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        vy = v_spoil_real
                    End If
                    Dim v As Double
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline2dbuildbilinearv(x, 3, y, 2, f, 1, s)

                    '  calculate S(0.25,0.50)
                    v = xalglib.spline2dcalc(s, vx, vy)
                    _TestResult = _TestResult And doc_test_real(v, 1.1250, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_bilinear")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_bicubic
            '      Bilinear spline interpolation
            '
            _TestResult = true
            For _spoil_scenario = -1 To 15
                Try
                    ' 
                    '  We use bilinear spline to interpolate f(x,y)=x^2+2*y^2 sampled 
                    '  at (x,y) from [0.0, 0.5, 1.0] X [0.0, 1.0].
                    ' 
                    Dim x() As Double = New Double(){0.0,0.5,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim f() As Double = New Double(){0.00,0.25,1.00,2.00,2.25,3.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim vx As Double = 0.25
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.PositiveInfinity
                        vx = v_spoil_real
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.NegativeInfinity
                        vx = v_spoil_real
                    End If
                    Dim vy As Double = 0.50
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.PositiveInfinity
                        vy = v_spoil_real
                    End If
                    If _spoil_scenario=15 Then
                        v_spoil_real = Double.NegativeInfinity
                        vy = v_spoil_real
                    End If
                    Dim v As Double
                    Dim dx As Double
                    Dim dy As Double
                    Dim dxy As Double
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline2dbuildbicubicv(x, 3, y, 2, f, 1, s)

                    '  calculate S(0.25,0.50)
                    v = xalglib.spline2dcalc(s, vx, vy)
                    _TestResult = _TestResult And doc_test_real(v, 1.0625, 0.00005)

                    '  calculate derivatives
                    xalglib.spline2ddiff(s, vx, vy, v, dx, dy, dxy)
                    _TestResult = _TestResult And doc_test_real(v, 1.0625, 0.00005)
                    _TestResult = _TestResult And doc_test_real(dx, 0.5000, 0.00005)
                    _TestResult = _TestResult And doc_test_real(dy, 2.0000, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_bicubic")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_fit_blocklls
            '      Fitting bicubic spline to irregular data
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    ' 
                    '  We use bicubic spline to reproduce f(x,y)=1/(1+x^2+2*y^2) sampled
                    '  at irregular points (x,y) from [-1,+1]*[-1,+1]
                    ' 
                    '  We have 5 such points, located approximately at corners of the area
                    '  and its center -  but not exactly at the grid. Thus, we have to FIT
                    '  the spline, i.e. to solve least squares problem
                    ' 
                    Dim xy(,) As Double = New Double(,){{-0.987,-0.902,0.359},{0.948,-0.992,0.347},{-1.000,1.000,0.333},{1.000,0.973,0.339},{0.017,0.180,0.968}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(xy)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(xy)
                    End If

                    ' 
                    '  First step is to create spline2dbuilder object and set its properties:
                    '  * d=1 means that we create vector-valued spline with 1 component
                    '  * we specify dataset xy
                    '  * we rely on automatic selection of interpolation area
                    '  * we tell builder that we want to use 5x5 grid for an underlying spline
                    '  * we choose least squares solver named BlockLLS and configure it by
                    '    telling that we want to apply zero nonlinearity penalty.
                    ' 
                    '  NOTE: you can specify non-zero lambdav if you want to make your spline
                    '        more "rigid", i.e. to penalize nonlinearity.
                    ' 
                    '  NOTE: ALGLIB has two solvers which fit bicubic splines to irregular data,
                    '        one of them is BlockLLS and another one is FastDDM. Former is
                    '        intended for moderately sized grids (up to 512x512 nodes, although
                    '        it may take up to few minutes); it is the most easy to use and
                    '        control spline fitting function in the library. Latter, FastDDM,
                    '        is intended for efficient solution of large-scale problems
                    '        (up to 100.000.000 nodes). Both solvers can be parallelized, but
                    '        FastDDM is much more efficient. See comments for more information.
                    ' 
                    Dim builder As spline2dbuilder = New XAlglib.spline2dbuilder() ' initializer can be dropped, but compiler will issue warning
                    Dim d As Integer = 1
                    Dim lambdav As Double = 0.000
                    xalglib.spline2dbuildercreate(d, builder)
                    xalglib.spline2dbuildersetpoints(builder, xy, 5)
                    xalglib.spline2dbuildersetgrid(builder, 5, 5)
                    xalglib.spline2dbuildersetalgoblocklls(builder, lambdav)

                    ' 
                    '  Now we are ready to fit and evaluate our results
                    ' 
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As spline2dfitreport = New XAlglib.spline2dfitreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.spline2dfit(builder, s, rep)

                    '  evaluate results - function value at the grid is reproduced exactly
                    Dim v As Double
                    v = xalglib.spline2dcalc(s, -1, 1)
                    _TestResult = _TestResult And doc_test_real(v, 0.333000, 0.005)

                    '  check maximum error - it must be nearly zero
                    _TestResult = _TestResult And doc_test_real(rep.maxerror, 0.000, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_fit_blocklls")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_unpack
            '      Unpacking bilinear spline
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    ' 
                    '  We build bilinear spline for f(x,y)=x+2*y+3*xy for (x,y) in [0,1].
                    '  Then we demonstrate how to unpack it.
                    ' 
                    Dim x() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim f() As Double = New Double(){0.00,1.00,2.00,6.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim c(,) As Double = New Double(,){{}}
                    Dim m As Integer
                    Dim n As Integer
                    Dim d As Integer
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning

                    '  build spline
                    xalglib.spline2dbuildbilinearv(x, 2, y, 2, f, 1, s)

                    '  unpack and test
                    xalglib.spline2dunpackv(s, m, n, d, c)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0,1,0,1,0,2,0,0,1,3,0,0,0,0,0,0,0,0,0,0}}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_unpack")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_copytrans
            '      Copy and transform
            '
            _TestResult = true
            For _spoil_scenario = -1 To 15
                Try
                    ' 
                    '  We build bilinear spline for f(x,y)=x+2*y for (x,y) in [0,1].
                    '  Then we apply several transformations to this spline.
                    ' 
                    Dim x() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim f() As Double = New Double(){0.00,1.00,2.00,3.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim snew As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim v As Double
                    xalglib.spline2dbuildbilinearv(x, 2, y, 2, f, 1, s)

                    '  copy spline, apply transformation x:=2*xnew, y:=4*ynew
                    '  evaluate at (xnew,ynew) = (0.25,0.25) - should be same as (x,y)=(0.5,1.0)
                    xalglib.spline2dcopy(s, snew)
                    xalglib.spline2dlintransxy(snew, 2.0, 0.0, 4.0, 0.0)
                    v = xalglib.spline2dcalc(snew, 0.25, 0.25)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.00005)

                    '  copy spline, apply transformation SNew:=2*S+3
                    xalglib.spline2dcopy(s, snew)
                    xalglib.spline2dlintransf(snew, 2.0, 3.0)
                    v = xalglib.spline2dcalc(snew, 0.5, 1.0)
                    _TestResult = _TestResult And doc_test_real(v, 8.000, 0.00005)

                    ' 
                    '  Same example, but for vector spline (f0,f1) = {x+2*y, 2*x+y}
                    ' 
                    Dim f2() As Double = New Double(){0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00}
                    If _spoil_scenario=12 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f2, v_spoil_real)
                    End If
                    If _spoil_scenario=13 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f2, v_spoil_real)
                    End If
                    If _spoil_scenario=14 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f2, v_spoil_real)
                    End If
                    If _spoil_scenario=15 Then
                        spoil_vector_by_deleting_element(f2)
                    End If
                    Dim vr() As Double = New Double(){}
                    xalglib.spline2dbuildbilinearv(x, 2, y, 2, f2, 2, s)

                    '  copy spline, apply transformation x:=2*xnew, y:=4*ynew
                    xalglib.spline2dcopy(s, snew)
                    xalglib.spline2dlintransxy(snew, 2.0, 0.0, 4.0, 0.0)
                    xalglib.spline2dcalcv(snew, 0.25, 0.25, vr)
                    _TestResult = _TestResult And doc_test_real_vector(vr, New Double(){2.500,2.000}, 0.00005)

                    '  copy spline, apply transformation SNew:=2*S+3
                    xalglib.spline2dcopy(s, snew)
                    xalglib.spline2dlintransf(snew, 2.0, 3.0)
                    xalglib.spline2dcalcv(snew, 0.5, 1.0, vr)
                    _TestResult = _TestResult And doc_test_real_vector(vr, New Double(){8.000,7.000}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_copytrans")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST spline2d_vector
            '      Copy and transform
            '
            _TestResult = true
            For _spoil_scenario = -1 To 11
                Try
                    ' 
                    '  We build bilinear vector-valued spline (f0,f1) = {x+2*y, 2*x+y}
                    '  Spline is built using function values at 2x2 grid: (x,y)=[0,1]*[0,1]
                    '  Then we perform evaluation at (x,y)=(0.1,0.3)
                    ' 
                    Dim x() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(x)
                    End If
                    Dim y() As Double = New Double(){0.0,1.0}
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=6 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(y, v_spoil_real)
                    End If
                    If _spoil_scenario=7 Then
                        spoil_vector_by_deleting_element(y)
                    End If
                    Dim f() As Double = New Double(){0.00,0.00,1.00,2.00,2.00,1.00,3.00,3.00}
                    If _spoil_scenario=8 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=9 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=10 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(f, v_spoil_real)
                    End If
                    If _spoil_scenario=11 Then
                        spoil_vector_by_deleting_element(f)
                    End If
                    Dim s As spline2dinterpolant = New XAlglib.spline2dinterpolant() ' initializer can be dropped, but compiler will issue warning
                    Dim vr() As Double = New Double(){}
                    xalglib.spline2dbuildbilinearv(x, 2, y, 2, f, 2, s)
                    xalglib.spline2dcalcv(s, 0.1, 0.3, vr)
                    _TestResult = _TestResult And doc_test_real_vector(vr, New Double(){0.700,0.500}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "spline2d_vector")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST rbf_d_hrbf
            '      Simple model built with HRBF algorithm
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example illustrates basic concepts of the RBF models: creation, modification,
                    '  evaluation.
                    '  
                    '  Suppose that we have set of 2-dimensional points with associated
                    '  scalar function values, and we want to build a RBF model using
                    '  our data.
                    '  
                    '  NOTE: we can work with 3D models too :)
                    '  
                    '  Typical sequence of steps is given below:
                    '  1. we create RBF model object
                    '  2. we attach our dataset to the RBF model and tune algorithm settings
                    '  3. we rebuild RBF model using QNN algorithm on new data
                    '  4. we use RBF model (evaluate, serialize, etc.)
                    ' 
                    Dim v As Double

                    ' 
                    '  Step 1: RBF model creation.
                    ' 
                    '  We have to specify dimensionality of the space (2 or 3) and
                    '  dimensionality of the function (scalar or vector).
                    ' 
                    '  New model is empty - it can be evaluated,
                    '  but we just get zero value at any point.
                    ' 
                    Dim model As rbfmodel = New XAlglib.rbfmodel() ' initializer can be dropped, but compiler will issue warning
                    xalglib.rbfcreate(2, 1, model)

                    v = xalglib.rbfcalc2(model, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.000, 0.005)

                    ' 
                    '  Step 2: we add dataset.
                    ' 
                    '  XY contains two points - x0=(-1,0) and x1=(+1,0) -
                    '  and two function values f(x0)=2, f(x1)=3.
                    ' 
                    '  We added points, but model was not rebuild yet.
                    '  If we call rbfcalc2(), we still will get 0.0 as result.
                    ' 
                    Dim xy(,) As Double = New Double(,){{-1,0,2},{+1,0,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.rbfsetpoints(model, xy)

                    v = xalglib.rbfcalc2(model, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.000, 0.005)

                    ' 
                    '  Step 3: rebuild model
                    ' 
                    '  After we've configured model, we should rebuild it -
                    '  it will change coefficients stored internally in the
                    '  rbfmodel structure.
                    ' 
                    '  We use hierarchical RBF algorithm with following parameters:
                    '  * RBase - set to 1.0
                    '  * NLayers - three layers are used (although such simple problem
                    '    does not need more than 1 layer)
                    '  * LambdaReg - is set to zero value, no smoothing is required
                    ' 
                    Dim rep As rbfreport = New XAlglib.rbfreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0)
                    xalglib.rbfbuildmodel(model, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)

                    ' 
                    '  Step 4: model was built
                    ' 
                    '  After call of rbfbuildmodel(), rbfcalc2() will return
                    '  value of the new model.
                    ' 
                    v = xalglib.rbfcalc2(model, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_hrbf")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST rbf_d_vector
            '      Working with vector functions
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    ' 
                    '  Suppose that we have set of 2-dimensional points with associated VECTOR
                    '  function values, and we want to build a RBF model using our data.
                    '  
                    '  Typical sequence of steps is given below:
                    '  1. we create RBF model object
                    '  2. we attach our dataset to the RBF model and tune algorithm settings
                    '  3. we rebuild RBF model using new data
                    '  4. we use RBF model (evaluate, serialize, etc.)
                    ' 
                    Dim x() As Double = New Double(){}
                    Dim y() As Double = New Double(){}

                    ' 
                    '  Step 1: RBF model creation.
                    ' 
                    '  We have to specify dimensionality of the space (equal to 2) and
                    '  dimensionality of the function (2-dimensional vector function).
                    ' 
                    '  New model is empty - it can be evaluated,
                    '  but we just get zero value at any point.
                    ' 
                    Dim model As rbfmodel = New XAlglib.rbfmodel() ' initializer can be dropped, but compiler will issue warning
                    xalglib.rbfcreate(2, 2, model)

                    x = New Double(){+1,+1}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(x, v_spoil_real)
                    End If
                    xalglib.rbfcalc(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,0.000}, 0.005)

                    ' 
                    '  Step 2: we add dataset.
                    ' 
                    '  XY arrays containt four points:
                    '  * (x0,y0) = (+1,+1), f(x0,y0)=(0,-1)
                    '  * (x1,y1) = (+1,-1), f(x1,y1)=(-1,0)
                    '  * (x2,y2) = (-1,-1), f(x2,y2)=(0,+1)
                    '  * (x3,y3) = (-1,+1), f(x3,y3)=(+1,0)
                    ' 
                    Dim xy(,) As Double = New Double(,){{+1,+1,0,-1},{+1,-1,-1,0},{-1,-1,0,+1},{-1,+1,+1,0}}
                    If _spoil_scenario=3 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    xalglib.rbfsetpoints(model, xy)

                    '  We added points, but model was not rebuild yet.
                    '  If we call rbfcalc(), we still will get 0.0 as result.
                    xalglib.rbfcalc(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,0.000}, 0.005)

                    ' 
                    '  Step 3: rebuild model
                    ' 
                    '  We use hierarchical RBF algorithm with following parameters:
                    '  * RBase - set to 1.0
                    '  * NLayers - three layers are used (although such simple problem
                    '    does not need more than 1 layer)
                    '  * LambdaReg - is set to zero value, no smoothing is required
                    ' 
                    '  After we've configured model, we should rebuild it -
                    '  it will change coefficients stored internally in the
                    '  rbfmodel structure.
                    ' 
                    Dim rep As rbfreport = New XAlglib.rbfreport() ' initializer can be dropped, but compiler will issue warning
                    xalglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0)
                    xalglib.rbfbuildmodel(model, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)

                    ' 
                    '  Step 4: model was built
                    ' 
                    '  After call of rbfbuildmodel(), rbfcalc() will return
                    '  value of the new model.
                    ' 
                    xalglib.rbfcalc(model, x, y)
                    _TestResult = _TestResult And doc_test_real_vector(y, New Double(){0.000,-1.000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_vector")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST rbf_d_polterm
            '      RBF models - working with polynomial term
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example show how to work with polynomial term
                    '  
                    '  Suppose that we have set of 2-dimensional points with associated
                    '  scalar function values, and we want to build a RBF model using
                    '  our data.
                    ' 
                    '  We use hierarchical RBF algorithm with following parameters:
                    '  * RBase - set to 1.0
                    '  * NLayers - three layers are used (although such simple problem
                    '    does not need more than 1 layer)
                    '  * LambdaReg - is set to zero value, no smoothing is required
                    ' 
                    Dim v As Double
                    Dim model As rbfmodel = New XAlglib.rbfmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{-1,0,2},{+1,0,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    Dim rep As rbfreport = New XAlglib.rbfreport() ' initializer can be dropped, but compiler will issue warning

                    xalglib.rbfcreate(2, 1, model)
                    xalglib.rbfsetpoints(model, xy)
                    xalglib.rbfsetalgohierarchical(model, 1.0, 3, 0.0)

                    ' 
                    '  By default, RBF model uses linear term. It means that model
                    '  looks like
                    '      f(x,y) = SUM(RBF[i]) + a*x + b*y + c
                    '  where RBF[i] is I-th radial basis function and a*x+by+c is a
                    '  linear term. Having linear terms in a model gives us:
                    '  (1) improved extrapolation properties
                    '  (2) linearity of the model when data can be perfectly fitted
                    '      by the linear function
                    '  (3) linear asymptotic behavior
                    ' 
                    '  Our simple dataset can be modelled by the linear function
                    '      f(x,y) = 0.5*x + 2.5
                    '  and rbfbuildmodel() with default settings should preserve this
                    '  linearity.
                    ' 
                    Dim nx As Integer
                    Dim ny As Integer
                    Dim nc As Integer
                    Dim modelversion As Integer
                    Dim xwr(,) As Double = New Double(,){{}}
                    Dim c(,) As Double = New Double(,){{}}
                    xalglib.rbfbuildmodel(model, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)
                    xalglib.rbfunpack(model, nx, ny, xwr, nc, c, modelversion)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0.500,0.000,2.500}}, 0.005)

                    '  asymptotic behavior of our function is linear
                    v = xalglib.rbfcalc2(model, 1000.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 502.50, 0.05)

                    ' 
                    '  Instead of linear term we can use constant term. In this case
                    '  we will get model which has form
                    '      f(x,y) = SUM(RBF[i]) + c
                    '  where RBF[i] is I-th radial basis function and c is a constant,
                    '  which is equal to the average function value on the dataset.
                    ' 
                    '  Because we've already attached dataset to the model the only
                    '  thing we have to do is to call rbfsetconstterm() and then
                    '  rebuild model with rbfbuildmodel().
                    ' 
                    xalglib.rbfsetconstterm(model)
                    xalglib.rbfbuildmodel(model, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)
                    xalglib.rbfunpack(model, nx, ny, xwr, nc, c, modelversion)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0.000,0.000,2.500}}, 0.005)

                    '  asymptotic behavior of our function is constant
                    v = xalglib.rbfcalc2(model, 1000.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.005)

                    ' 
                    '  Finally, we can use zero term. Just plain RBF without polynomial
                    '  part:
                    '      f(x,y) = SUM(RBF[i])
                    '  where RBF[i] is I-th radial basis function.
                    ' 
                    xalglib.rbfsetzeroterm(model)
                    xalglib.rbfbuildmodel(model, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)
                    xalglib.rbfunpack(model, nx, ny, xwr, nc, c, modelversion)
                    _TestResult = _TestResult And doc_test_real_matrix(c, New Double(,){{0.000,0.000,0.000}}, 0.005)

                    '  asymptotic behavior of our function is just zero constant
                    v = xalglib.rbfcalc2(model, 1000.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.000, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_polterm")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST rbf_d_serialize
            '      Serialization/unserialization
            '
            _TestResult = true
            For _spoil_scenario = -1 To 2
                Try
                    ' 
                    '  This example show how to serialize and unserialize RBF model
                    '  
                    '  Suppose that we have set of 2-dimensional points with associated
                    '  scalar function values, and we want to build a RBF model using
                    '  our data. Then we want to serialize it to string and to unserialize
                    '  from string, loading to another instance of RBF model.
                    ' 
                    '  Here we assume that you already know how to create RBF models.
                    ' 
                    Dim s As String = "" ' initializer can be dropped, but compiler will issue warning 
                    Dim v As Double
                    Dim model0 As rbfmodel = New XAlglib.rbfmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim model1 As rbfmodel = New XAlglib.rbfmodel() ' initializer can be dropped, but compiler will issue warning
                    Dim xy(,) As Double = New Double(,){{-1,0,2},{+1,0,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(xy, v_spoil_real)
                    End If
                    Dim rep As rbfreport = New XAlglib.rbfreport() ' initializer can be dropped, but compiler will issue warning

                    '  model initialization
                    xalglib.rbfcreate(2, 1, model0)
                    xalglib.rbfsetpoints(model0, xy)
                    xalglib.rbfsetalgohierarchical(model0, 1.0, 3, 0.0)
                    xalglib.rbfbuildmodel(model0, rep)
                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)

                    ' 
                    '  Serialization - it looks easy,
                    '  but you should carefully read next section.
                    ' 
                    xalglib.rbfserialize(model0, s)
                    xalglib.rbfunserialize(s, model1)

                    '  both models return same value
                    v = xalglib.rbfcalc2(model0, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.005)
                    v = xalglib.rbfcalc2(model1, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.005)

                    ' 
                    '  Previous section shows that model state is saved/restored during
                    '  serialization. However, some properties are NOT serialized.
                    ' 
                    '  Serialization saves/restores RBF model, but it does NOT saves/restores
                    '  settings which were used to build current model. In particular, dataset
                    '  which was used to build model, is not preserved.
                    ' 
                    '  What does it mean in for us?
                    ' 
                    '  Do you remember this sequence: rbfcreate-rbfsetpoints-rbfbuildmodel?
                    '  First step creates model, second step adds dataset and tunes model
                    '  settings, third step builds model using current dataset and model
                    '  construction settings.
                    ' 
                    '  If you call rbfbuildmodel() without calling rbfsetpoints() first, you
                    '  will get empty (zero) RBF model. In our example, model0 contains
                    '  dataset which was added by rbfsetpoints() call. However, model1 does
                    '  NOT contain dataset - because dataset is NOT serialized.
                    ' 
                    '  This, if we call rbfbuildmodel(model0,rep), we will get same model,
                    '  which returns 2.5 at (x,y)=(0,0). However, after same call model1 will
                    '  return zero - because it contains RBF model (coefficients), but does NOT
                    '  contain dataset which was used to build this model.
                    ' 
                    '  Basically, it means that:
                    '  * serialization of the RBF model preserves anything related to the model
                    '    EVALUATION
                    '  * but it does NOT creates perfect copy of the original object.
                    ' 
                    xalglib.rbfbuildmodel(model0, rep)
                    v = xalglib.rbfcalc2(model0, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 2.500, 0.005)

                    xalglib.rbfbuildmodel(model1, rep)
                    v = xalglib.rbfcalc2(model1, 0.0, 0.0)
                    _TestResult = _TestResult And doc_test_real(v, 0.000, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "rbf_d_serialize")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_d_1
            '      Determinant calculation, real matrix, short form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim b(,) As Double = New Double(,){{1,2},{2,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_row(b, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_col(b, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim a As Double
                    a = xalglib.rmatrixdet(b)
                    _TestResult = _TestResult And doc_test_real(a, -3, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_d_2
            '      Determinant calculation, real matrix, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim b(,) As Double = New Double(,){{5,4},{4,5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim a As Double
                    a = xalglib.rmatrixdet(b, 2)
                    _TestResult = _TestResult And doc_test_real(a, 9, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_d_3
            '      Determinant calculation, complex matrix, short form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim b(,) As alglib.complex = New alglib.complex(,){{New alglib.complex(1,+1),2},{2,New alglib.complex(1,-1)}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_row(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_col(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim a As alglib.complex
                    a = xalglib.cmatrixdet(b)
                    _TestResult = _TestResult And doc_test_complex(a, -2, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_3")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_d_4
            '      Determinant calculation, complex matrix, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim a As alglib.complex
                    Dim b(,) As alglib.complex = New alglib.complex(,){{New alglib.complex(0,5),4},{New alglib.complex(0,4),5}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    a = xalglib.cmatrixdet(b, 2)
                    _TestResult = _TestResult And doc_test_complex(a, New alglib.complex(0,9), 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_4")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_d_5
            '      Determinant calculation, complex matrix with zero imaginary part, short form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 6
                Try
                    Dim a As alglib.complex
                    Dim b(,) As alglib.complex = New alglib.complex(,){{9,1},{2,1}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_row(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_col(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    a = xalglib.cmatrixdet(b)
                    _TestResult = _TestResult And doc_test_complex(a, 7, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_d_5")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_0
            '      Determinant calculation, real matrix, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim a As Double
                    Dim b(,) As Double = New Double(,){{3,4},{-4,3}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    a = xalglib.rmatrixdet(b, 2)
                    _TestResult = _TestResult And doc_test_real(a, 25, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_0")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_1
            '      Determinant calculation, real matrix, LU, short form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    Dim a As Double
                    Dim b(,) As Double = New Double(,){{1,2},{2,5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_row(b, v_spoil_real)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_real = 0
                        spoil_matrix_by_adding_col(b, v_spoil_real)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim p() As Integer = New Integer(){1,1}
                    If _spoil_scenario=7 Then
                        v_spoil_int = 0
                        spoil_vector_by_adding_element(p, v_spoil_int)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(p)
                    End If
                    a = xalglib.rmatrixludet(b, p)
                    _TestResult = _TestResult And doc_test_real(a, -5, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_2
            '      Determinant calculation, real matrix, LU, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim a As Double
                    Dim b(,) As Double = New Double(,){{5,4},{4,5}}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim p() As Integer = New Integer(){0,1}
                    If _spoil_scenario=5 Then
                        spoil_vector_by_deleting_element(p)
                    End If
                    a = xalglib.rmatrixludet(b, p, 2)
                    _TestResult = _TestResult And doc_test_real(a, 25, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_2")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_3
            '      Determinant calculation, complex matrix, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 4
                Try
                    Dim a As alglib.complex
                    Dim b(,) As alglib.complex = New alglib.complex(,){{New alglib.complex(0,5),4},{-4,New alglib.complex(0,5)}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    a = xalglib.cmatrixdet(b, 2)
                    _TestResult = _TestResult And doc_test_complex(a, -9, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_3")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_4
            '      Determinant calculation, complex matrix, LU, short form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 8
                Try
                    Dim a As alglib.complex
                    Dim b(,) As alglib.complex = New alglib.complex(,){{1,2},{2,New alglib.complex(0,5)}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_row(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=4 Then
                        v_spoil_complex = 0
                        spoil_matrix_by_adding_col(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=5 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=6 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim p() As Integer = New Integer(){1,1}
                    If _spoil_scenario=7 Then
                        v_spoil_int = 0
                        spoil_vector_by_adding_element(p, v_spoil_int)
                    End If
                    If _spoil_scenario=8 Then
                        spoil_vector_by_deleting_element(p)
                    End If
                    a = xalglib.cmatrixludet(b, p)
                    _TestResult = _TestResult And doc_test_complex(a, New alglib.complex(0,-5), 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_4")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST matdet_t_5
            '      Determinant calculation, complex matrix, LU, full form
            '
            _TestResult = true
            For _spoil_scenario = -1 To 5
                Try
                    Dim a As alglib.complex
                    Dim b(,) As alglib.complex = New alglib.complex(,){{5,New alglib.complex(0,4)},{4,5}}
                    If _spoil_scenario=0 Then
                        v_spoil_complex = Double.NaN
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_complex = Double.PositiveInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_complex = Double.NegativeInfinity
                        spoil_matrix_by_value(b, v_spoil_complex)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_matrix_by_deleting_row(b)
                    End If
                    If _spoil_scenario=4 Then
                        spoil_matrix_by_deleting_col(b)
                    End If
                    Dim p() As Integer = New Integer(){0,1}
                    If _spoil_scenario=5 Then
                        spoil_vector_by_deleting_element(p)
                    End If
                    a = xalglib.cmatrixludet(b, p, 2)
                    _TestResult = _TestResult And doc_test_complex(a, 25, 0.0001)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "matdet_t_5")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST solvesks_d_1
            '      Solving positive definite sparse system using Skyline (SKS) solver
            '
            _TestResult = true
            For _spoil_scenario = -1 To 3
                Try
                    ' 
                    '  This example demonstrates creation/initialization of the sparse matrix
                    '  in the SKS (Skyline) storage format and solution using SKS-based direct
                    '  solver.
                    ' 
                    '  First, we have to create matrix and initialize it. Matrix is created
                    '  in the SKS format, using fixed bandwidth initialization function.
                    '  Several points should be noted:
                    ' 
                    '  1. SKS sparse storage format also allows variable bandwidth matrices;
                    '     we just do not want to overcomplicate this example.
                    ' 
                    '  2. SKS format requires you to specify matrix geometry prior to
                    '     initialization of its elements with sparseset(). If you specified
                    '     bandwidth=1, you can not change your mind afterwards and call
                    '     sparseset() for non-existent elements.
                    '  
                    '  3. Because SKS solver need just one triangle of SPD matrix, we can
                    '     omit initialization of the lower triangle of our matrix.
                    ' 
                    Dim n As Integer = 4
                    Dim bandwidth As Integer = 1
                    Dim s As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    xalglib.sparsecreatesksband(n, n, bandwidth, s)
                    xalglib.sparseset(s, 0, 0, 2.0)
                    xalglib.sparseset(s, 0, 1, 1.0)
                    xalglib.sparseset(s, 1, 1, 3.0)
                    xalglib.sparseset(s, 1, 2, 1.0)
                    xalglib.sparseset(s, 2, 2, 3.0)
                    xalglib.sparseset(s, 2, 3, 1.0)
                    xalglib.sparseset(s, 3, 3, 2.0)

                    ' 
                    '  Now we have symmetric positive definite 4x4 system width bandwidth=1:
                    ' 
                    '      [ 2 1     ]   [ x0]]   [  4 ]
                    '      [ 1 3 1   ]   [ x1 ]   [ 10 ]
                    '      [   1 3 1 ] * [ x2 ] = [ 15 ]
                    '      [     1 2 ]   [ x3 ]   [ 11 ]
                    ' 
                    '  After successful creation we can call SKS solver.
                    ' 
                    Dim b() As Double = New Double(){4,10,15,11}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(b)
                    End If
                    Dim rep As sparsesolverreport = New XAlglib.sparsesolverreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){}
                    Dim isuppertriangle As Boolean = true
                    xalglib.sparsesolvesks(s, n, isuppertriangle, b, rep, x)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){1.0000,2.0000,3.0000,4.0000}, 0.00005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "solvesks_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            '
            ' TEST lincg_d_1
            '      Solution of sparse linear systems with CG
            '
            System.Console.WriteLine("150/151")
            _TestResult = true
            For _spoil_scenario = -1 To 3
                Try
                    ' 
                    '  This example illustrates solution of sparse linear systems with
                    '  conjugate gradient method.
                    '  
                    '  Suppose that we have linear system A*x=b with sparse symmetric
                    '  positive definite A (represented by sparsematrix object)
                    '          [ 5 1       ]
                    '          [ 1 7 2     ]
                    '      A = [   2 8 1   ]
                    '          [     1 4 1 ]
                    '          [       1 4 ]
                    '  and right part b
                    '      [  7 ]
                    '      [ 17 ]
                    '  b = [ 14 ]
                    '      [ 10 ]
                    '      [  6 ]
                    '  and we want to solve this system using sparse linear CG. In order
                    '  to do so, we have to create left part (sparsematrix object) and
                    '  right part (dense array).
                    ' 
                    '  Initially, sparse matrix is created in the Hash-Table format,
                    '  which allows easy initialization, but do not allow matrix to be
                    '  used in the linear solvers. So after construction you should convert
                    '  sparse matrix to CRS format (one suited for linear operations).
                    ' 
                    '  It is important to note that in our example we initialize full
                    '  matrix A, both lower and upper triangles. However, it is symmetric
                    '  and sparse solver needs just one half of the matrix. So you may
                    '  save about half of the space by filling only one of the triangles.
                    ' 
                    Dim a As sparsematrix = New XAlglib.sparsematrix() ' initializer can be dropped, but compiler will issue warning
                    xalglib.sparsecreate(5, 5, a)
                    xalglib.sparseset(a, 0, 0, 5.0)
                    xalglib.sparseset(a, 0, 1, 1.0)
                    xalglib.sparseset(a, 1, 0, 1.0)
                    xalglib.sparseset(a, 1, 1, 7.0)
                    xalglib.sparseset(a, 1, 2, 2.0)
                    xalglib.sparseset(a, 2, 1, 2.0)
                    xalglib.sparseset(a, 2, 2, 8.0)
                    xalglib.sparseset(a, 2, 3, 1.0)
                    xalglib.sparseset(a, 3, 2, 1.0)
                    xalglib.sparseset(a, 3, 3, 4.0)
                    xalglib.sparseset(a, 3, 4, 1.0)
                    xalglib.sparseset(a, 4, 3, 1.0)
                    xalglib.sparseset(a, 4, 4, 4.0)

                    ' 
                    '  Now our matrix is fully initialized, but we have to do one more
                    '  step - convert it from Hash-Table format to CRS format (see
                    '  documentation on sparse matrices for more information about these
                    '  formats).
                    ' 
                    '  If you omit this call, ALGLIB will generate exception on the first
                    '  attempt to use A in linear operations. 
                    ' 
                    xalglib.sparseconverttocrs(a)

                    ' 
                    '  Initialization of the right part
                    ' 
                    Dim b() As Double = New Double(){7,17,14,10,6}
                    If _spoil_scenario=0 Then
                        v_spoil_real = Double.NaN
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=1 Then
                        v_spoil_real = Double.PositiveInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=2 Then
                        v_spoil_real = Double.NegativeInfinity
                        spoil_vector_by_value(b, v_spoil_real)
                    End If
                    If _spoil_scenario=3 Then
                        spoil_vector_by_deleting_element(b)
                    End If

                    ' 
                    '  Now we have to create linear solver object and to use it for the
                    '  solution of the linear system.
                    ' 
                    '  NOTE: lincgsolvesparse() accepts additional parameter which tells
                    '        what triangle of the symmetric matrix should be used - upper
                    '        or lower. Because we've filled both parts of the matrix, we
                    '        can use any part - upper or lower.
                    ' 
                    Dim s As lincgstate = New XAlglib.lincgstate() ' initializer can be dropped, but compiler will issue warning
                    Dim rep As lincgreport = New XAlglib.lincgreport() ' initializer can be dropped, but compiler will issue warning
                    Dim x() As Double = New Double(){}
                    xalglib.lincgcreate(5, s)
                    xalglib.lincgsolvesparse(s, a, true, b)
                    xalglib.lincgresults(s, x, rep)

                    _TestResult = _TestResult And doc_test_int(rep.terminationtype, 1)
                    _TestResult = _TestResult And doc_test_real_vector(x, New Double(){1.000,2.000,1.000,2.000,1.000}, 0.005)
                    _TestResult = _TestResult And (_spoil_scenario=-1)
                Catch E As AlglibException
                    _TestResult = _TestResult And (_spoil_scenario<>-1)
                Catch E As Exception
                    Throw
                End Try
            Next _spoil_scenario
            If Not _TestResult Then
                System.Console.WriteLine("{0,-32} FAILED", "lincg_d_1")
            End If
            _TotalResult = _TotalResult And _TestResult


            System.Console.WriteLine("151/151")
        Catch E As Exception
            System.Console.WriteLine("Unhandled exception was raised!")
            Environment.Exit(1)
        End Try
        If _TotalResult Then
            Environment.Exit(0)
        End If
        Environment.Exit(1)
    End Sub
End Module
