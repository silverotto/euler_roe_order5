    program order5
    use order5_module
    implicit none
    real(8)::x1_delta,x2_delta,t_delta,q1,q2,q3
    real(8),dimension(3)::err_array_1,err_array_2
    x1_delta=8.0D-2
    t_delta=1.0D-3
    x2_delta=4.0D-2
    err_array_1=order5_function(x1_delta,t_delta)
    err_array_2=order5_function(x2_delta,t_delta)
    q1=log10(err_array_1(1)/err_array_2(1))/log10(x1_delta/x2_delta);
    q2=log10(err_array_1(2)/err_array_2(2))/log10(x1_delta/x2_delta);
    q3=log10(err_array_1(3)/err_array_2(3))/log10(x1_delta/x2_delta);
    print*,"第一个算例的第一项误差:",err_array_1(1)
    print*,"第一个算例的第二项误差:",err_array_1(2)
    print*,"第一个算例的第三项误差:",err_array_1(3)
    print*,"第二个算例的第一项误差:",err_array_2(1)
    print*,"第二个算例的第二项误差:",err_array_2(2)
    print*,"第二个算例的第三项误差:",err_array_2(3)
    print*,"第一项格式精度:",q1
    print*,"第二项格式精度:",q2
    print*,"第三项格式精度:",q3
    end program order5