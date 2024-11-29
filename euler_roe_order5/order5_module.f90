module order5_module
    use, intrinsic :: iso_fortran_env ! 绘图模块初始化
    use pyplot_module
    use position_tranfer
    use conservation_to_flux
    implicit none
    public :: lagrange_interpolation
    contains
    function order5_function(x_delta, t_delta) result(err_array)
 
    real(kind=8), intent(in) :: x_delta, t_delta
    real(kind=8) :: pi
    real(kind=8) :: rho, u, p, E
    real(kind=8) :: gamma
    real(kind=8) :: err_array(3)
    real(kind=8) :: err1,err2,err3
    integer :: x_total, t_total, n, t_max, var_num
    integer :: ti, xi
    real(kind=8) :: x
 
    real(kind=8) :: u_l, u_r, u_delta
    real(kind=8), allocatable, dimension(:) :: bc_t0, bc_t1, bc_t2, bc_t3
    real(kind=8), allocatable, dimension(:,:) :: R_F0, R_F1, R_F2, R_F3
    real(kind=8), allocatable, dimension(:,:) :: F0, F1, F2, F3, Q0, Q1, Q2, Q3
    real(kind=8), allocatable, dimension(:,:,:) :: Q, F
    real(kind=8), allocatable, dimension(:) :: x_plt
    real(kind=8), allocatable, dimension(:,:) :: y_err, y_nu, y_ac
    character(len=10) :: temp_string
    type(pyplot) :: plt1, plt2, plt3
    x_total = 2
    t_total = 1
    n = x_total / x_delta
    t_max = t_total / t_delta
    var_num = 3 ! 守恒变量数，分别为rho , rhou, rhoE
    gamma = 1.4 ! 比热比
    pi = acos(-1.0D0) ! π
    allocate(bc_t0(var_num))
    allocate(bc_t1(var_num))
    allocate(bc_t2(var_num))
    allocate(bc_t3(var_num))
    allocate(F0(n, var_num))
    allocate(F1(n, var_num))
    allocate(F2(n, var_num))
    allocate(F3(n, var_num))
    allocate(R_F0(n, var_num))
    allocate(R_F1(n, var_num))
    allocate(R_F2(n, var_num))
    allocate(R_F3(n, var_num))
    allocate(Q(0:t_max, n, var_num)) ! 从0时刻推进至t_total时刻，即共有t_max+1个时间步
    allocate(F(0:t_max, n, var_num))
    allocate(x_plt(n))
    allocate(y_err(n, var_num))
    allocate(y_nu(n, var_num))
    allocate(y_ac(n, var_num))
    Q_initialization: do concurrent (xi = 1:n) ! 初始条件初始化
    !Q_initialization: do xi = 1,n ! 初始条件初始化
        x = x_delta * (xi - 1)
        rho = 1 + 0.5 * sin(pi * x)
        u = 1.0D0
        p = 1.0D0
        E = p / ((gamma - 1) * rho) + 0.5 * u ** 2
        Q(0, xi, 1) = rho
        Q(0, xi, 2) = rho * u
        Q(0, xi, 3) = rho * E
        F(0, xi, 1) = rho * u
        F(0, xi, 2) = rho * u * u + p
        F(0, xi, 3) = u * (rho * E + p)
    end do Q_initialization
    time_stepping: do ti = 1, t_max 
        F0 = F(ti - 1, :, :)
        Q0 = Q(ti - 1, :, :)
        bc_t0 = -sin((ti - 1) * t_delta) - sin(2 * (ti - 1) * t_delta)
        R_F0 = position_propagation_total(F0, bc_t0, x_delta, n)
        Q1 = Q0 + 0.5 * t_delta * R_F0 ! 这里得到推进的应该是守恒量，要转为通量的形式
        F1 = conservation_to_flux_function(Q1)
        bc_t1 = -sin((ti - 1 + 0.5_8) * t_delta) - sin(2 * (ti - 1 + 0.5_8) * t_delta)
        R_F1 = position_propagation_total(F1, bc_t1, x_delta, n)
        Q2 = Q0 + 0.5_8 * t_delta * R_F1
        F2 = conservation_to_flux_function(Q2)
        bc_t2 = -sin((ti - 1 + 0.5) * t_delta) - sin(2 * (ti - 1 + 0.5) * t_delta)
        R_F2 = position_propagation_total(F2, bc_t2, x_delta, n)
        Q3 = Q0 + t_delta * R_F2
        F3 = conservation_to_flux_function(Q3)
        bc_t3 = -sin((ti - 1 + 1) * t_delta) - sin(2 * (ti - 1 + 1) * t_delta)
        R_F3 = position_propagation_total(F3, bc_t3, x_delta, n)
    
        Q(ti, :, :) = Q0 + 1.0_8 / 6.0_8 * t_delta * (R_F0 + 2.0_8 * R_F1 + 2.0_8 * R_F2 + R_F3)
        F(ti, :, :) = conservation_to_flux_function(Q(ti, :, :))
    end do time_stepping
    
    err_initialization: do concurrent (xi = 1:n)
    !err_initialization: do xi = 1, n
        x_plt(xi) = (xi - 1) * x_delta
        y_nu(xi,1) = Q(t_max, xi,1)
        y_ac(xi,1) = 1+0.5 * sin(pi * (x_plt(xi)-t_total) )
        y_err(xi,1) = y_nu(xi,1)-y_ac(xi,1)
        y_nu(xi,2) = Q(t_max, xi,2)
        y_ac(xi,2) = 1+0.5 * sin(pi * (x_plt(xi)-t_total) )
        y_err(xi,2) = y_nu(xi,2)-y_ac(xi,2)
        y_nu(xi,3) = Q(t_max, xi,3)
        y_ac(xi,3) = 1/(gamma - 1) + 0.5 + 0.25_8 * sin(pi * (x_plt(xi)-t_total) )
        y_err(xi,3) = y_nu(xi,3)-y_ac(xi,3)
    end do err_initialization
    err1 = maxval(abs(y_err(:,1)) )
    err2 = maxval(abs(y_err(:,2)) )
    err3 = maxval(abs(y_err(:,3)) )
    !err1=norm2(y_err(:,1))
    !err2=norm2(y_err(:,2))
    !err3=norm2(y_err(:,3))
    err_array(1) = err1
    err_array(2) = err2
    err_array(3) = err3
    write(temp_string, '(F10.5)') x_delta
    
    call plt1%initialize(grid=.true., title="numerical result")
    call plt1%add_plot(x_plt, y_nu(:,3), label='res', linestyle='r-o', markersize=3, linewidth=1)
    call plt1%add_str('plt.show(block=False)') 
    call plt1%finish_ops()
    
    call plt2%initialize(grid=.true., title="accuracy result")
    call plt2%add_plot(x_plt, y_ac(:,3), label='res', linestyle='g-o', markersize=3, linewidth=1)
    call plt2%add_str('plt.show(block=False)') 
    call plt2%finish_ops()
    
    call plt3%initialize(grid=.true., title="error")
    call plt3%add_plot(x_plt, y_err(:,3), label='res', linestyle='b-o', markersize=3, linewidth=1)
    call plt3%add_str('plt.show(block=False)') 
    call plt3%finish_ops()
    
    !call plt1%savefig('x_delta='//temp_string//'numerical_result.png', pyfile='x_delta='//temp_string//'numerical_result.py')
    !call plt2%savefig('x_delta='//temp_string//'accuracy_result.png', pyfile='x_delta='//temp_string//'accuracy_result.py')
    !call plt3%savefig('x_delta='//temp_string//'error_result.png', pyfile='x_delta='//temp_string//'error_result.py')
    
    print*, 'done'
    deallocate(bc_t0)
    deallocate(bc_t1)
    deallocate(bc_t2)
    deallocate(bc_t3)
    deallocate(F0)
    deallocate(F1)
    deallocate(F2)
    deallocate(F3)
    deallocate(R_F0)
    deallocate(R_F1)
    deallocate(R_F2)
    deallocate(R_F3)
    deallocate(Q)
    deallocate(x_plt)
    deallocate(y_err)
    deallocate(y_nu)
    deallocate(y_ac)
    end function order5_function
end module order5_module