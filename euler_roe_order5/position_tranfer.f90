module position_tranfer
    use lagrange_method
    implicit none
    public :: position_propagation_total
    contains
    function position_propagation_total(F, bc_total, x_delta, n) result(R_F)
        integer, intent(in) :: n
        real(kind=8), dimension(:,:), intent(in) :: F
        real(kind=8), dimension(:), intent(in) :: bc_total
        real(kind=8), intent(in) :: x_delta
        real(kind=8) ::  bc
        real(kind=8), allocatable, dimension(:) :: F_single
        real(kind=8), allocatable, dimension(:,:) :: R_F
        real(kind=8) u_i_11, u_i_12, u_i_22, u_i_32, u_i_21, u_i_31
        real(kind=8), allocatable, dimension(:,:) :: flux_l, flux_r, flux_final
        real(kind=8) :: u_delta
        integer :: var_num, i , xi ,var_i
        var_num = size(F, 2)
        allocate(R_F(n, var_num))
        allocate(flux_l(n,var_num), flux_r(n,var_num),flux_final(n,var_num))
        allocate(F_single(n))
        flux_array_initialization: do var_i = 1, var_num
            F_single = F(:,var_i)
            flux_l(:,var_i) = lagrange_interpolation(F_single, bc,.true.)
            flux_r(:,var_i) = lagrange_interpolation(F_single, bc,.false.)
        end do flux_array_initialization
        ! 准备把lagrange_interpolation函数写成直接返回虚拟点数组的函数，不再一个一个点插值
        flux_final = flux_l
        var_iteration: do var_i = 1,var_num
            position_iteration: do i = 1, n 
            if (i==1) then
                u_i_31 = flux_final(n-1,var_i)
                u_i_21 = flux_final(n,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(i+1,var_i)
                u_i_22 = flux_final(i+2,var_i)
                u_i_32 = flux_final(i+3,var_i)
            elseif(i==2) then
                u_i_31 = flux_final(n,var_i)
                u_i_21 = flux_final(i-1,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(i+1,var_i)
                u_i_22 = flux_final(i+2,var_i)
                u_i_32 = flux_final(i+3,var_i)
            elseif(i==n) then
                u_i_31 = flux_final(i-2,var_i)
                u_i_21 = flux_final(i-1,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(1,var_i)
                u_i_22 = flux_final(2,var_i)
                u_i_32 = flux_final(3,var_i)
            elseif(i==n-1) then
                u_i_31 = flux_final(i-2,var_i)
                u_i_21 = flux_final(i-1,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(i+1,var_i)
                u_i_22 = flux_final(1,var_i)
                u_i_32 = flux_final(2,var_i)
            elseif(i==n-2) then
                u_i_31 = flux_final(i-2,var_i)
                u_i_21 = flux_final(i-1,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(i+1,var_i)
                u_i_22 = flux_final(i+2,var_i)
                u_i_32 = flux_final(1,var_i)
            else
                u_i_31 = flux_final(i-2,var_i)
                u_i_21 = flux_final(i-1,var_i)
                u_i_11 = flux_final(i,var_i)
                u_i_12 = flux_final(i+1,var_i)
                u_i_22 = flux_final(i+2,var_i)
                u_i_32 = flux_final(i+3,var_i)
            end if
            ! 由于只有内点，故都采用六阶中心差分格式
            u_delta = 75.0_8/64.0_8*(u_i_12 - u_i_11) & 
                - 25.0_8/384.0_8*(u_i_22 - u_i_21) &
                + 3.0_8/640.0_8*(u_i_32 - u_i_31)
            R_F(i,var_i) = -u_delta / x_delta
            !print*, -u_delta / x_delta
            end do position_iteration
        end do var_iteration
    end function position_propagation_total
end module position_tranfer