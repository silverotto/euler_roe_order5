    module lagrange_method
    implicit none
    public:: lagrange_interpolation
    contains
    function lagrange_interpolation(u,bc,is_direction_left) result(flux_array) !���lagrage��ֵ
        real(8),intent(in) :: bc  
        real(8),intent(in) :: u(:)  !u��ԭ�����趨Ϊ��0~n-1,����n����
        logical,intent(in) :: is_direction_left
        integer:: n,i
        real(8) :: coeffs(5),u_neighbor(5),u_i
        real(8), allocatable, dimension(:) :: flux_array
        n = size(u)
        allocate(flux_array(n))

          !���� �߽������������Ա߽磬��ȫ�����ڵ㣬����Ҫ�Ա߽���д���
        coeffs = (/ 0.0234375_8, -0.15625_8, 0.703125_8, 0.46875_8, -0.0390625_8 /)
        flux_array_generation: do i = 1, n
            ! ��ֵ��ֵ
            if(is_direction_left) then
                if (i == 1) then
                    u_neighbor = (/ u(n-2), u(n-1), u(n), u(i), u(i+1) /)
                ! ��߽���һ
                elseif (i == 2) then
                    u_neighbor = (/ u(n-1), u(n), u(i-1), u(i), u(i+1) /)
                ! ��߽��Ҷ�
                elseif (i == 3) then
                    u_neighbor = (/ u(n), u(i-2), u(i-1), u(i), u(i+1) /)
                ! �ұ߽�
                elseif (i == n) then
                    u_neighbor = (/ u(i-3), u(i-2), u(i-1), u(i), u(1) /)
                ! �м��
                else
                    u_neighbor = (/ u(i-3), u(i-2), u(i-1), u(i), u(i+1) /)
                end if
            ! ��ֵ��ֵ
            else
                if (i==1) then
                    u_neighbor = (/ u(n-1), u(n), u(i), u(i+1), u(i+2) /)
                ! ��߽���һ
                elseif (i == 2) then
                    u_neighbor = (/ u(n), u(i-1), u(i), u(i+1), u(i+2) /)
                ! �ұ߽�
                elseif (i == n) then
                    u_neighbor = (/ u(i-2), u(i-1), u(i), u(1), u(2) /)
                ! �ұ߽���һ
                elseif (i == n-1) then
                    u_neighbor = (/ u(i-2), u(i-1), u(i), u(i+1), u(1) /)
                ! �м��
                else
                    u_neighbor = (/ u(i-2), u(i-1), u(i), u(i+1), u(i+2) /)
                end if   
            end if
            u_i=dot_product(coeffs,u_neighbor)
            flux_array(i) = u_i
        end do flux_array_generation
        
        !print*,u_i
    end function lagrange_interpolation
    end module lagrange_method