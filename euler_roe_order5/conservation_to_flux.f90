    module conservation_to_flux
    implicit none
    public:: conservation_to_flux_function
    contains
    function conservation_to_flux_function(Q) result(F) !五阶lagrage插值
        real(kind=8), dimension(:,:), intent(in) :: Q
        integer:: n,var_num
        real(kind=8) :: gamma
        real(kind=8), allocatable, dimension(:) :: rho,u,p,E ! 若var_num不为3，则需要再修改该部分
        real(kind=8), allocatable, dimension(:) :: A,B,C
        real(kind=8), allocatable, dimension(:,:) :: F
        n=size(Q,1)
        var_num=size(Q,2)
        gamma=1.4
        allocate(rho(n))
        allocate(u(n))
        allocate(p(n))
        allocate(E(n))
        allocate(A(n))
        allocate(B(n))
        allocate(C(n))
        allocate(F(n,var_num))
        A=Q(:,1)
        B=Q(:,2)
        C=Q(:,3)
        rho=A
        u=B/A
        p=-(gamma-1)*(B*B-2*A*C)/(2*A)
        E=C/A
        F(:,1)=rho*u
        F(:,2)=rho*u*u+p
        F(:,3)=u * (rho * E + p)
    end function conservation_to_flux_function
    end module conservation_to_flux