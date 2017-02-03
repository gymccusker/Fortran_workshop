module CFL_mod

       use Types_mod

       implicit none

contains
        subroutine fd1d_heat_explicit_cfl( k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, cfl )
        real(kind=DP), intent(inout) :: cfl
        real(kind=DP) :: dx
        real(kind=DP) :: dt
        real(kind=DP), intent(in) :: k
        real(kind=DP), intent(in) :: t_max
        real(kind=DP), intent(in) :: t_min
        integer(kind=SI), intent(in) :: T_NUM
        real(kind=DP), intent(in) :: x_max
        real(kind=DP), intent(in) :: x_min
        integer(kind=SI), intent(in) :: X_NUM

        dx = ( x_max - x_min ) / dble( X_NUM - 1 )
        dt = ( t_max - t_min) / dble( T_NUM - 1 )

        cfl = k * dt / dx / dx

        write ( *, '(a)') ' '
        write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl

        end subroutine fd1d_heat_explicit_cfl





end module CFL_mod
