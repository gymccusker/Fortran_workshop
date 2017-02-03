module Solver_mod

      use Types_mod
      use RHS_mod

      implicit none


contains

subroutine fd1d_heat_explicit( x, t, dt, cfl, h, h_new )
      implicit none

      !integer(kind=SI), intent(in) :: X_NUM

      real(kind=DP), intent(in) :: cfl
      real(kind=DP), intent(in) :: dt
      !real(kind=DP), intent(in) :: h(X_NUM)
      real(kind=DP), dimension(:), intent(in) :: h
      !real(kind=DP), intent(inout) :: h_new(X_NUM)
      real(kind=DP), dimension(:), intent(inout) :: h_new
      integer(kind=SI) :: j
      real(kind=DP), intent(in) :: t
      !real(kind=DP), intent(in) :: x(X_NUM)
      real(kind=DP), dimension(:), intent(in) :: x
      !real(kind=DP) :: f(X_NUM)
      real(kind=DP) :: f(size(x))

      do j = 1, size(x)
        f(j) = func( j, x )
      end do

      h_new(1) = 0.0_DP

      do j = 2, size(x) - 1
        h_new(j) = h(j) + dt * f(j) + cfl * ( h(j-1) - 2.0_DP * h(j) + h(j+1) )
      end do

      ! set the boundary conditions again
      h_new(1) = 90.0_DP
      h_new(size(x)) = 70.0_DP
    end subroutine fd1d_heat_explicit




end module Solver_mod
