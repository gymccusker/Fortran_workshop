module RHS_mod


      use Types_mod

      implicit none
      
        

contains
  function func( j, x ) result ( d )
      !implicit none
      
      integer(kind=SI), intent(in) :: j!, X_NUM
      real(kind=DP) :: d
      !real(kind=DP), intent(in) :: x(X_NUM)
      real(kind=DP), dimension(:) :: x

      d = 0.0_DP

    end function func
end module RHS_mod
