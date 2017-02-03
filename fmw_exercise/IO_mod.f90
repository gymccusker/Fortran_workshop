module IO_mod

      use Types_mod

      implicit none


contains
      subroutine r8mat_write( output_filename, table )
      implicit none

      integer(kind=SI) :: m
      integer(kind=SI) :: n

      integer(kind=SI) :: j
      character(len=*), intent(in) :: output_filename
      integer(kind=SI) :: output_unit
      character(len=30) :: string 
      real(kind=DP), dimension(:,:), intent(in) :: table
 
      m = size( table(:,:), 1 )
      n = size( table(:,:), 2 )

      output_unit = 10
      open( unit = output_unit, file = output_filename, status = 'replace' )

      write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'

      do j = 1, n
        write ( output_unit, string ) table(1:m, j)
      end do

      close( unit = output_unit )
    end subroutine r8mat_write

    subroutine r8vec_linspace ( a_first, a_last, a )

      implicit none

      integer(kind=SI) :: n
      real(kind=DP), dimension(:), intent(inout) :: a
      real(kind=DP), intent(in) :: a_first
      real(kind=DP), intent(in) :: a_last
      integer(kind=SI) :: i

      n = size(a)

      do i = 1, n
        a(i) = ( dble( n - i ) * a_first + dble( i - 1 ) * a_last ) / dble( n - 1 )
      end do

    end subroutine r8vec_linspace

    subroutine r8vec_write ( output_filename, x )

      implicit none

      integer(kind=SI) :: m
      integer(kind=SI) :: n

      integer(kind=SI) :: j
      character(len=*), intent(in) :: output_filename
      integer(kind=SI) :: output_unit
      real(kind=DP), dimension(:), intent(in) :: x

      n=size(x)

      output_unit = 11
      open( unit = output_unit, file = output_filename, status = 'replace' )

      do j = 1, n
        write ( output_unit, '(2x,g24.16)' ) x(j)
      end do

      close ( unit = output_unit )
  end subroutine r8vec_write




end module IO_mod

