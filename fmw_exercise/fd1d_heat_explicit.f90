program fd1d_heat_explicit_prb

      use Types_mod

      implicit none

      integer(kind=SI), parameter :: T_NUM = 201
      integer(kind=SI), parameter :: X_NUM = 21
      
      real(kind=DP) :: cfl
      real(kind=DP) :: dt
      real(kind=DP), dimension(:), allocatable :: h
      real(kind=DP), dimension(:), allocatable :: h_new
      ! the "matrix" stores all x-values for all t-values
      ! remember Fortran is column major, meaning that rows are contiguous
      real(kind=DP), dimension(:,:), allocatable :: hmat
      integer(kind=SI) :: i
      integer(kind=SI) :: j
      real(kind=DP) :: k

      real(kind=DP), dimension(:), allocatable :: t
      real(kind=DP) :: t_max
      real(kind=DP) :: t_min
      real(kind=DP), dimension(:), allocatable :: x
      real(kind=DP) :: x_max
      real(kind=DP) :: x_min

      allocate( h(1:X_NUM), h_new(1:X_NUM), hmat(1:X_NUM,1:T_NUM), t(1:T_NUM), x(1:X_NUM) ) 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the FD1D_HEAT_EXPLICIT library.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_TEST01:'
      write ( *, '(a)' ) '  Compute an approximate solution to the time-dependent'
      write ( *, '(a)' ) '  one dimensional heat equation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    dH/dt - K * d2H/dx2 = f(x,t)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Run a simple test case.'

      ! heat coefficient
      k = 0.002_DP

      ! the x-range values
      x_min = 0.0_DP
      x_max = 1.0_DP
      ! X_NUM is the number of intervals in the x-direction
      call r8vec_linspace( x_min, x_max, x )

      ! the t-range values. integrate from t_min to t_max
      t_min = 0.0_DP
      t_max = 80.0_DP

      ! T_NUM is the number of intervals in the t-direction
      dt = ( t_max - t_min ) / dble( T_NUM - 1 )
      call r8vec_linspace( t_min, t_max, t )

      ! get the CFL coefficient
      call fd1d_heat_explicit_cfl( k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, cfl )

     if ( 0.5_DP <= cfl ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
        write ( *, '(a)' ) '  CFL condition failed.'
        write ( *, '(a)' ) '  0.5 <= K * dT / dX / dX = CFL.'
        stop
      end if

      ! set the initial condition
      do j = 1, X_NUM
        h(j) = 50.0_DP
      end do

      ! set the bounday condition
      h(1) = 90.0_DP
      h(X_NUM) = 70.0_DP

      ! initialise the matrix to the initial condition
      do i = 1, X_NUM
        hmat(i, 1) = h(i)
      end do

      ! the main time integration loop 
      do j = 2, T_NUM
        call fd1d_heat_explicit( x, t(j-1), dt, cfl, h, h_new )

        do i = 1, X_NUM
          hmat(i, j) = h_new(i)
          h(i) = h_new(i)
        end do
      end do

      ! write data to files
      call r8mat_write( 'h_test01.txt', hmat )
      call r8vec_write( 't_test01.txt', t )
      call r8vec_write( 'x_test01.txt', x )

      deallocate( h, h_new, hmat, t, x )

    contains

    function func( j, x ) result ( d )
      implicit none
      
      integer(kind=SI), intent(in) :: j!, X_NUM
      real(kind=DP) :: d
      !real(kind=DP), intent(in) :: x(X_NUM)
      real(kind=DP), dimension(:) :: x

      d = 0.0_DP
    end function func

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

    subroutine fd1d_heat_explicit_cfl( k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, cfl )

      implicit none

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
      dt = ( t_max - t_min ) / dble( T_NUM - 1 )

      cfl = k * dt / dx / dx

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl

    end subroutine fd1d_heat_explicit_cfl

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

end program fd1d_heat_explicit_prb

