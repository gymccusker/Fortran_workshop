module Types_mod

	use, intrinsic :: iso_fortran_env

	implicit none

	public :: SP, DP, SI, DI

	integer, parameter :: SP = REAL32
	integer, parameter :: DP = REAL64
	integer, parameter :: SI = INT32
	integer, parameter :: DI = INT64

contains

end module Types_mod
