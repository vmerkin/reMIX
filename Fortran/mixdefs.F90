! MIX definitions/constants

module mixdefs
  implicit none

  !Define variable precisions
  integer, parameter :: mix_single = kind(1.0)
  integer, parameter :: mix_double = kind(1.0D0)
  
  integer, parameter :: mix_real = mix_double
  integer, parameter :: mix_io_real  = mix_single !Precision for IO
  
  integer, parameter :: strLen = 100 !Default size for strings
  
  real (mix_real), parameter :: mix_pi = 4.0D0*atan(1.0D0)
  
  integer, parameter :: nVars = 4 ! change together wiht the enumerator below
  enum, bind(C)
     enumerator :: POT=1,FAC,SIGMAP,SIGMAH
  end enum
end module mixdefs
