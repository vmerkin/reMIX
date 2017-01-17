module mixtypes
  use mixdefs

  type Params_T
     ! place holder for model parameters
  end type Params_T

  type State_T
     real(mix_real), dimension(:,:,:),allocatable :: V
  end type State_T

  type Grid_T
     integer :: Np, Nt
     
     real(mix_real), dimension(:,:), allocatable :: x,y,t,p
  end type Grid_T

end module mixtypes
