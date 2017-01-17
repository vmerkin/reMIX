module mixgeom
  use mixdefs
  use mixtypes

  implicit none
  
  contains
    
    subroutine set_grid(G)
      type(Grid_T),intent(inout) :: G

      if (.not.allocated(G%dp)) then
         allocate(G%dp(G%Nt,G%Np-1))
      endif

      if (.not.allocated(G%dt)) then
         allocate(G%dt(G%Nt-1,G%Np))
      endif

      G%dp = G%p(:,2:G%Np)-G%p(:,1:G%Np-1)
      G%dt = G%t(2:G%Nt,:)-G%t(1:G%Nt-1,:)
    end subroutine set_grid
end module mixgeom
