module mixgeom
  use mixdefs
  use mixtypes

  implicit none
  
  contains
    
    subroutine set_grid(G)
      type(Grid_T),intent(inout) :: G

      ! note allocating everything with the same size
      ! careful with unused points
      allocate(G%dp(G%Nt,G%Np))    ! True size (Nt,Np-1)
      allocate(G%dt(G%Nt,G%Np))    ! True size (Nt-1,Np)

      ! things needed for the solver and dependent only on the grid, not conductances
      ! same caveats about the sizes
      allocate(G%ft(G%Nt,G%Np))
      allocate(G%fp(G%Nt,G%Np))
      allocate(G%dtdt(G%Nt,G%Np))
      allocate(G%dpdp(G%Nt,G%Np))

      ! note, keep explicit size on the LHS to avoid compiler-dependent
      ! problems down the road
      G%dp(:,1:G%Np-1) = G%p(:,2:G%Np)-G%p(:,1:G%Np-1)
      G%dp(:,G%Np) = mod(G%p(:,1)-G%p(:,G%Np),2*mix_pi)      ! fix up periodic
      G%dt(1:G%Nt-1,:) = G%t(2:G%Nt,:)-G%t(1:G%Nt-1,:)  

      ! note, unlike dp and dt above that are edge-centered, the things
      ! below are vortex centered; we just don't define them on the ends
      ! (e.g., ft(1,:) where we don't need them)
!      G%ft(2:G%Nt,:) = 1.0_mix_real/(G%dt(2:G%Nt,:)+G%dt(1:G%Nt-1,:))/&
!      sin(G%t(2:G%Nt,:))
!      G%fp(:,2:G%Np-1) = 1.0_mix_real/(G%dp(:,2:G%Np-1)+G%dp(:,1:G%Np-2))
!      G%fp(:,1) = 1.0_mix_real/(G%dp(:,1)+G%dp(:,G%Np))  ! fix up periodic

!      G%dtdt(2:G%Nt,:) = G%dt(2:G%Nt,:)/G%dt(1:G%Nt-1,:)-&
!      G%dt(1:G%Nt-1,:)/G%dt(2:G%Nt,:)

!      G%dpdp(:,2:G%Np) = G%dp(:,2:G%Np)/G%dp(:,1:G%Np-1)-&
!      G%dp(:,1:G%Np-1)/G%dp(:,2:G%Np)

!      G%dpdp(:,1) = G%dp(:,1)/G%dp(:,1:G%Np-1)-&
!      G%dp(:,1:G%Np-1)/G%dp(:,2:G%Np)

    end subroutine set_grid
end module mixgeom
