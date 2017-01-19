module mixgeom
  use mixdefs
  use mixtypes

  implicit none
  
  contains
    real(mix_real) function cosDipAngle(t)   ! FIXME: modify to allow both hemispheres
      real(mix_real),intent(in) :: t
      cosDipAngle = -2.0D0*cos(t)/sqrt(1.0D0+3.0D0*cos(t)**2)
    end function cosDipAngle

    subroutine init_grid_fromXY(G,x,y)
      type(Grid_T),intent(inout) :: G
      real(mix_real), dimension(:,:), intent(in) :: x,y
      integer, dimension(2) :: dims

      ! set grid size
      dims = shape(x); G%Nt = dims(1); G%Np = dims(2)

      allocate(G%x(G%Nt,G%Np))
      allocate(G%y(G%Nt,G%Np))

      G%x = x
      G%y = y
    end subroutine init_grid_fromXY

    subroutine set_grid(G)
      type(Grid_T),intent(inout) :: G
      integer :: i,j

      allocate(G%t(G%Nt,G%Np))
      allocate(G%p(G%Nt,G%Np))

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

      ! define spherical angular coordinates
      G%t = asin(sqrt(G%x**2+G%y**2))
      G%p = modulo((atan2(G%y,G%x)+2*mix_pi),(2*mix_pi)) 
      ! note, this mangles phi at theta=0; Fix it, although we don't
      ! need it, since pole coordinates are never used
      G%p(1,:)=G%p(2,:)

      ! note, keep explicit size on the LHS to avoid compiler-dependent
      ! problems down the road
      G%dt(1:G%Nt-1,:) = G%t(2:G%Nt,:)-G%t(1:G%Nt-1,:)  
      G%dp(:,1:G%Np-1) = G%p(:,2:G%Np)-G%p(:,1:G%Np-1)
      G%dp(:,G%Np) = modulo(G%p(:,1)-G%p(:,G%Np),2*mix_pi)      ! fix up periodic


      ! note, unlike dp and dt above that are edge-centered, the things
      ! below are vortex centered; we just don't define them on the ends
      ! (e.g., ft(1,:) where we don't need them)
      G%ft(2:G%Nt,:) = 1.0_mix_real/(G%dt(2:G%Nt,:)+G%dt(1:G%Nt-1,:))/sin(G%t(2:G%Nt,:))
      G%fp(:,2:G%Np) = 1.0_mix_real/(G%dp(:,2:G%Np)+G%dp(:,1:G%Np-1)) ! note, dp(Np) defined above
      G%fp(:,1) = 1.0_mix_real/(G%dp(:,1)+G%dp(:,G%Np))  ! fix up periodic

      G%dtdt(2:G%Nt,:) = G%dt(2:G%Nt,:)/G%dt(1:G%Nt-1,:)-G%dt(1:G%Nt-1,:)/G%dt(2:G%Nt,:)
      G%dpdp(:,2:G%Np) = G%dp(:,2:G%Np)/G%dp(:,1:G%Np-1)-G%dp(:,1:G%Np-1)/G%dp(:,2:G%Np)
      G%dpdp(:,1) = G%dp(:,1)/G%dp(:,G%Np)-G%dp(:,G%Np)/G%dp(:,1) ! fix up periodic

      ! dip angle: compute and store for access later
      allocate(G%cosd(G%Nt,G%Np))    
      do j=1,G%Np
         do i=1,G%Nt
            G%cosd(i,j) = cosDipAngle(G%t(i,j))
         enddo
      enddo
    end subroutine set_grid
end module mixgeom
