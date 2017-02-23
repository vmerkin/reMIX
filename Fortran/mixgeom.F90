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
      dims = shape(x); G%Nt = dims(2); G%Np = dims(1)

      allocate(G%x(G%Np,G%Nt))
      allocate(G%y(G%Np,G%Nt))

      G%x = x  ! note the transposes to conform to the new definition (Np,Nt)
      G%y = y
    end subroutine init_grid_fromXY

    subroutine set_grid(G)
      type(Grid_T),intent(inout) :: G
      integer :: i,j

      allocate(G%t(G%Np,G%Nt))
      allocate(G%p(G%Np,G%Nt))

      ! note allocating everything with the same size
      ! careful with unused points
      allocate(G%dp(G%Np,G%Nt))    ! True size (Np-1,Nt)
      allocate(G%dt(G%Np,G%Nt))    ! True size (Np,Nt-1)

      ! things needed for the solver and dependent only on the grid, not conductances
      ! same caveats about the sizes
      allocate(G%ft(G%Np,G%Nt))
      allocate(G%fp(G%Np,G%Nt))
      allocate(G%dtdt(G%Np,G%Nt))
      allocate(G%dpdp(G%Np,G%Nt))

      ! define spherical angular coordinates
      G%t = asin(sqrt(G%x**2+G%y**2))
      G%p = modulo((atan2(G%y,G%x)+2*mix_pi),(2*mix_pi)) 
      ! note, this mangles phi at theta=0; Fix it, although we don't
      ! need it, since pole coordinates are never used
      G%p(:,1)=G%p(:,2)

      ! FIXME: use fortran shift/merge functions to treat periodic boundaries

      ! note, keep explicit size on the LHS to avoid compiler-dependent
      ! problems down the road
      G%dt(:,1:G%Nt-1) = G%t(:,2:G%Nt)-G%t(:,1:G%Nt-1)  
      G%dp(1:G%Np-1,:) = G%p(2:G%Np,:)-G%p(1:G%Np-1,:)
      G%dp(G%Np,:) = modulo(G%p(1,:)-G%p(G%Np,:),2*mix_pi)      ! fix up periodic

      ! note, unlike dp and dt above that are edge-centered, the things
      ! below are vortex centered; we just don't define them on the ends
      ! (e.g., ft(1,:) where we don't need them)
      G%ft(:,2:G%Nt) = 1.0D0/(G%dt(:,2:G%Nt)+G%dt(:,1:G%Nt-1))/sin(G%t(:,2:G%Nt))
      G%fp(2:G%Np,:) = 1.0D0/(G%dp(2:G%Np,:)+G%dp(1:G%Np-1,:)) ! note, dp(Np) defined above
      G%fp(1,:) = 1.0D0/(G%dp(1,:)+G%dp(G%Np,:))  ! fix up periodic

      ! this should work but since we don't use the last element (G%Nt) let's use the next line instead to avoid possibly dividing by zero (G%dt(:,2:G%Nt))
!      G%dtdt(:,2:G%Nt) = G%dt(:,2:G%Nt)/G%dt(:,1:G%Nt-1)-G%dt(:,1:G%Nt-1)/G%dt(:,2:G%Nt)
      G%dtdt(:,2:G%Nt-1) = G%dt(:,2:G%Nt-1)/G%dt(:,1:G%Nt-2)-G%dt(:,1:G%Nt-2)/G%dt(:,2:G%Nt-1)
      G%dpdp(2:G%Np,:) = G%dp(2:G%Np,:)/G%dp(1:G%Np-1,:)-G%dp(1:G%Np-1,:)/G%dp(2:G%Np,:)
      G%dpdp(1,:) = G%dp(1,:)/G%dp(G%Np,:)-G%dp(G%Np,:)/G%dp(1,:) ! fix up periodic

      ! dip angle: compute and store for access later
      allocate(G%cosd(G%Np,G%Nt))    
      do i=1,G%Nt
         do j=1,G%Np
            G%cosd(j,i) = cosDipAngle(G%t(j,i))
         enddo
      enddo
    end subroutine set_grid
end module mixgeom
