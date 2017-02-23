! sets up and solves the stencil matrix
module mixsolver
  use mixdefs
  use mixtypes

  implicit none

  contains
    ! running index
    integer(kind=8) function K(j,i,Nj)
      integer, intent(in) :: i,j,Nj
      K=(i-1)*Nj+j
    end function K

    ! map 2 to -1, -2 to 1, and [-1,0,1] to themselves
    ! wt stands for wrap tripls
    integer function wt(q)
      integer,intent(in) :: q
      wt = nint(asin(sin(q*2*mix_pi/3))/asin(sin(2*mix_pi/3)))
    end function wt

    subroutine init_solver(P,G,S)
      type(Params_T), intent(in) :: P
      type(Grid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j

      ! this is the number of non-zeros in the matrix. Note, this
      ! depends on the stencil and boundary conditions
      S%nnz = G%Np*(G%Nt-2)*5 + G%Np + G%Np*(G%Np+1) 

      ! allocate RHS vector
      allocate(S%RHS(G%Nt*G%Np))
      ! allocate row index
      allocate(S%rowI(G%Nt*G%Np+1))
      ! Matrix 
      allocate(S%data(S%nnz))
      allocate(S%II(S%nnz))
      allocate(S%JJ(S%nnz))
      ! Equation terms
      allocate(S%F11(G%Np,G%Nt))
      allocate(S%F22(G%Np,G%Nt))
      allocate(S%F12(G%Np,G%Nt))
    end subroutine init_solver

    subroutine set_solver_terms(P,G,St,S)
      type(State_T), intent(in) :: St
      type(Params_T), intent(in) :: P  ! passing parameters just in case. Not used.
      type(Grid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      
      S%F11 = sin(G%t)*St%Vars(:,:,SIGMAP)/G%cosd**2
      S%F22 = St%Vars(:,:,SIGMAP)
      S%F12 =-St%Vars(:,:,SIGMAH)/G%cosd ! (assuming F21=-F12)
    end subroutine set_solver_terms

    subroutine set_solver_matrix_and_rhs(P,G,St,S,LLBC)
      type(State_T), intent(in) :: St
      type(Params_T), intent(in) :: P  ! passing parameters just in case. Not used.
      type(Grid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j,jm1,jp1,jj,a1,a2,q,u
      integer(kind=8) :: count,nextRow,nextRowI
      real(mix_real) :: dF12t,dF12p
      real(mix_real) :: d(-2:2),c(-2:2)
      real(mix_real),dimension(:),intent(in) :: LLBC

      ! init to zero. This is important, since we're only filling in
      ! non-zero elements in the matrix construction
      S%RHS = 0.0D0   
      S%data=0.0D0
      S%II = 0.0D0
      S%JJ = 0.0D0
      
      count=1
      nextRowI=1
      ! make a note about the fact that the order is corrects for the CSR format

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Pole boundary 
      do j=1,G%Np
         S%data(count) = 1.0D0
         S%II(count)   = K(j,1,G%Np)
         S%JJ(count)   = K(j,1,G%Np)
         count=count+1

         do jj=1,G%Np   ! with (Np,Nt) definition of the grid, this is a natural alignment
            S%data(count) = -G%dp(jj,2)/(2*mix_pi)
            S%II(count) = K(j,1,G%Np)
            S%JJ(count) = K(jj,2,G%Np)
            count=count+1
         enddo
         
         ! These are points on the pole boundary Thus, each row of the
         ! matrix has (Np+1) entries: for the point itself + Np
         ! entries for averaging the points at i=2 (see above).
         ! Therefore, the rows start with indices 1, Np+2, 2Np+3, etc.
         S%rowI(K(j,1,G%Np)) = nextRowI       !1+(j-1)*(G%Np+1)
         nextRowI=nextRowI+G%Np+1
      enddo
      ! note, not setting RHS because it's initializaed to zero anyway
      ! end pole boundary
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! inner block
      do i=2,G%Nt-1  ! excluding pole and low lat boundaries
         do j=1,G%Np
            jm1 = merge(G%Np,j-1,j.eq.1)   ! maps j-1=0 to j-1=Np, otherwise returns j-1
            jp1 = merge(1,j+1,j.eq.G%Np)   ! similar for j+1

            ! the above functions involve if statements, which I don't
            ! want to pack inside inner loop. Do this instead: these
            ! wonderful functions of my invention maps j-1=0 to Np,
            ! otherwise j-1 maps to itself. Similarly, j+1=Np is
            ! mapped to j+1=1, otherwise nothing's done
!            jm1 = modulo(j-1,G%Np)+G%Np*(1-int(ceiling(real(j-1)/G%Np)))  
!            jp1 = modulo(j+1,G%Np)+G%Np*(1-int(ceiling(real(modulo(j+1,G%Np))/G%Np)))

            ! derivatives for off diagonal conductance terms
            dF12p = G%fp(j,i)*( G%dp(jm1,i)/G%dp(j,i)*S%F12(jp1,i) + G%dpdp(j,i)*S%F12(j,i) - G%dp(j,i)/G%dp(jm1,i)*S%F12(jm1,i) )
            dF12t = G%ft(j,i)*( G%dt(j,i-1)/G%dt(j,i)*S%F12(j,i+1) + G%dtdt(j,i)*S%F12(j,i) - G%dt(j,i)/G%dt(j,i-1)*S%F12(j,i-1) )

            ! for periodic boundaries jm1>jp1 maybe the case, which
            ! breaks the order of things. So, let's figure it out

            a1 = (jp1-j)/abs(jp1-j)
            a2 = (jm1-j)/abs(jm1-j)
            q  = nint(0.5*(a1+a2))
            
            d(-2) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1)+&
                 dF12p*G%ft(j,i)*G%dt(j,i)/G%dt(j,i-1)
            c(-2) = K(j,i-1,G%Np)

            d(wt(-q-1)) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i)-&
                 dF12t*G%fp(j,i)*G%dp(j,i)/G%dp(jm1,i)
            c(wt(-q-1)) = K(jm1,i,G%Np)

            d(-q) = -G%ft(j,i)*( (S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)+(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1) ) - &
                 G%fp(j,i)/sin(G%t(j,i))**2*( (S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i) ) + &
                 dF12t*G%fp(j,i)*G%dpdp(j,i)-&
                 dF12p*G%ft(j,i)*G%dtdt(j,i) 
            c(-q) = K(j,i,G%Np)

            d(wt(-q+1)) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+&
                 dF12t*G%fp(j,i)*G%dp(jm1,i)/G%dp(j,i)
            c(wt(-q+1)) = K(jp1,i,G%Np)
            
            d(2) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)-&
                 dF12p*G%ft(j,i)*G%dt(j,i-1)/G%dt(j,i)
            c(2) = K(j,i+1,G%Np)

            do u=-2,2
               s%data(count) = d(u)
               S%JJ(count)   = c(u)
               S%II(count)   = K(j,i,G%Np)
               count=count+1
            enddo

            ! At this point, coming out of the pole boundary condition, we have: S%rowI(K(j,i,G%Np)-1) = G%Np*(G%Np+1)+1
!            S%rowI(K(j,i,G%Np)) = S%rowI(K(j,i,G%Np)-1)+G%Np+1 + 5*(K(j,i,G%Np)-1)
!            S%rowI(K(j,i,G%Np)) = G%Np*(G%Np+1)+1 + 5*(K(j,i,G%Np)-1)
            S%rowI(K(j,i,G%Np)) = nextRowI
            nextRowI = nextRowI+5
            S%RHS(K(j,i,G%Np)) = St%Vars(j,i,FAC)
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! low lat boundary
      do j=1,G%Np
         S%data(count)  = 1.0D0
         S%II(count)    = K(j,G%Nt,G%Np)
         S%JJ(count)    = K(j,G%Nt,G%Np)
         count=count+1

         S%rowI(K(j,G%Nt,G%Np)) = nextRowI 
         nextRowI=nextRowI+1
         S%RHS(K(j,G%Nt,G%Np)) = LLBC(j)
      enddo
      S%rowI(K(G%Np,G%Nt,G%Np)+1) = nextRowI   ! the dummy last element as required by the CSR storage
    end subroutine set_solver_matrix_and_rhs
    ! end low lat boundary
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mixsolver
