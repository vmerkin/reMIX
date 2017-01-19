! sets up and solves the stencil matrix

module mixsolver
  use mixdefs
  use mixtypes

  implicit none

  contains
    ! running index
    integer(kind=8) function K(i,j,Ni)
      integer, intent(in) :: i,j,Ni
      K=(j-1)*Ni+i
    end function K

    subroutine init_solver(P,G,St,S)
      type(State_T), intent(in) :: St
      type(Params_T), intent(in) :: P
      type(Grid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j

      ! this is the number of non-zeros in the matrix. Note, this
      ! depends on the stencil and boundary conditions
      S%nnz = G%Np*(G%Nt-2)*5 + G%Np + G%Np*(G%Np+1) 
      print *,S%nnz

      ! allocate RHS vector
      allocate(S%RHS(G%Nt*G%Np))
      ! Matrix 
      allocate(S%data(S%nnz))
      allocate(S%II(S%nnz))
      allocate(S%JJ(S%nnz))
      ! Equation terms
      allocate(S%F11(G%Nt,G%Np))
      allocate(S%F22(G%Nt,G%Np))
      allocate(S%F12(G%Nt,G%Np))
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
      integer :: i,j,jm1,jp1,jj
      integer(kind=8) :: count
      real(mix_real) :: dF12t,dF12p
      real(mix_real),dimension(:),intent(in) :: LLBC

      ! init to zero. This is important, since we're only filling in
      ! non-zero elements in the matrix construction
      S%RHS = 0.0_mix_real   
      S%data=0.0_mix_real
      S%II = 0.0_mix_real
      S%JJ = 0.0_mix_real
      
      count=1
      do j=1,G%Np
         jm1 = merge(G%Np,j-1,j.eq.1)   ! maps j-1=0 to j-1=Np, otherwise returns j-1
         jp1 = merge(1,j+1,j.eq.G%Np)   ! similar for j+1

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Pole boundary 
         S%data(count) = 1.0_mix_real
         S%II(count)   = K(1,j,G%Nt)
         S%JJ(count)   = K(1,j,G%Nt)
         count=count+1

         do jj=1,G%Np
            S%data(count) = -G%dp(2,jj)/(2*mix_pi)
            S%II(count) = K(1,j,G%Nt)
            S%JJ(count) = K(2,jj,G%Nt)
            count=count+1
         enddo
         ! note, not setting RHS because it's initializaed to zero anyway
         ! end pole boundary
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do i=2,G%Nt-1  ! excluding pole and low lat boundaries
            ! derivatives for off diagonal conductance terms
            dF12p = G%fp(i,j)*( G%dp(i,jm1)/G%dp(i,j)*S%F12(i,jp1) + G%dpdp(i,j)*S%F12(i,j) - G%dp(i,j)/G%dp(i,jm1)*S%F12(i,jm1) )
            dF12t = G%ft(i,j)*( G%dt(i-1,j)/G%dt(i,j)*S%F12(i+1,j) + G%dtdt(i,j)*S%F12(i,j) - G%dt(i,j)/G%dt(i-1,j)*S%F12(i-1,j) )

            S%data(count) = -G%ft(i,j)*( (S%F11(i,j)+S%F11(i+1,j))/G%dt(i,j)+(S%F11(i,j)+S%F11(i-1,j))/G%dt(i-1,j) ) - &
                 G%fp(i,j)/sin(G%t(i,j))**2*( (S%F22(i,j)+S%F22(i,jp1))/G%dp(i,j)+(S%F22(i,j)+S%F22(i,jm1))/G%dp(i,jm1) ) + &
                 dF12t*G%fp(i,j)*G%dpdp(i,j)-&
                 dF12p*G%ft(i,j)*G%dtdt(i,j) 
            S%II(count) = K(i,j,G%Nt)
            S%JJ(count) = K(i,j,G%Nt)
            count=count+1

            S%data(count) = G%ft(i,j)*(S%F11(i,j)+S%F11(i+1,j))/G%dt(i,j)-&
                 dF12p*G%ft(i,j)*G%dt(i-1,j)/G%dt(i,j)
            S%II(count)  = K(i,j,G%Nt)
            S%JJ(count)  = K(i+1,j,G%Nt)
            count=count+1

            S%data(count) = G%ft(i,j)*(S%F11(i,j)+S%F11(i-1,j))/G%dt(i-1,j)+&
                 dF12p*G%ft(i,j)*G%dt(i,j)/G%dt(i-1,j)
            S%II(count) = K(i,j,G%Nt)
            S%JJ(count)  = K(i-1,j,G%Nt)
            count=count+1

            S%data(count) = G%fp(i,j)/sin(G%t(i,j))**2*(S%F22(i,j)+S%F22(i,jp1))/G%dp(i,j)+&
                 dF12t*G%fp(i,j)*G%dp(i,jm1)/G%dp(i,j)
            S%II(count) = K(i,j,G%Nt)
            S%JJ(count) = K(i,jp1,G%Nt)
            count=count+1

            S%data(count) = G%fp(i,j)/sin(G%t(i,j))**2*(S%F22(i,j)+S%F22(i,jm1))/G%dp(i,jm1)-&
                 dF12t*G%fp(i,j)*G%dp(i,j)/G%dp(i,jm1)
            S%II(count) = K(i,j,G%Nt)
            S%JJ(count) = K(i,jm1,G%Nt)
            count=count+1

            S%RHS(K(i,j,G%Nt)) = St%Vars(i,j,FAC)
         enddo

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! low lat boundary
         S%data(count)  = 1.0_mix_real
         S%II(count)    = K(G%Nt,j,G%Nt)
         S%JJ(count)    = K(G%Nt,j,G%Nt)
         count=count+1

         S%RHS(K(G%Nt,j,G%Nt)) = LLBC(j)
         ! end low lat boundary
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      
    end subroutine set_solver_matrix_and_rhs

end module mixsolver
