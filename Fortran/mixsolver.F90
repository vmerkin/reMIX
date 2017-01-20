! sets up and solves the stencil matrix
#ifdef pardiso
  include 'mkl_pardiso.f90'
#elif dss
  include 'mkl_dss.f90'
#endif

module mixsolver
  use mixdefs
  use mixtypes
#ifdef pardiso
  use mkl_pardiso
#elif dss
  use mkl_dss
#endif

  implicit none

  contains
    ! running index
    integer(kind=8) function K(j,i,Nj)
      integer, intent(in) :: i,j,Nj
      K=(i-1)*Nj+j
    end function K

    subroutine init_solver(P,G,St,S)
      type(State_T), intent(in) :: St
      type(Params_T), intent(in) :: P
      type(Grid_T), intent(in) :: G
      type(Solver_T), intent(inout) :: S
      integer :: i,j

      TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
      INTEGER :: error

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! INITIAL STUFF FOR SOLVER
      ! ! Initialize the solver.
      ! error = DSS_CREATE( handle, MKL_DSS_DEFAULTS )
      ! if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error

      ! ! Define the non-zero structure of the matrix.
      ! error = DSS_DEFINE_STRUCTURE( handle, MKL_DSS_NON_SYMMETRIC, rowIndex, nRows, &
      !      & nCols, columns, nNonZeros )
      ! IF (error /= MKL_DSS_SUCCESS) GOTO 999


      ! ! Deallocate solver storage and various local arrays.
      ! error = DSS_DELETE( handle, MKL_DSS_DEFAULTS )
      ! if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error
      ! INITIAL STUFF FOR SOLVER
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! this is the number of non-zeros in the matrix. Note, this
      ! depends on the stencil and boundary conditions
      S%nnz = G%Np*(G%Nt-2)*5 + G%Np + G%Np*(G%Np+1) 

      ! allocate RHS vector
      allocate(S%RHS(G%Nt*G%Np))
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Pole boundary 
      do j=1,G%Np
         S%data(count) = 1.0_mix_real
         S%II(count)   = K(j,1,G%Np)
         S%JJ(count)   = K(j,1,G%Np)
         count=count+1

         do jj=1,G%Np   ! with (Np,Nt) definition of the grid, this is a natural alignment
            S%data(count) = -G%dp(jj,2)/(2*mix_pi)
            S%II(count) = K(j,1,G%Np)
            S%JJ(count) = K(jj,2,G%Np)
            count=count+1
         enddo
      enddo
      ! note, not setting RHS because it's initializaed to zero anyway
      ! end pole boundary
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! inner block
      do i=2,G%Nt-1  ! excluding pole and low lat boundaries
         do j=1,G%Np
!            jm1 = merge(G%Np,j-1,j.eq.1)   ! maps j-1=0 to j-1=Np, otherwise returns j-1
!            jp1 = merge(1,j+1,j.eq.G%Np)   ! similar for j+1

            ! the above functions involve if statements, which I don't
            ! want to pack inside inner loop. Do this instead: these
            ! wonderful functions of my invention maps j-1=0 to Np,
            ! otherwise j-1 maps to itself. Similarly, j+1=Np is
            ! mapped to j+1=1, otherwise nothing's done
            jm1 = modulo(j-1,G%Np)+G%Np*(1-int(ceiling(real(j-1)/G%Np)))  
            jp1 = modulo(j+1,G%Np)+G%Np*(1-int(ceiling(real(modulo(j+1,G%Np))/G%Np)))

            ! derivatives for off diagonal conductance terms
            dF12p = G%fp(j,i)*( G%dp(jm1,i)/G%dp(j,i)*S%F12(jp1,i) + G%dpdp(j,i)*S%F12(j,i) - G%dp(j,i)/G%dp(jm1,i)*S%F12(jm1,i) )
            dF12t = G%ft(j,i)*( G%dt(j,i-1)/G%dt(j,i)*S%F12(j,i+1) + G%dtdt(j,i)*S%F12(j,i) - G%dt(j,i)/G%dt(j,i-1)*S%F12(j,i-1) )

            S%data(count) = -G%ft(j,i)*( (S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)+(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1) ) - &
                 G%fp(j,i)/sin(G%t(j,i))**2*( (S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i) ) + &
                 dF12t*G%fp(j,i)*G%dpdp(j,i)-&
                 dF12p*G%ft(j,i)*G%dtdt(j,i) 
            S%II(count) = K(j,i,G%Np)
            S%JJ(count) = K(j,i,G%Np)
            count=count+1

            S%data(count) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i+1))/G%dt(j,i)-&
                 dF12p*G%ft(j,i)*G%dt(j,i-1)/G%dt(j,i)
            S%II(count)  = K(j,i,G%Np)
            S%JJ(count)  = K(j,i+1,G%Np)
            count=count+1

            S%data(count) = G%ft(j,i)*(S%F11(j,i)+S%F11(j,i-1))/G%dt(j,i-1)+&
                 dF12p*G%ft(j,i)*G%dt(j,i)/G%dt(j,i-1)
            S%II(count) = K(j,i,G%Np)
            S%JJ(count)  = K(j,i-1,G%Np)
            count=count+1

            S%data(count) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jp1,i))/G%dp(j,i)+&
                 dF12t*G%fp(j,i)*G%dp(jm1,i)/G%dp(j,i)
            S%II(count) = K(j,i,G%Np)
            S%JJ(count) = K(jp1,i,G%Np)
            count=count+1

            S%data(count) = G%fp(j,i)/sin(G%t(j,i))**2*(S%F22(j,i)+S%F22(jm1,i))/G%dp(jm1,i)-&
                 dF12t*G%fp(j,i)*G%dp(j,i)/G%dp(jm1,i)
            S%II(count) = K(j,i,G%Np)
            S%JJ(count) = K(jm1,i,G%Np)
            count=count+1

            S%RHS(K(j,i,G%Np)) = St%Vars(j,i,FAC)
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! low lat boundary
      do j=1,G%Np
         S%data(count)  = 1.0_mix_real
         S%II(count)    = K(j,G%Nt,G%Np)
         S%JJ(count)    = K(j,G%Nt,G%Np)
         count=count+1

         S%RHS(K(j,G%Nt,G%Np)) = LLBC(j)
      enddo
    end subroutine set_solver_matrix_and_rhs
    ! end low lat boundary
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mixsolver
