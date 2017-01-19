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

      allocate(S%RHS(G%Nt*G%Np))
      ! init to zero. This is important, since we're only filling in
      ! non-zero elements in the matrix construction
      S%RHS = 0.0_mix_real   

      ! this is the number of non-zeros in the matrix. Note, this
      ! depends on the stencil and boundary conditions
      S%nnz = G%Np*(G%Nt-2)*5 + G%Np + G%Np*(G%Np+1) 

      ! Matrix 
      allocate(S%data(S%nnz))
      allocate(S%II(S%nnz))
      allocate(S%JJ(S%nnz))
      S%data=0.0_mix_real
      S%II = 0.0_mix_real
      S%JJ = 0.0_mix_real

      ! unlike the Python solver, treat the periodic boundary separately
      ! to avoid introducing ghost cells (
      do j=1,G%Np
         do i=1,G%Nt
            print *,G%cosd(i,j)
         enddo
      enddo
    end subroutine init_solver

end module mixsolver
