! sets up and solves the stencil matrix

module mixsolver
  use mixdefs
  use mixtypes

  implicit none

  contains

    subroutine solver_init(P,G,S)
      type(State_T), intent(inout) :: S
      type(Params_T), intent(in) :: P
      type(Grid_T), intent(in) :: G

      integer :: i,j

      ! unlike the Python solver, treat the periodic boundary separately
      ! to avoid introducing ghost cells (
      do j=1,G%Np
         do i=2,G%Nt-1

         enddo
      enddo
    end subroutine solver_init
end module mixsolver
