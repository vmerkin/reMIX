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

      print *,G%dt
    end subroutine solver_init
end module mixsolver
