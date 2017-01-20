program MIX
  use mixdefs
  use mixio
  use mixtypes
  use mixgeom
  use mixsolver
  
  implicit none

  character(len=strLen) :: fname
  integer, dimension(6) :: simtime
  integer, dimension(2) :: dims
  integer :: Nt,Np
  integer :: u
  real(mix_io_real), dimension(:,:), allocatable :: x_in,y_in,&
      pot_in,j_in,sigmap_in,sigmah_in
  type(Grid_T) :: Grid
  type(State_T) :: State
  type(Params_T) :: Params  ! to be read from xml file, placeholder for now
  type(Solver_T) :: Solver
  real(mix_real),dimension(:),allocatable :: LLBC ! low latitude boundary condition


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TEMPORARY IO STUFF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call h5open_f(herror) !Setup H5 Fortran interface

  fname = "../data/Aug2010_mix_2010-08-04T00-00-00Z.h5"

  ! FIXME: pack everything into one 3D array eventually; also define grid class
  call getUT(fname,simtime)
  call getVar(fname,"Grid X",x_in)
  call getVar(fname,"Grid Y",y_in)
  call getVar(fname,"Potential North [kV]",pot_in)
  call getVar(fname,"FAC North [microAm2]",j_in)
  call getVar(fname,"Pedersen conductance North [S]",sigmap_in)
  call getVar(fname,"Hall conductance North [S]",sigmah_in)
  ! note, things are stored in the H5 file as (Nt,Np),
  ! while we want (Np,Nt); Note this for later.
  dims = shape(x_in); Nt = dims(1); Np = dims(2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initiate and set grid variables
  ! note, casting the arrays into whatever was defined as mix_real type (double by default)
  call init_grid_fromXY(Grid,real(x_in,mix_real),real(y_in,mix_real))
  call set_grid(Grid)

  ! FIXME: pack these into a module, like we did for grid above in mixgeom
  ! set up variable state
  allocate(State%Vars(Grid%Np,Grid%Nt,nVars))
  State%Vars(:,:,POT) = transpose(pot_in)
  State%Vars(:,:,FAC) = transpose(j_in)
  State%Vars(:,:,SIGMAP) = transpose(sigmap_in)
  State%Vars(:,:,SIGMAH) = transpose(sigmah_in)

  ! ok, done with the setup of the main data structures.
  ! now can do the solve
  call init_solver(Params,Grid,State,Solver)

  ! FIXME: low latitude boundary condition
  allocate(LLBC(Grid%Np))
  LLBC = 0.0_mix_real

  ! MAIN LOOP WILL START HERE
  call set_solver_terms(Params,Grid,State,Solver)
  call set_solver_matrix_and_rhs(Params,Grid,State,Solver,LLBC)

  open(newunit=u, file="data.dat", status="replace")
  write(u, *) Solver%data
  close(u)
  open(newunit=u, file="II.dat", status="replace")  
  write(u, *) Solver%II
  close(u)
  open(newunit=u, file="JJ.dat", status="replace")
  write(u, *) Solver%JJ
  close(u)
  open(newunit=u, file="RHS.dat", status="replace")  
  write(u, *) Solver%RHS
  close(u)

  open(newunit=u, file="x.dat", status="replace")  
  write(u, *) Grid%x
  close(u)
  open(newunit=u, file="y.dat", status="replace")  
  write(u, *) Grid%y
  close(u)

  call h5close_f(herror)  ! Close H5 Fortran interface
end program MIX
