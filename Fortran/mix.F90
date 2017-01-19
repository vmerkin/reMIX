include 'mkl_pardiso.f90'

program MIX
  use mkl_pardiso
  
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
  real(mix_io_real), dimension(:,:), allocatable :: x_in,y_in,&
      pot_in,j_in,sigmap_in,sigmah_in
  type(Grid_T) :: Grid
  type(State_T) :: State
  type(Params_T) :: Params  ! to be read from xml file, placeholder for now
  type(Solver_T) :: Solver

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



  ! initiate and set grid variables
  ! note, casting the arrays into whatever was defined as mix_real type (double by default)
  call init_grid_fromXY(Grid,real(x_in,mix_real),real(y_in,mix_real))
  call set_grid(Grid)

  ! FIXME: pack these into a module, like we did for grid above in mixgeom
  ! set up variable state
  dims = shape(x_in); Nt = dims(1); Np = dims(2)
  allocate(State%Vars(Grid%Nt,Grid%Np,nVars))
  State%Vars(:,:,POT) = pot_in
  State%Vars(:,:,FAC) = j_in
  State%Vars(:,:,SIGMAP) = sigmap_in
  State%Vars(:,:,SIGMAH) = sigmah_in

  ! ok, done with the setup of the main data structures.
  ! now can do the solve
  call init_solver(Params,Grid,State,Solver)

  call h5close_f(herror)  ! Close H5 Fortran interface
end program MIX
