#ifdef pardiso_solver
  include 'mkl_pardiso.f90'
#endif
program MIX
  use mixdefs
  use mixio
  use mixtypes
  use mixgeom
  use mixsolver
#ifdef pardiso_solver
  use mkl_pardiso
#endif
  
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
  real(mix_real),dimension(:),allocatable :: solution

#ifdef pardiso_solver
!  INTEGER ::  pt(64)
  TYPE(MKL_PARDISO_HANDLE) :: pt(64) 
  INTEGER :: error
  integer :: nrhs ! number of RHSs
  integer :: msglvl ! verbosity level
  integer :: mtype,iparm(64)  ! pardisoinit parameters
  integer, dimension(:), allocatable :: perm
  integer :: maxfct,mnum,phase ! pardiso parameters
#elif mgmres_solver
  integer :: maxitr
  integer :: mr
#endif

  ! timing 
  real :: start_time, finish_time


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TEMPORARY IO STUFF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call h5open_f(herror) !Setup H5 Fortran interface

!  fname = "../data/Aug2010_mix_2010-08-04T00-00-00Z.h5"
  fname = '../data/interp.h5'

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
  dims = shape(x_in); Nt = dims(2); Np = dims(1)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initiate and set grid variables
  ! note, casting the arrays into whatever was defined as mix_real type (double by default)
  call init_grid_fromXY(Grid,real(x_in,mix_real),real(y_in,mix_real))
  call set_grid(Grid)

  ! FIXME: pack these into a module, like we did for grid above in mixgeom
  ! set up variable state
  allocate(State%Vars(Grid%Np,Grid%Nt,nVars))
  State%Vars(:,:,POT) = pot_in
  State%Vars(:,:,FAC) = j_in
  State%Vars(:,:,SIGMAP) = sigmap_in
  State%Vars(:,:,SIGMAH) = sigmah_in

  ! ok, done with the setup of the main data structures.
  ! now can do the solve
  call init_solver(Params,Grid,State,Solver)

  ! FIXME: low latitude boundary condition
  allocate(LLBC(Grid%Np))
  LLBC = 0.0_mix_real

  ! MAIN LOOP WILL START HERE
  call set_solver_terms(Params,Grid,State,Solver)
  call cpu_time(start_time)
  call set_solver_matrix_and_rhs(Params,Grid,State,Solver,LLBC)
  call cpu_time(finish_time)
  print '("Time in matrix = ",f6.3," seconds.")',finish_time-start_time


#ifdef pardiso_solver
  ! see description of parameters here: https://software.intel.com/en-us/node/470284#E44B4021-701A-48DA-BA29-70CFA20766AA
  mtype = 11 ! real nonsymmetric matrix

  call pardisoinit (pt, mtype, iparm)

  allocate(perm(Grid%Np*Grid%Nt))
  allocate(solution(Grid%Np*Grid%Nt))
  iparm(1) = 0  ! all default parameters
  iparm(27)= 0  ! Matrix checker
  !iparm(28) =0; iparm (6) =0;
  maxfct = 1; mnum =1; 
  nrhs = 1
  msglvl = 0  ! no verbosity
  ! phase =11
  ! call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,nrhs,iparm,msglvl,Solver%RHS,solution,error)
  ! phase =33
  ! call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,nrhs,iparm,msglvl,Solver%RHS,solution,error)
  phase =13
  call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,nrhs,iparm,msglvl,Solver%RHS,solution,error)
  phase =-1  ! release
  call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,nrhs,iparm,0,Solver%RHS,solution,error)

#elif mgmres_solver  
  allocate(solution(Grid%Np*Grid%Nt))
  solution = 0.0_mix_real
  maxitr = 400
  mr = 30
!  call mgmres_st ( int(Grid%Np*Grid%Nt),int(Solver%nnz),int(Solver%II),int(Solver%JJ),Solver%data,solution,Solver%RHS,maxitr,mr,1.0D-3,1.0D-3)
  call cpu_time(start_time)
  call pmgmres_ilu_cr ( int(Grid%Np*Grid%Nt),int(Solver%nnz),int(Solver%rowI),int(Solver%JJ),Solver%data,solution,Solver%RHS,maxitr,mr,1.0D-8,1.0D-8)
  call cpu_time(finish_time)
  print '("Time in solve = ",f6.3," seconds.")',finish_time-start_time
#endif

  ! open(newunit=u, file="data.dat", status="replace")
  ! write(u, *) Solver%data
  ! close(u)
  ! open(newunit=u, file="II.dat", status="replace")  
  ! write(u, *) Solver%II
  ! close(u)
  ! open(newunit=u, file="JJ.dat", status="replace")
  ! write(u, *) Solver%JJ
  ! close(u)
  ! open(newunit=u, file="RHS.dat", status="replace")  
  ! write(u, *) Solver%RHS
  ! close(u)

  ! open(newunit=u, file="rhoind.dat", status="replace")  
  ! write(u, *) Solver%rowI
  ! close(u)

  open(newunit=u, file="solution.dat", status="replace")  
  write(u, *) solution
  close(u)

  call h5close_f(herror)  ! Close H5 Fortran interface
end program MIX
