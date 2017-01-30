#ifdef pardiso_solver
  include 'mkl_pardiso.f90'
#elif dss_solver
  include 'mkl_dss.f90'
#endif
program MIX
  use mixdefs
  use mixio
  use mixtypes
  use mixgeom
  use mixsolver
#ifdef pardiso_solver
  use mkl_pardiso
#elif dss_solver
  use mkl_dss
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

#ifdef dss_solver
  TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
  INTEGER :: error
!  integer, dimension(:), allocatable :: perm
  INTEGER perm(1)
#elif pardiso_solver
!  INTEGER ::  pt(64)
  TYPE(MKL_PARDISO_HANDLE) :: pt(64) 
  INTEGER :: error
  integer :: mtype,iparm(64)  ! pardisoinit parameters
  integer, dimension(:), allocatable :: perm
  integer :: maxfct,mnum,phase ! pardiso parameters

#endif


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
  call set_solver_matrix_and_rhs(Params,Grid,State,Solver,LLBC)


#ifdef dss_solver
  ! INITIAL STUFF FOR SOLVER
  ! Initialize the solver.
  error = DSS_CREATE( handle, MKL_DSS_DEFAULTS )
  if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error

  ! Define the non-zero structure of the matrix.
   error = DSS_DEFINE_STRUCTURE( handle, MKL_DSS_NON_SYMMETRIC,int(Solver%rowI),int(Grid%Np*Grid%Nt), &
        int(Grid%Np*Grid%Nt),int(Solver%JJ),int(Solver%nnz))
  if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error

  ! Reorder the matrix.
  !  allocate(perm(Grid%Np*Grid%Nt))
  perm(1) = 0
  error = DSS_REORDER( handle, MKL_DSS_AUTO_ORDER, perm )
  if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error

  ! Deallocate solver storage and various local arrays.
  error = DSS_DELETE( handle, MKL_DSS_DEFAULTS )
  if (error /=MKL_DSS_SUCCESS) WRITE(*,*) "Solver returned error code ", error
  ! INITIAL STUFF FOR SOLVER
#elif pardiso_solver
  ! see description of parameters here: https://software.intel.com/en-us/node/470284#E44B4021-701A-48DA-BA29-70CFA20766AA
  mtype = 11 ! real nonsymmetric matrix

  call pardisoinit (pt, mtype, iparm)

  allocate(perm(Grid%Np*Grid%Nt))
  allocate(solution(Grid%Np*Grid%Nt))
  iparm(1) = 0;! iparm(27)=1 !iparm(28) =0; iparm (6) =0;
  maxfct = 1; mnum =1; 
  ! phase =11
  ! call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,1,iparm,1,Solver%RHS,solution,error)
  ! phase =33
  ! call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,1,iparm,1,Solver%RHS,solution,error)
  phase =13
  call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,1,iparm,1,Solver%RHS,solution,error)
  phase =-1  ! release
  call pardiso(pt,maxfct,mnum,mtype,phase,int(Grid%Np*Grid%Nt),Solver%data,int(Solver%rowI),int(Solver%JJ),perm,1,iparm,0,Solver%RHS,solution,error)

  
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

  ! open(newunit=u, file="solution.dat", status="replace")  
  ! write(u, *) solution
  ! close(u)

  call h5close_f(herror)  ! Close H5 Fortran interface
end program MIX
