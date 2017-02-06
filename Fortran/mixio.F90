! Input/output for mix

module mixio
  use hdf5
  use mixdefs  

  implicit none

  !Necessary for HDF5 routines
  integer :: herror

  contains

    ! get UT: int array (Year, Month, Day, Hr, Min, Sec)
    subroutine getUT(fname,simtime)
      character(len=*),intent(in) :: fname
      integer, dimension(6), intent(out) :: simtime
      
      integer(HID_T) :: h5fId,attrId
#ifdef NO2003_HDF5
      integer(HSIZE_T), dimension(1) :: dims
#else
      TYPE(C_PTR) :: f_ptr
#endif

      call checkFile(fname)
      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      call h5aopen_f(h5fId,"SimTime",attrId,herror)
#ifdef NO2003_HDF5
      call h5aread_f(attrId,H5T_NATIVE_INTEGER,simtime,dims,herror)
#else
      f_ptr = C_LOC(simtime)
      call h5aread_f(attrId,H5T_NATIVE_INTEGER,f_ptr,herror)
#endif
      
      call h5aclose_f(attrId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine getUT

    subroutine getVar(fname,varname,var)
      character(len=*),intent(in) :: fname
      character(len=*),intent(in) :: varname
      
      integer(HID_T) :: h5fId,dsId,dspaceId
      integer(HSIZE_T), dimension(2) :: dims,maxdims

#ifndef NO2003_HDF5
      TYPE(C_PTR) :: f_ptr
#endif

      real, dimension(:,:), allocatable :: var
!      real, dimension(27,180) :: var

      call checkFile(fname)
      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      call h5dopen_f(h5fId,varname,dsId,herror)
      call h5dget_space_f(dsId,dspaceId,herror);
      call h5sget_simple_extent_dims_f(dspaceId, dims, maxdims, herror) 

      if (.not. allocated(var)) then 
         allocate(var(1:dims(1),1:dims(2)))
      end if

#ifdef NO2003_HDF5
      call h5dread_f(dsId,H5T_NATIVE_REAL,var,dims,herror)
#else
      f_ptr = C_LOC(var)
      call h5dread_f(dsId,H5T_NATIVE_REAL,f_ptr,herror)
#endif
      
      call h5dclose_f(dsId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine getVar

    
    ! check H5 file existence
    subroutine checkFile(fname)
      character(len=*),intent(in) :: fname
      logical :: fExist

      inquire(file=fname,exist=fExist)
      if (.not.(fExist)) then
         write(*,"(a,a)") "Cannot read file ", fname
         stop
      endif
    end subroutine checkFile
      
end module mixio
