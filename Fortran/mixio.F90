! Input/output for mix

module mixio
  use hdf5
  use mixdefs  
  use mixtypes

  implicit none

  !Necessary for HDF5 routines
  integer :: herror

  contains

    ! get UT: int array (Year, Month, Day, Hr, Min, Sec)
    subroutine getUT(fname,simtime)
      character(len=*),intent(in) :: fname
      integer, dimension(6), target, intent(out) :: simtime
      
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

    subroutine readVar(fname,varname,var)
      character(len=*),intent(in) :: fname
      character(len=*),intent(in) :: varname
      
      integer(HID_T) :: h5fId,dsId,dspaceId
      integer(HSIZE_T), dimension(2) :: dims,maxdims

#ifndef NO2003_HDF5
      TYPE(C_PTR) :: f_ptr
#endif

      real, dimension(:,:), target, allocatable :: var

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
    end subroutine readVar

    subroutine writeState(fname,St)
      character(len=*),intent(in) :: fname
      type(State_T), intent(in) :: St
      integer(HID_T) :: h5fId
      integer :: v

      call h5fcreate_f(fname,H5F_ACC_EXCL_F, h5fId, herror) 
      call h5fopen_f(fname,H5F_ACC_RDWR_F, h5fId, herror)
      
      do v=1,nVars
         select case (v)
            case (POT)
               call writeVar(h5fId,St%Vars(:,:,v),"Potential","kV")
            case (FAC)
               call writeVar(h5fId,St%Vars(:,:,v),"Field-aligned current","muA/m**2")
            case (SIGMAP)
               call writeVar(h5fId,St%Vars(:,:,v),"Pedersen conductance","S")
            case (SIGMAH)
               call writeVar(h5fId,St%Vars(:,:,v),"Hall conductance","S")
         end select
      enddo

      call h5fclose_f(h5fId, herror)
    end subroutine writeState

    subroutine writeGrid(fname,G)
      character(len=*),intent(in) :: fname
      type(Grid_T), intent(in) :: G
      integer(HID_T) :: h5fId

      call h5fopen_f(fname,H5F_ACC_RDWR_F, h5fId, herror)
      call writeVar(h5fId,G%x,'X',"Ri")
      call writeVar(h5fId,G%y,'Y',"Ri")
      call h5fclose_f(h5fId, herror)
    end subroutine writeGrid

    subroutine writeVar(fId,var,varName,units)
      integer(HID_T), intent(in) :: fId
      real(mix_real), dimension(:,:), intent(in) :: var
      character(len=*),intent(in) :: varName
      character(len=*),intent(in) :: units
      
      integer(HID_T) :: dsId,dspaceId,aspaceId,attrId,atypeId
      integer(HSIZE_T), dimension(2) :: dims
      integer(HSIZE_T), dimension(1) :: adims = [1]
      integer(SIZE_T) :: attrlen=100
      integer :: rank = 2, arank = 1
      
      dims = shape(var)
      ! create data space
      call h5screate_simple_f(rank, dims,dspaceId,herror)
      !Create data set and write data
      call h5dcreate_f(fId,varName,H5T_NATIVE_REAL, dspaceId,dsId, herror)

      ! data set attributes
      call h5screate_simple_f(arank, adims, aspaceId,herror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER, atypeId, herror)
      call h5tset_size_f(atypeId,attrlen, herror)
      call h5acreate_f (dsId,"Units",atypeId,aspaceId,attrId,herror)
      call h5awrite_f (attrId,atypeId,trim(units),adims,herror) 
      call h5aclose_f (attrId,herror)

      ! note type conversion to io_real (single, typically)
      call h5dwrite_f(dsId, H5T_NATIVE_REAL,real(var,mix_io_real),dims, herror)

      !Close up shop
      call h5dclose_f(dsId,herror)
      call h5sclose_f(dspaceId,herror)
    end subroutine writeVar

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
