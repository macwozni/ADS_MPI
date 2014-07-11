
      module utils

      include "mpif.h"
      contains


c Calculates the range of the direction that is assigned to processor
c with specified rank.

c rank     rank of the current process
c nrproc   # of processors for this direction
c n        size of the problem
c nrcpp    # of columns per processor
c ibeg     index of first assigned slice
c iend     index of last assigned slice

      subroutine ComputeEndpoints(rank,nrproc,n,nrcpp,ibeg,iend)
      implicit none
      integer :: rank, nrproc, n
      integer, intent(out) :: nrcpp, ibeg, iend

      nrcpp = (n+1+1) / nrproc
      if(rank .eq. nrproc-1)then
        ibeg = nrcpp*(nrproc-1)+1
        iend = n+1
      else
        ibeg = nrcpp*(nrproc-2)+1
        iend = nrcpp*(nrproc-1)
      endif

      end subroutine ComputeEndpoints


c Calculates sizes and ranges of slices for each processor in
c given direction.
c
c dims      sizes of 'slices' of dimension
c shifts    offsets of slices
c nrcpp     # of columns per processor
c stride    size of full, 3d slice
c n         size of the problem
c nrproc    # of processors for this direction

      subroutine FillDimVector(dims,shifts,nrcpp,stride,n,nrproc)
      implicit none
      integer, allocatable :: dims(:), shifts(:)
      integer :: nrcpp, stride, n, nrproc
      integer:: i

      allocate(dims(nrproc))
      allocate(shifts(nrproc))

      shifts = 0
      dims = 0

      do i = 1,nrproc-1
        dims(i) = nrcpp * stride
        if (i .gt. 1) shifts(i) = shifts(i-1) + dims(i-1)
      enddo

      dims(nrproc) = ((n+1) - nrcpp * (nrproc-1)) * stride
      shifts(nrproc) = shifts(nrproc-1) + dims(nrproc-1)

      end subroutine FillDimVector


c Linearizes and transposes an array
c
c F       input rank-2 array
c F_lin   output rank-1 array
c elems   first dimension of F
c stride  second dimension of F

      subroutine Linearize(F,F_lin,elems,stride)
      implicit none
      integer :: elems, stride
      real (kind=8), intent(in)  :: F(elems, stride)
      real (kind=8), intent(out) :: F_lin(elems * stride)
      integer :: i, a, b

      do i = 1,elems
        a = (i-1) * stride + 1
        b = i * stride
        F_lin(a:b) = F(i,:)
      enddo

      end subroutine Linearize


c Delinearizes an array
c
c F_lin   input rank-1 array
c F       output rank-2 array
c elems   first dimension of F
c stride  second dimension of F

      subroutine Delinearize(F_lin,F,elems,stride)
      implicit none
      integer :: elems, stride
      real (kind=8), intent(in)  :: F_lin(elems * stride)
      real (kind=8), intent(out) :: F(elems, stride)
      integer :: i, a, b

      do i = 1,elems
        a = (i-1) * stride + 1
        b = i * stride
        F(i,:) = F_lin(a:b) 
      enddo

      end subroutine Delinearize


c Gathers data along one axis
c
c F       input array
c F_out   output array
c n       problem size (total length of the axis)
c elems   length of slice assigned to this process
c stride  total size of each slice layer
c dims    sizes of slices for all processors
c shifts  offsets (linearized) of slices for all processors
c comm    communicator of the axis
c ierr    error code output

      subroutine Gather(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
      implicit none
      integer :: n, elems, stride, comm
      real (kind=8), intent(in)  :: F(elems, stride)
      real (kind=8), intent(out) :: F_out(n+1, stride)
      integer :: dims(:), shifts(:)
      integer, intent(out) :: ierr
      real (kind=8) :: F_lin(elems * stride), F_out_lin((n+1) * stride)

      call Linearize(F,F_lin,elems,stride)

      call mpi_gatherv(F_lin,
     .  elems * stride,
     .  MPI_DOUBLE_PRECISION,
     .  F_out_lin,
     .  dims, shifts,
     .  MPI_DOUBLE_PRECISION,
     .  0, comm, ierr) 

      call Delinearize(F_out_lin,F_out,n+1,stride)

      end subroutine Gather


c Scatters computed partial solution along one axis, but does not
c delinearize the output (used at the very end of computation)
c
c F       data to scatter
c F_out   buffer to receive data
c n       problem size
c elems   length of received slice
c stride  total size of each slice layer
c dims    sizes of slices for all processors
c shifts  offsets (linearized) of slices for all processors
c comm    communicator of the axis
c ierr    error code output

      subroutine Scatter2(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
      implicit none
      integer :: n, elems, stride, comm
      real (kind=8), intent(in) :: F(n+1, stride)
      real (kind=8), intent(out) :: F_out(elems * stride)
      integer :: dims(:), shifts(:)
      integer, intent(out) :: ierr
      real (kind=8) :: F_lin((n+1) * stride)

      call Linearize(F, F_lin, n+1, stride)

      call mpi_scatterv(F_lin,
     .  dims, shifts,
     .  MPI_DOUBLE_PRECISION,
     .  F_out,
     .  elems * stride,
     .  MPI_DOUBLE_PRECISION,
     .  0, comm, ierr)

      end subroutine Scatter2


c Scatters computed partial solution along one axis
c
c F       data to scatter
c F_out   buffer to receive data
c n       problem size
c elems   length of received slice
c stride  total size of each slice layer
c dims    sizes of slices for all processors
c shifts  offsets (linearized) of slices for all processors
c comm    communicator of the axis
c ierr    error code output

      subroutine Scatter(F,F_out,n,elems,stride,dims,shifts,comm,ierr)
      implicit none
      integer :: n, elems, stride, comm
      real (kind=8), intent(in) :: F(n+1, stride)
      real (kind=8), intent(out) :: F_out(elems, stride)
      integer :: dims(:), shifts(:)
      integer, intent(out) :: ierr
      real (kind=8) :: F_out_lin(elems * stride)

      call Scatter2(F,F_out_lin,n,elems,stride,dims,shifts,comm,ierr)
      call Delinearize(F_out_lin, F_out, elems, stride)

      end subroutine Scatter


c Broadcasts computed partial solution along one axisutation)
c
c F       data to distribute
c F_out   buffer to receive data
c n       problem size
c elems   length of received slice
c stride  total size of each slice layer
c dims    sizes of slices for all processors
c shifts  offsets (linearized) of slices for all processors
c comm    communicator of the axis
c ierr    error code output

      subroutine AllGather(F,F_out,n,elems,stride,dims,shifts,comm)
      implicit none
      integer :: n, elems, stride, comm
      real (kind=8), intent(in) :: F(elems, stride)
      real (kind=8), intent(out) :: F_out(n+1, stride)
      integer :: dims(:), shifts(:)
      real (kind=8) :: F_lin(elems*stride), F_out_lin((n+1)*stride)
      integer :: ierr

      call Linearize(F,F_lin,elems,stride)

      call mpi_allgatherv(F_lin,
     .  elems * stride,
     .  MPI_DOUBLE_PRECISION,
     .  F_out_lin,
     .  dims, shifts,
     .  MPI_DOUBLE_PRECISION,
     .  comm, ierr) 

      call Delinearize(F_out_lin,F_out,n+1,stride)

      end subroutine AllGather


      end module

