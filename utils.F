      module utils

      contains

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

      end module
