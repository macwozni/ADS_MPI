module analysis

   implicit none

contains

   !!!!!!! udostepnic uzytkownikowi
   ! udostepnic wrapper
   ! utils?

   ! Ux,Uy,Uz  knots vectors
   ! px,py,py  orders
   ! nx,ny,nz  number of intervals (problems size is (nx+1)*(ny+1)*(nz+1)
   ! nelemx,nelemy,nelemz - number of elements
   ! F coefficients of the function
   subroutine NormL2(&
      Ux, Uy, Uz, &
      p, n, nelem, &
      ibeg, iend, nrank, nrp, &
      F)
      use basis, ONLY: BasisData
      use projection_engine, ONLY: global2local
      implicit none
      integer(kind = 4), dimension(3), intent(in) :: n, p, nelem
      integer(kind = 4), dimension(3), intent(in) :: ibeg
      integer(kind = 4), dimension(3), intent(in) :: iend
      integer(kind = 4), dimension(3), intent(in) :: nrank
      integer(kind = 4), dimension(3), intent(in) :: nrp
      real (kind = 8), dimension(0:n(1) + p(1) + 1), intent(in) :: Ux
      real (kind = 8), dimension(0:n(2) + p(2) + 1), intent(in) :: Uy
      real (kind = 8), dimension(0:n(3) + p(3) + 1), intent(in) :: Uz
      real (kind = 8), intent(out) :: F(0:(iend(1) - ibeg(1) + 1) - 1, &
      0:(iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1)
      integer(kind = 4) :: mx, my, mz, ngx, ngy, ngz, ex, ey, ez
      integer(kind = 4) :: kx, ky, kz, ax, ay, az, d
      integer(kind = 4), dimension(nelem(1)) :: Ox
      integer(kind = 4), dimension(nelem(2)) :: Oy
      integer(kind = 4), dimension(nelem(3)) :: Oz
      real (kind = 8), dimension(nelem(1)) :: Jx
      real (kind = 8), dimension(nelem(2)) :: Jy
      real (kind = 8), dimension(nelem(3)) :: Jz
      real (kind = 8), dimension(p(1) + 1) :: Wx
      real (kind = 8), dimension(p(2) + 1) :: Wy
      real (kind = 8), dimension(p(3) + 1) :: Wz
      real (kind = 8), dimension(p(1) + 1, nelem(1)) :: Xx
      real (kind = 8), dimension(p(2) + 1, nelem(2)) :: Xy
      real (kind = 8), dimension(p(3) + 1, nelem(3)) :: Xz
      real (kind = 8), dimension(0:0, 0:p(1), p(1) + 1, nelem(1)) :: NNx
      real (kind = 8), dimension(0:0, 0:p(2), p(2) + 1, nelem(2)) :: NNy
      real (kind = 8), dimension(0:0, 0:p(3), p(3) + 1, nelem(3)) :: NNz
      real (kind = 8) :: J, W, value
      integer(kind = 4) :: nreppx, nreppy, nreppz !# elements per proc along x,y,z
      integer(kind = 4) :: ind, ind1, ind23, ind23a, indx, indy, indz
      integer(kind = 4) :: iprint

      iprint = 0

      d = 0
      mx = n(1) + p(1) + 1
      ngx = p(1) + 1
      my = n(2) + p(2) + 1
      ngy = p(2) + 1
      mz = n(3) + p(3) + 1
      ngz = p(3) + 1

      call BasisData(p(1), mx, Ux, 0, ngx, nelem(1), Ox, Jx, Wx, Xx, NNx)
      call BasisData(p(2), my, Uy, 0, ngy, nelem(2), Oy, Jy, Wy, Xy, NNy)
      call BasisData(p(3), mz, Uz, 0, ngz, nelem(3), Oz, Jz, Wz, Xz, NNz)

      ! parallel number of elements per processors
      nreppx = nelem(1)/nrp(1)
      nreppy = nelem(2)/nrp(2)
      nreppz = nelem(3)/nrp(3)
      F = 0
      do ex = max(nreppx * nrank(1) - p(1) + 1, 1), min(nelem(1), nreppx * (nrank(1) + 1) + p(1))
         do ey = max(nreppy * nrank(2) - p(2) + 1, 1), min(nelem(2), nreppy * (nrank(2) + 1) + p(2))
            do ez = max(nreppz * nrank(3) - p(3) + 1, 1), min(nelem(3), nreppz * (nrank(3) + 1) + p(3))
               J = Jx(ex) * Jy(ey) * Jz(ez)
               do kx = 1, ngx
                  do ky = 1, ngy
                     do kz = 1, ngz
                        W = Wx(kx) * Wy(ky) * Wz(kz)
                        do ax = 0, p(1)
                           do ay = 0, p(2)
                              do az = 0, p(3)
                                 d = (Ox(ex) + ax)+(Oy(ey) + ay)*(n(1) + 1)+(Oz(ez) + az)*(n(2) + 1)*(n(1) + 1)
                                 call global2local(ind, [n(1), n(2), n(3)], indx, indy, indz)
                                 if (indx < ibeg(1) - 1 .or. indx > iend(1) - 1) cycle
                                 if (indy < ibeg(2) - 1 .or. indy > iend(2) - 1) cycle
                                 if (indz < ibeg(3) - 1 .or. indz > iend(3) - 1) cycle
                                 ind1 = indx - ibeg(1) + 1
                                 ind23 = (indy - ibeg(2) + 1)+(indz - ibeg(3) + 1)*(iend(2) - ibeg(2) + 1)

                                 if (ind1 < 0 .or. ind1 > (iend(1) - ibeg(1))) cycle
                                 if (ind23 < 0 .or. ind23 > (iend(2) - ibeg(2) + 1)*(iend(3) - ibeg(3) + 1) - 1) cycle

                                 ! parallel
                                 F(ind1, ind23) = F(ind1, ind23) + &
                                 NNx(0, ax, kx, ex) * NNy(0, ay, ky, ey) * NNz(0, az, kz, ez) * J * W*value
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

   end subroutine NormL2

end module analysis

