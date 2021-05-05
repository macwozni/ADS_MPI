module DEBUGG


    contains
    

! -------------------------------------------------------------------
! Prints debugging information about results of distributing
! data to neighbouring processes.
! -------------------------------------------------------------------
subroutine PrintDistributedData(ads, ads_data)
    use Setup, ONLY: ADS_Setup, ADS_compute_data
    use parallelism, ONLY: MYRANKX, MYRANKY, MYRANKZ, PRINTRANK, &
    NRPROCX, NRPROCY, NRPROCZ, ComputeEndpoints
    implicit none
    type (ADS_setup), intent(inout) :: ads
    type (ADS_compute_data), intent(in) :: ads_data
    integer(kind = 4) :: i, j, k
    integer(kind = 4) :: obegx, oendx, obegy, oendy, obegz, oendz
    integer(kind = 4) :: mine, maxe

    write(*, *) PRINTRANK, 'R:'

    do i = max(MYRANKX - 1, 0) + 1, min(MYRANKX + 1, NRPROCX - 1) + 1
        do j = max(MYRANKY - 1, 0) + 1, min(MYRANKY + 1, NRPROCY - 1) + 1
            do k = max(MYRANKZ - 1, 0) + 1, min(MYRANKZ + 1, NRPROCZ - 1) + 1
            write(*, *) '(i,j,k)=', i + 1, j + 1, k + 1

            call ComputeEndpoints(i - 1, NRPROCX, ads % n(1), ads % p(1), ads % nrcpp(1), obegx, oendx, mine, maxe)
            call ComputeEndpoints(j - 1, NRPROCY, ads % n(2), ads % p(2), ads % nrcpp(2), obegy, oendy, mine, maxe)
            call ComputeEndpoints(k - 1, NRPROCZ, ads % n(3), ads % p(3), ads % nrcpp(3), obegz, oendz, mine, maxe)

            write(*, *) reshape(&
            ads_data % R(:, i - MYRANKX + 1, j - MYRANKY + 1, k - MYRANKZ + 1), &
            (/ oendz - obegz + 1, oendx - obegx + 1, oendy - obegy + 1 /))
            enddo
        enddo
    enddo

    write(*, *) '----'

end subroutine PrintDistributedData



! -------------------------------------------------------------------
! Sanity-check of dimensions vector
! -------------------------------------------------------------------
subroutine ValidateDimensions(n, s, nrcpp, &
    dimensionsX, dimensionsY, dimensionsZ)
    use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK
    use mpi
    implicit none
    integer(kind = 4), intent(in), dimension(3) :: n
    integer(kind = 4), intent(in), dimension(3) :: s
    integer(kind = 4), intent(in), dimension(3) :: nrcpp
    integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsX
    integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsY
    integer(kind = 4), intent(in), allocatable, dimension(:) :: dimensionsZ

    integer(kind = 4) :: i, k

    k = 0
    do i = 1, NRPROCX
        k = k + dimensionsX(i)
    enddo
    if (k /= (n(1) + 1) * s(2) * s(3)) then
        write(*, *) PRINTRANK, 'problem with dimensionsX', dimensionsX
        write(*, *) PRINTRANK, 'nx+1', n(1) + 1
        write(*, *) PRINTRANK, 'sy', s(2)
        write(*, *) PRINTRANK, 'sz', s(3)
        write(*, *) PRINTRANK, 'nrcppx', nrcpp(1)
        stop
    endif

    k = 0
    do i = 1, NRPROCY
        k = k + dimensionsY(i)
    enddo
    if (k /= (n(2) + 1) * s(1) * s(3)) then
        write(*, *) PRINTRANK, 'problem with dimensionsY', dimensionsY
        write(*, *) PRINTRANK, 'n+1', n(2) + 1
        write(*, *) PRINTRANK, 'sx', s(1)
        write(*, *) PRINTRANK, 'sz', s(3)
        stop
    endif

    k = 0
    do i = 1, NRPROCZ
        k = k + dimensionsZ(i)
    enddo
    if (k /= (n(3) + 1) * s(1) * s(2)) then
        write(*, *) PRINTRANK, 'problem with dimensionsZ', dimensionsZ
        write(*, *) PRINTRANK, 'n+1', n(3) + 1
        write(*, *) PRINTRANK, 'sx', s(1)
        write(*, *) PRINTRANK, 'sy', s(2)
        stop
    endif

end subroutine ValidateDimensions



! -------------------------------------------------------------------
! Displays computed domain decomposition, for debugging.
! -------------------------------------------------------------------
subroutine PrintDecompositionInfo(n, nrcpp, ibeg, iend)
        use parallelism, ONLY: NRPROCX, NRPROCY, NRPROCZ, PRINTRANK, &
        MYRANKX, MYRANKY, MYRANKZ
        implicit none
        integer(kind = 4), intent(in), dimension(3) :: n
        integer(kind = 4), intent(in), dimension(3) :: nrcpp
        integer(kind = 4), intent(in), dimension(3) :: ibeg
        integer(kind = 4), intent(in), dimension(3) :: iend

        write(*, *) PRINTRANK, 'MYRANKX,MYRANKY,MYRANKZ', MYRANKX, MYRANKY, MYRANKZ
        write(*, *) PRINTRANK, 'NRPROCX,NRPROCY,NRPROCZ', NRPROCX, NRPROCY, NRPROCZ
        write(*, *) PRINTRANK, 'nx+1', n(1) + 1
        write(*, *) PRINTRANK, 'ny+1', n(2) + 1
        write(*, *) PRINTRANK, 'nz+1', n(3) + 1
        write(*, *) PRINTRANK, 'nrcppx,nrcppy,nrcppz', nrcpp(1), nrcpp(2), nrcpp(3)
        write(*, *) PRINTRANK, 'ibegx,iendx', ibeg(1), iend(1)
        write(*, *) PRINTRANK, 'ibegy,iendy', ibeg(2), iend(2)
        write(*, *) PRINTRANK, 'ibegz,iendz', ibeg(3), iend(3)
end subroutine PrintDecompositionInfo

endmodule DEBUGG