program legacy_abort_probe
   use mpi
   use sparse, only: sparse_matrix, initialize_sparse, clear_matrix, add
   use mumps_solver, only: SolveOneDirection
   use mumps_test_support, only: configure_mumps_stub
   implicit none

   type(sparse_matrix), pointer :: matrix
   real(kind=8) :: rhs(2, 1)
   integer(kind=4) :: rank, ierr

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   if (rank == 0) then
      call initialize_sparse(2, 2, matrix)
      call add(matrix, 0, 0, 1.d0)
      call add(matrix, 1, 1, 1.d0)
      rhs = 1.d0
      call configure_mumps_stub(1, -77)
      call SolveOneDirection(rhs, 1, 1, 1, matrix)
      call clear_matrix(matrix)
   else
      ! This rank represents callers which would otherwise wait forever
      ! after a rank-local STOP in the legacy interface.
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
   end if

   call MPI_Finalize(ierr)
   stop 99
end program legacy_abort_probe
