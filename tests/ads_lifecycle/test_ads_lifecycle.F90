program test_ads_lifecycle
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mpi
   use Setup, only: ADS_Setup, ADS_compute_data
   use parallelism, only: MYRANK, MYRANKX, MYRANKY, MYRANKZ, &
      NRPROC, NRPROCX, NRPROCY, NRPROCZ, InitializeParallelism, Cleanup_Parallelism
   use ads_lifecycle, only: initialize, initialize_setup, &
      ComputeDecomposition, AllocateADSdata, AllocateADS, &
      Cleanup_data, Cleanup_ADS
   implicit none

   real(kind=8), parameter :: TOL = 4096.d0*epsilon(1.d0)
   integer(kind=4) :: checks, failures, ierr
   integer(kind=4) :: process_grid(3)

   call read_process_grid(process_grid)
   call InitializeParallelism(process_grid(1), process_grid(2), &
                              process_grid(3), ierr)

   checks = 0
   failures = 0

   call test_allocate_ads(failures)
   call test_compute_decomposition(failures)
   call test_allocate_data_and_cleanup(failures)
   call test_initialize_setup(failures)
   call test_full_initialize_and_reinitialize(failures)

   if (MYRANK == 0) then
      if (failures == 0) then
         write (*, '(A,I0,A)') 'OK (', checks, ' ADS lifecycle checks)'
      else
         write (*, '(A,I0,A,I0,A)') 'FAILED (', failures, ' of ', checks, &
                                    ' ADS lifecycle checks)'
      end if
   end if

   call Cleanup_Parallelism(ierr)
   if (failures /= 0) stop 1

contains

   subroutine test_allocate_ads(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      type(ADS_Setup) :: ads
      integer(kind=4), parameter :: n(3) = (/2, 3, 4/)
      integer(kind=4), parameter :: nelem(3) = (/2, 3, 2/)
      integer(kind=4), parameter :: p(3) = (/1, 2, 1/)
      integer(kind=4), parameter :: ng(3) = (/2, 3, 4/)
      integer(kind=4) :: cleanup_ierr
      logical :: local_ok

      call AllocateADS(n, nelem, p, ng, ads)
      local_ok = allocate_ads_storage_matches(ads, n, nelem, p, ng)
      call assert_global('AllocateADS allocates exact anisotropic storage', &
                         local_ok, failure_count)

      call Cleanup_ADS(ads, cleanup_ierr)
      call assert_global('Cleanup_ADS releases a partial setup', &
                         cleanup_ierr == 0 .and. setup_is_clean(ads), &
                         failure_count)

      cleanup_ierr = -1
      call Cleanup_ADS(ads, cleanup_ierr)
      call assert_global('Cleanup_ADS is idempotent', &
                         cleanup_ierr == 0 .and. setup_is_clean(ads), &
                         failure_count)
   end subroutine test_allocate_ads


   subroutine test_compute_decomposition(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      type(ADS_Setup) :: ads
      integer(kind=4) :: cleanup_ierr

      ads%n = (/5, 7, 6/)
      ads%p = (/2, 1, 3/)
      call ComputeDecomposition(ads)

      call assert_global('ComputeDecomposition wires ownership, gathers, and neighbours', &
                         decomposition_matches(ads), failure_count)

      call Cleanup_ADS(ads, cleanup_ierr)
      call assert_global('Cleanup_ADS releases decomposition vectors', &
                         cleanup_ierr == 0 .and. setup_is_clean(ads), &
                         failure_count)
   end subroutine test_compute_decomposition


   subroutine test_allocate_data_and_cleanup(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      type(ADS_Setup) :: test_space, trial_space
      type(ADS_compute_data) :: data
      integer(kind=4), parameter :: nelem(3) = (/8, 7, 9/)
      integer(kind=4), parameter :: ptest(3) = (/3, 1, 2/)
      integer(kind=4), parameter :: ptrial(3) = (/1, 3, 1/)
      integer(kind=4), parameter :: test_ng(3) = (/4, 2, 3/)
      integer(kind=4), parameter :: trial_ng(3) = (/2, 4, 2/)
      integer(kind=4) :: cleanup_ierr, setup_ierr

      call initialize_setup(nelem + ptest - 1, ptest, (/0, 0, 0/), &
                            test_ng, test_space, setup_ierr)
      call initialize_setup(nelem + ptrial - 1, ptrial, (/0, 0, 0/), &
                            trial_ng, trial_space, setup_ierr)

      call AllocateADSdata(test_space, trial_space, data)
      call assert_global('AllocateADSdata uses union ranges and maximum quadrature', &
                         data_core_matches(data, &
                            min(test_space%mine, trial_space%mine), &
                            max(test_space%maxe, trial_space%maxe), &
                            max(test_ng, trial_ng)), failure_count)
      call assert_global('AllocateADSdata derives an exact reusable halo plan', &
                         halo_plan_matches(data, trial_space), failure_count)

      call allocate_optional_rhs(data)
      call Cleanup_data(data, cleanup_ierr)
      call assert_global('Cleanup_data releases every runtime buffer', &
                         cleanup_ierr == 0 .and. data_is_clean(data), &
                         failure_count)

      cleanup_ierr = -1
      call Cleanup_data(data, cleanup_ierr)
      call assert_global('Cleanup_data is idempotent', &
                         cleanup_ierr == 0 .and. data_is_clean(data), &
                         failure_count)

      call Cleanup_ADS(test_space, cleanup_ierr)
      call Cleanup_ADS(trial_space, cleanup_ierr)
   end subroutine test_allocate_data_and_cleanup


   subroutine test_initialize_setup(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      type(ADS_Setup) :: ads
      integer(kind=4), parameter :: n(3) = (/4, 5, 6/)
      integer(kind=4), parameter :: p(3) = (/1, 1, 1/)
      integer(kind=4), parameter :: ng(3) = (/1, 1, 1/)
      integer(kind=4), parameter :: nelem(3) = (/4, 5, 6/)
      integer(kind=4) :: cleanup_ierr, setup_ierr

      call initialize_setup(n, p, (/0, 0, 0/), ng, ads, setup_ierr)

      call assert_global('initialize_setup stores metadata and exact shapes', &
                         setup_ierr == 0 .and. &
                         initialized_storage_matches(ads, n, p, ng, nelem), &
                         failure_count)
      call assert_global('initialize_setup computes MPI decomposition', &
                         decomposition_matches(ads) .and. &
                         all(ads%lnelem == ads%maxe - ads%mine + 1), &
                         failure_count)
      call assert_global('initialize_setup builds open knots and basis tables', &
                         numerical_setup_matches(ads), failure_count)
      call assert_global('initialize_setup has exact linear midpoint tables', &
                         linear_midpoint_setup_matches(ads), failure_count)

      call Cleanup_ADS(ads, cleanup_ierr)
      call assert_global('Cleanup_ADS releases a fully initialized setup', &
                         cleanup_ierr == 0 .and. setup_is_clean(ads), &
                         failure_count)
   end subroutine test_initialize_setup


   subroutine test_full_initialize_and_reinitialize(failure_count)
      integer(kind=4), intent(inout) :: failure_count
      type(ADS_Setup) :: test_space, trial_space
      type(ADS_compute_data) :: data
      integer(kind=4), parameter :: nelem1(3) = (/4, 5, 6/)
      integer(kind=4), parameter :: ptest1(3) = (/2, 1, 3/)
      integer(kind=4), parameter :: ptrial1(3) = (/1, 2, 1/)
      integer(kind=4), parameter :: test_ng1(3) = (/3, 2, 4/)
      integer(kind=4), parameter :: trial_ng1(3) = (/3, 3, 4/)
      integer(kind=4), parameter :: nelem2(3) = (/4, 4, 5/)
      integer(kind=4), parameter :: ptest2(3) = (/1, 2, 1/)
      integer(kind=4), parameter :: ptrial2(3) = (/2, 1, 2/)
      integer(kind=4), parameter :: test_ng2(3) = (/2, 3, 2/)
      integer(kind=4), parameter :: trial_ng2(3) = (/3, 3, 3/)
      integer(kind=4) :: lifecycle_ierr
      logical :: local_ok

      call initialize(nelem1, ptest1, ptrial1, ptrial1 - 1, &
                      test_space, trial_space, data, lifecycle_ierr)

      local_ok = lifecycle_ierr == 0
      local_ok = local_ok .and. initialized_storage_matches( &
         test_space, nelem1 + ptest1 - 1, ptest1, test_ng1, nelem1)
      local_ok = local_ok .and. initialized_storage_matches( &
         trial_space, nelem1 + ptrial1 - 1, ptrial1, trial_ng1, nelem1)
      local_ok = local_ok .and. decomposition_matches(test_space)
      local_ok = local_ok .and. decomposition_matches(trial_space)
      local_ok = local_ok .and. numerical_setup_matches(test_space)
      local_ok = local_ok .and. numerical_setup_matches(trial_space)
      call assert_global('initialize builds test and trial spaces with promoted quadrature', &
                         local_ok, failure_count)

      call assert_global('initialize allocates zeroed shared runtime state', &
                         data_core_matches(data, &
                            min(test_space%mine, trial_space%mine), &
                            max(test_space%maxe, trial_space%maxe), &
                            max(test_ng1, trial_ng1)), failure_count)

      call allocate_optional_rhs(data)
      call initialize(nelem2, ptest2, ptrial2, ptrial2 - 1, &
                      test_space, trial_space, data, lifecycle_ierr)

      local_ok = lifecycle_ierr == 0
      local_ok = local_ok .and. initialized_storage_matches( &
         test_space, nelem2 + ptest2 - 1, ptest2, test_ng2, nelem2)
      local_ok = local_ok .and. initialized_storage_matches( &
         trial_space, nelem2 + ptrial2 - 1, ptrial2, trial_ng2, nelem2)
      local_ok = local_ok .and. decomposition_matches(test_space)
      local_ok = local_ok .and. decomposition_matches(trial_space)
      local_ok = local_ok .and. numerical_setup_matches(test_space)
      local_ok = local_ok .and. numerical_setup_matches(trial_space)
      local_ok = local_ok .and. data_core_matches(data, &
         min(test_space%mine, trial_space%mine), &
         max(test_space%maxe, trial_space%maxe), &
         max(test_ng2, trial_ng2))
      call assert_global('initialize safely replaces existing allocatable state', &
                         local_ok, failure_count)

      call Cleanup_data(data, lifecycle_ierr)
      local_ok = lifecycle_ierr == 0 .and. data_is_clean(data)
      call Cleanup_ADS(test_space, lifecycle_ierr)
      local_ok = local_ok .and. lifecycle_ierr == 0 .and. &
                 setup_is_clean(test_space)
      call Cleanup_ADS(trial_space, lifecycle_ierr)
      local_ok = local_ok .and. lifecycle_ierr == 0 .and. &
                 setup_is_clean(trial_space)
      call assert_global('full lifecycle cleanup leaves no allocation behind', &
                         local_ok, failure_count)
   end subroutine test_full_initialize_and_reinitialize


   logical function allocate_ads_storage_matches(ads, n, nelem, p, ng) &
      result(matches)
      type(ADS_Setup), intent(in) :: ads
      integer(kind=4), intent(in) :: n(3), nelem(3), p(3), ng(3)

      matches = allocated(ads%Ox) .and. allocated(ads%Oy) .and. &
         allocated(ads%Oz) .and. allocated(ads%Jx) .and. &
         allocated(ads%Jy) .and. allocated(ads%Jz) .and. &
         allocated(ads%Xx) .and. allocated(ads%Xy) .and. &
         allocated(ads%Xz) .and. allocated(ads%NNx) .and. &
         allocated(ads%NNy) .and. allocated(ads%NNz) .and. &
         allocated(ads%Wx) .and. allocated(ads%Wy) .and. &
         allocated(ads%Wz) .and. allocated(ads%IPIVx) .and. &
         allocated(ads%IPIVy) .and. allocated(ads%IPIVz)
      if (.not. matches) return

      matches = size(ads%Ox) == nelem(1) .and. &
         size(ads%Oy) == nelem(2) .and. size(ads%Oz) == nelem(3)
      matches = matches .and. size(ads%Jx) == nelem(1) .and. &
         size(ads%Jy) == nelem(2) .and. size(ads%Jz) == nelem(3)
      matches = matches .and. all(shape(ads%Xx) == (/ng(1), nelem(1)/))
      matches = matches .and. all(shape(ads%Xy) == (/ng(2), nelem(2)/))
      matches = matches .and. all(shape(ads%Xz) == (/ng(3), nelem(3)/))
      matches = matches .and. &
         all(shape(ads%NNx) == (/2, p(1) + 1, ng(1), nelem(1)/))
      matches = matches .and. &
         all(shape(ads%NNy) == (/2, p(2) + 1, ng(2), nelem(2)/))
      matches = matches .and. &
         all(shape(ads%NNz) == (/2, p(3) + 1, ng(3), nelem(3)/))
      matches = matches .and. &
         all(lbound(ads%NNx) == (/0, 0, 1, 1/))
      matches = matches .and. &
         all(lbound(ads%NNy) == (/0, 0, 1, 1/))
      matches = matches .and. &
         all(lbound(ads%NNz) == (/0, 0, 1, 1/))
      matches = matches .and. all(lbound(ads%Xx) == 1) .and. &
         all(lbound(ads%Xy) == 1) .and. all(lbound(ads%Xz) == 1)
      matches = matches .and. lbound(ads%Ox, 1) == 1 .and. &
         lbound(ads%Oy, 1) == 1 .and. lbound(ads%Oz, 1) == 1
      matches = matches .and. lbound(ads%Jx, 1) == 1 .and. &
         lbound(ads%Jy, 1) == 1 .and. lbound(ads%Jz, 1) == 1
      matches = matches .and. size(ads%Wx) == ng(1) .and. &
         size(ads%Wy) == ng(2) .and. size(ads%Wz) == ng(3)
      matches = matches .and. lbound(ads%Wx, 1) == 1 .and. &
         lbound(ads%Wy, 1) == 1 .and. lbound(ads%Wz, 1) == 1
      matches = matches .and. size(ads%IPIVx) == n(1) + 1 .and. &
         size(ads%IPIVy) == n(2) + 1 .and. size(ads%IPIVz) == n(3) + 1
      matches = matches .and. lbound(ads%IPIVx, 1) == 1 .and. &
         lbound(ads%IPIVy, 1) == 1 .and. lbound(ads%IPIVz, 1) == 1
      matches = matches .and. .not. allocated(ads%Ux) .and. &
         .not. allocated(ads%Uy) .and. .not. allocated(ads%Uz)
      matches = matches .and. .not. allocated(ads%dimensionsX) .and. &
         .not. allocated(ads%dimensionsY) .and. &
         .not. allocated(ads%dimensionsZ)
      matches = matches .and. .not. allocated(ads%shiftsX) .and. &
         .not. allocated(ads%shiftsY) .and. .not. allocated(ads%shiftsZ)
   end function allocate_ads_storage_matches


   logical function initialized_storage_matches(ads, n, p, ng, nelem) &
      result(matches)
      type(ADS_Setup), intent(in) :: ads
      integer(kind=4), intent(in) :: n(3), p(3), ng(3), nelem(3)

      matches = allocated(ads%Ux) .and. allocated(ads%Uy) .and. &
         allocated(ads%Uz) .and. allocated(ads%dimensionsX) .and. &
         allocated(ads%dimensionsY) .and. allocated(ads%dimensionsZ) .and. &
         allocated(ads%shiftsX) .and. allocated(ads%shiftsY) .and. &
         allocated(ads%shiftsZ)
      if (.not. matches) return
      matches = allocate_ads_storage_matches_initialized(ads, n, nelem, p, ng)
      if (.not. matches) return

      matches = all(ads%n == n) .and. all(ads%p == p) .and. &
         all(ads%ng == ng) .and. all(ads%nelem == nelem) .and. &
         all(ads%m == n + p + 1) .and. ads%tau == 0.d0
      matches = matches .and. size(ads%Ux) == n(1) + p(1) + 2 .and. &
         size(ads%Uy) == n(2) + p(2) + 2 .and. &
         size(ads%Uz) == n(3) + p(3) + 2
      matches = matches .and. lbound(ads%Ux, 1) == 1 .and. &
         lbound(ads%Uy, 1) == 1 .and. lbound(ads%Uz, 1) == 1
      matches = matches .and. size(ads%dimensionsX) == NRPROCX .and. &
         size(ads%dimensionsY) == NRPROCY .and. &
         size(ads%dimensionsZ) == NRPROCZ
      matches = matches .and. size(ads%shiftsX) == NRPROCX .and. &
         size(ads%shiftsY) == NRPROCY .and. size(ads%shiftsZ) == NRPROCZ
      matches = matches .and. lbound(ads%dimensionsX, 1) == 1 .and. &
         lbound(ads%dimensionsY, 1) == 1 .and. &
         lbound(ads%dimensionsZ, 1) == 1
      matches = matches .and. lbound(ads%shiftsX, 1) == 1 .and. &
         lbound(ads%shiftsY, 1) == 1 .and. lbound(ads%shiftsZ, 1) == 1
   end function initialized_storage_matches


   logical function allocate_ads_storage_matches_initialized(ads, n, nelem, p, ng) &
      result(matches)
      type(ADS_Setup), intent(in) :: ads
      integer(kind=4), intent(in) :: n(3), nelem(3), p(3), ng(3)

      matches = allocated(ads%Ox) .and. allocated(ads%Oy) .and. &
         allocated(ads%Oz) .and. allocated(ads%Jx) .and. &
         allocated(ads%Jy) .and. allocated(ads%Jz) .and. &
         allocated(ads%Xx) .and. allocated(ads%Xy) .and. &
         allocated(ads%Xz) .and. allocated(ads%NNx) .and. &
         allocated(ads%NNy) .and. allocated(ads%NNz) .and. &
         allocated(ads%Wx) .and. allocated(ads%Wy) .and. &
         allocated(ads%Wz) .and. allocated(ads%IPIVx) .and. &
         allocated(ads%IPIVy) .and. allocated(ads%IPIVz)
      if (.not. matches) return

      matches = size(ads%Ox) == nelem(1) .and. &
         size(ads%Oy) == nelem(2) .and. size(ads%Oz) == nelem(3)
      matches = matches .and. size(ads%Jx) == nelem(1) .and. &
         size(ads%Jy) == nelem(2) .and. size(ads%Jz) == nelem(3)
      matches = matches .and. all(shape(ads%Xx) == (/ng(1), nelem(1)/))
      matches = matches .and. all(shape(ads%Xy) == (/ng(2), nelem(2)/))
      matches = matches .and. all(shape(ads%Xz) == (/ng(3), nelem(3)/))
      matches = matches .and. &
         all(shape(ads%NNx) == (/2, p(1) + 1, ng(1), nelem(1)/))
      matches = matches .and. &
         all(shape(ads%NNy) == (/2, p(2) + 1, ng(2), nelem(2)/))
      matches = matches .and. &
         all(shape(ads%NNz) == (/2, p(3) + 1, ng(3), nelem(3)/))
      matches = matches .and. &
         all(lbound(ads%NNx) == (/0, 0, 1, 1/))
      matches = matches .and. &
         all(lbound(ads%NNy) == (/0, 0, 1, 1/))
      matches = matches .and. &
         all(lbound(ads%NNz) == (/0, 0, 1, 1/))
      matches = matches .and. all(lbound(ads%Xx) == 1) .and. &
         all(lbound(ads%Xy) == 1) .and. all(lbound(ads%Xz) == 1)
      matches = matches .and. lbound(ads%Ox, 1) == 1 .and. &
         lbound(ads%Oy, 1) == 1 .and. lbound(ads%Oz, 1) == 1
      matches = matches .and. lbound(ads%Jx, 1) == 1 .and. &
         lbound(ads%Jy, 1) == 1 .and. lbound(ads%Jz, 1) == 1
      matches = matches .and. size(ads%Wx) == ng(1) .and. &
         size(ads%Wy) == ng(2) .and. size(ads%Wz) == ng(3)
      matches = matches .and. lbound(ads%Wx, 1) == 1 .and. &
         lbound(ads%Wy, 1) == 1 .and. lbound(ads%Wz, 1) == 1
      matches = matches .and. size(ads%IPIVx) == n(1) + 1 .and. &
         size(ads%IPIVy) == n(2) + 1 .and. size(ads%IPIVz) == n(3) + 1
      matches = matches .and. lbound(ads%IPIVx, 1) == 1 .and. &
         lbound(ads%IPIVy, 1) == 1 .and. lbound(ads%IPIVz, 1) == 1
   end function allocate_ads_storage_matches_initialized


   logical function decomposition_matches(ads) result(matches)
      type(ADS_Setup), intent(in) :: ads
      integer(kind=4) :: axis, expected_ibeg(3), expected_iend(3)
      integer(kind=4) :: expected_mine(3), expected_maxe(3)
      integer(kind=4) :: expected_nrcpp(3), expected_s(3)
      integer(kind=4) :: coords(3), procs(3), stride

      coords = (/MYRANKX, MYRANKY, MYRANKZ/)
      procs = (/NRPROCX, NRPROCY, NRPROCZ/)
      do axis = 1, 3
         expected_nrcpp(axis) = (ads%n(axis) + procs(axis))/procs(axis)
         expected_ibeg(axis) = expected_nrcpp(axis)*coords(axis) + 1
         expected_iend(axis) = min(expected_nrcpp(axis)*(coords(axis) + 1), &
                                   ads%n(axis) + 1)
         expected_mine(axis) = max(expected_ibeg(axis) - ads%p(axis) - 1, 1)
         expected_maxe(axis) = min(expected_iend(axis), &
                                   ads%n(axis) + 1 - ads%p(axis))
         expected_s(axis) = expected_iend(axis) - expected_ibeg(axis) + 1
      end do

      matches = all(ads%nrcpp == expected_nrcpp) .and. &
         all(ads%ibeg == expected_ibeg) .and. all(ads%iend == expected_iend)
      matches = matches .and. all(ads%mine == expected_mine) .and. &
         all(ads%maxe == expected_maxe) .and. all(ads%s == expected_s)
      matches = matches .and. all(expected_s > 0) .and. &
         all(expected_mine <= expected_maxe)
      if (.not. matches) return

      stride = expected_s(2)*expected_s(3)
      matches = dim_vector_matches(ads%dimensionsX, ads%shiftsX, &
         ads%n(1), expected_nrcpp(1), NRPROCX, stride)
      stride = expected_s(1)*expected_s(3)
      matches = matches .and. dim_vector_matches( &
         ads%dimensionsY, ads%shiftsY, ads%n(2), expected_nrcpp(2), &
         NRPROCY, stride)
      stride = expected_s(1)*expected_s(2)
      matches = matches .and. dim_vector_matches( &
         ads%dimensionsZ, ads%shiftsZ, ads%n(3), expected_nrcpp(3), &
         NRPROCZ, stride)

      matches = matches .and. neighbour_vector_matches( &
         ads%ibegsx, ads%iendsx, ads%n(1), NRPROCX, MYRANKX)
      matches = matches .and. neighbour_vector_matches( &
         ads%ibegsy, ads%iendsy, ads%n(2), NRPROCY, MYRANKY)
      matches = matches .and. neighbour_vector_matches( &
         ads%ibegsz, ads%iendsz, ads%n(3), NRPROCZ, MYRANKZ)
   end function decomposition_matches


   logical function dim_vector_matches(dims, shifts, n, nrcpp, nrproc, stride) &
      result(matches)
      integer(kind=4), allocatable, intent(in) :: dims(:), shifts(:)
      integer(kind=4), intent(in) :: n, nrcpp, nrproc, stride
      integer(kind=4) :: expected_count, expected_shift, i

      matches = allocated(dims) .and. allocated(shifts)
      if (.not. matches) return
      matches = size(dims) == nrproc .and. size(shifts) == nrproc
      if (.not. matches) return

      expected_shift = 0
      do i = 1, nrproc
         expected_count = min(nrcpp, n + 1 - nrcpp*(i - 1))*stride
         matches = matches .and. expected_count > 0
         matches = matches .and. dims(i) == expected_count
         matches = matches .and. shifts(i) == expected_shift
         expected_shift = expected_shift + expected_count
      end do
      matches = matches .and. expected_shift == (n + 1)*stride
   end function dim_vector_matches


   logical function neighbour_vector_matches(begins, ends, n, nrproc, rank) &
      result(matches)
      integer(kind=4), intent(in) :: begins(3), ends(3)
      integer(kind=4), intent(in) :: n, nrproc, rank
      integer(kind=4) :: expected_begins(3), expected_ends(3)
      integer(kind=4) :: neighbour, nrcpp, slot

      expected_begins = -1
      expected_ends = -1
      nrcpp = (n + nrproc)/nrproc
      do neighbour = max(rank - 1, 0), min(rank + 1, nrproc - 1)
         slot = neighbour - rank + 2
         expected_begins(slot) = nrcpp*neighbour + 1
         expected_ends(slot) = min(nrcpp*(neighbour + 1), n + 1)
      end do
      matches = all(begins == expected_begins) .and. &
                all(ends == expected_ends)
   end function neighbour_vector_matches


   logical function numerical_setup_matches(ads) result(matches)
      type(ADS_Setup), intent(in) :: ads

      matches = allocated(ads%Ux) .and. allocated(ads%Uy) .and. &
         allocated(ads%Uz) .and. allocated(ads%Ox) .and. &
         allocated(ads%Oy) .and. allocated(ads%Oz) .and. &
         allocated(ads%Jx) .and. allocated(ads%Jy) .and. &
         allocated(ads%Jz) .and. allocated(ads%Wx) .and. &
         allocated(ads%Wy) .and. allocated(ads%Wz) .and. &
         allocated(ads%Xx) .and. allocated(ads%Xy) .and. &
         allocated(ads%Xz) .and. allocated(ads%NNx) .and. &
         allocated(ads%NNy) .and. allocated(ads%NNz)
      if (.not. matches) return

      matches = open_knot_matches(ads%Ux, ads%n(1), ads%p(1), ads%nelem(1))
      matches = matches .and. open_knot_matches( &
         ads%Uy, ads%n(2), ads%p(2), ads%nelem(2))
      matches = matches .and. open_knot_matches( &
         ads%Uz, ads%n(3), ads%p(3), ads%nelem(3))
      matches = matches .and. basis_axis_matches( &
         ads%Ox, ads%Jx, ads%Wx, ads%Xx, ads%NNx, &
         ads%nelem(1), ads%p(1), ads%ng(1))
      matches = matches .and. basis_axis_matches( &
         ads%Oy, ads%Jy, ads%Wy, ads%Xy, ads%NNy, &
         ads%nelem(2), ads%p(2), ads%ng(2))
      matches = matches .and. basis_axis_matches( &
         ads%Oz, ads%Jz, ads%Wz, ads%Xz, ads%NNz, &
         ads%nelem(3), ads%p(3), ads%ng(3))
   end function numerical_setup_matches


   logical function open_knot_matches(U, n, p, nelem) result(matches)
      real(kind=8), intent(in) :: U(:)
      integer(kind=4), intent(in) :: n, p, nelem
      integer(kind=4) :: e, index

      matches = size(U) == n + p + 2
      if (.not. matches) return
      matches = all(U(1:p + 1) == 0.d0)
      do e = 1, nelem - 1
         index = p + 1 + e
         matches = matches .and. &
            abs(U(index) - real(e, kind=8)/real(nelem, kind=8)) <= TOL
      end do
      matches = matches .and. &
         all(U(size(U) - p:size(U)) == 1.d0)
      matches = matches .and. all(U(2:) >= U(:size(U) - 1))
   end function open_knot_matches


   logical function linear_midpoint_setup_matches(ads) result(matches)
      type(ADS_Setup), intent(in) :: ads

      matches = allocated(ads%Ox) .and. allocated(ads%Oy) .and. &
         allocated(ads%Oz) .and. allocated(ads%Jx) .and. &
         allocated(ads%Jy) .and. allocated(ads%Jz) .and. &
         allocated(ads%Wx) .and. allocated(ads%Wy) .and. &
         allocated(ads%Wz) .and. allocated(ads%Xx) .and. &
         allocated(ads%Xy) .and. allocated(ads%Xz) .and. &
         allocated(ads%NNx) .and. allocated(ads%NNy) .and. &
         allocated(ads%NNz)
      if (.not. matches) return
      matches = all(ads%p == 1) .and. all(ads%ng == 1)
      if (.not. matches) return

      matches = linear_midpoint_axis_matches( &
         ads%Ox, ads%Jx, ads%Wx, ads%Xx, ads%NNx, ads%nelem(1))
      matches = matches .and. linear_midpoint_axis_matches( &
         ads%Oy, ads%Jy, ads%Wy, ads%Xy, ads%NNy, ads%nelem(2))
      matches = matches .and. linear_midpoint_axis_matches( &
         ads%Oz, ads%Jz, ads%Wz, ads%Xz, ads%NNz, ads%nelem(3))
   end function linear_midpoint_setup_matches


   logical function linear_midpoint_axis_matches(O, J, W, X, N, nelem) &
      result(matches)
      integer(kind=4), intent(in) :: O(:)
      real(kind=8), intent(in) :: J(:), W(:), X(:, :)
      real(kind=8), intent(in) :: N(0:, 0:, :, :)
      integer(kind=4), intent(in) :: nelem
      integer(kind=4) :: e
      real(kind=8) :: expected_x

      matches = size(O) == nelem .and. size(J) == nelem .and. &
         size(W) == 1 .and. all(shape(X) == (/1, nelem/)) .and. &
         all(shape(N) == (/2, 2, 1, nelem/))
      if (.not. matches) return
      matches = W(1) == 2.d0
      do e = 1, nelem
         expected_x = (real(e, kind=8) - 0.5d0)/real(nelem, kind=8)
         matches = matches .and. O(e) == e - 1
         matches = matches .and. &
            abs(J(e) - 0.5d0/real(nelem, kind=8)) <= TOL
         matches = matches .and. abs(X(1, e) - expected_x) <= TOL
         matches = matches .and. &
            all(abs(N(0, :, 1, e) - (/0.5d0, 0.5d0/)) <= TOL)
         matches = matches .and. &
            all(abs(N(1, :, 1, e) - &
               (/-real(nelem, kind=8), real(nelem, kind=8)/)) <= TOL)
      end do
   end function linear_midpoint_axis_matches


   logical function basis_axis_matches(O, J, W, X, N, nelem, p, ng) &
      result(matches)
      integer(kind=4), intent(in) :: O(:)
      real(kind=8), intent(in) :: J(:), W(:), X(:, :)
      real(kind=8), intent(in) :: N(0:, 0:, :, :)
      integer(kind=4), intent(in) :: nelem, p, ng
      integer(kind=4) :: e, k
      real(kind=8) :: left, right, scale

      matches = all(ieee_is_finite(J)) .and. all(ieee_is_finite(W)) .and. &
         all(ieee_is_finite(X)) .and. all(ieee_is_finite(N))
      if (.not. matches) return
      matches = size(O) == nelem .and. size(J) == nelem .and. &
         size(W) == ng .and. all(shape(X) == (/ng, nelem/))
      matches = matches .and. &
         all(shape(N) == (/2, p + 1, ng, nelem/))
      if (.not. matches) return

      scale = max(1.d0, real(nelem, kind=8))
      matches = all(O == (/(e - 1, e = 1, nelem)/))
      matches = matches .and. &
         all(abs(J - 0.5d0/real(nelem, kind=8)) <= TOL)
      matches = matches .and. all(W > 0.d0) .and. &
         abs(sum(W) - 2.d0) <= TOL
      do e = 1, nelem
         left = real(e - 1, kind=8)/real(nelem, kind=8)
         right = real(e, kind=8)/real(nelem, kind=8)
         matches = matches .and. all(X(:, e) > left) .and. &
            all(X(:, e) < right)
         do k = 1, ng
            matches = matches .and. &
               abs(sum(N(0, :, k, e)) - 1.d0) <= TOL
            matches = matches .and. &
               abs(sum(N(1, :, k, e))) <= TOL*scale
         end do
      end do
   end function basis_axis_matches


   logical function data_core_matches(data, expected_mine, expected_maxe, &
                                      expected_ng) &
      result(matches)
      type(ADS_compute_data), intent(in) :: data
      integer(kind=4), intent(in) :: expected_mine(3), expected_maxe(3)
      integer(kind=4), intent(in) :: expected_ng(3)
      integer(kind=4) :: extents(3)

      matches = allocated(data%Un) .and. allocated(data%Un13) .and. &
         allocated(data%Un23) .and. allocated(data%dUn) .and. &
         allocated(data%dUn0) .and. allocated(data%dUn13) .and. &
         allocated(data%dUn23) .and. allocated(data%R)
      if (.not. matches) return

      extents = expected_maxe - expected_mine + 1
      matches = all(data%state_mine == expected_mine) .and. &
         all(data%state_maxe == expected_maxe) .and. &
         all(data%state_lnelem == extents)
      matches = matches .and. all(shape(data%Un) == &
         (/extents, expected_ng/))
      matches = matches .and. all(shape(data%Un13) == &
         (/extents, expected_ng/))
      matches = matches .and. all(shape(data%Un23) == &
         (/extents, expected_ng/))
      matches = matches .and. all(shape(data%dUn) == &
         (/extents, expected_ng, 3/))
      matches = matches .and. all(shape(data%dUn0) == &
         (/extents, expected_ng, 3/))
      matches = matches .and. all(shape(data%dUn13) == &
         (/extents, expected_ng, 3/))
      matches = matches .and. all(shape(data%dUn23) == &
         (/extents, expected_ng, 3/))
      matches = matches .and. all(data%halo_end >= data%halo_begin)
      matches = matches .and. all(shape(data%R) == &
         (/data%halo_end - data%halo_begin + 1, 1/))
      matches = matches .and. allocated(data%halo_send_begin) .and. &
         allocated(data%halo_send_end) .and. &
         allocated(data%halo_recv_begin) .and. &
         allocated(data%halo_recv_end) .and. &
         allocated(data%halo_send_count) .and. &
         allocated(data%halo_send_displ) .and. &
         allocated(data%halo_recv_count) .and. &
         allocated(data%halo_recv_displ) .and. &
         allocated(data%halo_send_buffer) .and. &
         allocated(data%halo_recv_buffer) .and. &
         allocated(data%halo_requests) .and. &
         allocated(data%halo_statuses)
      if (.not. matches) return
      matches = all(shape(data%halo_send_begin) == (/3, NRPROC/)) .and. &
         all(shape(data%halo_send_end) == (/3, NRPROC/)) .and. &
         all(shape(data%halo_recv_begin) == (/3, NRPROC/)) .and. &
         all(shape(data%halo_recv_end) == (/3, NRPROC/))
      matches = matches .and. size(data%halo_send_count) == NRPROC .and. &
         size(data%halo_send_displ) == NRPROC .and. &
         size(data%halo_recv_count) == NRPROC .and. &
         size(data%halo_recv_displ) == NRPROC
      matches = matches .and. all(data%halo_send_count >= 0) .and. &
         all(data%halo_recv_count >= 0)
      matches = matches .and. size(data%halo_send_buffer) >= &
         max(sum(data%halo_send_count), 1)
      matches = matches .and. size(data%halo_recv_buffer) >= &
         max(sum(data%halo_recv_count), 1)
      matches = matches .and. size(data%halo_requests) >= &
         max(2*(NRPROC - 1), 1)
      matches = matches .and. all(shape(data%halo_statuses) == &
         (/MPI_STATUS_SIZE, max(2*(NRPROC - 1), 1)/))
      matches = matches .and. all(lbound(data%Un) == 1) .and. &
         all(lbound(data%Un13) == 1) .and. &
         all(lbound(data%Un23) == 1)
      matches = matches .and. all(lbound(data%dUn) == 1) .and. &
         all(lbound(data%dUn0) == 1) .and. &
         all(lbound(data%dUn13) == 1) .and. &
         all(lbound(data%dUn23) == 1) .and. all(lbound(data%R) == 1)
      matches = matches .and. all(data%Un == 0.d0) .and. &
         all(data%Un13 == 0.d0) .and. all(data%Un23 == 0.d0)
      matches = matches .and. all(data%dUn == 0.d0) .and. &
         all(data%dUn0 == 0.d0) .and. all(data%dUn13 == 0.d0) .and. &
         all(data%dUn23 == 0.d0) .and. all(data%R == 0.d0)
      matches = matches .and. all(data%rhs_du_state == 0)
      matches = matches .and. optional_rhs_is_unallocated(data)
   end function data_core_matches


   logical function halo_plan_matches(data, trial_space) result(matches)
      type(ADS_compute_data), intent(in) :: data
      type(ADS_Setup), intent(in) :: trial_space
      integer(kind=4) :: expected_begin(3), expected_end(3)

      expected_begin = (/ &
         trial_space%Ox(data%state_mine(1)), &
         trial_space%Oy(data%state_mine(2)), &
         trial_space%Oz(data%state_mine(3))/)
      expected_end = (/ &
         trial_space%Ox(data%state_maxe(1)) + trial_space%p(1), &
         trial_space%Oy(data%state_maxe(2)) + trial_space%p(2), &
         trial_space%Oz(data%state_maxe(3)) + trial_space%p(3)/)

      matches = all(data%halo_begin == expected_begin) .and. &
         all(data%halo_end == expected_end)
      matches = matches .and. all(shape(data%R) == &
         (/expected_end - expected_begin + 1, 1/))
   end function halo_plan_matches


   logical function optional_rhs_is_unallocated(data) result(matches)
      type(ADS_compute_data), intent(in) :: data

      matches = .not. allocated(data%F) .and. .not. allocated(data%F2) .and. &
         .not. allocated(data%F3) .and. .not. allocated(data%FF) .and. &
         .not. allocated(data%Ft) .and. .not. allocated(data%Ft2) .and. &
         .not. allocated(data%Ft3) .and. .not. allocated(data%FFt)
   end function optional_rhs_is_unallocated


   subroutine allocate_optional_rhs(data)
      type(ADS_compute_data), intent(inout) :: data

      allocate(data%F(2, 3), data%F2(3, 2), data%F3(1, 4), data%FF(4, 1))
      allocate(data%Ft(3, 4), data%Ft2(4, 3), data%Ft3(2, 5), data%FFt(5, 2))
      data%F = 1.d0
      data%F2 = 2.d0
      data%F3 = 3.d0
      data%FF = 4.d0
      data%Ft = 5.d0
      data%Ft2 = 6.d0
      data%Ft3 = 7.d0
      data%FFt = 8.d0
   end subroutine allocate_optional_rhs


   logical function setup_is_clean(ads) result(clean)
      type(ADS_Setup), intent(in) :: ads

      clean = .not. allocated(ads%Ux) .and. .not. allocated(ads%Uy) .and. &
         .not. allocated(ads%Uz) .and. &
         .not. allocated(ads%dimensionsX) .and. &
         .not. allocated(ads%dimensionsY) .and. &
         .not. allocated(ads%dimensionsZ) .and. &
         .not. allocated(ads%shiftsX) .and. &
         .not. allocated(ads%shiftsY) .and. &
         .not. allocated(ads%shiftsZ) .and. &
         .not. allocated(ads%IPIVx) .and. &
         .not. allocated(ads%IPIVy) .and. &
         .not. allocated(ads%IPIVz) .and. &
         .not. allocated(ads%Ox) .and. .not. allocated(ads%Oy) .and. &
         .not. allocated(ads%Oz) .and. .not. allocated(ads%Jx) .and. &
         .not. allocated(ads%Jy) .and. .not. allocated(ads%Jz) .and. &
         .not. allocated(ads%Xx) .and. .not. allocated(ads%Xy) .and. &
         .not. allocated(ads%Xz) .and. .not. allocated(ads%NNx) .and. &
         .not. allocated(ads%NNy) .and. .not. allocated(ads%NNz) .and. &
         .not. allocated(ads%Wx) .and. .not. allocated(ads%Wy) .and. &
         .not. allocated(ads%Wz)
   end function setup_is_clean


   logical function data_is_clean(data) result(clean)
      type(ADS_compute_data), intent(in) :: data

      clean = optional_rhs_is_unallocated(data) .and. &
         .not. allocated(data%R) .and. .not. allocated(data%Un) .and. &
         .not. allocated(data%Un13) .and. .not. allocated(data%Un23) .and. &
         .not. allocated(data%dUn) .and. .not. allocated(data%dUn0) .and. &
         .not. allocated(data%dUn13) .and. .not. allocated(data%dUn23) .and. &
         .not. allocated(data%halo_send_begin) .and. &
         .not. allocated(data%halo_send_end) .and. &
         .not. allocated(data%halo_recv_begin) .and. &
         .not. allocated(data%halo_recv_end) .and. &
         .not. allocated(data%halo_send_count) .and. &
         .not. allocated(data%halo_send_displ) .and. &
         .not. allocated(data%halo_recv_count) .and. &
         .not. allocated(data%halo_recv_displ) .and. &
         .not. allocated(data%halo_send_buffer) .and. &
         .not. allocated(data%halo_recv_buffer) .and. &
         .not. allocated(data%halo_requests) .and. &
         .not. allocated(data%halo_statuses)
   end function data_is_clean


   subroutine assert_global(label, local_ok, failure_count)
      character(len=*), intent(in) :: label
      logical, intent(in) :: local_ok
      integer(kind=4), intent(inout) :: failure_count
      integer(kind=4) :: global_flag, local_flag, reduce_ierr
      logical :: global_ok

      local_flag = merge(1, 0, local_ok)
      call MPI_Allreduce(local_flag, global_flag, 1, MPI_INTEGER, MPI_MIN, &
                         MPI_COMM_WORLD, reduce_ierr)
      global_ok = reduce_ierr == MPI_SUCCESS .and. global_flag == 1

      checks = checks + 1
      if (.not. global_ok) failure_count = failure_count + 1
      if (MYRANK == 0) then
         if (global_ok) then
            write (*, '(A,A)') 'PASS ', trim(label)
         else
            write (*, '(A,A)') 'FAIL ', trim(label)
         end if
      end if
   end subroutine assert_global


   subroutine read_process_grid(grid)
      integer(kind=4), intent(out) :: grid(3)
      character(len=32) :: argument
      integer(kind=4) :: i, read_status

      if (command_argument_count() /= 3) then
         write (*, '(A)') 'Usage: test_ads_lifecycle <procx> <procy> <procz>'
         stop 2
      end if
      do i = 1, 3
         call get_command_argument(i, argument)
         read (argument, *, iostat=read_status) grid(i)
         if (read_status /= 0 .or. grid(i) <= 0) then
            write (*, '(A,I0,A,A)') 'Invalid process-grid argument ', i, &
                                    ': ', trim(argument)
            stop 2
         end if
      end do
   end subroutine read_process_grid

end program test_ads_lifecycle
