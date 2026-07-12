!------------------------------------------------------------------------------
!
! MODULE: solution_output
!
! DESCRIPTION:
!> @file solution_output.F90
!> @brief Solution gathering and visualization-output orchestration.
!>
!> @details
!> This module gathers distributed spline coefficients, broadcasts the
!> full coefficient field needed for distributed sampling, and delegates
!> the actual file serialization to the plotting/output backends.
!>
!> The module separates visualization coordination from the time-stepping
!> workflow while preserving the historical \ref PrintSolution entry
!> point re-exported by \ref ADSS. The current default output path writes
!> VTK image-data files through \ref VtkOutput.
!>
!> MPI-aware sampling remains in the plotting layer. This module is the
!> bridge between ADS solution storage and those backend-neutral output
!> routines.
!
!------------------------------------------------------------------------------
module solution_output

   implicit none

   private
   public :: PrintSolution

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Gathers the global solution and exports it for visualization.
!>
!> @details
!> This routine reconstructs the full spline solution on rank zero using
!> `GatherFullSolution`, broadcasts the coefficient tensor to all ranks,
!> prepares a plotting-parameter structure, and delegates MPI-distributed
!> sampling/export to `SaveSplinePlotMPI`, which writes VTK output via
!> \ref VtkOutput on rank zero.
!>
!> The generated filename has the form `stepNNNN...` derived from the
!> iteration counter.
!
! Input:
! ------
!> @param[in] iter
!> Iteration or time-step index used in the output filename.
!>
!> @param[in] ads
!> ADS setup structure containing the spline-space definition.
!>
!> @param[in] part
!> Local piece of the solution to be gathered.
!
!---------------------------------------------------------------------------
subroutine PrintSolution(iter, ads, part)
      use Setup, ONLY: ADS_Setup
      use parallelism, ONLY: MYRANK
      use plot, ONLY: SaveSplinePlotMPI, PlotParams
      use vtk, ONLY: VtkOutput
      use my_mpi, ONLY: GatherFullSolution
      use mpi
      implicit none
!> @brief Local piece of the solution.
      real(kind=8), dimension(:,:), intent(in) :: part
!> @brief ADS setup structure used for spline evaluation and export.
      type(ADS_setup), intent(in) :: ads
!> @brief Iteration or time-step index.
      integer(kind=4), intent(in) :: iter
      real(kind=8), allocatable :: solution(:, :, :)
      type(PlotParams) :: params
      character(len=20) :: filename
      integer(kind=4) :: coefficient_count, ierr

      call GatherFullSolution(0, part, solution, &
                              ads%n, ads%p, ads%s)

      coefficient_count = (ads%n(1) + 1)*(ads%n(2) + 1)*(ads%n(3) + 1)
      if (MYRANK /= 0) then
            allocate (solution(0:ads%n(1), 0:ads%n(2), 0:ads%n(3)))
      end if
      call MPI_Bcast(solution, coefficient_count, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      if (MYRANK == 0) then
            write (filename, '(I10)') iter
            filename = 'step'//adjustl(filename)
            ! filename = trim(filename) // '_'
      end if

      call MPI_Bcast(filename, len(filename), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      params = PlotParams(0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 31, 31, 31)
      call SaveSplinePlotMPI(trim(filename), &
                           ads%Ux, ads%p(1), ads%n(1), ads%nelem(1), &
                           ads%Uy, ads%p(2), ads%n(2), ads%nelem(2), &
                           ads%Uz, ads%p(3), ads%n(3), ads%nelem(3), &
                           ! solution, GnuPlotOutput, params, 0, MPI_COMM_WORLD)
                           solution, VtkOutput, params, 0, MPI_COMM_WORLD)

end subroutine PrintSolution

end module solution_output
