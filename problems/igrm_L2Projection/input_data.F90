!------------------------------------------------------------------------------
!
! MODULE: input_data
!
! DESCRIPTION:
!> @file input_data.F90
!> @brief Input parameters for the iGRM L2 projection example.
!>
!> @details
!> The module stores the command-line grid dimensions, test/trial
!> polynomial orders, MPI process-grid dimensions, and selected iGRM time
!> scheme used by the iGRM L2 projection driver. The actual forcing
!> callback is defined in \ref RHS_fun for this example.
!
!------------------------------------------------------------------------------
module input_data

   use argument_parser, ONLY: ARG_TIME_SCHEME_BE => TIME_SCHEME_BE, &
                              ARG_TIME_SCHEME_DG => TIME_SCHEME_DG, &
                              ARG_TIME_SCHEME_PR => TIME_SCHEME_PR, &
                              ReadIntegerArgument, ReadRealArgument, &
                              ReadTimeSchemeArgument, SelectTimeScheme

   implicit none

!> @brief Numbers of elements in the three parametric directions.
   integer(kind = 4), dimension(3) :: nelem

!> @brief Polynomial orders for the enriched iGRM test space.
   integer(kind = 4), dimension(3) :: p_test

!> @brief Polynomial orders for the trial space.
   integer(kind = 4), dimension(3) :: p_trial

!> @brief Numbers of MPI processes in the three process-grid directions.
   integer(kind = 4) :: procx, procy, procz

!> @brief Supported iGRM time-scheme identifiers.
   integer(kind = 4), parameter :: TIME_SCHEME_DG = ARG_TIME_SCHEME_DG
   integer(kind = 4), parameter :: TIME_SCHEME_PR = ARG_TIME_SCHEME_PR
   integer(kind = 4), parameter :: TIME_SCHEME_BE = ARG_TIME_SCHEME_BE

!> @brief Time-step length used to configure the selected iGRM scheme.
   real(kind = 8) :: scheme_tau = 1.d0

!> @brief Time scheme selected for the iGRM multi-step driver.
   character(len = 32) :: time_scheme = "dg"

!> @brief Normalized iGRM time-scheme identifier selected from arguments.
   integer(kind = 4) :: selected_time_scheme = TIME_SCHEME_DG


contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Reads iGRM L2 projection parameters from command-line arguments.
!>
!> @details
!> The accepted argument lists are:
!> `<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz>`,
!> `<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> <ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <scheme>`,
!> or
!> the same list followed by `<tau> <scheme>`.
!
!---------------------------------------------------------------------------
   subroutine InitializeParameters
      implicit none
      integer(kind = 4) :: nargs

      nargs = COMMAND_ARGUMENT_COUNT()

      if (nargs .NE. 12 .AND. nargs .NE. 13 .AND. nargs .NE. 14) then
         write(*,*) "proper usage with arguments: ", &
         "<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> ", &
         "<ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> [scheme]"
         write(*,*) "or: ", &
         "<nelem_x> <nelem_y> <nelem_z> <ptest_x> <ptest_y> <ptest_z> ", &
         "<ptrial_x> <ptrial_y> <ptrial_z> <procx> <procy> <procz> <tau> <scheme>"
         STOP 5
      end if

      call ReadIntegerArgument(1, nelem(1))
      call ReadIntegerArgument(2, nelem(2))
      call ReadIntegerArgument(3, nelem(3))
      call ReadIntegerArgument(4, p_test(1))
      call ReadIntegerArgument(5, p_test(2))
      call ReadIntegerArgument(6, p_test(3))
      call ReadIntegerArgument(7, p_trial(1))
      call ReadIntegerArgument(8, p_trial(2))
      call ReadIntegerArgument(9, p_trial(3))
      call ReadIntegerArgument(10, procx)
      call ReadIntegerArgument(11, procy)
      call ReadIntegerArgument(12, procz)
      scheme_tau = 1.d0
      call SelectTimeScheme("dg", selected_time_scheme, time_scheme)
      if (nargs .EQ. 14) then
         call ReadRealArgument(13, scheme_tau)
         call ReadTimeSchemeArgument(14, selected_time_scheme, time_scheme)
      else if (nargs .EQ. 13) then
         call ReadTimeSchemeArgument(13, selected_time_scheme, time_scheme)
      end if

   end subroutine InitializeParameters



end module input_data
