!------------------------------------------------------------------------------
!
! MODULE: gnuplot
!
! DESCRIPTION:
!> @file gnuplot.F90
!> @brief Module providing routines for exporting scalar fields in a
!> plain-text format suitable for gnuplot.
!>
!> @details
!> This module groups output procedures used to serialize sampled scalar
!> data into text files that can be post-processed by gnuplot or similar
!> plotting tools.
!>
!> The provided functionality includes:
!> - export of a full three-dimensional stack of layers through
!>   \ref GnuPlotOutput,
!> - construction of per-layer filenames through \ref BuildFileName,
!> - writing of a single two-dimensional layer through \ref OutputLayer.
!>
!> In the current design, each \f$z\f$-layer of a three-dimensional field
!> is written to a separate `.plot` file. The resulting files contain
!> triples of the form \f$(x, y, value)\f$ arranged in a grid-like
!> plain-text layout.
!
!------------------------------------------------------------------------------
module gnuplot

   implicit none

contains

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Writes a three-dimensional scalar field as a sequence of
!> gnuplot-compatible layer files.
!>
!> @details
!> This routine traverses the scalar field \p vals along the third grid
!> direction and writes each two-dimensional slice to a separate output
!> file by calling \ref OutputLayer.
!>
!> The generated filenames follow the pattern
!>
!> \f[
!> \texttt{filename}\{\texttt{layer}\}\texttt{.plot}
!> \f]
!>
!> where `layer` is the zero-based index of the exported slice.
!
! Input:
! ------
!> @param[in] filename
!> Base filename pattern used to construct per-layer output filenames.
!>
!> @param[in] vals
!> Three-dimensional scalar field sampled on a structured grid.
!>
!> @param[in] params
!> Plot-parameter structure defining the grid resolution.
!
! Notes:
! ------
!> @note
!> One output file is generated for each layer in the third direction.
!
!---------------------------------------------------------------------------
subroutine GnuPlotOutput(filename, vals, params)
   use plot, ONLY: PlotParams
   implicit none
   character(len=*), intent(in) :: filename
!> @brief Plot parameters defining the grid resolution.
   type(PlotParams), intent(in) :: params
!> @brief Scalar field values written layer by layer.
   real(kind=8), intent(in) :: vals(params%resx, params%resy, params%resz)
!> @brief Loop counter over layers in the third direction.
   integer :: z

#ifdef IPRINT
   write (*, *) 'Starting GNUPLOT output...'
#endif
   do z = 1, params%resz
      call OutputLayer(filename, z - 1, vals(:, :, z), params)
   end do
#ifdef IPRINT
   write (*, *) 'Done with output.'
#endif

end subroutine GnuPlotOutput

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Builds the filename associated with one output layer.
!>
!> @details
!> This routine constructs an output filename from a base pattern and a
!> layer index. The resulting filename has the form
!>
!> \f[
!> \texttt{pattern}\,\texttt{layer}\texttt{.plot}.
!> \f]
!>
!> The integer layer index is first written to a temporary character
!> buffer, left-adjusted, trimmed, and then concatenated with the trimmed
!> base pattern and the `.plot` extension.
!
! Input:
! ------
!> @param[in] pattern
!> Base filename pattern.
!>
!> @param[in] layer
!> Layer index appended to the base pattern.
!
! Output:
! -------
!> @param[out] filename
!> Constructed filename of the output layer file.
!
! Notes:
! ------
!> @note
!> The layer index is written using fixed-width integer formatting and
!> then left-adjusted and trimmed before concatenation.
!
!---------------------------------------------------------------------------
subroutine BuildFileName(pattern, layer, filename)
   implicit none
!> @brief Base filename pattern.
   character(len=*), intent(in) :: pattern
!> @brief Layer index appended to the filename.
   integer, intent(in) :: layer
!> @brief Constructed output filename.
   character(len=*), intent(out) :: filename
!> @brief Temporary buffer holding the textual layer index.
   character(len=32) :: buffer

   write (buffer, '(I32)') layer
   buffer = adjustl(buffer)
   filename = trim(pattern)//trim(buffer)//'.plot'

end subroutine BuildFileName

!---------------------------------------------------------------------------
!
! DESCRIPTION:
!> @brief Writes one two-dimensional layer of scalar data to a
!> gnuplot-compatible text file.
!>
!> @details
!> This routine exports a single layer of scalar data to a plain-text
!> file whose name is constructed by \ref BuildFileName.
!>
!> For each grid point, the output file contains one line with:
!> - the physical coordinate in the first direction,
!> - the physical coordinate in the second direction,
!> - the scalar value at that point.
!>
!> A blank line is inserted after each sweep in the first direction,
!> which is convenient for grid-based plotting in gnuplot.
!
! Input:
! ------
!> @param[in] pattern
!> Base filename pattern used to construct the output filename.
!>
!> @param[in] zlayer
!> Index of the exported layer.
!>
!> @param[in] vals
!> Two-dimensional scalar field corresponding to the selected layer.
!>
!> @param[in] params
!> Plot-parameter structure defining the grid resolution.
!
! Notes:
! ------
!> @note
!> The file is written in formatted sequential text mode.
!
!---------------------------------------------------------------------------
subroutine OutputLayer(pattern, zlayer, vals, params)
   use plot, ONLY: PlotParams
   implicit none
!> @brief Base filename pattern.
   character(len=*), intent(in) :: pattern
!> @brief Index of the layer being exported.
   integer, intent(in) :: zlayer
!> @brief Plot parameters defining the grid resolution.
   type(PlotParams), intent(in) :: params
!> @brief Two-dimensional scalar field of the selected layer.
   real(kind=8), intent(in) :: vals(params%resx, params%resy)

!> @brief Output filename.
   character(len=512) :: filename
!> @brief Output unit number used for file writing.
   integer :: outFile = 57 ! random value, Grothendieck's prime
!> @brief Loop counters over the two-dimensional grid.
   integer :: ix, iy
!> @brief Physical coordinates of one grid sample.
   real(kind=8) :: x, y

#ifdef IPRINT
   write (*, *) 'Layer', zlayer
#endif
   call BuildFileName(pattern, zlayer, filename)

   open (unit=outFile, file=trim(filename), &
         form='formatted', access='sequential', status='replace', &
         action='write')

   do ix = 1, params%resx
      if (params%resx > 1) then
         x = params%startx + (params%endx - params%startx)*dble(ix - 1)/dble(params%resx - 1)
      else
         x = params%startx
      end if
      do iy = 1, params%resy
         if (params%resy > 1) then
            y = params%starty + (params%endy - params%starty)*dble(iy - 1)/dble(params%resy - 1)
         else
            y = params%starty
         end if
         write (outFile, '(3(ES24.16,1X))') x, y, vals(ix, iy)
      end do
      write (outFile, *) ! blank line
   end do

   close (outFile)

end subroutine OutputLayer

end module gnuplot
