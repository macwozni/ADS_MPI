module communicators

use parallelism, only: NRPROCX, NRPROCY, NRPROCZ

integer :: processors(NRPROCX,NRPROCY,NRPROCZ)

! group involving processors along X, Y, Z
integer :: GROUPX(NRPROCY,NRPROCZ)
integer :: GROUPY(NRPROCX,NRPROCZ)
integer :: GROUPZ(NRPROCX,NRPROCY)

! communicatorx along X, Y, Z
integer :: COMMXALL(NRPROCY,NRPROCZ)
integer :: COMMYALL(NRPROCX,NRPROCZ)
integer :: COMMZALL(NRPROCX,NRPROCY)

! local communicators
integer :: COMMX,COMMY,COMMZ

! corresponding contexts for SCALAPACK calls
integer :: CONTEXTX
integer :: CONTEXTY
integer :: CONTEXTZ

end
