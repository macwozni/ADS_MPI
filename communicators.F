      module communicators

      use parallelism, only :NRPROCX, NRPROCY, NRPROCZ

      integer :: processors(NRPROCX,NRPROCY,NRPROCZ)

c  ...group involving processors along X
      integer :: GROUPX(NRPROCY,NRPROCZ)
c  ...group involving processors along Y
      integer :: GROUPY(NRPROCX,NRPROCZ)
c  ...group involving processors along Z
      integer :: GROUPZ(NRPROCX,NRPROCY)
c  ...communicatorx along X
      integer :: COMMXALL(NRPROCY,NRPROCZ)
c  ...communicatorx along Y
      integer :: COMMYALL(NRPROCX,NRPROCZ)
c  ...communicatorx along Z
      integer :: COMMZALL(NRPROCX,NRPROCY)
c  ...local communicators
      integer :: COMMX,COMMY,COMMZ
c  ...corresponding contexts for SCALAPACK calls
      integer :: CONTEXTX
      integer :: CONTEXTY
      integer :: CONTEXTZ
      end
