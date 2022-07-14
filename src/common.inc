c****************************************************************************
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
c written by: Mike Kleeman (Feb 1997)
c             California Institute of Technology
c             Pasadena CA 91125
c
c The purpose of this include file is to pull together all the remaining
c common blocks to one location so that the code is cleaner to read.
c
c****************************************************************************

      INTEGER NALIM, NGLIM
      PARAMETER (NALIM=0, NGLIM=306)

      INTEGER nd2
      parameter (nd2 = maxact)

      REAL CONST(MAXCON)
      REAL EMT(MAXACT)
      INTEGER NSPECS

      COMMON / PARAM / CONST, NSPECS
      COMMON / TRAJT / EMT

      COMMON / integ / flag1,flag4,print4,print1,print2,print3,
     +                 nstiff
C     ,icount1,icount2,icount3,icount4
      COMMON / othcon / PSSA, PSSA2
      COMMON / othcon2 / PSSA3

      REAL PSSA(MAXSTS), PSSA2(MAXSTS)
      REAL*8 PSSA3(MAXSTS)

      INTEGER nstiff
      
      logical flag1,flag4,print4,print1,print2,print3

