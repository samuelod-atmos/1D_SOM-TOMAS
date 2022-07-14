    
C     **************************************************
C     *  report                                        *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Used by box model for testing aerosol microphysics code.  This
C     routine outputs the simulation status.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE report(time)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

CAliA      include 'BB192SM9.COM'
CAliA      include 'BT263box.COM'
      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision time

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer n,k,j

      character*120 runname, datname, diagname
      CHARACTER*120 NCONCFILE, GCONCFILE, AECONCFILE
      common /FLNAMES/ datname, diagname,
     &                 NCONCFILE, GCONCFILE, AECONCFILE

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

 1    format(A15,40E15.6)
 4    FORMAT(41E15.6)
 5    FORMAT(459E15.6)
C-----CODE--------------------------------------------------------------

C Time

      write(*,1) 'Time (s): ', time
      write(*,1) 'Time (h): ', time/3600.

C Bulk species

C      DO N=1,ICOMP-1
C         WRITE(*,1) PARNAME(N),Gc(N)
C      END DO

C Aerosol number distribution

C      WRITE(*,1) 'A#: ',(Nk(k),k=1,IBINS)

C      print*,'ICOMP=',ICOMP
C Aerosol mass distribution
C      print*,'shape of Mk =',SHAPE(Mk)
C      print*,'Mk =',Mk(:,457)
C      DO J=1,ICOMP
C         write(*,1) PARNAME(J), (Mk(k,J),k=1,IBINS)
C         print*,'j=',j
C      END DO

C      print*,'All good 2'

      OPEN(UNIT=34, FILE=NCONCFILE, STATUS='OLD', ACCESS='APPEND')
      WRITE(34,4) time/3600., (Nk(K)/boxvol,K=1,IBINS)
      CLOSE(34)

C      print*,'All good 3'
      
      OPEN(UNIT=35, FILE=GCONCFILE, STATUS='OLD', ACCESS='APPEND')
      WRITE(35,5) time/3600., (Gc(N),N=1,ICOMP-1)
      CLOSE(35)

C      print*,'All good 4'

      OPEN(UNIT=36, FILE=AECONCFILE, STATUS='OLD', ACCESS='APPEND')
      DO J=1,ICOMP
         WRITE(36,4) time/3600., (Mk(N,J)*1.0E9 / 
     &                boxvol * 1.0E6, N=1, IBINS)
      END DO
      CLOSE(36)

      RETURN
      END
