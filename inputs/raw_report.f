    
C     **************************************************
C     *  report                                        *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Used by box model for testing aerosol microphysics code.  This
C     routine outputs the simulation status.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE report(time)

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

 1    format(A15,36E15.6)
 4    FORMAT(37E15.6)
 5    FORMAT(204E15.6)
C-----CODE--------------------------------------------------------------

C Time

      write(*,1) 'Time (s): ', time
      write(*,1) 'Time (h): ', time/3600.

C Bulk species

      DO N=1,ICOMP-1
         WRITE(*,1) PARNAME(N), Gc(N)
      END DO

C Aerosol number distribution

      WRITE(*,1) 'A#: ',(Nk(k),k=1,IBINS)

C Aerosol mass distribution

      DO J=1,ICOMP
         write(*,1) PARNAME(J), (Mk(k,J),k=1,IBINS)
      END DO

      OPEN(UNIT=34, FILE=NCONCFILE, STATUS='OLD', ACCESS='APPEND')
      WRITE(34,4) time/3600., (Nk(K),K=1,IBINS)
      CLOSE(34)
      
      OPEN(UNIT=35, FILE=GCONCFILE, STATUS='OLD', ACCESS='APPEND')
      WRITE(35,5) time/3600., (Gc(N),N=1,ICOMP-1)
      CLOSE(35)

      OPEN(UNIT=36, FILE=AECONCFILE, STATUS='OLD', ACCESS='APPEND')
      DO J=1,ICOMP
         WRITE(36,4) time/3600., (Mk(N,J), N=1,IBINS)
      END DO
      CLOSE(36)

      RETURN
      END
