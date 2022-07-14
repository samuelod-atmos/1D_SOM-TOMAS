
C     **************************************************
C     *  boxinit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, July 2000

C     Initializes the box model for testing size-resolved microphysics.
C     Starts either with zero concentration or with a restart file if
C     one exists.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE boxinit(time,No1,Dp1,sigma1,No2,Dp2,sigma2)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision time
      double precision No1, No2
      double precision Dp1, Dp2
      double precision sigma1, sigma2
      INTEGER INDSP

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,l,n

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

CAliA 2    format(200E15.6)

C-----CODE--------------------------------------------------------------

C Initializations

      !Initialize everything to zero at first
      DO N=1,IBINS
         Nk(N) = 0.0
         DO J=1,IORG
            Mk(N,J) = 0.0 
         END DO
      END DO

      call initbounds()
      call loginit(No1,Dp1,sigma1,No2,Dp2,sigma2)
C       call loginit_2mode(No,Dp,sigma,IDTSO4) !for 2 mode input

C initallize to ammonium sulfate
CCC      DO N=1,IBINS
c         Mk(N,strnh4) = Mk(N,strnh4)/96.0*18.0*2.0
CCC         Mk(N,strnh4) = 0.0
CC      END DO
      
      RETURN
      END
