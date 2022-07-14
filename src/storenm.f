
C     **************************************************
C     *  storenm                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000

C     Stores values of Nk and Mk into Nkd and Mkd for diagnostic
C     purposes.  Also do gas phase concentrations.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE storenm()

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer j,k

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      do j=1,icomp-1
         Gcd(j)=Gc(j)
      enddo
      do k=1,ibins
         Nkd(k)=Nk(k)
         do j=1,icomp
            Mkd(k,j)=Mk(k,j)
         enddo
      enddo

      RETURN
      END
