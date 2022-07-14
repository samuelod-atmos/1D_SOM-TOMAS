
C     **************************************************
C     *  loginit                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, November 1999

C     Initializes aerosol number and mass distributions with sulfate
C     only and a lognormal size distribution.

C-----INPUTS------------------------------------------------------------

C     No - total number concentration (#/cm3)
C     Dp - lognormal mean diameter (microns)
C     sigma - width parameter

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE loginit(N1, Dp1, sigma1, N2, Dp2, sigma2)

      IMPLICIT NONE

C     INCLUDE FILES

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision N1, Dp1, sigma1
      double precision N2, Dp2, sigma2
      integer idtcomp    !indicates chemical component to be added

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer b,n     !bin and species counters
      integer I
      double precision Dl, Dh   !lower and upper bounds of bin (microns)
      double precision Dk       !mean diameter of current bin (microns)
      double precision np       !number of particles in size bin in GCM cell
      double precision pi
      double precision No1
      double precision No2      

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

C Convert No from #/cm3 to #/box
      No1=N1*boxvol
      No2=N2*boxvol      

      !Loop over number of size bins
      DO N=1,IBINS

         !Calculate diameter of this size bin
         Dl=1.0e+6*((6.*xk(n))/(SEED_DENS*pi))**(1.0/3.0)
         Dh=1.0e+6*((6.*xk(n+1))/(SEED_DENS*pi))**(1.0/3.0)
         Dk=sqrt(Dl*Dh)
c         write(*,*) 'Dk(',n,')= ',Dk,' microns'
         print*,'Dl =',Dl,'Dh =',Dh,'Dk =',Dk
         !Calculate number concentration
         np=(No1/(sqrt(2.*pi)*Dk*log(sigma1))*exp(-( (log(Dk/Dp1))**2./
     &         (2.*(log(sigma1))**2.) )) * (Dh-Dl)) +
     &      (No2/(sqrt(2.*pi)*Dk*log(sigma2))*exp(-( (log(Dk/Dp2))**2./
     &         (2.*(log(sigma2))**2.) )) * (Dh-Dl))

         Nk(n) = Nk(n) + np
c         print*,'N=',Nk(n)

         !Calculate component mass
CCC      AliA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
         DO I=1,ICOMP
            Mk(n,I) = Mk(n,I) +
     &                np*sqrt(xk(n))*sqrt(xk(n+1))*P_FRAC(I)
         END DO
CCC      ////////////////////////////////////////////////////////////////// AliA

      END DO

      RETURN
      END
