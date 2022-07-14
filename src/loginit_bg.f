
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

      SUBROUTINE loginit_bg(N_bg1, Dp_bg1, sigma_bg1,
     &     N_bg2, Dp_bg2, sigma_bg2, N_bg3, Dp_bg3, sigma_bg3,
     &     density_bg, orgfrac_bg2,
     &     boxvol_old, time_check)


      IMPLICIT NONE

C     INCLUDE FILES

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision N_bg1, Dp_bg1, sigma_bg1
      double precision N_bg2, Dp_bg2, sigma_bg2
      double precision N_bg3, Dp_bg3, sigma_bg3
      double precision density_bg
      double precision boxvol_old ! initial box size
      double precision time_check
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
      double precision No3
      DOUBLE PRECISION orgfrac_bg2

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

C Convert No from #/cm3 to #/box
      No1=N_bg1*boxvol
      No2=N_bg2*boxvol
      No3=N_bg3*boxvol

      !Loop over number of size bins
      DO N=1,IBINS

         !Calculate diameter of this size bin
         Dl=1.0e+6*((6.*xk(n))/(density_bg*pi))**(1.0/3.0)
         Dh=1.0e+6*((6.*xk(n+1))/(density_bg*pi))**(1.0/3.0)
         Dk=sqrt(Dl*Dh)
c         write(*,*) 'Dk(',n,')= ',Dk,' microns'

         !Calculate number concentration
         np=(No1/(sqrt(2.0*pi)*Dk*log(sigma_bg1))*
     &        exp(-((log(Dk/Dp_bg1))**2.0/
     &        (2.0*(log(sigma_bg1))**2.0))) * (Dh-Dl)) +
     &        (No2/(sqrt(2.*pi)*Dk*log(sigma_bg2))*
     &        exp(-((log(Dk/Dp_bg2))**2.0/
     &        (2.0*(log(sigma_bg2))**2.0))) * (Dh-Dl)) +
     &        (No3/(sqrt(2.*pi)*Dk*log(sigma_bg3))*
     &        exp(-((log(Dk/Dp_bg3))**2.0/
     &        (2.0*(log(sigma_bg3))**2.0))) * (Dh-Dl))

         Nk(N) = Nk(N) + np
C         print*,'N=',Nk(n)
C         print*,'sigma2',sigma2
         !Calculate component mass
         Mk(N,SRTORGLAST) = Mk(n,SRTORGLAST) +
     &        np*sqrt(xk(n))*sqrt(xk(n+1))*orgfrac_bg2
         Mk(N,SRTSO4) = Mk(n,SRTSO4) +
     &        np*sqrt(xk(n))*sqrt(xk(n+1))*(1-orgfrac_bg2)         

      END DO

C      IF (time_check>0) THEN
CC     Convert No from #/cm3 to #/box
C         No1=N_bg1*boxvol_old
C         No2=N_bg2*boxvol_old
C
C      !Loop over number of size bins
C         DO N=1,IBINS
C
C         !Calculate diameter of this size bin
C            Dl=1.0e+6*((6.*xk(n))/(density_bg*pi))**(1.0/3.0)
C            Dh=1.0e+6*((6.*xk(n+1))/(density_bg*pi))**(1.0/3.0)
C            Dk=sqrt(Dl*Dh)
Cc         write(*,*) 'Dk(',n,')= ',Dk,' microns'
C
C         !Calculate number concentration
C            np=(No1/(sqrt(2.0*pi)*Dk*log(sigma_bg1))*
C     &           exp(-((log(Dk/Dp_bg1))**2.0/
C     &           (2.0*(log(sigma_bg1))**2.0))) * (Dh-Dl)) +
C     &           (No2/(sqrt(2.*pi)*Dk*log(sigma_bg2))*
C     &           exp(-((log(Dk/Dp_bg2))**2.0/
C     &           (2.0*(log(sigma_bg2))**2.0))) * (Dh-Dl))
C            
C            Nk(N) = Nk(N) - np
C            print*,'N=', N ,Nk(N)
C
C         !Calculate component mass
C            Mk(N,SRTORGLAST) = Mk(n,SRTORGLAST) -
C     &           np*sqrt(xk(n))*sqrt(xk(n+1))*orgfrac_bg2
C            Mk(N,SRTSO4) = Mk(n,SRTORGLAST) -
C     &           np*sqrt(xk(n))*sqrt(xk(n+1))*(1-orgfrac_bg2)
C            
C         END DO
C      END IF
      
      RETURN
      END
