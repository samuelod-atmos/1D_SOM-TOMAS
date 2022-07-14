
C     **************************************************
C     *  waternacl                                     *
C     **************************************************

C     WRITTEN BY Peter Adams, November 2001

C     This function uses the current RH to calculate how much water is 
C     in equilibrium with the seasalt.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

C-----Comments - SamO --------------------------------------------------

C     SamO changed this to calculate the wet ratio for organic species
C     based on molecular weight and rh

C----------------------------------------------------------------------

      double precision FUNCTION waterorg(rhe,ji,kappa,ORG_O2C3)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------
      include 'sizecode.COM'
C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision rhe   !relative humidity (0-100 scale)
      DOUBLE PRECISION A,R,pi
      DOUBLE PRECISION Dl,Dh,Dk
      DOUBLE PRECISION ORG_O2C3,k,kappa
      INTEGER ji
C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C     waternacl is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate waternacl based on the current humidity.
C     The curve fit is based on ISORROPIA results for sodium sulfate
C     at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(R=8.314159)
      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

C      Dl=1.0e+6*((6.*xk(ki))/(SEED_DENS*pi))**(1.0/3.0)
C      Dh=1.0e+6*((6.*xk(ki+1))/(SEED_DENS*pi))**(1.0/3.0)
C      Dk=sqrt(Dl*Dh)/10000000.0
      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.
C      A = (4.0*0.072*18.015)/(R*temp*1000.0*1000.0)
C      k = (4.0*(A**3))/(27.0*(Dk**3)*(log(rhe)**2))
      k = 1000000.0*(0.018015/1000.0)*(1/(17-0.85*ORG_O2C3))* !Kuang et al. 2020 
     &    (ORG_O2C3+0.28)/(7-log10(CSTAR(ji))/1.7)
C      k = 0.15
      waterorg = (rhe/100.-rhe/100.*k-1.0)/(rhe/100.-1.0)
      kappa = k

         !check for error
         if (waterorg .gt. 90.) then
            write(*,*) 'ERROR in waterorg'
            write(*,*) rhe,waterorg,k
            STOP
         endif

      RETURN
      END
