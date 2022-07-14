
C     **************************************************
C     *  waterso4                                      *
C     **************************************************

C     Adaptation of ezwatereqm used in size-resolved sulfate only sim
C     November, 2001
C     ezwatereqm WRITTEN BY Peter Adams, March 2000

C     This function uses the current RH to calculate how much water is 
C     in equilibrium with the sulfate.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      double precision FUNCTION waterso4(rhe,ki)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision rhe   !relative humidity (0-100 scale)
      DOUBLE PRECISION k,A,R,pi
      DOUBLE PRECISION Dl,Dh,Dk
      INTEGER ji, ki

C-----VARIABLE DECLARATIONS---------------------------------------------

C     VARIABLE COMMENTS...

C     waterso4 is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fit to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(R=8.314159)
      parameter(pi=3.14159)

C-----CODE--------------------------------------------------------------

      if (rhe .gt. 99.) rhe=99.
      if (rhe .lt. 1.) rhe=1.
C      print*,'rhe/100.=',rhe/100.
C      Dl=1.0e+6*((6.*xk(ki))/(SEED_DENS*pi))**(1.0/3.0)
C      Dh=1.0e+6*((6.*xk(ki+1))/(SEED_DENS*pi))**(1.0/3.0)
C      Dk=sqrt(Dl*Dh)/10000000.0
C      A = (4.0*0.072*18.015)/(R*temp*1000.0*1000.0)
C      print*,'A =',A
C      k = (4.0*(A**3))/(27.0*(Dk**3)*(log(rhe)**2))
      k = 0.5
C      print*,'k =',k
      waterso4 = (rhe/100.-rhe/100.*k-1.0)/(rhe/100.-1.0)
C      print*,'waterso4 =',waterso4

         !check for error
         if (waterso4 .gt. 30.) then
            write(*,*) 'ERROR in waterso4'
            write(*,*) rhe,waterso4
            STOP
         endif

      RETURN
      END
