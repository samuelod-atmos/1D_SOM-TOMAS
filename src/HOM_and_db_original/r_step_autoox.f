C     =================================================================
C     This subroutine steps forward for auto-oxidation of the precursor
C     =================================================================
      
      SUBROUTINE STEP_AUTOOX(AutoFrac,OH,Temp,Pres,BoxVol,time,dt)

      IMPLICIT NONE

      INCLUDE 'sizecode.COM'
      INCLUDE 'modlspc.h'
      INCLUDE 'common.inc'
      INCLUDE 'gaskin.h'

C      INCLUDE '../headers/inc_tomas_arrays.f'
C      INCLUDE '../headers/inc_saprc_arrays.f'
      
      DOUBLE PRECISION AutoFrac
      DOUBLE PRECISION dt
      DOUBLE PRECISION OH
C      DOUBLE PRECISION Temp
C      DOUBLE PRECISION Pres
C      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION time
      
      DOUBLE PRECISION PREC1,PREC2
      DOUBLE PRECISION dPREC1,dPREC2
      
      DOUBLE PRECISION NA
      PARAMETER(NA = 6.02E23)

      DOUBLE PRECISION R
      PARAMETER(R = 8.314)
      
      DOUBLE PRECISION nPRECS
      PARAMETER(nPRECS = 2)
      
      DOUBLE PRECISION kOH1
      PARAMETER(kOH1 = 1.63e-11)

      DOUBLE PRECISION kOH2
      PARAMETER(kOH2 = 1.32e-10)

      DOUBLE PRECISION MW1
      PARAMETER(MW1 = 184.37)

      DOUBLE PRECISION MW2
      PARAMETER(MW2 = 96.13)

      IF (time.LT.900.0) GOTO 100
      
      PREC1 = SAPRCGC(xPREC_FIRST) ! find index and specify
C      PREC2 = SAPRCGC(xPREC_FIRST+1)

      dPREC1 = (kOH1*7.34e15/300.)*PREC1*OH*dt*AutoFrac
C      dPREC2 = (kOH2*7.34e15/300.)*PREC2*OH*dt*AutoFrac


      SAPRCGC(xPREC_FIRST) = SAPRCGC(xPREC_FIRST)-dPREC1
C      SAPRCGC(xPREC_FIRST+1) = SAPRCGC(xPREC_FIRST+1)-dPREC2

      GC(iCOMP-4) = GC(iCOMP-4) +
     &           dPREC1*1e-6*(Pres/Temp/R)*MW1*0.001*(BoxVol*1e-6)
C     &           dPREC2*1e-6*(Pres/Temp/R)*MW2*0.001*(BoxVol*1e-6) 



 100  CONTINUE
      
      RETURN
      END SUBROUTINE
