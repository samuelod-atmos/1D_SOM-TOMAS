C     =======================================================================
C       THIS SUBROUTINE STEPS FORWARD FOR AUTO-OXIDATION OF THE PRECURSORS
C
C     AUTHOR: CHARLES HE
C     =======================================================================
      
      SUBROUTINE STEP_AUTOOX(f_HOMS,PRECS_kOH,
     &                       OH,Temp,Pres,BoxVol,timer,dt)

      IMPLICIT NONE

      INCLUDE '../headers/h_tomas_arrays.f'
      INCLUDE '../headers/h_som_arrays.f'
      
      DOUBLE PRECISION f_HOMS(nSOMPRECS)
      DOUBLE PRECISION PRECS_kOH(nSOMPRECS)
      DOUBLE PRECISION dt
      DOUBLE PRECISION OH
      DOUBLE PRECISION Temp
      DOUBLE PRECISION Pres
      DOUBLE PRECISION BoxVol
      DOUBLE PRECISION timer
            
      DOUBLE PRECISION NA
      PARAMETER(NA = 6.02E23)

      DOUBLE PRECISION R
      PARAMETER(R = 8.314)
      
      INTEGER i,j,k

      DOUBLE PRECISION PREC,dPREC
      INTEGER INDEX,xHOMS(nSOMCLASS)
      INTEGER INDEX1,INDEX2

      IF (timer.LT.0.0) GOTO 100
      
C     STEP FOR AUTO-OXIDATION:
C     =======================================================================
      k = 1
      
      DO i = 1,nSOMCLASS
         DO j = 1,uSOMPRECS(i)
C           PRECURSOR INDEX:
            INDEX = xPREC_FIRST(i)+j-1

C           PRECURSOR REACTED:
            PREC  = SOMGC(INDEX)
            dPREC = PRECS_kOH(k)*7.34e15/300.*PREC*OH*dt*f_HOMS(k)

C           DEDUCT FROM PRECURSOR:
            SOMGC(INDEX) = SOMGC(INDEX) - dPREC      
            k = k + 1
         END DO

C        ADD TO HOMS:         
         INDEX1 = xGRID_LAST(i)-xGRID_FIRST(1)+xORG_1ST
         INDEX2 = xGRID_LAST(i)-xGRID_FIRST(1)+1

         GC(INDEX1) = GC(INDEX1) + ! will have to substitute 
     &   dPREC*1e-6*(Pres/Temp/R)*ORG_MW(INDEX2)*0.001*(BoxVol*1e-6)

      END DO
      
 100  CONTINUE
      RETURN
      END SUBROUTINE
