c CONDENSATION
c Based on Tzivion, Feingold, Levin, JAS 1989 and 
c Stevens, Feingold, Cotton, JAS 1996
c--------------------------------------------------------------------

c TAU(k) ......... Forcing for diffusion = (2/3)*CPT*ETA_BAR*DELTA_T
c X(K) ........ Array of bin limits in mass space
c AMKD(K,J) ... Input array of mass moments
c ANKD(K) ..... Input array of number moments
c AMK(K,J) .... Output array of mass moments
c ANK(K) ...... Output array of number moments
c CSPECIES .... Index of chemical species that is condensing
c
c The supersaturation is calculated outside of the routine and assumed
c to be constant at its average value over the timestep.
c 
c The method has three basic components:
c (1) first a top hat representation of the distribution is construced
c     in each bin and these are translated according to the analytic
c     solutions
c (2) The translated tophats are then remapped to bins.  Here if a 
c     top hat entirely or in part lies below the lowest bin it is 
c     not counted.
c 

Cpja Additional notes (Peter Adams)

C     I have changed the routine to handle multicomponent aerosols.  The
C     arrays of mass moments are now two dimensional (size and species).
C     Only a single component (CSPECIES) is allowed to condense during
C     a given call to this routine.  Multicomponent condensation/evaporation
C     is accomplished via multiple calls.  Variables YLC and YUC are
C     similar to YL and YU except that they refer to the mass of the 
C     condensing species, rather than total aerosol mass.

C     I have removed ventilation variables (VSW/VNTF) from the subroutine
C     call.  They still exist internally within this subroutine, but
C     are initialized such that they do nothing.

C     I have created a new variable, AMKDRY, which is the total mass in
C     a size bin (sum of all chemical components excluding water).  I
C     have also created WR, which is the ratio of total wet mass to 
C     total dry mass in a size bin.

C     AMKC(k,j) is the total amount of mass after condensation of species
C     j in particles that BEGAN in bin k.  It is used as a diagnostic
C     for tracking down numerical errors.

Cpja End of my additional notes

      SUBROUTINE TMCOND(TAU,X,AMKD,ANKD,AMK,ANK,CSPECIES,moxd)

      IMPLICIT NONE

C      include '../headers/inc_tomas_arrays.f'
      include 'sizecode.COM' 

      INTEGER L,I,J,K,IMN,CSPECIES,kk
      double precision DN,DM,DYI,TAU(ibins),XL,XU,YL,YLC,YU,YUC
      double precision TEPS,NEPS,EX2,ZERO
      double precision XI,XX,XP,YM,WTH,W1,W2,WW,AVG
      double precision VSW,VNTF(ibins)
      double precision TAU_L, maxtau
      double precision X(ibins+1),AMKD(ibins,icomp),ANKD(ibins)
      double precision AMK(ibins,icomp),ANK(ibins)
      double precision AMKDRY(ibins), AMKWET(ibins), WR(ibins)
      double precision DMDT_INT
      double precision AMKD_tot
      PARAMETER (TEPS=1.0d-40,NEPS=1.0d-20)
      PARAMETER (EX2=2.d0/3.d0,ZERO=0.0d0)
      double precision moxd(ibins) !moxid/Nact (win, 5/25/06)
      double precision c1, c2 !correction factor (win, 5/25/06)
      double precision xk_hi,tmpvar,xk_lo,frac_lo_n,frac_lo_m

      external DMDT_INT

      DOUBLE PRECISION QQQ,AAA,BBB,CCC
      DOUBLE PRECISION RATIO
      DOUBLE PRECISION EFF,PSI
      DOUBLE PRECISION ATOT,A1,A2,HH
      DOUBLE PRECISION DIV

 3    format(I4,200E20.11)

C-----CODE-----------------------------------------------------------

      WW=0.5d0
      
C-----CALCULATE THE BOUND MULTIPLYING RATIO--------------------------

      RATIO = Xk(2)/Xk(1)

C If any ANKD are zero, set them to a small value to avoid division by zero
      do k=1,ibins
         if (ANKD(k) .lt. NEPS) then
            ANKD(k)=NEPS
            AMKD(k,SRTSO4)=NEPS*1.4*X(k) !make the added particles SO4
            do j=1,icomp
               if (j.ne.SRTSO4)then
                  AMKD(k,j)=0.d0
               endif
            enddo
         endif
      enddo

Cpja Sometimes, after repeated condensation calls, the average bin mass
Cpja can be just above the bin boundary - in that case, transfer a some
Cpja to the next highest bin
      do k=1,ibins-1
         
        AMKD_tot=0.0d0
        do kk=1,icomp-idiag
           AMKD_tot=AMKD_tot+AMKD(k,kk)
           
        enddo
C        print*,'AMKD_tot=',AMKD_tot,'bin=',k 
         if (AMKD_tot/ANKD(k).gt.X(k+1)) then
            do j=1,icomp
               AMKD(k+1,j)=AMKD(k+1,j)+0.1d0*AMKD(k,j)
               AMKD(k,j)=AMKD(k,j)*0.9d0
            enddo
            ANKD(k+1)=ANKD(k+1)+0.1d0*ANKD(k)
            ANKD(k)=ANKD(k)*0.9d0
         endif
      enddo

Cpja Initialize ventilation variables so they don't do anything
      VSW=0.0d0
      DO L=1,ibins
         VNTF(L)=0.0d0
      ENDDO

Cpja Initialize AMKDRY and WR
      DO L=1,ibins
         AMKDRY(L)=0.0d0
         AMKWET(L)=0.0d0
         DO J=1,icomp-idiag
            AMKDRY(L)=AMKDRY(L)+AMKD(L,J)
         ENDDO
         DO J=1,icomp
            AMKWET(L)=AMKWET(L)+AMKD(L,J)
         ENDDO
         WR(L)=AMKWET(L)/AMKDRY(L)
      ENDDO

Cpja Initialize X() array of particle masses based on xk()
      DO L=1,ibins
         X(L)=xk(L) ! These are the the lower limits of each bin
      ENDDO

C-----ONLY SOLVE WHEN SIGNIFICANT FORCING IS AVAILABLE---------------
      
      MAXTAU = 0.0 
      DO i=1,ibins
         MAXTAU = MAX(MAXTAU,ABS(TAU(i)))
      END DO

      
      IF (ABS(MAXTAU).LT.TEPS) THEN
         
         DO i = 1,iBINS
            DO j = 1,iCOMP
               AMK(i,j) = AMKD(i,j)
            ENDDO
            ANK(i) = ANKD(i)
         ENDDO
         
         GOTO 300
      END IF

         

      DO L=1,ibins
         DO J=1,icomp
            AMK(L,J)=0.d0
         ENDDO
         ANK(L)=0.d0
      ENDDO


c identify tophats and do lagrangian growth
c
      DO L = 1,iBINS

C--------SKIP IF NO PARTICLES IN BIN---------------------------------         
         
         IF (ANKD(L).EQ.0.) GOTO 200

C--------SKIP IF NO COND./EVAP. FORCING IN BIN-----------------------
         
         IF (TAU(L).EQ.0.) THEN
            ANK(L) = ANK(L) + ANKD(L)
            DO j = 1,iCOMP
               AMK(L,j) = AMK(L,j) + AMKD(L,j)
            ENDDO
            GOTO 200
         ENDIF

C--------CALCULATE AVE. MASS AND LINEAR APPROX. COEFFICIENTS---------
         
         AVG = AMKDRY(L)/ANKD(L)
         DIV = AVG/X(L)
         
         IF (DIV.GT.RATIO) DIV = RATIO
         IF (DIV.LT.1.0  ) DIV = 1.0
         
         EFF = 2.0*ANKD(L)/X(L)*(RATIO-DIV)*(1.0/(RATIO-1.0)**2.0)
         PSI = 2.0*ANKD(L)/X(L)*(DIV  -1.0)*(1.0/(RATIO-1.0)**2.0)
          
C         print*,'------------------------------------------'
C         print*,'AVG=',AVG
C         print*,'DIV=',DIV
C         print*,'EFF=',EFF
C         print*,'PSI=',PSI
C         print*,'------------------------------------------'

C--------DETERMINE THE INITIAL GROWTH WINDOW-------------------------
         
         XX=X(L)/AVG

         XI=.5d0 + XX*(1.5d0 - XX)
         
         if (XI .LT. 1.d0) then
            !W1 will have sqrt of negative number
            write(*,*)'ERROR: tmcond - XI<1 for bin: ',L
            write(*,*)'lower limit is',X(L)
            write(*,*)'AVG is ',AVG
            write(*,*)'Nk is ', ANKD(L)
            write(*,*)'Mk are ', (AMKD(k,j),j=1,icomp)
            write(*,*)'Initial N and M are: ',ANKD(L),AMKDRY(L)
            STOP
         endif

         W1 = SQRT(12.d0*(XI-1.d0))*AVG
         W2 = MIN(X(L+1)-AVG,AVG-X(L))
         WTH=W1*WW+W2*(1.d0-WW)

         IF(WTH.GT.1.) then
            write(*,*)'WTH>1 in cond.f, bin #',L
            STOP
         ENDIF

C         XU=AVG+WTH*.25d0
C         XL=AVG-WTH*.25d0

         XU = X(L+1)
         XL = X(L)

C--------GROW THE PARTICLES------------------------------------------
         
         TAU_L = TAU(L)*MAX(1.d0,VNTF(L)*VSW)

C         print*,'TAU_L=',TAU_L,'XU=',XU,'XL=',XL

         IF(TAU_L/TAU(l).GT. 6.) THEN
            PRINT *,'TAU..>6.',TAU(l),TAU_L,VSW,L
         ENDIF
         IF(TAU_L.GT.TAU(l)) THEN 
            PRINT *,'TAU...',TAU(l),TAU_L,VSW,L
         ENDIF
         
C         print*,'DMDT_INT(XU,TAU_L,WR(L))=',DMDT_INT(XU,TAU_L,WR(L))
C         print*,'DMDT_INT(XL,TAU_L,WR(L))=',DMDT_INT(XL,TAU_L,WR(L))
         YU = DMDT_INT(XU,TAU_L,WR(L))
         YL = DMDT_INT(XL,TAU_L,WR(L))
         
         if(moxd(L).eq.0d0) then
            c1=1.d0          !for so4cond call, without correction factor.
         else
            c1 = moxd(L)*2.d0/(YU+YL-XU-XL)
         endif

         c2 = c1 - (c1-1.d0)*(XU+XL)/(YU+YL)
         YU = YU*c2
         YL = YL*c2

         YUC = XU*AMKD(L,CSPECIES)/AMKDRY(L) + YU - XU
         YLC = XL*AMKD(L,CSPECIES)/AMKDRY(L) + YL - XL

         IF (YU .GT. X(ibins+1) ) THEN
            YUC = YUC*X(ibins+1)/YU
            YU  = X(ibins+1)
         ENDIF

C--------CALCULATE TOTAL TRAPEZOIDAL AREA UNDER CURVE----------------
         
         ATOT = (EFF + PSI)*(YU - YL)*0.5

C         print*,'ATOT=',ATOT
C         print*,'EFF=',EFF
C         print*,'PSI=',PSI
C         print*,'YU=',YU
C         print*,'YL=',YL
C--------REMAP BIN NUMBER AND MASS-----------------------------------
         
         DYI=1.d0/(YU-YL)
         
C--------IF UPPER HAT IS LESS THAN THE LOWEST BOUND------------------
         
         IF (YU.LT.X(1)) GOTO 200
         
C--------IF LOWER HAT IS LESS THAN THE LOWEST BOUND------------------
         
         IF (YL.LT.X(1)) THEN
            
            ANK(L) = ANK(L) + ANKD(L)
            
            do j=1,icomp
               if (J.EQ.CSPECIES) then
                  AMK(L,J)=AMK(L,J) + (YUC+YLC)*.5d0 *ANKD(L) ! Mass * number of particles
               else
                  AMK(L,J)=AMK(L,J)+AMKD(L,J)
               endif
            enddo
            
            GOTO 200
            
         ENDIF
         
         
c Begin remapping (start search at present location if condensation)
         
         IMN = 1
         
         IF(TAU(l).GT.0.)IMN=L
         
         ! 0.5*(YUC+YLC) seems to be the mass of the particle after condensation/evaporation
         
         DO I = IMN,ibins

C-----------FIND WHICH BIN THE LOWER HAT IS IN-----------------------

            IF (YL.LT.X(I+1)) THEN

C--------------IF THE UPPER HAT IS ALSO IN THIS BIN------------------

               IF (YU.LE.X(I+1)) THEN
                  DN=ANKD(L)

                  do j=1,icomp
                     DM=AMKD(L,J)

                     IF (J.EQ.CSPECIES) THEN
                        if ((YUC+YLC)*.5d0*DN .lt. 0.d0) then
                           !print*,'S O A C O N D T E S T  1!!!!'
                           !print*,'L=',L,'I=',I
                           !print*,'DN=',DN
                           !print*,'YLC=',YLC
                           !print*,'YUC=',YUC
                           !print*,'val=',(YUC+YLC)*.5d0*DN
                           AMK(I,J)=AMK(I,J)
                        else
                           AMK(I,J)=(YUC+YLC)*.5d0*DN+AMK(I,J)
                        endif
                     ELSE
                        AMK(I,J)=AMK(I,J)+DM
                     ENDIF

                  enddo


                  ANK(I) = ANK(I) + DN

C--------------IF THE UPPER HAT IS IN THE NEXT BIN-------------------

               ELSE

C-----------------CALCULATE THE TRAPEZOIDAL AREA IN THE CURRRENT BIN-                  
                  
                  HH = EFF + (PSI - EFF)*(X(I+1) - YL)/(YU - YL)
                  A1 = (HH + EFF)*(X(I+1)-YL)*0.5

C                  DN=ANKD(L)*(X(I+1)-YL)*DYI  ! THE AMOUNT TO STAY IN THE CURRENT BIN
                  
                  DN=ANKD(L)*A1/ATOT

                  do j=1,icomp
C                     DM=AMKD(L,J)*(X(I+1)-YL)*DYI
                     DM=AMKD(L,J)*A1/ATOT
                     
                     IF (J.EQ.CSPECIES) THEN
                        XP=DMDT_INT(X(I+1),-1.0d0*TAU_L,WR(L))
                        YM=XP*AMKD(L,J)/AMKDRY(L)+X(I+1)-XP
                        if ((YM+YLC)*.5d0*DN .lt. 0.d0) then
                           !print*,'S O A C O N D T E S T  2!!!!'
                           !print*,'L=',L,'I=',I
                           !print*,'DN=',DN
                           !print*,'YM=',YM
                           !print*,'YLC=',YLC
                           !print*,'YUC=',YUC
                           !print*,'val=',(YM+YLC)*.5d0*DN
                           AMK(I,J)=AMK(I,J)
                        else
                           AMK(I,J)=DN*(YM+YLC)*0.5d0+AMK(I,J)
                        endif
                     ELSE
                        AMK(I,J)=AMK(I,J)+DM
                     ENDIF
                  enddo


                  ANK(I) = ANK(I) + DN

                  DO K=I+1,ibins
                     IF(YU.LE.X(K+1))GOTO 100 ! This section is almost always skipped
                     DN=ANKD(L)*(X(K+1)-X(K))*DYI ! THE AMOUNT TO ADD TO THE NEXT BINS
                     
                     do j=1,icomp
                        DM=AMKD(L,J)*(X(K+1)-X(K))*DYI
                        IF (J.EQ.CSPECIES) THEN
                           XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                           YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                           AMK(K,J)=DN*1.5d0*YM+AMK(K,J)
                        ELSE
                           AMK(K,J)=AMK(K,J)+DM
                        ENDIF
                     enddo
                     ANK(K)=ANK(K)+DN
                  ENDDO
                  STOP 'Trying to put stuff in bin ibins+1'
100               CONTINUE
                  

                  A2 = ATOT - A1


                  !DN=ANKD(L)*(YU-X(K))*DYI ! THE AMOUNT TO ADD TO THE NEXT BINS

                  DN=ANKD(L)*A2/ATOT ! THE AMOUNT TO ADD TO THE NEXT BINS
                  
                  
                  do j=1,icomp
C                     DM=AMKD(L,J)*(YU-X(K))*DYI
                     DM=AMKD(L,J)*A2/ATOT

                     IF (J.EQ.CSPECIES) THEN
                        XP=DMDT_INT(X(K),-1.0d0*TAU_L,WR(L))
                        YM=XP*AMKD(L,J)/AMKDRY(L)+X(K)-XP
                        AMK(K,J)=DN*(YUC+YM)*0.5d0+AMK(K,J)
                     ELSE
                        AMK(K,J)=AMK(K,J)+DM
                     ENDIF


                  enddo




                  ANK(K)=ANK(K)+DN
                  
               ENDIF  !YU.LE.X(I+1)

               GOTO 200
               
            ENDIF   !YL.LT.X(I+1)
         ENDDO !I loop

200   CONTINUE  ! START NEXT BIN

      ENDDO                  !L loop
C      STOP
     
 300  RETURN
      END
