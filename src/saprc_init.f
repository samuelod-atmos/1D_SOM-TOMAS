      SUBROUTINE SAPRC_INIT(SAPRCGC2, NSOMPREC2, PREC_ID2,  
     &                      NSPEMISS2, EMISS_SPNAME2, EMISS_IPPM2,  
     &                      SEED_FRAC2, DLVP2, MWT_ACTSP2,
     &                      SOM_1STINDX2, SO2_INDX2, SO4_INDX2, time2,
     &                      ORG_O2C,SAPRC_SPNAME,INDEX_VOC2,kOH2)
      
C ====================================================================================================
C     This subroutine will get the input variables and read saprc files and intialize the gas phase
C     concentration for SOM-TOMAS (COSO)
C
C     Written by Ali Akherati - January 2019 - Colorado State University - MECH-E
C ====================================================================================================

      IMPLICIT NONE
      
C --- INCLUDE FILES ----------------------------------------------------------------------------------

      INCLUDE 'sizecode.COM'
      INCLUDE 'modlspc.h'
      INCLUDE 'common.inc'
      INCLUDE 'gaskin.h'

C --- VARIABLE LIST ----------------------------------------------------------------------------------
      
      INTEGER IERR
      INTEGER OFFSET
      INTEGER I, J
      INTEGER INDEX1, INDEX2
      INTEGER CNUMBER(IORG) ! Carbon number for organics
      INTEGER ONUMBER(IORG) ! Oxygen number for organics
      INTEGER NSOMPREC2      ! Number of SOM precursors
C      INTEGER POA_CNO_ST
C      INTEGER POA_CNO
C      INTEGER POA_1ST_LENNAME2   ! The length of the first POA name
C      INTEGER POA_1STINDX2
      INTEGER SOM_1STINDX2
      INTEGER SO2_INDX2
      INTEGER SO4_INDX2
      INTEGER NSPEMISS2 ! Number of species that has emissions
      
      REAL TEMP1
      REAL*4 SLA
      REAL EMISS_IPPM2(NSPEMISS2)
      REAL MWT_ACTSP2(NACT)
      REAL SAPRCGC2(NACT)

      DOUBLE PRECISION DLVP2(NSOMPREC2) ! DLVP for SOM precursors
      DOUBLE PRECISION DLVP_SP(IORG) ! DLVP for all organic species
      DOUBLE PRECISION MWORG_CB(IORG) ! Molecular weight of carbon back bone
      DOUBLE PRECISION SEED_FRAC2(NSPEMISS2) ! Fraction of species in initial seed aerosol
      DOUBLE PRECISION R ! gas constant
      DOUBLE PRECISION time2 ! seconds
      DOUBLE PRECISION ORG_O2C(iORG)
      DOUBLE PRECISION kOH2(6)

      CHARACTER*16 MODL ! SAPRC Model name
      CHARACTER*16 SP ! dummy variable for species names
C      CHARACTER*16 POA_ID ! POA name in SAPRC
C      CHARACTER*16 POA_1STNAME2  ! POA 1st name
      CHARACTER*16 SAPRC_SPNAME(NACT) ! SAPRC active species name
      CHARACTER*6 SOMSP_ID(IORG)
      CHARACTER*6 PREC_ID2(NSOMPREC2)
      CHARACTER*16 EMISS_SPNAME2(NSPEMISS2)
C      CHARACTER*16 POA_SPNAME
      CHARACTER*120 DATNAME, DIAGNAME
      CHARACTER*120 NCONCFILE, GCONCFILE, AECONCFILE, SPECFILE
      CHARACTER*120 SAPRCGCFILE, OH_FILE

      INTEGER INDEX_VOC2(6)
      PARAMETER(R=8.314)
      
      DATA  MODL / 'SAPRC14_REV1' /

      COMMON /FLNAMES/ DATNAME,   DIAGNAME,
     &                 NCONCFILE, GCONCFILE, AECONCFILE,
     &                 SPECFILE,  OH_FILE

C --- Species name in TOMAS --------------------------------------------------------------------------
      PARNAME(1) = 'SO4'
      PARNAME(ICOMP-1) = 'NH4'
      PARNAME(ICOMP) = 'H2O'
    
C --- Testing the shape
C ------------------------------------------------------------------------------------------------
c      print*,'Here are the shapes: SAPRCCG2, NSOMPREC2, PREC_ID2,
c     &                   EMISS_IPPM'
c      
c      print*,shape(SAPRCGC2),shape(NSOMPREC2),shape(PREC_ID2),
c     &                   shape(EMISS_IPPM2)
C ------------------------------------------------------------------------------------------------ 
C --- Main program to initialize SPARC and TOMAS gas phase -------------------------------------------

      ! zero emission
c      EMT(:) = 0.0
c      
c      ! SAPRC Initialization
c      OPEN(UNIT=1, FILE='saprc14_rev1.mod', STATUS='old')
c      OPEN(UNIT=3, FILE='saprc14_rev1.doc', STATUS='unknown')
c
c      print*, '44444'
c      
c      TREF = 298.0
c      TEMP1 = temp
c      CALL RDMOD(MODL, 1, 3, .true., TEMP1, 0, IERR,
c     +     NS, NACT, NSC, NSB, NDUM,
c     +     NRXN, LENKBF, NCOC1, NCOC2,
c     +     NCOEFV, NPHOTK, TREF, SLA,
c     +     DECLN, TSTART, CZMIN, MWT,
c     +     CONC0, CNO, NNO, SNO,
c     +     ONO, XNO, NRHC1, INDRHC,
c     +     RK, KPBUF, COEF, ABSQY,
c     +     DPHOTK, LKBUF, IPHRK, IRXTYP,
c     +     NAME, RXNLBL, COEFNM,
c     +     PHOTNM, TITLE) 
c      call flush(6) 
c      CLOSE(UNIT=1)
c      CLOSE(UNIT=3)

      print*, '55555'
      
C ######################################################################################################
C     DEFINE OFFSET BETWEEN ALL SPECIES AND ACTIVE SPECIES
      OFFSET = 13               ! ************ WARNING - CAN BE CHANGED BASED ON SAPRC FILE ***********
C ######################################################################################################

      print*, 'NACT =', NACT
      
C --- Initialize SAPRC species -------------------------------------------------------------------------
      SAPRC_SPNAME(:) = NAME(1+OFFSET:NACT+OFFSET) ! SAVE SAPRC gas phase Names
      SAPRCGC2(:) = 0.0
      
      print*, '66666'
      
C --- Initialize the POA species -----------------------------------------------------------------------
C      READ(POA_1STNAME2(1:5), '(A5)') POA_ID
C      READ(POA_1STNAME2(6:7), '(I2)') POA_CNO_ST

C      DO I=SRTORG1,SRTORG1+POA_NSP-1
C         POA_CNO = POA_CNO_ST+I-SRTORG1
C         WRITE(POA_SPNAME, '(A5,I2)') POA_ID, POA_CNO
C         PARNAME(I) = TRIM(POA_SPNAME)
C         CNUMBER(I-SRTORG1+1) = POA_CNO
C         ONUMBER(I-SRTORG1+1) = 0
C      END DO

      print*, '77777'
      
C      DO I=1,NACT
C         READ(SAPRC_SPNAME(I)(1:POA_1ST_LENNAME2), '(A7)') SP ! read the name of the species
C         IF (SP.EQ.POA_1STNAME2) THEN
C            POA_1STINDX2 = I ! Find the POA 1st name
C            GOTO 200
C         END IF
C      END DO
C 200  CONTINUE

      print*, '88888'

C --- Find SO2 index -----------------------------------------------------------------------------------
      DO I=1,NACT
         READ(SAPRC_SPNAME(I)(1:3), '(A3)') SP ! read the name of the species
         IF (SP.EQ.'SO2') THEN
            SO2_INDX2 = I ! Find the SO2 in active species
C            GOTO 201
         END IF
         IF(SP.EQ.'BEN') THEN
           INDEX_VOC2(1)=24
           kOH2(1) = 1.22E-12
         END IF
         IF(SP.EQ.'TRI') THEN
           INDEX_VOC2(2)=25
           kOH2(2) = 3.27e-11
         END IF
         IF(SP.EQ.'TOL') THEN
           INDEX_VOC2(3)=26
           kOH2(3) = 5.63e-12
         END IF
         IF(SP.EQ.'XYL') THEN
           INDEX_VOC2(4)=27
           kOH2(4) = 2.31E-11
         END IF
         IF(SP.EQ.'ISO') THEN
           INDEX_VOC2(5)=28
           kOH2(5) = 1.00e-10 
         END IF
         IF(SP.EQ.'TER') THEN
           INDEX_VOC2(6)=29
           kOH2(6) = 5.23e-11 
         END IF
      END DO
 201  CONTINUE

      print*,'INDEX_VOC2 =',INDEX_VOC2
      print*,'kOH2 =',kOH2   
   
C --- Find SO4 index -----------------------------------------------------------------------------------
      DO I=1,NACT
         READ(SAPRC_SPNAME(I)(1:4), '(A4)') SP ! read the name of the species
         IF (SP.EQ.'SULF') THEN
            SO4_INDX2 = I ! Find the SO4 in active species
            GOTO 202
         END IF
      END DO
 202  CONTINUE
      print*, 'SO4_indx2', SO4_indx2
C --- Initialize the oxygenated species (SOM Species) --------------------------------------------------
      DO I=1,NACT
         READ(SAPRC_SPNAME(I)(4:7), '(A7)') SP ! read the name of the species
         IF (SP.EQ.'SOMG') THEN
            SOM_1STINDX2 = I ! Find the 1st oxygenated SOM species
            GOTO 203
         END IF
      END DO
 203  CONTINUE

      print*, '99999'
      
      PARNAME(SRTORG1:SRTORG1+IORG-1-1) = ! Writing SOM species name in PARNAME
     &     SAPRC_SPNAME(SOM_1STINDX2:NACT)

      PARNAME(SRTORGLAST) = 'NONVOLATILE'
      
      SOMSP_ID(:) = PARNAME(SRTORG1:SRTORG1+IORG-1)(1:3) ! SOM species precursor ID
      
C      print*, 'POA_NSP, SRTORG1, IORG'
C      print*, POA_NSP, SRTORG1, IORG
      DO I=SRTORG1,IORG+SRTORG1-2
         print*, i, parname(i)
         READ(PARNAME(I)( 9:10), '(I2)') CNUMBER(I-SRTORG1+1)
         READ(PARNAME(I)(12:13), '(I2)') ONUMBER(I-SRTORG1+1)
      END DO
      CNUMBER(IORG) = 30
      ONUMBER(IORG) = 7

      DO I=1,iORG
        ORG_O2C(I) = DBLE(ONUMBER(I))/DBLE(CNUMBER(I))
      ENDDO
C      print*,'ORG_O2C =',ORG_O2C
C      print*,'ONUMBER =',ONUMBER
C      print*,'CNUMBER =',CNUMBER  

      print*, '10 10 10'
      
C --- Define dlvp for all the organic species ----------------------------------------------------------
      DO I=1,NSOMPREC2
         WHERE(SOMSP_ID==PREC_ID2(I))
            DLVP_SP = DLVP2(I)
         END WHERE
      END DO
      
      print*, '11 11 11'
      
C --- Calculate the molecular weight of organic species ------------------------------------------------
      MWORG = CNUMBER*12.0107 + ONUMBER*15.999
     &           + ((CNUMBER*2) + 2 - ONUMBER)*1.00794
      MWORG(IORG) = 1000.0
      
C --- Calculate the molecular weight for Cstar caluclation - MW of carbon back bone (CB) ---------------
      MWORG_CB = CNUMBER*12.0107 + ((CNUMBER*2)+2)*1.00794
      MWORG_CB(IORG) = 1000.0
      
C --- Calculate the Cstar values -----------------------------------------------------------------------
      CSTAR = 10.0**(-0.0337*MWORG_CB + 11.56 - ONUMBER*DLVP_SP) ! [ug m-3]
      CSTAR(IORG) = 1.0d-30 ! [ug m-3]
      
C *** Making POA non-volatile **************************************************************************
C      CSTAR(1:21) = 1.0d-20
C *** Making POA non-volatile **************************************************************************
      HVAP = -11.0*LOG10(CSTAR) + 131.0 ! HVAP - Epstein et al. 2010 ACP
      PSTAR = CSTAR*R*T1/MWORG/1.0e6 ! [pa]

      print*, '12 12 12'
      print*, 'NSPEMISS2=',NSPEMISS2 
C --- Implement SAPRC gas-phase conc., TOMAS gas-phase conc., and TOMAS particle-phase frac. -----------
      DO I=1,NSPEMISS2
         WHERE(SAPRC_SPNAME==EMISS_SPNAME2(I))
            SAPRCGC2 = EMISS_IPPM2(I) ! saprc gas-phase
         END WHERE
         WHERE(PARNAME==EMISS_SPNAME2(I))
            P_FRAC = SEED_FRAC2(I) ! particle-phase fraction of the given species
         END WHERE
      END DO

      print*, '13 13 13'

C      DO I=1,POA_NSP
C         GC(I+SRTORG1-1) = DBLE(APRCGC2(I+POA_1STINDX2-1))*1.0E-6 ! TOMAS gas-phase POA
C     &                     * PRES/R/TEMP * MWORG(I) *
C     &                     1.0E-3 * BOXVOL * 1.0e-6 ! ppm -> g/m3 -> kg/bag
C      END DO
      DO I=1,IORG-1
         GC(I+SRTORG1-1)=DBLE(SAPRCGC2(I-1+SOM_1STINDX2))*1.0E-6 ! TOMAS gas-phase SOM
     &                    * PRES/R/TEMP * MWORG(I) *
     &                    1.0E-3 * BOXVOL * 1.0e-6 ! ppm -> g/m3 -> kg/bag
      END DO
      GC(SRTORGLAST) = 0.0

      print*, '14 14 14'

!      print*,shape(SAPRCGC2),shape(NSOMPREC2),shape(SAPRC_SPNAME)
!      print*, 'name(i), SAPRCGC2(i)'
!      do i=1,NACT
!         print*, i , '-',name(i+offset),SAPRCGC2(i),SAPRC_SPNAME(i)
!         print*,i,'-',SAPRCGC2(i)
!      end do 
      
C --- Molecular weight of active species ---------------------------------------------------------------
      MWT_ACTSP2 = MWT(OFFSET+1:OFFSET+NACT)
      
C --- Save the specification ---------------------------------------------------------------------------
      OPEN(UNIT=37, FILE=SPECFILE, STATUS='old', ACCESS='APPEND')
 900  FORMAT(A20, 556A20)
 901  FORMAT(A20, 456A20)
 902  FORMAT(A20, 456I5)
 903  FORMAT(A20, 456F10.3)
 904  FORMAT(A20, 456F40.20)
 905  FORMAT(F8.1, E8.1)

      print*, '15 15 15'
      
C      WRITE(37,903) 'MWT of SOM:', MWT(STRINX+OFFSET:NACT+OFFSET)
      WRITE(37,900) 'ACTIVE GAS SP:', NAME(1+OFFSET:NACT+OFFSET)
      WRITE(37,900) 'SAPRC_SPNAME:', SAPRC_SPNAME(1:NACT)
      WRITE(37,901) 'SOM SP:', PARNAME(SRTORG1:SRTORG1+IORG-1)
      WRITE(37,902) 'SOM CARBON NO:', CNUMBER(:)
      WRITE(37,902) 'SOM OXYGEN NO:', ONUMBER(:)
      WRITE(37,903) 'MWORG of SOM:', MWORG(:)
      WRITE(37,903) 'MWORG_CB of SOM:', MWORG_CB(:)
      WRITE(37,904) 'CSTAR:', CSTAR(:)
      WRITE(37,904) 'HVAP:', HVAP(:)
      WRITE(37,904) 'PSTAR:', PSTAR(:)

      CLOSE(37)

      print*, '16 16 16'
      
C --- Copy rates of rxns to reactions array ------------------------------------------------------------
      DO I=1,MAXRXN
         RKZ(I) = RK(I)
      END DO
c      print*, 'RKZ=', RKZ
c      print*, 'CONST=', CONST

C --- Initilize the constant species concentration -----------------------------------------------------      
      NSPECS = NACT  !nact      
      DO I=1,6
         CONST(I) = CONC0(I)
      END DO
      CONST(6) = 0.0
C      CONST(6) = ( a_oh*exp(-1.0*ax_oh*0.0) +
C     &             b_oh*exp(-1.0*bx_oh*0.0) ) *
C     &           1.0E6/6.022E23/PRES*R*TEMP*1.0E6 ! OH Concentration - molecules/cm3 -> ppm

C --- Print the OH concentration -----------------------------------------------------------------------
C      OPEN(UNIT=39, FILE=OH_FILE, STATUS='old', ACCESS='append')
C      WRITE(39,905) time2, CONST(6)
C      CLOSE(39)

      print*, 'NSOMPREC=', NSOMPREC2
      print*, 'PREC_ID=', PREC_ID2
      print*, 'NSPEMISS=', NSPEMISS2
      print*, 'DLVP=', DLVP2
C      print*, 'POA_1STNAME=', POA_1STNAME2
C      print*, 'POA_1ST_LENNAME=', POA_1ST_LENNAME2
C      print*, 'POA_1STINDX=', POA_1STINDX2
      print*, 'SOM_1STINDX=', SOM_1STINDX2
      print*, 'time=', time2
      print*, 'ORG_O2C(1:10) init =',ORG_O2C(1:10)
C      print*, 'Name, Carbon, Oxygen, Cstar, MWorg, MWORG_CB'
C      do i=1,iorg
C         print*, PARNAME(i-1+srtorg1), CNUMBER(i), ONUMBER(i), CSTAR(i),
C     &        MWORG(i), MWORG_CB(i)
C      end do
     
      RETURN

      END SUBROUTINE
