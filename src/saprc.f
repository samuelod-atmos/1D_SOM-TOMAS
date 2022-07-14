      SUBROUTINE RDMOD (MOD,IUNIT,IOUT,DOCOUT,TEMP,IRXP,IERR,
     +                  NS, NACT, NSC, NSB, NDUM, NRXN,
     +                  LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                  TREF, LAT, DEC, TSTART, CZMIN, 
     +                  MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                  NRHC1, INDRHC,
     +                  RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                  LKBUF, IPHRK, IRXTYP,
     +                  NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c*****************************************************************************
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
c  created by: William P. Carter
c
c  modified by: Mike Kleeman (June 1995)
c         Environmental Engineering Science 138-78
c         California Institute of Technology
c         Pasadena, California, 91125
c
c  The subroutine was modified so that it could be used to read more than
c  one kinetic mechanism per model run.  Variables formerly passed in common
c  blocks are now passed explicitly.
c
c  modified by: Mike Kleeman (Feb 1999)
c         Civil and Environmental Engineering
c         University of California
c         Davis, California, 95616
c
c  The subroutine "NEWPHK" was modified to become "CITPHK".  The section of
c  code to calculate the solar zenith angle "cosz" was removed.  It is assumed
c  that cosz is determined by a previous call to the subroutine "SOLAR" and
c  is passed into "CITPHK" as an input variable.
c
c*****************************************************************************
c
c  inputs:
c    mod          name of the kinetic model being read [character*16]
c    iunit        unit number of the model compilation output file
c    iout         unit number for the output file from this subroutine
c    docout       true if model listing is to be produced [logical]
c    temp         default temperature
c    irxp         unit number to read mechanistic and kinetic parameters
c                 from if different from model default.  If non zero,
c                 subroutine RDPRM is called for this purpose.  If zero,
c                 model default parameters are used.
c    lat       lattitude used for calculation of solar intensity [ ]
c    dec       declination used for calculation of solar intensity [ ]
c    tstart    difference between time and true solar time [ ]
c
c  outputs:
c   ierr         return code passed to calling program
c                 IERR=0:   SUCCESSFUL COMPLETION
c                 IERR<>0:  ERROR.  (IF ERROR IN RDPRM, RDPRM ERROR CODE
c                           (NEGATIVE) USED.  IF ERROR IN RDMOD, IERR=1)
c    ns        number of species in the kinetic model
c    nact      number of active species in the kinetic model
c    nsc       number of constant species in the kinetic model
c    nsb       number of buildup species in the kinetic model
c    ndum      number of dummy species in the kinetic model
c    nrxn      number of reactions in the kinetic model
c    lenkbf    max number of parameters describing rate constants
c    ncoc1     index number for first constant coefficient in coef array
c    ncoc2     index number for last constant coefficient in coef array
c    ncoefv    number of variable coefficients in the model
c    nphotk    number of photolysis absorption profiles in the model
c    tref      reference temperature [K]
c    czmin     minimum dark value for which light intensity is proportional
c              to cos(z)
c    mwt       molecular weight of species in the model [g/mole]
c    conc0     initial concentrations of species in the model
c    cno       number of carbons in each species
c    nno       number of nitrogens in each species
c    nrhc1     number of emitted reactive hydrocarbon species
c    indrhc    array of emitted rhc species index #'s used to get species #
c              from rhc # used in emissions input
c    rk        array of rate constants
c    kpbuf     array of kinetic parameters for calculate of rate constants
c    coef      array of variable and constant coefficient values
c    absqy     array of quantum yield data used to calculate photolysis rates
c    dphotk    array of "dark" photolysis rates
c    lkbuf     location of parameters for reaction in array of kinetic params
c    iphrk     index where rate constant is stored in rate constant array
c    irxtyp    kinetic parameter type index for reaction
c    name      array for species names
c    rxnlbl    array of reaction labels
c    coefnm    array for variable coefficient names
c    photnm    array of labels for photolysis reactions
c    title     model title

C      implicit real (a-h,o-z)
      IMPLICIT NONE
      
      CHARACTER*16 MOD
      LOGICAL DOCOUT

      INCLUDE 'modlspc.h'
C
C       SPECIFICATIONS FOR LOCAL VARIABLES AND ARRAYS
C
      CHARACTER CHARXN*78
      LOGICAL*1 NXTNEW
      CHARACTER*8 ACTIVE,CONSTN,BLDP,DUMY,FORM,FORMOK
      CHARACTER*12 STST,EMRG

      INTEGER    NS, NSC, NACT, NSB, NDUM
      INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
      INTEGER    SNO(MAXNS)
      INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
      INTEGER    IRXTYP(MAXRXN)
      INTEGER    IUNIT, IOUT, IRXP, IERR
      INTEGER    NCOC1, NCOC2, NCOEFV
      INTEGER    I, I1, J, N1, N2, NSS
      INTEGER    IND, IPHOT, IRXN, ISOLWL
      INTEGER    LOCB

      REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
      REAL*4     XNO(MAXNS), ONO(MAXNS)
      REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
      REAL*4     LAT, DEC, TSTART, CZMIN
      REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF

      REAL       TEMP, TEMPX, X
      
      CHARACTER*16 NAME(MAXNS)
      CHARACTER*6  RXNLBL(MAXRXN)
      CHARACTER*8  COEFNM(MAXCOV)
      CHARACTER*8  PHOTNM(MAXPHK)
      CHARACTER*80 TITLE
C
C       FOLLOWING ARE WAVELENGTHS FOR WHICH ABS*QY DATA ARE
C       STORED.  SOME (HIGHER) WAVELENGTH INTERVALS LUMPED
C       TOGETHER, AND SDIST SUMMED FOR THESE.
C
      REAL SOLWL(NSOLWL)
C
      DATA ACTIVE/'ACTIVE  '/,CONSTN/'CONSTANT'/,BLDP /'BUILDUP'/,
     &     STST/'STEADY STATE'/,DUMY/'DUMMY   '/,FORMOK/'UCD:REV9'/,
     &     EMRG/'EMITTED ROG '/

      DATA SOLWL /     0.295,   0.300,    0.305,    0.310,    0.315,
     &     0.320,   0.325,   0.330,    0.335,    0.340,    0.345,
     &     0.350,            0.360,              0.370,
     &     0.380,            0.390,              0.400,
     &     0.410,            0.420,    0.430,    0.440,    0.450,
     &     0.460,   0.470,   0.480,    0.490,    0.500/
C               (LAST INTERVAL TO 1.000)
C
C
C       READ FORMAT CHARACTERS, AND CHECK IF OK WITH THIS VERSION
C
      READ (IUNIT,'(A8)') FORM
      IF (FORM.NE.FORMOK) THEN
         WRITE (IOUT,7001) FORM,FORMOK
 7001    FORMAT (' WRONG MODEL FILE FORMAT.  FORMAT = ''',
     &        A8,''', SHOULD BE ''',A8,'''.')
         CALL FLUSH(IOUT)
         STOP 'WRONG MODEL FILE FORMAT'
      ENDIF
C
      READ (IUNIT,7002) TITLE,NS,NACT,NSC,NSB,NDUM,NRXN,NCOEFV,
     &     NCOC1,NCOC2,NPHOTK,LENKBF,ISOLWL,TEMPX,TREF
 7002 FORMAT (A80/12I6/2F10.2)
c        print*,'rdmod marker1'
      call flush(6)
C
C       SET DEFAULT TEMPERATURE IF NONE SET ALREADY
C
      IF (TEMP.LE.0.0) TEMP=TEMPX
C
C       CHECK IF NUMBER OF SOLAR WL'S COMPATIBLE.
C
      IF (NSOLWL.NE.ISOLWL) THEN
         WRITE (IOUT,*) 'NSOLWL= ',NSOLWL,' MODEL INPUT =',
     &        ISOLWL,' NOT THE SAME'
         STOP 'NSOLWL INCOMPATIBLE'
      ENDIF
C
C       ABORT IF ANY DIMENSIONS EXCEEDED BY THIS MODEL
C
      IF (NS.GT.MAXNS) GOTO 9900
      IF (NACT.GT.MAXACT) GOTO 9900
      IF (NSC.GT.MAXCON) GOTO 9900
      IF (NSB.GT.MAXBO) GOTO 9900
      IF (NDUM.GT.MAXDUM) GOTO 9900
      IF (NRXN.GT.MAXRXN) GOTO 9900
      IF (NCOEFV.GT.MAXCOV) GOTO 9900
      IF (NCOC2.GT.MAXCOE) GOTO 9900
      IF (LENKBF.GT.MAXKBF) GOTO 9900
      IF (NPHOTK.GT.MAXPHK) GOTO 9900
C
C       SPECIES NAMES DEFAULT MOLECULAR WT, NO. CARBONS AND NITROGENS
C       AND INITIAL CONCENCENTRATIONS FOR NON-STEADY-STATE SPECIES.
C       FOR STEADY STATE SPECIES, READ NAME ONLY.
C
      I1=NSC+NACT+NSB+NDUM
      NSS=NS-I1
      READ (IUNIT,'(A16,0PF7.2,F6.2,2I3,2F6.2,1PE11.3)') 
     &     (NAME(I),MWT(I),
     &     CNO(I),NNO(I),SNO(I),ONO(I),XNO(I),CONC0(I),I=1,I1)
      IF (NSS.NE.0) READ (IUNIT,'(5A16)') (NAME(I),I=I1+1,NS)
C
C       READ VARIABLE COEFFICIENT NAMES AND VALUES
C
      IF (NCOEFV.GT.0) READ (IUNIT,'((4(A8,1PE11.3,1X)))')
     &     (COEFNM(I),COEF(I),I=1,NCOEFV)
C
C       READ CONSTANT COEFFICIENT VALUES
C
      IF (NCOC2.GE.NCOC1) READ (IUNIT,'(6E12.3)')
     &     (COEF(I),I=NCOC1,NCOC2)
C
C       READ ABSORPTION COEFFICIENT AND QUANTUM YIELD DATA
C
      IF (NPHOTK.GT.0) THEN
         DO IPHOT=1,NPHOTK
            READ (IUNIT,'(A8,I5,7X,1P6E10.3/(8E10.3))')
     &           PHOTNM(IPHOT),IPHRK(IPHOT),(ABSQY(I,IPHOT),I=1,NSOLWL)
         END DO
      ENDIF
C
C       READ IN KINETIC PARMS DATA
C
      READ (IUNIT,'(1P8E10.3)') (KPBUF(I),I=1,LENKBF)
      READ (IUNIT,'(16I5)') (LKBUF(IRXN),IRXN=1,NRXN)
      READ (IUNIT,'(16I5)') (IRXTYP(IRXN),IRXN=1,NRXN)
C
C       READ REACTION LABELS
C
c        READ (IUNIT,'(10A8)') (RXNLBL(IRXN),IRXN=1,NRXN)
      READ (IUNIT,'(10(A6,2x))') (RXNLBL(IRXN),IRXN=1,NRXN)
C
C
C       AT THIS POINT, ALL DATA REQUIRED BY MODEL ARE READ IN.  INPUT
C       FILE STILL HAS REACTION LISTINGS, WHICH CAN BE READ IF
C       DOCUMENTATION OUTPUT IS REQUESTED.
C
C
C       INITIALIZE RATE CONSTANTS
C
C       INITIALIZE RATE CONSTANTS FOR REACTIONS WITH CONSTANT RATE
C       CONSTANTS.  OTHERS INITIALIZED TO 0.
C
c        print*,'rdmod marker2'
      call flush(6)
      DO IRXN=1,NRXN
         IF (IRXTYP(IRXN).EQ.3) THEN
            RK(IRXN)=KPBUF(LKBUF(IRXN))
         ELSE
            RK(IRXN)=0.0
         ENDIF
      END DO
C
C       INITIALIZE TEMPERATURE DEPENDENT RATE CONSTANTS FOR
C       DEFAULT TEMPERATURE OF MODEL.
C
      CALL NEWRK (TEMP,1.0e6,
     +            NS, NACT, NSC, NSB, NDUM, NRXN,
     +            LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +            TREF, LAT, DEC, TSTART, CZMIN, 
     +            MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +            NRHC1, INDRHC,
     +            RK, KPBUF, COEF, ABSQY, DPHOTK,
     +            LKBUF, IPHRK, IRXTYP,
     +            NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
C
C       INITIALIZE PARAMETERS FOR PHOTOLYSIS REACTIONS.  ALSO CALCULATES
C       Z=0 PHOTKS IF MODEL LISTING IS REQUIRED.  IT IS ASSUMED THAT
C       THE LATITUDE, DECLINATION, AND TIME OFFSET (LAT, DEC, AND TSTART)
C       HAVE ALREADY BEEN DEFINED.
C
      CALL PHKINI (DOCOUT,
     +             NS, NACT, NSC, NSB, NDUM, NRXN,
     +             LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +             TREF, LAT, DEC, TSTART, CZMIN, 
     +             MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +             NRHC1, INDRHC,
     +             RK, KPBUF, COEF, ABSQY, DPHOTK,
     +             LKBUF, IPHRK, IRXTYP,
     +             NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
C
C
C       CALL RDPRM TO INPUT MODIFIED PARAMETERS IF NON ZERO IRXP
C       SPECIFIED.  NONZERO IERR INDICATES AN ERROR OCCURRED, AND
C       NO FURTHER PROCESSING DONE.
C
      IF (IRXP.GT.0) THEN
C
         CALL RDPRM (IRXP,IOUT,IERR,
     +               NS, NACT, NSC, NSB, NDUM, NRXN,
     +               LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +               TREF, LAT, DEC, TSTART, CZMIN, 
     +               MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +               NRHC1, INDRHC,
     +               RK, KPBUF, COEF, ABSQY, DPHOTK,
     +               LKBUF, IPHRK, IRXTYP,
     +               NAME, RXNLBL, COEFNM, PHOTNM, TITLE)

         IF (IERR.LT.0) RETURN
C
      ELSE
C
C               RDPRM IS USED TO DEFINE INDEX NO'S FOR EMITTED RHC'S.
C               SINCE IT IS NOT CALLED, SET  NO. OF RHC'S TO ZERO
C
         NRHC1=0
C
      ENDIF
C
C
C       AT THIS POINT, INITIALIZATION IS COMPLETE.  RETURN IF NO MODEL
C       LISTING OUTPUT IS REQUIRED.
C
      IERR=0
      IF (.NOT.DOCOUT) RETURN
C
C
C       THE CODE FROM HERE ON IS TO PRODUCE LISTING OF MODEL.
C
C
C       START NEW PAGE, AND OUTPUT TITLE OF MODEL.
C
      WRITE (IOUT,1030) TITLE,MOD
 1030 FORMAT ('1',130('*')/32X,'L I S T I N G   O F   C H E M I C A '
     &   ,'L   M E C H A N I S M'//30X,A80//41X,'MECHANISM PREP FILE NA'
     &   ,'ME = ',A8/1X/1X,130('*')//1X,'LISTING OF SPECIES AND THEIR '
     &   ,'DEFAULT MOLECULAR WEIGHTS, NO OF CARBONS AND NITROGENS, AND'
     &   ,' INITIAL CONCENTRATIONS.')
C
C       LIST SPECIES NAMES AND DEFAULTS.
C
      N1=NSC+NSB+NDUM
      N2=N1+NACT
      IF (NSC.NE.0) WRITE (IOUT,1032) CONSTN,NSC,(I,NAME(I),MWT(I),
     &     CNO(I),NNO(I),CONC0(I),I=1,NSC)
      IF (NSB.NE.0) WRITE (IOUT,1032) BLDP,NSB,(I,NAME(I),MWT(I),
     &     CNO(I),NNO(I),CONC0(I),I=NSC+1,NSC+NSB)
      IF (NDUM.NE.0) WRITE (IOUT,1032) DUMY,NDUM,(I,NAME(I),MWT(I),
     &     CNO(I),NNO(I),CONC0(I),I=NSC+NSB+1,N1)
      WRITE (IOUT,1032) ACTIVE,NACT,(I,NAME(I),MWT(I),CNO(I),NNO(I),
     &     CONC0(I),I=N1+1,N2)
 1032 FORMAT (/1X,A8,' SPECIES (',I3,'):'//7X,'NO. NAME             '
     &   ,'    MWT    C''S   N''S    CONC(0)'/7X,'--- ----------------'
     &   ,'  -------  -----  ---  -----------'/(6X,I3,2X,A16,0PF9.2,
     &   F7.2,I4,1PG13.3))
C
      IF (N2.LT.NS) WRITE (IOUT,1033) STST,NS-N2,(I,NAME(I),I=N2+1,NS)
 1033 FORMAT (/1X,A12,' SPECIES (',I3,'):'//7X,'NO. NAME             '
     &   ,'   NO. NAME                NO. NAME                NO. NAME'
     &   ,'                NO. NAME'/7X,'--- ----------------    --- -'
     &   ,'---------------    --- ----------------    --- ------------'
     &   ,'----    --- ----------------'/(6X,I3,2X,A16,3X,I3,2X,A16
     &   ,3X,I3,2X,A16,3X,I3,2X,A16,3X,I3,2X,A16))
C
C       LIST ROG EMITTED SPECIES, IF DEFINED BY RDPRM
C
      IF (NRHC1.NE.0) THEN
         WRITE (IOUT,1033) EMRG,NRHC1,(I,NAME(INDRHC(I)),I=1,NRHC1)
      ENDIF
C
C       LIST VARIABLE COEFFICIENT NAMES AND VALUES
C
      IF (NCOEFV.NE.0) WRITE (IOUT,1034) NCOEFV,(COEFNM(I),COEF(I),
     &     I=1,NCOEFV)
 1034 FORMAT (/1X,I3,' VARIABLE COEFFICIENTS:'//(' ',5(4X,A8,
     &     ' =',1PG11.4)))
C
C       LIST CONSTANT COEFFICIENT VALUES
C
      IF (NCOC1.LE.NCOC2) WRITE (IOUT,1036) NCOC2-NCOC1+1,(COEF(I),
     &     I=NCOC1,NCOC2)
 1036 FORMAT (/1X,I3,' CONSTANT COEFFICIENTS:'//(1X,1P10G13.3))
C
C       LIST PHOTOLYSIS DATA
C
      IF (NPHOTK.GT.0) THEN
         WRITE (IOUT,1040) NPHOTK,SOLWL
 1040    FORMAT (/1X,I3,' TYPES OF ABSORBTION COEFFICIENTS-QUA'
     &        ,'NTUM YIELDS  READ IN.  VALUES GIVEN FOR FOLLOWING '
     &        ,'WAVELENGTHS (IN NM):'//(32X,8F11.3))
         WRITE (IOUT,1042)
 1042    FORMAT (/1X,'PHOTNAME   K( Z=0 )   K( MIN )      ABSO'
     &        ,'RPTION COEFFICIENTS * QUANTUM YIELDS (CM2, BASE E)')
C
         X=CZMIN/CZLOW
         DO IPHOT=1,NPHOTK
            WRITE (IOUT,1046) PHOTNM(IPHOT),RK(IPHRK(IPHOT)),
     &           X*DPHOTK(IPHOT),(ABSQY(I,IPHOT),I=1,NSOLWL)
         END DO
 1046    FORMAT (/1X,A8,1P2E11.3,3X,1P8E11.3/(34X,8E11.3))
      ENDIF
C
C
C       WRITE REACTIONS.  TURN PAGE AND WRITE HEADER.
C
      WRITE (IOUT,63) TITLE,TEMP,TREF
 63   FORMAT ('1',A80/'0LIST OF REACTIONS AND RATE CONSTANTS.'//
     &   7X,'THERMAL RATE CONSTNATS GIVEN FOR DEFAULT TEMPERATURE ='
     &   ,F8.2/7X,'PHOTOLYSIS RATE CONSTANTS GIVEN FOR SOLAR Z = 0.'/
     &   7X,'ARRHENIUS RATE EXPRESSION USED IS K = A * EXP (EA/0.001'
     &   ,'9872 * TEMP) * (TEMP/TREF)**B,  WHERE TREF = ',F8.2//1X,
     &   'NO.  LBL      K       (     A       EA     B   )  REACTION'
     &   /1X)
C
C       LOOP OVER REACTIONS.
C
      DO IRXN=1,NRXN
C
C       READ FIRST RECORD FOR LISTING FOR REACTION.  LISTING MAY TAKE
C       MORE THAN ONE RECORD.  'NXTNEW' TELLS WHETHER THE NEXT RECORD
C       IS FOR A NEW REACTION.  SUBSEQUENT RECORDS WILL BE READ LATER,
C       AFTER KINETIC PARMS OUTPUT.
C
         READ (IUNIT,'(L1,1X,A)') NXTNEW,CHARXN
C
C       DETERMINE TYPE OF KINETIC PARMS FOR THIS REACTION.  THIS WILL
C       AFFECT THE OUTPUT FORMAT.
C
         IND=IRXTYP(IRXN)
         LOCB=LKBUF(IRXN)
C
         IF (IND.EQ.0) THEN
C           SAME K AS ANOTHER REACTION
            WRITE (IOUT,1052) IRXN,RXNLBL(IRXN),RK(LOCB),
     &           RXNLBL(LOCB),CHARXN
 1052       FORMAT (1X,I4,2X,A6,1PE12.3,5X,'K SAME AS RXN ',A6,6X,A)
         ELSEIF (IND.EQ.1) THEN
C               NOT USED
            GOTO 8010
         ELSEIF (IND.EQ.2) THEN
C               NOT USED
            GOTO 8010
         ELSEIF (IND.EQ.3) THEN
C               CONSTANT K
            WRITE (IOUT,1054) IRXN,RXNLBL(IRXN),RK(IRXN),CHARXN
 1054       FORMAT (1X,I3,2X,A6,1PE12.3,29X,A)
         ELSEIF (IND.EQ.4) THEN
C               K CALCULATED USING ARRHENIUS EXPRESSION
            WRITE (IOUT,1056) IRXN,RXNLBL(IRXN),RK(IRXN),
     &           KPBUF(LOCB+1),KPBUF(LOCB+2),KPBUF(LOCB+3),CHARXN
 1056       FORMAT (1X,I3,2X,A6,1PE12.3,' (',E11.3,0PF6.2,F7.3,')',
     &           2X,A)
         ELSEIF (IND.EQ.5) THEN
C               K CALCULATED USING FALLOFF EXPRESSION
            WRITE (IOUT,1058) IRXN,RXNLBL(IRXN),RK(IRXN),
     &           (KPBUF(LOCB+J),J=0,7),CHARXN
 1058       FORMAT (1X,I3,2X,A6,1PE12.3,' 0(',E10.3,0PF6.2,F7.3,')'
     &           /23X,'A(',1PE10.3,0PF6.2,F7.3,')'/25X,'F =',F6.3,4X,
     &           'N =',F6.3,4X,A)
         ELSEIF (IND.EQ.6) THEN
C               NOT USED
            GOTO 8010
         ELSEIF (IND.EQ.7) THEN
C               PHOTOLYSIS REACTION
            WRITE (IOUT,1060) IRXN,RXNLBL(IRXN),RK(IPHRK(LOCB)),
     &           PHOTNM(LOCB),CHARXN
 1060       FORMAT (1X,I3,2X,A6,1PE12.3,8X,'PHOT=',A8,8X,A)
         ELSE
            GOTO 8010
         ENDIF
C
C       READ AND WRITE SUBSEQUENT RECORDS GIVING REACTION LISTING, IF
C       ANY
C
 1062    IF (.NOT.NXTNEW) THEN
            READ (IUNIT,'(L1,1X,A)') NXTNEW,CHARXN
            WRITE (IOUT,'(51X,A)') CHARXN
            GOTO 1062
         ENDIF
C
C       END LOOP OVER REACTIONS
C
      END DO
C
C       MODEL DOCUMENTATION DONE.  RESET PHOTOLYSIS RATE CONSTANTS
C       TO ZERO
C
      DO I=1,NPHOTK
         RK(IPHRK(I))=0.0
      END DO
C
      RETURN
C
C       ERROR SECTION
C
 8000 WRITE (IOUT,*) 'RDMOD:  ERROR IN OPENING MODEL DATA SET.'
      GOTO 8999
C
 8010 WRITE (IOUT,8011) IND,RXNLBL(IRXN)
 8011 FORMAT (' RDMOD:  MODEL COMPATIBILITY ERROR.'/' ILLEGAL IND =',
     &     I8,' FOR REACTION = ',A6)
C
 8999 IERR=1
      RETURN
C
 9900 WRITE (IOUT,9902) NS,MAXNS,NACT,MAXACT,NSC,MAXCON,NSB,MAXBO,
     &     NDUM,MAXDUM,NRXN,MAXRXN,NCOEFV,MAXCOV,NCOC2,MAXCOE,LENKBF,
     &     MAXKBF,NPHOTK,MAXPHK
      print *,NS,MAXNS,NACT,MAXACT,NSC,MAXCON,NSB,MAXBO,
     &     NDUM,MAXDUM,NRXN,MAXRXN,NCOEFV,MAXCOV,NCOC2,MAXCOE,LENKBF,
     &     MAXKBF,NPHOTK,MAXPHK
 9902 FORMAT (/' RDMOD:  PROGRAM ABORTED BECAUSE MECHANISM EXCEEDS M'
     &   ,'AXIMUM LIMITS FOR THIS VERSION'/'0GIVEN BELOW IS CURRENT AND'
     &   ,' MAXIMUM VALUES OF RELEVANT PARAMETERS.'
     &   /' TOTAL SPECIES:',2I6/' ACTIVE SPECIES:',2I6
     &   /' CONSTANT SPECIES:',2I6,/' BUILDUP SPECIES:',2I6
     &   /' DUMMY SPECIES:',2I6/' REACTIONS:',2I6
     &   /' VARIABLE COEFFICIENTS:',2I6
     &   /' HIGHEST INDEX IN COEF ARRAY (CONSTANT COEFS):',2I6
     &   /' KINETIC PARAMETER BUFFER LENGTH:',2I6
     &   /' PHOTOLYSIS REACTIONS:',2I6)
      call flush(iout)
      STOP 'PIG RDMPD:  MODEL EXCEEDS PROGRAM LIMITS.  SEE OUTPUT FILE.'
      END










C       FILE:  RDPRM.FTN                CHEMICAL MECH SUBROUTINE FOR
C                                       AIRSHED MODEL
C
C       LAST UPDATE:    W.P.L. CARTER  3/6/88
C
C
      SUBROUTINE RDPRM (IUNIT,IOUT,IERR,
     +                  NS, NACT, NSC, NSB, NDUM, NRXN,
     +                  LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                  TREF, LAT, DEC, TSTART, CZMIN, 
     +                  MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                  NRHC1, INDRHC,
     +                  RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                  LKBUF, IPHRK, IRXTYP,
     +                  NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
C
C       SUBROUTINE TO READ MECHANISTIC PARAMETERS FOR LUMPED SPECIES
C       IN THE MECHANISM IN THE FORMAT OUTPUT BY THE "LUMPCOMP" PROGRAM.
C       SUITABLE FOR "CUSTOMIZING" MECHANISTIC PARAMETERS FOR LUMPED
C       SPECIES FOR SPECIFIC APPLICATIONS, BASED ON COMPOSITIONS OF
C       EMITTED SPECIES.
C
C       ALSO USED TO DETERMINE ORDERING OF EMITTED RHC SPECIES, AND
C       TO DERIVE INDEX ARRAY GIVING INTERNAL SPECIES NUMGER FROM
C       THE ORDERING OF THE EMITTED SPECIES.
C
C       THIS IS OPTIONALLY CALLED BY RDMOD WHEN THE MECHANISM IS
C       INITIATED, OR COULD BE CALLED SEPERATELY.  HOWEVER, IT SHOULD
C       NOT BE CALLED AFTER THE SIMULATION BEGINS, SINCE THE MECHANISTIC
C       PARAMETERS SHOULD BE CONSTANT DURING A SIMULATION.
C
C       ARGUMENTS:
C
C       IUNIT = (INPUT)  INPUT UNIT FOR FILE WITH MECHANISTIC PARAMETERS.
C               THE PROGRAM ASSUMES THAT THE THIS FILE IS OPEN.  INPUT
C               FROM THIS UNIT IS READ UNTIL END-OF-FILE, OR AN "END" RECORD
C               IS INPUT.  EXPECTS DATA IN THE FORMAT OF MECHANISTIC
C               PARAMETERS, SUCH AS PRODUCED BY THE "LUMPCOMP" PROGRAM.
C               (**NOTE** IF CALLED BY RDMOD, THIS SHOULD NOT BE
C               THE SAME UNIT AS USED BY RDMOD READING OUTPUT OF
C               THE MECHANISM PREPARATION PROGRAM (THE .MOD FILES),
C               SINCE RDPRM IS CALLED BEFORE ALL THE DATA IN THE .MOD
C               FILE IS READ.)
C
C               INPUT FORMAT OF DATA ON UNIT IUNIT IS ASSUMED TO CONSIST OF
C               EITHER COMMENT RECORDS (WITH A "!" IN COLUMN 1) OR RECORDS
C               WITH THE FOLLOWING FORMAT:
C                       FOR MECHANISTIC COEFFICIENT PARAMETERS:
C                               COL 1-8:   PARAMETER NAME (A4,4X)
C                               COL 9-20:  PARAMETER VALUE (G12.0)
C                       FOR RATE CONSTANTS:
C                               COL 1-2:   "K("
C                               COL 3-6:   RATE CONSTANT LABEL
C                               COL 7:     ")"
C                               COL 8-20:  RATE CONSTANT (PPM-MIN UNITS)
C                                          (G13.0)
C                       FOR EMITTED SPECIES ORDER AND PARAMETERS:
C                               COL 1-4:   "RHC("
C                               COL 5-8:   FIRST 4 CHARACTERS OF SPECIES
C                                          NAME
C                               COL 9:     ")"
C                               COL 10-15:  NO. CARBONS (F6.0)
C                               COL 16-24: MOLECULAR WEIGHT (F9.2)
C
C       IOUT =   OUTPUT UNIT.  USED FOR ERROR MESSAGES ONLY.  MUST BE OPEN.
C
C       IERR =  STATUS INDICATOR:
C                0 = SUCCESSFUL CONCLUSION.  INPUT ON IUNIT READ UNTIL
C                    END OF FILE (BUT FILE IS STILL OPEN).
C                1 = SUCCESSFUL CONCLUSION.  INPUT ON IUNIT READ UNTIL
C                    AN "END" RECORD WAS ENCOUNTERED.
C               -1 = ERROR ENCOUNTERED.  LAST RECORD READ ON IUNIT WAS
C                    ONE THAT COULD NOT BE SUCCESSFULLY PROCESSSED.
C               -2 = RDPRM CALLED WITH ZERO OR NEGATIVE IUNIT
C
      IMPLICIT NONE
      
      INTEGER IUNIT,IOUT,IERR
C
C
C       SPECIFICATION FOR MODEL-DEPENDENT DATA, INCLUDING KINETIC
C       AND MECHANISTIC PARAMETERS UPDATED BY THIS SUBROUTINE,
C       AND THE ARRAY OF EMITTED SPECIES INDEX NO.S PRODUCED BY
C       THIS SUBROUTINE.
C
        INCLUDE 'modlspc.h'
C
C
C       LOCAL SPECIFICATIONS
C
        CHARACTER  INLIN*100, SUBSTR*100
        INTEGER    NS, NSC, NACT, NSB, NDUM
        INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
        INTEGER    SNO(MAXNS)
        INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
        INTEGER    IRXTYP(MAXRXN)
        INTEGER    NCOC1, NCOC2, NCOEFV
        INTEGER    I, L1, L2, N1, NBLANK, LOCB

        REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
        REAL*4     ONO(MAXNS), XNO(MAXNS)
        REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
        REAL*4     LAT, DEC, TSTART, CZMIN
        REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF

        REAL       A, B, EA
        
        CHARACTER*16 NAME(MAXNS)
        CHARACTER*6  RXNLBL(MAXRXN)
        CHARACTER*8  COEFNM(MAXCOV)
        CHARACTER*8  PHOTNM(MAXPHK)
        CHARACTER*80 TITLE
C
C
C       ABORT WITH NO ACTION IF ZERO OR NEGATIVE INPUT UNIT
C       SPECIFIED.
C
        IF (IUNIT.LE.0) THEN
                IERR=-2
                RETURN
        ENDIF
C
C       NUMBER OF NON STEADY-STATE SPECIES
C
        N1=NSC+NSB+NDUM+NACT
C
C
C       NUMBER OF EMITTED SPECIES SET (OR RESTET) TO ZERO
C
        NRHC1=0
C
C       READ AND PROCESS NEXT INPUT.
C
  100   READ (IUNIT,'(A)',END=900) INLIN
C
        IF (INLIN(1:1).EQ.'!'.OR.INLIN(1:1).EQ.' ') THEN
C
C               RECORDS WITH A "!" OR BLANK IN FIRST POSITION ARE IGNORED.
C
                GOTO 100
C
        ELSEIF (INLIN(1:4).EQ.'END') THEN
C
C               END ... END OF INPUT.  DONE WITH PROCESSING.
C
                IERR=1
                RETURN
C
        ELSEIF (INLIN(1:2).EQ.'K(') THEN
C
C               K(rxnlbl) ... UPDATES RATE CONSTANTS FOR REACTION
C                             LABELED rxnlbl
C
                L1 = INDEX(INLIN,' ')-1
                IF (INLIN(L1:L1).NE.')') THEN
                        WRITE (IOUT,*) 'RDPRM:  ILLEGAL RATE CONSTANT',
     &                   ' SPECIFICATION INPUT: ',INLIN(1:L1)
                        IERR=-1
                        RETURN
                ENDIF
C
C               FIND WHICH REACTION THIS IS, THEN UPDATE RATE CONSTANT
C               ACCORDINGLY.  THE KINETIC PARAMETER TYPE FOR THIS REACTION
C               IS CHANGED TO "CONSTANT", SO IT WON'T BE OVERWRITTEN
C               BY NEWRK.  NOTE THAT RATE CONSTANTS MODIFIED THROUGH
C               THIS INPUT BECOME TEMPERATURE INDEPENDENT.  "SAMEK" AND
C               PHOTOLYSIS TYPE OF RATE CONSTANTS CANNOT BE UPDATED
C               IN THIS WAY.
C
                DO 110 I=1,NRXN
                L2 = NBLANK(RXNLBL(I))
                IF (RXNLBL(I)(:L2).EQ.INLIN(3:L1-1)) THEN
                 IF (IRXTYP(I).EQ.0) THEN
                         WRITE (IOUT,*) 'RDPRM:  ATTEMPT TO ',
     &                    'CHANGE RATE CONSTANT FOR ''SAMEK'''
     &                    ,'REACTION = ',INLIN(3:L1-1)
                         IERR=-1
                         RETURN
                 ELSEIF (IRXTYP(I).EQ.7) THEN
                         WRITE (IOUT,*) 'RDPRM:  ATTEMPT TO ',
     &                    'CHANGE RATE CONSTANT FOR PHOTOLY',
     &                    'SIS REACTION = ',INLIN(3:L1-1)
                         IERR=-1
                         RETURN
                 ENDIF
c                 READ (INLIN(L1+1:),'(G50.0)',ERR=800) RK(I)
c                 IRXTYP(I)=3
c first assume that we change to constant RK
                 IRXTYP(I)=3
                 SUBSTR=ADJUSTL(INLIN(L1+1:))
                 L1=INDEX(SUBSTR,' ')-1
                 READ (SUBSTR(:L1),'(G50.0)',ERR=800) RK(I)
c then check for Arrehnius RK=A(T/TREF)**B EXP(-Ea/0.0019872/T)
                 SUBSTR=ADJUSTL(SUBSTR(L1+1:))
                 L1=INDEX(SUBSTR,' ')-1
                 IF(L1.GT.0)THEN
                  IRXTYP(I)=4 
                  EA=0.0
                  B=0.0
                  READ(SUBSTR(:L1),'(G50.0)',ERR=800) A
                  SUBSTR=ADJUSTL(SUBSTR(L1+1:))
                  L1=INDEX(SUBSTR,' ')-1
                  IF(L1.GT.0)READ(SUBSTR(:L1),'(G50.0)',ERR=800) EA
                  SUBSTR=ADJUSTL(SUBSTR(L1+1:))
                  L1=INDEX(SUBSTR,' ')-1
                  IF(L1.GT.0)READ(SUBSTR(:L1),'(G50.0)',ERR=800) B
                  LOCB=LKBUF(I)
                  KPBUF(LOCB+1)=A
                  KPBUF(LOCB+2)=EA
                  KPBUF(LOCB+3)=B
                 ENDIF ! L1.GT.0
                 GOTO 100
                ENDIF ! RXNLBL.EQ.INLIN
  110           CONTINUE
C
C               SPECIFICATIONS FOR UNDEFINED REACTIONS  RESULTS IN ERROR.
C
                WRITE (IOUT,*) 'RDPRM:  UNRECOGNIZED REACTION LABEL = ',
     &           INLIN(3:L1-1),iunit
               
                IERR=-1
                RETURN
C
C
        ELSEIF (INLIN(1:4).EQ.'RHC(') THEN
C
C               W(species) ... UPDATES MOLECULAR WEIGHT AND NO. CARBONS
C                              FOR SPECIES WHOSE FIRST 4 CHARACTERS
C                              IN NAME ARE species.
C                              ALSO, ADDS THIS SPECIES TO THE LIST OF
C                              EMITTED SPECIES.  THESE ARE ASSUMED TO
C                              BE INPUT IN THE EXACT SAME ORDER AS
C                              EMISSIONS DATA.
C
                IF (INLIN(9:9).NE.')') THEN
                        WRITE (IOUT,*) 'RDPRM:  ILLEGAL SPECIES',
     &                   ' SPECIFICATION INPUT: ',INLIN(1:10)
                        IERR=-1
                        RETURN
                ENDIF
C
C               FIND WHICH SPECIES THIS IS, AND THEN UPDATE
C               PARAMETERS ACCORDINGLY.  DOES NOT SEARCH THROUGH
C               STEADY STATE SPECIES NAMES.
C
                DO 120 I=1,N1
                IF (NAME(I)(1:4).EQ.INLIN(5:8)) THEN
C
C                       NAME FOUND.  INDEX NO. = I.  UPDATE PARAMETERS
C
                        READ (INLIN,'(10X,F6.0,F9.2)',ERR=800) CNO(I),
     &                   MWT(I)
C
C                       ADD THIS TO LIST OF EMITTED SPECIES
C
                        IF (NRHC1.EQ.MAXRHC)
     &                   STOP 'TOO MANY EMITTED SPECIES'
                        NRHC1=NRHC1+1
                        INDRHC(NRHC1)=I
C
                        GOTO 100
                ENDIF
  120           CONTINUE
C
C               SPECIES NAME ON LIST.  OK OF NAME = "NROG", "LOST",
C               OR "UNKN".   ERROR OTHERWISE.
C
                IF (INLIN(5:8).EQ.'NROG') GOTO 100
                IF (INLIN(5:8).EQ.'LOST') GOTO 100
                IF (INLIN(5:8).EQ.'UNKN') GOTO 100
C
                WRITE (IOUT,*) 'RDPRM:  UNRECOGNIZED SPECIES = ',
     &           INLIN(5:8)
                IERR=-1
                RETURN
C
C
        ELSE
C
C               DETERMINE WHICH MECHANISTIC PARAMETER THIS IS.
C
                DO 130 I=1,NCOEFV
                IF (COEFNM(I).EQ.INLIN(1:8)) THEN
                        READ (INLIN,'(8X,G14.0)',ERR=800) COEF(I)
                        GOTO 100
                ENDIF
  130           CONTINUE
C
C               PARAMETER NAME NOT RECOGNIZED.  THIS IS NOT AN ERROR
C               IF THE "PARAMETER" IS ACTUALLY A SPECIES NAME.
C
                DO 132 I=1,N1
                IF (NAME(I)(1:8).EQ.INLIN(1:8)) GOTO 100
  132           CONTINUE
C
C               ALSO NOT AN ERROR IF "PARAMETER" IS SPECIES NAME
C               "NROG", "LOST", OR "UNKN"
C
                IF (INLIN(1:4).EQ.'NROG') GOTO 100
                IF (INLIN(1:4).EQ.'LOST') GOTO 100
                IF (INLIN(1:4).EQ.'UNKN') GOTO 100
C
C               ALSO NOT AN ERROR IF THE PARAMETER CODE HAS
C               THE SUFFIX "NC"
C
                L1=INDEX(INLIN,' ')-1
                IF (INLIN(L1-1:L1).EQ.'NC') GOTO 100
C
C               UNRECOGNIZED OR ILLEGAL PARAMETER INPUT
C
                WRITE (IOUT,*) 'RDPRM:  UNRECOGNIZED PARAMETER = ',
     &           '**'//INLIN(1:8)//'**',ncoefv
                do i=1,ncoefv
                 write(iout,*)i,'**'//coefnm(i)//'**'
                enddo
                IERR=-1
                RETURN
C
        ENDIF
C
C
C       ERROR DECODING INPUT
C
  800   WRITE (IOUT,*) 'RDPRM:  ERROR DECODING INPUT.  LINE = ',
     &   INLIN(1:22)
        IERR=-1
        RETURN
C
C
C       END OF FILE.  NORMAL TERMINATION
C
  900   IERR=0
        RETURN
        END



        SUBROUTINE CITPHK (cosz,FAC,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c------------------------------------------------------------------------
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
c  created by: Mike Kleeman (Feb 1999)
c         Civil and Environmental Engineering
c         University of California
c         Davis, California, 95616
c
c  This subroutine is based on the original SAPRC subroutine "NEWPHK".
c  The parameter "cosz" is now passed into the subroutine instead
c  of being calculated based on the simulation time.  These changes were
c  made to mirror changes made by Rob Harley during development of the
c  CIT airshed model.
c
c-------------------------------------------------------------------------
C
C       SUBROUTINE TO CALCULATE SOLAR PHOTOLYSIS RATE CONSTANTS
C       FOR TIME OF DAY, LATITUDE, AND DECLINATION INDICATED.  TIME
C       IS ASSUMED TO BE MINUTES AFTER MIDNIGHT.
C
C       THIS MAY BE CALLED BY RDMOD IF LISTING OF MECHANISM IS REQUIRED,
C       SO IT CAN GIVE CALCULATED PHOTOLYSIS RATES FOR Z=0 AND "DARK"
C       RATES.
C
C       THIS MUST BE CALLED PERIODICALLY AS TIME OF DAY CHANGES, SO
C       PHOTOLYSIS RATES WILL BE CURRENT.
C
C       PHKINI MUST BE CALLED BEFORE THIS IS CALLED.
C
C       ARGUMENTS:
C
C               cosine of solar zenith angle
        IMPLICIT NONE

        REAL*4 cosz
C
C               TOTAL LIGHT INTENSITY FACTOR.
        REAL*4 FAC
C
C
C       SPECIFICATION FOR MODEL-DEPENDENT DATA, INCLUDING KINETIC
C       PARAMETERS UPDATED BY THIS SUBROUTINE.
C
        INCLUDE 'modlspc.h'

        INTEGER    NS, NSC, NACT, NSB, NDUM
        INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
        INTEGER    SNO(MAXNS)
        INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
        INTEGER    IRXTYP(MAXRXN)
        INTEGER    NCOC1, NCOC2, NCOEFV
        INTEGER    I
        
        REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
        REAL*4     ONO(MAXNS), XNO(MAXNS)
        REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
        REAL*4     LAT, DEC, TSTART, CZMIN
        REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF
        REAL       FX
        
        CHARACTER*16 NAME(MAXNS)
        CHARACTER*6  RXNLBL(MAXRXN)
        CHARACTER*8  COEFNM(MAXCOV)
        CHARACTER*8  PHOTNM(MAXPHK)
        CHARACTER*80 TITLE
C
C
C       THE PARAMETER CZLOW IS THE MINIMUM COS(Z) FOR WHICH THE SOLAR
C       LIGHT INTENSITIES ARE CALCULATED USING FITS TO PETERSON'S
C       DATA.  FOR COS(Z)'S LESS THAN THAT, THEY ARE ASSUMED TO BE
C       PROPORTIONAL TO COS(Z).  THIS HAS TO BE ASSUMED, BECAUSE THE
C       FIT FUNCTION FAILS TO GIVE REASONABLE RESULTS FOR LOW VALUES
C       OF COS(Z).  THE VALUE OF CZLOW USED IS 7.0E-2, WHICH CORRESPONDS
C       TO Z=~86 DEGREES, WHICH IS THE HIGHEST Z FOR WHICH PETERSON
C       GIVES DATA.  *NOTE THAT THIS PARAMETER MUST HAVE THE SAME
C       VALUE IN THE PHKINI SUBROUTINE.
C       NOTE ALSO, THE "DPHOTK" VALUES ARE CALCULATED FOR COS(Z)=CZLOW.
C       THE "DARK" PHOTKS ARE THEN GIVEN BY (CZMIN/CZLOW) * THE VALUES
C       IN "DPHOTK"
C       CZLOW IS SPECIFIED IN MODLSPC.INC
C
C
c        print*,'cosz from citphk:',cosz
        IF (COSZ.LE.CZLOW) THEN
C
C               IF COSZ BETWEEN CZLOW AND CZMIN, ASSUME PHOTKS ARE
C               PROPORTIONAL TO COSZ, AND USE PHOTKS FOR COSZ=CZLOW
C               WHICH ARE IN DPHOTK.  IF COSZ IS LESS THAN CZMIN, USE
C               VALES CALCULATED FOR CZMIN.  NOTE THAT PHKINI DOES
C               NOT ALLOW CZMIN TO BE GREATER THAN CZLOW.
C
                IF (COSZ.GT.CZMIN) THEN
                        FX=FAC*COSZ/CZLOW
                ELSE
                        FX=FAC*CZMIN/CZLOW
                ENDIF
                DO 100 I=1,NPHOTK
  100           RK(IPHRK(I))=DPHOTK(I)*FX
C
        ELSE
C
C               COSZ POSITIVE AND GREATER GHAN CZMIN.  CALCULATE
C               SOLAR PHOTOLYSIS RATES.
C
                CALL SOLPHK (COSZ,FAC,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
C
        ENDIF
        RETURN
        END
C
C-----------------------------------------------------------------------
        SUBROUTINE PHKINI (CALCK,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
C
C       SUBROUTINE TO INITIALIZE PHOTOLYSIS RATES
C
C       CALLED BY RDMOD.  DOES NOT NEED TO BE CALLED AGAIN UNLESS THERE
C       IS A CHANGE IN THE LATITUDE, DECLINATION, ABSORBTION COEFFICEINTS
C       OR QUANTUM YIELDS, OR THE VALUES OF CZLOW OR CZMIN
C
C       ARGUMENTS:
C
C               CALCK:  IF TRUE, THE CALCULATE PHOTOLYSIS RATES FOR
C               Z=0.
        IMPLICIT NONE

        LOGICAL CALCK
C
C
C       SPECIFICATION FOR MODEL-DEPENDENT DATA, INCLUDING KINETIC
C       PARAMETERS UPDATED BY THIS SUBROUTINE.
C
        INCLUDE 'modlspc.h'
C
        REAL   PID180, RADDEG

        PARAMETER (PID180=1.745329E-2,RADDEG=4.363323E-3)
        INTEGER    NS, NSC, NACT, NSB, NDUM
        INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
        INTEGER    SNO(MAXNS)
        INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
        INTEGER    IRXTYP(MAXRXN)
        INTEGER    NCOC1, NCOC2, NCOEFV
        INTEGER    I

        REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
        REAL*4     ONO(MAXNS), XNO(MAXNS)
        REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
        REAL*4     LAT, DEC, TSTART, CZMIN
        REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF

        CHARACTER*16 NAME(MAXNS)
        CHARACTER*6  RXNLBL(MAXRXN)
        CHARACTER*8  COEFNM(MAXCOV)
        CHARACTER*8  PHOTNM(MAXPHK)
        CHARACTER*80 TITLE
C
C
C       CZLOW IS THE COS(Z) FOR WHICH PHOTKS ARE CALCULATED TO GO INTO
C       THE DPHOTK ARRAY, FOR DETERMINING LIGHT INTENSITIES AT HIGH
C       ZENITH ANGLES.  CZMIN IS THE EFFECTIVE ZENITH ANGLE FOR "DARK"
C       PHOTOLISIS RATES.  A NONZERO VALUE MAY IMPROVE INTEGRATION
C       PERFORMANCE FOR SOME SOLVERS.  CZMIN MUST BE NO GREATER THAN
C       CZLOW; IF IT IS, IT IS SET TO CZLOW.
C
        IF (CZMIN.GT.CZLOW) CZMIN=CZLOW
C
C       CALCULATE PHOTOLYSIS RATES FOR CZLOW
C
        CALL SOLPHK (CZLOW,1.0,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
        DO 102 I=1,NPHOTK
  102   DPHOTK(I)=RK(IPHRK(I))
C
C       CALCULATE Z=0 PHOTKS, IF CALCK SET.  OTHERWISE, SET THEM TO
C       ZERO.
C
        IF (CALCK) THEN
                CALL SOLPHK (1.0,1.0,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
        ELSE
                DO 104 I=1,NPHOTK
  104           RK(IPHRK(I))=0.0
        ENDIF
C
        RETURN
        END
C
C --------------------------------------------------------------------
C

        REAL FUNCTION SOLSD(AI, BI, EI, COSZ2)
        IMPLICIT NONE
        REAL AI, BI, EI, COSZ2
        SOLSD = AI*EXP(-1.0*BI*SQRT((1.0+EI)/(COSZ2+EI)))
        RETURN
        END

        SUBROUTINE SOLPHK (COSZ,FAC,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
C
C       SUBROUTINE TO CALCULATE SOLAR PHOTOLYSIS RATE FOR COS(Z) = COSZ
C       ALL MULTIPLIED BY FACTOR = FAC
C
C       CALLED BY NEWPHK OR PHKINI.  NORMALLY NOT CALLED DIRECTLY BY OTHER
C       ROUTINES.
C
C       ARGUMENTS:
C
C               COSINE OF ZENITH ANGLE
        IMPLICIT NONE

        REAL*4 COSZ
C
C               LIGHT INTENSITY FACTOR
        REAL*4 FAC
C
C
C       SPECIFICATION FOR MODEL-DEPENDENT DATA, INCLUDING KINETIC
C       PARAMETERS UPDATED BY THIS SUBROUTINE.
C
        INCLUDE 'modlspc.h'

        INTEGER    NS, NSC, NACT, NSB, NDUM
        INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
        INTEGER    SNO(MAXNS)
        INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
        INTEGER    IRXTYP(MAXRXN)
        INTEGER    NCOC1, NCOC2, NCOEFV, NSOLSD
        INTEGER    I, I1, J
        
        REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
        REAL*4     ONO(MAXNS), XNO(MAXNS)
        REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
        REAL*4     LAT, DEC, TSTART, CZMIN
        REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF
        REAL       RKI
        
        CHARACTER*16 NAME(MAXNS)
        CHARACTER*6  RXNLBL(MAXRXN)
        CHARACTER*8  COEFNM(MAXCOV)
        CHARACTER*8  PHOTNM(MAXPHK)
        CHARACTER*80 TITLE
C
C
C       PARAMETERS USED TO CALCULATE SOLAR PHOTKS:
C
C       PARAMETERS USED WERE ADJUSTED TO FIT ACTINIC FLUXES AS FUNCTION OF
C       ZENITH ANGLE CALCULATED BY PETERSON (EPA-600/4-76-025) USING
C       'BEST ESTIMATE' SURFACE ALBEDOS (SEE TABLE 4, P. 21 AND APPENDIX B
C       ,P.51).   CALCULATED FLUXES WERE FOR O3STRT=0.295  CM-ATM, AND FLU
C       XES FOR DIFFERENT O3STRT CAN'T BE CALCULATED IN THIS VERSION.  FLUXES
C       FOR WL=290-700 CALCULATED EXPLICITLY.  SDIST DEFINED FOR 290-1000,
C       WITH 700-1000 ASSUMED TO BE CONSTANT.
C
        PARAMETER (NSOLSD=48)
C
C       NSOLWL=27.  NSOLWL SPECIFIED IN IN RDMOD.FTN
C
        REAL SSDIST(NSOLWL),APARM(NSOLSD),BPARM(NSOLSD),EPARM(NSOLSD)
        REAL COSZ2, SOLSD
C
C       FOLLOWING ARE WL'S FOR WHICH PARMS ARE GIVEN
C
C       DATA SOLWL /0.295,   0.300,    0.305,    0.310,    0.315,
C     &    0.320,   0.325,   0.330,    0.335,    0.340,    0.345,
C     &    0.350,   0.355,   0.360,    0.365,    0.370,    0.375,
C     &    0.380,   0.385,   0.390,    0.395,    0.400,    0.405,
C     &    0.410,   0.415,   0.420,    0.430,    0.440,    0.450,
C     &    0.460,   0.470,   0.480,    0.490,    0.500,    0.510,
C     &    0.520,   0.530,   0.540,    0.550,    0.560,    0.570,
C     &    0.580,   0.600,   0.620,    0.640,    0.660,    0.680,
C     &    0.700,   1.000/
C
C       FOLLOWING ARE WAVELENGTHS FOR WHICH ABS*QY DATA ARE STORED.
C       SOME (HIGHER) WAVELENGTH INTERVALS LUMPED TOGETHER, AND SDIST
C       SUMMED FOR THESE.  *NOTE* THIS ARRAY IS LOCAL TO RDMOD, BUT
C       MUST BE AS SHOWN BELOW:
C
C       DATA SOLWL /     0.295,   0.300,    0.305,    0.310,    0.315,
C     &         0.320,   0.325,   0.330,    0.335,    0.340,    0.345,
C     &         0.350,            0.360,              0.370,
C     &         0.380,            0.390,              0.400,
C     &         0.410,            0.420,    0.430,    0.440,    0.450,
C     &         0.460,   0.470,   0.480,    0.490,    0.500,
C     &                  1.000/
C
        DATA APARM /
     &            4.71E16,  4.98E16,  6.13E16,  7.33E16,  7.54E16,
     &   8.06E16, 1.04E17,  1.03E17,  9.68E16,  1.00E17,  9.75E16,
     &   1.08E17, 9.69E16,  1.05E17,  1.26E17,  1.12E17,  1.22E17,
     &   9.96E16, 1.08E17,  1.10E17,  1.31E17,  1.56E17,  1.73E17,
     &   1.81E17, 1.81E17,  3.50E17,  3.54E17,  4.06E17,  4.60E17,
     &   4.66E17, 4.72E17,  4.64E17,  4.62E17,  4.80E17,  4.69E17,
     &   4.75E17, 4.75E17,  4.66E17,  4.71E17,  4.75E17,  4.85E17,
     &   9.80E17, 9.94E17,  9.83E17,  1.00E18,  1.02E18,  1.00E18,
     &   1.51E19/
C
        DATA BPARM /
     &            5.253,    3.037,    1.981,    1.354,    1.057,
     &   0.896,   0.791,    0.722,    0.671,    0.634,    0.605,
     &   0.578,   0.556,    0.533,    0.513,    0.495,    0.476,
     &   0.458,   0.443,    0.429,    0.414,    0.409,    0.394,
     &   0.384,   0.372,    0.357,    0.339,    0.326,    0.324,
     &   0.312,   0.301,    0.292,    0.285,    0.290,    0.285,
     &   0.282,   0.281,    0.278,    0.280,    0.279,    0.280,
     &   0.273,   0.273,    0.258,    0.251,    0.245,    0.236,
     &   0.236/
C
        DATA EPARM /
     &            0.046,    0.044,    0.043,    0.035,    0.034,
     &    0.034,  0.032,    0.030,    0.027,    0.025,    0.023,
     &    0.020,  0.018,    0.016,    0.015,    0.013,    0.011,
     &    0.010,  0.009,    0.008,    0.007,    0.007,    0.005,
     &    0.005,  0.004,    0.003,    0.002,    0.002,    0.002,
     &    0.001,  0.001,    0.000,    0.000,    0.000,    0.000,
     &    0.000,  0.000,    0.000,    0.000,    0.000,    0.000,
     &    0.000,  0.000,    0.000,    0.000,    0.000,    0.000,
     &    0.000/
C
C       THE FOLLOWING FUNCTION IS USED TO DEFINE SPECTRAL DIST
C        SOLSD(AI,BI,EI)=AI*EXP(-BI*SQRT((1.0+EI)/(COSZ2+EI)))
C----------------------------------------------------------------------
C
        COSZ2=COSZ**2
C
        DO 100 I=1,11
  100   SSDIST(I)=SOLSD(APARM(I),BPARM(I),EPARM(I),COSZ2)
        I1=11
        DO 102 J=12,18
        I=I1+1
        I1=I+1
  102   SSDIST(J)=SOLSD(APARM(I),BPARM(I),EPARM(I),COSZ2)
     &   +SOLSD(APARM(I1),BPARM(I1),EPARM(I1),COSZ2)
        DO 104 I=26,33
  104   SSDIST(I-7)=SOLSD(APARM(I),BPARM(I),EPARM(I),COSZ2)
        RKI=0.0
        DO 106 I=34,48
  106   RKI=RKI+SOLSD(APARM(I),BPARM(I),EPARM(I),COSZ2)
        SSDIST(27)=RKI
C
C       CALCULATE PHOTKS
C
        DO 110 J=1,NPHOTK
        RKI=0.0
        DO 112 I=1,NSOLWL
  112   RKI=RKI+SSDIST(I)*ABSQY(I,J)
        RK(IPHRK(J))=RKI*FAC
  110   CONTINUE
C
        RETURN
        END


        REAL FUNCTION ARREN(A,EA,B,TEMPR,TREF,RT)
        REAL A, TEMPR, TREF, B, EA, RT
        ARREN = A*((TEMPR/TREF)**B)*EXP(-EA/RT)
        RETURN
        END


C       FILE:  NEWRK.FTN                CHEMICAL MECH SUBROUTINE FOR
C                                       AIRSHED MODEL
C
C       LAST UPDATE:    W.P.L. CARTER  6/20/87
C
C
        SUBROUTINE NEWRK (TEMPR,ppmatm,
     +                    NS, NACT, NSC, NSB, NDUM, NRXN,
     +                    LENKBF, NCOC1, NCOC2, NCOEFV, NPHOTK,
     +                    TREF, LAT, DEC, TSTART, CZMIN, 
     +                    MWT, CONC0, CNO, NNO, SNO, ONO, XNO,
     +                    NRHC1, INDRHC,
     +                    RK, KPBUF, COEF, ABSQY, DPHOTK,
     +                    LKBUF, IPHRK, IRXTYP,
     +                    NAME, RXNLBL, COEFNM, PHOTNM, TITLE)
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
C
C       SUBROUTINE TO CALCULATE TEMPERATURE DEPENDENT RATE CONSTANTS
C       FOR TEMPERATURE INDICATED.  TEMPR = NEW TEMPERATURE FOR RATE
C       CONSTANTS.
C
C       THIS IS CALLED BY RDMOD WHEN MECHANISM DATA ARE INITIALIZED,
C       AND SHOULD BE CALLED EVERY TIME THE TEMPERATURE CHANGES.
C
C
C       UPDATES RK'S FOR REACTION TYPES 3 (ARRHENIUS) AND 4 (FALLOFF)
C       ONLY.  OTHER RK'S UNCHANGED.
C
C       THIS VERSION CALCULATES FALLOFF K'S FOR PRESSURE = 1 ATM
C       (M = 1.0E+6 PPM) ONLY.  TO GIVE THIS PROGRAM A VARIABLE PRESSURE
C       CAPABILITY, CHANGE THE PARAMETER PPMATM, WHICH IS THE TOTAL
C       PPM PRESENT, TO A PROGRAM VARIABLE
C
c        PARAMETER (PPMATM=1.0E+6)
C
C       SPECIFICATION FOR MODEL-DEPENDENT DATA, INCLUDING KINETIC
C       PARAMETERS UPDATED BY THIS SUBROUTINE.
C
        IMPLICIT NONE
        
        INCLUDE 'modlspc.h'
C
C       LOCAL SPECIFICATIONS
C
        REAL K0, KINF
        INTEGER    NS, NSC, NACT, NSB, NDUM
        INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
        INTEGER    SNO(MAXNS)
        INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
        INTEGER    IRXTYP(MAXRXN)
        INTEGER    NCOC1, NCOC2, NCOEFV
        INTEGER    IRXN, LOCB
        
        REAL*4     MWT(MAXNS), CNO(MAXNS), CONC0(MAXNS), RK(MAXRXN)
        REAL*4     ONO(MAXNS), XNO(MAXNS)
        REAL*4     KPBUF(MAXKBF), COEF(MAXCOE)
        REAL*4     LAT, DEC, TSTART, CZMIN
        REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF
        REAL       TEMPR, PPMATM, RT, RATIO
        REAL       ARREN
        
        CHARACTER*16 NAME(MAXNS)
        CHARACTER*6  RXNLBL(MAXRXN)
        CHARACTER*8  COEFNM(MAXCOV)
        CHARACTER*8  PHOTNM(MAXPHK)
        CHARACTER*80 TITLE
C
C       LOCAL STATEMENT FUNCTION DEFINITION
C
C        ARREN(A,EA,B)=A*((TEMPR/TREF)**B)*EXP(-EA/RT)
C
C       
C        print*,'RXNLBL',RXNLBL
        RT=0.0019872*TEMPR
C
C       LOOP OVER REACTIONS INVOLVED
C
        DO 1000 IRXN=1,NRXN
C
C       DETERMINE KINETICS TYPE, AND CALCULATE K ACCORDINGLY
C       ONLY TYPES 4 AND 5 INVOLVED.
C
        IF (IRXTYP(IRXN).EQ.4) THEN
C
C               CALCULATE RK USING ARRHENIUS EXPRESSION
                LOCB=LKBUF(IRXN)
c MJK corrects the pressure dependance 5/25/12
                RK(IRXN)=(PPMATM/1.e6)**nint(KPBUF(LOCB))
     +               *ARREN(KPBUF(LOCB+1),KPBUF(LOCB+2),KPBUF(LOCB+3),
     &               TEMPR, TREF, RT)
c MJK corrects the pressure dependance 5/25/12
C
        ELSEIF (IRXTYP(IRXN).EQ.5) THEN
C
C               CALCULATES RATE CONSTANT USING TROE FALLOFF EXPERSSION
C               FOR TEMPERATURE = TEMPR
C               KPBUF(0-2) = A,EA,B FOR P=0
C               KPBUF(3-5) = A,EA,B FOR P=INF.
C               KPBUF(6,7) = F AND N PARAMETERS
C
                LOCB=LKBUF(IRXN)
C                print*,'LOCB =',LOCB,'RXNLBL(LOCB) = ',RXNLBL(LOCB),
C     &           'IRXN = ',IRXN
c MJK and JH corrects the pressure dependance 5/29/12
c                K0=ARREN(KPBUF(LOCB),KPBUF(LOCB+1),KPBUF(LOCB+2))*PPMATM
c                RATIO=K0/ARREN(KPBUF(LOCB+3),KPBUF(LOCB+4),KPBUF(LOCB+5
c     &           ))
                K0=ARREN(KPBUF(LOCB),KPBUF(LOCB+1),KPBUF(LOCB+2),
     &               TEMPR,TREF,RT)*PPMATM**2/1.e6
                KINF=ARREN(KPBUF(LOCB+3),KPBUF(LOCB+4),KPBUF(LOCB+5),
     &              TEMPR,TREF,RT)*PPMATM/1.e6
                RATIO=K0/KINF
c MJK and JH corrects the pressure dependance 5/29/12
                RK(IRXN)=(K0/(1.0+RATIO))*KPBUF(LOCB+6)**(1.0/(1.0+
     &           (0.43429*ALOG(RATIO)/KPBUF(LOCB+7))**2))
C 
C Find reaction label = 44 and print rate constant. 
                IF (RXNLBL(IRXN) .EQ. '44') THEN
                  print*,'RK(IRXN)= ', RK(IRXN)
                  print*,'RXNLBL(IRXN)= ',RXNLBL(IRXN)
                ENDIF
C
        ENDIF
C
 1000   CONTINUE
C
        RETURN
        END

      function nblank(buf)
c*****************************************************************************
c
c written by: Mike Kleeman (May 2002)
c             Dept of Civil and Env Eng
c             UC Davis 95616
c
c The purpose of this function is to return the rightmost position of
c non blank characters in the input buffer.
c
c*******************************************************************************
      IMPLICIT NONE
      
      integer    nblank, ilen
      INTEGER    I
      character* (*) buf

      nblank = 1
      ilen = len(buf)
      do i = ilen,1,-1
       if(buf(i:i).ne.' ')then
        nblank = i
        return
       endif
      enddo
      return
      end

