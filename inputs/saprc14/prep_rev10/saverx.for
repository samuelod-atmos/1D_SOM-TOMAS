C       FILE:  SAVERX.FTN               PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  1/4/88
C		updated MJM to include SNO(I) 11/9/03
C
C
        SUBROUTINE SAVERX
C
C       WRITES DATA FOR PREPARED MODEL ONTO ".MOD" FILE.
C       DATA SET IS WRITTEN IN ASCII FORMAT.
C
C       THIS INCLUDES PHOTOLYSIS ABS QY'S, WHICH ARE READ FROM TEMPORARY
C       DATA SETS, ALLIGNED IN STANDARD WAVELENGTHS FOR SOLAR PHOTK CALC,
C       AND THEN SAVED.
C
C       CALLED BY:      PREP MAIN
C
C       SUBROUTINES CALLED:
C
C               "NEWSUBS" (SAPRC) LIBRARY UTILITY ROUTINES
C
C       FILNAM  ... ROUTINE TO CREATE FILE NAME
C       ALIN16  ... SPLITS UP INPUT STRING INTO GROUPS OF 16 CHARS
C
C               PREP SUBROUTINES
C
C       PFALIN  ... ALLIGNS ABS QY'S INTO FORM REQUIRED FOR SOLAR
C                   PHOTK CALCULATIONS
C
C
C       SPECIFICATION OF PREP PARAMETERS, VARIABLES AND ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL VARIABLES AND PARAMETERS
C
        CHARACTER*1 RXNSTR(160)
        CHARACTER*9 FORM
        LOGICAL   NXTNEW
        REAL*4 WL(MAXAQY),ABSQY(MAXAQY)
C
C       FOLLOWING ARE WAVELENGTHS FOR WHICH ABS*QY DATA ARE
C       STORED.  SOME (HIGHER) WAVELENGTH INTERVALS LUMPED
C       TOGETHER, AND SDIST SUMMED FOR THESE.
C
        PARAMETER (NSOLWL=27,NSOLW1=28)
        REAL SOLWL(NSOLW1),AQ(NSOLWL)
C
        DATA SOLWL /     0.295,   0.300,    0.305,    0.310,    0.315,
     &          0.320,   0.325,   0.330,    0.335,    0.340,    0.345,
     &          0.350,            0.360,              0.370,
     &          0.380,            0.390,              0.400,
     &          0.410,            0.420,    0.430,    0.440,    0.450,
     &          0.460,   0.470,   0.480,    0.490,    0.500,
     &                   1.000/
C
C       'MOD' FILE FORMAT INDICATOR
C
        DATA FORM /'UCD:REV9'/
C
C
C      NOTE: DON'T SAVE DELETABLE SPECIES
       NS=NS-NDEL
C
C
        IF (NOPRP) GOTO 1000
C
C       OPEN MODEL FILE
C
        IOB160=' '
        CALL FILNAM(I,IOB160,MODUIC,MODFIL,'.mod ')
        OPEN (UNIT=UNIT1,NAME=IOB160,TYPE='UNKNOWN')
c     &   CARRIAGECONTROL='LIST')
C
C       WRITE GENERAL PARMS
C
        I=0
        WRITE (UNIT1,'(A9,'' = MODEL FORMAT'')') FORM
        WRITE (UNIT1,101) TITLE,NS,NDY,NSC,NB,NDUM,NRXN,NCOEFV,
     &   NCOC1,NCOC2,NPHOTK,LOCKBF,NSOLWL,TEMPR,TREF
 101    FORMAT (A80/12I6/2F10.2)
C
C
C      WRITE SPECIES NAMES AND CHARACTERISTICS
C
       WRITE (UNIT1,'(A16,0PF7.2,F6.2,2I3,2F6.2,1PE11.3)') (NAME(I),
     &  MWT(I),CNO(I),NNO(I),SNO(I),ONO(I),XNO(I),CONC0(I),I=1,NA2)
       IF (NSS.NE.0) WRITE (UNIT1,'(5A16)') (NAME(I),I=NS1,NS)
C
C      WRITE COEF'S
C
       IF (NCOEFV.GT.0) WRITE (UNIT1,'((4(A8,1PE11.3,1X)))')
     &   (COEFNM(I),COEF(I),I=1,NCOEFV)
       IF (NCOC2.GE.NCOC1) WRITE (UNIT1,'(1P6E12.3)')
     &   (COEF(I),I=NCOC1,NCOC2)
C
C      WRITE ABSQY'S
C
 1000   IF (NPHOTK.GT.0) THEN
C
        WRITE (OUT,81) NPHOTK
   81   FORMAT (/1X,I2,' TYPES OF PHOTOLYSIS REACTIONS:')
C
        DO 80 IPHOT=1,NPHOTK
C
        IOBF16(1:8)=PHOTNM(IPHOT)
        IOBF16(9:9)=' '
        IOB160=' '
        IF (PFTYPE(IPHOT)) THEN
                CALL FILNAM (I,IOB160,TMPUIC,IOBF16,'.PTM ')
                OPEN (UNIT=UNIT2,NAME=IOB160,TYPE='OLD',ERR=69)
        ELSE
                CALL FILNAM (I,IOB160,PHFUIC,IOBF16,'.PHF ')
                OPEN (UNIT=UNIT2,NAME=IOB160,READONLY,TYPE='OLD',ERR=69)
        ENDIF
        READ (UNIT2,50,END=70,ERR=70) IOB160
   50   FORMAT (A80)
        WRITE (OUT,82) PHOTNM(IPHOT),IOB160,RXNLBL(IPHOTR(IPHOT))
   82   FORMAT (/1X,5X,'FILE = ',A8,6X,'LABEL =',A80/' ',22X,
     &   'FIRST REACTION = ',A6)
C       WRITE (OUT,82) PHOTNM(IPHOT),IOB160,(RXNLBL(IPHOTR(I,IPHOT)),
C     &    I=1,N)
C   82  FORMAT (/1X,5X,'FILE = ',A8,6X,'LABEL =',A80/' ',22X,
C     &  'REACTIONS =',15(1X,A6),:/(' ',33X,15(1X,A6)))
        N=0
        F=1.0
   52   READ (UNIT2,50,ERR=70,END=75) IOB160
        IF (IOB160(1:1).EQ.'!') GOTO 52
        IF (IOB160(1:8).EQ.'        ') GOTO 75
        IF (IOB160(1:2).EQ.'FA') GOTO 54
C
        N=N+1
        CALL ALIN16 (3,IOB160)
        DECODE (48,53,IOB160,ERR=73) WL(N),ABSQY(N),X
   53   FORMAT (3G16.0)
        IF (IOB160(41:48).EQ.'        ') X=1.0
        ABSQY(N)=ABSQY(N)*X*F
        GOTO 52
C
  54    CALL ALIN16 (2,IOB160)

        DECODE (16,53,IOB160(17:80),ERR=73) F
        IF (F.LE.0.0) GOTO 73
        GOTO 52
C
   69   WRITE (OUT,71) PHOTNM(IPHOT),RXNLBL(IPHOTR(I))
   71   FORMAT ('0ERROR: FILE = ',A8,' COULD NOT BE OPENED.'/' ',22X,
     &   'FIRST REACTION = ',A6)
C   69  WRITE (OUT,71) PHOTNM(IPHOT),(RXNLBL(IPHOTR(I,IPHOT)),I=1,N)
C   71  FORMAT ('0ERROR: FILE = ',A8,' COULD NOT BE OPENED.'/' ',22X,
C     &  'REACTIONS =',15(1X,A6),:/(' ',33X,15(1X,A6)))
        GOTO 72
C
   70   N=N+1
        WRITE (OUT,701) N
  701   FORMAT (' ',5X,'ERROR IN READING FILE, LINE NO. =',I4)
        GOTO 72
C
   73   WRITE (OUT,702) IOB160(1:16),IOB160(17:32),IOB160(33:48)
  702   FORMAT (' ',5X,'ERROR IN DECODING FILE INPUT. INPUT = ''',
     &   A16,''', ''',A16,''', ''',A16,'''')
C
c   72   WRITE (ICRT,74) PHOTNM(IPHOT)
   72   WRITE (ICRT,74) IOB160
   74   FORMAT (' ERROR IN READING PHOTOLYIS FILE = ''',(A),'''')
        NERR=NERR+1
        N=0
        WL(1)=0.0
        ABSQY(1)=0.0
   75   IF (PFTYPE(IPHOT)) THEN
C               (DELETE IF TEMPORARY FILE)
                CLOSE (UNIT=UNIT2,STATUS='DELETE')
        ELSE
                CLOSE (UNIT=UNIT2)
        ENDIF
C
C       ALLIGN TO SOLAR WAVELENGTHS
C
        CALL PFALIN (AQ,N,WL,ABSQY,NSOLWL,SOLWL)
C
        WRITE (OUT,83) (SOLWL(I),AQ(I),I=1,NSOLWL)
   83   FORMAT (' ',20X,'(WL,ABS*QY) =',5('  (',0PF5.3,1PE10.3,')',:)/
     &   (' ',33X,5('  (',0PF5.3,1PE10.3,')',:)))
        IF (YESPRP) THEN
                WRITE (UNIT1,'(A8,I5,7X,1P6E10.3/(8E10.3))')
     &           PHOTNM(IPHOT),IPHOTR(IPHOT),(AQ(I),I=1,NSOLWL)
        ENDIF
   80   CONTINUE
C
C       ENDIF FOR NPHOTK.GT.0
C
        ENDIF
C
        IF (NOPRP) GOTO 2000
C
C       WRITE KINETIC  PARMS DATA
C
        WRITE (UNIT1,'(1P8E10.3)') (KPBUF(I),I=1,LOCKBF)
        WRITE (UNIT1,'(16I5)') (LKBUF(IRXN),IRXN=1,NRXN)
        WRITE (UNIT1,'(16I5)') (RXTYP(IRXN),IRXN=1,NRXN)
C
C       WRITE REACTION LABELS
C
c        WRITE (UNIT1,'(10A8)') (RXNLBL(IRXN),IRXN=1,NRXN)
        WRITE (UNIT1,'(10(A6,2x))') (RXNLBL(IRXN),IRXN=1,NRXN)
C
C       WRITE REACTIONS.  REWIND TEMPORARY FILE UNIT WITH REACTION LISTS
C
        IOB160=' '
        CALL FILNAM (NLEN,IOB160,TMPUIC,'PREP.TMP ',' ')
        OPEN (UNIT=UNIT3,NAME=IOB160,TYPE='OLD',FORM='UNFORMATTED')
C
 100    READ (UNIT3,END=110) NXTNEW,L,(BYT160(I),I=1,L)
        IF (L.GT.78) THEN
                WRITE (OUT,102) 
 102            FORMAT (' WARNING:  STRING FOR RXN TRUNCATED')
                write(out,'(80a1)')(byt160(i),i=1,L)
                L=78
        ENDIF
        WRITE (UNIT1,92) NXTNEW,(BYT160(I),I=1,L)
  92    FORMAT (L1,1X,78A1)
        GOTO 100
 110    CONTINUE
C
        CLOSE (UNIT=UNIT3,STATUS='DELETE')
C
C       DONE
C
        CLOSE (UNIT=UNIT1)
C
C       SUMMARIZE LENGTH OF MODEL
C
 2000   NCOCF=NCOC2-NCOC1+1
        WRITE (OUT,1003) TEMPR,NRXN,NS,NCOEFV,NCOCF,LOCKBF
 1003   FORMAT ('0DEFAULT TEMPERATURE =',F8.2/' NO. REACTIONS ='
     &   I5,/' NO. SPECIES =',I5/' NO. VARIABLE COEFFICIENTS =',I5
     &   /' NO. CONSTANT COEFFICIENTS =',I5
     &   /' KINETIC PARM BUFFER LENGTH =',I5)
        IF (YESPRP) THEN
                WRITE (OUT,1005) NPRODS(NRXN+1)-1, MODFIL
 1005           FORMAT (' PRODUCT ARRAY LENGTH =',I7/' '/
     &           ' MODEL = ''',A16,''' STORED.')
        ELSE
                WRITE (OUT,1006) MODFIL
 1006           FORMAT (' MODEL = ''',A16,''' NOT STORED.')
        ENDIF
        RETURN
        END
