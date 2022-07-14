C       FILE:  RDRXN.FTN                PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  1/4/88
C
C	MJM added reading in of SNO (# of sulfur) 11/8/03
C
C
        SUBROUTINE RDRXN
C
C
C       SUBROUTINE TO READ IN MECHANISM-SPECIFIC DATA.  INITIALIZES MECHANISM
C       PARAMETERS, READS REACTIONS AND KINETIC PARAMETERS, READS IN INPUT
C       SPECIFYING REACTANT TYPES, ETC.  PROCESSES ALL INPUT IN MECHANISM
C       INPUT FILE EXCEPT CONTROL PARAMETERS ALREADY READ BY PREP MAIN.
C       IN MODETERMINE WHICH SPECIES
C
C
C       CALLED BY:      PREP MAIN
C
C       SUBROUTINES CALLED:
C
C               "NEWSUBS" (SAPRC) LIBRARY UTILITY ROUTINES
C
C       FILNAM  ... ROUTINE TO CREATE FILE NAME
C       INBUF   ... GET LINE FROM INPUT FILE.  GIVES PROMPT IF INPUT
C                   FROM USER REQUIRED.
C       ALIN16  ... SPLITS UP INPUT STRING INTO GROUPS OF 16 CHARS
C       ALIN8   ... SPLITS UP INPUT STRING INTO GROUPS OF 8 CHARS
C       MOVLFT  ... MOVES CHARACTERS IN 16-CHAR STRING TO LEFT
C
C               PREP SUBROUTINES (ALL IN FILE RXLST1.FTN)
C
C       RXLST1  ... PROCESSES A CHARACTER STRING DESCRIBING REACTANTS,
C                   PRODUCTS AND COEFFICIENTS FOR A REACTION.  CREATS
C                   NEW SPECIES AND COEFFICIENTS AS NEEDED.  CALLS
C                   SUBROUTINES
C       RXLST2  ... PROCESSES A LIST OF SPECIES, GIVING THEM A SPECIES
C                   TYPE AS INDICATED BY CURRENT VALUE OF "IND" VARIABLE.
C                   CREATES NEW SPECIES AS NEEDED.
C       SPCNAM  ... DETERMINES SPECIES NUMBER FROM NAME.  DEFINES TYPE
C                   BASED ON CURRENT VALUE OF "IND".  CREATES NEW SPECIES
C                   IF NEEDED.
C
C
C       INCLUDES SPECIFICATIONS OF COMMON VARIABLES, PARMS, AND ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C       SPECIFICATIONS OF LOCAL PARMS, VARIABLES AND ARRAYS
C
        PARAMETER (MXCOND=30,MXDEFD=30)
        LOGICAL   BLANKC,DUMRXN,KEQ,RXOPN,NEWRXN
        INTEGER   NDEFD,NCOND
        CHARACTER*8 LBCOND(MXCOND),LBDEFD(MXDEFD),GNAME
        CHARACTER*6 LBL,LBL2
        CHARACTER RKBUF*80,TMPRXN*160,SNAME*16
        CHARACTER*1 RKBF1(80)
        EQUIVALENCE (RKBUF,RKBF1(1))
C
C
C
C       INITIALIZE
C
        PPM=.FALSE.
        NEWRXN=.TRUE.
        IF (NOPRP) THEN
                RXOPN=.TRUE.
                VCODFD=.TRUE.
        ELSE
                RXOPN=.FALSE.
                VCODFD=.FALSE.
        ENDIF
        DO 5 I=1,8
    5   VCOCAL(I)=.FALSE.
C
C       DEFAULT SPECIES TYPE = ACTIVE
C
        IND=4
C
C
C       CREATE SPECIAL COEFFICIENT "LITTLE" WHICH IS MINIMUM CONSUMPTION
C       RATE FOR STEADY STATE SPECIES.  IT SHOULD BE SMALL TO HAVE NO
C       EFFECT ON CALCULATION, BUT BE NON-ZERO TO AVOID DIVIDE-BY-ZERO
C       AT THE START OF CALCULATIONS IF NO REACTIONS FORM OR CONSUME
C       THE SPECIES.
C       NO. VARIABLE COEFS INITIALIZED TO 1.
C
        INDLIL=1
        COEFNM(1)='LITTLE  '
        COEF(1)=LITTLE
        NCOEFV=1
C
C       INITIALIZE NO. SPECIES, REACTIONS, KINETIC PARMS, CONSTANT COEFS
C
        NS=0
        NRXN=0
        LOCKBF=0
        NCOC2=NCOC1-1
C
        do i=1,maxns
         TATOM(i)=" "
        enddo
C
C
C      MAIN COMMAND SECTION:
C
C      'COMMAND' INPUT CONSISTS OF LINES WITH '.' IN COL. 1.
C
C      COMMAND INPUT PROCESSING:
C
C      .  (ALONE)      -DELINIATES SECTIONS  (SEE BELOW).  ALL NON-COMMAND
c                       INPUT FOLLOWING IT UNTIL NEXT COMMAND IGNORED
C      . label         -DEFINES CONDITION = label AS TRUE (SEE BELOW)
C      . IF label      -FOLLOWING INPUT IGNORED UNTIL A BLANK COMMAND
C                       INPUT IS READ UNLESS CONDITION = label HAS BEEN
C                       DEFINED AS TRUE.  THIS INCLUDES NON-BLANK COMMAND
C                       INPUT.
C      . FOR label     -FOLLOWING INPUT IGNORED UNTIL A BLANK COMMAND
C                       INPUT IS READ IF A PREVIOUS 'FOR' COMMAND CARD
C                       WITH THE SAME LABEL = label HAS BEEN PROCESSED
C      .cmd            -FOLLOWING NON-COMMAND INPUT (E.G. REACTIONS,
C                       SPECIES DEFS, ETC) PROCESSED ACCORDING TO WHAT
C                       COMMAND cmd IS UNTIL ANOTHER COMMAND INPUT
C                       IS INPUT.  SEE CODE BELOW FOR LEGAL COMMANDS.
C      .cmd IF label   -FOLLOWING INPUT IGNORED UNTIL ANOTHER COMMAND
C                       IS INPUT UNLESS CONDITION = label HAS BEEN
C                       DEFINED AS TRUE.
C      .cmd FOR label  -FOLLOWING INPUT IGNORED UNTIL ANOTHER COMMAND
C                       IS INPUT IF A PREVIOUS 'FOR' COMMAND CARD
C                       WITH SAME LABEL HAS BEEN PROCESSED.
C
C       WHEN THIS SUBROUTINE IS FIRST CALLED, THE FIRST INPUT LINE
C       HAS ALREADY BEEN READ BY PREP MAIN, SO SKIP PAST INPUT PART.
C
       GOTO 1
C
  100  CALL INBUF ('RDRXN  :',8)
       IF (IORET.EQ.2) GOTO 1999
       IF (IOLEN.EQ.0) GOTO 100
       IF (IOB160(1:1).NE.'.') GOTO 100
C
    1  CALL ALIN8(3,IOB160)
       IF (IOB160(1:8).EQ.'  .UNITS') GOTO 150
       IF (IOB160(1:8).EQ.'   .PHOT') GOTO 1500
       IF (IOB160(1:8).EQ.'    .INS') GOTO 300
       BLANKC=IOB160(7:8).EQ.' .'
       IF (IOB160(9:16).NE.'        ') GOTO 110
       IF (BLANKC) GOTO 100
       GOTO 180
  110  IF (IOB160(9:16).EQ.'      IF') GOTO 120
       IF (IOB160(9:16).EQ.'     FOR') GOTO 130
       IF (BLANKC) GOTO 140
       WRITE (OUT,111) IOB160(9:16),IOB160(5:8)
  111  FORMAT (' UNRECOGNIZED COMMAND QUALIFIER = ''',A8,
     &  '''.  COMMAND = ''',A4,''' PROCESSED ANYWAY.')
       WRITE (ICRT,113) IOB160(1:48)
  113  FORMAT (' UNRECOGNIZED COMMAND QUALIFIER.',
     &  '  COMMAND PROCESSED ANYWAY.'/' DECODED LINE =  ''',
     &  A48,'''')
       NERR=NERR+1
       GOTO 180
C
C      'IF' PROCESSING.  SEE IF LABEL MATCHES LIST IN LBCOND
C
  120  IF (NCOND.EQ.0) GOTO 123
       DO 122 I=1,NCOND
       IF (LBCOND(I).EQ.IOB160(17:24)) GOTO 125
  122  CONTINUE
C  -     DOESN'T MATCH
  123  IF (BLANKC) GOTO 170
       GOTO 100
  125  IF (BLANKC) GOTO 100
       GOTO 180
C
C      'FOR' PROCESSING.
C
  130  IF (NDEFD.EQ.0) GOTO 133
       DO 132 I=1,NDEFD
       IF (LBDEFD(I).EQ.IOB160(17:24)) GOTO 135
  132  CONTINUE
C  -     ADD NEW 'FOR' LABEL, AND PROCESS COMMAND
  133  NDEFD=NDEFD+1
       LBDEFD(NDEFD)=IOB160(17:24)
       IF (BLANKC) GOTO 100
       GOTO 180
C  -     MATCHES PREVIOUS 'FOR' LABEL.  BYPASS PROCESSING
  135  IF (BLANKC) GOTO 170
       GOTO 100
C
C      ADD NEW CONDITION, AFTER CHECKING TO SEE NOT DUPLICATE
C
  140  IF (NCOND.EQ.0) GOTO 144
       DO 142 I=1,NCOND
       IF (IOB160(9:16).EQ.LBCOND(I)) GOTO 100
  142  CONTINUE
  144  NCOND=NCOND+1
       LBCOND(NCOND)=IOB160(9:16)
       GOTO 100
C
C      '.UNITS' PROCESSING.  FOR REACTION INPUT, 'UNITS=PPM'
C              CONVERTS UNITS, ASSUMED TO BE CMS UPON INPUT
C              INTO PPM-MIN UNITS.  '.UNITS=OK' CANCELS
C              '.UNITS=PPM' COMMAND
C
  150  IF (IOB160(13:16).EQ.' PPM') GOTO 154
       IF (IOB160(13:16).EQ.'  OK') GOTO 152
       WRITE (OUT,151) IOB160(9:16)
  151  FORMAT (' UNRECOGNIZED ''.UNITS'' SUBFIELD = ''',A8,'''')
       NERR=NERR+1
       WRITE (ICRT,151) IOB160(9:16)
  152  PPM=.FALSE.
       GOTO 100
  154  PPM=.TRUE.
       GOTO 100
C
C      '. ' FIELD TO BE BYPASSED
C
  170  CALL INBUF ('BYPASS  :',8)
       IF (IORET.EQ.2) GOTO 1999
       IF (IOLEN.EQ.0) GOTO 170
       IF (IOB160(1:2).EQ.'. ') GOTO 1
       GOTO 170
C
C      DETERMINE cmd
C
  180   IF (IOB160(5:8).EQ.'.CON') THEN
                IND=1
                GOTO 200
        ELSEIF (IOB160(5:8).EQ.'.BLD') THEN
                IND=2
                GOTO 200
        ELSEIF (IOB160(5:8).EQ.'.DUM') THEN
                IND=3
                GOTO 200
        ELSEIF (IOB160(5:8).EQ.'.STS') THEN
                IND=5
                GOTO 200
        ELSEIF (IOB160(5:8).EQ.'.ACT') THEN
                IND=4
                GOTO 200
C       reading in the type-o-matic info
        ELSEIF (IOB160(5:8).EQ.'.TYP') THEN
                GOTO 240
C       reading in the SOM defintions
        ELSEIF (IOB160(5:8).EQ.'.SOM') THEN
                GOTO 250
C
        ELSEIF (IOB160(5:8).EQ.'.RXN') THEN
                GOTO 1000
        ELSEIF (IOB160(5:8).EQ.'.COE') THEN
                GOTO 1200
        ELSEIF (IOB160(4:8).EQ.'.COEF') THEN
                GOTO 1200
C
        ELSEIF (IOB160(5:8).EQ.'.#C1') THEN
                IV=1
                GOTO 400
        ELSEIF (IOB160(5:8).EQ.'.#D1') THEN
                IV=2
                GOTO 400
        ENDIF
C
        IF (IOB160(5:8).EQ.'.END') GOTO 2000
C
C       COMMAND NOT RECOGNIZED
C
        WRITE (OUT,189) IOB160(5:8)
  189   FORMAT (' UNRECOGNIZED RDRXN COMMAND = ''',A4,'''')
        NERR=NERR+1
        IF (OUT.NE.1) WRITE (ICRT,189) IOB160(5:8)
        GOTO 100
C
C
C       END COMMAND PROCESSING SECTION.  PROCESS SPECIFIC COMMANDS
C
C
C       SPECIES DEFINITION ... READ NEXT SPECIES
C
  200   CALL INBUF ('SPECIE :',8)
        IF (IORET.EQ.2) GOTO 1999
        IF (IOLEN.EQ.0) GOTO 200
        IF (IOB160(1:1).EQ.'!') GOTO 200
        IF (IOB160(1:1).EQ.'.') GOTO 1
        IF (IOB160.EQ.' ') GO TO 200
        IF (IOB160(1:1).EQ.'='.OR.IOB160(1:2).EQ.' =') GO TO 220
C
        CALL ALIN16(8,IOB160)
        SNAME = IOB160(:16)
        CALL MOVLFT(SNAME)
        CALL SPCNAM(I2,SNAME,.FALSE.)
        READ (IOB160,'(16X,2F16.0)',ERR=201) CONC0(I2),MWT(I2)
C
C CODE WAS ADDED HERE TO ALLOW THE USER TO INPUT BOTH THE NUMBER
C OF TRACED ATOMS IN EACH SPECIES, BUT ALSO THE NUMBER OF SOURCES
C THAT THOSE ATOMS COULD HAVE ORIGINATED FROM.  THIS IS USEFUL WHEN
C LOOKING AT LARGE SPECIES THAT HAVE MANY TRACED ATOMS BUT ONLY 1
C OR 2 POTENTIAL SOURCES.  MJK 9/05
C
c        READ (IOB160,'(32X,F16.0,2I16,F16.0)',ERR=201)
c     +  CNO(I2),NNO(I2),SNO(I2),XNO(I2)
c--carbon--
        ib = 49
        j1 = index(iob160(ib:ib+15),'(')-2
        j2 = index(iob160(ib:ib+15),')')-2
        IF(j1.eq.-2)THEN
         READ(IOB160(ib:ib+15),'(F16.0)')CNO(I2)
         CSNO(I2) = CNO(I2)
        ELSE
         READ(IOB160(ib:ib+j1),'(F16.0)')CNO(I2)
         READ(IOB160(ib+j1+2:ib+j2),'(I16)')CSNO(I2)
        ENDIF
c--nitrogen--
        ib = 65
        j1 = index(iob160(ib:ib+15),'(')-2
        j2 = index(iob160(ib:ib+15),')')-2
        IF(j1.eq.-2)THEN
         READ(IOB160(ib:ib+15),'(I16)')NNO(I2)
         NSNO(I2) = NNO(I2)
        ELSE
         READ(IOB160(ib:ib+j1),'(I16)')NNO(I2)
         READ(IOB160(ib+j1+2:ib+j2),'(I16)')NSNO(I2)
        ENDIF
c--sulfur--
        ib = 81
        j1 = index(iob160(ib:ib+15),'(')-2
        j2 = index(iob160(ib:ib+15),')')-2
        IF(j1.eq.-2)THEN
         READ(IOB160(ib:ib+15),'(I16)')SNO(I2)
         SSNO(I2) = SNO(I2)
        ELSE
         READ(IOB160(ib:ib+j1),'(I16)')SNO(I2)
         READ(IOB160(ib+j1+2:ib+j2),'(I16)')SSNO(I2)
        ENDIF
c--oxygen--
        ib = 97
        j1 = index(iob160(ib:ib+15),'(')-2
        j2 = index(iob160(ib:ib+15),')')-2
        IF(j1.eq.-2)THEN
         READ(IOB160(ib:ib+15),'(F16.0)')ONO(I2)
         OSNO(I2) = ONO(I2)
        ELSE
         READ(IOB160(ib:ib+j1),'(F16.0)')ONO(I2)
         READ(IOB160(ib+j1+2:ib+j2),'(I16)')OSNO(I2)
        ENDIF
c--other--
        ib = 113
        j1 = index(iob160(ib:ib+15),'(')-2
        j2 = index(iob160(ib:ib+15),')')-2
        IF(j1.eq.-2)THEN
         READ(IOB160(ib:ib+15),'(F16.0)')XNO(I2)
         XSNO(I2) = XNO(I2)
        ELSE
         READ(IOB160(ib:ib+j1),'(F16.0)')XNO(I2)
         READ(IOB160(ib+j1+2:ib+j2),'(I16)')XSNO(I2)
        ENDIF
C
C	MJM DEBUG
C
c        write (*,'(A16,3F16.0,5I16,2(F16.0,I16))',ERR=201) NAME(I2),
c     &   CONC0(I2),MWT(I2),CNO(I2),CSNO(I2),NNO(I2),NSNO(I2),
c     &   SNO(I2),SSNO(I2),ONO(I2),OSNO(I2),XNO(I2),XSNO(I2)

        IF (NAME(I2).EQ.'M') INDM=I2
        GOTO 200
C
  201   WRITE (OUT,203) NAME(I2),IOB160(17:96)
  203   FORMAT (' ERROR IN READING PARMS FOR SPECIES = ',A16
     &   /'  ALLIGNED INPUT LINE FOLLOWS:'/A80)
        NERR=NERR+1
        WRITE (ICRT,203) NAME(I2),IOB160(17:96)
        GOTO 200
C
C  -     SPECIES INPUT AS 'PRODUCT' LIST
  220  I2=0
       CALL RXLST2 (IOB160,I2)
       IF (I2.GT.0) GOTO 200
       IF (IORET.EQ.2) GOTO 1999
       WRITE (ICRT,228) IOB160
  228  FORMAT (' ERROR IN READING REACTANT LIST GIVEN BELOW.'
     &  /' ',A80)
       GOTO 200
C
c      READ IN THE TYPE-O-MATIC INFO
c
  240   CALL INBUF ('SPECIE :',8)
        IF (IORET.EQ.2) GOTO 1999
        IF (IOLEN.EQ.0) GOTO 240
        IF (IOB160(1:1).EQ.'!') GOTO 240
        IF (IOB160(1:1).EQ.'.') GOTO 1
        IF (IOB160.EQ.' ') GO TO 240
        IF (IOB160(1:1).EQ.'='.OR.IOB160(1:2).EQ.' =')THEN
C  -     ATTEMPTED SPECIES INPUT AS 'PRODUCT' LIST
         I2=0
         WRITE(ICRT,248) IOB160
  248    FORMAT (' ERROR IN READING TYPEOMATIC LIST GIVEN BELOW.'
     &           /' ',A80)
         GOTO 240
        ENDIF
C
        CALL ALIN16(3,IOB160)
        CALL MOVLFT(IOB160)
        CALL TSPCNAM(I2,IOB160(1:16),.FALSE.)
        CALL TRDATOM(I2,iob160(17:32))
        READ (IOB160,'(32X,I16)',ERR=201) TTYPES(I2) 
        goto 240
C
c      READ IN THE SOM DEFINITIONS
c
  250   CALL INBUF ('SPECIE :',8)
        IF (IORET.EQ.2) GOTO 1999
        IF (IOLEN.EQ.0) GOTO 250
        IF (IOB160(1:1).EQ.'!') GOTO 250
        IF (IOB160(1:1).EQ.'.') GOTO 1
        IF (IOB160.EQ.' ') GO TO 250
        IF (IOB160(1:1).EQ.'='.OR.IOB160(1:2).EQ.' =')THEN
C  -     ATTEMPTED SPECIES INPUT AS 'PRODUCT' LIST
         I2=0
         WRITE(ICRT,258) IOB160
  258    FORMAT (' ERROR IN READING TYPEOMATIC LIST GIVEN BELOW.'
     &           /' ',A80)
         GOTO 250
        ENDIF
C
        CALL ALIN16(10,IOB160)
        CALL MOVLFT(IOB160)
        CALL SOMSPCNAM(I3,IOB160(1:16),.FALSE.)
        READ (IOB160,252)(SOMGRID(II,I3),II=1,4),(PFUNC(II,I3),II=1,4),
     &                    CFRAG(I3)
 252    FORMAT(16X,4(I16),4(F16.0),F16.0)
        print*,'SOM SPECIES:',SOMNAME(I3),(SOMGRID(II,I3),II=1,4),
     &          (PFUNC(II,I3),II=1,4),CFRAG(I3)
        if(somgrid(1,i3).lt.1 .or. somgrid(2,i3).gt.maxsomc)then
         print*,'bad som carbon grid specifications for :',somname(i3)
         print*,'entered:',somgrid(1,i3),somgrid(2,i3)
         print*,'min, max:',1,maxsomc
         print*,'reset maxsomc in pspecs.inc'
         stop 'rdrxn som marker1'
        endif
        if(somgrid(3,i3).lt.0 .or. somgrid(4,i3).gt.maxsomo)then
         print*,'bad som oxygen grid specifications for :',somname(i3)
         print*,'entered:',somgrid(3,i3),somgrid(4,i3)
         print*,'min, max:',1,maxsomo
         print*,'reset maxsomo in pspecs.inc'
         stop 'rdrxn som marker2'
        endif
        goto 250
C
C      DETERMINE TYPE OF INSERTED CODE
C
  300  IF (IOB160(9:16).EQ.'    INIT') THEN
             IV=4
       ELSEIF (IOB160(9:16).EQ.'    DIFF') THEN
             IV=5
       ELSE
             WRITE (OUT,311) IOB160(1:16)
  311        FORMAT (' UNRECOGNIZED TYPE FOR INSERTED CODE. INPUT=',
     &        A16,'.  COMMANDS FOR INSERTED CODE WILL BE IGNORED.')
             IF (OUT.NE.1) WRITE (ICRT,311) IOB160(1:16)
             IV=-1
             NERR=NERR+1
       ENDIF
C
C      READ IN VCOEF COMMANDS (FOR INSERTED CODE) AND SAVE IN
C      TEMPORARY FILE  (TYPE = IV)
C
  400  IF (IV.GT.0) THEN
           VCOCAL(IV)=.TRUE.
           IF (.NOT.VCODFD) THEN
               IOB160=' '
               CALL FILNAM (NLEN,IOB160,TMPUIC,'VCOPRP.TMP ',' ')
               OPEN (UNIT=UNIT1,NAME=IOB160,STATUS='UNKNOWN')
               VCODFD=.TRUE.
           ENDIF
       ENDIF
  410  CALL INBUF ('V.COEF :',8)
       IF (IORET.EQ.2) GOTO 1999
       IF (IOLEN.EQ.0) GOTO 410
       IF (IOB160(1:1).EQ.'!') GOTO 410
       IF (IOB160(1:1).EQ.'.') GOTO 1
       IF (IOB160(1:4).EQ.'    ') GOTO 410
       IF (NOPRP.OR.IV.EQ.-1) GOTO 410
       WRITE (UNIT1,411) IV,(BYT160(I),I=1,IOLEN)
  411  FORMAT (I1,128A1,128A1)
       GOTO 410
C
C
C      INPUT REACTION LIST (PROCESS .RXN COMMAND)
C
 1000  IND=6
       NPRODS(1)=1
       IF (.NOT.RXOPN) THEN
                IOB160=' '
                CALL FILNAM (NLEN,IOB160,TMPUIC,'PREP.TMP ',' ')
                OPEN (UNIT=UNIT3,NAME=IOB160,TYPE='UNKNOWN',
     &           FORM='UNFORMATTED')
                RXOPN=.TRUE.
                ENDIF
C
 1001  CALL INBUF ('REACTN :',8)
       IF (IORET.EQ.2) GOTO 1999
       IF (IOLEN.EQ.0) GOTO 1001
       IF (IOB160(1:1).EQ.'!') GOTO 1001
       IF (IOB160(1:1).EQ.'.') GOTO 1
C
        NRXN=NRXN+1
        NRXNLN=1
C       (NRXNLN GIVES NO. OF RXN LINES FOR CURRENT RXN)
        IF (NRXN.GT.MAXRXN) THEN
                I=MAXRXN
                WRITE (OUT,*) 'TOO MANY REACTIONS.  MAX =',I
                STOP 'TOO MANY REACTIONS'
                ENDIF
C
C  -     DECODE RXN PARMS
C        EXPECTED INPUT:  cccc) a,ea,b ;rxn
C                     OR  cccc: a,ea,b
C        WHERE cccc = REACTION LABEL, a = ARRHENIUS A FACTOR,
C        ea = ACTIVATION ENERGY (OR BLANK), b = B FACTOR FOR
C        TEMPERATURE DEPENDENCE (OR BLANK), rxn = REACTION STRING.
C        LEADING BLANKS OK.
C
        DUMRXN=.FALSE.
        RKBUF=' '
        J=0
C       GET PAST LEADING BLANKS
        DO 1003 I=1,IOLEN
        IF (BYT160(I).NE.' ') GOTO 1004
 1003   J=J+1
C       (RECORD IS TOTALLY BLANK)
        GOTO 1001
 1004   JP=J+1
C       (FIRST NON-BLANK CHAR)
c  the ";" separates the reaction rate info from the reaction  (ie reactants and products)
c
        DO 1005 L=JP,IOLEN
        IF (BYT160(L).EQ.';') GOTO 1006
 1005   RKBF1(L-J)=BYT160(L)
C       (COPY STRING UP TO ';' INTO RKBUF)
        DUMRXN=.TRUE.
C       (IF NO ';', THEN THIS IS 'DUMMY')
        NRTOSR(NRXN)=-1
c
c  L should be pointing to the ";" increment to get to the first reactant
c
 1006   L=L+1
        CALL ALIN16 (4,RKBUF)
        IF (RKBF1(16).NE.')')goto 1090
        RKBF1(16)=' '
C       (GET RID OF ')')
        CALL MOVLFT (RKBUF(1:16))
        LBL=RKBUF(1:6)
        RXNLBL(NRXN)=LBL
C
C       DETERMINE TYPE OF KINETICS
C
        IF (RKBUF(28:32).EQ.'CONST') THEN
C               CONSTANT REACTION 'CONST=k'
                ITYP=3
                LOCKBF=LOCKBF+1
                LKBUF(NRXN)=LOCKBF
                DECODE (48,1101,RKBUF,ERR=1090) KPBUF(LOCKBF)
 1101           FORMAT (32X,G16.0)
C
        ELSEIF (RKBUF(29:32).EQ.'  PF' .OR. RKBUF(29:32).EQ.'PHOT')
     &   THEN
C               PHOTOLYSIS REACTION 'PF=filename' OR  'PHOT=name'
                ITYP=7
                CALL MOVLFT (RKBUF(33:40))
                IF (NPHOTK.NE.0) THEN
                       DO 1011 I=1,NPHOTK
                       IF (RKBUF(33:40).EQ.PHOTNM(I)) THEN
                          LKBUF(NRXN)=I
                          GOTO 1014
                       ENDIF
 1011                  CONTINUE
                ENDIF
                NPHOTK=NPHOTK+1
                PHOTNM(NPHOTK)=RKBUF(33:40)
C               (IPHOTR IS FIRST RXN. NO OF THIS TYPE, WHICH IS THIS RXN)
                IPHOTR(NPHOTK)=NRXN
                LKBUF(NRXN)=NPHOTK
                PFTYPE(NPHOTK)=(RKBUF(29:32).EQ.'PHOT')
 1014           CONTINUE
C
        ELSEIF (RKBUF(28:32).EQ.'SAMEK') THEN
C               'SAMEK=rxnlbl' :  SAME RATE CONST AS ANOTHER RXN
                ITYP=0
                CALL MOVLFT (RKBUF(33:48))
                LBL2=RKBUF(33:38)
C
C               GET REACTION NO (I) FROM REACTION LABEL (LBL2)
C
                IF (NRXN.LE.1) GOTO 5011
                DO 5010 I=1,NRXN-1
                IF (LBL2.EQ.RXNLBL(I)) GOTO 5012
 5010           CONTINUE
 5011           WRITE (OUT,*) 'REACTION LABEL = ',LBL2,' NOT FOUND'
                GOTO 1090
 5012           IF (RXTYP(I).LT.3) THEN
                    WRITE (OUT,*) 'REF''D REACTION OF THE WRONG TYPE'
                    GOTO 1090
                ENDIF
                LKBUF(NRXN)=I
C
        ELSEIF (RKBUF(26:32).EQ.'FALLOFF') THEN
C               'FALLOFF' : TROE-TYPE FALLOFF USED TO DEFINE RATE CON.
C                       NEED TO READ PARMS LATER
                ITYP=5
                LKBUF(NRXN)=LOCKBF+1
                LOCKBF=LOCKBF+8
C
        ELSEIF (RKBUF(25:32).EQ.'K1+K2[M]') THEN
                ITYP=9
                LKBUF(NRXN)=LOCKBF+1
                LOCKBF=LOCKBF+7
                KPBUF(LOCKBF-6)=0 ! default rxn order for the P dependence
C
        ELSEIF (RKBUF(18:32).EQ.'K0+K3M/1+K3M/K2') THEN
                ITYP=10
                LKBUF(NRXN)=LOCKBF+1
                LOCKBF=LOCKBF+10
                KPBUF(LOCKBF-9)=0 ! default rxn order for the P dependence
C
        ELSE
C               DEFAULT:  READ A, EA, B
                ITYP=4
                LKBUF(NRXN)=LOCKBF+1
                LOCKBF=LOCKBF+4
                DECODE (64,1007,RKBUF,ERR=1090) KPBUF(LOCKBF-2),
     &           KPBUF(LOCKBF-1),KPBUF(LOCKBF)
                KPBUF(LOCKBF-3)=0 ! default rxn order for the P dependence
 1007           FORMAT (16X,4G16.0)
       ENDIF
C
c  now that up through the reaction has been processed, read in reactions and products
c     remember that in ASCII, [tab] = 9
c
 1019  IF (.NOT.DUMRXN) THEN
                J=0
                DO 1018 I=L,IOLEN
                ILCHAR=ICHAR(BYT160(I))
                IF (ILCHAR.NE.9) THEN
                        J=J+1
                        TMPRXN(J:J)=BYT160(I)
                ELSE
                        K=8-MOD(J+1,8)
C                       (ASSUMES ';' AT TAB POS)
                        TMPRXN(J+1:J+K)='        '
C                       (CONVERT TABS TO SPACES)
                        J=J+K
                ENDIF
 1018           CONTINUE
                IOB160=TMPRXN
                IOLEN=J
c
c   unit3 is a .TMP file used to develop the .mod file erased before prep exits
c
                IF (YESPRP) WRITE (UNIT3) NEWRXN,J,IOB160(1:IOLEN)
                OUTBUF=' '
                WRITE (OUTBUF,'(a)') LBL//' '//IOB160(1:iolen)
                LOUTBF=J+8
                IRTOSR(2,NRXN)=0
                IRTOSR(3,NRXN)=0
                INDEXP=NPRODS(NRXN)-1
C                 (LAST INDEX USED IN IPRODS)
                CALL RXLST1 (BYT160,NRTOSR(NRXN),INDEXP,NRXNLN)
                NPRODS(NRXN+1)=INDEXP+1
                EOF=IORET.EQ.2
                IF (NRTOSR(NRXN).EQ.0) GOTO 1092
        ELSE
C               (WRITE BLANK RXN STRING IF DUMMY RXN)
                IF (YESPRP) WRITE (UNIT3) NEWRXN,1,' '
                OUTBUF=' '
                WRITE (OUTBUF,1022) LBL
 1022           FORMAT (A6,') (NO REACTION)')
                LOUTBF=19
                IRTOSR(1,NRXN)=0
                IRTOSR(2,NRXN)=0
                IRTOSR(3,NRXN)=0
                NPRODS(NRXN+1)=NPRODS(NRXN)
                IORDR=0
        ENDIF
        RXTYP(NRXN)=ITYP
C
C       READ IN PARMS FOR 'FALLOFF' REACTION TYPE.
C
        IF (ITYP.EQ.5) THEN
                CALL INBUF ('FALLOFF:  K0 = ',15)
                IF (IORET.EQ.2) GOTO 1056
                CALL ALIN16 (3,IOB160)
                DECODE (48,204,IOB160,ERR=1057) KPBUF(LOCKBF-7),
     &           KPBUF(LOCKBF-6),KPBUF(LOCKBF-5)
  204           FORMAT (5G16.0)
                CALL INBUF ('FALLOFF:  KI = ',15)
                IF (IORET.EQ.2) GOTO 1056
                CALL ALIN16 (3,IOB160)
                DECODE (48,204,IOB160,ERR=1057) KPBUF(LOCKBF-4),
     &           KPBUF(LOCKBF-3),KPBUF(LOCKBF-2)
                CALL INBUF ('FALLOFF: F,N = ',15)
                IF (IORET.EQ.2) GOTO 1056
                CALL ALIN16 (2,IOB160)
                DECODE (32,204,IOB160,ERR=1057) KPBUF(LOCKBF-1),
     &           KPBUF(LOCKBF)
                IF (KPBUF(LOCKBF).EQ.0.0) KPBUF(LOCKBF)=1.0
        ELSEIF (ITYP.EQ.9) THEN
                CALL INBUF ('K1+K2[M]: K1 = ',15)
                CALL ALIN16 (4,IOB160)
                DECODE (48,204,IOB160,ERR=1090) KPBUF(LOCKBF-5),
     &           KPBUF(LOCKBF-4),KPBUF(LOCKBF-3)
                CALL INBUF ('K1+K2[M]: K2 = ',15)
                CALL ALIN16 (4,IOB160)
                DECODE (48,204,IOB160,ERR=1090) KPBUF(LOCKBF-2),
     &           KPBUF(LOCKBF-1),KPBUF(LOCKBF)
        ELSEIF (ITYP.EQ.10) THEN
                CALL INBUF ('K1+K3M/1+K3M/K2: K0 = ',22)
                CALL ALIN16 (4,IOB160)
                DECODE (48,204,IOB160,ERR=1090) KPBUF(LOCKBF-8),
     &           KPBUF(LOCKBF-7),KPBUF(LOCKBF-6)
                CALL INBUF ('K1+K3M/1+K3M/K2: K2 = ',22)
                CALL ALIN16 (4,IOB160)
                DECODE (48,204,IOB160,ERR=1090) KPBUF(LOCKBF-5),
     &           KPBUF(LOCKBF-4),KPBUF(LOCKBF-3)
                CALL INBUF ('K1+K3M/1+K3M/K2: K3 = ',22)
                CALL ALIN16 (4,IOB160)
                DECODE (48,204,IOB160,ERR=1090) KPBUF(LOCKBF-2),
     &           KPBUF(LOCKBF-1),KPBUF(LOCKBF)
        ENDIF
        GOTO 1060
C
C       ERROR READING ADDED PARMS
 1056   WRITE (OUT,*) '<EOF> WHEN ADDED KINETIC PARMS EXPECTED'
        GOTO 1059
 1057   WRITE (OUT,*) 'INPUT ERROR ADDED KINETIC PARMS EXPECTED'
 1059   NERR=NERR+1
C
C       DETERMINE IORDR = NUMBER OF NON-COEFFICIENT REACTANTS
C       AND CONVERT UNITS, IF NECESSARY
C
 1060   IORDR=0
        KEQ=.FALSE.
        IF (.NOT.PPM) GOTO 1030
        IF (ITYP.LT.3.OR.(ITYP.GT.5.AND.ITYP.LT.9)) GOTO 1030
        NRI=NRTOSR(NRXN)
        IF (NRI.LE.0) GOTO 1030
        DO 1032 J=1,NRI
        IR=IRTOSR(J,NRXN)
        IF (IR.GT.0) THEN
                IORDR=IORDR+1
        ELSEIF (-IR.GT.MAXMAX) THEN
C               (REFERENCES K OF ANOTHER REACTION, SO SUBTRACT OFF ORDER
C                OF THAT REACTION)
                KEQ=.TRUE.
                JRXN=-(IR+MAXMAX)
                IF (RXTYP(JRXN).EQ.7) THEN
                        IORDR=IORDR-1
C                       ( ALL PHOT RXNS ORDER=1)
                ELSE
                        NRJ=NRTOSR(JRXN)
                        DO 1034 I=1,NRJ
                        IR=IRTOSR(I,JRXN)
                        IF (IR.GT.0) THEN
                                IORDR=IORDR-1
                        ELSEIF (-IR.GT.MAXMAX) THEN
                                WRITE (OUT,*) 'CAN''T CONVERT UNITS'
     &                           ,'FOR THIS REACTION.'
                                WRITE (OUT,*) 'REF''S REACTION WHICH'
     &                           ,'REF''S ANOTHER'
                                GOTO 1090
                        ENDIF
 1034                   CONTINUE
                ENDIF
        ENDIF
 1032   CONTINUE
C  -     *NOTE*  CONVERSION FACTORS ARE FOR CONVERTING CMS UNITS INPUT
C                TO PPM-MIN UNITS.  VALID FOR 1 ATM PRESSURE ONLY.
C       CONV(1)=60.0
C       CONV(2)=4.40373E+17/TREF
C       CONV(3)=3.23206E+33/(TREF*TREF)
C
        IF (IORDR.EQ.0) GOTO 1030
        IF (KEQ) THEN
C               (EQUILIBRIUM CONSTANT)
                FAC=1.0
        ELSE
C               (RATE CONSTANT.  TAKE INTO ACCOUNT TIME UNITS AND CONC
C                UNITS = ORDER-1)
                FAC=60.0
                IORDR=IORDR-1
        ENDIF
        IF (IORDR.NE.0) FAC=FAC*(7.3395E+15/TREF)**(IORDR)
C
        IF (ITYP.EQ.3) THEN
                KPBUF(LKBUF(NRXN))=FAC*KPBUF(LKBUF(NRXN))
        ELSEIF (ITYP.EQ.4) THEN
                KPBUF(LOCKBF-2)=FAC*KPBUF(LOCKBF-2)
                KPBUF(LOCKBF)=KPBUF(LOCKBF)-FLOAT(IORDR)
                KPBUF(LOCKBF-3)=IORDR
        ELSEIF (ITYP.EQ.5) THEN
                KPBUF(LOCKBF-7)=FAC*(7.3395E+15/TREF)*KPBUF(LOCKBF-7)
                KPBUF(LOCKBF-5)=KPBUF(LOCKBF-5)-FLOAT(IORDR+1)
                KPBUF(LOCKBF-4)=FAC*KPBUF(LOCKBF-4)
                KPBUF(LOCKBF-2)=KPBUF(LOCKBF-2)-FLOAT(IORDR)
        ELSEIF (ITYP.EQ.9) THEN
                IF (KPBUF(LOCKBF-5).GT.0)THEN
                KPBUF(LOCKBF-5)=FAC*KPBUF(LOCKBF-5)
                KPBUF(LOCKBF-3)=KPBUF(LOCKBF-3)-FLOAT(IORDR)
                ENDIF
                IF (KPBUF(LOCKBF-2).GT.0)THEN
                KPBUF(LOCKBF-2)=FAC*(7.3395E+15/TREF)*KPBUF(LOCKBF-2)
                KPBUF(LOCKBF)=KPBUF(LOCKBF)-FLOAT(IORDR+1)
                ENDIF
                KPBUF(LOCKBF-6)=IORDR
        ELSEIF (ITYP.EQ.10) THEN
                IF (KPBUF(LOCKBF-8).GT.0)THEN
                KPBUF(LOCKBF-8)=FAC*KPBUF(LOCKBF-8)
                KPBUF(LOCKBF-6)=KPBUF(LOCKBF-6)-FLOAT(IORDR)
                ENDIF
                IF (KPBUF(LOCKBF-5).GT.0)THEN
                KPBUF(LOCKBF-5)=FAC*KPBUF(LOCKBF-5)
                KPBUF(LOCKBF-3)=KPBUF(LOCKBF-3)-FLOAT(IORDR)
                ENDIF
                IF (KPBUF(LOCKBF-2).GT.0)THEN
                KPBUF(LOCKBF-2)=FAC*(7.3395E+15/TREF)*KPBUF(LOCKBF-2)
                KPBUF(LOCKBF)=KPBUF(LOCKBF)-FLOAT(IORDR+1)
                ENDIF
                KPBUF(LOCKBF-9)=IORDR
        ENDIF
C
C       DONE WITH KINETIC PARM INPUT AND PROCESSING.  WRITE THEM OUT.
C
 1030   IF (ITYP.EQ.3) THEN
                WRITE (OUTBUF(73: ),1040) KPBUF(LKBUF(NRXN))
 1040           FORMAT (1PE11.3)
        ELSEIF (ITYP.EQ.7) THEN
                WRITE (OUTBUF(73: ),1041) RKBUF(33:40)
 1041           FORMAT (7X,'PHOT. = ',A8)
        ELSEIF (ITYP.EQ.0) THEN
                WRITE (OUTBUF(73: ),1042) RXNLBL(LKBUF(NRXN))
 1042           FORMAT (7X,'SAME K AS ',A6)
        ELSEIF (ITYP.EQ.1) THEN
                WRITE (OUTBUF(73: ),1043) RXNLBL(LKBUF(NRXN))
 1043           FORMAT (7X,'K = RATE OF ',A6)
        ELSEIF (ITYP.EQ.4) THEN
                WRITE (OUTBUF(73: ),1044) KPBUF(LOCKBF-2),
     &           KPBUF(LOCKBF-1),KPBUF(LOCKBF)
 1044           FORMAT (1PE11.3,0PF7.2,F8.3)
        ELSEIF (ITYP.EQ.5) THEN
                WRITE (OUTBUF(73: ),1045) KPBUF(LOCKBF-1),
     &           KPBUF(LOCKBF)
 1045           FORMAT (' FALLOFF F=',F6.3,', N=',F6.3)
        ELSEIF (ITYP.EQ.9) THEN
                WRITE (OUTBUF(73: ),1046) 
 1046           FORMAT (' K1+K2[M]')
        ELSEIF (ITYP.EQ.10) THEN
                WRITE (OUTBUF(73: ),1047)
 1047           FORMAT (' K0+K3[M]/(1+K3M/K2)')
        ELSE
                WRITE (OUT,*) 'PGM ERR.  BAD ITYP =',ITYP
        ENDIF
        LOUTBF=100
        IF (.NOT.DUMRXN) THEN
              NP=NPRODS(NRXN+1)-NPRODS(NRXN)
              IF (NP.LE.4) THEN
                  WRITE (OUTBUF(LOUTBF+1: ),1048) (IRTOSR(I,NRXN),I=1,3)
     &             ,(IPRODS(I),I=NPRODS(NRXN),NPRODS(NRXN+1)-1)
 1048             FORMAT (3I6,' =',4I6)
                  WRITE (OUT,1051) OUTBUF
 1051             FORMAT (' ',A160)
              ELSE
                  WRITE (OUTBUF(LOUTBF+1: ),1048) (IRTOSR(I,NRXN),I=1,3)
     &             ,(IPRODS(I),I=NPRODS(NRXN),NPRODS(NRXN)+3)
                  WRITE (OUT,1051) OUTBUF
                  WRITE (OUT,1050) (IPRODS(I),I=NPRODS(NRXN)+4,
     &              NPRODS(NRXN+1)-1)
 1050             FORMAT (115X,4I6)
              ENDIF
        ENDIF
        IF (ITYP.EQ.5) THEN
                WRITE (OUT,1049) (KPBUF(LOCKBF-I),I=7,2,-1)
 1049           FORMAT (70X,'K0:',1PE11.3,0PF7.2,F8.3/70X,'KI:',
     &           1PE11.3,0PF7.2,F8.3)
        ELSEIF (ITYP.EQ.9) THEN
                WRITE (OUT,1052) (KPBUF(LOCKBF-I),I=5,0,-1)
 1052           FORMAT (70X,'K1:',1PE11.3,0PF7.2,F8.3/
     &                  70X,'K2:',1PE11.3,0PF7.2,F8.3)
        ELSEIF (ITYP.EQ.10) THEN
                WRITE (OUT,1053) (KPBUF(LOCKBF-I),I=8,0,-1)
 1053           FORMAT (70X,'K0:',1PE11.3,0PF7.2,F8.3/
     &                  70X,'K2:',1PE11.3,0PF7.2,F8.3/
     &                  70X,'K3:',1PE11.3,0PF7.2,F8.3)
        ENDIF
        IF (LOCKBF.GT.MAXKBF) THEN
                WRITE (OUT,*) 'KINETIC PARM BUFFER LENGTH =',LOCKBF
     &           ,' EXCEEDS PROGRAM DIMENSIONS OF ',MAXKBF
                STOP 'TOO MANY RATE PARMS'
        ENDIF
        GO TO 1001
C
C  -     ERROR
 1090  WRITE (OUT,1091) IOB160(1:IOLEN)
 1091  FORMAT (' REACTION PARM INPUT ERROR IN LINE GIVEN BELOW.'
     &  /' ',A)
       NERR=NERR+1
       WRITE (ICRT,1091) IOB160(1:IOLEN)
       GOTO 1095
 1092  WRITE (ICRT,1094) IOB160(1:IOLEN)
 1094  FORMAT (' REACTION LIST DECODING ERROR.  INPUT LINE BELOW:'
     &  /' ',A)
       IF (YESPRP) THEN
               DO 1096 I=1,NRXNLN
 1096          BACKSPACE UNIT3
       ENDIF
 1095  NRXN=NRXN-1
       IF (EOF) GOTO 1999
       GOTO 1001
C
C
C     PROCESS VARIABLE COEFFICIENT INPUT DATA (.COE COMMAND)
C     VALUE OF VARIABLE COEFFICIENT READ IN    (DEFAULT VALUE = 1.0)
C
 1200  CALL INBUF ('COEFNT :',8)
       IF (IORET.EQ.2) GOTO 1999
       IF (IOLEN.EQ.0) GOTO 1200
       IF (IOB160(1:1).EQ.'!') GOTO 1200
       IF (IOB160(1:1).EQ.'.') GOTO 1
       IF (IOB160(1:8).EQ.'        ') GOTO 1200
       CALL ALIN16(2,IOB160)
       CALL MOVLFT(IOB160)
       DECODE (16,204,IOB160(17:32),ERR=1290) X
       IF (NCOEFV.EQ.0) GOTO 1206
       DO 1204 I=1,NCOEFV
       IF (IOB160(1:8).EQ.COEFNM(I)) GOTO 1205
 1204  CONTINUE
 1206  NCOEFV=NCOEFV+1
       IF (NCOEFV.GT.MAXCOV) THEN
                I=MAXCOV
                WRITE (OUT,*) 'TOO MANY VARIABLE COEFFICIENTS.  MAX =',I
                STOP 'TOO MANY VARIABLE COEFFICIENTS'
                ENDIF
       COEFNM(NCOEFV)=IOB160(1:8)
       COEF(NCOEFV)=X
       GOTO 1200
 1205  COEF(I)=X
       GOTO 1200
 1290  WRITE (OUT,1291) IOB160
 1291  FORMAT (' COEF INPUT ERROR IN LINE GIVEN BELOW.'/' ',A80)
       NERR=NERR+1
       WRITE (ICRT,1291) IOB160
       GOTO 1200
C
C
C      PROCESS PHOTOLYSIS ABSORPTION COEFFICIENT, QUANTUM YIELD INPUT
C      (.PHOT COMMAND).  THESE DATA ARE WRITTEN TO SEPERATE TEMPORARY
C      DATA SETS, WHICH ARE READ LATER BY SAVERX AND PROCESSED THERE
C      BEFORE BEING OUTPUT TO '.MOD' FILE.
C
 1500  IOBF16=' '
C      (GET NAME TO ASSOCIATE WITH THIS SET OF ABS QY'S)
       IOBF16=IOB160(9:16)
       CALL MOVLFT (IOBF16)
       IOB160=' '
       CALL FILNAM (I,IOB160,TMPUIC,IOBF16,'.PTM ')
       OPEN (UNIT=UNIT2,NAME=IOB160,TYPE='UNKNOWN')
 1501  CALL INBUF ('PHOT   :',8)
       IF (IORET.EQ.2) THEN
              CLOSE (UNIT=UNIT2)
              GOTO 1999
       ENDIF
       IF (IOLEN.EQ.0 .OR. IOB160.EQ.' ') THEN
               CLOSE (UNIT=UNIT2)
               GOTO 100
       ENDIF
       IF (IOB160(1:1).EQ.'.') THEN
               CLOSE (UNIT=UNIT2)
               GOTO 1
       ENDIF
       WRITE (UNIT2,1510) IOB160
 1510  FORMAT (A80)
       GOTO 1501
C
C
C      END OF INPUT:  .END STATEMENT OR END-OF-FILE ENCOUNTERED.
C      THIS VERSION REQUIRES AN ".END" STATEMENT FOR NORMAL END OF
C      PROCESSING, SO IF END-OF-FILE ENCOUNTERED FIRST, IT IS TREATED
C      AS AN ERROR.  CONTROL GOES TO 1999 IF END-OF-FILE ENCOUNTERED,
C      OR TO 2000 IF ".END" RECORD PROCESSED.
C
 1999  EOF=.TRUE.
       WRITE (OUT,*) ' ERROR -- UNEXPECTED END-OF-FILE.'
       NERR=NERR+1
C
 2000  IF (.NOT.RXOPN) STOP 'PREP ABORTED.  NO REACTIONS.'
C
C       CLOSE AND REWIND TEMPORARY DATA SETS
C
        IF (YESPRP) THEN
                IF (VCODFD) CLOSE (UNIT=UNIT1,STATUS='KEEP')
                WRITE (UNIT3) NEWRXN,1,' '
C               (NEED DUMMY RXN LINE AT END FOR SAVERX TO READ FILE OK)
                CLOSE (UNIT=UNIT3,STATUS='KEEP')
        ENDIF
C
C       WRITE NUMBERS OF REACTIONS, SPECIES, COEFS, ETC.  ABORT IF
C       BOUNDS EXCEEDED
C
        WRITE (OUT,2001) NRXN,NS,NCOEFV,NCOC2-NCOC1+1,NPHOTK
 2001   FORMAT ('0',I4,' REACTIONS,',I5,' SPECIES,',I5,' VARIABLE',
     &   ' COEFFICIENTS,',I4,' CONSTANT COEFFICIENTS, AND',I4,
     &   ' PHOTOLYSIS FILES.')
C
        IF (NRXN.GT.MAXRXN) STOP 'TOO MANY REACTIONS'
        IF (NS.GT.MAXNS) STOP 'TOO MANY SPECIES'
        IF (NCOC2.GT.MAXCO) STOP 'TOO MANY CONSTANT COEFFICIENTS rdrxn'
        IF (NCOEFV.GT.MAXCOV) STOP 'TOO MANY VARIABLE COEFFICIENTS'
        IF (NPHOTK.GT.MAXPHK) STOP 'TOO MANY PHOT. FILES'
C
        RETURN
C
        END
