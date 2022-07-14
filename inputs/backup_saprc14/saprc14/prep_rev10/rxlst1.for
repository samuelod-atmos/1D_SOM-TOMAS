C       FILE RXLST1.FTN                 PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATED:   W.P.L. CARTER 6/20/87
C
C
        SUBROUTINE RXLST1 (CHARXN,IR,IP,NLINES)
C
C       PROCESSES A CHARACTER STRING DESCRIBING REACTANTS, PRODUCTS AND
C       COEFFICIENTS FOR A REACTION.  CREATS NEW SPECIES AND COEFFICIENTS
C       AS NEEDED.  DOES EXTENSIVE ERROR CHECKING FOR APPROPRIATE REACTION
C       INPUT FORMAT.
C
C       CALLED BY:  RDRXN
C
C       SUBROUTINES CALLED (ALL IN FILE RXSUBS.FTN)
C
C       COFNAM  ... CREATES NEW COEFFICIENT GIVEN COEFFICIENT NAME.
C                   DETERMINES WHETHER IT IS VARIABLE OR CONSTANT,
C                   DEPENDING ON WHETHER ITS NAME CAN BE DECODED AS A
C                   NUMBER.
C       CONTLN  ... READS IN CONTINUATION LINE FOR REACTION
C       SPCNAM  ... DETERMINES SPECIES NUMBER FROM NAME.  DEFINES TYPE
C                   BASED ON CURRENT VALUE OF "IND".  CREATES NEW SPECIES
C                   IF NEEDED.
C
C
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL VARIABLES.  IR IS USED FOR NO. OF REACTANTS, IP AS INDEX
C       INTO PRODUCT ARRAY.
C
        INTEGER   IR,IP
        LOGICAL   COIND,CNEW,RCTNT,EOLN,RX
        CHARACTER*1 CHARXN(1),LCHAR
C
C
C       INITIALIZE
C
       I=0
       IR=0
       IP1=IP
       RX=.TRUE.
       CNEW=.FALSE.
       EOLN=.FALSE.
       LASTCO=0
C
C      NEW SPECIES OR COEF IN REACTANT LIST
C
 1000  J=0
       COIND=.FALSE.
       IOBF16=' '
C  -     EXAMINE NEXT CHAR
 1001  I=I+1
       IF (I.GT.IOLEN) GOTO 9001
C                                       !ERROR <EOL> ENCOUNTERED BEFORE '='
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
       IF (LCHAR.EQ.' ') GOTO 1001
 1005  IF (LCHAR.EQ.'+') GOTO 9002
C                                       !ERROR WHEN RCTNT OR COEF EXPECTED
       IF (LCHAR.EQ.'=') GOTO 9002
C                                       !ERROR WHEN RCTNT OR COEF EXPECTED
       IF (LCHAR.EQ.'"') GOTO 9002
C                                       !ERROR WHEN RCTNT OR COEF EXPECTED
       IF (LCHAR.EQ.'!') GOTO 9002
C                                       !ERROR WHEN RCTNT OR COEF EXPECTED
C
       IF (LCHAR.EQ.'#') THEN
                COIND=.TRUE.
                GOTO 1201
C                                       !GET NAME OF VARIABLE COEFFICIENT
       ELSEIF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1005
       ENDIF
C
C===================================
C       FINISH DECODING NAME OF A REACTANT OR VARIABLE COEFFICIENT
C
 1200  IF (J.GE.16) GOTO 9020
C                                               !ERROR. SPECIES NAME TOO LONG
       IF (J.GE.8 .AND. COIND) GOTO 9010
C                                               !ERROR. COEF NAME TOO LONG
       J=J+1
       BYT16(J)=LCHAR
C  -     EXAMINE NEXT CHAR
 1201  I=I+1
       IF (I.GT.IOLEN) GOTO 9001
C                                       !ERROR <EOL> ENCOUNTERED IN RCTNT LIST
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
 1205  IF (LCHAR.EQ.'"') GOTO 9007
C                                       !ERROR WHEN DECODING RCTNT OR COEF NAME
       IF (LCHAR.EQ.'!') GOTO 9007
C                                       !ERROR WHEN DECODING RCTNT OR COEF NAME
       IF (LCHAR.EQ.'#') GOTO 9007
C                                       !ERROR WHEN DECODING RCTNT OR COEF NAME
       IF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1205
       ENDIF
       IF (LCHAR.NE.'+'.AND.LCHAR.NE.'='.AND.LCHAR.NE.' ') GOTO 1200
C                                       !CONTINUE DECODING SPECIES OR COEF NAME
C  -   PROCESS SPECIES OR COEF
        IF (COIND) THEN
                IF (J.EQ.0) GOTO 9004
C                                               !MISSING COEF NAME AFTER '#'
                RCTNT=.TRUE.
                CALL COFNAM(IS,RCTNT)
                IF (IS.EQ.0) GOTO 9990
C                                               !ERROR IN COFNAM
        ELSE
                CALL SPCNAM(IS,IOBF16,RX)
        ENDIF
C
C===================================
C     ADD SPECIES OR COEF AS A REACTANT
C
        IR=IR+1
        IF (IR.GT.MAXRCT) GOTO 9003
C                                               !ERROR. TOO MANY REACTANTS
        IRTOSR(IR,NRXN)=IS
        IF (LCHAR.EQ.'+') GOTO 1000
C                                       !GET NEXT REACTANT OR COEFFICIENT
        IF (LCHAR.EQ.'=') GOTO 1499
C                                       !NOW PROCESS PRODUCT LIST
C
C===================================
C      LCHAR=' ' AFTER SPECIES OR COEF IN RCTNT LIST.
C      NOW LOOK FOR '+' OR '='.
C
 1210  I=I+1
       IF (I.GT.IOLEN) GOTO 9001
C                                       !ERROR. <EOL> ENCOUNTERED BEFORE '='
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
       IF (LCHAR.EQ.' ') GOTO 1210
 1215  IF (LCHAR.EQ.'+') GOTO 1000
C                                               !GET NEXT REACTANT OR COEF
       IF (LCHAR.EQ.'=') GOTO 1499
C                                               !NOW PROCESS PRODUCT LIST
       IF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1215
       ELSE
                GOTO 9008
C                                       !ERROR WHEN '=' OR '+' EXPECTED
       ENDIF
C
C===================================
C      NEW SPECIES OR COEF IN PRODUCT LIST
C
 1499  RCTNT=.FALSE.
 1500  J=0
       COIND=.FALSE.
       IOBF16=' '
C  -     EXAMINE NEXT CHAR
 1501  I=I+1
       IF (I.GT.IOLEN) THEN
                IF (LASTCO.NE.0) THEN
                        GOTO 9011
C                                       !ERROR. <EOL> WITHOUT CLOSING '"'
                ELSEIF (IP.EQ.IP1) THEN
                        RETURN
                ELSE
                        GOTO 9005
C                                       !ERROR. PRODUCT OR COEF EXPECTED
                ENDIF
       ENDIF
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
       IF (LCHAR.EQ.' ') GOTO 1501
 1505  IF (LCHAR.EQ.'+') GOTO 9012
C                                       !ERROR WHEN PRODUCT OR COEF EXPECTED
       IF (LCHAR.EQ.'=') GOTO 9012
C                                       !ERROR WHEN PRODUCT OR COEF EXPECTED
       IF (LCHAR.EQ.'"') GOTO 9012
C                                       !ERROR WHEN PRODUCT OR COEF EXPECTED
       IF (LCHAR.EQ.'!') GOTO 9012
C                                       !ERROR WHEN PRODUCT OR COEF EXPECTED
C
       IF (LCHAR.EQ.'#') THEN
                COIND=.TRUE.
                GOTO 1801
C                                       !GET NAME OF VARIABLE COEFFICIENT
       ELSEIF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1505
       ENDIF
C
C===================================
C      FINISH DECODING NAME OF A PRODUCT SPECIE
C
 1600  IF (J.GE.16) GOTO 9020
C                                       !ERROR. SPECIES NAME TOO LONG
       J=J+1
       BYT16(J)=LCHAR
C  -     EXAMINE NEXT CHAR
 1601  I=I+1
       IF (I.GT.IOLEN) THEN
                IF (LASTCO.NE.0) GOTO 9011
C                                               !<EOL> WITHOUT CLOSING '"'
                EOLN=.TRUE.
                GOTO 1700
C                                               !PROCESS PRODUCT
       ENDIF
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
 1605  IF (LCHAR.EQ.'#') GOTO 9013
C                                       !ERROR WHEN DECODING PRODUCT NAME
       IF (LCHAR.EQ.'=') GOTO 9013
C                                       !ERROR WHEN DECODING PRODUCT NAME
       IF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1605
       ENDIF
       IF (LCHAR.EQ.'"') THEN
C                                               !THIS SHOULD BE A CLOSING '"'
                IF (LASTCO.EQ.0) GOTO 9014
C                                               !ERROR. MISPLACED '"'
                GOTO 1700
       ENDIF
       IF (LCHAR.EQ.'+'.OR.LCHAR.EQ.'!') GOTO 1700
       IF (LCHAR.EQ.' ') GOTO 1700
       GOTO 1600
C                                               !CONTINUE DECODING NAME
C
C===================================
C      ADD SPECIES TO PRODUCT LIST. (ADD COEF FIRST IF APPLICABLE.)
C
 1700  CALL SPCNAM(IS,IOBF16,RX)
       IF (LASTCO.NE.0) THEN
                IP=IP+1
                IF (IP.GT.MAXPRD) GOTO 9017
C                                               !TOO MANY PRODUCTS IN MODEL
                IPRODS(IP)=LASTCO
       ENDIF
       IP=IP+1
       IF (IP.GT.MAXPRD) GOTO 9017
C                                               !TOO MANY PRODUCTS IN MODEL
       IPRODS(IP)=IS
       IF (EOLN) THEN
                RETURN
       ENDIF
       IF (LCHAR.EQ.'+') GOTO 1500
C                                               !GET NEXT SPECIES OR COEF
       IF (LCHAR.EQ.'!') THEN
                IF (LASTCO.NE.0) THEN
                        GOTO 9018
C                                               !ERROR. MISSING CLOSING '"'
                ELSE
                        RETURN
                ENDIF
       ELSEIF (LCHAR.EQ.'"') THEN
                LASTCO=0
       ENDIF
C
C===================================
C      LCHAR='"' OR ' ' AFTER SPECIES NAME IN PROD LIST.
C      NOW  LOOK FOR '+', '!' OR <EOL>
C
 1710  I=I+1
       IF (I.GT.IOLEN) THEN
                RETURN
       ENDIF
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
       IF (LCHAR.EQ.' ') GOTO 1710
 1715  IF (LCHAR.EQ.'+') GOTO 1500
C                                               !GET NEXT PRODUCT OR COEF
       IF (LCHAR.EQ.'!') THEN
                RETURN
       ENDIF
       IF (LCHAR.EQ.'"') THEN
                IF (LASTCO.NE.0) THEN
                        LASTCO=0
                        GOTO 1710
                ELSE
                        GOTO 9014
C                                               !ERROR. MISPLACED '"'
                ENDIF
       ELSEIF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1715
       ELSE
                GOTO 9021
C                                       !UNEXPECTED CHAR AFTER PRODUCT SPECIE
       ENDIF
C
C===================================
C       FINISH DECODING NAME OF A VARIABLE COEF IN PRODUCT LIST
C
C  -     EXAMINE NEXT CHAR
 1801  I=I+1
       IF (I.GT.IOLEN) GOTO 9019
C                                       !<EOL> DECODING COEF NAME
       IF (CNEW) THEN
               LCHAR=BYT160(I)
       ELSE
               LCHAR=CHARXN(I)
       ENDIF
1805   IF (LCHAR.EQ.'!') GOTO 9015
C                                       !MISSING SPECIES AFTER VAR COEF
       IF (LCHAR.EQ.'#') GOTO 9016
C                                       !ERROR DECODING COEF NAME IN PROD LIST
       IF (LCHAR.EQ.'+') GOTO 9015
C                                       !MISSING SPECIES AFTER VAR COEF
       IF (LCHAR.EQ.'&') THEN
                CALL CONTLN(I,LCHAR,NLINES)
C                                               !READ IN CONTINUATION LINE
                IF (I.EQ.0) GOTO 9998
C                                               !ERROR CONDITION
                CNEW=.TRUE.
                GOTO 1805
       ENDIF
       IF (LCHAR.NE.' '.AND.LCHAR.NE.'"') THEN
                IF (J.GE.8) GOTO 9010
C                                               !COEF NAME TOO LONG
                J=J+1
                BYT16(J)=LCHAR
                GOTO 1801
       ENDIF
       IF (J.EQ.0) GOTO 9004
C                                               !MISSING COEF NAME AFTER '#'
C
C===================================
C      LCHAR='"' AFTER COEF NAME.  ENCLOSES LIST OF SPECIES WITH SAME COEF.
C        (I.E. #C"A + B + C" + D   IS SAME AS #C A +#C B+#C C + D)
C      SET LASTCO TO COEF #.  DON'T ADD COEF TO PRODUCT LIST YET.
C
 1810  IF (LCHAR.EQ.'"') THEN
                IF (LASTCO.NE.0) THEN
                        GOTO 9014
C                                               !ERROR. MISPLACED '"'
                ELSE
                        CALL COFNAM(IS,RCTNT)
                        IF (IS.EQ.0) GOTO 9990
C                                               !ERROR IN COFNAM
                        LASTCO=IS
                        GOTO 1500
C                                               !GET NEXT SPECIES OR COEF
                ENDIF
       ENDIF
C
C===================================
C      LCHAR=' ' AFTER COEF NAME IN PROD LIST.
C      NOW LOOK FOR '"' OR BEGINNING OF NEXT SPECIES OR COEF NAME.
C      IF NOT '"', ADD COEF TO PRODUCT LIST.
C
       IF (LCHAR.EQ.' ') THEN
 1901           I=I+1
                IF (I.GT.IOLEN) GOTO 9019
C                                               !<EOL> DECODING COEF NAME
                IF (CNEW) THEN
                        LCHAR=BYT160(I)
                ELSE
                        LCHAR=CHARXN(I)
                ENDIF
                IF (LCHAR.EQ.' ') GOTO 1901
 1905           IF (LCHAR.EQ.'!') GOTO 9015
C                                               !MISSING SPECIES AFTER VAR COEF
                IF (LCHAR.EQ.'+') GOTO 9015
C                                               !MISSING SPECIES AFTER VAR COEF
                IF (LCHAR.EQ.'&') THEN
                        CALL CONTLN(I,LCHAR,NLINES)
C                                                       !READ IN CONT. LINE
                        IF (I.EQ.0) GOTO 9998
C                                                       !ERROR CONDITION
                        CNEW=.TRUE.
                        GOTO 1905
                ENDIF
                IF (LCHAR.EQ.'"') THEN
                        GOTO 1810
                ELSEIF (LCHAR.EQ.'#') THEN
                        CALL COFNAM(IS,RCTNT)
                        IF (IS.EQ.0) GOTO 9990
C                                                       !ERROR IN COFNAM
                        IP=IP+1
                        IF (IP.GT.MAXPRD) GOTO 9017
C                                                       !TOO MANY PRODUCTS
                        IPRODS(IP)=IS
                        J=0
                        COIND=.TRUE.
                        IOBF16=' '
                        GOTO 1801
C                                                       !DECODE NEXT COEF NAME
                ELSE
                        CALL COFNAM(IS,RCTNT)
                        IF (IS.EQ.0) GOTO 9990
C                                                       !ERROR IN COFNAM
                        IP=IP+1
                        IF (IP.GT.MAXPRD) GOTO 9017
C                                                       !TOO MANY PRODUCTS
                        IPRODS(IP)=IS
C                  ALREADY HAVE FIRST CHAR OF SPECIES NAME IN LCHAR
                        J=0
                        COIND=.FALSE.
                        IOBF16=' '
                        GOTO 1600
C                                               !CONTINUE DECODING SPECIES NAME
                ENDIF
       ENDIF
C
C====================================
C      ERRORS
C
 9001  WRITE (OUT,9101)
 9101  FORMAT (' ERROR:  <EOL> ENCOUNTERED IN RCTNT LIST BEFORE ''=''')
       GOTO 9998
C
 9002  WRITE (OUT,9102)
 9102  FORMAT (' ERROR:  ILLEGAL CHAR WHEN REACTANT OR COEF EXPECTED')
       GOTO 9990
C
 9003  WRITE (OUT,9103) MAXRCT
 9103  FORMAT (' ERROR: TOO MANY REACTANTS.  MAX =',I2)
       GOTO 9990
C
 9004  WRITE (OUT,9104)
 9104  FORMAT (' ERROR:  MISSING VAR COEF NAME AFTER ''#'' IN RXN',
     &    ' LIST')
       GOTO 9990
C
 9005  WRITE (OUT,9105)
 9105  FORMAT (' ERROR:  <EOL> WHEN PRODUCT SPECIES OR COEF EXPECTED')
       GOTO 9998
C
 9007  IF (COIND) GOTO 9016
       WRITE (OUT,9107) IOBF16
 9107  FORMAT (' ERROR:  ILLEGAL CHAR IN REACTANT NAME.  FIRST PART=',
     &   A16)
       GOTO 9990
C
 9008  WRITE (OUT,9108) IOBF16
 9108  FORMAT (' ERROR:  ''='' OR ''+'' EXPECTED AFTER BLANK.  ',
     &   'LAST SPECIES OR COEF=',A16)
       GOTO 9990
C
 9010  WRITE (OUT,9110) IOBF16(1:8)
 9110  FORMAT (' ERROR:  VAR COEF NAME TOO LONG.  FIRST PART=',A8)
       GOTO 9990
C
 9011  WRITE (OUT,9111)
 9111  FORMAT (' ERROR:  <EOL> WITHOUT CLOSING ''"'' IN PRODUCT LIST')
       GOTO 9998
C
 9012  WRITE (OUT,9112)
 9112  FORMAT (' ERROR:  ILLEGAL CHAR WHEN PRODUCT OR COEF EXPECTED')
       GOTO 9990
C
 9013  WRITE (OUT,9113) IOBF16
 9113  FORMAT (' ERROR:  ILLEGAL CHAR IN PRODUCT NAME.  FIRST PART=',
     &   A16)
       GOTO 9990
C
 9014  WRITE (OUT,9114)
 9114  FORMAT (' ERROR:  MISPLACED ''"''')
       GOTO 9990
C
 9015  WRITE (OUT,9115) IOBF16(1:8)
 9115  FORMAT (' ERROR:  MISSING SPECIES AFTER VAR COEF #',A8)
       GOTO 9990
C
 9016  WRITE (OUT,9116) IOBF16(1:8)
 9116  FORMAT (' ERROR:  ILLEGAL CHAR IN COEF NAME.  FIRST PART=#',A8)
       GOTO 9990
C
 9017  WRITE (OUT,9117) MAXPRD
 9117  FORMAT (' ERROR:  TOO MANY PRODUCTS IN MODEL.  MAX=',I4)
       STOP 'TOO MANY PRODUCT SPECIES IN MODEL'
C
 9018  WRITE (OUT,9118)
 9118  FORMAT (' ERROR:  MISSING CLOSING ''"''')
       GOTO 9990
C
 9019  WRITE (OUT,9119)
 9119  FORMAT (' ERROR:  <EOL> WHEN EXPECTING SPECIES AFTER VAR COEF')
       GOTO 9998
C
 9020  WRITE (OUT,9120) IOBF16
 9120  FORMAT (' ERROR:  SPECIES NAME TOO LONG.  FIRST PART=',A16)
       GOTO 9990
C
 9021  WRITE (OUT,9121) IOBF16
 9121  FORMAT (' ERROR:  UNEXPECTED CHAR AFTER PRODUCT SPECIES ',A16)
       GOTO 9990
C
 9990  WRITE (OUT,9991) I,LCHAR
 9991  FORMAT (' ABOVE ERROR DISCOVERED AT CHARACTER #',I4,
     &   '.  CHAR = ',A1)
 9998  WRITE (OUT,9999) OUTBUF(1:80)
 9999  FORMAT (' REACTION WITH ERROR:'/' ',A80)
       IR=0
       IP=IP1
       NERR=NERR+1
       RETURN
       END
