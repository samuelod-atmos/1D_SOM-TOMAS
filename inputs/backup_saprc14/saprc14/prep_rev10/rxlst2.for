C       FILE: RXLST2.FTN                        PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATED:  W.P.L. CARTER  6/10/87
C
C
        SUBROUTINE RXLST2 (CHARXN,IR)
C
C       PROCESSES A LIST OF SPECIES, GIVING THEM A SPECIES TYPE AS INDICATED
C       BY CURRENT VALUE OF "IND" VARIABLE.   CREATES NEW SPECIES AS NEEDED.
C       DOES SOME ERROR CHECKING
C
C      INPUT EXPECTED TO HAVE THE FORMAT:
C        '=(species 1) + (species 2) + etc.'
C
C
C       CALLED BY:  RDRXN
C
C       SUBROUTINES CALLED:
C
C       SPCNAM  ... DETERMINES SPECIES NUMBER FROM NAME.  DEFINES TYPE
C                   BASED ON CURRENT VALUE OF "IND".  CREATES NEW SPECIES
C                   IF NEEDED.
C
C
C       INCLUDES SPECIFICATIONS OF COMMON VARIABLES, PARMS, AND ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL SPECIFICATIONS
C
        CHARACTER*1 CHARXN(1),LCHAR
        LOGICAL   RX
C
C
C       INITIALIZE
C
       I=0
       IR=0
       RX=.FALSE.
C
C  -   LOOK FOR '=' FIRST
C
 100   I=I+1
       IF (I.GT.IOLEN) GOTO 9001
C                                    !BLANK LINE WHEN LIST OF SPECIES EXPECTED
       LCHAR=CHARXN(I)
       ILCHAR=ICHAR(LCHAR)
       IF (LCHAR.EQ.' '.OR.ILCHAR.EQ.9) THEN
C                                               !BLANK OR TAB
                GOTO 100
       ELSEIF (LCHAR.EQ.'=') THEN
                GOTO 1000
       ELSE
                J=0
                IOBF16=' '
                GOTO 1005
C                               !LCHAR SHOULD BE 1ST CHAR OF NAME IF NOT '='
       ENDIF
C
C      NEW SPECIES IN LIST
C
 1000  J=0
       IOBF16=' '
C  -     EXAMINE NEXT CHAR
 1001  I=I+1
       IF (I.GT.IOLEN) GOTO 9003
C                                       !ERR: <EOL> WHEN SPECIES NAME EXPECTED
       LCHAR=CHARXN(I)
       ILCHAR=ICHAR(LCHAR)
       IF (LCHAR.EQ.' '.OR.ILCHAR.EQ.9) GOTO 1001
 1005  IF (LCHAR.EQ.'+') GOTO 9004
C                                       !ERROR. SPECIES NAME EXPECTED
       IF (LCHAR.EQ.'!') GOTO 9004
C                                       !ERROR. SPECIES NAME EXPECTED
       IF (LCHAR.EQ.'=') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'"') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'#') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'&') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
C
C       VALID CHARACTER.  FINISH DECODING NAME OF A REACTANT
C
 1200  IF (J.GE.16) GOTO 9006
C                                       !ERROR. SPECIES NAME TOO LONG
       J=J+1
       BYT16(J)=LCHAR
C  -     EXAMINE NEXT CHAR
 1201  I=I+1
       IF (I.GT.IOLEN) THEN
                CALL SPCNAM(IDUM,IOBF16,RX)
                IR=IR+1
                RETURN
       ENDIF
       LCHAR=CHARXN(I)
       ILCHAR=ICHAR(LCHAR)
 1205  IF (LCHAR.EQ.'"') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'#') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'&') GOTO 9005
C                                       !ERROR. ILLEGAL CHARACTER
       IF (LCHAR.EQ.'+') THEN
                CALL SPCNAM(IDUM,IOBF16,RX)
                IR=IR+1
                GOTO 1000
C                                       !GET NEXT SPECIES NAME
       ELSEIF (LCHAR.EQ.' '.OR.ILCHAR.EQ.9) THEN
                CALL SPCNAM(IDUM,IOBF16,RX)
                IR=IR+1
C  -         LOOK FOR '+' OR '!' OR <EOL> NEXT
 1210           I=I+1
                IF (I.GT.IOLEN) RETURN
                LCHAR=CHARXN(I)
                ILCHAR=ICHAR(LCHAR)
                IF (LCHAR.EQ.' '.OR.ILCHAR.EQ.9) THEN
                        GOTO 1210
                ELSEIF (LCHAR.EQ.'!') THEN
                        RETURN
                ELSEIF (LCHAR.EQ.'+') THEN
                        GOTO 1000
C                                               !GET NEXT SPECIES NAME
                ELSE
                        GOTO 9007
C                                               !ERROR. '+' EXPECTED
                ENDIF
       ELSEIF (LCHAR.EQ.'!') THEN
                CALL SPCNAM(IDUM,IOBF16,RX)
                IR=IR+1
                RETURN
       ELSE
                GOTO 1200
C                                       !CONTINUE DECODING SPECIES NAME
       ENDIF
C
C  ------ ERRORS ------
C
 9001  WRITE (OUT,9101)
 9101  FORMAT (' ERROR:  BLANK LINE WHEN LIST OF REACTANTS EXPECTED')
       GOTO 9998
C
 9002  WRITE (OUT,9102)
 9102  FORMAT (' ERROR:  ''='' EXPECTED AT BEGINNING OF REACTANT LIST')
       GOTO 9990
C
 9003  WRITE (OUT,9103) IOB160
 9103  FORMAT (' ERROR:  <EOL> ENCOUNTERED WHEN SPECIES NAME EXPECTED.',
     &   '  INPUT LINE BELOW:'/' ',A80)
       GOTO 9998
C
 9004  WRITE (OUT,9104)
 9104  FORMAT (' ERROR:  ILLEGAL CHAR WHEN SPECIES NAME EXPECTED')
       GOTO 9990
C
 9005  WRITE (OUT,9105)
 9105  FORMAT (' ERROR:  ILLEGAL CHAR IN REACTANT LIST')
       GOTO 9990
C
 9006  WRITE (OUT,9106) IOBF16
 9106  FORMAT (' ERROR:  SPECIES NAME TOO LONG.  FIRST PART=',A16)
       GOTO 9990
C
 9007  WRITE (OUT,9107) IOBF16
 9107  FORMAT (' ERROR:  ''+'' EXPECTED AFTER BLANK.  ',
     &   'LAST SPECIES NAME=',A16)
       GOTO 9990
C
 9990  WRITE (OUT,9991) I,LCHAR,IOB160
 9991  FORMAT (' ABOVE ERROR DISCOVERED AT CHARACTER #',I4,
     &   '.  CHAR = ',A1,'.  INPUT LINE BELOW:'/' ',A80)
 9998  NERR=NERR+1
       IR=0
       RETURN
       END
