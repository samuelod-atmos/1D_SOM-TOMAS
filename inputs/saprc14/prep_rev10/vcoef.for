C       FILE:  VCOEF.FTN                        PREP SUBROUTINE
C
C       FIRST WRITTEN BY W.P.L. CARTER
C       MODIFIED BY SUSAN E. HEFFRON
C       RE-MODIFIED BY W.P.L. CARTER
C       LAST UPDATE:  W.P.L. CARTER   6/26/87
C
C
        SUBROUTINE VCOEF(LEN,ISUB)
C
C       PRODUCES CODE TO INSERT INTO SUBROUTINES, BASED ON EXPLICIT
C       "INSERTED" INSTRUCTIONS.  THIS PROCESSES "OLD" FORMAT CODE, WHICH
C       IS NOT AS EASY TO UNDERSTAND AS THE "NEW" VERSION, BUT WHICH
C       IS COMPATABLE BOTH WITH VERSIONS OF PREP WHICH PRODUCE MODEL-
C       DEPENDENT SUBROUTINES IN FORTRAN, AND THOSE WHICH PRODUCE THEM
C       IN MACHINE CODE.  IT HAS THE SAME FUNCTION AS "INSCOD", BUT IS
C       INCLUDED FOR COMPATIBLITY WITH PREP INPUT FILES USED FOR RESEARCH,
C       WHICH GENERALLY USE THE MACHINE CODE VERSIONS.  NOTE THE MACHINE
C       CODE VERSIONS OF PREP AND INT ARE AVAILABLE ONLY FOR DEC PDP-11
C       AND VAX COMPUTERS.
C
C      TAKES INPUT FROM IOB160.  ONLY PROCESSES ONE LINE.  LEN IS LENGTH OF
C      STRING IN IOB160.
C
C     INPUT FORMAT:  LIST OF MACHINE INSTRUCTIONS SEPARATED BY SPACE
C     OR IGNORED CHARACTER.
C     INSTRUCTION FORMAT:  INSTRUCTION CHARACTER | 1ST OPERAND CHAR |
C     2ND OPERAND CHAR(S)
C     INSTRUCTION (OR FORMAT) CHAR.....MEANING
C          *        SET REGISTERS.  IGNORED IN VAX VERSION
C          C        CLEAR REG  (NO 2ND OPERAND)
C          N        NEGATE REG   (NO 2ND OPERAND)
C          L        LOAD REG
C          A        ADD TO REG
C          S        SUBTRACT FROM REG
C          M        MULTIPLY
C          D        DIVIDE
C          E        STORE  CONTENTS OF REG
C          &   (FMT CHAR)   LIST CONTINUES ON NEXT CARD
C        BLANK  (FMT CHAR)  END OF THIS INSTRUCTION LIST
C
C     1ST OPERAND:  LEGAL CHARS=0,1,2...REFER TO FP REGS 9,10,11.
C     2ND OPERAND:  0,1,2...FP REGS (FOR R-R INST).   FOR R-M INST:
C        TYPE           BASE ASSUMED
C     C'specie name'        0          ANY SPECIES
C     K'reaction no.'       1          RATE CONSTANT (ADJ. FOR CONST.)
C     R'reaction no.'       2          REACTION RATE **
C     #'coef name'          3          VARIABLE COEF.
C     D'specie name'        4          ACTIVE OR B.O. DERIVATIVE  **
C     K(RXN LABEL), #(COEF NAME)...REFER TO SPECIES CONC (ACTIVE ONLY),
C     RATE CONSTANT (INCL CONST REACTANTS), AND PRESENT VAR COEF VALUE.
C     BLANK IN COL.1 MEANS END OF INPUT FOR VCOEF.
C
C     ** NOTE:  Rate or Derivitive reference disabled in this version.
C               Code to input it is retained, but program aborts at time
C               index no. is determined.
C
C     ISUB=1 for inserted Initial code (CONSTR routine)
C     Legal operand types: C (active or constant species only), K, #
C
C     ISUB=2 for inserted DIFFUN code
C     Legal operand types: C (active or constant species only), #, K
C
C     ISUB=3 for inserted Print/Save code (SXTRA routine)
C     Legal operand types: C, #, K
C
C
C       CALLED BY:  PNHRXN
C
C       SUBROUTINES CALLED:
C
C               "NEWSUBS" (SAPRC) LIBRARY UTILITY ROUTINES
C
C       ALIN16  ... SPLITS UP INPUT STRING INTO GROUPS OF 16 CHARS
C       MOVLFT  ... MOVES CHARACTERS IN 16-CHAR STRING TO LEFT
C       MOVRIT  ... MOVES CHARACTERS IN STRING TO RIGHT
C
C               PREP SUBROUTINES
C
C       ASMD    ... CONTAINS ENTRY POINTS FOR PRODUCING SPECIFIC MACHINE
C                   LANGUAGE INSTRUCTIONS.  ENTRY POINTS USED BY THIS
C                   ROUTINE ARE:  ADD, ADDR , DIV, DIVR, LOAD ,LOADR, MULT,
C                   MULTR, NEGR, STORE, SUBT, SUBTR
C
C
C       SPECIFICATIONS FOR PREP PARMS, VARIABLES, AND ARRAYS IN COMMON.
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL SPECIFICATIONS
C
        CHARACTER*1 IC1,INC(8),INUMC(9)
        INTEGER   BASENO(6)
        DATA INC/'C','N','L','A','S','M','D','E'/,
     1   INUMC/'0','1','2','C','#','K','R','D','A'/
C
C                   CBASE COBASE KBASE RBASE DYBASE ABASE
        DATA BASENO / 0,    3,     1,    2,    4,    7 /
C
C
C      IGNORE  LEADING BLANKS AND TAB CHARACTERS
C
       IC=0
  100  IC=IC+1
       IF (IC.GT.LEN) RETURN
       IC1=BYT160(IC)
       ICH1=ICHAR(IC1)
       IF (IC1.EQ.' ') GOTO 100
       IF (ICH1.EQ.9) GOTO 100
C                                     !TAB
       GOTO 201
C
C      LOOK AT 1ST CHAR OF NEW INSTRUCTION
C
  200 IC=IC+1
      IF (IC.GT.LEN) RETURN
      IC1=BYT160(IC)
      ICH1=ICHAR(IC1)
      IF (ICH1.EQ.9) GOTO 200
C                                     !TAB
      IF (IC1.EQ.' ') GOTO 200
C
C     DETERMINE INSTRUCTION CHAR
C
  201 IF (IC1.EQ.'*') GOTO 400
      DO 203 INST=1,8
      IF (IC1.EQ.INC(INST)) GOTO 204
  203 CONTINUE
        WRITE (OUT,*) 'ILLEGAL INSTRUCTION CHAR = ''',IC1,''''
        GOTO 800
C
C       IGNORE CHARACTERS AFTER '*' UNTIL BLANK OR ',' ENCOUNTERED.
C
  400   IC=IC+1
        IF (IC.GT.LEN) RETURN
        IC1=BYT160(IC)
        IF (IC1.EQ.' ') GOTO 200
        IF (IC1.EQ.',') GOTO 200
        GOTO 400
C
C     DETERMINE 1ST OPERAND
C
  204 IC=IC+1
      IC1=BYT160(IC)
      DO 205 INUM=1,3
      IF (INUMC(INUM).EQ.IC1) GOTO 206
  205 CONTINUE
        WRITE (OUT,*) 'ILLEGAL FIRST OPERAND = ''',IC1,''''
        GOTO 800
  206 IC=IC+1
      INUM=INUM-1
      IF (INST-2) 210,211,220
C
C      ONE REGISTER INSTRUCTION
C
  210 CALL SUBTR(INUM,INUM)
      GO TO 200
  211 CALL NEGR (INUM)
      GO TO 200
C
C     DETERMINE 2ND OPERAND TYPE
C
  220 IC1=BYT160(IC)
      INST=INST-2
      DO 221 INUM2=1,8
      IF (INUMC(INUM2).EQ.IC1) GO TO 222
  221 CONTINUE
        WRITE (OUT,*) 'ILLEGAL 2ND OPERAND CHAR = ''',IC1,''''
        GO TO 800
  222 IC=IC+1
C
C       DETERMINE IF OPERAND TYPE IS APPROPRIATE FOR THE TYPE
C       OF SUBROUTINE BEING PRODUCED.  NO RESTRICTION FOR TYPE 3
C       TYPE 1 AND 2 DOESN'T SUPPORT R OR D.
C
        IF (ISUB.EQ.1) THEN
                IF (INUM2.GE.7) GOTO 850
        ELSEIF (ISUB.EQ.2) THEN
                IF (INUM2.GE.7) GOTO 860
        ENDIF
C
C
C       INUM2 OF 3 OR LESS IS REGISTER-REGISTER INST.
C
        IF (INUM2.GT.3) GOTO 230
C
        INUM2=INUM2-1
        IF (INST.EQ.1) THEN
                CALL LOADR (INUM,INUM2)
        ELSEIF (INST.EQ.2) THEN
                CALL ADDR (INUM,INUM2)
        ELSEIF (INST.EQ.3) THEN
                CALL SUBTR (INUM,INUM2)
        ELSEIF (INST.EQ.4) THEN
                CALL MULTR (INUM,INUM2)
        ELSEIF (INST.EQ.5) THEN
                CALL DIVR (INUM,INUM2)
        ELSEIF (INST.EQ.6) THEN
                CALL LOADR (INUM2,INUM)
        ELSE
                WRITE (OUT,*) 'PGM ER BY 230.  INST=',INST
                GOTO 890
        ENDIF
        GOTO 200
C
C     REGISTER-MEMORY INST...DETERMINE ADDRESS
C
  230 INUM2=INUM2-3
      IC1=BYT160(IC)
        IF (IC1.NE.'''') THEN
                WRITE (OUT,*) 'FORMAT ERR IN DECODING LABEL.  EXPECTED'
     &           ,' QUOTE, GOT ''',IC1,''' (BY 230)'
                GO TO 800
        ENDIF
      IOBF16=' '
      IC=IC+1
      DO 231 I=1,16
      IC1=BYT160(IC)
      IC=IC+1
      IF (IC1.EQ.'''') GO TO 232
  231 BYT16(I)=IC1
      IC1=BYT160(IC)
        IF (IC1.NE.'''') THEN
                WRITE (OUT,*) 'FORMAT ERR IN DECODING LABEL.  EXPECTED'
     &           ,' QUOTE, GOT ''',IC1,''' (BY 231)'
                GO TO 800
        ENDIF
      IC=IC+1
C
C      FIND ARRAY INDEX NO.
C
C        INUM2:     1   2   3   4   5    6
C        TYPE:      C   #   K   R   D    A
C
  232   IF (INUM2.EQ.1) THEN
C
C               CONCENTRATION.  FIND SPECIES NO.
C
                DO 241 I=1,NS
                IF (NAME(I).EQ.IOBF16) GOTO 242
  241           CONTINUE
                GO TO 810
  242           IF (I.GE.NA1.AND.I.LE.NA2) THEN
C   -                   ACTIVE SPECIES.  USE 'A', OR INUM2=6
                        INUM2=6
                        I=I-N0
                ELSEIF (I.GT.NSC) THEN
C   -                   ONLY ACTIVE OR CONSTANT SPECIES ALLOWED
                        GOTO 870
                ENDIF
C   -                   NO CHANGE IN I OR INUM2 FOR CONSTANT SPECIES
C
        ELSEIF (INUM2.EQ.2) THEN
C
C               COEFFICIENT.  FIND COEF NO, CREATE NEW ONE IF NECC.
C
                CALL MOVRIT(16,IOBF16)
                DECODE (16,248,IOBF16,ERR=249) X
  248           FORMAT (F16.0)
C  -            DECODED, SO THIS IS A CONSTANT COEF
                IF (NCOC2.GE.NCOC1) THEN
                        DO 321 I=NCOC1,NCOC2
                        IF (X.EQ.COEF(I)) GOTO 250
  321                   CONTINUE
                ENDIF
C  -            ADD NEW CONSTANT COEF
                NCOC2=NCOC2+1
                IF (NCOC2.GT.MAXCO) THEN
                        J=MAXCO-MAXCOV
                        WRITE (OUT,*) 'TOO MANY CONSTANT COEFFICIENT',
     &                   'S.  MAX =',J
                        STOP 'TOO MANY CONSTANT COEFFICIENTS'
                ENDIF
                COEF(NCOC2)=X
                I=NCOC2
                GOTO 250
C
C  -            DIDN'T DECODE TO A NUMBER, SO THIS IS A VARIABLE COEF
  249           CALL MOVLFT (IOBF16)
                IF (NCOEFV.GT.0) THEN
                        DO 323 I=1,NCOEFV
                        IF (IOBF16(1:8).EQ.COEFNM(I)) GO TO 250
  323                   CONTINUE
                ENDIF
C  -            ADD NEW VARIABLE COEF
                NCOEFV=NCOEFV+1
                IF (NCOEFV.GT.MAXCOV) THEN
                        J=MAXCOV
                        WRITE (OUT,*) 'TOO MANY VARIABLE COEFFICIENT',
     &                   'S.  MAX =',J
                        STOP 'TOO MANY VARIABLE COEFFICIENTS'
                ENDIF
                COEFNM(NCOEFV)=IOBF16(1:8)
                COEF(NCOEFV)=0.0
                I = NCOEFV
C
        ELSEIF (INUM2.EQ.3) THEN
C
C               RATE CONSTANT.  FIND REACTION LABEL.
C
                CALL ALIN16(1,IOBF16)
                CALL MOVLFT(IOBF16(1:8))
                DO 245 I=1,NRXN
                IF (IOBF16(1:4).EQ.RXNLBL(I)) GOTO 246
  245           CONTINUE
                GOTO 830
  246           IF (RXTYP(I).EQ.2.OR.RXTYP(I).EQ.3) GOTO 840
        ELSEIF (INUM2.EQ.4) THEN
                WRITE (OUT,*) 'REFERENCE TO REACTION RATE NOT SUPPORTED'
                GOTO 890
        ELSEIF (INUM2.EQ.5) THEN
                WRITE (OUT,*) 'REFERENCE TO DERIVITIVE NOT SUPPORTED'
                GOTO 890
        ELSE
                WRITE (OUT,*) 'VCOEF PGM ER BY 250.  INUM2=',INUM2
                GOTO 890
        ENDIF
C
C       MEMORY REFERENCE DEFINED.  PRODUCE INSTRUCTION
C
  250   INUM2=BASENO(INUM2)
        IF (I.LE.0) THEN
                WRITE (OUT,*) 'PGM ERR.  BAD MEMORY REF INDEX = ',I
                GOTO 800
        ENDIF
        IF (INST.EQ.1) THEN
                CALL LOAD (INUM,I,INUM2)
        ELSEIF (INST.EQ.2) THEN
                CALL ADD (INUM,I,INUM2)
        ELSEIF (INST.EQ.3) THEN
                CALL SUBT (INUM,I,INUM2)
        ELSEIF (INST.EQ.4) THEN
                CALL MULT (INUM,I,INUM2)
        ELSEIF (INST.EQ.5) THEN
                CALL DIV (INUM,I,INUM2)
        ELSEIF (INST.EQ.6) THEN
                CALL STORE (INUM,I,INUM2)
        ELSE
                WRITE (OUT,*) 'PGM ER AT 250.  INST=',INST
                GOTO 890
        ENDIF
        GOTO 200
C
C      ERRORS
C
  800  WRITE (OUT,801) IC,IC1
  801  FORMAT ('$ INPUT ERROR AT CHARACTER #',I4,' (CHAR = ''',
     &  A1,''')')
       GOTO 890
  810  WRITE (OUT,811) INUMC(INUM2),IOBF16
  811  FORMAT ('$ INPUT ERROR: ',A1,'''',A16,'''',' NOT FOUND')
       GOTO 890
  830  WRITE (OUT,831) IOBF16(1:4)
  831  FORMAT ('$ INPUT ERROR.  REACTION NO. ',A4,' NOT FOUND')
       GOTO 890
  840  WRITE (OUT,841) IOBF16(1:4)
  841  FORMAT ('$ INPUT ERROR.  REACTION NO. ',A4,' HAS NO K')
       GOTO 890
  850  WRITE (OUT,851) INUMC(INUM2)
  851  FORMAT ('$ INPUT ERROR.  ILLEGAL OPERAND FOR INSERTED ',
     &   'CONSTR CODE. OPERAND TYPE=',A1)
       GOTO 890
  860  WRITE (OUT,861) INUMC(INUM2)
  861  FORMAT ('$ INPUT ERROR.  ILLEGAL OPERAND FOR INSERTED ',
     &   'DIFFUN CODE. OPERAND TYPE=',A1)
       GOTO 890
  870  WRITE (OUT,871) IOBF16
  871  FORMAT ('$ INPUT ERROR.  ILLEGAL SPECIES BYPE.  ',
     &   'CODE MUST BE ACTIVE OR CONSTANT.  ILLEGAL SPECIES =',A16)
  890  WRITE (OUT,891)
  891  FORMAT ('+.   REMAINING INSTRUCTIONS IGNORED.')
       NERR=NERR+1
       WRITE (ICRT,892) IOB160(1:LEN)
  892  FORMAT (' VCOEF ERROR.  INPUT LINE GIVEN BELOW:'
     &  /' ',A)
C
       RETURN
       END
