C       FILE:  INSCOD.FTN                       PREP SUBROUTINE
C
C       FIRST WRITTEN BY SUSAN E. HEFFRON
C       MODIFIED BY W.P.L. CARTER
C       LAST UPDATE:  W.P.L. CARTER   6/26/87
C
C
        SUBROUTINE INSCOD(LEN,ISUB)
C
C       PRODUCES CODE TO INSERT INTO SUBROUTINES, BASED ON EXPLICIT
C       "INSERTED" INSTRUCTIONS.  THIS PROCESSES "NEW" FORMAT CODE, WHICH
C       IS COMPATABLE WITH VERSIONS OF PREP WHICH PRODUCE MODEL-DEPENDENT
C       SUBROUTINES IN FORTRAN.  (IT IS NOT COMPATIBLE WITH VERSIONS OF
C       PREP WHICH PRODUCE MODEL-DEPENDENT SUBROUTINES IN MACHINE CODE)
C
C       TAKES INPUT FROM IOB160.  ONLY PROCESSES ONE LINE.  LEN IS LENGTH
C       OF STRING IN IOB160.
C
C        TYPE           BASE ASSUMED
C     "C species name"       0          ANY SPECIES
C     "K reaction label"     1          RATE CONSTANT (ADJ. FOR CONST.)
C     "R reaction label"     2          REACTION RATE
C     "# coef name"          3          VARIABLE COEF.
C     "D species name"       4          ACTIVE OR B.O. DERIVATIVE  **
C     [species name]         0          SAME AS "C species name"
C     'coef name'            3          SAME AS "# coef name"
C
C     ** Derivitive reference disabled in this version
C
C     ISUB=1 for inserted Initial code (CONSTR routine)
C     Legal operand types: C (active or constant species only), #, K
C
C     ISUB=2 for inserted DIFFUN code
C     Legal operand types: C (active or constant species only), #, K
C
C     ISUB=3 for inserted Print/Save code (SXTRA routine)
C     Legal operand types: C, #, K, R, D
C
C
C       CALLED BY:  PNHRXN
C
C       SUBROUTINES CALLED: (ALL "NEWSUBS" UTILITY ROUTINES)
C
C       MOVLFT  ... MOVES CHARACTERS IN 16-CHAR STRING TO LEFT
C       MOVRIT  ... MOVES CHARACTERS IN STRING TO RIGHT
C
C
C       SPECIFICATIONS FOR PREP PARMS, VARIABLES, AND ARRAYS IN COMMON.
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL SPECIFICATIONS
C
        CHARACTER*1 IC1,OPCHR(5),ENDCHR,CMARK
        CHARACTER*2 ARRNAM(6)
        DATA OPCHR/'C','#','K','R','D'/
        DATA ARRNAM /' C','CO','RK',' R','DY',' A'/
C
C
C      IGNORE  LEADING BLANKS AND TAB CHARACTERS
C
       IC=0
       OUTBUF=' '
       CMARK=' '
       LOUTBF=0
  100  IC=IC+1
       IF (IC.GT.LEN) RETURN
       IC1=BYT160(IC)
       ICH1=ICHAR(IC1)
       IF (IC1.EQ.' ') GOTO 100
       IF (ICH1.EQ.9) GOTO 100
C                                     !TAB
       GOTO 210
C
C      LOOK AT NEXT CHAR OF INSTRUCTION
C
  200 IC=IC+1
C
      IF(LOUTBF.GT.56)THEN
       WRITE(UNIT3,901)CMARK//OUTBUF(:LOUTBF)
       CMARK='&'
       LOUTBF=0
      endif
C      
      IF (IC.GT.LEN) GOTO 900
C                                     ! OUTPUT OUTBUF
      IC1=BYT160(IC)
  210 IF (IC1.EQ.'''') THEN
              ENDCHR=''''
              IOP=2
              GOTO 500
C                                     ! PROCESS COEF
      ELSEIF (IC1.EQ.'"') THEN
              GOTO 300
C                                     ! GET TYPE AND OPERAND
      ELSEIF (IC1.EQ.'[') THEN
              IOP=1
              ENDCHR=']'
C                                     ! PROCESS SPECIES NAME
              GOTO 400
      ELSE
C             COPY CHARACTER INTO OUTPUT BUFFER
              LOUTBF=LOUTBF+1
              OUTBUF(LOUTBF:LOUTBF)=IC1
              GOTO 200
      ENDIF
C
C     GET OPERAND OF THE FORM:
C        "C species"               IOP= 1
C        "# coefficient"                2
C        "K rxnlabel"                   3
C        "R rxnlabel"                   4
C        "D species"                    5
C
  300 IC=IC+1
      IF (IC.GT.LEN) GOTO 9901
C                                    ! ERROR
      IC1=BYT160(IC)
      DO 305 IOP=1,5
      IF (OPCHR(IOP).EQ.IC1) GOTO 310
  305 CONTINUE
      GOTO 9902
C                                    ! ERROR: UNRECOGNIZED OPERAND CHAR
C
C       DETERMINE IF ARRAY TYPE OK FOR SUBROUTINE USED.  NO RESTRICTION
C       FOR ISUB=3.  'D' AND 'R' NOT LEGAL FOR TYPES 1,2.
C
  310   IF (ISUB.EQ.1) THEN
                IF (IOP.GE.4) GOTO 9903
        ELSEIF (ISUB.EQ.2) THEN
                IF (IOP.GE.4) GOTO 9904
        ENDIF
C
      IC=IC+1
      IF (IC.GT.LEN) GOTO 9901
C                                        ! ERROR
      IF (BYT160(IC).NE.' ') GOTO 9905
C                                        ! ERROR: BLANK EXPECTED AFTER OP CHAR
      IF (IOP.EQ.1 .OR. IOP.EQ.5) THEN
C                                       ! C OR D
              ENDCHR='"'
              GOTO 400
C                                     ! PROCESS SPECIES NAME
C
      ELSEIF (IOP.EQ.2) THEN
              ENDCHR='"'
              GOTO 500
C                                     ! PROCESS COEFFICIENT
C
      ELSE
C                         ! (IOP.EQ.3 .OR. IOP.EQ.4) FOR K OR R
C             PROCESS RXN LABEL
              IOBF16=' '
              DO 320 J=1,6
              IC=IC+1
              IF (IC.GT.LEN) GOTO 9901
              IC1=BYT160(IC)
              IF (IC1.EQ.'"') GOTO 330
              IOBF16(J:J)=IC1
 320          CONTINUE
              IC=IC+1
              IF (IC.GT.LEN) GOTO 9901
              IC1=BYT160(IC)
              IF (IC1.NE.'"') GOTO 9906
C                                           ! ERROR: RXN LABEL TOO LONG
C
C             DETERMINE WHICH RXN
C
 330          IF (IOBF16.EQ.' ') GOTO 9907
C                                           ! ERROR: NO RXN LABEL AFTER K OR R
              DO 340 IX=1,NRXN
              IF (IOBF16.EQ.RXNLBL(IX)) THEN
                    IF (IOP.EQ.3) THEN
                           IF (RXTYP(IX).EQ.2 .OR. RXTYP(IX).EQ.3)
     &                       GOTO 9908
C                                           ! ERROR: RXN HAS NO K
                    ENDIF
                    WRITE (OUTBUF(LOUTBF+1: ),332) ARRNAM(IOP),IX
 332                FORMAT (A2,'(',I4,')')
                    LOUTBF=LOUTBF+8
                    GOTO 200
              ENDIF
 340          CONTINUE
              GOTO 9909
C                                            ! ERROR: RXN NOT FOUND
       ENDIF
C
C
C      PROCESS SPECIES NAME (OPERAND TYPE C OR D)
C       **NOTE** TYPE "D" DISABLED IN THIS VERSION
C
 400    IF (IOP.NE.1) THEN
                WRITE (OUT,401) ARRNAM(IOP)
 401            FORMAT (' USE OF ARRAY TYPE = '',A1,'' DISABLED.')
                GOTO 9990
        ENDIF
C
       IOBF16=' '
       DO 420 J=1,16
       IC=IC+1
       IF (IC.GT.LEN) GOTO 9901
       IC1=BYT160(IC)
       IF (IC1.EQ.ENDCHR) GOTO 430
       IOBF16(J:J)=IC1
 420   CONTINUE
       IC=IC+1
       IF (IC.GT.LEN) GOTO 9901
       IC1=BYT160(IC)
       IF (IC1.NE.ENDCHR) GOTO 9912
C                                       ! ERROR: SPECIES NAME TOO LONG
C
C  -     WHICH SPECIES NO.?
  430   IF (IOBF16.EQ.' ') GOTO 9915
C                                        ! ERROR: MISSING SPECIES NAME
        DO 431 I=1,NS
        IF (NAME(I).EQ.IOBF16) GOTO 433
  431   CONTINUE
        GOTO 9917
  433   IF (I.LE.NSC) THEN
C               CONSTANT SPECIES
                WRITE (OUTBUF(LOUTBF+1: ),432) ARRNAM(1),I
  432           FORMAT (A2,'(',I4,')')
                LOUTBF=LOUTBF+8
        ELSEIF (I.GE.NA1.AND.I.LE.NA2) THEN
C               ACTIVE SPECIES
                WRITE (OUTBUF(LOUTBF+1: ),432) ARRNAM(6),I-N0
                LOUTBF=LOUTBF+8
        ELSE
                GOTO 9913
        ENDIF
        GOTO 200
C
C      PROCESS COEFFICIENT
C
 500   IOBF16=' '
       DO 520 J=1,8
       IC=IC+1
       IF (IC.GT.LEN) GOTO 9901
       IC1=BYT160(IC)
       IF (IC1.EQ.ENDCHR) GOTO 530
       IOBF16(J:J)=IC1
 520   CONTINUE
       IC=IC+1
       IF (IC.GT.LEN) GOTO 9901
       IC1=BYT160(IC)
       IF (IC1.NE.ENDCHR) GOTO 9914
C                                       ! ERROR: COEF NAME TOO LONG
C
C  -     WHICH COEF. NO.?
  530  IF (IOBF16.EQ.' ') GOTO 9916
C                                       ! ERROR: COEF NAME MISSING
       CALL MOVRIT(16,IOBF16)
       DECODE (16,531,IOBF16,ERR=540) X
  531  FORMAT (F16.0)
C  -     DECODED, SO THIS IS A CONSTANT COEF
       IF (NCOC2.GE.NCOC1) THEN
                DO 532 I=NCOC1,NCOC2
                IF (X.EQ.COEF(I)) GOTO 550
  532           CONTINUE
       ENDIF
C  -     ADD NEW CONSTANT COEF
       NCOC2=NCOC2+1
       IF (NCOC2.GT.MAXCO) THEN
                J=MAXCO-MAXCOV
                WRITE (OUT,*) 'TOO MANY CONSTANT COEFFICIENTS.  MAX =',J
                STOP 'TOO MANY CONSTANT COEFFICIENTS IN inscod'
       ENDIF
       COEF(NCOC2)=X
       I=NCOC2
       GOTO 550
C
C  -     DIDN'T DECODE TO A NUMBER, SO THIS IS A VARIABLE COEF
  540  CALL MOVLFT (IOBF16)
       IF (NCOEFV.GT.0) THEN
                DO 541 I=1,NCOEFV
                IF (IOBF16(1:8).EQ.COEFNM(I)) GOTO 550
  541           CONTINUE
       ENDIF
C  -     ADD NEW VARIABLE COEF
       NCOEFV=NCOEFV+1
       IF (NCOEFV.GT.MAXCOV) THEN
                J=MAXCOV
                WRITE (OUT,*) 'TOO MANY VARIABLE COEFFICIENTS.  MAX =',J
                STOP 'TOO MANY VARIABLE COEFFICIENTS'
       ENDIF
       COEFNM(NCOEFV)=IOBF16(1:8)
       COEF(NCOEFV)=0.0
       I = NCOEFV
C
 550   WRITE (OUTBUF(LOUTBF+1: ),552) ARRNAM(IOP),I
 552   FORMAT (A2,'(',I4,')')
       LOUTBF=LOUTBF+8
       GOTO 200
C
C      OUTPUT INSTRUCTION TO FILE OF FORTRAN CODE
C
  900  WRITE (UNIT3,901) CMARK//OUTBUF(1:LOUTBF)
  901  FORMAT (5X,A)
       CMARK=' '
       RETURN
C
C      ERRORS
C
 9901  WRITE (OUT,8801)
 8801  FORMAT (' INPUT ERROR.  EOL WHEN DECODING OPERAND FOR',
     &   ' INSERTED CODE.')
       GOTO 9990
 9902  WRITE (OUT,8802) IC1
 8802  FORMAT (' INPUT ERROR.  UNRECOGNIZED OPERAND CHARACTER: ',A1)
       GOTO 9990
 9903  WRITE (OUT,8803) IC1
 8803  FORMAT (' INPUT ERROR.  ILLEGAL OPERAND FOR INSERTED ',
     &   'CONSTR CODE. OPERAND TYPE=',A1)
       GOTO 9990
 9904  WRITE (OUT,8804) IC1
 8804  FORMAT (' INPUT ERROR.  ILLEGAL OPERAND FOR INSERTED ',
     &   'DIFFUN CODE. OPERAND TYPE=',A1)
       GOTO 9990
 9905  WRITE (OUT,8805) IC1
 8805  FORMAT (' INPUT ERROR.  BLANK EXPECTED AFTER OPERAND CHAR: ',A1)
       GOTO 9990
 9906  WRITE (OUT,8806) IOBF16(1:4)
 8806  FORMAT (' INPUT ERROR.  RXN LABEL TOO LONG OR CLOSING " MISSING'
     &   ,'.  FIRST PART OF LABEL=',A4)
       GOTO 9990
 9907  WRITE (OUT,8807) OPCHR(IOP)
 8807  FORMAT (' INPUT ERROR.  NO RXN LABEL AFTER OPERAND CHAR ',A1)
       GOTO 9990
 9908  WRITE (OUT,8808) IOBF16(1:4)
 8808  FORMAT (' INPUT ERROR.  REACTION NO. ',A4,' HAS NO K')
       GOTO 9990
 9909  WRITE (OUT,8809) IOBF16(1:4)
 8809  FORMAT (' INPUT ERROR.  REACTION NO. ',A4,' NOT FOUND')
       GOTO 9990
 9912  WRITE (OUT,8812) ENDCHR,IOBF16
 8812  FORMAT (' INPUT ERROR.  SPECIES NAME TOO LONG OR CLOSING ',A1,
     &   ' MISSING.  FIRST PART OF NAME=',A16)
       GOTO 9990
 9913  WRITE (OUT,8813) IOBF16
 8813  FORMAT (' INPUT ERROR.  SPECIES MUST BE ACTIVE OR CONSTANT.',
     &   '  ILLEGAL SPECIES =',A16)
       GOTO 9990
 9914  WRITE (OUT,8814) ENDCHR,IOBF16(1:8)
 8814  FORMAT (' INPUT ERROR.  COEF NAME TOO LONG OR CLOSING ',A1,
     &   ' MISSING.  FIRST PART OF NAME=',A8)
       GOTO 9990
 9915  IF (ENDCHR.EQ.'"') THEN
             WRITE (OUT,8815) OPCHR(IOP)
 8815        FORMAT (' INPUT ERROR.  MISSING SPECIES NAME AFTER OPERAND'
     &       ,' CHAR ',A1)
       ELSE
             WRITE (OUT,8715)
 8715        FORMAT (' INPUT ERROR.  MISSING SPECIES NAME AFTER ''[''.')
       ENDIF
       GOTO 9990
 9916  IF (ENDCHR.EQ.'"') THEN
             WRITE (OUT,8816) OPCHR(IOP)
 8816        FORMAT (' INPUT ERROR.  MISSING COEF NAME AFTER OPERAND ',
     &        'CHAR ',A1)
       ELSE
             WRITE (OUT,8716)
 8716        FORMAT (' INPUT ERROR.  MISSING COEF NAME AFTER ''.')
       ENDIF
       GOTO 9990
 9917  WRITE (OUT,8817) IOBF16
 8817  FORMAT (' INPUT ERROR.  SPECIES ',A16,' NOT FOUND')
C
 9990  NERR=NERR+1
        WRITE (ICRT,8892) IOB160(1:LEN)
 8892  FORMAT (' ERROR IN INSERTED CODE.  INPUT LINE GIVEN BELOW:'
     &  /' ',A)
C
       RETURN
       END
