C       FILE: NEWSUBS.FTN       UTILITY SUBROUTINES USED BY SAPRC
C                               MODELING PROGRAMS.
C
C       WRITTEN BY W.P.L. CARTER, M.P. POE, AND S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  6/87
C
C==== IN THIS FILE ARE SUBROUTINES:
C       ALIN16  TMARK   FILNAM  MOVRIT  ALIN8   MOVLFT  MOVLF8
C       TNOU    TNOUA
C==== IN THIS FILE ARE FUNCTIONS:
C      LENFNM  NEQ16    LENSTR
C
C ============================================
       SUBROUTINE TNOU (C,L)
C
        PARAMETER ICRT=6
C                               ! FOR PDP ONLY
C
       LOGICAL   C(L)
C
C ***  CALL TNOU
C
       WRITE (ICRT,1) C
 1     FORMAT (' ',80A1)
       RETURN
       END
C
C *** CALL TNOUA
C
       SUBROUTINE TNOUA (C,L)
C
        PARAMETER ICRT=6
C                               ! FOR PDP ONLY
C
       LOGICAL   C(L)
       WRITE (ICRT,2) C
 2     FORMAT ('$',80A1)
       RETURN
C
C *** CALL TNOUB
C
       ENTRY TNOUB (C,L)
       WRITE (ICRT,3) C
 3     FORMAT ('+',80A1)
       RETURN
C
C *** CALL TNOUC
C
       ENTRY TNOUC (C,L)
       WRITE (ICRT,4) C
 4     FORMAT ('+',80A1,$)
       RETURN
       END
C ============================================
       SUBROUTINE ALIN16(L,STR)
C
C      ALIGNS CONTENTS OF STRING INTO MAX. OF L SECTIONS
C      OF 16 CHARACTERS EACH, PADDED TO THE LEFT WITH BLANKS
C      IF L>10 OR <=0, L=10 IS ASSUMED.
C      SECTONS DELINIATED BY ','. '=', OR BLANK.
C
C	MJK modified to have up to 10 16 character segments
C
C      *****************************************************
C      PRIME FORM      DECODE(80,*,BUF) I,X,J
C
C      DEC FORM        CALL ALIN16(3,BUF)
C                      DECODE(80,1000,BUF) I,X,J
C                 1000 FORMAT (I16,G16.0,I16)
C      *****************************************************
       CHARACTER*1 STR(160),OUTSTR(160),A
       LOGICAL   COMMA
       CHARACTER*8 OUT8(20)
       EQUIVALENCE (OUT8(1),OUTSTR(1))
C
C      OUTPUT STRING INITIALIZED TO BLANK
C
       DO 1 I=1,20
 1     OUT8(I)='        '
C
       IF (L.LE.0 .OR.L.GT.10) L=10
       N=16*L+1
       COMMA=.FALSE.
       M=0
C
C      NEW SECTION.  FIND FIRST NO-BLANK CHAR
C
 100   N=N-1
       IF (N.EQ.0) THEN
              IF (COMMA) M=M+1
              GOTO 200
       ENDIF
       A=STR(N)
       IA=ICHAR(A)
       IF (A.EQ.' ') GOTO 100
       IF (IA.EQ.9) GOTO 100
C                                  !TAB
       IF (.NOT.COMMA.AND.(A.EQ.','.OR.A.EQ.'=')) THEN
                COMMA=.TRUE.
                GOTO 100
       ENDIF
       M=M+1
       IF (M.GT.L) GOTO 800
       I2=16*(L-M)
       I=I2+16
       IF (A.NE.',' .AND. A.NE.'=') GOTO 122
       COMMA=.TRUE.
       GOTO 100
C
C      ADD CHARACTERS FOR SECTION.  BLANK, COMMA, OR BEGIN-OF-STRING
C      TERMINATES IT
C
  121   N=N-1
       IF (N.EQ.0) GOTO 200
       A=STR(N)
       IA=ICHAR(A)
       IF (A.EQ.' ' .OR. IA.EQ.9) THEN
C                                             !SPACE OR TAB
             COMMA=.FALSE.
             GOTO 100
       ENDIF
       IF (A.EQ.',' .OR. A.EQ.'=') THEN
             COMMA=.TRUE.
             GOTO 100
       ENDIF
       IF (I.LE.I2)  GOTO 810
 122   OUTSTR(I)=A
       I=I-1
       GOTO 121
C
C      END OF INPUT.  COPY OUTPUT STRING INTO INPUT ONE, AND RETURN
C
 200   IF (M.GT.L) GOTO 800
C                                 !MORE THAN MAX NO. SECTIONS FOUND
       L2=2*L
       M2=2*(L-M)
       M1=M2+1
       IF (M.EQ.0) THEN
C                                !INPUT STRING HAD ALL BLANK SECTIONS
            ENCODE (16*L,245,STR) ('        ',I=1,L2)
       ELSEIF (M.EQ.L) THEN
            ENCODE (16*L,245,STR) (OUT8(I),I=1,L2)
       ELSE
C              !(M.LT.L)
            ENCODE (16*L,245,STR) (OUT8(I),I=M1,L2),('        ',I=1,M2)
       ENDIF
  245  FORMAT (20A8)
       RETURN
C
C      ERROR SECTION
C
 800   ENCODE (16,802,STR) L
 802   FORMAT ('ERROR   NSEC>',I1,'  ')
       RETURN
C
 810   ENCODE (16,812,STR)
 812   FORMAT ('ERROR   SEC>16  ')
       RETURN
       END
C =============================================================
       FUNCTION LENFNM (NAME,N)
       CHARACTER*1 NAME(1)
C --- FIND FIRST BLANK
       DO 100 I=1,N
       IF (NAME(I).EQ.' ') GOTO 110
 100   CONTINUE
       LENFNM=N
       RETURN
 110   NAME(I)=' '
       LENFNM=I-1
       RETURN
       END
C =============================================================
       FUNCTION LENSTR (NAME,N)
       CHARACTER*1 NAME(*)
C ---  FIND FIRST NON-BLANK, STARTING FROM END
       DO 100 I=N,1,-1
         IF (NAME(I).NE.' ') GOTO 110
 100   CONTINUE
       LENSTR=0
       RETURN
 110   LENSTR=I
       RETURN
       END
C ============================================================
       SUBROUTINE MOVRIT (N,BUF)
C
C      MOVES CHARACTERS SO PADDED TO RIGHT
C
       CHARACTER*1 BUF(N)
C
       DO 10 L=N,1,-1
       IF (BUF(L).NE.' ') GOTO 20
   10  CONTINUE
       RETURN
C
   20  IF (L.EQ.N) RETURN
       NMOV=N-L
       DO 30 I=L,1,-1
   30  BUF(I+NMOV)=BUF(I)
C
       DO 40 I=1,NMOV
   40  BUF(I)=' '
       RETURN
       END
C=================================================================
      FUNCTION NEQ16 (C,N,CS)
      CHARACTER*8 C(2),CS(2,N)
      DO 1 I=1,N
      IF (C(1).EQ.CS(1,I).AND.C(2).EQ.CS(2,I)) GO TO 2
    1 CONTINUE
      NEQ16=0
      RETURN
    2 NEQ16=I
      RETURN
      END                                                               NEQ16
C====================================================================
       SUBROUTINE TMARK (BUF)
       CHARACTER BUF*16,BUF2*8
       CALL DATE (BUF)
       CALL TIME (BUF2)
       BUF(10:10)=' '
       BUF(11:15)=BUF2(1:5)
       BUF(16:16)=' '
       RETURN
       END
C=====================================================================
       SUBROUTINE FILNAM (N,OUTSTR,STR1,STR2,STR3)
C
C      FORMS FILE NAME IN OUTSTR FROM STR1+STR2+STR3
C      **NOTE: OUTSTR SHOULD BE INITIALIZED TO BLANKS BEFORE CALLING THIS
C        SUBROUTINE, OR JUNK WILL BE LEFT AT THE END OF THE STRING
C
       CHARACTER*1 OUTSTR(1),STR1(1),STR2(1),STR3(1)
C
       N=0
       M=0
   11  M=M+1
       IF (STR1(M).EQ.' ') GOTO 20
       N=N+1
       OUTSTR(N)=STR1(M)
       GOTO 11
   20  M=0
D      WRITE (5,15) N,(OUTSTR(I),I=1,N)
C                                           !NOTE: UNIT 5 FOR PDP ONLY
D 15   FORMAT(' STR1:  LEN= ',I4,'NAME= ',<N>A1)
   21  M=M+1
       IF (STR2(M).EQ.' ') GOTO 30
       N=N+1
       OUTSTR(N)=STR2(M)
       GOTO 21
   30  M=0
D      WRITE (5,35) N,(OUTSTR(I),I=1,N)
C                                            !NOTE: UNIT 5 FOR PDP ONLY
D 35   FORMAT(' STR2:  LEN= ',I4,'NAME= ',<N>A1)
   31  M=M+1
       IF (STR3(M).EQ.' ') GOTO 40
       N=N+1
       OUTSTR(N)=STR3(M)
       GOTO 31
   40  CONTINUE
D      WRITE (5,50) (OUTSTR(I),I=1,N)
C                                          !NOTE: UNIT 5 FOR PDP ONLY
D 50   FORMAT (' FILNAM : ',<N>A1)
       RETURN
       END
C===========================================================
       SUBROUTINE MOVLFT (BUF)
C
C      MOVES CHARACTERS ORIGINALY PADDED TO LEFT (FROM ALIN16) WITH
C      BLANKS SO THAT BLANKS ARE PADDED TO THE RIGHT INSTEAD.
C
       CHARACTER*1 BUF(16)
C
       DO 10 N=1,16
       IF (BUF(N).NE.' ') GOTO 20
   10  CONTINUE
       RETURN
C
   20  IF (N.EQ.1) RETURN
       N=N-1
       L=16-N
       DO 30 I=1,L
   30  BUF(I)=BUF(I+N)
C
       N=L+1
       DO 40 I=N,16
   40  BUF(I)=' '
       RETURN
C
C-------------------------------------------------------------------------
       ENTRY MOVLF8 (BUF)
C
       DO 110 N=1,8
       IF (BUF(N).NE.' ') GOTO 120
  110  CONTINUE
       RETURN
C
  120  IF (N.EQ.1) RETURN
       N=N-1
       L=8-N
       DO 130 I=1,L
  130  BUF(I)=BUF(I+N)
C
       N=L+1
       DO 140 I=N,8
  140  BUF(I)=' '
       RETURN
       END
C=================================================
C
       SUBROUTINE ALIN8(LEN,STR)
C
C      ALIGNS CONTENTS OF STRING INTO MAX. OF L SECTIONS
C      OF 8 CHARACTERS EACH, PADDED TO THE LEFT WITH BLANKS
C      IF L>10 OR <=0, L=10 IS ASSUMED.
C      SECTONS DELINIATED BY ','. '=', OR BLANK.
C      (TAB [CHR 9] TREATED AS SINGLE BLANK)
C
C      *****************************************************
C      PRIME FORM      DECODE(80,*,BUF) I,X,J
C
C      DEC FORM        CALL ALIN8(3,BUF)
C                      DECODE(80,1000,BUF) I.X.J)
C                 1000 FORMAT (I8,G8.0,I8)
C      *****************************************************
       CHARACTER*1 STR(80),OUTSTR(88),A
       LOGICAL   COMMA
       CHARACTER*8 OUT8(11)
       EQUIVALENCE (OUT8(1),OUTSTR(1))
C
       L=LEN
C
C      OUTPUT STRING INITIALIZED TO BLANK
C
       IF (L.LE.0 .OR.L.GT.10) L=10
       N=8*L+1
C
       DO 1 I=1,L
 1     OUT8(I)='        '
C
       COMMA=.FALSE.
       M=0
C
C      NEW SECTION.  FIND FIRST NO-BLANK CHAR
C
 100   N=N-1
       IF (N.EQ.0) THEN
              IF (COMMA) M=M+1
              GOTO 200
       ENDIF
       A=STR(N)
       IA=ICHAR(A)
       IF (A.EQ.' '.OR.IA.EQ.9) GOTO 100
C                                             !SPACE OR TAB
       IF (.NOT.COMMA.AND.(A.EQ.','.OR.A.EQ.'=')) THEN
                COMMA=.TRUE.
                GOTO 100
       ENDIF
       M=M+1
       IF (M.GT.L) GOTO 800
       I2=8*(L-M)
       I=I2+8
       IF (A.NE.','.AND.A.NE.'=') GOTO 122
       COMMA=.TRUE.
       GOTO 100
C
C      ADD CHARACTERS FOR SECTION.  BLANK, COMMA, OR BEGIN-OF-STRING
C      TERMINATES IT
C
  121   N=N-1
       IF (N.EQ.0) GOTO 200
       A=STR(N)
       IA=ICHAR(A)
       IF (A.EQ.' '.OR.IA.EQ.9) THEN
C                                         !SPACE OR TAB
             COMMA=.FALSE.
             GOTO 100
       ENDIF
       IF (A.EQ.','.OR.A.EQ.'=') THEN
             COMMA=.TRUE.
             GOTO 100
       ENDIF
       IF (I.LE.I2)  GOTO 800
 122   OUTSTR(I)=A
       I=I-1
       GOTO 121
C
C      END OF INPUT.  COPY OUTPUT STRING INTO INPUT ONE, AND RETURN
C
 200   IF (M.GT.L) GOTO 800
C                                !MORE THAN L SECTIONS FOUND
       M2=L-M
       M1=M2+1
       IF (M.EQ.0) THEN
C                                !INPUT STRING HAD ALL BLANK SECTIONS
             ENCODE (8*L,245,STR) ('        ',I=1,L)
       ELSEIF (M.EQ.L) THEN
             ENCODE (8*L,245,STR) (OUT8(I),I=1,L)
       ELSE
C             ! (M.LT.L)
             ENCODE (8*L,245,STR) (OUT8(I),I=M1,L),('        ',I=1,M2)
       ENDIF
  245  FORMAT (10A8)
       RETURN
C
C      ERROR SECTION
C
 800   ENCODE (8,801,STR)
 801   FORMAT ('ERROR   ')
       RETURN
       END
C
C================================================================
