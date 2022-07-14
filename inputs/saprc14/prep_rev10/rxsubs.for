C       FILE:  RXSUBS.FTN       PREP PROGRAM SUBROUTINES
C
C       CONTAINS SOURCE CODES FOR MISC. SUBS CALLED BY RXLST1 AND
C       RXLST2, WHICH IN TURN ARE CALLED BY RDRXN.
C       MODULES INCLUDE:
C
C       SPCNAM  ... DETERMINES SPECIES NUMBER FROM NAME.  DEFINES TYPE
C                   BASED ON CURRENT VALUE OF "IND".  CREATES NEW SPECIES
C                   IF NEEDED.
C       COFNAM  ... CREATES NEW COEFFICIENT GIVEN COEFFICIENT NAME.
C                   DETERMINES WHETHER IT IS VARIABLE OR CONSTANT,
C                   DEPENDING ON WHETHER ITS NAME CAN BE DECODED AS A
C                   NUMBER.
C       CONTLN  ... READS IN CONTINUATION LINE FOR REACTION
C
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  6/20/87
C
C -----------------------------------------------------------------------------
        SUBROUTINE SPCNAM(ISP,SPNAME,RX)
C
C       DETERMINE SPECIES NO. FOR NAME IN SPNAME.  RETURNS SPECIES NO. IN
C       ISP.  ADDS SPECIES TO NAME ARRAY IF NEW, AND DO OTHER UPDATING.
C
C       CALLED FROM:  RDRXN, RXLST2 (RX=.FALSE.)
C                     RXLST1 (RX=.TRUE.)
C
C
C       ARGUMENTS
C
        LOGICAL   RX
        CHARACTER*16 SPNAME
C
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C
       IF (NS.EQ.0) GOTO 1125
       DO 1122 ISP=1,NS
       IF (NAME(ISP).NE.SPNAME) GOTO 1122
D      WRITE (OUT,*) 'SPCNAM: (OLD) ',NAME(ISP),'  SPC NO.=',ISP
       IF (.NOT.RX) THEN
                IF (IND.GT.0) THEN
                        SPTYP(ISP)=IND
D                       WRITE (OUT,*) 'SPTYP = ',SPTYP(ISP)
                        CONC0(ISP)=0.0
                        MWT(ISP)=0.0
                        CNO(ISP)=0.0
                        NNO(ISP)=0
			SNO(ISP)=0
                        ONO(ISP) = 0.0
                        XNO(ISP) = 0.0
                ELSE
                        WRITE (OUT,*) 'PGM ER. IND<0 AT SPCNAM'
                ENDIF
       ENDIF
       RETURN
 1122  CONTINUE
C  -     NAME NOT FOUND - NEW SPECIES ADDED TO LIST
 1125  NS=NS+1
       IF (NS.GT.MAXNS) THEN
                I=MAXNS
                IF (RX) WRITE (OUT,1127) OUTBUF
 1127           FORMAT (' ',A80)
                WRITE (OUT,*) 'TOO MANY SPECIES.  MAX =',I
                STOP 'TOO MANY SPECIES'
                ENDIF
       NAME(NS)=SPNAME
       IF (IND.GT.0) THEN
                SPTYP(NS)=IND
                CONC0(NS)=0.0
                MWT(NS)=0.0
                CNO(NS)=0.0
                NNO(NS)=0
		SNO(NS)=0
                ONO(NS) = 0.0
                XNO(NS) = 0.0
       ELSE
                WRITE (OUT,*) 'PGM ER. IND<0 AT SPCNAM'
       ENDIF
       ISP=NS
       RETURN
       END

C -----------------------------------------------------------------------------
        SUBROUTINE SPCNAM2(ISP,SPNAME,RX)
C
C       DETERMINE SPECIES NO. FOR NAME IN SPNAME.  RETURNS SPECIES NO. IN
C       ISP.  
C
C
C       ARGUMENTS
C
        LOGICAL   RX
        CHARACTER*16 SPNAME
C
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C
       ISP=-1
       IF (NS.EQ.0) RETURN
       DO 1122 ISP2=1,NS
       IF (NAME(ISP2).EQ.SPNAME)THEN
        ISP=ISP2
        RETURN
       ENDIF
 1122  CONTINUE
       RETURN
       END
C
C -----------------------------------------------------------------------------
        SUBROUTINE COFNAM(ISP,REAC)
C
C       PROCESS VAR COEF INPUT.  EXPECTS NAME TO BE IN IOBF16.  SETS IS TO BE
C       COEFFICIENT NO.  REAC INDICATES REACTANT OR PRODUCT SIDE OF RXN LIST.
C       RETURNS COEFFICIENT NO. IN ISP.
C
C       CALLED BY: RXLST1:
C
C       SUBROUTINES CALLED (ALL "NEWSUBS" UTILITY ROUTINES)
C
C       ALIN16  ... SPLITS UP INPUT STRING INTO GROUPS OF 16 CHARS
C       MOVLFT  ... MOVES CHARACTERS IN 16-CHAR STRING TO LEFT
C       MOVRIT  ... MOVES CHARACTERS IN STRING TO RIGHT
C
C       SPECIFICATION OF ARGUMENTS
C
        LOGICAL   REAC
        CHARACTER*6 RXNUM
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C
       CALL MOVRIT(16,IOBF16)
       DECODE (16,1145,IOBF16,ERR=1146) X
 1145  FORMAT (F16.0)
C  -     DECODED, SO THIS IS A CONSTANT COEF
       IF (NCOC2.GE.NCOC1) THEN
                DO 1147 IC=NCOC1,NCOC2
                IF (X.EQ.COEF(IC)) THEN
                        ISP=-IC
                        RETURN
                ENDIF
 1147           CONTINUE
       ENDIF
C  -     ADD NEW CONSTANT COEF
       NCOC2=NCOC2+1
       IF (NCOC2.GT.MAXCO) THEN
                IC=MAXCO-MAXCOV
                WRITE (OUT,1127) OUTBUF
 1127           FORMAT (' ',A80)
                WRITE (OUT,*) 'TOO MANY CONSTANT COEFFICIENTS. MAX =',IC
                STOP 'TOO MANY CONSTANT COEFFICIENTS'
       ENDIF
       COEF(NCOC2)=X
D      WRITE (OUT,1171) NCOC2,COEF(NCOC2)
D1171  FORMAT (' IN RXLST1: NCOC2=',I4,'  COEF= ',F16.0)
       ISP=-NCOC2
       RETURN
C
C  -     DIDNT DECODE TO A NUMBER, SO THIS IS A VARIABLE COEF
 1146  CALL MOVLFT(IOBF16)
       IF (NCOEFV.GT.0) THEN
                DO 1143 IC=1,NCOEFV
                IF (IOBF16(1:8).EQ.COEFNM(IC)) THEN
                        ISP=-IC
                        RETURN
                ENDIF
 1143           CONTINUE
       ENDIF
C  -     ADD NEW VARIABLE COEF
 1142  IF (IOBF16(1:4).EQ.'RCON') GOTO 1150
       NCOEFV=NCOEFV+1
       IF (NCOEFV.GT.MAXCOV) THEN
                IC=MAXCOV
                WRITE (OUT,1127) OUTBUF
                WRITE (OUT,*) 'TOO MANY VARIABLE COEFFICIENTS. MAX =',IC
                STOP 'TOO MANY VARIABLE COEFFICIENTS'
       ENDIF
       COEFNM(NCOEFV)=IOBF16(1:8)
       COEF(NCOEFV)=0.0
       ISP = -NCOEFV
       RETURN
C
C  -     'RCONnnnn' MEANS INSTEAD OF USING COEF, USE RATE CONSTANT OF RXN nnnn.
C               ONLY VALID IF REACTANT
 1150   IF (.NOT.REAC) GOTO 1195
C                                       !ERROR
        IOBF16(1:4)='    '
        CALL ALIN16(1,IOBF16)
        CALL MOVLFT(IOBF16)
        RXNUM=IOBF16(1:6)
        DO 1154 KK=1,NRXN
        IF (RXNUM.EQ.RXNLBL(KK)) GOTO 1156
 1154   CONTINUE
C       (ERROR)
        GOTO 1190
 1156   ISP=-(KK+MAXMAX)
        RETURN
C
C  -    ERRORS
C
 1190   WRITE (OUT,1191) RXNUM
 1191   FORMAT (' ERROR:  REFERENCE TO UNDEFINED REACTION NO. =',A6)
        ISP=0
        RETURN
C
 1195   WRITE (OUT,1196)
 1196   FORMAT (' ERROR:  RCONnn NOT VALID IN LIST OF PRODUCTS')
        ISP=0
        RETURN
        END
C
C-----------------------------------------------------------------------------
        SUBROUTINE CONTLN(IC,LCHR,NLIN)
C
C       READS IN CONTINUATION LINE.  PUTS FIRST NON-BLANK, NON-TAB CHARACTER
C       IN LCHR.  CONVERTS LATER TABS TO SPACES AND SHIFTS RXN STRING SO
C       THAT LCHR IS FIRST CHAR IN IOB160.
C       IC = NO. CHARACTERS IN LINE.  NLIN IS NUMBER OF LINES FOR THIS
C       REACTION, WHICH IS INCREMENTED BY ONE BY THIS SUBROUTINE.
C
C
C       SPECIFICATIONS FOR ARGUMENTS
C
        CHARACTER*1 LCHR
C
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C       LOCAL VARIABLES
C
        CHARACTER*160 TMPRXN
        LOGICAL   NEWRXN
C
C
       WRITE (OUT,1127) OUTBUF
 1127  FORMAT (' ',A80)
       IPOSN=0
       IC=0
       CALL INBUF('MORE?  :',8)
       IF (IORET.EQ.2) GOTO 8005
C                                       !ERROR.
C  -     IGNORE LEADING BLANKS AND TAB CHAR'S.
 8001  IC=IC+1
       IPOSN=IPOSN+1
       IF (IC.GT.IOLEN) GOTO 8006
C                                       !ERROR.
       LCHR=BYT160(IC)
       ILCHAR=ICHAR(LCHR)
       IF (LCHR.EQ.' ') GOTO 8001
       IF (ILCHAR.EQ.9) THEN
C                                       !TAB
                IPOSN=IPOSN+7
                GOTO 8001
       ENDIF
C
       NLIN=NLIN+1
C                                       ! NO. OF LINES READ FOR THIS RXN
C
C  -    CONVERT TABS TO SPACES
       J=0
       DO 1200 I=IC,IOLEN
       ILCHAR=ICHAR(BYT160(I))
       IF (ILCHAR.NE.9) THEN
                J=J+1
                TMPRXN(J:J)=BYT160(I)
       ELSE
                K=J+IPOSN-1
                K=8-MOD(K,8)
                TMPRXN(J+1:J+K)='        '
C                                               !CONVERT TABS TO SPACES
                J=J+K
       ENDIF
 1200  CONTINUE
       IOB160=TMPRXN
       IOLEN=J
       IC=1
       OUTBUF=' '
       WRITE (OUTBUF,'(a)') '       '//iob160(1:iolen)
       LOUTBF=73
       IF (YESPRP) THEN
               NEWRXN=.FALSE.
C                                       !TO INDICATE NOT FIRST LINE OF RXN
               WRITE (UNIT3) NEWRXN,J,IOB160(1:IOLEN)
D              WRITE (ICRT,*) 'UNIT3: ',NEWRXN,J,IOB160(1:IOLEN)
       ENDIF
       RETURN
C
C  --   ERRORS
C
 8005  WRITE (OUT,8015)
 8015  FORMAT (' ERROR:  END-OF-FILE WHEN CONTINUATION LINE OF RXN',
     &   ' LIST EXPECTED')
       IC=0
       RETURN
C
 8006  WRITE (OUT,8016)
 8016  FORMAT (' ERROR:  BLANK LINE ENCOUNTERED WHEN CONTINUATION',
     &   ' LINE OF RXN LIST EXPECTED')
       IC=0
       RETURN
       END
C
