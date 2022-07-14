C       FILE:  BLDUP.FTN                PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       MAJOR MODIFICATIONS TO PRODUCE FORTRAN CODE BY S.E. HEFFRON
C       LAST UPDATED:  W.P.L. CARTER  6/20/87
C
C
        SUBROUTINE BLDUP
C
C       PRODUCES CODE TO CALCULATE RATES OF FORMATION OF "BUILDUP-ONLY"
C       SPECIES, GIVEN ALREADY CALCULATED RATES OF REACTIONS.
C
C       CALLED BY:  PNHRXN
C
C
C       SPECIFICATIONS OF PREP VARIABLES, PARMS, AND ARRAYS IN COMMON
C
        INCLUDE 'pspecs.inc'
C
C       THIS COMMON IS USED TO COMMUNICATE WITH SUBROUTINES CALLED
C       BY PNHRXN WHICH PRODUCE THE FORTRAN CODE.
C
        INTEGER   DIFNUM,DIFTOT,CODEU
        COMMON /CODE/ NC,NC1,DIFNUM,CODEU
C
C       LOCAL SPECIFICATIONS
C
        CHARACTER*80 FTNSTR
        CHARACTER*16 CNUM
        LOGICAL   FIRST
C
C
C       CALC DERIVITIVES.  LOOP OVER BO SPECIES
C
        DO 5000 IR= NB1,NB2
        FIRST=.TRUE.
C
        WRITE (CODEU,105) NAME(IR)
 105    FORMAT ('C'/'C',8X,A)
        FTNSTR=' '
        WRITE(CNUM,'(I16)')IR-NB0
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR,110) CNUM(:IW)
 110    FORMAT (6X,'DY(',(A),') = ')
        LS=13+IW
        LCOUNT=1
C                                  ! LENGTH OF FTNSTR
C
C       LOOP OVER RXNS TO FIND RXN FORMING IR
C
        DO 5010 IRXN=1,NRXN
        IP1=NPRODS(IRXN)
        IP2=NPRODS(IRXN+1)-1
        IF (IP2.LT.IP1) GOTO 5010
        DO 5015 I=IP1,IP2
        IF (IPRODS(I).EQ.IR) THEN
C                                       ! THIS REACTION FORMS IR
                J=I-1
C MJK changed to (LS.GT.63)                IF (LS.GT.56) THEN
                IF (LS.GT.63) THEN
C                                                   ! CHECK LENGTH OF FTNSTR
                       WRITE (CODEU,'(A)') FTNSTR(:LS)
                       NC=NC+1
C                                                   ! NUMBER OF LINES OF CODE
                       FTNSTR='     &'
C                                                   ! START OF CONTINUATION LINE
                       LS=15
c--check for more than 20 lines--
              LCOUNT=LCOUNT+1
              if(ftnstr(:25).eq.'     &                   ' .and.
     +             LCOUNT.gt.20)then
               WRITE(CNUM,'(I16)')IR-NB0
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               write(ftnstr,286)CNUM(:IW),CNUM(:IW)
 286           format(6x,"DY(",(a),") = DY(",(a),")")
               LS = 17+IW+IW
               LCOUNT = 1
              endif 
c--check for more than 20 lines--
C                                                   ! SET LS TO LINE UP WITH
C                                                     LINE ABOVE
                ENDIF
                IF (FIRST) THEN
                        WRITE(CNUM,'(I16)')IRXN
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),121) CNUM(:IW)
 121                    FORMAT ('R(',(a),')')
                        LS=LS+3+IW
                        LCOUNT=1
                        FIRST=.FALSE.
                ELSE
                        WRITE(CNUM,'(I16)')IRXN
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),120) CNUM(:IW)
 120                    FORMAT ('+R(',(a),')')
                        LS=LS+4+IW
                ENDIF
                IF (I.GT.IP1.AND.IPRODS(J).LT.0) THEN
C                                                       ! MULT RATE BY COEFS
C MJK changed to (LS.GT.62) 5014                   IF (LS.GT.56) THEN
 5014                   IF (LS.GT.62) THEN
C                                                           ! CHECK LEN OF STR
                               WRITE (CODEU,'(A)') FTNSTR(:LS)
                               NC=NC+1
                               FTNSTR='     &'
C                                                           ! START OF CONT. LN.
                               LS=15
c--check for more than 20 lines--
              LCOUNT=LCOUNT+1
              if(ftnstr(:25).eq.'     &                   ' .and.
     +             LCOUNT.gt.20)then
               WRITE(CNUM,'(I16)')IR-NB0
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               write(ftnstr,286)CNUM(:IW),CNUM(:IW)
               LS = 17+IW+IW
               LCOUNT = 1
              endif 
c--check for more than 20 lines--
                        ENDIF
                        WRITE(CNUM,'(I16)')-IPRODS(J)
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),130) CNUM(:IW) 
 130                    FORMAT ('*CO(',(a),')')
                        LS=LS+5+IW
                        J=J-1
                        IF (J.GE.IP1.AND.IPRODS(J).LT.0) GOTO 5014
                ENDIF
        ENDIF
 5015   CONTINUE
C                       ! NOTE THAT IF IR IS PRODUCT TWICE, RATE ADDED TWICE
C
 5010   CONTINUE
C                       ! END LOOP OVER REACTIONS
C
        IF (FIRST) THEN
                WRITE (FTNSTR(LS+1:),'(''0'')')
                LS=LS+1
        ENDIF
C
        WRITE (CODEU,'(A)') FTNSTR(:LS)
        NC=NC+1
 5000   CONTINUE
C                       ! END LOOP OVER BO SPECIES
        NC1=NC
C                       ! NUMBER OF LINES OF CODE IN CURRENT CODE FILE
        RETURN
        END
