C       FILE:  NEWINBUF.FTN                     SAPRC UTILITY ROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       MODIFIED FOR CHARACTER VARIABLES BY M.P. POE AND S.E.HEFFRON
C       LAST UPDATE:  W.P.L. CARTER 6/20/87
C
C
       SUBROUTINE INBUF (PROMPT,L)
C
C      READS IN IOBUF FROM INPUT FILE UNIT.
C      INPUT FILE UNIT IS INITIALLY IN=5 (KEYBOARD),
C      BUT WHEN @(filename) IS INPUT, TAKES INPUT
C      FROM SPECIFIED FILE.  UP TO 4 LEVELS ALLOWED,
C      PROVIDED TASK CAN HAVE AT LEAST THAT MANY
C      ACTIVE FILES.
C
C      CAN READ UP TO 160 CHARACTERS.  READS INTO BYT160.  NO. CHARACTERS
C      INPUT STORED AS IOLEN.  IF LESS THAN 80 CHARS INPUT, CHARACTERS
C      IOLEN+1 THRU 80 ARE BLANK.
C
C      IORET=2 MEANS END-OF-FILE
C      IORET=1 MEANS ERROR
C
       INCLUDE 'newiobuf.inc'
C
       CHARACTER*1 PROMPT(L)
       character*3 chari
       CHARACTER*25 TEXT
       IORET=0
C
C       WRITE PROMPT IF NECESSARY
C
    1   IF (L.GT.0.AND.(IN.EQ.IKEY.OR.SHOW))THEN
          write(text,'(a1)') (prompt(i),i=1,L)
          WRITE (ICRT,'(a1)') '$'//text(1:L)
         ENDIF
C
C       INITIALIZE BUFFER
C
        IOB160=' '
C
C       READ STRING FROM INPUT UNIT
C
        READ (IN,3,END=10) IOLEN,(BYT160(I),I=1,MIN0(IOLEN,MAXIBF))
    3   FORMAT (Q,160A1)
        IOLEN=MIN0(IOLEN,MAXIBF)
C
        IF (SHOW.AND.IN.NE.IKEY) WRITE (ICRT,7) (BYT160(I),I=1,IOLEN)
    7   FORMAT ('+',80A1,:/(' ',131A1))
C
        IF (BYT160(1).EQ.'!') GOTO 1
C                                         !IGNORE STRINGS STARTING WITH '!'
C
C       ANY STRING WITHOUT ! OR @ RETURNED TO CALLING PGM
        IF (BYT160(1).NE.'@') RETURN
C
C------------------------------------------------------------------
C
C       STRING STARTS WITH '@'.  MEANS SPECIAL COMMAND TO INBUF
C
C      '@filename': INPUT FROM NEXT LEVEL.  INCREMENT IN AND OPEN FILE
C
        IF (IOB160(1:2).NE.'@#') THEN
                IF (IN.EQ.MAXUNT) THEN
                        WRITE (ICRT,22) MAXUNT
   22                   FORMAT (' INDIRECT NESTING OF MORE THAN',I2,
     &                   ' NOT ALLOWED')
                        DO 20 I=1,IN
   20                   CLOSE (UNIT=I)
                        IN=IKEY
                        IORET=2
                        RETURN
                ELSE
                        IOB160(1:IOLEN-1)=IOB160(2:IOLEN)
                        IOB160(IOLEN:IOLEN)=' '
                        IF (IN.EQ.IKEY) THEN
                                IN=1
                        ELSE
                                IN=IN+1
                        ENDIF
                        OPEN (UNIT=IN,NAME=IOB160,READONLY,TYPE='OLD'
     &                   ,ERR=90)
                        if(show)then
                         write(chari,'(i3)')IN
                         WRITE(icrt,'(a)') 
     &                   '0---------- INPUT FROM FILE = '//
     &                   iob160(1:iolen)//'   UNIT = '//
     &                   chari//'   ----------'
                        endif
                        GOTO 1
                ENDIF
C
C       '@#END' FORCES END OF ALL INPUT (ALSO TOO MANY INDIRECTS)
C
        ELSEIF (IOB160(1:5).EQ.'@#END') THEN
                IF (IN.LE.MAXUNT) THEN
                        DO 21 I=1,IN
                        IF (SHOW) WRITE (ICRT,28) I
   28                   FORMAT ('0',11X,'UNIT =',I2,' CLOSED.'/' ')
   21                   CLOSE (UNIT=I)
                ENDIF
                IN=IKEY
                IORET=2
                RETURN
C
C       @#SH ... SHOW INDIRECT INPUT LINES
C       @#NS ... DON'T SHOW INDIRECT INPUT
C
        ELSEIF (IOB160(1:4).EQ.'@#SH') THEN
                SHOW=.TRUE.
                GOTO 1
C
        ELSEIF (IOB160(1:4).EQ.'@#NS') THEN
                SHOW=.FALSE.
                GOTO 1
C
C       UNRECOGNIZED COMMAND STARTING WITH '@#'
C
        ELSE
                WRITE (ICRT,94) IOB160(1:32)
   94           FORMAT (' UNRECOGNIZED INBUF COMMAND = ''',A4,'''')
                GOTO 1
C
        ENDIF
C--------------------------------------------------------------
C
C       END OF FILE.  INDICATE DONE WITH ALL INPUTS IF UNIT=5,
C       OTHERWISE, CLOSE FILE AND DEINCREMENT UNIT NO.
C
   10   IF (IN.EQ.IKEY) THEN
                IORET=2
                RETURN
        ELSE
                IF (SHOW) WRITE (ICRT,28) IN
                CLOSE (UNIT=IN)
                IF (IN.EQ.1) THEN
                        IN=IKEY
                        IF (NOKEY) THEN
                                IORET=2
                                RETURN
                        ELSE
                                GOTO 1
                        ENDIF
                ELSE
                        IN=IN-1
                        GOTO 1
                ENDIF
        ENDIF
C-------------------------------------------------------------------
C
C       ERROR IN OPENING FILE
C
   90   WRITE (ICRT,91) IOB160(1:32)
   91   FORMAT ('0ERROR IN OPENING FILE FOR INPUT.'/
     &   ' FILE NAME = ',A)
        IF (IN.GT.1) THEN
                IN=IN-1
                DO 93 I=1,IN
   93           CLOSE (UNIT=IN)
                IORET=2
                IN=IKEY
                RETURN
        ELSE
                IN=IKEY
                WRITE (ICRT,92)
   92           FORMAT ('$TRY AGAIN : ')
                GOTO 1
        ENDIF
C
C--------------------------------------------------------------
C      ENTRY POINT TO INITIALIZE AND GIVE MAX NO INPUT FILES
C
        ENTRY INITIN (L)
        MAXUNT=MIN0(L,4)
        SHOW=.FALSE.
        NOKEY=.FALSE.
        IF (IN.EQ.IKEY) RETURN
        DO 200 I=1,IN
  200   CLOSE (UNIT=I)
        IN=IKEY
        RETURN
C
C--------------------------------------------------------------------------
C      ENTRY POINT TO CLOSE ALL FILES, AND START INPUT FROM FILE
C      WHOSE NAME IS IN IOB160 ARRAY
C
        ENTRY FILEIN
C
        IORET=0
        IF (IN.LE.MAXUNT) THEN
                DO 302 I=1,IN
                IF (SHOW) WRITE (ICRT,28) I
  302           CLOSE (UNIT=I)
        ENDIF
        IN=1
        NOKEY=.TRUE.
        OPEN (UNIT=IN,NAME=IOB160,READONLY,STATUS='OLD',ERR=320)
        IF (SHOW) WRITE (ICRT,324)
  324   FORMAT ('0---------- INPUT FROM FILE GIVEN AT CALL TO FILEIN'
     &   ,' ----------')
        RETURN
C
C       ERROR.  INDICATE SO.
C
  320   IORET=1
        RETURN
C
C------------------------------------------------------------------------
C      ENTRY POINT TO START THE NEXT LEVEL OF INPUT FROM FILE
C      WHOSE NAME IS IN IOB160 ARRAY.  (SAME EFFECT AS @(FILENAME)).
C
        ENTRY FILNXT
C
        IF (IN.EQ.MAXUNT) THEN
                WRITE (ICRT,52) MAXUNT
   52           FORMAT (' INDIRECT NESTING OF MORE THAN',I2,
     &           ' NOT ALLOWED')
                DO 50 I=1,IN
   50           CLOSE (UNIT=I)
                IN=IKEY
                IORET=2
        ELSE
                IF (IN.EQ.IKEY) THEN
                        IN=1
                ELSE
                        IN=IN+1
                ENDIF
                OPEN (UNIT=IN,NAME=IOB160,READONLY,TYPE='OLD'
     &           ,ERR=90)
                IF (SHOW) THEN
                        IOLEN=LENSTR(IOB160,80)
                         write(chari,'(i3)')IN
                         WRITE(icrt,'(a)') 
     &                   '0---------- INPUT FROM FILE = '//
     &                   iob160(1:iolen)//'   UNIT = '//
     &                   chari//'   ----------'
                ENDIF
        ENDIF
        RETURN
        END
C=================================================================
       BLOCK DATA
C
C      INITIALIZE I/O BUFFER ARRAY FOR INBUF.
C
        INCLUDE 'newiobuf.inc'
C
        CHARACTER*4 BUF(44)
        EQUIVALENCE (BYT16(1),BUF(1))
C
       DATA IN/5/,NOKEY/.FALSE./,SHOW/.FALSE./,BUF/44*' '/,MAXUNT/4/
       END
