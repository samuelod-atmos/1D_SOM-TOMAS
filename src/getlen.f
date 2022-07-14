
C     **************************************************
C     *  GETLEN                                        *
C     **************************************************

C     Taken from "Problem Solving and Structured Programming in
C     FORTRAN 77", 4th edition, by Koffman and Friedman

C     Returns the length of any string excluding any blank "padding"
C     characters.

      INTEGER FUNCTION GETLEN(STRING)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS------------------------------------------
      CHARACTER *(*) STRING
C-----VARIABLE DECLARATIONS------------------------------------------
      CHARACTER *1 BLANK
      PARAMETER (BLANK = ' ')
      INTEGER NEXT
C-----CODE-----------------------------------------------------------
C     Start with the last character and find the first nonblank
      DO 10 NEXT = LEN(STRING),1,-1
         IF (STRING(NEXT:NEXT).NE.BLANK) THEN
            GETLEN=NEXT
            RETURN
         ENDIF
 10   CONTINUE
C     All characters are blanks
      GETLEN=0
      RETURN
      END
