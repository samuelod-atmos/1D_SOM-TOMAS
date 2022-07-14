C       FILE:  LISTSC.FTN               PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER   6/20/87
C
C	MJM added SNO functionality 11/8/03
C
C
        SUBROUTINE LISTSC
C
C       LISTS REACTANTS AND COEFFICIENTS AFTER THEY HAVE ALREADY
C       BEEN INPUT, AND THE REACTANTS ORDERED AND CLASSIFIED.
C
C       CALLED BY:  PREP MAIN
C
C
C       SPECIFICATIONS OF PREP VARIABLES, PARMS, AND ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C
        WRITE (OUT,600) TITLE
  600   FORMAT ('1',10X,A80/' ')
C
        WRITE (OUT,602) NDY,NSC,NB,NSS,NDUM,NDEL
  602   FORMAT (' '/' LIST OF REACTANTS:  (',I4,' ACTIVE,',I4,
     &   ' CONSTANT,',I4,' BUILDUP,',I4,' STEADY STATE,',I4,
     &   ' DUMMY, AND',I4,' DELETED'/' '/' TYPE           NO. NAME',
     &   '                 MWT    C''S   N''S  S''S     O''S   X''S '
     &   '  CONC(0)'/,' ------',
     &   '------   --- ----------------  -------  -----  ---  ',
     &   '---  -----  -----  ----',
     &   '-------')

C
        WRITE (OUT,603) (I,NAME(I),MWT(I),CNO(I),NNO(I),SNO(I),ONO(I),
     &   XNO(I),CONC0(I),I=NA1,NA2)
  603   FORMAT (' '/' ACTIVE        ',I4,2X,A16,F9.2,F7.2,I4,I5,F7.2,
     &           F7.2,1PG13.3/
     &         (15X,I4,2X,A16,0PF9.2,F7.2,I4,I5,F7.2,F7.2,1PG13.3))
C
        IF (NSC.GT.0) WRITE (OUT,605) (I,NAME(I),MWT(I),CNO(I),
     &   NNO(I),SNO(I),ONO(I),XNO(I),CONC0(I),I=NSC1,NSC2)
  605   FORMAT (' '/' CONSTANT      ',I4,2X,A16,F9.2,F7.2,I4,I5,F7.2,
     &           F7.2,1PG13.3/
     &         (15X,I4,2X,A16,0PF9.2,F7.2,I4,I5,F7.2,F7.2,1PG13.3))
C
        IF (NB.GT.0) WRITE (OUT,607) (I,NAME(I),MWT(I),CNO(I),
     &   NNO(I),SNO(I),ONO(I),XNO(I),CONC0(I),I=NB1,NB2)
  607   FORMAT (' '/' BUILDUP-ONLY  ',I4,2X,A16,F9.2,F7.2,I4,I5,F7.2,
     &             F7.2,1PG13.3/
     &         (15X,I4,2X,A16,0PF9.2,F7.2,I4,I5,F7.2,F7.2,1PG13.3))
C
        IF (NDUM.GT.0) WRITE (OUT,609) (I,NAME(I),MWT(I),CNO(I),
     &   NNO(I),SNO(I),ONO(I),XNO(I),CONC0(I),I=ND1,ND2)
  609   FORMAT (' '/' DUMMY         ',I4,2X,A16,F9.2,F7.2,I4,I5,F7.2,
     &            F7.2,1PG13.3/
     &         (15X,I4,2X,A16,0PF9.2,F7.2,I4,I5,F7.2,F7.2,1PG13.3))
C
        IF (NSS.GT.0) WRITE (OUT,611) (I,NAME(I),I=NS1,NS2)
  611   FORMAT (' '/' STEADY STATE  ',I4,2X,A16/(15X,I4,2X,A16))
C
        IF (NDEL.GT.0) WRITE (OUT,613) (I,NAME(I),I=NS2+1,NS2+NDEL)
  613   FORMAT (' '/' DELETED       ',I4,2X,A16/(15X,I4,2X,A16))
C
       IF (NCOEFV.GT.0) THEN
              WRITE (OUT,616) NCOEFV,(I,COEFNM(I),COEF(I),I=1,NCOEFV)
  616         FORMAT (/1X,I4,' VARIABLE COEFFICIENTS:'/('   ',5(I5,') '
     &         ,A8,1PG10.3)))
       ENDIF
       IF (NCOC2.GE.NCOC1) THEN
              NCOCF=NCOC2-NCOC1+1
              WRITE (OUT,622) NCOCF,(I,COEF(I),I=NCOC1,NCOC2)
  622         FORMAT (/1X,I4,' CONSTANT COEFFICIENTS:'/('   ',5(I5,') '
     &         ,1PG10.3)))
      ENDIF
C
      RETURN
      END
