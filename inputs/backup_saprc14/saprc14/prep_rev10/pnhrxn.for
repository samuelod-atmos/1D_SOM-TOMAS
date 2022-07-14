C       FILE:  PNHRXN.FTN               PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  6/20/87
C
        SUBROUTINE PNHRXN
C
C       CONTROLS PRODUCTION OF FORTRAN CODE FOR MODEL DEPENDENT
C       SUBROUTINES.  OPENS model.FTN FILE, WRITES INITIALIZATION
C       STATEMENTS, CALLS ROUTINES TO WRITE MODEL-SPECIFIC CODE,
C       CLOSES FILES.
C
C       CALLED BY:  PREP MAIN
C
C               "NEWSUBS" (SAPRC) LIBRARY UTILITY ROUTINES
C
C       FILNAM  ... ROUTINE TO CREATE FILE NAME
C
C               PREP SUBROUTINES
C
C       BLDUP   ... PRODUCES CODE TO CALCULATE RATES OF FORMATION
C                   OF "BUILDUP-ONLY" SPECIES, GIVEN ALREADY
C                   CALCULATED RATES OF REACTIONS.
C       DIFFUN  ... PRODUCES CODE TO CALCULATE RATES OF REACTIONS
C                   AND DERIVITIES, RATES OF FORMATION, AND PSEUDO-
C                   UNIMOLECULAR RATES OF REACTION FOR ACTIVE SPECIES
C                   ALSO CALCULATES STEADY STATE SPECIES CONC'S.
C       VCOEF   ... PRODUCES CODE TO INSERT INTO SUBROUTINES, BASED
C                   ON EXPLICIT "INSERTED" INSTRUCTIONS.  THIS PROCESSES
C                   "OLD" FORMAT CODE, WHICH IS COMPATIBLE WITH MACHINE
C                   INSTRUCTION VERSIONS OF PREP FOR PDP-11 AND VAX
C                   COMPUTERS.  THIS "OLD" FORMAT IS NOT AS EASY
C                   TO UNDERSTAND AS THE "NEW" FORMAT, BUT IS RETAINED
C                   BECAUSE MACHINE-INSTRUCTION VERSIONS OF PREP ARE
C                   USED IN MODEL DEVELOPMENT RESEARCH.
C       INSCOD  ... SAME AS VCOEF, EXCEPT PROCESSES INPUT GIVEN IN THE
C                   "NEW" FORTRAN-LIKE FORMAT.
C
C
C      SPECIFICATIONS FOR PREP PARMS, VARIABLES, AND ARRAYS IN COMMON.
C
       INCLUDE 'pspecs.inc'
C
C       THIS COMMON IS USED TO COMMUNICATE WITH SUBROUTINES CALLED
C       BY PNHRXN WHICH PRODUCE THE FORTRAN CODE.
C
        COMMON /CODE/ NC,NC1,DIFNUM,CODEU
C
C       LOCAL SPECIFICATIONS
C
        LOGICAL   U2OPEN
        INTEGER   DIFNUM,DIFTOT,CODEU
        CHARACTER CNAME*16
C
C
        U2OPEN=.FALSE.
C
        IOB160=' '
        CALL FILNAM (NLEN,IOB160,MODUIC,MODFIL,'.ftn ')
        CODEU=UNIT3
        OPEN (UNIT=CODEU,NAME=IOB160,STATUS='UNKNOWN')
c     &        CARRIAGECONTROL='LIST')
C
C
C       PRODUCE 'CONSTR' SUBROUTINE.  THIS IS CALLED AT START OF
C       INTEGRATION ONLY.  IN THIS VERSION, IT INCLUDES ONLY EXPLICITLY
C       "INCLUDED" STATEMENTS, SUCH AS DEFINING A VALUE OF COEFFICIENTS
C       WHICH DEPEND ON OTHER COEF'S, INITIAL CONCENTRATIONS WHICH DEPEND
C       ON OTHER COEFFICIENTS, ETC.  SUCH CODE IS NOT AUTOMATICALLY
C       PRODUCED, AND IF NO SUCH INSTRUCTIONS ARE GIVEN, THEN SUBROUTINE
C       DOES NOTHING.
C
        WRITE (CODEU,700)
 700    FORMAT ('C'/'C'/6X,'SUBROUTINE CONSTR (C,A,RK,CO)'/'C'/
     &   'C',5X,'MODEL DEPENDENT SUBROUTINE FOR INITIALIZATION'/
     &   'C'/6X,'REAL*4 C(*),A(*),RK(*),CO(*)'/'C')
        NC=0
        NC1=0
        DIFNUM=-1
C
        IF (VCOCAL(1).OR.VCOCAL(4)) THEN
C               ('.#C1' OR 'INS INIT)
                WRITE (OUT,1001)
 1001           FORMAT ('0INITIAL CONSTR CODE INSERTED:')
                WRITE (CODEU,1002)
 1002           FORMAT ('C',5X,'INSERTED INSTRUCTIONS'/'C')
                IF (VCOCAL(1)) THEN
C                       ('.#C1')
                       ASSIGN 1004 TO IGO
                       IND=1
                       GOTO 2000
                ENDIF
 1004           IF (VCOCAL(4)) THEN
C                       ('.INS INIT')
                       ASSIGN 1005 TO IGO
                       IND=4
                       GOTO 2000
                ENDIF
        ENDIF
 1005   WRITE (CODEU,'(6X,''RETURN''/6X,''END'')')
C
C
C       PRODUCE 'BLDUP' SUBROUTINE TO CALCULATE RATES OF FORMATION
C       OF UNREACTING SPECIES.
C
        WRITE (CODEU,701)
 701    FORMAT ('C'/'C'/6X,'SUBROUTINE BLDUP (R,CO,DY)'/'C'/
     &   'C',5X,'MODEL DEPENDENT SUBROUTINE TO CALCULATE FORMATION'/
     &   'C',5X,'RATES OF NON-REACTING SPECIES'/
     &   'C'/6X,'REAL*4 R(*),CO(*),DY(*)'/'C')
C
        IF (NB.GT.0) CALL BLDUP
C
        WRITE (CODEU,'(6X,''RETURN''/6X,''END'')')
C
C
C       PRODUCE 'DIFUN' SUBROUTINE TO CALCULATE DERIVITES OF ACTIVE
C       SPECIES, ETC.  IF PRE-DIFFUN INSERTED CODE IS SPECIFIED, THE
C       ROUTINES TO PROCESS THESE ARE CALLED FIRST.  THIS WOULD
C       INCLUDE INSTRUCTIONS TO DEFINE VALUES OF COEFFICIENTS
C       WHICH DEPEND ON TIME-VARYING CONCENTRATIONS, RATE CONSTANTS,
C       ETC.
C
        WRITE (CODEU,702)
 702    FORMAT ('C'/'C'/6X,'SUBROUTINE DIFUN (C,A,S,RK,R,CO,FR,LR)'/
     &   'C'/'C',5X,'MODEL DEPENDENT SUBROUTINE TO CALCULATE RATES OF'/
     &   'C',5X,'FORMATION AND LOSS FOR ACTIVE SPECIES'/'C'/
     &   6X,'REAL*4 C(*),A(*),S(*),RK(*),R(*),CO(*),FR(*),LR(*)'/'C')
C
        IF (VCOCAL(2).OR.VCOCAL(5)) THEN
                WRITE (OUT,1201)
 1201           FORMAT ('0INITIAL DIFFUN CODE INSERTED:')
                WRITE (CODEU,1002)
                IF (VCOCAL(2)) THEN
C                                                     ! '.#D1'
                        ASSIGN 1204 TO IGO
                        IND=2
                        print*,'igo:',igo
                        GOTO 2000
                ENDIF
 1204           IF (VCOCAL(5)) THEN
C                                                     ! '.INS DIFF'
                        ASSIGN 1205 TO IGO
                        IND=5
                        print*,'igo:',igo
                        GOTO 2000
                ENDIF
        ENDIF
 1205   DIFNUM=0
        CALL DIFFUN
        WRITE (CODEU,'(6X,''RETURN''/6X,''END'')')
C
C       MODIFIED BY MIKE KLEEMAN (JUNE 2018)
C       PRODUCE 'CHKMOD' SUBROUTINE TO CHECK THE VALUES OF
C       NS, NACT, NSC, ...etc IN THE MOD FILE FOR COMPATIBILITY
C
        WRITE (CODEU,600)
 600    FORMAT ('C',/,'C',/,6X,'SUBROUTINE CHKMOD (',
     &              'NS,NACT,NSC,NSB,NDUM,NRXN,NCOEFV,'/
     &   5X,'&',19x,'NCOC1,NCOC2,NPHOTK,LENKBF)'/
     &   'C'/
     &   'C',5X,'COMPATIBILTY CHECK FOR MOD AND FORTRAN FILES'/
     &   'C')

        write(cname,'(i16)')ns-ndel
        call movlft(cname)
        write(codeu,710)'NS',cname(:nblank(cname)),'NS'

        write(cname,'(i16)')ndy
        call movlft(cname)
        write(codeu,710)'NACT',cname(:nblank(cname)),'NACT'

        write(cname,'(i16)')nsc
        call movlft(cname)
        write(codeu,710)'NSC',cname(:nblank(cname)),'NSC'

        write(cname,'(i16)')nb
        call movlft(cname)
        write(codeu,710)'NSB',cname(:nblank(cname)),'NSB'

        write(cname,'(i16)')ndum
        call movlft(cname)
        write(codeu,710)'NDUM',cname(:nblank(cname)),'NDUM'

        write(cname,'(i16)')nrxn
        call movlft(cname)
        write(codeu,710)'NRXN',cname(:nblank(cname)),'NRXN'

        write(cname,'(i16)')ncoefv
        call movlft(cname)
        write(codeu,710)'NCOEFV',cname(:nblank(cname)),'NCOEFV'

        write(cname,'(i16)')ncoc1
        call movlft(cname)
        write(codeu,710)'NCOC1',cname(:nblank(cname)),'NCOC1'

        write(cname,'(i16)')ncoc2
        call movlft(cname)
        write(codeu,710)'NCOC2',cname(:nblank(cname)),'NCOC2'

        write(cname,'(i16)')nphotk
        call movlft(cname)
        write(codeu,710)'NPHOTK',cname(:nblank(cname)),'NPHOTK'

        write(cname,'(i16)')lockbf
        call movlft(cname)
        write(codeu,710)'LENKBF',cname(:nblank(cname)),'LENKBF'

 710    format(6x,'IF(',(a),'.NE.',
     &        (a),')STOP "MOD FILE MISMATCH FOR ',(a),'"')

        write(codeu,'(a)') ' '
        WRITE (CODEU,'(6X,''RETURN''/6X,''END'')')


        CLOSE (UNIT=CODEU)
C                                        ! CLOSE (UNIT=UNIT1)
        CODEU=UNIT3
C
C       THIS CODE IS FOR VERSIONS OF PREP WHICH PREPARED OVERLAYED
C       "DIFFUN" SUBROUTINES.  THIS IS NOT DONE IN THIS VERSION,
C       BUT IS RETAINED FOR COMPATIBILITY.
        DIFTOT=DIFNUM
        DIFNUM=-1
C
C      DONE
C
 1900  NSUBS=3
       IF (U2OPEN) CLOSE (UNIT=UNIT2,STATUS='DELETE')
C
       DIFNUM=DIFTOT+1
       WRITE (OUT,6) NC
   6   FORMAT ('0CODE PREPARED.  LENGTH=',I7,' STATEMENTS')
C
       RETURN
C-------------------------------------------------------------------
C      INTERNAL SUBROUTINE TO READ COMMANDS FOR VCOEF OR INSCOD.
C
 2000  IF (.NOT.U2OPEN) THEN
              IOB160=' '
              CALL FILNAM (NLEN,IOB160,TMPUIC,'VCOPRP.TMP ',' ')
              OPEN (UNIT=UNIT2,NAME=IOB160,TYPE='OLD')
              U2OPEN=.TRUE.
       ELSE
              REWIND UNIT2
       ENDIF
C
 2010  IOB160=' '
       READ (UNIT2,2011,END=2090) J,L,(BYT160(I),I=1,L)
 2011  FORMAT (I1,Q,128A1,128A1)
       IF (J.NE.IND) GOTO 2010
       WRITE (OUT,2021) (BYT160(I),I=1,L)
 2021  FORMAT (' ',6X,80A1,:/(' ',4X,'& ',80A1))
       print*,'ind:',ind
       IF (IND.LE.3) THEN
              CALL VCOEF(L,IND)
       ELSE
              CALL INSCOD(L,IND-3)
       ENDIF
       GOTO 2010
 2090  GOTO IGO,(1004,1005,1204,1205)
C-------------------------------------------------------------------
       END
