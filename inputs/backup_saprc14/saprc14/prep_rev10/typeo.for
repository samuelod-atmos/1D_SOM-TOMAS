C	TYPEO.FOR - the type-o-matic subroutine
C
C
       subroutine typeo
C
C
C
C       OUTPUT FILES ARE:         input-name.typ    TYPED MECHANISM LISTING
C
C
C       SUBROUTINES CALLED:
C       TYPRXN  ... processes selected rection for typing
C
C               SYSTEM ROUTINES:
C
C               "NEWSUBS" (SAPRC) LIBRARY UTILITY ROUTINES
C
C       FILEIN  ... ROUTINE TO OPEN INPUT FILE
C       FILNAM  ... ROUTINE TO CREATE FILE NAME
C       INBUF   ... GET LINE FROM INPUT FILE.  GIVES PROMPT IF INPUT
C                   FROM USER REQUIRED.
C
       include 'pspecs.inc'
C
C
       INTEGER k,irxn,i,itype,orgns,L,L2
       CHARACTER*1 RXNSTR(80)
       character*16 aname
       logical rxopn,newrxn,nxtnew
c
       newrxn = .true.
C
C******************************************************************************
C
C  Perform error checking to ensure the user specified species containing
C     the specified atoms to track 
C
C******************************************************************************
C
       inquire(unit=unit3,opened=rxopn)
       IF (.NOT.RXOPN) THEN
                IOB160=' '
                CALL FILNAM (NLEN,IOB160,TMPUIC,'PREP.TMP ',' ')
                OPEN (UNIT=UNIT3,NAME=IOB160,TYPE='OLD',
     &           FORM='UNFORMATTED')
                RXOPN=.TRUE.
       ELSE
        rewind (unit3)
       ENDIF
c
c "fast-forward" through the PREP.TMP file
c
 100    READ (UNIT3,END=110) NXTNEW,L2,(RXNSTR(I),I=1,L2)
        GOTO 100
 110    CONTINUE
C
C
C******************************************************************************
C
C  perform the typing: loop over all the species in tname so that all
C     the specified atoms in tatom are typed.  master loop k=1,nst
C
C******************************************************************************
C
c keep track of the original # of species
c
       orgns = ns
C
C OPEN OUTPUT FILE
C
       IOB160=' '
       CALL FILNAM (LENAM,IOB160,LOGUIC,MODFIL,'.typ ')
       OPEN (UNIT=TYP1,NAME=IOB160,status='UNKNOWN')
c
c print header to .pro file
c
       write(out,4400)
 4400  format(/,/,'These are the expanded type reactions',/)
c
c initialize the trace reaction label to be a blank string
c
       do i = 1,maxrxn
        typ_rxnlbl(i) = ' '
       enddo
c
c search for traced species in all the reactions.  note that the 
c traced list can be updated in TYPRXN.
c
       k = 0
       typ_nrxnlbl = 0
       do 1000 while (k.lt.nst)
        k=k+1
c
c write information to typing file and screen
c
        write(icrt,3007) tname(k),tatom(k),k,nst,nrxn
 3007   format ('Tracing species ',a16,' for ',a1,
     +           i3,' of ',i3,' with nrxn= ',i8)
        write(typ1,4200) tname(k),tatom(k),ttypes(k)
 4200   format('Species selected for typing: name, atom, sources:',
     +          2x,a16,1x,a1,4x,i4)
        write(typ1,4210)
 4210   format('Reactions that must be typed')
c
c find the length of the target trace species name
c
        L = nblank(tname(k))
c
c loop through the reactions to see if tname is in products of irxn.  note
c that we skip reactions with the sequence "NT" in the rxnlbl.
c
        do 1100 irxn=1,nrxn
c         print*,'irxn,nrxn:',irxn,nrxn
c skip reactions with the sequence "NT" in the rxnlbl
         if(index(rxnlbl(irxn),'NT').ne.0)goto 1100
c check the products (skip coeficients)
         do 1200 i=NPRODS(irxn),NPRODS(irxn+1)-1
          if(iprods(i).gt.0)then
           if(name(iprods(i))(:L).eq.tname(k)(:L) .and. 
     +       (name(iprods(i))(L+1:L+1).eq.' ' .or. 
     +        name(iprods(i))(L+1:L+1).eq.'_'))then
c current reaction requires typing!!!
            write(TYP1,3010) rxnlbl(irxn),tname(k),name(iprods(i)),
     +                       typ_rxnlbl(irxn)
 3010       format (a8,' has ',a16,' as a product ',a16, a8)
            CALL TYPRXN(k, irxn) 
            goto 1100 
           endif
          endif
 1200    continue
c check the reactants (skip coeficients)
         do 1300 i=1,nrtosr(irxn)
          if(irtosr(i,irxn).gt.0)then
           if(name(irtosr(i,irxn))(:L).eq.tname(k)(:L) .and. 
     +       (name(irtosr(i,irxn))(L+1:L+1).eq.' ' .or. 
     +        name(irtosr(i,irxn))(L+1:L+1).eq.'_' .or.
     +        name(irtosr(i,irxn))(L+1:L+1).eq.'['))then
c current reaction requires typing!!!
            write(TYP1,3020) rxnlbl(irxn),tname(k),name(irtosr(i,irxn)),
     +                       typ_rxnlbl(irxn)
 3020       format (a8,' has ',a16,' as a reactant ',a16, a8) 
            CALL TYPRXN(k, irxn) 
            goto 1100
           endif
          endif
 1300    continue
 1100   continue

 1000  end do ! bottom of the typing loop
c
       write(typ1,9500)
 9500  format(/,'These are the expanded species:')
         write(typ1,9600) (name(i), i=orgns+1,ns)
 9600  format(4(4x,a16))
c
c
        WRITE (typ1,2002) (I,SPTYP(I),NAME(I),I=1+orgns,NS)
 2002   FORMAT ('0LIST OF REACTANTS INPUT'//'         DECLARED TYPES:'
     &   ,'  1 = CONSTANT,  2 = BUILDUP,  3 = DUMMY, 4 = ACTIVE,  5 ='
     &   ,' ST. STATE,  6 = UNSPECIFIED'//(' ',4(I6,') [',I1,'] ',A16)))
c
c check for reaction templates.  note that some reactions may be
c deleted or changed to vanilla versions.
c
      CALL CHECK_FOR_TEMPLATES
c
c close the .typ file
c
      CLOSE(TYP1)
c
c close the temp file
c
      CLOSE (UNIT=UNIT3,STATUS='KEEP')
c
c return to the calling program
c
       return
       end
C
C end of file typeo.for
