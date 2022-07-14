C	SOM.FOR - the SOM subroutine
C
C
       subroutine som
C
C
C       OUTPUT FILES ARE:        
C
C       SUBROUTINES CALLED:
C
       include 'pspecs.inc'
       CHARACTER*1 RXNSTR(80)
       logical rxopn,newrxn,nxtnew
c
C
C
C******************************************************************************
C
C  create the new SOM species if they were used in a reaction
C
C******************************************************************************
C
C
C OPEN OUTPUT FILE
C
       IOB160=' '
       CALL FILNAM (LENAM,IOB160,LOGUIC,MODFIL,'.som ')
       OPEN (UNIT=TYP1,NAME=IOB160,status='UNKNOWN')
c
c print header to .pro file
c
       write(out,4400)
 4400  format(/,/,'These are the som reactions',/)
c
c initialize the trace reaction label to be a blank string
c
       do i = 1,maxrxn
        som_rxnlbl(i) = ' '
       enddo
c
c initialize somcheck to 0
c
       do k=1,maxns
       do no=0,maxsomo
       do nc=1,maxsomc
        somcheck(nc,no,k)=0
       enddo
       enddo
       enddo
c
c search for SOM species in all the reactions. 
c
       k = 0
       som_nrxnlbl = 0
       do 1000 while (k.lt.nsom)
        k=k+1
c
c find the length of the target SOM species name
c
        L = nblank(somname(k))
c
c loop over all rxns
c
        do 1100 irxn=1,nrxn
c skip reactions with the sequence "NSOM" in the rxnlbl
         if(index(rxnlbl(irxn),'NSOM').ne.0)goto 1100
c check the products (skip coeficients)
         do 1200 i=NPRODS(irxn),NPRODS(irxn+1)-1
          if(iprods(i).gt.0)then
           if(name(iprods(i))(:L).eq.somname(k)(:L) .and. 
     +       (name(iprods(i))(L+1:L+1).eq.' ' .or. 
     +        name(iprods(i))(L+1:L+1).eq.'_'))then
c current reaction requires SOM species!!!
            write(TYP1,3010) rxnlbl(irxn),somname(k),name(iprods(i)),
     +                       som_rxnlbl(irxn)
 3010       format (a8,' has ',a16,' as a product ',a16,(a))
            CALL SOMRXN(k, irxn, i) 
            goto 1100 
           endif
          endif ! iprods(i).gt.0
 1200    continue
 1100   continue

 1000  continue ! k.lt.nsom

c
c close the .som file
c
      CLOSE(TYP1)
c
c******************************************************************************
c 
c  rewrite the temp file used for the construction of the .mod file
c
c******************************************************************************
       IF (RXOPN) CLOSE(UNIT3)
       IOB160=' '
       CALL FILNAM (NLEN,IOB160,TMPUIC,'PREP.TMP ',' ')
       OPEN (UNIT=UNIT3,NAME=IOB160,FORM='UNFORMATTED')
       RXOPN=.TRUE.
c
        DO 100 IRXN=1,NRXN
         call write_rxn(irxn,.false.,.true.)
  100   CONTINUE
c
       CLOSE (UNIT=UNIT3,STATUS='KEEP')
       RXOPN=.FALSE.
c
c return to the calling subroutine
c
       return
       end
