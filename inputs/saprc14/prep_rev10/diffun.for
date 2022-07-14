C       FILE:  DIFFUN.FTN                       PREP SUBROUTINE
C
C       WRITTEN BY W.P.L. CARTER
C       UPDATED AND MAINTAINED BY S.E. HEFFRON
C       MAJOR REVISIONS TO PRODUCE FORTRAN CODE WRITTEN BY S.E. HEFFRON
C       LAST UPDATE:  W.P.L. CARTER  6/20/87
C
C
        SUBROUTINE DIFFUN
C
C       PRODUCES CODE TO CALCULATE RATES OF REACTIONS AND DERIVITIES,
C       RATES OF FORMATION, AND PSEUDO-UNIMOLECULAR RATES OF REACTION
C       FOR ACTIVE SPECIES ALSO CALCULATES STEADY STATE SPECIES CONC'S.
C
C       CALLED BY:  PNHRXN
C
C
C       SPECIFICATIONS FOR PREP PARMS, VARIABLES, AND ARRAYS IN COMMON.
C
        INCLUDE 'pspecs.inc'
C
C
C       THIS COMMON IS USED TO COMMUNICATE WITH SUBROUTINES CALLED
C       BY PNHRXN WHICH PRODUCE THE FORTRAN CODE.
C
        INTEGER   DIFNUM,DIFTOT,CODEU
        COMMON /CODE/ NC,NC1,DIFNUM,CODEU
C
C       SPECIFICATIONS FOR LOCAL VARIABLES AND ARRAYS.
C
        PARAMETER (MAXIRP=MAXRXN*4,MAXSS1=MAXSS+1)
        INTEGER   LRXBUF(MAXSS1),RXBUF(MAXIRP)
        INTEGER   DFSS(MAXSS1*MAXSS1),LDFSS(MAXSS1),NDFSS(MAXSS1)
C
        LOGICAL   LX,CHANGE,SSSET(MAXSS),REACTS,NEWSUM,COUPLED
        logical   react_and_prod
        BYTE DFRATE(MAXRXN)
C
        CHARACTER*1 PLUSCH
        CHARACTER*80  FTNSTR,TEMPSTR
        CHARACTER*16 CNUM,CNUM2
C
C
C
C       DEFINE RATES OF REACTIONS NOT INVOLVING STEADY STATE SPECIES,
C       AND INDICATE THESE BY DFRATE=1.
C
C       NOTE THAT REACTIONS WHOSE RATE IS USED AS THE RATE CONSTANT
C       FOR ANOTHER MUST NOT HAVE STEADY STATE REACTANTS.
C
C       FOR REACTIONS WITH STEADY STATE REACTANTS, DEFINE PARTIAL RATE
C       IF IT HAS OTHER REACTANTS, AND INDICATE THESE BY DFRATE=0.  IF
C       THE STEADY STATE SPECIE IS THE ONLY REACTANT, DO NOTHING YET,
C       EXCEPT INDICATE THIS BY DFRATE=-1.
C
        WRITE (CODEU,190)
 190    FORMAT ('C'/'C       DEFINE RATES OF REACTIONS, PARTIAL RATES',
     &             /'C       FOR REACTIONS WITH STEADY STATE REACTANTS'
     &             /'C')
        DO 110 IRXN = 1,NRXN
        NR = NRTOSR (IRXN)
        IF (NR.LT.0) THEN
C                                       ! DUMMY REACTION. NO RATE TO DEFINE.
                DFRATE(IRXN)=1
                GOTO 110
        ELSEIF (NR.EQ.1) THEN
C                                       ! SEE IF ONLY REACTANT SS
                IF (IRTOSR(1,IRXN).GE.NS1.AND.IRTOSR(1,IRXN).LE.NS2)
     &           THEN
                        DFRATE(IRXN)=-1
                        GOTO 110
                ENDIF
        ENDIF
        FTNSTR=' '
        WRITE(CNUM,'(I16)')IRXN
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR,200) CNUM(:IW)
 200    FORMAT (6X,'R(',(a),') = ')
        LS=12+IW
C                                               ! LENGTH OF STRING
        LKB=LKBUF(IRXN)
C                                               ! MULT K BY CONC'S
        IF (RXTYP(IRXN).EQ.0) THEN
                IF (RXTYP(LKB).EQ.7) THEN
                        WRITE(CNUM,'(I16)') IPHOTR(LKBUF(LKB))
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),210) CNUM(:IW)
 210                    FORMAT ('RK(',(A),')')
                ELSE
                        WRITE(CNUM,'(I16)') LKB
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),210) CNUM(:IW)
                ENDIF
        ELSEIF (RXTYP(IRXN).EQ.1) THEN
                IF (LKB.GT.IRXN.OR.DFRATE(LKB).NE.1) THEN
                        WRITE (OUT,*) 'ERROR.  RATE OF RXN'
     &                   ,LKB,' UNDEFINED.  CAN''T DEFINE '
     &                   ,'K FOR RXN',IRXN
                        NERR=NERR+1
                ENDIF
                WRITE(CNUM,'(I16)') LKB
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),215) CNUM(:IW)
C                                                      ! R(LKB)
 215            FORMAT ('R(',(a),')')
        ELSEIF (RXTYP(IRXN).EQ.7) THEN
                WRITE(CNUM,'(I16)') IPHOTR(LKB)
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),210) CNUM(:IW)
C                                                      ! RK(IPHOTR(LKB))
        ELSE
                WRITE(CNUM,'(I16)') IRXN
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),210) CNUM(:IW)
C                                                      ! RK(IRXN)
        ENDIF
        LS=LS+4+IW
        DFRATE(IRXN)=1
C                                       ! DFRATE=1 IF NO SS REACTANT FOUND
        DO 104 I=1,NR
        J=IRTOSR(I,IRXN)
        IF (J.GE.NS1.AND.J.LE.NS2) THEN
                IF (DFRATE(IRXN).EQ.0) THEN
C                                               ! SS FOUND.
                        WRITE (OUT,*) 'ERROR.  SS + SS REACTION',IRXN
                        NERR=NERR+1
                ELSE
                        DFRATE(IRXN)=0
                ENDIF
        ELSE
                IF (J.GT.0) THEN
C                                                       ! MULT BY CONC
                        IF (J.GT.N0) THEN
                                WRITE(CNUM,'(I16)') J-N0
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                                WRITE (FTNSTR(LS+1:),220) CNUM(:IW)
 220                            FORMAT ('*A(',(A),')')
                        ELSE
                                WRITE(CNUM,'(I16)') J
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                               WRITE (FTNSTR(LS+1:),221)CNUM(:IW)
 221                           FORMAT ('*C(',(a),')')
                        ENDIF
                        LS=LS+4+IW
                ELSEIF (-J.LE.MAXMAX) THEN
C                                                       ! OR MULT BY COEF
                        WRITE(CNUM,'(I16)') -J
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),222) CNUM(:IW)
 222                    FORMAT ('*CO(',(a),')')
                        LS=LS+5+IW
                ELSE
C                                                       ! OR BY RATE CONSTANT
                        JRXN=-J-MAXMAX
                        IF (RXTYP(JRXN).EQ.7) THEN
                                WRITE(CNUM,'(I16)') IPHOTR(LKBUF(JRXN))
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                                WRITE (FTNSTR(LS+1:),224) CNUM(:IW)
 224                            FORMAT ('*RK(',(a),')')
                                LS=LS+5+IW
                        ELSE
                                WRITE(CNUM,'(I16)') JRXN
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                                WRITE (FTNSTR(LS+1:),224) CNUM(:IW)
                                LS=LS+5+IW
                        ENDIF
                ENDIF
        ENDIF
  104   CONTINUE
        WRITE (CODEU,'(A)') FTNSTR(:LS)
        NC=NC+1
        NC1=NC1+1
  110   CONTINUE
C
C       PRODUCE CODE TO OBTAIN STEADY STATE CONCENTRATIONS.
C       LOOP THRU SS SPECIES.  PRODUCE CODE FOR SPECIE I IF
C       IT IS FORMED ONLY BY REACTION OF NON-SS SPECIES, OR OF
C       SS SPECIES FOR WHICH CODE HAVE ALREADY BEEN PRODUCED.
C       REACTIONS WHICH FORM SS SPECIES HAVE RATE DEFINED, DFRATE=1
C       MARKS THIS AS HAVING BEEN DEFINED.
C
        IF (NSS.EQ.0) GOTO 1050
C
C       DEFINE LRXBUF AND RXBUF, WHICH TELL WHICH REACTIONS
C       FORM A GIVEN STEADY STATE SPECIE
C
C       IN THIS VERSION, REACTIONS WHERE A SS SPECIES BOTH REACTS
C       AND IS FORMED IS ALLOWED, PROVIDED THAT THE STEADY STATE SPECIE IS
C       BOTH THE FIRST REACTANT AND THE FIRST PRODUCT.  IN THIS CASE, THE
C       CONCENTRATION OF THE SS SPECIE IS NOT DEPENDENT ON THE RATE OF
C       THE REACTION, SO IT IS NOT INCLUDED IN THE RXBUF ARRAY.
C
        WRITE (CODEU,240)
 240    FORMAT ('C'/'C      CALCULATE STEADY-STATE CONCENTRATIONS AND'
     &             /'C      RATES OF REACTIONS WITH SS SPECIES')
        NRX1=NRXN+1
        NSP=1
        ndp=1  ! next position to write in the ldfss vector 
        DO 1100 IX=1,NSS
        SSSET(IX)=.FALSE.
C                               ! INITIALIZE SS AS NOT 'SET'
        IR=IX+NA2
        print*,'mjk find rxns that produce:',name(ir)
        LRXBUF(IX)=NSP
        ndfss(ix)=0
        DO 1110 J=1,NRXN
        IP1=NPRODS(J)
        IP2=NPRODS(J+1)-1
        IF (IP2.LT.IP1) GOTO 1110
        REACTS=.FALSE.
        DO 1115 I=IP1,IP2
        IF (IPRODS(I).EQ.IR) THEN
           print*,'mjk rxn:',j,' ('//rxnlbl(j)//') produces '//name(ir)
           coef_react=0.0
           do ii=1,NRTOSR(J) ! loop over all the reactants
c check for reactions that also list the ss product as a reactant
            if(irtosr(ii,j).eq.ir)then ! ss product also a reactant
              reacts=.true.
              c1=1.0          ! default reactant stoich coef 
              if(ii.gt.1)then ! check for reactant stoich coeficient
                k=irtosr(ii-1,j)
                if(k.lt.0.and.-k.lt.maxcov)stop 'react vcoef notallowed'
                if(k.lt.0.and.-k.lt.maxco)c1=coef(-k)
              endif
              coef_react=coef_react+c1
c build a list of other ss reactants that form this ss product
            elseif(irtosr(ii,j).ge.ns1 .and. irtosr(ii,j).le.ns2)then
              do k=1,ndfss(ix)
               if(irtosr(ii,j).eq.dfss(ldfss(ix)+k-1))goto 555
              enddo
              print*,name(irtosr(ii,j)),' reacts to form ',name(ir)
              ndfss(ix)=ndfss(ix)+1
              if(ndfss(ix).eq.1)ldfss(ix)=ndp
              dfss(ndp)=irtosr(ii,j)
              ndp=ndp+1
 555          continue
            endif
           enddo ! ii=1,nrtosr(j)
c calc the net consumption for stoich, and do some common-sense checks
           if(reacts)then
             print*,'spec also a reactant with coef_react',coef_react
             do ii=i+1,ip2
               if(iprods(ii).eq.ir)then
                print*,'ss spec appears in react list and multiple'
                print*,'times in product list.  this challenges the'
                print*,'mech compiler too much.  please simplify the'
                print*,'product list for rxn ('//rxnlbl(j)//')'
                stop 'diffun marker mjk1'
               endif
             enddo
             tempstr=' '
             ltstr=1
             c2=1.0          ! default product stoich coef
             do kk=i-1,ip1,-1 ! check for product stoich coefficient
               k=iprods(kk)
               if(k.gt.0)goto 556 ! exit loop if normal spec found
c               if(k.lt.0.and.-k.lt.maxcov)stop 'prod vcoef notallowed'
               if(k.lt.0.and.-k.lt.maxcov)then
                WRITE(CNUM,'(I16)')-k
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                if(ltstr.eq.1)then
                 write(tempstr(ltstr:),500)cnum
 500             format("-CO(",(a),")")
                 ltstr=ltstr+5+iw
                else
                 write(tempstr(ltstr:),501)cnum
 501             format("*CO(",(a),")")
                 ltstr=ltstr+5+iw
                endif
               elseif(k.lt.0.and.-k.lt.maxco)then
                c2=c2*coef(-k)
               endif !if(k.lt.0.and.-k.lt.maxcov)
             enddo ! kk=i-1,ip1,-1
 556         continue
             if(c2.lt.0.0)then
              print*,'ss spec appears in react list and has neg'
              print*,'coef in the product list.  this looks like'
              print*,'an error.  please check rxn ('//rxnlbl(j)//')'
              stop 'diffun marker mjk2'
             endif
             if(ltstr.le.1)then
              coef_react=coef_react-c2
              if(coef_react.lt.0)then 
               print*,'net stoich coef:',coef_react
               print*,'ss spec appears in react list and has higher'
               print*,'production rate than consumption rate.  this '
               print*,'looks like an error.  please check rxn '
               print*,'('//rxnlbl(j)//')'
               stop 'diffun marker mjk3'
              endif
             else
              write(codeu,776)coef_react,tempstr(:ltstr),c2
 776          format(6x,"IF(",e14.7,(a),"*",e14.7,".LT.0.0)THEN")
              iiii=index(name(ir),' ')-1
              write(codeu,777)rxnlbl(j),name(ir)(:iiii)
 777          format(8x,'PRINT*,"rxn ',(a),' LISTS ',(a),
     +                  ' AS REACT AND PROD"')
              write(codeu,778)coef_react,tempstr(:ltstr),c2
 778          format(8x,'PRINT*,"rxn coef must be +ve: "',/
     +               8x,'PRINT*,',e14.7,(a),"*",e14.7)
              write(codeu,779)rxnlbl(j)
 779          format(8x,"STOP 'MJK MARKER ",(a),"'")
              write(codeu,'(6x,"ENDIF")')
             endif ! ltstr.gt.1
c if we get here, the species reacts and is formed, and all prods have
c been checked.  the species undergoes net consumption and is therefore
c not really produced by this reaction.  so skip to the next reaction.
c note that much of this will need to be redone when working with the 
c consumption rate for the species below.
             goto 1110
           ELSEIF (I.EQ.IP1) THEN
                   RXBUF(NSP)=J
           ELSEIF (IPRODS(I-1).GE.0) THEN
                   RXBUF(NSP)=J
           ELSE
                   RXBUF(NSP)=J+NRX1*(I-IP1)
           ENDIF
           NSP=NSP+1
        ENDIF
 1115   CONTINUE
C                       ! END LOOP OVER PRODUCTS
 1110   CONTINUE
C                       ! END LOOP OVER REACTIONS
        if(ndfss(ix).ge.1)
     +  print*,'SS species ',name(ir),' is formed by ',
     +        (name(dfss(ldfss(ix)+j-1)),j=1,ndfss(ix))
 1100   CONTINUE
C                       ! END LOOP OVER SS SPECIES
        LRXBUF(NSS+1)=NSP
C                               ! NSP USED LATER
        NSP1=NSP+1
C
        NSLOOP=0                ! NO SS ITERATIONS
        NSSET=0
C                               ! NO SS SPECIES SET
C
C       START OF ITERATIVE LOOP
C
 1002   NSLOOP=NSLOOP+1
        CHANGE=.FALSE.
C                               ! T IF SOME SS DEFINED DURING LOOP
C
        DO 1020 IX=1,NSS
        IF (SSSET(IX)) GOTO 1020
        I=IX+NA2
C
C       FIND OUT IF A REACTION FORMING I INVOLVES AN UN SSSET SS REACTANT
C
        coupled=.false.
        I1=LRXBUF(IX)
        I2=LRXBUF(IX+1)-1
C
        DO 1004 IP=I1,I2
        IRXN=RXBUF(IP)
        IF (IRXN.GT.NRXN) IRXN=IRXN-(IRXN/NRX1)*NRX1
        IF (DFRATE(IRXN).NE.1)THEN
c check for the special case of two coupled SS species. recall that no
c reaction is allowed to use two ss reactants, so all coupled cases will
c involve linear sets of equations with the form [A]{SS}={B} where [A]
c and {B} involve expressions where all the terms can be calcuated.  we
c will only deal with the two by two system since matrix inversion is 
c likely too expensive to be worthwhile.  so the equation set will be
c   a*SS1=b*SS2+c    with solution    SS1=(c*e+f*b)/(a*e-b*d)
c   e*SS2=d*SS1+f                     SS2=(a*f+d*c)/(a*e-b*d)
c where SS1 and SS2 are the two ss species, respectively, and a-f are
c the deterministic coefficients that can be calculated from the rate
c equations.
c
c count the number of unset ss spec that react to form ix
         notset1=0
         do k=1,ndfss(ix)
          jj=dfss(ldfss(ix)+k-1) ! this ss spec produces ix
          jx=jj-na2
          if(.not. ssset(jx))then
           notset1=notset1+1
           jnotset1=jx
          endif
         enddo ! k=1,ndfss(ix)
         if(notset1.le.0 .or. notset1.gt.1)goto 1020 ! abort this process
c count the number of unset ss spec that react to form jnotset1
         notset2=0
         do k=1,ndfss(jnotset1)
          jj=dfss(ldfss(jnotset1)+k-1) ! this ss spec produces jnotset1
          jx=jj-na2
          if(.not. ssset(jx))then
           notset2=notset2+1
           jnotset2=jx
          endif
         enddo ! k=1,ndfss(jx)
         if(notset2.le.0 .or. notset2.gt.1)goto 1020 ! abort this process
c check to see if ix and jnotset1 are coupled in a two by two system
         if(jnotset2.eq.ix)then
          coupled=.true.
          jj=jnotset1+na2
          print*,name(i),' and ',name(jj),' are coupled'
         else
          GOTO 1020   ! REACTION FORMING I NOT DEFINED.  QUIT
         endif
        ENDIF        ! DFRATE(IRXN).NE.1
C                                                    
 1004   CONTINUE
c
c if ix and jnotset1 are coupled, then treat them appropriately 
        if(coupled)then
         write(codeu,243)'deterministic solution terms for '//
     +                   name(i)//' and '//name(jj)
         ftnstr=' '
c         write(ftnstr,'(6x,"A1=CO(",I4,")")')indlil
         write(ftnstr,'(6x,"A1=0.")')
         ls=15
         NSSRXN=NSP
         nc_orig=nc
         DO 2220 IRXN=1,NRXN
         REACTS=.FALSE.
         IF (NRTOSR(IRXN).EQ.0) GOTO 2220
         DO 2225 J=1,NRTOSR(IRXN)
         IF (IRTOSR(J,IRXN).NE.I) GOTO 2225
         IF (REACTS) THEN
            WRITE (OUT,*) 'SS+SS REACTION NO.',IRXN
            NERR=NERR+1
            GOTO 2220
         ENDIF
         NSSRXN=NSSRXN+1
         RXBUF(NSSRXN)=IRXN
         REACTS=.TRUE.
         call check_react_and_prod_rxn(irxn,i,ftnstr,ls,react_and_prod)
         if(react_and_prod)goto 2220
C                                                
         lst=63
         if(nc-nc_orig.gt.20)lst=32
         call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
         if(ftnstr(:16).eq.'     &          '.and. 
     &      nc-nc_orig.gt.20)then
          ftnstr='      A1=A1 '
          ls=12
          nc_orig=nc
         endif
         IF (DFRATE(IRXN).EQ.-1) THEN ! IRXN IS A REACTION OF SS SPECIE I
            IF (RXTYP(IRXN).EQ.0) THEN
               JRXN=LKBUF(IRXN)
               IF (RXTYP(JRXN).EQ.7)JRXN=IPHOTR(LKBUF(JRXN)) 
               WRITE(CNUM,'(I16)')jrxn
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(...)
               LS=LS+5+iw
            ELSEIF (RXTYP(IRXN).EQ.1) THEN
               WRITE(CNUM,'(I16)')LKBUF(IRXN)
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),250) cnum(:iw) ! +R(...)
               LS=LS+4+iw
            ELSEIF (RXTYP(IRXN).EQ.7) THEN
               WRITE(CNUM,'(I16)')IPHOTR(LKBUF(IRXN))
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(IPHOTR(.))
                    LS=LS+5+iw
            ELSE
               WRITE(CNUM,'(I16)')IRXN
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(IRXN)
               LS=LS+5+iw
            ENDIF
         ELSE
            WRITE(CNUM,'(I16)')IRXN
            CALL MOVLFT(CNUM)
            IW=NBLANK(CNUM)
            WRITE (FTNSTR(LS+1:),250) cnum(:iw) ! +R(IRXN)
            LS=LS+4+iw
         ENDIF
 2225    CONTINUE ! END LOOP OVER REACTANTS
 2220    CONTINUE ! END LOOP OVER REACTIONS
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

         print*,'building B1'
         print*,'looking for reactions forming :',name(i)
         ftnstr=' '
         write(ftnstr,'(6x,"B1=0.")')
         LS=11
         nc_orig=nc
         do ip=i1,i2 ! all rxns forming ss spec i
          irxn=rxbuf(ip)
          if(irxn.gt.nrxn)then
           ICOEF=IRXN/NRX1
           IRXN=IRXN-ICOEF*NRX1
          endif
          print*,'dfrate for rxn:',irxn,rxnlbl(irxn),dfrate(irxn)
          if(dfrate(irxn).lt.1)then
           lst=64
           if(nc-nc_orig.gt.20)lst=32
           call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
           if(ftnstr(:16).eq.'     &          '.and. 
     &        nc-nc_orig.gt.20)then
            ftnstr='      B1=B1 '
            ls=12
            nc_orig=nc
           endif
           WRITE(CNUM,'(I16)')IRXN
           CALL MOVLFT(CNUM)
           IW=NBLANK(CNUM)
           write(ftnstr(ls+1:),250)cnum(:iw)  ! "+R(I4)"
           LS=LS+4+iw
           IP1=NPRODS(irxn)
           IP2=NPRODS(irxn+1)-1
           do ii=ip1,ip2
            iii=iprods(ii)
            if(iii.eq.i)then
             do kk=ii-1,ip1,-1 ! check for product stoich coefficient
              k=iprods(kk)
              print*,ii,k
              if(k.gt.0)goto 557 ! exit loop if normal spec found
              call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
              WRITE(CNUM,'(I16)')-K
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              WRITE (FTNSTR(LS+1:),222) cnum(:iw) ! "*CO(I4)"
              LS=LS+5+iw
             enddo ! kk=ii-1,ip1,-1
 557         continue
            endif ! ii.eq.i
           enddo ! ii=ip1,ip2
          endif ! dfrate(irxn).lt.1
         enddo ! ip=i1,i2
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

         print*,'building C1'
         print*,'looking for reactions forming :',name(i)
         ftnstr=' '
         write(ftnstr,'(6x,"C1=0.")')
         LS=11
         nc_orig=nc
         do ip=i1,i2 ! all rxns forming ss spec i
          irxn=rxbuf(ip)
          if(irxn.gt.nrxn)then
           ICOEF=IRXN/NRX1
           IRXN=IRXN-ICOEF*NRX1
          endif
          print*,'dfrate for rxn:',irxn,rxnlbl(irxn),dfrate(irxn)
          if(dfrate(irxn).eq.1)then
           lst=64
           if(nc-nc_orig.gt.20)lst=32
           call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
           if(ftnstr(:16).eq.'     &          '.and. 
     &        nc-nc_orig.gt.20)then
            ftnstr='      C1=C1 '
            ls=12
            nc_orig=nc
           endif
           WRITE(CNUM,'(I16)')irxn
           CALL MOVLFT(CNUM)
           IW=NBLANK(CNUM)
           write(ftnstr(ls+1:),250)cnum(:iw)  ! "+R(I4)"
           LS=LS+4+iw
           IP1=NPRODS(irxn)
           IP2=NPRODS(irxn+1)-1
           do ii=ip1,ip2
            iii=iprods(ii)
            if(iii.eq.i)then
             do kk=ii-1,ip1,-1 ! check for product stoich coefficient
              k=iprods(kk)
              print*,ii,k
              if(k.gt.0)goto 559 ! exit loop if normal spec found
              call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
              WRITE(CNUM,'(I16)')-K
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              WRITE (FTNSTR(LS+1:),222) cnum(:iw) ! "*CO(I4)"
              LS=LS+5+iw
             enddo ! kk=ii-1,ip1,-1
 559         continue
            endif ! ii.eq.i
           enddo ! ii=ip1,ip2
          endif ! dfrate(irxn).eq.1
         enddo ! ip=i1,i2
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

         print*,'building D1'
         ftnstr=' '
         write(ftnstr,'(6x,"D1=0.")')
         LS=11
         nc_orig=nc
         II1=LRXBUF(jnotset1)
         II2=LRXBUF(jnotset1+1)-1
         jj=jnotset1+na2
         print*,'looking for reactions forming :',name(jj)
         do ip=ii1,ii2 ! all rxns forming ss spec jnotset1
          irxn=rxbuf(ip)
          if(irxn.gt.nrxn)then
           ICOEF=IRXN/NRX1
           IRXN=IRXN-ICOEF*NRX1
          endif
          print*,'dfrate for rxn:',irxn,rxnlbl(irxn),dfrate(irxn)
          if(dfrate(irxn).lt.1)then
           lst=64
           if(nc-nc_orig.gt.20)lst=32
           call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
           if(ftnstr(:16).eq.'     &          '.and. 
     &        nc-nc_orig.gt.20)then
            ftnstr='      D1=D1 '
            ls=12
            nc_orig=nc
           endif
           WRITE(CNUM,'(I16)')irxn
           CALL MOVLFT(CNUM)
           IW=NBLANK(CNUM)
           write(ftnstr(ls+1:),250)cnum(:iw)  ! "+R(I4)"
           LS=LS+4+iw
           IP1=NPRODS(irxn)
           IP2=NPRODS(irxn+1)-1
           do ii=ip1,ip2
            iii=iprods(ii)
            if(iii.eq.jj)then
             do kk=ii-1,ip1,-1 ! check for product stoich coefficient
              k=iprods(kk)
              print*,ii,k
              if(k.gt.0)goto 561 ! exit loop if normal spec found
              call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
              WRITE(CNUM,'(I16)')-k
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              WRITE (FTNSTR(LS+1:),222) cnum(:iw) ! "*CO(I4)"
              LS=LS+5+iw
             enddo ! kk=ii-1,ip1,-1
 561         continue
            endif ! ii.eq.i
           enddo ! ii=ip1,ip2
          endif ! dfrate(irxn).eq.1
         enddo ! ip=ii1,ii2
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

c         write(ftnstr,'(6x,"E1=CO(",I4,")")')indlil
         write(ftnstr,'(6x,"E1=0.")')
         ls=17
         nc_orig=nc
         DO 2320 IRXN=1,NRXN
         REACTS=.FALSE.
         IF (NRTOSR(IRXN).EQ.0) GOTO 2320
         DO 2325 J=1,NRTOSR(IRXN)
         IF (IRTOSR(J,IRXN).NE.JJ) GOTO 2325
         IF (REACTS) THEN
            WRITE (OUT,*) 'SS+SS REACTION NO.',IRXN
            NERR=NERR+1
            GOTO 2320
         ENDIF
         REACTS=.TRUE.
         call check_react_and_prod_rxn(irxn,jj,ftnstr,ls,react_and_prod)
         if(react_and_prod)goto 2320
C                                                
         lst=63
         if(nc-nc_orig.gt.20)lst=32
         call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
         if(ftnstr(:16).eq.'     &          '.and. 
     &      nc-nc_orig.gt.20)then
          ftnstr='      E1=E1 '
          ls=12
          nc_orig=nc
         endif
         IF (DFRATE(IRXN).EQ.-1) THEN ! IRXN IS A REACTION OF SS SPECIE JJ
            IF (RXTYP(IRXN).EQ.0) THEN
               JRXN=LKBUF(IRXN)
               IF (RXTYP(JRXN).EQ.7)JRXN=IPHOTR(LKBUF(JRXN)) 
               WRITE(CNUM,'(I16)')JRXN
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(...)
               LS=LS+5+iw
            ELSEIF (RXTYP(IRXN).EQ.1) THEN
               WRITE(CNUM,'(I16)')LKBUF(IRXN)
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),250)cnum(:iw)  ! +R(...)
               LS=LS+4+iw
            ELSEIF (RXTYP(IRXN).EQ.7) THEN
               WRITE(CNUM,'(I16)')IPHOTR(LKBUF(IRXN))
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265)cnum(:iw) ! +RK(IPHOTR(.))
                    LS=LS+5+iw
            ELSE
               WRITE(CNUM,'(I16)')IRXN
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(IRXN)
               LS=LS+5+iw
            ENDIF
         ELSE
            WRITE(CNUM,'(I16)')IRXN
            CALL MOVLFT(CNUM)
            IW=NBLANK(CNUM)
            WRITE (FTNSTR(LS+1:),250) cnum(:iw) ! +R(IRXN)
            LS=LS+4+iw
         ENDIF
 2325    CONTINUE ! END LOOP OVER REACTANTS
 2320    CONTINUE ! END LOOP OVER REACTIONS
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

         print*,'building F1'
         print*,'looking for reactions forming :',name(jj)
         ftnstr=' '
         write(ftnstr,'(6x,"F1=0.")')
         LS=11
         nc_orig=nc
         II1=LRXBUF(jnotset1)
         II2=LRXBUF(jnotset1+1)-1
         jj=jnotset1+na2
         do ip=ii1,ii2 ! all rxns forming ss spec jnotset1
          irxn=rxbuf(ip)
          if(irxn.gt.nrxn)then
           ICOEF=IRXN/NRX1
           IRXN=IRXN-ICOEF*NRX1
          endif
          print*,'dfrate for rxn:',irxn,rxnlbl(irxn),dfrate(irxn)
          if(dfrate(irxn).eq.1)then
           lst=64
           if(nc-nc_orig.gt.20)lst=32
           call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
           if(ftnstr(:16).eq.'     &          '.and. 
     &        nc-nc_orig.gt.20)then
            ftnstr='      F1=F1 '
            ls=12
            nc_orig=nc
           endif
           WRITE(CNUM,'(I16)')IRXN
           CALL MOVLFT(CNUM)
           IW=NBLANK(CNUM)
           write(ftnstr(ls+1:),250)cnum(:iw)  ! "+R(I4)"
           LS=LS+4+iw
           IP1=NPRODS(irxn)
           IP2=NPRODS(irxn+1)-1
           do ii=ip1,ip2
            iii=iprods(ii)
            if(iii.eq.jj)then
             do kk=ii-1,ip1,-1 ! check for product stoich coefficient
              k=iprods(kk)
              print*,ii,k
              if(k.gt.0)goto 560 ! exit loop if normal spec found
              call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
              WRITE(CNUM,'(I16)')-K
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              WRITE (FTNSTR(LS+1:),222) cnum(:iw) ! "*CO(I4)"
              LS=LS+5+iw
             enddo ! kk=ii-1,ip1,-1
 560         continue
            endif ! ii.eq.i
           enddo ! ii=ip1,ip2
          endif ! dfrate(irxn).eq.1
         enddo ! ip=ii1,ii2
         call flushbuf(codeu,ls,ftnstr,6,nc,nc1)

c implement the solution for SS1=(c*e+f*b)/(a*e-b*d)
         write(codeu,'(6x,a)')'IF(C1*E1+F1*B1.GT.0.0)THEN'
         write(codeu,'(8x,a)')'IF(A1*E1-B1*D1.LE.0)THEN'
         write(codeu,'(10x,a)')'PRINT*,"coupled ss solution failed"'
         write(codeu,'(10x,a)')'PRINT*,"A1=",A1'
         write(codeu,'(10x,a)')'PRINT*,"B1=",B1'
         write(codeu,'(10x,a)')'PRINT*,"C1=",C1'
         write(codeu,'(10x,a)')'PRINT*,"D1=",D1'
         write(codeu,'(10x,a)')'PRINT*,"E1=",E1'
         write(codeu,'(10x,a)')'PRINT*,"F1=",F1'
         write(codeu,'(10x,a)')'stop "mjk marker1"'
         write(codeu,'(8x,a)')'ENDIF' 
         WRITE(CNUM,'(I16)')i-na2
         CALL MOVLFT(CNUM)
         IW=NBLANK(CNUM)
         WRITE(CNUM2,'(I16)')indlil
         CALL MOVLFT(CNUM2)
         IW2=NBLANK(CNUM2)
         write(codeu,246)cnum(:iw),cnum2(:iw2)
 246     FORMAT (8X,'S(',(a),')=(C1*E1+F1*B1)/(CO(',(a),
     +           ')+A1*E1-B1*D1)')
         write(codeu,'(6x,a)')'ELSE' 
         WRITE(CNUM,'(I16)')i-na2
         CALL MOVLFT(CNUM)
         IW=NBLANK(CNUM)
         write(codeu,247)cnum(:iw)
 247     FORMAT (8X,'S(',(a),')=0.0')
         write(codeu,'(6x,a)')'ENDIF' 
         goto 1121
        endif ! coupled
C
C       OK TO DEFINE SS SPECIE I, BECAUSE ALL REACTIONS FORMING IT
C       HAVE DFRATE=1.
C       CALCULATE SS CONC = FORMATION RATE/DESTRUCTION RATE
C
        WRITE (CODEU,243) NAME(I)
 243    FORMAT ('C'/'C',8X,A)
C                                          ! COMMENT LINE WITH NAME OF SPECIES
        FTNSTR=' '
        WRITE(CNUM,'(I16)')i-na2
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR,245) cnum(:iw)
C  MJM changed to FORMAT (6X,'S(',I4,') = (0.') from FORMAT (6X,'S(',I4,') = (')
 245    FORMAT (6X,'S(',(a),') = 0.')
C  MJM changed to LS=19 from LS=16 and NEWSUM=.FALSE. from NEWSUM=.TRUE.
C  marker1
        LS=14+iw
        NEWSUM=.FALSE.
C
        nc_orig=nc
        DO 1061 IP=I1,I2
        IRXN=RXBUF(IP)
        IF (IRXN.LE.NRXN) THEN
                WRITE(CNUM,'(I16)')irxn
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                IF (.NOT.NEWSUM) THEN
                        lst=63
                        if(nc-nc_orig.gt.20)lst=32
                        call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
                        if(ftnstr(:16).eq.'     &          '.and. 
     &                     nc-nc_orig.gt.20)then
                           WRITE(CNUM,'(I16)')i-na2
                           CALL MOVLFT(CNUM)
                           IW=NBLANK(CNUM)
                           WRITE (FTNSTR,2460) cnum(:iw),cnum(:iw)
 2460                      FORMAT (6X,'S(',(a),') = S(',(a),')')
                           ls=15+iw+iw
                           nc_orig=nc
                        endif
                        WRITE(CNUM,'(I16)')irxn
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),250)cnum(:iw)
 250                    FORMAT ('+R(',(a),')')
                        LS=LS+4+iw
                ELSE
                        WRITE(CNUM,'(I16)')irxn
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),252) cnum(:iw)
 252                    FORMAT ('R(',(a),')')
                        LS=LS+3+iw
                        NEWSUM=.FALSE.
                ENDIF
        ELSE
                ICOEF=IRXN/NRX1
                IRXN=IRXN-ICOEF*NRX1
C                                               ! MULTIPLY RATE BY COEF, IF
C                                               ! APPLICABLE.
C
                lst=63
                if(nc-nc_orig.gt.20)lst=32
                call flushbuf(codeu,ls,ftnstr,lst,nc,nc1)
                if(ftnstr(:16).eq.'     &          '.and. 
     &             nc-nc_orig.gt.20)then
                   WRITE(CNUM,'(I16)')i-na2
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                   WRITE (FTNSTR,2460) cnum(:iw),cnum(:iw)
                   ls=15+iw+iw
                   nc_orig=nc
                   newsum=.false.
                endif
                WRITE(CNUM,'(I16)')irxn
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                IF (.NOT.NEWSUM) THEN
C                                                    ! SUM ALREADY STARTED
                        WRITE (FTNSTR(LS+1:),250) cnum(:iw)
C                                                            ! +R(IRXN)
                        LS=LS+4+iw
                ELSE
                        WRITE (FTNSTR(LS+1:),252) cnum(:iw)
C                                                            ! R(IRXN)
                        LS=LS+3+iw
                ENDIF
c mjk - not sure if this code works correctly?
                DO 1066 J=ICOEF,1,-1
                KCOEF=IPRODS(NPRODS(IRXN)+J-1)
C                                                       ! NONZERO ICOEF IS
                IF (KCOEF.GT.0) GOTO 1067
C                                                       ! PREVIOUS REACTANT.
C                                                       ! IF <0, IT IS COEF
                call flushbuf(codeu,ls,ftnstr,62,nc,nc1)
                WRITE(CNUM,'(I16)')-KCOEF
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),222) cnum(:iw)
C                                                       ! MULT BY COEF
C  MJM changed to ...+9                LS=LS+8
                LS=LS+5+iw
 1066           CONTINUE
c mjk - not sure if this code works correctly?
 1067           NEWSUM=.FALSE.
        ENDIF
 1061   CONTINUE
        call flushbuf(codeu,ls,ftnstr,6,nc,nc1)
! END OF FORMATION RATE SUM
        WRITE(CNUM,'(I16)')i-na2
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR,2450) cnum(:iw),cnum(:iw)
 2450   FORMAT (6X,'S(',(a),') = S(',(a),')/')
        LS=16+iw+iw
C
C       NOW, FOR EACH REACTION DESTROYING I, FORM PSEUDO-UNIMOLECULAR
C       DESTRUCTION RATE.
C       NOTE THAT DILUTION ASSUMED NEGLIGABLE FOR STEADY STATE SPECIES.
C       BUT EACH HAS A MINIMUM DECAY RATE OF 'LITTLE' (COEF #INDLIL), TO
C       AVOID DIVIDE-BY-0 ERRORS WHICH MAY OCCUR IF IT IS NOT
C       INVOLVED IN ANY REACTIONS WITH NON-ZERO RATES.
C
C       NOTE THAT ANY REACTION DESTROYING I HAS DFRATE=0 OR -1
C
        call flushbuf(codeu,ls,ftnstr,62,nc,nc1)
        WRITE(CNUM,'(I16)')INDLIL
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR(LS+1:),260) cnum(:iw)
 260    FORMAT ('(CO(',(a),')')
        LS=LS+5+iw
C
C       LOOP OVER REACTIONS TO FIND REACTIONS OF I
C
        NSSRXN=NSP
C                       ! COUNTS HOW MANY REACTIONS OF I, STORED IN RXBUF
        DO 1120 IRXN=1,NRXN
c        print*,'checking rxn ',rxnlbl(irxn),' for react ',name(i)
        IR=NRTOSR(IRXN)
        REACTS=.FALSE.
        IF (IR.EQ.0) GOTO 1120
        DO 1125 J=1,IR
        IF (IRTOSR(J,IRXN).NE.I) GOTO 1125
        IF (REACTS) THEN
                WRITE (OUT,*) 'SS+SS REACTION NO.',IRXN
                NERR=NERR+1
                GOTO 1120
        ENDIF
        REACTS=.TRUE.
        NSSRXN=NSSRXN+1
        RXBUF(NSSRXN)=IRXN
c mjk check to see if this ss spec also appears in prod list.  this
c would mean that reaction is not handled in the conventional way.
c it would never have been marked as one producting the spec.
        call check_react_and_prod_rxn(irxn,i,ftnstr,ls,react_and_prod)
        if(react_and_prod)goto 1120
        call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
C                                               
        IF (DFRATE(IRXN).EQ.-1) THEN ! IRXN IS A REACTION OF SS SPECIE I
                IF (RXTYP(IRXN).EQ.0) THEN
                   JRXN=LKBUF(IRXN)
                   IF (RXTYP(JRXN).EQ.7)JRXN=IPHOTR(LKBUF(JRXN)) 
                   WRITE(CNUM,'(I16)')JRXN
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                   WRITE (FTNSTR(LS+1:),265) cnum(:iw)
 265               FORMAT ('+RK(',(a),')')
                   LS=LS+5+iw
                ELSEIF (RXTYP(IRXN).EQ.1) THEN
                   WRITE(CNUM,'(I16)')LKBUF(IRXN)
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                   WRITE (FTNSTR(LS+1:),250) cnum(:iw) ! +R(...)
                   LS=LS+4+iw
                ELSEIF (RXTYP(IRXN).EQ.7) THEN
                   WRITE(CNUM,'(I16)')IPHOTR(LKBUF(IRXN))
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                   WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(IPHOTR(.))
                        LS=LS+5+iw
                ELSE
                   WRITE(CNUM,'(I16)')IRXN
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                   WRITE (FTNSTR(LS+1:),265) cnum(:iw) ! +RK(IRXN)
                   LS=LS+5+iw
                ENDIF
        ELSE
                   WRITE(CNUM,'(I16)')IRXN
                   CALL MOVLFT(CNUM)
                   IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),250) cnum(:iw) ! +R(IRXN)
                LS=LS+4+iw
        ENDIF
 1125   CONTINUE ! END LOOP OVER REACTANTS
 1120   CONTINUE ! END LOOP OVER REACTIONS

! END OF DESTRUCTION RATE SUM
        WRITE (FTNSTR(LS+1:),'(A)')  ')'
        LS=LS+1
        WRITE (CODEU,'(A)') FTNSTR(:LS)
C                                            ! CONC OF SS SPECIES
        NC=NC+1
        NC1=NC1+1
C
C       SEE IF RXBUF ARRAY BOUNDS EXCEEDED, AND ABORT IF SO
C
 1121   CONTINUE
        IF (NSSRXN.GT.MAXIRP) THEN
                I=MAXIRP
                WRITE (OUT,*) 'MAXIMUM DIMENSION FOR RXBUF (',I,
     &           ') EXCEEDED IN DIFFUN.'
                STOP 'MODEL TOO BIG FOR RXBUF'
        ENDIF
C
C       INDICATE SPECIES AS 'SSSET', AND THAT A NEW ONE WAS PROCESSED.
C
        CHANGE=.TRUE.
        NSSET=NSSET+1
        SSSET(I-NA2)=.TRUE.
C
C       SS CONC NOW DEFINED AND IN "S"
C       DEFINE AND STORE RATES OF ALL ITS REACTIONS.
C
        DO 1035 IP=NSP1,NSSRXN
        IRXN=RXBUF(IP)
        FTNSTR=' '
        WRITE(CNUM,'(I16)')IRXN
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE(CNUM2,'(I16)')I-NA2
        CALL MOVLFT(CNUM2)
        IW2=NBLANK(CNUM2)
        WRITE (FTNSTR,280) cnum(:iw),cnum2(:iw2)
 280    FORMAT (6X,'R(',(a),') = S(',(a),') *')
C                                                   !"R(IRXN)=S(ISS)*"
C MJM changed to LS=25        LS=23
        LS=17+iw+iw2
        IF (DFRATE(IRXN).EQ.-1) THEN
                IF (RXTYP(IRXN).EQ.0) THEN
                        JRXN=LKBUF(IRXN)
                        IF (RXTYP(JRXN).EQ.7) THEN
                                WRITE(CNUM,'(I16)')IPHOTR(LKBUF(JRXN))
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                                WRITE (FTNSTR(LS+1:),210) cnum(:iw)
C                                                              ! RK(IPHOTR(...))
                        ELSE
                                WRITE(CNUM,'(I16)')JRXN
                                CALL MOVLFT(CNUM)
                                IW=NBLANK(CNUM)
                                WRITE (FTNSTR(LS+1:),210) cnum(:iw)
C                                                                     ! RK(JRXN)
                        ENDIF
                        LS=LS+4+IW
                ELSEIF (RXTYP(IRXN).EQ.1) THEN
                        WRITE(CNUM,'(I16)')LKBUF(IRXN)
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),215) cnum(:iw)
                        LS=LS+3+IW
C                                                                     ! R(...)
                ELSEIF (RXTYP(IRXN).EQ.7) THEN
                        WRITE(CNUM,'(I16)')IPHOTR(LKBUF(IRXN))
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),210)cnum(:iw) 
                        LS=LS+4+IW
C                                                                     ! RK(...)
                ELSE
                        WRITE(CNUM,'(I16)')IRXN
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),210) cnum(:iw)
                        LS=LS+4+IW
C                                                                     ! RK(IRXN)
                ENDIF
        ELSE
                WRITE(CNUM,'(I16)')IRXN
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),215) cnum(:iw)
                LS=LS+3+IW
C                                                                     ! R(IRXN)
        ENDIF
        WRITE (CODEU,'(A)') FTNSTR(:LS)
        NC=NC+1
        NC1=NC1+1
        DFRATE(IRXN)=1
 1035   CONTINUE
 1020   CONTINUE
C
C       END LOOP OVER SS SPECIES.  SEE IF ANOTHER LOOP NEEDED, OR IF
C       A LOOP WENT THRU W/O CHANGE, WHICH MEANS THAT THE SS SPECIES
C       CAN'T BE ORDERED.
C
        IF (NSSET.LT.NSS) THEN
                IF (CHANGE) THEN
                        GOTO 1002
                ELSE
                        WRITE (OUT,*) 'ERROR.  SOME SS SPECIES CAN''T  '
     &                   ,'BE ORDERED.  NO. ORDERED =',NSSET,'.  TOTAL'
     &                   ,' NO SS SPECIES =',NSS
                        NERR=NERR+1
                        DO IX=1,NSS
                           I=IX+NA2
                           IF (SSSET(IX))then
                              write(out,*)ix,i,name(i),' was set'
                           ELSE
                              write(out,*)ix,i,name(i),' was not set'
                           ENDIF
                        ENDDO
                ENDIF
        ENDIF
        WRITE (OUT,*) NSLOOP,' ITERATIONS REQUIRED TO ORDER '
     &   ,'STEADY STATE SPECIES.'
C
 1050  CONTINUE
C
C     END SS SECTION.  ALL SS CONC'S AND RATES OF REACTIONS FORMING
C     SS SPECIES DEFINED.  RATES OF ALL REACTIONS SHOULD BE DEFINED.
C       CHECK AND MAKE SURE THIS IS SO.
C
        DO 1900 IRXN = 1,NRXN
        IF (DFRATE(IRXN).NE.1) THEN
                WRITE (OUT,*) 'ERROR.  RATE OF REACTION',IRXN
     &           ,' UNDEFINED AFTER END OF SS LOOP. LABEL:'
     &            //rxnlbl(irxn)
                NERR=NERR+1
        ENDIF
 1900   CONTINUE
C
C
C     ALL RATES DEFINED.  DEFINE DERIVITIVES FOR ACTIVE SPECIES.
C
C       DEFINE TOTAL FORMATION RATE AND SAVE IN FR ARRAY
C
        WRITE (CODEU,282)
 282    FORMAT ('C'/'C       DEFINE FORMATION RATES FOR ACTIVE SPECIES')
C
        DO 2000 I=NA1,NA2
        NEWSUM=.TRUE.
        WRITE (CODEU,243) NAME(I)
C                                       ! COMMENT LINE WITH NAME OF SPECIES
        FTNSTR=' '
        WRITE(CNUM,'(I16)')I-N0
        CALL MOVLFT(CNUM)
        IW=NBLANK(CNUM)
        WRITE (FTNSTR,285) cnum(:iw)
 285    FORMAT (6X,'FR(',(a),') = ')
        LS=13+iw
C
C
C       LOOP OVER ALL REACTIONS.  ADD RATE IF IT IS A PRODUCT.
C
        LCOUNT = 1
        DO 2120 IRXN=1,NRXN
        IP1=NPRODS(IRXN)
        IP2=NPRODS(IRXN+1)-1
        IF (IP2.GE.IP1) THEN
             DO 2130 IX=IP1,IP2
             IF (IPRODS(IX).EQ.I) THEN         ! THIS REACTION FORMS IR
c--calculate the additional length to be added to the current line--
              WRITE(CNUM,'(I16)')IRXN
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              lsadd = 4+iw
              jx = ix
 2131         jx=jx-1   ! potential index location for a leading coefficient
              if(jx.ge.ip1.and.iprods(max(jx,ip1)).lt.0)then
               WRITE(CNUM,'(I16)')-iprods(max(jx,ip1))
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               lsadd = lsadd + 5 + iw
               goto 2131
              endif
c--dump the line and start over if needed--
              if(ls+lsadd.gt.71)then
               write(codeu,'(a)')ftnstr(:ls)
               nc = nc+1
               nc1=nc1+1
               lcount=lcount+1
               ftnstr='     &                   '
               ls=15
              endif
c--check for more than 20 lines--
              if(ftnstr(:25).eq.'     &                   ' .and.
     +                                           lcount.gt.20)then
               WRITE(CNUM,'(I16)')i-n0
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               write(ftnstr,286)cnum(:iw),cnum(:iw)
 286           format(6x,"FR(",(a),") = FR(",(a),")")
               ls = 17+iw+iw
               lcount = 1
              endif 
c--add a "+" to the current line if this is not the first entry--
              if(.not.newsum)then
               write(ftnstr(ls+1:),'(a)')'+'
               ls = ls+1
              else
               newsum=.false.
              endif
c--add the reaction to the current line--
              WRITE(CNUM,'(I16)')irxn
              CALL MOVLFT(CNUM)
              IW=NBLANK(CNUM)
              write(ftnstr(ls+1:),215)cnum(:iw)
              ls = ls + 3 + iw
              jx = ix
 2132         jx = jx-1 ! potential index location for a leading coefficient
              if(jx.ge.ip1.and.iprods(max(jx,ip1)).lt.0)then
               WRITE(CNUM,'(I16)')-iprods(max(jx,ip1))
               CALL MOVLFT(CNUM)
               IW=NBLANK(CNUM)
               if(ls+5+iw.gt.80)then
                print*,'error for rxn ',rxnlbl(irxn)
                print*,'output line exceeds 80 char in diffun.for'
                stop 'line too long in diffun.for'
               endif
               write(ftnstr(ls+1:),222)cnum(:iw)
               ls = ls + 5+iw
               goto 2132
              endif
             ENDIF
 2130        CONTINUE
        ENDIF
 2120   CONTINUE
C                                       ! END LOOP OVER REACTIONS
C
C       STORE FORMATION RATE
C
        IF (NEWSUM)  WRITE (FTNSTR(LS+1:),'(A)') '0'
        LS=LS+1
        WRITE (CODEU,'(A)') FTNSTR(:LS)
               if(ls.gt.72)stop 'marker 2 wrote bad line'
        NC=NC+1
        NC1=NC1+1
 2000   CONTINUE
C                                       ! END LOOP OVER SPECIES FOR FORMATION
C
C
C       LOOP OVER SPECIES FOR LOSS RATES.
C
        WRITE (CODEU,340)
 340    FORMAT ('C'/'C       DEFINE LOSS RATES OF ACTIVE SPECIES')
C
        DO 2001 I=NA1,NA2
C
C       SUM UP PSEUDO-UNIMOLECULAR LOSS RATES IN LR ARRAY
C       (IN THIS VERSION, FORM OUTPUT STRING IN FTNSTR)
C
        WRITE (CODEU,243) NAME(I)
C                                        ! COMMENT LINE WITH NAME OF SPECIES
C
        LX=.TRUE.
        PLUSCH=' '
C
C       LOOP OVER ALL REACTIONS.  SUBTRACT RATE IF I IS A REACTANT,
C
        DO 2200 IRXN=1,NRXN
        NR=NRTOSR(IRXN)
        IF (NR.GT.0) THEN
                NIX=0
                DO 2210 IX=1,NR
                IF (IRTOSR(IX,IRXN).EQ.I) THEN
                        NIX=NIX+1
                        IXX=IX
                ENDIF
 2210           CONTINUE
                IF (NIX.EQ.0) GOTO 2200
        ELSE
                GOTO 2200
        ENDIF
C
C       SPECIES I REACTS IN REACTION = IRXN.  CALCULATE PU DESTRUCTION
C       RATE.  IXX = LAST SPECIES I INDEX.  NIX = NO. TIMES SPECIE IN
C       REACTION.  FIRST WRITE START IF NONE WRITTEN.
C
        IF (LX) THEN
                LX=.FALSE.
                FTNSTR=' '
                WRITE(CNUM,'(I16)')I-N0
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR,345) cnum(:iw)
C                                        ! "LR(I-N0) = "
 345            FORMAT (6X,'LR(',(a),') = ')
                LS=17
        ENDIF
C
        LKB=LKBUF(IRXN)
C                                                       ! MULT K BY CONC'S
        call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
        IF (RXTYP(IRXN).EQ.0) THEN
                IF (RXTYP(LKB).EQ.7) THEN
                        WRITE(CNUM,'(I16)')IPHOTR(LKBUF(LKB))
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),267) PLUSCH,cnum(:iw)
C                                                                     ! +RK(...)
                ELSE
                        WRITE(CNUM,'(I16)')LKB
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),267) PLUSCH,cnum(:iw)
C                                                                     ! +RK(LKB)
                ENDIF
        ELSEIF (RXTYP(IRXN).EQ.1) THEN
                STOP 'K = RATE NOT SUPPORTED IN THIS VERSION'
        ELSEIF (RXTYP(IRXN).EQ.7) THEN
                WRITE(CNUM,'(I16)')IPHOTR(LKB)
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),267) PLUSCH,cnum(:iw)
C                                                                    ! +RK(...)
        ELSE
                WRITE(CNUM,'(I16)')IRXN
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),267) PLUSCH,cnum(:iw)
C                                                                    ! +RK(IRXN)
        ENDIF
 267    FORMAT (A1,'RK(',(a),')')
        PLUSCH='+'
        LS=LS+5+iw
        DO 2230 JX=1,NR
        IF (JX.EQ.IXX) GOTO 2230
        call flushbuf(codeu,ls,ftnstr,63,nc,nc1)
        J=IRTOSR(JX,IRXN)
        IF (J.GT.0) THEN
C                                                     ! MULT BY CONC
                IF (J.LT.NA1) THEN
                        WRITE(CNUM,'(I16)')J
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),221) cnum(:iw)
C 221                   FORMAT ('*C(',I4,')')
                        LS=LS+4+iw
                ELSEIF (J.GT.NA2) THEN
                        WRITE(CNUM,'(I16)')J-NA2
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),223)cnum(:iw)
 223                    FORMAT ('*S(',(a),')')
                        LS=LS+4+iw
                ELSE
                        WRITE(CNUM,'(I16)')J-N0
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),220)cnum(:iw)
C 220                   FORMAT ('*A(',I4,')')
                        LS=LS+4+iw
                ENDIF
        ELSEIF (-J.LE.MAXMAX) THEN
C                                                     ! OR MULT BY COEF
                WRITE(CNUM,'(I16)')-J
                CALL MOVLFT(CNUM)
                IW=NBLANK(CNUM)
                WRITE (FTNSTR(LS+1:),222) cnum(:iw)
C 222           FORMAT ('*CO(',I4,')')                  ! (DEFINED ABOVE)
                LS=LS+5+iw
        ELSE
C                                                     ! OR BY RATE CONSTANT
                JRXN=-J-MAXMAX
                IF (RXTYP(JRXN).EQ.7) THEN
                        WRITE(CNUM,'(I16)') IPHOTR(LKBUF(JRXN))
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),224)cnum(:iw)
C 224                   FORMAT ('*RK(',I4,')')          ! (DEFINED ABOVE)
                ELSE
                        WRITE(CNUM,'(I16)') JRXN
                        CALL MOVLFT(CNUM)
                        IW=NBLANK(CNUM)
                        WRITE (FTNSTR(LS+1:),224)cnum(:iw)
                ENDIF
                LS=LS+5+iw
        ENDIF
 2230   CONTINUE
        IF (NIX.GT.1) THEN
               IF (LS.GT.70) THEN
                       WRITE (CODEU,'(A)') FTNSTR(:LS)
                       NC=NC+1
                       NC1=NC1+1
                       FTNSTR='     &  '
                       LS=15
                ENDIF
                WRITE (FTNSTR(LS+1:),360) NIX
 360            FORMAT ('*',I1)
                LS=LS+2
        ENDIF
C
 2200   CONTINUE
C                       ! END LOOP OVER REACTIONS
C
C       STORE PU DESTRUCTION RATE
C
        IF (LX) THEN
                WRITE (CODEU,'(6X,''LR('',I4,'') = 0.0'')') I-N0
        ELSE
                WRITE (CODEU,'(A)') FTNSTR(:LS)
        ENDIF
C
        NC=NC+1
        NC1=NC1+1
C
 2001   CONTINUE
C                       ! END LOOP OVER ACTIVE SPECIES FOR LOSS
C
        RETURN
        END

        subroutine flushbuf(iounit,ls,ftnstr,lstarg,nc,nc1)
        integer iounit,ls,lstart
        character* (*) ftnstr
        if(ls.gt.lstarg)then
         write(iounit,'(a)')ftnstr(:ls)
         nc=nc+1
         nc1=nc1+1
         ftnstr='     &'
         ls=14
        endif
        return
        end

        subroutine check_react_and_prod_rxn(irxn,i,ftnstr,ls,
     +                                      react_and_prod)
c****************************************************************
c  written by: Mike Kleeman (April 2011)
c              UCD CEE
c  The purpose of this subroutine is to check for equations that
c  list species "i" as both a reactant and a product.  We want to 
c  add terms to output unit "codeu" to correctly describe this.
c
c  inputs:  irxn - reaction to check
c           i    - species to look for
c
c  outputs: codeu - output unit (passed via common block)
c           react_and_prod - true if "i" is a react and prod, 
c                            otherwise false
c***************************************************************
        include 'pspecs.inc'
        COMMON /CODE/ NC,NC1,DIFNUM,CODEU
        logical react_and_prod
        CHARACTER*80  FTNSTR, TEMPSTR

c start with the assumption that this reaction doesn't list the 
c target species "i" as both a reactant and product
        react_and_prod=.false.

c mjk check to see if this ss spec also appears in prod list.  this
c would mean that reaction is not handled in the conventional way.
c it would never have been marked as one producting the spec.
        IP1=NPRODS(IRXN)
        IP2=NPRODS(IRXN+1)-1
        DO II=IP1,IP2
         if(iprods(ii).eq.i)then ! ss spec is also a prod
c           print*,'reaction ',rxnlbl(irxn),' lists ',name(i),
c     +            'as a reactant and product'
c scan all the reactants
           coef_react=0.0        
           c1=1.0                ! default reactant stoich coef 
           do jj=1,nrtosr(irxn)  ! loop over all the reactants
            if(irtosr(jj,irxn).eq.i)then ! this index is target ss spec
              if(jj.gt.1)then ! check for reactant stoich coeficient
                k=irtosr(jj-1,irxn)
                if(k.lt.0.and.-k.lt.maxcov)stop 'react vcoef notallowed'
                if(k.lt.0.and.-k.lt.maxco)c1=coef(-k)
              endif
              coef_react=coef_react+c1
            endif ! irtosr(jj,irxn).eq.i
           enddo ! jj=1,ir
c scan the remaining products
           do jj=ii+1,ip2
             if(iprods(jj).eq.i)then
              print*,'ss spec appears in react list and multiple'
              print*,'times in product list.  this challenges the'
              print*,'mech compiler too much.  please simplify the'
              print*,'product list for rxn ('//rxnlbl(irxn)//')'
              stop 'diffun marker mjk5'
             endif
           enddo
           tempstr=' '
           ltstr=1
           c2=1.0              ! default product stoich coef
           do kk=ii-1,ip1,-1 ! check for product stoich coefficient
             k=iprods(kk)
             if(k.gt.0)goto 558 ! exit loop if normal spec found
c             if(k.lt.0.and.-k.lt.maxcov)stop 'prod vcoef notallowed'
             if(k.lt.0.and.-k.lt.maxcov)then
              if(ltstr.eq.1)then
               write(tempstr(ltstr:),500)-k
 500           format("-CO(",I4,")")
               ltstr=ltstr+9
              else
               write(tempstr(ltstr:),501)-k
 501           format("*CO(",I4,")")
               ltstr=ltstr+9
              endif
             elseif(k.lt.0.and.-k.lt.maxco)then
              c2=c2*coef(-k)
             endif !if(k.lt.0.and.-k.lt.maxcov)
           enddo ! kk=ii-1,ip1,-1
 558       continue
           if(c2.lt.0.0)then
            print*,'ss spec appears in react list and has neg'
            print*,'coef in the product list.  this looks like'
            print*,'an error.  please check rxn ('//rxnlbl(irxn)//')'
            stop 'diffun marker mjk6'
           endif
           if(ltstr.le.1)then
             coef_react=coef_react-c2
             if(coef_react.lt.0)then 
              print*,'ss spec appears in react list and has higher'
              print*,'production rate than consumption rate.  this '
              print*,'looks like an error.  please check rxn '
              print*,'('//rxnlbl(irxn)//')'
              print*,'net coef_react=',coef_react
              stop 'diffun marker mjk7'
             endif
           endif
c write the appropriate terms to the denominator for consumption rate
          print*,'mjk consumption rate for '//name(i)//' is special'
          call flushbuf(codeu,ls,ftnstr,62,nc,nc1)
          jrxn=irxn
          if(rxtyp(jrxn).eq.0)jrxn=lkbuf(jrxn)
          if(rxtyp(jrxn).eq.7)then
           write(ftnstr(ls+1:),265)iphotr(lkbuf(jrxn))
          else
           write(ftnstr(ls+1:),265)jrxn
          endif
 265      FORMAT ('+RK(',I4,')')
          ls=ls+9
          if(ltstr.le.1)then
            call flushbuf(codeu,ls,ftnstr,55,nc,nc1)
            write(ftnstr(ls+1:),266)coef_react
 266        format('*',e15.8)
            ls=ls+16
          else
            print*,'tempstr**'//tempstr(:ltstr)//'**'
            call flushbuf(codeu,ls,ftnstr,39-ltstr,nc,nc1)
            write(ftnstr(ls+1:),267)coef_react,tempstr(:ltstr),c2
 267        format('*(',e15.8,(a),'*',e14.7,')')
            ls=ls+18+ltstr+14
          endif!(ltstr.le.1)
          do jj=1,NRTOSR(IRXN)  ! loop over all the reactants
           k=IRTOSR(JJ,IRXN)
           if(k.ne.i)then ! this index is target ss spec
            call flushbuf(codeu,ls,ftnstr,62,nc,nc1)
            IF (k.GT.0) THEN                          ! MULT BY CONC
                    IF (k.LT.NA1) THEN
                            WRITE (FTNSTR(LS+1:),221) k
 221                       FORMAT ('*C(',I4,')')
                            LS=LS+8
                    ELSEIF (k.GT.NA2) THEN
                            WRITE (FTNSTR(LS+1:),223) k-NA2
 223                        FORMAT ('*S(',I4,')')
                            LS=LS+8
                    ELSE
                            WRITE (FTNSTR(LS+1:),220) k-N0
 220                       FORMAT ('*A(',I4,')')
                            LS=LS+8
                    ENDIF
            ELSEIF (-k.LE.MAXMAX) THEN                ! OR MULT BY COEF
                    WRITE (FTNSTR(LS+1:),222) -k
 222               FORMAT ('*CO(',I4,')')   
                    LS=LS+9
            ELSE                                      ! OR BY RATE CONSTANT
                    JRXN=-k-MAXMAX
                    IF (RXTYP(JRXN).EQ.7) THEN
                          WRITE (FTNSTR(LS+1:),224) IPHOTR(LKBUF(JRXN))
 224                     FORMAT ('*RK(',I4,')')  
                    ELSE
                            WRITE (FTNSTR(LS+1:),224) JRXN
                    ENDIF
                    LS=LS+9
            ENDIF ! k.gt.0 etc
           endif ! k.ne.i
          enddo ! jj=1,ir
          react_and_prod=.true.
          goto 1120 ! done with this reaction
         endif ! iprods(ii).eq.i
        ENDDO ! II=IP1,IP2
 1120   continue
c mjk end of code insert to check for ss spec in prod list.

c return to the calling subroutine
        return
        end
