C	TYPEOSUBS.FOR
C
C
C       Subroutines included to support the TYPEO.for procedures
C
C		TSPCNAM - map typed specied to array index, add if new
C		TPNTRXN - printout the typed reactions to name.pro2 file
c		TYPRXN  - type the reaction
C
C
C
C -----------------------------------------------------------------------------
c
c      Check to see if reaction irxn contains species tname(k), hence typing
C
       subroutine TYPRXN(k,irxn)
       include 'pspecs.inc'
       character aname*16, sname*21, suffix*30, sufgen*30
       character a1*1, suffix2*30, suffix3*30, sufext*30
       integer i,j,k,irxn
       integer isrc_p(0:9)
       real    tno(0:100),tsno(0:100)
       logical newrxn, use_src_template
c
c copy the traced atom name to the working variable
c
       write(aname,'(a16)') tatom(k)
       call movlft(aname)
       a1=aname(:1)
c
c check if this reaction has already been typed for the target atom
c
      if(index(typ_rxnlbl(irxn),aname(:1)).ne.0)return
c
c check if this reaction uses a source template map and initialize
      call check_src_template(irxn,aname,icnt_src,use_src_template)
c
c calculate the number of traced atoms in the reactant list, and count
c the species with traced atoms.  tr is the number of traced atoms, while
c trs is the maximum number of sources that those atoms could have come
c from.
c
      tr = 0.0
      trs = 0.0
      nr = 0
      c1 = 1.0
      icount = 0
      do i = 1,nrtosr(irxn)
       icount = icount + 1
       if(irtosr(i,irxn).gt.0)then
c contribution from this reactant
        if(aname(:1).eq.'N')then
         tr = tr + c1*nno(irtosr(i,irxn))
         trs = trs + c1*nsno(irtosr(i,irxn))
         tno(icount) = nno(irtosr(i,irxn))
         tsno(icount) = nsno(irtosr(i,irxn))
        elseif(aname(:1).eq.'S')then
         tr = tr + c1*sno(irtosr(i,irxn))
         trs = trs + c1*ssno(irtosr(i,irxn))
         tno(icount) = sno(irtosr(i,irxn))
         tsno(icount) = ssno(irtosr(i,irxn))
        elseif(aname(:1).eq.'C')then
         tr = tr + c1*cno(irtosr(i,irxn))
         trs = trs + c1*csno(irtosr(i,irxn))
         tno(icount) = cno(irtosr(i,irxn))
         tsno(icount) = csno(irtosr(i,irxn))
        elseif(aname(:1).eq.'O')then
         tr = tr + c1*ono(irtosr(i,irxn))
         trs = trs + c1*osno(irtosr(i,irxn))
         tno(icount) = ono(irtosr(i,irxn))
         tsno(icount) = osno(irtosr(i,irxn))
        elseif(aname(:1).eq.'X')then
         tr = tr + c1*xno(irtosr(i,irxn))
         trs = trs + c1*xsno(irtosr(i,irxn))
         tno(icount) = xno(irtosr(i,irxn))
         tsno(icount) = xsno(irtosr(i,irxn))
        endif
c count the number of reactants with traced atoms
        if(tno(icount).gt.0)nr = nr + 1
c fractional coefficients indicate that this is not an elementary rxn
c this rxn will need to be written as elementary steps for typing 
        if((amod(c1,1.0).ne.0.0.or.c1.eq.0) .and. tno(icount).gt.0)then
         if(.not. use_src_template)then
         print*,'non-integer coefficient for rxnlbl:',rxnlbl(irxn)
         print*,'coef,species:',c1,name(irtosr(i,irxn))
         print*,'rewrite rxn as elementrary steps for typing'
         print*,'only CONSTANT int coefs are allowed for typed reac'
         stop 'marker1'
         endif ! .not. use_src_template
        endif
        c1 = 1.0
       else
c possible leading coefficient
        if(-irtosr(i,irxn).lt.maxco)c1 = c1*coef(-irtosr(i,irxn))
        tno(icount) = 0
       endif
      enddo
c
c calculate the number of traced atoms in the product list.  also
c detect the case where we put a reactant in the product list with 
c a negative coefficient.  tp is the number of traced atoms, while 
c tps is the maximum number of sources that those traced atoms could
c have come from.
c
      tp = 0.0
      tps = 0.0
      c1 = 1.0
      do i = nprods(irxn),nprods(irxn+1)-1
       icount = icount + 1
c check the dimensions of icount
       if(icount.gt.100)then
        print*,'icount has exceeded 100 for reaction:',irxn,
     +      ' '//rxnlbl(irxn)
        print*,'reactants:'
        do j = 1,nrtosr(irxn)
         print*,j,name(irtosr(j,irxn))
        enddo
        print*,'products:'
        do j = nprods(irxn),nprods(irxn+1)-1
         jj = iprods(j)
         if(jj.lt.0 .and.-jj.le.maxcov)then
          print*,'coeff preceding next spec:',coefnm(-jj),coef(-jj)
         elseif(jj.lt.0 .and. -jj.le.maxco)then
          print*,'coeff preceding next spec:',coef(-jj)
         else
          print*,name(jj)
         endif
        enddo
       endif
       if(iprods(i).gt.0)then
c contribution from this product. note that a coefficient of -1.0 is
c allowed to signify that this is actually a reactant.  coefficients
c less than -1.0 are not allowed.
        if(c1.lt.0)then
         print*,'negative product coefficient detected for rxnlbl:',
     +           rxnlbl(irxn),c1
         if(c1.ne.-1)then
          if(index(rxnlbl(irxn),'-'//aname(:1)).ne.0
     +        .and. int(tr).eq.1)then
           print*,'tag from first typed reactant passed to all products'
           print*,'because reaction label contains: -'//aname(:1)
          else
           nnn=1 ! assume species with -ve coef being tracked
           if(aname(:1).eq.'N')then
            nnn=nno(iprods(i))
           elseif(aname(:1).eq.'S')then
            nnn=sno(iprods(i))
           elseif(aname(:1).eq.'C')then
            nnn=cno(iprods(i))
           elseif(aname(:1).eq.'O')then
            nnn=ono(iprods(i))
           elseif(aname(:1).eq.'X')then
            nnn=xno(iprods(i))
           else
            print*,'dont know how to check for atom: '//aname(:1)
           endif
           if(nnn.gt.0)then
            print*,'the only negative value allowed is -1.0'
            stop 'typeosubs marker 0'
           endif ! nnn.gt.0
          endif ! check for -X flag in rxnlbl and tr=1
         endif ! c1.ne.-1
        endif ! c1.lt.0
        if(aname(:1).eq.'N')then
         if(c1.gt.0)then
          tp = tp + c1*nno(iprods(i))
          tps = tps + c1*nsno(iprods(i))
         else
          tr = tr - c1*nno(iprods(i))
          trs = trs - c1*nsno(iprods(i))
         endif
         tno(icount) = nno(iprods(i))
         tsno(icount) = nsno(iprods(i))
        elseif(aname(:1).eq.'S')then
         if(c1.gt.0)then
          tp = tp + c1*sno(iprods(i))
          tps = tps + c1*ssno(iprods(i))
         else
          tr = tr - c1*sno(iprods(i))
          trs = trs - c1*ssno(iprods(i))
         endif
         tno(icount) = sno(iprods(i))
         tsno(icount) = ssno(iprods(i))
        elseif(aname(:1).eq.'C')then
         if(c1.gt.0)then
          tp = tp + c1*cno(iprods(i))
          tps = tps + c1*csno(iprods(i))
         else
          tr = tr - c1*cno(iprods(i))
          trs = trs - c1*csno(iprods(i))
         endif
         tno(icount) = cno(iprods(i))
         tsno(icount) = csno(iprods(i))
        elseif(aname(:1).eq.'O')then
         if(c1.gt.0)then
          tp = tp + c1*ono(iprods(i))
          tps = tps + c1*osno(iprods(i))
         else
          tr = tr - c1*ono(iprods(i))
          trs = trs - c1*osno(iprods(i))
         endif
         tno(icount) = ono(iprods(i))
         tsno(icount) = osno(iprods(i))
        elseif(aname(:1).eq.'X')then
         if(c1.gt.0)then
          tp = tp + c1*xno(iprods(i))
          tps = tps + c1*xsno(iprods(i))
         else
          tr = tr - c1*xno(iprods(i))
          trs = trs - c1*xsno(iprods(i))
         endif
         tno(icount) = xno(iprods(i))
         tsno(icount) = xsno(iprods(i))
        endif
c fractional coefficients indicate that this is not an elementary rxn
c this rxn will need to be written as elementary steps for typing unless
c we only have a single typed reactant 
        if(.not. use_src_template)then
        if((amod(c1,1.0).ne.0.0.or.c1.eq.0) .and. 
     +                     tno(icount).gt.0 .and. nr.gt.1)then
         print*,'non-integer coefficient for rxnlbl:',rxnlbl(irxn),nr
         print*,'coef,name:',c1,name(iprods(i))
         if(index(rxnlbl(irxn),'-'//aname(:1)).ne.0
     +   .and. amod(tr,1.0).eq.0 .and. int(tr).eq.1 .and. tp.gt.tr)then
          print*,'tag from first typed reactant passed to all products'
          print*,'because reaction label contains: -'//aname(:1)
         else
          print*,'rewrite rxn as elementrary steps for typing'
          stop 'marker2'
         endif
        endif
        endif ! .not. use_src_template
        c1 = 1.0
       else
c possible leading coefficient
        if(-iprods(i).le.maxco .and. -iprods(i).gt.maxcov)then
         c1 = c1*coef(-iprods(i))
        else
         c1 = 0
        endif
        tno(icount) = 0
        tsno(icount) = 0
       endif
      enddo
c
c make sure the reactant and product sums are equal
c
      print*,'trs:',trs
      call flush(6)
      if(tr.le.0.0 .or. amod(tr,1.0).ne.0.0 .or. tp.gt.tr)then
       print*,'mismatch for reaction tracing:',rxnlbl(irxn),
     +        ' '//aname(:1)
       print*,'tr, tp:',tr,tp
c check for -N, -S, -O, -X in reaction label to indicate that tag for
c first reactant should be passed to all products
       if(index(rxnlbl(irxn),'-'//aname(:1)).ne.0
     +   .and. amod(tr,1.0).eq.0 .and. int(tr).eq.1 .and. tp.gt.tr)then
         print*,'tag from first typed reactant passed to all products'
         print*,'because reaction label contains: -'//aname(:1)
c check for +N, +S, +O, +X in reaction label and [_N [_S [_O [_X in the 
c reactant or product list to indicate that custom src template map is used
       elseif(use_src_template)then
         print*,'ok because custom tag passing based on template rxn'
         print*,'with names R[_'//a1//'] + R[_'//a1//'] = ...'
       else
         print*,'this may be caused by variable coefs used as '//
     +          'product stoichiometric values'
         print*,'only constant integer values can be stoichiometric '//
     +          'values'
         print*,'reactant, tno:'
         icount = 0
         do i = 1,nrtosr(irxn)
          icount = icount + 1
          j = irtosr(i,irxn)
          if(j.gt.0)then
           print*,name(j),tno(icount)
          elseif(-j.le.maxcov)then
           print*,'constant to modify K',coefnm(-j)
          elseif(-j.le.maxco)then
           print*,'constant to modify K',coef(-j)
          endif
         enddo
         print*,'product, tno:'
         do i = nprods(irxn),nprods(irxn+1)-1
          icount = icount + 1
          j = iprods(i)
          if(j.lt.0 .and.-j.le.maxcov)then
           print*,'coeff preceding next spec:',coefnm(-j),coef(-j)
          elseif(j.lt.0 .and. -j.le.maxco)then
           print*,'coeff preceding next spec:',coef(-j)
          else
           print*,name(j),tno(icount)
          endif
         enddo
         print*,'if there is 1 reactant tag and you want to pass it '
         print*,'to all products then put -'//aname(:1)//
     +          ' in the rxn label '//rxnlbl(irxn)
         stop 'marker1'
       endif
      endif
c write a warning if this reaction will discard traced information
      if(int(tp).le.0)then
       print*,'warning - traced information discarded for rxn:',
     +         rxnlbl(irxn), aname(:1)
      endif
c
c check the typed reactants to make sure they exist; add them if not.
c note - we add the typed species to the regular list, and the 
c base species to the traced list
c
      icount = 0
      do i = 1,nrtosr(irxn)
       icount = icount + 1
       if(irtosr(i,irxn).gt.0)then
        call chktspc2(aname,tsno(icount),irtosr(i,irxn),k)
        if(tsno(icount).gt.0)then
         call tspcnam(isp,name(irtosr(i,irxn)),.true.)
         tatom(isp) = tatom(k)
         ttypes(isp) = ttypes(k)
        endif
       endif
      enddo
c
c check the typed products to make sure they exist; add them if not.
c note - we add the typed species to the regular list, and the 
c base species to the traced list
c
      do i = nprods(irxn),nprods(irxn+1)-1
       icount = icount + 1
       if(iprods(i).gt.0)then
        call chktspc2(aname,tsno(icount),iprods(i),k)
        if(tno(icount).gt.0)then
         call tspcnam(isp,name(iprods(i)),.true.)
         tatom(isp) = tatom(k)
         ttypes(isp) = ttypes(k)
        endif
       endif
      enddo
c
c go through and type this reaction. note that we skip the first
c index because this is the vanilla reaction.  nt1 is the total
c number of sources to type.  nt is the total number of reactions
c generated.
c
      nt1 = ttypes(k) + 1
      nt = nt1**int(trs)
      print*,'typing reaction :',irxn,tr,trs,nt
      call flush(6)
      do 10 i = 2,nt ! top of the added reaction loop
c
c index the reaction count and check to see if it exceeds maxrxn
c
       nrxn = nrxn + 1
       if(nrxn.gt.maxrxn)then
        write(out,*)'TOO MAN REACTIONS.  MAX=',maxrxn
        stop 'too many reactions in typrxn. increase maxrxn'
       endif
c
c set the new reaction label as a counter followed by the atom symbol
c note that the primary rxnlbl has the added reaction count followed by
c the traced atom.  the secondary typ_rxnlbl only contains the traced atom
c so that we can quickly determine if the reaction has been typed.
c
c primary rxn label for new rxn
       typ_nrxnlbl = typ_nrxnlbl + 1
       write(rxnlbl(nrxn),'(i5,a1)')typ_nrxnlbl,aname(:1)
       call movlft2(rxnlbl(nrxn))
c       print*,'new rxn lbl:**'//rxnlbl(nrxn)//'**'
c secondary rxn label for new rxn
       i1 = max(nblank(typ_rxnlbl(nrxn)),1)
       if(i1.gt.0)then
        write(typ_rxnlbl(nrxn),'(a)')typ_rxnlbl(nrxn)(:i1)//aname(:1)
       else
        write(typ_rxnlbl(nrxn),'(a)')aname(:1)
       endif
       call movlft2(typ_rxnlbl(nrxn))
c       print*,'new srxn lbl:**'//typ_rxnlbl(nrxn)//'**'
c
c set the new reaction type and check to see if irxn is a samek
c
       rxtyp(nrxn) = 0
       if(rxtyp(irxn).eq.0)then
        lkbuf(nrxn) = lkbuf(irxn)
       else
        lkbuf(nrxn) = irxn
       endif
c
c set the number of reactants
c
       nrtosr(nrxn) = nrtosr(irxn)
c
c set the number of products.  this is tricky because we allow integer
c coefficients in the product list.
c
c note: in the first revision of typed soam we considered a special
c case where n2o5_N1_N2 -> 0.5 no2_N1 + 0.5 no2_N2 + 0.5 no3_N1 + 0.5 no3_N2
c but now we simplify to n2o5_N1_N2 -> no2_N1 + no3_N2 because this seems
c more generally correct for soa source apportionment (and a lot easier to
c program)
c
       icount = nrtosr(irxn)
       ip = 0
c       print*,'# base products:',nprods(irxn+1)-nprods(irxn)
       do j = nprods(irxn),nprods(irxn+1)-1
        icount = icount + 1
        if(iprods(j).gt.0 .and. tno(icount).gt.0)then
c count any leading coefficients and calculate their product
         c1 = 1.0
         do j2 = j-1,nprods(irxn),-1
          if(iprods(j2).gt.0)goto 50
          if(-iprods(j2).le.maxco)c1 = c1 * coef(-iprods(j2))
         enddo
 50      continue
         nlc = j - 1 - j2 ! number of leading coefficients
c if we have an int leading coefficient, then expand the products.  otherwise
c check to make sure we have a single typed reactant, and continue
         if(amod(c1,1.0).eq.0.0 .and. c1.gt.0)then
          ip = ip - nlc + int(c1) 
         elseif(nr.le.1)then
          ip = ip + 1
c note that we allow a leading coefficient of -1 to indicate that the
c product is actually a reactant.
         elseif(c1.eq.-1)then
          ip = ip + 1
         else
          if(.not. use_src_template)then
          print*,'somehow a fractional product coeficient exists for ',
     +       'a reaction with more than 1 typed reactant',rxnlbl(irxn)
          stop 'marker1'
          endif ! .not. use_src_template
         endif
        else
         ip = ip + 1
        endif
       enddo
       nprods(nrxn+1) = nprods(nrxn) + ip 
       if(nprods(nrxn+1).gt.maxprd)then
        write(icrt,'(a)')"too many products in expansions"
        stop 'too many products in typrxn'
       endif
c       print*,'number of products:',ip
c
c generate the suffix for this index location.  note: this is the 
c combined suffix for all reactants or all products.  we will need
c to take pieces of this suffix for each species.
c
       suffix = sufgen(i,int(trs),nt1,aname(:1))
       print*,'rxn suffix**'//suffix(:nblank(suffix))//'**'
       icount = 0
       outbuf = ' '
c
c set the reactant numbers
c
       i1 = 1
       do j = 1,nrtosr(nrxn)
        icount = icount + 1
        if(tno(icount).gt.0)then
c find the proper substring of the suffix; i1 is the start of the current
c substring, and i2 is the start of the next substring
         isum = 0
         do i2 = i1,len(suffix)
          if(suffix(i2:i2+1).eq.'_'//aname(:1) .or. 
     +       suffix(i2:i2).eq.' ')isum=isum+1
          if(isum.gt.tsno(icount))goto 20
         enddo
 20      continue 
c now check the reactant names
         irtosr(j,nrxn)=ifindspc(name(irtosr(j,irxn)),suffix(i1:i2-1))
         if(irtosr(j,nrxn).le.0)then
          print*,'couldnt find expected reactant',tr,tp,aname(:1)
          print*,'basename,suffix:**'//name(irtosr(j,irxn))//'**'//
     +            suffix(i1:i2-1)//'**',i1,i2-1
          print*,'reactant list:',rxnlbl(irxn),nrtosr(irxn)
          do j2 = 1,nrtosr(irxn)
           if(irtosr(j2,irxn).gt.0)print*,name(irtosr(j2,irxn))
          enddo
          stop 'marker 1'
         endif
         i1 = i2
        else
         irtosr(j,nrxn) = irtosr(j,irxn)
        endif
       enddo
c
c remember any remaining portion of the suffix - this indicates that
c we have products with a coefficient of -1
c
       i1react = i1
c       if(i1react.lt.len(suffix))
c     +         print*,'leftover suffix:**'//suffix(i1react:)
c
c set the product numbers and coeficients.  this is tricky because we
c allow integer coefficients in the product list. jbase will be the 
c index count in the new rxn, while j is the index count from the 
c original rxn.
c
       i1 = 1
       jbase = nprods(nrxn)
       do j = nprods(irxn),nprods(irxn+1)-1
        icount = icount + 1
        if(iprods(j).gt.0)then
          print*,'# tracked atoms in:',name(iprods(j)),tno(icount)
        endif
        if(tno(icount).le.0)then ! this product has no tracked atoms
         iprods(jbase) = iprods(j)
         jbase = jbase + 1
        else ! this product has tracked atoms 
c count any leading coefficients and calculate their product
         c1 = 1.0
         do j2 = j-1,nprods(irxn),-1
          if(iprods(j2).gt.0)goto 40
          if(-iprods(j2).le.maxco)c1 = c1 * coef(-iprods(j2))
         enddo
 40      continue
c if we have integer leading coeficients then expand them.  check again
c to make sure fractional coeficients only occur when we have 1 typed
c reactant.  np is the number of products produced by current base product.
         if(amod(c1,1.0).eq.0.0 .and. c1.gt.0)then
          nlc = j - 1 - j2 
          jbase = jbase - nlc
          np = int(c1)
         elseif(nr.eq.1)then
          np = 1
          i1 = 1
c note that we allow a leading coefficient of -1 to indicate that the
c product is actually a reactant.
         elseif(c1.eq.-1)then
          np = 1
         else
          if(.not. use_src_template)then
          print*,'somehow a fractional product coeficient exists for ',
     +       'a reaction with more than 1 typed reactant',rxnlbl(irxn)
          stop 'marker2'
          endif ! .not. use_src_template
          np=1
         endif
c loop over the products, expanding as necessary.  np is the number of
c products derived from the current product.  note the change to handle
c the special case of products with a coefficient of -1.
         do 60 j2 = 1,np
c find the proper substring of the suffix; i1 is the start of the current
c substring, and i2 is the start of the next substring.  note that we 
c finish adding portions of the reaction suffix if this product has a
c leading coefficient of -1.  the current product suffix will be stored.
          if(c1.eq.-1)then
           i1prod = i1
           i1 = i1react
          endif
c use indicated src tempalte pattern if this is a src template rxn
          if(use_src_template)then
           ii=iprods(j)
           call get_srcs(irxn,name(ii),a1,isrc_p,icnt_src,0)
           do ii=1,isrc_p(0)
            suffix3=sufext(suffix,isrc_p(ii))
            if(ii.eq.1)then
             suffix2=suffix3
            else
             jj=nblank(suffix2)
             write(suffix2(jj+1:),'((a))')suffix3
            endif
           enddo ! ii=1,isrc_p(0)
c normal automatic type assignments to products
          else
           print*,'j2,suffix:',j2,suffix
           print*,'tp,tr,aname:',tp,tr,aname(:1)
           isum = 0
           do i2 = i1,len(suffix)
            if(suffix(i2:i2+1).eq.'_'//aname(:1) .or. 
     +         suffix(i2:i2).eq.' ')isum=isum+1
            if(isum.gt.tsno(icount))goto 30
           enddo ! i2=i1,len(suffix)
 30        continue 
           write(suffix2,'(a)')suffix(i1:i2-1)
          endif
c now check the name
          i4=nblank(suffix2)
          iprods(jbase) = ifindspc(name(iprods(j)),suffix2(:i4))
          if(iprods(jbase).le.0)then
           i3=nblank(name(iprods(j)))
           print*,'could not find expected product '//
     +             name(iprods(j))(:i3)//suffix2(:i4)
           print*,'name:**'//name(iprods(j))//'**'
           print*,'suffix2: **'//suffix2//'**'
           stop 'could not find expected product'
          endif
          if(tp.gt.tr)then
           if(index(rxnlbl(irxn),'-'//aname(:1)).ne.0)then
            print*,'keep 1st source tag for product '//
     +             name(iprods(jbase))
           elseif(use_src_template)then
            print*,'rxn template used to map product '//
     +             name(iprods(jbase))
           else
            print*,'marker1 error in typeosubs.  program should not'
            print*,'have tp>tr.  tp=',tp,' tr=',tr
            stop 'typeosubs marker1'
           endif
          else
           i1 = i2
          endif
          jbase = jbase + 1
 60      continue
c remember the index if we just added a portion of the reaction suffix to
c this product because it has a leading coefficient of -1.  also resume
c using the product suffix.
         if(c1.eq.-1)then
          i1react = i1
          i1 = i1prod
         endif
        endif ! tno(icount).le.0
       enddo ! j = 1,nrtosr(nrxn)
c
c write the products to the output buffer
c
       call write_rxn(nrxn,.true.,.true.)

 10   continue ! bottom of the added reaction loop
c
c modify the base secondary reaction label to include the typed atom
c
      i1 = max(nblank(typ_rxnlbl(irxn)),1)
      if(i1.gt.0)then
       write(typ_rxnlbl(irxn),'(a)')typ_rxnlbl(irxn)(:i1)//aname(:1)
      else
       write(typ_rxnlbl(irxn),'(a)')aname(:1)
      endif
      call movlft2(typ_rxnlbl(irxn))
c
c return to the calling subroutine
c
      return
       end

      subroutine chktspc2(atom,tno,i1,k)
c********************************************************************
c written by: Mike Kleeman (Nov 2004)
c             UC Davis CEE
c
c The purpose of this subroutine is to check and see if the species
c i1 has been typed for "atom" tracking base species k.  It will be 
c added to the species array of it is missing.  
c
c The suffix added to "basename" will be _X?_X?_X?... where X is
c the traced molecule (N, S, C, X) and ? is some number 1<?<Xno.  The 
c number of iterations is set by tno
c
c Note that the include file pspecs.inc has common blocks that pass 
c most of the i/o.
c
c Inputs:
c  atom  - atom being tracked (N, S, or C)
c  i1    - index number for species of interest
c  k     - base tracked species
c
c Outputs:
c
c********************************************************************

      include 'pspecs.inc'
      character atom*1,suffix*30,sufgen*30,name2*16

c--determine the number of atoms in the target species--
      n = int(tno)

c--check the number of atoms for limiting cases--
      if(n.eq.0)return

c--get the length of the basename and the number of typed species--
      l1 = index(name(i1),'[_'//atom)-1
      if(l1.lt.1)l1 = nblank(name(i1))
      nt1 = ttypes(k) + 1
      nt = nt1**n

c--loop over all the index values and check to see if the species is
c  in the active list.  note: vanilla species 1 is skipped--
      do 10 j1 = 2,nt
       suffix = sufgen(j1,n,nt1,atom)
c hunt through the existing names for basename plus the target suffix
       i = ifindspc(name(i1),suffix)
       if(i.gt.0)goto 10
c if we reach this point, the typed species must be added
       do 20 i = 1,ns
c expand any species that matches the basename
        if(name(i)(:l1).eq.name(i1)(:l1) .and. 
     +     index(name(i),'_'//atom).eq.0 .and.
     +     (name(i)(l1+1:l1+1).eq.' ' .or. 
     +                name(i)(l1+1:l1+1).eq.'_') )then
c expand this species; note that spcnam will change ns, but this 
c doesn't affect the loop for i after it is initialized
         if(ns+1.gt.maxns) stop 'chktspc2 marker3: exceeded ns'
         write(name2,'(a)')name(i)(:nblank(name(i)))//
     +                     suffix(:nblank(suffix))
c         print*,'adding species 1**'//name2//'**'
         call spcnam(isp,name2,.true.)
         sptyp(isp) = sptyp(i)
         mwt(isp) = mwt(i)
         conc0(isp) = conc0(i)
         nno(isp) = nno(i)
         nsno(isp) = nsno(i)
         cno(isp) = cno(i)
         csno(isp) = csno(i)
         ono(isp) = ono(i)
         osno(isp) = osno(i)
         sno(isp) = sno(i)
         ssno(isp) = ssno(i)
         xno(isp) = xno(i)
         xsno(isp) = xsno(i)
         write(typ1, 1005) name(isp),sptyp(isp),mwt(isp),conc0(isp),
     +                     nno(isp),cno(isp),sno(isp),ono(isp)
 1005    format(2x,'Added species ',a16,i4,f8.3,f8.3,i4,f8.3,i4,f8.3)
        endif
 20    continue
 10   continue

c--return to the calling subroutine--
      return
      end
C
C -----------------------------------------------------------------------------
C
       subroutine TRDATOM(ISP,ANAME)
c
c       Read in the atom to be typed and check to see if it's okay
c         make sure that the atom is either N,C,S, or -
c
       character*16 aname
C
        INCLUDE 'pspecs.inc'
C
       call movlft(aname)
c       write(icrt,'(a16)') aname
       if(aname(2:2) .eq. ' ') then
         if(aname(1:1) .eq. 'n' .or. aname(1:1) .eq. 'N') then
           aname(1:1) = 'N'
         elseif(aname(1:1) .eq. 'c' .or. aname(1:1) .eq. 'C') then
           aname(1:1) = 'C'
         elseif(aname(1:1) .eq. 's' .or. aname(1:1) .eq. 'S') then
           aname(1:1) = 'S'
         elseif(aname(1:1) .eq. 'o' .or. aname(1:1) .eq. 'O') then
           aname(1:1) = 'O'
         elseif(aname(1:1) .eq. 'x' .or. aname(1:1) .eq. 'X') then
           aname(1:1) = 'X'
         else
c          none of the above and that's a problem what to do
           aname(1:1) = '*'
         endif
       else
c        is it written out as "Nitrogren" ?
         if(aname(1:8).eq.'nitrogen'.or.aname(1:8).eq.'Nitrogen') then
           aname(1:1) = 'N'
         elseif(aname(1:6).eq.'carbon'.or.aname(1:6).eq.'Carbon') then
           aname(1:1) = 'C'
         elseif(aname(1:6).eq.'sulfur'.or.aname(1:6).eq.'Sulfur') then
           aname(1:1) = 'S'
         elseif(aname(1:6).eq.'oxygen'.or.aname(1:6).eq.'Oxygen') then
           aname(1:1) = 'O'
         elseif(aname(1:4).eq.'none'.or.aname(1:4).eq.'None'
     &                      .or.aname(1:4).eq.'NONE') then
           aname(1:1) = 'X'
         else
c          none of the above and that's a problem what to do
           aname(1:1) = '*'
         endif
       endif
c
       tatom(ISP)=aname(1:1)
c
c        finished so return
c
       return
       end
c
C -----------------------------------------------------------------------------
        SUBROUTINE TSPCNAM(ISP,SPNAME,RX)
C
C       DETERMINE SPECIES NO. FOR NAME IN TSPNAME.  RETURNS SPECIES NO. IN
C       ISP.  ADDS SPECIES TO NAME to type ARRAY IF NEW, AND DO OTHER UPDATING.
C
C       CALLED FROM:  RDRXN, TYPEO
C
C
C       ARGUMENTS
C
        LOGICAL   RX
        integer i
        CHARACTER*16 SPNAME
C
C       SPECIFICATIONS FOR PREPARATION PROGRAM VARIABLES, PARAMETERS, AND
C       ARRAYS
C
        INCLUDE 'pspecs.inc'
C
C
       IF (NST.EQ.0) GOTO 1125
       DO 1122 ISP=1,NST
       IF (TNAME(ISP).NE.SPNAME) GOTO 1122
D      WRITE (OUT,*) 'SPCNAM: (OLD) ',TNAME(ISP),'  SPC NO.=',ISP
       IF (.NOT.RX) THEN
C if there needs to be any processing of an already specified typeo species
C put it here.
       ENDIF
       RETURN
 1122  CONTINUE
C  -     NAME NOT FOUND - NEW SPECIES ADDED TO LIST
 1125  NST=NST+1
       IF (NST.GT.MAXNS) THEN
                I=MAXNS
                IF (RX) WRITE (OUT,1127) OUTBUF
 1127           FORMAT (' ',A80)
                WRITE (OUT,*) 'TOO MANY SPECIES.  MAX =',I
                STOP 'TOO MANY SPECIES'
       ENDIF
       TNAME(NST)=SPNAME
       IF (IND.GT.0) THEN
c        write in the default atom and number of sources
       ELSE
                WRITE (OUT,*) 'PGM ER. IND<0 AT TSPCNAM'
       ENDIF
       ISP=NST
       RETURN
       END


      function nblank(buf)
c*****************************************************************************
c
c written by: Mike Kleeman (May 2002)
c             Dept of Civil and Env Eng
c             UC Davis 95616
c
c The purpose of this function is to return the rightmost position of non
c blank characters in the input buffer.
c
c*******************************************************************************

      integer    nblank, ilen
      character* (*) buf
      nblank = -1
      ilen = len(buf)
      do i = ilen,1,-1
       if(buf(i:i).ne.' ')then
        nblank = i
        return
       endif
      enddo
      return
      end

      function sufgen(j1,n,nt1,atom)
c***************************************************************************
c  written by: Mike Kleeman (Nov 2004)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this function is to generate the suffix used in the 
c naming procedure for typed species.
c
c Inputs:
c  j1   - current count in all the possibilities
c  n    - number of target atoms 
c  nt1  - number of sources to be typed
c  atom - N, S, C, or - (name of the traced atom)
c
c Outputs:
c  sufgen - character* string with the required suffix
c***************************************************************************
      character* (*) sufgen
      character  atom*1
      integer i2(10)

c--check the size--
      if(n.gt.10)then 
       print*,'must increase the dimension of working vector '//
     +        'i2 in sufgen'
       stop 'marker2 in sufgen'
      endif

c--remember the original count--
      j2 = j1

c--basically convert the value of j1 to base nt1 and record in i2--
      do j3 = n,1,-1
       i2(j3) = (j2-1) / nt1**(j3-1)
       j2 = j2 - i2(j3)*nt1**(j3-1)
      enddo
c      print*,j1,(i2(j3),j3=n,1,-1)

c--write the results to the output string--
      sufgen = ' '
      do j3 = n,1,-1
       L = max(nblank(sufgen)+1,1)
       if(i2(j3).lt.10)then
        write(sufgen(L:),10)atom,i2(j3)
 10     format('_',a1,i1)
       elseif(i2(j3).lt.100)then
        write(sufgen(L:),20)atom,i2(j3)
 20     format('_',a1,i2)
       elseif(i2(j3).lt.1000)then
        write(sufgen(L:),30)atom,i2(j3)
 30     format('_',a1,i3)
       else
        print*,'more than 1000 sources doesnt seem reasonable'
        stop 'marker1 in sufgen'
       endif
      enddo

c--return to the calling subroutine--
      return
      end

      function ifindspc(basename,suffix)
c*************************************************************************
c written by: Mike Kleeman (Nov 2004)
c             UC Davis CEE
c             Davis CA 95616
c
c The purpose of this function is to find the index value of the species
c name that contains basename and suffix.  It is assumed that each suffix
c has the form _X?_X?_X? where X is the traced atom, and ? is the source 
c number.  Note that all vanilla species will have suffix _X0
c
c Modified by: Mike Kleeman (June 2018).  Changes made to deal with the
c [_X name fragment.  If suffix has this value, then return the base
c species index.  If suffix doesn't have this value but basename does
c have this value, then substitute suffix in place of [_X* and find
c that species
c
c**************************************************************************
      include 'pspecs.inc'
      character* (*) basename, suffix
      character a1
      integer    ifindspc

c--set default value to not found
      ifindspc = -1
      a1=' '

c      print*,'ifindspc basename, suffix ',basename,suffix
c--check for the string [_X in the target suffix--
      if(suffix(:2).eq.'[_'.and.(suffix(3:3).eq.'N' .or.
     +   suffix(3:3).eq.'S' .or. suffix(3:3).eq.'O' .or.
     +   suffix(3:3).eq.'C' .or. suffix(3:3).eq.'X') )then
c sanity checks
       a1=suffix(3:3)
       ii=index(basename,'[_'//a1)-1
       if(ii.lt.1)then
        print*,'ifindspc called with suffix [_'//a1
        print*,'but basename does not contain that string'
        print*,'basename **'//basename//'**'
        stop 'ifindspc marker1'
       endif
c search for the basename species
       do i=1,ns
        jj=nblank(name(i))
c        print*,'ifindspc name ',name(i)
        if(name(i)(:jj).eq.basename(:ii))then
         ifindspc = i
c         print*,'ifindspc found'
         return
        endif ! trim(name(i)).eq.basename(ii-1)
       enddo ! i=1,ns 
       return
      endif ! suffix(:2).eq.'[_'

c--check the suffix to see if it is vanilla.  this is guaranteed if
c  the sum of all _X? values is zero--
      l2 = nblank(suffix)
      isum = 0
      do i = 1,l2
       if(suffix(i:i).eq.'_')then
        a1=suffix(i+1:i+1)
        read(suffix(i+2:i+2),'(i1)')j
        isum = isum + j
       endif
      enddo

c--check to see if the suffix is _X0 which signifies vanilla target for
c  this traced atom--
      if(isum.eq.0)then
c hunt through the species names for exact match on basename
c       print*,'hunting for name0:**'//basename//'**'
       l1 = index(basename,'[_'//a1)-1
       if(l1.le.0)l1 = nblank(basename)
       do i = 1,ns
        jj=nblank(name(i))
        if(name(i)(:jj).eq.basename(:l1))then 
         ifindspc = i
         return
        endif
       enddo
      else
c hunt through the species names for basename_suffix
       l1 = index(basename,'[_'//a1)-1
       if(l1.le.0)l1 = index(basename,'_')-1
       if(l1.le.0)l1 = nblank(basename)
       l2 = nblank(suffix)
c       print*,'hunting for name1:**'//basename(:l1)//suffix(:l2)//'**'
       do i = 1,ns
        l3 = index(name(i),'_')-1
        if(l3.le.0)l3 = nblank(name(i))
        if(name(i)(:l3).eq.basename(:l1) .and.
     +     index(name(i),suffix(:l2)).gt.0)then
         ifindspc = i
         return
        endif
       enddo
      endif

c--species not found, so return empty--
      return
      end

      SUBROUTINE MOVLFT2 (BUF)
C
C     MOVES CHARACTERS ORIGINALY PADDED TO LEFT WITH
C     BLANKS SO THAT BLANKS ARE PADDED TO THE RIGHT INSTEAD.
C
      CHARACTER* (*) BUF
C
      ILEN = len(buf)
      DO 10 N=1,ILEN
      IF (BUF(N:N).NE.' ') GOTO 20
   10 CONTINUE
      RETURN
C
   20 IF (N.EQ.1) RETURN
      N=N-1
      L=ILEN-N
      DO 30 I=1,L
   30 BUF(I:I)=BUF(I+N:I+N)
C
      N=L+1
      DO 40 I=N,ILEN
   40 BUF(I:I)=' '
      RETURN
      END

      subroutine write_rxn(irxn,wout,wunit3)
c
c     writes a reaction out to unit=out and unit3=3 (set is pspecs.inc)
c
      include 'pspecs.inc'
      character sname*21
      logical newrxn,wout,wunit3
c
c write the reactants to the output buffer
c
       outbuf=' '
       do j = 1,nrtosr(irxn)
        if(irtosr(j,irxn).gt.0)then
         sname = name(irtosr(j,irxn))
        elseif(-irtosr(j,irxn).gt.maxmax)then
         write(sname,'(a)')rxnlbl(-irtosr(j,irxn)-maxmax)
         call movlft(sname)
         write(sname,'(a)')'#RCON'//sname(:nblank(sname))
        elseif(-irtosr(j,irxn).ge.maxcov)then
         write(sname,'(f16.3)')coef(-irtosr(j,irxn))
         call movlft(sname)
         write(sname,'(a)')'#'//sname(:nblank(sname))
        elseif(-irtosr(j,irxn).ne.0)then
         write(sname,'(a)')coefnm(-irtosr(j,irxn))
         call movlft(sname)
         write(sname,'(a)')'#'//sname(:nblank(sname))
        else
         write(icrt,*)'coefnm(0) is an error!'
         write(icrt,*)outbuf(:nblank(outbuf))
         write(icrt,*)(irtosr(i1,irxn),i1=1,nrtosr(irxn))
         stop 'marker1'
        endif
        call movlft(sname)
        jj = nblank(outbuf) + 2
        write(outbuf(jj:),'(a)')sname(:nblank(sname))
c add a plus sign for real species
        if(irtosr(j,irxn).gt.0)then
         if(j.lt.nrtosr(irxn))then
          jj = nblank(outbuf) + 2
          write(outbuf(jj:),'(a)')'+'
         endif
        endif
       enddo
c add an equal sign 
       jj = nblank(outbuf) + 2
       write(outbuf(jj:),'(a)')'='
c
c write the products to the output buffer
c
       newrxn = .true.
       do j = nprods(irxn),nprods(irxn+1)-1
        if(iprods(j).gt.0)then
         sname = name(iprods(j))
        elseif(-iprods(j).gt.maxmax)then
         write(sname,'(a)')rxnlbl(-iprods(j)-maxmax)
         call movlft(sname)
         write(sname,'(a)')'#RCON'//sname(:nblank(sname))
        elseif(-iprods(j).ge.maxcov)then
         write(sname,'(f16.3)')coef(-iprods(j))
         call movlft(sname)
         write(sname,'(a)')'#'//sname(:nblank(sname))
        elseif(-iprods(j).ne.0)then
         write(sname,'(a)')coefnm(-iprods(j))
         call movlft(sname)
         write(sname,'(a)')'#'//sname(:nblank(sname))
        else
         write(icrt,*)'coefnm(0) is an error!'
         write(icrt,*)outbuf(:nblank(outbuf))
         write(icrt,*)(iprods(i1),i1=nprods(irxn),nprods(irxn+1)-1)
         stop 'marker2'
        endif
        call movlft(sname)
        jj = nblank(outbuf) + 2
        write(outbuf(jj:),'(a)')sname(:nblank(sname))
        jj = nblank(outbuf) + 2
        if(iprods(j).gt.0 .and. j.lt.nprods(irxn+1)-1)then
         write(outbuf(jj:),'(a)')'+'
        endif
c intermediate dump of the output buffer if necessary
        jj = nblank(outbuf) + 2
        if(jj.ge.62 .and. j.lt.nprods(irxn+1)-1)then
         write(outbuf(jj:),'(a)')'&'
         if(wunit3)write(unit3)newrxn,jj,outbuf(:jj)
         if(wout)then
          if(newrxn)then
           write(out,'(1x,a6,1x,(a))')rxnlbl(irxn),outbuf(:jj)
          else
           write(out,'(10x,(a))')outbuf(:jj)
          endif
         endif
         if(newrxn)then
          newrxn = .false.
         endif
         outbuf = ' '
        endif
       enddo
c
c final dump of the output buffer
c
       jj = nblank(outbuf)
       if(wunit3)write(unit3)newrxn,jj,outbuf(:jj)
       if(wout)then
        if(newrxn)then
         write(out,'(1x,a6,1x,(a))')rxnlbl(irxn),outbuf(:jj)
        else
         write(out,'(10x,(a))')outbuf(:jj)
        endif
       endif
       if(newrxn)then
        newrxn = .false.
       endif
c       print*,'outbuf:'
c       print*,'**'//outbuf(:nblank(outbuf))//'**'
c       stop 'debug point 1'
       return
       end

      subroutine check_src_template(irxn,aname,icnt_src,
     +                              use_src_template)
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to determine if this is a 
c template rxn used to explicitly map reactants to products for
c cases where automatic typing is not possible.
c Inputs:
c   irxn      - reaction number being checked
c   aname     - character string with the target typed molecule 
c
c Outputs:
c   icnt_src         - total number of sources tags in reactants
c   use_src_template - true if this is a template rxn
c***************************************************************************
      include 'pspecs.inc'
      integer ireact(9),isrc_r(0:9),isrc_p(0:9)
      character * (*) aname
      character a1
      logical use_src_template

c set the default to not using src template map
      use_src_template = .false.
      a1=aname(:1)

c scan the reactant list for [_X and find the
c corresponding vanilla species
      icnt_src=0
      do i = 1,nrtosr(irxn)
       ii = irtosr(i,irxn)
       if(ii.gt.0)then
         jj = index(name(ii),'[_'//a1)
         if(jj.gt.1)then
          if(index(rxnlbl(irxn),'+'//a1).le.0)then
           print*,'[_'//a1//' detected in reactant names but +'//a1
           print*,'not found in rxnlbl :',rxnlbl(irxn)
           print*,'change the reaction label to use source template'
           stop 'check_src_template marker0'
          endif
          if(.not.use_src_template)print*,'src template used for rxn '
     +                             ,rxnlbl(irxn)
          use_src_template=.true.
          call get_srcs(irxn,name(ii),a1,isrc_r,icnt_src,1)
          ib = ifindspc(name(ii),'[_'//a1)
c base spec found
          if(ib.gt.0)then
           print*,'template species, base species:',name(ii),name(ib)
           do kk=1,isrc_r(0)
            ireact(icnt_src-kk+1)=ib ! remember the tracked reactants
           enddo
c sync the N counts
           if(a1.eq.'N' .and. isrc_r(0).ne.nno(ib))then
            print*,'base species tracked N count =',nno(ib)
            print*,'template species tracked N count =',isrc_r(0)
            print*,'these must match.'
            stop 'check_src_template marker 0N'
           endif
           kk=max(nno(ii),nno(ib))
           nno(ii)=kk
           nno(ib)=kk
           kk=max(nsno(ii),nsno(ib))
           nsno(ii)=kk
           nsno(ib)=kk
c sync the S counts
           if(a1.eq.'S' .and. isrc_r(0).ne.sno(ib))then
            print*,'base species tracked S count =',sno(ib)
            print*,'template species tracked S count =',isrc_r(0)
            print*,'these must match.'
            stop 'check_src_template marker 0S'
           endif
           kk=max(sno(ii),sno(ib))
           sno(ii)=kk
           sno(ib)=kk
           kk=max(ssno(ii),ssno(ib))
           ssno(ii)=kk
           ssno(ib)=kk
c sync the C counts
           if(a1.eq.'C' .and. isrc_r(0).ne.cno(ib))then
            print*,'base species tracked C count =',cno(ib)
            print*,'template species tracked C count =',isrc_r(0)
            print*,'these must match.'
            stop 'check_src_template marker 0C'
           endif
           kk=max(cno(ii),cno(ib))
           cno(ii)=kk
           cno(ib)=kk
           kk=max(csno(ii),csno(ib))
           csno(ii)=kk
           csno(ib)=kk
c sync the O counts
           if(a1.eq.'O' .and. isrc_r(0).ne.ono(ib))then
            print*,'base species tracked O count =',ono(ib)
            print*,'template species tracked O count =',isrc_r(0)
            print*,'these must match.'
            stop 'check_src_template marker 0O'
           endif
           kk=max(ono(ii),ono(ib))
           ono(ii)=kk
           ono(ib)=kk
           kk=max(osno(ii),osno(ib))
           osno(ii)=kk
           osno(ib)=kk
c sync the X counts
           if(a1.eq.'X' .and. isrc_r(0).ne.xno(ib))then
            print*,'base species tracked X count =',xno(ib)
            print*,'template species tracked X count =',isrc_r(0)
            print*,'these must match.'
            stop 'check_src_template marker 0X'
           endif
           kk=max(xno(ii),xno(ib))
           xno(ii)=kk
           xno(ib)=kk
           kk=max(xsno(ii),xsno(ib))
           xsno(ii)=kk
           xsno(ib)=kk
          else
c force user to add the base spec
           print*,'please add a vanilla species to match'
           print*,'the source template species **'//name(ii)//'**'
           stop 'check_src_template marker'
          endif ! ib.gt.0
         endif ! index(name(ii),'[_'//a1).gt.0
       endif !irtosr(i,irxn).gt.0
      enddo ! i=1,nrtosr(irxn)
      print*,'finished scanning reactants. # srcs found: ',icnt_src

c scan the product list for [_X and find the
c corresponding vanilla species
      do i = nprods(irxn),nprods(irxn+1)-1
       ii=iprods(i)
       if(ii.gt.0)then
c         print*,'i,ii,name(ii):',i,ii,name(ii)
         jj = index(name(ii),'[_'//a1)
         if(jj.gt.1)then
          if(index(rxnlbl(irxn),'+'//a1).le.0)then
           print*,'[_'//a1//' detected in product names but +'//a1
           print*,'not found in rxnlbl :',rxnlbl(irxn)
           print*,'change the reaction label to use source template'
           stop 'check_src_template marker0'
          endif
          call get_srcs(irxn,name(ii),a1,isrc_p,icnt_src,0)
          use_src_template = .true.
          ib = ifindspc(name(ii),'[_'//a1)
c base spec found
          if(ib.gt.0)then
           print*,'template prod,base prod:',name(ii),name(ib)
c sync the N counts
           if(a1.eq.'N' .and. isrc_p(0).ne.nno(ib))then
            print*,'base species tracked N count =',nno(ib)
            print*,'template species tracked N count =',isrc_p(0)
            print*,'these must match.'
            stop 'check_src_template marker 1N'
           endif
           kk=max(nno(ii),nno(ib))
           nno(ii)=kk
           nno(ib)=kk
           kk=max(nsno(ii),nsno(ib))
           nsno(ii)=kk
           nsno(ib)=kk
c sync the S counts
           if(a1.eq.'S' .and. isrc_p(0).ne.sno(ib))then
            print*,'base species tracked S count =',sno(ib)
            print*,'template species tracked S count =',isrc_p(0)
            print*,'these must match.'
            stop 'check_src_template marker 1S'
           endif
           kk=max(sno(ii),sno(ib))
           sno(ii)=kk
           sno(ib)=kk
           kk=max(ssno(ii),ssno(ib))
           ssno(ii)=kk
           ssno(ib)=kk
c sync the C counts
           if(a1.eq.'C' .and. isrc_p(0).ne.cno(ib))then
            print*,'base species tracked C count =',cno(ib)
            print*,'template species tracked C count =',isrc_p(0)
            print*,'these must match.'
            stop 'check_src_template marker 1C'
           endif
           kk=max(cno(ii),cno(ib))
           cno(ii)=kk
           cno(ib)=kk
           kk=max(csno(ii),csno(ib))
           csno(ii)=kk
           csno(ib)=kk
c sync the O counts
           if(a1.eq.'O' .and. isrc_p(0).ne.ono(ib))then
            print*,'base species tracked O count =',ono(ib)
            print*,'template species tracked O count =',isrc_p(0)
            print*,'these must match.'
            stop 'check_src_template marker 1O'
           endif
           kk=max(ono(ii),ono(ib))
           ono(ii)=kk
           ono(ib)=kk
           kk=max(osno(ii),osno(ib))
           osno(ii)=kk
           osno(ib)=kk
c sync the X counts
           if(a1.eq.'X' .and. isrc_p(0).ne.xno(ib))then
            print*,'base species tracked X count =',xno(ib)
            print*,'template species tracked X count =',isrc_p(0)
            print*,'these must match.'
            stop 'check_src_template marker 1X'
           endif
           kk=max(xno(ii),xno(ib))
           xno(ii)=kk
           xno(ib)=kk
           kk=max(xsno(ii),xsno(ib))
           xsno(ii)=kk
           xsno(ib)=kk
c print out reactants contributing to this product
           print*,'# of reactants contributing to product:',isrc_p(0)
           do kk=1,isrc_p(0)
            print*,kk,isrc_p(kk),name(ireact(isrc_p(kk)))
           enddo ! kk=1,isrc_p(0)
          else
c force user to add the base spec
           print*,'please add a vanilla species to match'
           print*,'the source template species **'//name(ii)//'**'
           stop 'check_src_template marker'
          endif ! ib.gt.0
         endif ! jj.gt.1
       endif ! iprods(i).gt.0
      enddo ! i=nprods(irxn),nprods(irxn+1)-1

c end of the subroutine check_src_template
      return
      end

      subroutine get_srcs(irxn,sname,a1,isrc,icnt_src,icheck)
c********************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c inputs:  sname    - R[_X1_X2...]
c          a1       - one of N S C O X tags
c          icheck   - 0 dont check order
c                     1 check order
c                   
c outputs: isrc     -vector of sources
c                    isrc(0) total number of srcs
c                    isrc(1) first src tracked
c                    isrc(2) second src tracked, etc
c          icnt_src - global source count from all names
c*********************************************************
      include 'pspecs.inc'
      character * (*) sname
      character a1
      integer isrc(0:9)
c set all counters to zero
      isrc(0)=0                ! number of sources tracked in sname
      ib=index(sname,'[_'//a1) ! base index in string sname
      if(ib.le.0)return        ! no sources in this species
      l2=nblank(sname)
 10   if(ib.ge.l2)return
c      print*,'reading **'//sname(ib:)//'**'
      jj = index(sname(ib:),'_'//a1)
      if(jj.le.0)return
c read the src number from this location in the species
      read(sname(ib-1+jj+2:ib-1+jj+2),'(i1)')kk
c check the sequence of source tags in reactants
      if(icheck.ne.0)then 
       if(kk.ne.icnt_src+1)then
        print*,'error detected in reaction label ',rxnlbl(irxn)
        print*,'_'//a1//' detected in reactant names but number'
        print*,'must increase monotonically in reactants from'
        print*,'left to right.  found reactant '//sname
        print'((a),i1,a)','expected reactant '//sname(:jj+2)
     +                 ,icnt_src+1,']'
        stop 'get_srcs marker0a'
       endif ! kk.ne.icnt_src+1
       icnt_src=kk
c or check that the product sources match to the reactant sources
      else 
       if(kk.gt.icnt_src)then
        print*,'reaction template product mismatch for lbl',
     +          rxnlbl(irxn)
        print*,'product ',sname,'doesnt have matching reactant'
        print*,'current product source tag number:',kk
        print*,'maximum reactant source tag count:',icnt_src
        stop 'get_srcs marker0b'
       endif ! isrc.gt.icnt_src
      endif ! kk.ne.icnt_src+1
c increase the number of tracked sources
      i=isrc(0)+1
      isrc(0)=i
      if(i.gt.9)then
       print*,'expand vector isrc in get_srcs ',sname
       stop 'subroutine get_srcs marker1'
      endif ! i.gt.9
c remember the sources that are tracked
      isrc(i)=kk
      ib=ib+jj+2 ! change if not using '_SRC'
      goto 10
      end 

      function sufext(suffix,i)
c***************************************************************************
c  written by: Mike Kleeman (Nov 2004)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this function is to generate the ith fragment of the suffix
c
c Inputs:
c  suffix - string like _X1_X3_X2 or similar
c  i1     - target fragment of the suffix delimited by '_'
c
c Outputs:
c  sufext - character string with the required fragment of the suffix
c***************************************************************************
      character*30 sufext,suffix
      sufext=' '
      icount=0  ! counter in the extracted sufext
      jfrag=0   ! counter for fragments in the original suffix
c      print*,' '
c      print*,'sufext fragments'
      do j=1,nblank(suffix)
       if(suffix(j:j).eq.'_')jfrag=jfrag+1
       if(jfrag.eq.i)then
        icount=icount+1
        if(icount.gt.29)stop 'sufext suffix too long'
        write(sufext(icount:icount),'(a)')suffix(j:j)
c        print*,icount,sufext(icount:icount),suffix(j:j)
       endif ! jfrag.eq.i
      enddo ! j=1,nblank(suffix)
      icount=icount+1
      write(sufext(icount:icount),'(a)')' '
      if(icount.eq.0)then
       print*,'warning - could not find fragment ',i
       print*,'in suffix **'//suffix//'**'
       stop 'function sufext marker 1'
      endif ! icount.eq.0
      return
      end

      subroutine check_for_templates
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to search for template rxns.  These
c should either be deleted or they should be converted to vanilla rxns
c if the vanilla version is not present.
c***************************************************************************
      include 'pspecs.inc'
      character a1
      integer itrack(maxrxn)
      do irxn=1,nrxn
       itrack(irxn)=0
c scan the reactants looking for template fragments
       do i = 1,nrtosr(irxn)
        ii = irtosr(i,irxn)
        if(ii.gt.0)then
          j = index(name(ii),'[_')
          if(j.gt.1)then
           read(name(ii)(j+2:j+2),'(a)') a1
           if(index(rxnlbl(irxn),'+'//a1).gt.0)then
c we have identified that this is a template reaction
            print*,rxnlbl(irxn),'is is a template reaction'
            call check_for_vanilla(irxn,irxn2)
            if(irxn2.gt.0 .and. irxn2.le.nrxn)then
c we have found a matching vanilla reaction
             print*,'vanilla reaction found :',rxnlbl(irxn2)
             if(index(rxnlbl(irxn2),'NT').lt.1)then
              print*,'Error - the vanilla reaction cant use'
              print*,'automatic typing.  Please add NT to the label'
              stop 'check_for_templates marker1'
             endif
             print*,'template rxn ',rxnlbl(irxn),' will be deleted'
             itrack(irxn)=1
c redirect samek rate constants to vanilla rxn
             call redirect_samek(irxn,irxn2)
            else
c we will convert this to a vanilla reaction
             print*,'no vanilla reaction found'
             print*,'rxn ',rxnlbl(irxn),' will be converted to vanilla'
             itrack(irxn)=2
            endif
            goto 10
           endif ! 
          endif ! j.gt.1 .and. ...
        endif ! ii.gt.0
       enddo ! i = 1,nrtosr(irxn)
 10    continue
      enddo ! irxn=1,nrxn

c make changes to the reaction mechanisms 
      call modify_rxn(itrack)

      end

      subroutine redirect_samek(irxn1,irxn2)
c**************************************************************************
c  written by: Mike Kleeman (July 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to search for samek reactions aimed
c at irxn1 and redirect them to irxn2.
c
c Inputs:
c  irxn1  - original samek target
c  irxn2  - new samek target
c***************************************************************************
C
      include 'pspecs.inc'
      
      do irxn=1,nrxn
       if(rxtyp(irxn).eq.0 .and. lkbuf(irxn).eq.irxn1)then
         lkbuf(irxn)=irxn2
       endif
      enddo ! irxn=1,nrxn

c--return to the calling program--
      return
      end

      subroutine check_for_vanilla(irxn1,irxn2)
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to search for vanilla rxns matching
c a src template reaction. Vanilla reactions have vanilla versions of 
c the reactants and products.
c
c Inputs:
c  irxn1  - original reaction
c
c Outputs:
c  irxn2  - duplicate reaction
c***************************************************************************
C
      include 'pspecs.inc'

      character a1
      integer ir1(100),ir2(100)  ! reactant index
      integer ip1(100),ip2(100)  ! product index
      real    cr1(100),cr2(100)  ! reactant leading coefficients
      real    cp1(100),cp2(100)  ! product leading coefficients

      irxn2=-1

c assemble the "fingerprint" of the original reaction
      call get_fingerprint(irxn1,nr1,ir1,cr1,np1,ip1,cp1,1)

c loop over all reactions to search for matches to the original reaction
      do 10 irxn=1,nrxn
       if(irxn.eq.irxn1)goto 10
       call get_fingerprint(irxn,nr2,ir2,cr2,np2,ip2,cp2,0)
c check reactants
       if(nr1.ne.nr2)goto 10
       do 20 i=1,nr1
        do j=1,nr2
         if(ir1(i).eq.ir2(j) .and. abs(cr1(i)-cr2(j)).lt.0.001)goto 20 
        enddo ! j=1,nr2
        goto 10 ! does not match
 20    continue
c check products
       if(np1.ne.np2)goto 10
       do 30 i=1,np1
        do j=1,np2
         if(ip1(i).eq.ip2(j) .and. abs(cp1(i)-cp2(j)).lt.0.001)goto 30 
        enddo ! j=1,np2
        goto 10 ! does not match
 30    continue
c if we make it here, the reactions match.  set irxn2 and return
       irxn2=irxn
c       print*,'matching reaction fingerprints',irxn1,irxn2
c       print*,'cr1,ir1:',(cr1(i),ir1(i),i=1,nr1)
c       print*,'cr2,ir2:',(cr2(i),ir2(i),i=1,nr2)
c       print*,'cp1,ip1:',(cp1(i),ip1(i),i=1,np1)
c       print*,'cp2,ip2:',(cp2(i),ip2(i),i=1,np2)
       return
 10   continue

c return to the calling program with no match found
      return
      end

      subroutine get_fingerprint(irxn,nr,ir,cr,np,ip,cp,iflag)
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to return the "fingerprint" of a
c rxn composed of the reactant index values, reactant leading
c coefficients, product index values, and product leading coefficients.
c
c Inputs:
c  irxn   - reaction index number
c  iflag  - 0=do not convert all template species to vanilla
c             otherwise convert
c
c Outputs:
c  nr     - number of reactants
c  ir     - reactant index values
c  cr     - reactant leading coefficients
c  np     - number of products
c  ip     - product index values
c  cp     - product leading coefficients
c***************************************************************************
      include 'pspecs.inc'
      integer ir(100),ip(100)
      real    cr(100),cp(100)
      character a1

c initialize
      nr=0
      np=0
      do i=1,100
       ir(i)=0
       cr(i)=1.0
       ip(i)=0
       cp(i)=1.0
      enddo

c scan the reactants 
      c1=1.0
      do i = 1,nrtosr(irxn)
       i1 = irtosr(i,irxn)  ! index for reactant
       if(i1.lt.0)then      ! negative values are constants
        c1=c1*coef(-i1)
       else                 ! positive values are active species
c check to see if this species has a template fragment in the name
        if(iflag.ne.0)then
         j = index(name(i1),'[_')
         if(j.gt.1)then
          read(name(i1)(j+2:j+2),'(a)') a1
          i1 = ifindspc(name(i1),'[_'//a1)
          if(i1.lt.1)then
           print*,'problem analyzing template rxn ',rxnlbl(irxn1)
           print*,'could not find vanilla version of reactant ',
     +             name(i1)
           stop 'check_for_vanilla marker1'
          endif ! i1.lt.1
         endif ! j.gt.1
        endif ! iflag.ne.0
c check to see if this reactant is already in the list
        do j=1,nr
         if(ir(j).eq.i1)then
          cr(j)=cr(j)+c1
          c1=1.0
          goto 10
         endif ! ir(j).eq.i1
        enddo ! j=1,nr
c add the reactant to the list if we make it to this point
        nr=nr+1
        ir(nr)=i1
        cr(nr)=c1
        c1=1.0
 10     continue
       endif ! i1.lt.0
      enddo ! i = 1,nrtosr(irxn)

c scan the products
      c1=1.0
      do i = nprods(irxn),nprods(irxn+1)-1
       i1 = iprods(i)
       if(i1.lt.0)then      ! negative values are constants
        c1=c1*coef(-i1)
       else                 ! positive values are active species
c check to see if this species has a template fragment in the name
        if(iflag.ne.0)then
         j = index(name(i1),'[_')
         if(j.gt.1)then
          read(name(i1)(j+2:j+2),'(a)') a1
          i1 = ifindspc(name(i1),'[_'//a1)
          if(i1.lt.1)then
           print*,'problem analyzing template rxn ',rxnlbl(irxn1)
           print*,'could not find vanilla version of product ',
     +             name(i1)
           stop 'check_for_vanilla marker2'
          endif ! i1.lt.1
         endif ! j.gt.1
        endif ! iflag.ne.0
c check to see if this reactant is already in the list
        do j=1,np
         if(ip(j).eq.i1)then
          cp(j)=cp(j)+c1
          c1=1.0
          goto 20
         endif ! ip(j).eq.i1
        enddo ! j=1,np
c add the reactant to the list if we make it to this point
        np=np+1
        ip(np)=i1
        cp(np)=c1
        c1=1.0
 20     continue
       endif ! i1.lt.0
      enddo ! i = nprods(irxn),nprods(irxn+1)-1

c return to the calling program
      return
      end

      subroutine modify_rxn(itrack)
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to modify rxns based on the input
c flag itrack.  Reactions can be deleted or converted from reaction
c templates to vanilla versions.
c
c Inputs:
c  iflag   = 0 no change
c            1 delete rxn
c            2 conver rxn to vanilla
c
c**************************************************************************
      include 'pspecs.inc'
      byte rxtyp2(maxrxn)
      logical rxopn
      integer itrack(maxrxn)
      INTEGER NRTOSR2(MAXRXN),NPRODS2(MAXRXN+1),LKBUF2(MAXRXN)
      INTEGER IRTOSR2(MAXRCT,MAXRXN),IPRODS2(MAXPRD)
      REAL*4  KPBUF2(MAXKBF)
      CHARACTER*1 a1
      CHARACTER*6 RXNLBL2(MAXRXN),typ_rxnlbl2(maxrxn)
 
c loop over all the reactions and modify as needed
      nrxn2=0
      LOCKBF=0
      do 1000 irxn=1,nrxn
c delete this rxn by not copying to working variables
       if(itrack(irxn).eq.1)then
        print*,'delete lbl '//rxnlbl(irxn),'irxn',irxn
        goto 1000
       endif
c if we are here, then this reaction will be used in some form
c copy the reaction parameters to the working variable
       nrxn2=nrxn2+1
       call tparm(irxn,rxnlbl,typ_rxnlbl,nrtosr,irtosr,
     +            nprods,iprods,rxtyp,lkbuf,kpbuf,
     +            nrxn2,rxnlbl2,typ_rxnlbl2,nrtosr2,irtosr2,
     +            nprods2,iprods2,rxtyp2,lkbuf2,kpbuf2,
     +            maxrxn,maxrct,maxprd,maxkbf,lockbf)
c convert this rxn to vanilla
       if(itrack(irxn).eq.2)then
        print*,'convert to vanilla lbl '//rxnlbl(irxn),'irxn',irxn
c        print*,'reactants:'
        do i=1,nrtosr2(nrxn2)
         ii=irtosr2(i,nrxn2)
         if(ii.gt.0)then
          j=index(name(ii),'[_')
          if(j.gt.1)then
           a1=name(ii)(j+2:j+2)
           ib = ifindspc(name(ii),'[_'//a1)
c           print*,name(ii),name(ib)
           irtosr2(i,nrxn2)=ib
          else
c           print*,name(ii)
          endif ! j.gt.1
         endif ! ii.gt.0
        enddo ! i=1,nrtosr2(nrxn2)
c        print*,'products:'
        do i=nprods2(nrxn2),nprods2(nrxn2+1)-1
         ii=iprods2(i)
         if(ii.gt.0)then
          j=index(name(ii),'[_')
          if(j.gt.1)then
           a1=name(ii)(j+2:j+2)
           ib = ifindspc(name(ii),'[_X')
c           print*,name(ii),name(ib)
           iprods2(i)=ib
          else
c           print*,name(ii)
          endif ! j.gt.1
         endif ! ii.gt.0
        enddo ! i=nprods2(nrxn2),nprods2(nrxn2+1)-1
       endif ! itrack(irxn).eq.2

 1000 continue ! irxn=1,nrxn

c loop over all reactions and copy the modified parameters back
c to original variables
      lockbf=0
      do irxn=1,nrxn2
       call tparm(irxn,rxnlbl2,typ_rxnlbl2,nrtosr2,irtosr2,
     +            nprods2,iprods2,rxtyp2,lkbuf2,kpbuf2,
     +            irxn,rxnlbl,typ_rxnlbl,nrtosr,irtosr,
     +            nprods,iprods,rxtyp,lkbuf,kpbuf,
     +            maxrxn,maxrct,maxprd,maxkbf,lockbf)
      enddo ! i=1,nrxn2
      nrxn=nrxn2

c done processing all reactions. 
c rewrite the temp file used for the construction of the .mod file
      inquire(unit=unit3,opened=rxopn)
      IF (RXOPN) CLOSE(UNIT3)
      IOB160=' '
      CALL FILNAM (NLEN,IOB160,TMPUIC,'PREP.TMP ',' ')
      OPEN (UNIT=UNIT3,NAME=IOB160,FORM='UNFORMATTED')
      RXOPN=.TRUE.
      DO IRXN=1,NRXN
       call write_rxn(irxn,.false.,.true.)
      ENDDO
      CLOSE (UNIT=UNIT3,STATUS='KEEP')
      RXOPN=.FALSE.

c return to the calling program
      return
      end

      subroutine tparm(irxn,rxnlbl,typ_rxnlbl,nrtosr,irtosr,
     +            nprods,iprods,rxtyp,lkbuf,kpbuf,
     +            irxn2,rxnlbl2,typ_rxnlbl2,nrtosr2,irtosr2,
     +            nprods2,iprods2,rxtyp2,lkbuf2,kpbuf2,
     +            maxrxn,maxrct,maxprd,maxkbf,lockbf)
c***************************************************************************
c  written by: Mike Kleeman (June 2018)
c              UC Davis CEE
c              Davis CA 95616
c
c The purpose of this subroutine is to transfer the parameters from 
c one irxn to the working variables at location irxn2.
c***************************************************************************
      byte rxtyp(maxrxn)
      INTEGER NRTOSR(MAXRXN),NPRODS(MAXRXN+1),LKBUF(MAXRXN)
      INTEGER IRTOSR(MAXRCT,MAXRXN),IPRODS(MAXPRD)
      REAL*4  KPBUF(MAXKBF)
      CHARACTER*6 RXNLBL(MAXRXN),typ_rxnlbl(maxrxn)
      byte rxtyp2(maxrxn)
      INTEGER NRTOSR2(MAXRXN),NPRODS2(MAXRXN+1),LKBUF2(MAXRXN)
      INTEGER IRTOSR2(MAXRCT,MAXRXN),IPRODS2(MAXPRD)
      REAL*4  KPBUF2(MAXKBF)
      CHARACTER*6 RXNLBL2(MAXRXN),typ_rxnlbl2(maxrxn)

      rxnlbl2(irxn2)=rxnlbl(irxn)
      typ_rxnlbl2(irxn2)=typ_rxnlbl(irxn)
      nrtosr2(irxn2)=nrtosr(irxn)
      do i=1,nrtosr(irxn)
       irtosr2(i,irxn2)=irtosr(i,irxn)
      enddo
      if(irxn2.eq.1)nprods2(irxn2)=1
      i2=nprods2(irxn2)
      do j = nprods(irxn),nprods(irxn+1)-1
       j2=j-nprods(irxn)
       iprods2(i2+j2) = iprods(j)
      enddo
      np=nprods(irxn+1)-nprods(irxn)
      nprods2(irxn2+1)=nprods2(irxn2)+np
      rxtyp2(irxn2)=rxtyp(irxn)
      if(rxtyp(irxn).eq.0)then
       lkbuf2(irxn2)=lkbuf(irxn)
      elseif(rxtyp(irxn).eq.3)then
       LOCKBF=LOCKBF+1
       LKBUF2(IRXN2)=LOCKBF
       KPBUF2(LOCKBF)=KPBUF(LKBUF(IRXN))
      elseif(rxtyp(irxn).eq.4)then
       LOCKBF=LOCKBF+1
       LKBUF2(IRXN2)=LOCKBF
       KPBUF2(LOCKBF)=KPBUF(LKBUF(IRXN))
       KPBUF2(LOCKBF+1)=KPBUF(LKBUF(IRXN)+1)
       KPBUF2(LOCKBF+2)=KPBUF(LKBUF(IRXN)+2)
       KPBUF2(LOCKBF+3)=KPBUF(LKBUF(IRXN)+3)
       LOCKBF=LOCKBF+3
      elseif(rxtyp(irxn).eq.5)then
       LOCKBF=LOCKBF+1
       LKBUF2(IRXN2)=LOCKBF
       KPBUF2(LOCKBF)=KPBUF(LKBUF(IRXN))
       KPBUF2(LOCKBF+1)=KPBUF(LKBUF(IRXN)+1)
       KPBUF2(LOCKBF+2)=KPBUF(LKBUF(IRXN)+2)
       KPBUF2(LOCKBF+3)=KPBUF(LKBUF(IRXN)+3)
       KPBUF2(LOCKBF+4)=KPBUF(LKBUF(IRXN)+4)
       KPBUF2(LOCKBF+5)=KPBUF(LKBUF(IRXN)+5)
       KPBUF2(LOCKBF+6)=KPBUF(LKBUF(IRXN)+6)
       KPBUF2(LOCKBF+7)=KPBUF(LKBUF(IRXN)+7)
       LOCKBF=LOCKBF+7
      elseif(rxtyp(irxn).eq.7)then
       LKBUF2(IRXN2)=LKBUF(IRXN)
      elseif(rxtyp(irxn).eq.9)then
       LOCKBF=LOCKBF+1
       LKBUF2(IRXN2)=LOCKBF
       KPBUF2(LOCKBF)=KPBUF(LKBUF(IRXN))
       KPBUF2(LOCKBF+1)=KPBUF(LKBUF(IRXN)+1)
       KPBUF2(LOCKBF+2)=KPBUF(LKBUF(IRXN)+2)
       KPBUF2(LOCKBF+3)=KPBUF(LKBUF(IRXN)+3)
       KPBUF2(LOCKBF+4)=KPBUF(LKBUF(IRXN)+4)
       KPBUF2(LOCKBF+5)=KPBUF(LKBUF(IRXN)+5)
       KPBUF2(LOCKBF+6)=KPBUF(LKBUF(IRXN)+6)
       LOCKBF=LOCKBF+6
      elseif(rxtyp(irxn).eq.10)then
       LOCKBF=LOCKBF+1
       LKBUF2(IRXN2)=LOCKBF
       KPBUF2(LOCKBF)=KPBUF(LKBUF(IRXN))
       KPBUF2(LOCKBF+1)=KPBUF(LKBUF(IRXN)+1)
       KPBUF2(LOCKBF+2)=KPBUF(LKBUF(IRXN)+2)
       KPBUF2(LOCKBF+3)=KPBUF(LKBUF(IRXN)+3)
       KPBUF2(LOCKBF+4)=KPBUF(LKBUF(IRXN)+4)
       KPBUF2(LOCKBF+5)=KPBUF(LKBUF(IRXN)+5)
       KPBUF2(LOCKBF+6)=KPBUF(LKBUF(IRXN)+6)
       KPBUF2(LOCKBF+7)=KPBUF(LKBUF(IRXN)+7)
       KPBUF2(LOCKBF+8)=KPBUF(LKBUF(IRXN)+8)
       KPBUF2(LOCKBF+9)=KPBUF(LKBUF(IRXN)+9)
       LOCKBF=LOCKBF+9
      else
       print*,'unknown rxtyp: ',rxtyp(irxn)
       stop 'modify_rxn marker1'
      endif ! rxtyp(irxn).eq. 0, ...

c return to the calling program
      return
      end

C end of typeosubs.for
