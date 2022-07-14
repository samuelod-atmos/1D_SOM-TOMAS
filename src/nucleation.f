
C     **************************************************
C     *  nucleation                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.
C     The number of particles added to the first size bin is calculated
C     by comparing the growth rate of the new particles to the coagulation
C     sink.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     dt - total model time step to be taken (s)
C     sfctemp - surface temp for hacking nucleation dependent
C               temperature 
C     fion - ion recombination coefficient 
C     cstar_nuc - volatility threshold for organics to participate in
C                 nucleation.
C     org_nuc - switch for organic nucleation 
C     dunne_nuc - switch for inorganic nucleation 
C     nh3molec_cm31 - concentration of NH3 for inorganic nucleation
C     
C
C      
C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values

!      SUBROUTINE nucleation(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,nucscale)
      SUBROUTINE nucleation(Nki, Mki, Gci, Nkf, Mkf, Gcf,sfctemp,
     &     dt, fion, cstar_nuc, org_nuc,dunne_nuc,nh3molec_cm31,MW_ORG,
     &      HOM)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision dt
      double precision nucscale
      double precision nh3molec_cm31  ! Added by SamO
      DOUBLE PRECISION MW_ORG(iorg)    ! Added by SamO 
      DOUBLE PRECISION HOM

C-----VARIABLE DECLARATIONS---------------------------------------------

!      double precision nh3ppt   ! gas phase ammonia in pptv
!      double precision nh3moleccm3   ! gas phase ammonia in molec cc-1
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision org      ! gas phase organics in molec cc-1
      double precision fn       ! nucleation rate cm-3 s-1
      double precision rnuc     ! critical nucleation radius [nm]
      double precision gtime    ! time to grow to first size bin [s]
      double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      double precision Mktot    ! total mass in bin
      double precision neps
      double precision meps
      double precision density  ! density of particle [kg/m3]
      double precision pi
      double precision frac     ! fraction of particles growing into first size bin
      double precision d1,d2    ! diameters of particles [m]
      double precision mp       ! mass of particle [kg]
      double precision mold     ! saved mass in first bin
      double precision mnuc     !mass of nucleation
      double precision Mair     ! concentration of air [molec cm-3]
      double precision fntemp1  ! nucleation rate from Ricco 
      double precision fntemp2  ! nucleation rate from Dunne
      DOUBLE PRECISION TOG_LVOC
      DOUBLE PRECISION CSTAR_NUC
      INTEGER org_nuc, dunne_nuc
      DOUBLE PRECISION fion
      DOUBLE PRECISION INX      ! Counter 
      DOUBLE PRECISION MW_TEMP,sfctemp  ! Average molecular weight of organics
      DOUBLE PRECISION Jbn2, Jtn2, Jbi2, Jti2       
      DOUBLE PRECISION HOM2
      integer nuc_bin

      CHARACTER*120 FN_FILE

      COMMON /FLNAMES/ FN_FILE
C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens
      external aerodens

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159)

C-----CODE--------------------------------------------------------------

C     AliA ------------------------------------------------------------
      TOG_LVOC = 0.d0
      INX = 0.d0
      MW_TEMP = 0.d0
      HOM2 = HOM*1000.0/300.0*6.022E23/BOXVOL

      do k=1,icomp-1
        Gcf(k)=Gci(k)
      enddo
      
      IF (CSTAR_NUC.GT.0.0) THEN
        DO k=1,IORG
           IF (CSTAR(k).LT.CSTAR_NUC) THEN
              TOG_LVOC = TOG_LVOC + Gci(k+SRTORG1-1)
              INX = INX + 1.0
              MW_TEMP = MW_TEMP + MW_ORG(k)
           END IF
        END DO
        MW_TEMP = MW_TEMP/INX
        TOG_LVOC = TOG_LVOC/BOXVOL/MW_TEMP*1000*6.022e2
        PRINT*,'MW_TEMP =',MW_TEMP,'INX =',INX
      ENDIF
C      TOG_LVOC = TOG_LVOC*1.0E+9/BOXVOL*1.0E+6 ! kg/bag -> ug/m3                                                                                                                                             
C      MW_TEMP = MW_TEMP/INX
C      TOG_LVOC = TOG_LVOC/BOXVOL/MW_TEMP*1000*6.022e2
C      PRINT*, 'TOG_LVOC preHOM =', TOG_LVOC, '[molec/cm^3]'
      TOG_LVOC = TOG_LVOC + HOM2*0.1
C      PRINT*, 'TOG_LVOC =', TOG_LVOC, '[molec/cm^3]'
      PRINT*, 'HOM2 = ',HOM2
C      NC_RATE = A_NUC * (TOG_LVOC**B_NUC)
C      PRINT*, 'NC_RATE =', NC_RATE
C      Nk(1) = Nk(1) + NC_RATE*BOXVOL*ddt
C      PRINT*, 'Nk(1) =', Nk(:1)
C      Mki(1,SRTORGLAST)=Mki(1,SRTORGLAST)+NC_RATE*BOXVOL*ddt*
C     &                sqrt(xk(1)*xk(2))

C     AliA -------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23 ! [molec cm-3]
!      org = Gci(srtorg)/boxvol*1000.d0/200.d0*6.022d23*forgnuc ! Changed by AliA
      print*, 'h2so4 [molec/cm^3]=',h2so4,'org=',TOG_LVOC !org ! Changed by AliA
!      print*, h2so4
!      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d1 ! Changed by SamO
!      nh3moleccm3 = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23 ! Changed

C SamO changed this for temperature hijacking [temp or sfctemp]
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Mair = 2.69E19*273.15/temp*pres/101325.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      fn = 0.d0
      fntemp1 = 0.d0
      fntemp2 = 0.d0
      rnuc = 0.d0
      gtime = 0.d0

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (org_nuc .eq. 1) then

C Sam changed [temp or sfctemp]
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         call ricco_nucl(temp,h2so4,TOG_LVOC,fntemp1,rnuc) ! Changed by AliA
         print*, 'rnuc =', rnuc
         print*, 'fntemp1 =', fntemp1
         fn = fn + fntemp1
      endif

C Sam changed [temp or sfctemp]
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (dunne_nuc .eq. 1) then

        call dunne_inorg_nucl(temp,fion,h2so4,nh3molec_cm31,Mair,
     &                        fntemp2,rnuc,Jbn2,Jtn2,Jbi2,Jti2)
      print*,'Jbn2=', Jbn2,'Jtn2=', Jtn2,'Jbi2=', Jbi2,'Jti2=', Jti2
        print*,'fntemp2',fntemp2,'rnuc',rnuc 
        fn = fn + fntemp2
      endif
      print*,'fn = ',fn

C Write nucleation rates to FN_FILE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(UNIT=40,FILE=FN_FILE,STATUS='old',ACCESS='append')
      WRITE(40,*) fn, fntemp1, fntemp2, Jbn2, Jtn2, Jbi2, Jti2
      CLOSE(40)
      fntemp1 = 0.d0
      fntemp2 = 0.d0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!seo -------------------------------------------------------------- 
!seo      if (h2so4.gt.1.d4) then
!seo        if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
c              print*, 'napari'
!seo            call napa_nucl(temp,rh,h2so4,nh3ppt,fntemp,rnuc) !ternary nuc
!seo            fn = fn + fntemp
!seo            fntemp = 0.d0
!seo        elseif (bin_nuc.eq.1) then
!seo            print*, 'vehk'
!seo            call vehk_nucl(temp,rh,h2so4,fntemp,rnuc) !binary nuc
!seo            fn = fn + fntemp
!seo            fntemp = 0.d0
!seo        endif
!seo      endif
!seo       !      print*,fn
!seo      fn=fn*nucscale 

C     if nucleation occured, see how many particles grow to join the first size
C     section
      if (fn.gt.0.d0) then
c         print*, 'entered fn'
c         print*, 'fn', fn, 'rnuc',rnuc
         ! get the time it takes to grow to first size bin
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         d2 = (6.d0*sqrt(xk(1)*xk(2))/1800.d0/pi)**(1.d0/3.d0)
Cjrp         print*, 'd1', d1, 'd2', d2
Cjrp         call getGrowthTime(d1,d2,Gc(srtso4),temp,
Cjrp     &        boxvol,gtime)
Cjrp         print*, 'gtime', gtime
Cjrp         ! get the first order loss rate of these particles due to coagulation
Cjrp         ! get loss rate for cluster size
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         call getCoagLoss(d1,ltc1)
Cjrp         print*, 'd1', d1, 'ltc1', ltc1
Cjrp         ! we want the loss rate for particles at the size where they are added
Cjrp         ! to bin 1
Cjrp         Mktot = 0.d0
Cjrp         do j=1,icomp
Cjrp            Mktot=Mktot+Mk(1,j)
Cjrp         enddo
Cjrp         Mktot=Mktot+2.d0*Mk(1,srtso4)/96.d0-Mk(1,srtnh4)/18.d0 ! account for h+
Cjrp
Cjrp         if (Mktot.gt.meps)then
Cjrp            density=aerodens(Mk(1,srtso4),0.d0,Mk(1,srtnh4),
Cjrp     &           0.d0,Mk(1,srth2o)) !assume bisulfate
Cjrp         else
Cjrp            density = 1400.
Cjrp         endif
Cjrp
Cjrp         if(Nk(1).gt.neps .and. Mktot.gt.meps)then
Cjrp            mp=Mktot/Nk(1)
Cjrp         else
Cjrp            mp=sqrt(xk(1)*xk(2))
Cjrp         endif
Cjrp
Cjrp         d1 = (6.d0*mp/density/pi)**(1.d0/3.d0)
Cjrp         call getCoagLoss(d1,ltc2)
Cjrp         print*, 'd1', d1, 'ltc2', ltc2
Cjrp
Cjrp         ! take average of loss rate constants
Cjrp         ltc = (ltc1 + ltc2)/2.d0
Cjrp         print*, 'ltc', ltc
Cjrp
Cjrp	 if (gtime.lt.0.d0)then
Cjrp	    gtime=0.d0
Cjrp         endif
Cjrp
Cjrp         frac = exp(-ltc*gtime)
Cjrp         print*,'frac', frac
Cjrpc         frac = 1.d0

         mnuc = (4.d0/3.d0*pi*(rnuc*1D-9)**3)*1350.d0
         if (mnuc.lt.xk(1))then
            print*,'mnuc < xk(1) in nucleation',mnuc
            stop
         endif
         nuc_bin = 1
         do while (mnuc .gt. xk(nuc_bin+1))
            nuc_bin = nuc_bin + 1
         enddo
c         print*,'nuc_bin',nuc_bin

         mold = Mki(nuc_bin,srtso4)
         Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4)+fn*mnuc*
     &        boxvol*dt
         Nkf(nuc_bin) = Nki(nuc_bin)+fn*boxvol*dt
         Gcf(srtso4) = Gci(srtso4) - (Mkf(nuc_bin,srtso4)-mold)

C there is a chance that Gcf will go less than zero because we are artificially growing
C particles into the first size bin.  don't let it go less than zero.
         if (Gcf(srtso4).lt.0.d0)then
            Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4) +
     &           Gci(srtso4)*96./98.
            Nkf(nuc_bin) = Nki(1) + Gci(srtso4)*96./98./mnuc
            Gcf(srtso4) = 0.d0
            print*,'used up h2so4 in nuc subroutine'
         endif

      endif

      do k=1,ibins
         if (k .ne. nuc_bin)then
            Nkf(k) = Nki(k)
            do i=1,icomp
               Mkf(k,i) = Mki(k,i)
            enddo
         endif
      enddo
C      print*,'Finished nucleation.f'


      RETURN
      END
