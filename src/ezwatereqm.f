
C     **************************************************
C     *  ezwatereqm                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, March 2000

C     This routine uses the current RH to calculate how much water is 
C     in equilibrium with the aerosol.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE ezwatereqm(Mke,ORG_O2C2,kappa_avg)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Mke(ibins,icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k,j
      double precision so4mass, naclmass, orgmass
      double precision wrso4, wrnacl, wrorg
      double precision rhe
      DOUBLE PRECISION ORG_O2C2(iorg)
      DOUBLE PRECISION kappa_avg(ibins),kappa
      double precision waterso4, waternacl, waterorg
      DOUBLE PRECISION tot_orgmass
      external waterso4, waternacl, waterorg

C     VARIABLE COMMENTS...

C     This version of the routine works for sulfate and sea salt
C     particles.  They are assumed to be externally mixed and their
C     associated water is added up to get total aerosol water.
C     wr is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fits to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K and sea salt at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      rhe=100.d0*rh
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0
C      print*,'ORG_O2C in ezwater =',ORG_O2C2
C      print*,'SUM(Mke(:,srth2o) before=',SUM(Mke(:,srth2o))
      do k=1,ibins
         kappa_avg(k) = 0.d0
         so4mass=Mke(k,srtso4)*1.2  !1.2 converts kg so4 to kg nh4hso4
         naclmass=0.d0
         orgmass=0.d0
         Mke(k,srth2o) = 0.d0
         tot_orgmass = SUM(Mk(k,srtorg1:srtorglast))
         do j=1,iorg
           orgmass=Mke(k,srtorg1+j-1)
           wrorg = waterorg(rhe,j,kappa,ORG_O2C2(j))
           Mke(k,srth2o) = Mke(k,srth2o) + 
     &       orgmass*(wrorg-1.d0)
           kappa_avg(k) = kappa_avg(k)+kappa*(orgmass/tot_orgmass)
         enddo
C         print*,'Bin =',k
C         print*,'Mass weighted kappa =',kappa_avg
         wrso4 = waterso4(rhe,k)
         wrnacl= 0.d0
         ! assume organics have same water uptake as ammonium bisulfate
         Mke(k,srth2o)= Mke(k,srth2o)+so4mass*(wrso4-1.d0)
     &                  +naclmass*(wrnacl-1.d0)

      enddo
C      print*,'SUM(Mke(:,srth2o) after=',SUM(Mke(:,srth2o))

      RETURN
      END

