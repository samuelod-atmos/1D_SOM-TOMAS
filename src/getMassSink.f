
C     **************************************************
C     *  getCondSink                                   *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine calculates the condensation sink (first order loss
C     rate of condensing gases) from the aerosol size distribution.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nk(ibins) - number of particles per size bin in grid cell
C     Nnuc - number of particles per size bin in grid cell
C     Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mk(ibins, icomp) - mass of a given species per size bin/grid cell
C     spec - number of the species we are finding the condensation sink for

C-----OUTPUTS-----------------------------------------------------------

C     CS - condensation sink [s^-1]
C     sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getMassSink(Nko,Mko,spec,CS,sinkfrac)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nko(ibins), Mko(ibins, icomp)
      double precision CS, sinkfrac(ibins)
      integer spec

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision mu                  !viscosity of air (kg/m s)
      double precision mfp                 !mean free path of air molecule (m)
      real Di       !diffusivity of gas in air (m2/s)
      double precision Neps     !tolerance for number
      real density  !density [kg m^-3]
      double precision mp       !mass per particle [kg]
      double precision Dpk(ibins) !diameter of particle [m]
      double precision Kn       !Knudson number
      double precision beta(ibins) !non-continuum correction factor
      double precision Mktot    !total mass in bin [kg]
      double precision tmass

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------
c      real aerodens, gasdiff
c      external aerodens, gasdiff

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(Neps=1.0d10)
!      double precision alpha(icomp) ! accomodation coef  
!      data alpha/1.0,1.0,1.0,1.0/
!      real Sv(icomp)         !parameter used for estimating diffusivity
!      data Sv /42.88,42.88,42.88,42.8/

C-----CODE--------------------------------------------------------------


C     get condensation sink
      CS = 0.d0
      do k=1,ibins
         do c=1,icomp ! partition to all mass now
            CS = CS + Mko(k,c)
         enddo
      enddo
      do k=1,ibins
         tmass = 0.d0
         do c=1,icomp
            tmass = tmass + Mko(k,c)
         enddo
         sinkfrac(k) = tmass/CS
      enddo
c      CS = 2.d0*pi*Di*CS/(boxvol*1D-6)      

      return
      end
