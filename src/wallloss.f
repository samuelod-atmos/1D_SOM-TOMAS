C =============================================================================================
C This subroutine has been written by Emily Bian at Atmospheric Science at CSU to simulate 
C both particle wall losses and vapor wall losses
C 
C The codes have been modified in Oct. 2018 by Ali Akherati at MECH-E at CSU     
C =============================================================================================

      subroutine wallloss(kei, kw0i, VWL, PWL,
     &                    Gci, Nki, Mki, dt, time, P_temp, P_num,
     &                    WalVapMas, tau_w, kw_HYBRID, aadt)

CCCCC DECLERATION...............................................................................
      IMPLICIT NONE

      INCLUDE 'sizecode.COM'

      DOUBLE PRECISION aadt   ! microphysics timestep
      DOUBLE PRECISION cdt    ! time step [seconds]
      DOUBLE PRECISION kB     ! Boltzmann constant [kg-m2/s2-K]
      DOUBLE PRECISION g      ! gravity [m/s2]
      DOUBLE PRECISION R      ! universal gas constant [kg-m3/s2-mol-K]
      DOUBLE PRECISION pi
      DOUBLE PRECISION lambda ! mean free path of particles in air [nm]
      DOUBLE PRECISION Na     ! Avogadro's number [molec/mol]
      DOUBLE PRECISION MWair  ! molecular weight of air, kg/mol
      INTEGER num
      
      PARAMETER(kB = 1.38e-23, g = 9.81, R = 8.314)
      PARAMETER(pi = 3.1415926, lambda = 65.2, Na = 6.0221409e23)
      PARAMETER(MWair = 0.02897, num = 10)
      
      INTEGER i, j, k, n
      INTEGER VWL, PWL
      
      DOUBLE PRECISION debye1
      EXTERNAL debye1

      DOUBLE PRECISION dt, time
      DOUBLE PRECISION kw0i, kei
      DOUBLE PRECISION WalVapMas(iorg)
      DOUBLE PRECISION tau_w(iorg)
      DOUBLE PRECISION rho, dil, rhoSI ! read in parameres **********
      DOUBLE PRECISION Mo, Dl, Dh, Dk(ibins), kw_turb(ibins)
      DOUBLE PRECISION Nki(ibins), Mki(ibins,icomp), Gci(icomp-1)
      DOUBLE PRECISION mu
      DOUBLE PRECISION Cc(ibins), B(ibins), D(ibins)
      DOUBLE PRECISION vT(ibins),x(ibins),D1(ibins)
      DOUBLE PRECISION V_bag, R_bag, StoV, S_bag, height, S1, S2
      DOUBLE PRECISION R_bag1,StoV1, S_bag1, factor ! for FSL volume test
      DOUBLE PRECISION gamma_p, alpha_w
      DOUBLE PRECISION Ra, morg(iorg), mair, dorg, dair
      DOUBLE PRECISION cmean(iorg), kw_gas(iorg), Dg(iorg)
      DOUBLE PRECISION wlc(iorg)
      DOUBLE PRECISION wlc_zhang
      DOUBLE PRECISION Gc_temp(iorg,num), P_temp(icomp)
      DOUBLE PRECISION WalVap(iorg,num)
      DOUBLE PRECISION kw_gas_off(iorg), before(icomp), after(icomp)
      DOUBLE PRECISION P_temp1(ibins,icomp), before1(ibins,icomp) ! for g-p-p on the wall
      DOUBLE PRECISION P_num(ibins), before2(ibins)	  
      DOUBLE PRECISION cstarT(IORG)
      DOUBLE PRECISION kw_HYBRID(IBINS)
      
      DATA alpha_w /1.0E-5/
      DATA wlc_zhang /10000.0/  ! in unit of ug/m3
CCCCC ..........................................................................................

      cdt = aadt/FLOAT(num)
      
      mu = 1.458E-6 * (TEMP**1.5) /(110.4 + TEMP) ! viscosity of air [kg/m-s]
      V_bag = BOXVOL/1.0e6 ! volume of filled bag [m3]
      R_bag = V_bag**(1./3.)
      S_bag = 6.*R_bag**2.
      StoV = S_bag / V_bag ! surface to volume ratio for smoke chamber
CCCCC ------for FSL chamber -------------------------------------------
C      R_bag1 = V_bag**(1./3.)
C      S_bag1 = 6.*R_bag1**2.
C      StoV1 = S_bag1/7. ! surface to volume ratio for smoke chamber
C      StoV = (12.4*12.4+12.4*19.6*2.)*2./V_bag	 !for FSL chamber
C      R_bag = (12.4*12.4*19.6)**(1./3.) !for FSL chamber
C      factor = StoV/StoV1
CCCCC -----for cylinder / rectangle reactor----------------------------
C      height = V_bag/(14.*14.*1e-4) !m
C      S1 = 14.0**2*1e-4 !m
C      S2 = height*14.0*1e-2 !m
C      StoV = (S1*2.+S2*4.)/V_bag ! surface to volume ratio for smoke chamber
	  
      factor = 1.0
      rhoSI = 1000*rho !density of particles, kg/m3
      gamma_p = 1.0

      do j=1,iorg
         cstarT(j)=cstar(j) * (298.0/TEMP) *
     &        EXP(Hvap(j)*1.0d3*(1.0d0/298.0-1.0d0/TEMP)/R)
C ------ Krechmer et al. 2016 ES&T ------------------------------------------
         if (cstarT(j).ge.1.0e4) then
            wlc(j) = 10.0e3 !30e3 ! [ug/m3] 
         else if ((cstarT(j).lt.1.0e4).and.(cstarT(j).gt.1.0)) then
            wlc(j) = 16.0*(cstarT(j)**0.6) ! [ug/m3]
         else
            wlc(j) = 16.0   ! [ug/m3]
         end if
C ------ Krechmer et al. 2016 ES&T ------------------------------------------
      end do	
	 
CCCCC Molecular weight of hypothetical compound, kg/mol, assumed from n-alkanes	 
      do j=1,IORG
         morg(j) = MWORG(j)/Na
      end do
	 
CCCCC Size calculation
      do i=1,ibins
         Dl=1.0e+6*((6.*xk(i))/(densorg(1)*pi))**(1./3.) ! ????? what is 1770?
         Dh=1.0e+6*((6.*xk(i+1))/(densorg(1)*pi))**(1./3.) ! ????  what is the unit? um?
         Dk(i)=sqrt(Dl*Dh)*1000. ! ??????? [nm]?
      end do
	  
CCCCC calculated aerosol properties
      do i=1,ibins
         kw_HYBRID(i) = 0.0
      end do
	 
      do i=1,ibins
         Cc(i) = 1. + (2*lambda/Dk(i)) *
     &        (1.257 + 0.4 * exp(-(0.55*Dk(i)/lambda)))
         B(i) = Cc(i) / (3*pi*mu*Dk(i)*1e-9) ! particle mobility, kg/s
         D(i) = kB * TEMP * B(i)       ! diffusivity, m2/s
         vT(i) = rhoSI * (Dk(i)*1e-9)**2 * g / (18.*mu)
         x(i) = pi * vT(i) / (2.*sqrt(kei*D(i))) !  skip the fitting
         D1(i) = debye1(x(i)) 
         kw_TURB(i) = 6.*sqrt(kei*D(i))/(pi*R_bag)*D1(i) +
     &        vT(i)/(4.*R_bag/3.) ! for cubic / spheric
c         kw_TURB(i) = 2.0*sqrt(kei*D(i))/pi*2*(S1+S2) +
c     &    vT(i) / tanh(pi*vT(i)/4./sqrt(kei*D(i))) ! for rectangle, cubic

         IF (PWL.eq.1) THEN
            kw_HYBRID(i) = kw0i*factor + kw_TURB(i)
         ELSE
            kw_HYBRID(i) = 0.0
         END IF
c         print*, StoV1,StoV,kw0i*factor*3600.,kw_TURB(i)*3600.,kw_HYBRID(i)*3600.
      end do


CCCCC loss particle to the wall
CCCCC for bulky non-organic species, use k_flat
      do j=1,icomp
         before(j) = 0.0
         after(j) = 0.0        
      end do	

      do n=1,ibins
         do j=1,icomp
	    before1(n,j) = 0.0
         end do
         before2(n) = 0.0
      end do
	
CCCCCCCCCCCCCCCCCCCCCC particle loss calculation, at 2014/12/11 CCCCCCCCCCCCCCCCCCCCCC
      do n=1,IBINS              !for size-resolved species, use k_turb
         before2(n) = Nki(n)
         Nki(n) = Nki(n)*exp(-dt*kw_HYBRID(n))
         P_num(n) = P_num(n) + before2(n)- Nki(n)
      end do
		
      do j=1,icomp
         do n = 1, IBINS        ! normal cases
!            if ((n .GE. 7) .AND. (n .LE. 22)) then ! for comparison with APE model
!            if (n .LE. 26) then !for coarse mode calculation
            before(j) = before(j) + Mki(n,j) !for normal case
!            before1(n,j) = Mki(n,j)  !! for g-p-p on the wall
!            endif
            Mki(n,j)= Mki(n,j)* exp(-dt*kw_HYBRID(n))
!            if ((n .GE. 7) .AND. (n .LE. 22)) then  ! for comparison with APE model
!            if (N .LE. 26) then ! for coarse mode calculation
            after(j) = after(j) + Mki(n,j)
!            P_temp1(n,j) = P_temp1(n,j) + before1(n,j) - Mki(n,j) ! for g-p-p on the wall
!	     endif
         end do
         P_temp(j) = P_temp(j)+before(j)-after(j)
      end do
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! Question: ke can be applied on both of the number and mass size distribution?
!   organic gas properties
      alpha_w = 1.e-3           ! mass accomodation coefficient of organics onto chamber wall
	   !(i.e., sorption; Matsunaga and Ziemann, 2010)
      
      mair = MWair/Na;
      dorg = 10e-10 ! diameter of organic molecule [m]
      dair = (0.79*1.09e-10 + 0.21*1.21e-10) ! diameter of air molecule, m

      IF (VWL .eq. 0) THEN
         kw_gas(j) = 0.0
         kw_gas_off(j) = 0.0
      ELSE
C      print*, '========================================================'
C         print*, 'j  -  cmean(j), kw_gas(j), Dg(j), kw_gas_off(j),'
C     &     //'cstarT(j)'
         do j=1,iorg
            cmean(j) = sqrt(3.*R*TEMP/MWORG(j)) ! root mean square speed of gas, m/s
            Dg(j) = (2./3.)*sqrt(kB*TEMP/pi*0.5*(1./morg(j)+1/mair))*1./ 
     &           (pi*(0.5*(dorg + dair))**2.)/Na*(R*TEMP/pres) ! diffusivity of organic vapors, m2/s
C         Dg(j) = diorg(j)*(TEMP/298.0)**ditemporg(j)
C            kw_gas(j) = StoV*(alpha_w*cmean(j)/4.)/
C     &           (1.+(pi/2.)*(alpha_w*cmean(j)/(4.*sqrt(kei*Dg(j))) ) ) ! Zhang et al., 2014
            kw_gas(j) = StoV*(2./pi)*(sqrt(kei*Dg(j))) ! Zhang et al., 2014

C         solubility(j,m,time) = (cmean(j,m)*sqrt(pi*time)/
C     &        (4*10^logHb(j,m)*Ra*TEMP(j)*sqrt(Dg(j,m))))^-1 ! inverse uptake based on Henry's law
         ! uptake accounting for mass accommodation and solubility, uptake coeffcient
C         gamma(j) = (1./alpha_w + cmean(j)*sqrt(pi*time)/
C     &        (4.*10.^logHb(j)*Ra*TEMP*sqrt(Dg(j))))^-1.
C         kw_gas(j) = StoV*(gamma(j)*cmean(j)/4.)/(1.+(pi/2.)*
C     &        (gamma(j)*cmean(j)/(4.*sqrt(kei*Dg(j))))) ! gas-phase wall loss rate, 1/s
C         kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc/200./gamma_p
C         kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc/(MWORG(j)*1000.0)/gamma_p
C         kw_gas_off(j) = kw_gas(j)*cstarT(j)/wlc(j)/(MWORG(j)*1000.0)/gamma_p
            kw_gas_off(j) = kw_gas(j) * cstarT(j) / wlc(j) !wlc_zhang
            tau_w(j) = 1./(kw_gas(j)+kw_gas_off(j))
C         write(*,*) 'vapor loss rate constant'
C         print *,'kw_gas=',kw_gas(j),kw_gas_off(j) !,Dg(j)	 
C            print*, j, '-', cmean(j), kw_gas(j), Dg(j), kw_gas_off(j),
C     &           cstarT(j)
         end do	
      END IF

C      print*, '========================================================'
CCCCC swap Gc array into Gc_temp array dependent on time step
      do j=1,iorg
         Gc_temp(j,1) =  Gci(j+1)
C         write(*,*) 'initial vapor concentration'
C         print*, Gc_temp(j,1)/7000000.*1.E15
C         Gc_temp2(j,1) = 0.0
         WalVap(j,1) = WalVapMas(j)
      end do
	 
CCCCC Euler's method
C      do k=2,num   ! time step
C         do j=1,iorg !species
C            if ((WalVap(j,k-1)+Gc_temp(j,k-1)*
C     &           kw_gas(j)-WalVap(j,k-1)*kw_gas_off(j))
C     &           .GE. 0.0) then ! added in vbs extension
C               Gc_temp(j,k) = Gc_temp(j,k-1)+WalVap(j,k-1)*kw_gas_off(j)-Gc_temp(j,k-1)*kw_gas(j)
C               WalVap(j,k) = WalVap(j,k-1)+Gc_temp(j,k-1)*kw_gas(j)-WalVap(j,k-1)*kw_gas_off(j)
C            else
C               Gc_temp(j,k) = Gc_temp(j,k-1) + WalVap(j,k-1)
C               WalVap(j,k) = 0.0
C            endif
C         end do
C      end do		

CCCCC integrate method
      do k=2,num   !time step
         do j=1,iorg !species
            Gc_temp(j,k) = Gc_temp(j,k-1)*exp(-kw_gas(j)*cdt) +
     &           WalVap(j,k-1)*(1-exp(-kw_gas_off(j)*cdt))
            WalVap(j,k) = WalVap(j,k-1)*exp(-kw_gas_off(j)*cdt) +
     &           Gc_temp(j,k-1)*(1-exp(-kw_gas(j)*cdt))
         end do
      end do
 

!   Swap Gc_temp back into Gc array and Swap out vapour wall loss
      do j=1,iorg
         Gci(j+1) = Gc_temp(j,num)
         if (WalVap(j,num) .LT. 1e-30) then !add in vbs extension
            WalVapMas(j) = 0.0
         else
	    WalVapMas(j) = WalVap(j,num)
         endif
      end do
	     
	
 100  Format(1X,E10.3)	
      END subroutine
	
 
	
!debye1 function, taken from MISFUN	
	
      function debye1 ( xvalue )

!*****************************************************************************80
!
!! DEBYE1 calculates the Debye function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the function.
!
!    Output, real ( kind = 8 ) DEBYE1, the value of the function.
!
      implicit none

      double precision adeb1(0:18)
      double precision cheval
      double precision debye1
      double precision, parameter :: eight = 8.0D+00
      double precision, parameter :: four = 4.0D+00
      double precision, parameter :: half = 0.5D+00
      integer ( kind = 4 ) i
      integer ( kind = 4 ) nexp
      integer ( kind = 4 ), parameter :: nterms = 15
      double precision, parameter :: one = 1.0D+00
      double precision, parameter :: quart = 0.25D+00
      double precision x
      double precision xvalue
      double precision, parameter :: zero = 0.0D+00

      double precision debinf, expmx,
     &     nine, rk, sum1, t, thirt6, xk, xlim, xlow,
     &     xupper
      data nine,thirt6 /9.0d0, 36.0d0 /
      data debinf/0.60792710185402662866d0/
      data adeb1/2.40065971903814101941d0,
     &           0.19372130421893600885d0,
     &          -0.623291245548957703d-2, 
     &           0.35111747702064800d-3, 
     &          -0.2282224667012310d-4, 
     &           0.158054678750300d-5, 
     &          -0.11353781970719d-6, 
     &           0.835833611875d-8, 
     &          -0.62644247872d-9, 
     &           0.4760334890d-10, 
     &          -0.365741540d-11, 
     &           0.28354310d-12, 
     &          -0.2214729d-13, 
     &           0.174092d-14, 
     &          -0.13759d-15, 
     &           0.1093d-16, 
     &          -0.87d-18,
     &           0.7d-19,
     &          -0.1d-19/
!
!   Machine-dependent constants (suitable for IEEE machines)
!
      data xlow,xupper,xlim/0.298023d-7,35.35051d0,708.39642d0/

      x = xvalue

      if ( x < zero ) then
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'DEBYE1 - Fatal error!'
         write ( *, '(a)' ) '  Argument X < 0.'
         debye1 = zero
      else if ( x < xlow ) then
         debye1 = ( ( x - nine ) * x + thirt6 ) / thirt6
      else if ( x <= four ) then
         t = ( ( x * x / eight ) - half ) - half
         debye1 = cheval ( nterms, adeb1, t ) - quart * x
      else

         debye1 = one / ( x * debinf )
         if ( x < xlim ) then
            expmx = exp ( -x )
            if ( xupper < x ) then
               debye1 = debye1 - expmx * ( one + one / x )
            else
               sum1 = zero
               rk = aint ( xlim / x )
               nexp = int ( rk )
               xk = rk * x
               do i = nexp, 1, -1
                  t = ( one + one / xk ) / rk
                  sum1 = sum1 * expmx + t
                  rk = rk - one
                  xk = xk - x
               end do
               debye1 = debye1 - sum1 * expmx
            end if
         end if
      end if
      
      return
      end


      subroutine debye1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DEBYE1_VALUES returns some values of Debye's function of order 1.
!
!  Discussion:
!
!    The function is defined by:
!
!      DEBYE1(x) = 1 / x * Integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
!
!    The data was reported by McLeod.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
      implicit none

      integer ( kind = 4 ), parameter :: n_max = 20
      
      double precision fx
      double precision, save, dimension ( n_max ) :: fx_vec = (/ 
     &     0.99951182471380889183D+00, 
     &     0.99221462647120597836D+00,
     &     0.96918395997895308324D+00,
     &     0.88192715679060552968D+00,
     &     0.77750463411224827642D+00,
     &     0.68614531078940204342D+00,
     &     0.60694728460981007205D+00,
     &     0.53878956907785587703D+00,
     &     0.48043521957304283829D+00,
     &     0.38814802129793784501D+00,
     &     0.36930802829242526815D+00,
     &     0.32087619770014612104D+00,
     &     0.29423996623154246701D+00,
     &     0.27126046678502189985D+00,
     &     0.20523930310221503723D+00,
     &     0.16444346567994602563D+00,
     &     0.10966194482735821276D+00,
     &     0.82246701178200016086D-01,
     &     0.54831135561510852445D-01,
     &     0.32898681336964528729D-01 /)
      integer ( kind = 4 ) n_data
      double precision x
      double precision, save, dimension ( n_max ) :: x_vec = (/ 
     &     0.0019531250D+00,
     &     0.0312500000D+00,
     &     0.1250000000D+00,
     &     0.5000000000D+00,
     &     1.0000000000D+00,
     &     1.5000000000D+00,
     &     2.0000000000D+00,
     &     2.5000000000D+00,
     &     3.0000000000D+00,
     &     4.0000000000D+00,
     &     4.2500000000D+00,
     &     5.0000000000D+00,
     &     5.5000000000D+00,
     &     6.0000000000D+00,
     &     8.0000000000D+00,
     &     10.0000000000D+00,
     &     15.0000000000D+00,
     &     20.0000000000D+00,
     &     30.0000000000D+00,
     &     50.0000000000D+00 /)

      if ( n_data < 0 ) then
         n_data = 0
      end if

      n_data = n_data + 1
      
      if ( n_max < n_data ) then
         n_data = 0
         x = 0.0D+00
         fx = 0.0D+00
      else
         x  = x_vec(n_data)
         fx = fx_vec(n_data)
      end if
      
      return
      
      end

      function cheval ( n, a, t )

!*****************************************************************************80
!
!! CHEVAL evaluates a Chebyshev series.
!
!  Discussion:
!
!    This function evaluates a Chebyshev series, using the
!    Clenshaw method with Reinsch modification, as analysed
!    in the paper by Oliver.
!
!  Modified:
!
!    07 August 2004
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!    J Oliver,
!    An error analysis of the modified Clenshaw method for
!    evaluating Chebyshev and Fourier series,
!    Journal of the IMA,
!    Volume 20, 1977, pages 379-391.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of terms in the sequence.
!
!    Input, real ( kind = 8 ) A(0:N), the coefficients of the Chebyshev series.
!
!    Input, real ( kind = 8 ) T, the value at which the series is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) CHEVAL, the value of the Chebyshev series at T.
!
      implicit none
      
      integer ( kind = 4 ) n
      
      double precision :: a(0:n)
      double precision :: cheval
      double precision :: d1
      double precision :: d2
      double precision, parameter :: half = 0.5D+00
      integer ( kind = 4 ) i
      double precision :: t
      double precision, parameter :: test = 0.6D+00
      double precision :: tt
      double precision, parameter :: two = 2.0D+00
      double precision :: u0
      double precision :: u1
      double precision :: u2
      double precision, parameter :: zero = 0.0D+00
      
      u1 = zero
C ---------------------------------------------------     
C      T <= -0.6, Reinsch modification.
C ---------------------------------------------------
      if ( t <= -test ) then         
         d1 = zero
         tt = ( t + half ) + half
         tt = tt + tt
         do i = n, 0, -1
            d2 = d1
            u2 = u1
            d1 = tt * u2 + a(i) - d2
            u1 = d1 - u2
         end do
         
         cheval = ( d1 - d2 ) / two
C ---------------------------------------------------
C     -0.6 < T < 0.6, Standard Clenshaw method.
C ---------------------------------------------------
      else if ( t < test ) then         
         u0 = zero
         tt = t + t         
         do i = n, 0, -1
            u2 = u1
            u1 = u0
            u0 = tt * u1 + a(i) - u2
         end do         
         cheval = ( u0 - u2 ) / two
C ---------------------------------------------------
C     0.6 <= T, Reinsch modification.
C ---------------------------------------------------
      else
         d1 = zero
         tt = ( t - half ) - half
         tt = tt + tt
         do i = n, 0, -1
            d2 = d1
            u2 = u1
            d1 = tt * u2 + a(i) + d2
            u1 = d1 + u2
         end do
         cheval = ( d1 + d2 ) / two
         
      end if
      
      return
      end
