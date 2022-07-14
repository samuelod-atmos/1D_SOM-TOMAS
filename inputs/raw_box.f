C     **************************************************
C     *  box                                           *
C     **************************************************

C     WRITTEN BY Peter Adams, October 1999

C     This program is a box model that is used to test the size
C     resolved aerosol code for the GISS GCM.

c      subroutine box(nucbase,H2SO4rate,life)a
      PROGRAM box

      IMPLICIT NONE
      
C-----INCLUDE FILES--------------------------------------------------

      include 'sizecode.COM'

      ! AliA ////////////////////////////////////////////////////////
      include 'modlspc.h'  
      include 'common.inc'
      include 'gaskin.h' 
      ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ AliA

      
C-----INPUTS---------------------------------------------------------

C-----VARIABLE DECLARATIONS------------------------------------------

      integer i,j,l,n,k, ttt      !counters
      integer iseed, getlen
      integer nsteps            !number of timesteps per hour
      integer hourint
      integer nintern,tstep ! number of internal mphys steps
      integer tvar 
      INTEGER INDOMFRACI

      logical dayflag, clearflag  !true/false for day/night or clear/cloudy

      external getlen, setoh

      double precision No, Dp, sigma !lognormal parameters (No in #/cm3,Dp in um)
      double precision Nc, Dpc,sigmac
      double precision time     !elapsed time (seconds)
      double precision hour  !elapsed time (hours)
      double precision adt,aadt,aadtn      !aerosol microphysics time step (seconds)
      double precision aadtp ! previous dynamic timestep
      double precision Dpk, vd  !particle diameter (um) and drydep velocity (cm/s)
      double precision mnucl    !mass nucleating in a time step (kg)
      double precision Ntot     !used to check number conservation by condensation
      double precision area     !xsect area of plume (m2)
      double precision boxvolp     !previous boxvol (m3)
      double precision depth    ! depth of plume being tracked (m)
      double precision ohc   !OH conc (molec cm-3)
      double precision cyclet   !how far into 12-hour diurnal cycle we are (s)
      double precision setoh
      double precision u10
      double precision focean
      double precision endtime
      double precision tot_n_1, tot_n_2, tot_num
      double precision tot_s_1, tot_s_2
      double precision test     !test criteria for if we need to call report
      double precision tau(IBINS)
      double precision Nkout(IBINS), Mkout(IBINS,icomp)
      double precision Gcout(icomp-1)
      double precision Nkt(IBINS), Mkt(IBINS,icomp)
      double precision Gct(icomp-1)
      double precision nucrate ! nucleation rate [cm-3 s-1]
      double precision dilt, dilratio ! dilution timescale [hr]
      double precision timerem
      double precision soarate ! soa formation rate [ug m-3 hr-1]
      double precision sulfrate ! sulf formation rate [ug m-3 hr-1]
      double precision massin ! mass to condense in timestep [kg]
      double precision tmass1,tmass2
      double precision tt, tn, tsu, toa ! temporary variables
      double precision timeall(500) ! all times from input [hr]
      double precision nucall(500) ! all nuc rates from input [cm-3 s-1]
      double precision sulfall(500) ! all sulfate conc from input [molec cm-3]
      double precision soaall(500) ! all soa rates from input [ug m-3 hr-1]
      double precision tempnum
      double precision CS, sinkfrac(IBINS)
      double precision tmp
      double precision R ! gas constant
      double precision orgfracinit ! initial fraction of mass that is organic
      double precision soascale ! scale factor for VBS mass
      double precision totamass ! sum to get total aerosol mass [ug m-3]
      double precision totmass ! sum to get total mass [ug m-3]
      double precision faero(iorg) ! fraction of POA initially in the aerosol phase
      double precision ampb(iorg) ! organic aerosol mass in each volatility bin
      double precision fampb(iorg) ! fraction of organic aerosol mass in each volatility bin
      double precision guess ! check for initializing mass	  
      double precision kw0,ke!dilt2, kw_turb(36) !add by Emily
      double precision WalVapMas(iorg), P_temp(icomp)
      double precision P_temp1(ibins,icomp), P_num(ibins)	  ! add by Emily
      double precision tau_p, tau_w(iorg),kw_HYBRID(ibins)
      double precision Nkout1(IBINS), Mkout1(IBINS,icomp),Gc_v(icomp-1)
      double precision Nk1(IBINS), Mk1(IBINS,icomp)	  
      double precision alpha_i(5), kOH(15), Gc_i(15,15)
      double precision OH_expo, Gc_ii(15) !edit by Emily for SOA chem, 11/30/2015
      DOUBLE PRECISION omfraci_val
      
      character*120 runname, datname, diagname, infilename, vapname
      CHARACTER*120 NCONCFILE, GCONCFILE, AECONCFILE, SPECFILE
      CHARACTER*120 SAPRCGCFILE
      CHARACTER*120 DIRECTORY
      character*120 tauname,kwname,cglname,walname
      character*120 binfilename

      common /FLNAMES/ datname, diagname,
     &                 NCONCFILE, GCONCFILE, AECONCFILE
	  
CCC   AliA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      INTEGER IERR
      REAL*4 SLA

      CHARACTER*16 MODL
      CHARACTER*16 SP
      CHARACTER*8 MCHNAME1, MCHNAME2
      
      CHARACTER*7, ALLOCATABLE :: PRECNAME(:)
      CHARACTER*16, ALLOCATABLE :: SOMPREC(:)
      
      INTEGER STRINX, ENDINX, NSOMPREC, TNSOMSP
      INTEGER CTR, CTRN, OFFSET, CTRW
      INTEGER CARBON, HYDROGEN, OXYGEN
      INTEGER NINPOXYSP, NINPPRECSP
      INTEGER OXYSPEMISS, PRECSPEMISS
      INTEGER NEQNT
      INTEGER VWL, PWL, COAG
      
      INTEGER, ALLOCATABLE :: PRECSTRINX(:)
      INTEGER, ALLOCATABLE :: NSOMSP(:)
      INTEGER, ALLOCATABLE :: SOMCARBON(:)
      INTEGER, ALLOCATABLE :: SOMOXYGEN(:)
      INTEGER, ALLOCATABLE :: SOMHYDROGEN(:)
      INTEGER, ALLOCATABLE :: INPC(:), INPO(:)
      INTEGER, ALLOCATABLE :: PRECLENNAME(:)
      
      REAL, ALLOCATABLE :: SAPRCGC(:)
      REAL, ALLOCATABLE :: IPPMPREC(:)
      REAL, ALLOCATABLE :: IPPMOXY(:)      
      REAL, ALLOCATABLE :: DLVP(:)
      REAL TEMP1, TIN, TOUT
      REAL GCREAL(IORG)
      REAL CALVALREAL
      DOUBLE PRECISION CALVALDBLE
      DOUBLE PRECISION MWORG_HC(IORG)

      DOUBLE PRECISION Dbk ! particle-phase diffusion coefficient [m2/s] 
      DOUBLE PRECISION Kc  ! first-order loss rate of species in the particle phase [1/s]

      REAL a_oh, ax_oh, b_oh, bx_oh ! OH parameters for changing
                                    ! the OH concentration over time
      
      DATA  MODL / 'SAPRC14_REV1' /

      EXTERNAL DIFF
      
CCC   ////////////////////////////////////////////////////////////////////////// AliA

C-----VARIABLE COMMENTS----------------------------------------------

C-----ADJUSTABLE PARAMETERS------------------------------------------

      parameter(R=8.314)
      !parameter(No=0.d0, Dp=0.1, sigma=2.0)

      !SPECIAL TREATMENT FOR PARMETER NAMES
      PARNAME(1) = 'SO4'
      PARNAME(ICOMP-1) = 'NH4'
      PARNAME(ICOMP) = 'H2O'

C-----CODE-----------------------------------------------------------

      write(*,*) 'Enter name of run'
      read(*,'(A120)') runname
C      runname = 'burn'

      datname='../outputs/'//runname(1:getlen(runname))//
     &     '.dat'
      diagname='../outputs/'//runname(1:getlen(runname))//
     &     '.diag'
      vapname='../outputs/'//runname(1:getlen(runname))//
     &     '_vl.dat'
      tauname='../outputs/'//runname(1:getlen(runname))//
     &     '_tau.dat'
      kwname='../outputs/'//runname(1:getlen(runname))//
     &     '_kw.dat'  
      walname='../outputs/'//runname(1:getlen(runname))//
     &     '_wal.dat'
      NCONCFILE='../outputs/'//RUNNAME(1:GETLEN(RUNNAME))//
     &     '_noconc.dat'
      GCONCFILE='../outputs/'//RUNNAME(1:GETLEN(RUNNAME))//
     &     '_gc.dat'
      SAPRCGCFILE='../outputs/'//RUNNAME(1:GETLEN(RUNNAME))//
     &     '_saprcgc.dat'
      AECONCFILE='../outputs/'//RUNNAME(1:GETLEN(RUNNAME))//
     & '_aemass.dat'
      SPECFILE='../outputs/'//RUNNAME(1:GETLEN(RUNNAME))//
     &     '_spec.dat'
      write(*,*) datname
      write(*,*) diagname
      write(*,*) vapname
      write(*,*) kwname
      write(*,*) walname
      WRITE(*,*) NCONCFILE
      WRITE(*,*) GCONCFILE
      WRITE(*,*) SAPRCGCFILE
      WRITE(*,*) AECONCFILE
      WRITE(*,*) SPECFILE
      open(unit=30,file=datname,status='new')
      close(30)
      open(unit=31,file=vapname,status='new')
      close(31)	 
      open(unit=32,file=tauname,status='new')
      close(32)
      open(unit=33,file=kwname,status='new')
      close(33)
      OPEN(UNIT=34, FILE=NCONCFILE, STATUS='new')
      CLOSE(34)
      OPEN(UNIT=35, FILE=GCONCFILE, STATUS='new')
      CLOSE(35)
      OPEN(UNIT=36, FILE=AECONCFILE, STATUS='new')
      CLOSE(36)
      OPEN(UNIT=37, FILE=SPECFILE, STATUS='new')
      CLOSE(37)
      OPEN(UNIT=38, FILE=SAPRCGCFILE, STATUS='new')
      CLOSE(38)
	 	
      WRITE(*,*) 'Enter 0/1 to swith Off/On coagulation'
      read(*,'(I1)') COAG
C      VWL = 1
      print*, 'COAG=', COAG

      WRITE(*,*) 'Enter 0/1 to swith Off/On vapor wall losses'
      read(*,'(I1)') VWL
C      VWL = 1
      print*, 'VWL=', VWL

      WRITE(*,*) 'Enter 0/1 to swith Off/On particle wall losses'
      read(*,'(I1)') PWL
C      PWL = 1
      print*, 'PWL=', PWL

      read(*,*) ke
C      ke = 1.074685063
      print*,'ke=', ke
      
      read(*,*) kw0
C      kw0 = 8.03E-05 
      print*, 'kw0=', kw0
      
      WRITE(*,*) 'Enter OH concentration parameters [molecules cm-3]'
      WRITE(*,*) 'a_oh*exp(-1.*ax_oh*t) + b_oh*exp(-1.*bx_oh*t)'
      WRITE(*,*) 't in hours'
      READ(*,*) a_oh
      PRINT*, 'a_oh =', a_oh
      READ(*,*) ax_oh
      PRINT*, 'ax_oh =', ax_oh
      READ(*,*) b_oh
      PRINT*, 'b_oh =', b_oh
      READ(*,*) bx_oh
      PRINT*, 'bx_oh =', bx_oh

      write(*,*) 'Enter nucleation rate [cm-3 s-1]'
      read(*,'(E5.1)') nucrate
C      nucrate = 0.0
      print*, 'nucrate=', nucrate

      write(*,*) 'Enter dilution ratio'
      read(*,'(E5.1)') dilt
C      dilt = 2.493
      print*, 'dilt=', dilt

      write(*,*) 'Enter the size of the bag'
      read(*,'(E15.1)') BOXVOL
C      boxvol = 7000000.0
      print*, 'boxvol=', BOXVOL
      
      write(*,*) 'Enter the hours of run time (make even hours)'
      read(*,'(F8.5)') endtime
C      endtime = 4.0
      print*, 'endtime=', endtime

      write(*,*) 'Enter the accommodation coefficient (alpha)'
      read(*,'(F8.5)') alpha
C      alpha = 1.0
      print*, 'alpha=', alpha

      write(*,*) 'particle-phase diffusion coefficient [m2/s] (Dbk)'
      read(*,'(F25.23)') Dbk
C      Dbk = 1.0d-10
      print*, 'Dbk =', Dbk      
      
      write(*,*) 'first-order loss rate of species in the'
     &     //' particle phase'
      read(*,'(F10.8)') kc
C      Dbk = 1.0d-10
      print*, 'kc =', kc     
      
write(*,*) 'Enter surface tension [N/m]'
      read(*,'(F8.5)') storg
C      storg = 0.025
      print*, 'storg=', storg

      write(*,*) 'Enter pressure [Pa]'
      read(*,'(F12.3)') pres
C      temp = 292.85
      print*,'pres=', pres

      write(*,*) 'Enter temperature [K]'
      read(*,'(F8.5)') temp
C      temp = 292.85
      print*,'temp=', temp
      
      write(*,*) 'Enter relative humidity'
      read(*,'(8F8.5)') rh
C      rh = 0.2
      print*,'rh=', rh

      write(*,*) 'Enter background number conc'
      read(*,'(F8.5)') No
C      No = 5843.0
      print*, 'No=', No

      write(*,*) 'Enter background median diameter [microns]'
      read(*,'(F8.5)') Dp
C      Dp = 0.157
      print*, 'Dp=', Dp

      write(*,*) 'Enter background sigma'
      read(*,'(F8.5)') sigma
C      sigma = 1.73
      print*, 'sigma=', sigma

CCC   AliA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      READ(*,'(I1)') PRECSPEMISS
C      PRECSPEMISS = 1
      print*, 'PRECSPEMISS=', PRECSPEMISS
      IF (PRECSPEMISS .EQ. 1) THEN
         READ(*,'(I3)') NINPPRECSP
C         NINPPRECSP = 1
         print*, 'NINPPRECSP=', NINPPRECSP
         
         ALLOCATE(IPPMPREC(NINPPRECSP),
     &            PRECNAME(NINPPRECSP),
     &            PRECLENNAME(NINPPRECSP))
         READ(*,*) PRECNAME
C         PRECNAME(1) = 'ALK5'
         print*, 'PRECNAME=', PRECNAME
         READ(*,*) PRECLENNAME
C         PRECLENNAME(1) = 4
         print*, 'PRECLENNAME=', PRECLENNAME
         READ(*,*) IPPMPREC
C         IPPMPREC(1) = 0.001
         print*, 'IPPMPREC=', IPPMPREC
      END IF

      
      READ(*,'(I1)') OXYSPEMISS
C      OXYSPEMISS = 0
      print*, 'OXYSPEMISS=', OXYSPEMISS
      IF (OXYSPEMISS .EQ. 1) THEN
         READ(*,'(I3)') NINPOXYSP

         ALLOCATE(INPC(NINPOXYSP),
     &            INPO(NINPOXYSP),
     &            IPPMOXY(NINPOXYSP))
         READ(*,*) INPC
         READ(*,*) INPO
         READ(*,*) IPPMOXY
      END IF

      READ(*,'(I2)') NSOMPREC
C      NSOMPREC = 1
      print*, 'NSOMPREC=', NSOMPREC
      ALLOCATE(SOMPREC(NSOMPREC),
     &         DLVP(NSOMPREC))
      READ(*,*) SOMPREC
C      SOMPREC(1) = 'ALKSOMG'
      print*, 'SOMPREC=', SOMPREC
      READ(*,*) DLVP
C      DLVP(1) = 1.21
      print*, 'DLVP=', DLVP
      
CCC   /////////////////////////////////////////////////////////////////////////// AliA
      
C-------initialize wall vapour concentration---------------------
      do j = 1,iorg
        WalVapMas(j) = 0.0
      enddo
      do j = 1,icomp
        P_temp(j) = 0.0
      enddo	
      do j = 1,iorg
         tau_w(j) = 0.0
      enddo	
      tau_p = 0.0
  
      do n = 1,ibins
         do j = 1,icomp			
            P_temp1(n,j) = 0.0
            Mk1(n,j) = 0.0
         enddo
         P_num(n) = 0.0
         Nk1(n) = 0.0
      enddo
      
      do j=1,icomp-1
         Gc_v(j) = 0.0
      end do

      do n=1,15
         do j=1,15
            Gc_i(n,j) = 0.0
         enddo
      enddo
      
C------------------------------------------------------------
   
      time=0.0
      adt=10.0*60.0 ! timestep for writing output [seconds]
      nintern=60  ! number of internal microphysics steps
      aadt=adt/nintern ! microphysics timestep
      print*,'adt',adt
      depth=1.          ! m (will be held fixed)

      boxmass=.0289*pres*boxvol*1.0e-6/R/temp
      hour=time/3600.

C AliA ///////////////////////////////////////////////////////////////////////////////////////////////
      ! zero emission
      EMT(:) = 0.0
      
      ! SAPRC Initialization
      OPEN(UNIT=1, FILE='saprc14_rev1.mod', STATUS='old')
      OPEN(UNIT=3, FILE='saprc14_rev1.doc', STATUS='unknown')
      
      TREF = 298.0
      TEMP1 = temp
      CALL RDMOD(MODL, 1, 3, .true., TEMP1, 0, IERR,
     +     NS, NACT, NSC, NSB, NDUM,
     +     NRXN, LENKBF, NCOC1, NCOC2,
     +     NCOEFV, NPHOTK, TREF, SLA,
     +     DECLN, TSTART, CZMIN, MWT,
     +     CONC0, CNO, NNO, SNO,
     +     ONO, XNO, NRHC1, INDRHC,
     +     RK, KPBUF, COEF, ABSQY,
     +     DPHOTK, LKBUF, IPHRK, IRXTYP,
     +     NAME, RXNLBL, COEFNM,
     +     PHOTNM, TITLE) 
      call flush(6) 
      CLOSE(UNIT=1)
      CLOSE(UNIT=3)

CCC   DEFINE OFFSET BETWEEN ALL SPECIES AND ACTIVE SPECIES
      OFFSET = 13               ! ************ WARNING - CAN BE CHANGED BASED ON SAPRC FILE ***********
         
CCC   FINDING THE STARTING INDEX FOR SOM SPECIES BASED ON DIFFERENT PRECURSORS
      ALLOCATE(PRECSTRINX(NSOMPREC), NSOMSP(NSOMPREC))
      
      PRECSTRINX(:) = 0         ! starting index of SOM species
      NSOMSP(:) = 0             ! how many species exist for that specific precursor
      DO N=1,NSOMPREC           ! Number of SOM precursores in input
         CTR = 0
         CTRN = 0
         DO I=1,NACT            ! loop through total species
            CTR = CTR + 1
            READ(NAME(I+OFFSET)(1:7), '(A7)') SP ! read the name of the species
            IF (SP .EQ. SOMPREC(N)) THEN ! check if it matches with the precursor name
               CTRN = CTRN + 1
               ENDINX = CTR
c               print*, 'SP=', SP, '  SOMPREC(N)=', SOMPREC(N)
c               print*, 'CTRN=', CTRN, 'ENDINX=', ENDINX
            END IF
         END DO
         NSOMSP(N) = CTRN
         PRECSTRINX(N) = ENDINX - CTRN + 1
      END DO

      OPEN(UNIT=37, FILE=SPECFILE, STATUS='old', ACCESS='APPEND')
 900  FORMAT(307A20)
      WRITE(37,900) 'ACTIVE GAS SP:', NAME(1+OFFSET:NACT+OFFSET)
      CLOSE(37)
      
CCC   FINDING THE FIRST SOM SPECIES
      CTR = 0
      DO I=1,NACT ! loop through total species
         CTR = CTR + 1
         READ(NAME(I+OFFSET)(4:6), '(A7)') SP ! read the name of the species 
         IF (SP .EQ. 'SOM') THEN ! check if it matches with the precursor name
            STRINX = CTR
C            print*, 'STRINX=', STRINX
            GOTO 200
         END IF
      END DO
 200  CONTINUE

      OPEN(UNIT=37, FILE=SPECFILE, STATUS='old', ACCESS='APPEND')      
 901  FORMAT(202A20)
      WRITE(37,901) 'SOM SP:', NAME(STRINX+OFFSET:NACT+OFFSET)
      CLOSE(37)
      
CCC   MAKE THE CARBON, OXYGEN, HYDROGEN ARRAY
      TNSOMSP = SUM(NSOMSP)     ! total number of SOM species
      print*, 'TNSOMSP=', TNSOMSP
      ALLOCATE(SOMCARBON(TNSOMSP), SOMOXYGEN(TNSOMSP),
     &     SOMHYDROGEN(TNSOMSP)) ! allocate the size of C, O, H arrays
      DO I=STRINX,NACT
         READ(NAME(I+OFFSET)(4:6), '(A7)') SP ! read the name of the species
         IF (SP .EQ. 'SOM') THEN
            READ(NAME(I+OFFSET)(9:10),'(I2)') CARBON ! read the carbon number
            READ(NAME(I+OFFSET)(12:13),'(I2)') OXYGEN ! read the oxygen number
            SOMCARBON(I-STRINX+1) = CARBON ! save the carbon number
            SOMOXYGEN(I-STRINX+1) = OXYGEN ! save the oxygen number
            SOMHYDROGEN(I-STRINX+1) = (CARBON*2)+2-OXYGEN ! calculate the hydrogen number
            !print*, NAME(I+OFFSET), CARBON, OXYGEN
            MWORG(I-STRINX+1) = CARBON*12.0107 + OXYGEN*15.999
     &           + ((CARBON*2)+2-OXYGEN)*1.00794
            MWORG_HC(I-STRINX+1) = CARBON*12.0107 !+ OXYGEN*15.999
     &           + ((CARBON*2)+2-OXYGEN)*1.00794
         END IF
      END DO

      OPEN(UNIT=37, FILE=SPECFILE, STATUS='old', ACCESS='APPEND')      
 902  FORMAT(A20,201I5)
 903  FORMAT(A20,201F10.3)
      WRITE(37,902) 'SOM CARBON NO:', SOMCARBON(1:TNSOMSP)
      WRITE(37,902) 'SOM OXYGEN NO:', SOMOXYGEN(1:TNSOMSP)
      WRITE(37,903) 'MWT of SOM:', MWT(STRINX+OFFSET:NACT+OFFSET)
      WRITE(37,903) 'MWORG of SOM:', MWORG(1:TNSOMSP)
      WRITE(37,903) 'MWORG_HC of SOM:', MWORG_HC(1:TNSOMSP)
      CLOSE(37)

CCC   MAKE AN ARRAY OF ACTIVE GAS PHASE ACTIVE SPECIES IN SAPRC
      ALLOCATE(SAPRCGC(NACT))   ! initialize the array size of total gas concentration
      SAPRCGC(:) = 0.0

      IF (PRECSPEMISS .EQ. 1) THEN
         print*, 'NINPPRECSP=', NINPPRECSP
         DO N=1,NINPPRECSP
            DO I=1,NACT
               print*, NAME(I+OFFSET), PRECNAME(N)
               MCHNAME1 = NAME(I+OFFSET)(1:PRECLENNAME(N))
               MCHNAME2 = PRECNAME(N)(1:PRECLENNAME(N))
               IF (MCHNAME1 .EQ. MCHNAME2) THEN
                  SAPRCGC(I) = IPPMPREC(N)
                  GOTO 201
               END IF
            END DO
 201        CONTINUE
         END DO
      END IF
         
      IF (OXYSPEMISS .EQ. 1) THEN
         DO N=1,NINPOXYSP
            DO I=1,TNSOMSP
               IF ((SOMCARBON(I).EQ.INPC(N)).AND.
     &              (SOMOXYGEN(I).EQ.INPO(N))) THEN
                  SAPRCGC(I+STRINX-1) = IPPMOXY(N)
                  GOTO 202
               END IF
            END DO
 202        CONTINUE
         END DO
      END IF

CCC   CSTAR CALCULATION BASED ON dlvp, OXYGEN No., and MW
      !ALLOCATE(CSTAR(TNSOMSP), HVAP(TNSOMSP))
      DO N=1,NSOMPREC
       ENDINX = PRECSTRINX(N)+NSOMSP(N)-1
       DO I=PRECSTRINX(N),ENDINX
        PARNAME(I-STRINX+1+1) = NAME(I+OFFSET)
        PRINT*, PARNAME(I-STRINX+1+1), NAME(I+OFFSET)
C        CSTAR(I-STRINX+1) =
C     &       10.0**(-.0337*MWT(I+OFFSET)+11.56-
C     &       SOMOXYGEN(I-STRINX+1)*DLVP(N))
        CSTAR(I-STRINX+1) =
     &       10.0**(-.0337*MWORG_HC(I-PRECSTRINX(1)+1)+11.56-
     &       SOMOXYGEN(I-STRINX+1)*DLVP(N))
        HVAP(I-STRINX+1)=-11.0*LOG10(CSTAR(I-STRINX+1))+131.0 ! HVAP - Epstein et al. 2010 ACP
       END DO
      END DO

      DO J=1,TNSOMSP
         PSTAR(J) = CSTAR(J)*R*T1/MWT(J+STRINX+OFFSET-1)/1.0e6 ! [pa]
C         PSTAR(J) = CSTAR(J)*R*T1/MWORG(J)/1.0e6 ! [pa]
C         print*, 'CSTAR=', CSTAR(j) , MWORG(j), pstar(j)
         !stop
      END DO

      OPEN(UNIT=37, FILE=SPECFILE, STATUS='old', ACCESS='APPEND')      
 904  FORMAT(A20,201F40.20)
      WRITE(37,904) 'CSTAR:', CSTAR(1:TNSOMSP)
      WRITE(37,904) 'HVAP:', HVAP(1:TNSOMSP)
      WRITE(37,904) 'PSTAR:', PSTAR(1:TNSOMSP)
      CLOSE(37)

      
      NSPECS = NACT !nact      
      DO I=1,6
         CONST(I) = CONC0(I)
      END DO

      CONST(6) = ( a_oh*exp(-1.0*ax_oh*0.0) +
     &             b_oh*exp(-1.0*bx_oh*0.0) ) *
     &           1.0E6/6.022E23/PRES*R*TEMP*1.0E6 ! OH Concentration - molecules/cm3 -> ppm

c     Copy rates of rxns to reactions array
      do I=1,MAXRXN
         RKZ(I) = RK(I)         
      enddo
c      print*, 'RKZ=', RKZ
c      print*, 'CONST=', CONST
      
CCC   COPY THE SOM CONCENTRATION FROM SAPRC-GC TO GC
CCC   *************************** CONVERT THE UNIT FROM ppm TO ug/m3 *********************************
      DO I=STRINX,NACT
         CALVALREAL = SAPRCGC(I) * 1.0E-6 * PRES/R/TEMP*
     &        MWT(I+OFFSET) * 1.0E-3 * BOXVOL * 1.0e-6 ! ppm -> g/m3 -> kg/bag
         GC(I-STRINX+SRTORG1) = DBLE(CALVALREAL)
      END DO

      print*, 'SP NAME, CARBON NO, OXYGEN NO, HYDROGEN NO, MW, CSTAR'
      do i=1,tnsomsp
         print*, NAME(i+STRINX+OFFSET-1), SOMCARBON(i),
     &        SOMOXYGEN(i), SOMHYDROGEN(I), MWT(i+STRINX+OFFSET-1),
     &        CSTAR(i)
      end do
      
C \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ AliA

      ! get organic aerosol psats at current temperature 
      ! modified cstar as the function of temp

      call boxinit(time,No,Dp,sigma,SRTSO4)

      ! read in background aerosol concentrations
      totamass=0.d0
      do n=1,ibins
         totamass=totamass+Nk(n)*sqrt(xk(n)*xk(n+1))
      enddo

      print*,'totamass [ug m-3]', totamass

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO N=1,IBINS
         IF (Nk(N) .LT. 1d-60) THEN
            Nk(N) = 0.d0
         END IF
         DO J=1,ICOMP
            IF (Mk(N,J) .LT. 1d-60) THEN
               Mk(N,J) = 0.d0
            END IF
         END DO
      END DO
      
      DO I=1,ICOMP-1
         IF (GC(I) .LT. 1d-60) THEN
            GC(I) = 0.d0
         END IF
      END DO
      DO I=1,NACT+OFFSET
         IF (SAPRCGC(I) .LT. 1d-60) THEN
            SAPRCGC(I) = 0.d0
         END IF
      END DO
      
      call report(time)
 999  FORMAT(306F40.20)
      OPEN(UNIT=38, FILE=SAPRCGCFILE, STATUS='old',ACCESS='APPEND')      
      WRITE(38,999) SAPRCGC(1:NACT)
      CLOSE(38)
      
      open(unit=31,file=vapname,status='old',access='append')
         write(31,2) time/3600.,(WalVapMas(n),P_temp(n+1), ! for g-p setting
     &           n=1,iorg)		 
      close(31)	  

      CTRW = 0

 10   continue  !beginning of new time step

      ! get organic aerosol psats at current temperature
      do j=1,iorg
         psatorg(j)=pstar(j)*exp(Hvap(j)*1.d3*(1.d0/T1-1.d0/temp)/R)
C         print*, 'psatorg', psatorg(j)
      enddo

C------------------------------------------------
C Now timestep through
      CTR = 0      
      do tstep=1,nintern
C C------------------------------------------------
	  

CCCCC Call Wall loss subroutine
C ----------------------------------------------------------------
         call wallloss(ke, kw0, VWL, PWL,
     &        Gc, Nk, Mk, aadt, time, P_temp, P_num,
     &        WalVapMas, tau_w, kw_HYBRID, aadt) !for normal case
	 
CCCCC add soa to gas phase
C ----------------------------------------------------------------         
         soarate = 0.0
         do j=1,iorg
            Gc(srtorg1+j-1)=Gc(srtorg1+j-1)
         enddo
         
         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 1'
               stop
            endif
            Nkout(k)=Nk(k)
            do j=1,icomp
               Mkout(k,j)=Mk(k,j)
            enddo
         enddo

CCCCC Coagulation
C ----------------------------------------------------------------         
         call storenm() 
         
         call multicoag(aadt)

         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 2'
               do n=1,ibins
                  print*,n,Nk(n),Nkout(n)
                  do j=1,icomp
                     print*,n,j,Mk(n,j),Mkout(n,j)
                  enddo
               enddo
               stop
            endif
         enddo
		 
CCCCC Do water eqm at appropriate times
C ----------------------------------------------------------------         
         call eznh3eqm(Gc,Mk)
         call ezwatereqm(Mk)

         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 3'
               stop
            endif
         enddo

CCCCC Nucleation
C ----------------------------------------------------------------         
         call storenm()
         Nk(1)=Nk(1)+nucrate*boxvol*aadt
         Mk(1,srtso4)=Mk(1,srtso4)+nucrate*boxvol*aadt*
     &                sqrt(xk(1)*xk(2))

         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 4'
               stop
            endif
         enddo
		 
         call ezwatereqm(Mk)
         call eznh3eqm(Gc,Mk)

         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 5'
               stop
            endif
         enddo

         call nucleation(Nk,Mk,Gc,Nkout,Mkout,Gcout,aadt,a_nuc,b_nuc,cstar_nuc)
         DO n=1,NBINS
            DO J=1,ICOMP
               Mk(n,j)=Mkout(n,j)
            ENDDO
            Nk(n)=Nkout(n)
         ENDDO
         Gc(srtso4)=Gcout(srtso4)
         Gc(srtnh4)=Gcout(srtnh4)

	 
CCCCC Sulfate condensation
C ----------------------------------------------------------------
         call storenm()
         call getCondSink(Nk,Mk,srtso4,CS,sinkfrac) 
         massin=0.0 !sulfconc/6.022D23*0.098d0*boxvol*CS*aadt
         call ezcond(Nk,Mk,massin,srtso4,Nkout,Mkout)
         tmass1=massin
         tmass2=0.d0
         do n=1,IBINS
            do j=1,icomp
               tmass1=tmass1+Mk(n,j)
               tmass2=tmass2+Mkout(n,j)
               Mk(n,j)=Mkout(n,j)
            enddo
            Nk(n)=Nkout(n)
         enddo
         if (abs(tmass1-tmass2)/massin.gt.1.E-3)then
            print*,'ERROR IN BOX.f'
            print*,'no conservation of mass in sulf cond'
            print*,tmass1,tmass2
            print*,massin
            stop
         endif

         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 6'
               stop
            endif
         enddo

CCCCC SOA condensation
C ----------------------------------------------------------------         
         call storenm()

         CTRW = CTRW + 1
         
CCCCC Gas phase chemistry
C ----------------------------------------------------------------         
         IF (CTRW.GT.1) THEN
            DO I=STRINX,NACT
               CALVALDBLE = GC(I-STRINX+SRTORG1) * 1.0E6 /
     &              PRES * R * TEMP / MWT(I+OFFSET)*1.0E3/BOXVOL*1.0e6 ! kg/bag -> g/m3 -> ppm
C               print*, I-STRINX+1, I, 'SAPRGC=',
C     &              SAPRCGC(I), 'CALVALDBLE=', CALVALDBLE
               SAPRCGC(I) = REAL(CALVALDBLE)
            END DO
         ELSE
            DO I=STRINX,NACT
               CALVALDBLE = GC(I-STRINX+SRTORG1) * 1.0E6 /
     &              PRES * R * TEMP / MWT(I+OFFSET)*1.0E3/BOXVOL*1.0e6 ! kg/bag -> g/m3 -> ppm
C               print*, I-STRINX+1, I, 'SAPRGC=',
C     &              SAPRCGC(I), 'CALVALDBLE=', CALVALDBLE
               SAPRCGC(I) = REAL(CALVALDBLE) + SAPRCGC(I)
            END DO
         END IF
         CTR = CTR + 1
         TOUT = time + tstep*aadt
         TIN = (TOUT - aadt)/60.0
         TOUT = TOUT/60.0
         NEQNT = NACT

         CONST(6) = ( a_oh*exp(-1.0*ax_oh*TIN/60.) +
     &                b_oh*exp(-1.0*bx_oh*TIN/60.) ) *
     &              1.0E6/6.022E23/PRES*R*TEMP*1.0E6 ! OH Concentration - molecules/cm3 -> ppm

         
         call INTEGR2(TIN,TOUT,NEQNT,SAPRCGC,DIFF,
     &        2.e-3, 1.0, 1.e-3, 0., 1.0E-30, 1.0e-09, 1, 5)
         
         DO I=STRINX,NACT
            CALVALREAL = SAPRCGC(I) * 1.0E-6 * PRES/R/TEMP*
     &           MWT(I+OFFSET) * 1.0E-3 * BOXVOL * 1.0e-6 ! ppm -> g/m3 -> kg/bag
            GC(I-STRINX+SRTORG1) = DBLE(CALVALREAL)
         END DO

CCCCC SOA condensation
C ----------------------------------------------------------------
         call soacond(Nk,Mk,Gc,Nkout,Mkout,Gcout,aadt,tau_p,Dbk,kc)
CAliA         call soacond(Nk,Mk,Gc,Nkout,Mkout,Gcout,aadt,tau_p)
         do n=1,IBINS
            do j=1,icomp
               Mk(n,j)=Mkout(n,j)
            enddo
            Nk(n)=Nkout(n)
         enddo

         do j=1,icomp-idiag
            Gc(j)=Gcout(j)
         enddo		
		 
  
         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 7'
               stop
            endif
         enddo

         call eznh3eqm(Gc,Mk)
         call ezwatereqm(Mk)
         call mnfix(Nk,Mk)
         
         do k=1,ibins
            if (isnan(Nk(k)))then
               print*,'NaN 8'
               stop
            endif
         enddo

         DO N=1,IBINS
            IF (Nk(N) .LT. 1d-60) THEN
               Nk(N) = 0.d0
            END IF
            DO J=1,ICOMP
               IF (Mk(N,J) .LT. 1d-60) THEN
                  Mk(N,J) = 0.d0
               END IF
            END DO
         END DO
            
         DO I=1,ICOMP-1
            IF (GC(I) .LT. 1d-60) THEN
               GC(I) = 0.d0
            END IF
         END DO
         DO I=1,NACT+OFFSET
            IF (SAPRCGC(I) .LT. 1d-60) THEN
               SAPRCGC(I) = 0.d0
            END IF
         END DO
         
      enddo     

C ***********************
C End of aerosol dynamics
C ***********************
      DO N=1,IBINS
         IF (Nk(N) .LT. 1d-60) THEN
            Nk(N) = 0.d0
         END IF
         DO J=1,ICOMP
            IF (Mk(N,J) .LT. 1d-60) THEN
               Mk(N,J) = 0.d0
            END IF
         END DO
      END DO

      DO I=1,ICOMP-1
         IF (GC(I) .LT. 1d-60) THEN
            GC(I) = 0.d0
         END IF
      END DO
      
C End of time step
      time=time+adt
      hour=time/3600.d0
      OPEN(UNIT=38, FILE=SAPRCGCFILE, STATUS='old',ACCESS='APPEND')      
      WRITE(38,999) SAPRCGC(1:NACT)
      CLOSE(38)

      call report(time)
		 
C----------------report vapor wall loss----------------	  
C----------------for g-p-p on the wall system-------------
C         do n=1,icomp
C   	     P_temp(n)=0.0
C         enddo
C   	  
C         do j=1,icomp
C         do n=1,ibins
C   	     P_temp(j)=P_temp(j)+P_temp1(n,j)
C         enddo
C         enddo	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 2    format(485E15.6)	
 2    format(710E15.6)
      open(unit=31,file=vapname,status='old',access='append')
      write(31,2) time/3600.,(WalVapMas(n),P_temp(n+1),n=1,iorg)

      close(31)
C----------------report taus for gp and gw partitioning---------------  

      open(unit=32, file=tauname, status='old',access='append')
      write(32,2) time/3600., tau_p,(tau_w(n),n=1,iorg)
      close(32)
C------------------report kw for each burns--------------------------------------		 
      open(unit=33, file=kwname, status='old',access='append')
      write(33,2) time/3600., (kw_HYBRID(n),n=1,ibins)
      close(33)

      if (hour .lt. endtime) goto 10
      totmass = 0.0
	    
      do n=1,ibins
         totmass=totmass+Nk(N)*sqrt(xk(n)*xk(n+1))
      enddo
      print*,'totmass [ug m-3]', totmass
      
      END   !of main


