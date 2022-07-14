      SUBROUTINE INTEGR2(tin,tout,nlast,c,dc,
     &     e1,e2,e3,e4,e5,e6,itmin,itmax)
C************************************************************************
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
C
C     THIS SUBROUTINE INTEGRATES A SET OF COUPLED DIFFERENTIAL EQNS
C     USING A HYBRID INTEGRATION TECHNIQUE.
C
C     THE FORM OF THE DIFFERENTIAL EXPRESSIONS ASSUMED BY THIS PROGRAM
C     IS:
C            DC/DT = A(C,T) - B(C,T)*C = F(C,T)
C
C     A SUBROUTINE OF THE FORM:
C
C            DC(N, C, A, B, F)
C
C     MUST BE PROVIDED FOR THE SYSTEM OF INTEREST TO COMPUTE THE A, B AND
C     F ARRAYS FROM THE CONCENTRATIONS, C.
C
C************************************************************************
C
C     REFERENCE: T.R. YOUNG AND J.P. BORIS, "A NUMERICAL TECHNIQUE
C     FOR SOLVING STIFF ORDINARY DIFFERENTIAL EQUATIONS ASSOCIATED
C     WITH THE CHEMICAL KINETICS OF REACTIVE-FLOW PROBLEMS," J. PHYS.
C     CHEM., 81, 2424-2427, 1977.
C
C     REVISED BY:- MIKE KLEEMAN (NOVEMBER 1994)
C            ENVIRONMENTAL ENGINEERING SCIENCE 138-78
C              CALIFORNIA INSTITUTE OF TECHNOLOGY
C                PASADENA, CALIFORNIA, 91125
C
C     REVISED SO THAT:
C        1. THE CODE COULD BE USED TO SOLVE BOTH THE PHOTOCHEMICAL
C           REACTION SYSTEM AND THE GAS - PARTICLE CONVERSION SYSTEM.
C        2. THE CODE WOULD RUN MORE EFFICIENTLY.
C        3. AN ERROR IN THE TIME STEP MANIPULATION WAS CORRECTED.
C
C     REVISED BY:- ROBERT HARLEY (AUGUST 1989)
C            ENVIRONMENTAL ENGINEERING SCIENCE 138-78
C              CALIFORNIA INSTITUTE OF TECHNOLOGY
C                PASADENA, CALIFORNIA, 91125
C
C     REVISED TO MAKE SMALL EFFICIENCY IMPROVEMENTS, AND TO MAKE THE
C     CODE CLEARER AND MORE READABLE.  THE NUMERICAL CONSTANTS HAVE BEEN
C     SET BACK TO VALUES RECOMMENDED BY YOUNG & BORIS.
C
C************************************************************************
      IMPLICIT NONE

      INCLUDE 'modlspc.h'
      INCLUDE 'common.inc'
                                                                                
      REAL
     &     C(nd2),                 ! SPECIES CONCENTRATIONS
     &     C2(nd2), C3(nd2),       ! ITERATIVE PREDICTIONS FOR C(T+DT)
     &     A(nd2), A2(nd2),        ! FORMATION RATES
     &     B(nd2), B2(nd2),        ! DEPLETION "RATE CONSTANTS"
     &     F(nd2), F2(nd2)         ! NET FORMATION, F = A - B*C
      
      LOGICAL
     &     STIFF(nd2),             ! FLAG FOR STIFF EQUATIONS
     &     convg                   ! flag for convergence
      
      REAL
     &     TIN, TOUT, DTOVR2,      ! TIMES; START/CURRENT, STOP AND DT/2
     &     DT, DTMIN,              ! TIME STEP VARIABLES
     &     Y2, Y3, SUMTAU,         ! WORK VARIABLES
     &     E1, E2, E3, E4, E5, E6, ! INTEGRATION PARAMETERS
     &     reps,eps,               ! relative error and max relative error
     &     XXX, TIN1, DT1, DTK2,
     &     ESUM1, ESUM2,
     &     E4R, EPS1
                                                                                
      INTEGER
     &     J,k1,k2,k3,             ! LOOP CONTROL VARIABLES
     &     K4, NLAST,
     &     ITMIN, JFAIL,
     &     ITMAX                   ! MAX ITERATIONS BEFORE MODIFY DT
                                                                                
C***Control parameters for hybrid solution technique.************************
C         E1 - Scaling factor to determine initial time step.
C         E2 - Determines whether particular ODE falls into stiff
C              or normal class (normal => Euler's method for integration,
C              stiff => asymptotic method of Young & Boris).
C         E3 - Convergenece criterion.  Now specified as an argument
C              to INTEGR when the subroutine is called.  A typical
C              value would be 1.0E-3.
C         E4 - The minimum value of EPS used for time step modification
C              after a converged step.  This parameter should be set equal
C              to 0 if the equations are sensitive (this sets no limit) but
C              mariginal gains in efficiency can be realized for more stable
C              equations by using E4 = 0.95 or lower (this forces dt to
C              increase by some amount).
C         E5 - The lower bound for species concentrations (must be greater
C              than zero but sould be very small.
C         E6 - for species with low concentrations, we automatically
C              assume convergence.  E6 is the concentration below which
C              convergence is assumed.
C
C     ITMAX  - If the integration has not converged after ITMAX corrections
C              the time step is reduced.
C
C     ITMIN  - Must iterate at least itmin times to get optimal convegence.
                                                                                
C     *** REMEMBER THAT EVERY TIME DC IS CALLED, STEADY STATE / EQUILIBRIUM
C         CONCENTRATIONS WILL BE UPDATED.***
C
C***INITIALIZE VARIABLES***
      k2 = 0
      k3 = 0
      k4 = 0
      dtk2 = 0.
      e4r = e4
      dtmin = 1.0e-35
      convg = .true.
      do j = 1,nlast
       IF (C(J) .LT. e5) C(J)=e5
      enddo
C***COMPUTE INITIAL ESTIMATE FOR THE ARRAYS A, B, AND F***
      if(print4)print*,'calling diff with convg',convg,dt
      dt = 1.e-10

      CALL DC(nlast, C, A, B, F, tin, dt, convg)

      k3 = k3 + 1
      if(flag1)return
                                                                                
C***COMPUTE INITIAL ESTIMATE FOR THE TIME STEP***
      DT   =  10.0
      XXX  =  10.0
      DO 10 J = 1,nlast
       IF (F(J) .NE. 0.0)THEN
        XXX = ABS(C(J)/F(J))
       ENDIF
       if(b(j).gt.0.1) xxx = amin1(abs(1./b(j)),xxx)
       IF (XXX .GT. 0.0)  DT = amin1(DT, XXX)
 10   CONTINUE
      DT = E1*DT
      IF (DT .LT. DTMIN)  DT = DTMIN
                                                                                
c***BEGINNING OF THE PREDICTOR/CORRECTOR LOOP***
 20   IF (tin+dt.gt.tout)dt=tout-tin
      if(dt.lt.dtmin)then
       print*,'stepsize too small in integr',
     +      k1,k2,k3,jfail,eps,e4r,dt,tin,tout
       write(*,'(5(1x,e15.8))')eps,e4r,dt,tin,tout
       do j = 1,nlast
        print*,j,c(j),a(j),b(j),b2(j),stiff(j)
       enddo
       print*,'calling dc for diagnostic purposes:'
       print4 = .true.
       CALL DC(nlast,C, A, B, F, tin, dt, convg)
       stop 'integr: marker1'
      endif
      DTOVR2 = DT*0.5
      esum1 = 0.
      if(print4)then
        print*,'first step started at',tin,'with dt',dt,nlast
      endif
c--get the a and b values to be used in the next step--
      if(print4)print*,'calling diff with convg',convg,dt
      CALL DC(nlast,C, A, B, F, tin, dt, convg)
      if(print4)print*,'returned from diff with convg',convg,dt
      if(flag1)return
      if(print4)then
       print*,'first estimate of c,a,b based on converged conc'
       do j = 1,nlast
        write(*,*)j,c(j),a(j),b(j)
       enddo
      endif

C***DETERMINE STIFF EQUATIONS AND PREDICT CONCENTRATION AT T+DT***
      convg = .false.
      nstiff = 0
      DO 40 J = 1,nlast
       STIFF(J) = (DT*B(J)).GT.E2
       IF (STIFF(J)) THEN
        xxx = exp(-dt*b(j))
        c2(j) = a(j)/b(j)*(1-xxx) + c(j)*xxx
c        TAU2  = 2.0/B(J)
c        C2(J) = (C(J)*(TAU2-DT) + TAU2*DT*A(J)) / (TAU2+DT)
        nstiff = nstiff + 1
       ELSE
        C2(J) = C(J) + F(J)*DT
       ENDIF
c--reset negative concentrations--
       IF (C2(J) .LT. e5) C2(J)=e5
 40   CONTINUE
                                                                                
C***ITERATIVELY CORRECT THE PREDICTION FOR T+DT***
      k1 = 1                ! loop counter (< itmax)
c--compute a, b and f using the latest values at t+dt--
 50   if(print4)print*,'calling diff with convg',convg,dt
      flag4 = .false.
      CALL DC(nlast,C2, A2, B2, F2, tin, dt, convg)
      k3 = k3+1
      if(flag1)return
      esum2 = esum1

c--correct the predictions for c(t+dt)--
       DO 60 J = 1,nlast
        IF (.NOT.STIFF(J).OR.B2(J).LE.1.0E-30) THEN
          C3(J) = C(J) + (F(J)+F2(J))*DTOVR2
        ELSE
          SUMTAU = 1.0/B(J) + 1.0/B2(J)
          xxx = exp(-2.*dt/sumtau)
          c3(j) = (a(j) + a2(j))*sumtau/4.*(1.-xxx) + c(j)*xxx
        ENDIF
c--reset negative concentrations--
        IF (C3(J) .LT. e5) C3(J)=e5
 60   CONTINUE
                                                                                
c***TEST FOR CONVERGENCE AND CALCULATE EPS***
      eps = 0.
      esum1 = 0.
      jfail = 0
      do 70 j = 1,nlast
        Y2 = C2(J)
        Y3 = C3(J)
        if(y3.gt.e6 .or. y2.gt.e6)then
         reps = abs(y3-y2)/y3
         esum1 = esum1 + reps
         if(reps.gt.eps)then
          jfail = j
          eps = reps
         endif
        endif
 70   continue
c--calculate the relative error and remember for time step manipulation--
      eps = eps/e3
      esum1 = esum1/e3
c      write(6,*)'total sum of error',esum1
      if(k1.eq.1)eps1 = eps
      if(eps.gt.1.01)then
       convg = .false.
      else
       convg = .true.
      endif
                                                                                
c--check the total number of convergence checks and start printing info--
      k4 = k4 + 1
      IF (K4.GE.1.5E6)then
        print *,'setting print4 = true'
        PRINT4 = .TRUE.
      ENDIF
      IF (K4.GE.2.E6)STOP 'INTEGR2 TOO MANY ITERATIONS'
c**SYSTEM DID NOT CONVERGE OR FAILURE FLAG SET********************************
      if(flag4)then
       if(print4)then
        print*,'flag4 failure recognized with time step',dt
       endif
       eps1 = amax1(eps1,1.5)
       call step2(eps1,dt)
       if(print4)then
        print*,'failure pattern recognized, time step reduced to',dt
       endif
       k2 = 0
       goto 20
      endif
      if(.not. convg)then
       if(print4)then
        print*,'iteration',k1,'using time step',dt,'neq=',nlast
        write(*,'(1x,a7,i4,a20,2e15.8,l1)')'species',jfail,
     +  'did not converge',eps,esum1,flag1
        print*,c2(jfail),C3(jfail),f(jfail),f2(jfail),stiff(jfail)
        print*,a(jfail),a2(jfail),b(jfail),b2(jfail)
       endif
       k2 = 0
c--attempt to recognize failure pattern (reduce dt if pattern found)--
       if(esum2.gt.0.)then
        xxx = esum1/esum2
        if(xxx.gt.0.9)then
         call step2(eps1,dt)
         if(print4)then
          print*,'failure pattern recognized, time step reduced to',dt
         endif
         goto 20
        endif
       endif
c--no failure pattern, check iteration counters and decide what to do--
       if(k1.le.itmax)then
         do 90 j = 1,nlast
          c2(j) = C3(j)
 90      continue
         k1 = k1 + 1
         goto 50
       else
         call step2(eps1,dt)
         if(print4)then
          print*,'it>itmax so reduced time step to',dt
         endif
         goto 20
       endif
c**SYSTEM DID CONVERGE*******************************************************
      ELSE
c--check minimum iterations and update the concentrations--
       if(k3.lt.itmin)then
        do j = 1,nlast
         c2(j) = C3(j)
        enddo
        k1 = k1 + 1
c        if(print4)then
c         print*,'system converged but it<itmint',k3,itmin
c        endif
        goto 50
       else
        do j = 1,nlast
         c(j) = C3(j)
        enddo
       endif
c--check if time is equal to the target time--
       tin = tin+dt
       if(print4)then
         write(*,'(1x,a24,1x,e15.7,i3,i6,i3,e15.7)')
     +       'converged->dt,it1,k2,it2,t',dt,k1,k2,k3,tin
       endif
        k3 = 0
        dt1 = tout - tin
        tin1 = tin + dt1
       if(tout-tin.lt.dtmin .or. tin1.le.tin)then
c changed return condition
c       if (abs((tout-tin)/tout).lt.1.e-6) then
        print4 = .false.
        tin = tout
        return
       endif
c--consider what to do if we have no failure but dt constant--
        k2 = k2 + 1
        if(k2.eq.1)dtk2 = dt
        if(k2.gt.100)then
         if(dt.le.1.01*dtk2)then
          if(e4r.gt.0.1)then
           e4r = amin1(e4r,0.95)
          else
           e4r = 0.95
          endif
c          if(print4)print*,'setting e4r based on k2',k2,e4r
         endif
         dtk2 = dt
        endif
c--modify the stepsize and goto top of predictor/corrector loop--
        eps = amax1(0.2,eps)
        if(e4r.gt.0.1)eps = amin1(e4r,eps)
        call step2(eps,dt)
        GO TO 20
       endif
C
C-----END OF THE PROGRAM----
C
      END

      SUBROUTINE STEP2(Eps,Dt)
c****************************************************************************
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
c    created by: Mike Kleeman (November 1994)
c       Environmental Engineering Science 138-78
c       California Institute Of Technology
c       Pasadena California 91125
c
c The purpose of this subroutine is to estimate the appropriate time step
c for the integration of the equations so that eps->1.
c
c inputs
c       eps    max value of abs((c(i+1)-c(i))/c(i))
c       dt     value of previous timestep
c
c outputs
c       dt     selected time step for integration
c
c****************************************************************************
 
      REAL xx1 , xx2 , yyy , Eps , Dt
      INTEGER icount
 
c--modify the value of eps when it is just > 1 to improve convergence--
      IF ( Eps.GT.1.0 .AND. Eps.LE.1.1 ) THEN
         Eps = Eps**3
      ELSEIF ( Eps.GT.1.1 .AND. Eps.LT.1.2 ) THEN
         Eps = Eps**2
      ENDIF
 
c--calculate new time step using Newton iteration--
      icount = 1
      xx1 = 0.5*(Eps+1.0)
 100  xx2 = 0.5*(xx1+Eps/xx1)
      yyy = ABS(xx2-xx1)/xx2
      IF ( yyy.GT.0.01 .AND. icount.LE.5 ) THEN
         xx1 = xx2
         icount = icount + 1
         GOTO 100
      ENDIF
      xx2 = MIN(xx2,1.E3)
      Dt = Dt*(1./xx2)
      Dt = MIN(Dt,5.)
      END
 
