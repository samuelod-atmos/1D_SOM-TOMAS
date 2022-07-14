
      SUBROUTINE DIFF(NEQN, C, A, B, F, TIN, DT, CONVG)
C
C     Program written by Armistead Russell
C     Revised version by Robert Harley (June 1992)
C     Revised for SAPRC mechanism compiler by Robert Harley (August 1992)
C     Revised to add elevated sources and vert adv by Mike Kleeman (March 2000)
C
C...PROGRAM VARIABLES...
C
C     A      - A COEFFICIENT FOR THE DIFFERENTIAL EQUATIONS
C     B      - B COEFFICIENT FOR THE DIFFERENTIAL EQUATIONS
C     C      - CONCENTRATION ARRAY AT TIME T
C     F      - DC/DT AT TIME T
C     TIN    - CURRENT TIME
C     DT     - CURRENT TIME STEP
C     CONVG  - LOGICAL FLAG TRUE IF LAST STEP CONVERGED
C
C...SUBPROGRAMS REQUIRED...
C
C     DIFUN - PROGRAM TO CALCULATE THE CHEMICAL REACTION RATES
C
      IMPLICIT NONE
      
      INCLUDE 'modlspc.h'
      INCLUDE 'gaskin.h'
      INCLUDE 'common.inc'
      
      INTEGER I, J
      INTEGER NEQN
      
      REAL    TIN, DT
      Real    C(NEQN), A(NEQN), B(NEQN), F(NEQN)
      Real    RXNRAT(maxrxn) !, PSSA(maxsts)

      logical convg
C
C     COMMON BLOCKS
C
CAliA      COMMON / PARAM / NSPECS, CONST(ncz)
CAliA      COMMON / TRAJT / EMT(LCOLSP)
C
C----------------------------------------------------------------------
C
C     STEP 1:- INITIALIZE THE DIFFERENTIAL TERM
C
C----------------------------------------------------------------------
C
C      print*, '2222222222222222222222222222222222222222222'
C      print*, 'CONST=', CONST
      DO 10 J=1, NEQN
        A(J) = 0.0
        B(J) = 0.0
 10   CONTINUE

c      print*, 'AAAAAAAAAAAAAAAAAAAA'
c      print*, A(1:NEQN)
c      print*, 'BBBBBBBBBBBBBBBBBBBB'
c      print*, B(1:NEQN)

C----------------------------------------------------------------------
C
C     STEP 2:- CALCULATE THE CONTRIBUTION DUE TO REACTION
C
C----------------------------------------------------------------------
C
c--  Write to file to check that information transferred
c        open(unit=1,file='rateofrxns.txt',status='unknown')
c        do ic = 1, NEQN
c           write(unit=1,*) ic, RKZ(ic), c(ic)
c        enddo        
C      print*, 'rhs.f before*******************************************'
c      print*, 'RKZ(1)=', RKZ(1)
c      print*, 'RKZ(2)=', RKZ(2)
c      print*, 'RKZ(3)=', RKZ(3)
c      print*, 'RKZ(6)=', RKZ(6)
c      print*, 'RKZ(7)=', RKZ(7)
c      print*, 'RKZ(8)=', RKZ(8)
c      print*, 'RKZ(18)=', RKZ(18)
c      print*, 'RKZ(19)=', RKZ(19)
c      print*, 'RKZ(30)=', RKZ(30)
c      print*, 'RKZ(36)=', RKZ(36)
c      print*, 'RKZ(188)=', RKZ(188)
c      print*, 'RKZ(191)=', RKZ(191)
c      print*, 'RKZ(196)=', RKZ(196)
c      print*, 'RKZ(200)=', RKZ(200)
c      print*, 'RKZ(256)=', RKZ(256)
c      print*, 'RKZ(260)=', RKZ(260)
c      print*, 'RKZ(264)=', RKZ(264)
c      print*, 'RKZ(268)=', RKZ(268)
c      print*, 'RKZ(277)=', RKZ(277)
c      print*, 'RKZ(281)=', RKZ(281)
c      print*, 'RKZ(291)=', RKZ(291)
c      print*, 'RKZ(303)=', RKZ(303)
c      print*, '=============================='
c      print*, 'SAPRCGC(1)=', C(1)
c      print*, 'SAPRCGC(2)=', C(2)
c      print*, 'SAPRCGC(3)=', C(3)
c      print*, 'SAPRCGC(4)=', C(4)
c      print*, 'SAPRCGC(35)=', C(35)
c      print*, 'SAPRCGC(37)=', C(37)
c      print*, 'SAPRCGC(39)=', C(39)
c      print*, 'SAPRCGC(41)=', C(41)
c      print*, 'SAPRCGC(67)=', C(67)
c      print*, 'SAPRCGC(74)=', C(74)
c      print*, 'SAPRCGC(78)=', C(78)
c      print*, 'SAPRCGC(79)=', C(79)
c      print*, 'SAPRCGC(84)=', C(84)
c      print*, 'SAPRCGC(85)=', C(85)
c      print*, 'SAPRCGC(90)=', C(90)
c      print*, 'SAPRCGC(91)=', C(91)
c      print*, 'SAPRCGC(93)=', C(93)
c      print*, '=============================='
c      print*, 'CONST(1)=', CONST(1)
c      print*, 'CONST(2)=', CONST(2)
c      print*, 'CONST(5)=', CONST(5)
c      print*, 'CONST(6)=', CONST(6)
c      print*, '=============================='
c      print*, 'COEF(1)=', COEF(1)
c      print*, '=============================='

      CALL DIFUN(CONST, C, PSSA, RKZ, RXNRAT,
     &     COEF, A, B)

C      print*, 'CONST=', CONST
C      print*, 'C=', C(1:NEQN)
C      print*, 'PSSA ******************************************'
C      print*, PSSA
c      print*, 'rhs.f after******************************************'
c      print*, 'RKZ(1)=', RKZ(1)
c      print*, 'RKZ(2)=', RKZ(2)
c      print*, 'RKZ(3)=', RKZ(3)
c      print*, 'RKZ(6)=', RKZ(6)
c      print*, 'RKZ(7)=', RKZ(7)
c      print*, 'RKZ(8)=', RKZ(8)
c      print*, 'RKZ(18)=', RKZ(18)
c      print*, 'RKZ(19)=', RKZ(19)
c      print*, 'RKZ(30)=', RKZ(30)
c      print*, 'RKZ(36)=', RKZ(36)
c      print*, 'RKZ(188)=', RKZ(188)
c      print*, 'RKZ(191)=', RKZ(191)
c      print*, 'RKZ(196)=', RKZ(196)
c      print*, 'RKZ(200)=', RKZ(200)
c      print*, 'RKZ(256)=', RKZ(256)
c      print*, 'RKZ(260)=', RKZ(260)
c      print*, 'RKZ(264)=', RKZ(264)
c      print*, 'RKZ(268)=', RKZ(268)
c      print*, 'RKZ(277)=', RKZ(277)
c      print*, 'RKZ(281)=', RKZ(281)
c      print*, 'RKZ(291)=', RKZ(291)
c      print*, 'RKZ(303)=', RKZ(303)
c      print*, '=============================='
c      print*, 'SAPRCGC(1)=', C(1)
c      print*, 'SAPRCGC(2)=', C(2)
c      print*, 'SAPRCGC(3)=', C(3)
c      print*, 'SAPRCGC(4)=', C(4)
c      print*, 'SAPRCGC(35)=', C(35)
c      print*, 'SAPRCGC(37)=', C(37)
c      print*, 'SAPRCGC(39)=', C(39)
c      print*, 'SAPRCGC(41)=', C(41)
c      print*, 'SAPRCGC(67)=', C(67)
c      print*, 'SAPRCGC(74)=', C(74)
c      print*, 'SAPRCGC(78)=', C(78)
c      print*, 'SAPRCGC(79)=', C(79)
c      print*, 'SAPRCGC(84)=', C(84)
c      print*, 'SAPRCGC(85)=', C(85)
c      print*, 'SAPRCGC(90)=', C(90)
c      print*, 'SAPRCGC(91)=', C(91)
c      print*, 'SAPRCGC(93)=', C(93)
c      print*, '=============================='
c      print*, 'CONST(1)=', CONST(1)
c      print*, 'CONST(2)=', CONST(2)
c      print*, 'CONST(5)=', CONST(5)
c      print*, 'CONST(6)=', CONST(6)
c      print*, '=============================='
c      print*, 'COEF(1)=', COEF(1)
c      print*, '=============================='
C      print*, 'RKZ ******************************************'
C      print*, RKZ
C      print*, 'RXNRAT=', RXNRAT
C      print*, 'COEF=', COEF
c      print*, 'A=', A(1:NEQN)
c      print*, 'B=', B(1:NEQN)
c      stop
C--------------------------------------------------------------------
C
C     PHASE 3 :- CALCULATE THE CONTRIBUTION FROM DIRECT SOURCES
C
C---------------------------------------------------------------------
C
CAliA      DO J=1, NEQN
CAliA         A(J) = A(J) + EMT(J)
CAliA      ENDDO
c
C
C-----------------------------------------------------------------------
C
C     STEP 3 :- compute the net rate of formation for each species
C
C-----------------------------------------------------------------------
C
      DO 200 I = 1, NEQN
        F(I) = A(I) - B(I)*C(I)
 200  CONTINUE
c      write(*,'(a12,5(1x,e8.2))'),'rates:',F(1),emt(1),a(1),b(1),c(1)
c
c
C-----END OF THE PROGRAM----
C
      RETURN
      END


