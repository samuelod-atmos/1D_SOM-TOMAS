c $Id: gaskin.h,v 1.1 2006/09/29 16:24:23 yingqi Exp yingqi $
c $Log: gaskin.h,v $
c Revision 1.1  2006/09/29 16:24:23  yingqi
c Initial revision
c
c Revision 1.1  2006/01/06 18:16:17  yingqi
c Initial revision
c
c*****************************************************************************
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
c  created by: Mike Kleeman (June 1995)
c         Environmental Engineering Science 138-78
c         California Institute of Technology
c         Pasedena, California, 91125
c
c  The purpose of this file is to define the common block containing the
c  information describing the gas phase kinetic model.
c
c  NOTE:  THIS FILE MUST BE PRECEDED BY THE SAPRC MECHANISMS INCLUDE FILE
c         "MODLSPC.H"
c
c*****************************************************************************
c
c  variables:
c    ns        number of species in the kinetic model
c    nact      number of active species in the kinetic model
c    nsc       number of constant species in the kinetic model
c    nsb       number of buildup species in the kinetic model
c    ndum      number of dummy species in the kinetic model
c
c For each of the next 5 arrays the first nsc are constant, next nsb are
c buildup, next ndum are dummy, next nact are active, rest are steady state.
c    name      array for species names
c    mwt       molecular weight of species in the model [g/mole]
c    cno       number of carbons in each species
c    conc0     initial concentrations of species in the model
c    nno       number of nitrogens in each species
c    sno       number of sulfurs in each species
c    xno       number of atom x in each species (user defined)
c
c    nrxn      number of reactions in the kinetic model
c    lenkbf    max number of parameters describing rate constants
c
c    rxnlbl    array of reaction labels
c    rk        array of rate constants
c    kpbuf     array of kinetic parameters for calculate of rate constants
c    lkbuf     location of parameters for reaction in array of kinetic params
c    irxtyp    kinetic parameter type index for reaction
c The index values and use of LKBUF are:
C               0 = SAME K AS ANOTHER REACTION.  LKBUF GIVES RXN NO.
C               1 = (NOT USED)
C               2 = (NOT USED)
C               3 = CONSTANT RATE CONSTANT.  T, LIGHT INDEPENDENT.
C                   LKBUF GIVES LOC. OF K IN KPBUF ARRAY.
C               4 = K CALCULATED USING ARRHENIUS EXPRESSION.  LKBUF
C                   GIVES LOC OF PARMS IN KPBUF ARRAY.  VALUES IN KPBUF
C                   ARE A, EA, AND B, WHERE
C                   K = A*EXP(EA/0.0019872*TEMP)*(TEMP/TREF)**B
C               5 = K CALCULATED USING FALLOFF EXPRESSION.  LKBUF GIVES
C                   LOC OF PARMS IN KPBUF ARRAY.  VALUES IN KPBUF ARE
C                   A0, E0, B0, AI, EI, BI, F AND N.
C               6 = (NOT USED)
C               7 = PHOTOLYSIS REACTION.  LKBUF GIVES INDEX FOR SET OF
C                   ABSORPTION COEFFICIENT/QUANTUM YIELDS.
C
c    ncoc1     index number for first constant coefficient in coef array
c    ncoc2     index number for last constant coefficient in coef array
c    ncoefv    number of variable coefficients in the model
c    coef      array of variable and constant coefficient values
c    nphotk    number of photolysis absorption profiles in the model
c    tref      reference temperature [K]
c    lat       lattitude used for calculation of solar intensity [ ]
c    dec       declination used for calculation of solar intensity [ ]
c    tstart    difference between time and true solar time [ ]
c    czmin     minimum dark value for which light intensity is proportional
c              to cos(z)
c    nrhc1     number of emitted reactive hydrocarbon species
c    indrhc    array of emitted rhc species index #s used to get species #
c              from rhc # used in emissions input
c    absqy     array of quantum yield data used to calculate photolysis rates
c    dphotk    array of "dark" photolysis rates
c    iphrk     index where rate constant is stored in rate constant array
c    coefnm    array for variable coefficient names
c    photnm    array of labels for photolysis reactions
c    title     model title
c
c--definition of variable types--
      INTEGER    NS, NSC, NACT, NSB, NDUM
      INTEGER    NNO(MAXNS), NRHC1,INDRHC(MAXRHC), NRXN, LENKBF
      INTEGER    SNO(MAXNS)
      INTEGER    LKBUF(MAXRXN), NPHOTK, IPHRK(MAXPHK)
      INTEGER    NCOC1, NCOC2, NCOEFV
      INTEGER    IRXTYP(MAXRXN)

      REAL*4     MWT(MAXNS), CNO(MAXNS), XNO(MAXNS), CONC0(MAXNS)
      REAL*4     ONO(MAXNS)
      REAL*4     KPBUF(MAXKBF), COEF(MAXCOE), RK(MAXRXN), RKZ(MAXRXN)
      REAL*4     DECLN, TSTART, CZMIN
      REAL*4     ABSQY(NSOLWL,MAXPHK), DPHOTK(MAXPHK), TREF

      CHARACTER*16 NAME(MAXNS)
      CHARACTER*6  RXNLBL(MAXRXN)
      CHARACTER*8  COEFNM(MAXCOV)
      CHARACTER*8  PHOTNM(MAXPHK)
      CHARACTER*80 TITLE


c--definition of common block for variables--
      COMMON /MODDAT/ NS, NACT, NSC, NSB, NDUM, NRXN, LENKBF, NCOC1,
     & NCOC2, NCOEFV, NPHOTK,
     & TREF, DECLN, TSTART, CZMIN,
     & MWT, CONC0, CNO, NNO, SNO, ONO, XNO, NRHC1, INDRHC,
     & RK, KPBUF, COEF, ABSQY, DPHOTK,
     & LKBUF, IPHRK, IRXTYP, RKZ

      COMMON /MODCHR/ NAME,RXNLBL, COEFNM, PHOTNM, TITLE
