c $Id: modlspc.h,v 1.1 2006/09/29 16:24:23 yingqi Exp yingqi $
c $Log: modlspc.h,v $
c Revision 1.1  2006/09/29 16:24:23  yingqi
c Initial revision
c
c Revision 1.1  2006/01/06 18:16:17  yingqi
c Initial revision
c
c
c  copyrighted (c) 1997 California Institute of Technology
c                  ALL RIGHTS RESERVED
c
C       MODLSPC.H:      SPECIFICATIONS CHEMICAL MECHANISM-DEPENDENT
C                       DATA.  TO BE USED WITH CHEMICAL MECHANISM
C                       MODULE OF ARB AIRSHED MODEL
C
C       CREATED BY:  WILLIAM P. L. CARTER  6/18/87.
C       LAST UPDATED: 11/16/88
C       modified by: Mike Kleeman (June 1995)
C       Modifications were made to remove the initialization and common
C       block of variables from this file and place this information
C       in the file "gaskin.h".  A similar file was creted to describe
C       an aqueous kinetic mechanism: "partkin.h".
C
C       SPECIFICATION OF ARRAY SIZES:  PARAMETERS GIVE MAXIMUM FOR
C       ARRAY DIMENSIONS.
C
C
C       VARIABLES AND PARAMETERS RELATING TO THE SPECIES, REACTIONS, AND
C       KINETIC PARAMETERS IN THE MECHANSIM.
C
C               TOTAL NUMBER OF SPECIES, AND NUMBERS OF ACTIVE, DUMMY,
C               UNREACTING PRODUCTS ("BUILDUP") AND CONSTANT SPECIES.
C

c changed by Qi Ying, Feburary 2003
c change to hold expanded gas phase chemistry

        INTEGER MAXNS, MAXCON, MAXACT, MAXBO, MAXDUM

c first row for sarpc99
c second row for sarpc07
        PARAMETER (MAXNS=2000,MAXCON=350,MAXACT=1800,MAXBO=10,MAXDUM=10)
c changed to hold expanded gas phase chemistry

C
C               NEW PARAMETER TO DEFINE MAXIMUM NUMBER OF STEADY STATE SPECIES
C
        INTEGER MAXSTS
        PARAMETER (MAXSTS=30)
C
C               NUMBER OF EMITTED RHC SPECIES
C
        INTEGER MAXRHC
        PARAMETER (MAXRHC=100)
C
C               NUMBER OF REACTIONS AND NUMBER OF KINETIC PARAMETERS
C               USED TO CALCULATE RATE CONSTANTS.
C
c changed to hold expanded chemistry
c modified by Qi Ying, Feburary 2003
        INTEGER MAXRXN, MAXKBF
        PARAMETER (MAXRXN=2000, MAXKBF=3*MAXRXN)
c modified by Qi Ying, Feburary 2003
C
C               STOICHIOMETRIC AND OTHER COEFFICIENTS.  "MAXCOV" CAN BE
C               VARIABLE (I.E., HAVE COEFFICIENT NAMES AND THUS BE CHANGED
C               AT INTEGRATION TIME), REST MUST BE CONSTANT.
C
        INTEGER MAXCOE, MAXCOV
        PARAMETER (MAXCOE=3500, MAXCOV=30)
C
C
C       PARAMETERS AND DATA USED IN SOLAR LIGHT INTENSITY CALCULATION.  FOR
C       ZENITH ANGLES BETWEEN 0 AND ~78 DEGREES, THE LIGHT INTENSITIES
C       ARE CALCULATED USING FUNCTIONS FIT TO PETERSON (EPA-600/4-76-025)S
C       ACTINIC IRRADIANCES FOR HIS "BEST ESTIMATE" SURFACE ALBITOS.  FOR
C       ZENITH ANGLES GREATER THAN ~78 DEGREES, LIGHT INTENSITIES ARE
C       ASSUMED TO BE PROPORTINAL TO COS(Z), DOWN TO A MINIMUM "DARK" VALUE
C       GIVEN BY CZMIN.  THE PARMAETER "CZLOW" IS COS(Z) FOR ~78 DEGREES.
C       THE ARRAY "DPHOTK" CONTAINS THE SOLAR LIGHT INTENSITIES FOR
C       COS(Z)=CZLOW.
C
        REAL CZLOW
        PARAMETER (CZLOW=7.0E-2)
C
C               NUMBER OF WAVELENGTHS FOR WHICH ABSORPTION COEFFICIENT
C               AND QUANTUM YIELD DATA ARE GIVEN, FOR THE PURPOSE OF
C               CALCULATING SOLAR PHOTOLYSIS RATES.  THE WAVELENGTHS
C               *MUST* BE AS INDICATED BELOW
C
        INTEGER NSOLWL
        PARAMETER (NSOLWL=27)
C       DATA SOLWL /     0.295,   0.300,    0.305,    0.310,    0.315,
C     &         0.320,   0.325,   0.330,    0.335,    0.340,    0.345,
C     &         0.350,            0.360,              0.370,
C     &         0.380,            0.390,              0.400,
C     &         0.410,            0.420,    0.430,    0.440,    0.450,
C     &         0.460,   0.470,   0.480,    0.490,    0.500,
C     &                  1.000/
C
C               NUMBER OF DIFFERENT PHOTOLYSIS REACTIONS
C
        INTEGER MAXPHK
        PARAMETER (MAXPHK=35)
