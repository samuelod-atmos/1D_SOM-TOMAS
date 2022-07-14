1**********************************************************************************************************************************

                              S A P R C   M O D E L   P R E P A R A T I O N   P R O G R A M

                                                     29-Nov-21 10:17 

                             TEMP 300.                                                       

 **********************************************************************************************************************************
0PARAMETERS SPECIFIED:
 

 REACTIONS INPUT:                                                              A        EA      B
                                                                               (IN PPM-MIN UNITS)
 1      NO2 + HV = NO + O3P                                                     PHOT. = NO2-06            3     8     0 =     2    51                            
 2      O3P + O2 + M = O3 + M                                              2.040E-05   0.00  -4.600      51     4     5 =     1     5                            
 3      O3P + O3 = #2 O2                                                   1.174E+04   4.09  -1.000      51     1     0 = -2001     4                            
 4      O3P + NO = NO2                                                    FALLOFF F= 0.600, N= 1.000     51     2     0 =     3                                  
                                                                      K0:  3.232E-03   0.00  -3.500
                                                                      KI:  4.404E+04   0.00  -1.000
 5      O3P + NO2 = NO + O2                                                8.073E+03  -0.37  -1.000      51     3     0 =     2     4                            
 6      O3P + NO2 = NO3                                                   FALLOFF F= 0.600, N= 1.000     51     3     0 =    10                                  
                                                                      K0:  8.978E-03   0.00  -3.800
                                                                      KI:  3.229E+04   0.00  -1.700
 7      O3 + NO = NO2 + O2                                                 4.404E+03   2.98  -1.000       1     2     0 =     3     4                            
 8      O3 + NO2 = O2 + NO3                                                2.055E+02   4.91  -1.000       1     3     0 =     4    10                            
 9      NO + NO3 = #2 NO2                                                  2.642E+04  -0.22  -1.000       2    10     0 = -2001     3                            
 10     NO + NO + O2 = #2 NO2                                              0.000E+00  -1.05  -2.000       2     2     4 = -2001     3                            
 11     NO2 + NO3 = N2O5                                                  FALLOFF F= 0.350, N= 1.330      3    10     0 =    11                                  
                                                                      K0:  1.293E-01   0.00  -6.100
                                                                      KI:  2.789E+03   0.00  -0.800
 12     N2O5 = NO2 + NO3                                                  FALLOFF F= 0.350, N= 1.330     11     0     0 =     3    10                            
                                                                      K0:  1.908E+12  21.86  -4.500
                                                                      KI:  5.820E+16  22.02   0.100
 13     N2O5 + H2O = #2 HNO3                                               3.670E-07   0.00  -1.000      11     6     0 = -2001    12                            
 14     N2O5 + H2O + H2O = #2 HNO3 + H2O                                   0.000E+00   0.00  -2.000      11     6     6 = -2001    12     6                      
 15     NO2 + NO3 = NO + NO2 + O2                                          6.606E+01   2.50  -1.000       3    10     0 =     2     3     4                      
 16     NO3 + HV = NO + O2                                                      PHOT. = NO3NO-06         10     8     0 =     2     4                            
 17     NO3 + HV = NO2 + O3P                                                    PHOT. = NO3NO2-6         10     8     0 =     3    51                            
 18     O3 + HV = O1D + O2                                                      PHOT. = O3O1D-06          1     8     0 =    52     4                            
 19     O3 + HV = O3P + O2                                                      PHOT. = O3O3P-06          1     8     0 =    51     4                            
 20     O1D + H2O = #2 OH                                                  2.393E+05  -0.12  -1.000      52     6     0 = -2001     9                            
 21     O1D + M = O3P + M                                                  3.494E+04  -0.19  -1.000      52     5     0 =    51     5                            
 22     OH + NO = HONO                                                    FALLOFF F= 0.600, N= 1.000      9     2     0 =    13                                  
                                                                      K0:  2.514E-02   0.00  -4.600
                                                                      KI:  5.284E+04   0.00  -1.100
 23     HONO + HV = OH + NO                                                     PHOT. = HONO-06          13     8     0 =     9     2                            
 24     OH + HONO = H2O + NO2                                              3.670E+03  -0.52  -1.000       9    13     0 =     6     3                            
 25     OH + NO2 = HNO3                                                   FALLOFF F= 0.600, N= 1.000      9     3     0 =    12                                  
                                                                      K0:  6.464E-02   0.00  -5.000
                                                                      KI:  4.110E+04   0.00  -1.000
 26     OH + NO3 = HO2 + NO2                                               2.936E+04   0.00  -1.000       9    10     0 =    88     3                            
 27     OH + HNO3 = H2O + NO3                                              3.523E+01  -0.91  -1.000       9    12     0 =     6    10                            
 28     HNO3 + HV = OH + NO2                                                    PHOT. = HNO3             12     8     0 =     9     3                            
 29     OH + CO = HO2 + CO2                                                2.114E+02   0.00  -1.000       9    89     0 =    88    90                            
 30     OH + O3 = HO2 + O2                                                 2.495E+03   1.87  -1.000       9     1     0 =    88     4                            
 31     HO2 + NO = OH + NO2                                                5.284E+03  -0.54  -1.000      88     2     0 =     9     3                            
 32     HO2 + NO2 = HNO4                                                  FALLOFF F= 0.600, N= 1.000     88     3     0 =    14                                  
                                                                      K0:  7.182E-03   0.00  -5.400
                                                                      KI:  4.257E+03   0.00  -2.100
 33     HNO4 = HO2 + NO2                                                  FALLOFF F= 0.600, N= 1.000     14     0     0 =    88     3                            
                                                                      K0:  5.461E+10  21.16  -3.400
                                                                      KI:  3.252E+17  22.20  -2.300
 34     HNO4 + HV = #.61 HO2 + #.61 NO2 + #.39 OH + #.39 NO3                    PHOT. = HNO4-06          14     8     0 = -2002    88 -2002     3                
                                                                                                                    -2003     9 -2003    10
 35     HNO4 + OH = H2O + NO2 + O2                                         1.908E+03  -0.75  -1.000      14     9     0 =     6     3     4                      
 36     HO2 + O3 = OH + #2 O2                                              2.980E-01  -1.38   3.570      88     1     0 =     9 -2001     4                      
 37     HO2 + HO2 = HO2H + O2                                              3.229E+02  -1.19  -1.000      88    88     0 =    91     4                            
 38     HO2 + HO2 + H2O = HO2H + O2 + H2O                                  1.106E-05  -5.56  -2.000      88    88     6 =    91     4     6                      
 39     NO3 + HO2 = #.8 OH + #.8 NO2 + #.8 O2 + #.2 HNO3 + #.2 O2          5.872E+03   0.00  -1.000      10    88     0 = -2004     9 -2004     3                
                                                                                                                    -2004     4 -2005    12
                                                                                                                    -2005     4
 40     NO3 + NO3 = #2 NO2 + O2                                            1.248E+03   4.87  -1.000      10    10     0 = -2001     3     4                      
 41     HO2H + HV = #2 OH                                                       PHOT. = H2O2             91     8     0 = -2001     9                            
 42     HO2H + OH = HO2 + H2O                                              2.642E+03   0.00  -1.000      91     9     0 =    88     6                            
 43     OH + HO2 = H2O + O2                                                7.046E+04  -0.50  -1.000       9    88     0 =     6     4                            
 44     OH + SO2 = HO2 + SULF + SO2                                       FALLOFF F= 0.600, N= 1.000      9    15     0 =    88    16    15                      
                                                                      K0:  1.185E-02   0.00  -6.300
                                                                      KI:  2.349E+03   0.00  -1.000
 45     OH + H2 = HO2 + H2O                                                1.130E+04   4.17  -1.000       9     7     0 =    88     6                            
 BR01   MEO2 + NO = NO2 + HCHO + HO2                                       3.376E+03  -0.71  -1.000      92     2     0 =     3    93    88                      
 BR02   MEO2 + HO2 = COOH + O2                                             5.079E+02  -1.55  -0.640      92    88     0 =    94     4                            
 BR03   MEO2 + HO2 = HCHO + O2 + H2O                                       4.903E+01  -1.55  -4.530      92    88     0 =    93     4     6                      
 BR04   MEO2 + NO3 = HCHO + HO2 + NO2                                      1.908E+03   0.00  -1.000      92    10     0 =    93    88     3                      
 BR05   MEO2 + MEO2 = MEOH + HCHO + O2                                     9.380E+01  -0.73  -2.800      92    92     0 =    95    93     4                      
 BR06   MEO2 + MEO2 = #2 HCHO + #2 HO2                                     1.086E+03   1.03  -1.000      92    92     0 = -2001    93 -2001    88                
 BR07   RO2C + NO = NO2                                                    3.817E+03  -0.75  -1.000      96     2     0 =     3                                  
 BR08   RO2C + HO2 =                                                       5.578E+02  -1.79  -1.000      96    88     0 =                                        
 BR09   RO2C + NO3 = NO2                                                   3.376E+03   0.00  -1.000      96    10     0 =     3                                  
 BR10   RO2C + MEO2 = #.5 HO2 + #.75 HCHO + #.25 MEOH                      2.936E+02   0.00  -1.000      96    92     0 = -2006    88 -2007    93                
                                                                                                                    -2008    95
 BR11   RO2C + RO2C =                                                      5.138E+01   0.00  -1.000      96    96     0 =                                        
 BR12   RO2XC + NO = XN                                                         SAME K AS BR07           97     2     0 =    18                                  
 BR13   RO2XC + HO2 =                                                           SAME K AS BR08           97    88     0 =                                        
 BR14   RO2XC + NO3 = NO2                                                       SAME K AS BR09           97    10     0 =     3                                  
 BR15   RO2XC + MEO2 = #.5 HO2 + #.75 HCHO + #.25 MEOH                          SAME K AS BR10           97    92     0 = -2006    88 -2007    93                
                                                                                                                    -2008    95
 BR16   RO2XC + RO2C =                                                          SAME K AS BR11           97    96     0 =                                        
 BR17   RO2XC + RO2XC =                                                         SAME K AS BR11           97    97     0 =                                        
 BR18   MECO3 + NO2 = PAN                                                 FALLOFF F= 0.300, N= 1.410     98     3     0 =    20                                  
                                                                      K0:  9.696E+00   0.00  -9.100
                                                                      KI:  1.776E+04   0.00  -1.900
 BR19   PAN = MECO3 + NO2                                                 FALLOFF F= 0.300, N= 1.410     20     0     0 =    98     3                            
                                                                      K0:  7.193E+12  24.05  -1.000
                                                                      KI:  2.400E+18  27.03   0.000
 BR20   PAN + HV = #.6 MECO3 + #.6 NO2 + #.4 MEO2 + &                            
        #.4 CO2 + #.4 NO3                                                       PHOT. = PAN              20     8     0 = -2009    98 -2009     3                
                                                                                                                    -2010    92 -2010    90
                                                                                                                    -2010    10
 BR21   MECO3 + NO = MEO2 + CO2 + NO2                                      1.101E+04  -0.58  -1.000      98     2     0 =    92    90     3                      
 BR22   MECO3 + HO2 = #.44 OH + #.44 MEO2 + #.44 CO2 + &                         
        #.41 CCOOH + #.15 O3 + #.15 CCOOH                                  7.633E+02  -1.95  -1.000      98    88     0 = -2011     9 -2011    92                
                                                                                                                    -2011    90 -2012    99
                                                                                                                    -2013     1 -2013    99
 BR23   MECO3 + NO3 = MEO2 + CO2 + NO2 + O2                                     SAME K AS BR09           98    10     0 =    92    90     3     4                
 BR24   MECO3 + MEO2 = #.1 CCOOH + #.1 HCHO + #.1 O2 + &                         
        #.9 HCHO + #.9 HO2 + #.9 MEO2 + #.9 CO2                            2.936E+03  -0.99  -1.000      98    92     0 = -2014    99 -2014    93                
                                                                                                                    -2014     4 -2015    93
                                                                                                                    -2015    88 -2015    92
                                                                                                                    -2015    90
 BR25   MECO3 + RO2C = MEO2 + CO2                                          6.459E+02  -2.13  -1.000      98    96     0 =    92    90                            
 BR26   MECO3 + RO2XC = MEO2 + CO2                                              SAME K AS BR25           98    97     0 =    92    90                            
 BR27   MECO3 + MECO3 = #2 MEO2 + #2 CO2 + O2                              4.257E+03  -0.99  -1.000      98    98     0 = -2001    92 -2001    90                
                                                                                                                        4
 BR28   RCO3 + NO2 = PAN2                                                  1.776E+04   0.00  -2.070     100     3     0 =    21                                  
 BR29   PAN2 = RCO3 + NO2                                                  4.980E+18  27.70   0.000      21     0     0 =   100     3                            
 BR30   PAN2 + HV = #.6 RCO3 + #.6 NO2 + #.4 RO2C + &                            
        #.4 xHO2 + #.4 yROOH + #.4 xCCHO + #.4 CO2 + &                           
        #.4 N3                                                                  PHOT. = PAN              21     8     0 = -2009   100 -2009     3                
                                                                                                                    -2010    96 -2010    56
                                                                                                                    -2010    82 -2010    66
                                                                                                                    -2010    90 -2010   101
 BR31   RCO3 + NO = NO2 + RO2C + xHO2 + yROOH + xCCHO + &                        
        CO2                                                                9.835E+03  -0.68  -1.000     100     2     0 =     3    96    56    82                
                                                                                                                       66    90
 BR32   RCO3 + HO2 = #.44 OH + #.44 RO2C + #.44 xHO2 + &                         
        #.44 xCCHO + #.44 yROOH + #.44 CO2 + &                                   
        #.41 RCOOH + #.15 O3 + #.15 RCOOH                                       SAME K AS BR22          100    88     0 = -2011     9 -2011    96                
                                                                                                                    -2011    56 -2011    66
                                                                                                                    -2011    82 -2011    90
                                                                                                                    -2012   102 -2013     1
                                                                                                                    -2013   102
 BR33   RCO3 + NO3 = NO2 + RO2C + xHO2 + yROOH + xCCHO + &                       
        CO2 + O2                                                                SAME K AS BR09          100    10     0 =     3    96    56    82                
                                                                                                                       66    90     4
 BR34   RCO3 + MEO2 = HCHO + HO2 + RO2C + xHO2 + xCCHO + &                       
        yROOH + CO2                                                             SAME K AS BR24          100    92     0 =    93    88    96    56                
                                                                                                                       66    82    90
 BR35   RCO3 + RO2C = RO2C + xHO2 + xCCHO + yROOH + CO2                         SAME K AS BR25          100    96     0 =    96    56    66    82                
                                                                                                                       90
 BR36   RCO3 + RO2XC = RO2C + xHO2 + xCCHO + yROOH + CO2                        SAME K AS BR25          100    97     0 =    96    56    66    82                
                                                                                                                       90
 BR37   RCO3 + MECO3 = #2 CO2 + MEO2 + RO2C + xHO2 + &                           
        yROOH + xCCHO + O2                                                      SAME K AS BR27          100    98     0 = -2001    90    92    96                
                                                                                                                       56    82    66     4
 BR38   RCO3 + RCO3 = #2 RO2C + #2 xHO2 + #2 xCCHO + &                           
        #2 yROOH + #2 CO2                                                       SAME K AS BR27          100   100     0 = -2001    96 -2001    56                
                                                                                                                    -2001    66 -2001    82
                                                                                                                    -2001    90
 BR39   BZCO3 + NO2 = PBZN                                                 2.011E+04   0.00  -1.000     103     3     0 =    22                                  
 BR40   PBZN = BZCO3 + NO2                                                 4.740E+18  27.82   0.000      22     0     0 =   103     3                            
 BR41   PBZN + HV = #.6 BZCO3 + #.6 NO2 + #.4 CO2 + &                            
        #.4 BZO + #.4 RO2C + #.4 NO3                                            PHOT. = PAN              22     8     0 = -2009   103 -2009     3                
                                                                                                                    -2010    90 -2010    54
                                                                                                                    -2010    96 -2010    10
 BR42   BZCO3 + NO = NO2 + CO2 + BZO + RO2C                                     SAME K AS BR31          103     2     0 =     3    90    54    96                
 BR43   BZCO3 + HO2 = #.44 OH + #.44 BZO + #.44 RO2C + &                         
        #.44 CO2 + #.41 RCOOH + &                                                
        #.15 O3 + #.44 RCOOH + #2.24 XC                                         SAME K AS BR22          103    88     0 = -2011     9 -2011    54                
                                                                                                                    -2011    96 -2011    90
                                                                                                                    -2012   102 -2013     1
                                                                                                                    -2011   102 -2016   104
 BR44   BZCO3 + NO3 = NO2 + CO2 + BZO + RO2C + O2                               SAME K AS BR09          103    10     0 =     3    90    54    96                
                                                                                                                        4
 BR45   BZCO3 + MEO2 = HCHO + HO2 + RO2C + BZO + CO2                            SAME K AS BR24          103    92     0 =    93    88    96    54                
                                                                                                                       90
 BR46   BZCO3 + RO2C = RO2C + BZO + CO2                                         SAME K AS BR25          103    96     0 =    96    54    90                      
 BR47   BZCO3 + RO2XC = RO2C + BZO + CO2                                        SAME K AS BR25          103    97     0 =    96    54    90                      
 BR48   BZCO3 + MECO3 = #2 CO2 + MEO2 + BZO + RO2C                              SAME K AS BR27          103    98     0 = -2001    90    92    54                
                                                                                                                       96
 BR49   BZCO3 + RCO3 = #2 CO2 + RO2C + xHO2 + yROOH + &                          
        xCCHO + BZO + RO2C                                                      SAME K AS BR27          103   100     0 = -2001    90    96    56                
                                                                                                                       82    66    54    96
 BR50   BZCO3 + BZCO3 = #2 BZO + #2 RO2C + #2 CO2                               SAME K AS BR27          103   103     0 = -2001    54 -2001    96                
                                                                                                                    -2001    90
 BR51   MACO3 + NO2 = MAPAN                                                     SAME K AS BR28          105     3     0 =    23                                  
 BR52   MAPAN = MACO3 + NO2                                                9.600E+17  26.80   0.000      23     0     0 =   105     3                            
 BR53   MAPAN + HV = #.6 MACO3 + #.6 NO2 + #.4 CO2 + &                           
        #.4 HCHO + #.4 MECO3 + #.4 NO3                                          PHOT. = PAN              23     8     0 = -2009   105 -2009     3                
                                                                                                                    -2010    90 -2010    93
                                                                                                                    -2010    98 -2010    10
 BR54   MACO3 + NO = NO2 + CO2 + HCHO + MECO3                                   SAME K AS BR31          105     2     0 =     3    90    93    98                
 BR55   MACO3 + HO2 = #.44 OH + #.44 HCHO + #.44 MECO3 + &                       
        #.44 CO2 + #.41 RCOOH + &                                                
        #.15 O3 + #.15 RCOOH + #.56 XC                                          SAME K AS BR22          105    88     0 = -2011     9 -2011    93                
                                                                                                                    -2011    98 -2011    90
                                                                                                                    -2012   102 -2013     1
                                                                                                                    -2013   102 -2017   104
 BR56   MACO3 + NO3 = NO2 + CO2 + HCHO + MECO3 + O2                             SAME K AS BR09          105    10     0 =     3    90    93    98                
                                                                                                                        4
 BR57   MACO3 + MEO2 = #2 HCHO + HO2 + CO2 + MECO3                              SAME K AS BR24          105    92     0 = -2001    93    88    90                
                                                                                                                       98
 BR58   MACO3 + RO2C = CO2 + HCHO + MECO3                                       SAME K AS BR25          105    96     0 =    90    93    98                      
 BR59   MACO3 + RO2XC = CO2 + HCHO + MECO3                                      SAME K AS BR25          105    97     0 =    90    93    98                      
 BR60   MACO3 + MECO3 = #2 CO2 + MEO2 + HCHO + MECO3 + O2                       SAME K AS BR27          105    98     0 = -2001    90    92    93                
                                                                                                                       98     4
 BR61   MACO3 + RCO3 = HCHO + MECO3 + RO2C + xHO2 + &                            
        yROOH + xCCHO + #2 CO2                                                  SAME K AS BR27          105   100     0 =    93    98    96    56                
                                                                                                                       82    66 -2001    90
 BR62   MACO3 + BZCO3 = HCHO + MECO3 + BZO + RO2C + &                            
        #2 CO2                                                                  SAME K AS BR27          105   103     0 =    93    98    54    96                
                                                                                                                    -2001    90
 BR63   MACO3 + MACO3 = #2 HCHO + #2 MECO3 + #2 CO2                             SAME K AS BR27          105   105     0 = -2001    93 -2001    98                
                                                                                                                    -2001    90
 BR64   TBUO + NO2 = RNO3 + #-1 XC + #-1 XC                                3.523E+04   0.00  -1.000      53     3     0 =    19 -2018   104 -2018                
                                                                                                                      104
 BR65   TBUO = ACET + MEO2                                                 4.500E+16  16.20   0.000      53     0     0 =   106    92                            
 BR66   BZO + NO2 = NPHE                                                   3.376E+04  -0.30  -1.000      54     3     0 =    17                                  
 BR67   BZO + HO2 = CRES + #-1 XC                                               SAME K AS BR08           54    88     0 =   107 -2018   104                      
 BR68   BZO = CRES + RO2C + xHO2 + #-1 XC                                  6.000E-02   0.00   0.000      54     0     0 =   107    96    56 -2018                
                                                                                                                      104
 RO01   xHO2 + #RO2RO = HO2                                                1.000E+00   0.00   0.000      56    -5     0 =    88                                  
 RO02   xHO2 + #RO2XRO =                                                   1.000E+00   0.00   0.000      56    -6     0 =                                        
 RO03   xOH + #RO2RO = OH                                                  1.000E+00   0.00   0.000      57    -5     0 =     9                                  
 RO04   xOH + #RO2XRO =                                                    1.000E+00   0.00   0.000      57    -6     0 =                                        
 RO05   xNO2 + #RO2RO = NO2                                                1.000E+00   0.00   0.000      58    -5     0 =     3                                  
 RO06   xNO2 + #RO2XRO = XN                                                1.000E+00   0.00   0.000      58    -6     0 =    18                                  
 RO07   xMEO2 + #RO2RO = MEO2                                              1.000E+00   0.00   0.000      59    -5     0 =    92                                  
 RO08   xMEO2 + #RO2XRO = XC                                               1.000E+00   0.00   0.000      59    -6     0 =   104                                  
 RO09   xMECO3 + #RO2RO = MECO3                                            1.000E+00   0.00   0.000      60    -5     0 =    98                                  
 RO10   xMECO3 + #RO2XRO = #2 XC                                           1.000E+00   0.00   0.000      60    -6     0 = -2001   104                            
 RO11   xRCO3 + #RO2RO = RCO3                                              1.000E+00   0.00   0.000      61    -5     0 =   100                                  
 RO12   xRCO3 + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      61    -6     0 = -2019   104                            
 RO13   xMACO3 + #RO2RO = MACO3                                            1.000E+00   0.00   0.000      62    -5     0 =   105                                  
 RO14   xMACO3 + #RO2XRO = #4 XC                                           1.000E+00   0.00   0.000      62    -6     0 = -2020   104                            
 RO15   xTBUO + #RO2RO = TBUO                                              1.000E+00   0.00   0.000      63    -5     0 =    53                                  
 RO16   xTBUO + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      63    -6     0 = -2020   104                            
 RO17   xCO + #RO2RO = CO                                                  1.000E+00   0.00   0.000      64    -5     0 =    89                                  
 RO18   xCO + #RO2XRO = XC                                                 1.000E+00   0.00   0.000      64    -6     0 =   104                                  
 BP01   HCHO + HV = #2 HO2 + CO                                                 PHOT. = HCHOR-06         93     8     0 = -2001    88    89                      
 BP02   HCHO + HV = H2 + CO                                                     PHOT. = HCHOM-06         93     8     0 =     7    89                            
 BP03   HCHO + OH = HO2 + CO + H2O                                         7.927E+03  -0.27  -1.000      93     9     0 =    88    89     6                      
 BP07   HCHO + NO3 = HNO3 + HO2 + CO                                       2.936E+03   4.83  -1.000      93    10     0 =    12    88    89                      
 BP08   CCHO + OH = MECO3 + H2O                                            6.459E+03  -0.73  -1.000     108     9     0 =    98     6                            
 BP09   CCHO + HV = CO + HO2 + MEO2                                             PHOT. = CCHO_R          108     8     0 =    89    88    92                      
 BP10   CCHO + NO3 = HNO3 + MECO3                                          2.055E+03   3.70  -1.000     108    10     0 =    12    98                            
 BP11   RCHO + OH = #.965 RCO3 + #.035 RO2C + &                                  
        #.035 xHO2 + #.035 xCO + #.035 xCCHO + #.035 yROOH                 7.486E+03  -0.81  -1.000     109     9     0 = -2021   100 -2022    96                
                                                                                                                    -2022    56 -2022    64
                                                                                                                    -2022    66 -2022    82
 BP12   RCHO + HV = RO2C + xHO2 + yROOH + xCCHO + CO + &                         
        HO2                                                                     PHOT. = C2CHO           109     8     0 =    96    56    82    66                
                                                                                                                       89    88
 BP13   RCHO + NO3 = HNO3 + RCO3                                           2.055E+03   3.18  -1.000     109    10     0 =    12   100                            
 BP14   ACET + OH = RO2C + xMECO3 + xHCHO + yROOH                          6.694E+01  -0.85   2.650     106     9     0 =    96    60    65    82                
 BP15   ACET + HV + #QY.ACET = #.62 MECO3 + #1.38 MEO2 + #.38 CO                PHOT. = ACET-06         106     8    -7 = -2023    98 -2024    92                
                                                                                                                    -2025    89
 BP16   MEK + OH = #.967 RO2C + #.039 RO2XC + #.039 zRNO3 + &                    
        #.376 xHO2 + #.51 xMECO3 + #.074 xRCO3 + &                               
        #.088 xHCHO + #.504 xCCHO + #.376 xRCHO + yROOH + &                      
        #.3 XC                                                             1.908E+03   0.05   1.000     110     9     0 = -2026    96 -2027    97                
                                                                                                                    -2027    81 -2028    56
                                                                                                                    -2029    60 -2030    61
                                                                                                                    -2031    65 -2032    66
                                                                                                                    -2028    67    82 -2033
                                                                                                                      104
 BP17   MEK + HV + #QY.MEK = MECO3 + RO2C + xHO2 + xCCHO + yROOH                PHOT. = MEK-06          110     8    -8 =    98    96    56    66                
                                                                                                                       82
 BP18   MEOH + OH = HCHO + HO2                                             4.184E+03   0.69  -1.000      95     9     0 =    93    88                            
 BP19   HCOOH + OH = HO2 + CO2                                             6.606E+02   0.00  -1.000     111     9     0 =    88    90                            
 BP20   CCOOH + OH = #.509 MEO2 + #.491 RO2C + #.509 CO2 + &                     
        #.491 xHO2 + #.491 xMGLY + #.491 yROOH + &                               
        #-0.491 XC                                                         6.165E+01  -1.70  -1.000      99     9     0 = -2034    92 -2035    96                
                                                                                                                    -2034    90 -2035    56
                                                                                                                    -2035    72 -2035    82
                                                                                                                    -2036   104
 BP21   RCOOH + OH = RO2C + xHO2 + #.143 CO2 + &                                 
        #.142 xCCHO + #.4 xRCHO + #.457 xBACL + yROOH + &                        
        #-0.455 XC                                                         1.761E+03   0.00  -1.000     102     9     0 =    96    56 -2037    90                
                                                                                                                    -2038    66 -2010    67
                                                                                                                    -2039    73    82 -2040
                                                                                                                      104
 BP22   COOH + OH = H2O + #.3 HCHO + #.3 OH + #.7 MEO2                     5.578E+03  -0.40  -1.000      94     9     0 =     6 -2033    93 -2033                
                                                                                                                        9 -2041    92
 BP23   COOH + HV = HCHO + HO2 + OH                                             PHOT. = COOH             94     8     0 =    93    88     9                      
 BP24   ROOH + OH = #.744 OH + #.251 RO2C + #.004 RO2XC + &                      
        #.004 zRNO3 + #.744 RCHO + #.239 xHO2 + &                                
        #.012 xOH + #.012 xHCHO + #.012 xCCHO + &                                
        #.205 xRCHO + #.034 xPROD2 + #.256 yROOH + &                             
        #-0.115 XC                                                         3.670E+04   0.00  -1.000     112     9     0 = -2042     9 -2043    96                
                                                                                                                    -2044    97 -2044    81
                                                                                                                    -2042   109 -2045    56
                                                                                                                    -2046    57 -2046    65
                                                                                                                    -2046    66 -2047    67
                                                                                                                    -2048    70 -2049    82
                                                                                                                    -2050   104
 BP25   ROOH + HV = RCHO + HO2 + OH                                             PHOT. = COOH            112     8     0 =   109    88     9                      
 BP26   R6OOH + OH = #.84 OH + #.222 RO2C + #.029 RO2XC + &                      
        #.029 zRNO3 + #.84 PROD2 + #.09 xHO2 + #.041 xOH + &                     
        #.02 xCCHO + #.075 xRCHO + #.084 xPROD2 + &                              
        #.16 yROOH + #.02 XC                                               8.220E+04   0.00  -1.000     113     9     0 = -2051     9 -2052    96                
                                                                                                                    -2053    97 -2053    81
                                                                                                                    -2051   114 -2054    56
                                                                                                                    -2055    57 -2056    66
                                                                                                                    -2057    67 -2058    70
                                                                                                                    -2059    82 -2056   104
 BP27   R6OOH + HV = OH + #.142 HO2 + #.782 RO2C + &                             
        #.077 RO2XC + #.077 zRNO3 + #.085 RCHO + &                               
        #.142 PROD2 + #.782 xHO2 + #.026 xCCHO + &                               
        #.058 xRCHO + #.698 xPROD2 + #.858 yR6OOH + &                            
        #.017 XC                                                                PHOT. = COOH            113     8     0 =     9 -2038    88 -2060                
                                                                                                                       96 -2061    97 -2061
                                                                                                                       81 -2062   109 -2038
                                                                                                                      114 -2060    56 -2063
                                                                                                                       66 -2064    67 -2065
                                                                                                                       70 -2066    83 -2067
                                                                                                                      104
 BP28   RAOOH + OH = #.139 OH + #.148 HO2 + #.589 RO2C + &                       
        #.124 RO2XC + #.124 zRNO3 + #.074 PROD2 + &                              
        #.147 MGLY + #.139 IPRD + #.565 xHO2 + #.024 xOH + &                     
        #.448 xRCHO + #.026 xGLY + #.030 xMEK + &                                
        #.252 xMGLY + #.073 xAFG1 + #.073 xAFG2 + &                              
        #.713 yR6OOH + #2.674 XC                                           2.070E+05   0.00  -1.000     115     9     0 = -2068     9 -2069    88                
                                                                                                                    -2070    96 -2071    97
                                                                                                                    -2071    81 -2030   114
                                                                                                                    -2072   116 -2068   117
                                                                                                                    -2073    56 -2074    57
                                                                                                                    -2075    67 -2063    71
                                                                                                                    -2076    69 -2077    72
                                                                                                                    -2078    75 -2078    76
                                                                                                                    -2079    83 -2080   104
 BP29   RAOOH + HV = OH + HO2 + #.5 GLY + #.5 MGLY + &                           
        #.5 AFG1 + #.5 AFG2 + &                                                  
        #.5 XC                                                                  PHOT. = COOH            115     8     0 =     9    88 -2006   118                
                                                                                                                    -2006   116 -2006   119
                                                                                                                    -2006   120 -2006   104
 BP30   GLY + HV = #2 CO + #2 HO2                                               PHOT. = GLY-07R         118     8     0 = -2001    89 -2001    88                
 BP31   GLY + HV = HCHO + CO                                                    PHOT. = GLY-07M         118     8     0 =    93    89                            
 BP32   GLY + OH = #.7 HO2 + #1.4 CO + #.3 HCOCO3                          4.550E+03  -0.68  -1.000     118     9     0 = -2041    88 -2081    89                
                                                                                                                    -2033   121
 BP33   GLY + NO3 = HNO3 + #.7 HO2 + #1.4 CO + #.3 HCOCO3                  4.110E+03   4.72  -1.000     118    10     0 =    12 -2041    88 -2081                
                                                                                                                       89 -2033   121
 BP80   HCOCO3 + NO = HO2 + CO + CO2 + NO2                                      SAME K AS BR31          121     2     0 =    88    89    90     3                
 BP81   HCOCO3 + NO2 = HO2 + CO + CO2 + NO3                                     SAME K AS BR28          121     3     0 =    88    89    90    10                
 BP82   HCOCO3 + HO2 =  #.44 HO2 + #.44 CO + &                                   
        #.44 CO2 + #.56 GLY + #.15 O3 + #.44 OH                                 SAME K AS BR22          121    88     0 = -2011    88 -2011    89                
                                                                                                                    -2011    90 -2017   118
                                                                                                                    -2013     1 -2011     9
 BP34   MGLY + HV = HO2 + CO + MECO3                                            PHOT. = MGLY-06         116     8     0 =    88    89    98                      
 BP35   MGLY + OH = CO + MECO3                                             2.202E+04   0.00  -1.000     116     9     0 =    89    98                            
 BP36   MGLY + NO3 = HNO3 + CO + MECO3                                     2.055E+03   3.77  -1.000     116    10     0 =    12    89    98                      
 BP37   BACL + HV = #2 MECO3                                                    PHOT. = BACL-07         122     8     0 = -2001    98                            
 BP38   CRES + OH = #.700 HO2 + #.100 BZO + #.170 xHO2 + &                       
        #.030 OH + #.170 RO2C + #.700 CATL + #.030 AFG3 + &                      
        #.085 xAFG1 + #.085 xAFG2 + #.085 xGLY + &                               
        #.085 xMGLY + #.100 xCNDPP + #.130 yRAOOH + &                            
        #.040 yRAOOHP + #.115 XC                                           2.349E+03  -1.93  -1.000     107     9     0 = -2041    88 -2014    54                
                                                                                                                    -2082    56 -2076     9
                                                                                                                    -2082    96 -2041   123
                                                                                                                    -2076   124 -2062    75
                                                                                                                    -2062    76 -2062    71
                                                                                                                    -2062    72 -2014   125
                                                                                                                    -2083    84 -2084   126
                                                                                                                    -2085   104
 BP39   CRES + NO3 = #.100 HNO3 + #.900 XN + #.700 HO2 + &                       
        #.100 BZO + #.170 xHO2 + #.030 OH + #.170 RO2C + &                       
        #.700 CATL + #.030 AFG3 + #.085 xAFG1 + &                                
        #.085 xAFG2 + #.085 xGLY + #.085 xMGLY + &                               
        #.130 yRAOOH + #.715 XC                                            2.055E+04   0.00  -1.000     107    10     0 = -2014    12 -2015    18                
                                                                                                                    -2041    88 -2014    54
                                                                                                                    -2082    56 -2076     9
                                                                                                                    -2082    96 -2041   123
                                                                                                                    -2076   124 -2062    75
                                                                                                                    -2062    76 -2062    71
                                                                                                                    -2062    72 -2083    84
                                                                                                                    -2086   104
 BP85   XYNL + OH = #.700 HO2 + #.070 BZO + #.230 xHO2 + &                       
        #.230 RO2C + #.700 CATL + #.115 xAFG1 + &                                
        #.115 xAFG2 + #.115 xGLY + #.115 xMGLY + &                               
        #.110 xCNDPP + #.020 yRAOOH + #.210 yRAOOHP + &                          
        #.995 XC                                                           1.083E+05   0.00  -1.000     127     9     0 = -2041    88 -2087    54                
                                                                                                                    -2088    56 -2088    96
                                                                                                                    -2041   123 -2085    75
                                                                                                                    -2085    76 -2085    71
                                                                                                                    -2085    72 -2089   125
                                                                                                                    -2056    84 -2090   126
                                                                                                                    -2091   104
 BP86   XYNL + NO3 = #.070 HNO3 + #.930 XN + #.700 HO2 + &                       
        #.070 BZO + #.230 xHO2 + #.230 RO2C + #.700 CATL + &                     
        #.115 xAFG1 + #.115 xAFG2 + #.115 xGLY + &                               
        #.115 xMGLY + #.020 yRAOOH + #1.655 XC                             4.492E+04   0.00  -1.000     127    10     0 = -2087    12 -2092    18                
                                                                                                                    -2041    88 -2087    54
                                                                                                                    -2088    56 -2088    96
                                                                                                                    -2041   123 -2085    75
                                                                                                                    -2085    76 -2085    71
                                                                                                                    -2085    72 -2056    84
                                                                                                                    -2093   104
 BP87   CATL + OH = #.400 HO2 + #.200 BZO + #.200 xHO2 + &                       
        #.200 OH + #.200 RO2C + #.200 AFG3 + #.100 xAFG1 + &                     
        #.100 xAFG2 + #.100 xGLY + #.100 xMGLY + &                               
        #.330 CNDPP + #.200 yRAOOH + #-.080 XC                             2.936E+05   0.00  -1.000     123     9     0 = -2010    88 -2005    54                
                                                                                                                    -2005    56 -2005     9
                                                                                                                    -2005    96 -2005   124
                                                                                                                    -2014    75 -2014    76
                                                                                                                    -2014    71 -2014    72
                                                                                                                    -2094   128 -2005    84
                                                                                                                    -2095   104
 BP88   CATL + NO3 = #.200 HNO3 + #.800 XN + #.400 HO2 + &                       
        #.200 BZO + #.200 xHO2 + #.200 OH + #.200 RO2C + &                       
        #.200 AFG3 + #.100 xAFG1 + #.100 xAFG2 + &                               
        #.100 xGLY + #.100 xMGLY + #.200 yRAOOH + &                              
        #1.900 XC                                                          2.495E+05   0.00  -1.000     123    10     0 = -2005    12 -2004    18                
                                                                                                                    -2010    88 -2005    54
                                                                                                                    -2005    56 -2005     9
                                                                                                                    -2005    96 -2005   124
                                                                                                                    -2014    75 -2014    76
                                                                                                                    -2014    71 -2014    72
                                                                                                                    -2005    84 -2096   104
 BP40   NPHE + OH = BZO + XN                                               5.138E+03   0.00  -1.000      17     9     0 =    54    18                            
 BP41   NPHE + HV + #QY.NO2A = HONO + #6 XC                                     PHOT. = NO2-06           17     8    -9 =    13 -2097   104                      
 BP42   NPHE + HV + #QY.NPHE = #6 XC + XN                                       PHOT. = NO2-06           17     8   -10 = -2097   104    18                      
 BP43   BALD + OH = BZCO3                                                  1.761E+04   0.00  -1.000     129     9     0 =   103                                  
 BP44   BALD + HV + #QY.BALD = #7 XC                                            PHOT. = BALD-06         129     8   -11 = -2098   104                            
 BP45   BALD + NO3 = HNO3 + BZCO3                                          1.967E+03   3.70  -1.000     129    10     0 =    12   103                            
 BP83   PHEN + OH = #.700 HO2 + #.100 BZO + #.095 xHO2 + &                       
        #.105 OH + #.095 RO2C + #.700 CATL + #.105 AFG3 + &                      
        #.048 xAFG1 + #.048 xAFG2 + #.095 xGLY + &                               
        #.065 yRAOOH + #.030 yRAOOHP + #-.205 XC                           6.899E+02  -2.42  -1.000     130     9     0 = -2041    88 -2014    54                
                                                                                                                    -2099    56 -2100     9
                                                                                                                    -2099    96 -2041   123
                                                                                                                    -2100   124 -2101    75
                                                                                                                    -2101    76 -2099    71
                                                                                                                    -2102    84 -2076   126
                                                                                                                    -2103   104
 BP84   PHEN + NO3 = #.100 HNO3 + #.900 XN + #.700 HO2 + &                       
        #.100 BZO + #.095 xHO2 + #.105 OH + #.095 RO2C + &                       
        #.700 CATL + #.105 AFG3 + #.048 xAFG1 + &                                
        #.048 xAFG2 + #.095 xGLY + #.065 yRAOOH + &                              
        #.030 yRAOOHP + #-.205 XC                                          5.578E+03   0.00  -1.000     130    10     0 = -2014    12 -2015    18                
                                                                                                                    -2041    88 -2014    54
                                                                                                                    -2099    56 -2100     9
                                                                                                                    -2099    96 -2041   123
                                                                                                                    -2100   124 -2101    75
                                                                                                                    -2101    76 -2099    71
                                                                                                                    -2102    84 -2076   126
                                                                                                                    -2103   104
 BP46   AFG1 + OH = #.217 MACO3 + #.723 RO2C + &                                 
        #.060 RO2XC + #.060 zRNO3 + #.521 xHO2 + &                               
        #.201 xMECO3 + #.334 xCO + #.407 xRCHO + &                               
        #.129 xMEK + #.107 xGLY + #.267 xMGLY + &                                
        #.783 yR6OOH + #.284 XC                                            1.086E+05   0.00  -1.000     119     9     0 = -2104   105 -2105    96                
                                                                                                                    -2106    97 -2106    81
                                                                                                                    -2107    56 -2108    60
                                                                                                                    -2109    64 -2110    67
                                                                                                                    -2111    69 -2112    71
                                                                                                                    -2113    72 -2114    83
                                                                                                                    -2115   104
 BP48   AFG1 + HV = #1.023 HO2 + #.173 MEO2 + &                                  
        #.305 MECO3 + #.500 MACO3 + #.695 CO + #.195 GLY + &                     
        #.305 MGLY + #.217 XC                                                   PHOT. = AFG1            119     8     0 = -2116    88 -2117    92                
                                                                                                                    -2118    98 -2006   105
                                                                                                                    -2119    89 -2120   118
                                                                                                                    -2118   116 -2104   104
 BP49   AFG2 + OH = #.217 MACO3 + #.723 RO2C + &                                 
        #.060 RO2XC + #.060 zRNO3 + #.521 xHO2 + &                               
        #.201 xMECO3 + #.334 xCO + #.407 xRCHO + &                               
        #.129 xMEK + #.107 xGLY + #.267 xMGLY + &                                
        #.783 yR6OOH + #.284 XC                                            1.086E+05   0.00  -1.000     120     9     0 = -2104   105 -2105    96                
                                                                                                                    -2106    97 -2106    81
                                                                                                                    -2107    56 -2108    60
                                                                                                                    -2109    64 -2110    67
                                                                                                                    -2111    69 -2112    71
                                                                                                                    -2113    72 -2114    83
                                                                                                                    -2115   104
 BP51   AFG2 + HV = PROD2 + #-1 XC                                              PHOT. = AFG1            120     8     0 =   114 -2018   104                      
 BP52   AFG3 + OH = #.206 MACO3 + #.733 RO2C + &                                 
        #.117 RO2XC + #.117 zRNO3 + #.561 xHO2 + &                               
        #.117 xMECO3 + #.114 xCO + #.274 xGLY + &                                
        #.153 xMGLY + #.019 xBACL + #.195 xAFG1 + &                              
        #.195 xAFG2 + #.231 xIPRD + #.794 yR6OOH + &                             
        #.938 XC                                                           1.372E+05   0.00  -1.000     124     9     0 = -2121   105 -2122    96                
                                                                                                                    -2123    97 -2123    81
                                                                                                                    -2124    56 -2123    60
                                                                                                                    -2125    64 -2126    71
                                                                                                                    -2127    72 -2128    73
                                                                                                                    -2120    75 -2120    76
                                                                                                                    -2129    79 -2130    83
                                                                                                                    -2131   104
 BP53   AFG3 + O3 = #.471 OH + #.554 HO2 + #.013 MECO3 + &                       
        #.258 RO2C + #.007 RO2XC + #.007 zRNO3 + #.580 CO + &                    
        #.190 CO2 + #.366 GLY + #.184 MGLY + #.350 AFG1 + &                      
        #.350 AFG2 + #.139 AFG3 + #.003 MACR + #.004 MVK + &                     
        #.003 IPRD + #.095 xHO2 + #.163 xRCO3 + &                                
        #.163 xHCHO + #.095 xMGLY + #.264 yR6OOH + &                             
        #-.575 XC                                                          2.099E-02   0.00  -1.000     124     1     0 = -2132     9 -2133    88                
                                                                                                                    -2134    98 -2135    96
                                                                                                                    -2136    97 -2136    81
                                                                                                                    -2137    89 -2138    90
                                                                                                                    -2139   118 -2140   116
                                                                                                                    -2141   119 -2141   120
                                                                                                                    -2068   124 -2142   131
                                                                                                                    -2044   132 -2142   117
                                                                                                                    -2099    56 -2143    61
                                                                                                                    -2143    65 -2099    72
                                                                                                                    -2144    83 -2145   104
 BP89   AFG4 + OH = #.902 RO2C + #.098 RO2XC + &                                 
        #.098 zRNO3 + #.902 xMECO3 + #.902 xRCHO + yROOH + &                     
        #.902 XC                                                           9.248E+04   0.00  -1.000     133     9     0 = -2146    96 -2147    97                
                                                                                                                    -2147    81 -2146    60
                                                                                                                    -2146    67    82 -2146
                                                                                                                      104
 BP54   MACR + OH = #.5 MACO3 + #.5 RO2C + #.5 xHO2 + &                          
        #.416 xCO + #.084 xHCHO + #.416 xMEK + &                                 
        #.084 xMGLY + #.5 yROOH + #-0.416 XC                               1.174E+04  -0.75  -1.000     131     9     0 = -2006   105 -2006    96                
                                                                                                                    -2006    56 -2148    64
                                                                                                                    -2058    65 -2148    69
                                                                                                                    -2058    72 -2006    82
                                                                                                                    -2149   104
 BP55   MACR + O3 = #.208 OH + #.108 HO2 + #.1 RO2C + &                          
        #.45 CO + #.117 CO2 + #.1 HCHO + #.9 MGLY + &                            
        #.333 HCOOH + #.1 xRCO3 + #.1 xHCHO + #.1 yROOH + &                      
        #-0.1 XC                                                           2.055E+00   4.17  -1.000     131     1     0 = -2150     9 -2151    88                
                                                                                                                    -2014    96 -2152    89
                                                                                                                    -2123    90 -2014    93
                                                                                                                    -2015   116 -2153   111
                                                                                                                    -2014    61 -2014    65
                                                                                                                    -2014    82 -2154   104
 BP56   MACR + NO3 = #.5 MACO3 + #.5 RO2C + #.5 HNO3 + &                         
        #.5 xHO2 + #.5 xCO + &                                                   
        #.5 yROOH + #1.5 XC + #.5 XN                                       2.202E+03   3.61  -1.000     131    10     0 = -2006   105 -2006    96                
                                                                                                                    -2006    12 -2006    56
                                                                                                                    -2006    64 -2006    82
                                                                                                                    -2155   104 -2006    18
 BP57   MACR + O3P = RCHO + XC                                             9.306E+03   0.00  -1.000     131    51     0 =   109   104                            
 BP58   MACR + HV = #.33 OH + #.67 HO2 + #.34 MECO3 + &                          
        #.33 MACO3 + #.33 RO2C + #.67 CO + #.34 HCHO + &                         
        #.33 xMECO3 + #.33 xHCHO + #.33 yROOH                                   PHOT. = MACR-06         131     8     0 = -2094     9 -2156    88                
                                                                                                                    -2157    98 -2094   105
                                                                                                                    -2094    96 -2156    89
                                                                                                                    -2157    93 -2094    60
                                                                                                                    -2094    65 -2094    82
 BP59   MVK + OH = #.975 RO2C + #.025 RO2XC + &                                  
        #.025 zRNO3 + &                                                          
        #.3 xHO2 + #.675 xMECO3 + #.3 xHCHO + &                                  
        #.675 xRCHO + #.3 xMGLY + yROOH + #-0.725 XC                       3.817E+03  -1.21  -1.000     132     9     0 = -2158    96 -2159    97                
                                                                                                                    -2159    81 -2033    56
                                                                                                                    -2160    60 -2033    65
                                                                                                                    -2160    67 -2033    72
                                                                                                                       82 -2161   104
 BP60   MVK + O3 = #.164 OH + #.064 HO2 + #.05 RO2C + &                          
        #.05 xHO2 + #.475 CO + #.124 CO2 + #.05 HCHO + &                         
        #.95 MGLY + #.351 HCOOH + #.05 xRCO3 + &                                 
        #.05 xHCHO + #.05 yROOH + #-0.05 XC                                1.248E+00   3.02  -1.000     132     1     0 = -2162     9 -2163    88                
                                                                                                                    -2164    96 -2164    56
                                                                                                                    -2165    89 -2071    90
                                                                                                                    -2164    93 -2166   116
                                                                                                                    -2167   111 -2164    61
                                                                                                                    -2164    65 -2164    82
                                                                                                                    -2168   104
 BP62   MVK + O3P = #.45 RCHO + #.55 MEK + #.45 XC                         6.341E+03   0.00  -1.000     132    51     0 = -2152   109 -2169   110                
                                                                                                                    -2152   104
 BP63   MVK + HV = #.4 MEO2 + #.6 CO + #.6 PROD2 + &                             
        #.4 MACO3 + #-2.2 XC                                                    PHOT. = MVK-06          132     8     0 = -2010    92 -2009    89                
                                                                                                                    -2009   114 -2010   105
                                                                                                                    -2170   104
 BP64   IPRD + OH = #.289 MACO3 + #.67 RO2C + #.67 xHO2 + &                      
        #.041 RO2XC + #.041 zRNO3 + #.336 xCO + #.055 xHCHO + &                  
        #.129 xCCHO + #.013 xRCHO + #.15 xMEK + &                                
        #.332 xPROD2 + #.15 xGLY + #.174 xMGLY + &                               
        #-0.504 XC + #.711 yR6OOH                                          9.086E+04   0.00  -1.000     117     9     0 = -2171   105 -2156    96                
                                                                                                                    -2156    56 -2055    97
                                                                                                                    -2055    81 -2172    64
                                                                                                                    -2173    65 -2111    66
                                                                                                                    -2134    67 -2013    69
                                                                                                                    -2174    70 -2013    71
                                                                                                                    -2175    72 -2176   104
                                                                                                                    -2177    83
 BP65   IPRD + O3 = #.285 OH + #.4 HO2 + #.048 RO2C + &                          
        #.048 xRCO3 + #.498 CO + #.14 CO2 + #.124 HCHO + &                       
        #.21 MEK + #.023 GLY + #.742 MGLY + #.1 HCOOH + &                        
        #.372 RCOOH + #.047 xCCHO + #.001 xHCHO + &                              
        #.048 yR6OOH + #-.329 XC                                           6.136E-03   0.00  -1.000     117     1     0 = -2178     9 -2010    88                
                                                                                                                    -2101    96 -2101    61
                                                                                                                    -2179    89 -2180    90
                                                                                                                    -2071    93 -2090   110
                                                                                                                    -2181   118 -2182   116
                                                                                                                    -2014   111 -2183   102
                                                                                                                    -2184    66 -2185    65
                                                                                                                    -2101    83 -2186   104
 BP66   IPRD + NO3 = #.15 MACO3 + #.15 HNO3 + #.799 RO2C + &                     
        #.799 xHO2 + #.051 RO2XC + #.051 zRNO3 + #.572 xCO + &                   
        #.227 xHCHO + #.218 xRCHO + #.008 xMGLY + &                              
        #.572 xRNO3 + #.85 yR6OOH + #.278 XN + #-1 XC + #.185 XC           1.468E+02   0.00  -1.000     117    10     0 = -2013   105 -2013    12                
                                                                                                                    -2187    96 -2187    56
                                                                                                                    -2188    97 -2188    81
                                                                                                                    -2189    64 -2190    65
                                                                                                                    -2191    67 -2192    72
                                                                                                                    -2189    80 -2193    83
                                                                                                                    -2194    18 -2018   104
                                                                                                                    -2195   104
 BP67   IPRD + HV = #1.233 HO2 + #.467 MECO3 + #.3 RCO3 + &                      
        #1.233 CO + #.3 HCHO + #.467 CCHO + #.233 MEK + &                        
        #-.233 XC                                                               PHOT. = MACR-06         117     8     0 = -2196    88 -2197    98                
                                                                                                                    -2033   100 -2196    89
                                                                                                                    -2033    93 -2197   108
                                                                                                                    -2198   110 -2199   104
 BP68   PROD2 + OH = #.472 HO2 + #.379 xHO2 + &                                  
        #.029 xMECO3 + #.049 xRCO3 + #.473 RO2C + &                              
        #.071 RO2XC + #.071 zRNO3 + #.002 HCHO + &                               
        #.211 xHCHO + #.001 CCHO + #.083 xCCHO + &                               
        #.143 RCHO + #.402 xRCHO + #.115 xMEK + &                                
        #.329 PROD2 + #.007 xPROD2 + #.528 yR6OOH + &                            
        #.877 XC                                                           2.275E+04   0.00  -1.000     114     9     0 = -2200    88 -2201    56                
                                                                                                                    -2053    60 -2202    61
                                                                                                                    -2203    96 -2204    97
                                                                                                                    -2204    81 -2205    93
                                                                                                                    -2206    65 -2185   108
                                                                                                                    -2207    66 -2037   109
                                                                                                                    -2208    67 -2085    69
                                                                                                                    -2209   114 -2136    70
                                                                                                                    -2210    83 -2211   104
 BP69   PROD2 + HV + #QY.PROD2 = #.913 xHO2 + #.4 MECO3 + #.6 RCO3 + &           
        #1.59 RO2C + #.087 RO2XC + #.087 zRNO3 + &                               
        #.303 xHCHO + #.163 xCCHO + #.78 xRCHO + yR6OOH + &                      
        #-.091 XC                                                               PHOT. = MEK-06          114     8   -12 = -2212    56 -2010    98                
                                                                                                                    -2009   100 -2213    96
                                                                                                                    -2214    97 -2214    81
                                                                                                                    -2215    65 -2143    66
                                                                                                                    -2216    67    83 -2217
                                                                                                                      104
 BP70   RNO3 + OH = #.189 HO2 + #.305 xHO2 + #.019 NO2 + &                       
        #.313 xNO2 + #.976 RO2C + #.175 RO2XC + &                                
        #.175 zRNO3 + #.011 xHCHO + #.429 xCCHO + &                              
        #.001 RCHO + #.036 xRCHO + #.004 xACET + &                               
        #.01 MEK + #.17 xMEK + #.008 PROD2 + &                                   
        #.031 xPROD2 + #.189 RNO3 + #.305 xRNO3 + &                              
        #.157 yROOH + #.636 yR6OOH + #.174 XN + #.04 XC                    1.057E+04   0.00  -1.000      19     9     0 = -2218    88 -2118    56                
                                                                                                                    -2128     3 -2219    58
                                                                                                                    -2220    96 -2221    97
                                                                                                                    -2221    81 -2222    65
                                                                                                                    -2223    66 -2185   109
                                                                                                                    -2224    67 -2044    68
                                                                                                                    -2225   110 -2082    69
                                                                                                                    -2192   114 -2226    70
                                                                                                                    -2218    19 -2118    80
                                                                                                                    -2227    82 -2228    83
                                                                                                                    -2175    18 -2084   104
 BP71   RNO3 + HV = #.344 HO2 + #.554 xHO2 + NO2 + &                             
        #.721 RO2C + #.102 RO2XC + #.102 zRNO3 + &                               
        #.074 HCHO + #.061 xHCHO + #.214 CCHO + &                                
        #.23 xCCHO + #.074 RCHO + #.063 xRCHO + &                                
        #.008 xACET + #.124 MEK + #.083 xMEK + &                                 
        #.19 PROD2 + #.261 xPROD2 + #.066 yROOH + &                              
        #.591 yR6OOH + #.396 XC                                                 PHOT. = IC3ONO2          19     8     0 = -2229    88 -2133    56                
                                                                                                                        3 -2230    96 -2231
                                                                                                                       97 -2231    81 -2030
                                                                                                                       93 -2232    65 -2233
                                                                                                                      108 -2088    66 -2030
                                                                                                                      109 -2234    67 -2192
                                                                                                                       68 -2071   110 -2207
                                                                                                                       69 -2138   114 -2235
                                                                                                                       70 -2236    82 -2237
                                                                                                                       83 -2238   104
 BP72   ETOH + OH = #.95 HO2 + #.05 xHO2 + #.05 RO2C + &                         
        #.081 xHCHO + #.95 CCHO + #.01 xCCHO + &                                 
        #.05 yROOH + #-.001 XC                                             8.059E+02  -1.05   1.000     134     9     0 = -2166    88 -2164    56                
                                                                                                                    -2164    96 -2239    65
                                                                                                                    -2166   108 -2225    66
                                                                                                                    -2164    82 -2240   104
 PO01   xHCHO + #RO2RO  = HCHO                                             1.000E+00   0.00   0.000      65    -5     0 =    93                                  
 PO02   xHCHO + #RO2XRO = XC                                               1.000E+00   0.00   0.000      65    -6     0 =   104                                  
 PO03   xCCHO + #RO2RO = CCHO                                              1.000E+00   0.00   0.000      66    -5     0 =   108                                  
 PO04   xCCHO + #RO2XRO = #2 XC                                            1.000E+00   0.00   0.000      66    -6     0 = -2001   104                            
 PO05   xRCHO + #RO2RO = RCHO                                              1.000E+00   0.00   0.000      67    -5     0 =   109                                  
 PO06   xRCHO + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      67    -6     0 = -2019   104                            
 PO07   xACET + #RO2RO = ACET                                              1.000E+00   0.00   0.000      68    -5     0 =   106                                  
 PO08   xACET + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      68    -6     0 = -2019   104                            
 PO09   xMEK + #RO2RO = MEK                                                1.000E+00   0.00   0.000      69    -5     0 =   110                                  
 PO10   xMEK + #RO2XRO = #4 XC                                             1.000E+00   0.00   0.000      69    -6     0 = -2020   104                            
 PO11   xPROD2 + #RO2RO = PROD2                                            1.000E+00   0.00   0.000      70    -5     0 =   114                                  
 PO12   xPROD2 + #RO2XRO = #6 XC                                           1.000E+00   0.00   0.000      70    -6     0 = -2097   104                            
 PO13   xGLY + #RO2RO = GLY                                                1.000E+00   0.00   0.000      71    -5     0 =   118                                  
 PO14   xGLY + #RO2XRO = #2 XC                                             1.000E+00   0.00   0.000      71    -6     0 = -2001   104                            
 PO15   xMGLY + #RO2RO = MGLY                                              1.000E+00   0.00   0.000      72    -5     0 =   116                                  
 PO16   xMGLY + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      72    -6     0 = -2019   104                            
 PO17   xBACL + #RO2RO = BACL                                              1.000E+00   0.00   0.000      73    -5     0 =   122                                  
 PO18   xBACL + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      73    -6     0 = -2020   104                            
 PO19   xBALD + #RO2RO = BALD                                              1.000E+00   0.00   0.000      74    -5     0 =   129                                  
 PO20   xBALD + #RO2XRO = #7 XC                                            1.000E+00   0.00   0.000      74    -6     0 = -2098   104                            
 PO21   xAFG1 + #RO2RO = AFG1                                              1.000E+00   0.00   0.000      75    -5     0 =   119                                  
 PO22   xAFG1 + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      75    -6     0 = -2241   104                            
 PO23   xAFG2 + #RO2RO = AFG2                                              1.000E+00   0.00   0.000      76    -5     0 =   120                                  
 PO24   xAFG2 + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      76    -6     0 = -2241   104                            
 PO50   xAFG4 + #RO2XRO = #6 XC                                            1.000E+00   0.00   0.000     135    -6     0 = -2097   104                            
 PO51   xAFG4 + #RO2RO = AFG4                                              1.000E+00   0.00   0.000     135    -5     0 =   133                                  
 PO27   xMACR + #RO2RO = MACR                                              1.000E+00   0.00   0.000      77    -5     0 =   131                                  
 PO28   xMACR + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      77    -6     0 = -2020   104                            
 PO29   xMVK + #RO2RO = MVK                                                1.000E+00   0.00   0.000      78    -5     0 =   132                                  
 PO30   xMVK + #RO2XRO = #4 XC                                             1.000E+00   0.00   0.000      78    -6     0 = -2020   104                            
 PO31   xIPRD + #RO2RO = IPRD                                              1.000E+00   0.00   0.000      79    -5     0 =   117                                  
 PO32   xIPRD + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      79    -6     0 = -2241   104                            
 PO33   xRNO3 + #RO2RO = RNO3                                              1.000E+00   0.00   0.000      80    -5     0 =    19                                  
 PO34   xRNO3 + #RO2XRO = #6 XC + XN                                       1.000E+00   0.00   0.000      80    -6     0 = -2097   104    18                      
 PO35   zRNO3 + #RO2NO = RNO3 + #-1 XN                                     1.000E+00   0.00   0.000      81   -13     0 =    19 -2018    18                      
 PO36   zRNO3 + #RO22NN = PROD2 + HO2                                      1.000E+00   0.00   0.000      81   -14     0 =   114    88                            
 PO37   zRNO3 + #RO2XRO = #6 XC                                            1.000E+00   0.00   0.000      81    -6     0 = -2097   104                            
 PO38   yROOH + #RO2HO2 = ROOH + #-3 XC                                    1.000E+00   0.00   0.000      82   -15     0 =   112 -2242   104                      
 PO39   yROOH + #RO2RO2M = MEK + #-4 XC                                    1.000E+00   0.00   0.000      82   -16     0 =   110 -2243   104                      
 PO40   yROOH + #RO2RO =                                                   1.000E+00   0.00   0.000      82    -5     0 =                                        
 PO41   yR6OOH + #RO2HO2 = R6OOH + #-6 XC                                  1.000E+00   0.00   0.000      83   -15     0 =   113 -2244   104                      
 PO42   yR6OOH + #RO2RO2M  = PROD2 + #-6 XC                                1.000E+00   0.00   0.000      83   -16     0 =   114 -2244   104                      
 PO43   yR6OOH + #RO2RO =                                                  1.000E+00   0.00   0.000      83    -5     0 =                                        
 PO44   yRAOOH + #RO2HO2 = RAOOH + #-8 XC                                  1.000E+00   0.00   0.000      84   -15     0 =   115 -2245   104                      
 PO45   yRAOOH + #RO2RO2M  = PROD2 + #-6 XC                                1.000E+00   0.00   0.000      84   -16     0 =   114 -2244   104                      
 PO46   yRAOOH + #RO2RO =                                                  1.000E+00   0.00   0.000      84    -5     0 =                                        
 BE01   CH4 + OH = H2O + MEO2                                              2.716E+03   3.36  -1.000     136     9     0 =     6    92                            
 BE02   ETHENE + OH = RO2C + xHO2 + #1.61 xHCHO + &                              
        #.195 xCCHO + yROOH                                               FALLOFF F= 0.600, N= 1.000    137     9     0 =    96    56 -2246    65                
                                                                                                                    -2120    66    82
                                                                      K0:  3.591E+00   0.00  -6.500
                                                                      KI:  1.292E+04   0.00  -1.850
 BE03   ETHENE + O3 = #.16 OH + #.16 HO2 + #.51 CO + &                           
        #.12 CO2 + HCHO + #.37 HCOOH                                       1.342E+01   5.13  -1.000     137     1     0 = -2059     9 -2059    88                
                                                                                                                    -2029    89 -2247    90
                                                                                                                       93 -2248   111
 BE04   ETHENE + NO3 = RO2C + xHO2 + xRCHO + yROOH + &                           
        #-1 XC + XN                                                        4.844E+03   5.72   1.000     137    10     0 =    96    56    67    82                
                                                                                                                    -2018   104    18
 BE05   ETHENE + O3P = #.8 HO2 + #.51 MEO2 + #.29 RO2C + &                       
        #.51 CO + #.1 CCHO + #.29 xHO2 + #.278 xCO + &                           
        #.278 xHCHO + #.012 xGLY + #.29 yROOH + #.2 XC                     1.571E+04   1.59  -1.000     137    51     0 = -2004    88 -2029    92                
                                                                                                                    -2249    96 -2029    89
                                                                                                                    -2014   108 -2249    56
                                                                                                                    -2194    64 -2194    65
                                                                                                                    -2046    71 -2249    82
                                                                                                                    -2005   104
 BE06   ISOPRN + OH = #.986 RO2C + #.093 RO2XC + &                               
        #.093 zRNO3 + #.907 xHO2 + #.624 xHCHO + #.23 xMACR + &                  
        #.32 xMVK + #.357 xIPRD + yR6OOH + &                                     
        #-1.000 XC + #0.833 XC                                             3.728E+04  -0.81  -1.000      37     9     0 = -2250    96 -2251    97                
                                                                                                                    -2251    81 -2252    56
                                                                                                                    -2253    65 -2088    77
                                                                                                                    -2254    78 -2255    79
                                                                                                                       83 -2018   104 -2256
                                                                                                                      104
 BE07   ISOPRN + O3 = #.266 OH + #.066 HO2 + &                                   
        #.192 RO2C + #.008 RO2XC + #.008 zRNO3 + #.275 CO + &                    
        #.122 CO2 + #.4 HCHO + #.1 PROD2 + #.39 MACR + &                         
        #.16 MVK + #.15 IPRD + #.204 HCOOH + &                                   
        #.192 xMACO3 + #.192 xHCHO + #.2 yR6OOH + &                              
        #-1.000 XC + #0.401 XC                                             1.154E+01   3.80  -1.000      37     1     0 = -2257     9 -2236    88                
                                                                                                                    -2258    96 -2192    97
                                                                                                                    -2192    81 -2259    89
                                                                                                                    -2260    90 -2010    93
                                                                                                                    -2014   114 -2003   131
                                                                                                                    -2059   132 -2013   117
                                                                                                                    -2261   111 -2258    62
                                                                                                                    -2258    65 -2005    83
                                                                                                                    -2018   104 -2262   104
 BE08   ISOPRN + NO3 = #.936 RO2C + #.064 RO2XC + &                              
        #.064 zRNO3 + #.749 xHO2 + #.187 xNO2 + #.936 xIPRD + &                  
        yR6OOH + #-1 XC + #.936 XC + #.813 XN                              4.448E+03   0.89  -1.000      37    10     0 = -2263    96 -2163    97                
                                                                                                                    -2163    81 -2264    56
                                                                                                                    -2265    58 -2263    79
                                                                                                                       83 -2018   104 -2263
                                                                                                                      104 -2266    18
 BE09   ISOPRN + O3P = #.25 MEO2 + #.24 RO2C + &                                 
        #.01 RO2XC + #.01 zRNO3 + #.75 PROD2 + #.24 xMACO3 + &                   
        #.24 xHCHO + #.25 yR6OOH + &                                             
        #-1.0 XC                                                           5.138E+04   0.00  -1.000      37    51     0 = -2008    92 -2267    96                
                                                                                                                    -2225    97 -2225    81
                                                                                                                    -2007   114 -2267    62
                                                                                                                    -2267    65 -2008    83
                                                                                                                    -2018   104
 BE10   ACETYLEN + OH = #.7 OH + #.3 HO2 + #.3 CO + &                            
        #.7 GLY + #.3 HCOOH                                               FALLOFF F= 0.600, N= 1.000    138     9     0 = -2041     9 -2033    88                
                                                                                                                    -2033    89 -2041   118
                                                                                                                    -2033   111
                                                                      K0:  1.975E-01   0.00  -2.000
                                                                      KI:  1.218E+03   0.00   1.000
 BE11   ACETYLEN + O3 = #.5 OH + #1.5 HO2 + #1.5 CO + &                          
        #.5 CO2                                                            1.468E+01   8.15  -1.000     138     1     0 = -2006     9 -2155    88                
                                                                                                                    -2155    89 -2006    90
 BE12   BENZ + OH = #.116 OH + #.29 RO2C + #.29 xHO2 + &                         
        #.024 RO2XC + #.024 zRNO3 + #.57 HO2 + #.57 CRES + &                     
        #.116 AFG3 + #.290 xGLY + #.029 xAFG1 + &                                
        #.261 xAFG2 + #.314 yRAOOH + &                                           
        #-1.0 XC + #0.024 XC + BNZSOMG                                     3.420E+03   0.38  -1.000      29     9     0 = -2268     9 -2249    96                
                                                                                                                    -2249    56 -2074    97
                                                                                                                    -2074    81 -2269    88
                                                                                                                    -2269   107 -2268   124
                                                                                                                    -2249    71 -2053    75
                                                                                                                    -2235    76 -2270    84
                                                                                                                    -2018   104 -2074   104
                                                                                                                       46
 BL01   ALK1 + OH = xHO2 + RO2C + xCCHO + yROOH                            1.967E+03   0.99   1.000     139     9     0 =    56    96    66    82                
 BL02   ALK2 + OH = #.965 xHO2 + #.965 RO2C + &                                  
        #.035 RO2XC + #.035 zRNO3 + #.261 xRCHO + &                              
        #.704 xACET + yROOH + #-.105 XC                                    2.187E+03   0.17   1.000     140     9     0 = -2021    56 -2021    96                
                                                                                                                    -2022    97 -2022    81
                                                                                                                    -2235    67 -2271    68
                                                                                                                       82 -2272   104
 BL03   ALK3 + OH = #.695 xHO2 + #.236 xTBUO + &                                 
        #1.253 RO2C + #.07 RO2XC + #.07 zRNO3 + &                                
        #.026 xHCHO + #.445 xCCHO + #.122 xRCHO + &                              
        #.024 xACET + #.332 xMEK + #.983 yROOH + &                               
        #.017 yR6OOH + #-.046 XC                                           2.217E+03  -0.25  -1.000     141     9     0 = -2119    56 -2273    63                
                                                                                                                    -2274    96 -2087    97
                                                                                                                    -2087    81 -2063    65
                                                                                                                    -2275    66 -2260    67
                                                                                                                    -2074    68 -2174    69
                                                                                                                    -2276    82 -2067    83
                                                                                                                    -2277   104
 BL04   ALK4 + OH = #.83 xHO2 + #.01 xMEO2 + &                                   
        #.011 xMECO3 + #1.763 RO2C + #.149 RO2XC + &                             
        #.149 zRNO3 + #.002 xCO + #.029 xHCHO + &                                
        #.438 xCCHO + #.236 xRCHO + #.426 xACET + &                              
        #.106 xMEK + #.146 xPROD2 + yR6OOH + #-.119 XC                     5.505E+03  -0.09  -1.000     142     9     0 = -2278    56 -2225    59                
                                                                                                                    -2222    60 -2279    96
                                                                                                                    -2280    97 -2280    81
                                                                                                                    -2205    64 -2053    65
                                                                                                                    -2281    66 -2273    67
                                                                                                                    -2282    68 -2283    69
                                                                                                                    -2284    70    83 -2285
                                                                                                                      104
 BL05   ALK5 + OH = #.647 xHO2 + #1.605 RO2C + &                                 
        #.353 RO2XC + #.353 zRNO3 + #.04 xHCHO + &                               
        #.106 xCCHO + #.209 xRCHO + #.071 xACET + &                              
        #.086 xMEK + #.407 xPROD2 + yR6OOH + #2.004 XC                     2.046E+04   0.00  -1.000      28     9     0 = -2286    56 -2287    96                
                                                                                                                    -2288    97 -2288    81
                                                                                                                    -2084    65 -2283    66
                                                                                                                    -2289    67 -2204    68
                                                                                                                    -2290    69 -2110    70
                                                                                                                       83 -2291   104
 BL06   OLE1 + OH = #.904 xHO2 + #.001 xMEO2 + &                                 
        #1.138 RO2C + #.095 RO2XC + #.095 zRNO3 + &                              
        #.7 xHCHO + #.301 xCCHO + #.47 xRCHO + &                                 
        #.005 xACET + #.026 xMACR + #.008 xMVK + &                               
        #.006 xIPRD + #.119 xPROD2 + #.413 yROOH + &                             
        #.587 yR6OOH + #.822 XC                                            9.072E+03  -1.00  -1.000     143     9     0 = -2292    56 -2185    59                
                                                                                                                    -2293    96 -2099    97
                                                                                                                    -2099    81 -2041    65
                                                                                                                    -2294    66 -2295    67
                                                                                                                    -2296    68 -2063    77
                                                                                                                    -2192    78 -2297    79
                                                                                                                    -2298    70 -2299    82
                                                                                                                    -2300    83 -2301   104
 BL07   OLE1 + O3 = #.116 HO2 + #.04 xHO2 + #.193 OH + &                         
        #.104 MEO2 + #.063 RO2C + #.004 RO2XC + &                                
        #.004 zRNO3 + #.368 CO + #.125 CO2 + #.5 HCHO + &                        
        #.147 CCHO + #.007 xCCHO + #.353 RCHO + &                                
        #.031 xRCHO + #.002 xACET + #.006 MEK + &                                
        #.185 HCOOH + #.022 CCOOH + #.112 RCOOH + &                              
        #.189 PROD2 + #.007 yROOH + #.037 yR6OOH + &                             
        #.69 XC                                                            4.624E+00   3.38  -1.000     143     1     0 = -2268    88 -2084    56                
                                                                                                                    -2302     9 -2303    92
                                                                                                                    -2234    96 -2044    97
                                                                                                                    -2044    81 -2304    89
                                                                                                                    -2305    90 -2006    93
                                                                                                                    -2072   108 -2136    66
                                                                                                                    -2288   109 -2226    67
                                                                                                                    -2205    68 -2297   110
                                                                                                                    -2195   111 -2306    99
                                                                                                                    -2307   102 -2218   114
                                                                                                                    -2136    82 -2308    83
                                                                                                                    -2309   104
 BL08   OLE1 + NO3 = #.824 xHO2 + #1.312 RO2C + &                                
        #.176 RO2XC + #.176 zRNO3 + #.009 xCCHO + &                              
        #.002 xRCHO + #.024 xACET + #.546 xRNO3 + &                              
        #.413 yROOH + #.587 yR6OOH + #.454 XN + #.572 XC                   6.943E+02   2.08  -1.000     143    10     0 = -2310    56 -2311    96                
                                                                                                                    -2312    97 -2312    81
                                                                                                                    -2313    66 -2205    67
                                                                                                                    -2074    68 -2314    80
                                                                                                                    -2299    82 -2300    83
                                                                                                                    -2315    18 -2189   104
 BL09   OLE1 + O3P = #.45 RCHO + #.437 MEK + #.113 PROD2 + &                     
        #1.224 XC                                                          2.187E+04   0.65  -1.000     143    51     0 = -2152   109 -2316   110                
                                                                                                                    -2317   114 -2318   104
 BL10   OLE2 + OH = #.914 xHO2 + #.966 RO2C + &                                  
        #.086 RO2XC + #.086 zRNO3 + #.209 xHCHO + &                              
        #.788 xCCHO + #.481 xRCHO + #.136 xACET + &                              
        #.076 xMEK + #.027 xMACR + #.002 xMVK + &                                
        #.037 xIPRD + #.022 xPROD2 + #.357 yROOH + &                             
        #.643 yR6OOH + #.111 XC                                            1.850E+04  -0.97  -1.000     144     9     0 = -2319    56 -2320    96                
                                                                                                                    -2290    97 -2290    81
                                                                                                                    -2289    65 -2321    66
                                                                                                                    -2322    67 -2323    68
                                                                                                                    -2324    69 -2325    77
                                                                                                                    -2205    78 -2308    79
                                                                                                                    -2306    70 -2255    82
                                                                                                                    -2326    83 -2327   104
 BL11   OLE2 + O3 = #.093 HO2 + #.039 xHO2 + #.423 OH + &                        
        #.29 MEO2 + #.147 xMECO3 + #.008 xRCO3 + &                               
        #.2 RO2C + #.003 RO2XC + #.003 zRNO3 + #.297 CO + &                      
        #.162 CO2 + #.152 HCHO + #.108 xHCHO + &                                 
        #.428 CCHO + #.067 xCCHO + #.315 RCHO + &                                
        #.018 xRCHO + #.048 ACET + #.031 MEK + &                                 
        #.001 xMEK + #.033 HCOOH + #.061 CCOOH + &                               
        #.222 RCOOH + #.028 MACR + #.021 MVK + &                                 
        #.042 PROD2 + #.069 yROOH + #.128 yR6OOH + &                             
        #.125 XC                                                           1.195E+01   2.49  -1.000     144     1     0 = -2251    88 -2027    56                
                                                                                                                    -2328     9 -2249    92
                                                                                                                    -2072    60 -2192    61
                                                                                                                    -2005    96 -2142    97
                                                                                                                    -2142    81 -2329    89
                                                                                                                    -2330    90 -2331    93
                                                                                                                    -2151    65 -2332   108
                                                                                                                    -2333    66 -2334   109
                                                                                                                    -2335    67 -2101   106
                                                                                                                    -2226   110 -2185    69
                                                                                                                    -2336   111 -2232    99
                                                                                                                    -2052   102 -2337   131
                                                                                                                    -2338   132 -2339   114
                                                                                                                    -2340    82 -2341    83
                                                                                                                    -2305   104
 BL12   OLE2 + NO3 = #.423 xHO2 + #.409 xNO2 + &                                 
        #.033 xMEO2 + #1.185 RO2C + #.136 RO2XC + &                              
        #.136 zRNO3 + #.074 xHCHO + #.546 xCCHO + &                              
        #.154 xRCHO + #.11 xACET + #.002 xMEK + &                                
        #.026 xMVK + #.007 xIPRD + #.322 xRNO3 + &                               
        #.357 yROOH + #.643 yR6OOH + #.269 XN + #.114 XC                   3.229E+02  -0.76  -1.000     144    10     0 = -2328    56 -2342    58                
                                                                                                                    -2336    59 -2343    96
                                                                                                                    -2323    97 -2323    81
                                                                                                                    -2030    65 -2314    66
                                                                                                                    -2344    67 -2089    68
                                                                                                                    -2205    69 -2063    78
                                                                                                                    -2136    79 -2345    80
                                                                                                                    -2255    82 -2326    83
                                                                                                                    -2346    18 -2125   104
 BL13   OLE2 + O3P = #.014 HO2 + #.007 xHO2 + &                                  
        #.007 xMACO3 + #.013 RO2C + #.001 RO2XC + &                              
        #.001 zRNO3 + #.006 xCO + #.074 RCHO + #.709 MEK + &                     
        #.006 xMACR + #.202 PROD2 + #.014 yROOH + &                              
        #.666 XC                                                           2.099E+04  -0.22  -1.000     144    51     0 = -2347    88 -2136    56                
                                                                                                                    -2136    62 -2134    96
                                                                                                                    -2185    97 -2185    81
                                                                                                                    -2297    64 -2030   109
                                                                                                                    -2348   110 -2297    77
                                                                                                                    -2349   114 -2347    82
                                                                                                                    -2350   104
 BL14   ARO1 + OH = #.166 HO2 + #.482 xHO2 + #.284 OH + &                        
        #.482 RO2C + #.068 RO2XC + #.068 zRNO3 + &                               
        #.218 xGLY + #.138 xMGLY + #.166 CRES + &                                
        #.049 xBALD + #.164 xAFG1 + #.193 xAFG2 + &                              
        #.284 AFG3 + #.077 xPROD2 + #.403 yRAOOH + &                             
        #.147 yR6OOH + #.002 XC                                            9.028E+03   0.00  -1.000      40     9     0 = -2351    88 -2352    56                
                                                                                                                    -2115     9 -2352    96
                                                                                                                    -2353    97 -2353    81
                                                                                                                    -2191    71 -2354    72
                                                                                                                    -2351   107 -2202    74
                                                                                                                    -2162    75 -2302    76
                                                                                                                    -2115   124 -2061    70
                                                                                                                    -2355    84 -2072    83
                                                                                                                    -2205   104
 BL15   ARO2 + OH = #.108 HO2 + #.58 xHO2 + #.202 OH + &                         
        #.58 RO2C + #.11 RO2XC + #.11 zRNO3 + #.116 xGLY + &                     
        #.286 xMGLY + #.104 xBACL + #.108 CRES + &                               
        #.039 xBALD + #.217 xAFG1 + #.21 xAFG2 + &                               
        #.202 AFG3 + #.035 xPROD2 + &                                            
        #.089 yR6OOH + #.601 yRAOOH + #1.486 XC                            3.773E+04   0.00  -1.000      41     9     0 = -2151    88 -2137    56                
                                                                                                                    -2349     9 -2137    96
                                                                                                                    -2089    97 -2089    81
                                                                                                                    -2268    71 -2356    72
                                                                                                                    -2303    73 -2151   107
                                                                                                                    -2027    74 -2104    75
                                                                                                                    -2090    76 -2349   124
                                                                                                                    -2022    70 -2357    83
                                                                                                                    -2358    84 -2359   104
 BL16   TERP + OH = #.759 xHO2 + #.042 xRCO3 + &                                 
        #1.147 RO2C + #.2 RO2XC + #.2 zRNO3 + #.001 xCO + &                      
        #.264 xHCHO + #.533 xRCHO + #.036 xACET + &                              
        #.005 xMEK + #.009 xMGLY + #.014 xBACL + &                               
        #.002 xMVK + #.001 xIPRD + #.255 xPROD2 + yR6OOH + &                     
        #5.056 XC                                                          2.745E+04  -0.86  -1.000      38     9     0 = -2360    56 -2339    61                
                                                                                                                    -2361    96 -2005    97
                                                                                                                    -2005    81 -2185    64
                                                                                                                    -2144    65 -2362    67
                                                                                                                    -2224    68 -2296    69
                                                                                                                    -2313    72 -2347    73
                                                                                                                    -2205    78 -2185    79
                                                                                                                    -2363    70    83 -2364
                                                                                                                      104
 BL17   TERP + O3 = #.052 HO2 + #.067 xHO2 + #.585 OH + &                        
        #.126 xMECO3 + #.149 xRCO3 + #.875 RO2C + &                              
        #.203 RO2XC + #.203 zRNO3 + #.166 CO + #.019 xCO + &                     
        #.045 CO2 + #.079 HCHO + #.15 xHCHO + #.22 xRCHO + &                     
        #.165 xACET + #.004 MEK + #.107 HCOOH + &                                
        #.043 RCOOH + #.001 xGLY + #.002 xMGLY + &                               
        #.055 xBACL + #.001 xMACR + #.001 xIPRD + &                              
        #.409 PROD2 + #.545 yR6OOH + #3.526 XC                             1.405E+00   1.56  -1.000      38     1     0 = -2365    88 -2333    56                
                                                                                                                    -2366     9 -2367    60
                                                                                                                    -2280    61 -2368    96
                                                                                                                    -2369    97 -2369    81
                                                                                                                    -2351    89 -2128    64
                                                                                                                    -2370    90 -2371    93
                                                                                                                    -2013    65 -2372    67
                                                                                                                    -2373    68 -2044   110
                                                                                                                    -2112   111 -2374   102
                                                                                                                    -2185    71 -2205    72
                                                                                                                    -2173    73 -2185    77
                                                                                                                    -2185    79 -2342   114
                                                                                                                    -2375    83 -2376   104
 BL18   TERP + NO3 = #.162 xHO2 + #.421 xNO2 + &                                 
        #.019 xRCO3 + #1.509 RO2C + #.397 RO2XC + &                              
        #.397 zRNO3 + #.01 xCO + #.017 xHCHO + &                                 
        #.001 xCCHO + #.509 xRCHO + #.175 xACET + &                              
        #.001 xMGLY + #.003 xMACR + #.001 xMVK + &                               
        #.002 xIPRD + #.163 xRNO3 + yR6OOH + #.416 XN + &                        
        #4.473 XC                                                          1.879E+03  -0.97  -1.000      38    10     0 = -2330    56 -2377    58                
                                                                                                                    -2128    61 -2378    96
                                                                                                                    -2379    97 -2379    81
                                                                                                                    -2225    64 -2067    65
                                                                                                                    -2185    66 -2034    67
                                                                                                                    -2221    68 -2185    72
                                                                                                                    -2142    77 -2185    78
                                                                                                                    -2205    79 -2143    80
                                                                                                                       83 -2148    18 -2380
                                                                                                                      104
 BL19   TERP + O3P = #.147 RCHO + #.853 PROD2 + #4.441 XC                  5.446E+04   0.00  -1.000      38    51     0 = -2072   109 -2381   114                
                                                                                                                    -2382   104
 BL20   BENZENE + OH = BENZENE                                             1.791E+03   0.00  -1.000      30     9     0 =    30                                  
 BL21   TRIMBNZ + OH = BNZSOMG + TRIMBNZ                                   4.800E+04   0.00  -1.000      31     9     0 =    46    31                            
 BL22   TOLUENE + OH = TOLSOMG + TOLUENE                                   8.264E+03   0.00  -1.000      32     9     0 =    44    32                            
 BL23   XYLENE  + OH = XYLSOMG + XYLENE                                    1.996E+04   0.00  -1.000      33     9     0 =    45    33                            
 BL24   TERPENE + OH = TRPSOMG + TERPENE                                   7.677E+04   0.00  -1.000      35     9     0 =    48    35                            
 BL25   SESQTRP + OH = TRPSOMG + SESQTRP                                   4.404E+05   0.00  -1.000      36     9     0 =    48    36                            
 BL26   ISOPRENE+ OH = ISPSOMG + ISOPRENE                                  1.468E+05   0.00  -1.000      34     9     0 =    47    34                            
 BL27   IVOC    + OH = IVOSOMG + IVOC                                      8.705E+04   0.00  -1.000      42     9     0 =    49    42                            
 BL28   SVOC    + OH = SVOSOMG + SVOC                                      1.527E+05   0.00  -1.000      43     9     0 =    50    43                            
 CI01   CL2 + HV = #2 CL                                                        PHOT. = CL2             145     8     0 = -2001   146                            
 CI02   CL + NO + M = CLNO + M                                             2.729E-03   0.00  -3.800     146     2     5 =    24     5                            
 CI03   CLNO + HV = CL + NO                                                     PHOT. = CLNO-06          24     8     0 =   146     2                            
 CI04   CL + NO2 = CLONO                                                  FALLOFF F= 0.600, N= 1.000    146     3     0 =    25                                  
                                                                      K0:  4.669E-02   0.00  -4.000
                                                                      KI:  1.468E+05   0.00  -2.000
 CI05   CL + NO2 = CLNO2                                                  FALLOFF F= 0.600, N= 1.000    146     3     0 =    26                                  
                                                                      K0:  6.464E-03   0.00  -4.000
                                                                      KI:  1.468E+05   0.00  -2.000
 CI06   CLONO + HV = CL + NO2                                                   PHOT. = CLONO            25     8     0 =   146     3                            
 CI07   CLNO2 + HV = CL + NO2                                                   PHOT. = CLNO2            26     8     0 =   146     3                            
 CI08   CL + HO2 = HCL + O2                                                5.050E+04   0.00  -1.560     146    88     0 =   147     4                            
 CI09   CL + HO2 = CLO + OH                                                1.381E+04   0.00   1.100     146    88     0 =   148     9                            
 CI10   CL + O3 = CLO + O2                                                 4.110E+04   0.50  -1.000     146     1     0 =   148     4                            
 CI11   CL + NO3 = CLO + NO2                                               3.523E+04   0.00  -1.000     146    10     0 =   148     3                            
 CI12   CLO + NO = CL + NO2                                                9.101E+03  -0.59  -1.000     148     2     0 =   146     3                            
 CI13   CLO + NO2 = CLONO2                                                FALLOFF F= 0.600, N= 1.000    148     3     0 =    27                                  
                                                                      K0:  6.464E-03   0.00  -5.400
                                                                      KI:  2.202E+04   0.00  -2.900
 CI14   CLONO2 + HV = CLO + NO2                                                 PHOT. = CLONO2-1         27     8     0 =   148     3                            
 CI15   CLONO2 + HV = CL + NO3                                                  PHOT. = CLONO2-2         27     8     0 =   146    10                            
 CI16   CLONO2 = CLO + NO2                                                FALLOFF F= 0.600, N= 1.000     27     0     0 =   148     3                            
                                                                      K0:  6.576E+10  24.90  -2.000
                                                                      KI:  2.226E+17  24.90   3.500
 CI17   CL + CLONO2 = CL2 + NO3                                            9.101E+03  -0.29  -1.000     146    27     0 =   145    10                            
 CI18   CLO + HO2 = HOCL + O2                                              3.229E+03  -0.68  -1.000     148    88     0 =   149     4                            
 CI19   HOCL + HV = OH + CL                                                     PHOT. = HOCL-06         149     8     0 =     9   146                            
 CI20   CLO + CLO = #.29 CL2 + #1.42 CL + O2                               1.835E+04   3.89  -1.000     148   148     0 = -2249   145 -2383   146                
                                                                                                                        4
 CI21   OH + HCL = H2O + CL                                                2.495E+03   0.46  -1.000       9   147     0 =     6   146                            
 CI22   CL + H2 = HCL + HO2                                                5.725E+04   4.59  -1.000     146     7     0 =   147    88                            
 CP01   HCHO + CL = HCL + HO2 + CO                                         1.189E+05   0.06  -1.000      93   146     0 =   147    88    89                      
 CP02   CCHO + CL = HCL + MECO3                                            1.174E+05   0.00  -1.000     108   146     0 =   147    98                            
 CP03   MEOH + CL = HCL + HCHO + HO2                                       8.073E+04   0.00  -1.000      95   146     0 =   147    93    88                      
 CP04   RCHO + CL = HCL + #.9 RCO3 + #.1 RO2C + &                                
        #.1 xCCHO + &                                                            
        #.1 xCO + #.1 xHO2 + #.1 yROOH                                     1.806E+05   0.00  -1.000     109   146     0 =   147 -2015   100 -2014                
                                                                                                                       96 -2014    66 -2014
                                                                                                                       64 -2014    56 -2014
                                                                                                                       82
 CP05   ACET + CL = HCL + RO2C + xHCHO + xMECO3 + yROOH                    1.130E+05   1.99  -1.000     106   146     0 =   147    96    65    60                
                                                                                                                       82
 CP06   MEK + CL = HCL + #.975 RO2C + #.039 RO2XC + &                            
        #.039 zRNO3 + #.84 xHO2 + #.085 xMECO3 + &                               
        #.036 xRCO3 + #.065 xHCHO + #.07 xCCHO + &                               
        #.84 xRCHO + yROOH + #.763 XC                                      5.284E+04   0.00  -1.000     110   146     0 =   147 -2158    96 -2027                
                                                                                                                       97 -2027    81 -2051
                                                                                                                       56 -2062    60 -2224
                                                                                                                       61 -2102    65 -2087
                                                                                                                       66 -2051    67    82
                                                                                                                    -2384   104
 CP07   RNO3 + CL = HCL + #.038 NO2 + #.055 HO2 + &                              
        #1.282 RO2C + #.202 RO2XC + #.202 zRNO3 + &                              
        #.009 RCHO + #.018 MEK + #.012 PROD2 + &                                 
        #.055 RNO3 + #.159 xNO2 + #.547 xHO2 + &                                 
        #.045 xHCHO + #.300 xCCHO + #.020 xRCHO + &                              
        #.003 xACET + #.041 xMEK + #.046 xPROD2 + &                              
        #.547 xRNO3 + #.908 yR6OOH + #.201 XN + #-1 XC + #.851 XC          2.818E+05   0.00  -1.000      19   146     0 =   147 -2385     3 -2173                
                                                                                                                       88 -2386    96 -2349
                                                                                                                       97 -2349    81 -2313
                                                                                                                      109 -2335   110 -2046
                                                                                                                      114 -2173    19 -2387
                                                                                                                       58 -2388    56 -2370
                                                                                                                       65 -2033    66 -2056
                                                                                                                       67 -2142    68 -2055
                                                                                                                       69 -2389    70 -2388
                                                                                                                       80 -2390    83 -2108
                                                                                                                       18 -2018   104 -2391
                                                                                                                      104
 CP08   PROD2 + CL = HCL + #.314 HO2 + #.680 RO2C + &                            
        #.116 RO2XC + #.116 zRNO3 + #.198 RCHO + &                               
        #.116 PROD2 + #.541 xHO2 + #.007 xMECO3 + &                              
        #.022 xRCO3 + #.237 xHCHO + #.109 xCCHO + &                              
        #.591 xRCHO + #.051 xMEK + #.040 xPROD2 + &                              
        #.686 yR6OOH + #1.262 XC                                           2.936E+05   0.00  -1.000     114   146     0 =   147 -2270    88 -2392                
                                                                                                                       96 -2268    97 -2268
                                                                                                                       81 -2393   109 -2268
                                                                                                                      114 -2394    56 -2136
                                                                                                                       60 -2306    61 -2395
                                                                                                                       65 -2396    66 -2237
                                                                                                                       67 -2188    69 -2084
                                                                                                                       70 -2397    83 -2398
                                                                                                                      104
 CP09   GLY + CL = HCL + #.63 HO2 + #1.26 CO + #.37 RCO3 + &                     
        #-.37 XC                                                           1.189E+05   0.06  -1.000     118   146     0 =   147 -2399    88 -2400                
                                                                                                                       89 -2248   100 -2401
                                                                                                                      104
 CP10   MGLY + CL = HCL + CO + MECO3                                       1.174E+05   0.00  -1.000     116   146     0 =   147    89    98                      
 CP11   CRES + CL = HCL + xHO2 + xBALD + yR6OOH                            9.101E+04   0.00  -1.000     107   146     0 =   147    56    74    83                
 CP12   BALD + CL = HCL + BZCO3                                            1.174E+05   0.00  -1.000     129   146     0 =   147   103                            
 CP13   ROOH + CL = HCL + #.414 OH + #.588 RO2C + &                              
        #.414 RCHO + #.104 xOH + #.482 xHO2 + &                                  
        #.106 xHCHO + #.104 xCCHO + #.197 xRCHO + &                              
        #.285 xMEK + #.586 yROOH + #-0.287 XC                              2.437E+05   0.00  -1.000     112   146     0 =   147 -2402     9 -2403                
                                                                                                                       96 -2402   109 -2303
                                                                                                                       57 -2352    56 -2283
                                                                                                                       65 -2303    66 -2404
                                                                                                                       67 -2178    69 -2405
                                                                                                                       82 -2406   104
 CP14   R6OOH + CL = HCL + #.145 OH + #1.078 RO2C + &                            
        #.117 RO2XC + #.117 zRNO3 + #.145 PROD2 + #.502 xOH + &                  
        #.237 xHO2 + #.186 xCCHO + #.676 xRCHO + &                               
        #.28 xPROD2 + #.855 yR6OOH + #.348 XC                              4.404E+05   0.00  -1.000     113   146     0 =   147 -2407     9 -2408                
                                                                                                                       96 -2123    97 -2123
                                                                                                                       81 -2407   114 -2409
                                                                                                                       57 -2395    56 -2410
                                                                                                                       66 -2411    67 -2412
                                                                                                                       70 -2413    83 -2414
                                                                                                                      104
 CP15   RAOOH + CL = #.404 HCL + #.139 OH + #.148 HO2 + &                        
        #.589 RO2C + #.124 RO2XC + #.124 zRNO3 + &                               
        #.074 PROD2 + #.147 MGLY + #.139 IPRD + &                                
        #.565 xHO2 + #.024 xOH + #.448 xRCHO + &                                 
        #.026 xGLY + #.030 xMEK + #.252 xMGLY + &                                
        #.073 xAFG1 + #.073 xAFG2 + #.713 yR6OOH + &                             
        #2.674 XC                                                          6.297E+05   0.00  -1.000     115   146     0 = -2415   147 -2068     9                
                                                                                                                    -2069    88 -2070    96
                                                                                                                    -2071    97 -2071    81
                                                                                                                    -2030   114 -2072   116
                                                                                                                    -2068   117 -2073    56
                                                                                                                    -2074    57 -2075    67
                                                                                                                    -2063    71 -2076    69
                                                                                                                    -2077    72 -2078    75
                                                                                                                    -2078    76 -2079    83
                                                                                                                    -2080   104
 CP16   MACR + CL = #.25 HCL + #.165 MACO3 + #.802 RO2C + &                      
        #.033 RO2XC + #.033 zRNO3 + #.802 xHO2 + &                               
        #.541 xCO + #.082 xIPRD + #.18 xCLCCHO + &                               
        #.541 xCLACET + #.835 yROOH + #.208 XC                             5.651E+05   0.00  -1.000     131   146     0 = -2008   147 -2373   105                
                                                                                                                    -2416    96 -2336    97
                                                                                                                    -2336    81 -2416    56
                                                                                                                    -2394    64 -2417    79
                                                                                                                    -2418    86 -2394    87
                                                                                                                    -2419    82 -2150   104
 CP17   MVK + CL = #1.283 RO2C + #.053 RO2XC + #.053 zRNO3 + &                   
        #.322 xHO2 + #.625 xMECO3 + #.947 xCLCCHO + &                            
        yROOH + #.538 XC                                                   3.406E+05   0.00  -1.000     132   146     0 = -2420    96 -2421    97                
                                                                                                                    -2421    81 -2345    56
                                                                                                                    -2422    60 -2423    86
                                                                                                                       82 -2424   104
 CP18   IPRD + CL = #.401 HCL + #.084 HO2 + #.154 MACO3 + &                      
        #.73 RO2C + #.051 RO2XC + #.051 zRNO3 + &                                
        #.042 AFG1 + #.042 AFG2 + #.712 xHO2 + #.498 xCO + &                     
        #.195 xHCHO + #.017 xMGLY + #.009 xAFG1 + &                              
        #.009 xAFG2 + #.115 xIPRD + #.14 xCLCCHO + &                             
        #.42 xCLACET + #.762 yR6OOH + #.709 XC                             6.048E+05   0.00  -1.000     117   146     0 = -2262   147 -2058    88                
                                                                                                                    -2344   105 -2425    96
                                                                                                                    -2188    97 -2188    81
                                                                                                                    -2339   119 -2339   120
                                                                                                                    -2426    56 -2179    64
                                                                                                                    -2120    65 -2067    72
                                                                                                                    -2313    75 -2313    76
                                                                                                                    -2085    79 -2180    86
                                                                                                                    -2427    87 -2428    83
                                                                                                                    -2348   104
 CP19   CLCCHO + HV = HO2 + CO + RO2C + xCL + xHCHO + &                          
        yROOH                                                                   PHOT. = CLCCHO          150     8     0 =    88    89    96    85                
                                                                                                                       65    82
 CP20   CLCCHO + OH = RCO3 + #-1 XC                                        4.550E+03   0.00  -1.000     150     9     0 =   100 -2018   104                      
 CP21   CLCCHO + CL = HCL + RCO3 + #-1 XC                                  1.894E+04   0.00  -1.000     150   146     0 =   147   100 -2018   104                
 CP22   CLACET + HV + #QY.CLAC = MECO3 + RO2C + xCL + xHCHO + yROOH             PHOT. = CLACET          151     8   -17 =    98    96    85    65                
                                                                                                                       82
 CP23   xCL + #RO2RO = CL                                                  1.000E+00   0.00   0.000      85    -5     0 =   146                                  
 CP24   xCL + #RO2XRO =                                                    1.000E+00   0.00   0.000      85    -6     0 =                                        
 CP25   xCLCCHO + #RO2RO = CLCCHO                                          1.000E+00   0.00   0.000      86    -5     0 =   150                                  
 CP26   xCLCCHO + #RO2XRO = #2 XC                                          1.000E+00   0.00   0.000      86    -6     0 = -2001   104                            
 CP27   xCLACET + #RO2RO = CLACET                                          1.000E+00   0.00   0.000      87    -5     0 =   151                                  
 CP28   xCLACET + #RO2XRO = #3 XC                                          1.000E+00   0.00   0.000      87    -6     0 = -2019   104                            
 CE01   CH4 + CL = HCL + MEO2                                              1.072E+04   2.54  -1.000     136   146     0 =   147    92                            
 CE02   ETHENE + CL = #2 RO2C + xHO2 + xHCHO + CLCHO                      FALLOFF F= 0.600, N= 1.000    137   146     0 = -2001    96    56    65                
                                                                                                                      152
                                                                      K0:  5.746E-01   0.00  -5.300
                                                                      KI:  4.550E+05   0.00  -2.000
 CE03   ISOPRN + CL = #.15 HCL + #1.168 RO2C + &                                 
        #.085 RO2XC + #.085 zRNO3 + #.738 xHO2 + &                               
        #.177 xCL + #.275 xHCHO + #.177 xMVK + &                                 
        #.671 xIPRD + #.067 xCLCCHO + yR6OOH + #.018 XC                    7.046E+05   0.00  -1.000      37   146     0 = -2013   147 -2429    96                
                                                                                                                    -2062    97 -2062    81
                                                                                                                    -2430    56 -2431    85
                                                                                                                    -2259    65 -2431    78
                                                                                                                    -2432    79 -2333    86
                                                                                                                       83 -2335   104
 CE04   ACETYLEN + CL = HO2 + CO + XC                                     FALLOFF F= 0.600, N= 1.000    138   146     0 =    88    89   104                      
                                                                      K0:  1.867E-01   0.00  -4.400
                                                                      KI:  3.229E+05   0.00  -1.000
 BC01   ALK1 + CL = xHO2 + RO2C + HCL + xCCHO + yROOH                      1.218E+05   0.20  -1.000     139   146     0 =    56    96   147    66                
                                                                                                                       82
 BC02   ALK2 + CL = #.97 xHO2 + #.97 RO2C + #.03 RO2XC + &                       
        #.03 zRNO3 + HCL + #.482 xRCHO + #.488 xACET + &                         
        yROOH + #-.09 XC                                                   1.761E+05  -0.08  -1.000     140   146     0 = -2433    56 -2433    96                
                                                                                                                    -2076    97 -2076    81
                                                                                                                      147 -2352    67 -2434
                                                                                                                       68    82 -2435   104
 BC03   ALK3 + CL = #.835 xHO2 + #.094 xTBUO + &                                 
        #1.361 RO2C + #.07 RO2XC + #.07 zRNO3 + HCL + &                          
        #.078 xHCHO + #.34 xCCHO + #.343 xRCHO + &                               
        #.075 xACET + #.253 xMEK + #.983 yROOH + &                               
        #.017 yR6OOH + #.18 XC                                             2.730E+05   0.00  -1.000     141   146     0 = -2419    56 -2436    63                
                                                                                                                    -2437    96 -2087    97
                                                                                                                    -2087    81   147 -2438
                                                                                                                       65 -2157    66 -2439
                                                                                                                       67 -2057    68 -2440
                                                                                                                       69 -2276    82 -2067
                                                                                                                       83 -2418   104
 BC04   ALK4 + CL = #.827 xHO2 + #.003 xMEO2 + &                                 
        #.004 xMECO3 + #1.737 RO2C + #.165 RO2XC + &                             
        #.165 zRNO3 + HCL + #.003 xCO + #.034 xHCHO + &                          
        #.287 xCCHO + #.412 xRCHO + #.247 xACET + &                              
        #.076 xMEK + #.13 xPROD2 + yR6OOH + #.327 XC                       3.861E+05   0.00  -1.000     142   146     0 = -2441    56 -2142    59                
                                                                                                                    -2044    60 -2442    96
                                                                                                                    -2373    97 -2373    81
                                                                                                                      147 -2142    64 -2048
                                                                                                                       65 -2443    66 -2444
                                                                                                                       67 -2445    68 -2324
                                                                                                                       69 -2083    70    83
                                                                                                                    -2446   104
 BC05   ALK5 + CL = #.647 xHO2 + #1.541 RO2C + &                                 
        #.352 RO2XC + #.352 zRNO3 + HCL + #.022 xHCHO + &                        
        #.08 xCCHO + #.258 xRCHO + #.044 xACET + &                               
        #.041 xMEK + #.378 xPROD2 + yR6OOH + #2.368 XC                     6.180E+05   0.00  -1.000      28   146     0 = -2286    56 -2447    96                
                                                                                                                    -2448    97 -2448    81
                                                                                                                      147 -2306    65 -2449
                                                                                                                       66 -2135    67 -2450
                                                                                                                       68 -2055    69 -2451
                                                                                                                       70    83 -2452   104
 BC06   OLE1 + CL = #.902 xHO2 + #1.42 RO2C + &                                  
        #.098 RO2XC + #.098 zRNO3 + #.308 HCL + &                                
        #.025 xHCHO + #.146 xCCHO + #.051 xRCHO + &                              
        #.188 xMACR + #.014 xMVK + #.027 xIPRD + &                               
        #.225 xCLCCHO + #.396 xCLACET + #.413 yROOH + &                          
        #.587 yR6OOH + #1.361 XC                                           5.211E+05   0.00  -1.000     143   146     0 = -2146    56 -2383    96                
                                                                                                                    -2147    97 -2147    81
                                                                                                                    -2453   147 -2159    65
                                                                                                                    -2284    66 -2188    67
                                                                                                                    -2454    77 -2347    78
                                                                                                                    -2325    79 -2455    86
                                                                                                                    -2238    87 -2299    82
                                                                                                                    -2300    83 -2437   104
 BC07   OLE2 + CL = #.447 xHO2 + #.448 xCL + #.001 xMEO2 + &                     
        #1.514 RO2C + #.104 RO2XC + #.104 zRNO3 + &                              
        #.263 HCL + #.228 xHCHO + #.361 xCCHO + &                                
        #.3 xRCHO + #.081 xACET + #.04 xMEK + &                                  
        #.049 xMACR + #.055 xMVK + #.179 xIPRD + &                               
        #.012 xCLCCHO + #.18 xCLACET + #.357 yROOH + &                           
        #.643 yR6OOH + #.247 XC                                            5.622E+05   0.00  -1.000     144   146     0 = -2456    56 -2075    85                
                                                                                                                    -2185    59 -2457    96
                                                                                                                    -2303    97 -2303    81
                                                                                                                    -2458   147 -2459    65
                                                                                                                    -2460    66 -2033    67
                                                                                                                    -2239    68 -2084    69
                                                                                                                    -2202    77 -2173    78
                                                                                                                    -2461    79 -2046    86
                                                                                                                    -2418    87 -2255    82
                                                                                                                    -2326    83 -2445   104
 BC08   ARO1 + CL = #.88 xHO2 + #.88 RO2C + #.12 RO2XC + &                       
        #.12 zRNO3 + #.671 xBALD + #.21 xPROD2 + #.323 XC                  1.468E+05   0.00  -1.000      40   146     0 = -2462    56 -2462    96                
                                                                                                                    -2247    97 -2247    81
                                                                                                                    -2432    74 -2090    70
                                                                                                                    -2463   104
 BC09   ARO2 + CL = #.842 xHO2 + #.842 RO2C + &                                  
        #.158 RO2XC + #.158 zRNO3 + #.618 xBALD + &                              
        #.224 xPROD2 + #2.382 XC                                           3.200E+05   0.00  -1.000      41   146     0 = -2464    56 -2464    96                
                                                                                                                    -2465    97 -2465    81
                                                                                                                    -2466    74 -2467    70
                                                                                                                    -2468   104
 BC10   TERP + CL = #.252 xHO2 + #.068 xCL + &                                   
        #.034 xMECO3 + #.05 xRCO3 + #.016 xMACO3 + &                             
        #2.258 RO2C + #.582 RO2XC + #.582 zRNO3 + &                              
        #.548 HCL + #.035 xCO + #.158 xHCHO + &                                  
        #.185 xRCHO + #.274 xACET + #.007 xGLY + &                               
        #.003 xBACL + #.003 xMVK + #.158 xIPRD + &                               
        #.006 xAFG1 + #.006 xAFG2 + &                                            
        #.109 xCLCCHO + yR6OOH + #3.543 XC                                 8.015E+05   0.00  -1.000      38   146     0 = -2077    56 -2353    85                
                                                                                                                    -2048    60 -2164    61
                                                                                                                    -2469    62 -2470    96
                                                                                                                    -2471    97 -2471    81
                                                                                                                    -2472   147 -2022    64
                                                                                                                    -2465    65 -2195    67
                                                                                                                    -2126    68 -2136    71
                                                                                                                    -2142    73 -2142    78
                                                                                                                    -2465    79 -2297    75
                                                                                                                    -2297    76 -2396    86
                                                                                                                       83 -2473   104
 D1       HCL = HCL                                                        1.000E-25   0.00   0.000     147     0     0 =   147                                  
 D2       NH3 = NH3                                                        1.000E-25   0.00   0.000     153     0     0 =   153                                  
 D3       GDENS = GDENS                                                    1.000E-25   0.00   0.000     154     0     0 =   154                                  
 D4       MMOLES = MMOLES                                                  1.000E-25   0.00   0.000     155     0     0 =   155                                  
 D5       XICBC = XICBC                                                    1.000E-25   0.00   0.000     156     0     0 =   156                                  
 D6       YICBC = YICBC                                                    1.000E-25   0.00   0.000     157     0     0 =   157                                  
 D7       ZICBC = ZICBC                                                    1.000E-25   0.00   0.000     158     0     0 =   158                                  
 TS01   CHCL3 + OH = xCL + RO2C + yROOH + XC                               8.323E+02   1.00   1.000     159     9     0 =    85    96    82   104                
 TS02   CHCL3 + CL = xCL + RO2C + HCL + yROOH + XC                         4.844E+03   1.97  -1.000     159   146     0 =    85    96   147    82                
                                                                                                                      104
0 366 REACTIONS,  159 SPECIES,   17 VARIABLE COEFFICIENTS, 473 CONSTANT COEFFICIENTS, AND  35 PHOTOLYSIS FILES.


These are the som reactions

 BE12   BENZ + OH = #0.116 OH + #0.290 RO2C + #0.290 xHO2 + #0.024 RO2XC + &
          #0.024 zRNO3 + #0.570 HO2 + #0.570 CRES + #0.116 AFG3 + #0.290 &
          xGLY + #0.029 xAFG1 + #0.261 xAFG2 + #0.314 yRAOOH + #-1.000 &
          XC + #0.024 XC + #1.000 #0.514 BNZSOMG_06_01 + #1.000 #0.087 &
          BNZSOMG_06_02 + #1.000 #0.207 BNZSOMG_06_03 + #1.000 #0.192 BNZSOMG_06_04
 S1.1   BNZSOMG_02_01 + OH = #0.001 #0.514 BNZSOMG_02_02 + #0.001 #0.087 &
          BNZSOMG_02_03 + #0.001 #0.399 BNZSOMG_02_04 + #0.999 #1.000 BNZSOMG_01_01 + &
          #0.999 #1.000 BNZSOMG_01_02
 S2.1   BNZSOMG_02_02 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #2.000 &
          BNZSOMG_01_02
 S3.1   BNZSOMG_02_03 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_01_02
 S4.1   BNZSOMG_02_04 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_01_02
 S5.1   BNZSOMG_03_01 + OH = #0.001 #0.514 BNZSOMG_03_02 + #0.001 #0.087 &
          BNZSOMG_03_03 + #0.001 #0.207 BNZSOMG_03_04 + #0.001 #0.192 BNZSOMG_03_05 + &
          #0.999 #0.500 BNZSOMG_01_01 + #0.999 #0.500 BNZSOMG_02_01 + #0.999 &
          #0.500 BNZSOMG_01_02 + #0.999 #0.500 BNZSOMG_02_02
 S6.1   BNZSOMG_03_02 + OH = #0.000 #0.514 BNZSOMG_03_03 + #0.000 #0.087 &
          BNZSOMG_03_04 + #0.000 #0.207 BNZSOMG_03_05 + #0.000 #0.192 BNZSOMG_03_06 + &
          #1.000 #0.500 BNZSOMG_01_01 + #1.000 #0.000 BNZSOMG_02_01 + #1.000 &
          #0.500 BNZSOMG_01_02 + #1.000 #0.500 BNZSOMG_02_02 + #1.000 #0.500 &
          BNZSOMG_02_03
 S7.1   BNZSOMG_03_03 + OH = #1.000 #0.500 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.500 BNZSOMG_01_02 + #1.000 #0.000 BNZSOMG_02_02 + &
          #1.000 #0.500 BNZSOMG_02_03 + #1.000 #0.500 BNZSOMG_02_04
 S8.1   BNZSOMG_03_04 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #1.000 BNZSOMG_01_02 + #1.000 #0.000 BNZSOMG_02_02 + &
          #1.000 #0.000 BNZSOMG_02_03 + #1.000 #1.000 BNZSOMG_02_04
 S9.1   BNZSOMG_03_05 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_01_02 + #1.000 #0.000 BNZSOMG_02_02 + &
          #1.000 #0.000 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_02_04
 S10.1  BNZSOMG_03_06 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_01_02 + #1.000 #0.000 BNZSOMG_02_02 + &
          #1.000 #0.000 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_02_04
 S11.1  BNZSOMG_04_01 + OH = #0.001 #0.514 BNZSOMG_04_02 + #0.001 #0.087 &
          BNZSOMG_04_03 + #0.001 #0.207 BNZSOMG_04_04 + #0.001 #0.192 BNZSOMG_04_05 + &
          #0.999 #0.333 BNZSOMG_01_01 + #0.999 #0.333 BNZSOMG_02_01 + #0.999 &
          #0.333 BNZSOMG_03_01 + #0.999 #0.333 BNZSOMG_01_02 + #0.999 #0.333 &
          BNZSOMG_02_02 + #0.999 #0.333 BNZSOMG_03_02
 S12.1  BNZSOMG_04_02 + OH = #0.001 #0.514 BNZSOMG_04_03 + #0.001 #0.087 &
          BNZSOMG_04_04 + #0.001 #0.207 BNZSOMG_04_05 + #0.001 #0.192 BNZSOMG_04_06 + &
          #0.999 #0.286 BNZSOMG_01_01 + #0.999 #0.286 BNZSOMG_02_01 + #0.999 &
          #0.000 BNZSOMG_03_01 + #0.999 #0.286 BNZSOMG_01_02 + #0.999 #0.286 &
          BNZSOMG_02_02 + #0.999 #0.286 BNZSOMG_03_02 + #0.999 #0.286 BNZSOMG_02_03 + &
          #0.999 #0.286 BNZSOMG_03_03
 S13.1  BNZSOMG_04_03 + OH = #0.000 #0.514 BNZSOMG_04_04 + #0.000 #0.087 &
          BNZSOMG_04_05 + #0.000 #0.207 BNZSOMG_04_06 + #0.000 #0.192 BNZSOMG_04_07 + &
          #1.000 #0.250 BNZSOMG_01_01 + #1.000 #0.250 BNZSOMG_02_01 + #1.000 &
          #0.000 BNZSOMG_03_01 + #1.000 #0.250 BNZSOMG_01_02 + #1.000 #0.250 &
          BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 #0.250 BNZSOMG_02_03 + &
          #1.000 #0.250 BNZSOMG_03_03 + #1.000 #0.250 BNZSOMG_02_04 + #1.000 &
          #0.250 BNZSOMG_03_04
 S14.1  BNZSOMG_04_04 + OH = #1.000 #0.286 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.286 BNZSOMG_01_02 + &
          #1.000 #0.286 BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 &
          #0.286 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.286 &
          BNZSOMG_02_04 + #1.000 #0.286 BNZSOMG_03_04 + #1.000 #0.286 BNZSOMG_03_05
 S15.1  BNZSOMG_04_05 + OH = #1.000 #0.333 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.333 BNZSOMG_01_02 + &
          #1.000 #0.000 BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 &
          #0.333 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.333 &
          BNZSOMG_02_04 + #1.000 #0.000 BNZSOMG_03_04 + #1.000 #0.333 BNZSOMG_03_05 + &
          #1.000 #0.333 BNZSOMG_03_06
 S16.1  BNZSOMG_04_06 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.667 BNZSOMG_01_02 + &
          #1.000 #0.000 BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 &
          #0.000 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.667 &
          BNZSOMG_02_04 + #1.000 #0.000 BNZSOMG_03_04 + #1.000 #0.000 BNZSOMG_03_05 + &
          #1.000 #0.667 BNZSOMG_03_06
 S17.1  BNZSOMG_04_07 + OH = #1.000 #0.000 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_01_02 + &
          #1.000 #0.000 BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 &
          #0.000 BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.000 &
          BNZSOMG_02_04 + #1.000 #0.000 BNZSOMG_03_04 + #1.000 #0.000 BNZSOMG_03_05 + &
          #1.000 #0.000 BNZSOMG_03_06
 S18.1  BNZSOMG_05_01 + OH = #0.002 #0.514 BNZSOMG_05_02 + #0.002 #0.087 &
          BNZSOMG_05_03 + #0.002 #0.207 BNZSOMG_05_04 + #0.002 #0.192 BNZSOMG_05_05 + &
          #0.998 #0.250 BNZSOMG_01_01 + #0.998 #0.250 BNZSOMG_02_01 + #0.998 &
          #0.250 BNZSOMG_03_01 + #0.998 #0.250 BNZSOMG_04_01 + #0.998 #0.250 &
          BNZSOMG_01_02 + #0.998 #0.250 BNZSOMG_02_02 + #0.998 #0.250 BNZSOMG_03_02 + &
          #0.998 #0.250 BNZSOMG_04_02
 S19.1  BNZSOMG_05_02 + OH = #0.001 #0.514 BNZSOMG_05_03 + #0.001 #0.087 &
          BNZSOMG_05_04 + #0.001 #0.207 BNZSOMG_05_05 + #0.001 #0.192 BNZSOMG_05_06 + &
          #0.999 #0.182 BNZSOMG_01_01 + #0.999 #0.182 BNZSOMG_02_01 + #0.999 &
          #0.182 BNZSOMG_03_01 + #0.999 #0.000 BNZSOMG_04_01 + #0.999 #0.182 &
          BNZSOMG_01_02 + #0.999 #0.182 BNZSOMG_02_02 + #0.999 #0.182 BNZSOMG_03_02 + &
          #0.999 #0.182 BNZSOMG_04_02 + #0.999 #0.182 BNZSOMG_02_03 + #0.999 &
          #0.182 BNZSOMG_03_03 + #0.999 #0.182 BNZSOMG_04_03
 S20.1  BNZSOMG_05_03 + OH = #0.001 #0.514 BNZSOMG_05_04 + #0.001 #0.087 &
          BNZSOMG_05_05 + #0.001 #0.207 BNZSOMG_05_06 + #0.001 #0.192 BNZSOMG_05_07 + &
          #0.999 #0.154 BNZSOMG_01_01 + #0.999 #0.154 BNZSOMG_02_01 + #0.999 &
          #0.154 BNZSOMG_03_01 + #0.999 #0.000 BNZSOMG_04_01 + #0.999 #0.154 &
          BNZSOMG_01_02 + #0.999 #0.154 BNZSOMG_02_02 + #0.999 #0.154 BNZSOMG_03_02 + &
          #0.999 #0.000 BNZSOMG_04_02 + #0.999 #0.154 BNZSOMG_02_03 + #0.999 &
          #0.154 BNZSOMG_03_03 + #0.999 #0.154 BNZSOMG_04_03 + #0.999 #0.154 &
          BNZSOMG_02_04 + #0.999 #0.154 BNZSOMG_03_04 + #0.999 #0.154 BNZSOMG_04_04
 S21.1  BNZSOMG_05_04 + OH = #0.000 #0.514 BNZSOMG_05_05 + #0.000 #0.087 &
          BNZSOMG_05_06 + #0.000 #0.399 BNZSOMG_05_07 + #1.000 #0.143 BNZSOMG_01_01 + &
          #1.000 #0.143 BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 &
          #0.000 BNZSOMG_04_01 + #1.000 #0.143 BNZSOMG_01_02 + #1.000 #0.143 &
          BNZSOMG_02_02 + #1.000 #0.143 BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + &
          #1.000 #0.143 BNZSOMG_02_03 + #1.000 #0.143 BNZSOMG_03_03 + #1.000 &
          #0.000 BNZSOMG_04_03 + #1.000 #0.143 BNZSOMG_02_04 + #1.000 #0.143 &
          BNZSOMG_03_04 + #1.000 #0.143 BNZSOMG_04_04 + #1.000 #0.143 BNZSOMG_03_05 + &
          #1.000 #0.143 BNZSOMG_04_05
 S22.1  BNZSOMG_05_05 + OH = #1.000 #0.143 BNZSOMG_01_01 + #1.000 #0.143 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.143 BNZSOMG_01_02 + #1.000 #0.143 BNZSOMG_02_02 + #1.000 &
          #0.000 BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + #1.000 #0.143 &
          BNZSOMG_02_03 + #1.000 #0.143 BNZSOMG_03_03 + #1.000 #0.000 BNZSOMG_04_03 + &
          #1.000 #0.143 BNZSOMG_02_04 + #1.000 #0.143 BNZSOMG_03_04 + #1.000 &
          #0.000 BNZSOMG_04_04 + #1.000 #0.143 BNZSOMG_03_05 + #1.000 #0.143 &
          BNZSOMG_04_05 + #1.000 #0.143 BNZSOMG_03_06 + #1.000 #0.143 BNZSOMG_04_06
 S23.1  BNZSOMG_05_06 + OH = #1.000 #0.154 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.154 BNZSOMG_01_02 + #1.000 #0.154 BNZSOMG_02_02 + #1.000 &
          #0.000 BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + #1.000 #0.154 &
          BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.000 BNZSOMG_04_03 + &
          #1.000 #0.154 BNZSOMG_02_04 + #1.000 #0.154 BNZSOMG_03_04 + #1.000 &
          #0.000 BNZSOMG_04_04 + #1.000 #0.154 BNZSOMG_03_05 + #1.000 #0.000 &
          BNZSOMG_04_05 + #1.000 #0.154 BNZSOMG_03_06 + #1.000 #0.154 BNZSOMG_04_06 + &
          #1.000 #0.154 BNZSOMG_04_07
 S24.1  BNZSOMG_05_07 + OH = #1.000 #0.200 BNZSOMG_01_01 + #1.000 #0.000 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.200 BNZSOMG_01_02 + #1.000 #0.000 BNZSOMG_02_02 + #1.000 &
          #0.000 BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + #1.000 #0.200 &
          BNZSOMG_02_03 + #1.000 #0.000 BNZSOMG_03_03 + #1.000 #0.000 BNZSOMG_04_03 + &
          #1.000 #0.200 BNZSOMG_02_04 + #1.000 #0.000 BNZSOMG_03_04 + #1.000 &
          #0.000 BNZSOMG_04_04 + #1.000 #0.200 BNZSOMG_03_05 + #1.000 #0.000 &
          BNZSOMG_04_05 + #1.000 #0.200 BNZSOMG_03_06 + #1.000 #0.000 BNZSOMG_04_06 + &
          #1.000 #0.200 BNZSOMG_04_07
 S25.1  BNZSOMG_06_01 + OH = #0.002 #0.514 BNZSOMG_06_02 + #0.002 #0.087 &
          BNZSOMG_06_03 + #0.002 #0.207 BNZSOMG_06_04 + #0.002 #0.192 BNZSOMG_06_05 + &
          #0.998 #0.200 BNZSOMG_01_01 + #0.998 #0.200 BNZSOMG_02_01 + #0.998 &
          #0.200 BNZSOMG_03_01 + #0.998 #0.200 BNZSOMG_04_01 + #0.998 #0.200 &
          BNZSOMG_05_01 + #0.998 #0.200 BNZSOMG_01_02 + #0.998 #0.200 BNZSOMG_02_02 + &
          #0.998 #0.200 BNZSOMG_03_02 + #0.998 #0.200 BNZSOMG_04_02 + #0.998 &
          #0.200 BNZSOMG_05_02
 S26.1  BNZSOMG_06_02 + OH = #0.001 #0.514 BNZSOMG_06_03 + #0.001 #0.087 &
          BNZSOMG_06_04 + #0.001 #0.207 BNZSOMG_06_05 + #0.001 #0.192 BNZSOMG_06_06 + &
          #0.999 #0.154 BNZSOMG_01_01 + #0.999 #0.154 BNZSOMG_02_01 + #0.999 &
          #0.154 BNZSOMG_03_01 + #0.999 #0.154 BNZSOMG_04_01 + #0.999 #0.000 &
          BNZSOMG_05_01 + #0.999 #0.154 BNZSOMG_01_02 + #0.999 #0.154 BNZSOMG_02_02 + &
          #0.999 #0.154 BNZSOMG_03_02 + #0.999 #0.154 BNZSOMG_04_02 + #0.999 &
          #0.154 BNZSOMG_05_02 + #0.999 #0.154 BNZSOMG_02_03 + #0.999 #0.154 &
          BNZSOMG_03_03 + #0.999 #0.154 BNZSOMG_04_03 + #0.999 #0.154 BNZSOMG_05_03
 S27.1  BNZSOMG_06_03 + OH = #0.001 #0.514 BNZSOMG_06_04 + #0.001 #0.087 &
          BNZSOMG_06_05 + #0.001 #0.207 BNZSOMG_06_06 + #0.001 #0.192 BNZSOMG_06_07 + &
          #0.999 #0.125 BNZSOMG_01_01 + #0.999 #0.125 BNZSOMG_02_01 + #0.999 &
          #0.125 BNZSOMG_03_01 + #0.999 #0.125 BNZSOMG_04_01 + #0.999 #0.000 &
          BNZSOMG_05_01 + #0.999 #0.125 BNZSOMG_01_02 + #0.999 #0.125 BNZSOMG_02_02 + &
          #0.999 #0.125 BNZSOMG_03_02 + #0.999 #0.125 BNZSOMG_04_02 + #0.999 &
          #0.000 BNZSOMG_05_02 + #0.999 #0.125 BNZSOMG_02_03 + #0.999 #0.125 &
          BNZSOMG_03_03 + #0.999 #0.125 BNZSOMG_04_03 + #0.999 #0.125 BNZSOMG_05_03 + &
          #0.999 #0.125 BNZSOMG_02_04 + #0.999 #0.125 BNZSOMG_03_04 + #0.999 &
          #0.125 BNZSOMG_04_04 + #0.999 #0.125 BNZSOMG_05_04
 S28.1  BNZSOMG_06_04 + OH = #0.000 #0.514 BNZSOMG_06_05 + #0.000 #0.087 &
          BNZSOMG_06_06 + #0.000 #0.399 BNZSOMG_06_07 + #1.000 #0.111 BNZSOMG_01_01 + &
          #1.000 #0.111 BNZSOMG_02_01 + #1.000 #0.111 BNZSOMG_03_01 + #1.000 &
          #0.000 BNZSOMG_04_01 + #1.000 #0.000 BNZSOMG_05_01 + #1.000 #0.111 &
          BNZSOMG_01_02 + #1.000 #0.111 BNZSOMG_02_02 + #1.000 #0.111 BNZSOMG_03_02 + &
          #1.000 #0.111 BNZSOMG_04_02 + #1.000 #0.000 BNZSOMG_05_02 + #1.000 &
          #0.111 BNZSOMG_02_03 + #1.000 #0.111 BNZSOMG_03_03 + #1.000 #0.111 &
          BNZSOMG_04_03 + #1.000 #0.000 BNZSOMG_05_03 + #1.000 #0.111 BNZSOMG_02_04 + &
          #1.000 #0.111 BNZSOMG_03_04 + #1.000 #0.111 BNZSOMG_04_04 + #1.000 &
          #0.111 BNZSOMG_05_04 + #1.000 #0.111 BNZSOMG_03_05 + #1.000 #0.111 &
          BNZSOMG_04_05 + #1.000 #0.111 BNZSOMG_05_05
 S29.1  BNZSOMG_06_05 + OH = #0.000 #0.514 BNZSOMG_06_06 + #0.000 #0.486 &
          BNZSOMG_06_07 + #1.000 #0.105 BNZSOMG_01_01 + #1.000 #0.105 BNZSOMG_02_01 + &
          #1.000 #0.105 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + #1.000 &
          #0.000 BNZSOMG_05_01 + #1.000 #0.105 BNZSOMG_01_02 + #1.000 #0.105 &
          BNZSOMG_02_02 + #1.000 #0.105 BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + &
          #1.000 #0.000 BNZSOMG_05_02 + #1.000 #0.105 BNZSOMG_02_03 + #1.000 &
          #0.105 BNZSOMG_03_03 + #1.000 #0.105 BNZSOMG_04_03 + #1.000 #0.000 &
          BNZSOMG_05_03 + #1.000 #0.105 BNZSOMG_02_04 + #1.000 #0.105 BNZSOMG_03_04 + &
          #1.000 #0.105 BNZSOMG_04_04 + #1.000 #0.000 BNZSOMG_05_04 + #1.000 &
          #0.105 BNZSOMG_03_05 + #1.000 #0.105 BNZSOMG_04_05 + #1.000 #0.105 &
          BNZSOMG_05_05 + #1.000 #0.105 BNZSOMG_03_06 + #1.000 #0.105 BNZSOMG_04_06 + &
          #1.000 #0.105 BNZSOMG_05_06
 S30.1  BNZSOMG_06_06 + OH = #1.000 #0.105 BNZSOMG_01_01 + #1.000 #0.105 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.000 BNZSOMG_05_01 + #1.000 #0.105 BNZSOMG_01_02 + #1.000 &
          #0.105 BNZSOMG_02_02 + #1.000 #0.105 BNZSOMG_03_02 + #1.000 #0.000 &
          BNZSOMG_04_02 + #1.000 #0.000 BNZSOMG_05_02 + #1.000 #0.105 BNZSOMG_02_03 + &
          #1.000 #0.105 BNZSOMG_03_03 + #1.000 #0.000 BNZSOMG_04_03 + #1.000 &
          #0.000 BNZSOMG_05_03 + #1.000 #0.105 BNZSOMG_02_04 + #1.000 #0.105 &
          BNZSOMG_03_04 + #1.000 #0.105 BNZSOMG_04_04 + #1.000 #0.000 BNZSOMG_05_04 + &
          #1.000 #0.105 BNZSOMG_03_05 + #1.000 #0.105 BNZSOMG_04_05 + #1.000 &
          #0.000 BNZSOMG_05_05 + #1.000 #0.105 BNZSOMG_03_06 + #1.000 #0.105 &
          BNZSOMG_04_06 + #1.000 #0.105 BNZSOMG_05_06 + #1.000 #0.105 BNZSOMG_04_07 + &
          #1.000 #0.105 BNZSOMG_05_07
 S31.1  BNZSOMG_06_07 + OH = #1.000 #0.125 BNZSOMG_01_01 + #1.000 #0.125 &
          BNZSOMG_02_01 + #1.000 #0.000 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.000 BNZSOMG_05_01 + #1.000 #0.125 BNZSOMG_01_02 + #1.000 &
          #0.125 BNZSOMG_02_02 + #1.000 #0.000 BNZSOMG_03_02 + #1.000 #0.000 &
          BNZSOMG_04_02 + #1.000 #0.000 BNZSOMG_05_02 + #1.000 #0.125 BNZSOMG_02_03 + &
          #1.000 #0.125 BNZSOMG_03_03 + #1.000 #0.000 BNZSOMG_04_03 + #1.000 &
          #0.000 BNZSOMG_05_03 + #1.000 #0.125 BNZSOMG_02_04 + #1.000 #0.125 &
          BNZSOMG_03_04 + #1.000 #0.000 BNZSOMG_04_04 + #1.000 #0.000 BNZSOMG_05_04 + &
          #1.000 #0.125 BNZSOMG_03_05 + #1.000 #0.125 BNZSOMG_04_05 + #1.000 &
          #0.000 BNZSOMG_05_05 + #1.000 #0.125 BNZSOMG_03_06 + #1.000 #0.125 &
          BNZSOMG_04_06 + #1.000 #0.000 BNZSOMG_05_06 + #1.000 #0.125 BNZSOMG_04_07 + &
          #1.000 #0.125 BNZSOMG_05_07
 S32.1  BNZSOMG_07_01 + OH = #0.002 #0.514 BNZSOMG_07_02 + #0.002 #0.087 &
          BNZSOMG_07_03 + #0.002 #0.207 BNZSOMG_07_04 + #0.002 #0.192 BNZSOMG_07_05 + &
          #0.998 #0.167 BNZSOMG_01_01 + #0.998 #0.167 BNZSOMG_02_01 + #0.998 &
          #0.167 BNZSOMG_03_01 + #0.998 #0.167 BNZSOMG_04_01 + #0.998 #0.167 &
          BNZSOMG_05_01 + #0.998 #0.167 BNZSOMG_06_01 + #0.998 #0.167 BNZSOMG_01_02 + &
          #0.998 #0.167 BNZSOMG_02_02 + #0.998 #0.167 BNZSOMG_03_02 + #0.998 &
          #0.167 BNZSOMG_04_02 + #0.998 #0.167 BNZSOMG_05_02 + #0.998 #0.167 &
          BNZSOMG_06_02
 S33.1  BNZSOMG_07_02 + OH = #0.001 #0.514 BNZSOMG_07_03 + #0.001 #0.087 &
          BNZSOMG_07_04 + #0.001 #0.207 BNZSOMG_07_05 + #0.001 #0.192 BNZSOMG_07_06 + &
          #0.999 #0.125 BNZSOMG_01_01 + #0.999 #0.125 BNZSOMG_02_01 + #0.999 &
          #0.125 BNZSOMG_03_01 + #0.999 #0.125 BNZSOMG_04_01 + #0.999 #0.125 &
          BNZSOMG_05_01 + #0.999 #0.000 BNZSOMG_06_01 + #0.999 #0.125 BNZSOMG_01_02 + &
          #0.999 #0.125 BNZSOMG_02_02 + #0.999 #0.125 BNZSOMG_03_02 + #0.999 &
          #0.125 BNZSOMG_04_02 + #0.999 #0.125 BNZSOMG_05_02 + #0.999 #0.125 &
          BNZSOMG_06_02 + #0.999 #0.125 BNZSOMG_02_03 + #0.999 #0.125 BNZSOMG_03_03 + &
          #0.999 #0.125 BNZSOMG_04_03 + #0.999 #0.125 BNZSOMG_05_03 + #0.999 &
          #0.125 BNZSOMG_06_03
 S34.1  BNZSOMG_07_03 + OH = #0.001 #0.514 BNZSOMG_07_04 + #0.001 #0.087 &
          BNZSOMG_07_05 + #0.001 #0.207 BNZSOMG_07_06 + #0.001 #0.192 BNZSOMG_07_07 + &
          #0.999 #0.100 BNZSOMG_01_01 + #0.999 #0.100 BNZSOMG_02_01 + #0.999 &
          #0.100 BNZSOMG_03_01 + #0.999 #0.100 BNZSOMG_04_01 + #0.999 #0.100 &
          BNZSOMG_05_01 + #0.999 #0.000 BNZSOMG_06_01 + #0.999 #0.100 BNZSOMG_01_02 + &
          #0.999 #0.100 BNZSOMG_02_02 + #0.999 #0.100 BNZSOMG_03_02 + #0.999 &
          #0.100 BNZSOMG_04_02 + #0.999 #0.100 BNZSOMG_05_02 + #0.999 #0.000 &
          BNZSOMG_06_02 + #0.999 #0.100 BNZSOMG_02_03 + #0.999 #0.100 BNZSOMG_03_03 + &
          #0.999 #0.100 BNZSOMG_04_03 + #0.999 #0.100 BNZSOMG_05_03 + #0.999 &
          #0.100 BNZSOMG_06_03 + #0.999 #0.100 BNZSOMG_02_04 + #0.999 #0.100 &
          BNZSOMG_03_04 + #0.999 #0.100 BNZSOMG_04_04 + #0.999 #0.100 BNZSOMG_05_04 + &
          #0.999 #0.100 BNZSOMG_06_04
 S35.1  BNZSOMG_07_04 + OH = #0.001 #0.514 BNZSOMG_07_05 + #0.001 #0.087 &
          BNZSOMG_07_06 + #0.001 #0.399 BNZSOMG_07_07 + #0.999 #0.091 BNZSOMG_01_01 + &
          #0.999 #0.091 BNZSOMG_02_01 + #0.999 #0.091 BNZSOMG_03_01 + #0.999 &
          #0.091 BNZSOMG_04_01 + #0.999 #0.000 BNZSOMG_05_01 + #0.999 #0.000 &
          BNZSOMG_06_01 + #0.999 #0.091 BNZSOMG_01_02 + #0.999 #0.091 BNZSOMG_02_02 + &
          #0.999 #0.091 BNZSOMG_03_02 + #0.999 #0.091 BNZSOMG_04_02 + #0.999 &
          #0.091 BNZSOMG_05_02 + #0.999 #0.000 BNZSOMG_06_02 + #0.999 #0.091 &
          BNZSOMG_02_03 + #0.999 #0.091 BNZSOMG_03_03 + #0.999 #0.091 BNZSOMG_04_03 + &
          #0.999 #0.091 BNZSOMG_05_03 + #0.999 #0.000 BNZSOMG_06_03 + #0.999 &
          #0.091 BNZSOMG_02_04 + #0.999 #0.091 BNZSOMG_03_04 + #0.999 #0.091 &
          BNZSOMG_04_04 + #0.999 #0.091 BNZSOMG_05_04 + #0.999 #0.091 BNZSOMG_06_04 + &
          #0.999 #0.091 BNZSOMG_03_05 + #0.999 #0.091 BNZSOMG_04_05 + #0.999 &
          #0.091 BNZSOMG_05_05 + #0.999 #0.091 BNZSOMG_06_05
 S36.1  BNZSOMG_07_05 + OH = #0.000 #0.514 BNZSOMG_07_06 + #0.000 #0.486 &
          BNZSOMG_07_07 + #1.000 #0.083 BNZSOMG_01_01 + #1.000 #0.083 BNZSOMG_02_01 + &
          #1.000 #0.083 BNZSOMG_03_01 + #1.000 #0.083 BNZSOMG_04_01 + #1.000 &
          #0.000 BNZSOMG_05_01 + #1.000 #0.000 BNZSOMG_06_01 + #1.000 #0.083 &
          BNZSOMG_01_02 + #1.000 #0.083 BNZSOMG_02_02 + #1.000 #0.083 BNZSOMG_03_02 + &
          #1.000 #0.083 BNZSOMG_04_02 + #1.000 #0.000 BNZSOMG_05_02 + #1.000 &
          #0.000 BNZSOMG_06_02 + #1.000 #0.083 BNZSOMG_02_03 + #1.000 #0.083 &
          BNZSOMG_03_03 + #1.000 #0.083 BNZSOMG_04_03 + #1.000 #0.083 BNZSOMG_05_03 + &
          #1.000 #0.000 BNZSOMG_06_03 + #1.000 #0.083 BNZSOMG_02_04 + #1.000 &
          #0.083 BNZSOMG_03_04 + #1.000 #0.083 BNZSOMG_04_04 + #1.000 #0.083 &
          BNZSOMG_05_04 + #1.000 #0.000 BNZSOMG_06_04 + #1.000 #0.083 BNZSOMG_03_05 + &
          #1.000 #0.083 BNZSOMG_04_05 + #1.000 #0.083 BNZSOMG_05_05 + #1.000 &
          #0.083 BNZSOMG_06_05 + #1.000 #0.083 BNZSOMG_03_06 + #1.000 #0.083 &
          BNZSOMG_04_06 + #1.000 #0.083 BNZSOMG_05_06 + #1.000 #0.083 BNZSOMG_06_06
 S37.1  BNZSOMG_07_06 + OH = #0.000 #1.000 BNZSOMG_07_07 + #1.000 #0.080 &
          BNZSOMG_01_01 + #1.000 #0.080 BNZSOMG_02_01 + #1.000 #0.080 BNZSOMG_03_01 + &
          #1.000 #0.000 BNZSOMG_04_01 + #1.000 #0.000 BNZSOMG_05_01 + #1.000 &
          #0.000 BNZSOMG_06_01 + #1.000 #0.080 BNZSOMG_01_02 + #1.000 #0.080 &
          BNZSOMG_02_02 + #1.000 #0.080 BNZSOMG_03_02 + #1.000 #0.080 BNZSOMG_04_02 + &
          #1.000 #0.000 BNZSOMG_05_02 + #1.000 #0.000 BNZSOMG_06_02 + #1.000 &
          #0.080 BNZSOMG_02_03 + #1.000 #0.080 BNZSOMG_03_03 + #1.000 #0.080 &
          BNZSOMG_04_03 + #1.000 #0.000 BNZSOMG_05_03 + #1.000 #0.000 BNZSOMG_06_03 + &
          #1.000 #0.080 BNZSOMG_02_04 + #1.000 #0.080 BNZSOMG_03_04 + #1.000 &
          #0.080 BNZSOMG_04_04 + #1.000 #0.080 BNZSOMG_05_04 + #1.000 #0.000 &
          BNZSOMG_06_04 + #1.000 #0.080 BNZSOMG_03_05 + #1.000 #0.080 BNZSOMG_04_05 + &
          #1.000 #0.080 BNZSOMG_05_05 + #1.000 #0.000 BNZSOMG_06_05 + #1.000 &
          #0.080 BNZSOMG_03_06 + #1.000 #0.080 BNZSOMG_04_06 + #1.000 #0.080 &
          BNZSOMG_05_06 + #1.000 #0.080 BNZSOMG_06_06 + #1.000 #0.080 BNZSOMG_04_07 + &
          #1.000 #0.080 BNZSOMG_05_07 + #1.000 #0.080 BNZSOMG_06_07
 S38.1  BNZSOMG_07_07 + OH = #1.000 #0.091 BNZSOMG_01_01 + #1.000 #0.091 &
          BNZSOMG_02_01 + #1.000 #0.091 BNZSOMG_03_01 + #1.000 #0.000 BNZSOMG_04_01 + &
          #1.000 #0.000 BNZSOMG_05_01 + #1.000 #0.000 BNZSOMG_06_01 + #1.000 &
          #0.091 BNZSOMG_01_02 + #1.000 #0.091 BNZSOMG_02_02 + #1.000 #0.091 &
          BNZSOMG_03_02 + #1.000 #0.000 BNZSOMG_04_02 + #1.000 #0.000 BNZSOMG_05_02 + &
          #1.000 #0.000 BNZSOMG_06_02 + #1.000 #0.091 BNZSOMG_02_03 + #1.000 &
          #0.091 BNZSOMG_03_03 + #1.000 #0.091 BNZSOMG_04_03 + #1.000 #0.000 &
          BNZSOMG_05_03 + #1.000 #0.000 BNZSOMG_06_03 + #1.000 #0.091 BNZSOMG_02_04 + &
          #1.000 #0.091 BNZSOMG_03_04 + #1.000 #0.091 BNZSOMG_04_04 + #1.000 &
          #0.000 BNZSOMG_05_04 + #1.000 #0.000 BNZSOMG_06_04 + #1.000 #0.091 &
          BNZSOMG_03_05 + #1.000 #0.091 BNZSOMG_04_05 + #1.000 #0.091 BNZSOMG_05_05 + &
          #1.000 #0.000 BNZSOMG_06_05 + #1.000 #0.091 BNZSOMG_03_06 + #1.000 &
          #0.091 BNZSOMG_04_06 + #1.000 #0.091 BNZSOMG_05_06 + #1.000 #0.000 &
          BNZSOMG_06_06 + #1.000 #0.091 BNZSOMG_04_07 + #1.000 #0.091 BNZSOMG_05_07 + &
          #1.000 #0.091 BNZSOMG_06_07
 S39.1  BNZSOMG_08_01 + OH = #0.002 #0.514 BNZSOMG_08_02 + #0.002 #0.087 &
          BNZSOMG_08_03 + #0.002 #0.207 BNZSOMG_08_04 + #0.002 #0.192 BNZSOMG_08_05 + &
          #0.998 #0.143 BNZSOMG_01_01 + #0.998 #0.143 BNZSOMG_02_01 + #0.998 &
          #0.143 BNZSOMG_03_01 + #0.998 #0.143 BNZSOMG_04_01 + #0.998 #0.143 &
          BNZSOMG_05_01 + #0.998 #0.143 BNZSOMG_06_01 + #0.998 #0.143 BNZSOMG_07_01 + &
          #0.998 #0.143 BNZSOMG_01_02 + #0.998 #0.143 BNZSOMG_02_02 + #0.998 &
          #0.143 BNZSOMG_03_02 + #0.998 #0.143 BNZSOMG_04_02 + #0.998 #0.143 &
          BNZSOMG_05_02 + #0.998 #0.143 BNZSOMG_06_02 + #0.998 #0.143 BNZSOMG_07_02
 S40.1  BNZSOMG_08_02 + OH = #0.001 #0.514 BNZSOMG_08_03 + #0.001 #0.087 &
          BNZSOMG_08_04 + #0.001 #0.207 BNZSOMG_08_05 + #0.001 #0.192 BNZSOMG_08_06 + &
          #0.999 #0.105 BNZSOMG_01_01 + #0.999 #0.105 BNZSOMG_02_01 + #0.999 &
          #0.105 BNZSOMG_03_01 + #0.999 #0.105 BNZSOMG_04_01 + #0.999 #0.105 &
          BNZSOMG_05_01 + #0.999 #0.105 BNZSOMG_06_01 + #0.999 #0.000 BNZSOMG_07_01 + &
          #0.999 #0.105 BNZSOMG_01_02 + #0.999 #0.105 BNZSOMG_02_02 + #0.999 &
          #0.105 BNZSOMG_03_02 + #0.999 #0.105 BNZSOMG_04_02 + #0.999 #0.105 &
          BNZSOMG_05_02 + #0.999 #0.105 BNZSOMG_06_02 + #0.999 #0.105 BNZSOMG_07_02 + &
          #0.999 #0.105 BNZSOMG_02_03 + #0.999 #0.105 BNZSOMG_03_03 + #0.999 &
          #0.105 BNZSOMG_04_03 + #0.999 #0.105 BNZSOMG_05_03 + #0.999 #0.105 &
          BNZSOMG_06_03 + #0.999 #0.105 BNZSOMG_07_03
 S41.1  BNZSOMG_08_03 + OH = #0.001 #0.514 BNZSOMG_08_04 + #0.001 #0.087 &
          BNZSOMG_08_05 + #0.001 #0.207 BNZSOMG_08_06 + #0.001 #0.192 BNZSOMG_08_07 + &
          #0.999 #0.083 BNZSOMG_01_01 + #0.999 #0.083 BNZSOMG_02_01 + #0.999 &
          #0.083 BNZSOMG_03_01 + #0.999 #0.083 BNZSOMG_04_01 + #0.999 #0.083 &
          BNZSOMG_05_01 + #0.999 #0.083 BNZSOMG_06_01 + #0.999 #0.000 BNZSOMG_07_01 + &
          #0.999 #0.083 BNZSOMG_01_02 + #0.999 #0.083 BNZSOMG_02_02 + #0.999 &
          #0.083 BNZSOMG_03_02 + #0.999 #0.083 BNZSOMG_04_02 + #0.999 #0.083 &
          BNZSOMG_05_02 + #0.999 #0.083 BNZSOMG_06_02 + #0.999 #0.000 BNZSOMG_07_02 + &
          #0.999 #0.083 BNZSOMG_02_03 + #0.999 #0.083 BNZSOMG_03_03 + #0.999 &
          #0.083 BNZSOMG_04_03 + #0.999 #0.083 BNZSOMG_05_03 + #0.999 #0.083 &
          BNZSOMG_06_03 + #0.999 #0.083 BNZSOMG_07_03 + #0.999 #0.083 BNZSOMG_02_04 + &
          #0.999 #0.083 BNZSOMG_03_04 + #0.999 #0.083 BNZSOMG_04_04 + #0.999 &
          #0.083 BNZSOMG_05_04 + #0.999 #0.083 BNZSOMG_06_04 + #0.999 #0.083 &
          BNZSOMG_07_04
 S42.1  BNZSOMG_08_04 + OH = #0.001 #0.514 BNZSOMG_08_05 + #0.001 #0.087 &
          BNZSOMG_08_06 + #0.001 #0.399 BNZSOMG_08_07 + #0.999 #0.074 BNZSOMG_01_01 + &
          #0.999 #0.074 BNZSOMG_02_01 + #0.999 #0.074 BNZSOMG_03_01 + #0.999 &
          #0.074 BNZSOMG_04_01 + #0.999 #0.074 BNZSOMG_05_01 + #0.999 #0.000 &
          BNZSOMG_06_01 + #0.999 #0.000 BNZSOMG_07_01 + #0.999 #0.074 BNZSOMG_01_02 + &
          #0.999 #0.074 BNZSOMG_02_02 + #0.999 #0.074 BNZSOMG_03_02 + #0.999 &
          #0.074 BNZSOMG_04_02 + #0.999 #0.074 BNZSOMG_05_02 + #0.999 #0.074 &
          BNZSOMG_06_02 + #0.999 #0.000 BNZSOMG_07_02 + #0.999 #0.074 BNZSOMG_02_03 + &
          #0.999 #0.074 BNZSOMG_03_03 + #0.999 #0.074 BNZSOMG_04_03 + #0.999 &
          #0.074 BNZSOMG_05_03 + #0.999 #0.074 BNZSOMG_06_03 + #0.999 #0.000 &
          BNZSOMG_07_03 + #0.999 #0.074 BNZSOMG_02_04 + #0.999 #0.074 BNZSOMG_03_04 + &
          #0.999 #0.074 BNZSOMG_04_04 + #0.999 #0.074 BNZSOMG_05_04 + #0.999 &
          #0.074 BNZSOMG_06_04 + #0.999 #0.074 BNZSOMG_07_04 + #0.999 #0.074 &
          BNZSOMG_03_05 + #0.999 #0.074 BNZSOMG_04_05 + #0.999 #0.074 BNZSOMG_05_05 + &
          #0.999 #0.074 BNZSOMG_06_05 + #0.999 #0.074 BNZSOMG_07_05
 S43.1  BNZSOMG_08_05 + OH = #0.000 #0.514 BNZSOMG_08_06 + #0.000 #0.486 &
          BNZSOMG_08_07 + #1.000 #0.067 BNZSOMG_01_01 + #1.000 #0.067 BNZSOMG_02_01 + &
          #1.000 #0.067 BNZSOMG_03_01 + #1.000 #0.067 BNZSOMG_04_01 + #1.000 &
          #0.067 BNZSOMG_05_01 + #1.000 #0.000 BNZSOMG_06_01 + #1.000 #0.000 &
          BNZSOMG_07_01 + #1.000 #0.067 BNZSOMG_01_02 + #1.000 #0.067 BNZSOMG_02_02 + &
          #1.000 #0.067 BNZSOMG_03_02 + #1.000 #0.067 BNZSOMG_04_02 + #1.000 &
          #0.067 BNZSOMG_05_02 + #1.000 #0.000 BNZSOMG_06_02 + #1.000 #0.000 &
          BNZSOMG_07_02 + #1.000 #0.067 BNZSOMG_02_03 + #1.000 #0.067 BNZSOMG_03_03 + &
          #1.000 #0.067 BNZSOMG_04_03 + #1.000 #0.067 BNZSOMG_05_03 + #1.000 &
          #0.067 BNZSOMG_06_03 + #1.000 #0.000 BNZSOMG_07_03 + #1.000 #0.067 &
          BNZSOMG_02_04 + #1.000 #0.067 BNZSOMG_03_04 + #1.000 #0.067 BNZSOMG_04_04 + &
          #1.000 #0.067 BNZSOMG_05_04 + #1.000 #0.067 BNZSOMG_06_04 + #1.000 &
          #0.000 BNZSOMG_07_04 + #1.000 #0.067 BNZSOMG_03_05 + #1.000 #0.067 &
          BNZSOMG_04_05 + #1.000 #0.067 BNZSOMG_05_05 + #1.000 #0.067 BNZSOMG_06_05 + &
          #1.000 #0.067 BNZSOMG_07_05 + #1.000 #0.067 BNZSOMG_03_06 + #1.000 &
          #0.067 BNZSOMG_04_06 + #1.000 #0.067 BNZSOMG_05_06 + #1.000 #0.067 &
          BNZSOMG_06_06 + #1.000 #0.067 BNZSOMG_07_06
 S44.1  BNZSOMG_08_06 + OH = #0.000 #1.000 BNZSOMG_08_07 + #1.000 #0.065 &
          BNZSOMG_01_01 + #1.000 #0.065 BNZSOMG_02_01 + #1.000 #0.065 BNZSOMG_03_01 + &
          #1.000 #0.065 BNZSOMG_04_01 + #1.000 #0.000 BNZSOMG_05_01 + #1.000 &
          #0.000 BNZSOMG_06_01 + #1.000 #0.000 BNZSOMG_07_01 + #1.000 #0.065 &
          BNZSOMG_01_02 + #1.000 #0.065 BNZSOMG_02_02 + #1.000 #0.065 BNZSOMG_03_02 + &
          #1.000 #0.065 BNZSOMG_04_02 + #1.000 #0.065 BNZSOMG_05_02 + #1.000 &
          #0.000 BNZSOMG_06_02 + #1.000 #0.000 BNZSOMG_07_02 + #1.000 #0.065 &
          BNZSOMG_02_03 + #1.000 #0.065 BNZSOMG_03_03 + #1.000 #0.065 BNZSOMG_04_03 + &
          #1.000 #0.065 BNZSOMG_05_03 + #1.000 #0.000 BNZSOMG_06_03 + #1.000 &
          #0.000 BNZSOMG_07_03 + #1.000 #0.065 BNZSOMG_02_04 + #1.000 #0.065 &
          BNZSOMG_03_04 + #1.000 #0.065 BNZSOMG_04_04 + #1.000 #0.065 BNZSOMG_05_04 + &
          #1.000 #0.065 BNZSOMG_06_04 + #1.000 #0.000 BNZSOMG_07_04 + #1.000 &
          #0.065 BNZSOMG_03_05 + #1.000 #0.065 BNZSOMG_04_05 + #1.000 #0.065 &
          BNZSOMG_05_05 + #1.000 #0.065 BNZSOMG_06_05 + #1.000 #0.000 BNZSOMG_07_05 + &
          #1.000 #0.065 BNZSOMG_03_06 + #1.000 #0.065 BNZSOMG_04_06 + #1.000 &
          #0.065 BNZSOMG_05_06 + #1.000 #0.065 BNZSOMG_06_06 + #1.000 #0.065 &
          BNZSOMG_07_06 + #1.000 #0.065 BNZSOMG_04_07 + #1.000 #0.065 BNZSOMG_05_07 + &
          #1.000 #0.065 BNZSOMG_06_07 + #1.000 #0.065 BNZSOMG_07_07
 S45.1  BNZSOMG_08_07 + OH = #1.000 #0.071 BNZSOMG_01_01 + #1.000 #0.071 &
          BNZSOMG_02_01 + #1.000 #0.071 BNZSOMG_03_01 + #1.000 #0.071 BNZSOMG_04_01 + &
          #1.000 #0.000 BNZSOMG_05_01 + #1.000 #0.000 BNZSOMG_06_01 + #1.000 &
          #0.000 BNZSOMG_07_01 + #1.000 #0.071 BNZSOMG_01_02 + #1.000 #0.071 &
          BNZSOMG_02_02 + #1.000 #0.071 BNZSOMG_03_02 + #1.000 #0.071 BNZSOMG_04_02 + &
          #1.000 #0.000 BNZSOMG_05_02 + #1.000 #0.000 BNZSOMG_06_02 + #1.000 &
          #0.000 BNZSOMG_07_02 + #1.000 #0.071 BNZSOMG_02_03 + #1.000 #0.071 &
          BNZSOMG_03_03 + #1.000 #0.071 BNZSOMG_04_03 + #1.000 #0.071 BNZSOMG_05_03 + &
          #1.000 #0.000 BNZSOMG_06_03 + #1.000 #0.000 BNZSOMG_07_03 + #1.000 &
          #0.071 BNZSOMG_02_04 + #1.000 #0.071 BNZSOMG_03_04 + #1.000 #0.071 &
          BNZSOMG_04_04 + #1.000 #0.071 BNZSOMG_05_04 + #1.000 #0.000 BNZSOMG_06_04 + &
          #1.000 #0.000 BNZSOMG_07_04 + #1.000 #0.071 BNZSOMG_03_05 + #1.000 &
          #0.071 BNZSOMG_04_05 + #1.000 #0.071 BNZSOMG_05_05 + #1.000 #0.071 &
          BNZSOMG_06_05 + #1.000 #0.000 BNZSOMG_07_05 + #1.000 #0.071 BNZSOMG_03_06 + &
          #1.000 #0.071 BNZSOMG_04_06 + #1.000 #0.071 BNZSOMG_05_06 + #1.000 &
          #0.071 BNZSOMG_06_06 + #1.000 #0.000 BNZSOMG_07_06 + #1.000 #0.071 &
          BNZSOMG_04_07 + #1.000 #0.071 BNZSOMG_05_07 + #1.000 #0.071 BNZSOMG_06_07 + &
          #1.000 #0.071 BNZSOMG_07_07
 S46.1  BNZSOMG_09_01 + OH = #0.002 #0.514 BNZSOMG_09_02 + #0.002 #0.087 &
          BNZSOMG_09_03 + #0.002 #0.207 BNZSOMG_09_04 + #0.002 #0.192 BNZSOMG_09_05 + &
          #0.998 #0.125 BNZSOMG_01_01 + #0.998 #0.125 BNZSOMG_02_01 + #0.998 &
          #0.125 BNZSOMG_03_01 + #0.998 #0.125 BNZSOMG_04_01 + #0.998 #0.125 &
          BNZSOMG_05_01 + #0.998 #0.125 BNZSOMG_06_01 + #0.998 #0.125 BNZSOMG_07_01 + &
          #0.998 #0.125 BNZSOMG_08_01 + #0.998 #0.125 BNZSOMG_01_02 + #0.998 &
          #0.125 BNZSOMG_02_02 + #0.998 #0.125 BNZSOMG_03_02 + #0.998 #0.125 &
          BNZSOMG_04_02 + #0.998 #0.125 BNZSOMG_05_02 + #0.998 #0.125 BNZSOMG_06_02 + &
          #0.998 #0.125 BNZSOMG_07_02 + #0.998 #0.125 BNZSOMG_08_02
 S47.1  BNZSOMG_09_02 + OH = #0.002 #0.514 BNZSOMG_09_03 + #0.002 #0.087 &
          BNZSOMG_09_04 + #0.002 #0.207 BNZSOMG_09_05 + #0.002 #0.192 BNZSOMG_09_06 + &
          #0.998 #0.087 BNZSOMG_01_01 + #0.998 #0.087 BNZSOMG_02_01 + #0.998 &
          #0.087 BNZSOMG_03_01 + #0.998 #0.087 BNZSOMG_04_01 + #0.998 #0.087 &
          BNZSOMG_05_01 + #0.998 #0.087 BNZSOMG_06_01 + #0.998 #0.087 BNZSOMG_07_01 + &
          #0.998 #0.000 BNZSOMG_08_01 + #0.998 #0.087 BNZSOMG_01_02 + #0.998 &
          #0.087 BNZSOMG_02_02 + #0.998 #0.087 BNZSOMG_03_02 + #0.998 #0.087 &
          BNZSOMG_04_02 + #0.998 #0.087 BNZSOMG_05_02 + #0.998 #0.087 BNZSOMG_06_02 + &
          #0.998 #0.087 BNZSOMG_07_02 + #0.998 #0.087 BNZSOMG_08_02 + #0.998 &
          #0.087 BNZSOMG_02_03 + #0.998 #0.087 BNZSOMG_03_03 + #0.998 #0.087 &
          BNZSOMG_04_03 + #0.998 #0.087 BNZSOMG_05_03 + #0.998 #0.087 BNZSOMG_06_03 + &
          #0.998 #0.087 BNZSOMG_07_03 + #0.998 #0.087 BNZSOMG_08_03
 S48.1  BNZSOMG_09_03 + OH = #0.001 #0.514 BNZSOMG_09_04 + #0.001 #0.087 &
          BNZSOMG_09_05 + #0.001 #0.207 BNZSOMG_09_06 + #0.001 #0.192 BNZSOMG_09_07 + &
          #0.999 #0.069 BNZSOMG_01_01 + #0.999 #0.069 BNZSOMG_02_01 + #0.999 &
          #0.069 BNZSOMG_03_01 + #0.999 #0.069 BNZSOMG_04_01 + #0.999 #0.069 &
          BNZSOMG_05_01 + #0.999 #0.069 BNZSOMG_06_01 + #0.999 #0.069 BNZSOMG_07_01 + &
          #0.999 #0.000 BNZSOMG_08_01 + #0.999 #0.069 BNZSOMG_01_02 + #0.999 &
          #0.069 BNZSOMG_02_02 + #0.999 #0.069 BNZSOMG_03_02 + #0.999 #0.069 &
          BNZSOMG_04_02 + #0.999 #0.069 BNZSOMG_05_02 + #0.999 #0.069 BNZSOMG_06_02 + &
          #0.999 #0.069 BNZSOMG_07_02 + #0.999 #0.000 BNZSOMG_08_02 + #0.999 &
          #0.069 BNZSOMG_02_03 + #0.999 #0.069 BNZSOMG_03_03 + #0.999 #0.069 &
          BNZSOMG_04_03 + #0.999 #0.069 BNZSOMG_05_03 + #0.999 #0.069 BNZSOMG_06_03 + &
          #0.999 #0.069 BNZSOMG_07_03 + #0.999 #0.069 BNZSOMG_08_03 + #0.999 &
          #0.069 BNZSOMG_02_04 + #0.999 #0.069 BNZSOMG_03_04 + #0.999 #0.069 &
          BNZSOMG_04_04 + #0.999 #0.069 BNZSOMG_05_04 + #0.999 #0.069 BNZSOMG_06_04 + &
          #0.999 #0.069 BNZSOMG_07_04 + #0.999 #0.069 BNZSOMG_08_04
 S49.1  BNZSOMG_09_04 + OH = #0.001 #0.514 BNZSOMG_09_05 + #0.001 #0.087 &
          BNZSOMG_09_06 + #0.001 #0.399 BNZSOMG_09_07 + #0.999 #0.059 BNZSOMG_01_01 + &
          #0.999 #0.059 BNZSOMG_02_01 + #0.999 #0.059 BNZSOMG_03_01 + #0.999 &
          #0.059 BNZSOMG_04_01 + #0.999 #0.059 BNZSOMG_05_01 + #0.999 #0.059 &
          BNZSOMG_06_01 + #0.999 #0.000 BNZSOMG_07_01 + #0.999 #0.000 BNZSOMG_08_01 + &
          #0.999 #0.059 BNZSOMG_01_02 + #0.999 #0.059 BNZSOMG_02_02 + #0.999 &
          #0.059 BNZSOMG_03_02 + #0.999 #0.059 BNZSOMG_04_02 + #0.999 #0.059 &
          BNZSOMG_05_02 + #0.999 #0.059 BNZSOMG_06_02 + #0.999 #0.059 BNZSOMG_07_02 + &
          #0.999 #0.000 BNZSOMG_08_02 + #0.999 #0.059 BNZSOMG_02_03 + #0.999 &
          #0.059 BNZSOMG_03_03 + #0.999 #0.059 BNZSOMG_04_03 + #0.999 #0.059 &
          BNZSOMG_05_03 + #0.999 #0.059 BNZSOMG_06_03 + #0.999 #0.059 BNZSOMG_07_03 + &
          #0.999 #0.000 BNZSOMG_08_03 + #0.999 #0.059 BNZSOMG_02_04 + #0.999 &
          #0.059 BNZSOMG_03_04 + #0.999 #0.059 BNZSOMG_04_04 + #0.999 #0.059 &
          BNZSOMG_05_04 + #0.999 #0.059 BNZSOMG_06_04 + #0.999 #0.059 BNZSOMG_07_04 + &
          #0.999 #0.059 BNZSOMG_08_04 + #0.999 #0.059 BNZSOMG_03_05 + #0.999 &
          #0.059 BNZSOMG_04_05 + #0.999 #0.059 BNZSOMG_05_05 + #0.999 #0.059 &
          BNZSOMG_06_05 + #0.999 #0.059 BNZSOMG_07_05 + #0.999 #0.059 BNZSOMG_08_05
 S50.1  BNZSOMG_09_05 + OH = #0.001 #0.514 BNZSOMG_09_06 + #0.001 #0.486 &
          BNZSOMG_09_07 + #0.999 #0.053 BNZSOMG_01_01 + #0.999 #0.053 BNZSOMG_02_01 + &
          #0.999 #0.053 BNZSOMG_03_01 + #0.999 #0.053 BNZSOMG_04_01 + #0.999 &
          #0.053 BNZSOMG_05_01 + #0.999 #0.053 BNZSOMG_06_01 + #0.999 #0.000 &
          BNZSOMG_07_01 + #0.999 #0.000 BNZSOMG_08_01 + #0.999 #0.053 BNZSOMG_01_02 + &
          #0.999 #0.053 BNZSOMG_02_02 + #0.999 #0.053 BNZSOMG_03_02 + #0.999 &
          #0.053 BNZSOMG_04_02 + #0.999 #0.053 BNZSOMG_05_02 + #0.999 #0.053 &
          BNZSOMG_06_02 + #0.999 #0.000 BNZSOMG_07_02 + #0.999 #0.000 BNZSOMG_08_02 + &
          #0.999 #0.053 BNZSOMG_02_03 + #0.999 #0.053 BNZSOMG_03_03 + #0.999 &
          #0.053 BNZSOMG_04_03 + #0.999 #0.053 BNZSOMG_05_03 + #0.999 #0.053 &
          BNZSOMG_06_03 + #0.999 #0.053 BNZSOMG_07_03 + #0.999 #0.000 BNZSOMG_08_03 + &
          #0.999 #0.053 BNZSOMG_02_04 + #0.999 #0.053 BNZSOMG_03_04 + #0.999 &
          #0.053 BNZSOMG_04_04 + #0.999 #0.053 BNZSOMG_05_04 + #0.999 #0.053 &
          BNZSOMG_06_04 + #0.999 #0.053 BNZSOMG_07_04 + #0.999 #0.000 BNZSOMG_08_04 + &
          #0.999 #0.053 BNZSOMG_03_05 + #0.999 #0.053 BNZSOMG_04_05 + #0.999 &
          #0.053 BNZSOMG_05_05 + #0.999 #0.053 BNZSOMG_06_05 + #0.999 #0.053 &
          BNZSOMG_07_05 + #0.999 #0.053 BNZSOMG_08_05 + #0.999 #0.053 BNZSOMG_03_06 + &
          #0.999 #0.053 BNZSOMG_04_06 + #0.999 #0.053 BNZSOMG_05_06 + #0.999 &
          #0.053 BNZSOMG_06_06 + #0.999 #0.053 BNZSOMG_07_06 + #0.999 #0.053 &
          BNZSOMG_08_06
 S51.1  BNZSOMG_09_06 + OH = #0.000 #1.000 BNZSOMG_09_07 + #1.000 #0.049 &
          BNZSOMG_01_01 + #1.000 #0.049 BNZSOMG_02_01 + #1.000 #0.049 BNZSOMG_03_01 + &
          #1.000 #0.049 BNZSOMG_04_01 + #1.000 #0.049 BNZSOMG_05_01 + #1.000 &
          #0.000 BNZSOMG_06_01 + #1.000 #0.000 BNZSOMG_07_01 + #1.000 #0.000 &
          BNZSOMG_08_01 + #1.000 #0.049 BNZSOMG_01_02 + #1.000 #0.049 BNZSOMG_02_02 + &
          #1.000 #0.049 BNZSOMG_03_02 + #1.000 #0.049 BNZSOMG_04_02 + #1.000 &
          #0.049 BNZSOMG_05_02 + #1.000 #0.049 BNZSOMG_06_02 + #1.000 #0.000 &
          BNZSOMG_07_02 + #1.000 #0.000 BNZSOMG_08_02 + #1.000 #0.049 BNZSOMG_02_03 + &
          #1.000 #0.049 BNZSOMG_03_03 + #1.000 #0.049 BNZSOMG_04_03 + #1.000 &
          #0.049 BNZSOMG_05_03 + #1.000 #0.049 BNZSOMG_06_03 + #1.000 #0.000 &
          BNZSOMG_07_03 + #1.000 #0.000 BNZSOMG_08_03 + #1.000 #0.049 BNZSOMG_02_04 + &
          #1.000 #0.049 BNZSOMG_03_04 + #1.000 #0.049 BNZSOMG_04_04 + #1.000 &
          #0.049 BNZSOMG_05_04 + #1.000 #0.049 BNZSOMG_06_04 + #1.000 #0.049 &
          BNZSOMG_07_04 + #1.000 #0.000 BNZSOMG_08_04 + #1.000 #0.049 BNZSOMG_03_05 + &
          #1.000 #0.049 BNZSOMG_04_05 + #1.000 #0.049 BNZSOMG_05_05 + #1.000 &
          #0.049 BNZSOMG_06_05 + #1.000 #0.049 BNZSOMG_07_05 + #1.000 #0.000 &
          BNZSOMG_08_05 + #1.000 #0.049 BNZSOMG_03_06 + #1.000 #0.049 BNZSOMG_04_06 + &
          #1.000 #0.049 BNZSOMG_05_06 + #1.000 #0.049 BNZSOMG_06_06 + #1.000 &
          #0.049 BNZSOMG_07_06 + #1.000 #0.049 BNZSOMG_08_06 + #1.000 #0.049 &
          BNZSOMG_04_07 + #1.000 #0.049 BNZSOMG_05_07 + #1.000 #0.049 BNZSOMG_06_07 + &
          #1.000 #0.049 BNZSOMG_07_07 + #1.000 #0.049 BNZSOMG_08_07
 S52.1  BNZSOMG_09_07 + OH = #1.000 #0.053 BNZSOMG_01_01 + #1.000 #0.053 &
          BNZSOMG_02_01 + #1.000 #0.053 BNZSOMG_03_01 + #1.000 #0.053 BNZSOMG_04_01 + &
          #1.000 #0.053 BNZSOMG_05_01 + #1.000 #0.000 BNZSOMG_06_01 + #1.000 &
          #0.000 BNZSOMG_07_01 + #1.000 #0.000 BNZSOMG_08_01 + #1.000 #0.053 &
          BNZSOMG_01_02 + #1.000 #0.053 BNZSOMG_02_02 + #1.000 #0.053 BNZSOMG_03_02 + &
          #1.000 #0.053 BNZSOMG_04_02 + #1.000 #0.053 BNZSOMG_05_02 + #1.000 &
          #0.000 BNZSOMG_06_02 + #1.000 #0.000 BNZSOMG_07_02 + #1.000 #0.000 &
          BNZSOMG_08_02 + #1.000 #0.053 BNZSOMG_02_03 + #1.000 #0.053 BNZSOMG_03_03 + &
          #1.000 #0.053 BNZSOMG_04_03 + #1.000 #0.053 BNZSOMG_05_03 + #1.000 &
          #0.053 BNZSOMG_06_03 + #1.000 #0.000 BNZSOMG_07_03 + #1.000 #0.000 &
          BNZSOMG_08_03 + #1.000 #0.053 BNZSOMG_02_04 + #1.000 #0.053 BNZSOMG_03_04 + &
          #1.000 #0.053 BNZSOMG_04_04 + #1.000 #0.053 BNZSOMG_05_04 + #1.000 &
          #0.053 BNZSOMG_06_04 + #1.000 #0.000 BNZSOMG_07_04 + #1.000 #0.000 &
          BNZSOMG_08_04 + #1.000 #0.053 BNZSOMG_03_05 + #1.000 #0.053 BNZSOMG_04_05 + &
          #1.000 #0.053 BNZSOMG_05_05 + #1.000 #0.053 BNZSOMG_06_05 + #1.000 &
          #0.053 BNZSOMG_07_05 + #1.000 #0.000 BNZSOMG_08_05 + #1.000 #0.053 &
          BNZSOMG_03_06 + #1.000 #0.053 BNZSOMG_04_06 + #1.000 #0.053 BNZSOMG_05_06 + &
          #1.000 #0.053 BNZSOMG_06_06 + #1.000 #0.053 BNZSOMG_07_06 + #1.000 &
          #0.000 BNZSOMG_08_06 + #1.000 #0.053 BNZSOMG_04_07 + #1.000 #0.053 &
          BNZSOMG_05_07 + #1.000 #0.053 BNZSOMG_06_07 + #1.000 #0.053 BNZSOMG_07_07 + &
          #1.000 #0.053 BNZSOMG_08_07
 BL21   TRIMBNZ + OH = #1.000 #0.514 BNZSOMG_09_01 + #1.000 #0.087 BNZSOMG_09_02 + &
          #1.000 #0.207 BNZSOMG_09_03 + #1.000 #0.192 BNZSOMG_09_04 + TRIMBNZ
 BL22   TOLUENE + OH = #1.000 #0.694 TOLSOMG_07_01 + #1.000 #0.031 TOLSOMG_07_02 + &
          #1.000 #0.002 TOLSOMG_07_03 + #1.000 #0.273 TOLSOMG_07_04 + TOLUENE
 S53.2  TOLSOMG_02_01 + OH = #0.001 #0.694 TOLSOMG_02_02 + #0.001 #0.031 &
          TOLSOMG_02_03 + #0.001 #0.275 TOLSOMG_02_04 + #0.999 #1.000 TOLSOMG_01_01 + &
          #0.999 #1.000 TOLSOMG_01_02
 S54.2  TOLSOMG_02_02 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #2.000 &
          TOLSOMG_01_02
 S55.2  TOLSOMG_02_03 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_01_02
 S56.2  TOLSOMG_02_04 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_01_02
 S57.2  TOLSOMG_03_01 + OH = #0.001 #0.694 TOLSOMG_03_02 + #0.001 #0.031 &
          TOLSOMG_03_03 + #0.001 #0.002 TOLSOMG_03_04 + #0.001 #0.273 TOLSOMG_03_05 + &
          #0.999 #0.500 TOLSOMG_01_01 + #0.999 #0.500 TOLSOMG_02_01 + #0.999 &
          #0.500 TOLSOMG_01_02 + #0.999 #0.500 TOLSOMG_02_02
 S58.2  TOLSOMG_03_02 + OH = #0.000 #0.694 TOLSOMG_03_03 + #0.000 #0.031 &
          TOLSOMG_03_04 + #0.000 #0.002 TOLSOMG_03_05 + #0.000 #0.273 TOLSOMG_03_06 + &
          #1.000 #0.500 TOLSOMG_01_01 + #1.000 #0.000 TOLSOMG_02_01 + #1.000 &
          #0.500 TOLSOMG_01_02 + #1.000 #0.500 TOLSOMG_02_02 + #1.000 #0.500 &
          TOLSOMG_02_03
 S59.2  TOLSOMG_03_03 + OH = #1.000 #0.500 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.500 TOLSOMG_01_02 + #1.000 #0.000 TOLSOMG_02_02 + &
          #1.000 #0.500 TOLSOMG_02_03 + #1.000 #0.500 TOLSOMG_02_04
 S60.2  TOLSOMG_03_04 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #1.000 TOLSOMG_01_02 + #1.000 #0.000 TOLSOMG_02_02 + &
          #1.000 #0.000 TOLSOMG_02_03 + #1.000 #1.000 TOLSOMG_02_04
 S61.2  TOLSOMG_03_05 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_01_02 + #1.000 #0.000 TOLSOMG_02_02 + &
          #1.000 #0.000 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_02_04
 S62.2  TOLSOMG_03_06 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_01_02 + #1.000 #0.000 TOLSOMG_02_02 + &
          #1.000 #0.000 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_02_04
 S63.2  TOLSOMG_04_01 + OH = #0.001 #0.694 TOLSOMG_04_02 + #0.001 #0.031 &
          TOLSOMG_04_03 + #0.001 #0.002 TOLSOMG_04_04 + #0.001 #0.273 TOLSOMG_04_05 + &
          #0.999 #0.333 TOLSOMG_01_01 + #0.999 #0.333 TOLSOMG_02_01 + #0.999 &
          #0.333 TOLSOMG_03_01 + #0.999 #0.333 TOLSOMG_01_02 + #0.999 #0.333 &
          TOLSOMG_02_02 + #0.999 #0.333 TOLSOMG_03_02
 S64.2  TOLSOMG_04_02 + OH = #0.001 #0.694 TOLSOMG_04_03 + #0.001 #0.031 &
          TOLSOMG_04_04 + #0.001 #0.002 TOLSOMG_04_05 + #0.001 #0.273 TOLSOMG_04_06 + &
          #0.999 #0.286 TOLSOMG_01_01 + #0.999 #0.286 TOLSOMG_02_01 + #0.999 &
          #0.000 TOLSOMG_03_01 + #0.999 #0.286 TOLSOMG_01_02 + #0.999 #0.286 &
          TOLSOMG_02_02 + #0.999 #0.286 TOLSOMG_03_02 + #0.999 #0.286 TOLSOMG_02_03 + &
          #0.999 #0.286 TOLSOMG_03_03
 S65.2  TOLSOMG_04_03 + OH = #0.000 #0.694 TOLSOMG_04_04 + #0.000 #0.031 &
          TOLSOMG_04_05 + #0.000 #0.002 TOLSOMG_04_06 + #0.000 #0.273 TOLSOMG_04_07 + &
          #1.000 #0.250 TOLSOMG_01_01 + #1.000 #0.250 TOLSOMG_02_01 + #1.000 &
          #0.000 TOLSOMG_03_01 + #1.000 #0.250 TOLSOMG_01_02 + #1.000 #0.250 &
          TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 #0.250 TOLSOMG_02_03 + &
          #1.000 #0.250 TOLSOMG_03_03 + #1.000 #0.250 TOLSOMG_02_04 + #1.000 &
          #0.250 TOLSOMG_03_04
 S66.2  TOLSOMG_04_04 + OH = #1.000 #0.286 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.286 TOLSOMG_01_02 + &
          #1.000 #0.286 TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 &
          #0.286 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.286 &
          TOLSOMG_02_04 + #1.000 #0.286 TOLSOMG_03_04 + #1.000 #0.286 TOLSOMG_03_05
 S67.2  TOLSOMG_04_05 + OH = #1.000 #0.333 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.333 TOLSOMG_01_02 + &
          #1.000 #0.000 TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 &
          #0.333 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.333 &
          TOLSOMG_02_04 + #1.000 #0.000 TOLSOMG_03_04 + #1.000 #0.333 TOLSOMG_03_05 + &
          #1.000 #0.333 TOLSOMG_03_06
 S68.2  TOLSOMG_04_06 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.667 TOLSOMG_01_02 + &
          #1.000 #0.000 TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 &
          #0.000 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.667 &
          TOLSOMG_02_04 + #1.000 #0.000 TOLSOMG_03_04 + #1.000 #0.000 TOLSOMG_03_05 + &
          #1.000 #0.667 TOLSOMG_03_06
 S69.2  TOLSOMG_04_07 + OH = #1.000 #0.000 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_01_02 + &
          #1.000 #0.000 TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 &
          #0.000 TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.000 &
          TOLSOMG_02_04 + #1.000 #0.000 TOLSOMG_03_04 + #1.000 #0.000 TOLSOMG_03_05 + &
          #1.000 #0.000 TOLSOMG_03_06
 S70.2  TOLSOMG_05_01 + OH = #0.002 #0.694 TOLSOMG_05_02 + #0.002 #0.031 &
          TOLSOMG_05_03 + #0.002 #0.002 TOLSOMG_05_04 + #0.002 #0.273 TOLSOMG_05_05 + &
          #0.998 #0.250 TOLSOMG_01_01 + #0.998 #0.250 TOLSOMG_02_01 + #0.998 &
          #0.250 TOLSOMG_03_01 + #0.998 #0.250 TOLSOMG_04_01 + #0.998 #0.250 &
          TOLSOMG_01_02 + #0.998 #0.250 TOLSOMG_02_02 + #0.998 #0.250 TOLSOMG_03_02 + &
          #0.998 #0.250 TOLSOMG_04_02
 S71.2  TOLSOMG_05_02 + OH = #0.001 #0.694 TOLSOMG_05_03 + #0.001 #0.031 &
          TOLSOMG_05_04 + #0.001 #0.002 TOLSOMG_05_05 + #0.001 #0.273 TOLSOMG_05_06 + &
          #0.999 #0.182 TOLSOMG_01_01 + #0.999 #0.182 TOLSOMG_02_01 + #0.999 &
          #0.182 TOLSOMG_03_01 + #0.999 #0.000 TOLSOMG_04_01 + #0.999 #0.182 &
          TOLSOMG_01_02 + #0.999 #0.182 TOLSOMG_02_02 + #0.999 #0.182 TOLSOMG_03_02 + &
          #0.999 #0.182 TOLSOMG_04_02 + #0.999 #0.182 TOLSOMG_02_03 + #0.999 &
          #0.182 TOLSOMG_03_03 + #0.999 #0.182 TOLSOMG_04_03
 S72.2  TOLSOMG_05_03 + OH = #0.001 #0.694 TOLSOMG_05_04 + #0.001 #0.031 &
          TOLSOMG_05_05 + #0.001 #0.002 TOLSOMG_05_06 + #0.001 #0.273 TOLSOMG_05_07 + &
          #0.999 #0.154 TOLSOMG_01_01 + #0.999 #0.154 TOLSOMG_02_01 + #0.999 &
          #0.154 TOLSOMG_03_01 + #0.999 #0.000 TOLSOMG_04_01 + #0.999 #0.154 &
          TOLSOMG_01_02 + #0.999 #0.154 TOLSOMG_02_02 + #0.999 #0.154 TOLSOMG_03_02 + &
          #0.999 #0.000 TOLSOMG_04_02 + #0.999 #0.154 TOLSOMG_02_03 + #0.999 &
          #0.154 TOLSOMG_03_03 + #0.999 #0.154 TOLSOMG_04_03 + #0.999 #0.154 &
          TOLSOMG_02_04 + #0.999 #0.154 TOLSOMG_03_04 + #0.999 #0.154 TOLSOMG_04_04
 S73.2  TOLSOMG_05_04 + OH = #0.000 #0.694 TOLSOMG_05_05 + #0.000 #0.031 &
          TOLSOMG_05_06 + #0.000 #0.275 TOLSOMG_05_07 + #1.000 #0.143 TOLSOMG_01_01 + &
          #1.000 #0.143 TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 &
          #0.000 TOLSOMG_04_01 + #1.000 #0.143 TOLSOMG_01_02 + #1.000 #0.143 &
          TOLSOMG_02_02 + #1.000 #0.143 TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + &
          #1.000 #0.143 TOLSOMG_02_03 + #1.000 #0.143 TOLSOMG_03_03 + #1.000 &
          #0.000 TOLSOMG_04_03 + #1.000 #0.143 TOLSOMG_02_04 + #1.000 #0.143 &
          TOLSOMG_03_04 + #1.000 #0.143 TOLSOMG_04_04 + #1.000 #0.143 TOLSOMG_03_05 + &
          #1.000 #0.143 TOLSOMG_04_05
 S74.2  TOLSOMG_05_05 + OH = #1.000 #0.143 TOLSOMG_01_01 + #1.000 #0.143 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.143 TOLSOMG_01_02 + #1.000 #0.143 TOLSOMG_02_02 + #1.000 &
          #0.000 TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + #1.000 #0.143 &
          TOLSOMG_02_03 + #1.000 #0.143 TOLSOMG_03_03 + #1.000 #0.000 TOLSOMG_04_03 + &
          #1.000 #0.143 TOLSOMG_02_04 + #1.000 #0.143 TOLSOMG_03_04 + #1.000 &
          #0.000 TOLSOMG_04_04 + #1.000 #0.143 TOLSOMG_03_05 + #1.000 #0.143 &
          TOLSOMG_04_05 + #1.000 #0.143 TOLSOMG_03_06 + #1.000 #0.143 TOLSOMG_04_06
 S75.2  TOLSOMG_05_06 + OH = #1.000 #0.154 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.154 TOLSOMG_01_02 + #1.000 #0.154 TOLSOMG_02_02 + #1.000 &
          #0.000 TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + #1.000 #0.154 &
          TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.000 TOLSOMG_04_03 + &
          #1.000 #0.154 TOLSOMG_02_04 + #1.000 #0.154 TOLSOMG_03_04 + #1.000 &
          #0.000 TOLSOMG_04_04 + #1.000 #0.154 TOLSOMG_03_05 + #1.000 #0.000 &
          TOLSOMG_04_05 + #1.000 #0.154 TOLSOMG_03_06 + #1.000 #0.154 TOLSOMG_04_06 + &
          #1.000 #0.154 TOLSOMG_04_07
 S76.2  TOLSOMG_05_07 + OH = #1.000 #0.200 TOLSOMG_01_01 + #1.000 #0.000 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.200 TOLSOMG_01_02 + #1.000 #0.000 TOLSOMG_02_02 + #1.000 &
          #0.000 TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + #1.000 #0.200 &
          TOLSOMG_02_03 + #1.000 #0.000 TOLSOMG_03_03 + #1.000 #0.000 TOLSOMG_04_03 + &
          #1.000 #0.200 TOLSOMG_02_04 + #1.000 #0.000 TOLSOMG_03_04 + #1.000 &
          #0.000 TOLSOMG_04_04 + #1.000 #0.200 TOLSOMG_03_05 + #1.000 #0.000 &
          TOLSOMG_04_05 + #1.000 #0.200 TOLSOMG_03_06 + #1.000 #0.000 TOLSOMG_04_06 + &
          #1.000 #0.200 TOLSOMG_04_07
 S77.2  TOLSOMG_06_01 + OH = #0.002 #0.694 TOLSOMG_06_02 + #0.002 #0.031 &
          TOLSOMG_06_03 + #0.002 #0.002 TOLSOMG_06_04 + #0.002 #0.273 TOLSOMG_06_05 + &
          #0.998 #0.200 TOLSOMG_01_01 + #0.998 #0.200 TOLSOMG_02_01 + #0.998 &
          #0.200 TOLSOMG_03_01 + #0.998 #0.200 TOLSOMG_04_01 + #0.998 #0.200 &
          TOLSOMG_05_01 + #0.998 #0.200 TOLSOMG_01_02 + #0.998 #0.200 TOLSOMG_02_02 + &
          #0.998 #0.200 TOLSOMG_03_02 + #0.998 #0.200 TOLSOMG_04_02 + #0.998 &
          #0.200 TOLSOMG_05_02
 S78.2  TOLSOMG_06_02 + OH = #0.001 #0.694 TOLSOMG_06_03 + #0.001 #0.031 &
          TOLSOMG_06_04 + #0.001 #0.002 TOLSOMG_06_05 + #0.001 #0.273 TOLSOMG_06_06 + &
          #0.999 #0.154 TOLSOMG_01_01 + #0.999 #0.154 TOLSOMG_02_01 + #0.999 &
          #0.154 TOLSOMG_03_01 + #0.999 #0.154 TOLSOMG_04_01 + #0.999 #0.000 &
          TOLSOMG_05_01 + #0.999 #0.154 TOLSOMG_01_02 + #0.999 #0.154 TOLSOMG_02_02 + &
          #0.999 #0.154 TOLSOMG_03_02 + #0.999 #0.154 TOLSOMG_04_02 + #0.999 &
          #0.154 TOLSOMG_05_02 + #0.999 #0.154 TOLSOMG_02_03 + #0.999 #0.154 &
          TOLSOMG_03_03 + #0.999 #0.154 TOLSOMG_04_03 + #0.999 #0.154 TOLSOMG_05_03
 S79.2  TOLSOMG_06_03 + OH = #0.001 #0.694 TOLSOMG_06_04 + #0.001 #0.031 &
          TOLSOMG_06_05 + #0.001 #0.002 TOLSOMG_06_06 + #0.001 #0.273 TOLSOMG_06_07 + &
          #0.999 #0.125 TOLSOMG_01_01 + #0.999 #0.125 TOLSOMG_02_01 + #0.999 &
          #0.125 TOLSOMG_03_01 + #0.999 #0.125 TOLSOMG_04_01 + #0.999 #0.000 &
          TOLSOMG_05_01 + #0.999 #0.125 TOLSOMG_01_02 + #0.999 #0.125 TOLSOMG_02_02 + &
          #0.999 #0.125 TOLSOMG_03_02 + #0.999 #0.125 TOLSOMG_04_02 + #0.999 &
          #0.000 TOLSOMG_05_02 + #0.999 #0.125 TOLSOMG_02_03 + #0.999 #0.125 &
          TOLSOMG_03_03 + #0.999 #0.125 TOLSOMG_04_03 + #0.999 #0.125 TOLSOMG_05_03 + &
          #0.999 #0.125 TOLSOMG_02_04 + #0.999 #0.125 TOLSOMG_03_04 + #0.999 &
          #0.125 TOLSOMG_04_04 + #0.999 #0.125 TOLSOMG_05_04
 S80.2  TOLSOMG_06_04 + OH = #0.000 #0.694 TOLSOMG_06_05 + #0.000 #0.031 &
          TOLSOMG_06_06 + #0.000 #0.275 TOLSOMG_06_07 + #1.000 #0.111 TOLSOMG_01_01 + &
          #1.000 #0.111 TOLSOMG_02_01 + #1.000 #0.111 TOLSOMG_03_01 + #1.000 &
          #0.000 TOLSOMG_04_01 + #1.000 #0.000 TOLSOMG_05_01 + #1.000 #0.111 &
          TOLSOMG_01_02 + #1.000 #0.111 TOLSOMG_02_02 + #1.000 #0.111 TOLSOMG_03_02 + &
          #1.000 #0.111 TOLSOMG_04_02 + #1.000 #0.000 TOLSOMG_05_02 + #1.000 &
          #0.111 TOLSOMG_02_03 + #1.000 #0.111 TOLSOMG_03_03 + #1.000 #0.111 &
          TOLSOMG_04_03 + #1.000 #0.000 TOLSOMG_05_03 + #1.000 #0.111 TOLSOMG_02_04 + &
          #1.000 #0.111 TOLSOMG_03_04 + #1.000 #0.111 TOLSOMG_04_04 + #1.000 &
          #0.111 TOLSOMG_05_04 + #1.000 #0.111 TOLSOMG_03_05 + #1.000 #0.111 &
          TOLSOMG_04_05 + #1.000 #0.111 TOLSOMG_05_05
 S81.2  TOLSOMG_06_05 + OH = #0.000 #0.694 TOLSOMG_06_06 + #0.000 #0.306 &
          TOLSOMG_06_07 + #1.000 #0.105 TOLSOMG_01_01 + #1.000 #0.105 TOLSOMG_02_01 + &
          #1.000 #0.105 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + #1.000 &
          #0.000 TOLSOMG_05_01 + #1.000 #0.105 TOLSOMG_01_02 + #1.000 #0.105 &
          TOLSOMG_02_02 + #1.000 #0.105 TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + &
          #1.000 #0.000 TOLSOMG_05_02 + #1.000 #0.105 TOLSOMG_02_03 + #1.000 &
          #0.105 TOLSOMG_03_03 + #1.000 #0.105 TOLSOMG_04_03 + #1.000 #0.000 &
          TOLSOMG_05_03 + #1.000 #0.105 TOLSOMG_02_04 + #1.000 #0.105 TOLSOMG_03_04 + &
          #1.000 #0.105 TOLSOMG_04_04 + #1.000 #0.000 TOLSOMG_05_04 + #1.000 &
          #0.105 TOLSOMG_03_05 + #1.000 #0.105 TOLSOMG_04_05 + #1.000 #0.105 &
          TOLSOMG_05_05 + #1.000 #0.105 TOLSOMG_03_06 + #1.000 #0.105 TOLSOMG_04_06 + &
          #1.000 #0.105 TOLSOMG_05_06
 S82.2  TOLSOMG_06_06 + OH = #1.000 #0.105 TOLSOMG_01_01 + #1.000 #0.105 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.000 TOLSOMG_05_01 + #1.000 #0.105 TOLSOMG_01_02 + #1.000 &
          #0.105 TOLSOMG_02_02 + #1.000 #0.105 TOLSOMG_03_02 + #1.000 #0.000 &
          TOLSOMG_04_02 + #1.000 #0.000 TOLSOMG_05_02 + #1.000 #0.105 TOLSOMG_02_03 + &
          #1.000 #0.105 TOLSOMG_03_03 + #1.000 #0.000 TOLSOMG_04_03 + #1.000 &
          #0.000 TOLSOMG_05_03 + #1.000 #0.105 TOLSOMG_02_04 + #1.000 #0.105 &
          TOLSOMG_03_04 + #1.000 #0.105 TOLSOMG_04_04 + #1.000 #0.000 TOLSOMG_05_04 + &
          #1.000 #0.105 TOLSOMG_03_05 + #1.000 #0.105 TOLSOMG_04_05 + #1.000 &
          #0.000 TOLSOMG_05_05 + #1.000 #0.105 TOLSOMG_03_06 + #1.000 #0.105 &
          TOLSOMG_04_06 + #1.000 #0.105 TOLSOMG_05_06 + #1.000 #0.105 TOLSOMG_04_07 + &
          #1.000 #0.105 TOLSOMG_05_07
 S83.2  TOLSOMG_06_07 + OH = #1.000 #0.125 TOLSOMG_01_01 + #1.000 #0.125 &
          TOLSOMG_02_01 + #1.000 #0.000 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.000 TOLSOMG_05_01 + #1.000 #0.125 TOLSOMG_01_02 + #1.000 &
          #0.125 TOLSOMG_02_02 + #1.000 #0.000 TOLSOMG_03_02 + #1.000 #0.000 &
          TOLSOMG_04_02 + #1.000 #0.000 TOLSOMG_05_02 + #1.000 #0.125 TOLSOMG_02_03 + &
          #1.000 #0.125 TOLSOMG_03_03 + #1.000 #0.000 TOLSOMG_04_03 + #1.000 &
          #0.000 TOLSOMG_05_03 + #1.000 #0.125 TOLSOMG_02_04 + #1.000 #0.125 &
          TOLSOMG_03_04 + #1.000 #0.000 TOLSOMG_04_04 + #1.000 #0.000 TOLSOMG_05_04 + &
          #1.000 #0.125 TOLSOMG_03_05 + #1.000 #0.125 TOLSOMG_04_05 + #1.000 &
          #0.000 TOLSOMG_05_05 + #1.000 #0.125 TOLSOMG_03_06 + #1.000 #0.125 &
          TOLSOMG_04_06 + #1.000 #0.000 TOLSOMG_05_06 + #1.000 #0.125 TOLSOMG_04_07 + &
          #1.000 #0.125 TOLSOMG_05_07
 S84.2  TOLSOMG_07_01 + OH = #0.002 #0.694 TOLSOMG_07_02 + #0.002 #0.031 &
          TOLSOMG_07_03 + #0.002 #0.002 TOLSOMG_07_04 + #0.002 #0.273 TOLSOMG_07_05 + &
          #0.998 #0.167 TOLSOMG_01_01 + #0.998 #0.167 TOLSOMG_02_01 + #0.998 &
          #0.167 TOLSOMG_03_01 + #0.998 #0.167 TOLSOMG_04_01 + #0.998 #0.167 &
          TOLSOMG_05_01 + #0.998 #0.167 TOLSOMG_06_01 + #0.998 #0.167 TOLSOMG_01_02 + &
          #0.998 #0.167 TOLSOMG_02_02 + #0.998 #0.167 TOLSOMG_03_02 + #0.998 &
          #0.167 TOLSOMG_04_02 + #0.998 #0.167 TOLSOMG_05_02 + #0.998 #0.167 &
          TOLSOMG_06_02
 S85.2  TOLSOMG_07_02 + OH = #0.001 #0.694 TOLSOMG_07_03 + #0.001 #0.031 &
          TOLSOMG_07_04 + #0.001 #0.002 TOLSOMG_07_05 + #0.001 #0.273 TOLSOMG_07_06 + &
          #0.999 #0.125 TOLSOMG_01_01 + #0.999 #0.125 TOLSOMG_02_01 + #0.999 &
          #0.125 TOLSOMG_03_01 + #0.999 #0.125 TOLSOMG_04_01 + #0.999 #0.125 &
          TOLSOMG_05_01 + #0.999 #0.000 TOLSOMG_06_01 + #0.999 #0.125 TOLSOMG_01_02 + &
          #0.999 #0.125 TOLSOMG_02_02 + #0.999 #0.125 TOLSOMG_03_02 + #0.999 &
          #0.125 TOLSOMG_04_02 + #0.999 #0.125 TOLSOMG_05_02 + #0.999 #0.125 &
          TOLSOMG_06_02 + #0.999 #0.125 TOLSOMG_02_03 + #0.999 #0.125 TOLSOMG_03_03 + &
          #0.999 #0.125 TOLSOMG_04_03 + #0.999 #0.125 TOLSOMG_05_03 + #0.999 &
          #0.125 TOLSOMG_06_03
 S86.2  TOLSOMG_07_03 + OH = #0.001 #0.694 TOLSOMG_07_04 + #0.001 #0.031 &
          TOLSOMG_07_05 + #0.001 #0.002 TOLSOMG_07_06 + #0.001 #0.273 TOLSOMG_07_07 + &
          #0.999 #0.100 TOLSOMG_01_01 + #0.999 #0.100 TOLSOMG_02_01 + #0.999 &
          #0.100 TOLSOMG_03_01 + #0.999 #0.100 TOLSOMG_04_01 + #0.999 #0.100 &
          TOLSOMG_05_01 + #0.999 #0.000 TOLSOMG_06_01 + #0.999 #0.100 TOLSOMG_01_02 + &
          #0.999 #0.100 TOLSOMG_02_02 + #0.999 #0.100 TOLSOMG_03_02 + #0.999 &
          #0.100 TOLSOMG_04_02 + #0.999 #0.100 TOLSOMG_05_02 + #0.999 #0.000 &
          TOLSOMG_06_02 + #0.999 #0.100 TOLSOMG_02_03 + #0.999 #0.100 TOLSOMG_03_03 + &
          #0.999 #0.100 TOLSOMG_04_03 + #0.999 #0.100 TOLSOMG_05_03 + #0.999 &
          #0.100 TOLSOMG_06_03 + #0.999 #0.100 TOLSOMG_02_04 + #0.999 #0.100 &
          TOLSOMG_03_04 + #0.999 #0.100 TOLSOMG_04_04 + #0.999 #0.100 TOLSOMG_05_04 + &
          #0.999 #0.100 TOLSOMG_06_04
 S87.2  TOLSOMG_07_04 + OH = #0.001 #0.694 TOLSOMG_07_05 + #0.001 #0.031 &
          TOLSOMG_07_06 + #0.001 #0.275 TOLSOMG_07_07 + #0.999 #0.091 TOLSOMG_01_01 + &
          #0.999 #0.091 TOLSOMG_02_01 + #0.999 #0.091 TOLSOMG_03_01 + #0.999 &
          #0.091 TOLSOMG_04_01 + #0.999 #0.000 TOLSOMG_05_01 + #0.999 #0.000 &
          TOLSOMG_06_01 + #0.999 #0.091 TOLSOMG_01_02 + #0.999 #0.091 TOLSOMG_02_02 + &
          #0.999 #0.091 TOLSOMG_03_02 + #0.999 #0.091 TOLSOMG_04_02 + #0.999 &
          #0.091 TOLSOMG_05_02 + #0.999 #0.000 TOLSOMG_06_02 + #0.999 #0.091 &
          TOLSOMG_02_03 + #0.999 #0.091 TOLSOMG_03_03 + #0.999 #0.091 TOLSOMG_04_03 + &
          #0.999 #0.091 TOLSOMG_05_03 + #0.999 #0.000 TOLSOMG_06_03 + #0.999 &
          #0.091 TOLSOMG_02_04 + #0.999 #0.091 TOLSOMG_03_04 + #0.999 #0.091 &
          TOLSOMG_04_04 + #0.999 #0.091 TOLSOMG_05_04 + #0.999 #0.091 TOLSOMG_06_04 + &
          #0.999 #0.091 TOLSOMG_03_05 + #0.999 #0.091 TOLSOMG_04_05 + #0.999 &
          #0.091 TOLSOMG_05_05 + #0.999 #0.091 TOLSOMG_06_05
 S88.2  TOLSOMG_07_05 + OH = #0.000 #0.694 TOLSOMG_07_06 + #0.000 #0.306 &
          TOLSOMG_07_07 + #1.000 #0.083 TOLSOMG_01_01 + #1.000 #0.083 TOLSOMG_02_01 + &
          #1.000 #0.083 TOLSOMG_03_01 + #1.000 #0.083 TOLSOMG_04_01 + #1.000 &
          #0.000 TOLSOMG_05_01 + #1.000 #0.000 TOLSOMG_06_01 + #1.000 #0.083 &
          TOLSOMG_01_02 + #1.000 #0.083 TOLSOMG_02_02 + #1.000 #0.083 TOLSOMG_03_02 + &
          #1.000 #0.083 TOLSOMG_04_02 + #1.000 #0.000 TOLSOMG_05_02 + #1.000 &
          #0.000 TOLSOMG_06_02 + #1.000 #0.083 TOLSOMG_02_03 + #1.000 #0.083 &
          TOLSOMG_03_03 + #1.000 #0.083 TOLSOMG_04_03 + #1.000 #0.083 TOLSOMG_05_03 + &
          #1.000 #0.000 TOLSOMG_06_03 + #1.000 #0.083 TOLSOMG_02_04 + #1.000 &
          #0.083 TOLSOMG_03_04 + #1.000 #0.083 TOLSOMG_04_04 + #1.000 #0.083 &
          TOLSOMG_05_04 + #1.000 #0.000 TOLSOMG_06_04 + #1.000 #0.083 TOLSOMG_03_05 + &
          #1.000 #0.083 TOLSOMG_04_05 + #1.000 #0.083 TOLSOMG_05_05 + #1.000 &
          #0.083 TOLSOMG_06_05 + #1.000 #0.083 TOLSOMG_03_06 + #1.000 #0.083 &
          TOLSOMG_04_06 + #1.000 #0.083 TOLSOMG_05_06 + #1.000 #0.083 TOLSOMG_06_06
 S89.2  TOLSOMG_07_06 + OH = #0.000 #1.000 TOLSOMG_07_07 + #1.000 #0.080 &
          TOLSOMG_01_01 + #1.000 #0.080 TOLSOMG_02_01 + #1.000 #0.080 TOLSOMG_03_01 + &
          #1.000 #0.000 TOLSOMG_04_01 + #1.000 #0.000 TOLSOMG_05_01 + #1.000 &
          #0.000 TOLSOMG_06_01 + #1.000 #0.080 TOLSOMG_01_02 + #1.000 #0.080 &
          TOLSOMG_02_02 + #1.000 #0.080 TOLSOMG_03_02 + #1.000 #0.080 TOLSOMG_04_02 + &
          #1.000 #0.000 TOLSOMG_05_02 + #1.000 #0.000 TOLSOMG_06_02 + #1.000 &
          #0.080 TOLSOMG_02_03 + #1.000 #0.080 TOLSOMG_03_03 + #1.000 #0.080 &
          TOLSOMG_04_03 + #1.000 #0.000 TOLSOMG_05_03 + #1.000 #0.000 TOLSOMG_06_03 + &
          #1.000 #0.080 TOLSOMG_02_04 + #1.000 #0.080 TOLSOMG_03_04 + #1.000 &
          #0.080 TOLSOMG_04_04 + #1.000 #0.080 TOLSOMG_05_04 + #1.000 #0.000 &
          TOLSOMG_06_04 + #1.000 #0.080 TOLSOMG_03_05 + #1.000 #0.080 TOLSOMG_04_05 + &
          #1.000 #0.080 TOLSOMG_05_05 + #1.000 #0.000 TOLSOMG_06_05 + #1.000 &
          #0.080 TOLSOMG_03_06 + #1.000 #0.080 TOLSOMG_04_06 + #1.000 #0.080 &
          TOLSOMG_05_06 + #1.000 #0.080 TOLSOMG_06_06 + #1.000 #0.080 TOLSOMG_04_07 + &
          #1.000 #0.080 TOLSOMG_05_07 + #1.000 #0.080 TOLSOMG_06_07
 S90.2  TOLSOMG_07_07 + OH = #1.000 #0.091 TOLSOMG_01_01 + #1.000 #0.091 &
          TOLSOMG_02_01 + #1.000 #0.091 TOLSOMG_03_01 + #1.000 #0.000 TOLSOMG_04_01 + &
          #1.000 #0.000 TOLSOMG_05_01 + #1.000 #0.000 TOLSOMG_06_01 + #1.000 &
          #0.091 TOLSOMG_01_02 + #1.000 #0.091 TOLSOMG_02_02 + #1.000 #0.091 &
          TOLSOMG_03_02 + #1.000 #0.000 TOLSOMG_04_02 + #1.000 #0.000 TOLSOMG_05_02 + &
          #1.000 #0.000 TOLSOMG_06_02 + #1.000 #0.091 TOLSOMG_02_03 + #1.000 &
          #0.091 TOLSOMG_03_03 + #1.000 #0.091 TOLSOMG_04_03 + #1.000 #0.000 &
          TOLSOMG_05_03 + #1.000 #0.000 TOLSOMG_06_03 + #1.000 #0.091 TOLSOMG_02_04 + &
          #1.000 #0.091 TOLSOMG_03_04 + #1.000 #0.091 TOLSOMG_04_04 + #1.000 &
          #0.000 TOLSOMG_05_04 + #1.000 #0.000 TOLSOMG_06_04 + #1.000 #0.091 &
          TOLSOMG_03_05 + #1.000 #0.091 TOLSOMG_04_05 + #1.000 #0.091 TOLSOMG_05_05 + &
          #1.000 #0.000 TOLSOMG_06_05 + #1.000 #0.091 TOLSOMG_03_06 + #1.000 &
          #0.091 TOLSOMG_04_06 + #1.000 #0.091 TOLSOMG_05_06 + #1.000 #0.000 &
          TOLSOMG_06_06 + #1.000 #0.091 TOLSOMG_04_07 + #1.000 #0.091 TOLSOMG_05_07 + &
          #1.000 #0.091 TOLSOMG_06_07
 BL23   XYLENE + OH = #1.000 #0.440 XYLSOMG_08_01 + #1.000 #0.128 XYLSOMG_08_02 + &
          #1.000 #0.016 XYLSOMG_08_03 + #1.000 #0.416 XYLSOMG_08_04 + XYLENE
 S91.3  XYLSOMG_02_01 + OH = #0.041 #0.440 XYLSOMG_02_02 + #0.041 #0.128 &
          XYLSOMG_02_03 + #0.041 #0.432 XYLSOMG_02_04 + #0.959 #1.000 XYLSOMG_01_01 + &
          #0.959 #1.000 XYLSOMG_01_02
 S92.3  XYLSOMG_02_02 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #2.000 &
          XYLSOMG_01_02
 S93.3  XYLSOMG_02_03 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_01_02
 S94.3  XYLSOMG_02_04 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_01_02
 S95.3  XYLSOMG_03_01 + OH = #0.065 #0.440 XYLSOMG_03_02 + #0.065 #0.128 &
          XYLSOMG_03_03 + #0.065 #0.016 XYLSOMG_03_04 + #0.065 #0.416 XYLSOMG_03_05 + &
          #0.935 #0.500 XYLSOMG_01_01 + #0.935 #0.500 XYLSOMG_02_01 + #0.935 &
          #0.500 XYLSOMG_01_02 + #0.935 #0.500 XYLSOMG_02_02
 S96.3  XYLSOMG_03_02 + OH = #0.024 #0.440 XYLSOMG_03_03 + #0.024 #0.128 &
          XYLSOMG_03_04 + #0.024 #0.016 XYLSOMG_03_05 + #0.024 #0.416 XYLSOMG_03_06 + &
          #0.976 #0.500 XYLSOMG_01_01 + #0.976 #0.000 XYLSOMG_02_01 + #0.976 &
          #0.500 XYLSOMG_01_02 + #0.976 #0.500 XYLSOMG_02_02 + #0.976 #0.500 &
          XYLSOMG_02_03
 S97.3  XYLSOMG_03_03 + OH = #1.000 #0.500 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.500 XYLSOMG_01_02 + #1.000 #0.000 XYLSOMG_02_02 + &
          #1.000 #0.500 XYLSOMG_02_03 + #1.000 #0.500 XYLSOMG_02_04
 S98.3  XYLSOMG_03_04 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #1.000 XYLSOMG_01_02 + #1.000 #0.000 XYLSOMG_02_02 + &
          #1.000 #0.000 XYLSOMG_02_03 + #1.000 #1.000 XYLSOMG_02_04
 S99.3  XYLSOMG_03_05 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_01_02 + #1.000 #0.000 XYLSOMG_02_02 + &
          #1.000 #0.000 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_02_04
 S100.3 XYLSOMG_03_06 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_01_02 + #1.000 #0.000 XYLSOMG_02_02 + &
          #1.000 #0.000 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_02_04
 S101.3 XYLSOMG_04_01 + OH = #0.081 #0.440 XYLSOMG_04_02 + #0.081 #0.128 &
          XYLSOMG_04_03 + #0.081 #0.016 XYLSOMG_04_04 + #0.081 #0.416 XYLSOMG_04_05 + &
          #0.919 #0.333 XYLSOMG_01_01 + #0.919 #0.333 XYLSOMG_02_01 + #0.919 &
          #0.333 XYLSOMG_03_01 + #0.919 #0.333 XYLSOMG_01_02 + #0.919 #0.333 &
          XYLSOMG_02_02 + #0.919 #0.333 XYLSOMG_03_02
 S102.3 XYLSOMG_04_02 + OH = #0.041 #0.440 XYLSOMG_04_03 + #0.041 #0.128 &
          XYLSOMG_04_04 + #0.041 #0.016 XYLSOMG_04_05 + #0.041 #0.416 XYLSOMG_04_06 + &
          #0.959 #0.286 XYLSOMG_01_01 + #0.959 #0.286 XYLSOMG_02_01 + #0.959 &
          #0.000 XYLSOMG_03_01 + #0.959 #0.286 XYLSOMG_01_02 + #0.959 #0.286 &
          XYLSOMG_02_02 + #0.959 #0.286 XYLSOMG_03_02 + #0.959 #0.286 XYLSOMG_02_03 + &
          #0.959 #0.286 XYLSOMG_03_03
 S103.3 XYLSOMG_04_03 + OH = #0.017 #0.440 XYLSOMG_04_04 + #0.017 #0.128 &
          XYLSOMG_04_05 + #0.017 #0.016 XYLSOMG_04_06 + #0.017 #0.416 XYLSOMG_04_07 + &
          #0.983 #0.250 XYLSOMG_01_01 + #0.983 #0.250 XYLSOMG_02_01 + #0.983 &
          #0.000 XYLSOMG_03_01 + #0.983 #0.250 XYLSOMG_01_02 + #0.983 #0.250 &
          XYLSOMG_02_02 + #0.983 #0.000 XYLSOMG_03_02 + #0.983 #0.250 XYLSOMG_02_03 + &
          #0.983 #0.250 XYLSOMG_03_03 + #0.983 #0.250 XYLSOMG_02_04 + #0.983 &
          #0.250 XYLSOMG_03_04
 S104.3 XYLSOMG_04_04 + OH = #1.000 #0.286 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.286 XYLSOMG_01_02 + &
          #1.000 #0.286 XYLSOMG_02_02 + #1.000 #0.000 XYLSOMG_03_02 + #1.000 &
          #0.286 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.286 &
          XYLSOMG_02_04 + #1.000 #0.286 XYLSOMG_03_04 + #1.000 #0.286 XYLSOMG_03_05
 S105.3 XYLSOMG_04_05 + OH = #1.000 #0.333 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.333 XYLSOMG_01_02 + &
          #1.000 #0.000 XYLSOMG_02_02 + #1.000 #0.000 XYLSOMG_03_02 + #1.000 &
          #0.333 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.333 &
          XYLSOMG_02_04 + #1.000 #0.000 XYLSOMG_03_04 + #1.000 #0.333 XYLSOMG_03_05 + &
          #1.000 #0.333 XYLSOMG_03_06
 S106.3 XYLSOMG_04_06 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.667 XYLSOMG_01_02 + &
          #1.000 #0.000 XYLSOMG_02_02 + #1.000 #0.000 XYLSOMG_03_02 + #1.000 &
          #0.000 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.667 &
          XYLSOMG_02_04 + #1.000 #0.000 XYLSOMG_03_04 + #1.000 #0.000 XYLSOMG_03_05 + &
          #1.000 #0.667 XYLSOMG_03_06
 S107.3 XYLSOMG_04_07 + OH = #1.000 #0.000 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_01_02 + &
          #1.000 #0.000 XYLSOMG_02_02 + #1.000 #0.000 XYLSOMG_03_02 + #1.000 &
          #0.000 XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.000 &
          XYLSOMG_02_04 + #1.000 #0.000 XYLSOMG_03_04 + #1.000 #0.000 XYLSOMG_03_05 + &
          #1.000 #0.000 XYLSOMG_03_06
 S108.3 XYLSOMG_05_01 + OH = #0.094 #0.440 XYLSOMG_05_02 + #0.094 #0.128 &
          XYLSOMG_05_03 + #0.094 #0.016 XYLSOMG_05_04 + #0.094 #0.416 XYLSOMG_05_05 + &
          #0.906 #0.250 XYLSOMG_01_01 + #0.906 #0.250 XYLSOMG_02_01 + #0.906 &
          #0.250 XYLSOMG_03_01 + #0.906 #0.250 XYLSOMG_04_01 + #0.906 #0.250 &
          XYLSOMG_01_02 + #0.906 #0.250 XYLSOMG_02_02 + #0.906 #0.250 XYLSOMG_03_02 + &
          #0.906 #0.250 XYLSOMG_04_02
 S109.3 XYLSOMG_05_02 + OH = #0.054 #0.440 XYLSOMG_05_03 + #0.054 #0.128 &
          XYLSOMG_05_04 + #0.054 #0.016 XYLSOMG_05_05 + #0.054 #0.416 XYLSOMG_05_06 + &
          #0.946 #0.182 XYLSOMG_01_01 + #0.946 #0.182 XYLSOMG_02_01 + #0.946 &
          #0.182 XYLSOMG_03_01 + #0.946 #0.000 XYLSOMG_04_01 + #0.946 #0.182 &
          XYLSOMG_01_02 + #0.946 #0.182 XYLSOMG_02_02 + #0.946 #0.182 XYLSOMG_03_02 + &
          #0.946 #0.182 XYLSOMG_04_02 + #0.946 #0.182 XYLSOMG_02_03 + #0.946 &
          #0.182 XYLSOMG_03_03 + #0.946 #0.182 XYLSOMG_04_03
 S110.3 XYLSOMG_05_03 + OH = #0.031 #0.440 XYLSOMG_05_04 + #0.031 #0.128 &
          XYLSOMG_05_05 + #0.031 #0.016 XYLSOMG_05_06 + #0.031 #0.416 XYLSOMG_05_07 + &
          #0.969 #0.154 XYLSOMG_01_01 + #0.969 #0.154 XYLSOMG_02_01 + #0.969 &
          #0.154 XYLSOMG_03_01 + #0.969 #0.000 XYLSOMG_04_01 + #0.969 #0.154 &
          XYLSOMG_01_02 + #0.969 #0.154 XYLSOMG_02_02 + #0.969 #0.154 XYLSOMG_03_02 + &
          #0.969 #0.000 XYLSOMG_04_02 + #0.969 #0.154 XYLSOMG_02_03 + #0.969 &
          #0.154 XYLSOMG_03_03 + #0.969 #0.154 XYLSOMG_04_03 + #0.969 #0.154 &
          XYLSOMG_02_04 + #0.969 #0.154 XYLSOMG_03_04 + #0.969 #0.154 XYLSOMG_04_04
 S111.3 XYLSOMG_05_04 + OH = #0.014 #0.440 XYLSOMG_05_05 + #0.014 #0.128 &
          XYLSOMG_05_06 + #0.014 #0.432 XYLSOMG_05_07 + #0.986 #0.143 XYLSOMG_01_01 + &
          #0.986 #0.143 XYLSOMG_02_01 + #0.986 #0.000 XYLSOMG_03_01 + #0.986 &
          #0.000 XYLSOMG_04_01 + #0.986 #0.143 XYLSOMG_01_02 + #0.986 #0.143 &
          XYLSOMG_02_02 + #0.986 #0.143 XYLSOMG_03_02 + #0.986 #0.000 XYLSOMG_04_02 + &
          #0.986 #0.143 XYLSOMG_02_03 + #0.986 #0.143 XYLSOMG_03_03 + #0.986 &
          #0.000 XYLSOMG_04_03 + #0.986 #0.143 XYLSOMG_02_04 + #0.986 #0.143 &
          XYLSOMG_03_04 + #0.986 #0.143 XYLSOMG_04_04 + #0.986 #0.143 XYLSOMG_03_05 + &
          #0.986 #0.143 XYLSOMG_04_05
 S112.3 XYLSOMG_05_05 + OH = #1.000 #0.143 XYLSOMG_01_01 + #1.000 #0.143 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.143 XYLSOMG_01_02 + #1.000 #0.143 XYLSOMG_02_02 + #1.000 &
          #0.000 XYLSOMG_03_02 + #1.000 #0.000 XYLSOMG_04_02 + #1.000 #0.143 &
          XYLSOMG_02_03 + #1.000 #0.143 XYLSOMG_03_03 + #1.000 #0.000 XYLSOMG_04_03 + &
          #1.000 #0.143 XYLSOMG_02_04 + #1.000 #0.143 XYLSOMG_03_04 + #1.000 &
          #0.000 XYLSOMG_04_04 + #1.000 #0.143 XYLSOMG_03_05 + #1.000 #0.143 &
          XYLSOMG_04_05 + #1.000 #0.143 XYLSOMG_03_06 + #1.000 #0.143 XYLSOMG_04_06
 S113.3 XYLSOMG_05_06 + OH = #1.000 #0.154 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.154 XYLSOMG_01_02 + #1.000 #0.154 XYLSOMG_02_02 + #1.000 &
          #0.000 XYLSOMG_03_02 + #1.000 #0.000 XYLSOMG_04_02 + #1.000 #0.154 &
          XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.000 XYLSOMG_04_03 + &
          #1.000 #0.154 XYLSOMG_02_04 + #1.000 #0.154 XYLSOMG_03_04 + #1.000 &
          #0.000 XYLSOMG_04_04 + #1.000 #0.154 XYLSOMG_03_05 + #1.000 #0.000 &
          XYLSOMG_04_05 + #1.000 #0.154 XYLSOMG_03_06 + #1.000 #0.154 XYLSOMG_04_06 + &
          #1.000 #0.154 XYLSOMG_04_07
 S114.3 XYLSOMG_05_07 + OH = #1.000 #0.200 XYLSOMG_01_01 + #1.000 #0.000 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.200 XYLSOMG_01_02 + #1.000 #0.000 XYLSOMG_02_02 + #1.000 &
          #0.000 XYLSOMG_03_02 + #1.000 #0.000 XYLSOMG_04_02 + #1.000 #0.200 &
          XYLSOMG_02_03 + #1.000 #0.000 XYLSOMG_03_03 + #1.000 #0.000 XYLSOMG_04_03 + &
          #1.000 #0.200 XYLSOMG_02_04 + #1.000 #0.000 XYLSOMG_03_04 + #1.000 &
          #0.000 XYLSOMG_04_04 + #1.000 #0.200 XYLSOMG_03_05 + #1.000 #0.000 &
          XYLSOMG_04_05 + #1.000 #0.200 XYLSOMG_03_06 + #1.000 #0.000 XYLSOMG_04_06 + &
          #1.000 #0.200 XYLSOMG_04_07
 S115.3 XYLSOMG_06_01 + OH = #0.104 #0.440 XYLSOMG_06_02 + #0.104 #0.128 &
          XYLSOMG_06_03 + #0.104 #0.016 XYLSOMG_06_04 + #0.104 #0.416 XYLSOMG_06_05 + &
          #0.896 #0.200 XYLSOMG_01_01 + #0.896 #0.200 XYLSOMG_02_01 + #0.896 &
          #0.200 XYLSOMG_03_01 + #0.896 #0.200 XYLSOMG_04_01 + #0.896 #0.200 &
          XYLSOMG_05_01 + #0.896 #0.200 XYLSOMG_01_02 + #0.896 #0.200 XYLSOMG_02_02 + &
          #0.896 #0.200 XYLSOMG_03_02 + #0.896 #0.200 XYLSOMG_04_02 + #0.896 &
          #0.200 XYLSOMG_05_02
 S116.3 XYLSOMG_06_02 + OH = #0.065 #0.440 XYLSOMG_06_03 + #0.065 #0.128 &
          XYLSOMG_06_04 + #0.065 #0.016 XYLSOMG_06_05 + #0.065 #0.416 XYLSOMG_06_06 + &
          #0.935 #0.154 XYLSOMG_01_01 + #0.935 #0.154 XYLSOMG_02_01 + #0.935 &
          #0.154 XYLSOMG_03_01 + #0.935 #0.154 XYLSOMG_04_01 + #0.935 #0.000 &
          XYLSOMG_05_01 + #0.935 #0.154 XYLSOMG_01_02 + #0.935 #0.154 XYLSOMG_02_02 + &
          #0.935 #0.154 XYLSOMG_03_02 + #0.935 #0.154 XYLSOMG_04_02 + #0.935 &
          #0.154 XYLSOMG_05_02 + #0.935 #0.154 XYLSOMG_02_03 + #0.935 #0.154 &
          XYLSOMG_03_03 + #0.935 #0.154 XYLSOMG_04_03 + #0.935 #0.154 XYLSOMG_05_03
 S117.3 XYLSOMG_06_03 + OH = #0.041 #0.440 XYLSOMG_06_04 + #0.041 #0.128 &
          XYLSOMG_06_05 + #0.041 #0.016 XYLSOMG_06_06 + #0.041 #0.416 XYLSOMG_06_07 + &
          #0.959 #0.125 XYLSOMG_01_01 + #0.959 #0.125 XYLSOMG_02_01 + #0.959 &
          #0.125 XYLSOMG_03_01 + #0.959 #0.125 XYLSOMG_04_01 + #0.959 #0.000 &
          XYLSOMG_05_01 + #0.959 #0.125 XYLSOMG_01_02 + #0.959 #0.125 XYLSOMG_02_02 + &
          #0.959 #0.125 XYLSOMG_03_02 + #0.959 #0.125 XYLSOMG_04_02 + #0.959 &
          #0.000 XYLSOMG_05_02 + #0.959 #0.125 XYLSOMG_02_03 + #0.959 #0.125 &
          XYLSOMG_03_03 + #0.959 #0.125 XYLSOMG_04_03 + #0.959 #0.125 XYLSOMG_05_03 + &
          #0.959 #0.125 XYLSOMG_02_04 + #0.959 #0.125 XYLSOMG_03_04 + #0.959 &
          #0.125 XYLSOMG_04_04 + #0.959 #0.125 XYLSOMG_05_04
 S118.3 XYLSOMG_06_04 + OH = #0.024 #0.440 XYLSOMG_06_05 + #0.024 #0.128 &
          XYLSOMG_06_06 + #0.024 #0.432 XYLSOMG_06_07 + #0.976 #0.111 XYLSOMG_01_01 + &
          #0.976 #0.111 XYLSOMG_02_01 + #0.976 #0.111 XYLSOMG_03_01 + #0.976 &
          #0.000 XYLSOMG_04_01 + #0.976 #0.000 XYLSOMG_05_01 + #0.976 #0.111 &
          XYLSOMG_01_02 + #0.976 #0.111 XYLSOMG_02_02 + #0.976 #0.111 XYLSOMG_03_02 + &
          #0.976 #0.111 XYLSOMG_04_02 + #0.976 #0.000 XYLSOMG_05_02 + #0.976 &
          #0.111 XYLSOMG_02_03 + #0.976 #0.111 XYLSOMG_03_03 + #0.976 #0.111 &
          XYLSOMG_04_03 + #0.976 #0.000 XYLSOMG_05_03 + #0.976 #0.111 XYLSOMG_02_04 + &
          #0.976 #0.111 XYLSOMG_03_04 + #0.976 #0.111 XYLSOMG_04_04 + #0.976 &
          #0.111 XYLSOMG_05_04 + #0.976 #0.111 XYLSOMG_03_05 + #0.976 #0.111 &
          XYLSOMG_04_05 + #0.976 #0.111 XYLSOMG_05_05
 S119.3 XYLSOMG_06_05 + OH = #0.011 #0.440 XYLSOMG_06_06 + #0.011 #0.560 &
          XYLSOMG_06_07 + #0.989 #0.105 XYLSOMG_01_01 + #0.989 #0.105 XYLSOMG_02_01 + &
          #0.989 #0.105 XYLSOMG_03_01 + #0.989 #0.000 XYLSOMG_04_01 + #0.989 &
          #0.000 XYLSOMG_05_01 + #0.989 #0.105 XYLSOMG_01_02 + #0.989 #0.105 &
          XYLSOMG_02_02 + #0.989 #0.105 XYLSOMG_03_02 + #0.989 #0.000 XYLSOMG_04_02 + &
          #0.989 #0.000 XYLSOMG_05_02 + #0.989 #0.105 XYLSOMG_02_03 + #0.989 &
          #0.105 XYLSOMG_03_03 + #0.989 #0.105 XYLSOMG_04_03 + #0.989 #0.000 &
          XYLSOMG_05_03 + #0.989 #0.105 XYLSOMG_02_04 + #0.989 #0.105 XYLSOMG_03_04 + &
          #0.989 #0.105 XYLSOMG_04_04 + #0.989 #0.000 XYLSOMG_05_04 + #0.989 &
          #0.105 XYLSOMG_03_05 + #0.989 #0.105 XYLSOMG_04_05 + #0.989 #0.105 &
          XYLSOMG_05_05 + #0.989 #0.105 XYLSOMG_03_06 + #0.989 #0.105 XYLSOMG_04_06 + &
          #0.989 #0.105 XYLSOMG_05_06
 S120.3 XYLSOMG_06_06 + OH = #1.000 #0.105 XYLSOMG_01_01 + #1.000 #0.105 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.000 XYLSOMG_05_01 + #1.000 #0.105 XYLSOMG_01_02 + #1.000 &
          #0.105 XYLSOMG_02_02 + #1.000 #0.105 XYLSOMG_03_02 + #1.000 #0.000 &
          XYLSOMG_04_02 + #1.000 #0.000 XYLSOMG_05_02 + #1.000 #0.105 XYLSOMG_02_03 + &
          #1.000 #0.105 XYLSOMG_03_03 + #1.000 #0.000 XYLSOMG_04_03 + #1.000 &
          #0.000 XYLSOMG_05_03 + #1.000 #0.105 XYLSOMG_02_04 + #1.000 #0.105 &
          XYLSOMG_03_04 + #1.000 #0.105 XYLSOMG_04_04 + #1.000 #0.000 XYLSOMG_05_04 + &
          #1.000 #0.105 XYLSOMG_03_05 + #1.000 #0.105 XYLSOMG_04_05 + #1.000 &
          #0.000 XYLSOMG_05_05 + #1.000 #0.105 XYLSOMG_03_06 + #1.000 #0.105 &
          XYLSOMG_04_06 + #1.000 #0.105 XYLSOMG_05_06 + #1.000 #0.105 XYLSOMG_04_07 + &
          #1.000 #0.105 XYLSOMG_05_07
 S121.3 XYLSOMG_06_07 + OH = #1.000 #0.125 XYLSOMG_01_01 + #1.000 #0.125 &
          XYLSOMG_02_01 + #1.000 #0.000 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.000 XYLSOMG_05_01 + #1.000 #0.125 XYLSOMG_01_02 + #1.000 &
          #0.125 XYLSOMG_02_02 + #1.000 #0.000 XYLSOMG_03_02 + #1.000 #0.000 &
          XYLSOMG_04_02 + #1.000 #0.000 XYLSOMG_05_02 + #1.000 #0.125 XYLSOMG_02_03 + &
          #1.000 #0.125 XYLSOMG_03_03 + #1.000 #0.000 XYLSOMG_04_03 + #1.000 &
          #0.000 XYLSOMG_05_03 + #1.000 #0.125 XYLSOMG_02_04 + #1.000 #0.125 &
          XYLSOMG_03_04 + #1.000 #0.000 XYLSOMG_04_04 + #1.000 #0.000 XYLSOMG_05_04 + &
          #1.000 #0.125 XYLSOMG_03_05 + #1.000 #0.125 XYLSOMG_04_05 + #1.000 &
          #0.000 XYLSOMG_05_05 + #1.000 #0.125 XYLSOMG_03_06 + #1.000 #0.125 &
          XYLSOMG_04_06 + #1.000 #0.000 XYLSOMG_05_06 + #1.000 #0.125 XYLSOMG_04_07 + &
          #1.000 #0.125 XYLSOMG_05_07
 S122.3 XYLSOMG_07_01 + OH = #0.112 #0.440 XYLSOMG_07_02 + #0.112 #0.128 &
          XYLSOMG_07_03 + #0.112 #0.016 XYLSOMG_07_04 + #0.112 #0.416 XYLSOMG_07_05 + &
          #0.888 #0.167 XYLSOMG_01_01 + #0.888 #0.167 XYLSOMG_02_01 + #0.888 &
          #0.167 XYLSOMG_03_01 + #0.888 #0.167 XYLSOMG_04_01 + #0.888 #0.167 &
          XYLSOMG_05_01 + #0.888 #0.167 XYLSOMG_06_01 + #0.888 #0.167 XYLSOMG_01_02 + &
          #0.888 #0.167 XYLSOMG_02_02 + #0.888 #0.167 XYLSOMG_03_02 + #0.888 &
          #0.167 XYLSOMG_04_02 + #0.888 #0.167 XYLSOMG_05_02 + #0.888 #0.167 &
          XYLSOMG_06_02
 S123.3 XYLSOMG_07_02 + OH = #0.074 #0.440 XYLSOMG_07_03 + #0.074 #0.128 &
          XYLSOMG_07_04 + #0.074 #0.016 XYLSOMG_07_05 + #0.074 #0.416 XYLSOMG_07_06 + &
          #0.926 #0.125 XYLSOMG_01_01 + #0.926 #0.125 XYLSOMG_02_01 + #0.926 &
          #0.125 XYLSOMG_03_01 + #0.926 #0.125 XYLSOMG_04_01 + #0.926 #0.125 &
          XYLSOMG_05_01 + #0.926 #0.000 XYLSOMG_06_01 + #0.926 #0.125 XYLSOMG_01_02 + &
          #0.926 #0.125 XYLSOMG_02_02 + #0.926 #0.125 XYLSOMG_03_02 + #0.926 &
          #0.125 XYLSOMG_04_02 + #0.926 #0.125 XYLSOMG_05_02 + #0.926 #0.125 &
          XYLSOMG_06_02 + #0.926 #0.125 XYLSOMG_02_03 + #0.926 #0.125 XYLSOMG_03_03 + &
          #0.926 #0.125 XYLSOMG_04_03 + #0.926 #0.125 XYLSOMG_05_03 + #0.926 &
          #0.125 XYLSOMG_06_03
 S124.3 XYLSOMG_07_03 + OH = #0.050 #0.440 XYLSOMG_07_04 + #0.050 #0.128 &
          XYLSOMG_07_05 + #0.050 #0.016 XYLSOMG_07_06 + #0.050 #0.416 XYLSOMG_07_07 + &
          #0.950 #0.100 XYLSOMG_01_01 + #0.950 #0.100 XYLSOMG_02_01 + #0.950 &
          #0.100 XYLSOMG_03_01 + #0.950 #0.100 XYLSOMG_04_01 + #0.950 #0.100 &
          XYLSOMG_05_01 + #0.950 #0.000 XYLSOMG_06_01 + #0.950 #0.100 XYLSOMG_01_02 + &
          #0.950 #0.100 XYLSOMG_02_02 + #0.950 #0.100 XYLSOMG_03_02 + #0.950 &
          #0.100 XYLSOMG_04_02 + #0.950 #0.100 XYLSOMG_05_02 + #0.950 #0.000 &
          XYLSOMG_06_02 + #0.950 #0.100 XYLSOMG_02_03 + #0.950 #0.100 XYLSOMG_03_03 + &
          #0.950 #0.100 XYLSOMG_04_03 + #0.950 #0.100 XYLSOMG_05_03 + #0.950 &
          #0.100 XYLSOMG_06_03 + #0.950 #0.100 XYLSOMG_02_04 + #0.950 #0.100 &
          XYLSOMG_03_04 + #0.950 #0.100 XYLSOMG_04_04 + #0.950 #0.100 XYLSOMG_05_04 + &
          #0.950 #0.100 XYLSOMG_06_04
 S125.3 XYLSOMG_07_04 + OH = #0.034 #0.440 XYLSOMG_07_05 + #0.034 #0.128 &
          XYLSOMG_07_06 + #0.034 #0.432 XYLSOMG_07_07 + #0.966 #0.091 XYLSOMG_01_01 + &
          #0.966 #0.091 XYLSOMG_02_01 + #0.966 #0.091 XYLSOMG_03_01 + #0.966 &
          #0.091 XYLSOMG_04_01 + #0.966 #0.000 XYLSOMG_05_01 + #0.966 #0.000 &
          XYLSOMG_06_01 + #0.966 #0.091 XYLSOMG_01_02 + #0.966 #0.091 XYLSOMG_02_02 + &
          #0.966 #0.091 XYLSOMG_03_02 + #0.966 #0.091 XYLSOMG_04_02 + #0.966 &
          #0.091 XYLSOMG_05_02 + #0.966 #0.000 XYLSOMG_06_02 + #0.966 #0.091 &
          XYLSOMG_02_03 + #0.966 #0.091 XYLSOMG_03_03 + #0.966 #0.091 XYLSOMG_04_03 + &
          #0.966 #0.091 XYLSOMG_05_03 + #0.966 #0.000 XYLSOMG_06_03 + #0.966 &
          #0.091 XYLSOMG_02_04 + #0.966 #0.091 XYLSOMG_03_04 + #0.966 #0.091 &
          XYLSOMG_04_04 + #0.966 #0.091 XYLSOMG_05_04 + #0.966 #0.091 XYLSOMG_06_04 + &
          #0.966 #0.091 XYLSOMG_03_05 + #0.966 #0.091 XYLSOMG_04_05 + #0.966 &
          #0.091 XYLSOMG_05_05 + #0.966 #0.091 XYLSOMG_06_05
 S126.3 XYLSOMG_07_05 + OH = #0.020 #0.440 XYLSOMG_07_06 + #0.020 #0.560 &
          XYLSOMG_07_07 + #0.980 #0.083 XYLSOMG_01_01 + #0.980 #0.083 XYLSOMG_02_01 + &
          #0.980 #0.083 XYLSOMG_03_01 + #0.980 #0.083 XYLSOMG_04_01 + #0.980 &
          #0.000 XYLSOMG_05_01 + #0.980 #0.000 XYLSOMG_06_01 + #0.980 #0.083 &
          XYLSOMG_01_02 + #0.980 #0.083 XYLSOMG_02_02 + #0.980 #0.083 XYLSOMG_03_02 + &
          #0.980 #0.083 XYLSOMG_04_02 + #0.980 #0.000 XYLSOMG_05_02 + #0.980 &
          #0.000 XYLSOMG_06_02 + #0.980 #0.083 XYLSOMG_02_03 + #0.980 #0.083 &
          XYLSOMG_03_03 + #0.980 #0.083 XYLSOMG_04_03 + #0.980 #0.083 XYLSOMG_05_03 + &
          #0.980 #0.000 XYLSOMG_06_03 + #0.980 #0.083 XYLSOMG_02_04 + #0.980 &
          #0.083 XYLSOMG_03_04 + #0.980 #0.083 XYLSOMG_04_04 + #0.980 #0.083 &
          XYLSOMG_05_04 + #0.980 #0.000 XYLSOMG_06_04 + #0.980 #0.083 XYLSOMG_03_05 + &
          #0.980 #0.083 XYLSOMG_04_05 + #0.980 #0.083 XYLSOMG_05_05 + #0.980 &
          #0.083 XYLSOMG_06_05 + #0.980 #0.083 XYLSOMG_03_06 + #0.980 #0.083 &
          XYLSOMG_04_06 + #0.980 #0.083 XYLSOMG_05_06 + #0.980 #0.083 XYLSOMG_06_06
 S127.3 XYLSOMG_07_06 + OH = #0.009 #1.000 XYLSOMG_07_07 + #0.991 #0.080 &
          XYLSOMG_01_01 + #0.991 #0.080 XYLSOMG_02_01 + #0.991 #0.080 XYLSOMG_03_01 + &
          #0.991 #0.000 XYLSOMG_04_01 + #0.991 #0.000 XYLSOMG_05_01 + #0.991 &
          #0.000 XYLSOMG_06_01 + #0.991 #0.080 XYLSOMG_01_02 + #0.991 #0.080 &
          XYLSOMG_02_02 + #0.991 #0.080 XYLSOMG_03_02 + #0.991 #0.080 XYLSOMG_04_02 + &
          #0.991 #0.000 XYLSOMG_05_02 + #0.991 #0.000 XYLSOMG_06_02 + #0.991 &
          #0.080 XYLSOMG_02_03 + #0.991 #0.080 XYLSOMG_03_03 + #0.991 #0.080 &
          XYLSOMG_04_03 + #0.991 #0.000 XYLSOMG_05_03 + #0.991 #0.000 XYLSOMG_06_03 + &
          #0.991 #0.080 XYLSOMG_02_04 + #0.991 #0.080 XYLSOMG_03_04 + #0.991 &
          #0.080 XYLSOMG_04_04 + #0.991 #0.080 XYLSOMG_05_04 + #0.991 #0.000 &
          XYLSOMG_06_04 + #0.991 #0.080 XYLSOMG_03_05 + #0.991 #0.080 XYLSOMG_04_05 + &
          #0.991 #0.080 XYLSOMG_05_05 + #0.991 #0.000 XYLSOMG_06_05 + #0.991 &
          #0.080 XYLSOMG_03_06 + #0.991 #0.080 XYLSOMG_04_06 + #0.991 #0.080 &
          XYLSOMG_05_06 + #0.991 #0.080 XYLSOMG_06_06 + #0.991 #0.080 XYLSOMG_04_07 + &
          #0.991 #0.080 XYLSOMG_05_07 + #0.991 #0.080 XYLSOMG_06_07
 S128.3 XYLSOMG_07_07 + OH = #1.000 #0.091 XYLSOMG_01_01 + #1.000 #0.091 &
          XYLSOMG_02_01 + #1.000 #0.091 XYLSOMG_03_01 + #1.000 #0.000 XYLSOMG_04_01 + &
          #1.000 #0.000 XYLSOMG_05_01 + #1.000 #0.000 XYLSOMG_06_01 + #1.000 &
          #0.091 XYLSOMG_01_02 + #1.000 #0.091 XYLSOMG_02_02 + #1.000 #0.091 &
          XYLSOMG_03_02 + #1.000 #0.000 XYLSOMG_04_02 + #1.000 #0.000 XYLSOMG_05_02 + &
          #1.000 #0.000 XYLSOMG_06_02 + #1.000 #0.091 XYLSOMG_02_03 + #1.000 &
          #0.091 XYLSOMG_03_03 + #1.000 #0.091 XYLSOMG_04_03 + #1.000 #0.000 &
          XYLSOMG_05_03 + #1.000 #0.000 XYLSOMG_06_03 + #1.000 #0.091 XYLSOMG_02_04 + &
          #1.000 #0.091 XYLSOMG_03_04 + #1.000 #0.091 XYLSOMG_04_04 + #1.000 &
          #0.000 XYLSOMG_05_04 + #1.000 #0.000 XYLSOMG_06_04 + #1.000 #0.091 &
          XYLSOMG_03_05 + #1.000 #0.091 XYLSOMG_04_05 + #1.000 #0.091 XYLSOMG_05_05 + &
          #1.000 #0.000 XYLSOMG_06_05 + #1.000 #0.091 XYLSOMG_03_06 + #1.000 &
          #0.091 XYLSOMG_04_06 + #1.000 #0.091 XYLSOMG_05_06 + #1.000 #0.000 &
          XYLSOMG_06_06 + #1.000 #0.091 XYLSOMG_04_07 + #1.000 #0.091 XYLSOMG_05_07 + &
          #1.000 #0.091 XYLSOMG_06_07
 S129.3 XYLSOMG_08_01 + OH = #0.119 #0.440 XYLSOMG_08_02 + #0.119 #0.128 &
          XYLSOMG_08_03 + #0.119 #0.016 XYLSOMG_08_04 + #0.119 #0.416 XYLSOMG_08_05 + &
          #0.881 #0.143 XYLSOMG_01_01 + #0.881 #0.143 XYLSOMG_02_01 + #0.881 &
          #0.143 XYLSOMG_03_01 + #0.881 #0.143 XYLSOMG_04_01 + #0.881 #0.143 &
          XYLSOMG_05_01 + #0.881 #0.143 XYLSOMG_06_01 + #0.881 #0.143 XYLSOMG_07_01 + &
          #0.881 #0.143 XYLSOMG_01_02 + #0.881 #0.143 XYLSOMG_02_02 + #0.881 &
          #0.143 XYLSOMG_03_02 + #0.881 #0.143 XYLSOMG_04_02 + #0.881 #0.143 &
          XYLSOMG_05_02 + #0.881 #0.143 XYLSOMG_06_02 + #0.881 #0.143 XYLSOMG_07_02
 S130.3 XYLSOMG_08_02 + OH = #0.081 #0.440 XYLSOMG_08_03 + #0.081 #0.128 &
          XYLSOMG_08_04 + #0.081 #0.016 XYLSOMG_08_05 + #0.081 #0.416 XYLSOMG_08_06 + &
          #0.919 #0.105 XYLSOMG_01_01 + #0.919 #0.105 XYLSOMG_02_01 + #0.919 &
          #0.105 XYLSOMG_03_01 + #0.919 #0.105 XYLSOMG_04_01 + #0.919 #0.105 &
          XYLSOMG_05_01 + #0.919 #0.105 XYLSOMG_06_01 + #0.919 #0.000 XYLSOMG_07_01 + &
          #0.919 #0.105 XYLSOMG_01_02 + #0.919 #0.105 XYLSOMG_02_02 + #0.919 &
          #0.105 XYLSOMG_03_02 + #0.919 #0.105 XYLSOMG_04_02 + #0.919 #0.105 &
          XYLSOMG_05_02 + #0.919 #0.105 XYLSOMG_06_02 + #0.919 #0.105 XYLSOMG_07_02 + &
          #0.919 #0.105 XYLSOMG_02_03 + #0.919 #0.105 XYLSOMG_03_03 + #0.919 &
          #0.105 XYLSOMG_04_03 + #0.919 #0.105 XYLSOMG_05_03 + #0.919 #0.105 &
          XYLSOMG_06_03 + #0.919 #0.105 XYLSOMG_07_03
 S131.3 XYLSOMG_08_03 + OH = #0.058 #0.440 XYLSOMG_08_04 + #0.058 #0.128 &
          XYLSOMG_08_05 + #0.058 #0.016 XYLSOMG_08_06 + #0.058 #0.416 XYLSOMG_08_07 + &
          #0.942 #0.083 XYLSOMG_01_01 + #0.942 #0.083 XYLSOMG_02_01 + #0.942 &
          #0.083 XYLSOMG_03_01 + #0.942 #0.083 XYLSOMG_04_01 + #0.942 #0.083 &
          XYLSOMG_05_01 + #0.942 #0.083 XYLSOMG_06_01 + #0.942 #0.000 XYLSOMG_07_01 + &
          #0.942 #0.083 XYLSOMG_01_02 + #0.942 #0.083 XYLSOMG_02_02 + #0.942 &
          #0.083 XYLSOMG_03_02 + #0.942 #0.083 XYLSOMG_04_02 + #0.942 #0.083 &
          XYLSOMG_05_02 + #0.942 #0.083 XYLSOMG_06_02 + #0.942 #0.000 XYLSOMG_07_02 + &
          #0.942 #0.083 XYLSOMG_02_03 + #0.942 #0.083 XYLSOMG_03_03 + #0.942 &
          #0.083 XYLSOMG_04_03 + #0.942 #0.083 XYLSOMG_05_03 + #0.942 #0.083 &
          XYLSOMG_06_03 + #0.942 #0.083 XYLSOMG_07_03 + #0.942 #0.083 XYLSOMG_02_04 + &
          #0.942 #0.083 XYLSOMG_03_04 + #0.942 #0.083 XYLSOMG_04_04 + #0.942 &
          #0.083 XYLSOMG_05_04 + #0.942 #0.083 XYLSOMG_06_04 + #0.942 #0.083 &
          XYLSOMG_07_04
 S132.3 XYLSOMG_08_04 + OH = #0.041 #0.440 XYLSOMG_08_05 + #0.041 #0.128 &
          XYLSOMG_08_06 + #0.041 #0.432 XYLSOMG_08_07 + #0.959 #0.074 XYLSOMG_01_01 + &
          #0.959 #0.074 XYLSOMG_02_01 + #0.959 #0.074 XYLSOMG_03_01 + #0.959 &
          #0.074 XYLSOMG_04_01 + #0.959 #0.074 XYLSOMG_05_01 + #0.959 #0.000 &
          XYLSOMG_06_01 + #0.959 #0.000 XYLSOMG_07_01 + #0.959 #0.074 XYLSOMG_01_02 + &
          #0.959 #0.074 XYLSOMG_02_02 + #0.959 #0.074 XYLSOMG_03_02 + #0.959 &
          #0.074 XYLSOMG_04_02 + #0.959 #0.074 XYLSOMG_05_02 + #0.959 #0.074 &
          XYLSOMG_06_02 + #0.959 #0.000 XYLSOMG_07_02 + #0.959 #0.074 XYLSOMG_02_03 + &
          #0.959 #0.074 XYLSOMG_03_03 + #0.959 #0.074 XYLSOMG_04_03 + #0.959 &
          #0.074 XYLSOMG_05_03 + #0.959 #0.074 XYLSOMG_06_03 + #0.959 #0.000 &
          XYLSOMG_07_03 + #0.959 #0.074 XYLSOMG_02_04 + #0.959 #0.074 XYLSOMG_03_04 + &
          #0.959 #0.074 XYLSOMG_04_04 + #0.959 #0.074 XYLSOMG_05_04 + #0.959 &
          #0.074 XYLSOMG_06_04 + #0.959 #0.074 XYLSOMG_07_04 + #0.959 #0.074 &
          XYLSOMG_03_05 + #0.959 #0.074 XYLSOMG_04_05 + #0.959 #0.074 XYLSOMG_05_05 + &
          #0.959 #0.074 XYLSOMG_06_05 + #0.959 #0.074 XYLSOMG_07_05
 S133.3 XYLSOMG_08_05 + OH = #0.028 #0.440 XYLSOMG_08_06 + #0.028 #0.560 &
          XYLSOMG_08_07 + #0.972 #0.067 XYLSOMG_01_01 + #0.972 #0.067 XYLSOMG_02_01 + &
          #0.972 #0.067 XYLSOMG_03_01 + #0.972 #0.067 XYLSOMG_04_01 + #0.972 &
          #0.067 XYLSOMG_05_01 + #0.972 #0.000 XYLSOMG_06_01 + #0.972 #0.000 &
          XYLSOMG_07_01 + #0.972 #0.067 XYLSOMG_01_02 + #0.972 #0.067 XYLSOMG_02_02 + &
          #0.972 #0.067 XYLSOMG_03_02 + #0.972 #0.067 XYLSOMG_04_02 + #0.972 &
          #0.067 XYLSOMG_05_02 + #0.972 #0.000 XYLSOMG_06_02 + #0.972 #0.000 &
          XYLSOMG_07_02 + #0.972 #0.067 XYLSOMG_02_03 + #0.972 #0.067 XYLSOMG_03_03 + &
          #0.972 #0.067 XYLSOMG_04_03 + #0.972 #0.067 XYLSOMG_05_03 + #0.972 &
          #0.067 XYLSOMG_06_03 + #0.972 #0.000 XYLSOMG_07_03 + #0.972 #0.067 &
          XYLSOMG_02_04 + #0.972 #0.067 XYLSOMG_03_04 + #0.972 #0.067 XYLSOMG_04_04 + &
          #0.972 #0.067 XYLSOMG_05_04 + #0.972 #0.067 XYLSOMG_06_04 + #0.972 &
          #0.000 XYLSOMG_07_04 + #0.972 #0.067 XYLSOMG_03_05 + #0.972 #0.067 &
          XYLSOMG_04_05 + #0.972 #0.067 XYLSOMG_05_05 + #0.972 #0.067 XYLSOMG_06_05 + &
          #0.972 #0.067 XYLSOMG_07_05 + #0.972 #0.067 XYLSOMG_03_06 + #0.972 &
          #0.067 XYLSOMG_04_06 + #0.972 #0.067 XYLSOMG_05_06 + #0.972 #0.067 &
          XYLSOMG_06_06 + #0.972 #0.067 XYLSOMG_07_06
 S134.3 XYLSOMG_08_06 + OH = #0.017 #1.000 XYLSOMG_08_07 + #0.983 #0.065 &
          XYLSOMG_01_01 + #0.983 #0.065 XYLSOMG_02_01 + #0.983 #0.065 XYLSOMG_03_01 + &
          #0.983 #0.065 XYLSOMG_04_01 + #0.983 #0.000 XYLSOMG_05_01 + #0.983 &
          #0.000 XYLSOMG_06_01 + #0.983 #0.000 XYLSOMG_07_01 + #0.983 #0.065 &
          XYLSOMG_01_02 + #0.983 #0.065 XYLSOMG_02_02 + #0.983 #0.065 XYLSOMG_03_02 + &
          #0.983 #0.065 XYLSOMG_04_02 + #0.983 #0.065 XYLSOMG_05_02 + #0.983 &
          #0.000 XYLSOMG_06_02 + #0.983 #0.000 XYLSOMG_07_02 + #0.983 #0.065 &
          XYLSOMG_02_03 + #0.983 #0.065 XYLSOMG_03_03 + #0.983 #0.065 XYLSOMG_04_03 + &
          #0.983 #0.065 XYLSOMG_05_03 + #0.983 #0.000 XYLSOMG_06_03 + #0.983 &
          #0.000 XYLSOMG_07_03 + #0.983 #0.065 XYLSOMG_02_04 + #0.983 #0.065 &
          XYLSOMG_03_04 + #0.983 #0.065 XYLSOMG_04_04 + #0.983 #0.065 XYLSOMG_05_04 + &
          #0.983 #0.065 XYLSOMG_06_04 + #0.983 #0.000 XYLSOMG_07_04 + #0.983 &
          #0.065 XYLSOMG_03_05 + #0.983 #0.065 XYLSOMG_04_05 + #0.983 #0.065 &
          XYLSOMG_05_05 + #0.983 #0.065 XYLSOMG_06_05 + #0.983 #0.000 XYLSOMG_07_05 + &
          #0.983 #0.065 XYLSOMG_03_06 + #0.983 #0.065 XYLSOMG_04_06 + #0.983 &
          #0.065 XYLSOMG_05_06 + #0.983 #0.065 XYLSOMG_06_06 + #0.983 #0.065 &
          XYLSOMG_07_06 + #0.983 #0.065 XYLSOMG_04_07 + #0.983 #0.065 XYLSOMG_05_07 + &
          #0.983 #0.065 XYLSOMG_06_07 + #0.983 #0.065 XYLSOMG_07_07
 S135.3 XYLSOMG_08_07 + OH = #1.000 #0.071 XYLSOMG_01_01 + #1.000 #0.071 &
          XYLSOMG_02_01 + #1.000 #0.071 XYLSOMG_03_01 + #1.000 #0.071 XYLSOMG_04_01 + &
          #1.000 #0.000 XYLSOMG_05_01 + #1.000 #0.000 XYLSOMG_06_01 + #1.000 &
          #0.000 XYLSOMG_07_01 + #1.000 #0.071 XYLSOMG_01_02 + #1.000 #0.071 &
          XYLSOMG_02_02 + #1.000 #0.071 XYLSOMG_03_02 + #1.000 #0.071 XYLSOMG_04_02 + &
          #1.000 #0.000 XYLSOMG_05_02 + #1.000 #0.000 XYLSOMG_06_02 + #1.000 &
          #0.000 XYLSOMG_07_02 + #1.000 #0.071 XYLSOMG_02_03 + #1.000 #0.071 &
          XYLSOMG_03_03 + #1.000 #0.071 XYLSOMG_04_03 + #1.000 #0.071 XYLSOMG_05_03 + &
          #1.000 #0.000 XYLSOMG_06_03 + #1.000 #0.000 XYLSOMG_07_03 + #1.000 &
          #0.071 XYLSOMG_02_04 + #1.000 #0.071 XYLSOMG_03_04 + #1.000 #0.071 &
          XYLSOMG_04_04 + #1.000 #0.071 XYLSOMG_05_04 + #1.000 #0.000 XYLSOMG_06_04 + &
          #1.000 #0.000 XYLSOMG_07_04 + #1.000 #0.071 XYLSOMG_03_05 + #1.000 &
          #0.071 XYLSOMG_04_05 + #1.000 #0.071 XYLSOMG_05_05 + #1.000 #0.071 &
          XYLSOMG_06_05 + #1.000 #0.000 XYLSOMG_07_05 + #1.000 #0.071 XYLSOMG_03_06 + &
          #1.000 #0.071 XYLSOMG_04_06 + #1.000 #0.071 XYLSOMG_05_06 + #1.000 &
          #0.071 XYLSOMG_06_06 + #1.000 #0.000 XYLSOMG_07_06 + #1.000 #0.071 &
          XYLSOMG_04_07 + #1.000 #0.071 XYLSOMG_05_07 + #1.000 #0.071 XYLSOMG_06_07 + &
          #1.000 #0.071 XYLSOMG_07_07
 BL26   ISOPRENE + OH = #1.000 #0.576 ISPSOMG_05_01 + #1.000 #0.053 ISPSOMG_05_02 + &
          #1.000 #0.356 ISPSOMG_05_03 + #1.000 #0.015 ISPSOMG_05_04 + ISOPRENE
 S136.4 ISPSOMG_02_01 + OH = #0.010 #0.576 ISPSOMG_02_02 + #0.010 #0.053 &
          ISPSOMG_02_03 + #0.010 #0.371 ISPSOMG_02_04 + #0.990 #1.000 ISPSOMG_01_01 + &
          #0.990 #1.000 ISPSOMG_01_02
 S137.4 ISPSOMG_02_02 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #2.000 &
          ISPSOMG_01_02
 S138.4 ISPSOMG_02_03 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_01_02
 S139.4 ISPSOMG_02_04 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_01_02
 S140.4 ISPSOMG_03_01 + OH = #0.016 #0.576 ISPSOMG_03_02 + #0.016 #0.053 &
          ISPSOMG_03_03 + #0.016 #0.356 ISPSOMG_03_04 + #0.016 #0.015 ISPSOMG_03_05 + &
          #0.984 #0.500 ISPSOMG_01_01 + #0.984 #0.500 ISPSOMG_02_01 + #0.984 &
          #0.500 ISPSOMG_01_02 + #0.984 #0.500 ISPSOMG_02_02
 S141.4 ISPSOMG_03_02 + OH = #0.006 #0.576 ISPSOMG_03_03 + #0.006 #0.053 &
          ISPSOMG_03_04 + #0.006 #0.356 ISPSOMG_03_05 + #0.006 #0.015 ISPSOMG_03_06 + &
          #0.994 #0.500 ISPSOMG_01_01 + #0.994 #0.000 ISPSOMG_02_01 + #0.994 &
          #0.500 ISPSOMG_01_02 + #0.994 #0.500 ISPSOMG_02_02 + #0.994 #0.500 &
          ISPSOMG_02_03
 S142.4 ISPSOMG_03_03 + OH = #1.000 #0.500 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.500 ISPSOMG_01_02 + #1.000 #0.000 ISPSOMG_02_02 + &
          #1.000 #0.500 ISPSOMG_02_03 + #1.000 #0.500 ISPSOMG_02_04
 S143.4 ISPSOMG_03_04 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #1.000 ISPSOMG_01_02 + #1.000 #0.000 ISPSOMG_02_02 + &
          #1.000 #0.000 ISPSOMG_02_03 + #1.000 #1.000 ISPSOMG_02_04
 S144.4 ISPSOMG_03_05 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_01_02 + #1.000 #0.000 ISPSOMG_02_02 + &
          #1.000 #0.000 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_02_04
 S145.4 ISPSOMG_03_06 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_01_02 + #1.000 #0.000 ISPSOMG_02_02 + &
          #1.000 #0.000 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_02_04
 S146.4 ISPSOMG_04_01 + OH = #0.021 #0.576 ISPSOMG_04_02 + #0.021 #0.053 &
          ISPSOMG_04_03 + #0.021 #0.356 ISPSOMG_04_04 + #0.021 #0.015 ISPSOMG_04_05 + &
          #0.979 #0.333 ISPSOMG_01_01 + #0.979 #0.333 ISPSOMG_02_01 + #0.979 &
          #0.333 ISPSOMG_03_01 + #0.979 #0.333 ISPSOMG_01_02 + #0.979 #0.333 &
          ISPSOMG_02_02 + #0.979 #0.333 ISPSOMG_03_02
 S147.4 ISPSOMG_04_02 + OH = #0.010 #0.576 ISPSOMG_04_03 + #0.010 #0.053 &
          ISPSOMG_04_04 + #0.010 #0.356 ISPSOMG_04_05 + #0.010 #0.015 ISPSOMG_04_06 + &
          #0.990 #0.286 ISPSOMG_01_01 + #0.990 #0.286 ISPSOMG_02_01 + #0.990 &
          #0.000 ISPSOMG_03_01 + #0.990 #0.286 ISPSOMG_01_02 + #0.990 #0.286 &
          ISPSOMG_02_02 + #0.990 #0.286 ISPSOMG_03_02 + #0.990 #0.286 ISPSOMG_02_03 + &
          #0.990 #0.286 ISPSOMG_03_03
 S148.4 ISPSOMG_04_03 + OH = #0.004 #0.576 ISPSOMG_04_04 + #0.004 #0.053 &
          ISPSOMG_04_05 + #0.004 #0.356 ISPSOMG_04_06 + #0.004 #0.015 ISPSOMG_04_07 + &
          #0.996 #0.250 ISPSOMG_01_01 + #0.996 #0.250 ISPSOMG_02_01 + #0.996 &
          #0.000 ISPSOMG_03_01 + #0.996 #0.250 ISPSOMG_01_02 + #0.996 #0.250 &
          ISPSOMG_02_02 + #0.996 #0.000 ISPSOMG_03_02 + #0.996 #0.250 ISPSOMG_02_03 + &
          #0.996 #0.250 ISPSOMG_03_03 + #0.996 #0.250 ISPSOMG_02_04 + #0.996 &
          #0.250 ISPSOMG_03_04
 S149.4 ISPSOMG_04_04 + OH = #1.000 #0.286 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.286 ISPSOMG_01_02 + &
          #1.000 #0.286 ISPSOMG_02_02 + #1.000 #0.000 ISPSOMG_03_02 + #1.000 &
          #0.286 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.286 &
          ISPSOMG_02_04 + #1.000 #0.286 ISPSOMG_03_04 + #1.000 #0.286 ISPSOMG_03_05
 S150.4 ISPSOMG_04_05 + OH = #1.000 #0.333 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.333 ISPSOMG_01_02 + &
          #1.000 #0.000 ISPSOMG_02_02 + #1.000 #0.000 ISPSOMG_03_02 + #1.000 &
          #0.333 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.333 &
          ISPSOMG_02_04 + #1.000 #0.000 ISPSOMG_03_04 + #1.000 #0.333 ISPSOMG_03_05 + &
          #1.000 #0.333 ISPSOMG_03_06
 S151.4 ISPSOMG_04_06 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.667 ISPSOMG_01_02 + &
          #1.000 #0.000 ISPSOMG_02_02 + #1.000 #0.000 ISPSOMG_03_02 + #1.000 &
          #0.000 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.667 &
          ISPSOMG_02_04 + #1.000 #0.000 ISPSOMG_03_04 + #1.000 #0.000 ISPSOMG_03_05 + &
          #1.000 #0.667 ISPSOMG_03_06
 S152.4 ISPSOMG_04_07 + OH = #1.000 #0.000 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.000 ISPSOMG_01_02 + &
          #1.000 #0.000 ISPSOMG_02_02 + #1.000 #0.000 ISPSOMG_03_02 + #1.000 &
          #0.000 ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.000 &
          ISPSOMG_02_04 + #1.000 #0.000 ISPSOMG_03_04 + #1.000 #0.000 ISPSOMG_03_05 + &
          #1.000 #0.000 ISPSOMG_03_06
 S153.4 ISPSOMG_05_01 + OH = #0.024 #0.576 ISPSOMG_05_02 + #0.024 #0.053 &
          ISPSOMG_05_03 + #0.024 #0.356 ISPSOMG_05_04 + #0.024 #0.015 ISPSOMG_05_05 + &
          #0.976 #0.250 ISPSOMG_01_01 + #0.976 #0.250 ISPSOMG_02_01 + #0.976 &
          #0.250 ISPSOMG_03_01 + #0.976 #0.250 ISPSOMG_04_01 + #0.976 #0.250 &
          ISPSOMG_01_02 + #0.976 #0.250 ISPSOMG_02_02 + #0.976 #0.250 ISPSOMG_03_02 + &
          #0.976 #0.250 ISPSOMG_04_02
 S154.4 ISPSOMG_05_02 + OH = #0.014 #0.576 ISPSOMG_05_03 + #0.014 #0.053 &
          ISPSOMG_05_04 + #0.014 #0.356 ISPSOMG_05_05 + #0.014 #0.015 ISPSOMG_05_06 + &
          #0.986 #0.182 ISPSOMG_01_01 + #0.986 #0.182 ISPSOMG_02_01 + #0.986 &
          #0.182 ISPSOMG_03_01 + #0.986 #0.000 ISPSOMG_04_01 + #0.986 #0.182 &
          ISPSOMG_01_02 + #0.986 #0.182 ISPSOMG_02_02 + #0.986 #0.182 ISPSOMG_03_02 + &
          #0.986 #0.182 ISPSOMG_04_02 + #0.986 #0.182 ISPSOMG_02_03 + #0.986 &
          #0.182 ISPSOMG_03_03 + #0.986 #0.182 ISPSOMG_04_03
 S155.4 ISPSOMG_05_03 + OH = #0.008 #0.576 ISPSOMG_05_04 + #0.008 #0.053 &
          ISPSOMG_05_05 + #0.008 #0.356 ISPSOMG_05_06 + #0.008 #0.015 ISPSOMG_05_07 + &
          #0.992 #0.154 ISPSOMG_01_01 + #0.992 #0.154 ISPSOMG_02_01 + #0.992 &
          #0.154 ISPSOMG_03_01 + #0.992 #0.000 ISPSOMG_04_01 + #0.992 #0.154 &
          ISPSOMG_01_02 + #0.992 #0.154 ISPSOMG_02_02 + #0.992 #0.154 ISPSOMG_03_02 + &
          #0.992 #0.000 ISPSOMG_04_02 + #0.992 #0.154 ISPSOMG_02_03 + #0.992 &
          #0.154 ISPSOMG_03_03 + #0.992 #0.154 ISPSOMG_04_03 + #0.992 #0.154 &
          ISPSOMG_02_04 + #0.992 #0.154 ISPSOMG_03_04 + #0.992 #0.154 ISPSOMG_04_04
 S156.4 ISPSOMG_05_04 + OH = #0.003 #0.576 ISPSOMG_05_05 + #0.003 #0.053 &
          ISPSOMG_05_06 + #0.003 #0.371 ISPSOMG_05_07 + #0.997 #0.143 ISPSOMG_01_01 + &
          #0.997 #0.143 ISPSOMG_02_01 + #0.997 #0.000 ISPSOMG_03_01 + #0.997 &
          #0.000 ISPSOMG_04_01 + #0.997 #0.143 ISPSOMG_01_02 + #0.997 #0.143 &
          ISPSOMG_02_02 + #0.997 #0.143 ISPSOMG_03_02 + #0.997 #0.000 ISPSOMG_04_02 + &
          #0.997 #0.143 ISPSOMG_02_03 + #0.997 #0.143 ISPSOMG_03_03 + #0.997 &
          #0.000 ISPSOMG_04_03 + #0.997 #0.143 ISPSOMG_02_04 + #0.997 #0.143 &
          ISPSOMG_03_04 + #0.997 #0.143 ISPSOMG_04_04 + #0.997 #0.143 ISPSOMG_03_05 + &
          #0.997 #0.143 ISPSOMG_04_05
 S157.4 ISPSOMG_05_05 + OH = #1.000 #0.143 ISPSOMG_01_01 + #1.000 #0.143 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.000 ISPSOMG_04_01 + &
          #1.000 #0.143 ISPSOMG_01_02 + #1.000 #0.143 ISPSOMG_02_02 + #1.000 &
          #0.000 ISPSOMG_03_02 + #1.000 #0.000 ISPSOMG_04_02 + #1.000 #0.143 &
          ISPSOMG_02_03 + #1.000 #0.143 ISPSOMG_03_03 + #1.000 #0.000 ISPSOMG_04_03 + &
          #1.000 #0.143 ISPSOMG_02_04 + #1.000 #0.143 ISPSOMG_03_04 + #1.000 &
          #0.000 ISPSOMG_04_04 + #1.000 #0.143 ISPSOMG_03_05 + #1.000 #0.143 &
          ISPSOMG_04_05 + #1.000 #0.143 ISPSOMG_03_06 + #1.000 #0.143 ISPSOMG_04_06
 S158.4 ISPSOMG_05_06 + OH = #1.000 #0.154 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.000 ISPSOMG_04_01 + &
          #1.000 #0.154 ISPSOMG_01_02 + #1.000 #0.154 ISPSOMG_02_02 + #1.000 &
          #0.000 ISPSOMG_03_02 + #1.000 #0.000 ISPSOMG_04_02 + #1.000 #0.154 &
          ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.000 ISPSOMG_04_03 + &
          #1.000 #0.154 ISPSOMG_02_04 + #1.000 #0.154 ISPSOMG_03_04 + #1.000 &
          #0.000 ISPSOMG_04_04 + #1.000 #0.154 ISPSOMG_03_05 + #1.000 #0.000 &
          ISPSOMG_04_05 + #1.000 #0.154 ISPSOMG_03_06 + #1.000 #0.154 ISPSOMG_04_06 + &
          #1.000 #0.154 ISPSOMG_04_07
 S159.4 ISPSOMG_05_07 + OH = #1.000 #0.200 ISPSOMG_01_01 + #1.000 #0.000 &
          ISPSOMG_02_01 + #1.000 #0.000 ISPSOMG_03_01 + #1.000 #0.000 ISPSOMG_04_01 + &
          #1.000 #0.200 ISPSOMG_01_02 + #1.000 #0.000 ISPSOMG_02_02 + #1.000 &
          #0.000 ISPSOMG_03_02 + #1.000 #0.000 ISPSOMG_04_02 + #1.000 #0.200 &
          ISPSOMG_02_03 + #1.000 #0.000 ISPSOMG_03_03 + #1.000 #0.000 ISPSOMG_04_03 + &
          #1.000 #0.200 ISPSOMG_02_04 + #1.000 #0.000 ISPSOMG_03_04 + #1.000 &
          #0.000 ISPSOMG_04_04 + #1.000 #0.200 ISPSOMG_03_05 + #1.000 #0.000 &
          ISPSOMG_04_05 + #1.000 #0.200 ISPSOMG_03_06 + #1.000 #0.000 ISPSOMG_04_06 + &
          #1.000 #0.200 ISPSOMG_04_07
 BL24   TERPENE + OH = #1.000 #0.361 TRPSOMG_10_01 + #1.000 #0.007 TRPSOMG_10_02 + &
          #1.000 #0.606 TRPSOMG_10_03 + #1.000 #0.026 TRPSOMG_10_04 + TERPENE
 S160.5 TRPSOMG_02_01 + OH = #0.059 #0.361 TRPSOMG_02_02 + #0.059 #0.007 &
          TRPSOMG_02_03 + #0.059 #0.632 TRPSOMG_02_04 + #0.941 #1.000 TRPSOMG_01_01 + &
          #0.941 #1.000 TRPSOMG_01_02
 S161.5 TRPSOMG_02_02 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #2.000 &
          TRPSOMG_01_02
 S162.5 TRPSOMG_02_03 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_01_02
 S163.5 TRPSOMG_02_04 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_01_02
 S164.5 TRPSOMG_03_01 + OH = #0.091 #0.361 TRPSOMG_03_02 + #0.091 #0.007 &
          TRPSOMG_03_03 + #0.091 #0.606 TRPSOMG_03_04 + #0.091 #0.026 TRPSOMG_03_05 + &
          #0.909 #0.500 TRPSOMG_01_01 + #0.909 #0.500 TRPSOMG_02_01 + #0.909 &
          #0.500 TRPSOMG_01_02 + #0.909 #0.500 TRPSOMG_02_02
 S165.5 TRPSOMG_03_02 + OH = #0.035 #0.361 TRPSOMG_03_03 + #0.035 #0.007 &
          TRPSOMG_03_04 + #0.035 #0.606 TRPSOMG_03_05 + #0.035 #0.026 TRPSOMG_03_06 + &
          #0.965 #0.500 TRPSOMG_01_01 + #0.965 #0.000 TRPSOMG_02_01 + #0.965 &
          #0.500 TRPSOMG_01_02 + #0.965 #0.500 TRPSOMG_02_02 + #0.965 #0.500 &
          TRPSOMG_02_03
 S166.5 TRPSOMG_03_03 + OH = #1.000 #0.500 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.500 TRPSOMG_01_02 + #1.000 #0.000 TRPSOMG_02_02 + &
          #1.000 #0.500 TRPSOMG_02_03 + #1.000 #0.500 TRPSOMG_02_04
 S167.5 TRPSOMG_03_04 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #1.000 TRPSOMG_01_02 + #1.000 #0.000 TRPSOMG_02_02 + &
          #1.000 #0.000 TRPSOMG_02_03 + #1.000 #1.000 TRPSOMG_02_04
 S168.5 TRPSOMG_03_05 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_01_02 + #1.000 #0.000 TRPSOMG_02_02 + &
          #1.000 #0.000 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_02_04
 S169.5 TRPSOMG_03_06 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_01_02 + #1.000 #0.000 TRPSOMG_02_02 + &
          #1.000 #0.000 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_02_04
 S170.5 TRPSOMG_04_01 + OH = #0.114 #0.361 TRPSOMG_04_02 + #0.114 #0.007 &
          TRPSOMG_04_03 + #0.114 #0.606 TRPSOMG_04_04 + #0.114 #0.026 TRPSOMG_04_05 + &
          #0.886 #0.333 TRPSOMG_01_01 + #0.886 #0.333 TRPSOMG_02_01 + #0.886 &
          #0.333 TRPSOMG_03_01 + #0.886 #0.333 TRPSOMG_01_02 + #0.886 #0.333 &
          TRPSOMG_02_02 + #0.886 #0.333 TRPSOMG_03_02
 S171.5 TRPSOMG_04_02 + OH = #0.059 #0.361 TRPSOMG_04_03 + #0.059 #0.007 &
          TRPSOMG_04_04 + #0.059 #0.606 TRPSOMG_04_05 + #0.059 #0.026 TRPSOMG_04_06 + &
          #0.941 #0.286 TRPSOMG_01_01 + #0.941 #0.286 TRPSOMG_02_01 + #0.941 &
          #0.000 TRPSOMG_03_01 + #0.941 #0.286 TRPSOMG_01_02 + #0.941 #0.286 &
          TRPSOMG_02_02 + #0.941 #0.286 TRPSOMG_03_02 + #0.941 #0.286 TRPSOMG_02_03 + &
          #0.941 #0.286 TRPSOMG_03_03
 S172.5 TRPSOMG_04_03 + OH = #0.025 #0.361 TRPSOMG_04_04 + #0.025 #0.007 &
          TRPSOMG_04_05 + #0.025 #0.606 TRPSOMG_04_06 + #0.025 #0.026 TRPSOMG_04_07 + &
          #0.975 #0.250 TRPSOMG_01_01 + #0.975 #0.250 TRPSOMG_02_01 + #0.975 &
          #0.000 TRPSOMG_03_01 + #0.975 #0.250 TRPSOMG_01_02 + #0.975 #0.250 &
          TRPSOMG_02_02 + #0.975 #0.000 TRPSOMG_03_02 + #0.975 #0.250 TRPSOMG_02_03 + &
          #0.975 #0.250 TRPSOMG_03_03 + #0.975 #0.250 TRPSOMG_02_04 + #0.975 &
          #0.250 TRPSOMG_03_04
 S173.5 TRPSOMG_04_04 + OH = #1.000 #0.286 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.286 TRPSOMG_01_02 + &
          #1.000 #0.286 TRPSOMG_02_02 + #1.000 #0.000 TRPSOMG_03_02 + #1.000 &
          #0.286 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.286 &
          TRPSOMG_02_04 + #1.000 #0.286 TRPSOMG_03_04 + #1.000 #0.286 TRPSOMG_03_05
 S174.5 TRPSOMG_04_05 + OH = #1.000 #0.333 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.333 TRPSOMG_01_02 + &
          #1.000 #0.000 TRPSOMG_02_02 + #1.000 #0.000 TRPSOMG_03_02 + #1.000 &
          #0.333 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.333 &
          TRPSOMG_02_04 + #1.000 #0.000 TRPSOMG_03_04 + #1.000 #0.333 TRPSOMG_03_05 + &
          #1.000 #0.333 TRPSOMG_03_06
 S175.5 TRPSOMG_04_06 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.667 TRPSOMG_01_02 + &
          #1.000 #0.000 TRPSOMG_02_02 + #1.000 #0.000 TRPSOMG_03_02 + #1.000 &
          #0.000 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.667 &
          TRPSOMG_02_04 + #1.000 #0.000 TRPSOMG_03_04 + #1.000 #0.000 TRPSOMG_03_05 + &
          #1.000 #0.667 TRPSOMG_03_06
 S176.5 TRPSOMG_04_07 + OH = #1.000 #0.000 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_01_02 + &
          #1.000 #0.000 TRPSOMG_02_02 + #1.000 #0.000 TRPSOMG_03_02 + #1.000 &
          #0.000 TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.000 &
          TRPSOMG_02_04 + #1.000 #0.000 TRPSOMG_03_04 + #1.000 #0.000 TRPSOMG_03_05 + &
          #1.000 #0.000 TRPSOMG_03_06
 S177.5 TRPSOMG_05_01 + OH = #0.131 #0.361 TRPSOMG_05_02 + #0.131 #0.007 &
          TRPSOMG_05_03 + #0.131 #0.606 TRPSOMG_05_04 + #0.131 #0.026 TRPSOMG_05_05 + &
          #0.869 #0.250 TRPSOMG_01_01 + #0.869 #0.250 TRPSOMG_02_01 + #0.869 &
          #0.250 TRPSOMG_03_01 + #0.869 #0.250 TRPSOMG_04_01 + #0.869 #0.250 &
          TRPSOMG_01_02 + #0.869 #0.250 TRPSOMG_02_02 + #0.869 #0.250 TRPSOMG_03_02 + &
          #0.869 #0.250 TRPSOMG_04_02
 S178.5 TRPSOMG_05_02 + OH = #0.077 #0.361 TRPSOMG_05_03 + #0.077 #0.007 &
          TRPSOMG_05_04 + #0.077 #0.606 TRPSOMG_05_05 + #0.077 #0.026 TRPSOMG_05_06 + &
          #0.923 #0.182 TRPSOMG_01_01 + #0.923 #0.182 TRPSOMG_02_01 + #0.923 &
          #0.182 TRPSOMG_03_01 + #0.923 #0.000 TRPSOMG_04_01 + #0.923 #0.182 &
          TRPSOMG_01_02 + #0.923 #0.182 TRPSOMG_02_02 + #0.923 #0.182 TRPSOMG_03_02 + &
          #0.923 #0.182 TRPSOMG_04_02 + #0.923 #0.182 TRPSOMG_02_03 + #0.923 &
          #0.182 TRPSOMG_03_03 + #0.923 #0.182 TRPSOMG_04_03
 S179.5 TRPSOMG_05_03 + OH = #0.043 #0.361 TRPSOMG_05_04 + #0.043 #0.007 &
          TRPSOMG_05_05 + #0.043 #0.606 TRPSOMG_05_06 + #0.043 #0.026 TRPSOMG_05_07 + &
          #0.957 #0.154 TRPSOMG_01_01 + #0.957 #0.154 TRPSOMG_02_01 + #0.957 &
          #0.154 TRPSOMG_03_01 + #0.957 #0.000 TRPSOMG_04_01 + #0.957 #0.154 &
          TRPSOMG_01_02 + #0.957 #0.154 TRPSOMG_02_02 + #0.957 #0.154 TRPSOMG_03_02 + &
          #0.957 #0.000 TRPSOMG_04_02 + #0.957 #0.154 TRPSOMG_02_03 + #0.957 &
          #0.154 TRPSOMG_03_03 + #0.957 #0.154 TRPSOMG_04_03 + #0.957 #0.154 &
          TRPSOMG_02_04 + #0.957 #0.154 TRPSOMG_03_04 + #0.957 #0.154 TRPSOMG_04_04
 S180.5 TRPSOMG_05_04 + OH = #0.019 #0.361 TRPSOMG_05_05 + #0.019 #0.007 &
          TRPSOMG_05_06 + #0.019 #0.632 TRPSOMG_05_07 + #0.981 #0.143 TRPSOMG_01_01 + &
          #0.981 #0.143 TRPSOMG_02_01 + #0.981 #0.000 TRPSOMG_03_01 + #0.981 &
          #0.000 TRPSOMG_04_01 + #0.981 #0.143 TRPSOMG_01_02 + #0.981 #0.143 &
          TRPSOMG_02_02 + #0.981 #0.143 TRPSOMG_03_02 + #0.981 #0.000 TRPSOMG_04_02 + &
          #0.981 #0.143 TRPSOMG_02_03 + #0.981 #0.143 TRPSOMG_03_03 + #0.981 &
          #0.000 TRPSOMG_04_03 + #0.981 #0.143 TRPSOMG_02_04 + #0.981 #0.143 &
          TRPSOMG_03_04 + #0.981 #0.143 TRPSOMG_04_04 + #0.981 #0.143 TRPSOMG_03_05 + &
          #0.981 #0.143 TRPSOMG_04_05
 S181.5 TRPSOMG_05_05 + OH = #1.000 #0.143 TRPSOMG_01_01 + #1.000 #0.143 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.143 TRPSOMG_01_02 + #1.000 #0.143 TRPSOMG_02_02 + #1.000 &
          #0.000 TRPSOMG_03_02 + #1.000 #0.000 TRPSOMG_04_02 + #1.000 #0.143 &
          TRPSOMG_02_03 + #1.000 #0.143 TRPSOMG_03_03 + #1.000 #0.000 TRPSOMG_04_03 + &
          #1.000 #0.143 TRPSOMG_02_04 + #1.000 #0.143 TRPSOMG_03_04 + #1.000 &
          #0.000 TRPSOMG_04_04 + #1.000 #0.143 TRPSOMG_03_05 + #1.000 #0.143 &
          TRPSOMG_04_05 + #1.000 #0.143 TRPSOMG_03_06 + #1.000 #0.143 TRPSOMG_04_06
 S182.5 TRPSOMG_05_06 + OH = #1.000 #0.154 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.154 TRPSOMG_01_02 + #1.000 #0.154 TRPSOMG_02_02 + #1.000 &
          #0.000 TRPSOMG_03_02 + #1.000 #0.000 TRPSOMG_04_02 + #1.000 #0.154 &
          TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.000 TRPSOMG_04_03 + &
          #1.000 #0.154 TRPSOMG_02_04 + #1.000 #0.154 TRPSOMG_03_04 + #1.000 &
          #0.000 TRPSOMG_04_04 + #1.000 #0.154 TRPSOMG_03_05 + #1.000 #0.000 &
          TRPSOMG_04_05 + #1.000 #0.154 TRPSOMG_03_06 + #1.000 #0.154 TRPSOMG_04_06 + &
          #1.000 #0.154 TRPSOMG_04_07
 S183.5 TRPSOMG_05_07 + OH = #1.000 #0.200 TRPSOMG_01_01 + #1.000 #0.000 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.200 TRPSOMG_01_02 + #1.000 #0.000 TRPSOMG_02_02 + #1.000 &
          #0.000 TRPSOMG_03_02 + #1.000 #0.000 TRPSOMG_04_02 + #1.000 #0.200 &
          TRPSOMG_02_03 + #1.000 #0.000 TRPSOMG_03_03 + #1.000 #0.000 TRPSOMG_04_03 + &
          #1.000 #0.200 TRPSOMG_02_04 + #1.000 #0.000 TRPSOMG_03_04 + #1.000 &
          #0.000 TRPSOMG_04_04 + #1.000 #0.200 TRPSOMG_03_05 + #1.000 #0.000 &
          TRPSOMG_04_05 + #1.000 #0.200 TRPSOMG_03_06 + #1.000 #0.000 TRPSOMG_04_06 + &
          #1.000 #0.200 TRPSOMG_04_07
 S184.5 TRPSOMG_06_01 + OH = #0.144 #0.361 TRPSOMG_06_02 + #0.144 #0.007 &
          TRPSOMG_06_03 + #0.144 #0.606 TRPSOMG_06_04 + #0.144 #0.026 TRPSOMG_06_05 + &
          #0.856 #0.200 TRPSOMG_01_01 + #0.856 #0.200 TRPSOMG_02_01 + #0.856 &
          #0.200 TRPSOMG_03_01 + #0.856 #0.200 TRPSOMG_04_01 + #0.856 #0.200 &
          TRPSOMG_05_01 + #0.856 #0.200 TRPSOMG_01_02 + #0.856 #0.200 TRPSOMG_02_02 + &
          #0.856 #0.200 TRPSOMG_03_02 + #0.856 #0.200 TRPSOMG_04_02 + #0.856 &
          #0.200 TRPSOMG_05_02
 S185.5 TRPSOMG_06_02 + OH = #0.091 #0.361 TRPSOMG_06_03 + #0.091 #0.007 &
          TRPSOMG_06_04 + #0.091 #0.606 TRPSOMG_06_05 + #0.091 #0.026 TRPSOMG_06_06 + &
          #0.909 #0.154 TRPSOMG_01_01 + #0.909 #0.154 TRPSOMG_02_01 + #0.909 &
          #0.154 TRPSOMG_03_01 + #0.909 #0.154 TRPSOMG_04_01 + #0.909 #0.000 &
          TRPSOMG_05_01 + #0.909 #0.154 TRPSOMG_01_02 + #0.909 #0.154 TRPSOMG_02_02 + &
          #0.909 #0.154 TRPSOMG_03_02 + #0.909 #0.154 TRPSOMG_04_02 + #0.909 &
          #0.154 TRPSOMG_05_02 + #0.909 #0.154 TRPSOMG_02_03 + #0.909 #0.154 &
          TRPSOMG_03_03 + #0.909 #0.154 TRPSOMG_04_03 + #0.909 #0.154 TRPSOMG_05_03
 S186.5 TRPSOMG_06_03 + OH = #0.059 #0.361 TRPSOMG_06_04 + #0.059 #0.007 &
          TRPSOMG_06_05 + #0.059 #0.606 TRPSOMG_06_06 + #0.059 #0.026 TRPSOMG_06_07 + &
          #0.941 #0.125 TRPSOMG_01_01 + #0.941 #0.125 TRPSOMG_02_01 + #0.941 &
          #0.125 TRPSOMG_03_01 + #0.941 #0.125 TRPSOMG_04_01 + #0.941 #0.000 &
          TRPSOMG_05_01 + #0.941 #0.125 TRPSOMG_01_02 + #0.941 #0.125 TRPSOMG_02_02 + &
          #0.941 #0.125 TRPSOMG_03_02 + #0.941 #0.125 TRPSOMG_04_02 + #0.941 &
          #0.000 TRPSOMG_05_02 + #0.941 #0.125 TRPSOMG_02_03 + #0.941 #0.125 &
          TRPSOMG_03_03 + #0.941 #0.125 TRPSOMG_04_03 + #0.941 #0.125 TRPSOMG_05_03 + &
          #0.941 #0.125 TRPSOMG_02_04 + #0.941 #0.125 TRPSOMG_03_04 + #0.941 &
          #0.125 TRPSOMG_04_04 + #0.941 #0.125 TRPSOMG_05_04
 S187.5 TRPSOMG_06_04 + OH = #0.035 #0.361 TRPSOMG_06_05 + #0.035 #0.007 &
          TRPSOMG_06_06 + #0.035 #0.632 TRPSOMG_06_07 + #0.965 #0.111 TRPSOMG_01_01 + &
          #0.965 #0.111 TRPSOMG_02_01 + #0.965 #0.111 TRPSOMG_03_01 + #0.965 &
          #0.000 TRPSOMG_04_01 + #0.965 #0.000 TRPSOMG_05_01 + #0.965 #0.111 &
          TRPSOMG_01_02 + #0.965 #0.111 TRPSOMG_02_02 + #0.965 #0.111 TRPSOMG_03_02 + &
          #0.965 #0.111 TRPSOMG_04_02 + #0.965 #0.000 TRPSOMG_05_02 + #0.965 &
          #0.111 TRPSOMG_02_03 + #0.965 #0.111 TRPSOMG_03_03 + #0.965 #0.111 &
          TRPSOMG_04_03 + #0.965 #0.000 TRPSOMG_05_03 + #0.965 #0.111 TRPSOMG_02_04 + &
          #0.965 #0.111 TRPSOMG_03_04 + #0.965 #0.111 TRPSOMG_04_04 + #0.965 &
          #0.111 TRPSOMG_05_04 + #0.965 #0.111 TRPSOMG_03_05 + #0.965 #0.111 &
          TRPSOMG_04_05 + #0.965 #0.111 TRPSOMG_05_05
 S188.5 TRPSOMG_06_05 + OH = #0.016 #0.361 TRPSOMG_06_06 + #0.016 #0.639 &
          TRPSOMG_06_07 + #0.984 #0.105 TRPSOMG_01_01 + #0.984 #0.105 TRPSOMG_02_01 + &
          #0.984 #0.105 TRPSOMG_03_01 + #0.984 #0.000 TRPSOMG_04_01 + #0.984 &
          #0.000 TRPSOMG_05_01 + #0.984 #0.105 TRPSOMG_01_02 + #0.984 #0.105 &
          TRPSOMG_02_02 + #0.984 #0.105 TRPSOMG_03_02 + #0.984 #0.000 TRPSOMG_04_02 + &
          #0.984 #0.000 TRPSOMG_05_02 + #0.984 #0.105 TRPSOMG_02_03 + #0.984 &
          #0.105 TRPSOMG_03_03 + #0.984 #0.105 TRPSOMG_04_03 + #0.984 #0.000 &
          TRPSOMG_05_03 + #0.984 #0.105 TRPSOMG_02_04 + #0.984 #0.105 TRPSOMG_03_04 + &
          #0.984 #0.105 TRPSOMG_04_04 + #0.984 #0.000 TRPSOMG_05_04 + #0.984 &
          #0.105 TRPSOMG_03_05 + #0.984 #0.105 TRPSOMG_04_05 + #0.984 #0.105 &
          TRPSOMG_05_05 + #0.984 #0.105 TRPSOMG_03_06 + #0.984 #0.105 TRPSOMG_04_06 + &
          #0.984 #0.105 TRPSOMG_05_06
 S189.5 TRPSOMG_06_06 + OH = #1.000 #0.105 TRPSOMG_01_01 + #1.000 #0.105 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.000 TRPSOMG_05_01 + #1.000 #0.105 TRPSOMG_01_02 + #1.000 &
          #0.105 TRPSOMG_02_02 + #1.000 #0.105 TRPSOMG_03_02 + #1.000 #0.000 &
          TRPSOMG_04_02 + #1.000 #0.000 TRPSOMG_05_02 + #1.000 #0.105 TRPSOMG_02_03 + &
          #1.000 #0.105 TRPSOMG_03_03 + #1.000 #0.000 TRPSOMG_04_03 + #1.000 &
          #0.000 TRPSOMG_05_03 + #1.000 #0.105 TRPSOMG_02_04 + #1.000 #0.105 &
          TRPSOMG_03_04 + #1.000 #0.105 TRPSOMG_04_04 + #1.000 #0.000 TRPSOMG_05_04 + &
          #1.000 #0.105 TRPSOMG_03_05 + #1.000 #0.105 TRPSOMG_04_05 + #1.000 &
          #0.000 TRPSOMG_05_05 + #1.000 #0.105 TRPSOMG_03_06 + #1.000 #0.105 &
          TRPSOMG_04_06 + #1.000 #0.105 TRPSOMG_05_06 + #1.000 #0.105 TRPSOMG_04_07 + &
          #1.000 #0.105 TRPSOMG_05_07
 S190.5 TRPSOMG_06_07 + OH = #1.000 #0.125 TRPSOMG_01_01 + #1.000 #0.125 &
          TRPSOMG_02_01 + #1.000 #0.000 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.000 TRPSOMG_05_01 + #1.000 #0.125 TRPSOMG_01_02 + #1.000 &
          #0.125 TRPSOMG_02_02 + #1.000 #0.000 TRPSOMG_03_02 + #1.000 #0.000 &
          TRPSOMG_04_02 + #1.000 #0.000 TRPSOMG_05_02 + #1.000 #0.125 TRPSOMG_02_03 + &
          #1.000 #0.125 TRPSOMG_03_03 + #1.000 #0.000 TRPSOMG_04_03 + #1.000 &
          #0.000 TRPSOMG_05_03 + #1.000 #0.125 TRPSOMG_02_04 + #1.000 #0.125 &
          TRPSOMG_03_04 + #1.000 #0.000 TRPSOMG_04_04 + #1.000 #0.000 TRPSOMG_05_04 + &
          #1.000 #0.125 TRPSOMG_03_05 + #1.000 #0.125 TRPSOMG_04_05 + #1.000 &
          #0.000 TRPSOMG_05_05 + #1.000 #0.125 TRPSOMG_03_06 + #1.000 #0.125 &
          TRPSOMG_04_06 + #1.000 #0.000 TRPSOMG_05_06 + #1.000 #0.125 TRPSOMG_04_07 + &
          #1.000 #0.125 TRPSOMG_05_07
 S191.5 TRPSOMG_07_01 + OH = #0.156 #0.361 TRPSOMG_07_02 + #0.156 #0.007 &
          TRPSOMG_07_03 + #0.156 #0.606 TRPSOMG_07_04 + #0.156 #0.026 TRPSOMG_07_05 + &
          #0.844 #0.167 TRPSOMG_01_01 + #0.844 #0.167 TRPSOMG_02_01 + #0.844 &
          #0.167 TRPSOMG_03_01 + #0.844 #0.167 TRPSOMG_04_01 + #0.844 #0.167 &
          TRPSOMG_05_01 + #0.844 #0.167 TRPSOMG_06_01 + #0.844 #0.167 TRPSOMG_01_02 + &
          #0.844 #0.167 TRPSOMG_02_02 + #0.844 #0.167 TRPSOMG_03_02 + #0.844 &
          #0.167 TRPSOMG_04_02 + #0.844 #0.167 TRPSOMG_05_02 + #0.844 #0.167 &
          TRPSOMG_06_02
 S192.5 TRPSOMG_07_02 + OH = #0.103 #0.361 TRPSOMG_07_03 + #0.103 #0.007 &
          TRPSOMG_07_04 + #0.103 #0.606 TRPSOMG_07_05 + #0.103 #0.026 TRPSOMG_07_06 + &
          #0.897 #0.125 TRPSOMG_01_01 + #0.897 #0.125 TRPSOMG_02_01 + #0.897 &
          #0.125 TRPSOMG_03_01 + #0.897 #0.125 TRPSOMG_04_01 + #0.897 #0.125 &
          TRPSOMG_05_01 + #0.897 #0.000 TRPSOMG_06_01 + #0.897 #0.125 TRPSOMG_01_02 + &
          #0.897 #0.125 TRPSOMG_02_02 + #0.897 #0.125 TRPSOMG_03_02 + #0.897 &
          #0.125 TRPSOMG_04_02 + #0.897 #0.125 TRPSOMG_05_02 + #0.897 #0.125 &
          TRPSOMG_06_02 + #0.897 #0.125 TRPSOMG_02_03 + #0.897 #0.125 TRPSOMG_03_03 + &
          #0.897 #0.125 TRPSOMG_04_03 + #0.897 #0.125 TRPSOMG_05_03 + #0.897 &
          #0.125 TRPSOMG_06_03
 S193.5 TRPSOMG_07_03 + OH = #0.071 #0.361 TRPSOMG_07_04 + #0.071 #0.007 &
          TRPSOMG_07_05 + #0.071 #0.606 TRPSOMG_07_06 + #0.071 #0.026 TRPSOMG_07_07 + &
          #0.929 #0.100 TRPSOMG_01_01 + #0.929 #0.100 TRPSOMG_02_01 + #0.929 &
          #0.100 TRPSOMG_03_01 + #0.929 #0.100 TRPSOMG_04_01 + #0.929 #0.100 &
          TRPSOMG_05_01 + #0.929 #0.000 TRPSOMG_06_01 + #0.929 #0.100 TRPSOMG_01_02 + &
          #0.929 #0.100 TRPSOMG_02_02 + #0.929 #0.100 TRPSOMG_03_02 + #0.929 &
          #0.100 TRPSOMG_04_02 + #0.929 #0.100 TRPSOMG_05_02 + #0.929 #0.000 &
          TRPSOMG_06_02 + #0.929 #0.100 TRPSOMG_02_03 + #0.929 #0.100 TRPSOMG_03_03 + &
          #0.929 #0.100 TRPSOMG_04_03 + #0.929 #0.100 TRPSOMG_05_03 + #0.929 &
          #0.100 TRPSOMG_06_03 + #0.929 #0.100 TRPSOMG_02_04 + #0.929 #0.100 &
          TRPSOMG_03_04 + #0.929 #0.100 TRPSOMG_04_04 + #0.929 #0.100 TRPSOMG_05_04 + &
          #0.929 #0.100 TRPSOMG_06_04
 S194.5 TRPSOMG_07_04 + OH = #0.048 #0.361 TRPSOMG_07_05 + #0.048 #0.007 &
          TRPSOMG_07_06 + #0.048 #0.632 TRPSOMG_07_07 + #0.952 #0.091 TRPSOMG_01_01 + &
          #0.952 #0.091 TRPSOMG_02_01 + #0.952 #0.091 TRPSOMG_03_01 + #0.952 &
          #0.091 TRPSOMG_04_01 + #0.952 #0.000 TRPSOMG_05_01 + #0.952 #0.000 &
          TRPSOMG_06_01 + #0.952 #0.091 TRPSOMG_01_02 + #0.952 #0.091 TRPSOMG_02_02 + &
          #0.952 #0.091 TRPSOMG_03_02 + #0.952 #0.091 TRPSOMG_04_02 + #0.952 &
          #0.091 TRPSOMG_05_02 + #0.952 #0.000 TRPSOMG_06_02 + #0.952 #0.091 &
          TRPSOMG_02_03 + #0.952 #0.091 TRPSOMG_03_03 + #0.952 #0.091 TRPSOMG_04_03 + &
          #0.952 #0.091 TRPSOMG_05_03 + #0.952 #0.000 TRPSOMG_06_03 + #0.952 &
          #0.091 TRPSOMG_02_04 + #0.952 #0.091 TRPSOMG_03_04 + #0.952 #0.091 &
          TRPSOMG_04_04 + #0.952 #0.091 TRPSOMG_05_04 + #0.952 #0.091 TRPSOMG_06_04 + &
          #0.952 #0.091 TRPSOMG_03_05 + #0.952 #0.091 TRPSOMG_04_05 + #0.952 &
          #0.091 TRPSOMG_05_05 + #0.952 #0.091 TRPSOMG_06_05
 S195.5 TRPSOMG_07_05 + OH = #0.029 #0.361 TRPSOMG_07_06 + #0.029 #0.639 &
          TRPSOMG_07_07 + #0.971 #0.083 TRPSOMG_01_01 + #0.971 #0.083 TRPSOMG_02_01 + &
          #0.971 #0.083 TRPSOMG_03_01 + #0.971 #0.083 TRPSOMG_04_01 + #0.971 &
          #0.000 TRPSOMG_05_01 + #0.971 #0.000 TRPSOMG_06_01 + #0.971 #0.083 &
          TRPSOMG_01_02 + #0.971 #0.083 TRPSOMG_02_02 + #0.971 #0.083 TRPSOMG_03_02 + &
          #0.971 #0.083 TRPSOMG_04_02 + #0.971 #0.000 TRPSOMG_05_02 + #0.971 &
          #0.000 TRPSOMG_06_02 + #0.971 #0.083 TRPSOMG_02_03 + #0.971 #0.083 &
          TRPSOMG_03_03 + #0.971 #0.083 TRPSOMG_04_03 + #0.971 #0.083 TRPSOMG_05_03 + &
          #0.971 #0.000 TRPSOMG_06_03 + #0.971 #0.083 TRPSOMG_02_04 + #0.971 &
          #0.083 TRPSOMG_03_04 + #0.971 #0.083 TRPSOMG_04_04 + #0.971 #0.083 &
          TRPSOMG_05_04 + #0.971 #0.000 TRPSOMG_06_04 + #0.971 #0.083 TRPSOMG_03_05 + &
          #0.971 #0.083 TRPSOMG_04_05 + #0.971 #0.083 TRPSOMG_05_05 + #0.971 &
          #0.083 TRPSOMG_06_05 + #0.971 #0.083 TRPSOMG_03_06 + #0.971 #0.083 &
          TRPSOMG_04_06 + #0.971 #0.083 TRPSOMG_05_06 + #0.971 #0.083 TRPSOMG_06_06
 S196.5 TRPSOMG_07_06 + OH = #0.013 #1.000 TRPSOMG_07_07 + #0.987 #0.080 &
          TRPSOMG_01_01 + #0.987 #0.080 TRPSOMG_02_01 + #0.987 #0.080 TRPSOMG_03_01 + &
          #0.987 #0.000 TRPSOMG_04_01 + #0.987 #0.000 TRPSOMG_05_01 + #0.987 &
          #0.000 TRPSOMG_06_01 + #0.987 #0.080 TRPSOMG_01_02 + #0.987 #0.080 &
          TRPSOMG_02_02 + #0.987 #0.080 TRPSOMG_03_02 + #0.987 #0.080 TRPSOMG_04_02 + &
          #0.987 #0.000 TRPSOMG_05_02 + #0.987 #0.000 TRPSOMG_06_02 + #0.987 &
          #0.080 TRPSOMG_02_03 + #0.987 #0.080 TRPSOMG_03_03 + #0.987 #0.080 &
          TRPSOMG_04_03 + #0.987 #0.000 TRPSOMG_05_03 + #0.987 #0.000 TRPSOMG_06_03 + &
          #0.987 #0.080 TRPSOMG_02_04 + #0.987 #0.080 TRPSOMG_03_04 + #0.987 &
          #0.080 TRPSOMG_04_04 + #0.987 #0.080 TRPSOMG_05_04 + #0.987 #0.000 &
          TRPSOMG_06_04 + #0.987 #0.080 TRPSOMG_03_05 + #0.987 #0.080 TRPSOMG_04_05 + &
          #0.987 #0.080 TRPSOMG_05_05 + #0.987 #0.000 TRPSOMG_06_05 + #0.987 &
          #0.080 TRPSOMG_03_06 + #0.987 #0.080 TRPSOMG_04_06 + #0.987 #0.080 &
          TRPSOMG_05_06 + #0.987 #0.080 TRPSOMG_06_06 + #0.987 #0.080 TRPSOMG_04_07 + &
          #0.987 #0.080 TRPSOMG_05_07 + #0.987 #0.080 TRPSOMG_06_07
 S197.5 TRPSOMG_07_07 + OH = #1.000 #0.091 TRPSOMG_01_01 + #1.000 #0.091 &
          TRPSOMG_02_01 + #1.000 #0.091 TRPSOMG_03_01 + #1.000 #0.000 TRPSOMG_04_01 + &
          #1.000 #0.000 TRPSOMG_05_01 + #1.000 #0.000 TRPSOMG_06_01 + #1.000 &
          #0.091 TRPSOMG_01_02 + #1.000 #0.091 TRPSOMG_02_02 + #1.000 #0.091 &
          TRPSOMG_03_02 + #1.000 #0.000 TRPSOMG_04_02 + #1.000 #0.000 TRPSOMG_05_02 + &
          #1.000 #0.000 TRPSOMG_06_02 + #1.000 #0.091 TRPSOMG_02_03 + #1.000 &
          #0.091 TRPSOMG_03_03 + #1.000 #0.091 TRPSOMG_04_03 + #1.000 #0.000 &
          TRPSOMG_05_03 + #1.000 #0.000 TRPSOMG_06_03 + #1.000 #0.091 TRPSOMG_02_04 + &
          #1.000 #0.091 TRPSOMG_03_04 + #1.000 #0.091 TRPSOMG_04_04 + #1.000 &
          #0.000 TRPSOMG_05_04 + #1.000 #0.000 TRPSOMG_06_04 + #1.000 #0.091 &
          TRPSOMG_03_05 + #1.000 #0.091 TRPSOMG_04_05 + #1.000 #0.091 TRPSOMG_05_05 + &
          #1.000 #0.000 TRPSOMG_06_05 + #1.000 #0.091 TRPSOMG_03_06 + #1.000 &
          #0.091 TRPSOMG_04_06 + #1.000 #0.091 TRPSOMG_05_06 + #1.000 #0.000 &
          TRPSOMG_06_06 + #1.000 #0.091 TRPSOMG_04_07 + #1.000 #0.091 TRPSOMG_05_07 + &
          #1.000 #0.091 TRPSOMG_06_07
 S198.5 TRPSOMG_08_01 + OH = #0.165 #0.361 TRPSOMG_08_02 + #0.165 #0.007 &
          TRPSOMG_08_03 + #0.165 #0.606 TRPSOMG_08_04 + #0.165 #0.026 TRPSOMG_08_05 + &
          #0.835 #0.143 TRPSOMG_01_01 + #0.835 #0.143 TRPSOMG_02_01 + #0.835 &
          #0.143 TRPSOMG_03_01 + #0.835 #0.143 TRPSOMG_04_01 + #0.835 #0.143 &
          TRPSOMG_05_01 + #0.835 #0.143 TRPSOMG_06_01 + #0.835 #0.143 TRPSOMG_07_01 + &
          #0.835 #0.143 TRPSOMG_01_02 + #0.835 #0.143 TRPSOMG_02_02 + #0.835 &
          #0.143 TRPSOMG_03_02 + #0.835 #0.143 TRPSOMG_04_02 + #0.835 #0.143 &
          TRPSOMG_05_02 + #0.835 #0.143 TRPSOMG_06_02 + #0.835 #0.143 TRPSOMG_07_02
 S199.5 TRPSOMG_08_02 + OH = #0.114 #0.361 TRPSOMG_08_03 + #0.114 #0.007 &
          TRPSOMG_08_04 + #0.114 #0.606 TRPSOMG_08_05 + #0.114 #0.026 TRPSOMG_08_06 + &
          #0.886 #0.105 TRPSOMG_01_01 + #0.886 #0.105 TRPSOMG_02_01 + #0.886 &
          #0.105 TRPSOMG_03_01 + #0.886 #0.105 TRPSOMG_04_01 + #0.886 #0.105 &
          TRPSOMG_05_01 + #0.886 #0.105 TRPSOMG_06_01 + #0.886 #0.000 TRPSOMG_07_01 + &
          #0.886 #0.105 TRPSOMG_01_02 + #0.886 #0.105 TRPSOMG_02_02 + #0.886 &
          #0.105 TRPSOMG_03_02 + #0.886 #0.105 TRPSOMG_04_02 + #0.886 #0.105 &
          TRPSOMG_05_02 + #0.886 #0.105 TRPSOMG_06_02 + #0.886 #0.105 TRPSOMG_07_02 + &
          #0.886 #0.105 TRPSOMG_02_03 + #0.886 #0.105 TRPSOMG_03_03 + #0.886 &
          #0.105 TRPSOMG_04_03 + #0.886 #0.105 TRPSOMG_05_03 + #0.886 #0.105 &
          TRPSOMG_06_03 + #0.886 #0.105 TRPSOMG_07_03
 S200.5 TRPSOMG_08_03 + OH = #0.082 #0.361 TRPSOMG_08_04 + #0.082 #0.007 &
          TRPSOMG_08_05 + #0.082 #0.606 TRPSOMG_08_06 + #0.082 #0.026 TRPSOMG_08_07 + &
          #0.918 #0.083 TRPSOMG_01_01 + #0.918 #0.083 TRPSOMG_02_01 + #0.918 &
          #0.083 TRPSOMG_03_01 + #0.918 #0.083 TRPSOMG_04_01 + #0.918 #0.083 &
          TRPSOMG_05_01 + #0.918 #0.083 TRPSOMG_06_01 + #0.918 #0.000 TRPSOMG_07_01 + &
          #0.918 #0.083 TRPSOMG_01_02 + #0.918 #0.083 TRPSOMG_02_02 + #0.918 &
          #0.083 TRPSOMG_03_02 + #0.918 #0.083 TRPSOMG_04_02 + #0.918 #0.083 &
          TRPSOMG_05_02 + #0.918 #0.083 TRPSOMG_06_02 + #0.918 #0.000 TRPSOMG_07_02 + &
          #0.918 #0.083 TRPSOMG_02_03 + #0.918 #0.083 TRPSOMG_03_03 + #0.918 &
          #0.083 TRPSOMG_04_03 + #0.918 #0.083 TRPSOMG_05_03 + #0.918 #0.083 &
          TRPSOMG_06_03 + #0.918 #0.083 TRPSOMG_07_03 + #0.918 #0.083 TRPSOMG_02_04 + &
          #0.918 #0.083 TRPSOMG_03_04 + #0.918 #0.083 TRPSOMG_04_04 + #0.918 &
          #0.083 TRPSOMG_05_04 + #0.918 #0.083 TRPSOMG_06_04 + #0.918 #0.083 &
          TRPSOMG_07_04
 S201.5 TRPSOMG_08_04 + OH = #0.059 #0.361 TRPSOMG_08_05 + #0.059 #0.007 &
          TRPSOMG_08_06 + #0.059 #0.632 TRPSOMG_08_07 + #0.941 #0.074 TRPSOMG_01_01 + &
          #0.941 #0.074 TRPSOMG_02_01 + #0.941 #0.074 TRPSOMG_03_01 + #0.941 &
          #0.074 TRPSOMG_04_01 + #0.941 #0.074 TRPSOMG_05_01 + #0.941 #0.000 &
          TRPSOMG_06_01 + #0.941 #0.000 TRPSOMG_07_01 + #0.941 #0.074 TRPSOMG_01_02 + &
          #0.941 #0.074 TRPSOMG_02_02 + #0.941 #0.074 TRPSOMG_03_02 + #0.941 &
          #0.074 TRPSOMG_04_02 + #0.941 #0.074 TRPSOMG_05_02 + #0.941 #0.074 &
          TRPSOMG_06_02 + #0.941 #0.000 TRPSOMG_07_02 + #0.941 #0.074 TRPSOMG_02_03 + &
          #0.941 #0.074 TRPSOMG_03_03 + #0.941 #0.074 TRPSOMG_04_03 + #0.941 &
          #0.074 TRPSOMG_05_03 + #0.941 #0.074 TRPSOMG_06_03 + #0.941 #0.000 &
          TRPSOMG_07_03 + #0.941 #0.074 TRPSOMG_02_04 + #0.941 #0.074 TRPSOMG_03_04 + &
          #0.941 #0.074 TRPSOMG_04_04 + #0.941 #0.074 TRPSOMG_05_04 + #0.941 &
          #0.074 TRPSOMG_06_04 + #0.941 #0.074 TRPSOMG_07_04 + #0.941 #0.074 &
          TRPSOMG_03_05 + #0.941 #0.074 TRPSOMG_04_05 + #0.941 #0.074 TRPSOMG_05_05 + &
          #0.941 #0.074 TRPSOMG_06_05 + #0.941 #0.074 TRPSOMG_07_05
 S202.5 TRPSOMG_08_05 + OH = #0.040 #0.361 TRPSOMG_08_06 + #0.040 #0.639 &
          TRPSOMG_08_07 + #0.960 #0.067 TRPSOMG_01_01 + #0.960 #0.067 TRPSOMG_02_01 + &
          #0.960 #0.067 TRPSOMG_03_01 + #0.960 #0.067 TRPSOMG_04_01 + #0.960 &
          #0.067 TRPSOMG_05_01 + #0.960 #0.000 TRPSOMG_06_01 + #0.960 #0.000 &
          TRPSOMG_07_01 + #0.960 #0.067 TRPSOMG_01_02 + #0.960 #0.067 TRPSOMG_02_02 + &
          #0.960 #0.067 TRPSOMG_03_02 + #0.960 #0.067 TRPSOMG_04_02 + #0.960 &
          #0.067 TRPSOMG_05_02 + #0.960 #0.000 TRPSOMG_06_02 + #0.960 #0.000 &
          TRPSOMG_07_02 + #0.960 #0.067 TRPSOMG_02_03 + #0.960 #0.067 TRPSOMG_03_03 + &
          #0.960 #0.067 TRPSOMG_04_03 + #0.960 #0.067 TRPSOMG_05_03 + #0.960 &
          #0.067 TRPSOMG_06_03 + #0.960 #0.000 TRPSOMG_07_03 + #0.960 #0.067 &
          TRPSOMG_02_04 + #0.960 #0.067 TRPSOMG_03_04 + #0.960 #0.067 TRPSOMG_04_04 + &
          #0.960 #0.067 TRPSOMG_05_04 + #0.960 #0.067 TRPSOMG_06_04 + #0.960 &
          #0.000 TRPSOMG_07_04 + #0.960 #0.067 TRPSOMG_03_05 + #0.960 #0.067 &
          TRPSOMG_04_05 + #0.960 #0.067 TRPSOMG_05_05 + #0.960 #0.067 TRPSOMG_06_05 + &
          #0.960 #0.067 TRPSOMG_07_05 + #0.960 #0.067 TRPSOMG_03_06 + #0.960 &
          #0.067 TRPSOMG_04_06 + #0.960 #0.067 TRPSOMG_05_06 + #0.960 #0.067 &
          TRPSOMG_06_06 + #0.960 #0.067 TRPSOMG_07_06
 S203.5 TRPSOMG_08_06 + OH = #0.025 #1.000 TRPSOMG_08_07 + #0.975 #0.065 &
          TRPSOMG_01_01 + #0.975 #0.065 TRPSOMG_02_01 + #0.975 #0.065 TRPSOMG_03_01 + &
          #0.975 #0.065 TRPSOMG_04_01 + #0.975 #0.000 TRPSOMG_05_01 + #0.975 &
          #0.000 TRPSOMG_06_01 + #0.975 #0.000 TRPSOMG_07_01 + #0.975 #0.065 &
          TRPSOMG_01_02 + #0.975 #0.065 TRPSOMG_02_02 + #0.975 #0.065 TRPSOMG_03_02 + &
          #0.975 #0.065 TRPSOMG_04_02 + #0.975 #0.065 TRPSOMG_05_02 + #0.975 &
          #0.000 TRPSOMG_06_02 + #0.975 #0.000 TRPSOMG_07_02 + #0.975 #0.065 &
          TRPSOMG_02_03 + #0.975 #0.065 TRPSOMG_03_03 + #0.975 #0.065 TRPSOMG_04_03 + &
          #0.975 #0.065 TRPSOMG_05_03 + #0.975 #0.000 TRPSOMG_06_03 + #0.975 &
          #0.000 TRPSOMG_07_03 + #0.975 #0.065 TRPSOMG_02_04 + #0.975 #0.065 &
          TRPSOMG_03_04 + #0.975 #0.065 TRPSOMG_04_04 + #0.975 #0.065 TRPSOMG_05_04 + &
          #0.975 #0.065 TRPSOMG_06_04 + #0.975 #0.000 TRPSOMG_07_04 + #0.975 &
          #0.065 TRPSOMG_03_05 + #0.975 #0.065 TRPSOMG_04_05 + #0.975 #0.065 &
          TRPSOMG_05_05 + #0.975 #0.065 TRPSOMG_06_05 + #0.975 #0.000 TRPSOMG_07_05 + &
          #0.975 #0.065 TRPSOMG_03_06 + #0.975 #0.065 TRPSOMG_04_06 + #0.975 &
          #0.065 TRPSOMG_05_06 + #0.975 #0.065 TRPSOMG_06_06 + #0.975 #0.065 &
          TRPSOMG_07_06 + #0.975 #0.065 TRPSOMG_04_07 + #0.975 #0.065 TRPSOMG_05_07 + &
          #0.975 #0.065 TRPSOMG_06_07 + #0.975 #0.065 TRPSOMG_07_07
 S204.5 TRPSOMG_08_07 + OH = #1.000 #0.071 TRPSOMG_01_01 + #1.000 #0.071 &
          TRPSOMG_02_01 + #1.000 #0.071 TRPSOMG_03_01 + #1.000 #0.071 TRPSOMG_04_01 + &
          #1.000 #0.000 TRPSOMG_05_01 + #1.000 #0.000 TRPSOMG_06_01 + #1.000 &
          #0.000 TRPSOMG_07_01 + #1.000 #0.071 TRPSOMG_01_02 + #1.000 #0.071 &
          TRPSOMG_02_02 + #1.000 #0.071 TRPSOMG_03_02 + #1.000 #0.071 TRPSOMG_04_02 + &
          #1.000 #0.000 TRPSOMG_05_02 + #1.000 #0.000 TRPSOMG_06_02 + #1.000 &
          #0.000 TRPSOMG_07_02 + #1.000 #0.071 TRPSOMG_02_03 + #1.000 #0.071 &
          TRPSOMG_03_03 + #1.000 #0.071 TRPSOMG_04_03 + #1.000 #0.071 TRPSOMG_05_03 + &
          #1.000 #0.000 TRPSOMG_06_03 + #1.000 #0.000 TRPSOMG_07_03 + #1.000 &
          #0.071 TRPSOMG_02_04 + #1.000 #0.071 TRPSOMG_03_04 + #1.000 #0.071 &
          TRPSOMG_04_04 + #1.000 #0.071 TRPSOMG_05_04 + #1.000 #0.000 TRPSOMG_06_04 + &
          #1.000 #0.000 TRPSOMG_07_04 + #1.000 #0.071 TRPSOMG_03_05 + #1.000 &
          #0.071 TRPSOMG_04_05 + #1.000 #0.071 TRPSOMG_05_05 + #1.000 #0.071 &
          TRPSOMG_06_05 + #1.000 #0.000 TRPSOMG_07_05 + #1.000 #0.071 TRPSOMG_03_06 + &
          #1.000 #0.071 TRPSOMG_04_06 + #1.000 #0.071 TRPSOMG_05_06 + #1.000 &
          #0.071 TRPSOMG_06_06 + #1.000 #0.000 TRPSOMG_07_06 + #1.000 #0.071 &
          TRPSOMG_04_07 + #1.000 #0.071 TRPSOMG_05_07 + #1.000 #0.071 TRPSOMG_06_07 + &
          #1.000 #0.071 TRPSOMG_07_07
 S205.5 TRPSOMG_09_01 + OH = #0.174 #0.361 TRPSOMG_09_02 + #0.174 #0.007 &
          TRPSOMG_09_03 + #0.174 #0.606 TRPSOMG_09_04 + #0.174 #0.026 TRPSOMG_09_05 + &
          #0.826 #0.125 TRPSOMG_01_01 + #0.826 #0.125 TRPSOMG_02_01 + #0.826 &
          #0.125 TRPSOMG_03_01 + #0.826 #0.125 TRPSOMG_04_01 + #0.826 #0.125 &
          TRPSOMG_05_01 + #0.826 #0.125 TRPSOMG_06_01 + #0.826 #0.125 TRPSOMG_07_01 + &
          #0.826 #0.125 TRPSOMG_08_01 + #0.826 #0.125 TRPSOMG_01_02 + #0.826 &
          #0.125 TRPSOMG_02_02 + #0.826 #0.125 TRPSOMG_03_02 + #0.826 #0.125 &
          TRPSOMG_04_02 + #0.826 #0.125 TRPSOMG_05_02 + #0.826 #0.125 TRPSOMG_06_02 + &
          #0.826 #0.125 TRPSOMG_07_02 + #0.826 #0.125 TRPSOMG_08_02
 S206.5 TRPSOMG_09_02 + OH = #0.123 #0.361 TRPSOMG_09_03 + #0.123 #0.007 &
          TRPSOMG_09_04 + #0.123 #0.606 TRPSOMG_09_05 + #0.123 #0.026 TRPSOMG_09_06 + &
          #0.877 #0.087 TRPSOMG_01_01 + #0.877 #0.087 TRPSOMG_02_01 + #0.877 &
          #0.087 TRPSOMG_03_01 + #0.877 #0.087 TRPSOMG_04_01 + #0.877 #0.087 &
          TRPSOMG_05_01 + #0.877 #0.087 TRPSOMG_06_01 + #0.877 #0.087 TRPSOMG_07_01 + &
          #0.877 #0.000 TRPSOMG_08_01 + #0.877 #0.087 TRPSOMG_01_02 + #0.877 &
          #0.087 TRPSOMG_02_02 + #0.877 #0.087 TRPSOMG_03_02 + #0.877 #0.087 &
          TRPSOMG_04_02 + #0.877 #0.087 TRPSOMG_05_02 + #0.877 #0.087 TRPSOMG_06_02 + &
          #0.877 #0.087 TRPSOMG_07_02 + #0.877 #0.087 TRPSOMG_08_02 + #0.877 &
          #0.087 TRPSOMG_02_03 + #0.877 #0.087 TRPSOMG_03_03 + #0.877 #0.087 &
          TRPSOMG_04_03 + #0.877 #0.087 TRPSOMG_05_03 + #0.877 #0.087 TRPSOMG_06_03 + &
          #0.877 #0.087 TRPSOMG_07_03 + #0.877 #0.087 TRPSOMG_08_03
 S207.5 TRPSOMG_09_03 + OH = #0.091 #0.361 TRPSOMG_09_04 + #0.091 #0.007 &
          TRPSOMG_09_05 + #0.091 #0.606 TRPSOMG_09_06 + #0.091 #0.026 TRPSOMG_09_07 + &
          #0.909 #0.069 TRPSOMG_01_01 + #0.909 #0.069 TRPSOMG_02_01 + #0.909 &
          #0.069 TRPSOMG_03_01 + #0.909 #0.069 TRPSOMG_04_01 + #0.909 #0.069 &
          TRPSOMG_05_01 + #0.909 #0.069 TRPSOMG_06_01 + #0.909 #0.069 TRPSOMG_07_01 + &
          #0.909 #0.000 TRPSOMG_08_01 + #0.909 #0.069 TRPSOMG_01_02 + #0.909 &
          #0.069 TRPSOMG_02_02 + #0.909 #0.069 TRPSOMG_03_02 + #0.909 #0.069 &
          TRPSOMG_04_02 + #0.909 #0.069 TRPSOMG_05_02 + #0.909 #0.069 TRPSOMG_06_02 + &
          #0.909 #0.069 TRPSOMG_07_02 + #0.909 #0.000 TRPSOMG_08_02 + #0.909 &
          #0.069 TRPSOMG_02_03 + #0.909 #0.069 TRPSOMG_03_03 + #0.909 #0.069 &
          TRPSOMG_04_03 + #0.909 #0.069 TRPSOMG_05_03 + #0.909 #0.069 TRPSOMG_06_03 + &
          #0.909 #0.069 TRPSOMG_07_03 + #0.909 #0.069 TRPSOMG_08_03 + #0.909 &
          #0.069 TRPSOMG_02_04 + #0.909 #0.069 TRPSOMG_03_04 + #0.909 #0.069 &
          TRPSOMG_04_04 + #0.909 #0.069 TRPSOMG_05_04 + #0.909 #0.069 TRPSOMG_06_04 + &
          #0.909 #0.069 TRPSOMG_07_04 + #0.909 #0.069 TRPSOMG_08_04
 S208.5 TRPSOMG_09_04 + OH = #0.068 #0.361 TRPSOMG_09_05 + #0.068 #0.007 &
          TRPSOMG_09_06 + #0.068 #0.632 TRPSOMG_09_07 + #0.932 #0.059 TRPSOMG_01_01 + &
          #0.932 #0.059 TRPSOMG_02_01 + #0.932 #0.059 TRPSOMG_03_01 + #0.932 &
          #0.059 TRPSOMG_04_01 + #0.932 #0.059 TRPSOMG_05_01 + #0.932 #0.059 &
          TRPSOMG_06_01 + #0.932 #0.000 TRPSOMG_07_01 + #0.932 #0.000 TRPSOMG_08_01 + &
          #0.932 #0.059 TRPSOMG_01_02 + #0.932 #0.059 TRPSOMG_02_02 + #0.932 &
          #0.059 TRPSOMG_03_02 + #0.932 #0.059 TRPSOMG_04_02 + #0.932 #0.059 &
          TRPSOMG_05_02 + #0.932 #0.059 TRPSOMG_06_02 + #0.932 #0.059 TRPSOMG_07_02 + &
          #0.932 #0.000 TRPSOMG_08_02 + #0.932 #0.059 TRPSOMG_02_03 + #0.932 &
          #0.059 TRPSOMG_03_03 + #0.932 #0.059 TRPSOMG_04_03 + #0.932 #0.059 &
          TRPSOMG_05_03 + #0.932 #0.059 TRPSOMG_06_03 + #0.932 #0.059 TRPSOMG_07_03 + &
          #0.932 #0.000 TRPSOMG_08_03 + #0.932 #0.059 TRPSOMG_02_04 + #0.932 &
          #0.059 TRPSOMG_03_04 + #0.932 #0.059 TRPSOMG_04_04 + #0.932 #0.059 &
          TRPSOMG_05_04 + #0.932 #0.059 TRPSOMG_06_04 + #0.932 #0.059 TRPSOMG_07_04 + &
          #0.932 #0.059 TRPSOMG_08_04 + #0.932 #0.059 TRPSOMG_03_05 + #0.932 &
          #0.059 TRPSOMG_04_05 + #0.932 #0.059 TRPSOMG_05_05 + #0.932 #0.059 &
          TRPSOMG_06_05 + #0.932 #0.059 TRPSOMG_07_05 + #0.932 #0.059 TRPSOMG_08_05
 S209.5 TRPSOMG_09_05 + OH = #0.050 #0.361 TRPSOMG_09_06 + #0.050 #0.639 &
          TRPSOMG_09_07 + #0.950 #0.053 TRPSOMG_01_01 + #0.950 #0.053 TRPSOMG_02_01 + &
          #0.950 #0.053 TRPSOMG_03_01 + #0.950 #0.053 TRPSOMG_04_01 + #0.950 &
          #0.053 TRPSOMG_05_01 + #0.950 #0.053 TRPSOMG_06_01 + #0.950 #0.000 &
          TRPSOMG_07_01 + #0.950 #0.000 TRPSOMG_08_01 + #0.950 #0.053 TRPSOMG_01_02 + &
          #0.950 #0.053 TRPSOMG_02_02 + #0.950 #0.053 TRPSOMG_03_02 + #0.950 &
          #0.053 TRPSOMG_04_02 + #0.950 #0.053 TRPSOMG_05_02 + #0.950 #0.053 &
          TRPSOMG_06_02 + #0.950 #0.000 TRPSOMG_07_02 + #0.950 #0.000 TRPSOMG_08_02 + &
          #0.950 #0.053 TRPSOMG_02_03 + #0.950 #0.053 TRPSOMG_03_03 + #0.950 &
          #0.053 TRPSOMG_04_03 + #0.950 #0.053 TRPSOMG_05_03 + #0.950 #0.053 &
          TRPSOMG_06_03 + #0.950 #0.053 TRPSOMG_07_03 + #0.950 #0.000 TRPSOMG_08_03 + &
          #0.950 #0.053 TRPSOMG_02_04 + #0.950 #0.053 TRPSOMG_03_04 + #0.950 &
          #0.053 TRPSOMG_04_04 + #0.950 #0.053 TRPSOMG_05_04 + #0.950 #0.053 &
          TRPSOMG_06_04 + #0.950 #0.053 TRPSOMG_07_04 + #0.950 #0.000 TRPSOMG_08_04 + &
          #0.950 #0.053 TRPSOMG_03_05 + #0.950 #0.053 TRPSOMG_04_05 + #0.950 &
          #0.053 TRPSOMG_05_05 + #0.950 #0.053 TRPSOMG_06_05 + #0.950 #0.053 &
          TRPSOMG_07_05 + #0.950 #0.053 TRPSOMG_08_05 + #0.950 #0.053 TRPSOMG_03_06 + &
          #0.950 #0.053 TRPSOMG_04_06 + #0.950 #0.053 TRPSOMG_05_06 + #0.950 &
          #0.053 TRPSOMG_06_06 + #0.950 #0.053 TRPSOMG_07_06 + #0.950 #0.053 &
          TRPSOMG_08_06
 S210.5 TRPSOMG_09_06 + OH = #0.035 #1.000 TRPSOMG_09_07 + #0.965 #0.049 &
          TRPSOMG_01_01 + #0.965 #0.049 TRPSOMG_02_01 + #0.965 #0.049 TRPSOMG_03_01 + &
          #0.965 #0.049 TRPSOMG_04_01 + #0.965 #0.049 TRPSOMG_05_01 + #0.965 &
          #0.000 TRPSOMG_06_01 + #0.965 #0.000 TRPSOMG_07_01 + #0.965 #0.000 &
          TRPSOMG_08_01 + #0.965 #0.049 TRPSOMG_01_02 + #0.965 #0.049 TRPSOMG_02_02 + &
          #0.965 #0.049 TRPSOMG_03_02 + #0.965 #0.049 TRPSOMG_04_02 + #0.965 &
          #0.049 TRPSOMG_05_02 + #0.965 #0.049 TRPSOMG_06_02 + #0.965 #0.000 &
          TRPSOMG_07_02 + #0.965 #0.000 TRPSOMG_08_02 + #0.965 #0.049 TRPSOMG_02_03 + &
          #0.965 #0.049 TRPSOMG_03_03 + #0.965 #0.049 TRPSOMG_04_03 + #0.965 &
          #0.049 TRPSOMG_05_03 + #0.965 #0.049 TRPSOMG_06_03 + #0.965 #0.000 &
          TRPSOMG_07_03 + #0.965 #0.000 TRPSOMG_08_03 + #0.965 #0.049 TRPSOMG_02_04 + &
          #0.965 #0.049 TRPSOMG_03_04 + #0.965 #0.049 TRPSOMG_04_04 + #0.965 &
          #0.049 TRPSOMG_05_04 + #0.965 #0.049 TRPSOMG_06_04 + #0.965 #0.049 &
          TRPSOMG_07_04 + #0.965 #0.000 TRPSOMG_08_04 + #0.965 #0.049 TRPSOMG_03_05 + &
          #0.965 #0.049 TRPSOMG_04_05 + #0.965 #0.049 TRPSOMG_05_05 + #0.965 &
          #0.049 TRPSOMG_06_05 + #0.965 #0.049 TRPSOMG_07_05 + #0.965 #0.000 &
          TRPSOMG_08_05 + #0.965 #0.049 TRPSOMG_03_06 + #0.965 #0.049 TRPSOMG_04_06 + &
          #0.965 #0.049 TRPSOMG_05_06 + #0.965 #0.049 TRPSOMG_06_06 + #0.965 &
          #0.049 TRPSOMG_07_06 + #0.965 #0.049 TRPSOMG_08_06 + #0.965 #0.049 &
          TRPSOMG_04_07 + #0.965 #0.049 TRPSOMG_05_07 + #0.965 #0.049 TRPSOMG_06_07 + &
          #0.965 #0.049 TRPSOMG_07_07 + #0.965 #0.049 TRPSOMG_08_07
 S211.5 TRPSOMG_09_07 + OH = #1.000 #0.053 TRPSOMG_01_01 + #1.000 #0.053 &
          TRPSOMG_02_01 + #1.000 #0.053 TRPSOMG_03_01 + #1.000 #0.053 TRPSOMG_04_01 + &
          #1.000 #0.053 TRPSOMG_05_01 + #1.000 #0.000 TRPSOMG_06_01 + #1.000 &
          #0.000 TRPSOMG_07_01 + #1.000 #0.000 TRPSOMG_08_01 + #1.000 #0.053 &
          TRPSOMG_01_02 + #1.000 #0.053 TRPSOMG_02_02 + #1.000 #0.053 TRPSOMG_03_02 + &
          #1.000 #0.053 TRPSOMG_04_02 + #1.000 #0.053 TRPSOMG_05_02 + #1.000 &
          #0.000 TRPSOMG_06_02 + #1.000 #0.000 TRPSOMG_07_02 + #1.000 #0.000 &
          TRPSOMG_08_02 + #1.000 #0.053 TRPSOMG_02_03 + #1.000 #0.053 TRPSOMG_03_03 + &
          #1.000 #0.053 TRPSOMG_04_03 + #1.000 #0.053 TRPSOMG_05_03 + #1.000 &
          #0.053 TRPSOMG_06_03 + #1.000 #0.000 TRPSOMG_07_03 + #1.000 #0.000 &
          TRPSOMG_08_03 + #1.000 #0.053 TRPSOMG_02_04 + #1.000 #0.053 TRPSOMG_03_04 + &
          #1.000 #0.053 TRPSOMG_04_04 + #1.000 #0.053 TRPSOMG_05_04 + #1.000 &
          #0.053 TRPSOMG_06_04 + #1.000 #0.000 TRPSOMG_07_04 + #1.000 #0.000 &
          TRPSOMG_08_04 + #1.000 #0.053 TRPSOMG_03_05 + #1.000 #0.053 TRPSOMG_04_05 + &
          #1.000 #0.053 TRPSOMG_05_05 + #1.000 #0.053 TRPSOMG_06_05 + #1.000 &
          #0.053 TRPSOMG_07_05 + #1.000 #0.000 TRPSOMG_08_05 + #1.000 #0.053 &
          TRPSOMG_03_06 + #1.000 #0.053 TRPSOMG_04_06 + #1.000 #0.053 TRPSOMG_05_06 + &
          #1.000 #0.053 TRPSOMG_06_06 + #1.000 #0.053 TRPSOMG_07_06 + #1.000 &
          #0.000 TRPSOMG_08_06 + #1.000 #0.053 TRPSOMG_04_07 + #1.000 #0.053 &
          TRPSOMG_05_07 + #1.000 #0.053 TRPSOMG_06_07 + #1.000 #0.053 TRPSOMG_07_07 + &
          #1.000 #0.053 TRPSOMG_08_07
 S212.5 TRPSOMG_10_01 + OH = #0.182 #0.361 TRPSOMG_10_02 + #0.182 #0.007 &
          TRPSOMG_10_03 + #0.182 #0.606 TRPSOMG_10_04 + #0.182 #0.026 TRPSOMG_10_05 + &
          #0.818 #0.111 TRPSOMG_01_01 + #0.818 #0.111 TRPSOMG_02_01 + #0.818 &
          #0.111 TRPSOMG_03_01 + #0.818 #0.111 TRPSOMG_04_01 + #0.818 #0.111 &
          TRPSOMG_05_01 + #0.818 #0.111 TRPSOMG_06_01 + #0.818 #0.111 TRPSOMG_07_01 + &
          #0.818 #0.111 TRPSOMG_08_01 + #0.818 #0.111 TRPSOMG_09_01 + #0.818 &
          #0.111 TRPSOMG_01_02 + #0.818 #0.111 TRPSOMG_02_02 + #0.818 #0.111 &
          TRPSOMG_03_02 + #0.818 #0.111 TRPSOMG_04_02 + #0.818 #0.111 TRPSOMG_05_02 + &
          #0.818 #0.111 TRPSOMG_06_02 + #0.818 #0.111 TRPSOMG_07_02 + #0.818 &
          #0.111 TRPSOMG_08_02 + #0.818 #0.111 TRPSOMG_09_02
 S213.5 TRPSOMG_10_02 + OH = #0.131 #0.361 TRPSOMG_10_03 + #0.131 #0.007 &
          TRPSOMG_10_04 + #0.131 #0.606 TRPSOMG_10_05 + #0.131 #0.026 TRPSOMG_10_06 + &
          #0.869 #0.080 TRPSOMG_01_01 + #0.869 #0.080 TRPSOMG_02_01 + #0.869 &
          #0.080 TRPSOMG_03_01 + #0.869 #0.080 TRPSOMG_04_01 + #0.869 #0.080 &
          TRPSOMG_05_01 + #0.869 #0.080 TRPSOMG_06_01 + #0.869 #0.080 TRPSOMG_07_01 + &
          #0.869 #0.080 TRPSOMG_08_01 + #0.869 #0.000 TRPSOMG_09_01 + #0.869 &
          #0.080 TRPSOMG_01_02 + #0.869 #0.080 TRPSOMG_02_02 + #0.869 #0.080 &
          TRPSOMG_03_02 + #0.869 #0.080 TRPSOMG_04_02 + #0.869 #0.080 TRPSOMG_05_02 + &
          #0.869 #0.080 TRPSOMG_06_02 + #0.869 #0.080 TRPSOMG_07_02 + #0.869 &
          #0.080 TRPSOMG_08_02 + #0.869 #0.080 TRPSOMG_09_02 + #0.869 #0.080 &
          TRPSOMG_02_03 + #0.869 #0.080 TRPSOMG_03_03 + #0.869 #0.080 TRPSOMG_04_03 + &
          #0.869 #0.080 TRPSOMG_05_03 + #0.869 #0.080 TRPSOMG_06_03 + #0.869 &
          #0.080 TRPSOMG_07_03 + #0.869 #0.080 TRPSOMG_08_03 + #0.869 #0.080 &
          TRPSOMG_09_03
 S214.5 TRPSOMG_10_03 + OH = #0.099 #0.361 TRPSOMG_10_04 + #0.099 #0.007 &
          TRPSOMG_10_05 + #0.099 #0.606 TRPSOMG_10_06 + #0.099 #0.026 TRPSOMG_10_07 + &
          #0.901 #0.062 TRPSOMG_01_01 + #0.901 #0.062 TRPSOMG_02_01 + #0.901 &
          #0.062 TRPSOMG_03_01 + #0.901 #0.062 TRPSOMG_04_01 + #0.901 #0.062 &
          TRPSOMG_05_01 + #0.901 #0.062 TRPSOMG_06_01 + #0.901 #0.062 TRPSOMG_07_01 + &
          #0.901 #0.062 TRPSOMG_08_01 + #0.901 #0.000 TRPSOMG_09_01 + #0.901 &
          #0.062 TRPSOMG_01_02 + #0.901 #0.062 TRPSOMG_02_02 + #0.901 #0.062 &
          TRPSOMG_03_02 + #0.901 #0.062 TRPSOMG_04_02 + #0.901 #0.062 TRPSOMG_05_02 + &
          #0.901 #0.062 TRPSOMG_06_02 + #0.901 #0.062 TRPSOMG_07_02 + #0.901 &
          #0.062 TRPSOMG_08_02 + #0.901 #0.000 TRPSOMG_09_02 + #0.901 #0.062 &
          TRPSOMG_02_03 + #0.901 #0.062 TRPSOMG_03_03 + #0.901 #0.062 TRPSOMG_04_03 + &
          #0.901 #0.062 TRPSOMG_05_03 + #0.901 #0.062 TRPSOMG_06_03 + #0.901 &
          #0.062 TRPSOMG_07_03 + #0.901 #0.062 TRPSOMG_08_03 + #0.901 #0.062 &
          TRPSOMG_09_03 + #0.901 #0.062 TRPSOMG_02_04 + #0.901 #0.062 TRPSOMG_03_04 + &
          #0.901 #0.062 TRPSOMG_04_04 + #0.901 #0.062 TRPSOMG_05_04 + #0.901 &
          #0.062 TRPSOMG_06_04 + #0.901 #0.062 TRPSOMG_07_04 + #0.901 #0.062 &
          TRPSOMG_08_04 + #0.901 #0.062 TRPSOMG_09_04
 S215.5 TRPSOMG_10_04 + OH = #0.077 #0.361 TRPSOMG_10_05 + #0.077 #0.007 &
          TRPSOMG_10_06 + #0.077 #0.632 TRPSOMG_10_07 + #0.923 #0.053 TRPSOMG_01_01 + &
          #0.923 #0.053 TRPSOMG_02_01 + #0.923 #0.053 TRPSOMG_03_01 + #0.923 &
          #0.053 TRPSOMG_04_01 + #0.923 #0.053 TRPSOMG_05_01 + #0.923 #0.053 &
          TRPSOMG_06_01 + #0.923 #0.053 TRPSOMG_07_01 + #0.923 #0.000 TRPSOMG_08_01 + &
          #0.923 #0.000 TRPSOMG_09_01 + #0.923 #0.053 TRPSOMG_01_02 + #0.923 &
          #0.053 TRPSOMG_02_02 + #0.923 #0.053 TRPSOMG_03_02 + #0.923 #0.053 &
          TRPSOMG_04_02 + #0.923 #0.053 TRPSOMG_05_02 + #0.923 #0.053 TRPSOMG_06_02 + &
          #0.923 #0.053 TRPSOMG_07_02 + #0.923 #0.053 TRPSOMG_08_02 + #0.923 &
          #0.000 TRPSOMG_09_02 + #0.923 #0.053 TRPSOMG_02_03 + #0.923 #0.053 &
          TRPSOMG_03_03 + #0.923 #0.053 TRPSOMG_04_03 + #0.923 #0.053 TRPSOMG_05_03 + &
          #0.923 #0.053 TRPSOMG_06_03 + #0.923 #0.053 TRPSOMG_07_03 + #0.923 &
          #0.053 TRPSOMG_08_03 + #0.923 #0.000 TRPSOMG_09_03 + #0.923 #0.053 &
          TRPSOMG_02_04 + #0.923 #0.053 TRPSOMG_03_04 + #0.923 #0.053 TRPSOMG_04_04 + &
          #0.923 #0.053 TRPSOMG_05_04 + #0.923 #0.053 TRPSOMG_06_04 + #0.923 &
          #0.053 TRPSOMG_07_04 + #0.923 #0.053 TRPSOMG_08_04 + #0.923 #0.053 &
          TRPSOMG_09_04 + #0.923 #0.053 TRPSOMG_03_05 + #0.923 #0.053 TRPSOMG_04_05 + &
          #0.923 #0.053 TRPSOMG_05_05 + #0.923 #0.053 TRPSOMG_06_05 + #0.923 &
          #0.053 TRPSOMG_07_05 + #0.923 #0.053 TRPSOMG_08_05 + #0.923 #0.053 &
          TRPSOMG_09_05
 S216.5 TRPSOMG_10_05 + OH = #0.059 #0.361 TRPSOMG_10_06 + #0.059 #0.639 &
          TRPSOMG_10_07 + #0.941 #0.047 TRPSOMG_01_01 + #0.941 #0.047 TRPSOMG_02_01 + &
          #0.941 #0.047 TRPSOMG_03_01 + #0.941 #0.047 TRPSOMG_04_01 + #0.941 &
          #0.047 TRPSOMG_05_01 + #0.941 #0.047 TRPSOMG_06_01 + #0.941 #0.047 &
          TRPSOMG_07_01 + #0.941 #0.000 TRPSOMG_08_01 + #0.941 #0.000 TRPSOMG_09_01 + &
          #0.941 #0.047 TRPSOMG_01_02 + #0.941 #0.047 TRPSOMG_02_02 + #0.941 &
          #0.047 TRPSOMG_03_02 + #0.941 #0.047 TRPSOMG_04_02 + #0.941 #0.047 &
          TRPSOMG_05_02 + #0.941 #0.047 TRPSOMG_06_02 + #0.941 #0.047 TRPSOMG_07_02 + &
          #0.941 #0.000 TRPSOMG_08_02 + #0.941 #0.000 TRPSOMG_09_02 + #0.941 &
          #0.047 TRPSOMG_02_03 + #0.941 #0.047 TRPSOMG_03_03 + #0.941 #0.047 &
          TRPSOMG_04_03 + #0.941 #0.047 TRPSOMG_05_03 + #0.941 #0.047 TRPSOMG_06_03 + &
          #0.941 #0.047 TRPSOMG_07_03 + #0.941 #0.047 TRPSOMG_08_03 + #0.941 &
          #0.000 TRPSOMG_09_03 + #0.941 #0.047 TRPSOMG_02_04 + #0.941 #0.047 &
          TRPSOMG_03_04 + #0.941 #0.047 TRPSOMG_04_04 + #0.941 #0.047 TRPSOMG_05_04 + &
          #0.941 #0.047 TRPSOMG_06_04 + #0.941 #0.047 TRPSOMG_07_04 + #0.941 &
          #0.047 TRPSOMG_08_04 + #0.941 #0.000 TRPSOMG_09_04 + #0.941 #0.047 &
          TRPSOMG_03_05 + #0.941 #0.047 TRPSOMG_04_05 + #0.941 #0.047 TRPSOMG_05_05 + &
          #0.941 #0.047 TRPSOMG_06_05 + #0.941 #0.047 TRPSOMG_07_05 + #0.941 &
          #0.047 TRPSOMG_08_05 + #0.941 #0.047 TRPSOMG_09_05 + #0.941 #0.047 &
          TRPSOMG_03_06 + #0.941 #0.047 TRPSOMG_04_06 + #0.941 #0.047 TRPSOMG_05_06 + &
          #0.941 #0.047 TRPSOMG_06_06 + #0.941 #0.047 TRPSOMG_07_06 + #0.941 &
          #0.047 TRPSOMG_08_06 + #0.941 #0.047 TRPSOMG_09_06
 S217.5 TRPSOMG_10_06 + OH = #0.043 #1.000 TRPSOMG_10_07 + #0.957 #0.043 &
          TRPSOMG_01_01 + #0.957 #0.043 TRPSOMG_02_01 + #0.957 #0.043 TRPSOMG_03_01 + &
          #0.957 #0.043 TRPSOMG_04_01 + #0.957 #0.043 TRPSOMG_05_01 + #0.957 &
          #0.043 TRPSOMG_06_01 + #0.957 #0.000 TRPSOMG_07_01 + #0.957 #0.000 &
          TRPSOMG_08_01 + #0.957 #0.000 TRPSOMG_09_01 + #0.957 #0.043 TRPSOMG_01_02 + &
          #0.957 #0.043 TRPSOMG_02_02 + #0.957 #0.043 TRPSOMG_03_02 + #0.957 &
          #0.043 TRPSOMG_04_02 + #0.957 #0.043 TRPSOMG_05_02 + #0.957 #0.043 &
          TRPSOMG_06_02 + #0.957 #0.043 TRPSOMG_07_02 + #0.957 #0.000 TRPSOMG_08_02 + &
          #0.957 #0.000 TRPSOMG_09_02 + #0.957 #0.043 TRPSOMG_02_03 + #0.957 &
          #0.043 TRPSOMG_03_03 + #0.957 #0.043 TRPSOMG_04_03 + #0.957 #0.043 &
          TRPSOMG_05_03 + #0.957 #0.043 TRPSOMG_06_03 + #0.957 #0.043 TRPSOMG_07_03 + &
          #0.957 #0.000 TRPSOMG_08_03 + #0.957 #0.000 TRPSOMG_09_03 + #0.957 &
          #0.043 TRPSOMG_02_04 + #0.957 #0.043 TRPSOMG_03_04 + #0.957 #0.043 &
          TRPSOMG_04_04 + #0.957 #0.043 TRPSOMG_05_04 + #0.957 #0.043 TRPSOMG_06_04 + &
          #0.957 #0.043 TRPSOMG_07_04 + #0.957 #0.043 TRPSOMG_08_04 + #0.957 &
          #0.000 TRPSOMG_09_04 + #0.957 #0.043 TRPSOMG_03_05 + #0.957 #0.043 &
          TRPSOMG_04_05 + #0.957 #0.043 TRPSOMG_05_05 + #0.957 #0.043 TRPSOMG_06_05 + &
          #0.957 #0.043 TRPSOMG_07_05 + #0.957 #0.043 TRPSOMG_08_05 + #0.957 &
          #0.000 TRPSOMG_09_05 + #0.957 #0.043 TRPSOMG_03_06 + #0.957 #0.043 &
          TRPSOMG_04_06 + #0.957 #0.043 TRPSOMG_05_06 + #0.957 #0.043 TRPSOMG_06_06 + &
          #0.957 #0.043 TRPSOMG_07_06 + #0.957 #0.043 TRPSOMG_08_06 + #0.957 &
          #0.043 TRPSOMG_09_06 + #0.957 #0.043 TRPSOMG_04_07 + #0.957 #0.043 &
          TRPSOMG_05_07 + #0.957 #0.043 TRPSOMG_06_07 + #0.957 #0.043 TRPSOMG_07_07 + &
          #0.957 #0.043 TRPSOMG_08_07 + #0.957 #0.043 TRPSOMG_09_07
 S218.5 TRPSOMG_10_07 + OH = #1.000 #0.045 TRPSOMG_01_01 + #1.000 #0.045 &
          TRPSOMG_02_01 + #1.000 #0.045 TRPSOMG_03_01 + #1.000 #0.045 TRPSOMG_04_01 + &
          #1.000 #0.045 TRPSOMG_05_01 + #1.000 #0.045 TRPSOMG_06_01 + #1.000 &
          #0.000 TRPSOMG_07_01 + #1.000 #0.000 TRPSOMG_08_01 + #1.000 #0.000 &
          TRPSOMG_09_01 + #1.000 #0.045 TRPSOMG_01_02 + #1.000 #0.045 TRPSOMG_02_02 + &
          #1.000 #0.045 TRPSOMG_03_02 + #1.000 #0.045 TRPSOMG_04_02 + #1.000 &
          #0.045 TRPSOMG_05_02 + #1.000 #0.045 TRPSOMG_06_02 + #1.000 #0.000 &
          TRPSOMG_07_02 + #1.000 #0.000 TRPSOMG_08_02 + #1.000 #0.000 TRPSOMG_09_02 + &
          #1.000 #0.045 TRPSOMG_02_03 + #1.000 #0.045 TRPSOMG_03_03 + #1.000 &
          #0.045 TRPSOMG_04_03 + #1.000 #0.045 TRPSOMG_05_03 + #1.000 #0.045 &
          TRPSOMG_06_03 + #1.000 #0.045 TRPSOMG_07_03 + #1.000 #0.000 TRPSOMG_08_03 + &
          #1.000 #0.000 TRPSOMG_09_03 + #1.000 #0.045 TRPSOMG_02_04 + #1.000 &
          #0.045 TRPSOMG_03_04 + #1.000 #0.045 TRPSOMG_04_04 + #1.000 #0.045 &
          TRPSOMG_05_04 + #1.000 #0.045 TRPSOMG_06_04 + #1.000 #0.045 TRPSOMG_07_04 + &
          #1.000 #0.000 TRPSOMG_08_04 + #1.000 #0.000 TRPSOMG_09_04 + #1.000 &
          #0.045 TRPSOMG_03_05 + #1.000 #0.045 TRPSOMG_04_05 + #1.000 #0.045 &
          TRPSOMG_05_05 + #1.000 #0.045 TRPSOMG_06_05 + #1.000 #0.045 TRPSOMG_07_05 + &
          #1.000 #0.045 TRPSOMG_08_05 + #1.000 #0.000 TRPSOMG_09_05 + #1.000 &
          #0.045 TRPSOMG_03_06 + #1.000 #0.045 TRPSOMG_04_06 + #1.000 #0.045 &
          TRPSOMG_05_06 + #1.000 #0.045 TRPSOMG_06_06 + #1.000 #0.045 TRPSOMG_07_06 + &
          #1.000 #0.045 TRPSOMG_08_06 + #1.000 #0.000 TRPSOMG_09_06 + #1.000 &
          #0.045 TRPSOMG_04_07 + #1.000 #0.045 TRPSOMG_05_07 + #1.000 #0.045 &
          TRPSOMG_06_07 + #1.000 #0.045 TRPSOMG_07_07 + #1.000 #0.045 TRPSOMG_08_07 + &
          #1.000 #0.045 TRPSOMG_09_07
 S219.5 TRPSOMG_11_01 + OH = #0.188 #0.361 TRPSOMG_11_02 + #0.188 #0.007 &
          TRPSOMG_11_03 + #0.188 #0.606 TRPSOMG_11_04 + #0.188 #0.026 TRPSOMG_11_05 + &
          #0.812 #0.100 TRPSOMG_01_01 + #0.812 #0.100 TRPSOMG_02_01 + #0.812 &
          #0.100 TRPSOMG_03_01 + #0.812 #0.100 TRPSOMG_04_01 + #0.812 #0.100 &
          TRPSOMG_05_01 + #0.812 #0.100 TRPSOMG_06_01 + #0.812 #0.100 TRPSOMG_07_01 + &
          #0.812 #0.100 TRPSOMG_08_01 + #0.812 #0.100 TRPSOMG_09_01 + #0.812 &
          #0.100 TRPSOMG_10_01 + #0.812 #0.100 TRPSOMG_01_02 + #0.812 #0.100 &
          TRPSOMG_02_02 + #0.812 #0.100 TRPSOMG_03_02 + #0.812 #0.100 TRPSOMG_04_02 + &
          #0.812 #0.100 TRPSOMG_05_02 + #0.812 #0.100 TRPSOMG_06_02 + #0.812 &
          #0.100 TRPSOMG_07_02 + #0.812 #0.100 TRPSOMG_08_02 + #0.812 #0.100 &
          TRPSOMG_09_02 + #0.812 #0.100 TRPSOMG_10_02
 S220.5 TRPSOMG_11_02 + OH = #0.138 #0.361 TRPSOMG_11_03 + #0.138 #0.007 &
          TRPSOMG_11_04 + #0.138 #0.606 TRPSOMG_11_05 + #0.138 #0.026 TRPSOMG_11_06 + &
          #0.862 #0.071 TRPSOMG_01_01 + #0.862 #0.071 TRPSOMG_02_01 + #0.862 &
          #0.071 TRPSOMG_03_01 + #0.862 #0.071 TRPSOMG_04_01 + #0.862 #0.071 &
          TRPSOMG_05_01 + #0.862 #0.071 TRPSOMG_06_01 + #0.862 #0.071 TRPSOMG_07_01 + &
          #0.862 #0.071 TRPSOMG_08_01 + #0.862 #0.071 TRPSOMG_09_01 + #0.862 &
          #0.000 TRPSOMG_10_01 + #0.862 #0.071 TRPSOMG_01_02 + #0.862 #0.071 &
          TRPSOMG_02_02 + #0.862 #0.071 TRPSOMG_03_02 + #0.862 #0.071 TRPSOMG_04_02 + &
          #0.862 #0.071 TRPSOMG_05_02 + #0.862 #0.071 TRPSOMG_06_02 + #0.862 &
          #0.071 TRPSOMG_07_02 + #0.862 #0.071 TRPSOMG_08_02 + #0.862 #0.071 &
          TRPSOMG_09_02 + #0.862 #0.071 TRPSOMG_10_02 + #0.862 #0.071 TRPSOMG_02_03 + &
          #0.862 #0.071 TRPSOMG_03_03 + #0.862 #0.071 TRPSOMG_04_03 + #0.862 &
          #0.071 TRPSOMG_05_03 + #0.862 #0.071 TRPSOMG_06_03 + #0.862 #0.071 &
          TRPSOMG_07_03 + #0.862 #0.071 TRPSOMG_08_03 + #0.862 #0.071 TRPSOMG_09_03 + &
          #0.862 #0.071 TRPSOMG_10_03
 S221.5 TRPSOMG_11_03 + OH = #0.107 #0.361 TRPSOMG_11_04 + #0.107 #0.007 &
          TRPSOMG_11_05 + #0.107 #0.606 TRPSOMG_11_06 + #0.107 #0.026 TRPSOMG_11_07 + &
          #0.893 #0.056 TRPSOMG_01_01 + #0.893 #0.056 TRPSOMG_02_01 + #0.893 &
          #0.056 TRPSOMG_03_01 + #0.893 #0.056 TRPSOMG_04_01 + #0.893 #0.056 &
          TRPSOMG_05_01 + #0.893 #0.056 TRPSOMG_06_01 + #0.893 #0.056 TRPSOMG_07_01 + &
          #0.893 #0.056 TRPSOMG_08_01 + #0.893 #0.056 TRPSOMG_09_01 + #0.893 &
          #0.000 TRPSOMG_10_01 + #0.893 #0.056 TRPSOMG_01_02 + #0.893 #0.056 &
          TRPSOMG_02_02 + #0.893 #0.056 TRPSOMG_03_02 + #0.893 #0.056 TRPSOMG_04_02 + &
          #0.893 #0.056 TRPSOMG_05_02 + #0.893 #0.056 TRPSOMG_06_02 + #0.893 &
          #0.056 TRPSOMG_07_02 + #0.893 #0.056 TRPSOMG_08_02 + #0.893 #0.056 &
          TRPSOMG_09_02 + #0.893 #0.000 TRPSOMG_10_02 + #0.893 #0.056 TRPSOMG_02_03 + &
          #0.893 #0.056 TRPSOMG_03_03 + #0.893 #0.056 TRPSOMG_04_03 + #0.893 &
          #0.056 TRPSOMG_05_03 + #0.893 #0.056 TRPSOMG_06_03 + #0.893 #0.056 &
          TRPSOMG_07_03 + #0.893 #0.056 TRPSOMG_08_03 + #0.893 #0.056 TRPSOMG_09_03 + &
          #0.893 #0.056 TRPSOMG_10_03 + #0.893 #0.056 TRPSOMG_02_04 + #0.893 &
          #0.056 TRPSOMG_03_04 + #0.893 #0.056 TRPSOMG_04_04 + #0.893 #0.056 &
          TRPSOMG_05_04 + #0.893 #0.056 TRPSOMG_06_04 + #0.893 #0.056 TRPSOMG_07_04 + &
          #0.893 #0.056 TRPSOMG_08_04 + #0.893 #0.056 TRPSOMG_09_04 + #0.893 &
          #0.056 TRPSOMG_10_04
 S222.5 TRPSOMG_11_04 + OH = #0.084 #0.361 TRPSOMG_11_05 + #0.084 #0.007 &
          TRPSOMG_11_06 + #0.084 #0.632 TRPSOMG_11_07 + #0.916 #0.048 TRPSOMG_01_01 + &
          #0.916 #0.048 TRPSOMG_02_01 + #0.916 #0.048 TRPSOMG_03_01 + #0.916 &
          #0.048 TRPSOMG_04_01 + #0.916 #0.048 TRPSOMG_05_01 + #0.916 #0.048 &
          TRPSOMG_06_01 + #0.916 #0.048 TRPSOMG_07_01 + #0.916 #0.048 TRPSOMG_08_01 + &
          #0.916 #0.000 TRPSOMG_09_01 + #0.916 #0.000 TRPSOMG_10_01 + #0.916 &
          #0.048 TRPSOMG_01_02 + #0.916 #0.048 TRPSOMG_02_02 + #0.916 #0.048 &
          TRPSOMG_03_02 + #0.916 #0.048 TRPSOMG_04_02 + #0.916 #0.048 TRPSOMG_05_02 + &
          #0.916 #0.048 TRPSOMG_06_02 + #0.916 #0.048 TRPSOMG_07_02 + #0.916 &
          #0.048 TRPSOMG_08_02 + #0.916 #0.048 TRPSOMG_09_02 + #0.916 #0.000 &
          TRPSOMG_10_02 + #0.916 #0.048 TRPSOMG_02_03 + #0.916 #0.048 TRPSOMG_03_03 + &
          #0.916 #0.048 TRPSOMG_04_03 + #0.916 #0.048 TRPSOMG_05_03 + #0.916 &
          #0.048 TRPSOMG_06_03 + #0.916 #0.048 TRPSOMG_07_03 + #0.916 #0.048 &
          TRPSOMG_08_03 + #0.916 #0.048 TRPSOMG_09_03 + #0.916 #0.000 TRPSOMG_10_03 + &
          #0.916 #0.048 TRPSOMG_02_04 + #0.916 #0.048 TRPSOMG_03_04 + #0.916 &
          #0.048 TRPSOMG_04_04 + #0.916 #0.048 TRPSOMG_05_04 + #0.916 #0.048 &
          TRPSOMG_06_04 + #0.916 #0.048 TRPSOMG_07_04 + #0.916 #0.048 TRPSOMG_08_04 + &
          #0.916 #0.048 TRPSOMG_09_04 + #0.916 #0.048 TRPSOMG_10_04 + #0.916 &
          #0.048 TRPSOMG_03_05 + #0.916 #0.048 TRPSOMG_04_05 + #0.916 #0.048 &
          TRPSOMG_05_05 + #0.916 #0.048 TRPSOMG_06_05 + #0.916 #0.048 TRPSOMG_07_05 + &
          #0.916 #0.048 TRPSOMG_08_05 + #0.916 #0.048 TRPSOMG_09_05 + #0.916 &
          #0.048 TRPSOMG_10_05
 S223.5 TRPSOMG_11_05 + OH = #0.066 #0.361 TRPSOMG_11_06 + #0.066 #0.639 &
          TRPSOMG_11_07 + #0.934 #0.042 TRPSOMG_01_01 + #0.934 #0.042 TRPSOMG_02_01 + &
          #0.934 #0.042 TRPSOMG_03_01 + #0.934 #0.042 TRPSOMG_04_01 + #0.934 &
          #0.042 TRPSOMG_05_01 + #0.934 #0.042 TRPSOMG_06_01 + #0.934 #0.042 &
          TRPSOMG_07_01 + #0.934 #0.042 TRPSOMG_08_01 + #0.934 #0.000 TRPSOMG_09_01 + &
          #0.934 #0.000 TRPSOMG_10_01 + #0.934 #0.042 TRPSOMG_01_02 + #0.934 &
          #0.042 TRPSOMG_02_02 + #0.934 #0.042 TRPSOMG_03_02 + #0.934 #0.042 &
          TRPSOMG_04_02 + #0.934 #0.042 TRPSOMG_05_02 + #0.934 #0.042 TRPSOMG_06_02 + &
          #0.934 #0.042 TRPSOMG_07_02 + #0.934 #0.042 TRPSOMG_08_02 + #0.934 &
          #0.000 TRPSOMG_09_02 + #0.934 #0.000 TRPSOMG_10_02 + #0.934 #0.042 &
          TRPSOMG_02_03 + #0.934 #0.042 TRPSOMG_03_03 + #0.934 #0.042 TRPSOMG_04_03 + &
          #0.934 #0.042 TRPSOMG_05_03 + #0.934 #0.042 TRPSOMG_06_03 + #0.934 &
          #0.042 TRPSOMG_07_03 + #0.934 #0.042 TRPSOMG_08_03 + #0.934 #0.042 &
          TRPSOMG_09_03 + #0.934 #0.000 TRPSOMG_10_03 + #0.934 #0.042 TRPSOMG_02_04 + &
          #0.934 #0.042 TRPSOMG_03_04 + #0.934 #0.042 TRPSOMG_04_04 + #0.934 &
          #0.042 TRPSOMG_05_04 + #0.934 #0.042 TRPSOMG_06_04 + #0.934 #0.042 &
          TRPSOMG_07_04 + #0.934 #0.042 TRPSOMG_08_04 + #0.934 #0.042 TRPSOMG_09_04 + &
          #0.934 #0.000 TRPSOMG_10_04 + #0.934 #0.042 TRPSOMG_03_05 + #0.934 &
          #0.042 TRPSOMG_04_05 + #0.934 #0.042 TRPSOMG_05_05 + #0.934 #0.042 &
          TRPSOMG_06_05 + #0.934 #0.042 TRPSOMG_07_05 + #0.934 #0.042 TRPSOMG_08_05 + &
          #0.934 #0.042 TRPSOMG_09_05 + #0.934 #0.042 TRPSOMG_10_05 + #0.934 &
          #0.042 TRPSOMG_03_06 + #0.934 #0.042 TRPSOMG_04_06 + #0.934 #0.042 &
          TRPSOMG_05_06 + #0.934 #0.042 TRPSOMG_06_06 + #0.934 #0.042 TRPSOMG_07_06 + &
          #0.934 #0.042 TRPSOMG_08_06 + #0.934 #0.042 TRPSOMG_09_06 + #0.934 &
          #0.042 TRPSOMG_10_06
 S224.5 TRPSOMG_11_06 + OH = #0.051 #1.000 TRPSOMG_11_07 + #0.949 #0.038 &
          TRPSOMG_01_01 + #0.949 #0.038 TRPSOMG_02_01 + #0.949 #0.038 TRPSOMG_03_01 + &
          #0.949 #0.038 TRPSOMG_04_01 + #0.949 #0.038 TRPSOMG_05_01 + #0.949 &
          #0.038 TRPSOMG_06_01 + #0.949 #0.038 TRPSOMG_07_01 + #0.949 #0.000 &
          TRPSOMG_08_01 + #0.949 #0.000 TRPSOMG_09_01 + #0.949 #0.000 TRPSOMG_10_01 + &
          #0.949 #0.038 TRPSOMG_01_02 + #0.949 #0.038 TRPSOMG_02_02 + #0.949 &
          #0.038 TRPSOMG_03_02 + #0.949 #0.038 TRPSOMG_04_02 + #0.949 #0.038 &
          TRPSOMG_05_02 + #0.949 #0.038 TRPSOMG_06_02 + #0.949 #0.038 TRPSOMG_07_02 + &
          #0.949 #0.038 TRPSOMG_08_02 + #0.949 #0.000 TRPSOMG_09_02 + #0.949 &
          #0.000 TRPSOMG_10_02 + #0.949 #0.038 TRPSOMG_02_03 + #0.949 #0.038 &
          TRPSOMG_03_03 + #0.949 #0.038 TRPSOMG_04_03 + #0.949 #0.038 TRPSOMG_05_03 + &
          #0.949 #0.038 TRPSOMG_06_03 + #0.949 #0.038 TRPSOMG_07_03 + #0.949 &
          #0.038 TRPSOMG_08_03 + #0.949 #0.000 TRPSOMG_09_03 + #0.949 #0.000 &
          TRPSOMG_10_03 + #0.949 #0.038 TRPSOMG_02_04 + #0.949 #0.038 TRPSOMG_03_04 + &
          #0.949 #0.038 TRPSOMG_04_04 + #0.949 #0.038 TRPSOMG_05_04 + #0.949 &
          #0.038 TRPSOMG_06_04 + #0.949 #0.038 TRPSOMG_07_04 + #0.949 #0.038 &
          TRPSOMG_08_04 + #0.949 #0.038 TRPSOMG_09_04 + #0.949 #0.000 TRPSOMG_10_04 + &
          #0.949 #0.038 TRPSOMG_03_05 + #0.949 #0.038 TRPSOMG_04_05 + #0.949 &
          #0.038 TRPSOMG_05_05 + #0.949 #0.038 TRPSOMG_06_05 + #0.949 #0.038 &
          TRPSOMG_07_05 + #0.949 #0.038 TRPSOMG_08_05 + #0.949 #0.038 TRPSOMG_09_05 + &
          #0.949 #0.000 TRPSOMG_10_05 + #0.949 #0.038 TRPSOMG_03_06 + #0.949 &
          #0.038 TRPSOMG_04_06 + #0.949 #0.038 TRPSOMG_05_06 + #0.949 #0.038 &
          TRPSOMG_06_06 + #0.949 #0.038 TRPSOMG_07_06 + #0.949 #0.038 TRPSOMG_08_06 + &
          #0.949 #0.038 TRPSOMG_09_06 + #0.949 #0.038 TRPSOMG_10_06 + #0.949 &
          #0.038 TRPSOMG_04_07 + #0.949 #0.038 TRPSOMG_05_07 + #0.949 #0.038 &
          TRPSOMG_06_07 + #0.949 #0.038 TRPSOMG_07_07 + #0.949 #0.038 TRPSOMG_08_07 + &
          #0.949 #0.038 TRPSOMG_09_07 + #0.949 #0.038 TRPSOMG_10_07
 S225.5 TRPSOMG_11_07 + OH = #1.000 #0.040 TRPSOMG_01_01 + #1.000 #0.040 &
          TRPSOMG_02_01 + #1.000 #0.040 TRPSOMG_03_01 + #1.000 #0.040 TRPSOMG_04_01 + &
          #1.000 #0.040 TRPSOMG_05_01 + #1.000 #0.040 TRPSOMG_06_01 + #1.000 &
          #0.040 TRPSOMG_07_01 + #1.000 #0.000 TRPSOMG_08_01 + #1.000 #0.000 &
          TRPSOMG_09_01 + #1.000 #0.000 TRPSOMG_10_01 + #1.000 #0.040 TRPSOMG_01_02 + &
          #1.000 #0.040 TRPSOMG_02_02 + #1.000 #0.040 TRPSOMG_03_02 + #1.000 &
          #0.040 TRPSOMG_04_02 + #1.000 #0.040 TRPSOMG_05_02 + #1.000 #0.040 &
          TRPSOMG_06_02 + #1.000 #0.040 TRPSOMG_07_02 + #1.000 #0.000 TRPSOMG_08_02 + &
          #1.000 #0.000 TRPSOMG_09_02 + #1.000 #0.000 TRPSOMG_10_02 + #1.000 &
          #0.040 TRPSOMG_02_03 + #1.000 #0.040 TRPSOMG_03_03 + #1.000 #0.040 &
          TRPSOMG_04_03 + #1.000 #0.040 TRPSOMG_05_03 + #1.000 #0.040 TRPSOMG_06_03 + &
          #1.000 #0.040 TRPSOMG_07_03 + #1.000 #0.040 TRPSOMG_08_03 + #1.000 &
          #0.000 TRPSOMG_09_03 + #1.000 #0.000 TRPSOMG_10_03 + #1.000 #0.040 &
          TRPSOMG_02_04 + #1.000 #0.040 TRPSOMG_03_04 + #1.000 #0.040 TRPSOMG_04_04 + &
          #1.000 #0.040 TRPSOMG_05_04 + #1.000 #0.040 TRPSOMG_06_04 + #1.000 &
          #0.040 TRPSOMG_07_04 + #1.000 #0.040 TRPSOMG_08_04 + #1.000 #0.000 &
          TRPSOMG_09_04 + #1.000 #0.000 TRPSOMG_10_04 + #1.000 #0.040 TRPSOMG_03_05 + &
          #1.000 #0.040 TRPSOMG_04_05 + #1.000 #0.040 TRPSOMG_05_05 + #1.000 &
          #0.040 TRPSOMG_06_05 + #1.000 #0.040 TRPSOMG_07_05 + #1.000 #0.040 &
          TRPSOMG_08_05 + #1.000 #0.040 TRPSOMG_09_05 + #1.000 #0.000 TRPSOMG_10_05 + &
          #1.000 #0.040 TRPSOMG_03_06 + #1.000 #0.040 TRPSOMG_04_06 + #1.000 &
          #0.040 TRPSOMG_05_06 + #1.000 #0.040 TRPSOMG_06_06 + #1.000 #0.040 &
          TRPSOMG_07_06 + #1.000 #0.040 TRPSOMG_08_06 + #1.000 #0.040 TRPSOMG_09_06 + &
          #1.000 #0.000 TRPSOMG_10_06 + #1.000 #0.040 TRPSOMG_04_07 + #1.000 &
          #0.040 TRPSOMG_05_07 + #1.000 #0.040 TRPSOMG_06_07 + #1.000 #0.040 &
          TRPSOMG_07_07 + #1.000 #0.040 TRPSOMG_08_07 + #1.000 #0.040 TRPSOMG_09_07 + &
          #1.000 #0.040 TRPSOMG_10_07
 S226.5 TRPSOMG_12_01 + OH = #0.194 #0.361 TRPSOMG_12_02 + #0.194 #0.007 &
          TRPSOMG_12_03 + #0.194 #0.606 TRPSOMG_12_04 + #0.194 #0.026 TRPSOMG_12_05 + &
          #0.806 #0.091 TRPSOMG_01_01 + #0.806 #0.091 TRPSOMG_02_01 + #0.806 &
          #0.091 TRPSOMG_03_01 + #0.806 #0.091 TRPSOMG_04_01 + #0.806 #0.091 &
          TRPSOMG_05_01 + #0.806 #0.091 TRPSOMG_06_01 + #0.806 #0.091 TRPSOMG_07_01 + &
          #0.806 #0.091 TRPSOMG_08_01 + #0.806 #0.091 TRPSOMG_09_01 + #0.806 &
          #0.091 TRPSOMG_10_01 + #0.806 #0.091 TRPSOMG_11_01 + #0.806 #0.091 &
          TRPSOMG_01_02 + #0.806 #0.091 TRPSOMG_02_02 + #0.806 #0.091 TRPSOMG_03_02 + &
          #0.806 #0.091 TRPSOMG_04_02 + #0.806 #0.091 TRPSOMG_05_02 + #0.806 &
          #0.091 TRPSOMG_06_02 + #0.806 #0.091 TRPSOMG_07_02 + #0.806 #0.091 &
          TRPSOMG_08_02 + #0.806 #0.091 TRPSOMG_09_02 + #0.806 #0.091 TRPSOMG_10_02 + &
          #0.806 #0.091 TRPSOMG_11_02
 S227.5 TRPSOMG_12_02 + OH = #0.144 #0.361 TRPSOMG_12_03 + #0.144 #0.007 &
          TRPSOMG_12_04 + #0.144 #0.606 TRPSOMG_12_05 + #0.144 #0.026 TRPSOMG_12_06 + &
          #0.856 #0.065 TRPSOMG_01_01 + #0.856 #0.065 TRPSOMG_02_01 + #0.856 &
          #0.065 TRPSOMG_03_01 + #0.856 #0.065 TRPSOMG_04_01 + #0.856 #0.065 &
          TRPSOMG_05_01 + #0.856 #0.065 TRPSOMG_06_01 + #0.856 #0.065 TRPSOMG_07_01 + &
          #0.856 #0.065 TRPSOMG_08_01 + #0.856 #0.065 TRPSOMG_09_01 + #0.856 &
          #0.065 TRPSOMG_10_01 + #0.856 #0.000 TRPSOMG_11_01 + #0.856 #0.065 &
          TRPSOMG_01_02 + #0.856 #0.065 TRPSOMG_02_02 + #0.856 #0.065 TRPSOMG_03_02 + &
          #0.856 #0.065 TRPSOMG_04_02 + #0.856 #0.065 TRPSOMG_05_02 + #0.856 &
          #0.065 TRPSOMG_06_02 + #0.856 #0.065 TRPSOMG_07_02 + #0.856 #0.065 &
          TRPSOMG_08_02 + #0.856 #0.065 TRPSOMG_09_02 + #0.856 #0.065 TRPSOMG_10_02 + &
          #0.856 #0.065 TRPSOMG_11_02 + #0.856 #0.065 TRPSOMG_02_03 + #0.856 &
          #0.065 TRPSOMG_03_03 + #0.856 #0.065 TRPSOMG_04_03 + #0.856 #0.065 &
          TRPSOMG_05_03 + #0.856 #0.065 TRPSOMG_06_03 + #0.856 #0.065 TRPSOMG_07_03 + &
          #0.856 #0.065 TRPSOMG_08_03 + #0.856 #0.065 TRPSOMG_09_03 + #0.856 &
          #0.065 TRPSOMG_10_03 + #0.856 #0.065 TRPSOMG_11_03
 S228.5 TRPSOMG_12_03 + OH = #0.114 #0.361 TRPSOMG_12_04 + #0.114 #0.007 &
          TRPSOMG_12_05 + #0.114 #0.606 TRPSOMG_12_06 + #0.114 #0.026 TRPSOMG_12_07 + &
          #0.886 #0.050 TRPSOMG_01_01 + #0.886 #0.050 TRPSOMG_02_01 + #0.886 &
          #0.050 TRPSOMG_03_01 + #0.886 #0.050 TRPSOMG_04_01 + #0.886 #0.050 &
          TRPSOMG_05_01 + #0.886 #0.050 TRPSOMG_06_01 + #0.886 #0.050 TRPSOMG_07_01 + &
          #0.886 #0.050 TRPSOMG_08_01 + #0.886 #0.050 TRPSOMG_09_01 + #0.886 &
          #0.050 TRPSOMG_10_01 + #0.886 #0.000 TRPSOMG_11_01 + #0.886 #0.050 &
          TRPSOMG_01_02 + #0.886 #0.050 TRPSOMG_02_02 + #0.886 #0.050 TRPSOMG_03_02 + &
          #0.886 #0.050 TRPSOMG_04_02 + #0.886 #0.050 TRPSOMG_05_02 + #0.886 &
          #0.050 TRPSOMG_06_02 + #0.886 #0.050 TRPSOMG_07_02 + #0.886 #0.050 &
          TRPSOMG_08_02 + #0.886 #0.050 TRPSOMG_09_02 + #0.886 #0.050 TRPSOMG_10_02 + &
          #0.886 #0.000 TRPSOMG_11_02 + #0.886 #0.050 TRPSOMG_02_03 + #0.886 &
          #0.050 TRPSOMG_03_03 + #0.886 #0.050 TRPSOMG_04_03 + #0.886 #0.050 &
          TRPSOMG_05_03 + #0.886 #0.050 TRPSOMG_06_03 + #0.886 #0.050 TRPSOMG_07_03 + &
          #0.886 #0.050 TRPSOMG_08_03 + #0.886 #0.050 TRPSOMG_09_03 + #0.886 &
          #0.050 TRPSOMG_10_03 + #0.886 #0.050 TRPSOMG_11_03 + #0.886 #0.050 &
          TRPSOMG_02_04 + #0.886 #0.050 TRPSOMG_03_04 + #0.886 #0.050 TRPSOMG_04_04 + &
          #0.886 #0.050 TRPSOMG_05_04 + #0.886 #0.050 TRPSOMG_06_04 + #0.886 &
          #0.050 TRPSOMG_07_04 + #0.886 #0.050 TRPSOMG_08_04 + #0.886 #0.050 &
          TRPSOMG_09_04 + #0.886 #0.050 TRPSOMG_10_04 + #0.886 #0.050 TRPSOMG_11_04
 S229.5 TRPSOMG_12_04 + OH = #0.091 #0.361 TRPSOMG_12_05 + #0.091 #0.007 &
          TRPSOMG_12_06 + #0.091 #0.632 TRPSOMG_12_07 + #0.909 #0.043 TRPSOMG_01_01 + &
          #0.909 #0.043 TRPSOMG_02_01 + #0.909 #0.043 TRPSOMG_03_01 + #0.909 &
          #0.043 TRPSOMG_04_01 + #0.909 #0.043 TRPSOMG_05_01 + #0.909 #0.043 &
          TRPSOMG_06_01 + #0.909 #0.043 TRPSOMG_07_01 + #0.909 #0.043 TRPSOMG_08_01 + &
          #0.909 #0.043 TRPSOMG_09_01 + #0.909 #0.000 TRPSOMG_10_01 + #0.909 &
          #0.000 TRPSOMG_11_01 + #0.909 #0.043 TRPSOMG_01_02 + #0.909 #0.043 &
          TRPSOMG_02_02 + #0.909 #0.043 TRPSOMG_03_02 + #0.909 #0.043 TRPSOMG_04_02 + &
          #0.909 #0.043 TRPSOMG_05_02 + #0.909 #0.043 TRPSOMG_06_02 + #0.909 &
          #0.043 TRPSOMG_07_02 + #0.909 #0.043 TRPSOMG_08_02 + #0.909 #0.043 &
          TRPSOMG_09_02 + #0.909 #0.043 TRPSOMG_10_02 + #0.909 #0.000 TRPSOMG_11_02 + &
          #0.909 #0.043 TRPSOMG_02_03 + #0.909 #0.043 TRPSOMG_03_03 + #0.909 &
          #0.043 TRPSOMG_04_03 + #0.909 #0.043 TRPSOMG_05_03 + #0.909 #0.043 &
          TRPSOMG_06_03 + #0.909 #0.043 TRPSOMG_07_03 + #0.909 #0.043 TRPSOMG_08_03 + &
          #0.909 #0.043 TRPSOMG_09_03 + #0.909 #0.043 TRPSOMG_10_03 + #0.909 &
          #0.000 TRPSOMG_11_03 + #0.909 #0.043 TRPSOMG_02_04 + #0.909 #0.043 &
          TRPSOMG_03_04 + #0.909 #0.043 TRPSOMG_04_04 + #0.909 #0.043 TRPSOMG_05_04 + &
          #0.909 #0.043 TRPSOMG_06_04 + #0.909 #0.043 TRPSOMG_07_04 + #0.909 &
          #0.043 TRPSOMG_08_04 + #0.909 #0.043 TRPSOMG_09_04 + #0.909 #0.043 &
          TRPSOMG_10_04 + #0.909 #0.043 TRPSOMG_11_04 + #0.909 #0.043 TRPSOMG_03_05 + &
          #0.909 #0.043 TRPSOMG_04_05 + #0.909 #0.043 TRPSOMG_05_05 + #0.909 &
          #0.043 TRPSOMG_06_05 + #0.909 #0.043 TRPSOMG_07_05 + #0.909 #0.043 &
          TRPSOMG_08_05 + #0.909 #0.043 TRPSOMG_09_05 + #0.909 #0.043 TRPSOMG_10_05 + &
          #0.909 #0.043 TRPSOMG_11_05
 S230.5 TRPSOMG_12_05 + OH = #0.073 #0.361 TRPSOMG_12_06 + #0.073 #0.639 &
          TRPSOMG_12_07 + #0.927 #0.037 TRPSOMG_01_01 + #0.927 #0.037 TRPSOMG_02_01 + &
          #0.927 #0.037 TRPSOMG_03_01 + #0.927 #0.037 TRPSOMG_04_01 + #0.927 &
          #0.037 TRPSOMG_05_01 + #0.927 #0.037 TRPSOMG_06_01 + #0.927 #0.037 &
          TRPSOMG_07_01 + #0.927 #0.037 TRPSOMG_08_01 + #0.927 #0.037 TRPSOMG_09_01 + &
          #0.927 #0.000 TRPSOMG_10_01 + #0.927 #0.000 TRPSOMG_11_01 + #0.927 &
          #0.037 TRPSOMG_01_02 + #0.927 #0.037 TRPSOMG_02_02 + #0.927 #0.037 &
          TRPSOMG_03_02 + #0.927 #0.037 TRPSOMG_04_02 + #0.927 #0.037 TRPSOMG_05_02 + &
          #0.927 #0.037 TRPSOMG_06_02 + #0.927 #0.037 TRPSOMG_07_02 + #0.927 &
          #0.037 TRPSOMG_08_02 + #0.927 #0.037 TRPSOMG_09_02 + #0.927 #0.000 &
          TRPSOMG_10_02 + #0.927 #0.000 TRPSOMG_11_02 + #0.927 #0.037 TRPSOMG_02_03 + &
          #0.927 #0.037 TRPSOMG_03_03 + #0.927 #0.037 TRPSOMG_04_03 + #0.927 &
          #0.037 TRPSOMG_05_03 + #0.927 #0.037 TRPSOMG_06_03 + #0.927 #0.037 &
          TRPSOMG_07_03 + #0.927 #0.037 TRPSOMG_08_03 + #0.927 #0.037 TRPSOMG_09_03 + &
          #0.927 #0.037 TRPSOMG_10_03 + #0.927 #0.000 TRPSOMG_11_03 + #0.927 &
          #0.037 TRPSOMG_02_04 + #0.927 #0.037 TRPSOMG_03_04 + #0.927 #0.037 &
          TRPSOMG_04_04 + #0.927 #0.037 TRPSOMG_05_04 + #0.927 #0.037 TRPSOMG_06_04 + &
          #0.927 #0.037 TRPSOMG_07_04 + #0.927 #0.037 TRPSOMG_08_04 + #0.927 &
          #0.037 TRPSOMG_09_04 + #0.927 #0.037 TRPSOMG_10_04 + #0.927 #0.000 &
          TRPSOMG_11_04 + #0.927 #0.037 TRPSOMG_03_05 + #0.927 #0.037 TRPSOMG_04_05 + &
          #0.927 #0.037 TRPSOMG_05_05 + #0.927 #0.037 TRPSOMG_06_05 + #0.927 &
          #0.037 TRPSOMG_07_05 + #0.927 #0.037 TRPSOMG_08_05 + #0.927 #0.037 &
          TRPSOMG_09_05 + #0.927 #0.037 TRPSOMG_10_05 + #0.927 #0.037 TRPSOMG_11_05 + &
          #0.927 #0.037 TRPSOMG_03_06 + #0.927 #0.037 TRPSOMG_04_06 + #0.927 &
          #0.037 TRPSOMG_05_06 + #0.927 #0.037 TRPSOMG_06_06 + #0.927 #0.037 &
          TRPSOMG_07_06 + #0.927 #0.037 TRPSOMG_08_06 + #0.927 #0.037 TRPSOMG_09_06 + &
          #0.927 #0.037 TRPSOMG_10_06 + #0.927 #0.037 TRPSOMG_11_06
 S231.5 TRPSOMG_12_06 + OH = #0.059 #1.000 TRPSOMG_12_07 + #0.941 #0.034 &
          TRPSOMG_01_01 + #0.941 #0.034 TRPSOMG_02_01 + #0.941 #0.034 TRPSOMG_03_01 + &
          #0.941 #0.034 TRPSOMG_04_01 + #0.941 #0.034 TRPSOMG_05_01 + #0.941 &
          #0.034 TRPSOMG_06_01 + #0.941 #0.034 TRPSOMG_07_01 + #0.941 #0.034 &
          TRPSOMG_08_01 + #0.941 #0.000 TRPSOMG_09_01 + #0.941 #0.000 TRPSOMG_10_01 + &
          #0.941 #0.000 TRPSOMG_11_01 + #0.941 #0.034 TRPSOMG_01_02 + #0.941 &
          #0.034 TRPSOMG_02_02 + #0.941 #0.034 TRPSOMG_03_02 + #0.941 #0.034 &
          TRPSOMG_04_02 + #0.941 #0.034 TRPSOMG_05_02 + #0.941 #0.034 TRPSOMG_06_02 + &
          #0.941 #0.034 TRPSOMG_07_02 + #0.941 #0.034 TRPSOMG_08_02 + #0.941 &
          #0.034 TRPSOMG_09_02 + #0.941 #0.000 TRPSOMG_10_02 + #0.941 #0.000 &
          TRPSOMG_11_02 + #0.941 #0.034 TRPSOMG_02_03 + #0.941 #0.034 TRPSOMG_03_03 + &
          #0.941 #0.034 TRPSOMG_04_03 + #0.941 #0.034 TRPSOMG_05_03 + #0.941 &
          #0.034 TRPSOMG_06_03 + #0.941 #0.034 TRPSOMG_07_03 + #0.941 #0.034 &
          TRPSOMG_08_03 + #0.941 #0.034 TRPSOMG_09_03 + #0.941 #0.000 TRPSOMG_10_03 + &
          #0.941 #0.000 TRPSOMG_11_03 + #0.941 #0.034 TRPSOMG_02_04 + #0.941 &
          #0.034 TRPSOMG_03_04 + #0.941 #0.034 TRPSOMG_04_04 + #0.941 #0.034 &
          TRPSOMG_05_04 + #0.941 #0.034 TRPSOMG_06_04 + #0.941 #0.034 TRPSOMG_07_04 + &
          #0.941 #0.034 TRPSOMG_08_04 + #0.941 #0.034 TRPSOMG_09_04 + #0.941 &
          #0.034 TRPSOMG_10_04 + #0.941 #0.000 TRPSOMG_11_04 + #0.941 #0.034 &
          TRPSOMG_03_05 + #0.941 #0.034 TRPSOMG_04_05 + #0.941 #0.034 TRPSOMG_05_05 + &
          #0.941 #0.034 TRPSOMG_06_05 + #0.941 #0.034 TRPSOMG_07_05 + #0.941 &
          #0.034 TRPSOMG_08_05 + #0.941 #0.034 TRPSOMG_09_05 + #0.941 #0.034 &
          TRPSOMG_10_05 + #0.941 #0.000 TRPSOMG_11_05 + #0.941 #0.034 TRPSOMG_03_06 + &
          #0.941 #0.034 TRPSOMG_04_06 + #0.941 #0.034 TRPSOMG_05_06 + #0.941 &
          #0.034 TRPSOMG_06_06 + #0.941 #0.034 TRPSOMG_07_06 + #0.941 #0.034 &
          TRPSOMG_08_06 + #0.941 #0.034 TRPSOMG_09_06 + #0.941 #0.034 TRPSOMG_10_06 + &
          #0.941 #0.034 TRPSOMG_11_06 + #0.941 #0.034 TRPSOMG_04_07 + #0.941 &
          #0.034 TRPSOMG_05_07 + #0.941 #0.034 TRPSOMG_06_07 + #0.941 #0.034 &
          TRPSOMG_07_07 + #0.941 #0.034 TRPSOMG_08_07 + #0.941 #0.034 TRPSOMG_09_07 + &
          #0.941 #0.034 TRPSOMG_10_07 + #0.941 #0.034 TRPSOMG_11_07
 S232.5 TRPSOMG_12_07 + OH = #1.000 #0.036 TRPSOMG_01_01 + #1.000 #0.036 &
          TRPSOMG_02_01 + #1.000 #0.036 TRPSOMG_03_01 + #1.000 #0.036 TRPSOMG_04_01 + &
          #1.000 #0.036 TRPSOMG_05_01 + #1.000 #0.036 TRPSOMG_06_01 + #1.000 &
          #0.036 TRPSOMG_07_01 + #1.000 #0.036 TRPSOMG_08_01 + #1.000 #0.000 &
          TRPSOMG_09_01 + #1.000 #0.000 TRPSOMG_10_01 + #1.000 #0.000 TRPSOMG_11_01 + &
          #1.000 #0.036 TRPSOMG_01_02 + #1.000 #0.036 TRPSOMG_02_02 + #1.000 &
          #0.036 TRPSOMG_03_02 + #1.000 #0.036 TRPSOMG_04_02 + #1.000 #0.036 &
          TRPSOMG_05_02 + #1.000 #0.036 TRPSOMG_06_02 + #1.000 #0.036 TRPSOMG_07_02 + &
          #1.000 #0.036 TRPSOMG_08_02 + #1.000 #0.000 TRPSOMG_09_02 + #1.000 &
          #0.000 TRPSOMG_10_02 + #1.000 #0.000 TRPSOMG_11_02 + #1.000 #0.036 &
          TRPSOMG_02_03 + #1.000 #0.036 TRPSOMG_03_03 + #1.000 #0.036 TRPSOMG_04_03 + &
          #1.000 #0.036 TRPSOMG_05_03 + #1.000 #0.036 TRPSOMG_06_03 + #1.000 &
          #0.036 TRPSOMG_07_03 + #1.000 #0.036 TRPSOMG_08_03 + #1.000 #0.036 &
          TRPSOMG_09_03 + #1.000 #0.000 TRPSOMG_10_03 + #1.000 #0.000 TRPSOMG_11_03 + &
          #1.000 #0.036 TRPSOMG_02_04 + #1.000 #0.036 TRPSOMG_03_04 + #1.000 &
          #0.036 TRPSOMG_04_04 + #1.000 #0.036 TRPSOMG_05_04 + #1.000 #0.036 &
          TRPSOMG_06_04 + #1.000 #0.036 TRPSOMG_07_04 + #1.000 #0.036 TRPSOMG_08_04 + &
          #1.000 #0.036 TRPSOMG_09_04 + #1.000 #0.000 TRPSOMG_10_04 + #1.000 &
          #0.000 TRPSOMG_11_04 + #1.000 #0.036 TRPSOMG_03_05 + #1.000 #0.036 &
          TRPSOMG_04_05 + #1.000 #0.036 TRPSOMG_05_05 + #1.000 #0.036 TRPSOMG_06_05 + &
          #1.000 #0.036 TRPSOMG_07_05 + #1.000 #0.036 TRPSOMG_08_05 + #1.000 &
          #0.036 TRPSOMG_09_05 + #1.000 #0.036 TRPSOMG_10_05 + #1.000 #0.000 &
          TRPSOMG_11_05 + #1.000 #0.036 TRPSOMG_03_06 + #1.000 #0.036 TRPSOMG_04_06 + &
          #1.000 #0.036 TRPSOMG_05_06 + #1.000 #0.036 TRPSOMG_06_06 + #1.000 &
          #0.036 TRPSOMG_07_06 + #1.000 #0.036 TRPSOMG_08_06 + #1.000 #0.036 &
          TRPSOMG_09_06 + #1.000 #0.036 TRPSOMG_10_06 + #1.000 #0.000 TRPSOMG_11_06 + &
          #1.000 #0.036 TRPSOMG_04_07 + #1.000 #0.036 TRPSOMG_05_07 + #1.000 &
          #0.036 TRPSOMG_06_07 + #1.000 #0.036 TRPSOMG_07_07 + #1.000 #0.036 &
          TRPSOMG_08_07 + #1.000 #0.036 TRPSOMG_09_07 + #1.000 #0.036 TRPSOMG_10_07 + &
          #1.000 #0.036 TRPSOMG_11_07
 S233.5 TRPSOMG_13_01 + OH = #0.200 #0.361 TRPSOMG_13_02 + #0.200 #0.007 &
          TRPSOMG_13_03 + #0.200 #0.606 TRPSOMG_13_04 + #0.200 #0.026 TRPSOMG_13_05 + &
          #0.800 #0.083 TRPSOMG_01_01 + #0.800 #0.083 TRPSOMG_02_01 + #0.800 &
          #0.083 TRPSOMG_03_01 + #0.800 #0.083 TRPSOMG_04_01 + #0.800 #0.083 &
          TRPSOMG_05_01 + #0.800 #0.083 TRPSOMG_06_01 + #0.800 #0.083 TRPSOMG_07_01 + &
          #0.800 #0.083 TRPSOMG_08_01 + #0.800 #0.083 TRPSOMG_09_01 + #0.800 &
          #0.083 TRPSOMG_10_01 + #0.800 #0.083 TRPSOMG_11_01 + #0.800 #0.083 &
          TRPSOMG_12_01 + #0.800 #0.083 TRPSOMG_01_02 + #0.800 #0.083 TRPSOMG_02_02 + &
          #0.800 #0.083 TRPSOMG_03_02 + #0.800 #0.083 TRPSOMG_04_02 + #0.800 &
          #0.083 TRPSOMG_05_02 + #0.800 #0.083 TRPSOMG_06_02 + #0.800 #0.083 &
          TRPSOMG_07_02 + #0.800 #0.083 TRPSOMG_08_02 + #0.800 #0.083 TRPSOMG_09_02 + &
          #0.800 #0.083 TRPSOMG_10_02 + #0.800 #0.083 TRPSOMG_11_02 + #0.800 &
          #0.083 TRPSOMG_12_02
 S234.5 TRPSOMG_13_02 + OH = #0.150 #0.361 TRPSOMG_13_03 + #0.150 #0.007 &
          TRPSOMG_13_04 + #0.150 #0.606 TRPSOMG_13_05 + #0.150 #0.026 TRPSOMG_13_06 + &
          #0.850 #0.057 TRPSOMG_01_01 + #0.850 #0.057 TRPSOMG_02_01 + #0.850 &
          #0.057 TRPSOMG_03_01 + #0.850 #0.057 TRPSOMG_04_01 + #0.850 #0.057 &
          TRPSOMG_05_01 + #0.850 #0.057 TRPSOMG_06_01 + #0.850 #0.057 TRPSOMG_07_01 + &
          #0.850 #0.057 TRPSOMG_08_01 + #0.850 #0.057 TRPSOMG_09_01 + #0.850 &
          #0.057 TRPSOMG_10_01 + #0.850 #0.057 TRPSOMG_11_01 + #0.850 #0.000 &
          TRPSOMG_12_01 + #0.850 #0.057 TRPSOMG_01_02 + #0.850 #0.057 TRPSOMG_02_02 + &
          #0.850 #0.057 TRPSOMG_03_02 + #0.850 #0.057 TRPSOMG_04_02 + #0.850 &
          #0.057 TRPSOMG_05_02 + #0.850 #0.057 TRPSOMG_06_02 + #0.850 #0.057 &
          TRPSOMG_07_02 + #0.850 #0.057 TRPSOMG_08_02 + #0.850 #0.057 TRPSOMG_09_02 + &
          #0.850 #0.057 TRPSOMG_10_02 + #0.850 #0.057 TRPSOMG_11_02 + #0.850 &
          #0.057 TRPSOMG_12_02 + #0.850 #0.057 TRPSOMG_02_03 + #0.850 #0.057 &
          TRPSOMG_03_03 + #0.850 #0.057 TRPSOMG_04_03 + #0.850 #0.057 TRPSOMG_05_03 + &
          #0.850 #0.057 TRPSOMG_06_03 + #0.850 #0.057 TRPSOMG_07_03 + #0.850 &
          #0.057 TRPSOMG_08_03 + #0.850 #0.057 TRPSOMG_09_03 + #0.850 #0.057 &
          TRPSOMG_10_03 + #0.850 #0.057 TRPSOMG_11_03 + #0.850 #0.057 TRPSOMG_12_03
 S235.5 TRPSOMG_13_03 + OH = #0.120 #0.361 TRPSOMG_13_04 + #0.120 #0.007 &
          TRPSOMG_13_05 + #0.120 #0.606 TRPSOMG_13_06 + #0.120 #0.026 TRPSOMG_13_07 + &
          #0.880 #0.044 TRPSOMG_01_01 + #0.880 #0.044 TRPSOMG_02_01 + #0.880 &
          #0.044 TRPSOMG_03_01 + #0.880 #0.044 TRPSOMG_04_01 + #0.880 #0.044 &
          TRPSOMG_05_01 + #0.880 #0.044 TRPSOMG_06_01 + #0.880 #0.044 TRPSOMG_07_01 + &
          #0.880 #0.044 TRPSOMG_08_01 + #0.880 #0.044 TRPSOMG_09_01 + #0.880 &
          #0.044 TRPSOMG_10_01 + #0.880 #0.044 TRPSOMG_11_01 + #0.880 #0.000 &
          TRPSOMG_12_01 + #0.880 #0.044 TRPSOMG_01_02 + #0.880 #0.044 TRPSOMG_02_02 + &
          #0.880 #0.044 TRPSOMG_03_02 + #0.880 #0.044 TRPSOMG_04_02 + #0.880 &
          #0.044 TRPSOMG_05_02 + #0.880 #0.044 TRPSOMG_06_02 + #0.880 #0.044 &
          TRPSOMG_07_02 + #0.880 #0.044 TRPSOMG_08_02 + #0.880 #0.044 TRPSOMG_09_02 + &
          #0.880 #0.044 TRPSOMG_10_02 + #0.880 #0.044 TRPSOMG_11_02 + #0.880 &
          #0.000 TRPSOMG_12_02 + #0.880 #0.044 TRPSOMG_02_03 + #0.880 #0.044 &
          TRPSOMG_03_03 + #0.880 #0.044 TRPSOMG_04_03 + #0.880 #0.044 TRPSOMG_05_03 + &
          #0.880 #0.044 TRPSOMG_06_03 + #0.880 #0.044 TRPSOMG_07_03 + #0.880 &
          #0.044 TRPSOMG_08_03 + #0.880 #0.044 TRPSOMG_09_03 + #0.880 #0.044 &
          TRPSOMG_10_03 + #0.880 #0.044 TRPSOMG_11_03 + #0.880 #0.044 TRPSOMG_12_03 + &
          #0.880 #0.044 TRPSOMG_02_04 + #0.880 #0.044 TRPSOMG_03_04 + #0.880 &
          #0.044 TRPSOMG_04_04 + #0.880 #0.044 TRPSOMG_05_04 + #0.880 #0.044 &
          TRPSOMG_06_04 + #0.880 #0.044 TRPSOMG_07_04 + #0.880 #0.044 TRPSOMG_08_04 + &
          #0.880 #0.044 TRPSOMG_09_04 + #0.880 #0.044 TRPSOMG_10_04 + #0.880 &
          #0.044 TRPSOMG_11_04 + #0.880 #0.044 TRPSOMG_12_04
 S236.5 TRPSOMG_13_04 + OH = #0.097 #0.361 TRPSOMG_13_05 + #0.097 #0.007 &
          TRPSOMG_13_06 + #0.097 #0.632 TRPSOMG_13_07 + #0.903 #0.037 TRPSOMG_01_01 + &
          #0.903 #0.037 TRPSOMG_02_01 + #0.903 #0.037 TRPSOMG_03_01 + #0.903 &
          #0.037 TRPSOMG_04_01 + #0.903 #0.037 TRPSOMG_05_01 + #0.903 #0.037 &
          TRPSOMG_06_01 + #0.903 #0.037 TRPSOMG_07_01 + #0.903 #0.037 TRPSOMG_08_01 + &
          #0.903 #0.037 TRPSOMG_09_01 + #0.903 #0.037 TRPSOMG_10_01 + #0.903 &
          #0.000 TRPSOMG_11_01 + #0.903 #0.000 TRPSOMG_12_01 + #0.903 #0.037 &
          TRPSOMG_01_02 + #0.903 #0.037 TRPSOMG_02_02 + #0.903 #0.037 TRPSOMG_03_02 + &
          #0.903 #0.037 TRPSOMG_04_02 + #0.903 #0.037 TRPSOMG_05_02 + #0.903 &
          #0.037 TRPSOMG_06_02 + #0.903 #0.037 TRPSOMG_07_02 + #0.903 #0.037 &
          TRPSOMG_08_02 + #0.903 #0.037 TRPSOMG_09_02 + #0.903 #0.037 TRPSOMG_10_02 + &
          #0.903 #0.037 TRPSOMG_11_02 + #0.903 #0.000 TRPSOMG_12_02 + #0.903 &
          #0.037 TRPSOMG_02_03 + #0.903 #0.037 TRPSOMG_03_03 + #0.903 #0.037 &
          TRPSOMG_04_03 + #0.903 #0.037 TRPSOMG_05_03 + #0.903 #0.037 TRPSOMG_06_03 + &
          #0.903 #0.037 TRPSOMG_07_03 + #0.903 #0.037 TRPSOMG_08_03 + #0.903 &
          #0.037 TRPSOMG_09_03 + #0.903 #0.037 TRPSOMG_10_03 + #0.903 #0.037 &
          TRPSOMG_11_03 + #0.903 #0.000 TRPSOMG_12_03 + #0.903 #0.037 TRPSOMG_02_04 + &
          #0.903 #0.037 TRPSOMG_03_04 + #0.903 #0.037 TRPSOMG_04_04 + #0.903 &
          #0.037 TRPSOMG_05_04 + #0.903 #0.037 TRPSOMG_06_04 + #0.903 #0.037 &
          TRPSOMG_07_04 + #0.903 #0.037 TRPSOMG_08_04 + #0.903 #0.037 TRPSOMG_09_04 + &
          #0.903 #0.037 TRPSOMG_10_04 + #0.903 #0.037 TRPSOMG_11_04 + #0.903 &
          #0.037 TRPSOMG_12_04 + #0.903 #0.037 TRPSOMG_03_05 + #0.903 #0.037 &
          TRPSOMG_04_05 + #0.903 #0.037 TRPSOMG_05_05 + #0.903 #0.037 TRPSOMG_06_05 + &
          #0.903 #0.037 TRPSOMG_07_05 + #0.903 #0.037 TRPSOMG_08_05 + #0.903 &
          #0.037 TRPSOMG_09_05 + #0.903 #0.037 TRPSOMG_10_05 + #0.903 #0.037 &
          TRPSOMG_11_05 + #0.903 #0.037 TRPSOMG_12_05
 S237.5 TRPSOMG_13_05 + OH = #0.080 #0.361 TRPSOMG_13_06 + #0.080 #0.639 &
          TRPSOMG_13_07 + #0.920 #0.032 TRPSOMG_01_01 + #0.920 #0.032 TRPSOMG_02_01 + &
          #0.920 #0.032 TRPSOMG_03_01 + #0.920 #0.032 TRPSOMG_04_01 + #0.920 &
          #0.032 TRPSOMG_05_01 + #0.920 #0.032 TRPSOMG_06_01 + #0.920 #0.032 &
          TRPSOMG_07_01 + #0.920 #0.032 TRPSOMG_08_01 + #0.920 #0.032 TRPSOMG_09_01 + &
          #0.920 #0.032 TRPSOMG_10_01 + #0.920 #0.000 TRPSOMG_11_01 + #0.920 &
          #0.000 TRPSOMG_12_01 + #0.920 #0.032 TRPSOMG_01_02 + #0.920 #0.032 &
          TRPSOMG_02_02 + #0.920 #0.032 TRPSOMG_03_02 + #0.920 #0.032 TRPSOMG_04_02 + &
          #0.920 #0.032 TRPSOMG_05_02 + #0.920 #0.032 TRPSOMG_06_02 + #0.920 &
          #0.032 TRPSOMG_07_02 + #0.920 #0.032 TRPSOMG_08_02 + #0.920 #0.032 &
          TRPSOMG_09_02 + #0.920 #0.032 TRPSOMG_10_02 + #0.920 #0.000 TRPSOMG_11_02 + &
          #0.920 #0.000 TRPSOMG_12_02 + #0.920 #0.032 TRPSOMG_02_03 + #0.920 &
          #0.032 TRPSOMG_03_03 + #0.920 #0.032 TRPSOMG_04_03 + #0.920 #0.032 &
          TRPSOMG_05_03 + #0.920 #0.032 TRPSOMG_06_03 + #0.920 #0.032 TRPSOMG_07_03 + &
          #0.920 #0.032 TRPSOMG_08_03 + #0.920 #0.032 TRPSOMG_09_03 + #0.920 &
          #0.032 TRPSOMG_10_03 + #0.920 #0.032 TRPSOMG_11_03 + #0.920 #0.000 &
          TRPSOMG_12_03 + #0.920 #0.032 TRPSOMG_02_04 + #0.920 #0.032 TRPSOMG_03_04 + &
          #0.920 #0.032 TRPSOMG_04_04 + #0.920 #0.032 TRPSOMG_05_04 + #0.920 &
          #0.032 TRPSOMG_06_04 + #0.920 #0.032 TRPSOMG_07_04 + #0.920 #0.032 &
          TRPSOMG_08_04 + #0.920 #0.032 TRPSOMG_09_04 + #0.920 #0.032 TRPSOMG_10_04 + &
          #0.920 #0.032 TRPSOMG_11_04 + #0.920 #0.000 TRPSOMG_12_04 + #0.920 &
          #0.032 TRPSOMG_03_05 + #0.920 #0.032 TRPSOMG_04_05 + #0.920 #0.032 &
          TRPSOMG_05_05 + #0.920 #0.032 TRPSOMG_06_05 + #0.920 #0.032 TRPSOMG_07_05 + &
          #0.920 #0.032 TRPSOMG_08_05 + #0.920 #0.032 TRPSOMG_09_05 + #0.920 &
          #0.032 TRPSOMG_10_05 + #0.920 #0.032 TRPSOMG_11_05 + #0.920 #0.032 &
          TRPSOMG_12_05 + #0.920 #0.032 TRPSOMG_03_06 + #0.920 #0.032 TRPSOMG_04_06 + &
          #0.920 #0.032 TRPSOMG_05_06 + #0.920 #0.032 TRPSOMG_06_06 + #0.920 &
          #0.032 TRPSOMG_07_06 + #0.920 #0.032 TRPSOMG_08_06 + #0.920 #0.032 &
          TRPSOMG_09_06 + #0.920 #0.032 TRPSOMG_10_06 + #0.920 #0.032 TRPSOMG_11_06 + &
          #0.920 #0.032 TRPSOMG_12_06
 S238.5 TRPSOMG_13_06 + OH = #0.065 #1.000 TRPSOMG_13_07 + #0.935 #0.029 &
          TRPSOMG_01_01 + #0.935 #0.029 TRPSOMG_02_01 + #0.935 #0.029 TRPSOMG_03_01 + &
          #0.935 #0.029 TRPSOMG_04_01 + #0.935 #0.029 TRPSOMG_05_01 + #0.935 &
          #0.029 TRPSOMG_06_01 + #0.935 #0.029 TRPSOMG_07_01 + #0.935 #0.029 &
          TRPSOMG_08_01 + #0.935 #0.029 TRPSOMG_09_01 + #0.935 #0.000 TRPSOMG_10_01 + &
          #0.935 #0.000 TRPSOMG_11_01 + #0.935 #0.000 TRPSOMG_12_01 + #0.935 &
          #0.029 TRPSOMG_01_02 + #0.935 #0.029 TRPSOMG_02_02 + #0.935 #0.029 &
          TRPSOMG_03_02 + #0.935 #0.029 TRPSOMG_04_02 + #0.935 #0.029 TRPSOMG_05_02 + &
          #0.935 #0.029 TRPSOMG_06_02 + #0.935 #0.029 TRPSOMG_07_02 + #0.935 &
          #0.029 TRPSOMG_08_02 + #0.935 #0.029 TRPSOMG_09_02 + #0.935 #0.029 &
          TRPSOMG_10_02 + #0.935 #0.000 TRPSOMG_11_02 + #0.935 #0.000 TRPSOMG_12_02 + &
          #0.935 #0.029 TRPSOMG_02_03 + #0.935 #0.029 TRPSOMG_03_03 + #0.935 &
          #0.029 TRPSOMG_04_03 + #0.935 #0.029 TRPSOMG_05_03 + #0.935 #0.029 &
          TRPSOMG_06_03 + #0.935 #0.029 TRPSOMG_07_03 + #0.935 #0.029 TRPSOMG_08_03 + &
          #0.935 #0.029 TRPSOMG_09_03 + #0.935 #0.029 TRPSOMG_10_03 + #0.935 &
          #0.000 TRPSOMG_11_03 + #0.935 #0.000 TRPSOMG_12_03 + #0.935 #0.029 &
          TRPSOMG_02_04 + #0.935 #0.029 TRPSOMG_03_04 + #0.935 #0.029 TRPSOMG_04_04 + &
          #0.935 #0.029 TRPSOMG_05_04 + #0.935 #0.029 TRPSOMG_06_04 + #0.935 &
          #0.029 TRPSOMG_07_04 + #0.935 #0.029 TRPSOMG_08_04 + #0.935 #0.029 &
          TRPSOMG_09_04 + #0.935 #0.029 TRPSOMG_10_04 + #0.935 #0.029 TRPSOMG_11_04 + &
          #0.935 #0.000 TRPSOMG_12_04 + #0.935 #0.029 TRPSOMG_03_05 + #0.935 &
          #0.029 TRPSOMG_04_05 + #0.935 #0.029 TRPSOMG_05_05 + #0.935 #0.029 &
          TRPSOMG_06_05 + #0.935 #0.029 TRPSOMG_07_05 + #0.935 #0.029 TRPSOMG_08_05 + &
          #0.935 #0.029 TRPSOMG_09_05 + #0.935 #0.029 TRPSOMG_10_05 + #0.935 &
          #0.029 TRPSOMG_11_05 + #0.935 #0.000 TRPSOMG_12_05 + #0.935 #0.029 &
          TRPSOMG_03_06 + #0.935 #0.029 TRPSOMG_04_06 + #0.935 #0.029 TRPSOMG_05_06 + &
          #0.935 #0.029 TRPSOMG_06_06 + #0.935 #0.029 TRPSOMG_07_06 + #0.935 &
          #0.029 TRPSOMG_08_06 + #0.935 #0.029 TRPSOMG_09_06 + #0.935 #0.029 &
          TRPSOMG_10_06 + #0.935 #0.029 TRPSOMG_11_06 + #0.935 #0.029 TRPSOMG_12_06 + &
          #0.935 #0.029 TRPSOMG_04_07 + #0.935 #0.029 TRPSOMG_05_07 + #0.935 &
          #0.029 TRPSOMG_06_07 + #0.935 #0.029 TRPSOMG_07_07 + #0.935 #0.029 &
          TRPSOMG_08_07 + #0.935 #0.029 TRPSOMG_09_07 + #0.935 #0.029 TRPSOMG_10_07 + &
          #0.935 #0.029 TRPSOMG_11_07 + #0.935 #0.029 TRPSOMG_12_07
 S239.5 TRPSOMG_13_07 + OH = #1.000 #0.030 TRPSOMG_01_01 + #1.000 #0.030 &
          TRPSOMG_02_01 + #1.000 #0.030 TRPSOMG_03_01 + #1.000 #0.030 TRPSOMG_04_01 + &
          #1.000 #0.030 TRPSOMG_05_01 + #1.000 #0.030 TRPSOMG_06_01 + #1.000 &
          #0.030 TRPSOMG_07_01 + #1.000 #0.030 TRPSOMG_08_01 + #1.000 #0.030 &
          TRPSOMG_09_01 + #1.000 #0.000 TRPSOMG_10_01 + #1.000 #0.000 TRPSOMG_11_01 + &
          #1.000 #0.000 TRPSOMG_12_01 + #1.000 #0.030 TRPSOMG_01_02 + #1.000 &
          #0.030 TRPSOMG_02_02 + #1.000 #0.030 TRPSOMG_03_02 + #1.000 #0.030 &
          TRPSOMG_04_02 + #1.000 #0.030 TRPSOMG_05_02 + #1.000 #0.030 TRPSOMG_06_02 + &
          #1.000 #0.030 TRPSOMG_07_02 + #1.000 #0.030 TRPSOMG_08_02 + #1.000 &
          #0.030 TRPSOMG_09_02 + #1.000 #0.000 TRPSOMG_10_02 + #1.000 #0.000 &
          TRPSOMG_11_02 + #1.000 #0.000 TRPSOMG_12_02 + #1.000 #0.030 TRPSOMG_02_03 + &
          #1.000 #0.030 TRPSOMG_03_03 + #1.000 #0.030 TRPSOMG_04_03 + #1.000 &
          #0.030 TRPSOMG_05_03 + #1.000 #0.030 TRPSOMG_06_03 + #1.000 #0.030 &
          TRPSOMG_07_03 + #1.000 #0.030 TRPSOMG_08_03 + #1.000 #0.030 TRPSOMG_09_03 + &
          #1.000 #0.030 TRPSOMG_10_03 + #1.000 #0.000 TRPSOMG_11_03 + #1.000 &
          #0.000 TRPSOMG_12_03 + #1.000 #0.030 TRPSOMG_02_04 + #1.000 #0.030 &
          TRPSOMG_03_04 + #1.000 #0.030 TRPSOMG_04_04 + #1.000 #0.030 TRPSOMG_05_04 + &
          #1.000 #0.030 TRPSOMG_06_04 + #1.000 #0.030 TRPSOMG_07_04 + #1.000 &
          #0.030 TRPSOMG_08_04 + #1.000 #0.030 TRPSOMG_09_04 + #1.000 #0.030 &
          TRPSOMG_10_04 + #1.000 #0.000 TRPSOMG_11_04 + #1.000 #0.000 TRPSOMG_12_04 + &
          #1.000 #0.030 TRPSOMG_03_05 + #1.000 #0.030 TRPSOMG_04_05 + #1.000 &
          #0.030 TRPSOMG_05_05 + #1.000 #0.030 TRPSOMG_06_05 + #1.000 #0.030 &
          TRPSOMG_07_05 + #1.000 #0.030 TRPSOMG_08_05 + #1.000 #0.030 TRPSOMG_09_05 + &
          #1.000 #0.030 TRPSOMG_10_05 + #1.000 #0.030 TRPSOMG_11_05 + #1.000 &
          #0.000 TRPSOMG_12_05 + #1.000 #0.030 TRPSOMG_03_06 + #1.000 #0.030 &
          TRPSOMG_04_06 + #1.000 #0.030 TRPSOMG_05_06 + #1.000 #0.030 TRPSOMG_06_06 + &
          #1.000 #0.030 TRPSOMG_07_06 + #1.000 #0.030 TRPSOMG_08_06 + #1.000 &
          #0.030 TRPSOMG_09_06 + #1.000 #0.030 TRPSOMG_10_06 + #1.000 #0.030 &
          TRPSOMG_11_06 + #1.000 #0.000 TRPSOMG_12_06 + #1.000 #0.030 TRPSOMG_04_07 + &
          #1.000 #0.030 TRPSOMG_05_07 + #1.000 #0.030 TRPSOMG_06_07 + #1.000 &
          #0.030 TRPSOMG_07_07 + #1.000 #0.030 TRPSOMG_08_07 + #1.000 #0.030 &
          TRPSOMG_09_07 + #1.000 #0.030 TRPSOMG_10_07 + #1.000 #0.030 TRPSOMG_11_07 + &
          #1.000 #0.030 TRPSOMG_12_07
 S240.5 TRPSOMG_14_01 + OH = #0.205 #0.361 TRPSOMG_14_02 + #0.205 #0.007 &
          TRPSOMG_14_03 + #0.205 #0.606 TRPSOMG_14_04 + #0.205 #0.026 TRPSOMG_14_05 + &
          #0.795 #0.077 TRPSOMG_01_01 + #0.795 #0.077 TRPSOMG_02_01 + #0.795 &
          #0.077 TRPSOMG_03_01 + #0.795 #0.077 TRPSOMG_04_01 + #0.795 #0.077 &
          TRPSOMG_05_01 + #0.795 #0.077 TRPSOMG_06_01 + #0.795 #0.077 TRPSOMG_07_01 + &
          #0.795 #0.077 TRPSOMG_08_01 + #0.795 #0.077 TRPSOMG_09_01 + #0.795 &
          #0.077 TRPSOMG_10_01 + #0.795 #0.077 TRPSOMG_11_01 + #0.795 #0.077 &
          TRPSOMG_12_01 + #0.795 #0.077 TRPSOMG_13_01 + #0.795 #0.077 TRPSOMG_01_02 + &
          #0.795 #0.077 TRPSOMG_02_02 + #0.795 #0.077 TRPSOMG_03_02 + #0.795 &
          #0.077 TRPSOMG_04_02 + #0.795 #0.077 TRPSOMG_05_02 + #0.795 #0.077 &
          TRPSOMG_06_02 + #0.795 #0.077 TRPSOMG_07_02 + #0.795 #0.077 TRPSOMG_08_02 + &
          #0.795 #0.077 TRPSOMG_09_02 + #0.795 #0.077 TRPSOMG_10_02 + #0.795 &
          #0.077 TRPSOMG_11_02 + #0.795 #0.077 TRPSOMG_12_02 + #0.795 #0.077 &
          TRPSOMG_13_02
 S241.5 TRPSOMG_14_02 + OH = #0.156 #0.361 TRPSOMG_14_03 + #0.156 #0.007 &
          TRPSOMG_14_04 + #0.156 #0.606 TRPSOMG_14_05 + #0.156 #0.026 TRPSOMG_14_06 + &
          #0.844 #0.054 TRPSOMG_01_01 + #0.844 #0.054 TRPSOMG_02_01 + #0.844 &
          #0.054 TRPSOMG_03_01 + #0.844 #0.054 TRPSOMG_04_01 + #0.844 #0.054 &
          TRPSOMG_05_01 + #0.844 #0.054 TRPSOMG_06_01 + #0.844 #0.054 TRPSOMG_07_01 + &
          #0.844 #0.054 TRPSOMG_08_01 + #0.844 #0.054 TRPSOMG_09_01 + #0.844 &
          #0.054 TRPSOMG_10_01 + #0.844 #0.054 TRPSOMG_11_01 + #0.844 #0.054 &
          TRPSOMG_12_01 + #0.844 #0.000 TRPSOMG_13_01 + #0.844 #0.054 TRPSOMG_01_02 + &
          #0.844 #0.054 TRPSOMG_02_02 + #0.844 #0.054 TRPSOMG_03_02 + #0.844 &
          #0.054 TRPSOMG_04_02 + #0.844 #0.054 TRPSOMG_05_02 + #0.844 #0.054 &
          TRPSOMG_06_02 + #0.844 #0.054 TRPSOMG_07_02 + #0.844 #0.054 TRPSOMG_08_02 + &
          #0.844 #0.054 TRPSOMG_09_02 + #0.844 #0.054 TRPSOMG_10_02 + #0.844 &
          #0.054 TRPSOMG_11_02 + #0.844 #0.054 TRPSOMG_12_02 + #0.844 #0.054 &
          TRPSOMG_13_02 + #0.844 #0.054 TRPSOMG_02_03 + #0.844 #0.054 TRPSOMG_03_03 + &
          #0.844 #0.054 TRPSOMG_04_03 + #0.844 #0.054 TRPSOMG_05_03 + #0.844 &
          #0.054 TRPSOMG_06_03 + #0.844 #0.054 TRPSOMG_07_03 + #0.844 #0.054 &
          TRPSOMG_08_03 + #0.844 #0.054 TRPSOMG_09_03 + #0.844 #0.054 TRPSOMG_10_03 + &
          #0.844 #0.054 TRPSOMG_11_03 + #0.844 #0.054 TRPSOMG_12_03 + #0.844 &
          #0.054 TRPSOMG_13_03
 S242.5 TRPSOMG_14_03 + OH = #0.125 #0.361 TRPSOMG_14_04 + #0.125 #0.007 &
          TRPSOMG_14_05 + #0.125 #0.606 TRPSOMG_14_06 + #0.125 #0.026 TRPSOMG_14_07 + &
          #0.875 #0.042 TRPSOMG_01_01 + #0.875 #0.042 TRPSOMG_02_01 + #0.875 &
          #0.042 TRPSOMG_03_01 + #0.875 #0.042 TRPSOMG_04_01 + #0.875 #0.042 &
          TRPSOMG_05_01 + #0.875 #0.042 TRPSOMG_06_01 + #0.875 #0.042 TRPSOMG_07_01 + &
          #0.875 #0.042 TRPSOMG_08_01 + #0.875 #0.042 TRPSOMG_09_01 + #0.875 &
          #0.042 TRPSOMG_10_01 + #0.875 #0.042 TRPSOMG_11_01 + #0.875 #0.042 &
          TRPSOMG_12_01 + #0.875 #0.000 TRPSOMG_13_01 + #0.875 #0.042 TRPSOMG_01_02 + &
          #0.875 #0.042 TRPSOMG_02_02 + #0.875 #0.042 TRPSOMG_03_02 + #0.875 &
          #0.042 TRPSOMG_04_02 + #0.875 #0.042 TRPSOMG_05_02 + #0.875 #0.042 &
          TRPSOMG_06_02 + #0.875 #0.042 TRPSOMG_07_02 + #0.875 #0.042 TRPSOMG_08_02 + &
          #0.875 #0.042 TRPSOMG_09_02 + #0.875 #0.042 TRPSOMG_10_02 + #0.875 &
          #0.042 TRPSOMG_11_02 + #0.875 #0.042 TRPSOMG_12_02 + #0.875 #0.000 &
          TRPSOMG_13_02 + #0.875 #0.042 TRPSOMG_02_03 + #0.875 #0.042 TRPSOMG_03_03 + &
          #0.875 #0.042 TRPSOMG_04_03 + #0.875 #0.042 TRPSOMG_05_03 + #0.875 &
          #0.042 TRPSOMG_06_03 + #0.875 #0.042 TRPSOMG_07_03 + #0.875 #0.042 &
          TRPSOMG_08_03 + #0.875 #0.042 TRPSOMG_09_03 + #0.875 #0.042 TRPSOMG_10_03 + &
          #0.875 #0.042 TRPSOMG_11_03 + #0.875 #0.042 TRPSOMG_12_03 + #0.875 &
          #0.042 TRPSOMG_13_03 + #0.875 #0.042 TRPSOMG_02_04 + #0.875 #0.042 &
          TRPSOMG_03_04 + #0.875 #0.042 TRPSOMG_04_04 + #0.875 #0.042 TRPSOMG_05_04 + &
          #0.875 #0.042 TRPSOMG_06_04 + #0.875 #0.042 TRPSOMG_07_04 + #0.875 &
          #0.042 TRPSOMG_08_04 + #0.875 #0.042 TRPSOMG_09_04 + #0.875 #0.042 &
          TRPSOMG_10_04 + #0.875 #0.042 TRPSOMG_11_04 + #0.875 #0.042 TRPSOMG_12_04 + &
          #0.875 #0.042 TRPSOMG_13_04
 S243.5 TRPSOMG_14_04 + OH = #0.103 #0.361 TRPSOMG_14_05 + #0.103 #0.007 &
          TRPSOMG_14_06 + #0.103 #0.632 TRPSOMG_14_07 + #0.897 #0.034 TRPSOMG_01_01 + &
          #0.897 #0.034 TRPSOMG_02_01 + #0.897 #0.034 TRPSOMG_03_01 + #0.897 &
          #0.034 TRPSOMG_04_01 + #0.897 #0.034 TRPSOMG_05_01 + #0.897 #0.034 &
          TRPSOMG_06_01 + #0.897 #0.034 TRPSOMG_07_01 + #0.897 #0.034 TRPSOMG_08_01 + &
          #0.897 #0.034 TRPSOMG_09_01 + #0.897 #0.034 TRPSOMG_10_01 + #0.897 &
          #0.034 TRPSOMG_11_01 + #0.897 #0.000 TRPSOMG_12_01 + #0.897 #0.000 &
          TRPSOMG_13_01 + #0.897 #0.034 TRPSOMG_01_02 + #0.897 #0.034 TRPSOMG_02_02 + &
          #0.897 #0.034 TRPSOMG_03_02 + #0.897 #0.034 TRPSOMG_04_02 + #0.897 &
          #0.034 TRPSOMG_05_02 + #0.897 #0.034 TRPSOMG_06_02 + #0.897 #0.034 &
          TRPSOMG_07_02 + #0.897 #0.034 TRPSOMG_08_02 + #0.897 #0.034 TRPSOMG_09_02 + &
          #0.897 #0.034 TRPSOMG_10_02 + #0.897 #0.034 TRPSOMG_11_02 + #0.897 &
          #0.034 TRPSOMG_12_02 + #0.897 #0.000 TRPSOMG_13_02 + #0.897 #0.034 &
          TRPSOMG_02_03 + #0.897 #0.034 TRPSOMG_03_03 + #0.897 #0.034 TRPSOMG_04_03 + &
          #0.897 #0.034 TRPSOMG_05_03 + #0.897 #0.034 TRPSOMG_06_03 + #0.897 &
          #0.034 TRPSOMG_07_03 + #0.897 #0.034 TRPSOMG_08_03 + #0.897 #0.034 &
          TRPSOMG_09_03 + #0.897 #0.034 TRPSOMG_10_03 + #0.897 #0.034 TRPSOMG_11_03 + &
          #0.897 #0.034 TRPSOMG_12_03 + #0.897 #0.000 TRPSOMG_13_03 + #0.897 &
          #0.034 TRPSOMG_02_04 + #0.897 #0.034 TRPSOMG_03_04 + #0.897 #0.034 &
          TRPSOMG_04_04 + #0.897 #0.034 TRPSOMG_05_04 + #0.897 #0.034 TRPSOMG_06_04 + &
          #0.897 #0.034 TRPSOMG_07_04 + #0.897 #0.034 TRPSOMG_08_04 + #0.897 &
          #0.034 TRPSOMG_09_04 + #0.897 #0.034 TRPSOMG_10_04 + #0.897 #0.034 &
          TRPSOMG_11_04 + #0.897 #0.034 TRPSOMG_12_04 + #0.897 #0.034 TRPSOMG_13_04 + &
          #0.897 #0.034 TRPSOMG_03_05 + #0.897 #0.034 TRPSOMG_04_05 + #0.897 &
          #0.034 TRPSOMG_05_05 + #0.897 #0.034 TRPSOMG_06_05 + #0.897 #0.034 &
          TRPSOMG_07_05 + #0.897 #0.034 TRPSOMG_08_05 + #0.897 #0.034 TRPSOMG_09_05 + &
          #0.897 #0.034 TRPSOMG_10_05 + #0.897 #0.034 TRPSOMG_11_05 + #0.897 &
          #0.034 TRPSOMG_12_05 + #0.897 #0.034 TRPSOMG_13_05
 S244.5 TRPSOMG_14_05 + OH = #0.086 #0.361 TRPSOMG_14_06 + #0.086 #0.639 &
          TRPSOMG_14_07 + #0.914 #0.030 TRPSOMG_01_01 + #0.914 #0.030 TRPSOMG_02_01 + &
          #0.914 #0.030 TRPSOMG_03_01 + #0.914 #0.030 TRPSOMG_04_01 + #0.914 &
          #0.030 TRPSOMG_05_01 + #0.914 #0.030 TRPSOMG_06_01 + #0.914 #0.030 &
          TRPSOMG_07_01 + #0.914 #0.030 TRPSOMG_08_01 + #0.914 #0.030 TRPSOMG_09_01 + &
          #0.914 #0.030 TRPSOMG_10_01 + #0.914 #0.030 TRPSOMG_11_01 + #0.914 &
          #0.000 TRPSOMG_12_01 + #0.914 #0.000 TRPSOMG_13_01 + #0.914 #0.030 &
          TRPSOMG_01_02 + #0.914 #0.030 TRPSOMG_02_02 + #0.914 #0.030 TRPSOMG_03_02 + &
          #0.914 #0.030 TRPSOMG_04_02 + #0.914 #0.030 TRPSOMG_05_02 + #0.914 &
          #0.030 TRPSOMG_06_02 + #0.914 #0.030 TRPSOMG_07_02 + #0.914 #0.030 &
          TRPSOMG_08_02 + #0.914 #0.030 TRPSOMG_09_02 + #0.914 #0.030 TRPSOMG_10_02 + &
          #0.914 #0.030 TRPSOMG_11_02 + #0.914 #0.000 TRPSOMG_12_02 + #0.914 &
          #0.000 TRPSOMG_13_02 + #0.914 #0.030 TRPSOMG_02_03 + #0.914 #0.030 &
          TRPSOMG_03_03 + #0.914 #0.030 TRPSOMG_04_03 + #0.914 #0.030 TRPSOMG_05_03 + &
          #0.914 #0.030 TRPSOMG_06_03 + #0.914 #0.030 TRPSOMG_07_03 + #0.914 &
          #0.030 TRPSOMG_08_03 + #0.914 #0.030 TRPSOMG_09_03 + #0.914 #0.030 &
          TRPSOMG_10_03 + #0.914 #0.030 TRPSOMG_11_03 + #0.914 #0.030 TRPSOMG_12_03 + &
          #0.914 #0.000 TRPSOMG_13_03 + #0.914 #0.030 TRPSOMG_02_04 + #0.914 &
          #0.030 TRPSOMG_03_04 + #0.914 #0.030 TRPSOMG_04_04 + #0.914 #0.030 &
          TRPSOMG_05_04 + #0.914 #0.030 TRPSOMG_06_04 + #0.914 #0.030 TRPSOMG_07_04 + &
          #0.914 #0.030 TRPSOMG_08_04 + #0.914 #0.030 TRPSOMG_09_04 + #0.914 &
          #0.030 TRPSOMG_10_04 + #0.914 #0.030 TRPSOMG_11_04 + #0.914 #0.030 &
          TRPSOMG_12_04 + #0.914 #0.000 TRPSOMG_13_04 + #0.914 #0.030 TRPSOMG_03_05 + &
          #0.914 #0.030 TRPSOMG_04_05 + #0.914 #0.030 TRPSOMG_05_05 + #0.914 &
          #0.030 TRPSOMG_06_05 + #0.914 #0.030 TRPSOMG_07_05 + #0.914 #0.030 &
          TRPSOMG_08_05 + #0.914 #0.030 TRPSOMG_09_05 + #0.914 #0.030 TRPSOMG_10_05 + &
          #0.914 #0.030 TRPSOMG_11_05 + #0.914 #0.030 TRPSOMG_12_05 + #0.914 &
          #0.030 TRPSOMG_13_05 + #0.914 #0.030 TRPSOMG_03_06 + #0.914 #0.030 &
          TRPSOMG_04_06 + #0.914 #0.030 TRPSOMG_05_06 + #0.914 #0.030 TRPSOMG_06_06 + &
          #0.914 #0.030 TRPSOMG_07_06 + #0.914 #0.030 TRPSOMG_08_06 + #0.914 &
          #0.030 TRPSOMG_09_06 + #0.914 #0.030 TRPSOMG_10_06 + #0.914 #0.030 &
          TRPSOMG_11_06 + #0.914 #0.030 TRPSOMG_12_06 + #0.914 #0.030 TRPSOMG_13_06
 S245.5 TRPSOMG_14_06 + OH = #0.071 #1.000 TRPSOMG_14_07 + #0.929 #0.027 &
          TRPSOMG_01_01 + #0.929 #0.027 TRPSOMG_02_01 + #0.929 #0.027 TRPSOMG_03_01 + &
          #0.929 #0.027 TRPSOMG_04_01 + #0.929 #0.027 TRPSOMG_05_01 + #0.929 &
          #0.027 TRPSOMG_06_01 + #0.929 #0.027 TRPSOMG_07_01 + #0.929 #0.027 &
          TRPSOMG_08_01 + #0.929 #0.027 TRPSOMG_09_01 + #0.929 #0.027 TRPSOMG_10_01 + &
          #0.929 #0.000 TRPSOMG_11_01 + #0.929 #0.000 TRPSOMG_12_01 + #0.929 &
          #0.000 TRPSOMG_13_01 + #0.929 #0.027 TRPSOMG_01_02 + #0.929 #0.027 &
          TRPSOMG_02_02 + #0.929 #0.027 TRPSOMG_03_02 + #0.929 #0.027 TRPSOMG_04_02 + &
          #0.929 #0.027 TRPSOMG_05_02 + #0.929 #0.027 TRPSOMG_06_02 + #0.929 &
          #0.027 TRPSOMG_07_02 + #0.929 #0.027 TRPSOMG_08_02 + #0.929 #0.027 &
          TRPSOMG_09_02 + #0.929 #0.027 TRPSOMG_10_02 + #0.929 #0.027 TRPSOMG_11_02 + &
          #0.929 #0.000 TRPSOMG_12_02 + #0.929 #0.000 TRPSOMG_13_02 + #0.929 &
          #0.027 TRPSOMG_02_03 + #0.929 #0.027 TRPSOMG_03_03 + #0.929 #0.027 &
          TRPSOMG_04_03 + #0.929 #0.027 TRPSOMG_05_03 + #0.929 #0.027 TRPSOMG_06_03 + &
          #0.929 #0.027 TRPSOMG_07_03 + #0.929 #0.027 TRPSOMG_08_03 + #0.929 &
          #0.027 TRPSOMG_09_03 + #0.929 #0.027 TRPSOMG_10_03 + #0.929 #0.027 &
          TRPSOMG_11_03 + #0.929 #0.000 TRPSOMG_12_03 + #0.929 #0.000 TRPSOMG_13_03 + &
          #0.929 #0.027 TRPSOMG_02_04 + #0.929 #0.027 TRPSOMG_03_04 + #0.929 &
          #0.027 TRPSOMG_04_04 + #0.929 #0.027 TRPSOMG_05_04 + #0.929 #0.027 &
          TRPSOMG_06_04 + #0.929 #0.027 TRPSOMG_07_04 + #0.929 #0.027 TRPSOMG_08_04 + &
          #0.929 #0.027 TRPSOMG_09_04 + #0.929 #0.027 TRPSOMG_10_04 + #0.929 &
          #0.027 TRPSOMG_11_04 + #0.929 #0.027 TRPSOMG_12_04 + #0.929 #0.000 &
          TRPSOMG_13_04 + #0.929 #0.027 TRPSOMG_03_05 + #0.929 #0.027 TRPSOMG_04_05 + &
          #0.929 #0.027 TRPSOMG_05_05 + #0.929 #0.027 TRPSOMG_06_05 + #0.929 &
          #0.027 TRPSOMG_07_05 + #0.929 #0.027 TRPSOMG_08_05 + #0.929 #0.027 &
          TRPSOMG_09_05 + #0.929 #0.027 TRPSOMG_10_05 + #0.929 #0.027 TRPSOMG_11_05 + &
          #0.929 #0.027 TRPSOMG_12_05 + #0.929 #0.000 TRPSOMG_13_05 + #0.929 &
          #0.027 TRPSOMG_03_06 + #0.929 #0.027 TRPSOMG_04_06 + #0.929 #0.027 &
          TRPSOMG_05_06 + #0.929 #0.027 TRPSOMG_06_06 + #0.929 #0.027 TRPSOMG_07_06 + &
          #0.929 #0.027 TRPSOMG_08_06 + #0.929 #0.027 TRPSOMG_09_06 + #0.929 &
          #0.027 TRPSOMG_10_06 + #0.929 #0.027 TRPSOMG_11_06 + #0.929 #0.027 &
          TRPSOMG_12_06 + #0.929 #0.027 TRPSOMG_13_06 + #0.929 #0.027 TRPSOMG_04_07 + &
          #0.929 #0.027 TRPSOMG_05_07 + #0.929 #0.027 TRPSOMG_06_07 + #0.929 &
          #0.027 TRPSOMG_07_07 + #0.929 #0.027 TRPSOMG_08_07 + #0.929 #0.027 &
          TRPSOMG_09_07 + #0.929 #0.027 TRPSOMG_10_07 + #0.929 #0.027 TRPSOMG_11_07 + &
          #0.929 #0.027 TRPSOMG_12_07 + #0.929 #0.027 TRPSOMG_13_07
 S246.5 TRPSOMG_14_07 + OH = #1.000 #0.028 TRPSOMG_01_01 + #1.000 #0.028 &
          TRPSOMG_02_01 + #1.000 #0.028 TRPSOMG_03_01 + #1.000 #0.028 TRPSOMG_04_01 + &
          #1.000 #0.028 TRPSOMG_05_01 + #1.000 #0.028 TRPSOMG_06_01 + #1.000 &
          #0.028 TRPSOMG_07_01 + #1.000 #0.028 TRPSOMG_08_01 + #1.000 #0.028 &
          TRPSOMG_09_01 + #1.000 #0.028 TRPSOMG_10_01 + #1.000 #0.000 TRPSOMG_11_01 + &
          #1.000 #0.000 TRPSOMG_12_01 + #1.000 #0.000 TRPSOMG_13_01 + #1.000 &
          #0.028 TRPSOMG_01_02 + #1.000 #0.028 TRPSOMG_02_02 + #1.000 #0.028 &
          TRPSOMG_03_02 + #1.000 #0.028 TRPSOMG_04_02 + #1.000 #0.028 TRPSOMG_05_02 + &
          #1.000 #0.028 TRPSOMG_06_02 + #1.000 #0.028 TRPSOMG_07_02 + #1.000 &
          #0.028 TRPSOMG_08_02 + #1.000 #0.028 TRPSOMG_09_02 + #1.000 #0.028 &
          TRPSOMG_10_02 + #1.000 #0.000 TRPSOMG_11_02 + #1.000 #0.000 TRPSOMG_12_02 + &
          #1.000 #0.000 TRPSOMG_13_02 + #1.000 #0.028 TRPSOMG_02_03 + #1.000 &
          #0.028 TRPSOMG_03_03 + #1.000 #0.028 TRPSOMG_04_03 + #1.000 #0.028 &
          TRPSOMG_05_03 + #1.000 #0.028 TRPSOMG_06_03 + #1.000 #0.028 TRPSOMG_07_03 + &
          #1.000 #0.028 TRPSOMG_08_03 + #1.000 #0.028 TRPSOMG_09_03 + #1.000 &
          #0.028 TRPSOMG_10_03 + #1.000 #0.028 TRPSOMG_11_03 + #1.000 #0.000 &
          TRPSOMG_12_03 + #1.000 #0.000 TRPSOMG_13_03 + #1.000 #0.028 TRPSOMG_02_04 + &
          #1.000 #0.028 TRPSOMG_03_04 + #1.000 #0.028 TRPSOMG_04_04 + #1.000 &
          #0.028 TRPSOMG_05_04 + #1.000 #0.028 TRPSOMG_06_04 + #1.000 #0.028 &
          TRPSOMG_07_04 + #1.000 #0.028 TRPSOMG_08_04 + #1.000 #0.028 TRPSOMG_09_04 + &
          #1.000 #0.028 TRPSOMG_10_04 + #1.000 #0.028 TRPSOMG_11_04 + #1.000 &
          #0.000 TRPSOMG_12_04 + #1.000 #0.000 TRPSOMG_13_04 + #1.000 #0.028 &
          TRPSOMG_03_05 + #1.000 #0.028 TRPSOMG_04_05 + #1.000 #0.028 TRPSOMG_05_05 + &
          #1.000 #0.028 TRPSOMG_06_05 + #1.000 #0.028 TRPSOMG_07_05 + #1.000 &
          #0.028 TRPSOMG_08_05 + #1.000 #0.028 TRPSOMG_09_05 + #1.000 #0.028 &
          TRPSOMG_10_05 + #1.000 #0.028 TRPSOMG_11_05 + #1.000 #0.028 TRPSOMG_12_05 + &
          #1.000 #0.000 TRPSOMG_13_05 + #1.000 #0.028 TRPSOMG_03_06 + #1.000 &
          #0.028 TRPSOMG_04_06 + #1.000 #0.028 TRPSOMG_05_06 + #1.000 #0.028 &
          TRPSOMG_06_06 + #1.000 #0.028 TRPSOMG_07_06 + #1.000 #0.028 TRPSOMG_08_06 + &
          #1.000 #0.028 TRPSOMG_09_06 + #1.000 #0.028 TRPSOMG_10_06 + #1.000 &
          #0.028 TRPSOMG_11_06 + #1.000 #0.028 TRPSOMG_12_06 + #1.000 #0.000 &
          TRPSOMG_13_06 + #1.000 #0.028 TRPSOMG_04_07 + #1.000 #0.028 TRPSOMG_05_07 + &
          #1.000 #0.028 TRPSOMG_06_07 + #1.000 #0.028 TRPSOMG_07_07 + #1.000 &
          #0.028 TRPSOMG_08_07 + #1.000 #0.028 TRPSOMG_09_07 + #1.000 #0.028 &
          TRPSOMG_10_07 + #1.000 #0.028 TRPSOMG_11_07 + #1.000 #0.028 TRPSOMG_12_07 + &
          #1.000 #0.028 TRPSOMG_13_07
 S247.5 TRPSOMG_15_01 + OH = #0.210 #0.361 TRPSOMG_15_02 + #0.210 #0.007 &
          TRPSOMG_15_03 + #0.210 #0.606 TRPSOMG_15_04 + #0.210 #0.026 TRPSOMG_15_05 + &
          #0.790 #0.071 TRPSOMG_01_01 + #0.790 #0.071 TRPSOMG_02_01 + #0.790 &
          #0.071 TRPSOMG_03_01 + #0.790 #0.071 TRPSOMG_04_01 + #0.790 #0.071 &
          TRPSOMG_05_01 + #0.790 #0.071 TRPSOMG_06_01 + #0.790 #0.071 TRPSOMG_07_01 + &
          #0.790 #0.071 TRPSOMG_08_01 + #0.790 #0.071 TRPSOMG_09_01 + #0.790 &
          #0.071 TRPSOMG_10_01 + #0.790 #0.071 TRPSOMG_11_01 + #0.790 #0.071 &
          TRPSOMG_12_01 + #0.790 #0.071 TRPSOMG_13_01 + #0.790 #0.071 TRPSOMG_14_01 + &
          #0.790 #0.071 TRPSOMG_01_02 + #0.790 #0.071 TRPSOMG_02_02 + #0.790 &
          #0.071 TRPSOMG_03_02 + #0.790 #0.071 TRPSOMG_04_02 + #0.790 #0.071 &
          TRPSOMG_05_02 + #0.790 #0.071 TRPSOMG_06_02 + #0.790 #0.071 TRPSOMG_07_02 + &
          #0.790 #0.071 TRPSOMG_08_02 + #0.790 #0.071 TRPSOMG_09_02 + #0.790 &
          #0.071 TRPSOMG_10_02 + #0.790 #0.071 TRPSOMG_11_02 + #0.790 #0.071 &
          TRPSOMG_12_02 + #0.790 #0.071 TRPSOMG_13_02 + #0.790 #0.071 TRPSOMG_14_02
 S248.5 TRPSOMG_15_02 + OH = #0.161 #0.361 TRPSOMG_15_03 + #0.161 #0.007 &
          TRPSOMG_15_04 + #0.161 #0.606 TRPSOMG_15_05 + #0.161 #0.026 TRPSOMG_15_06 + &
          #0.839 #0.050 TRPSOMG_01_01 + #0.839 #0.050 TRPSOMG_02_01 + #0.839 &
          #0.050 TRPSOMG_03_01 + #0.839 #0.050 TRPSOMG_04_01 + #0.839 #0.050 &
          TRPSOMG_05_01 + #0.839 #0.050 TRPSOMG_06_01 + #0.839 #0.050 TRPSOMG_07_01 + &
          #0.839 #0.050 TRPSOMG_08_01 + #0.839 #0.050 TRPSOMG_09_01 + #0.839 &
          #0.050 TRPSOMG_10_01 + #0.839 #0.050 TRPSOMG_11_01 + #0.839 #0.050 &
          TRPSOMG_12_01 + #0.839 #0.050 TRPSOMG_13_01 + #0.839 #0.000 TRPSOMG_14_01 + &
          #0.839 #0.050 TRPSOMG_01_02 + #0.839 #0.050 TRPSOMG_02_02 + #0.839 &
          #0.050 TRPSOMG_03_02 + #0.839 #0.050 TRPSOMG_04_02 + #0.839 #0.050 &
          TRPSOMG_05_02 + #0.839 #0.050 TRPSOMG_06_02 + #0.839 #0.050 TRPSOMG_07_02 + &
          #0.839 #0.050 TRPSOMG_08_02 + #0.839 #0.050 TRPSOMG_09_02 + #0.839 &
          #0.050 TRPSOMG_10_02 + #0.839 #0.050 TRPSOMG_11_02 + #0.839 #0.050 &
          TRPSOMG_12_02 + #0.839 #0.050 TRPSOMG_13_02 + #0.839 #0.050 TRPSOMG_14_02 + &
          #0.839 #0.050 TRPSOMG_02_03 + #0.839 #0.050 TRPSOMG_03_03 + #0.839 &
          #0.050 TRPSOMG_04_03 + #0.839 #0.050 TRPSOMG_05_03 + #0.839 #0.050 &
          TRPSOMG_06_03 + #0.839 #0.050 TRPSOMG_07_03 + #0.839 #0.050 TRPSOMG_08_03 + &
          #0.839 #0.050 TRPSOMG_09_03 + #0.839 #0.050 TRPSOMG_10_03 + #0.839 &
          #0.050 TRPSOMG_11_03 + #0.839 #0.050 TRPSOMG_12_03 + #0.839 #0.050 &
          TRPSOMG_13_03 + #0.839 #0.050 TRPSOMG_14_03
 S249.5 TRPSOMG_15_03 + OH = #0.131 #0.361 TRPSOMG_15_04 + #0.131 #0.007 &
          TRPSOMG_15_05 + #0.131 #0.606 TRPSOMG_15_06 + #0.131 #0.026 TRPSOMG_15_07 + &
          #0.869 #0.038 TRPSOMG_01_01 + #0.869 #0.038 TRPSOMG_02_01 + #0.869 &
          #0.038 TRPSOMG_03_01 + #0.869 #0.038 TRPSOMG_04_01 + #0.869 #0.038 &
          TRPSOMG_05_01 + #0.869 #0.038 TRPSOMG_06_01 + #0.869 #0.038 TRPSOMG_07_01 + &
          #0.869 #0.038 TRPSOMG_08_01 + #0.869 #0.038 TRPSOMG_09_01 + #0.869 &
          #0.038 TRPSOMG_10_01 + #0.869 #0.038 TRPSOMG_11_01 + #0.869 #0.038 &
          TRPSOMG_12_01 + #0.869 #0.038 TRPSOMG_13_01 + #0.869 #0.000 TRPSOMG_14_01 + &
          #0.869 #0.038 TRPSOMG_01_02 + #0.869 #0.038 TRPSOMG_02_02 + #0.869 &
          #0.038 TRPSOMG_03_02 + #0.869 #0.038 TRPSOMG_04_02 + #0.869 #0.038 &
          TRPSOMG_05_02 + #0.869 #0.038 TRPSOMG_06_02 + #0.869 #0.038 TRPSOMG_07_02 + &
          #0.869 #0.038 TRPSOMG_08_02 + #0.869 #0.038 TRPSOMG_09_02 + #0.869 &
          #0.038 TRPSOMG_10_02 + #0.869 #0.038 TRPSOMG_11_02 + #0.869 #0.038 &
          TRPSOMG_12_02 + #0.869 #0.038 TRPSOMG_13_02 + #0.869 #0.000 TRPSOMG_14_02 + &
          #0.869 #0.038 TRPSOMG_02_03 + #0.869 #0.038 TRPSOMG_03_03 + #0.869 &
          #0.038 TRPSOMG_04_03 + #0.869 #0.038 TRPSOMG_05_03 + #0.869 #0.038 &
          TRPSOMG_06_03 + #0.869 #0.038 TRPSOMG_07_03 + #0.869 #0.038 TRPSOMG_08_03 + &
          #0.869 #0.038 TRPSOMG_09_03 + #0.869 #0.038 TRPSOMG_10_03 + #0.869 &
          #0.038 TRPSOMG_11_03 + #0.869 #0.038 TRPSOMG_12_03 + #0.869 #0.038 &
          TRPSOMG_13_03 + #0.869 #0.038 TRPSOMG_14_03 + #0.869 #0.038 TRPSOMG_02_04 + &
          #0.869 #0.038 TRPSOMG_03_04 + #0.869 #0.038 TRPSOMG_04_04 + #0.869 &
          #0.038 TRPSOMG_05_04 + #0.869 #0.038 TRPSOMG_06_04 + #0.869 #0.038 &
          TRPSOMG_07_04 + #0.869 #0.038 TRPSOMG_08_04 + #0.869 #0.038 TRPSOMG_09_04 + &
          #0.869 #0.038 TRPSOMG_10_04 + #0.869 #0.038 TRPSOMG_11_04 + #0.869 &
          #0.038 TRPSOMG_12_04 + #0.869 #0.038 TRPSOMG_13_04 + #0.869 #0.038 &
          TRPSOMG_14_04
 S250.5 TRPSOMG_15_04 + OH = #0.109 #0.361 TRPSOMG_15_05 + #0.109 #0.007 &
          TRPSOMG_15_06 + #0.109 #0.632 TRPSOMG_15_07 + #0.891 #0.032 TRPSOMG_01_01 + &
          #0.891 #0.032 TRPSOMG_02_01 + #0.891 #0.032 TRPSOMG_03_01 + #0.891 &
          #0.032 TRPSOMG_04_01 + #0.891 #0.032 TRPSOMG_05_01 + #0.891 #0.032 &
          TRPSOMG_06_01 + #0.891 #0.032 TRPSOMG_07_01 + #0.891 #0.032 TRPSOMG_08_01 + &
          #0.891 #0.032 TRPSOMG_09_01 + #0.891 #0.032 TRPSOMG_10_01 + #0.891 &
          #0.032 TRPSOMG_11_01 + #0.891 #0.032 TRPSOMG_12_01 + #0.891 #0.000 &
          TRPSOMG_13_01 + #0.891 #0.000 TRPSOMG_14_01 + #0.891 #0.032 TRPSOMG_01_02 + &
          #0.891 #0.032 TRPSOMG_02_02 + #0.891 #0.032 TRPSOMG_03_02 + #0.891 &
          #0.032 TRPSOMG_04_02 + #0.891 #0.032 TRPSOMG_05_02 + #0.891 #0.032 &
          TRPSOMG_06_02 + #0.891 #0.032 TRPSOMG_07_02 + #0.891 #0.032 TRPSOMG_08_02 + &
          #0.891 #0.032 TRPSOMG_09_02 + #0.891 #0.032 TRPSOMG_10_02 + #0.891 &
          #0.032 TRPSOMG_11_02 + #0.891 #0.032 TRPSOMG_12_02 + #0.891 #0.032 &
          TRPSOMG_13_02 + #0.891 #0.000 TRPSOMG_14_02 + #0.891 #0.032 TRPSOMG_02_03 + &
          #0.891 #0.032 TRPSOMG_03_03 + #0.891 #0.032 TRPSOMG_04_03 + #0.891 &
          #0.032 TRPSOMG_05_03 + #0.891 #0.032 TRPSOMG_06_03 + #0.891 #0.032 &
          TRPSOMG_07_03 + #0.891 #0.032 TRPSOMG_08_03 + #0.891 #0.032 TRPSOMG_09_03 + &
          #0.891 #0.032 TRPSOMG_10_03 + #0.891 #0.032 TRPSOMG_11_03 + #0.891 &
          #0.032 TRPSOMG_12_03 + #0.891 #0.032 TRPSOMG_13_03 + #0.891 #0.000 &
          TRPSOMG_14_03 + #0.891 #0.032 TRPSOMG_02_04 + #0.891 #0.032 TRPSOMG_03_04 + &
          #0.891 #0.032 TRPSOMG_04_04 + #0.891 #0.032 TRPSOMG_05_04 + #0.891 &
          #0.032 TRPSOMG_06_04 + #0.891 #0.032 TRPSOMG_07_04 + #0.891 #0.032 &
          TRPSOMG_08_04 + #0.891 #0.032 TRPSOMG_09_04 + #0.891 #0.032 TRPSOMG_10_04 + &
          #0.891 #0.032 TRPSOMG_11_04 + #0.891 #0.032 TRPSOMG_12_04 + #0.891 &
          #0.032 TRPSOMG_13_04 + #0.891 #0.032 TRPSOMG_14_04 + #0.891 #0.032 &
          TRPSOMG_03_05 + #0.891 #0.032 TRPSOMG_04_05 + #0.891 #0.032 TRPSOMG_05_05 + &
          #0.891 #0.032 TRPSOMG_06_05 + #0.891 #0.032 TRPSOMG_07_05 + #0.891 &
          #0.032 TRPSOMG_08_05 + #0.891 #0.032 TRPSOMG_09_05 + #0.891 #0.032 &
          TRPSOMG_10_05 + #0.891 #0.032 TRPSOMG_11_05 + #0.891 #0.032 TRPSOMG_12_05 + &
          #0.891 #0.032 TRPSOMG_13_05 + #0.891 #0.032 TRPSOMG_14_05
 S251.5 TRPSOMG_15_05 + OH = #0.091 #0.361 TRPSOMG_15_06 + #0.091 #0.639 &
          TRPSOMG_15_07 + #0.909 #0.028 TRPSOMG_01_01 + #0.909 #0.028 TRPSOMG_02_01 + &
          #0.909 #0.028 TRPSOMG_03_01 + #0.909 #0.028 TRPSOMG_04_01 + #0.909 &
          #0.028 TRPSOMG_05_01 + #0.909 #0.028 TRPSOMG_06_01 + #0.909 #0.028 &
          TRPSOMG_07_01 + #0.909 #0.028 TRPSOMG_08_01 + #0.909 #0.028 TRPSOMG_09_01 + &
          #0.909 #0.028 TRPSOMG_10_01 + #0.909 #0.028 TRPSOMG_11_01 + #0.909 &
          #0.028 TRPSOMG_12_01 + #0.909 #0.000 TRPSOMG_13_01 + #0.909 #0.000 &
          TRPSOMG_14_01 + #0.909 #0.028 TRPSOMG_01_02 + #0.909 #0.028 TRPSOMG_02_02 + &
          #0.909 #0.028 TRPSOMG_03_02 + #0.909 #0.028 TRPSOMG_04_02 + #0.909 &
          #0.028 TRPSOMG_05_02 + #0.909 #0.028 TRPSOMG_06_02 + #0.909 #0.028 &
          TRPSOMG_07_02 + #0.909 #0.028 TRPSOMG_08_02 + #0.909 #0.028 TRPSOMG_09_02 + &
          #0.909 #0.028 TRPSOMG_10_02 + #0.909 #0.028 TRPSOMG_11_02 + #0.909 &
          #0.028 TRPSOMG_12_02 + #0.909 #0.000 TRPSOMG_13_02 + #0.909 #0.000 &
          TRPSOMG_14_02 + #0.909 #0.028 TRPSOMG_02_03 + #0.909 #0.028 TRPSOMG_03_03 + &
          #0.909 #0.028 TRPSOMG_04_03 + #0.909 #0.028 TRPSOMG_05_03 + #0.909 &
          #0.028 TRPSOMG_06_03 + #0.909 #0.028 TRPSOMG_07_03 + #0.909 #0.028 &
          TRPSOMG_08_03 + #0.909 #0.028 TRPSOMG_09_03 + #0.909 #0.028 TRPSOMG_10_03 + &
          #0.909 #0.028 TRPSOMG_11_03 + #0.909 #0.028 TRPSOMG_12_03 + #0.909 &
          #0.028 TRPSOMG_13_03 + #0.909 #0.000 TRPSOMG_14_03 + #0.909 #0.028 &
          TRPSOMG_02_04 + #0.909 #0.028 TRPSOMG_03_04 + #0.909 #0.028 TRPSOMG_04_04 + &
          #0.909 #0.028 TRPSOMG_05_04 + #0.909 #0.028 TRPSOMG_06_04 + #0.909 &
          #0.028 TRPSOMG_07_04 + #0.909 #0.028 TRPSOMG_08_04 + #0.909 #0.028 &
          TRPSOMG_09_04 + #0.909 #0.028 TRPSOMG_10_04 + #0.909 #0.028 TRPSOMG_11_04 + &
          #0.909 #0.028 TRPSOMG_12_04 + #0.909 #0.028 TRPSOMG_13_04 + #0.909 &
          #0.000 TRPSOMG_14_04 + #0.909 #0.028 TRPSOMG_03_05 + #0.909 #0.028 &
          TRPSOMG_04_05 + #0.909 #0.028 TRPSOMG_05_05 + #0.909 #0.028 TRPSOMG_06_05 + &
          #0.909 #0.028 TRPSOMG_07_05 + #0.909 #0.028 TRPSOMG_08_05 + #0.909 &
          #0.028 TRPSOMG_09_05 + #0.909 #0.028 TRPSOMG_10_05 + #0.909 #0.028 &
          TRPSOMG_11_05 + #0.909 #0.028 TRPSOMG_12_05 + #0.909 #0.028 TRPSOMG_13_05 + &
          #0.909 #0.028 TRPSOMG_14_05 + #0.909 #0.028 TRPSOMG_03_06 + #0.909 &
          #0.028 TRPSOMG_04_06 + #0.909 #0.028 TRPSOMG_05_06 + #0.909 #0.028 &
          TRPSOMG_06_06 + #0.909 #0.028 TRPSOMG_07_06 + #0.909 #0.028 TRPSOMG_08_06 + &
          #0.909 #0.028 TRPSOMG_09_06 + #0.909 #0.028 TRPSOMG_10_06 + #0.909 &
          #0.028 TRPSOMG_11_06 + #0.909 #0.028 TRPSOMG_12_06 + #0.909 #0.028 &
          TRPSOMG_13_06 + #0.909 #0.028 TRPSOMG_14_06
 S252.5 TRPSOMG_15_06 + OH = #0.077 #1.000 TRPSOMG_15_07 + #0.923 #0.025 &
          TRPSOMG_01_01 + #0.923 #0.025 TRPSOMG_02_01 + #0.923 #0.025 TRPSOMG_03_01 + &
          #0.923 #0.025 TRPSOMG_04_01 + #0.923 #0.025 TRPSOMG_05_01 + #0.923 &
          #0.025 TRPSOMG_06_01 + #0.923 #0.025 TRPSOMG_07_01 + #0.923 #0.025 &
          TRPSOMG_08_01 + #0.923 #0.025 TRPSOMG_09_01 + #0.923 #0.025 TRPSOMG_10_01 + &
          #0.923 #0.025 TRPSOMG_11_01 + #0.923 #0.000 TRPSOMG_12_01 + #0.923 &
          #0.000 TRPSOMG_13_01 + #0.923 #0.000 TRPSOMG_14_01 + #0.923 #0.025 &
          TRPSOMG_01_02 + #0.923 #0.025 TRPSOMG_02_02 + #0.923 #0.025 TRPSOMG_03_02 + &
          #0.923 #0.025 TRPSOMG_04_02 + #0.923 #0.025 TRPSOMG_05_02 + #0.923 &
          #0.025 TRPSOMG_06_02 + #0.923 #0.025 TRPSOMG_07_02 + #0.923 #0.025 &
          TRPSOMG_08_02 + #0.923 #0.025 TRPSOMG_09_02 + #0.923 #0.025 TRPSOMG_10_02 + &
          #0.923 #0.025 TRPSOMG_11_02 + #0.923 #0.025 TRPSOMG_12_02 + #0.923 &
          #0.000 TRPSOMG_13_02 + #0.923 #0.000 TRPSOMG_14_02 + #0.923 #0.025 &
          TRPSOMG_02_03 + #0.923 #0.025 TRPSOMG_03_03 + #0.923 #0.025 TRPSOMG_04_03 + &
          #0.923 #0.025 TRPSOMG_05_03 + #0.923 #0.025 TRPSOMG_06_03 + #0.923 &
          #0.025 TRPSOMG_07_03 + #0.923 #0.025 TRPSOMG_08_03 + #0.923 #0.025 &
          TRPSOMG_09_03 + #0.923 #0.025 TRPSOMG_10_03 + #0.923 #0.025 TRPSOMG_11_03 + &
          #0.923 #0.025 TRPSOMG_12_03 + #0.923 #0.000 TRPSOMG_13_03 + #0.923 &
          #0.000 TRPSOMG_14_03 + #0.923 #0.025 TRPSOMG_02_04 + #0.923 #0.025 &
          TRPSOMG_03_04 + #0.923 #0.025 TRPSOMG_04_04 + #0.923 #0.025 TRPSOMG_05_04 + &
          #0.923 #0.025 TRPSOMG_06_04 + #0.923 #0.025 TRPSOMG_07_04 + #0.923 &
          #0.025 TRPSOMG_08_04 + #0.923 #0.025 TRPSOMG_09_04 + #0.923 #0.025 &
          TRPSOMG_10_04 + #0.923 #0.025 TRPSOMG_11_04 + #0.923 #0.025 TRPSOMG_12_04 + &
          #0.923 #0.025 TRPSOMG_13_04 + #0.923 #0.000 TRPSOMG_14_04 + #0.923 &
          #0.025 TRPSOMG_03_05 + #0.923 #0.025 TRPSOMG_04_05 + #0.923 #0.025 &
          TRPSOMG_05_05 + #0.923 #0.025 TRPSOMG_06_05 + #0.923 #0.025 TRPSOMG_07_05 + &
          #0.923 #0.025 TRPSOMG_08_05 + #0.923 #0.025 TRPSOMG_09_05 + #0.923 &
          #0.025 TRPSOMG_10_05 + #0.923 #0.025 TRPSOMG_11_05 + #0.923 #0.025 &
          TRPSOMG_12_05 + #0.923 #0.025 TRPSOMG_13_05 + #0.923 #0.000 TRPSOMG_14_05 + &
          #0.923 #0.025 TRPSOMG_03_06 + #0.923 #0.025 TRPSOMG_04_06 + #0.923 &
          #0.025 TRPSOMG_05_06 + #0.923 #0.025 TRPSOMG_06_06 + #0.923 #0.025 &
          TRPSOMG_07_06 + #0.923 #0.025 TRPSOMG_08_06 + #0.923 #0.025 TRPSOMG_09_06 + &
          #0.923 #0.025 TRPSOMG_10_06 + #0.923 #0.025 TRPSOMG_11_06 + #0.923 &
          #0.025 TRPSOMG_12_06 + #0.923 #0.025 TRPSOMG_13_06 + #0.923 #0.025 &
          TRPSOMG_14_06 + #0.923 #0.025 TRPSOMG_04_07 + #0.923 #0.025 TRPSOMG_05_07 + &
          #0.923 #0.025 TRPSOMG_06_07 + #0.923 #0.025 TRPSOMG_07_07 + #0.923 &
          #0.025 TRPSOMG_08_07 + #0.923 #0.025 TRPSOMG_09_07 + #0.923 #0.025 &
          TRPSOMG_10_07 + #0.923 #0.025 TRPSOMG_11_07 + #0.923 #0.025 TRPSOMG_12_07 + &
          #0.923 #0.025 TRPSOMG_13_07 + #0.923 #0.025 TRPSOMG_14_07
 S253.5 TRPSOMG_15_07 + OH = #1.000 #0.026 TRPSOMG_01_01 + #1.000 #0.026 &
          TRPSOMG_02_01 + #1.000 #0.026 TRPSOMG_03_01 + #1.000 #0.026 TRPSOMG_04_01 + &
          #1.000 #0.026 TRPSOMG_05_01 + #1.000 #0.026 TRPSOMG_06_01 + #1.000 &
          #0.026 TRPSOMG_07_01 + #1.000 #0.026 TRPSOMG_08_01 + #1.000 #0.026 &
          TRPSOMG_09_01 + #1.000 #0.026 TRPSOMG_10_01 + #1.000 #0.026 TRPSOMG_11_01 + &
          #1.000 #0.000 TRPSOMG_12_01 + #1.000 #0.000 TRPSOMG_13_01 + #1.000 &
          #0.000 TRPSOMG_14_01 + #1.000 #0.026 TRPSOMG_01_02 + #1.000 #0.026 &
          TRPSOMG_02_02 + #1.000 #0.026 TRPSOMG_03_02 + #1.000 #0.026 TRPSOMG_04_02 + &
          #1.000 #0.026 TRPSOMG_05_02 + #1.000 #0.026 TRPSOMG_06_02 + #1.000 &
          #0.026 TRPSOMG_07_02 + #1.000 #0.026 TRPSOMG_08_02 + #1.000 #0.026 &
          TRPSOMG_09_02 + #1.000 #0.026 TRPSOMG_10_02 + #1.000 #0.026 TRPSOMG_11_02 + &
          #1.000 #0.000 TRPSOMG_12_02 + #1.000 #0.000 TRPSOMG_13_02 + #1.000 &
          #0.000 TRPSOMG_14_02 + #1.000 #0.026 TRPSOMG_02_03 + #1.000 #0.026 &
          TRPSOMG_03_03 + #1.000 #0.026 TRPSOMG_04_03 + #1.000 #0.026 TRPSOMG_05_03 + &
          #1.000 #0.026 TRPSOMG_06_03 + #1.000 #0.026 TRPSOMG_07_03 + #1.000 &
          #0.026 TRPSOMG_08_03 + #1.000 #0.026 TRPSOMG_09_03 + #1.000 #0.026 &
          TRPSOMG_10_03 + #1.000 #0.026 TRPSOMG_11_03 + #1.000 #0.026 TRPSOMG_12_03 + &
          #1.000 #0.000 TRPSOMG_13_03 + #1.000 #0.000 TRPSOMG_14_03 + #1.000 &
          #0.026 TRPSOMG_02_04 + #1.000 #0.026 TRPSOMG_03_04 + #1.000 #0.026 &
          TRPSOMG_04_04 + #1.000 #0.026 TRPSOMG_05_04 + #1.000 #0.026 TRPSOMG_06_04 + &
          #1.000 #0.026 TRPSOMG_07_04 + #1.000 #0.026 TRPSOMG_08_04 + #1.000 &
          #0.026 TRPSOMG_09_04 + #1.000 #0.026 TRPSOMG_10_04 + #1.000 #0.026 &
          TRPSOMG_11_04 + #1.000 #0.026 TRPSOMG_12_04 + #1.000 #0.000 TRPSOMG_13_04 + &
          #1.000 #0.000 TRPSOMG_14_04 + #1.000 #0.026 TRPSOMG_03_05 + #1.000 &
          #0.026 TRPSOMG_04_05 + #1.000 #0.026 TRPSOMG_05_05 + #1.000 #0.026 &
          TRPSOMG_06_05 + #1.000 #0.026 TRPSOMG_07_05 + #1.000 #0.026 TRPSOMG_08_05 + &
          #1.000 #0.026 TRPSOMG_09_05 + #1.000 #0.026 TRPSOMG_10_05 + #1.000 &
          #0.026 TRPSOMG_11_05 + #1.000 #0.026 TRPSOMG_12_05 + #1.000 #0.026 &
          TRPSOMG_13_05 + #1.000 #0.000 TRPSOMG_14_05 + #1.000 #0.026 TRPSOMG_03_06 + &
          #1.000 #0.026 TRPSOMG_04_06 + #1.000 #0.026 TRPSOMG_05_06 + #1.000 &
          #0.026 TRPSOMG_06_06 + #1.000 #0.026 TRPSOMG_07_06 + #1.000 #0.026 &
          TRPSOMG_08_06 + #1.000 #0.026 TRPSOMG_09_06 + #1.000 #0.026 TRPSOMG_10_06 + &
          #1.000 #0.026 TRPSOMG_11_06 + #1.000 #0.026 TRPSOMG_12_06 + #1.000 &
          #0.026 TRPSOMG_13_06 + #1.000 #0.000 TRPSOMG_14_06 + #1.000 #0.026 &
          TRPSOMG_04_07 + #1.000 #0.026 TRPSOMG_05_07 + #1.000 #0.026 TRPSOMG_06_07 + &
          #1.000 #0.026 TRPSOMG_07_07 + #1.000 #0.026 TRPSOMG_08_07 + #1.000 &
          #0.026 TRPSOMG_09_07 + #1.000 #0.026 TRPSOMG_10_07 + #1.000 #0.026 &
          TRPSOMG_11_07 + #1.000 #0.026 TRPSOMG_12_07 + #1.000 #0.026 TRPSOMG_13_07 + &
          #1.000 #0.026 TRPSOMG_14_07
 BL25   SESQTRP + OH = #1.000 #0.361 TRPSOMG_15_01 + #1.000 #0.007 TRPSOMG_15_02 + &
          #1.000 #0.606 TRPSOMG_15_03 + #1.000 #0.026 TRPSOMG_15_04 + SESQTRP
 BL27   IVOC + OH = #0.000 #0.361 IVOSOMG_03_04 + #0.000 #0.007 IVOSOMG_03_05 + &
          #0.000 #0.632 IVOSOMG_03_06 + #1.000 #0.167 IVOSOMG_01_01 + #1.000 &
          #0.167 IVOSOMG_02_01 + #1.000 #0.167 IVOSOMG_01_02 + #1.000 #0.167 &
          IVOSOMG_02_02 + #1.000 #0.167 IVOSOMG_02_03 + #1.000 #0.167 IVOSOMG_02_04 + &
          IVOC
 S254.6 IVOSOMG_02_01 + OH = #0.059 #0.361 IVOSOMG_02_02 + #0.059 #0.007 &
          IVOSOMG_02_03 + #0.059 #0.632 IVOSOMG_02_04 + #0.941 #1.000 IVOSOMG_01_01 + &
          #0.941 #1.000 IVOSOMG_01_02
 S255.6 IVOSOMG_02_02 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #2.000 &
          IVOSOMG_01_02
 S256.6 IVOSOMG_02_03 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_01_02
 S257.6 IVOSOMG_02_04 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_01_02
 S258.6 IVOSOMG_03_01 + OH = #0.091 #0.361 IVOSOMG_03_02 + #0.091 #0.007 &
          IVOSOMG_03_03 + #0.091 #0.606 IVOSOMG_03_04 + #0.091 #0.026 IVOSOMG_03_05 + &
          #0.909 #0.500 IVOSOMG_01_01 + #0.909 #0.500 IVOSOMG_02_01 + #0.909 &
          #0.500 IVOSOMG_01_02 + #0.909 #0.500 IVOSOMG_02_02
 S259.6 IVOSOMG_03_02 + OH = #0.035 #0.361 IVOSOMG_03_03 + #0.035 #0.007 &
          IVOSOMG_03_04 + #0.035 #0.606 IVOSOMG_03_05 + #0.035 #0.026 IVOSOMG_03_06 + &
          #0.965 #0.500 IVOSOMG_01_01 + #0.965 #0.000 IVOSOMG_02_01 + #0.965 &
          #0.500 IVOSOMG_01_02 + #0.965 #0.500 IVOSOMG_02_02 + #0.965 #0.500 &
          IVOSOMG_02_03
 S260.6 IVOSOMG_03_04 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #1.000 IVOSOMG_01_02 + #1.000 #0.000 IVOSOMG_02_02 + &
          #1.000 #0.000 IVOSOMG_02_03 + #1.000 #1.000 IVOSOMG_02_04
 S261.6 IVOSOMG_03_05 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_01_02 + #1.000 #0.000 IVOSOMG_02_02 + &
          #1.000 #0.000 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_02_04
 S262.6 IVOSOMG_03_06 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_01_02 + #1.000 #0.000 IVOSOMG_02_02 + &
          #1.000 #0.000 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_02_04
 S263.6 IVOSOMG_04_01 + OH = #0.114 #0.361 IVOSOMG_04_02 + #0.114 #0.007 &
          IVOSOMG_04_03 + #0.114 #0.606 IVOSOMG_04_04 + #0.114 #0.026 IVOSOMG_04_05 + &
          #0.886 #0.333 IVOSOMG_01_01 + #0.886 #0.333 IVOSOMG_02_01 + #0.886 &
          #0.333 IVOSOMG_03_01 + #0.886 #0.333 IVOSOMG_01_02 + #0.886 #0.333 &
          IVOSOMG_02_02 + #0.886 #0.333 IVOSOMG_03_02
 S264.6 IVOSOMG_04_02 + OH = #0.059 #0.361 IVOSOMG_04_03 + #0.059 #0.007 &
          IVOSOMG_04_04 + #0.059 #0.606 IVOSOMG_04_05 + #0.059 #0.026 IVOSOMG_04_06 + &
          #0.941 #0.286 IVOSOMG_01_01 + #0.941 #0.286 IVOSOMG_02_01 + #0.941 &
          #0.000 IVOSOMG_03_01 + #0.941 #0.286 IVOSOMG_01_02 + #0.941 #0.286 &
          IVOSOMG_02_02 + #0.941 #0.286 IVOSOMG_03_02 + #0.941 #0.286 IVOSOMG_02_03 + &
          #0.941 #0.286 IVOSOMG_03_03
 S265.6 IVOSOMG_04_03 + OH = #0.025 #0.361 IVOSOMG_04_04 + #0.025 #0.007 &
          IVOSOMG_04_05 + #0.025 #0.606 IVOSOMG_04_06 + #0.025 #0.026 IVOSOMG_04_07 + &
          #0.975 #0.250 IVOSOMG_01_01 + #0.975 #0.250 IVOSOMG_02_01 + #0.975 &
          #0.000 IVOSOMG_03_01 + #0.975 #0.250 IVOSOMG_01_02 + #0.975 #0.250 &
          IVOSOMG_02_02 + #0.975 #0.000 IVOSOMG_03_02 + #0.975 #0.250 IVOSOMG_02_03 + &
          #0.975 #0.250 IVOSOMG_03_03 + #0.975 #0.250 IVOSOMG_02_04 + #0.975 &
          #0.250 IVOSOMG_03_04
 S266.6 IVOSOMG_04_04 + OH = #1.000 #0.286 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.286 IVOSOMG_01_02 + &
          #1.000 #0.286 IVOSOMG_02_02 + #1.000 #0.000 IVOSOMG_03_02 + #1.000 &
          #0.286 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.286 &
          IVOSOMG_02_04 + #1.000 #0.286 IVOSOMG_03_04 + #1.000 #0.286 IVOSOMG_03_05
 S267.6 IVOSOMG_04_05 + OH = #1.000 #0.333 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.333 IVOSOMG_01_02 + &
          #1.000 #0.000 IVOSOMG_02_02 + #1.000 #0.000 IVOSOMG_03_02 + #1.000 &
          #0.333 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.333 &
          IVOSOMG_02_04 + #1.000 #0.000 IVOSOMG_03_04 + #1.000 #0.333 IVOSOMG_03_05 + &
          #1.000 #0.333 IVOSOMG_03_06
 S268.6 IVOSOMG_04_06 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.667 IVOSOMG_01_02 + &
          #1.000 #0.000 IVOSOMG_02_02 + #1.000 #0.000 IVOSOMG_03_02 + #1.000 &
          #0.000 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.667 &
          IVOSOMG_02_04 + #1.000 #0.000 IVOSOMG_03_04 + #1.000 #0.000 IVOSOMG_03_05 + &
          #1.000 #0.667 IVOSOMG_03_06
 S269.6 IVOSOMG_04_07 + OH = #1.000 #0.000 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_01_02 + &
          #1.000 #0.000 IVOSOMG_02_02 + #1.000 #0.000 IVOSOMG_03_02 + #1.000 &
          #0.000 IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.000 &
          IVOSOMG_02_04 + #1.000 #0.000 IVOSOMG_03_04 + #1.000 #0.000 IVOSOMG_03_05 + &
          #1.000 #0.000 IVOSOMG_03_06
 S270.6 IVOSOMG_05_01 + OH = #0.131 #0.361 IVOSOMG_05_02 + #0.131 #0.007 &
          IVOSOMG_05_03 + #0.131 #0.606 IVOSOMG_05_04 + #0.131 #0.026 IVOSOMG_05_05 + &
          #0.869 #0.250 IVOSOMG_01_01 + #0.869 #0.250 IVOSOMG_02_01 + #0.869 &
          #0.250 IVOSOMG_03_01 + #0.869 #0.250 IVOSOMG_04_01 + #0.869 #0.250 &
          IVOSOMG_01_02 + #0.869 #0.250 IVOSOMG_02_02 + #0.869 #0.250 IVOSOMG_03_02 + &
          #0.869 #0.250 IVOSOMG_04_02
 S271.6 IVOSOMG_05_02 + OH = #0.077 #0.361 IVOSOMG_05_03 + #0.077 #0.007 &
          IVOSOMG_05_04 + #0.077 #0.606 IVOSOMG_05_05 + #0.077 #0.026 IVOSOMG_05_06 + &
          #0.923 #0.182 IVOSOMG_01_01 + #0.923 #0.182 IVOSOMG_02_01 + #0.923 &
          #0.182 IVOSOMG_03_01 + #0.923 #0.000 IVOSOMG_04_01 + #0.923 #0.182 &
          IVOSOMG_01_02 + #0.923 #0.182 IVOSOMG_02_02 + #0.923 #0.182 IVOSOMG_03_02 + &
          #0.923 #0.182 IVOSOMG_04_02 + #0.923 #0.182 IVOSOMG_02_03 + #0.923 &
          #0.182 IVOSOMG_03_03 + #0.923 #0.182 IVOSOMG_04_03
 S272.6 IVOSOMG_05_03 + OH = #0.043 #0.361 IVOSOMG_05_04 + #0.043 #0.007 &
          IVOSOMG_05_05 + #0.043 #0.606 IVOSOMG_05_06 + #0.043 #0.026 IVOSOMG_05_07 + &
          #0.957 #0.154 IVOSOMG_01_01 + #0.957 #0.154 IVOSOMG_02_01 + #0.957 &
          #0.154 IVOSOMG_03_01 + #0.957 #0.000 IVOSOMG_04_01 + #0.957 #0.154 &
          IVOSOMG_01_02 + #0.957 #0.154 IVOSOMG_02_02 + #0.957 #0.154 IVOSOMG_03_02 + &
          #0.957 #0.000 IVOSOMG_04_02 + #0.957 #0.154 IVOSOMG_02_03 + #0.957 &
          #0.154 IVOSOMG_03_03 + #0.957 #0.154 IVOSOMG_04_03 + #0.957 #0.154 &
          IVOSOMG_02_04 + #0.957 #0.154 IVOSOMG_03_04 + #0.957 #0.154 IVOSOMG_04_04
 S273.6 IVOSOMG_05_04 + OH = #0.019 #0.361 IVOSOMG_05_05 + #0.019 #0.007 &
          IVOSOMG_05_06 + #0.019 #0.632 IVOSOMG_05_07 + #0.981 #0.143 IVOSOMG_01_01 + &
          #0.981 #0.143 IVOSOMG_02_01 + #0.981 #0.000 IVOSOMG_03_01 + #0.981 &
          #0.000 IVOSOMG_04_01 + #0.981 #0.143 IVOSOMG_01_02 + #0.981 #0.143 &
          IVOSOMG_02_02 + #0.981 #0.143 IVOSOMG_03_02 + #0.981 #0.000 IVOSOMG_04_02 + &
          #0.981 #0.143 IVOSOMG_02_03 + #0.981 #0.143 IVOSOMG_03_03 + #0.981 &
          #0.000 IVOSOMG_04_03 + #0.981 #0.143 IVOSOMG_02_04 + #0.981 #0.143 &
          IVOSOMG_03_04 + #0.981 #0.143 IVOSOMG_04_04 + #0.981 #0.143 IVOSOMG_03_05 + &
          #0.981 #0.143 IVOSOMG_04_05
 S274.6 IVOSOMG_05_05 + OH = #1.000 #0.143 IVOSOMG_01_01 + #1.000 #0.143 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.143 IVOSOMG_01_02 + #1.000 #0.143 IVOSOMG_02_02 + #1.000 &
          #0.000 IVOSOMG_03_02 + #1.000 #0.000 IVOSOMG_04_02 + #1.000 #0.143 &
          IVOSOMG_02_03 + #1.000 #0.143 IVOSOMG_03_03 + #1.000 #0.000 IVOSOMG_04_03 + &
          #1.000 #0.143 IVOSOMG_02_04 + #1.000 #0.143 IVOSOMG_03_04 + #1.000 &
          #0.000 IVOSOMG_04_04 + #1.000 #0.143 IVOSOMG_03_05 + #1.000 #0.143 &
          IVOSOMG_04_05 + #1.000 #0.143 IVOSOMG_03_06 + #1.000 #0.143 IVOSOMG_04_06
 S275.6 IVOSOMG_05_06 + OH = #1.000 #0.154 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.154 IVOSOMG_01_02 + #1.000 #0.154 IVOSOMG_02_02 + #1.000 &
          #0.000 IVOSOMG_03_02 + #1.000 #0.000 IVOSOMG_04_02 + #1.000 #0.154 &
          IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.000 IVOSOMG_04_03 + &
          #1.000 #0.154 IVOSOMG_02_04 + #1.000 #0.154 IVOSOMG_03_04 + #1.000 &
          #0.000 IVOSOMG_04_04 + #1.000 #0.154 IVOSOMG_03_05 + #1.000 #0.000 &
          IVOSOMG_04_05 + #1.000 #0.154 IVOSOMG_03_06 + #1.000 #0.154 IVOSOMG_04_06 + &
          #1.000 #0.154 IVOSOMG_04_07
 S276.6 IVOSOMG_05_07 + OH = #1.000 #0.200 IVOSOMG_01_01 + #1.000 #0.000 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.200 IVOSOMG_01_02 + #1.000 #0.000 IVOSOMG_02_02 + #1.000 &
          #0.000 IVOSOMG_03_02 + #1.000 #0.000 IVOSOMG_04_02 + #1.000 #0.200 &
          IVOSOMG_02_03 + #1.000 #0.000 IVOSOMG_03_03 + #1.000 #0.000 IVOSOMG_04_03 + &
          #1.000 #0.200 IVOSOMG_02_04 + #1.000 #0.000 IVOSOMG_03_04 + #1.000 &
          #0.000 IVOSOMG_04_04 + #1.000 #0.200 IVOSOMG_03_05 + #1.000 #0.000 &
          IVOSOMG_04_05 + #1.000 #0.200 IVOSOMG_03_06 + #1.000 #0.000 IVOSOMG_04_06 + &
          #1.000 #0.200 IVOSOMG_04_07
 S277.6 IVOSOMG_06_01 + OH = #0.144 #0.361 IVOSOMG_06_02 + #0.144 #0.007 &
          IVOSOMG_06_03 + #0.144 #0.606 IVOSOMG_06_04 + #0.144 #0.026 IVOSOMG_06_05 + &
          #0.856 #0.200 IVOSOMG_01_01 + #0.856 #0.200 IVOSOMG_02_01 + #0.856 &
          #0.200 IVOSOMG_03_01 + #0.856 #0.200 IVOSOMG_04_01 + #0.856 #0.200 &
          IVOSOMG_05_01 + #0.856 #0.200 IVOSOMG_01_02 + #0.856 #0.200 IVOSOMG_02_02 + &
          #0.856 #0.200 IVOSOMG_03_02 + #0.856 #0.200 IVOSOMG_04_02 + #0.856 &
          #0.200 IVOSOMG_05_02
 S278.6 IVOSOMG_06_02 + OH = #0.091 #0.361 IVOSOMG_06_03 + #0.091 #0.007 &
          IVOSOMG_06_04 + #0.091 #0.606 IVOSOMG_06_05 + #0.091 #0.026 IVOSOMG_06_06 + &
          #0.909 #0.154 IVOSOMG_01_01 + #0.909 #0.154 IVOSOMG_02_01 + #0.909 &
          #0.154 IVOSOMG_03_01 + #0.909 #0.154 IVOSOMG_04_01 + #0.909 #0.000 &
          IVOSOMG_05_01 + #0.909 #0.154 IVOSOMG_01_02 + #0.909 #0.154 IVOSOMG_02_02 + &
          #0.909 #0.154 IVOSOMG_03_02 + #0.909 #0.154 IVOSOMG_04_02 + #0.909 &
          #0.154 IVOSOMG_05_02 + #0.909 #0.154 IVOSOMG_02_03 + #0.909 #0.154 &
          IVOSOMG_03_03 + #0.909 #0.154 IVOSOMG_04_03 + #0.909 #0.154 IVOSOMG_05_03
 S279.6 IVOSOMG_06_03 + OH = #0.059 #0.361 IVOSOMG_06_04 + #0.059 #0.007 &
          IVOSOMG_06_05 + #0.059 #0.606 IVOSOMG_06_06 + #0.059 #0.026 IVOSOMG_06_07 + &
          #0.941 #0.125 IVOSOMG_01_01 + #0.941 #0.125 IVOSOMG_02_01 + #0.941 &
          #0.125 IVOSOMG_03_01 + #0.941 #0.125 IVOSOMG_04_01 + #0.941 #0.000 &
          IVOSOMG_05_01 + #0.941 #0.125 IVOSOMG_01_02 + #0.941 #0.125 IVOSOMG_02_02 + &
          #0.941 #0.125 IVOSOMG_03_02 + #0.941 #0.125 IVOSOMG_04_02 + #0.941 &
          #0.000 IVOSOMG_05_02 + #0.941 #0.125 IVOSOMG_02_03 + #0.941 #0.125 &
          IVOSOMG_03_03 + #0.941 #0.125 IVOSOMG_04_03 + #0.941 #0.125 IVOSOMG_05_03 + &
          #0.941 #0.125 IVOSOMG_02_04 + #0.941 #0.125 IVOSOMG_03_04 + #0.941 &
          #0.125 IVOSOMG_04_04 + #0.941 #0.125 IVOSOMG_05_04
 S280.6 IVOSOMG_06_04 + OH = #0.035 #0.361 IVOSOMG_06_05 + #0.035 #0.007 &
          IVOSOMG_06_06 + #0.035 #0.632 IVOSOMG_06_07 + #0.965 #0.111 IVOSOMG_01_01 + &
          #0.965 #0.111 IVOSOMG_02_01 + #0.965 #0.111 IVOSOMG_03_01 + #0.965 &
          #0.000 IVOSOMG_04_01 + #0.965 #0.000 IVOSOMG_05_01 + #0.965 #0.111 &
          IVOSOMG_01_02 + #0.965 #0.111 IVOSOMG_02_02 + #0.965 #0.111 IVOSOMG_03_02 + &
          #0.965 #0.111 IVOSOMG_04_02 + #0.965 #0.000 IVOSOMG_05_02 + #0.965 &
          #0.111 IVOSOMG_02_03 + #0.965 #0.111 IVOSOMG_03_03 + #0.965 #0.111 &
          IVOSOMG_04_03 + #0.965 #0.000 IVOSOMG_05_03 + #0.965 #0.111 IVOSOMG_02_04 + &
          #0.965 #0.111 IVOSOMG_03_04 + #0.965 #0.111 IVOSOMG_04_04 + #0.965 &
          #0.111 IVOSOMG_05_04 + #0.965 #0.111 IVOSOMG_03_05 + #0.965 #0.111 &
          IVOSOMG_04_05 + #0.965 #0.111 IVOSOMG_05_05
 S281.6 IVOSOMG_06_05 + OH = #0.016 #0.361 IVOSOMG_06_06 + #0.016 #0.639 &
          IVOSOMG_06_07 + #0.984 #0.105 IVOSOMG_01_01 + #0.984 #0.105 IVOSOMG_02_01 + &
          #0.984 #0.105 IVOSOMG_03_01 + #0.984 #0.000 IVOSOMG_04_01 + #0.984 &
          #0.000 IVOSOMG_05_01 + #0.984 #0.105 IVOSOMG_01_02 + #0.984 #0.105 &
          IVOSOMG_02_02 + #0.984 #0.105 IVOSOMG_03_02 + #0.984 #0.000 IVOSOMG_04_02 + &
          #0.984 #0.000 IVOSOMG_05_02 + #0.984 #0.105 IVOSOMG_02_03 + #0.984 &
          #0.105 IVOSOMG_03_03 + #0.984 #0.105 IVOSOMG_04_03 + #0.984 #0.000 &
          IVOSOMG_05_03 + #0.984 #0.105 IVOSOMG_02_04 + #0.984 #0.105 IVOSOMG_03_04 + &
          #0.984 #0.105 IVOSOMG_04_04 + #0.984 #0.000 IVOSOMG_05_04 + #0.984 &
          #0.105 IVOSOMG_03_05 + #0.984 #0.105 IVOSOMG_04_05 + #0.984 #0.105 &
          IVOSOMG_05_05 + #0.984 #0.105 IVOSOMG_03_06 + #0.984 #0.105 IVOSOMG_04_06 + &
          #0.984 #0.105 IVOSOMG_05_06
 S282.6 IVOSOMG_06_06 + OH = #1.000 #0.105 IVOSOMG_01_01 + #1.000 #0.105 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.000 IVOSOMG_05_01 + #1.000 #0.105 IVOSOMG_01_02 + #1.000 &
          #0.105 IVOSOMG_02_02 + #1.000 #0.105 IVOSOMG_03_02 + #1.000 #0.000 &
          IVOSOMG_04_02 + #1.000 #0.000 IVOSOMG_05_02 + #1.000 #0.105 IVOSOMG_02_03 + &
          #1.000 #0.105 IVOSOMG_03_03 + #1.000 #0.000 IVOSOMG_04_03 + #1.000 &
          #0.000 IVOSOMG_05_03 + #1.000 #0.105 IVOSOMG_02_04 + #1.000 #0.105 &
          IVOSOMG_03_04 + #1.000 #0.105 IVOSOMG_04_04 + #1.000 #0.000 IVOSOMG_05_04 + &
          #1.000 #0.105 IVOSOMG_03_05 + #1.000 #0.105 IVOSOMG_04_05 + #1.000 &
          #0.000 IVOSOMG_05_05 + #1.000 #0.105 IVOSOMG_03_06 + #1.000 #0.105 &
          IVOSOMG_04_06 + #1.000 #0.105 IVOSOMG_05_06 + #1.000 #0.105 IVOSOMG_04_07 + &
          #1.000 #0.105 IVOSOMG_05_07
 S283.6 IVOSOMG_06_07 + OH = #1.000 #0.125 IVOSOMG_01_01 + #1.000 #0.125 &
          IVOSOMG_02_01 + #1.000 #0.000 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.000 IVOSOMG_05_01 + #1.000 #0.125 IVOSOMG_01_02 + #1.000 &
          #0.125 IVOSOMG_02_02 + #1.000 #0.000 IVOSOMG_03_02 + #1.000 #0.000 &
          IVOSOMG_04_02 + #1.000 #0.000 IVOSOMG_05_02 + #1.000 #0.125 IVOSOMG_02_03 + &
          #1.000 #0.125 IVOSOMG_03_03 + #1.000 #0.000 IVOSOMG_04_03 + #1.000 &
          #0.000 IVOSOMG_05_03 + #1.000 #0.125 IVOSOMG_02_04 + #1.000 #0.125 &
          IVOSOMG_03_04 + #1.000 #0.000 IVOSOMG_04_04 + #1.000 #0.000 IVOSOMG_05_04 + &
          #1.000 #0.125 IVOSOMG_03_05 + #1.000 #0.125 IVOSOMG_04_05 + #1.000 &
          #0.000 IVOSOMG_05_05 + #1.000 #0.125 IVOSOMG_03_06 + #1.000 #0.125 &
          IVOSOMG_04_06 + #1.000 #0.000 IVOSOMG_05_06 + #1.000 #0.125 IVOSOMG_04_07 + &
          #1.000 #0.125 IVOSOMG_05_07
 S284.6 IVOSOMG_07_01 + OH = #0.156 #0.361 IVOSOMG_07_02 + #0.156 #0.007 &
          IVOSOMG_07_03 + #0.156 #0.606 IVOSOMG_07_04 + #0.156 #0.026 IVOSOMG_07_05 + &
          #0.844 #0.167 IVOSOMG_01_01 + #0.844 #0.167 IVOSOMG_02_01 + #0.844 &
          #0.167 IVOSOMG_03_01 + #0.844 #0.167 IVOSOMG_04_01 + #0.844 #0.167 &
          IVOSOMG_05_01 + #0.844 #0.167 IVOSOMG_06_01 + #0.844 #0.167 IVOSOMG_01_02 + &
          #0.844 #0.167 IVOSOMG_02_02 + #0.844 #0.167 IVOSOMG_03_02 + #0.844 &
          #0.167 IVOSOMG_04_02 + #0.844 #0.167 IVOSOMG_05_02 + #0.844 #0.167 &
          IVOSOMG_06_02
 S285.6 IVOSOMG_07_02 + OH = #0.103 #0.361 IVOSOMG_07_03 + #0.103 #0.007 &
          IVOSOMG_07_04 + #0.103 #0.606 IVOSOMG_07_05 + #0.103 #0.026 IVOSOMG_07_06 + &
          #0.897 #0.125 IVOSOMG_01_01 + #0.897 #0.125 IVOSOMG_02_01 + #0.897 &
          #0.125 IVOSOMG_03_01 + #0.897 #0.125 IVOSOMG_04_01 + #0.897 #0.125 &
          IVOSOMG_05_01 + #0.897 #0.000 IVOSOMG_06_01 + #0.897 #0.125 IVOSOMG_01_02 + &
          #0.897 #0.125 IVOSOMG_02_02 + #0.897 #0.125 IVOSOMG_03_02 + #0.897 &
          #0.125 IVOSOMG_04_02 + #0.897 #0.125 IVOSOMG_05_02 + #0.897 #0.125 &
          IVOSOMG_06_02 + #0.897 #0.125 IVOSOMG_02_03 + #0.897 #0.125 IVOSOMG_03_03 + &
          #0.897 #0.125 IVOSOMG_04_03 + #0.897 #0.125 IVOSOMG_05_03 + #0.897 &
          #0.125 IVOSOMG_06_03
 S286.6 IVOSOMG_07_03 + OH = #0.071 #0.361 IVOSOMG_07_04 + #0.071 #0.007 &
          IVOSOMG_07_05 + #0.071 #0.606 IVOSOMG_07_06 + #0.071 #0.026 IVOSOMG_07_07 + &
          #0.929 #0.100 IVOSOMG_01_01 + #0.929 #0.100 IVOSOMG_02_01 + #0.929 &
          #0.100 IVOSOMG_03_01 + #0.929 #0.100 IVOSOMG_04_01 + #0.929 #0.100 &
          IVOSOMG_05_01 + #0.929 #0.000 IVOSOMG_06_01 + #0.929 #0.100 IVOSOMG_01_02 + &
          #0.929 #0.100 IVOSOMG_02_02 + #0.929 #0.100 IVOSOMG_03_02 + #0.929 &
          #0.100 IVOSOMG_04_02 + #0.929 #0.100 IVOSOMG_05_02 + #0.929 #0.000 &
          IVOSOMG_06_02 + #0.929 #0.100 IVOSOMG_02_03 + #0.929 #0.100 IVOSOMG_03_03 + &
          #0.929 #0.100 IVOSOMG_04_03 + #0.929 #0.100 IVOSOMG_05_03 + #0.929 &
          #0.100 IVOSOMG_06_03 + #0.929 #0.100 IVOSOMG_02_04 + #0.929 #0.100 &
          IVOSOMG_03_04 + #0.929 #0.100 IVOSOMG_04_04 + #0.929 #0.100 IVOSOMG_05_04 + &
          #0.929 #0.100 IVOSOMG_06_04
 S287.6 IVOSOMG_07_04 + OH = #0.048 #0.361 IVOSOMG_07_05 + #0.048 #0.007 &
          IVOSOMG_07_06 + #0.048 #0.632 IVOSOMG_07_07 + #0.952 #0.091 IVOSOMG_01_01 + &
          #0.952 #0.091 IVOSOMG_02_01 + #0.952 #0.091 IVOSOMG_03_01 + #0.952 &
          #0.091 IVOSOMG_04_01 + #0.952 #0.000 IVOSOMG_05_01 + #0.952 #0.000 &
          IVOSOMG_06_01 + #0.952 #0.091 IVOSOMG_01_02 + #0.952 #0.091 IVOSOMG_02_02 + &
          #0.952 #0.091 IVOSOMG_03_02 + #0.952 #0.091 IVOSOMG_04_02 + #0.952 &
          #0.091 IVOSOMG_05_02 + #0.952 #0.000 IVOSOMG_06_02 + #0.952 #0.091 &
          IVOSOMG_02_03 + #0.952 #0.091 IVOSOMG_03_03 + #0.952 #0.091 IVOSOMG_04_03 + &
          #0.952 #0.091 IVOSOMG_05_03 + #0.952 #0.000 IVOSOMG_06_03 + #0.952 &
          #0.091 IVOSOMG_02_04 + #0.952 #0.091 IVOSOMG_03_04 + #0.952 #0.091 &
          IVOSOMG_04_04 + #0.952 #0.091 IVOSOMG_05_04 + #0.952 #0.091 IVOSOMG_06_04 + &
          #0.952 #0.091 IVOSOMG_03_05 + #0.952 #0.091 IVOSOMG_04_05 + #0.952 &
          #0.091 IVOSOMG_05_05 + #0.952 #0.091 IVOSOMG_06_05
 S288.6 IVOSOMG_07_05 + OH = #0.029 #0.361 IVOSOMG_07_06 + #0.029 #0.639 &
          IVOSOMG_07_07 + #0.971 #0.083 IVOSOMG_01_01 + #0.971 #0.083 IVOSOMG_02_01 + &
          #0.971 #0.083 IVOSOMG_03_01 + #0.971 #0.083 IVOSOMG_04_01 + #0.971 &
          #0.000 IVOSOMG_05_01 + #0.971 #0.000 IVOSOMG_06_01 + #0.971 #0.083 &
          IVOSOMG_01_02 + #0.971 #0.083 IVOSOMG_02_02 + #0.971 #0.083 IVOSOMG_03_02 + &
          #0.971 #0.083 IVOSOMG_04_02 + #0.971 #0.000 IVOSOMG_05_02 + #0.971 &
          #0.000 IVOSOMG_06_02 + #0.971 #0.083 IVOSOMG_02_03 + #0.971 #0.083 &
          IVOSOMG_03_03 + #0.971 #0.083 IVOSOMG_04_03 + #0.971 #0.083 IVOSOMG_05_03 + &
          #0.971 #0.000 IVOSOMG_06_03 + #0.971 #0.083 IVOSOMG_02_04 + #0.971 &
          #0.083 IVOSOMG_03_04 + #0.971 #0.083 IVOSOMG_04_04 + #0.971 #0.083 &
          IVOSOMG_05_04 + #0.971 #0.000 IVOSOMG_06_04 + #0.971 #0.083 IVOSOMG_03_05 + &
          #0.971 #0.083 IVOSOMG_04_05 + #0.971 #0.083 IVOSOMG_05_05 + #0.971 &
          #0.083 IVOSOMG_06_05 + #0.971 #0.083 IVOSOMG_03_06 + #0.971 #0.083 &
          IVOSOMG_04_06 + #0.971 #0.083 IVOSOMG_05_06 + #0.971 #0.083 IVOSOMG_06_06
 S289.6 IVOSOMG_07_06 + OH = #0.013 #1.000 IVOSOMG_07_07 + #0.987 #0.080 &
          IVOSOMG_01_01 + #0.987 #0.080 IVOSOMG_02_01 + #0.987 #0.080 IVOSOMG_03_01 + &
          #0.987 #0.000 IVOSOMG_04_01 + #0.987 #0.000 IVOSOMG_05_01 + #0.987 &
          #0.000 IVOSOMG_06_01 + #0.987 #0.080 IVOSOMG_01_02 + #0.987 #0.080 &
          IVOSOMG_02_02 + #0.987 #0.080 IVOSOMG_03_02 + #0.987 #0.080 IVOSOMG_04_02 + &
          #0.987 #0.000 IVOSOMG_05_02 + #0.987 #0.000 IVOSOMG_06_02 + #0.987 &
          #0.080 IVOSOMG_02_03 + #0.987 #0.080 IVOSOMG_03_03 + #0.987 #0.080 &
          IVOSOMG_04_03 + #0.987 #0.000 IVOSOMG_05_03 + #0.987 #0.000 IVOSOMG_06_03 + &
          #0.987 #0.080 IVOSOMG_02_04 + #0.987 #0.080 IVOSOMG_03_04 + #0.987 &
          #0.080 IVOSOMG_04_04 + #0.987 #0.080 IVOSOMG_05_04 + #0.987 #0.000 &
          IVOSOMG_06_04 + #0.987 #0.080 IVOSOMG_03_05 + #0.987 #0.080 IVOSOMG_04_05 + &
          #0.987 #0.080 IVOSOMG_05_05 + #0.987 #0.000 IVOSOMG_06_05 + #0.987 &
          #0.080 IVOSOMG_03_06 + #0.987 #0.080 IVOSOMG_04_06 + #0.987 #0.080 &
          IVOSOMG_05_06 + #0.987 #0.080 IVOSOMG_06_06 + #0.987 #0.080 IVOSOMG_04_07 + &
          #0.987 #0.080 IVOSOMG_05_07 + #0.987 #0.080 IVOSOMG_06_07
 S290.6 IVOSOMG_07_07 + OH = #1.000 #0.091 IVOSOMG_01_01 + #1.000 #0.091 &
          IVOSOMG_02_01 + #1.000 #0.091 IVOSOMG_03_01 + #1.000 #0.000 IVOSOMG_04_01 + &
          #1.000 #0.000 IVOSOMG_05_01 + #1.000 #0.000 IVOSOMG_06_01 + #1.000 &
          #0.091 IVOSOMG_01_02 + #1.000 #0.091 IVOSOMG_02_02 + #1.000 #0.091 &
          IVOSOMG_03_02 + #1.000 #0.000 IVOSOMG_04_02 + #1.000 #0.000 IVOSOMG_05_02 + &
          #1.000 #0.000 IVOSOMG_06_02 + #1.000 #0.091 IVOSOMG_02_03 + #1.000 &
          #0.091 IVOSOMG_03_03 + #1.000 #0.091 IVOSOMG_04_03 + #1.000 #0.000 &
          IVOSOMG_05_03 + #1.000 #0.000 IVOSOMG_06_03 + #1.000 #0.091 IVOSOMG_02_04 + &
          #1.000 #0.091 IVOSOMG_03_04 + #1.000 #0.091 IVOSOMG_04_04 + #1.000 &
          #0.000 IVOSOMG_05_04 + #1.000 #0.000 IVOSOMG_06_04 + #1.000 #0.091 &
          IVOSOMG_03_05 + #1.000 #0.091 IVOSOMG_04_05 + #1.000 #0.091 IVOSOMG_05_05 + &
          #1.000 #0.000 IVOSOMG_06_05 + #1.000 #0.091 IVOSOMG_03_06 + #1.000 &
          #0.091 IVOSOMG_04_06 + #1.000 #0.091 IVOSOMG_05_06 + #1.000 #0.000 &
          IVOSOMG_06_06 + #1.000 #0.091 IVOSOMG_04_07 + #1.000 #0.091 IVOSOMG_05_07 + &
          #1.000 #0.091 IVOSOMG_06_07
 S291.6 IVOSOMG_08_01 + OH = #0.165 #0.361 IVOSOMG_08_02 + #0.165 #0.007 &
          IVOSOMG_08_03 + #0.165 #0.606 IVOSOMG_08_04 + #0.165 #0.026 IVOSOMG_08_05 + &
          #0.835 #0.143 IVOSOMG_01_01 + #0.835 #0.143 IVOSOMG_02_01 + #0.835 &
          #0.143 IVOSOMG_03_01 + #0.835 #0.143 IVOSOMG_04_01 + #0.835 #0.143 &
          IVOSOMG_05_01 + #0.835 #0.143 IVOSOMG_06_01 + #0.835 #0.143 IVOSOMG_07_01 + &
          #0.835 #0.143 IVOSOMG_01_02 + #0.835 #0.143 IVOSOMG_02_02 + #0.835 &
          #0.143 IVOSOMG_03_02 + #0.835 #0.143 IVOSOMG_04_02 + #0.835 #0.143 &
          IVOSOMG_05_02 + #0.835 #0.143 IVOSOMG_06_02 + #0.835 #0.143 IVOSOMG_07_02
 S292.6 IVOSOMG_08_02 + OH = #0.114 #0.361 IVOSOMG_08_03 + #0.114 #0.007 &
          IVOSOMG_08_04 + #0.114 #0.606 IVOSOMG_08_05 + #0.114 #0.026 IVOSOMG_08_06 + &
          #0.886 #0.105 IVOSOMG_01_01 + #0.886 #0.105 IVOSOMG_02_01 + #0.886 &
          #0.105 IVOSOMG_03_01 + #0.886 #0.105 IVOSOMG_04_01 + #0.886 #0.105 &
          IVOSOMG_05_01 + #0.886 #0.105 IVOSOMG_06_01 + #0.886 #0.000 IVOSOMG_07_01 + &
          #0.886 #0.105 IVOSOMG_01_02 + #0.886 #0.105 IVOSOMG_02_02 + #0.886 &
          #0.105 IVOSOMG_03_02 + #0.886 #0.105 IVOSOMG_04_02 + #0.886 #0.105 &
          IVOSOMG_05_02 + #0.886 #0.105 IVOSOMG_06_02 + #0.886 #0.105 IVOSOMG_07_02 + &
          #0.886 #0.105 IVOSOMG_02_03 + #0.886 #0.105 IVOSOMG_03_03 + #0.886 &
          #0.105 IVOSOMG_04_03 + #0.886 #0.105 IVOSOMG_05_03 + #0.886 #0.105 &
          IVOSOMG_06_03 + #0.886 #0.105 IVOSOMG_07_03
 S293.6 IVOSOMG_08_03 + OH = #0.082 #0.361 IVOSOMG_08_04 + #0.082 #0.007 &
          IVOSOMG_08_05 + #0.082 #0.606 IVOSOMG_08_06 + #0.082 #0.026 IVOSOMG_08_07 + &
          #0.918 #0.083 IVOSOMG_01_01 + #0.918 #0.083 IVOSOMG_02_01 + #0.918 &
          #0.083 IVOSOMG_03_01 + #0.918 #0.083 IVOSOMG_04_01 + #0.918 #0.083 &
          IVOSOMG_05_01 + #0.918 #0.083 IVOSOMG_06_01 + #0.918 #0.000 IVOSOMG_07_01 + &
          #0.918 #0.083 IVOSOMG_01_02 + #0.918 #0.083 IVOSOMG_02_02 + #0.918 &
          #0.083 IVOSOMG_03_02 + #0.918 #0.083 IVOSOMG_04_02 + #0.918 #0.083 &
          IVOSOMG_05_02 + #0.918 #0.083 IVOSOMG_06_02 + #0.918 #0.000 IVOSOMG_07_02 + &
          #0.918 #0.083 IVOSOMG_02_03 + #0.918 #0.083 IVOSOMG_03_03 + #0.918 &
          #0.083 IVOSOMG_04_03 + #0.918 #0.083 IVOSOMG_05_03 + #0.918 #0.083 &
          IVOSOMG_06_03 + #0.918 #0.083 IVOSOMG_07_03 + #0.918 #0.083 IVOSOMG_02_04 + &
          #0.918 #0.083 IVOSOMG_03_04 + #0.918 #0.083 IVOSOMG_04_04 + #0.918 &
          #0.083 IVOSOMG_05_04 + #0.918 #0.083 IVOSOMG_06_04 + #0.918 #0.083 &
          IVOSOMG_07_04
 S294.6 IVOSOMG_08_04 + OH = #0.059 #0.361 IVOSOMG_08_05 + #0.059 #0.007 &
          IVOSOMG_08_06 + #0.059 #0.632 IVOSOMG_08_07 + #0.941 #0.074 IVOSOMG_01_01 + &
          #0.941 #0.074 IVOSOMG_02_01 + #0.941 #0.074 IVOSOMG_03_01 + #0.941 &
          #0.074 IVOSOMG_04_01 + #0.941 #0.074 IVOSOMG_05_01 + #0.941 #0.000 &
          IVOSOMG_06_01 + #0.941 #0.000 IVOSOMG_07_01 + #0.941 #0.074 IVOSOMG_01_02 + &
          #0.941 #0.074 IVOSOMG_02_02 + #0.941 #0.074 IVOSOMG_03_02 + #0.941 &
          #0.074 IVOSOMG_04_02 + #0.941 #0.074 IVOSOMG_05_02 + #0.941 #0.074 &
          IVOSOMG_06_02 + #0.941 #0.000 IVOSOMG_07_02 + #0.941 #0.074 IVOSOMG_02_03 + &
          #0.941 #0.074 IVOSOMG_03_03 + #0.941 #0.074 IVOSOMG_04_03 + #0.941 &
          #0.074 IVOSOMG_05_03 + #0.941 #0.074 IVOSOMG_06_03 + #0.941 #0.000 &
          IVOSOMG_07_03 + #0.941 #0.074 IVOSOMG_02_04 + #0.941 #0.074 IVOSOMG_03_04 + &
          #0.941 #0.074 IVOSOMG_04_04 + #0.941 #0.074 IVOSOMG_05_04 + #0.941 &
          #0.074 IVOSOMG_06_04 + #0.941 #0.074 IVOSOMG_07_04 + #0.941 #0.074 &
          IVOSOMG_03_05 + #0.941 #0.074 IVOSOMG_04_05 + #0.941 #0.074 IVOSOMG_05_05 + &
          #0.941 #0.074 IVOSOMG_06_05 + #0.941 #0.074 IVOSOMG_07_05
 S295.6 IVOSOMG_08_05 + OH = #0.040 #0.361 IVOSOMG_08_06 + #0.040 #0.639 &
          IVOSOMG_08_07 + #0.960 #0.067 IVOSOMG_01_01 + #0.960 #0.067 IVOSOMG_02_01 + &
          #0.960 #0.067 IVOSOMG_03_01 + #0.960 #0.067 IVOSOMG_04_01 + #0.960 &
          #0.067 IVOSOMG_05_01 + #0.960 #0.000 IVOSOMG_06_01 + #0.960 #0.000 &
          IVOSOMG_07_01 + #0.960 #0.067 IVOSOMG_01_02 + #0.960 #0.067 IVOSOMG_02_02 + &
          #0.960 #0.067 IVOSOMG_03_02 + #0.960 #0.067 IVOSOMG_04_02 + #0.960 &
          #0.067 IVOSOMG_05_02 + #0.960 #0.000 IVOSOMG_06_02 + #0.960 #0.000 &
          IVOSOMG_07_02 + #0.960 #0.067 IVOSOMG_02_03 + #0.960 #0.067 IVOSOMG_03_03 + &
          #0.960 #0.067 IVOSOMG_04_03 + #0.960 #0.067 IVOSOMG_05_03 + #0.960 &
          #0.067 IVOSOMG_06_03 + #0.960 #0.000 IVOSOMG_07_03 + #0.960 #0.067 &
          IVOSOMG_02_04 + #0.960 #0.067 IVOSOMG_03_04 + #0.960 #0.067 IVOSOMG_04_04 + &
          #0.960 #0.067 IVOSOMG_05_04 + #0.960 #0.067 IVOSOMG_06_04 + #0.960 &
          #0.000 IVOSOMG_07_04 + #0.960 #0.067 IVOSOMG_03_05 + #0.960 #0.067 &
          IVOSOMG_04_05 + #0.960 #0.067 IVOSOMG_05_05 + #0.960 #0.067 IVOSOMG_06_05 + &
          #0.960 #0.067 IVOSOMG_07_05 + #0.960 #0.067 IVOSOMG_03_06 + #0.960 &
          #0.067 IVOSOMG_04_06 + #0.960 #0.067 IVOSOMG_05_06 + #0.960 #0.067 &
          IVOSOMG_06_06 + #0.960 #0.067 IVOSOMG_07_06
 S296.6 IVOSOMG_08_06 + OH = #0.025 #1.000 IVOSOMG_08_07 + #0.975 #0.065 &
          IVOSOMG_01_01 + #0.975 #0.065 IVOSOMG_02_01 + #0.975 #0.065 IVOSOMG_03_01 + &
          #0.975 #0.065 IVOSOMG_04_01 + #0.975 #0.000 IVOSOMG_05_01 + #0.975 &
          #0.000 IVOSOMG_06_01 + #0.975 #0.000 IVOSOMG_07_01 + #0.975 #0.065 &
          IVOSOMG_01_02 + #0.975 #0.065 IVOSOMG_02_02 + #0.975 #0.065 IVOSOMG_03_02 + &
          #0.975 #0.065 IVOSOMG_04_02 + #0.975 #0.065 IVOSOMG_05_02 + #0.975 &
          #0.000 IVOSOMG_06_02 + #0.975 #0.000 IVOSOMG_07_02 + #0.975 #0.065 &
          IVOSOMG_02_03 + #0.975 #0.065 IVOSOMG_03_03 + #0.975 #0.065 IVOSOMG_04_03 + &
          #0.975 #0.065 IVOSOMG_05_03 + #0.975 #0.000 IVOSOMG_06_03 + #0.975 &
          #0.000 IVOSOMG_07_03 + #0.975 #0.065 IVOSOMG_02_04 + #0.975 #0.065 &
          IVOSOMG_03_04 + #0.975 #0.065 IVOSOMG_04_04 + #0.975 #0.065 IVOSOMG_05_04 + &
          #0.975 #0.065 IVOSOMG_06_04 + #0.975 #0.000 IVOSOMG_07_04 + #0.975 &
          #0.065 IVOSOMG_03_05 + #0.975 #0.065 IVOSOMG_04_05 + #0.975 #0.065 &
          IVOSOMG_05_05 + #0.975 #0.065 IVOSOMG_06_05 + #0.975 #0.000 IVOSOMG_07_05 + &
          #0.975 #0.065 IVOSOMG_03_06 + #0.975 #0.065 IVOSOMG_04_06 + #0.975 &
          #0.065 IVOSOMG_05_06 + #0.975 #0.065 IVOSOMG_06_06 + #0.975 #0.065 &
          IVOSOMG_07_06 + #0.975 #0.065 IVOSOMG_04_07 + #0.975 #0.065 IVOSOMG_05_07 + &
          #0.975 #0.065 IVOSOMG_06_07 + #0.975 #0.065 IVOSOMG_07_07
 S297.6 IVOSOMG_08_07 + OH = #1.000 #0.071 IVOSOMG_01_01 + #1.000 #0.071 &
          IVOSOMG_02_01 + #1.000 #0.071 IVOSOMG_03_01 + #1.000 #0.071 IVOSOMG_04_01 + &
          #1.000 #0.000 IVOSOMG_05_01 + #1.000 #0.000 IVOSOMG_06_01 + #1.000 &
          #0.000 IVOSOMG_07_01 + #1.000 #0.071 IVOSOMG_01_02 + #1.000 #0.071 &
          IVOSOMG_02_02 + #1.000 #0.071 IVOSOMG_03_02 + #1.000 #0.071 IVOSOMG_04_02 + &
          #1.000 #0.000 IVOSOMG_05_02 + #1.000 #0.000 IVOSOMG_06_02 + #1.000 &
          #0.000 IVOSOMG_07_02 + #1.000 #0.071 IVOSOMG_02_03 + #1.000 #0.071 &
          IVOSOMG_03_03 + #1.000 #0.071 IVOSOMG_04_03 + #1.000 #0.071 IVOSOMG_05_03 + &
          #1.000 #0.000 IVOSOMG_06_03 + #1.000 #0.000 IVOSOMG_07_03 + #1.000 &
          #0.071 IVOSOMG_02_04 + #1.000 #0.071 IVOSOMG_03_04 + #1.000 #0.071 &
          IVOSOMG_04_04 + #1.000 #0.071 IVOSOMG_05_04 + #1.000 #0.000 IVOSOMG_06_04 + &
          #1.000 #0.000 IVOSOMG_07_04 + #1.000 #0.071 IVOSOMG_03_05 + #1.000 &
          #0.071 IVOSOMG_04_05 + #1.000 #0.071 IVOSOMG_05_05 + #1.000 #0.071 &
          IVOSOMG_06_05 + #1.000 #0.000 IVOSOMG_07_05 + #1.000 #0.071 IVOSOMG_03_06 + &
          #1.000 #0.071 IVOSOMG_04_06 + #1.000 #0.071 IVOSOMG_05_06 + #1.000 &
          #0.071 IVOSOMG_06_06 + #1.000 #0.000 IVOSOMG_07_06 + #1.000 #0.071 &
          IVOSOMG_04_07 + #1.000 #0.071 IVOSOMG_05_07 + #1.000 #0.071 IVOSOMG_06_07 + &
          #1.000 #0.071 IVOSOMG_07_07
 S298.6 IVOSOMG_09_01 + OH = #0.174 #0.361 IVOSOMG_09_02 + #0.174 #0.007 &
          IVOSOMG_09_03 + #0.174 #0.606 IVOSOMG_09_04 + #0.174 #0.026 IVOSOMG_09_05 + &
          #0.826 #0.125 IVOSOMG_01_01 + #0.826 #0.125 IVOSOMG_02_01 + #0.826 &
          #0.125 IVOSOMG_03_01 + #0.826 #0.125 IVOSOMG_04_01 + #0.826 #0.125 &
          IVOSOMG_05_01 + #0.826 #0.125 IVOSOMG_06_01 + #0.826 #0.125 IVOSOMG_07_01 + &
          #0.826 #0.125 IVOSOMG_08_01 + #0.826 #0.125 IVOSOMG_01_02 + #0.826 &
          #0.125 IVOSOMG_02_02 + #0.826 #0.125 IVOSOMG_03_02 + #0.826 #0.125 &
          IVOSOMG_04_02 + #0.826 #0.125 IVOSOMG_05_02 + #0.826 #0.125 IVOSOMG_06_02 + &
          #0.826 #0.125 IVOSOMG_07_02 + #0.826 #0.125 IVOSOMG_08_02
 S299.6 IVOSOMG_09_02 + OH = #0.123 #0.361 IVOSOMG_09_03 + #0.123 #0.007 &
          IVOSOMG_09_04 + #0.123 #0.606 IVOSOMG_09_05 + #0.123 #0.026 IVOSOMG_09_06 + &
          #0.877 #0.087 IVOSOMG_01_01 + #0.877 #0.087 IVOSOMG_02_01 + #0.877 &
          #0.087 IVOSOMG_03_01 + #0.877 #0.087 IVOSOMG_04_01 + #0.877 #0.087 &
          IVOSOMG_05_01 + #0.877 #0.087 IVOSOMG_06_01 + #0.877 #0.087 IVOSOMG_07_01 + &
          #0.877 #0.000 IVOSOMG_08_01 + #0.877 #0.087 IVOSOMG_01_02 + #0.877 &
          #0.087 IVOSOMG_02_02 + #0.877 #0.087 IVOSOMG_03_02 + #0.877 #0.087 &
          IVOSOMG_04_02 + #0.877 #0.087 IVOSOMG_05_02 + #0.877 #0.087 IVOSOMG_06_02 + &
          #0.877 #0.087 IVOSOMG_07_02 + #0.877 #0.087 IVOSOMG_08_02 + #0.877 &
          #0.087 IVOSOMG_02_03 + #0.877 #0.087 IVOSOMG_03_03 + #0.877 #0.087 &
          IVOSOMG_04_03 + #0.877 #0.087 IVOSOMG_05_03 + #0.877 #0.087 IVOSOMG_06_03 + &
          #0.877 #0.087 IVOSOMG_07_03 + #0.877 #0.087 IVOSOMG_08_03
 S300.6 IVOSOMG_09_03 + OH = #0.091 #0.361 IVOSOMG_09_04 + #0.091 #0.007 &
          IVOSOMG_09_05 + #0.091 #0.606 IVOSOMG_09_06 + #0.091 #0.026 IVOSOMG_09_07 + &
          #0.909 #0.069 IVOSOMG_01_01 + #0.909 #0.069 IVOSOMG_02_01 + #0.909 &
          #0.069 IVOSOMG_03_01 + #0.909 #0.069 IVOSOMG_04_01 + #0.909 #0.069 &
          IVOSOMG_05_01 + #0.909 #0.069 IVOSOMG_06_01 + #0.909 #0.069 IVOSOMG_07_01 + &
          #0.909 #0.000 IVOSOMG_08_01 + #0.909 #0.069 IVOSOMG_01_02 + #0.909 &
          #0.069 IVOSOMG_02_02 + #0.909 #0.069 IVOSOMG_03_02 + #0.909 #0.069 &
          IVOSOMG_04_02 + #0.909 #0.069 IVOSOMG_05_02 + #0.909 #0.069 IVOSOMG_06_02 + &
          #0.909 #0.069 IVOSOMG_07_02 + #0.909 #0.000 IVOSOMG_08_02 + #0.909 &
          #0.069 IVOSOMG_02_03 + #0.909 #0.069 IVOSOMG_03_03 + #0.909 #0.069 &
          IVOSOMG_04_03 + #0.909 #0.069 IVOSOMG_05_03 + #0.909 #0.069 IVOSOMG_06_03 + &
          #0.909 #0.069 IVOSOMG_07_03 + #0.909 #0.069 IVOSOMG_08_03 + #0.909 &
          #0.069 IVOSOMG_02_04 + #0.909 #0.069 IVOSOMG_03_04 + #0.909 #0.069 &
          IVOSOMG_04_04 + #0.909 #0.069 IVOSOMG_05_04 + #0.909 #0.069 IVOSOMG_06_04 + &
          #0.909 #0.069 IVOSOMG_07_04 + #0.909 #0.069 IVOSOMG_08_04
 S301.6 IVOSOMG_09_04 + OH = #0.068 #0.361 IVOSOMG_09_05 + #0.068 #0.007 &
          IVOSOMG_09_06 + #0.068 #0.632 IVOSOMG_09_07 + #0.932 #0.059 IVOSOMG_01_01 + &
          #0.932 #0.059 IVOSOMG_02_01 + #0.932 #0.059 IVOSOMG_03_01 + #0.932 &
          #0.059 IVOSOMG_04_01 + #0.932 #0.059 IVOSOMG_05_01 + #0.932 #0.059 &
          IVOSOMG_06_01 + #0.932 #0.000 IVOSOMG_07_01 + #0.932 #0.000 IVOSOMG_08_01 + &
          #0.932 #0.059 IVOSOMG_01_02 + #0.932 #0.059 IVOSOMG_02_02 + #0.932 &
          #0.059 IVOSOMG_03_02 + #0.932 #0.059 IVOSOMG_04_02 + #0.932 #0.059 &
          IVOSOMG_05_02 + #0.932 #0.059 IVOSOMG_06_02 + #0.932 #0.059 IVOSOMG_07_02 + &
          #0.932 #0.000 IVOSOMG_08_02 + #0.932 #0.059 IVOSOMG_02_03 + #0.932 &
          #0.059 IVOSOMG_03_03 + #0.932 #0.059 IVOSOMG_04_03 + #0.932 #0.059 &
          IVOSOMG_05_03 + #0.932 #0.059 IVOSOMG_06_03 + #0.932 #0.059 IVOSOMG_07_03 + &
          #0.932 #0.000 IVOSOMG_08_03 + #0.932 #0.059 IVOSOMG_02_04 + #0.932 &
          #0.059 IVOSOMG_03_04 + #0.932 #0.059 IVOSOMG_04_04 + #0.932 #0.059 &
          IVOSOMG_05_04 + #0.932 #0.059 IVOSOMG_06_04 + #0.932 #0.059 IVOSOMG_07_04 + &
          #0.932 #0.059 IVOSOMG_08_04 + #0.932 #0.059 IVOSOMG_03_05 + #0.932 &
          #0.059 IVOSOMG_04_05 + #0.932 #0.059 IVOSOMG_05_05 + #0.932 #0.059 &
          IVOSOMG_06_05 + #0.932 #0.059 IVOSOMG_07_05 + #0.932 #0.059 IVOSOMG_08_05
 S302.6 IVOSOMG_09_05 + OH = #0.050 #0.361 IVOSOMG_09_06 + #0.050 #0.639 &
          IVOSOMG_09_07 + #0.950 #0.053 IVOSOMG_01_01 + #0.950 #0.053 IVOSOMG_02_01 + &
          #0.950 #0.053 IVOSOMG_03_01 + #0.950 #0.053 IVOSOMG_04_01 + #0.950 &
          #0.053 IVOSOMG_05_01 + #0.950 #0.053 IVOSOMG_06_01 + #0.950 #0.000 &
          IVOSOMG_07_01 + #0.950 #0.000 IVOSOMG_08_01 + #0.950 #0.053 IVOSOMG_01_02 + &
          #0.950 #0.053 IVOSOMG_02_02 + #0.950 #0.053 IVOSOMG_03_02 + #0.950 &
          #0.053 IVOSOMG_04_02 + #0.950 #0.053 IVOSOMG_05_02 + #0.950 #0.053 &
          IVOSOMG_06_02 + #0.950 #0.000 IVOSOMG_07_02 + #0.950 #0.000 IVOSOMG_08_02 + &
          #0.950 #0.053 IVOSOMG_02_03 + #0.950 #0.053 IVOSOMG_03_03 + #0.950 &
          #0.053 IVOSOMG_04_03 + #0.950 #0.053 IVOSOMG_05_03 + #0.950 #0.053 &
          IVOSOMG_06_03 + #0.950 #0.053 IVOSOMG_07_03 + #0.950 #0.000 IVOSOMG_08_03 + &
          #0.950 #0.053 IVOSOMG_02_04 + #0.950 #0.053 IVOSOMG_03_04 + #0.950 &
          #0.053 IVOSOMG_04_04 + #0.950 #0.053 IVOSOMG_05_04 + #0.950 #0.053 &
          IVOSOMG_06_04 + #0.950 #0.053 IVOSOMG_07_04 + #0.950 #0.000 IVOSOMG_08_04 + &
          #0.950 #0.053 IVOSOMG_03_05 + #0.950 #0.053 IVOSOMG_04_05 + #0.950 &
          #0.053 IVOSOMG_05_05 + #0.950 #0.053 IVOSOMG_06_05 + #0.950 #0.053 &
          IVOSOMG_07_05 + #0.950 #0.053 IVOSOMG_08_05 + #0.950 #0.053 IVOSOMG_03_06 + &
          #0.950 #0.053 IVOSOMG_04_06 + #0.950 #0.053 IVOSOMG_05_06 + #0.950 &
          #0.053 IVOSOMG_06_06 + #0.950 #0.053 IVOSOMG_07_06 + #0.950 #0.053 &
          IVOSOMG_08_06
 S303.6 IVOSOMG_09_06 + OH = #0.035 #1.000 IVOSOMG_09_07 + #0.965 #0.049 &
          IVOSOMG_01_01 + #0.965 #0.049 IVOSOMG_02_01 + #0.965 #0.049 IVOSOMG_03_01 + &
          #0.965 #0.049 IVOSOMG_04_01 + #0.965 #0.049 IVOSOMG_05_01 + #0.965 &
          #0.000 IVOSOMG_06_01 + #0.965 #0.000 IVOSOMG_07_01 + #0.965 #0.000 &
          IVOSOMG_08_01 + #0.965 #0.049 IVOSOMG_01_02 + #0.965 #0.049 IVOSOMG_02_02 + &
          #0.965 #0.049 IVOSOMG_03_02 + #0.965 #0.049 IVOSOMG_04_02 + #0.965 &
          #0.049 IVOSOMG_05_02 + #0.965 #0.049 IVOSOMG_06_02 + #0.965 #0.000 &
          IVOSOMG_07_02 + #0.965 #0.000 IVOSOMG_08_02 + #0.965 #0.049 IVOSOMG_02_03 + &
          #0.965 #0.049 IVOSOMG_03_03 + #0.965 #0.049 IVOSOMG_04_03 + #0.965 &
          #0.049 IVOSOMG_05_03 + #0.965 #0.049 IVOSOMG_06_03 + #0.965 #0.000 &
          IVOSOMG_07_03 + #0.965 #0.000 IVOSOMG_08_03 + #0.965 #0.049 IVOSOMG_02_04 + &
          #0.965 #0.049 IVOSOMG_03_04 + #0.965 #0.049 IVOSOMG_04_04 + #0.965 &
          #0.049 IVOSOMG_05_04 + #0.965 #0.049 IVOSOMG_06_04 + #0.965 #0.049 &
          IVOSOMG_07_04 + #0.965 #0.000 IVOSOMG_08_04 + #0.965 #0.049 IVOSOMG_03_05 + &
          #0.965 #0.049 IVOSOMG_04_05 + #0.965 #0.049 IVOSOMG_05_05 + #0.965 &
          #0.049 IVOSOMG_06_05 + #0.965 #0.049 IVOSOMG_07_05 + #0.965 #0.000 &
          IVOSOMG_08_05 + #0.965 #0.049 IVOSOMG_03_06 + #0.965 #0.049 IVOSOMG_04_06 + &
          #0.965 #0.049 IVOSOMG_05_06 + #0.965 #0.049 IVOSOMG_06_06 + #0.965 &
          #0.049 IVOSOMG_07_06 + #0.965 #0.049 IVOSOMG_08_06 + #0.965 #0.049 &
          IVOSOMG_04_07 + #0.965 #0.049 IVOSOMG_05_07 + #0.965 #0.049 IVOSOMG_06_07 + &
          #0.965 #0.049 IVOSOMG_07_07 + #0.965 #0.049 IVOSOMG_08_07
 S304.6 IVOSOMG_09_07 + OH = #1.000 #0.053 IVOSOMG_01_01 + #1.000 #0.053 &
          IVOSOMG_02_01 + #1.000 #0.053 IVOSOMG_03_01 + #1.000 #0.053 IVOSOMG_04_01 + &
          #1.000 #0.053 IVOSOMG_05_01 + #1.000 #0.000 IVOSOMG_06_01 + #1.000 &
          #0.000 IVOSOMG_07_01 + #1.000 #0.000 IVOSOMG_08_01 + #1.000 #0.053 &
          IVOSOMG_01_02 + #1.000 #0.053 IVOSOMG_02_02 + #1.000 #0.053 IVOSOMG_03_02 + &
          #1.000 #0.053 IVOSOMG_04_02 + #1.000 #0.053 IVOSOMG_05_02 + #1.000 &
          #0.000 IVOSOMG_06_02 + #1.000 #0.000 IVOSOMG_07_02 + #1.000 #0.000 &
          IVOSOMG_08_02 + #1.000 #0.053 IVOSOMG_02_03 + #1.000 #0.053 IVOSOMG_03_03 + &
          #1.000 #0.053 IVOSOMG_04_03 + #1.000 #0.053 IVOSOMG_05_03 + #1.000 &
          #0.053 IVOSOMG_06_03 + #1.000 #0.000 IVOSOMG_07_03 + #1.000 #0.000 &
          IVOSOMG_08_03 + #1.000 #0.053 IVOSOMG_02_04 + #1.000 #0.053 IVOSOMG_03_04 + &
          #1.000 #0.053 IVOSOMG_04_04 + #1.000 #0.053 IVOSOMG_05_04 + #1.000 &
          #0.053 IVOSOMG_06_04 + #1.000 #0.000 IVOSOMG_07_04 + #1.000 #0.000 &
          IVOSOMG_08_04 + #1.000 #0.053 IVOSOMG_03_05 + #1.000 #0.053 IVOSOMG_04_05 + &
          #1.000 #0.053 IVOSOMG_05_05 + #1.000 #0.053 IVOSOMG_06_05 + #1.000 &
          #0.053 IVOSOMG_07_05 + #1.000 #0.000 IVOSOMG_08_05 + #1.000 #0.053 &
          IVOSOMG_03_06 + #1.000 #0.053 IVOSOMG_04_06 + #1.000 #0.053 IVOSOMG_05_06 + &
          #1.000 #0.053 IVOSOMG_06_06 + #1.000 #0.053 IVOSOMG_07_06 + #1.000 &
          #0.000 IVOSOMG_08_06 + #1.000 #0.053 IVOSOMG_04_07 + #1.000 #0.053 &
          IVOSOMG_05_07 + #1.000 #0.053 IVOSOMG_06_07 + #1.000 #0.053 IVOSOMG_07_07 + &
          #1.000 #0.053 IVOSOMG_08_07
 S305.6 IVOSOMG_10_01 + OH = #0.182 #0.361 IVOSOMG_10_02 + #0.182 #0.007 &
          IVOSOMG_10_03 + #0.182 #0.606 IVOSOMG_10_04 + #0.182 #0.026 IVOSOMG_10_05 + &
          #0.818 #0.111 IVOSOMG_01_01 + #0.818 #0.111 IVOSOMG_02_01 + #0.818 &
          #0.111 IVOSOMG_03_01 + #0.818 #0.111 IVOSOMG_04_01 + #0.818 #0.111 &
          IVOSOMG_05_01 + #0.818 #0.111 IVOSOMG_06_01 + #0.818 #0.111 IVOSOMG_07_01 + &
          #0.818 #0.111 IVOSOMG_08_01 + #0.818 #0.111 IVOSOMG_09_01 + #0.818 &
          #0.111 IVOSOMG_01_02 + #0.818 #0.111 IVOSOMG_02_02 + #0.818 #0.111 &
          IVOSOMG_03_02 + #0.818 #0.111 IVOSOMG_04_02 + #0.818 #0.111 IVOSOMG_05_02 + &
          #0.818 #0.111 IVOSOMG_06_02 + #0.818 #0.111 IVOSOMG_07_02 + #0.818 &
          #0.111 IVOSOMG_08_02 + #0.818 #0.111 IVOSOMG_09_02
 S306.6 IVOSOMG_10_02 + OH = #0.131 #0.361 IVOSOMG_10_03 + #0.131 #0.007 &
          IVOSOMG_10_04 + #0.131 #0.606 IVOSOMG_10_05 + #0.131 #0.026 IVOSOMG_10_06 + &
          #0.869 #0.080 IVOSOMG_01_01 + #0.869 #0.080 IVOSOMG_02_01 + #0.869 &
          #0.080 IVOSOMG_03_01 + #0.869 #0.080 IVOSOMG_04_01 + #0.869 #0.080 &
          IVOSOMG_05_01 + #0.869 #0.080 IVOSOMG_06_01 + #0.869 #0.080 IVOSOMG_07_01 + &
          #0.869 #0.080 IVOSOMG_08_01 + #0.869 #0.000 IVOSOMG_09_01 + #0.869 &
          #0.080 IVOSOMG_01_02 + #0.869 #0.080 IVOSOMG_02_02 + #0.869 #0.080 &
          IVOSOMG_03_02 + #0.869 #0.080 IVOSOMG_04_02 + #0.869 #0.080 IVOSOMG_05_02 + &
          #0.869 #0.080 IVOSOMG_06_02 + #0.869 #0.080 IVOSOMG_07_02 + #0.869 &
          #0.080 IVOSOMG_08_02 + #0.869 #0.080 IVOSOMG_09_02 + #0.869 #0.080 &
          IVOSOMG_02_03 + #0.869 #0.080 IVOSOMG_03_03 + #0.869 #0.080 IVOSOMG_04_03 + &
          #0.869 #0.080 IVOSOMG_05_03 + #0.869 #0.080 IVOSOMG_06_03 + #0.869 &
          #0.080 IVOSOMG_07_03 + #0.869 #0.080 IVOSOMG_08_03 + #0.869 #0.080 &
          IVOSOMG_09_03
 S307.6 IVOSOMG_10_03 + OH = #0.099 #0.361 IVOSOMG_10_04 + #0.099 #0.007 &
          IVOSOMG_10_05 + #0.099 #0.606 IVOSOMG_10_06 + #0.099 #0.026 IVOSOMG_10_07 + &
          #0.901 #0.062 IVOSOMG_01_01 + #0.901 #0.062 IVOSOMG_02_01 + #0.901 &
          #0.062 IVOSOMG_03_01 + #0.901 #0.062 IVOSOMG_04_01 + #0.901 #0.062 &
          IVOSOMG_05_01 + #0.901 #0.062 IVOSOMG_06_01 + #0.901 #0.062 IVOSOMG_07_01 + &
          #0.901 #0.062 IVOSOMG_08_01 + #0.901 #0.000 IVOSOMG_09_01 + #0.901 &
          #0.062 IVOSOMG_01_02 + #0.901 #0.062 IVOSOMG_02_02 + #0.901 #0.062 &
          IVOSOMG_03_02 + #0.901 #0.062 IVOSOMG_04_02 + #0.901 #0.062 IVOSOMG_05_02 + &
          #0.901 #0.062 IVOSOMG_06_02 + #0.901 #0.062 IVOSOMG_07_02 + #0.901 &
          #0.062 IVOSOMG_08_02 + #0.901 #0.000 IVOSOMG_09_02 + #0.901 #0.062 &
          IVOSOMG_02_03 + #0.901 #0.062 IVOSOMG_03_03 + #0.901 #0.062 IVOSOMG_04_03 + &
          #0.901 #0.062 IVOSOMG_05_03 + #0.901 #0.062 IVOSOMG_06_03 + #0.901 &
          #0.062 IVOSOMG_07_03 + #0.901 #0.062 IVOSOMG_08_03 + #0.901 #0.062 &
          IVOSOMG_09_03 + #0.901 #0.062 IVOSOMG_02_04 + #0.901 #0.062 IVOSOMG_03_04 + &
          #0.901 #0.062 IVOSOMG_04_04 + #0.901 #0.062 IVOSOMG_05_04 + #0.901 &
          #0.062 IVOSOMG_06_04 + #0.901 #0.062 IVOSOMG_07_04 + #0.901 #0.062 &
          IVOSOMG_08_04 + #0.901 #0.062 IVOSOMG_09_04
 S308.6 IVOSOMG_10_04 + OH = #0.077 #0.361 IVOSOMG_10_05 + #0.077 #0.007 &
          IVOSOMG_10_06 + #0.077 #0.632 IVOSOMG_10_07 + #0.923 #0.053 IVOSOMG_01_01 + &
          #0.923 #0.053 IVOSOMG_02_01 + #0.923 #0.053 IVOSOMG_03_01 + #0.923 &
          #0.053 IVOSOMG_04_01 + #0.923 #0.053 IVOSOMG_05_01 + #0.923 #0.053 &
          IVOSOMG_06_01 + #0.923 #0.053 IVOSOMG_07_01 + #0.923 #0.000 IVOSOMG_08_01 + &
          #0.923 #0.000 IVOSOMG_09_01 + #0.923 #0.053 IVOSOMG_01_02 + #0.923 &
          #0.053 IVOSOMG_02_02 + #0.923 #0.053 IVOSOMG_03_02 + #0.923 #0.053 &
          IVOSOMG_04_02 + #0.923 #0.053 IVOSOMG_05_02 + #0.923 #0.053 IVOSOMG_06_02 + &
          #0.923 #0.053 IVOSOMG_07_02 + #0.923 #0.053 IVOSOMG_08_02 + #0.923 &
          #0.000 IVOSOMG_09_02 + #0.923 #0.053 IVOSOMG_02_03 + #0.923 #0.053 &
          IVOSOMG_03_03 + #0.923 #0.053 IVOSOMG_04_03 + #0.923 #0.053 IVOSOMG_05_03 + &
          #0.923 #0.053 IVOSOMG_06_03 + #0.923 #0.053 IVOSOMG_07_03 + #0.923 &
          #0.053 IVOSOMG_08_03 + #0.923 #0.000 IVOSOMG_09_03 + #0.923 #0.053 &
          IVOSOMG_02_04 + #0.923 #0.053 IVOSOMG_03_04 + #0.923 #0.053 IVOSOMG_04_04 + &
          #0.923 #0.053 IVOSOMG_05_04 + #0.923 #0.053 IVOSOMG_06_04 + #0.923 &
          #0.053 IVOSOMG_07_04 + #0.923 #0.053 IVOSOMG_08_04 + #0.923 #0.053 &
          IVOSOMG_09_04 + #0.923 #0.053 IVOSOMG_03_05 + #0.923 #0.053 IVOSOMG_04_05 + &
          #0.923 #0.053 IVOSOMG_05_05 + #0.923 #0.053 IVOSOMG_06_05 + #0.923 &
          #0.053 IVOSOMG_07_05 + #0.923 #0.053 IVOSOMG_08_05 + #0.923 #0.053 &
          IVOSOMG_09_05
 S309.6 IVOSOMG_10_05 + OH = #0.059 #0.361 IVOSOMG_10_06 + #0.059 #0.639 &
          IVOSOMG_10_07 + #0.941 #0.047 IVOSOMG_01_01 + #0.941 #0.047 IVOSOMG_02_01 + &
          #0.941 #0.047 IVOSOMG_03_01 + #0.941 #0.047 IVOSOMG_04_01 + #0.941 &
          #0.047 IVOSOMG_05_01 + #0.941 #0.047 IVOSOMG_06_01 + #0.941 #0.047 &
          IVOSOMG_07_01 + #0.941 #0.000 IVOSOMG_08_01 + #0.941 #0.000 IVOSOMG_09_01 + &
          #0.941 #0.047 IVOSOMG_01_02 + #0.941 #0.047 IVOSOMG_02_02 + #0.941 &
          #0.047 IVOSOMG_03_02 + #0.941 #0.047 IVOSOMG_04_02 + #0.941 #0.047 &
          IVOSOMG_05_02 + #0.941 #0.047 IVOSOMG_06_02 + #0.941 #0.047 IVOSOMG_07_02 + &
          #0.941 #0.000 IVOSOMG_08_02 + #0.941 #0.000 IVOSOMG_09_02 + #0.941 &
          #0.047 IVOSOMG_02_03 + #0.941 #0.047 IVOSOMG_03_03 + #0.941 #0.047 &
          IVOSOMG_04_03 + #0.941 #0.047 IVOSOMG_05_03 + #0.941 #0.047 IVOSOMG_06_03 + &
          #0.941 #0.047 IVOSOMG_07_03 + #0.941 #0.047 IVOSOMG_08_03 + #0.941 &
          #0.000 IVOSOMG_09_03 + #0.941 #0.047 IVOSOMG_02_04 + #0.941 #0.047 &
          IVOSOMG_03_04 + #0.941 #0.047 IVOSOMG_04_04 + #0.941 #0.047 IVOSOMG_05_04 + &
          #0.941 #0.047 IVOSOMG_06_04 + #0.941 #0.047 IVOSOMG_07_04 + #0.941 &
          #0.047 IVOSOMG_08_04 + #0.941 #0.000 IVOSOMG_09_04 + #0.941 #0.047 &
          IVOSOMG_03_05 + #0.941 #0.047 IVOSOMG_04_05 + #0.941 #0.047 IVOSOMG_05_05 + &
          #0.941 #0.047 IVOSOMG_06_05 + #0.941 #0.047 IVOSOMG_07_05 + #0.941 &
          #0.047 IVOSOMG_08_05 + #0.941 #0.047 IVOSOMG_09_05 + #0.941 #0.047 &
          IVOSOMG_03_06 + #0.941 #0.047 IVOSOMG_04_06 + #0.941 #0.047 IVOSOMG_05_06 + &
          #0.941 #0.047 IVOSOMG_06_06 + #0.941 #0.047 IVOSOMG_07_06 + #0.941 &
          #0.047 IVOSOMG_08_06 + #0.941 #0.047 IVOSOMG_09_06
 S310.6 IVOSOMG_10_06 + OH = #0.043 #1.000 IVOSOMG_10_07 + #0.957 #0.043 &
          IVOSOMG_01_01 + #0.957 #0.043 IVOSOMG_02_01 + #0.957 #0.043 IVOSOMG_03_01 + &
          #0.957 #0.043 IVOSOMG_04_01 + #0.957 #0.043 IVOSOMG_05_01 + #0.957 &
          #0.043 IVOSOMG_06_01 + #0.957 #0.000 IVOSOMG_07_01 + #0.957 #0.000 &
          IVOSOMG_08_01 + #0.957 #0.000 IVOSOMG_09_01 + #0.957 #0.043 IVOSOMG_01_02 + &
          #0.957 #0.043 IVOSOMG_02_02 + #0.957 #0.043 IVOSOMG_03_02 + #0.957 &
          #0.043 IVOSOMG_04_02 + #0.957 #0.043 IVOSOMG_05_02 + #0.957 #0.043 &
          IVOSOMG_06_02 + #0.957 #0.043 IVOSOMG_07_02 + #0.957 #0.000 IVOSOMG_08_02 + &
          #0.957 #0.000 IVOSOMG_09_02 + #0.957 #0.043 IVOSOMG_02_03 + #0.957 &
          #0.043 IVOSOMG_03_03 + #0.957 #0.043 IVOSOMG_04_03 + #0.957 #0.043 &
          IVOSOMG_05_03 + #0.957 #0.043 IVOSOMG_06_03 + #0.957 #0.043 IVOSOMG_07_03 + &
          #0.957 #0.000 IVOSOMG_08_03 + #0.957 #0.000 IVOSOMG_09_03 + #0.957 &
          #0.043 IVOSOMG_02_04 + #0.957 #0.043 IVOSOMG_03_04 + #0.957 #0.043 &
          IVOSOMG_04_04 + #0.957 #0.043 IVOSOMG_05_04 + #0.957 #0.043 IVOSOMG_06_04 + &
          #0.957 #0.043 IVOSOMG_07_04 + #0.957 #0.043 IVOSOMG_08_04 + #0.957 &
          #0.000 IVOSOMG_09_04 + #0.957 #0.043 IVOSOMG_03_05 + #0.957 #0.043 &
          IVOSOMG_04_05 + #0.957 #0.043 IVOSOMG_05_05 + #0.957 #0.043 IVOSOMG_06_05 + &
          #0.957 #0.043 IVOSOMG_07_05 + #0.957 #0.043 IVOSOMG_08_05 + #0.957 &
          #0.000 IVOSOMG_09_05 + #0.957 #0.043 IVOSOMG_03_06 + #0.957 #0.043 &
          IVOSOMG_04_06 + #0.957 #0.043 IVOSOMG_05_06 + #0.957 #0.043 IVOSOMG_06_06 + &
          #0.957 #0.043 IVOSOMG_07_06 + #0.957 #0.043 IVOSOMG_08_06 + #0.957 &
          #0.043 IVOSOMG_09_06 + #0.957 #0.043 IVOSOMG_04_07 + #0.957 #0.043 &
          IVOSOMG_05_07 + #0.957 #0.043 IVOSOMG_06_07 + #0.957 #0.043 IVOSOMG_07_07 + &
          #0.957 #0.043 IVOSOMG_08_07 + #0.957 #0.043 IVOSOMG_09_07
 S311.6 IVOSOMG_10_07 + OH = #1.000 #0.045 IVOSOMG_01_01 + #1.000 #0.045 &
          IVOSOMG_02_01 + #1.000 #0.045 IVOSOMG_03_01 + #1.000 #0.045 IVOSOMG_04_01 + &
          #1.000 #0.045 IVOSOMG_05_01 + #1.000 #0.045 IVOSOMG_06_01 + #1.000 &
          #0.000 IVOSOMG_07_01 + #1.000 #0.000 IVOSOMG_08_01 + #1.000 #0.000 &
          IVOSOMG_09_01 + #1.000 #0.045 IVOSOMG_01_02 + #1.000 #0.045 IVOSOMG_02_02 + &
          #1.000 #0.045 IVOSOMG_03_02 + #1.000 #0.045 IVOSOMG_04_02 + #1.000 &
          #0.045 IVOSOMG_05_02 + #1.000 #0.045 IVOSOMG_06_02 + #1.000 #0.000 &
          IVOSOMG_07_02 + #1.000 #0.000 IVOSOMG_08_02 + #1.000 #0.000 IVOSOMG_09_02 + &
          #1.000 #0.045 IVOSOMG_02_03 + #1.000 #0.045 IVOSOMG_03_03 + #1.000 &
          #0.045 IVOSOMG_04_03 + #1.000 #0.045 IVOSOMG_05_03 + #1.000 #0.045 &
          IVOSOMG_06_03 + #1.000 #0.045 IVOSOMG_07_03 + #1.000 #0.000 IVOSOMG_08_03 + &
          #1.000 #0.000 IVOSOMG_09_03 + #1.000 #0.045 IVOSOMG_02_04 + #1.000 &
          #0.045 IVOSOMG_03_04 + #1.000 #0.045 IVOSOMG_04_04 + #1.000 #0.045 &
          IVOSOMG_05_04 + #1.000 #0.045 IVOSOMG_06_04 + #1.000 #0.045 IVOSOMG_07_04 + &
          #1.000 #0.000 IVOSOMG_08_04 + #1.000 #0.000 IVOSOMG_09_04 + #1.000 &
          #0.045 IVOSOMG_03_05 + #1.000 #0.045 IVOSOMG_04_05 + #1.000 #0.045 &
          IVOSOMG_05_05 + #1.000 #0.045 IVOSOMG_06_05 + #1.000 #0.045 IVOSOMG_07_05 + &
          #1.000 #0.045 IVOSOMG_08_05 + #1.000 #0.000 IVOSOMG_09_05 + #1.000 &
          #0.045 IVOSOMG_03_06 + #1.000 #0.045 IVOSOMG_04_06 + #1.000 #0.045 &
          IVOSOMG_05_06 + #1.000 #0.045 IVOSOMG_06_06 + #1.000 #0.045 IVOSOMG_07_06 + &
          #1.000 #0.045 IVOSOMG_08_06 + #1.000 #0.000 IVOSOMG_09_06 + #1.000 &
          #0.045 IVOSOMG_04_07 + #1.000 #0.045 IVOSOMG_05_07 + #1.000 #0.045 &
          IVOSOMG_06_07 + #1.000 #0.045 IVOSOMG_07_07 + #1.000 #0.045 IVOSOMG_08_07 + &
          #1.000 #0.045 IVOSOMG_09_07
 S312.6 IVOSOMG_11_01 + OH = #0.188 #0.361 IVOSOMG_11_02 + #0.188 #0.007 &
          IVOSOMG_11_03 + #0.188 #0.606 IVOSOMG_11_04 + #0.188 #0.026 IVOSOMG_11_05 + &
          #0.812 #0.100 IVOSOMG_01_01 + #0.812 #0.100 IVOSOMG_02_01 + #0.812 &
          #0.100 IVOSOMG_03_01 + #0.812 #0.100 IVOSOMG_04_01 + #0.812 #0.100 &
          IVOSOMG_05_01 + #0.812 #0.100 IVOSOMG_06_01 + #0.812 #0.100 IVOSOMG_07_01 + &
          #0.812 #0.100 IVOSOMG_08_01 + #0.812 #0.100 IVOSOMG_09_01 + #0.812 &
          #0.100 IVOSOMG_10_01 + #0.812 #0.100 IVOSOMG_01_02 + #0.812 #0.100 &
          IVOSOMG_02_02 + #0.812 #0.100 IVOSOMG_03_02 + #0.812 #0.100 IVOSOMG_04_02 + &
          #0.812 #0.100 IVOSOMG_05_02 + #0.812 #0.100 IVOSOMG_06_02 + #0.812 &
          #0.100 IVOSOMG_07_02 + #0.812 #0.100 IVOSOMG_08_02 + #0.812 #0.100 &
          IVOSOMG_09_02 + #0.812 #0.100 IVOSOMG_10_02
 S313.6 IVOSOMG_11_02 + OH = #0.138 #0.361 IVOSOMG_11_03 + #0.138 #0.007 &
          IVOSOMG_11_04 + #0.138 #0.606 IVOSOMG_11_05 + #0.138 #0.026 IVOSOMG_11_06 + &
          #0.862 #0.071 IVOSOMG_01_01 + #0.862 #0.071 IVOSOMG_02_01 + #0.862 &
          #0.071 IVOSOMG_03_01 + #0.862 #0.071 IVOSOMG_04_01 + #0.862 #0.071 &
          IVOSOMG_05_01 + #0.862 #0.071 IVOSOMG_06_01 + #0.862 #0.071 IVOSOMG_07_01 + &
          #0.862 #0.071 IVOSOMG_08_01 + #0.862 #0.071 IVOSOMG_09_01 + #0.862 &
          #0.000 IVOSOMG_10_01 + #0.862 #0.071 IVOSOMG_01_02 + #0.862 #0.071 &
          IVOSOMG_02_02 + #0.862 #0.071 IVOSOMG_03_02 + #0.862 #0.071 IVOSOMG_04_02 + &
          #0.862 #0.071 IVOSOMG_05_02 + #0.862 #0.071 IVOSOMG_06_02 + #0.862 &
          #0.071 IVOSOMG_07_02 + #0.862 #0.071 IVOSOMG_08_02 + #0.862 #0.071 &
          IVOSOMG_09_02 + #0.862 #0.071 IVOSOMG_10_02 + #0.862 #0.071 IVOSOMG_02_03 + &
          #0.862 #0.071 IVOSOMG_03_03 + #0.862 #0.071 IVOSOMG_04_03 + #0.862 &
          #0.071 IVOSOMG_05_03 + #0.862 #0.071 IVOSOMG_06_03 + #0.862 #0.071 &
          IVOSOMG_07_03 + #0.862 #0.071 IVOSOMG_08_03 + #0.862 #0.071 IVOSOMG_09_03 + &
          #0.862 #0.071 IVOSOMG_10_03
 S314.6 IVOSOMG_11_03 + OH = #0.107 #0.361 IVOSOMG_11_04 + #0.107 #0.007 &
          IVOSOMG_11_05 + #0.107 #0.606 IVOSOMG_11_06 + #0.107 #0.026 IVOSOMG_11_07 + &
          #0.893 #0.056 IVOSOMG_01_01 + #0.893 #0.056 IVOSOMG_02_01 + #0.893 &
          #0.056 IVOSOMG_03_01 + #0.893 #0.056 IVOSOMG_04_01 + #0.893 #0.056 &
          IVOSOMG_05_01 + #0.893 #0.056 IVOSOMG_06_01 + #0.893 #0.056 IVOSOMG_07_01 + &
          #0.893 #0.056 IVOSOMG_08_01 + #0.893 #0.056 IVOSOMG_09_01 + #0.893 &
          #0.000 IVOSOMG_10_01 + #0.893 #0.056 IVOSOMG_01_02 + #0.893 #0.056 &
          IVOSOMG_02_02 + #0.893 #0.056 IVOSOMG_03_02 + #0.893 #0.056 IVOSOMG_04_02 + &
          #0.893 #0.056 IVOSOMG_05_02 + #0.893 #0.056 IVOSOMG_06_02 + #0.893 &
          #0.056 IVOSOMG_07_02 + #0.893 #0.056 IVOSOMG_08_02 + #0.893 #0.056 &
          IVOSOMG_09_02 + #0.893 #0.000 IVOSOMG_10_02 + #0.893 #0.056 IVOSOMG_02_03 + &
          #0.893 #0.056 IVOSOMG_03_03 + #0.893 #0.056 IVOSOMG_04_03 + #0.893 &
          #0.056 IVOSOMG_05_03 + #0.893 #0.056 IVOSOMG_06_03 + #0.893 #0.056 &
          IVOSOMG_07_03 + #0.893 #0.056 IVOSOMG_08_03 + #0.893 #0.056 IVOSOMG_09_03 + &
          #0.893 #0.056 IVOSOMG_10_03 + #0.893 #0.056 IVOSOMG_02_04 + #0.893 &
          #0.056 IVOSOMG_03_04 + #0.893 #0.056 IVOSOMG_04_04 + #0.893 #0.056 &
          IVOSOMG_05_04 + #0.893 #0.056 IVOSOMG_06_04 + #0.893 #0.056 IVOSOMG_07_04 + &
          #0.893 #0.056 IVOSOMG_08_04 + #0.893 #0.056 IVOSOMG_09_04 + #0.893 &
          #0.056 IVOSOMG_10_04
 S315.6 IVOSOMG_11_04 + OH = #0.084 #0.361 IVOSOMG_11_05 + #0.084 #0.007 &
          IVOSOMG_11_06 + #0.084 #0.632 IVOSOMG_11_07 + #0.916 #0.048 IVOSOMG_01_01 + &
          #0.916 #0.048 IVOSOMG_02_01 + #0.916 #0.048 IVOSOMG_03_01 + #0.916 &
          #0.048 IVOSOMG_04_01 + #0.916 #0.048 IVOSOMG_05_01 + #0.916 #0.048 &
          IVOSOMG_06_01 + #0.916 #0.048 IVOSOMG_07_01 + #0.916 #0.048 IVOSOMG_08_01 + &
          #0.916 #0.000 IVOSOMG_09_01 + #0.916 #0.000 IVOSOMG_10_01 + #0.916 &
          #0.048 IVOSOMG_01_02 + #0.916 #0.048 IVOSOMG_02_02 + #0.916 #0.048 &
          IVOSOMG_03_02 + #0.916 #0.048 IVOSOMG_04_02 + #0.916 #0.048 IVOSOMG_05_02 + &
          #0.916 #0.048 IVOSOMG_06_02 + #0.916 #0.048 IVOSOMG_07_02 + #0.916 &
          #0.048 IVOSOMG_08_02 + #0.916 #0.048 IVOSOMG_09_02 + #0.916 #0.000 &
          IVOSOMG_10_02 + #0.916 #0.048 IVOSOMG_02_03 + #0.916 #0.048 IVOSOMG_03_03 + &
          #0.916 #0.048 IVOSOMG_04_03 + #0.916 #0.048 IVOSOMG_05_03 + #0.916 &
          #0.048 IVOSOMG_06_03 + #0.916 #0.048 IVOSOMG_07_03 + #0.916 #0.048 &
          IVOSOMG_08_03 + #0.916 #0.048 IVOSOMG_09_03 + #0.916 #0.000 IVOSOMG_10_03 + &
          #0.916 #0.048 IVOSOMG_02_04 + #0.916 #0.048 IVOSOMG_03_04 + #0.916 &
          #0.048 IVOSOMG_04_04 + #0.916 #0.048 IVOSOMG_05_04 + #0.916 #0.048 &
          IVOSOMG_06_04 + #0.916 #0.048 IVOSOMG_07_04 + #0.916 #0.048 IVOSOMG_08_04 + &
          #0.916 #0.048 IVOSOMG_09_04 + #0.916 #0.048 IVOSOMG_10_04 + #0.916 &
          #0.048 IVOSOMG_03_05 + #0.916 #0.048 IVOSOMG_04_05 + #0.916 #0.048 &
          IVOSOMG_05_05 + #0.916 #0.048 IVOSOMG_06_05 + #0.916 #0.048 IVOSOMG_07_05 + &
          #0.916 #0.048 IVOSOMG_08_05 + #0.916 #0.048 IVOSOMG_09_05 + #0.916 &
          #0.048 IVOSOMG_10_05
 S316.6 IVOSOMG_11_05 + OH = #0.066 #0.361 IVOSOMG_11_06 + #0.066 #0.639 &
          IVOSOMG_11_07 + #0.934 #0.042 IVOSOMG_01_01 + #0.934 #0.042 IVOSOMG_02_01 + &
          #0.934 #0.042 IVOSOMG_03_01 + #0.934 #0.042 IVOSOMG_04_01 + #0.934 &
          #0.042 IVOSOMG_05_01 + #0.934 #0.042 IVOSOMG_06_01 + #0.934 #0.042 &
          IVOSOMG_07_01 + #0.934 #0.042 IVOSOMG_08_01 + #0.934 #0.000 IVOSOMG_09_01 + &
          #0.934 #0.000 IVOSOMG_10_01 + #0.934 #0.042 IVOSOMG_01_02 + #0.934 &
          #0.042 IVOSOMG_02_02 + #0.934 #0.042 IVOSOMG_03_02 + #0.934 #0.042 &
          IVOSOMG_04_02 + #0.934 #0.042 IVOSOMG_05_02 + #0.934 #0.042 IVOSOMG_06_02 + &
          #0.934 #0.042 IVOSOMG_07_02 + #0.934 #0.042 IVOSOMG_08_02 + #0.934 &
          #0.000 IVOSOMG_09_02 + #0.934 #0.000 IVOSOMG_10_02 + #0.934 #0.042 &
          IVOSOMG_02_03 + #0.934 #0.042 IVOSOMG_03_03 + #0.934 #0.042 IVOSOMG_04_03 + &
          #0.934 #0.042 IVOSOMG_05_03 + #0.934 #0.042 IVOSOMG_06_03 + #0.934 &
          #0.042 IVOSOMG_07_03 + #0.934 #0.042 IVOSOMG_08_03 + #0.934 #0.042 &
          IVOSOMG_09_03 + #0.934 #0.000 IVOSOMG_10_03 + #0.934 #0.042 IVOSOMG_02_04 + &
          #0.934 #0.042 IVOSOMG_03_04 + #0.934 #0.042 IVOSOMG_04_04 + #0.934 &
          #0.042 IVOSOMG_05_04 + #0.934 #0.042 IVOSOMG_06_04 + #0.934 #0.042 &
          IVOSOMG_07_04 + #0.934 #0.042 IVOSOMG_08_04 + #0.934 #0.042 IVOSOMG_09_04 + &
          #0.934 #0.000 IVOSOMG_10_04 + #0.934 #0.042 IVOSOMG_03_05 + #0.934 &
          #0.042 IVOSOMG_04_05 + #0.934 #0.042 IVOSOMG_05_05 + #0.934 #0.042 &
          IVOSOMG_06_05 + #0.934 #0.042 IVOSOMG_07_05 + #0.934 #0.042 IVOSOMG_08_05 + &
          #0.934 #0.042 IVOSOMG_09_05 + #0.934 #0.042 IVOSOMG_10_05 + #0.934 &
          #0.042 IVOSOMG_03_06 + #0.934 #0.042 IVOSOMG_04_06 + #0.934 #0.042 &
          IVOSOMG_05_06 + #0.934 #0.042 IVOSOMG_06_06 + #0.934 #0.042 IVOSOMG_07_06 + &
          #0.934 #0.042 IVOSOMG_08_06 + #0.934 #0.042 IVOSOMG_09_06 + #0.934 &
          #0.042 IVOSOMG_10_06
 S317.6 IVOSOMG_11_06 + OH = #0.051 #1.000 IVOSOMG_11_07 + #0.949 #0.038 &
          IVOSOMG_01_01 + #0.949 #0.038 IVOSOMG_02_01 + #0.949 #0.038 IVOSOMG_03_01 + &
          #0.949 #0.038 IVOSOMG_04_01 + #0.949 #0.038 IVOSOMG_05_01 + #0.949 &
          #0.038 IVOSOMG_06_01 + #0.949 #0.038 IVOSOMG_07_01 + #0.949 #0.000 &
          IVOSOMG_08_01 + #0.949 #0.000 IVOSOMG_09_01 + #0.949 #0.000 IVOSOMG_10_01 + &
          #0.949 #0.038 IVOSOMG_01_02 + #0.949 #0.038 IVOSOMG_02_02 + #0.949 &
          #0.038 IVOSOMG_03_02 + #0.949 #0.038 IVOSOMG_04_02 + #0.949 #0.038 &
          IVOSOMG_05_02 + #0.949 #0.038 IVOSOMG_06_02 + #0.949 #0.038 IVOSOMG_07_02 + &
          #0.949 #0.038 IVOSOMG_08_02 + #0.949 #0.000 IVOSOMG_09_02 + #0.949 &
          #0.000 IVOSOMG_10_02 + #0.949 #0.038 IVOSOMG_02_03 + #0.949 #0.038 &
          IVOSOMG_03_03 + #0.949 #0.038 IVOSOMG_04_03 + #0.949 #0.038 IVOSOMG_05_03 + &
          #0.949 #0.038 IVOSOMG_06_03 + #0.949 #0.038 IVOSOMG_07_03 + #0.949 &
          #0.038 IVOSOMG_08_03 + #0.949 #0.000 IVOSOMG_09_03 + #0.949 #0.000 &
          IVOSOMG_10_03 + #0.949 #0.038 IVOSOMG_02_04 + #0.949 #0.038 IVOSOMG_03_04 + &
          #0.949 #0.038 IVOSOMG_04_04 + #0.949 #0.038 IVOSOMG_05_04 + #0.949 &
          #0.038 IVOSOMG_06_04 + #0.949 #0.038 IVOSOMG_07_04 + #0.949 #0.038 &
          IVOSOMG_08_04 + #0.949 #0.038 IVOSOMG_09_04 + #0.949 #0.000 IVOSOMG_10_04 + &
          #0.949 #0.038 IVOSOMG_03_05 + #0.949 #0.038 IVOSOMG_04_05 + #0.949 &
          #0.038 IVOSOMG_05_05 + #0.949 #0.038 IVOSOMG_06_05 + #0.949 #0.038 &
          IVOSOMG_07_05 + #0.949 #0.038 IVOSOMG_08_05 + #0.949 #0.038 IVOSOMG_09_05 + &
          #0.949 #0.000 IVOSOMG_10_05 + #0.949 #0.038 IVOSOMG_03_06 + #0.949 &
          #0.038 IVOSOMG_04_06 + #0.949 #0.038 IVOSOMG_05_06 + #0.949 #0.038 &
          IVOSOMG_06_06 + #0.949 #0.038 IVOSOMG_07_06 + #0.949 #0.038 IVOSOMG_08_06 + &
          #0.949 #0.038 IVOSOMG_09_06 + #0.949 #0.038 IVOSOMG_10_06 + #0.949 &
          #0.038 IVOSOMG_04_07 + #0.949 #0.038 IVOSOMG_05_07 + #0.949 #0.038 &
          IVOSOMG_06_07 + #0.949 #0.038 IVOSOMG_07_07 + #0.949 #0.038 IVOSOMG_08_07 + &
          #0.949 #0.038 IVOSOMG_09_07 + #0.949 #0.038 IVOSOMG_10_07
 S318.6 IVOSOMG_11_07 + OH = #1.000 #0.040 IVOSOMG_01_01 + #1.000 #0.040 &
          IVOSOMG_02_01 + #1.000 #0.040 IVOSOMG_03_01 + #1.000 #0.040 IVOSOMG_04_01 + &
          #1.000 #0.040 IVOSOMG_05_01 + #1.000 #0.040 IVOSOMG_06_01 + #1.000 &
          #0.040 IVOSOMG_07_01 + #1.000 #0.000 IVOSOMG_08_01 + #1.000 #0.000 &
          IVOSOMG_09_01 + #1.000 #0.000 IVOSOMG_10_01 + #1.000 #0.040 IVOSOMG_01_02 + &
          #1.000 #0.040 IVOSOMG_02_02 + #1.000 #0.040 IVOSOMG_03_02 + #1.000 &
          #0.040 IVOSOMG_04_02 + #1.000 #0.040 IVOSOMG_05_02 + #1.000 #0.040 &
          IVOSOMG_06_02 + #1.000 #0.040 IVOSOMG_07_02 + #1.000 #0.000 IVOSOMG_08_02 + &
          #1.000 #0.000 IVOSOMG_09_02 + #1.000 #0.000 IVOSOMG_10_02 + #1.000 &
          #0.040 IVOSOMG_02_03 + #1.000 #0.040 IVOSOMG_03_03 + #1.000 #0.040 &
          IVOSOMG_04_03 + #1.000 #0.040 IVOSOMG_05_03 + #1.000 #0.040 IVOSOMG_06_03 + &
          #1.000 #0.040 IVOSOMG_07_03 + #1.000 #0.040 IVOSOMG_08_03 + #1.000 &
          #0.000 IVOSOMG_09_03 + #1.000 #0.000 IVOSOMG_10_03 + #1.000 #0.040 &
          IVOSOMG_02_04 + #1.000 #0.040 IVOSOMG_03_04 + #1.000 #0.040 IVOSOMG_04_04 + &
          #1.000 #0.040 IVOSOMG_05_04 + #1.000 #0.040 IVOSOMG_06_04 + #1.000 &
          #0.040 IVOSOMG_07_04 + #1.000 #0.040 IVOSOMG_08_04 + #1.000 #0.000 &
          IVOSOMG_09_04 + #1.000 #0.000 IVOSOMG_10_04 + #1.000 #0.040 IVOSOMG_03_05 + &
          #1.000 #0.040 IVOSOMG_04_05 + #1.000 #0.040 IVOSOMG_05_05 + #1.000 &
          #0.040 IVOSOMG_06_05 + #1.000 #0.040 IVOSOMG_07_05 + #1.000 #0.040 &
          IVOSOMG_08_05 + #1.000 #0.040 IVOSOMG_09_05 + #1.000 #0.000 IVOSOMG_10_05 + &
          #1.000 #0.040 IVOSOMG_03_06 + #1.000 #0.040 IVOSOMG_04_06 + #1.000 &
          #0.040 IVOSOMG_05_06 + #1.000 #0.040 IVOSOMG_06_06 + #1.000 #0.040 &
          IVOSOMG_07_06 + #1.000 #0.040 IVOSOMG_08_06 + #1.000 #0.040 IVOSOMG_09_06 + &
          #1.000 #0.000 IVOSOMG_10_06 + #1.000 #0.040 IVOSOMG_04_07 + #1.000 &
          #0.040 IVOSOMG_05_07 + #1.000 #0.040 IVOSOMG_06_07 + #1.000 #0.040 &
          IVOSOMG_07_07 + #1.000 #0.040 IVOSOMG_08_07 + #1.000 #0.040 IVOSOMG_09_07 + &
          #1.000 #0.040 IVOSOMG_10_07
 S319.6 IVOSOMG_12_01 + OH = #0.194 #0.361 IVOSOMG_12_02 + #0.194 #0.007 &
          IVOSOMG_12_03 + #0.194 #0.606 IVOSOMG_12_04 + #0.194 #0.026 IVOSOMG_12_05 + &
          #0.806 #0.091 IVOSOMG_01_01 + #0.806 #0.091 IVOSOMG_02_01 + #0.806 &
          #0.091 IVOSOMG_03_01 + #0.806 #0.091 IVOSOMG_04_01 + #0.806 #0.091 &
          IVOSOMG_05_01 + #0.806 #0.091 IVOSOMG_06_01 + #0.806 #0.091 IVOSOMG_07_01 + &
          #0.806 #0.091 IVOSOMG_08_01 + #0.806 #0.091 IVOSOMG_09_01 + #0.806 &
          #0.091 IVOSOMG_10_01 + #0.806 #0.091 IVOSOMG_11_01 + #0.806 #0.091 &
          IVOSOMG_01_02 + #0.806 #0.091 IVOSOMG_02_02 + #0.806 #0.091 IVOSOMG_03_02 + &
          #0.806 #0.091 IVOSOMG_04_02 + #0.806 #0.091 IVOSOMG_05_02 + #0.806 &
          #0.091 IVOSOMG_06_02 + #0.806 #0.091 IVOSOMG_07_02 + #0.806 #0.091 &
          IVOSOMG_08_02 + #0.806 #0.091 IVOSOMG_09_02 + #0.806 #0.091 IVOSOMG_10_02 + &
          #0.806 #0.091 IVOSOMG_11_02
 S320.6 IVOSOMG_12_02 + OH = #0.144 #0.361 IVOSOMG_12_03 + #0.144 #0.007 &
          IVOSOMG_12_04 + #0.144 #0.606 IVOSOMG_12_05 + #0.144 #0.026 IVOSOMG_12_06 + &
          #0.856 #0.065 IVOSOMG_01_01 + #0.856 #0.065 IVOSOMG_02_01 + #0.856 &
          #0.065 IVOSOMG_03_01 + #0.856 #0.065 IVOSOMG_04_01 + #0.856 #0.065 &
          IVOSOMG_05_01 + #0.856 #0.065 IVOSOMG_06_01 + #0.856 #0.065 IVOSOMG_07_01 + &
          #0.856 #0.065 IVOSOMG_08_01 + #0.856 #0.065 IVOSOMG_09_01 + #0.856 &
          #0.065 IVOSOMG_10_01 + #0.856 #0.000 IVOSOMG_11_01 + #0.856 #0.065 &
          IVOSOMG_01_02 + #0.856 #0.065 IVOSOMG_02_02 + #0.856 #0.065 IVOSOMG_03_02 + &
          #0.856 #0.065 IVOSOMG_04_02 + #0.856 #0.065 IVOSOMG_05_02 + #0.856 &
          #0.065 IVOSOMG_06_02 + #0.856 #0.065 IVOSOMG_07_02 + #0.856 #0.065 &
          IVOSOMG_08_02 + #0.856 #0.065 IVOSOMG_09_02 + #0.856 #0.065 IVOSOMG_10_02 + &
          #0.856 #0.065 IVOSOMG_11_02 + #0.856 #0.065 IVOSOMG_02_03 + #0.856 &
          #0.065 IVOSOMG_03_03 + #0.856 #0.065 IVOSOMG_04_03 + #0.856 #0.065 &
          IVOSOMG_05_03 + #0.856 #0.065 IVOSOMG_06_03 + #0.856 #0.065 IVOSOMG_07_03 + &
          #0.856 #0.065 IVOSOMG_08_03 + #0.856 #0.065 IVOSOMG_09_03 + #0.856 &
          #0.065 IVOSOMG_10_03 + #0.856 #0.065 IVOSOMG_11_03
 S321.6 IVOSOMG_12_03 + OH = #0.114 #0.361 IVOSOMG_12_04 + #0.114 #0.007 &
          IVOSOMG_12_05 + #0.114 #0.606 IVOSOMG_12_06 + #0.114 #0.026 IVOSOMG_12_07 + &
          #0.886 #0.050 IVOSOMG_01_01 + #0.886 #0.050 IVOSOMG_02_01 + #0.886 &
          #0.050 IVOSOMG_03_01 + #0.886 #0.050 IVOSOMG_04_01 + #0.886 #0.050 &
          IVOSOMG_05_01 + #0.886 #0.050 IVOSOMG_06_01 + #0.886 #0.050 IVOSOMG_07_01 + &
          #0.886 #0.050 IVOSOMG_08_01 + #0.886 #0.050 IVOSOMG_09_01 + #0.886 &
          #0.050 IVOSOMG_10_01 + #0.886 #0.000 IVOSOMG_11_01 + #0.886 #0.050 &
          IVOSOMG_01_02 + #0.886 #0.050 IVOSOMG_02_02 + #0.886 #0.050 IVOSOMG_03_02 + &
          #0.886 #0.050 IVOSOMG_04_02 + #0.886 #0.050 IVOSOMG_05_02 + #0.886 &
          #0.050 IVOSOMG_06_02 + #0.886 #0.050 IVOSOMG_07_02 + #0.886 #0.050 &
          IVOSOMG_08_02 + #0.886 #0.050 IVOSOMG_09_02 + #0.886 #0.050 IVOSOMG_10_02 + &
          #0.886 #0.000 IVOSOMG_11_02 + #0.886 #0.050 IVOSOMG_02_03 + #0.886 &
          #0.050 IVOSOMG_03_03 + #0.886 #0.050 IVOSOMG_04_03 + #0.886 #0.050 &
          IVOSOMG_05_03 + #0.886 #0.050 IVOSOMG_06_03 + #0.886 #0.050 IVOSOMG_07_03 + &
          #0.886 #0.050 IVOSOMG_08_03 + #0.886 #0.050 IVOSOMG_09_03 + #0.886 &
          #0.050 IVOSOMG_10_03 + #0.886 #0.050 IVOSOMG_11_03 + #0.886 #0.050 &
          IVOSOMG_02_04 + #0.886 #0.050 IVOSOMG_03_04 + #0.886 #0.050 IVOSOMG_04_04 + &
          #0.886 #0.050 IVOSOMG_05_04 + #0.886 #0.050 IVOSOMG_06_04 + #0.886 &
          #0.050 IVOSOMG_07_04 + #0.886 #0.050 IVOSOMG_08_04 + #0.886 #0.050 &
          IVOSOMG_09_04 + #0.886 #0.050 IVOSOMG_10_04 + #0.886 #0.050 IVOSOMG_11_04
 S322.6 IVOSOMG_12_04 + OH = #0.091 #0.361 IVOSOMG_12_05 + #0.091 #0.007 &
          IVOSOMG_12_06 + #0.091 #0.632 IVOSOMG_12_07 + #0.909 #0.043 IVOSOMG_01_01 + &
          #0.909 #0.043 IVOSOMG_02_01 + #0.909 #0.043 IVOSOMG_03_01 + #0.909 &
          #0.043 IVOSOMG_04_01 + #0.909 #0.043 IVOSOMG_05_01 + #0.909 #0.043 &
          IVOSOMG_06_01 + #0.909 #0.043 IVOSOMG_07_01 + #0.909 #0.043 IVOSOMG_08_01 + &
          #0.909 #0.043 IVOSOMG_09_01 + #0.909 #0.000 IVOSOMG_10_01 + #0.909 &
          #0.000 IVOSOMG_11_01 + #0.909 #0.043 IVOSOMG_01_02 + #0.909 #0.043 &
          IVOSOMG_02_02 + #0.909 #0.043 IVOSOMG_03_02 + #0.909 #0.043 IVOSOMG_04_02 + &
          #0.909 #0.043 IVOSOMG_05_02 + #0.909 #0.043 IVOSOMG_06_02 + #0.909 &
          #0.043 IVOSOMG_07_02 + #0.909 #0.043 IVOSOMG_08_02 + #0.909 #0.043 &
          IVOSOMG_09_02 + #0.909 #0.043 IVOSOMG_10_02 + #0.909 #0.000 IVOSOMG_11_02 + &
          #0.909 #0.043 IVOSOMG_02_03 + #0.909 #0.043 IVOSOMG_03_03 + #0.909 &
          #0.043 IVOSOMG_04_03 + #0.909 #0.043 IVOSOMG_05_03 + #0.909 #0.043 &
          IVOSOMG_06_03 + #0.909 #0.043 IVOSOMG_07_03 + #0.909 #0.043 IVOSOMG_08_03 + &
          #0.909 #0.043 IVOSOMG_09_03 + #0.909 #0.043 IVOSOMG_10_03 + #0.909 &
          #0.000 IVOSOMG_11_03 + #0.909 #0.043 IVOSOMG_02_04 + #0.909 #0.043 &
          IVOSOMG_03_04 + #0.909 #0.043 IVOSOMG_04_04 + #0.909 #0.043 IVOSOMG_05_04 + &
          #0.909 #0.043 IVOSOMG_06_04 + #0.909 #0.043 IVOSOMG_07_04 + #0.909 &
          #0.043 IVOSOMG_08_04 + #0.909 #0.043 IVOSOMG_09_04 + #0.909 #0.043 &
          IVOSOMG_10_04 + #0.909 #0.043 IVOSOMG_11_04 + #0.909 #0.043 IVOSOMG_03_05 + &
          #0.909 #0.043 IVOSOMG_04_05 + #0.909 #0.043 IVOSOMG_05_05 + #0.909 &
          #0.043 IVOSOMG_06_05 + #0.909 #0.043 IVOSOMG_07_05 + #0.909 #0.043 &
          IVOSOMG_08_05 + #0.909 #0.043 IVOSOMG_09_05 + #0.909 #0.043 IVOSOMG_10_05 + &
          #0.909 #0.043 IVOSOMG_11_05
 S323.6 IVOSOMG_12_05 + OH = #0.073 #0.361 IVOSOMG_12_06 + #0.073 #0.639 &
          IVOSOMG_12_07 + #0.927 #0.037 IVOSOMG_01_01 + #0.927 #0.037 IVOSOMG_02_01 + &
          #0.927 #0.037 IVOSOMG_03_01 + #0.927 #0.037 IVOSOMG_04_01 + #0.927 &
          #0.037 IVOSOMG_05_01 + #0.927 #0.037 IVOSOMG_06_01 + #0.927 #0.037 &
          IVOSOMG_07_01 + #0.927 #0.037 IVOSOMG_08_01 + #0.927 #0.037 IVOSOMG_09_01 + &
          #0.927 #0.000 IVOSOMG_10_01 + #0.927 #0.000 IVOSOMG_11_01 + #0.927 &
          #0.037 IVOSOMG_01_02 + #0.927 #0.037 IVOSOMG_02_02 + #0.927 #0.037 &
          IVOSOMG_03_02 + #0.927 #0.037 IVOSOMG_04_02 + #0.927 #0.037 IVOSOMG_05_02 + &
          #0.927 #0.037 IVOSOMG_06_02 + #0.927 #0.037 IVOSOMG_07_02 + #0.927 &
          #0.037 IVOSOMG_08_02 + #0.927 #0.037 IVOSOMG_09_02 + #0.927 #0.000 &
          IVOSOMG_10_02 + #0.927 #0.000 IVOSOMG_11_02 + #0.927 #0.037 IVOSOMG_02_03 + &
          #0.927 #0.037 IVOSOMG_03_03 + #0.927 #0.037 IVOSOMG_04_03 + #0.927 &
          #0.037 IVOSOMG_05_03 + #0.927 #0.037 IVOSOMG_06_03 + #0.927 #0.037 &
          IVOSOMG_07_03 + #0.927 #0.037 IVOSOMG_08_03 + #0.927 #0.037 IVOSOMG_09_03 + &
          #0.927 #0.037 IVOSOMG_10_03 + #0.927 #0.000 IVOSOMG_11_03 + #0.927 &
          #0.037 IVOSOMG_02_04 + #0.927 #0.037 IVOSOMG_03_04 + #0.927 #0.037 &
          IVOSOMG_04_04 + #0.927 #0.037 IVOSOMG_05_04 + #0.927 #0.037 IVOSOMG_06_04 + &
          #0.927 #0.037 IVOSOMG_07_04 + #0.927 #0.037 IVOSOMG_08_04 + #0.927 &
          #0.037 IVOSOMG_09_04 + #0.927 #0.037 IVOSOMG_10_04 + #0.927 #0.000 &
          IVOSOMG_11_04 + #0.927 #0.037 IVOSOMG_03_05 + #0.927 #0.037 IVOSOMG_04_05 + &
          #0.927 #0.037 IVOSOMG_05_05 + #0.927 #0.037 IVOSOMG_06_05 + #0.927 &
          #0.037 IVOSOMG_07_05 + #0.927 #0.037 IVOSOMG_08_05 + #0.927 #0.037 &
          IVOSOMG_09_05 + #0.927 #0.037 IVOSOMG_10_05 + #0.927 #0.037 IVOSOMG_11_05 + &
          #0.927 #0.037 IVOSOMG_03_06 + #0.927 #0.037 IVOSOMG_04_06 + #0.927 &
          #0.037 IVOSOMG_05_06 + #0.927 #0.037 IVOSOMG_06_06 + #0.927 #0.037 &
          IVOSOMG_07_06 + #0.927 #0.037 IVOSOMG_08_06 + #0.927 #0.037 IVOSOMG_09_06 + &
          #0.927 #0.037 IVOSOMG_10_06 + #0.927 #0.037 IVOSOMG_11_06
 S324.6 IVOSOMG_12_06 + OH = #0.059 #1.000 IVOSOMG_12_07 + #0.941 #0.034 &
          IVOSOMG_01_01 + #0.941 #0.034 IVOSOMG_02_01 + #0.941 #0.034 IVOSOMG_03_01 + &
          #0.941 #0.034 IVOSOMG_04_01 + #0.941 #0.034 IVOSOMG_05_01 + #0.941 &
          #0.034 IVOSOMG_06_01 + #0.941 #0.034 IVOSOMG_07_01 + #0.941 #0.034 &
          IVOSOMG_08_01 + #0.941 #0.000 IVOSOMG_09_01 + #0.941 #0.000 IVOSOMG_10_01 + &
          #0.941 #0.000 IVOSOMG_11_01 + #0.941 #0.034 IVOSOMG_01_02 + #0.941 &
          #0.034 IVOSOMG_02_02 + #0.941 #0.034 IVOSOMG_03_02 + #0.941 #0.034 &
          IVOSOMG_04_02 + #0.941 #0.034 IVOSOMG_05_02 + #0.941 #0.034 IVOSOMG_06_02 + &
          #0.941 #0.034 IVOSOMG_07_02 + #0.941 #0.034 IVOSOMG_08_02 + #0.941 &
          #0.034 IVOSOMG_09_02 + #0.941 #0.000 IVOSOMG_10_02 + #0.941 #0.000 &
          IVOSOMG_11_02 + #0.941 #0.034 IVOSOMG_02_03 + #0.941 #0.034 IVOSOMG_03_03 + &
          #0.941 #0.034 IVOSOMG_04_03 + #0.941 #0.034 IVOSOMG_05_03 + #0.941 &
          #0.034 IVOSOMG_06_03 + #0.941 #0.034 IVOSOMG_07_03 + #0.941 #0.034 &
          IVOSOMG_08_03 + #0.941 #0.034 IVOSOMG_09_03 + #0.941 #0.000 IVOSOMG_10_03 + &
          #0.941 #0.000 IVOSOMG_11_03 + #0.941 #0.034 IVOSOMG_02_04 + #0.941 &
          #0.034 IVOSOMG_03_04 + #0.941 #0.034 IVOSOMG_04_04 + #0.941 #0.034 &
          IVOSOMG_05_04 + #0.941 #0.034 IVOSOMG_06_04 + #0.941 #0.034 IVOSOMG_07_04 + &
          #0.941 #0.034 IVOSOMG_08_04 + #0.941 #0.034 IVOSOMG_09_04 + #0.941 &
          #0.034 IVOSOMG_10_04 + #0.941 #0.000 IVOSOMG_11_04 + #0.941 #0.034 &
          IVOSOMG_03_05 + #0.941 #0.034 IVOSOMG_04_05 + #0.941 #0.034 IVOSOMG_05_05 + &
          #0.941 #0.034 IVOSOMG_06_05 + #0.941 #0.034 IVOSOMG_07_05 + #0.941 &
          #0.034 IVOSOMG_08_05 + #0.941 #0.034 IVOSOMG_09_05 + #0.941 #0.034 &
          IVOSOMG_10_05 + #0.941 #0.000 IVOSOMG_11_05 + #0.941 #0.034 IVOSOMG_03_06 + &
          #0.941 #0.034 IVOSOMG_04_06 + #0.941 #0.034 IVOSOMG_05_06 + #0.941 &
          #0.034 IVOSOMG_06_06 + #0.941 #0.034 IVOSOMG_07_06 + #0.941 #0.034 &
          IVOSOMG_08_06 + #0.941 #0.034 IVOSOMG_09_06 + #0.941 #0.034 IVOSOMG_10_06 + &
          #0.941 #0.034 IVOSOMG_11_06 + #0.941 #0.034 IVOSOMG_04_07 + #0.941 &
          #0.034 IVOSOMG_05_07 + #0.941 #0.034 IVOSOMG_06_07 + #0.941 #0.034 &
          IVOSOMG_07_07 + #0.941 #0.034 IVOSOMG_08_07 + #0.941 #0.034 IVOSOMG_09_07 + &
          #0.941 #0.034 IVOSOMG_10_07 + #0.941 #0.034 IVOSOMG_11_07
 S325.6 IVOSOMG_12_07 + OH = #1.000 #0.036 IVOSOMG_01_01 + #1.000 #0.036 &
          IVOSOMG_02_01 + #1.000 #0.036 IVOSOMG_03_01 + #1.000 #0.036 IVOSOMG_04_01 + &
          #1.000 #0.036 IVOSOMG_05_01 + #1.000 #0.036 IVOSOMG_06_01 + #1.000 &
          #0.036 IVOSOMG_07_01 + #1.000 #0.036 IVOSOMG_08_01 + #1.000 #0.000 &
          IVOSOMG_09_01 + #1.000 #0.000 IVOSOMG_10_01 + #1.000 #0.000 IVOSOMG_11_01 + &
          #1.000 #0.036 IVOSOMG_01_02 + #1.000 #0.036 IVOSOMG_02_02 + #1.000 &
          #0.036 IVOSOMG_03_02 + #1.000 #0.036 IVOSOMG_04_02 + #1.000 #0.036 &
          IVOSOMG_05_02 + #1.000 #0.036 IVOSOMG_06_02 + #1.000 #0.036 IVOSOMG_07_02 + &
          #1.000 #0.036 IVOSOMG_08_02 + #1.000 #0.000 IVOSOMG_09_02 + #1.000 &
          #0.000 IVOSOMG_10_02 + #1.000 #0.000 IVOSOMG_11_02 + #1.000 #0.036 &
          IVOSOMG_02_03 + #1.000 #0.036 IVOSOMG_03_03 + #1.000 #0.036 IVOSOMG_04_03 + &
          #1.000 #0.036 IVOSOMG_05_03 + #1.000 #0.036 IVOSOMG_06_03 + #1.000 &
          #0.036 IVOSOMG_07_03 + #1.000 #0.036 IVOSOMG_08_03 + #1.000 #0.036 &
          IVOSOMG_09_03 + #1.000 #0.000 IVOSOMG_10_03 + #1.000 #0.000 IVOSOMG_11_03 + &
          #1.000 #0.036 IVOSOMG_02_04 + #1.000 #0.036 IVOSOMG_03_04 + #1.000 &
          #0.036 IVOSOMG_04_04 + #1.000 #0.036 IVOSOMG_05_04 + #1.000 #0.036 &
          IVOSOMG_06_04 + #1.000 #0.036 IVOSOMG_07_04 + #1.000 #0.036 IVOSOMG_08_04 + &
          #1.000 #0.036 IVOSOMG_09_04 + #1.000 #0.000 IVOSOMG_10_04 + #1.000 &
          #0.000 IVOSOMG_11_04 + #1.000 #0.036 IVOSOMG_03_05 + #1.000 #0.036 &
          IVOSOMG_04_05 + #1.000 #0.036 IVOSOMG_05_05 + #1.000 #0.036 IVOSOMG_06_05 + &
          #1.000 #0.036 IVOSOMG_07_05 + #1.000 #0.036 IVOSOMG_08_05 + #1.000 &
          #0.036 IVOSOMG_09_05 + #1.000 #0.036 IVOSOMG_10_05 + #1.000 #0.000 &
          IVOSOMG_11_05 + #1.000 #0.036 IVOSOMG_03_06 + #1.000 #0.036 IVOSOMG_04_06 + &
          #1.000 #0.036 IVOSOMG_05_06 + #1.000 #0.036 IVOSOMG_06_06 + #1.000 &
          #0.036 IVOSOMG_07_06 + #1.000 #0.036 IVOSOMG_08_06 + #1.000 #0.036 &
          IVOSOMG_09_06 + #1.000 #0.036 IVOSOMG_10_06 + #1.000 #0.000 IVOSOMG_11_06 + &
          #1.000 #0.036 IVOSOMG_04_07 + #1.000 #0.036 IVOSOMG_05_07 + #1.000 &
          #0.036 IVOSOMG_06_07 + #1.000 #0.036 IVOSOMG_07_07 + #1.000 #0.036 &
          IVOSOMG_08_07 + #1.000 #0.036 IVOSOMG_09_07 + #1.000 #0.036 IVOSOMG_10_07 + &
          #1.000 #0.036 IVOSOMG_11_07
 S326.6 IVOSOMG_13_01 + OH = #0.200 #0.361 IVOSOMG_13_02 + #0.200 #0.007 &
          IVOSOMG_13_03 + #0.200 #0.606 IVOSOMG_13_04 + #0.200 #0.026 IVOSOMG_13_05 + &
          #0.800 #0.083 IVOSOMG_01_01 + #0.800 #0.083 IVOSOMG_02_01 + #0.800 &
          #0.083 IVOSOMG_03_01 + #0.800 #0.083 IVOSOMG_04_01 + #0.800 #0.083 &
          IVOSOMG_05_01 + #0.800 #0.083 IVOSOMG_06_01 + #0.800 #0.083 IVOSOMG_07_01 + &
          #0.800 #0.083 IVOSOMG_08_01 + #0.800 #0.083 IVOSOMG_09_01 + #0.800 &
          #0.083 IVOSOMG_10_01 + #0.800 #0.083 IVOSOMG_11_01 + #0.800 #0.083 &
          IVOSOMG_12_01 + #0.800 #0.083 IVOSOMG_01_02 + #0.800 #0.083 IVOSOMG_02_02 + &
          #0.800 #0.083 IVOSOMG_03_02 + #0.800 #0.083 IVOSOMG_04_02 + #0.800 &
          #0.083 IVOSOMG_05_02 + #0.800 #0.083 IVOSOMG_06_02 + #0.800 #0.083 &
          IVOSOMG_07_02 + #0.800 #0.083 IVOSOMG_08_02 + #0.800 #0.083 IVOSOMG_09_02 + &
          #0.800 #0.083 IVOSOMG_10_02 + #0.800 #0.083 IVOSOMG_11_02 + #0.800 &
          #0.083 IVOSOMG_12_02
 S327.6 IVOSOMG_13_02 + OH = #0.150 #0.361 IVOSOMG_13_03 + #0.150 #0.007 &
          IVOSOMG_13_04 + #0.150 #0.606 IVOSOMG_13_05 + #0.150 #0.026 IVOSOMG_13_06 + &
          #0.850 #0.057 IVOSOMG_01_01 + #0.850 #0.057 IVOSOMG_02_01 + #0.850 &
          #0.057 IVOSOMG_03_01 + #0.850 #0.057 IVOSOMG_04_01 + #0.850 #0.057 &
          IVOSOMG_05_01 + #0.850 #0.057 IVOSOMG_06_01 + #0.850 #0.057 IVOSOMG_07_01 + &
          #0.850 #0.057 IVOSOMG_08_01 + #0.850 #0.057 IVOSOMG_09_01 + #0.850 &
          #0.057 IVOSOMG_10_01 + #0.850 #0.057 IVOSOMG_11_01 + #0.850 #0.000 &
          IVOSOMG_12_01 + #0.850 #0.057 IVOSOMG_01_02 + #0.850 #0.057 IVOSOMG_02_02 + &
          #0.850 #0.057 IVOSOMG_03_02 + #0.850 #0.057 IVOSOMG_04_02 + #0.850 &
          #0.057 IVOSOMG_05_02 + #0.850 #0.057 IVOSOMG_06_02 + #0.850 #0.057 &
          IVOSOMG_07_02 + #0.850 #0.057 IVOSOMG_08_02 + #0.850 #0.057 IVOSOMG_09_02 + &
          #0.850 #0.057 IVOSOMG_10_02 + #0.850 #0.057 IVOSOMG_11_02 + #0.850 &
          #0.057 IVOSOMG_12_02 + #0.850 #0.057 IVOSOMG_02_03 + #0.850 #0.057 &
          IVOSOMG_03_03 + #0.850 #0.057 IVOSOMG_04_03 + #0.850 #0.057 IVOSOMG_05_03 + &
          #0.850 #0.057 IVOSOMG_06_03 + #0.850 #0.057 IVOSOMG_07_03 + #0.850 &
          #0.057 IVOSOMG_08_03 + #0.850 #0.057 IVOSOMG_09_03 + #0.850 #0.057 &
          IVOSOMG_10_03 + #0.850 #0.057 IVOSOMG_11_03 + #0.850 #0.057 IVOSOMG_12_03
 S328.6 IVOSOMG_13_03 + OH = #0.120 #0.361 IVOSOMG_13_04 + #0.120 #0.007 &
          IVOSOMG_13_05 + #0.120 #0.606 IVOSOMG_13_06 + #0.120 #0.026 IVOSOMG_13_07 + &
          #0.880 #0.044 IVOSOMG_01_01 + #0.880 #0.044 IVOSOMG_02_01 + #0.880 &
          #0.044 IVOSOMG_03_01 + #0.880 #0.044 IVOSOMG_04_01 + #0.880 #0.044 &
          IVOSOMG_05_01 + #0.880 #0.044 IVOSOMG_06_01 + #0.880 #0.044 IVOSOMG_07_01 + &
          #0.880 #0.044 IVOSOMG_08_01 + #0.880 #0.044 IVOSOMG_09_01 + #0.880 &
          #0.044 IVOSOMG_10_01 + #0.880 #0.044 IVOSOMG_11_01 + #0.880 #0.000 &
          IVOSOMG_12_01 + #0.880 #0.044 IVOSOMG_01_02 + #0.880 #0.044 IVOSOMG_02_02 + &
          #0.880 #0.044 IVOSOMG_03_02 + #0.880 #0.044 IVOSOMG_04_02 + #0.880 &
          #0.044 IVOSOMG_05_02 + #0.880 #0.044 IVOSOMG_06_02 + #0.880 #0.044 &
          IVOSOMG_07_02 + #0.880 #0.044 IVOSOMG_08_02 + #0.880 #0.044 IVOSOMG_09_02 + &
          #0.880 #0.044 IVOSOMG_10_02 + #0.880 #0.044 IVOSOMG_11_02 + #0.880 &
          #0.000 IVOSOMG_12_02 + #0.880 #0.044 IVOSOMG_02_03 + #0.880 #0.044 &
          IVOSOMG_03_03 + #0.880 #0.044 IVOSOMG_04_03 + #0.880 #0.044 IVOSOMG_05_03 + &
          #0.880 #0.044 IVOSOMG_06_03 + #0.880 #0.044 IVOSOMG_07_03 + #0.880 &
          #0.044 IVOSOMG_08_03 + #0.880 #0.044 IVOSOMG_09_03 + #0.880 #0.044 &
          IVOSOMG_10_03 + #0.880 #0.044 IVOSOMG_11_03 + #0.880 #0.044 IVOSOMG_12_03 + &
          #0.880 #0.044 IVOSOMG_02_04 + #0.880 #0.044 IVOSOMG_03_04 + #0.880 &
          #0.044 IVOSOMG_04_04 + #0.880 #0.044 IVOSOMG_05_04 + #0.880 #0.044 &
          IVOSOMG_06_04 + #0.880 #0.044 IVOSOMG_07_04 + #0.880 #0.044 IVOSOMG_08_04 + &
          #0.880 #0.044 IVOSOMG_09_04 + #0.880 #0.044 IVOSOMG_10_04 + #0.880 &
          #0.044 IVOSOMG_11_04 + #0.880 #0.044 IVOSOMG_12_04
 S329.6 IVOSOMG_13_04 + OH = #0.097 #0.361 IVOSOMG_13_05 + #0.097 #0.007 &
          IVOSOMG_13_06 + #0.097 #0.632 IVOSOMG_13_07 + #0.903 #0.037 IVOSOMG_01_01 + &
          #0.903 #0.037 IVOSOMG_02_01 + #0.903 #0.037 IVOSOMG_03_01 + #0.903 &
          #0.037 IVOSOMG_04_01 + #0.903 #0.037 IVOSOMG_05_01 + #0.903 #0.037 &
          IVOSOMG_06_01 + #0.903 #0.037 IVOSOMG_07_01 + #0.903 #0.037 IVOSOMG_08_01 + &
          #0.903 #0.037 IVOSOMG_09_01 + #0.903 #0.037 IVOSOMG_10_01 + #0.903 &
          #0.000 IVOSOMG_11_01 + #0.903 #0.000 IVOSOMG_12_01 + #0.903 #0.037 &
          IVOSOMG_01_02 + #0.903 #0.037 IVOSOMG_02_02 + #0.903 #0.037 IVOSOMG_03_02 + &
          #0.903 #0.037 IVOSOMG_04_02 + #0.903 #0.037 IVOSOMG_05_02 + #0.903 &
          #0.037 IVOSOMG_06_02 + #0.903 #0.037 IVOSOMG_07_02 + #0.903 #0.037 &
          IVOSOMG_08_02 + #0.903 #0.037 IVOSOMG_09_02 + #0.903 #0.037 IVOSOMG_10_02 + &
          #0.903 #0.037 IVOSOMG_11_02 + #0.903 #0.000 IVOSOMG_12_02 + #0.903 &
          #0.037 IVOSOMG_02_03 + #0.903 #0.037 IVOSOMG_03_03 + #0.903 #0.037 &
          IVOSOMG_04_03 + #0.903 #0.037 IVOSOMG_05_03 + #0.903 #0.037 IVOSOMG_06_03 + &
          #0.903 #0.037 IVOSOMG_07_03 + #0.903 #0.037 IVOSOMG_08_03 + #0.903 &
          #0.037 IVOSOMG_09_03 + #0.903 #0.037 IVOSOMG_10_03 + #0.903 #0.037 &
          IVOSOMG_11_03 + #0.903 #0.000 IVOSOMG_12_03 + #0.903 #0.037 IVOSOMG_02_04 + &
          #0.903 #0.037 IVOSOMG_03_04 + #0.903 #0.037 IVOSOMG_04_04 + #0.903 &
          #0.037 IVOSOMG_05_04 + #0.903 #0.037 IVOSOMG_06_04 + #0.903 #0.037 &
          IVOSOMG_07_04 + #0.903 #0.037 IVOSOMG_08_04 + #0.903 #0.037 IVOSOMG_09_04 + &
          #0.903 #0.037 IVOSOMG_10_04 + #0.903 #0.037 IVOSOMG_11_04 + #0.903 &
          #0.037 IVOSOMG_12_04 + #0.903 #0.037 IVOSOMG_03_05 + #0.903 #0.037 &
          IVOSOMG_04_05 + #0.903 #0.037 IVOSOMG_05_05 + #0.903 #0.037 IVOSOMG_06_05 + &
          #0.903 #0.037 IVOSOMG_07_05 + #0.903 #0.037 IVOSOMG_08_05 + #0.903 &
          #0.037 IVOSOMG_09_05 + #0.903 #0.037 IVOSOMG_10_05 + #0.903 #0.037 &
          IVOSOMG_11_05 + #0.903 #0.037 IVOSOMG_12_05
 S330.6 IVOSOMG_13_05 + OH = #0.080 #0.361 IVOSOMG_13_06 + #0.080 #0.639 &
          IVOSOMG_13_07 + #0.920 #0.032 IVOSOMG_01_01 + #0.920 #0.032 IVOSOMG_02_01 + &
          #0.920 #0.032 IVOSOMG_03_01 + #0.920 #0.032 IVOSOMG_04_01 + #0.920 &
          #0.032 IVOSOMG_05_01 + #0.920 #0.032 IVOSOMG_06_01 + #0.920 #0.032 &
          IVOSOMG_07_01 + #0.920 #0.032 IVOSOMG_08_01 + #0.920 #0.032 IVOSOMG_09_01 + &
          #0.920 #0.032 IVOSOMG_10_01 + #0.920 #0.000 IVOSOMG_11_01 + #0.920 &
          #0.000 IVOSOMG_12_01 + #0.920 #0.032 IVOSOMG_01_02 + #0.920 #0.032 &
          IVOSOMG_02_02 + #0.920 #0.032 IVOSOMG_03_02 + #0.920 #0.032 IVOSOMG_04_02 + &
          #0.920 #0.032 IVOSOMG_05_02 + #0.920 #0.032 IVOSOMG_06_02 + #0.920 &
          #0.032 IVOSOMG_07_02 + #0.920 #0.032 IVOSOMG_08_02 + #0.920 #0.032 &
          IVOSOMG_09_02 + #0.920 #0.032 IVOSOMG_10_02 + #0.920 #0.000 IVOSOMG_11_02 + &
          #0.920 #0.000 IVOSOMG_12_02 + #0.920 #0.032 IVOSOMG_02_03 + #0.920 &
          #0.032 IVOSOMG_03_03 + #0.920 #0.032 IVOSOMG_04_03 + #0.920 #0.032 &
          IVOSOMG_05_03 + #0.920 #0.032 IVOSOMG_06_03 + #0.920 #0.032 IVOSOMG_07_03 + &
          #0.920 #0.032 IVOSOMG_08_03 + #0.920 #0.032 IVOSOMG_09_03 + #0.920 &
          #0.032 IVOSOMG_10_03 + #0.920 #0.032 IVOSOMG_11_03 + #0.920 #0.000 &
          IVOSOMG_12_03 + #0.920 #0.032 IVOSOMG_02_04 + #0.920 #0.032 IVOSOMG_03_04 + &
          #0.920 #0.032 IVOSOMG_04_04 + #0.920 #0.032 IVOSOMG_05_04 + #0.920 &
          #0.032 IVOSOMG_06_04 + #0.920 #0.032 IVOSOMG_07_04 + #0.920 #0.032 &
          IVOSOMG_08_04 + #0.920 #0.032 IVOSOMG_09_04 + #0.920 #0.032 IVOSOMG_10_04 + &
          #0.920 #0.032 IVOSOMG_11_04 + #0.920 #0.000 IVOSOMG_12_04 + #0.920 &
          #0.032 IVOSOMG_03_05 + #0.920 #0.032 IVOSOMG_04_05 + #0.920 #0.032 &
          IVOSOMG_05_05 + #0.920 #0.032 IVOSOMG_06_05 + #0.920 #0.032 IVOSOMG_07_05 + &
          #0.920 #0.032 IVOSOMG_08_05 + #0.920 #0.032 IVOSOMG_09_05 + #0.920 &
          #0.032 IVOSOMG_10_05 + #0.920 #0.032 IVOSOMG_11_05 + #0.920 #0.032 &
          IVOSOMG_12_05 + #0.920 #0.032 IVOSOMG_03_06 + #0.920 #0.032 IVOSOMG_04_06 + &
          #0.920 #0.032 IVOSOMG_05_06 + #0.920 #0.032 IVOSOMG_06_06 + #0.920 &
          #0.032 IVOSOMG_07_06 + #0.920 #0.032 IVOSOMG_08_06 + #0.920 #0.032 &
          IVOSOMG_09_06 + #0.920 #0.032 IVOSOMG_10_06 + #0.920 #0.032 IVOSOMG_11_06 + &
          #0.920 #0.032 IVOSOMG_12_06
 S331.6 IVOSOMG_13_06 + OH = #0.065 #1.000 IVOSOMG_13_07 + #0.935 #0.029 &
          IVOSOMG_01_01 + #0.935 #0.029 IVOSOMG_02_01 + #0.935 #0.029 IVOSOMG_03_01 + &
          #0.935 #0.029 IVOSOMG_04_01 + #0.935 #0.029 IVOSOMG_05_01 + #0.935 &
          #0.029 IVOSOMG_06_01 + #0.935 #0.029 IVOSOMG_07_01 + #0.935 #0.029 &
          IVOSOMG_08_01 + #0.935 #0.029 IVOSOMG_09_01 + #0.935 #0.000 IVOSOMG_10_01 + &
          #0.935 #0.000 IVOSOMG_11_01 + #0.935 #0.000 IVOSOMG_12_01 + #0.935 &
          #0.029 IVOSOMG_01_02 + #0.935 #0.029 IVOSOMG_02_02 + #0.935 #0.029 &
          IVOSOMG_03_02 + #0.935 #0.029 IVOSOMG_04_02 + #0.935 #0.029 IVOSOMG_05_02 + &
          #0.935 #0.029 IVOSOMG_06_02 + #0.935 #0.029 IVOSOMG_07_02 + #0.935 &
          #0.029 IVOSOMG_08_02 + #0.935 #0.029 IVOSOMG_09_02 + #0.935 #0.029 &
          IVOSOMG_10_02 + #0.935 #0.000 IVOSOMG_11_02 + #0.935 #0.000 IVOSOMG_12_02 + &
          #0.935 #0.029 IVOSOMG_02_03 + #0.935 #0.029 IVOSOMG_03_03 + #0.935 &
          #0.029 IVOSOMG_04_03 + #0.935 #0.029 IVOSOMG_05_03 + #0.935 #0.029 &
          IVOSOMG_06_03 + #0.935 #0.029 IVOSOMG_07_03 + #0.935 #0.029 IVOSOMG_08_03 + &
          #0.935 #0.029 IVOSOMG_09_03 + #0.935 #0.029 IVOSOMG_10_03 + #0.935 &
          #0.000 IVOSOMG_11_03 + #0.935 #0.000 IVOSOMG_12_03 + #0.935 #0.029 &
          IVOSOMG_02_04 + #0.935 #0.029 IVOSOMG_03_04 + #0.935 #0.029 IVOSOMG_04_04 + &
          #0.935 #0.029 IVOSOMG_05_04 + #0.935 #0.029 IVOSOMG_06_04 + #0.935 &
          #0.029 IVOSOMG_07_04 + #0.935 #0.029 IVOSOMG_08_04 + #0.935 #0.029 &
          IVOSOMG_09_04 + #0.935 #0.029 IVOSOMG_10_04 + #0.935 #0.029 IVOSOMG_11_04 + &
          #0.935 #0.000 IVOSOMG_12_04 + #0.935 #0.029 IVOSOMG_03_05 + #0.935 &
          #0.029 IVOSOMG_04_05 + #0.935 #0.029 IVOSOMG_05_05 + #0.935 #0.029 &
          IVOSOMG_06_05 + #0.935 #0.029 IVOSOMG_07_05 + #0.935 #0.029 IVOSOMG_08_05 + &
          #0.935 #0.029 IVOSOMG_09_05 + #0.935 #0.029 IVOSOMG_10_05 + #0.935 &
          #0.029 IVOSOMG_11_05 + #0.935 #0.000 IVOSOMG_12_05 + #0.935 #0.029 &
          IVOSOMG_03_06 + #0.935 #0.029 IVOSOMG_04_06 + #0.935 #0.029 IVOSOMG_05_06 + &
          #0.935 #0.029 IVOSOMG_06_06 + #0.935 #0.029 IVOSOMG_07_06 + #0.935 &
          #0.029 IVOSOMG_08_06 + #0.935 #0.029 IVOSOMG_09_06 + #0.935 #0.029 &
          IVOSOMG_10_06 + #0.935 #0.029 IVOSOMG_11_06 + #0.935 #0.029 IVOSOMG_12_06 + &
          #0.935 #0.029 IVOSOMG_04_07 + #0.935 #0.029 IVOSOMG_05_07 + #0.935 &
          #0.029 IVOSOMG_06_07 + #0.935 #0.029 IVOSOMG_07_07 + #0.935 #0.029 &
          IVOSOMG_08_07 + #0.935 #0.029 IVOSOMG_09_07 + #0.935 #0.029 IVOSOMG_10_07 + &
          #0.935 #0.029 IVOSOMG_11_07 + #0.935 #0.029 IVOSOMG_12_07
 S332.6 IVOSOMG_13_07 + OH = #1.000 #0.030 IVOSOMG_01_01 + #1.000 #0.030 &
          IVOSOMG_02_01 + #1.000 #0.030 IVOSOMG_03_01 + #1.000 #0.030 IVOSOMG_04_01 + &
          #1.000 #0.030 IVOSOMG_05_01 + #1.000 #0.030 IVOSOMG_06_01 + #1.000 &
          #0.030 IVOSOMG_07_01 + #1.000 #0.030 IVOSOMG_08_01 + #1.000 #0.030 &
          IVOSOMG_09_01 + #1.000 #0.000 IVOSOMG_10_01 + #1.000 #0.000 IVOSOMG_11_01 + &
          #1.000 #0.000 IVOSOMG_12_01 + #1.000 #0.030 IVOSOMG_01_02 + #1.000 &
          #0.030 IVOSOMG_02_02 + #1.000 #0.030 IVOSOMG_03_02 + #1.000 #0.030 &
          IVOSOMG_04_02 + #1.000 #0.030 IVOSOMG_05_02 + #1.000 #0.030 IVOSOMG_06_02 + &
          #1.000 #0.030 IVOSOMG_07_02 + #1.000 #0.030 IVOSOMG_08_02 + #1.000 &
          #0.030 IVOSOMG_09_02 + #1.000 #0.000 IVOSOMG_10_02 + #1.000 #0.000 &
          IVOSOMG_11_02 + #1.000 #0.000 IVOSOMG_12_02 + #1.000 #0.030 IVOSOMG_02_03 + &
          #1.000 #0.030 IVOSOMG_03_03 + #1.000 #0.030 IVOSOMG_04_03 + #1.000 &
          #0.030 IVOSOMG_05_03 + #1.000 #0.030 IVOSOMG_06_03 + #1.000 #0.030 &
          IVOSOMG_07_03 + #1.000 #0.030 IVOSOMG_08_03 + #1.000 #0.030 IVOSOMG_09_03 + &
          #1.000 #0.030 IVOSOMG_10_03 + #1.000 #0.000 IVOSOMG_11_03 + #1.000 &
          #0.000 IVOSOMG_12_03 + #1.000 #0.030 IVOSOMG_02_04 + #1.000 #0.030 &
          IVOSOMG_03_04 + #1.000 #0.030 IVOSOMG_04_04 + #1.000 #0.030 IVOSOMG_05_04 + &
          #1.000 #0.030 IVOSOMG_06_04 + #1.000 #0.030 IVOSOMG_07_04 + #1.000 &
          #0.030 IVOSOMG_08_04 + #1.000 #0.030 IVOSOMG_09_04 + #1.000 #0.030 &
          IVOSOMG_10_04 + #1.000 #0.000 IVOSOMG_11_04 + #1.000 #0.000 IVOSOMG_12_04 + &
          #1.000 #0.030 IVOSOMG_03_05 + #1.000 #0.030 IVOSOMG_04_05 + #1.000 &
          #0.030 IVOSOMG_05_05 + #1.000 #0.030 IVOSOMG_06_05 + #1.000 #0.030 &
          IVOSOMG_07_05 + #1.000 #0.030 IVOSOMG_08_05 + #1.000 #0.030 IVOSOMG_09_05 + &
          #1.000 #0.030 IVOSOMG_10_05 + #1.000 #0.030 IVOSOMG_11_05 + #1.000 &
          #0.000 IVOSOMG_12_05 + #1.000 #0.030 IVOSOMG_03_06 + #1.000 #0.030 &
          IVOSOMG_04_06 + #1.000 #0.030 IVOSOMG_05_06 + #1.000 #0.030 IVOSOMG_06_06 + &
          #1.000 #0.030 IVOSOMG_07_06 + #1.000 #0.030 IVOSOMG_08_06 + #1.000 &
          #0.030 IVOSOMG_09_06 + #1.000 #0.030 IVOSOMG_10_06 + #1.000 #0.030 &
          IVOSOMG_11_06 + #1.000 #0.000 IVOSOMG_12_06 + #1.000 #0.030 IVOSOMG_04_07 + &
          #1.000 #0.030 IVOSOMG_05_07 + #1.000 #0.030 IVOSOMG_06_07 + #1.000 &
          #0.030 IVOSOMG_07_07 + #1.000 #0.030 IVOSOMG_08_07 + #1.000 #0.030 &
          IVOSOMG_09_07 + #1.000 #0.030 IVOSOMG_10_07 + #1.000 #0.030 IVOSOMG_11_07 + &
          #1.000 #0.030 IVOSOMG_12_07
 S333.6 IVOSOMG_14_01 + OH = #0.205 #0.361 IVOSOMG_14_02 + #0.205 #0.007 &
          IVOSOMG_14_03 + #0.205 #0.606 IVOSOMG_14_04 + #0.205 #0.026 IVOSOMG_14_05 + &
          #0.795 #0.077 IVOSOMG_01_01 + #0.795 #0.077 IVOSOMG_02_01 + #0.795 &
          #0.077 IVOSOMG_03_01 + #0.795 #0.077 IVOSOMG_04_01 + #0.795 #0.077 &
          IVOSOMG_05_01 + #0.795 #0.077 IVOSOMG_06_01 + #0.795 #0.077 IVOSOMG_07_01 + &
          #0.795 #0.077 IVOSOMG_08_01 + #0.795 #0.077 IVOSOMG_09_01 + #0.795 &
          #0.077 IVOSOMG_10_01 + #0.795 #0.077 IVOSOMG_11_01 + #0.795 #0.077 &
          IVOSOMG_12_01 + #0.795 #0.077 IVOSOMG_13_01 + #0.795 #0.077 IVOSOMG_01_02 + &
          #0.795 #0.077 IVOSOMG_02_02 + #0.795 #0.077 IVOSOMG_03_02 + #0.795 &
          #0.077 IVOSOMG_04_02 + #0.795 #0.077 IVOSOMG_05_02 + #0.795 #0.077 &
          IVOSOMG_06_02 + #0.795 #0.077 IVOSOMG_07_02 + #0.795 #0.077 IVOSOMG_08_02 + &
          #0.795 #0.077 IVOSOMG_09_02 + #0.795 #0.077 IVOSOMG_10_02 + #0.795 &
          #0.077 IVOSOMG_11_02 + #0.795 #0.077 IVOSOMG_12_02 + #0.795 #0.077 &
          IVOSOMG_13_02
 S334.6 IVOSOMG_14_02 + OH = #0.156 #0.361 IVOSOMG_14_03 + #0.156 #0.007 &
          IVOSOMG_14_04 + #0.156 #0.606 IVOSOMG_14_05 + #0.156 #0.026 IVOSOMG_14_06 + &
          #0.844 #0.054 IVOSOMG_01_01 + #0.844 #0.054 IVOSOMG_02_01 + #0.844 &
          #0.054 IVOSOMG_03_01 + #0.844 #0.054 IVOSOMG_04_01 + #0.844 #0.054 &
          IVOSOMG_05_01 + #0.844 #0.054 IVOSOMG_06_01 + #0.844 #0.054 IVOSOMG_07_01 + &
          #0.844 #0.054 IVOSOMG_08_01 + #0.844 #0.054 IVOSOMG_09_01 + #0.844 &
          #0.054 IVOSOMG_10_01 + #0.844 #0.054 IVOSOMG_11_01 + #0.844 #0.054 &
          IVOSOMG_12_01 + #0.844 #0.000 IVOSOMG_13_01 + #0.844 #0.054 IVOSOMG_01_02 + &
          #0.844 #0.054 IVOSOMG_02_02 + #0.844 #0.054 IVOSOMG_03_02 + #0.844 &
          #0.054 IVOSOMG_04_02 + #0.844 #0.054 IVOSOMG_05_02 + #0.844 #0.054 &
          IVOSOMG_06_02 + #0.844 #0.054 IVOSOMG_07_02 + #0.844 #0.054 IVOSOMG_08_02 + &
          #0.844 #0.054 IVOSOMG_09_02 + #0.844 #0.054 IVOSOMG_10_02 + #0.844 &
          #0.054 IVOSOMG_11_02 + #0.844 #0.054 IVOSOMG_12_02 + #0.844 #0.054 &
          IVOSOMG_13_02 + #0.844 #0.054 IVOSOMG_02_03 + #0.844 #0.054 IVOSOMG_03_03 + &
          #0.844 #0.054 IVOSOMG_04_03 + #0.844 #0.054 IVOSOMG_05_03 + #0.844 &
          #0.054 IVOSOMG_06_03 + #0.844 #0.054 IVOSOMG_07_03 + #0.844 #0.054 &
          IVOSOMG_08_03 + #0.844 #0.054 IVOSOMG_09_03 + #0.844 #0.054 IVOSOMG_10_03 + &
          #0.844 #0.054 IVOSOMG_11_03 + #0.844 #0.054 IVOSOMG_12_03 + #0.844 &
          #0.054 IVOSOMG_13_03
 S335.6 IVOSOMG_14_03 + OH = #0.125 #0.361 IVOSOMG_14_04 + #0.125 #0.007 &
          IVOSOMG_14_05 + #0.125 #0.606 IVOSOMG_14_06 + #0.125 #0.026 IVOSOMG_14_07 + &
          #0.875 #0.042 IVOSOMG_01_01 + #0.875 #0.042 IVOSOMG_02_01 + #0.875 &
          #0.042 IVOSOMG_03_01 + #0.875 #0.042 IVOSOMG_04_01 + #0.875 #0.042 &
          IVOSOMG_05_01 + #0.875 #0.042 IVOSOMG_06_01 + #0.875 #0.042 IVOSOMG_07_01 + &
          #0.875 #0.042 IVOSOMG_08_01 + #0.875 #0.042 IVOSOMG_09_01 + #0.875 &
          #0.042 IVOSOMG_10_01 + #0.875 #0.042 IVOSOMG_11_01 + #0.875 #0.042 &
          IVOSOMG_12_01 + #0.875 #0.000 IVOSOMG_13_01 + #0.875 #0.042 IVOSOMG_01_02 + &
          #0.875 #0.042 IVOSOMG_02_02 + #0.875 #0.042 IVOSOMG_03_02 + #0.875 &
          #0.042 IVOSOMG_04_02 + #0.875 #0.042 IVOSOMG_05_02 + #0.875 #0.042 &
          IVOSOMG_06_02 + #0.875 #0.042 IVOSOMG_07_02 + #0.875 #0.042 IVOSOMG_08_02 + &
          #0.875 #0.042 IVOSOMG_09_02 + #0.875 #0.042 IVOSOMG_10_02 + #0.875 &
          #0.042 IVOSOMG_11_02 + #0.875 #0.042 IVOSOMG_12_02 + #0.875 #0.000 &
          IVOSOMG_13_02 + #0.875 #0.042 IVOSOMG_02_03 + #0.875 #0.042 IVOSOMG_03_03 + &
          #0.875 #0.042 IVOSOMG_04_03 + #0.875 #0.042 IVOSOMG_05_03 + #0.875 &
          #0.042 IVOSOMG_06_03 + #0.875 #0.042 IVOSOMG_07_03 + #0.875 #0.042 &
          IVOSOMG_08_03 + #0.875 #0.042 IVOSOMG_09_03 + #0.875 #0.042 IVOSOMG_10_03 + &
          #0.875 #0.042 IVOSOMG_11_03 + #0.875 #0.042 IVOSOMG_12_03 + #0.875 &
          #0.042 IVOSOMG_13_03 + #0.875 #0.042 IVOSOMG_02_04 + #0.875 #0.042 &
          IVOSOMG_03_04 + #0.875 #0.042 IVOSOMG_04_04 + #0.875 #0.042 IVOSOMG_05_04 + &
          #0.875 #0.042 IVOSOMG_06_04 + #0.875 #0.042 IVOSOMG_07_04 + #0.875 &
          #0.042 IVOSOMG_08_04 + #0.875 #0.042 IVOSOMG_09_04 + #0.875 #0.042 &
          IVOSOMG_10_04 + #0.875 #0.042 IVOSOMG_11_04 + #0.875 #0.042 IVOSOMG_12_04 + &
          #0.875 #0.042 IVOSOMG_13_04
 S336.6 IVOSOMG_14_04 + OH = #0.103 #0.361 IVOSOMG_14_05 + #0.103 #0.007 &
          IVOSOMG_14_06 + #0.103 #0.632 IVOSOMG_14_07 + #0.897 #0.034 IVOSOMG_01_01 + &
          #0.897 #0.034 IVOSOMG_02_01 + #0.897 #0.034 IVOSOMG_03_01 + #0.897 &
          #0.034 IVOSOMG_04_01 + #0.897 #0.034 IVOSOMG_05_01 + #0.897 #0.034 &
          IVOSOMG_06_01 + #0.897 #0.034 IVOSOMG_07_01 + #0.897 #0.034 IVOSOMG_08_01 + &
          #0.897 #0.034 IVOSOMG_09_01 + #0.897 #0.034 IVOSOMG_10_01 + #0.897 &
          #0.034 IVOSOMG_11_01 + #0.897 #0.000 IVOSOMG_12_01 + #0.897 #0.000 &
          IVOSOMG_13_01 + #0.897 #0.034 IVOSOMG_01_02 + #0.897 #0.034 IVOSOMG_02_02 + &
          #0.897 #0.034 IVOSOMG_03_02 + #0.897 #0.034 IVOSOMG_04_02 + #0.897 &
          #0.034 IVOSOMG_05_02 + #0.897 #0.034 IVOSOMG_06_02 + #0.897 #0.034 &
          IVOSOMG_07_02 + #0.897 #0.034 IVOSOMG_08_02 + #0.897 #0.034 IVOSOMG_09_02 + &
          #0.897 #0.034 IVOSOMG_10_02 + #0.897 #0.034 IVOSOMG_11_02 + #0.897 &
          #0.034 IVOSOMG_12_02 + #0.897 #0.000 IVOSOMG_13_02 + #0.897 #0.034 &
          IVOSOMG_02_03 + #0.897 #0.034 IVOSOMG_03_03 + #0.897 #0.034 IVOSOMG_04_03 + &
          #0.897 #0.034 IVOSOMG_05_03 + #0.897 #0.034 IVOSOMG_06_03 + #0.897 &
          #0.034 IVOSOMG_07_03 + #0.897 #0.034 IVOSOMG_08_03 + #0.897 #0.034 &
          IVOSOMG_09_03 + #0.897 #0.034 IVOSOMG_10_03 + #0.897 #0.034 IVOSOMG_11_03 + &
          #0.897 #0.034 IVOSOMG_12_03 + #0.897 #0.000 IVOSOMG_13_03 + #0.897 &
          #0.034 IVOSOMG_02_04 + #0.897 #0.034 IVOSOMG_03_04 + #0.897 #0.034 &
          IVOSOMG_04_04 + #0.897 #0.034 IVOSOMG_05_04 + #0.897 #0.034 IVOSOMG_06_04 + &
          #0.897 #0.034 IVOSOMG_07_04 + #0.897 #0.034 IVOSOMG_08_04 + #0.897 &
          #0.034 IVOSOMG_09_04 + #0.897 #0.034 IVOSOMG_10_04 + #0.897 #0.034 &
          IVOSOMG_11_04 + #0.897 #0.034 IVOSOMG_12_04 + #0.897 #0.034 IVOSOMG_13_04 + &
          #0.897 #0.034 IVOSOMG_03_05 + #0.897 #0.034 IVOSOMG_04_05 + #0.897 &
          #0.034 IVOSOMG_05_05 + #0.897 #0.034 IVOSOMG_06_05 + #0.897 #0.034 &
          IVOSOMG_07_05 + #0.897 #0.034 IVOSOMG_08_05 + #0.897 #0.034 IVOSOMG_09_05 + &
          #0.897 #0.034 IVOSOMG_10_05 + #0.897 #0.034 IVOSOMG_11_05 + #0.897 &
          #0.034 IVOSOMG_12_05 + #0.897 #0.034 IVOSOMG_13_05
 S337.6 IVOSOMG_14_05 + OH = #0.086 #0.361 IVOSOMG_14_06 + #0.086 #0.639 &
          IVOSOMG_14_07 + #0.914 #0.030 IVOSOMG_01_01 + #0.914 #0.030 IVOSOMG_02_01 + &
          #0.914 #0.030 IVOSOMG_03_01 + #0.914 #0.030 IVOSOMG_04_01 + #0.914 &
          #0.030 IVOSOMG_05_01 + #0.914 #0.030 IVOSOMG_06_01 + #0.914 #0.030 &
          IVOSOMG_07_01 + #0.914 #0.030 IVOSOMG_08_01 + #0.914 #0.030 IVOSOMG_09_01 + &
          #0.914 #0.030 IVOSOMG_10_01 + #0.914 #0.030 IVOSOMG_11_01 + #0.914 &
          #0.000 IVOSOMG_12_01 + #0.914 #0.000 IVOSOMG_13_01 + #0.914 #0.030 &
          IVOSOMG_01_02 + #0.914 #0.030 IVOSOMG_02_02 + #0.914 #0.030 IVOSOMG_03_02 + &
          #0.914 #0.030 IVOSOMG_04_02 + #0.914 #0.030 IVOSOMG_05_02 + #0.914 &
          #0.030 IVOSOMG_06_02 + #0.914 #0.030 IVOSOMG_07_02 + #0.914 #0.030 &
          IVOSOMG_08_02 + #0.914 #0.030 IVOSOMG_09_02 + #0.914 #0.030 IVOSOMG_10_02 + &
          #0.914 #0.030 IVOSOMG_11_02 + #0.914 #0.000 IVOSOMG_12_02 + #0.914 &
          #0.000 IVOSOMG_13_02 + #0.914 #0.030 IVOSOMG_02_03 + #0.914 #0.030 &
          IVOSOMG_03_03 + #0.914 #0.030 IVOSOMG_04_03 + #0.914 #0.030 IVOSOMG_05_03 + &
          #0.914 #0.030 IVOSOMG_06_03 + #0.914 #0.030 IVOSOMG_07_03 + #0.914 &
          #0.030 IVOSOMG_08_03 + #0.914 #0.030 IVOSOMG_09_03 + #0.914 #0.030 &
          IVOSOMG_10_03 + #0.914 #0.030 IVOSOMG_11_03 + #0.914 #0.030 IVOSOMG_12_03 + &
          #0.914 #0.000 IVOSOMG_13_03 + #0.914 #0.030 IVOSOMG_02_04 + #0.914 &
          #0.030 IVOSOMG_03_04 + #0.914 #0.030 IVOSOMG_04_04 + #0.914 #0.030 &
          IVOSOMG_05_04 + #0.914 #0.030 IVOSOMG_06_04 + #0.914 #0.030 IVOSOMG_07_04 + &
          #0.914 #0.030 IVOSOMG_08_04 + #0.914 #0.030 IVOSOMG_09_04 + #0.914 &
          #0.030 IVOSOMG_10_04 + #0.914 #0.030 IVOSOMG_11_04 + #0.914 #0.030 &
          IVOSOMG_12_04 + #0.914 #0.000 IVOSOMG_13_04 + #0.914 #0.030 IVOSOMG_03_05 + &
          #0.914 #0.030 IVOSOMG_04_05 + #0.914 #0.030 IVOSOMG_05_05 + #0.914 &
          #0.030 IVOSOMG_06_05 + #0.914 #0.030 IVOSOMG_07_05 + #0.914 #0.030 &
          IVOSOMG_08_05 + #0.914 #0.030 IVOSOMG_09_05 + #0.914 #0.030 IVOSOMG_10_05 + &
          #0.914 #0.030 IVOSOMG_11_05 + #0.914 #0.030 IVOSOMG_12_05 + #0.914 &
          #0.030 IVOSOMG_13_05 + #0.914 #0.030 IVOSOMG_03_06 + #0.914 #0.030 &
          IVOSOMG_04_06 + #0.914 #0.030 IVOSOMG_05_06 + #0.914 #0.030 IVOSOMG_06_06 + &
          #0.914 #0.030 IVOSOMG_07_06 + #0.914 #0.030 IVOSOMG_08_06 + #0.914 &
          #0.030 IVOSOMG_09_06 + #0.914 #0.030 IVOSOMG_10_06 + #0.914 #0.030 &
          IVOSOMG_11_06 + #0.914 #0.030 IVOSOMG_12_06 + #0.914 #0.030 IVOSOMG_13_06
 S338.6 IVOSOMG_14_06 + OH = #0.071 #1.000 IVOSOMG_14_07 + #0.929 #0.027 &
          IVOSOMG_01_01 + #0.929 #0.027 IVOSOMG_02_01 + #0.929 #0.027 IVOSOMG_03_01 + &
          #0.929 #0.027 IVOSOMG_04_01 + #0.929 #0.027 IVOSOMG_05_01 + #0.929 &
          #0.027 IVOSOMG_06_01 + #0.929 #0.027 IVOSOMG_07_01 + #0.929 #0.027 &
          IVOSOMG_08_01 + #0.929 #0.027 IVOSOMG_09_01 + #0.929 #0.027 IVOSOMG_10_01 + &
          #0.929 #0.000 IVOSOMG_11_01 + #0.929 #0.000 IVOSOMG_12_01 + #0.929 &
          #0.000 IVOSOMG_13_01 + #0.929 #0.027 IVOSOMG_01_02 + #0.929 #0.027 &
          IVOSOMG_02_02 + #0.929 #0.027 IVOSOMG_03_02 + #0.929 #0.027 IVOSOMG_04_02 + &
          #0.929 #0.027 IVOSOMG_05_02 + #0.929 #0.027 IVOSOMG_06_02 + #0.929 &
          #0.027 IVOSOMG_07_02 + #0.929 #0.027 IVOSOMG_08_02 + #0.929 #0.027 &
          IVOSOMG_09_02 + #0.929 #0.027 IVOSOMG_10_02 + #0.929 #0.027 IVOSOMG_11_02 + &
          #0.929 #0.000 IVOSOMG_12_02 + #0.929 #0.000 IVOSOMG_13_02 + #0.929 &
          #0.027 IVOSOMG_02_03 + #0.929 #0.027 IVOSOMG_03_03 + #0.929 #0.027 &
          IVOSOMG_04_03 + #0.929 #0.027 IVOSOMG_05_03 + #0.929 #0.027 IVOSOMG_06_03 + &
          #0.929 #0.027 IVOSOMG_07_03 + #0.929 #0.027 IVOSOMG_08_03 + #0.929 &
          #0.027 IVOSOMG_09_03 + #0.929 #0.027 IVOSOMG_10_03 + #0.929 #0.027 &
          IVOSOMG_11_03 + #0.929 #0.000 IVOSOMG_12_03 + #0.929 #0.000 IVOSOMG_13_03 + &
          #0.929 #0.027 IVOSOMG_02_04 + #0.929 #0.027 IVOSOMG_03_04 + #0.929 &
          #0.027 IVOSOMG_04_04 + #0.929 #0.027 IVOSOMG_05_04 + #0.929 #0.027 &
          IVOSOMG_06_04 + #0.929 #0.027 IVOSOMG_07_04 + #0.929 #0.027 IVOSOMG_08_04 + &
          #0.929 #0.027 IVOSOMG_09_04 + #0.929 #0.027 IVOSOMG_10_04 + #0.929 &
          #0.027 IVOSOMG_11_04 + #0.929 #0.027 IVOSOMG_12_04 + #0.929 #0.000 &
          IVOSOMG_13_04 + #0.929 #0.027 IVOSOMG_03_05 + #0.929 #0.027 IVOSOMG_04_05 + &
          #0.929 #0.027 IVOSOMG_05_05 + #0.929 #0.027 IVOSOMG_06_05 + #0.929 &
          #0.027 IVOSOMG_07_05 + #0.929 #0.027 IVOSOMG_08_05 + #0.929 #0.027 &
          IVOSOMG_09_05 + #0.929 #0.027 IVOSOMG_10_05 + #0.929 #0.027 IVOSOMG_11_05 + &
          #0.929 #0.027 IVOSOMG_12_05 + #0.929 #0.000 IVOSOMG_13_05 + #0.929 &
          #0.027 IVOSOMG_03_06 + #0.929 #0.027 IVOSOMG_04_06 + #0.929 #0.027 &
          IVOSOMG_05_06 + #0.929 #0.027 IVOSOMG_06_06 + #0.929 #0.027 IVOSOMG_07_06 + &
          #0.929 #0.027 IVOSOMG_08_06 + #0.929 #0.027 IVOSOMG_09_06 + #0.929 &
          #0.027 IVOSOMG_10_06 + #0.929 #0.027 IVOSOMG_11_06 + #0.929 #0.027 &
          IVOSOMG_12_06 + #0.929 #0.027 IVOSOMG_13_06 + #0.929 #0.027 IVOSOMG_04_07 + &
          #0.929 #0.027 IVOSOMG_05_07 + #0.929 #0.027 IVOSOMG_06_07 + #0.929 &
          #0.027 IVOSOMG_07_07 + #0.929 #0.027 IVOSOMG_08_07 + #0.929 #0.027 &
          IVOSOMG_09_07 + #0.929 #0.027 IVOSOMG_10_07 + #0.929 #0.027 IVOSOMG_11_07 + &
          #0.929 #0.027 IVOSOMG_12_07 + #0.929 #0.027 IVOSOMG_13_07
 S339.6 IVOSOMG_14_07 + OH = #1.000 #0.028 IVOSOMG_01_01 + #1.000 #0.028 &
          IVOSOMG_02_01 + #1.000 #0.028 IVOSOMG_03_01 + #1.000 #0.028 IVOSOMG_04_01 + &
          #1.000 #0.028 IVOSOMG_05_01 + #1.000 #0.028 IVOSOMG_06_01 + #1.000 &
          #0.028 IVOSOMG_07_01 + #1.000 #0.028 IVOSOMG_08_01 + #1.000 #0.028 &
          IVOSOMG_09_01 + #1.000 #0.028 IVOSOMG_10_01 + #1.000 #0.000 IVOSOMG_11_01 + &
          #1.000 #0.000 IVOSOMG_12_01 + #1.000 #0.000 IVOSOMG_13_01 + #1.000 &
          #0.028 IVOSOMG_01_02 + #1.000 #0.028 IVOSOMG_02_02 + #1.000 #0.028 &
          IVOSOMG_03_02 + #1.000 #0.028 IVOSOMG_04_02 + #1.000 #0.028 IVOSOMG_05_02 + &
          #1.000 #0.028 IVOSOMG_06_02 + #1.000 #0.028 IVOSOMG_07_02 + #1.000 &
          #0.028 IVOSOMG_08_02 + #1.000 #0.028 IVOSOMG_09_02 + #1.000 #0.028 &
          IVOSOMG_10_02 + #1.000 #0.000 IVOSOMG_11_02 + #1.000 #0.000 IVOSOMG_12_02 + &
          #1.000 #0.000 IVOSOMG_13_02 + #1.000 #0.028 IVOSOMG_02_03 + #1.000 &
          #0.028 IVOSOMG_03_03 + #1.000 #0.028 IVOSOMG_04_03 + #1.000 #0.028 &
          IVOSOMG_05_03 + #1.000 #0.028 IVOSOMG_06_03 + #1.000 #0.028 IVOSOMG_07_03 + &
          #1.000 #0.028 IVOSOMG_08_03 + #1.000 #0.028 IVOSOMG_09_03 + #1.000 &
          #0.028 IVOSOMG_10_03 + #1.000 #0.028 IVOSOMG_11_03 + #1.000 #0.000 &
          IVOSOMG_12_03 + #1.000 #0.000 IVOSOMG_13_03 + #1.000 #0.028 IVOSOMG_02_04 + &
          #1.000 #0.028 IVOSOMG_03_04 + #1.000 #0.028 IVOSOMG_04_04 + #1.000 &
          #0.028 IVOSOMG_05_04 + #1.000 #0.028 IVOSOMG_06_04 + #1.000 #0.028 &
          IVOSOMG_07_04 + #1.000 #0.028 IVOSOMG_08_04 + #1.000 #0.028 IVOSOMG_09_04 + &
          #1.000 #0.028 IVOSOMG_10_04 + #1.000 #0.028 IVOSOMG_11_04 + #1.000 &
          #0.000 IVOSOMG_12_04 + #1.000 #0.000 IVOSOMG_13_04 + #1.000 #0.028 &
          IVOSOMG_03_05 + #1.000 #0.028 IVOSOMG_04_05 + #1.000 #0.028 IVOSOMG_05_05 + &
          #1.000 #0.028 IVOSOMG_06_05 + #1.000 #0.028 IVOSOMG_07_05 + #1.000 &
          #0.028 IVOSOMG_08_05 + #1.000 #0.028 IVOSOMG_09_05 + #1.000 #0.028 &
          IVOSOMG_10_05 + #1.000 #0.028 IVOSOMG_11_05 + #1.000 #0.028 IVOSOMG_12_05 + &
          #1.000 #0.000 IVOSOMG_13_05 + #1.000 #0.028 IVOSOMG_03_06 + #1.000 &
          #0.028 IVOSOMG_04_06 + #1.000 #0.028 IVOSOMG_05_06 + #1.000 #0.028 &
          IVOSOMG_06_06 + #1.000 #0.028 IVOSOMG_07_06 + #1.000 #0.028 IVOSOMG_08_06 + &
          #1.000 #0.028 IVOSOMG_09_06 + #1.000 #0.028 IVOSOMG_10_06 + #1.000 &
          #0.028 IVOSOMG_11_06 + #1.000 #0.028 IVOSOMG_12_06 + #1.000 #0.000 &
          IVOSOMG_13_06 + #1.000 #0.028 IVOSOMG_04_07 + #1.000 #0.028 IVOSOMG_05_07 + &
          #1.000 #0.028 IVOSOMG_06_07 + #1.000 #0.028 IVOSOMG_07_07 + #1.000 &
          #0.028 IVOSOMG_08_07 + #1.000 #0.028 IVOSOMG_09_07 + #1.000 #0.028 &
          IVOSOMG_10_07 + #1.000 #0.028 IVOSOMG_11_07 + #1.000 #0.028 IVOSOMG_12_07 + &
          #1.000 #0.028 IVOSOMG_13_07
 S340.6 IVOSOMG_15_01 + OH = #0.210 #0.361 IVOSOMG_15_02 + #0.210 #0.007 &
          IVOSOMG_15_03 + #0.210 #0.606 IVOSOMG_15_04 + #0.210 #0.026 IVOSOMG_15_05 + &
          #0.790 #0.071 IVOSOMG_01_01 + #0.790 #0.071 IVOSOMG_02_01 + #0.790 &
          #0.071 IVOSOMG_03_01 + #0.790 #0.071 IVOSOMG_04_01 + #0.790 #0.071 &
          IVOSOMG_05_01 + #0.790 #0.071 IVOSOMG_06_01 + #0.790 #0.071 IVOSOMG_07_01 + &
          #0.790 #0.071 IVOSOMG_08_01 + #0.790 #0.071 IVOSOMG_09_01 + #0.790 &
          #0.071 IVOSOMG_10_01 + #0.790 #0.071 IVOSOMG_11_01 + #0.790 #0.071 &
          IVOSOMG_12_01 + #0.790 #0.071 IVOSOMG_13_01 + #0.790 #0.071 IVOSOMG_14_01 + &
          #0.790 #0.071 IVOSOMG_01_02 + #0.790 #0.071 IVOSOMG_02_02 + #0.790 &
          #0.071 IVOSOMG_03_02 + #0.790 #0.071 IVOSOMG_04_02 + #0.790 #0.071 &
          IVOSOMG_05_02 + #0.790 #0.071 IVOSOMG_06_02 + #0.790 #0.071 IVOSOMG_07_02 + &
          #0.790 #0.071 IVOSOMG_08_02 + #0.790 #0.071 IVOSOMG_09_02 + #0.790 &
          #0.071 IVOSOMG_10_02 + #0.790 #0.071 IVOSOMG_11_02 + #0.790 #0.071 &
          IVOSOMG_12_02 + #0.790 #0.071 IVOSOMG_13_02 + #0.790 #0.071 IVOSOMG_14_02
 S341.6 IVOSOMG_15_02 + OH = #0.161 #0.361 IVOSOMG_15_03 + #0.161 #0.007 &
          IVOSOMG_15_04 + #0.161 #0.606 IVOSOMG_15_05 + #0.161 #0.026 IVOSOMG_15_06 + &
          #0.839 #0.050 IVOSOMG_01_01 + #0.839 #0.050 IVOSOMG_02_01 + #0.839 &
          #0.050 IVOSOMG_03_01 + #0.839 #0.050 IVOSOMG_04_01 + #0.839 #0.050 &
          IVOSOMG_05_01 + #0.839 #0.050 IVOSOMG_06_01 + #0.839 #0.050 IVOSOMG_07_01 + &
          #0.839 #0.050 IVOSOMG_08_01 + #0.839 #0.050 IVOSOMG_09_01 + #0.839 &
          #0.050 IVOSOMG_10_01 + #0.839 #0.050 IVOSOMG_11_01 + #0.839 #0.050 &
          IVOSOMG_12_01 + #0.839 #0.050 IVOSOMG_13_01 + #0.839 #0.000 IVOSOMG_14_01 + &
          #0.839 #0.050 IVOSOMG_01_02 + #0.839 #0.050 IVOSOMG_02_02 + #0.839 &
          #0.050 IVOSOMG_03_02 + #0.839 #0.050 IVOSOMG_04_02 + #0.839 #0.050 &
          IVOSOMG_05_02 + #0.839 #0.050 IVOSOMG_06_02 + #0.839 #0.050 IVOSOMG_07_02 + &
          #0.839 #0.050 IVOSOMG_08_02 + #0.839 #0.050 IVOSOMG_09_02 + #0.839 &
          #0.050 IVOSOMG_10_02 + #0.839 #0.050 IVOSOMG_11_02 + #0.839 #0.050 &
          IVOSOMG_12_02 + #0.839 #0.050 IVOSOMG_13_02 + #0.839 #0.050 IVOSOMG_14_02 + &
          #0.839 #0.050 IVOSOMG_02_03 + #0.839 #0.050 IVOSOMG_03_03 + #0.839 &
          #0.050 IVOSOMG_04_03 + #0.839 #0.050 IVOSOMG_05_03 + #0.839 #0.050 &
          IVOSOMG_06_03 + #0.839 #0.050 IVOSOMG_07_03 + #0.839 #0.050 IVOSOMG_08_03 + &
          #0.839 #0.050 IVOSOMG_09_03 + #0.839 #0.050 IVOSOMG_10_03 + #0.839 &
          #0.050 IVOSOMG_11_03 + #0.839 #0.050 IVOSOMG_12_03 + #0.839 #0.050 &
          IVOSOMG_13_03 + #0.839 #0.050 IVOSOMG_14_03
 S342.6 IVOSOMG_15_03 + OH = #0.131 #0.361 IVOSOMG_15_04 + #0.131 #0.007 &
          IVOSOMG_15_05 + #0.131 #0.606 IVOSOMG_15_06 + #0.131 #0.026 IVOSOMG_15_07 + &
          #0.869 #0.038 IVOSOMG_01_01 + #0.869 #0.038 IVOSOMG_02_01 + #0.869 &
          #0.038 IVOSOMG_03_01 + #0.869 #0.038 IVOSOMG_04_01 + #0.869 #0.038 &
          IVOSOMG_05_01 + #0.869 #0.038 IVOSOMG_06_01 + #0.869 #0.038 IVOSOMG_07_01 + &
          #0.869 #0.038 IVOSOMG_08_01 + #0.869 #0.038 IVOSOMG_09_01 + #0.869 &
          #0.038 IVOSOMG_10_01 + #0.869 #0.038 IVOSOMG_11_01 + #0.869 #0.038 &
          IVOSOMG_12_01 + #0.869 #0.038 IVOSOMG_13_01 + #0.869 #0.000 IVOSOMG_14_01 + &
          #0.869 #0.038 IVOSOMG_01_02 + #0.869 #0.038 IVOSOMG_02_02 + #0.869 &
          #0.038 IVOSOMG_03_02 + #0.869 #0.038 IVOSOMG_04_02 + #0.869 #0.038 &
          IVOSOMG_05_02 + #0.869 #0.038 IVOSOMG_06_02 + #0.869 #0.038 IVOSOMG_07_02 + &
          #0.869 #0.038 IVOSOMG_08_02 + #0.869 #0.038 IVOSOMG_09_02 + #0.869 &
          #0.038 IVOSOMG_10_02 + #0.869 #0.038 IVOSOMG_11_02 + #0.869 #0.038 &
          IVOSOMG_12_02 + #0.869 #0.038 IVOSOMG_13_02 + #0.869 #0.000 IVOSOMG_14_02 + &
          #0.869 #0.038 IVOSOMG_02_03 + #0.869 #0.038 IVOSOMG_03_03 + #0.869 &
          #0.038 IVOSOMG_04_03 + #0.869 #0.038 IVOSOMG_05_03 + #0.869 #0.038 &
          IVOSOMG_06_03 + #0.869 #0.038 IVOSOMG_07_03 + #0.869 #0.038 IVOSOMG_08_03 + &
          #0.869 #0.038 IVOSOMG_09_03 + #0.869 #0.038 IVOSOMG_10_03 + #0.869 &
          #0.038 IVOSOMG_11_03 + #0.869 #0.038 IVOSOMG_12_03 + #0.869 #0.038 &
          IVOSOMG_13_03 + #0.869 #0.038 IVOSOMG_14_03 + #0.869 #0.038 IVOSOMG_02_04 + &
          #0.869 #0.038 IVOSOMG_03_04 + #0.869 #0.038 IVOSOMG_04_04 + #0.869 &
          #0.038 IVOSOMG_05_04 + #0.869 #0.038 IVOSOMG_06_04 + #0.869 #0.038 &
          IVOSOMG_07_04 + #0.869 #0.038 IVOSOMG_08_04 + #0.869 #0.038 IVOSOMG_09_04 + &
          #0.869 #0.038 IVOSOMG_10_04 + #0.869 #0.038 IVOSOMG_11_04 + #0.869 &
          #0.038 IVOSOMG_12_04 + #0.869 #0.038 IVOSOMG_13_04 + #0.869 #0.038 &
          IVOSOMG_14_04
 S343.6 IVOSOMG_15_04 + OH = #0.109 #0.361 IVOSOMG_15_05 + #0.109 #0.007 &
          IVOSOMG_15_06 + #0.109 #0.632 IVOSOMG_15_07 + #0.891 #0.032 IVOSOMG_01_01 + &
          #0.891 #0.032 IVOSOMG_02_01 + #0.891 #0.032 IVOSOMG_03_01 + #0.891 &
          #0.032 IVOSOMG_04_01 + #0.891 #0.032 IVOSOMG_05_01 + #0.891 #0.032 &
          IVOSOMG_06_01 + #0.891 #0.032 IVOSOMG_07_01 + #0.891 #0.032 IVOSOMG_08_01 + &
          #0.891 #0.032 IVOSOMG_09_01 + #0.891 #0.032 IVOSOMG_10_01 + #0.891 &
          #0.032 IVOSOMG_11_01 + #0.891 #0.032 IVOSOMG_12_01 + #0.891 #0.000 &
          IVOSOMG_13_01 + #0.891 #0.000 IVOSOMG_14_01 + #0.891 #0.032 IVOSOMG_01_02 + &
          #0.891 #0.032 IVOSOMG_02_02 + #0.891 #0.032 IVOSOMG_03_02 + #0.891 &
          #0.032 IVOSOMG_04_02 + #0.891 #0.032 IVOSOMG_05_02 + #0.891 #0.032 &
          IVOSOMG_06_02 + #0.891 #0.032 IVOSOMG_07_02 + #0.891 #0.032 IVOSOMG_08_02 + &
          #0.891 #0.032 IVOSOMG_09_02 + #0.891 #0.032 IVOSOMG_10_02 + #0.891 &
          #0.032 IVOSOMG_11_02 + #0.891 #0.032 IVOSOMG_12_02 + #0.891 #0.032 &
          IVOSOMG_13_02 + #0.891 #0.000 IVOSOMG_14_02 + #0.891 #0.032 IVOSOMG_02_03 + &
          #0.891 #0.032 IVOSOMG_03_03 + #0.891 #0.032 IVOSOMG_04_03 + #0.891 &
          #0.032 IVOSOMG_05_03 + #0.891 #0.032 IVOSOMG_06_03 + #0.891 #0.032 &
          IVOSOMG_07_03 + #0.891 #0.032 IVOSOMG_08_03 + #0.891 #0.032 IVOSOMG_09_03 + &
          #0.891 #0.032 IVOSOMG_10_03 + #0.891 #0.032 IVOSOMG_11_03 + #0.891 &
          #0.032 IVOSOMG_12_03 + #0.891 #0.032 IVOSOMG_13_03 + #0.891 #0.000 &
          IVOSOMG_14_03 + #0.891 #0.032 IVOSOMG_02_04 + #0.891 #0.032 IVOSOMG_03_04 + &
          #0.891 #0.032 IVOSOMG_04_04 + #0.891 #0.032 IVOSOMG_05_04 + #0.891 &
          #0.032 IVOSOMG_06_04 + #0.891 #0.032 IVOSOMG_07_04 + #0.891 #0.032 &
          IVOSOMG_08_04 + #0.891 #0.032 IVOSOMG_09_04 + #0.891 #0.032 IVOSOMG_10_04 + &
          #0.891 #0.032 IVOSOMG_11_04 + #0.891 #0.032 IVOSOMG_12_04 + #0.891 &
          #0.032 IVOSOMG_13_04 + #0.891 #0.032 IVOSOMG_14_04 + #0.891 #0.032 &
          IVOSOMG_03_05 + #0.891 #0.032 IVOSOMG_04_05 + #0.891 #0.032 IVOSOMG_05_05 + &
          #0.891 #0.032 IVOSOMG_06_05 + #0.891 #0.032 IVOSOMG_07_05 + #0.891 &
          #0.032 IVOSOMG_08_05 + #0.891 #0.032 IVOSOMG_09_05 + #0.891 #0.032 &
          IVOSOMG_10_05 + #0.891 #0.032 IVOSOMG_11_05 + #0.891 #0.032 IVOSOMG_12_05 + &
          #0.891 #0.032 IVOSOMG_13_05 + #0.891 #0.032 IVOSOMG_14_05
 S344.6 IVOSOMG_15_05 + OH = #0.091 #0.361 IVOSOMG_15_06 + #0.091 #0.639 &
          IVOSOMG_15_07 + #0.909 #0.028 IVOSOMG_01_01 + #0.909 #0.028 IVOSOMG_02_01 + &
          #0.909 #0.028 IVOSOMG_03_01 + #0.909 #0.028 IVOSOMG_04_01 + #0.909 &
          #0.028 IVOSOMG_05_01 + #0.909 #0.028 IVOSOMG_06_01 + #0.909 #0.028 &
          IVOSOMG_07_01 + #0.909 #0.028 IVOSOMG_08_01 + #0.909 #0.028 IVOSOMG_09_01 + &
          #0.909 #0.028 IVOSOMG_10_01 + #0.909 #0.028 IVOSOMG_11_01 + #0.909 &
          #0.028 IVOSOMG_12_01 + #0.909 #0.000 IVOSOMG_13_01 + #0.909 #0.000 &
          IVOSOMG_14_01 + #0.909 #0.028 IVOSOMG_01_02 + #0.909 #0.028 IVOSOMG_02_02 + &
          #0.909 #0.028 IVOSOMG_03_02 + #0.909 #0.028 IVOSOMG_04_02 + #0.909 &
          #0.028 IVOSOMG_05_02 + #0.909 #0.028 IVOSOMG_06_02 + #0.909 #0.028 &
          IVOSOMG_07_02 + #0.909 #0.028 IVOSOMG_08_02 + #0.909 #0.028 IVOSOMG_09_02 + &
          #0.909 #0.028 IVOSOMG_10_02 + #0.909 #0.028 IVOSOMG_11_02 + #0.909 &
          #0.028 IVOSOMG_12_02 + #0.909 #0.000 IVOSOMG_13_02 + #0.909 #0.000 &
          IVOSOMG_14_02 + #0.909 #0.028 IVOSOMG_02_03 + #0.909 #0.028 IVOSOMG_03_03 + &
          #0.909 #0.028 IVOSOMG_04_03 + #0.909 #0.028 IVOSOMG_05_03 + #0.909 &
          #0.028 IVOSOMG_06_03 + #0.909 #0.028 IVOSOMG_07_03 + #0.909 #0.028 &
          IVOSOMG_08_03 + #0.909 #0.028 IVOSOMG_09_03 + #0.909 #0.028 IVOSOMG_10_03 + &
          #0.909 #0.028 IVOSOMG_11_03 + #0.909 #0.028 IVOSOMG_12_03 + #0.909 &
          #0.028 IVOSOMG_13_03 + #0.909 #0.000 IVOSOMG_14_03 + #0.909 #0.028 &
          IVOSOMG_02_04 + #0.909 #0.028 IVOSOMG_03_04 + #0.909 #0.028 IVOSOMG_04_04 + &
          #0.909 #0.028 IVOSOMG_05_04 + #0.909 #0.028 IVOSOMG_06_04 + #0.909 &
          #0.028 IVOSOMG_07_04 + #0.909 #0.028 IVOSOMG_08_04 + #0.909 #0.028 &
          IVOSOMG_09_04 + #0.909 #0.028 IVOSOMG_10_04 + #0.909 #0.028 IVOSOMG_11_04 + &
          #0.909 #0.028 IVOSOMG_12_04 + #0.909 #0.028 IVOSOMG_13_04 + #0.909 &
          #0.000 IVOSOMG_14_04 + #0.909 #0.028 IVOSOMG_03_05 + #0.909 #0.028 &
          IVOSOMG_04_05 + #0.909 #0.028 IVOSOMG_05_05 + #0.909 #0.028 IVOSOMG_06_05 + &
          #0.909 #0.028 IVOSOMG_07_05 + #0.909 #0.028 IVOSOMG_08_05 + #0.909 &
          #0.028 IVOSOMG_09_05 + #0.909 #0.028 IVOSOMG_10_05 + #0.909 #0.028 &
          IVOSOMG_11_05 + #0.909 #0.028 IVOSOMG_12_05 + #0.909 #0.028 IVOSOMG_13_05 + &
          #0.909 #0.028 IVOSOMG_14_05 + #0.909 #0.028 IVOSOMG_03_06 + #0.909 &
          #0.028 IVOSOMG_04_06 + #0.909 #0.028 IVOSOMG_05_06 + #0.909 #0.028 &
          IVOSOMG_06_06 + #0.909 #0.028 IVOSOMG_07_06 + #0.909 #0.028 IVOSOMG_08_06 + &
          #0.909 #0.028 IVOSOMG_09_06 + #0.909 #0.028 IVOSOMG_10_06 + #0.909 &
          #0.028 IVOSOMG_11_06 + #0.909 #0.028 IVOSOMG_12_06 + #0.909 #0.028 &
          IVOSOMG_13_06 + #0.909 #0.028 IVOSOMG_14_06
 S345.6 IVOSOMG_15_06 + OH = #0.077 #1.000 IVOSOMG_15_07 + #0.923 #0.025 &
          IVOSOMG_01_01 + #0.923 #0.025 IVOSOMG_02_01 + #0.923 #0.025 IVOSOMG_03_01 + &
          #0.923 #0.025 IVOSOMG_04_01 + #0.923 #0.025 IVOSOMG_05_01 + #0.923 &
          #0.025 IVOSOMG_06_01 + #0.923 #0.025 IVOSOMG_07_01 + #0.923 #0.025 &
          IVOSOMG_08_01 + #0.923 #0.025 IVOSOMG_09_01 + #0.923 #0.025 IVOSOMG_10_01 + &
          #0.923 #0.025 IVOSOMG_11_01 + #0.923 #0.000 IVOSOMG_12_01 + #0.923 &
          #0.000 IVOSOMG_13_01 + #0.923 #0.000 IVOSOMG_14_01 + #0.923 #0.025 &
          IVOSOMG_01_02 + #0.923 #0.025 IVOSOMG_02_02 + #0.923 #0.025 IVOSOMG_03_02 + &
          #0.923 #0.025 IVOSOMG_04_02 + #0.923 #0.025 IVOSOMG_05_02 + #0.923 &
          #0.025 IVOSOMG_06_02 + #0.923 #0.025 IVOSOMG_07_02 + #0.923 #0.025 &
          IVOSOMG_08_02 + #0.923 #0.025 IVOSOMG_09_02 + #0.923 #0.025 IVOSOMG_10_02 + &
          #0.923 #0.025 IVOSOMG_11_02 + #0.923 #0.025 IVOSOMG_12_02 + #0.923 &
          #0.000 IVOSOMG_13_02 + #0.923 #0.000 IVOSOMG_14_02 + #0.923 #0.025 &
          IVOSOMG_02_03 + #0.923 #0.025 IVOSOMG_03_03 + #0.923 #0.025 IVOSOMG_04_03 + &
          #0.923 #0.025 IVOSOMG_05_03 + #0.923 #0.025 IVOSOMG_06_03 + #0.923 &
          #0.025 IVOSOMG_07_03 + #0.923 #0.025 IVOSOMG_08_03 + #0.923 #0.025 &
          IVOSOMG_09_03 + #0.923 #0.025 IVOSOMG_10_03 + #0.923 #0.025 IVOSOMG_11_03 + &
          #0.923 #0.025 IVOSOMG_12_03 + #0.923 #0.000 IVOSOMG_13_03 + #0.923 &
          #0.000 IVOSOMG_14_03 + #0.923 #0.025 IVOSOMG_02_04 + #0.923 #0.025 &
          IVOSOMG_03_04 + #0.923 #0.025 IVOSOMG_04_04 + #0.923 #0.025 IVOSOMG_05_04 + &
          #0.923 #0.025 IVOSOMG_06_04 + #0.923 #0.025 IVOSOMG_07_04 + #0.923 &
          #0.025 IVOSOMG_08_04 + #0.923 #0.025 IVOSOMG_09_04 + #0.923 #0.025 &
          IVOSOMG_10_04 + #0.923 #0.025 IVOSOMG_11_04 + #0.923 #0.025 IVOSOMG_12_04 + &
          #0.923 #0.025 IVOSOMG_13_04 + #0.923 #0.000 IVOSOMG_14_04 + #0.923 &
          #0.025 IVOSOMG_03_05 + #0.923 #0.025 IVOSOMG_04_05 + #0.923 #0.025 &
          IVOSOMG_05_05 + #0.923 #0.025 IVOSOMG_06_05 + #0.923 #0.025 IVOSOMG_07_05 + &
          #0.923 #0.025 IVOSOMG_08_05 + #0.923 #0.025 IVOSOMG_09_05 + #0.923 &
          #0.025 IVOSOMG_10_05 + #0.923 #0.025 IVOSOMG_11_05 + #0.923 #0.025 &
          IVOSOMG_12_05 + #0.923 #0.025 IVOSOMG_13_05 + #0.923 #0.000 IVOSOMG_14_05 + &
          #0.923 #0.025 IVOSOMG_03_06 + #0.923 #0.025 IVOSOMG_04_06 + #0.923 &
          #0.025 IVOSOMG_05_06 + #0.923 #0.025 IVOSOMG_06_06 + #0.923 #0.025 &
          IVOSOMG_07_06 + #0.923 #0.025 IVOSOMG_08_06 + #0.923 #0.025 IVOSOMG_09_06 + &
          #0.923 #0.025 IVOSOMG_10_06 + #0.923 #0.025 IVOSOMG_11_06 + #0.923 &
          #0.025 IVOSOMG_12_06 + #0.923 #0.025 IVOSOMG_13_06 + #0.923 #0.025 &
          IVOSOMG_14_06 + #0.923 #0.025 IVOSOMG_04_07 + #0.923 #0.025 IVOSOMG_05_07 + &
          #0.923 #0.025 IVOSOMG_06_07 + #0.923 #0.025 IVOSOMG_07_07 + #0.923 &
          #0.025 IVOSOMG_08_07 + #0.923 #0.025 IVOSOMG_09_07 + #0.923 #0.025 &
          IVOSOMG_10_07 + #0.923 #0.025 IVOSOMG_11_07 + #0.923 #0.025 IVOSOMG_12_07 + &
          #0.923 #0.025 IVOSOMG_13_07 + #0.923 #0.025 IVOSOMG_14_07
 S346.6 IVOSOMG_15_07 + OH = #1.000 #0.026 IVOSOMG_01_01 + #1.000 #0.026 &
          IVOSOMG_02_01 + #1.000 #0.026 IVOSOMG_03_01 + #1.000 #0.026 IVOSOMG_04_01 + &
          #1.000 #0.026 IVOSOMG_05_01 + #1.000 #0.026 IVOSOMG_06_01 + #1.000 &
          #0.026 IVOSOMG_07_01 + #1.000 #0.026 IVOSOMG_08_01 + #1.000 #0.026 &
          IVOSOMG_09_01 + #1.000 #0.026 IVOSOMG_10_01 + #1.000 #0.026 IVOSOMG_11_01 + &
          #1.000 #0.000 IVOSOMG_12_01 + #1.000 #0.000 IVOSOMG_13_01 + #1.000 &
          #0.000 IVOSOMG_14_01 + #1.000 #0.026 IVOSOMG_01_02 + #1.000 #0.026 &
          IVOSOMG_02_02 + #1.000 #0.026 IVOSOMG_03_02 + #1.000 #0.026 IVOSOMG_04_02 + &
          #1.000 #0.026 IVOSOMG_05_02 + #1.000 #0.026 IVOSOMG_06_02 + #1.000 &
          #0.026 IVOSOMG_07_02 + #1.000 #0.026 IVOSOMG_08_02 + #1.000 #0.026 &
          IVOSOMG_09_02 + #1.000 #0.026 IVOSOMG_10_02 + #1.000 #0.026 IVOSOMG_11_02 + &
          #1.000 #0.000 IVOSOMG_12_02 + #1.000 #0.000 IVOSOMG_13_02 + #1.000 &
          #0.000 IVOSOMG_14_02 + #1.000 #0.026 IVOSOMG_02_03 + #1.000 #0.026 &
          IVOSOMG_03_03 + #1.000 #0.026 IVOSOMG_04_03 + #1.000 #0.026 IVOSOMG_05_03 + &
          #1.000 #0.026 IVOSOMG_06_03 + #1.000 #0.026 IVOSOMG_07_03 + #1.000 &
          #0.026 IVOSOMG_08_03 + #1.000 #0.026 IVOSOMG_09_03 + #1.000 #0.026 &
          IVOSOMG_10_03 + #1.000 #0.026 IVOSOMG_11_03 + #1.000 #0.026 IVOSOMG_12_03 + &
          #1.000 #0.000 IVOSOMG_13_03 + #1.000 #0.000 IVOSOMG_14_03 + #1.000 &
          #0.026 IVOSOMG_02_04 + #1.000 #0.026 IVOSOMG_03_04 + #1.000 #0.026 &
          IVOSOMG_04_04 + #1.000 #0.026 IVOSOMG_05_04 + #1.000 #0.026 IVOSOMG_06_04 + &
          #1.000 #0.026 IVOSOMG_07_04 + #1.000 #0.026 IVOSOMG_08_04 + #1.000 &
          #0.026 IVOSOMG_09_04 + #1.000 #0.026 IVOSOMG_10_04 + #1.000 #0.026 &
          IVOSOMG_11_04 + #1.000 #0.026 IVOSOMG_12_04 + #1.000 #0.000 IVOSOMG_13_04 + &
          #1.000 #0.000 IVOSOMG_14_04 + #1.000 #0.026 IVOSOMG_03_05 + #1.000 &
          #0.026 IVOSOMG_04_05 + #1.000 #0.026 IVOSOMG_05_05 + #1.000 #0.026 &
          IVOSOMG_06_05 + #1.000 #0.026 IVOSOMG_07_05 + #1.000 #0.026 IVOSOMG_08_05 + &
          #1.000 #0.026 IVOSOMG_09_05 + #1.000 #0.026 IVOSOMG_10_05 + #1.000 &
          #0.026 IVOSOMG_11_05 + #1.000 #0.026 IVOSOMG_12_05 + #1.000 #0.026 &
          IVOSOMG_13_05 + #1.000 #0.000 IVOSOMG_14_05 + #1.000 #0.026 IVOSOMG_03_06 + &
          #1.000 #0.026 IVOSOMG_04_06 + #1.000 #0.026 IVOSOMG_05_06 + #1.000 &
          #0.026 IVOSOMG_06_06 + #1.000 #0.026 IVOSOMG_07_06 + #1.000 #0.026 &
          IVOSOMG_08_06 + #1.000 #0.026 IVOSOMG_09_06 + #1.000 #0.026 IVOSOMG_10_06 + &
          #1.000 #0.026 IVOSOMG_11_06 + #1.000 #0.026 IVOSOMG_12_06 + #1.000 &
          #0.026 IVOSOMG_13_06 + #1.000 #0.000 IVOSOMG_14_06 + #1.000 #0.026 &
          IVOSOMG_04_07 + #1.000 #0.026 IVOSOMG_05_07 + #1.000 #0.026 IVOSOMG_06_07 + &
          #1.000 #0.026 IVOSOMG_07_07 + #1.000 #0.026 IVOSOMG_08_07 + #1.000 &
          #0.026 IVOSOMG_09_07 + #1.000 #0.026 IVOSOMG_10_07 + #1.000 #0.026 &
          IVOSOMG_11_07 + #1.000 #0.026 IVOSOMG_12_07 + #1.000 #0.026 IVOSOMG_13_07 + &
          #1.000 #0.026 IVOSOMG_14_07
 BL28   SVOC + OH = #0.035 #0.361 SVOSOMG_06_05 + #0.035 #0.007 SVOSOMG_06_06 + &
          #0.035 #0.632 SVOSOMG_06_07 + #0.965 #0.048 SVOSOMG_01_01 + #0.965 &
          #0.048 SVOSOMG_02_01 + #0.965 #0.048 SVOSOMG_03_01 + #0.965 #0.048 &
          SVOSOMG_04_01 + #0.965 #0.048 SVOSOMG_05_01 + #0.965 #0.048 SVOSOMG_01_02 + &
          #0.965 #0.048 SVOSOMG_02_02 + #0.965 #0.048 SVOSOMG_03_02 + #0.965 &
          #0.048 SVOSOMG_04_02 + #0.965 #0.048 SVOSOMG_05_02 + #0.965 #0.048 &
          SVOSOMG_02_03 + #0.965 #0.048 SVOSOMG_03_03 + #0.965 #0.048 SVOSOMG_04_03 + &
          #0.965 #0.048 SVOSOMG_05_03 + #0.965 #0.048 SVOSOMG_02_04 + #0.965 &
          #0.048 SVOSOMG_03_04 + #0.965 #0.048 SVOSOMG_04_04 + #0.965 #0.048 &
          SVOSOMG_05_04 + #0.965 #0.048 SVOSOMG_03_05 + #0.965 #0.048 SVOSOMG_04_05 + &
          #0.965 #0.048 SVOSOMG_05_05 + SVOC
 S347.7 SVOSOMG_02_01 + OH = #0.059 #0.361 SVOSOMG_02_02 + #0.059 #0.007 &
          SVOSOMG_02_03 + #0.059 #0.632 SVOSOMG_02_04 + #0.941 #1.000 SVOSOMG_01_01 + &
          #0.941 #1.000 SVOSOMG_01_02
 S348.7 SVOSOMG_02_02 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #2.000 &
          SVOSOMG_01_02
 S349.7 SVOSOMG_02_03 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_01_02
 S350.7 SVOSOMG_02_04 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_01_02
 S351.7 SVOSOMG_03_01 + OH = #0.091 #0.361 SVOSOMG_03_02 + #0.091 #0.007 &
          SVOSOMG_03_03 + #0.091 #0.606 SVOSOMG_03_04 + #0.091 #0.026 SVOSOMG_03_05 + &
          #0.909 #0.500 SVOSOMG_01_01 + #0.909 #0.500 SVOSOMG_02_01 + #0.909 &
          #0.500 SVOSOMG_01_02 + #0.909 #0.500 SVOSOMG_02_02
 S352.7 SVOSOMG_03_02 + OH = #0.035 #0.361 SVOSOMG_03_03 + #0.035 #0.007 &
          SVOSOMG_03_04 + #0.035 #0.606 SVOSOMG_03_05 + #0.035 #0.026 SVOSOMG_03_06 + &
          #0.965 #0.500 SVOSOMG_01_01 + #0.965 #0.000 SVOSOMG_02_01 + #0.965 &
          #0.500 SVOSOMG_01_02 + #0.965 #0.500 SVOSOMG_02_02 + #0.965 #0.500 &
          SVOSOMG_02_03
 S353.7 SVOSOMG_03_03 + OH = #1.000 #0.500 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.500 SVOSOMG_01_02 + #1.000 #0.000 SVOSOMG_02_02 + &
          #1.000 #0.500 SVOSOMG_02_03 + #1.000 #0.500 SVOSOMG_02_04
 S354.7 SVOSOMG_03_04 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #1.000 SVOSOMG_01_02 + #1.000 #0.000 SVOSOMG_02_02 + &
          #1.000 #0.000 SVOSOMG_02_03 + #1.000 #1.000 SVOSOMG_02_04
 S355.7 SVOSOMG_03_05 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_01_02 + #1.000 #0.000 SVOSOMG_02_02 + &
          #1.000 #0.000 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_02_04
 S356.7 SVOSOMG_03_06 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_01_02 + #1.000 #0.000 SVOSOMG_02_02 + &
          #1.000 #0.000 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_02_04
 S357.7 SVOSOMG_04_01 + OH = #0.114 #0.361 SVOSOMG_04_02 + #0.114 #0.007 &
          SVOSOMG_04_03 + #0.114 #0.606 SVOSOMG_04_04 + #0.114 #0.026 SVOSOMG_04_05 + &
          #0.886 #0.333 SVOSOMG_01_01 + #0.886 #0.333 SVOSOMG_02_01 + #0.886 &
          #0.333 SVOSOMG_03_01 + #0.886 #0.333 SVOSOMG_01_02 + #0.886 #0.333 &
          SVOSOMG_02_02 + #0.886 #0.333 SVOSOMG_03_02
 S358.7 SVOSOMG_04_02 + OH = #0.059 #0.361 SVOSOMG_04_03 + #0.059 #0.007 &
          SVOSOMG_04_04 + #0.059 #0.606 SVOSOMG_04_05 + #0.059 #0.026 SVOSOMG_04_06 + &
          #0.941 #0.286 SVOSOMG_01_01 + #0.941 #0.286 SVOSOMG_02_01 + #0.941 &
          #0.000 SVOSOMG_03_01 + #0.941 #0.286 SVOSOMG_01_02 + #0.941 #0.286 &
          SVOSOMG_02_02 + #0.941 #0.286 SVOSOMG_03_02 + #0.941 #0.286 SVOSOMG_02_03 + &
          #0.941 #0.286 SVOSOMG_03_03
 S359.7 SVOSOMG_04_03 + OH = #0.025 #0.361 SVOSOMG_04_04 + #0.025 #0.007 &
          SVOSOMG_04_05 + #0.025 #0.606 SVOSOMG_04_06 + #0.025 #0.026 SVOSOMG_04_07 + &
          #0.975 #0.250 SVOSOMG_01_01 + #0.975 #0.250 SVOSOMG_02_01 + #0.975 &
          #0.000 SVOSOMG_03_01 + #0.975 #0.250 SVOSOMG_01_02 + #0.975 #0.250 &
          SVOSOMG_02_02 + #0.975 #0.000 SVOSOMG_03_02 + #0.975 #0.250 SVOSOMG_02_03 + &
          #0.975 #0.250 SVOSOMG_03_03 + #0.975 #0.250 SVOSOMG_02_04 + #0.975 &
          #0.250 SVOSOMG_03_04
 S360.7 SVOSOMG_04_04 + OH = #1.000 #0.286 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.286 SVOSOMG_01_02 + &
          #1.000 #0.286 SVOSOMG_02_02 + #1.000 #0.000 SVOSOMG_03_02 + #1.000 &
          #0.286 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.286 &
          SVOSOMG_02_04 + #1.000 #0.286 SVOSOMG_03_04 + #1.000 #0.286 SVOSOMG_03_05
 S361.7 SVOSOMG_04_05 + OH = #1.000 #0.333 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.333 SVOSOMG_01_02 + &
          #1.000 #0.000 SVOSOMG_02_02 + #1.000 #0.000 SVOSOMG_03_02 + #1.000 &
          #0.333 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.333 &
          SVOSOMG_02_04 + #1.000 #0.000 SVOSOMG_03_04 + #1.000 #0.333 SVOSOMG_03_05 + &
          #1.000 #0.333 SVOSOMG_03_06
 S362.7 SVOSOMG_04_06 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.667 SVOSOMG_01_02 + &
          #1.000 #0.000 SVOSOMG_02_02 + #1.000 #0.000 SVOSOMG_03_02 + #1.000 &
          #0.000 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.667 &
          SVOSOMG_02_04 + #1.000 #0.000 SVOSOMG_03_04 + #1.000 #0.000 SVOSOMG_03_05 + &
          #1.000 #0.667 SVOSOMG_03_06
 S363.7 SVOSOMG_04_07 + OH = #1.000 #0.000 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_01_02 + &
          #1.000 #0.000 SVOSOMG_02_02 + #1.000 #0.000 SVOSOMG_03_02 + #1.000 &
          #0.000 SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.000 &
          SVOSOMG_02_04 + #1.000 #0.000 SVOSOMG_03_04 + #1.000 #0.000 SVOSOMG_03_05 + &
          #1.000 #0.000 SVOSOMG_03_06
 S364.7 SVOSOMG_05_01 + OH = #0.131 #0.361 SVOSOMG_05_02 + #0.131 #0.007 &
          SVOSOMG_05_03 + #0.131 #0.606 SVOSOMG_05_04 + #0.131 #0.026 SVOSOMG_05_05 + &
          #0.869 #0.250 SVOSOMG_01_01 + #0.869 #0.250 SVOSOMG_02_01 + #0.869 &
          #0.250 SVOSOMG_03_01 + #0.869 #0.250 SVOSOMG_04_01 + #0.869 #0.250 &
          SVOSOMG_01_02 + #0.869 #0.250 SVOSOMG_02_02 + #0.869 #0.250 SVOSOMG_03_02 + &
          #0.869 #0.250 SVOSOMG_04_02
 S365.7 SVOSOMG_05_02 + OH = #0.077 #0.361 SVOSOMG_05_03 + #0.077 #0.007 &
          SVOSOMG_05_04 + #0.077 #0.606 SVOSOMG_05_05 + #0.077 #0.026 SVOSOMG_05_06 + &
          #0.923 #0.182 SVOSOMG_01_01 + #0.923 #0.182 SVOSOMG_02_01 + #0.923 &
          #0.182 SVOSOMG_03_01 + #0.923 #0.000 SVOSOMG_04_01 + #0.923 #0.182 &
          SVOSOMG_01_02 + #0.923 #0.182 SVOSOMG_02_02 + #0.923 #0.182 SVOSOMG_03_02 + &
          #0.923 #0.182 SVOSOMG_04_02 + #0.923 #0.182 SVOSOMG_02_03 + #0.923 &
          #0.182 SVOSOMG_03_03 + #0.923 #0.182 SVOSOMG_04_03
 S366.7 SVOSOMG_05_03 + OH = #0.043 #0.361 SVOSOMG_05_04 + #0.043 #0.007 &
          SVOSOMG_05_05 + #0.043 #0.606 SVOSOMG_05_06 + #0.043 #0.026 SVOSOMG_05_07 + &
          #0.957 #0.154 SVOSOMG_01_01 + #0.957 #0.154 SVOSOMG_02_01 + #0.957 &
          #0.154 SVOSOMG_03_01 + #0.957 #0.000 SVOSOMG_04_01 + #0.957 #0.154 &
          SVOSOMG_01_02 + #0.957 #0.154 SVOSOMG_02_02 + #0.957 #0.154 SVOSOMG_03_02 + &
          #0.957 #0.000 SVOSOMG_04_02 + #0.957 #0.154 SVOSOMG_02_03 + #0.957 &
          #0.154 SVOSOMG_03_03 + #0.957 #0.154 SVOSOMG_04_03 + #0.957 #0.154 &
          SVOSOMG_02_04 + #0.957 #0.154 SVOSOMG_03_04 + #0.957 #0.154 SVOSOMG_04_04
 S367.7 SVOSOMG_05_04 + OH = #0.019 #0.361 SVOSOMG_05_05 + #0.019 #0.007 &
          SVOSOMG_05_06 + #0.019 #0.632 SVOSOMG_05_07 + #0.981 #0.143 SVOSOMG_01_01 + &
          #0.981 #0.143 SVOSOMG_02_01 + #0.981 #0.000 SVOSOMG_03_01 + #0.981 &
          #0.000 SVOSOMG_04_01 + #0.981 #0.143 SVOSOMG_01_02 + #0.981 #0.143 &
          SVOSOMG_02_02 + #0.981 #0.143 SVOSOMG_03_02 + #0.981 #0.000 SVOSOMG_04_02 + &
          #0.981 #0.143 SVOSOMG_02_03 + #0.981 #0.143 SVOSOMG_03_03 + #0.981 &
          #0.000 SVOSOMG_04_03 + #0.981 #0.143 SVOSOMG_02_04 + #0.981 #0.143 &
          SVOSOMG_03_04 + #0.981 #0.143 SVOSOMG_04_04 + #0.981 #0.143 SVOSOMG_03_05 + &
          #0.981 #0.143 SVOSOMG_04_05
 S368.7 SVOSOMG_05_05 + OH = #1.000 #0.143 SVOSOMG_01_01 + #1.000 #0.143 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.143 SVOSOMG_01_02 + #1.000 #0.143 SVOSOMG_02_02 + #1.000 &
          #0.000 SVOSOMG_03_02 + #1.000 #0.000 SVOSOMG_04_02 + #1.000 #0.143 &
          SVOSOMG_02_03 + #1.000 #0.143 SVOSOMG_03_03 + #1.000 #0.000 SVOSOMG_04_03 + &
          #1.000 #0.143 SVOSOMG_02_04 + #1.000 #0.143 SVOSOMG_03_04 + #1.000 &
          #0.000 SVOSOMG_04_04 + #1.000 #0.143 SVOSOMG_03_05 + #1.000 #0.143 &
          SVOSOMG_04_05 + #1.000 #0.143 SVOSOMG_03_06 + #1.000 #0.143 SVOSOMG_04_06
 S369.7 SVOSOMG_05_06 + OH = #1.000 #0.154 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.154 SVOSOMG_01_02 + #1.000 #0.154 SVOSOMG_02_02 + #1.000 &
          #0.000 SVOSOMG_03_02 + #1.000 #0.000 SVOSOMG_04_02 + #1.000 #0.154 &
          SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.000 SVOSOMG_04_03 + &
          #1.000 #0.154 SVOSOMG_02_04 + #1.000 #0.154 SVOSOMG_03_04 + #1.000 &
          #0.000 SVOSOMG_04_04 + #1.000 #0.154 SVOSOMG_03_05 + #1.000 #0.000 &
          SVOSOMG_04_05 + #1.000 #0.154 SVOSOMG_03_06 + #1.000 #0.154 SVOSOMG_04_06 + &
          #1.000 #0.154 SVOSOMG_04_07
 S370.7 SVOSOMG_05_07 + OH = #1.000 #0.200 SVOSOMG_01_01 + #1.000 #0.000 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.200 SVOSOMG_01_02 + #1.000 #0.000 SVOSOMG_02_02 + #1.000 &
          #0.000 SVOSOMG_03_02 + #1.000 #0.000 SVOSOMG_04_02 + #1.000 #0.200 &
          SVOSOMG_02_03 + #1.000 #0.000 SVOSOMG_03_03 + #1.000 #0.000 SVOSOMG_04_03 + &
          #1.000 #0.200 SVOSOMG_02_04 + #1.000 #0.000 SVOSOMG_03_04 + #1.000 &
          #0.000 SVOSOMG_04_04 + #1.000 #0.200 SVOSOMG_03_05 + #1.000 #0.000 &
          SVOSOMG_04_05 + #1.000 #0.200 SVOSOMG_03_06 + #1.000 #0.000 SVOSOMG_04_06 + &
          #1.000 #0.200 SVOSOMG_04_07
 S371.7 SVOSOMG_06_01 + OH = #0.144 #0.361 SVOSOMG_06_02 + #0.144 #0.007 &
          SVOSOMG_06_03 + #0.144 #0.606 SVOSOMG_06_04 + #0.144 #0.026 SVOSOMG_06_05 + &
          #0.856 #0.200 SVOSOMG_01_01 + #0.856 #0.200 SVOSOMG_02_01 + #0.856 &
          #0.200 SVOSOMG_03_01 + #0.856 #0.200 SVOSOMG_04_01 + #0.856 #0.200 &
          SVOSOMG_05_01 + #0.856 #0.200 SVOSOMG_01_02 + #0.856 #0.200 SVOSOMG_02_02 + &
          #0.856 #0.200 SVOSOMG_03_02 + #0.856 #0.200 SVOSOMG_04_02 + #0.856 &
          #0.200 SVOSOMG_05_02
 S372.7 SVOSOMG_06_02 + OH = #0.091 #0.361 SVOSOMG_06_03 + #0.091 #0.007 &
          SVOSOMG_06_04 + #0.091 #0.606 SVOSOMG_06_05 + #0.091 #0.026 SVOSOMG_06_06 + &
          #0.909 #0.154 SVOSOMG_01_01 + #0.909 #0.154 SVOSOMG_02_01 + #0.909 &
          #0.154 SVOSOMG_03_01 + #0.909 #0.154 SVOSOMG_04_01 + #0.909 #0.000 &
          SVOSOMG_05_01 + #0.909 #0.154 SVOSOMG_01_02 + #0.909 #0.154 SVOSOMG_02_02 + &
          #0.909 #0.154 SVOSOMG_03_02 + #0.909 #0.154 SVOSOMG_04_02 + #0.909 &
          #0.154 SVOSOMG_05_02 + #0.909 #0.154 SVOSOMG_02_03 + #0.909 #0.154 &
          SVOSOMG_03_03 + #0.909 #0.154 SVOSOMG_04_03 + #0.909 #0.154 SVOSOMG_05_03
 S373.7 SVOSOMG_06_03 + OH = #0.059 #0.361 SVOSOMG_06_04 + #0.059 #0.007 &
          SVOSOMG_06_05 + #0.059 #0.606 SVOSOMG_06_06 + #0.059 #0.026 SVOSOMG_06_07 + &
          #0.941 #0.125 SVOSOMG_01_01 + #0.941 #0.125 SVOSOMG_02_01 + #0.941 &
          #0.125 SVOSOMG_03_01 + #0.941 #0.125 SVOSOMG_04_01 + #0.941 #0.000 &
          SVOSOMG_05_01 + #0.941 #0.125 SVOSOMG_01_02 + #0.941 #0.125 SVOSOMG_02_02 + &
          #0.941 #0.125 SVOSOMG_03_02 + #0.941 #0.125 SVOSOMG_04_02 + #0.941 &
          #0.000 SVOSOMG_05_02 + #0.941 #0.125 SVOSOMG_02_03 + #0.941 #0.125 &
          SVOSOMG_03_03 + #0.941 #0.125 SVOSOMG_04_03 + #0.941 #0.125 SVOSOMG_05_03 + &
          #0.941 #0.125 SVOSOMG_02_04 + #0.941 #0.125 SVOSOMG_03_04 + #0.941 &
          #0.125 SVOSOMG_04_04 + #0.941 #0.125 SVOSOMG_05_04
 S374.7 SVOSOMG_06_05 + OH = #0.016 #0.361 SVOSOMG_06_06 + #0.016 #0.639 &
          SVOSOMG_06_07 + #0.984 #0.105 SVOSOMG_01_01 + #0.984 #0.105 SVOSOMG_02_01 + &
          #0.984 #0.105 SVOSOMG_03_01 + #0.984 #0.000 SVOSOMG_04_01 + #0.984 &
          #0.000 SVOSOMG_05_01 + #0.984 #0.105 SVOSOMG_01_02 + #0.984 #0.105 &
          SVOSOMG_02_02 + #0.984 #0.105 SVOSOMG_03_02 + #0.984 #0.000 SVOSOMG_04_02 + &
          #0.984 #0.000 SVOSOMG_05_02 + #0.984 #0.105 SVOSOMG_02_03 + #0.984 &
          #0.105 SVOSOMG_03_03 + #0.984 #0.105 SVOSOMG_04_03 + #0.984 #0.000 &
          SVOSOMG_05_03 + #0.984 #0.105 SVOSOMG_02_04 + #0.984 #0.105 SVOSOMG_03_04 + &
          #0.984 #0.105 SVOSOMG_04_04 + #0.984 #0.000 SVOSOMG_05_04 + #0.984 &
          #0.105 SVOSOMG_03_05 + #0.984 #0.105 SVOSOMG_04_05 + #0.984 #0.105 &
          SVOSOMG_05_05 + #0.984 #0.105 SVOSOMG_03_06 + #0.984 #0.105 SVOSOMG_04_06 + &
          #0.984 #0.105 SVOSOMG_05_06
 S375.7 SVOSOMG_06_06 + OH = #1.000 #0.105 SVOSOMG_01_01 + #1.000 #0.105 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.000 SVOSOMG_05_01 + #1.000 #0.105 SVOSOMG_01_02 + #1.000 &
          #0.105 SVOSOMG_02_02 + #1.000 #0.105 SVOSOMG_03_02 + #1.000 #0.000 &
          SVOSOMG_04_02 + #1.000 #0.000 SVOSOMG_05_02 + #1.000 #0.105 SVOSOMG_02_03 + &
          #1.000 #0.105 SVOSOMG_03_03 + #1.000 #0.000 SVOSOMG_04_03 + #1.000 &
          #0.000 SVOSOMG_05_03 + #1.000 #0.105 SVOSOMG_02_04 + #1.000 #0.105 &
          SVOSOMG_03_04 + #1.000 #0.105 SVOSOMG_04_04 + #1.000 #0.000 SVOSOMG_05_04 + &
          #1.000 #0.105 SVOSOMG_03_05 + #1.000 #0.105 SVOSOMG_04_05 + #1.000 &
          #0.000 SVOSOMG_05_05 + #1.000 #0.105 SVOSOMG_03_06 + #1.000 #0.105 &
          SVOSOMG_04_06 + #1.000 #0.105 SVOSOMG_05_06 + #1.000 #0.105 SVOSOMG_04_07 + &
          #1.000 #0.105 SVOSOMG_05_07
 S376.7 SVOSOMG_06_07 + OH = #1.000 #0.125 SVOSOMG_01_01 + #1.000 #0.125 &
          SVOSOMG_02_01 + #1.000 #0.000 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.000 SVOSOMG_05_01 + #1.000 #0.125 SVOSOMG_01_02 + #1.000 &
          #0.125 SVOSOMG_02_02 + #1.000 #0.000 SVOSOMG_03_02 + #1.000 #0.000 &
          SVOSOMG_04_02 + #1.000 #0.000 SVOSOMG_05_02 + #1.000 #0.125 SVOSOMG_02_03 + &
          #1.000 #0.125 SVOSOMG_03_03 + #1.000 #0.000 SVOSOMG_04_03 + #1.000 &
          #0.000 SVOSOMG_05_03 + #1.000 #0.125 SVOSOMG_02_04 + #1.000 #0.125 &
          SVOSOMG_03_04 + #1.000 #0.000 SVOSOMG_04_04 + #1.000 #0.000 SVOSOMG_05_04 + &
          #1.000 #0.125 SVOSOMG_03_05 + #1.000 #0.125 SVOSOMG_04_05 + #1.000 &
          #0.000 SVOSOMG_05_05 + #1.000 #0.125 SVOSOMG_03_06 + #1.000 #0.125 &
          SVOSOMG_04_06 + #1.000 #0.000 SVOSOMG_05_06 + #1.000 #0.125 SVOSOMG_04_07 + &
          #1.000 #0.125 SVOSOMG_05_07
 S377.7 SVOSOMG_07_01 + OH = #0.156 #0.361 SVOSOMG_07_02 + #0.156 #0.007 &
          SVOSOMG_07_03 + #0.156 #0.606 SVOSOMG_07_04 + #0.156 #0.026 SVOSOMG_07_05 + &
          #0.844 #0.167 SVOSOMG_01_01 + #0.844 #0.167 SVOSOMG_02_01 + #0.844 &
          #0.167 SVOSOMG_03_01 + #0.844 #0.167 SVOSOMG_04_01 + #0.844 #0.167 &
          SVOSOMG_05_01 + #0.844 #0.167 SVOSOMG_06_01 + #0.844 #0.167 SVOSOMG_01_02 + &
          #0.844 #0.167 SVOSOMG_02_02 + #0.844 #0.167 SVOSOMG_03_02 + #0.844 &
          #0.167 SVOSOMG_04_02 + #0.844 #0.167 SVOSOMG_05_02 + #0.844 #0.167 &
          SVOSOMG_06_02
 S378.7 SVOSOMG_07_02 + OH = #0.103 #0.361 SVOSOMG_07_03 + #0.103 #0.007 &
          SVOSOMG_07_04 + #0.103 #0.606 SVOSOMG_07_05 + #0.103 #0.026 SVOSOMG_07_06 + &
          #0.897 #0.125 SVOSOMG_01_01 + #0.897 #0.125 SVOSOMG_02_01 + #0.897 &
          #0.125 SVOSOMG_03_01 + #0.897 #0.125 SVOSOMG_04_01 + #0.897 #0.125 &
          SVOSOMG_05_01 + #0.897 #0.000 SVOSOMG_06_01 + #0.897 #0.125 SVOSOMG_01_02 + &
          #0.897 #0.125 SVOSOMG_02_02 + #0.897 #0.125 SVOSOMG_03_02 + #0.897 &
          #0.125 SVOSOMG_04_02 + #0.897 #0.125 SVOSOMG_05_02 + #0.897 #0.125 &
          SVOSOMG_06_02 + #0.897 #0.125 SVOSOMG_02_03 + #0.897 #0.125 SVOSOMG_03_03 + &
          #0.897 #0.125 SVOSOMG_04_03 + #0.897 #0.125 SVOSOMG_05_03 + #0.897 &
          #0.125 SVOSOMG_06_03
 S379.7 SVOSOMG_07_03 + OH = #0.071 #0.361 SVOSOMG_07_04 + #0.071 #0.007 &
          SVOSOMG_07_05 + #0.071 #0.606 SVOSOMG_07_06 + #0.071 #0.026 SVOSOMG_07_07 + &
          #0.929 #0.100 SVOSOMG_01_01 + #0.929 #0.100 SVOSOMG_02_01 + #0.929 &
          #0.100 SVOSOMG_03_01 + #0.929 #0.100 SVOSOMG_04_01 + #0.929 #0.100 &
          SVOSOMG_05_01 + #0.929 #0.000 SVOSOMG_06_01 + #0.929 #0.100 SVOSOMG_01_02 + &
          #0.929 #0.100 SVOSOMG_02_02 + #0.929 #0.100 SVOSOMG_03_02 + #0.929 &
          #0.100 SVOSOMG_04_02 + #0.929 #0.100 SVOSOMG_05_02 + #0.929 #0.000 &
          SVOSOMG_06_02 + #0.929 #0.100 SVOSOMG_02_03 + #0.929 #0.100 SVOSOMG_03_03 + &
          #0.929 #0.100 SVOSOMG_04_03 + #0.929 #0.100 SVOSOMG_05_03 + #0.929 &
          #0.100 SVOSOMG_06_03 + #0.929 #0.100 SVOSOMG_02_04 + #0.929 #0.100 &
          SVOSOMG_03_04 + #0.929 #0.100 SVOSOMG_04_04 + #0.929 #0.100 SVOSOMG_05_04 + &
          #0.929 #0.100 SVOSOMG_06_04
 S380.7 SVOSOMG_07_04 + OH = #0.048 #0.361 SVOSOMG_07_05 + #0.048 #0.007 &
          SVOSOMG_07_06 + #0.048 #0.632 SVOSOMG_07_07 + #0.952 #0.091 SVOSOMG_01_01 + &
          #0.952 #0.091 SVOSOMG_02_01 + #0.952 #0.091 SVOSOMG_03_01 + #0.952 &
          #0.091 SVOSOMG_04_01 + #0.952 #0.000 SVOSOMG_05_01 + #0.952 #0.000 &
          SVOSOMG_06_01 + #0.952 #0.091 SVOSOMG_01_02 + #0.952 #0.091 SVOSOMG_02_02 + &
          #0.952 #0.091 SVOSOMG_03_02 + #0.952 #0.091 SVOSOMG_04_02 + #0.952 &
          #0.091 SVOSOMG_05_02 + #0.952 #0.000 SVOSOMG_06_02 + #0.952 #0.091 &
          SVOSOMG_02_03 + #0.952 #0.091 SVOSOMG_03_03 + #0.952 #0.091 SVOSOMG_04_03 + &
          #0.952 #0.091 SVOSOMG_05_03 + #0.952 #0.000 SVOSOMG_06_03 + #0.952 &
          #0.091 SVOSOMG_02_04 + #0.952 #0.091 SVOSOMG_03_04 + #0.952 #0.091 &
          SVOSOMG_04_04 + #0.952 #0.091 SVOSOMG_05_04 + #0.952 #0.091 SVOSOMG_06_04 + &
          #0.952 #0.091 SVOSOMG_03_05 + #0.952 #0.091 SVOSOMG_04_05 + #0.952 &
          #0.091 SVOSOMG_05_05 + #0.952 #0.091 SVOSOMG_06_05
 S381.7 SVOSOMG_07_05 + OH = #0.029 #0.361 SVOSOMG_07_06 + #0.029 #0.639 &
          SVOSOMG_07_07 + #0.971 #0.083 SVOSOMG_01_01 + #0.971 #0.083 SVOSOMG_02_01 + &
          #0.971 #0.083 SVOSOMG_03_01 + #0.971 #0.083 SVOSOMG_04_01 + #0.971 &
          #0.000 SVOSOMG_05_01 + #0.971 #0.000 SVOSOMG_06_01 + #0.971 #0.083 &
          SVOSOMG_01_02 + #0.971 #0.083 SVOSOMG_02_02 + #0.971 #0.083 SVOSOMG_03_02 + &
          #0.971 #0.083 SVOSOMG_04_02 + #0.971 #0.000 SVOSOMG_05_02 + #0.971 &
          #0.000 SVOSOMG_06_02 + #0.971 #0.083 SVOSOMG_02_03 + #0.971 #0.083 &
          SVOSOMG_03_03 + #0.971 #0.083 SVOSOMG_04_03 + #0.971 #0.083 SVOSOMG_05_03 + &
          #0.971 #0.000 SVOSOMG_06_03 + #0.971 #0.083 SVOSOMG_02_04 + #0.971 &
          #0.083 SVOSOMG_03_04 + #0.971 #0.083 SVOSOMG_04_04 + #0.971 #0.083 &
          SVOSOMG_05_04 + #0.971 #0.000 SVOSOMG_06_04 + #0.971 #0.083 SVOSOMG_03_05 + &
          #0.971 #0.083 SVOSOMG_04_05 + #0.971 #0.083 SVOSOMG_05_05 + #0.971 &
          #0.083 SVOSOMG_06_05 + #0.971 #0.083 SVOSOMG_03_06 + #0.971 #0.083 &
          SVOSOMG_04_06 + #0.971 #0.083 SVOSOMG_05_06 + #0.971 #0.083 SVOSOMG_06_06
 S382.7 SVOSOMG_07_06 + OH = #0.013 #1.000 SVOSOMG_07_07 + #0.987 #0.080 &
          SVOSOMG_01_01 + #0.987 #0.080 SVOSOMG_02_01 + #0.987 #0.080 SVOSOMG_03_01 + &
          #0.987 #0.000 SVOSOMG_04_01 + #0.987 #0.000 SVOSOMG_05_01 + #0.987 &
          #0.000 SVOSOMG_06_01 + #0.987 #0.080 SVOSOMG_01_02 + #0.987 #0.080 &
          SVOSOMG_02_02 + #0.987 #0.080 SVOSOMG_03_02 + #0.987 #0.080 SVOSOMG_04_02 + &
          #0.987 #0.000 SVOSOMG_05_02 + #0.987 #0.000 SVOSOMG_06_02 + #0.987 &
          #0.080 SVOSOMG_02_03 + #0.987 #0.080 SVOSOMG_03_03 + #0.987 #0.080 &
          SVOSOMG_04_03 + #0.987 #0.000 SVOSOMG_05_03 + #0.987 #0.000 SVOSOMG_06_03 + &
          #0.987 #0.080 SVOSOMG_02_04 + #0.987 #0.080 SVOSOMG_03_04 + #0.987 &
          #0.080 SVOSOMG_04_04 + #0.987 #0.080 SVOSOMG_05_04 + #0.987 #0.000 &
          SVOSOMG_06_04 + #0.987 #0.080 SVOSOMG_03_05 + #0.987 #0.080 SVOSOMG_04_05 + &
          #0.987 #0.080 SVOSOMG_05_05 + #0.987 #0.000 SVOSOMG_06_05 + #0.987 &
          #0.080 SVOSOMG_03_06 + #0.987 #0.080 SVOSOMG_04_06 + #0.987 #0.080 &
          SVOSOMG_05_06 + #0.987 #0.080 SVOSOMG_06_06 + #0.987 #0.080 SVOSOMG_04_07 + &
          #0.987 #0.080 SVOSOMG_05_07 + #0.987 #0.080 SVOSOMG_06_07
 S383.7 SVOSOMG_07_07 + OH = #1.000 #0.091 SVOSOMG_01_01 + #1.000 #0.091 &
          SVOSOMG_02_01 + #1.000 #0.091 SVOSOMG_03_01 + #1.000 #0.000 SVOSOMG_04_01 + &
          #1.000 #0.000 SVOSOMG_05_01 + #1.000 #0.000 SVOSOMG_06_01 + #1.000 &
          #0.091 SVOSOMG_01_02 + #1.000 #0.091 SVOSOMG_02_02 + #1.000 #0.091 &
          SVOSOMG_03_02 + #1.000 #0.000 SVOSOMG_04_02 + #1.000 #0.000 SVOSOMG_05_02 + &
          #1.000 #0.000 SVOSOMG_06_02 + #1.000 #0.091 SVOSOMG_02_03 + #1.000 &
          #0.091 SVOSOMG_03_03 + #1.000 #0.091 SVOSOMG_04_03 + #1.000 #0.000 &
          SVOSOMG_05_03 + #1.000 #0.000 SVOSOMG_06_03 + #1.000 #0.091 SVOSOMG_02_04 + &
          #1.000 #0.091 SVOSOMG_03_04 + #1.000 #0.091 SVOSOMG_04_04 + #1.000 &
          #0.000 SVOSOMG_05_04 + #1.000 #0.000 SVOSOMG_06_04 + #1.000 #0.091 &
          SVOSOMG_03_05 + #1.000 #0.091 SVOSOMG_04_05 + #1.000 #0.091 SVOSOMG_05_05 + &
          #1.000 #0.000 SVOSOMG_06_05 + #1.000 #0.091 SVOSOMG_03_06 + #1.000 &
          #0.091 SVOSOMG_04_06 + #1.000 #0.091 SVOSOMG_05_06 + #1.000 #0.000 &
          SVOSOMG_06_06 + #1.000 #0.091 SVOSOMG_04_07 + #1.000 #0.091 SVOSOMG_05_07 + &
          #1.000 #0.091 SVOSOMG_06_07
 S384.7 SVOSOMG_08_01 + OH = #0.165 #0.361 SVOSOMG_08_02 + #0.165 #0.007 &
          SVOSOMG_08_03 + #0.165 #0.606 SVOSOMG_08_04 + #0.165 #0.026 SVOSOMG_08_05 + &
          #0.835 #0.143 SVOSOMG_01_01 + #0.835 #0.143 SVOSOMG_02_01 + #0.835 &
          #0.143 SVOSOMG_03_01 + #0.835 #0.143 SVOSOMG_04_01 + #0.835 #0.143 &
          SVOSOMG_05_01 + #0.835 #0.143 SVOSOMG_06_01 + #0.835 #0.143 SVOSOMG_07_01 + &
          #0.835 #0.143 SVOSOMG_01_02 + #0.835 #0.143 SVOSOMG_02_02 + #0.835 &
          #0.143 SVOSOMG_03_02 + #0.835 #0.143 SVOSOMG_04_02 + #0.835 #0.143 &
          SVOSOMG_05_02 + #0.835 #0.143 SVOSOMG_06_02 + #0.835 #0.143 SVOSOMG_07_02
 S385.7 SVOSOMG_08_02 + OH = #0.114 #0.361 SVOSOMG_08_03 + #0.114 #0.007 &
          SVOSOMG_08_04 + #0.114 #0.606 SVOSOMG_08_05 + #0.114 #0.026 SVOSOMG_08_06 + &
          #0.886 #0.105 SVOSOMG_01_01 + #0.886 #0.105 SVOSOMG_02_01 + #0.886 &
          #0.105 SVOSOMG_03_01 + #0.886 #0.105 SVOSOMG_04_01 + #0.886 #0.105 &
          SVOSOMG_05_01 + #0.886 #0.105 SVOSOMG_06_01 + #0.886 #0.000 SVOSOMG_07_01 + &
          #0.886 #0.105 SVOSOMG_01_02 + #0.886 #0.105 SVOSOMG_02_02 + #0.886 &
          #0.105 SVOSOMG_03_02 + #0.886 #0.105 SVOSOMG_04_02 + #0.886 #0.105 &
          SVOSOMG_05_02 + #0.886 #0.105 SVOSOMG_06_02 + #0.886 #0.105 SVOSOMG_07_02 + &
          #0.886 #0.105 SVOSOMG_02_03 + #0.886 #0.105 SVOSOMG_03_03 + #0.886 &
          #0.105 SVOSOMG_04_03 + #0.886 #0.105 SVOSOMG_05_03 + #0.886 #0.105 &
          SVOSOMG_06_03 + #0.886 #0.105 SVOSOMG_07_03
 S386.7 SVOSOMG_08_03 + OH = #0.082 #0.361 SVOSOMG_08_04 + #0.082 #0.007 &
          SVOSOMG_08_05 + #0.082 #0.606 SVOSOMG_08_06 + #0.082 #0.026 SVOSOMG_08_07 + &
          #0.918 #0.083 SVOSOMG_01_01 + #0.918 #0.083 SVOSOMG_02_01 + #0.918 &
          #0.083 SVOSOMG_03_01 + #0.918 #0.083 SVOSOMG_04_01 + #0.918 #0.083 &
          SVOSOMG_05_01 + #0.918 #0.083 SVOSOMG_06_01 + #0.918 #0.000 SVOSOMG_07_01 + &
          #0.918 #0.083 SVOSOMG_01_02 + #0.918 #0.083 SVOSOMG_02_02 + #0.918 &
          #0.083 SVOSOMG_03_02 + #0.918 #0.083 SVOSOMG_04_02 + #0.918 #0.083 &
          SVOSOMG_05_02 + #0.918 #0.083 SVOSOMG_06_02 + #0.918 #0.000 SVOSOMG_07_02 + &
          #0.918 #0.083 SVOSOMG_02_03 + #0.918 #0.083 SVOSOMG_03_03 + #0.918 &
          #0.083 SVOSOMG_04_03 + #0.918 #0.083 SVOSOMG_05_03 + #0.918 #0.083 &
          SVOSOMG_06_03 + #0.918 #0.083 SVOSOMG_07_03 + #0.918 #0.083 SVOSOMG_02_04 + &
          #0.918 #0.083 SVOSOMG_03_04 + #0.918 #0.083 SVOSOMG_04_04 + #0.918 &
          #0.083 SVOSOMG_05_04 + #0.918 #0.083 SVOSOMG_06_04 + #0.918 #0.083 &
          SVOSOMG_07_04
 S387.7 SVOSOMG_08_04 + OH = #0.059 #0.361 SVOSOMG_08_05 + #0.059 #0.007 &
          SVOSOMG_08_06 + #0.059 #0.632 SVOSOMG_08_07 + #0.941 #0.074 SVOSOMG_01_01 + &
          #0.941 #0.074 SVOSOMG_02_01 + #0.941 #0.074 SVOSOMG_03_01 + #0.941 &
          #0.074 SVOSOMG_04_01 + #0.941 #0.074 SVOSOMG_05_01 + #0.941 #0.000 &
          SVOSOMG_06_01 + #0.941 #0.000 SVOSOMG_07_01 + #0.941 #0.074 SVOSOMG_01_02 + &
          #0.941 #0.074 SVOSOMG_02_02 + #0.941 #0.074 SVOSOMG_03_02 + #0.941 &
          #0.074 SVOSOMG_04_02 + #0.941 #0.074 SVOSOMG_05_02 + #0.941 #0.074 &
          SVOSOMG_06_02 + #0.941 #0.000 SVOSOMG_07_02 + #0.941 #0.074 SVOSOMG_02_03 + &
          #0.941 #0.074 SVOSOMG_03_03 + #0.941 #0.074 SVOSOMG_04_03 + #0.941 &
          #0.074 SVOSOMG_05_03 + #0.941 #0.074 SVOSOMG_06_03 + #0.941 #0.000 &
          SVOSOMG_07_03 + #0.941 #0.074 SVOSOMG_02_04 + #0.941 #0.074 SVOSOMG_03_04 + &
          #0.941 #0.074 SVOSOMG_04_04 + #0.941 #0.074 SVOSOMG_05_04 + #0.941 &
          #0.074 SVOSOMG_06_04 + #0.941 #0.074 SVOSOMG_07_04 + #0.941 #0.074 &
          SVOSOMG_03_05 + #0.941 #0.074 SVOSOMG_04_05 + #0.941 #0.074 SVOSOMG_05_05 + &
          #0.941 #0.074 SVOSOMG_06_05 + #0.941 #0.074 SVOSOMG_07_05
 S388.7 SVOSOMG_08_05 + OH = #0.040 #0.361 SVOSOMG_08_06 + #0.040 #0.639 &
          SVOSOMG_08_07 + #0.960 #0.067 SVOSOMG_01_01 + #0.960 #0.067 SVOSOMG_02_01 + &
          #0.960 #0.067 SVOSOMG_03_01 + #0.960 #0.067 SVOSOMG_04_01 + #0.960 &
          #0.067 SVOSOMG_05_01 + #0.960 #0.000 SVOSOMG_06_01 + #0.960 #0.000 &
          SVOSOMG_07_01 + #0.960 #0.067 SVOSOMG_01_02 + #0.960 #0.067 SVOSOMG_02_02 + &
          #0.960 #0.067 SVOSOMG_03_02 + #0.960 #0.067 SVOSOMG_04_02 + #0.960 &
          #0.067 SVOSOMG_05_02 + #0.960 #0.000 SVOSOMG_06_02 + #0.960 #0.000 &
          SVOSOMG_07_02 + #0.960 #0.067 SVOSOMG_02_03 + #0.960 #0.067 SVOSOMG_03_03 + &
          #0.960 #0.067 SVOSOMG_04_03 + #0.960 #0.067 SVOSOMG_05_03 + #0.960 &
          #0.067 SVOSOMG_06_03 + #0.960 #0.000 SVOSOMG_07_03 + #0.960 #0.067 &
          SVOSOMG_02_04 + #0.960 #0.067 SVOSOMG_03_04 + #0.960 #0.067 SVOSOMG_04_04 + &
          #0.960 #0.067 SVOSOMG_05_04 + #0.960 #0.067 SVOSOMG_06_04 + #0.960 &
          #0.000 SVOSOMG_07_04 + #0.960 #0.067 SVOSOMG_03_05 + #0.960 #0.067 &
          SVOSOMG_04_05 + #0.960 #0.067 SVOSOMG_05_05 + #0.960 #0.067 SVOSOMG_06_05 + &
          #0.960 #0.067 SVOSOMG_07_05 + #0.960 #0.067 SVOSOMG_03_06 + #0.960 &
          #0.067 SVOSOMG_04_06 + #0.960 #0.067 SVOSOMG_05_06 + #0.960 #0.067 &
          SVOSOMG_06_06 + #0.960 #0.067 SVOSOMG_07_06
 S389.7 SVOSOMG_08_06 + OH = #0.025 #1.000 SVOSOMG_08_07 + #0.975 #0.065 &
          SVOSOMG_01_01 + #0.975 #0.065 SVOSOMG_02_01 + #0.975 #0.065 SVOSOMG_03_01 + &
          #0.975 #0.065 SVOSOMG_04_01 + #0.975 #0.000 SVOSOMG_05_01 + #0.975 &
          #0.000 SVOSOMG_06_01 + #0.975 #0.000 SVOSOMG_07_01 + #0.975 #0.065 &
          SVOSOMG_01_02 + #0.975 #0.065 SVOSOMG_02_02 + #0.975 #0.065 SVOSOMG_03_02 + &
          #0.975 #0.065 SVOSOMG_04_02 + #0.975 #0.065 SVOSOMG_05_02 + #0.975 &
          #0.000 SVOSOMG_06_02 + #0.975 #0.000 SVOSOMG_07_02 + #0.975 #0.065 &
          SVOSOMG_02_03 + #0.975 #0.065 SVOSOMG_03_03 + #0.975 #0.065 SVOSOMG_04_03 + &
          #0.975 #0.065 SVOSOMG_05_03 + #0.975 #0.000 SVOSOMG_06_03 + #0.975 &
          #0.000 SVOSOMG_07_03 + #0.975 #0.065 SVOSOMG_02_04 + #0.975 #0.065 &
          SVOSOMG_03_04 + #0.975 #0.065 SVOSOMG_04_04 + #0.975 #0.065 SVOSOMG_05_04 + &
          #0.975 #0.065 SVOSOMG_06_04 + #0.975 #0.000 SVOSOMG_07_04 + #0.975 &
          #0.065 SVOSOMG_03_05 + #0.975 #0.065 SVOSOMG_04_05 + #0.975 #0.065 &
          SVOSOMG_05_05 + #0.975 #0.065 SVOSOMG_06_05 + #0.975 #0.000 SVOSOMG_07_05 + &
          #0.975 #0.065 SVOSOMG_03_06 + #0.975 #0.065 SVOSOMG_04_06 + #0.975 &
          #0.065 SVOSOMG_05_06 + #0.975 #0.065 SVOSOMG_06_06 + #0.975 #0.065 &
          SVOSOMG_07_06 + #0.975 #0.065 SVOSOMG_04_07 + #0.975 #0.065 SVOSOMG_05_07 + &
          #0.975 #0.065 SVOSOMG_06_07 + #0.975 #0.065 SVOSOMG_07_07
 S390.7 SVOSOMG_08_07 + OH = #1.000 #0.071 SVOSOMG_01_01 + #1.000 #0.071 &
          SVOSOMG_02_01 + #1.000 #0.071 SVOSOMG_03_01 + #1.000 #0.071 SVOSOMG_04_01 + &
          #1.000 #0.000 SVOSOMG_05_01 + #1.000 #0.000 SVOSOMG_06_01 + #1.000 &
          #0.000 SVOSOMG_07_01 + #1.000 #0.071 SVOSOMG_01_02 + #1.000 #0.071 &
          SVOSOMG_02_02 + #1.000 #0.071 SVOSOMG_03_02 + #1.000 #0.071 SVOSOMG_04_02 + &
          #1.000 #0.000 SVOSOMG_05_02 + #1.000 #0.000 SVOSOMG_06_02 + #1.000 &
          #0.000 SVOSOMG_07_02 + #1.000 #0.071 SVOSOMG_02_03 + #1.000 #0.071 &
          SVOSOMG_03_03 + #1.000 #0.071 SVOSOMG_04_03 + #1.000 #0.071 SVOSOMG_05_03 + &
          #1.000 #0.000 SVOSOMG_06_03 + #1.000 #0.000 SVOSOMG_07_03 + #1.000 &
          #0.071 SVOSOMG_02_04 + #1.000 #0.071 SVOSOMG_03_04 + #1.000 #0.071 &
          SVOSOMG_04_04 + #1.000 #0.071 SVOSOMG_05_04 + #1.000 #0.000 SVOSOMG_06_04 + &
          #1.000 #0.000 SVOSOMG_07_04 + #1.000 #0.071 SVOSOMG_03_05 + #1.000 &
          #0.071 SVOSOMG_04_05 + #1.000 #0.071 SVOSOMG_05_05 + #1.000 #0.071 &
          SVOSOMG_06_05 + #1.000 #0.000 SVOSOMG_07_05 + #1.000 #0.071 SVOSOMG_03_06 + &
          #1.000 #0.071 SVOSOMG_04_06 + #1.000 #0.071 SVOSOMG_05_06 + #1.000 &
          #0.071 SVOSOMG_06_06 + #1.000 #0.000 SVOSOMG_07_06 + #1.000 #0.071 &
          SVOSOMG_04_07 + #1.000 #0.071 SVOSOMG_05_07 + #1.000 #0.071 SVOSOMG_06_07 + &
          #1.000 #0.071 SVOSOMG_07_07
 S391.7 SVOSOMG_09_01 + OH = #0.174 #0.361 SVOSOMG_09_02 + #0.174 #0.007 &
          SVOSOMG_09_03 + #0.174 #0.606 SVOSOMG_09_04 + #0.174 #0.026 SVOSOMG_09_05 + &
          #0.826 #0.125 SVOSOMG_01_01 + #0.826 #0.125 SVOSOMG_02_01 + #0.826 &
          #0.125 SVOSOMG_03_01 + #0.826 #0.125 SVOSOMG_04_01 + #0.826 #0.125 &
          SVOSOMG_05_01 + #0.826 #0.125 SVOSOMG_06_01 + #0.826 #0.125 SVOSOMG_07_01 + &
          #0.826 #0.125 SVOSOMG_08_01 + #0.826 #0.125 SVOSOMG_01_02 + #0.826 &
          #0.125 SVOSOMG_02_02 + #0.826 #0.125 SVOSOMG_03_02 + #0.826 #0.125 &
          SVOSOMG_04_02 + #0.826 #0.125 SVOSOMG_05_02 + #0.826 #0.125 SVOSOMG_06_02 + &
          #0.826 #0.125 SVOSOMG_07_02 + #0.826 #0.125 SVOSOMG_08_02
 S392.7 SVOSOMG_09_02 + OH = #0.123 #0.361 SVOSOMG_09_03 + #0.123 #0.007 &
          SVOSOMG_09_04 + #0.123 #0.606 SVOSOMG_09_05 + #0.123 #0.026 SVOSOMG_09_06 + &
          #0.877 #0.087 SVOSOMG_01_01 + #0.877 #0.087 SVOSOMG_02_01 + #0.877 &
          #0.087 SVOSOMG_03_01 + #0.877 #0.087 SVOSOMG_04_01 + #0.877 #0.087 &
          SVOSOMG_05_01 + #0.877 #0.087 SVOSOMG_06_01 + #0.877 #0.087 SVOSOMG_07_01 + &
          #0.877 #0.000 SVOSOMG_08_01 + #0.877 #0.087 SVOSOMG_01_02 + #0.877 &
          #0.087 SVOSOMG_02_02 + #0.877 #0.087 SVOSOMG_03_02 + #0.877 #0.087 &
          SVOSOMG_04_02 + #0.877 #0.087 SVOSOMG_05_02 + #0.877 #0.087 SVOSOMG_06_02 + &
          #0.877 #0.087 SVOSOMG_07_02 + #0.877 #0.087 SVOSOMG_08_02 + #0.877 &
          #0.087 SVOSOMG_02_03 + #0.877 #0.087 SVOSOMG_03_03 + #0.877 #0.087 &
          SVOSOMG_04_03 + #0.877 #0.087 SVOSOMG_05_03 + #0.877 #0.087 SVOSOMG_06_03 + &
          #0.877 #0.087 SVOSOMG_07_03 + #0.877 #0.087 SVOSOMG_08_03
 S393.7 SVOSOMG_09_03 + OH = #0.091 #0.361 SVOSOMG_09_04 + #0.091 #0.007 &
          SVOSOMG_09_05 + #0.091 #0.606 SVOSOMG_09_06 + #0.091 #0.026 SVOSOMG_09_07 + &
          #0.909 #0.069 SVOSOMG_01_01 + #0.909 #0.069 SVOSOMG_02_01 + #0.909 &
          #0.069 SVOSOMG_03_01 + #0.909 #0.069 SVOSOMG_04_01 + #0.909 #0.069 &
          SVOSOMG_05_01 + #0.909 #0.069 SVOSOMG_06_01 + #0.909 #0.069 SVOSOMG_07_01 + &
          #0.909 #0.000 SVOSOMG_08_01 + #0.909 #0.069 SVOSOMG_01_02 + #0.909 &
          #0.069 SVOSOMG_02_02 + #0.909 #0.069 SVOSOMG_03_02 + #0.909 #0.069 &
          SVOSOMG_04_02 + #0.909 #0.069 SVOSOMG_05_02 + #0.909 #0.069 SVOSOMG_06_02 + &
          #0.909 #0.069 SVOSOMG_07_02 + #0.909 #0.000 SVOSOMG_08_02 + #0.909 &
          #0.069 SVOSOMG_02_03 + #0.909 #0.069 SVOSOMG_03_03 + #0.909 #0.069 &
          SVOSOMG_04_03 + #0.909 #0.069 SVOSOMG_05_03 + #0.909 #0.069 SVOSOMG_06_03 + &
          #0.909 #0.069 SVOSOMG_07_03 + #0.909 #0.069 SVOSOMG_08_03 + #0.909 &
          #0.069 SVOSOMG_02_04 + #0.909 #0.069 SVOSOMG_03_04 + #0.909 #0.069 &
          SVOSOMG_04_04 + #0.909 #0.069 SVOSOMG_05_04 + #0.909 #0.069 SVOSOMG_06_04 + &
          #0.909 #0.069 SVOSOMG_07_04 + #0.909 #0.069 SVOSOMG_08_04
 S394.7 SVOSOMG_09_04 + OH = #0.068 #0.361 SVOSOMG_09_05 + #0.068 #0.007 &
          SVOSOMG_09_06 + #0.068 #0.632 SVOSOMG_09_07 + #0.932 #0.059 SVOSOMG_01_01 + &
          #0.932 #0.059 SVOSOMG_02_01 + #0.932 #0.059 SVOSOMG_03_01 + #0.932 &
          #0.059 SVOSOMG_04_01 + #0.932 #0.059 SVOSOMG_05_01 + #0.932 #0.059 &
          SVOSOMG_06_01 + #0.932 #0.000 SVOSOMG_07_01 + #0.932 #0.000 SVOSOMG_08_01 + &
          #0.932 #0.059 SVOSOMG_01_02 + #0.932 #0.059 SVOSOMG_02_02 + #0.932 &
          #0.059 SVOSOMG_03_02 + #0.932 #0.059 SVOSOMG_04_02 + #0.932 #0.059 &
          SVOSOMG_05_02 + #0.932 #0.059 SVOSOMG_06_02 + #0.932 #0.059 SVOSOMG_07_02 + &
          #0.932 #0.000 SVOSOMG_08_02 + #0.932 #0.059 SVOSOMG_02_03 + #0.932 &
          #0.059 SVOSOMG_03_03 + #0.932 #0.059 SVOSOMG_04_03 + #0.932 #0.059 &
          SVOSOMG_05_03 + #0.932 #0.059 SVOSOMG_06_03 + #0.932 #0.059 SVOSOMG_07_03 + &
          #0.932 #0.000 SVOSOMG_08_03 + #0.932 #0.059 SVOSOMG_02_04 + #0.932 &
          #0.059 SVOSOMG_03_04 + #0.932 #0.059 SVOSOMG_04_04 + #0.932 #0.059 &
          SVOSOMG_05_04 + #0.932 #0.059 SVOSOMG_06_04 + #0.932 #0.059 SVOSOMG_07_04 + &
          #0.932 #0.059 SVOSOMG_08_04 + #0.932 #0.059 SVOSOMG_03_05 + #0.932 &
          #0.059 SVOSOMG_04_05 + #0.932 #0.059 SVOSOMG_05_05 + #0.932 #0.059 &
          SVOSOMG_06_05 + #0.932 #0.059 SVOSOMG_07_05 + #0.932 #0.059 SVOSOMG_08_05
 S395.7 SVOSOMG_09_05 + OH = #0.050 #0.361 SVOSOMG_09_06 + #0.050 #0.639 &
          SVOSOMG_09_07 + #0.950 #0.053 SVOSOMG_01_01 + #0.950 #0.053 SVOSOMG_02_01 + &
          #0.950 #0.053 SVOSOMG_03_01 + #0.950 #0.053 SVOSOMG_04_01 + #0.950 &
          #0.053 SVOSOMG_05_01 + #0.950 #0.053 SVOSOMG_06_01 + #0.950 #0.000 &
          SVOSOMG_07_01 + #0.950 #0.000 SVOSOMG_08_01 + #0.950 #0.053 SVOSOMG_01_02 + &
          #0.950 #0.053 SVOSOMG_02_02 + #0.950 #0.053 SVOSOMG_03_02 + #0.950 &
          #0.053 SVOSOMG_04_02 + #0.950 #0.053 SVOSOMG_05_02 + #0.950 #0.053 &
          SVOSOMG_06_02 + #0.950 #0.000 SVOSOMG_07_02 + #0.950 #0.000 SVOSOMG_08_02 + &
          #0.950 #0.053 SVOSOMG_02_03 + #0.950 #0.053 SVOSOMG_03_03 + #0.950 &
          #0.053 SVOSOMG_04_03 + #0.950 #0.053 SVOSOMG_05_03 + #0.950 #0.053 &
          SVOSOMG_06_03 + #0.950 #0.053 SVOSOMG_07_03 + #0.950 #0.000 SVOSOMG_08_03 + &
          #0.950 #0.053 SVOSOMG_02_04 + #0.950 #0.053 SVOSOMG_03_04 + #0.950 &
          #0.053 SVOSOMG_04_04 + #0.950 #0.053 SVOSOMG_05_04 + #0.950 #0.053 &
          SVOSOMG_06_04 + #0.950 #0.053 SVOSOMG_07_04 + #0.950 #0.000 SVOSOMG_08_04 + &
          #0.950 #0.053 SVOSOMG_03_05 + #0.950 #0.053 SVOSOMG_04_05 + #0.950 &
          #0.053 SVOSOMG_05_05 + #0.950 #0.053 SVOSOMG_06_05 + #0.950 #0.053 &
          SVOSOMG_07_05 + #0.950 #0.053 SVOSOMG_08_05 + #0.950 #0.053 SVOSOMG_03_06 + &
          #0.950 #0.053 SVOSOMG_04_06 + #0.950 #0.053 SVOSOMG_05_06 + #0.950 &
          #0.053 SVOSOMG_06_06 + #0.950 #0.053 SVOSOMG_07_06 + #0.950 #0.053 &
          SVOSOMG_08_06
 S396.7 SVOSOMG_09_06 + OH = #0.035 #1.000 SVOSOMG_09_07 + #0.965 #0.049 &
          SVOSOMG_01_01 + #0.965 #0.049 SVOSOMG_02_01 + #0.965 #0.049 SVOSOMG_03_01 + &
          #0.965 #0.049 SVOSOMG_04_01 + #0.965 #0.049 SVOSOMG_05_01 + #0.965 &
          #0.000 SVOSOMG_06_01 + #0.965 #0.000 SVOSOMG_07_01 + #0.965 #0.000 &
          SVOSOMG_08_01 + #0.965 #0.049 SVOSOMG_01_02 + #0.965 #0.049 SVOSOMG_02_02 + &
          #0.965 #0.049 SVOSOMG_03_02 + #0.965 #0.049 SVOSOMG_04_02 + #0.965 &
          #0.049 SVOSOMG_05_02 + #0.965 #0.049 SVOSOMG_06_02 + #0.965 #0.000 &
          SVOSOMG_07_02 + #0.965 #0.000 SVOSOMG_08_02 + #0.965 #0.049 SVOSOMG_02_03 + &
          #0.965 #0.049 SVOSOMG_03_03 + #0.965 #0.049 SVOSOMG_04_03 + #0.965 &
          #0.049 SVOSOMG_05_03 + #0.965 #0.049 SVOSOMG_06_03 + #0.965 #0.000 &
          SVOSOMG_07_03 + #0.965 #0.000 SVOSOMG_08_03 + #0.965 #0.049 SVOSOMG_02_04 + &
          #0.965 #0.049 SVOSOMG_03_04 + #0.965 #0.049 SVOSOMG_04_04 + #0.965 &
          #0.049 SVOSOMG_05_04 + #0.965 #0.049 SVOSOMG_06_04 + #0.965 #0.049 &
          SVOSOMG_07_04 + #0.965 #0.000 SVOSOMG_08_04 + #0.965 #0.049 SVOSOMG_03_05 + &
          #0.965 #0.049 SVOSOMG_04_05 + #0.965 #0.049 SVOSOMG_05_05 + #0.965 &
          #0.049 SVOSOMG_06_05 + #0.965 #0.049 SVOSOMG_07_05 + #0.965 #0.000 &
          SVOSOMG_08_05 + #0.965 #0.049 SVOSOMG_03_06 + #0.965 #0.049 SVOSOMG_04_06 + &
          #0.965 #0.049 SVOSOMG_05_06 + #0.965 #0.049 SVOSOMG_06_06 + #0.965 &
          #0.049 SVOSOMG_07_06 + #0.965 #0.049 SVOSOMG_08_06 + #0.965 #0.049 &
          SVOSOMG_04_07 + #0.965 #0.049 SVOSOMG_05_07 + #0.965 #0.049 SVOSOMG_06_07 + &
          #0.965 #0.049 SVOSOMG_07_07 + #0.965 #0.049 SVOSOMG_08_07
 S397.7 SVOSOMG_09_07 + OH = #1.000 #0.053 SVOSOMG_01_01 + #1.000 #0.053 &
          SVOSOMG_02_01 + #1.000 #0.053 SVOSOMG_03_01 + #1.000 #0.053 SVOSOMG_04_01 + &
          #1.000 #0.053 SVOSOMG_05_01 + #1.000 #0.000 SVOSOMG_06_01 + #1.000 &
          #0.000 SVOSOMG_07_01 + #1.000 #0.000 SVOSOMG_08_01 + #1.000 #0.053 &
          SVOSOMG_01_02 + #1.000 #0.053 SVOSOMG_02_02 + #1.000 #0.053 SVOSOMG_03_02 + &
          #1.000 #0.053 SVOSOMG_04_02 + #1.000 #0.053 SVOSOMG_05_02 + #1.000 &
          #0.000 SVOSOMG_06_02 + #1.000 #0.000 SVOSOMG_07_02 + #1.000 #0.000 &
          SVOSOMG_08_02 + #1.000 #0.053 SVOSOMG_02_03 + #1.000 #0.053 SVOSOMG_03_03 + &
          #1.000 #0.053 SVOSOMG_04_03 + #1.000 #0.053 SVOSOMG_05_03 + #1.000 &
          #0.053 SVOSOMG_06_03 + #1.000 #0.000 SVOSOMG_07_03 + #1.000 #0.000 &
          SVOSOMG_08_03 + #1.000 #0.053 SVOSOMG_02_04 + #1.000 #0.053 SVOSOMG_03_04 + &
          #1.000 #0.053 SVOSOMG_04_04 + #1.000 #0.053 SVOSOMG_05_04 + #1.000 &
          #0.053 SVOSOMG_06_04 + #1.000 #0.000 SVOSOMG_07_04 + #1.000 #0.000 &
          SVOSOMG_08_04 + #1.000 #0.053 SVOSOMG_03_05 + #1.000 #0.053 SVOSOMG_04_05 + &
          #1.000 #0.053 SVOSOMG_05_05 + #1.000 #0.053 SVOSOMG_06_05 + #1.000 &
          #0.053 SVOSOMG_07_05 + #1.000 #0.000 SVOSOMG_08_05 + #1.000 #0.053 &
          SVOSOMG_03_06 + #1.000 #0.053 SVOSOMG_04_06 + #1.000 #0.053 SVOSOMG_05_06 + &
          #1.000 #0.053 SVOSOMG_06_06 + #1.000 #0.053 SVOSOMG_07_06 + #1.000 &
          #0.000 SVOSOMG_08_06 + #1.000 #0.053 SVOSOMG_04_07 + #1.000 #0.053 &
          SVOSOMG_05_07 + #1.000 #0.053 SVOSOMG_06_07 + #1.000 #0.053 SVOSOMG_07_07 + &
          #1.000 #0.053 SVOSOMG_08_07
 S398.7 SVOSOMG_10_01 + OH = #0.182 #0.361 SVOSOMG_10_02 + #0.182 #0.007 &
          SVOSOMG_10_03 + #0.182 #0.606 SVOSOMG_10_04 + #0.182 #0.026 SVOSOMG_10_05 + &
          #0.818 #0.111 SVOSOMG_01_01 + #0.818 #0.111 SVOSOMG_02_01 + #0.818 &
          #0.111 SVOSOMG_03_01 + #0.818 #0.111 SVOSOMG_04_01 + #0.818 #0.111 &
          SVOSOMG_05_01 + #0.818 #0.111 SVOSOMG_06_01 + #0.818 #0.111 SVOSOMG_07_01 + &
          #0.818 #0.111 SVOSOMG_08_01 + #0.818 #0.111 SVOSOMG_09_01 + #0.818 &
          #0.111 SVOSOMG_01_02 + #0.818 #0.111 SVOSOMG_02_02 + #0.818 #0.111 &
          SVOSOMG_03_02 + #0.818 #0.111 SVOSOMG_04_02 + #0.818 #0.111 SVOSOMG_05_02 + &
          #0.818 #0.111 SVOSOMG_06_02 + #0.818 #0.111 SVOSOMG_07_02 + #0.818 &
          #0.111 SVOSOMG_08_02 + #0.818 #0.111 SVOSOMG_09_02
 S399.7 SVOSOMG_10_02 + OH = #0.131 #0.361 SVOSOMG_10_03 + #0.131 #0.007 &
          SVOSOMG_10_04 + #0.131 #0.606 SVOSOMG_10_05 + #0.131 #0.026 SVOSOMG_10_06 + &
          #0.869 #0.080 SVOSOMG_01_01 + #0.869 #0.080 SVOSOMG_02_01 + #0.869 &
          #0.080 SVOSOMG_03_01 + #0.869 #0.080 SVOSOMG_04_01 + #0.869 #0.080 &
          SVOSOMG_05_01 + #0.869 #0.080 SVOSOMG_06_01 + #0.869 #0.080 SVOSOMG_07_01 + &
          #0.869 #0.080 SVOSOMG_08_01 + #0.869 #0.000 SVOSOMG_09_01 + #0.869 &
          #0.080 SVOSOMG_01_02 + #0.869 #0.080 SVOSOMG_02_02 + #0.869 #0.080 &
          SVOSOMG_03_02 + #0.869 #0.080 SVOSOMG_04_02 + #0.869 #0.080 SVOSOMG_05_02 + &
          #0.869 #0.080 SVOSOMG_06_02 + #0.869 #0.080 SVOSOMG_07_02 + #0.869 &
          #0.080 SVOSOMG_08_02 + #0.869 #0.080 SVOSOMG_09_02 + #0.869 #0.080 &
          SVOSOMG_02_03 + #0.869 #0.080 SVOSOMG_03_03 + #0.869 #0.080 SVOSOMG_04_03 + &
          #0.869 #0.080 SVOSOMG_05_03 + #0.869 #0.080 SVOSOMG_06_03 + #0.869 &
          #0.080 SVOSOMG_07_03 + #0.869 #0.080 SVOSOMG_08_03 + #0.869 #0.080 &
          SVOSOMG_09_03
 S400.7 SVOSOMG_10_03 + OH = #0.099 #0.361 SVOSOMG_10_04 + #0.099 #0.007 &
          SVOSOMG_10_05 + #0.099 #0.606 SVOSOMG_10_06 + #0.099 #0.026 SVOSOMG_10_07 + &
          #0.901 #0.062 SVOSOMG_01_01 + #0.901 #0.062 SVOSOMG_02_01 + #0.901 &
          #0.062 SVOSOMG_03_01 + #0.901 #0.062 SVOSOMG_04_01 + #0.901 #0.062 &
          SVOSOMG_05_01 + #0.901 #0.062 SVOSOMG_06_01 + #0.901 #0.062 SVOSOMG_07_01 + &
          #0.901 #0.062 SVOSOMG_08_01 + #0.901 #0.000 SVOSOMG_09_01 + #0.901 &
          #0.062 SVOSOMG_01_02 + #0.901 #0.062 SVOSOMG_02_02 + #0.901 #0.062 &
          SVOSOMG_03_02 + #0.901 #0.062 SVOSOMG_04_02 + #0.901 #0.062 SVOSOMG_05_02 + &
          #0.901 #0.062 SVOSOMG_06_02 + #0.901 #0.062 SVOSOMG_07_02 + #0.901 &
          #0.062 SVOSOMG_08_02 + #0.901 #0.000 SVOSOMG_09_02 + #0.901 #0.062 &
          SVOSOMG_02_03 + #0.901 #0.062 SVOSOMG_03_03 + #0.901 #0.062 SVOSOMG_04_03 + &
          #0.901 #0.062 SVOSOMG_05_03 + #0.901 #0.062 SVOSOMG_06_03 + #0.901 &
          #0.062 SVOSOMG_07_03 + #0.901 #0.062 SVOSOMG_08_03 + #0.901 #0.062 &
          SVOSOMG_09_03 + #0.901 #0.062 SVOSOMG_02_04 + #0.901 #0.062 SVOSOMG_03_04 + &
          #0.901 #0.062 SVOSOMG_04_04 + #0.901 #0.062 SVOSOMG_05_04 + #0.901 &
          #0.062 SVOSOMG_06_04 + #0.901 #0.062 SVOSOMG_07_04 + #0.901 #0.062 &
          SVOSOMG_08_04 + #0.901 #0.062 SVOSOMG_09_04
 S401.7 SVOSOMG_10_04 + OH = #0.077 #0.361 SVOSOMG_10_05 + #0.077 #0.007 &
          SVOSOMG_10_06 + #0.077 #0.632 SVOSOMG_10_07 + #0.923 #0.053 SVOSOMG_01_01 + &
          #0.923 #0.053 SVOSOMG_02_01 + #0.923 #0.053 SVOSOMG_03_01 + #0.923 &
          #0.053 SVOSOMG_04_01 + #0.923 #0.053 SVOSOMG_05_01 + #0.923 #0.053 &
          SVOSOMG_06_01 + #0.923 #0.053 SVOSOMG_07_01 + #0.923 #0.000 SVOSOMG_08_01 + &
          #0.923 #0.000 SVOSOMG_09_01 + #0.923 #0.053 SVOSOMG_01_02 + #0.923 &
          #0.053 SVOSOMG_02_02 + #0.923 #0.053 SVOSOMG_03_02 + #0.923 #0.053 &
          SVOSOMG_04_02 + #0.923 #0.053 SVOSOMG_05_02 + #0.923 #0.053 SVOSOMG_06_02 + &
          #0.923 #0.053 SVOSOMG_07_02 + #0.923 #0.053 SVOSOMG_08_02 + #0.923 &
          #0.000 SVOSOMG_09_02 + #0.923 #0.053 SVOSOMG_02_03 + #0.923 #0.053 &
          SVOSOMG_03_03 + #0.923 #0.053 SVOSOMG_04_03 + #0.923 #0.053 SVOSOMG_05_03 + &
          #0.923 #0.053 SVOSOMG_06_03 + #0.923 #0.053 SVOSOMG_07_03 + #0.923 &
          #0.053 SVOSOMG_08_03 + #0.923 #0.000 SVOSOMG_09_03 + #0.923 #0.053 &
          SVOSOMG_02_04 + #0.923 #0.053 SVOSOMG_03_04 + #0.923 #0.053 SVOSOMG_04_04 + &
          #0.923 #0.053 SVOSOMG_05_04 + #0.923 #0.053 SVOSOMG_06_04 + #0.923 &
          #0.053 SVOSOMG_07_04 + #0.923 #0.053 SVOSOMG_08_04 + #0.923 #0.053 &
          SVOSOMG_09_04 + #0.923 #0.053 SVOSOMG_03_05 + #0.923 #0.053 SVOSOMG_04_05 + &
          #0.923 #0.053 SVOSOMG_05_05 + #0.923 #0.053 SVOSOMG_06_05 + #0.923 &
          #0.053 SVOSOMG_07_05 + #0.923 #0.053 SVOSOMG_08_05 + #0.923 #0.053 &
          SVOSOMG_09_05
 S402.7 SVOSOMG_10_05 + OH = #0.059 #0.361 SVOSOMG_10_06 + #0.059 #0.639 &
          SVOSOMG_10_07 + #0.941 #0.047 SVOSOMG_01_01 + #0.941 #0.047 SVOSOMG_02_01 + &
          #0.941 #0.047 SVOSOMG_03_01 + #0.941 #0.047 SVOSOMG_04_01 + #0.941 &
          #0.047 SVOSOMG_05_01 + #0.941 #0.047 SVOSOMG_06_01 + #0.941 #0.047 &
          SVOSOMG_07_01 + #0.941 #0.000 SVOSOMG_08_01 + #0.941 #0.000 SVOSOMG_09_01 + &
          #0.941 #0.047 SVOSOMG_01_02 + #0.941 #0.047 SVOSOMG_02_02 + #0.941 &
          #0.047 SVOSOMG_03_02 + #0.941 #0.047 SVOSOMG_04_02 + #0.941 #0.047 &
          SVOSOMG_05_02 + #0.941 #0.047 SVOSOMG_06_02 + #0.941 #0.047 SVOSOMG_07_02 + &
          #0.941 #0.000 SVOSOMG_08_02 + #0.941 #0.000 SVOSOMG_09_02 + #0.941 &
          #0.047 SVOSOMG_02_03 + #0.941 #0.047 SVOSOMG_03_03 + #0.941 #0.047 &
          SVOSOMG_04_03 + #0.941 #0.047 SVOSOMG_05_03 + #0.941 #0.047 SVOSOMG_06_03 + &
          #0.941 #0.047 SVOSOMG_07_03 + #0.941 #0.047 SVOSOMG_08_03 + #0.941 &
          #0.000 SVOSOMG_09_03 + #0.941 #0.047 SVOSOMG_02_04 + #0.941 #0.047 &
          SVOSOMG_03_04 + #0.941 #0.047 SVOSOMG_04_04 + #0.941 #0.047 SVOSOMG_05_04 + &
          #0.941 #0.047 SVOSOMG_06_04 + #0.941 #0.047 SVOSOMG_07_04 + #0.941 &
          #0.047 SVOSOMG_08_04 + #0.941 #0.000 SVOSOMG_09_04 + #0.941 #0.047 &
          SVOSOMG_03_05 + #0.941 #0.047 SVOSOMG_04_05 + #0.941 #0.047 SVOSOMG_05_05 + &
          #0.941 #0.047 SVOSOMG_06_05 + #0.941 #0.047 SVOSOMG_07_05 + #0.941 &
          #0.047 SVOSOMG_08_05 + #0.941 #0.047 SVOSOMG_09_05 + #0.941 #0.047 &
          SVOSOMG_03_06 + #0.941 #0.047 SVOSOMG_04_06 + #0.941 #0.047 SVOSOMG_05_06 + &
          #0.941 #0.047 SVOSOMG_06_06 + #0.941 #0.047 SVOSOMG_07_06 + #0.941 &
          #0.047 SVOSOMG_08_06 + #0.941 #0.047 SVOSOMG_09_06
 S403.7 SVOSOMG_10_06 + OH = #0.043 #1.000 SVOSOMG_10_07 + #0.957 #0.043 &
          SVOSOMG_01_01 + #0.957 #0.043 SVOSOMG_02_01 + #0.957 #0.043 SVOSOMG_03_01 + &
          #0.957 #0.043 SVOSOMG_04_01 + #0.957 #0.043 SVOSOMG_05_01 + #0.957 &
          #0.043 SVOSOMG_06_01 + #0.957 #0.000 SVOSOMG_07_01 + #0.957 #0.000 &
          SVOSOMG_08_01 + #0.957 #0.000 SVOSOMG_09_01 + #0.957 #0.043 SVOSOMG_01_02 + &
          #0.957 #0.043 SVOSOMG_02_02 + #0.957 #0.043 SVOSOMG_03_02 + #0.957 &
          #0.043 SVOSOMG_04_02 + #0.957 #0.043 SVOSOMG_05_02 + #0.957 #0.043 &
          SVOSOMG_06_02 + #0.957 #0.043 SVOSOMG_07_02 + #0.957 #0.000 SVOSOMG_08_02 + &
          #0.957 #0.000 SVOSOMG_09_02 + #0.957 #0.043 SVOSOMG_02_03 + #0.957 &
          #0.043 SVOSOMG_03_03 + #0.957 #0.043 SVOSOMG_04_03 + #0.957 #0.043 &
          SVOSOMG_05_03 + #0.957 #0.043 SVOSOMG_06_03 + #0.957 #0.043 SVOSOMG_07_03 + &
          #0.957 #0.000 SVOSOMG_08_03 + #0.957 #0.000 SVOSOMG_09_03 + #0.957 &
          #0.043 SVOSOMG_02_04 + #0.957 #0.043 SVOSOMG_03_04 + #0.957 #0.043 &
          SVOSOMG_04_04 + #0.957 #0.043 SVOSOMG_05_04 + #0.957 #0.043 SVOSOMG_06_04 + &
          #0.957 #0.043 SVOSOMG_07_04 + #0.957 #0.043 SVOSOMG_08_04 + #0.957 &
          #0.000 SVOSOMG_09_04 + #0.957 #0.043 SVOSOMG_03_05 + #0.957 #0.043 &
          SVOSOMG_04_05 + #0.957 #0.043 SVOSOMG_05_05 + #0.957 #0.043 SVOSOMG_06_05 + &
          #0.957 #0.043 SVOSOMG_07_05 + #0.957 #0.043 SVOSOMG_08_05 + #0.957 &
          #0.000 SVOSOMG_09_05 + #0.957 #0.043 SVOSOMG_03_06 + #0.957 #0.043 &
          SVOSOMG_04_06 + #0.957 #0.043 SVOSOMG_05_06 + #0.957 #0.043 SVOSOMG_06_06 + &
          #0.957 #0.043 SVOSOMG_07_06 + #0.957 #0.043 SVOSOMG_08_06 + #0.957 &
          #0.043 SVOSOMG_09_06 + #0.957 #0.043 SVOSOMG_04_07 + #0.957 #0.043 &
          SVOSOMG_05_07 + #0.957 #0.043 SVOSOMG_06_07 + #0.957 #0.043 SVOSOMG_07_07 + &
          #0.957 #0.043 SVOSOMG_08_07 + #0.957 #0.043 SVOSOMG_09_07
 S404.7 SVOSOMG_10_07 + OH = #1.000 #0.045 SVOSOMG_01_01 + #1.000 #0.045 &
          SVOSOMG_02_01 + #1.000 #0.045 SVOSOMG_03_01 + #1.000 #0.045 SVOSOMG_04_01 + &
          #1.000 #0.045 SVOSOMG_05_01 + #1.000 #0.045 SVOSOMG_06_01 + #1.000 &
          #0.000 SVOSOMG_07_01 + #1.000 #0.000 SVOSOMG_08_01 + #1.000 #0.000 &
          SVOSOMG_09_01 + #1.000 #0.045 SVOSOMG_01_02 + #1.000 #0.045 SVOSOMG_02_02 + &
          #1.000 #0.045 SVOSOMG_03_02 + #1.000 #0.045 SVOSOMG_04_02 + #1.000 &
          #0.045 SVOSOMG_05_02 + #1.000 #0.045 SVOSOMG_06_02 + #1.000 #0.000 &
          SVOSOMG_07_02 + #1.000 #0.000 SVOSOMG_08_02 + #1.000 #0.000 SVOSOMG_09_02 + &
          #1.000 #0.045 SVOSOMG_02_03 + #1.000 #0.045 SVOSOMG_03_03 + #1.000 &
          #0.045 SVOSOMG_04_03 + #1.000 #0.045 SVOSOMG_05_03 + #1.000 #0.045 &
          SVOSOMG_06_03 + #1.000 #0.045 SVOSOMG_07_03 + #1.000 #0.000 SVOSOMG_08_03 + &
          #1.000 #0.000 SVOSOMG_09_03 + #1.000 #0.045 SVOSOMG_02_04 + #1.000 &
          #0.045 SVOSOMG_03_04 + #1.000 #0.045 SVOSOMG_04_04 + #1.000 #0.045 &
          SVOSOMG_05_04 + #1.000 #0.045 SVOSOMG_06_04 + #1.000 #0.045 SVOSOMG_07_04 + &
          #1.000 #0.000 SVOSOMG_08_04 + #1.000 #0.000 SVOSOMG_09_04 + #1.000 &
          #0.045 SVOSOMG_03_05 + #1.000 #0.045 SVOSOMG_04_05 + #1.000 #0.045 &
          SVOSOMG_05_05 + #1.000 #0.045 SVOSOMG_06_05 + #1.000 #0.045 SVOSOMG_07_05 + &
          #1.000 #0.045 SVOSOMG_08_05 + #1.000 #0.000 SVOSOMG_09_05 + #1.000 &
          #0.045 SVOSOMG_03_06 + #1.000 #0.045 SVOSOMG_04_06 + #1.000 #0.045 &
          SVOSOMG_05_06 + #1.000 #0.045 SVOSOMG_06_06 + #1.000 #0.045 SVOSOMG_07_06 + &
          #1.000 #0.045 SVOSOMG_08_06 + #1.000 #0.000 SVOSOMG_09_06 + #1.000 &
          #0.045 SVOSOMG_04_07 + #1.000 #0.045 SVOSOMG_05_07 + #1.000 #0.045 &
          SVOSOMG_06_07 + #1.000 #0.045 SVOSOMG_07_07 + #1.000 #0.045 SVOSOMG_08_07 + &
          #1.000 #0.045 SVOSOMG_09_07
 S405.7 SVOSOMG_11_01 + OH = #0.188 #0.361 SVOSOMG_11_02 + #0.188 #0.007 &
          SVOSOMG_11_03 + #0.188 #0.606 SVOSOMG_11_04 + #0.188 #0.026 SVOSOMG_11_05 + &
          #0.812 #0.100 SVOSOMG_01_01 + #0.812 #0.100 SVOSOMG_02_01 + #0.812 &
          #0.100 SVOSOMG_03_01 + #0.812 #0.100 SVOSOMG_04_01 + #0.812 #0.100 &
          SVOSOMG_05_01 + #0.812 #0.100 SVOSOMG_06_01 + #0.812 #0.100 SVOSOMG_07_01 + &
          #0.812 #0.100 SVOSOMG_08_01 + #0.812 #0.100 SVOSOMG_09_01 + #0.812 &
          #0.100 SVOSOMG_10_01 + #0.812 #0.100 SVOSOMG_01_02 + #0.812 #0.100 &
          SVOSOMG_02_02 + #0.812 #0.100 SVOSOMG_03_02 + #0.812 #0.100 SVOSOMG_04_02 + &
          #0.812 #0.100 SVOSOMG_05_02 + #0.812 #0.100 SVOSOMG_06_02 + #0.812 &
          #0.100 SVOSOMG_07_02 + #0.812 #0.100 SVOSOMG_08_02 + #0.812 #0.100 &
          SVOSOMG_09_02 + #0.812 #0.100 SVOSOMG_10_02
 S406.7 SVOSOMG_11_02 + OH = #0.138 #0.361 SVOSOMG_11_03 + #0.138 #0.007 &
          SVOSOMG_11_04 + #0.138 #0.606 SVOSOMG_11_05 + #0.138 #0.026 SVOSOMG_11_06 + &
          #0.862 #0.071 SVOSOMG_01_01 + #0.862 #0.071 SVOSOMG_02_01 + #0.862 &
          #0.071 SVOSOMG_03_01 + #0.862 #0.071 SVOSOMG_04_01 + #0.862 #0.071 &
          SVOSOMG_05_01 + #0.862 #0.071 SVOSOMG_06_01 + #0.862 #0.071 SVOSOMG_07_01 + &
          #0.862 #0.071 SVOSOMG_08_01 + #0.862 #0.071 SVOSOMG_09_01 + #0.862 &
          #0.000 SVOSOMG_10_01 + #0.862 #0.071 SVOSOMG_01_02 + #0.862 #0.071 &
          SVOSOMG_02_02 + #0.862 #0.071 SVOSOMG_03_02 + #0.862 #0.071 SVOSOMG_04_02 + &
          #0.862 #0.071 SVOSOMG_05_02 + #0.862 #0.071 SVOSOMG_06_02 + #0.862 &
          #0.071 SVOSOMG_07_02 + #0.862 #0.071 SVOSOMG_08_02 + #0.862 #0.071 &
          SVOSOMG_09_02 + #0.862 #0.071 SVOSOMG_10_02 + #0.862 #0.071 SVOSOMG_02_03 + &
          #0.862 #0.071 SVOSOMG_03_03 + #0.862 #0.071 SVOSOMG_04_03 + #0.862 &
          #0.071 SVOSOMG_05_03 + #0.862 #0.071 SVOSOMG_06_03 + #0.862 #0.071 &
          SVOSOMG_07_03 + #0.862 #0.071 SVOSOMG_08_03 + #0.862 #0.071 SVOSOMG_09_03 + &
          #0.862 #0.071 SVOSOMG_10_03
 S407.7 SVOSOMG_11_03 + OH = #0.107 #0.361 SVOSOMG_11_04 + #0.107 #0.007 &
          SVOSOMG_11_05 + #0.107 #0.606 SVOSOMG_11_06 + #0.107 #0.026 SVOSOMG_11_07 + &
          #0.893 #0.056 SVOSOMG_01_01 + #0.893 #0.056 SVOSOMG_02_01 + #0.893 &
          #0.056 SVOSOMG_03_01 + #0.893 #0.056 SVOSOMG_04_01 + #0.893 #0.056 &
          SVOSOMG_05_01 + #0.893 #0.056 SVOSOMG_06_01 + #0.893 #0.056 SVOSOMG_07_01 + &
          #0.893 #0.056 SVOSOMG_08_01 + #0.893 #0.056 SVOSOMG_09_01 + #0.893 &
          #0.000 SVOSOMG_10_01 + #0.893 #0.056 SVOSOMG_01_02 + #0.893 #0.056 &
          SVOSOMG_02_02 + #0.893 #0.056 SVOSOMG_03_02 + #0.893 #0.056 SVOSOMG_04_02 + &
          #0.893 #0.056 SVOSOMG_05_02 + #0.893 #0.056 SVOSOMG_06_02 + #0.893 &
          #0.056 SVOSOMG_07_02 + #0.893 #0.056 SVOSOMG_08_02 + #0.893 #0.056 &
          SVOSOMG_09_02 + #0.893 #0.000 SVOSOMG_10_02 + #0.893 #0.056 SVOSOMG_02_03 + &
          #0.893 #0.056 SVOSOMG_03_03 + #0.893 #0.056 SVOSOMG_04_03 + #0.893 &
          #0.056 SVOSOMG_05_03 + #0.893 #0.056 SVOSOMG_06_03 + #0.893 #0.056 &
          SVOSOMG_07_03 + #0.893 #0.056 SVOSOMG_08_03 + #0.893 #0.056 SVOSOMG_09_03 + &
          #0.893 #0.056 SVOSOMG_10_03 + #0.893 #0.056 SVOSOMG_02_04 + #0.893 &
          #0.056 SVOSOMG_03_04 + #0.893 #0.056 SVOSOMG_04_04 + #0.893 #0.056 &
          SVOSOMG_05_04 + #0.893 #0.056 SVOSOMG_06_04 + #0.893 #0.056 SVOSOMG_07_04 + &
          #0.893 #0.056 SVOSOMG_08_04 + #0.893 #0.056 SVOSOMG_09_04 + #0.893 &
          #0.056 SVOSOMG_10_04
 S408.7 SVOSOMG_11_04 + OH = #0.084 #0.361 SVOSOMG_11_05 + #0.084 #0.007 &
          SVOSOMG_11_06 + #0.084 #0.632 SVOSOMG_11_07 + #0.916 #0.048 SVOSOMG_01_01 + &
          #0.916 #0.048 SVOSOMG_02_01 + #0.916 #0.048 SVOSOMG_03_01 + #0.916 &
          #0.048 SVOSOMG_04_01 + #0.916 #0.048 SVOSOMG_05_01 + #0.916 #0.048 &
          SVOSOMG_06_01 + #0.916 #0.048 SVOSOMG_07_01 + #0.916 #0.048 SVOSOMG_08_01 + &
          #0.916 #0.000 SVOSOMG_09_01 + #0.916 #0.000 SVOSOMG_10_01 + #0.916 &
          #0.048 SVOSOMG_01_02 + #0.916 #0.048 SVOSOMG_02_02 + #0.916 #0.048 &
          SVOSOMG_03_02 + #0.916 #0.048 SVOSOMG_04_02 + #0.916 #0.048 SVOSOMG_05_02 + &
          #0.916 #0.048 SVOSOMG_06_02 + #0.916 #0.048 SVOSOMG_07_02 + #0.916 &
          #0.048 SVOSOMG_08_02 + #0.916 #0.048 SVOSOMG_09_02 + #0.916 #0.000 &
          SVOSOMG_10_02 + #0.916 #0.048 SVOSOMG_02_03 + #0.916 #0.048 SVOSOMG_03_03 + &
          #0.916 #0.048 SVOSOMG_04_03 + #0.916 #0.048 SVOSOMG_05_03 + #0.916 &
          #0.048 SVOSOMG_06_03 + #0.916 #0.048 SVOSOMG_07_03 + #0.916 #0.048 &
          SVOSOMG_08_03 + #0.916 #0.048 SVOSOMG_09_03 + #0.916 #0.000 SVOSOMG_10_03 + &
          #0.916 #0.048 SVOSOMG_02_04 + #0.916 #0.048 SVOSOMG_03_04 + #0.916 &
          #0.048 SVOSOMG_04_04 + #0.916 #0.048 SVOSOMG_05_04 + #0.916 #0.048 &
          SVOSOMG_06_04 + #0.916 #0.048 SVOSOMG_07_04 + #0.916 #0.048 SVOSOMG_08_04 + &
          #0.916 #0.048 SVOSOMG_09_04 + #0.916 #0.048 SVOSOMG_10_04 + #0.916 &
          #0.048 SVOSOMG_03_05 + #0.916 #0.048 SVOSOMG_04_05 + #0.916 #0.048 &
          SVOSOMG_05_05 + #0.916 #0.048 SVOSOMG_06_05 + #0.916 #0.048 SVOSOMG_07_05 + &
          #0.916 #0.048 SVOSOMG_08_05 + #0.916 #0.048 SVOSOMG_09_05 + #0.916 &
          #0.048 SVOSOMG_10_05
 S409.7 SVOSOMG_11_05 + OH = #0.066 #0.361 SVOSOMG_11_06 + #0.066 #0.639 &
          SVOSOMG_11_07 + #0.934 #0.042 SVOSOMG_01_01 + #0.934 #0.042 SVOSOMG_02_01 + &
          #0.934 #0.042 SVOSOMG_03_01 + #0.934 #0.042 SVOSOMG_04_01 + #0.934 &
          #0.042 SVOSOMG_05_01 + #0.934 #0.042 SVOSOMG_06_01 + #0.934 #0.042 &
          SVOSOMG_07_01 + #0.934 #0.042 SVOSOMG_08_01 + #0.934 #0.000 SVOSOMG_09_01 + &
          #0.934 #0.000 SVOSOMG_10_01 + #0.934 #0.042 SVOSOMG_01_02 + #0.934 &
          #0.042 SVOSOMG_02_02 + #0.934 #0.042 SVOSOMG_03_02 + #0.934 #0.042 &
          SVOSOMG_04_02 + #0.934 #0.042 SVOSOMG_05_02 + #0.934 #0.042 SVOSOMG_06_02 + &
          #0.934 #0.042 SVOSOMG_07_02 + #0.934 #0.042 SVOSOMG_08_02 + #0.934 &
          #0.000 SVOSOMG_09_02 + #0.934 #0.000 SVOSOMG_10_02 + #0.934 #0.042 &
          SVOSOMG_02_03 + #0.934 #0.042 SVOSOMG_03_03 + #0.934 #0.042 SVOSOMG_04_03 + &
          #0.934 #0.042 SVOSOMG_05_03 + #0.934 #0.042 SVOSOMG_06_03 + #0.934 &
          #0.042 SVOSOMG_07_03 + #0.934 #0.042 SVOSOMG_08_03 + #0.934 #0.042 &
          SVOSOMG_09_03 + #0.934 #0.000 SVOSOMG_10_03 + #0.934 #0.042 SVOSOMG_02_04 + &
          #0.934 #0.042 SVOSOMG_03_04 + #0.934 #0.042 SVOSOMG_04_04 + #0.934 &
          #0.042 SVOSOMG_05_04 + #0.934 #0.042 SVOSOMG_06_04 + #0.934 #0.042 &
          SVOSOMG_07_04 + #0.934 #0.042 SVOSOMG_08_04 + #0.934 #0.042 SVOSOMG_09_04 + &
          #0.934 #0.000 SVOSOMG_10_04 + #0.934 #0.042 SVOSOMG_03_05 + #0.934 &
          #0.042 SVOSOMG_04_05 + #0.934 #0.042 SVOSOMG_05_05 + #0.934 #0.042 &
          SVOSOMG_06_05 + #0.934 #0.042 SVOSOMG_07_05 + #0.934 #0.042 SVOSOMG_08_05 + &
          #0.934 #0.042 SVOSOMG_09_05 + #0.934 #0.042 SVOSOMG_10_05 + #0.934 &
          #0.042 SVOSOMG_03_06 + #0.934 #0.042 SVOSOMG_04_06 + #0.934 #0.042 &
          SVOSOMG_05_06 + #0.934 #0.042 SVOSOMG_06_06 + #0.934 #0.042 SVOSOMG_07_06 + &
          #0.934 #0.042 SVOSOMG_08_06 + #0.934 #0.042 SVOSOMG_09_06 + #0.934 &
          #0.042 SVOSOMG_10_06
 S410.7 SVOSOMG_11_06 + OH = #0.051 #1.000 SVOSOMG_11_07 + #0.949 #0.038 &
          SVOSOMG_01_01 + #0.949 #0.038 SVOSOMG_02_01 + #0.949 #0.038 SVOSOMG_03_01 + &
          #0.949 #0.038 SVOSOMG_04_01 + #0.949 #0.038 SVOSOMG_05_01 + #0.949 &
          #0.038 SVOSOMG_06_01 + #0.949 #0.038 SVOSOMG_07_01 + #0.949 #0.000 &
          SVOSOMG_08_01 + #0.949 #0.000 SVOSOMG_09_01 + #0.949 #0.000 SVOSOMG_10_01 + &
          #0.949 #0.038 SVOSOMG_01_02 + #0.949 #0.038 SVOSOMG_02_02 + #0.949 &
          #0.038 SVOSOMG_03_02 + #0.949 #0.038 SVOSOMG_04_02 + #0.949 #0.038 &
          SVOSOMG_05_02 + #0.949 #0.038 SVOSOMG_06_02 + #0.949 #0.038 SVOSOMG_07_02 + &
          #0.949 #0.038 SVOSOMG_08_02 + #0.949 #0.000 SVOSOMG_09_02 + #0.949 &
          #0.000 SVOSOMG_10_02 + #0.949 #0.038 SVOSOMG_02_03 + #0.949 #0.038 &
          SVOSOMG_03_03 + #0.949 #0.038 SVOSOMG_04_03 + #0.949 #0.038 SVOSOMG_05_03 + &
          #0.949 #0.038 SVOSOMG_06_03 + #0.949 #0.038 SVOSOMG_07_03 + #0.949 &
          #0.038 SVOSOMG_08_03 + #0.949 #0.000 SVOSOMG_09_03 + #0.949 #0.000 &
          SVOSOMG_10_03 + #0.949 #0.038 SVOSOMG_02_04 + #0.949 #0.038 SVOSOMG_03_04 + &
          #0.949 #0.038 SVOSOMG_04_04 + #0.949 #0.038 SVOSOMG_05_04 + #0.949 &
          #0.038 SVOSOMG_06_04 + #0.949 #0.038 SVOSOMG_07_04 + #0.949 #0.038 &
          SVOSOMG_08_04 + #0.949 #0.038 SVOSOMG_09_04 + #0.949 #0.000 SVOSOMG_10_04 + &
          #0.949 #0.038 SVOSOMG_03_05 + #0.949 #0.038 SVOSOMG_04_05 + #0.949 &
          #0.038 SVOSOMG_05_05 + #0.949 #0.038 SVOSOMG_06_05 + #0.949 #0.038 &
          SVOSOMG_07_05 + #0.949 #0.038 SVOSOMG_08_05 + #0.949 #0.038 SVOSOMG_09_05 + &
          #0.949 #0.000 SVOSOMG_10_05 + #0.949 #0.038 SVOSOMG_03_06 + #0.949 &
          #0.038 SVOSOMG_04_06 + #0.949 #0.038 SVOSOMG_05_06 + #0.949 #0.038 &
          SVOSOMG_06_06 + #0.949 #0.038 SVOSOMG_07_06 + #0.949 #0.038 SVOSOMG_08_06 + &
          #0.949 #0.038 SVOSOMG_09_06 + #0.949 #0.038 SVOSOMG_10_06 + #0.949 &
          #0.038 SVOSOMG_04_07 + #0.949 #0.038 SVOSOMG_05_07 + #0.949 #0.038 &
          SVOSOMG_06_07 + #0.949 #0.038 SVOSOMG_07_07 + #0.949 #0.038 SVOSOMG_08_07 + &
          #0.949 #0.038 SVOSOMG_09_07 + #0.949 #0.038 SVOSOMG_10_07
 S411.7 SVOSOMG_11_07 + OH = #1.000 #0.040 SVOSOMG_01_01 + #1.000 #0.040 &
          SVOSOMG_02_01 + #1.000 #0.040 SVOSOMG_03_01 + #1.000 #0.040 SVOSOMG_04_01 + &
          #1.000 #0.040 SVOSOMG_05_01 + #1.000 #0.040 SVOSOMG_06_01 + #1.000 &
          #0.040 SVOSOMG_07_01 + #1.000 #0.000 SVOSOMG_08_01 + #1.000 #0.000 &
          SVOSOMG_09_01 + #1.000 #0.000 SVOSOMG_10_01 + #1.000 #0.040 SVOSOMG_01_02 + &
          #1.000 #0.040 SVOSOMG_02_02 + #1.000 #0.040 SVOSOMG_03_02 + #1.000 &
          #0.040 SVOSOMG_04_02 + #1.000 #0.040 SVOSOMG_05_02 + #1.000 #0.040 &
          SVOSOMG_06_02 + #1.000 #0.040 SVOSOMG_07_02 + #1.000 #0.000 SVOSOMG_08_02 + &
          #1.000 #0.000 SVOSOMG_09_02 + #1.000 #0.000 SVOSOMG_10_02 + #1.000 &
          #0.040 SVOSOMG_02_03 + #1.000 #0.040 SVOSOMG_03_03 + #1.000 #0.040 &
          SVOSOMG_04_03 + #1.000 #0.040 SVOSOMG_05_03 + #1.000 #0.040 SVOSOMG_06_03 + &
          #1.000 #0.040 SVOSOMG_07_03 + #1.000 #0.040 SVOSOMG_08_03 + #1.000 &
          #0.000 SVOSOMG_09_03 + #1.000 #0.000 SVOSOMG_10_03 + #1.000 #0.040 &
          SVOSOMG_02_04 + #1.000 #0.040 SVOSOMG_03_04 + #1.000 #0.040 SVOSOMG_04_04 + &
          #1.000 #0.040 SVOSOMG_05_04 + #1.000 #0.040 SVOSOMG_06_04 + #1.000 &
          #0.040 SVOSOMG_07_04 + #1.000 #0.040 SVOSOMG_08_04 + #1.000 #0.000 &
          SVOSOMG_09_04 + #1.000 #0.000 SVOSOMG_10_04 + #1.000 #0.040 SVOSOMG_03_05 + &
          #1.000 #0.040 SVOSOMG_04_05 + #1.000 #0.040 SVOSOMG_05_05 + #1.000 &
          #0.040 SVOSOMG_06_05 + #1.000 #0.040 SVOSOMG_07_05 + #1.000 #0.040 &
          SVOSOMG_08_05 + #1.000 #0.040 SVOSOMG_09_05 + #1.000 #0.000 SVOSOMG_10_05 + &
          #1.000 #0.040 SVOSOMG_03_06 + #1.000 #0.040 SVOSOMG_04_06 + #1.000 &
          #0.040 SVOSOMG_05_06 + #1.000 #0.040 SVOSOMG_06_06 + #1.000 #0.040 &
          SVOSOMG_07_06 + #1.000 #0.040 SVOSOMG_08_06 + #1.000 #0.040 SVOSOMG_09_06 + &
          #1.000 #0.000 SVOSOMG_10_06 + #1.000 #0.040 SVOSOMG_04_07 + #1.000 &
          #0.040 SVOSOMG_05_07 + #1.000 #0.040 SVOSOMG_06_07 + #1.000 #0.040 &
          SVOSOMG_07_07 + #1.000 #0.040 SVOSOMG_08_07 + #1.000 #0.040 SVOSOMG_09_07 + &
          #1.000 #0.040 SVOSOMG_10_07
 S412.7 SVOSOMG_12_01 + OH = #0.194 #0.361 SVOSOMG_12_02 + #0.194 #0.007 &
          SVOSOMG_12_03 + #0.194 #0.606 SVOSOMG_12_04 + #0.194 #0.026 SVOSOMG_12_05 + &
          #0.806 #0.091 SVOSOMG_01_01 + #0.806 #0.091 SVOSOMG_02_01 + #0.806 &
          #0.091 SVOSOMG_03_01 + #0.806 #0.091 SVOSOMG_04_01 + #0.806 #0.091 &
          SVOSOMG_05_01 + #0.806 #0.091 SVOSOMG_06_01 + #0.806 #0.091 SVOSOMG_07_01 + &
          #0.806 #0.091 SVOSOMG_08_01 + #0.806 #0.091 SVOSOMG_09_01 + #0.806 &
          #0.091 SVOSOMG_10_01 + #0.806 #0.091 SVOSOMG_11_01 + #0.806 #0.091 &
          SVOSOMG_01_02 + #0.806 #0.091 SVOSOMG_02_02 + #0.806 #0.091 SVOSOMG_03_02 + &
          #0.806 #0.091 SVOSOMG_04_02 + #0.806 #0.091 SVOSOMG_05_02 + #0.806 &
          #0.091 SVOSOMG_06_02 + #0.806 #0.091 SVOSOMG_07_02 + #0.806 #0.091 &
          SVOSOMG_08_02 + #0.806 #0.091 SVOSOMG_09_02 + #0.806 #0.091 SVOSOMG_10_02 + &
          #0.806 #0.091 SVOSOMG_11_02
 S413.7 SVOSOMG_12_02 + OH = #0.144 #0.361 SVOSOMG_12_03 + #0.144 #0.007 &
          SVOSOMG_12_04 + #0.144 #0.606 SVOSOMG_12_05 + #0.144 #0.026 SVOSOMG_12_06 + &
          #0.856 #0.065 SVOSOMG_01_01 + #0.856 #0.065 SVOSOMG_02_01 + #0.856 &
          #0.065 SVOSOMG_03_01 + #0.856 #0.065 SVOSOMG_04_01 + #0.856 #0.065 &
          SVOSOMG_05_01 + #0.856 #0.065 SVOSOMG_06_01 + #0.856 #0.065 SVOSOMG_07_01 + &
          #0.856 #0.065 SVOSOMG_08_01 + #0.856 #0.065 SVOSOMG_09_01 + #0.856 &
          #0.065 SVOSOMG_10_01 + #0.856 #0.000 SVOSOMG_11_01 + #0.856 #0.065 &
          SVOSOMG_01_02 + #0.856 #0.065 SVOSOMG_02_02 + #0.856 #0.065 SVOSOMG_03_02 + &
          #0.856 #0.065 SVOSOMG_04_02 + #0.856 #0.065 SVOSOMG_05_02 + #0.856 &
          #0.065 SVOSOMG_06_02 + #0.856 #0.065 SVOSOMG_07_02 + #0.856 #0.065 &
          SVOSOMG_08_02 + #0.856 #0.065 SVOSOMG_09_02 + #0.856 #0.065 SVOSOMG_10_02 + &
          #0.856 #0.065 SVOSOMG_11_02 + #0.856 #0.065 SVOSOMG_02_03 + #0.856 &
          #0.065 SVOSOMG_03_03 + #0.856 #0.065 SVOSOMG_04_03 + #0.856 #0.065 &
          SVOSOMG_05_03 + #0.856 #0.065 SVOSOMG_06_03 + #0.856 #0.065 SVOSOMG_07_03 + &
          #0.856 #0.065 SVOSOMG_08_03 + #0.856 #0.065 SVOSOMG_09_03 + #0.856 &
          #0.065 SVOSOMG_10_03 + #0.856 #0.065 SVOSOMG_11_03
 S414.7 SVOSOMG_12_03 + OH = #0.114 #0.361 SVOSOMG_12_04 + #0.114 #0.007 &
          SVOSOMG_12_05 + #0.114 #0.606 SVOSOMG_12_06 + #0.114 #0.026 SVOSOMG_12_07 + &
          #0.886 #0.050 SVOSOMG_01_01 + #0.886 #0.050 SVOSOMG_02_01 + #0.886 &
          #0.050 SVOSOMG_03_01 + #0.886 #0.050 SVOSOMG_04_01 + #0.886 #0.050 &
          SVOSOMG_05_01 + #0.886 #0.050 SVOSOMG_06_01 + #0.886 #0.050 SVOSOMG_07_01 + &
          #0.886 #0.050 SVOSOMG_08_01 + #0.886 #0.050 SVOSOMG_09_01 + #0.886 &
          #0.050 SVOSOMG_10_01 + #0.886 #0.000 SVOSOMG_11_01 + #0.886 #0.050 &
          SVOSOMG_01_02 + #0.886 #0.050 SVOSOMG_02_02 + #0.886 #0.050 SVOSOMG_03_02 + &
          #0.886 #0.050 SVOSOMG_04_02 + #0.886 #0.050 SVOSOMG_05_02 + #0.886 &
          #0.050 SVOSOMG_06_02 + #0.886 #0.050 SVOSOMG_07_02 + #0.886 #0.050 &
          SVOSOMG_08_02 + #0.886 #0.050 SVOSOMG_09_02 + #0.886 #0.050 SVOSOMG_10_02 + &
          #0.886 #0.000 SVOSOMG_11_02 + #0.886 #0.050 SVOSOMG_02_03 + #0.886 &
          #0.050 SVOSOMG_03_03 + #0.886 #0.050 SVOSOMG_04_03 + #0.886 #0.050 &
          SVOSOMG_05_03 + #0.886 #0.050 SVOSOMG_06_03 + #0.886 #0.050 SVOSOMG_07_03 + &
          #0.886 #0.050 SVOSOMG_08_03 + #0.886 #0.050 SVOSOMG_09_03 + #0.886 &
          #0.050 SVOSOMG_10_03 + #0.886 #0.050 SVOSOMG_11_03 + #0.886 #0.050 &
          SVOSOMG_02_04 + #0.886 #0.050 SVOSOMG_03_04 + #0.886 #0.050 SVOSOMG_04_04 + &
          #0.886 #0.050 SVOSOMG_05_04 + #0.886 #0.050 SVOSOMG_06_04 + #0.886 &
          #0.050 SVOSOMG_07_04 + #0.886 #0.050 SVOSOMG_08_04 + #0.886 #0.050 &
          SVOSOMG_09_04 + #0.886 #0.050 SVOSOMG_10_04 + #0.886 #0.050 SVOSOMG_11_04
 S415.7 SVOSOMG_12_04 + OH = #0.091 #0.361 SVOSOMG_12_05 + #0.091 #0.007 &
          SVOSOMG_12_06 + #0.091 #0.632 SVOSOMG_12_07 + #0.909 #0.043 SVOSOMG_01_01 + &
          #0.909 #0.043 SVOSOMG_02_01 + #0.909 #0.043 SVOSOMG_03_01 + #0.909 &
          #0.043 SVOSOMG_04_01 + #0.909 #0.043 SVOSOMG_05_01 + #0.909 #0.043 &
          SVOSOMG_06_01 + #0.909 #0.043 SVOSOMG_07_01 + #0.909 #0.043 SVOSOMG_08_01 + &
          #0.909 #0.043 SVOSOMG_09_01 + #0.909 #0.000 SVOSOMG_10_01 + #0.909 &
          #0.000 SVOSOMG_11_01 + #0.909 #0.043 SVOSOMG_01_02 + #0.909 #0.043 &
          SVOSOMG_02_02 + #0.909 #0.043 SVOSOMG_03_02 + #0.909 #0.043 SVOSOMG_04_02 + &
          #0.909 #0.043 SVOSOMG_05_02 + #0.909 #0.043 SVOSOMG_06_02 + #0.909 &
          #0.043 SVOSOMG_07_02 + #0.909 #0.043 SVOSOMG_08_02 + #0.909 #0.043 &
          SVOSOMG_09_02 + #0.909 #0.043 SVOSOMG_10_02 + #0.909 #0.000 SVOSOMG_11_02 + &
          #0.909 #0.043 SVOSOMG_02_03 + #0.909 #0.043 SVOSOMG_03_03 + #0.909 &
          #0.043 SVOSOMG_04_03 + #0.909 #0.043 SVOSOMG_05_03 + #0.909 #0.043 &
          SVOSOMG_06_03 + #0.909 #0.043 SVOSOMG_07_03 + #0.909 #0.043 SVOSOMG_08_03 + &
          #0.909 #0.043 SVOSOMG_09_03 + #0.909 #0.043 SVOSOMG_10_03 + #0.909 &
          #0.000 SVOSOMG_11_03 + #0.909 #0.043 SVOSOMG_02_04 + #0.909 #0.043 &
          SVOSOMG_03_04 + #0.909 #0.043 SVOSOMG_04_04 + #0.909 #0.043 SVOSOMG_05_04 + &
          #0.909 #0.043 SVOSOMG_06_04 + #0.909 #0.043 SVOSOMG_07_04 + #0.909 &
          #0.043 SVOSOMG_08_04 + #0.909 #0.043 SVOSOMG_09_04 + #0.909 #0.043 &
          SVOSOMG_10_04 + #0.909 #0.043 SVOSOMG_11_04 + #0.909 #0.043 SVOSOMG_03_05 + &
          #0.909 #0.043 SVOSOMG_04_05 + #0.909 #0.043 SVOSOMG_05_05 + #0.909 &
          #0.043 SVOSOMG_06_05 + #0.909 #0.043 SVOSOMG_07_05 + #0.909 #0.043 &
          SVOSOMG_08_05 + #0.909 #0.043 SVOSOMG_09_05 + #0.909 #0.043 SVOSOMG_10_05 + &
          #0.909 #0.043 SVOSOMG_11_05
 S416.7 SVOSOMG_12_05 + OH = #0.073 #0.361 SVOSOMG_12_06 + #0.073 #0.639 &
          SVOSOMG_12_07 + #0.927 #0.037 SVOSOMG_01_01 + #0.927 #0.037 SVOSOMG_02_01 + &
          #0.927 #0.037 SVOSOMG_03_01 + #0.927 #0.037 SVOSOMG_04_01 + #0.927 &
          #0.037 SVOSOMG_05_01 + #0.927 #0.037 SVOSOMG_06_01 + #0.927 #0.037 &
          SVOSOMG_07_01 + #0.927 #0.037 SVOSOMG_08_01 + #0.927 #0.037 SVOSOMG_09_01 + &
          #0.927 #0.000 SVOSOMG_10_01 + #0.927 #0.000 SVOSOMG_11_01 + #0.927 &
          #0.037 SVOSOMG_01_02 + #0.927 #0.037 SVOSOMG_02_02 + #0.927 #0.037 &
          SVOSOMG_03_02 + #0.927 #0.037 SVOSOMG_04_02 + #0.927 #0.037 SVOSOMG_05_02 + &
          #0.927 #0.037 SVOSOMG_06_02 + #0.927 #0.037 SVOSOMG_07_02 + #0.927 &
          #0.037 SVOSOMG_08_02 + #0.927 #0.037 SVOSOMG_09_02 + #0.927 #0.000 &
          SVOSOMG_10_02 + #0.927 #0.000 SVOSOMG_11_02 + #0.927 #0.037 SVOSOMG_02_03 + &
          #0.927 #0.037 SVOSOMG_03_03 + #0.927 #0.037 SVOSOMG_04_03 + #0.927 &
          #0.037 SVOSOMG_05_03 + #0.927 #0.037 SVOSOMG_06_03 + #0.927 #0.037 &
          SVOSOMG_07_03 + #0.927 #0.037 SVOSOMG_08_03 + #0.927 #0.037 SVOSOMG_09_03 + &
          #0.927 #0.037 SVOSOMG_10_03 + #0.927 #0.000 SVOSOMG_11_03 + #0.927 &
          #0.037 SVOSOMG_02_04 + #0.927 #0.037 SVOSOMG_03_04 + #0.927 #0.037 &
          SVOSOMG_04_04 + #0.927 #0.037 SVOSOMG_05_04 + #0.927 #0.037 SVOSOMG_06_04 + &
          #0.927 #0.037 SVOSOMG_07_04 + #0.927 #0.037 SVOSOMG_08_04 + #0.927 &
          #0.037 SVOSOMG_09_04 + #0.927 #0.037 SVOSOMG_10_04 + #0.927 #0.000 &
          SVOSOMG_11_04 + #0.927 #0.037 SVOSOMG_03_05 + #0.927 #0.037 SVOSOMG_04_05 + &
          #0.927 #0.037 SVOSOMG_05_05 + #0.927 #0.037 SVOSOMG_06_05 + #0.927 &
          #0.037 SVOSOMG_07_05 + #0.927 #0.037 SVOSOMG_08_05 + #0.927 #0.037 &
          SVOSOMG_09_05 + #0.927 #0.037 SVOSOMG_10_05 + #0.927 #0.037 SVOSOMG_11_05 + &
          #0.927 #0.037 SVOSOMG_03_06 + #0.927 #0.037 SVOSOMG_04_06 + #0.927 &
          #0.037 SVOSOMG_05_06 + #0.927 #0.037 SVOSOMG_06_06 + #0.927 #0.037 &
          SVOSOMG_07_06 + #0.927 #0.037 SVOSOMG_08_06 + #0.927 #0.037 SVOSOMG_09_06 + &
          #0.927 #0.037 SVOSOMG_10_06 + #0.927 #0.037 SVOSOMG_11_06
 S417.7 SVOSOMG_12_06 + OH = #0.059 #1.000 SVOSOMG_12_07 + #0.941 #0.034 &
          SVOSOMG_01_01 + #0.941 #0.034 SVOSOMG_02_01 + #0.941 #0.034 SVOSOMG_03_01 + &
          #0.941 #0.034 SVOSOMG_04_01 + #0.941 #0.034 SVOSOMG_05_01 + #0.941 &
          #0.034 SVOSOMG_06_01 + #0.941 #0.034 SVOSOMG_07_01 + #0.941 #0.034 &
          SVOSOMG_08_01 + #0.941 #0.000 SVOSOMG_09_01 + #0.941 #0.000 SVOSOMG_10_01 + &
          #0.941 #0.000 SVOSOMG_11_01 + #0.941 #0.034 SVOSOMG_01_02 + #0.941 &
          #0.034 SVOSOMG_02_02 + #0.941 #0.034 SVOSOMG_03_02 + #0.941 #0.034 &
          SVOSOMG_04_02 + #0.941 #0.034 SVOSOMG_05_02 + #0.941 #0.034 SVOSOMG_06_02 + &
          #0.941 #0.034 SVOSOMG_07_02 + #0.941 #0.034 SVOSOMG_08_02 + #0.941 &
          #0.034 SVOSOMG_09_02 + #0.941 #0.000 SVOSOMG_10_02 + #0.941 #0.000 &
          SVOSOMG_11_02 + #0.941 #0.034 SVOSOMG_02_03 + #0.941 #0.034 SVOSOMG_03_03 + &
          #0.941 #0.034 SVOSOMG_04_03 + #0.941 #0.034 SVOSOMG_05_03 + #0.941 &
          #0.034 SVOSOMG_06_03 + #0.941 #0.034 SVOSOMG_07_03 + #0.941 #0.034 &
          SVOSOMG_08_03 + #0.941 #0.034 SVOSOMG_09_03 + #0.941 #0.000 SVOSOMG_10_03 + &
          #0.941 #0.000 SVOSOMG_11_03 + #0.941 #0.034 SVOSOMG_02_04 + #0.941 &
          #0.034 SVOSOMG_03_04 + #0.941 #0.034 SVOSOMG_04_04 + #0.941 #0.034 &
          SVOSOMG_05_04 + #0.941 #0.034 SVOSOMG_06_04 + #0.941 #0.034 SVOSOMG_07_04 + &
          #0.941 #0.034 SVOSOMG_08_04 + #0.941 #0.034 SVOSOMG_09_04 + #0.941 &
          #0.034 SVOSOMG_10_04 + #0.941 #0.000 SVOSOMG_11_04 + #0.941 #0.034 &
          SVOSOMG_03_05 + #0.941 #0.034 SVOSOMG_04_05 + #0.941 #0.034 SVOSOMG_05_05 + &
          #0.941 #0.034 SVOSOMG_06_05 + #0.941 #0.034 SVOSOMG_07_05 + #0.941 &
          #0.034 SVOSOMG_08_05 + #0.941 #0.034 SVOSOMG_09_05 + #0.941 #0.034 &
          SVOSOMG_10_05 + #0.941 #0.000 SVOSOMG_11_05 + #0.941 #0.034 SVOSOMG_03_06 + &
          #0.941 #0.034 SVOSOMG_04_06 + #0.941 #0.034 SVOSOMG_05_06 + #0.941 &
          #0.034 SVOSOMG_06_06 + #0.941 #0.034 SVOSOMG_07_06 + #0.941 #0.034 &
          SVOSOMG_08_06 + #0.941 #0.034 SVOSOMG_09_06 + #0.941 #0.034 SVOSOMG_10_06 + &
          #0.941 #0.034 SVOSOMG_11_06 + #0.941 #0.034 SVOSOMG_04_07 + #0.941 &
          #0.034 SVOSOMG_05_07 + #0.941 #0.034 SVOSOMG_06_07 + #0.941 #0.034 &
          SVOSOMG_07_07 + #0.941 #0.034 SVOSOMG_08_07 + #0.941 #0.034 SVOSOMG_09_07 + &
          #0.941 #0.034 SVOSOMG_10_07 + #0.941 #0.034 SVOSOMG_11_07
 S418.7 SVOSOMG_12_07 + OH = #1.000 #0.036 SVOSOMG_01_01 + #1.000 #0.036 &
          SVOSOMG_02_01 + #1.000 #0.036 SVOSOMG_03_01 + #1.000 #0.036 SVOSOMG_04_01 + &
          #1.000 #0.036 SVOSOMG_05_01 + #1.000 #0.036 SVOSOMG_06_01 + #1.000 &
          #0.036 SVOSOMG_07_01 + #1.000 #0.036 SVOSOMG_08_01 + #1.000 #0.000 &
          SVOSOMG_09_01 + #1.000 #0.000 SVOSOMG_10_01 + #1.000 #0.000 SVOSOMG_11_01 + &
          #1.000 #0.036 SVOSOMG_01_02 + #1.000 #0.036 SVOSOMG_02_02 + #1.000 &
          #0.036 SVOSOMG_03_02 + #1.000 #0.036 SVOSOMG_04_02 + #1.000 #0.036 &
          SVOSOMG_05_02 + #1.000 #0.036 SVOSOMG_06_02 + #1.000 #0.036 SVOSOMG_07_02 + &
          #1.000 #0.036 SVOSOMG_08_02 + #1.000 #0.000 SVOSOMG_09_02 + #1.000 &
          #0.000 SVOSOMG_10_02 + #1.000 #0.000 SVOSOMG_11_02 + #1.000 #0.036 &
          SVOSOMG_02_03 + #1.000 #0.036 SVOSOMG_03_03 + #1.000 #0.036 SVOSOMG_04_03 + &
          #1.000 #0.036 SVOSOMG_05_03 + #1.000 #0.036 SVOSOMG_06_03 + #1.000 &
          #0.036 SVOSOMG_07_03 + #1.000 #0.036 SVOSOMG_08_03 + #1.000 #0.036 &
          SVOSOMG_09_03 + #1.000 #0.000 SVOSOMG_10_03 + #1.000 #0.000 SVOSOMG_11_03 + &
          #1.000 #0.036 SVOSOMG_02_04 + #1.000 #0.036 SVOSOMG_03_04 + #1.000 &
          #0.036 SVOSOMG_04_04 + #1.000 #0.036 SVOSOMG_05_04 + #1.000 #0.036 &
          SVOSOMG_06_04 + #1.000 #0.036 SVOSOMG_07_04 + #1.000 #0.036 SVOSOMG_08_04 + &
          #1.000 #0.036 SVOSOMG_09_04 + #1.000 #0.000 SVOSOMG_10_04 + #1.000 &
          #0.000 SVOSOMG_11_04 + #1.000 #0.036 SVOSOMG_03_05 + #1.000 #0.036 &
          SVOSOMG_04_05 + #1.000 #0.036 SVOSOMG_05_05 + #1.000 #0.036 SVOSOMG_06_05 + &
          #1.000 #0.036 SVOSOMG_07_05 + #1.000 #0.036 SVOSOMG_08_05 + #1.000 &
          #0.036 SVOSOMG_09_05 + #1.000 #0.036 SVOSOMG_10_05 + #1.000 #0.000 &
          SVOSOMG_11_05 + #1.000 #0.036 SVOSOMG_03_06 + #1.000 #0.036 SVOSOMG_04_06 + &
          #1.000 #0.036 SVOSOMG_05_06 + #1.000 #0.036 SVOSOMG_06_06 + #1.000 &
          #0.036 SVOSOMG_07_06 + #1.000 #0.036 SVOSOMG_08_06 + #1.000 #0.036 &
          SVOSOMG_09_06 + #1.000 #0.036 SVOSOMG_10_06 + #1.000 #0.000 SVOSOMG_11_06 + &
          #1.000 #0.036 SVOSOMG_04_07 + #1.000 #0.036 SVOSOMG_05_07 + #1.000 &
          #0.036 SVOSOMG_06_07 + #1.000 #0.036 SVOSOMG_07_07 + #1.000 #0.036 &
          SVOSOMG_08_07 + #1.000 #0.036 SVOSOMG_09_07 + #1.000 #0.036 SVOSOMG_10_07 + &
          #1.000 #0.036 SVOSOMG_11_07
 S419.7 SVOSOMG_13_01 + OH = #0.200 #0.361 SVOSOMG_13_02 + #0.200 #0.007 &
          SVOSOMG_13_03 + #0.200 #0.606 SVOSOMG_13_04 + #0.200 #0.026 SVOSOMG_13_05 + &
          #0.800 #0.083 SVOSOMG_01_01 + #0.800 #0.083 SVOSOMG_02_01 + #0.800 &
          #0.083 SVOSOMG_03_01 + #0.800 #0.083 SVOSOMG_04_01 + #0.800 #0.083 &
          SVOSOMG_05_01 + #0.800 #0.083 SVOSOMG_06_01 + #0.800 #0.083 SVOSOMG_07_01 + &
          #0.800 #0.083 SVOSOMG_08_01 + #0.800 #0.083 SVOSOMG_09_01 + #0.800 &
          #0.083 SVOSOMG_10_01 + #0.800 #0.083 SVOSOMG_11_01 + #0.800 #0.083 &
          SVOSOMG_12_01 + #0.800 #0.083 SVOSOMG_01_02 + #0.800 #0.083 SVOSOMG_02_02 + &
          #0.800 #0.083 SVOSOMG_03_02 + #0.800 #0.083 SVOSOMG_04_02 + #0.800 &
          #0.083 SVOSOMG_05_02 + #0.800 #0.083 SVOSOMG_06_02 + #0.800 #0.083 &
          SVOSOMG_07_02 + #0.800 #0.083 SVOSOMG_08_02 + #0.800 #0.083 SVOSOMG_09_02 + &
          #0.800 #0.083 SVOSOMG_10_02 + #0.800 #0.083 SVOSOMG_11_02 + #0.800 &
          #0.083 SVOSOMG_12_02
 S420.7 SVOSOMG_13_02 + OH = #0.150 #0.361 SVOSOMG_13_03 + #0.150 #0.007 &
          SVOSOMG_13_04 + #0.150 #0.606 SVOSOMG_13_05 + #0.150 #0.026 SVOSOMG_13_06 + &
          #0.850 #0.057 SVOSOMG_01_01 + #0.850 #0.057 SVOSOMG_02_01 + #0.850 &
          #0.057 SVOSOMG_03_01 + #0.850 #0.057 SVOSOMG_04_01 + #0.850 #0.057 &
          SVOSOMG_05_01 + #0.850 #0.057 SVOSOMG_06_01 + #0.850 #0.057 SVOSOMG_07_01 + &
          #0.850 #0.057 SVOSOMG_08_01 + #0.850 #0.057 SVOSOMG_09_01 + #0.850 &
          #0.057 SVOSOMG_10_01 + #0.850 #0.057 SVOSOMG_11_01 + #0.850 #0.000 &
          SVOSOMG_12_01 + #0.850 #0.057 SVOSOMG_01_02 + #0.850 #0.057 SVOSOMG_02_02 + &
          #0.850 #0.057 SVOSOMG_03_02 + #0.850 #0.057 SVOSOMG_04_02 + #0.850 &
          #0.057 SVOSOMG_05_02 + #0.850 #0.057 SVOSOMG_06_02 + #0.850 #0.057 &
          SVOSOMG_07_02 + #0.850 #0.057 SVOSOMG_08_02 + #0.850 #0.057 SVOSOMG_09_02 + &
          #0.850 #0.057 SVOSOMG_10_02 + #0.850 #0.057 SVOSOMG_11_02 + #0.850 &
          #0.057 SVOSOMG_12_02 + #0.850 #0.057 SVOSOMG_02_03 + #0.850 #0.057 &
          SVOSOMG_03_03 + #0.850 #0.057 SVOSOMG_04_03 + #0.850 #0.057 SVOSOMG_05_03 + &
          #0.850 #0.057 SVOSOMG_06_03 + #0.850 #0.057 SVOSOMG_07_03 + #0.850 &
          #0.057 SVOSOMG_08_03 + #0.850 #0.057 SVOSOMG_09_03 + #0.850 #0.057 &
          SVOSOMG_10_03 + #0.850 #0.057 SVOSOMG_11_03 + #0.850 #0.057 SVOSOMG_12_03
 S421.7 SVOSOMG_13_03 + OH = #0.120 #0.361 SVOSOMG_13_04 + #0.120 #0.007 &
          SVOSOMG_13_05 + #0.120 #0.606 SVOSOMG_13_06 + #0.120 #0.026 SVOSOMG_13_07 + &
          #0.880 #0.044 SVOSOMG_01_01 + #0.880 #0.044 SVOSOMG_02_01 + #0.880 &
          #0.044 SVOSOMG_03_01 + #0.880 #0.044 SVOSOMG_04_01 + #0.880 #0.044 &
          SVOSOMG_05_01 + #0.880 #0.044 SVOSOMG_06_01 + #0.880 #0.044 SVOSOMG_07_01 + &
          #0.880 #0.044 SVOSOMG_08_01 + #0.880 #0.044 SVOSOMG_09_01 + #0.880 &
          #0.044 SVOSOMG_10_01 + #0.880 #0.044 SVOSOMG_11_01 + #0.880 #0.000 &
          SVOSOMG_12_01 + #0.880 #0.044 SVOSOMG_01_02 + #0.880 #0.044 SVOSOMG_02_02 + &
          #0.880 #0.044 SVOSOMG_03_02 + #0.880 #0.044 SVOSOMG_04_02 + #0.880 &
          #0.044 SVOSOMG_05_02 + #0.880 #0.044 SVOSOMG_06_02 + #0.880 #0.044 &
          SVOSOMG_07_02 + #0.880 #0.044 SVOSOMG_08_02 + #0.880 #0.044 SVOSOMG_09_02 + &
          #0.880 #0.044 SVOSOMG_10_02 + #0.880 #0.044 SVOSOMG_11_02 + #0.880 &
          #0.000 SVOSOMG_12_02 + #0.880 #0.044 SVOSOMG_02_03 + #0.880 #0.044 &
          SVOSOMG_03_03 + #0.880 #0.044 SVOSOMG_04_03 + #0.880 #0.044 SVOSOMG_05_03 + &
          #0.880 #0.044 SVOSOMG_06_03 + #0.880 #0.044 SVOSOMG_07_03 + #0.880 &
          #0.044 SVOSOMG_08_03 + #0.880 #0.044 SVOSOMG_09_03 + #0.880 #0.044 &
          SVOSOMG_10_03 + #0.880 #0.044 SVOSOMG_11_03 + #0.880 #0.044 SVOSOMG_12_03 + &
          #0.880 #0.044 SVOSOMG_02_04 + #0.880 #0.044 SVOSOMG_03_04 + #0.880 &
          #0.044 SVOSOMG_04_04 + #0.880 #0.044 SVOSOMG_05_04 + #0.880 #0.044 &
          SVOSOMG_06_04 + #0.880 #0.044 SVOSOMG_07_04 + #0.880 #0.044 SVOSOMG_08_04 + &
          #0.880 #0.044 SVOSOMG_09_04 + #0.880 #0.044 SVOSOMG_10_04 + #0.880 &
          #0.044 SVOSOMG_11_04 + #0.880 #0.044 SVOSOMG_12_04
 S422.7 SVOSOMG_13_04 + OH = #0.097 #0.361 SVOSOMG_13_05 + #0.097 #0.007 &
          SVOSOMG_13_06 + #0.097 #0.632 SVOSOMG_13_07 + #0.903 #0.037 SVOSOMG_01_01 + &
          #0.903 #0.037 SVOSOMG_02_01 + #0.903 #0.037 SVOSOMG_03_01 + #0.903 &
          #0.037 SVOSOMG_04_01 + #0.903 #0.037 SVOSOMG_05_01 + #0.903 #0.037 &
          SVOSOMG_06_01 + #0.903 #0.037 SVOSOMG_07_01 + #0.903 #0.037 SVOSOMG_08_01 + &
          #0.903 #0.037 SVOSOMG_09_01 + #0.903 #0.037 SVOSOMG_10_01 + #0.903 &
          #0.000 SVOSOMG_11_01 + #0.903 #0.000 SVOSOMG_12_01 + #0.903 #0.037 &
          SVOSOMG_01_02 + #0.903 #0.037 SVOSOMG_02_02 + #0.903 #0.037 SVOSOMG_03_02 + &
          #0.903 #0.037 SVOSOMG_04_02 + #0.903 #0.037 SVOSOMG_05_02 + #0.903 &
          #0.037 SVOSOMG_06_02 + #0.903 #0.037 SVOSOMG_07_02 + #0.903 #0.037 &
          SVOSOMG_08_02 + #0.903 #0.037 SVOSOMG_09_02 + #0.903 #0.037 SVOSOMG_10_02 + &
          #0.903 #0.037 SVOSOMG_11_02 + #0.903 #0.000 SVOSOMG_12_02 + #0.903 &
          #0.037 SVOSOMG_02_03 + #0.903 #0.037 SVOSOMG_03_03 + #0.903 #0.037 &
          SVOSOMG_04_03 + #0.903 #0.037 SVOSOMG_05_03 + #0.903 #0.037 SVOSOMG_06_03 + &
          #0.903 #0.037 SVOSOMG_07_03 + #0.903 #0.037 SVOSOMG_08_03 + #0.903 &
          #0.037 SVOSOMG_09_03 + #0.903 #0.037 SVOSOMG_10_03 + #0.903 #0.037 &
          SVOSOMG_11_03 + #0.903 #0.000 SVOSOMG_12_03 + #0.903 #0.037 SVOSOMG_02_04 + &
          #0.903 #0.037 SVOSOMG_03_04 + #0.903 #0.037 SVOSOMG_04_04 + #0.903 &
          #0.037 SVOSOMG_05_04 + #0.903 #0.037 SVOSOMG_06_04 + #0.903 #0.037 &
          SVOSOMG_07_04 + #0.903 #0.037 SVOSOMG_08_04 + #0.903 #0.037 SVOSOMG_09_04 + &
          #0.903 #0.037 SVOSOMG_10_04 + #0.903 #0.037 SVOSOMG_11_04 + #0.903 &
          #0.037 SVOSOMG_12_04 + #0.903 #0.037 SVOSOMG_03_05 + #0.903 #0.037 &
          SVOSOMG_04_05 + #0.903 #0.037 SVOSOMG_05_05 + #0.903 #0.037 SVOSOMG_06_05 + &
          #0.903 #0.037 SVOSOMG_07_05 + #0.903 #0.037 SVOSOMG_08_05 + #0.903 &
          #0.037 SVOSOMG_09_05 + #0.903 #0.037 SVOSOMG_10_05 + #0.903 #0.037 &
          SVOSOMG_11_05 + #0.903 #0.037 SVOSOMG_12_05
 S423.7 SVOSOMG_13_05 + OH = #0.080 #0.361 SVOSOMG_13_06 + #0.080 #0.639 &
          SVOSOMG_13_07 + #0.920 #0.032 SVOSOMG_01_01 + #0.920 #0.032 SVOSOMG_02_01 + &
          #0.920 #0.032 SVOSOMG_03_01 + #0.920 #0.032 SVOSOMG_04_01 + #0.920 &
          #0.032 SVOSOMG_05_01 + #0.920 #0.032 SVOSOMG_06_01 + #0.920 #0.032 &
          SVOSOMG_07_01 + #0.920 #0.032 SVOSOMG_08_01 + #0.920 #0.032 SVOSOMG_09_01 + &
          #0.920 #0.032 SVOSOMG_10_01 + #0.920 #0.000 SVOSOMG_11_01 + #0.920 &
          #0.000 SVOSOMG_12_01 + #0.920 #0.032 SVOSOMG_01_02 + #0.920 #0.032 &
          SVOSOMG_02_02 + #0.920 #0.032 SVOSOMG_03_02 + #0.920 #0.032 SVOSOMG_04_02 + &
          #0.920 #0.032 SVOSOMG_05_02 + #0.920 #0.032 SVOSOMG_06_02 + #0.920 &
          #0.032 SVOSOMG_07_02 + #0.920 #0.032 SVOSOMG_08_02 + #0.920 #0.032 &
          SVOSOMG_09_02 + #0.920 #0.032 SVOSOMG_10_02 + #0.920 #0.000 SVOSOMG_11_02 + &
          #0.920 #0.000 SVOSOMG_12_02 + #0.920 #0.032 SVOSOMG_02_03 + #0.920 &
          #0.032 SVOSOMG_03_03 + #0.920 #0.032 SVOSOMG_04_03 + #0.920 #0.032 &
          SVOSOMG_05_03 + #0.920 #0.032 SVOSOMG_06_03 + #0.920 #0.032 SVOSOMG_07_03 + &
          #0.920 #0.032 SVOSOMG_08_03 + #0.920 #0.032 SVOSOMG_09_03 + #0.920 &
          #0.032 SVOSOMG_10_03 + #0.920 #0.032 SVOSOMG_11_03 + #0.920 #0.000 &
          SVOSOMG_12_03 + #0.920 #0.032 SVOSOMG_02_04 + #0.920 #0.032 SVOSOMG_03_04 + &
          #0.920 #0.032 SVOSOMG_04_04 + #0.920 #0.032 SVOSOMG_05_04 + #0.920 &
          #0.032 SVOSOMG_06_04 + #0.920 #0.032 SVOSOMG_07_04 + #0.920 #0.032 &
          SVOSOMG_08_04 + #0.920 #0.032 SVOSOMG_09_04 + #0.920 #0.032 SVOSOMG_10_04 + &
          #0.920 #0.032 SVOSOMG_11_04 + #0.920 #0.000 SVOSOMG_12_04 + #0.920 &
          #0.032 SVOSOMG_03_05 + #0.920 #0.032 SVOSOMG_04_05 + #0.920 #0.032 &
          SVOSOMG_05_05 + #0.920 #0.032 SVOSOMG_06_05 + #0.920 #0.032 SVOSOMG_07_05 + &
          #0.920 #0.032 SVOSOMG_08_05 + #0.920 #0.032 SVOSOMG_09_05 + #0.920 &
          #0.032 SVOSOMG_10_05 + #0.920 #0.032 SVOSOMG_11_05 + #0.920 #0.032 &
          SVOSOMG_12_05 + #0.920 #0.032 SVOSOMG_03_06 + #0.920 #0.032 SVOSOMG_04_06 + &
          #0.920 #0.032 SVOSOMG_05_06 + #0.920 #0.032 SVOSOMG_06_06 + #0.920 &
          #0.032 SVOSOMG_07_06 + #0.920 #0.032 SVOSOMG_08_06 + #0.920 #0.032 &
          SVOSOMG_09_06 + #0.920 #0.032 SVOSOMG_10_06 + #0.920 #0.032 SVOSOMG_11_06 + &
          #0.920 #0.032 SVOSOMG_12_06
 S424.7 SVOSOMG_13_06 + OH = #0.065 #1.000 SVOSOMG_13_07 + #0.935 #0.029 &
          SVOSOMG_01_01 + #0.935 #0.029 SVOSOMG_02_01 + #0.935 #0.029 SVOSOMG_03_01 + &
          #0.935 #0.029 SVOSOMG_04_01 + #0.935 #0.029 SVOSOMG_05_01 + #0.935 &
          #0.029 SVOSOMG_06_01 + #0.935 #0.029 SVOSOMG_07_01 + #0.935 #0.029 &
          SVOSOMG_08_01 + #0.935 #0.029 SVOSOMG_09_01 + #0.935 #0.000 SVOSOMG_10_01 + &
          #0.935 #0.000 SVOSOMG_11_01 + #0.935 #0.000 SVOSOMG_12_01 + #0.935 &
          #0.029 SVOSOMG_01_02 + #0.935 #0.029 SVOSOMG_02_02 + #0.935 #0.029 &
          SVOSOMG_03_02 + #0.935 #0.029 SVOSOMG_04_02 + #0.935 #0.029 SVOSOMG_05_02 + &
          #0.935 #0.029 SVOSOMG_06_02 + #0.935 #0.029 SVOSOMG_07_02 + #0.935 &
          #0.029 SVOSOMG_08_02 + #0.935 #0.029 SVOSOMG_09_02 + #0.935 #0.029 &
          SVOSOMG_10_02 + #0.935 #0.000 SVOSOMG_11_02 + #0.935 #0.000 SVOSOMG_12_02 + &
          #0.935 #0.029 SVOSOMG_02_03 + #0.935 #0.029 SVOSOMG_03_03 + #0.935 &
          #0.029 SVOSOMG_04_03 + #0.935 #0.029 SVOSOMG_05_03 + #0.935 #0.029 &
          SVOSOMG_06_03 + #0.935 #0.029 SVOSOMG_07_03 + #0.935 #0.029 SVOSOMG_08_03 + &
          #0.935 #0.029 SVOSOMG_09_03 + #0.935 #0.029 SVOSOMG_10_03 + #0.935 &
          #0.000 SVOSOMG_11_03 + #0.935 #0.000 SVOSOMG_12_03 + #0.935 #0.029 &
          SVOSOMG_02_04 + #0.935 #0.029 SVOSOMG_03_04 + #0.935 #0.029 SVOSOMG_04_04 + &
          #0.935 #0.029 SVOSOMG_05_04 + #0.935 #0.029 SVOSOMG_06_04 + #0.935 &
          #0.029 SVOSOMG_07_04 + #0.935 #0.029 SVOSOMG_08_04 + #0.935 #0.029 &
          SVOSOMG_09_04 + #0.935 #0.029 SVOSOMG_10_04 + #0.935 #0.029 SVOSOMG_11_04 + &
          #0.935 #0.000 SVOSOMG_12_04 + #0.935 #0.029 SVOSOMG_03_05 + #0.935 &
          #0.029 SVOSOMG_04_05 + #0.935 #0.029 SVOSOMG_05_05 + #0.935 #0.029 &
          SVOSOMG_06_05 + #0.935 #0.029 SVOSOMG_07_05 + #0.935 #0.029 SVOSOMG_08_05 + &
          #0.935 #0.029 SVOSOMG_09_05 + #0.935 #0.029 SVOSOMG_10_05 + #0.935 &
          #0.029 SVOSOMG_11_05 + #0.935 #0.000 SVOSOMG_12_05 + #0.935 #0.029 &
          SVOSOMG_03_06 + #0.935 #0.029 SVOSOMG_04_06 + #0.935 #0.029 SVOSOMG_05_06 + &
          #0.935 #0.029 SVOSOMG_06_06 + #0.935 #0.029 SVOSOMG_07_06 + #0.935 &
          #0.029 SVOSOMG_08_06 + #0.935 #0.029 SVOSOMG_09_06 + #0.935 #0.029 &
          SVOSOMG_10_06 + #0.935 #0.029 SVOSOMG_11_06 + #0.935 #0.029 SVOSOMG_12_06 + &
          #0.935 #0.029 SVOSOMG_04_07 + #0.935 #0.029 SVOSOMG_05_07 + #0.935 &
          #0.029 SVOSOMG_06_07 + #0.935 #0.029 SVOSOMG_07_07 + #0.935 #0.029 &
          SVOSOMG_08_07 + #0.935 #0.029 SVOSOMG_09_07 + #0.935 #0.029 SVOSOMG_10_07 + &
          #0.935 #0.029 SVOSOMG_11_07 + #0.935 #0.029 SVOSOMG_12_07
 S425.7 SVOSOMG_13_07 + OH = #1.000 #0.030 SVOSOMG_01_01 + #1.000 #0.030 &
          SVOSOMG_02_01 + #1.000 #0.030 SVOSOMG_03_01 + #1.000 #0.030 SVOSOMG_04_01 + &
          #1.000 #0.030 SVOSOMG_05_01 + #1.000 #0.030 SVOSOMG_06_01 + #1.000 &
          #0.030 SVOSOMG_07_01 + #1.000 #0.030 SVOSOMG_08_01 + #1.000 #0.030 &
          SVOSOMG_09_01 + #1.000 #0.000 SVOSOMG_10_01 + #1.000 #0.000 SVOSOMG_11_01 + &
          #1.000 #0.000 SVOSOMG_12_01 + #1.000 #0.030 SVOSOMG_01_02 + #1.000 &
          #0.030 SVOSOMG_02_02 + #1.000 #0.030 SVOSOMG_03_02 + #1.000 #0.030 &
          SVOSOMG_04_02 + #1.000 #0.030 SVOSOMG_05_02 + #1.000 #0.030 SVOSOMG_06_02 + &
          #1.000 #0.030 SVOSOMG_07_02 + #1.000 #0.030 SVOSOMG_08_02 + #1.000 &
          #0.030 SVOSOMG_09_02 + #1.000 #0.000 SVOSOMG_10_02 + #1.000 #0.000 &
          SVOSOMG_11_02 + #1.000 #0.000 SVOSOMG_12_02 + #1.000 #0.030 SVOSOMG_02_03 + &
          #1.000 #0.030 SVOSOMG_03_03 + #1.000 #0.030 SVOSOMG_04_03 + #1.000 &
          #0.030 SVOSOMG_05_03 + #1.000 #0.030 SVOSOMG_06_03 + #1.000 #0.030 &
          SVOSOMG_07_03 + #1.000 #0.030 SVOSOMG_08_03 + #1.000 #0.030 SVOSOMG_09_03 + &
          #1.000 #0.030 SVOSOMG_10_03 + #1.000 #0.000 SVOSOMG_11_03 + #1.000 &
          #0.000 SVOSOMG_12_03 + #1.000 #0.030 SVOSOMG_02_04 + #1.000 #0.030 &
          SVOSOMG_03_04 + #1.000 #0.030 SVOSOMG_04_04 + #1.000 #0.030 SVOSOMG_05_04 + &
          #1.000 #0.030 SVOSOMG_06_04 + #1.000 #0.030 SVOSOMG_07_04 + #1.000 &
          #0.030 SVOSOMG_08_04 + #1.000 #0.030 SVOSOMG_09_04 + #1.000 #0.030 &
          SVOSOMG_10_04 + #1.000 #0.000 SVOSOMG_11_04 + #1.000 #0.000 SVOSOMG_12_04 + &
          #1.000 #0.030 SVOSOMG_03_05 + #1.000 #0.030 SVOSOMG_04_05 + #1.000 &
          #0.030 SVOSOMG_05_05 + #1.000 #0.030 SVOSOMG_06_05 + #1.000 #0.030 &
          SVOSOMG_07_05 + #1.000 #0.030 SVOSOMG_08_05 + #1.000 #0.030 SVOSOMG_09_05 + &
          #1.000 #0.030 SVOSOMG_10_05 + #1.000 #0.030 SVOSOMG_11_05 + #1.000 &
          #0.000 SVOSOMG_12_05 + #1.000 #0.030 SVOSOMG_03_06 + #1.000 #0.030 &
          SVOSOMG_04_06 + #1.000 #0.030 SVOSOMG_05_06 + #1.000 #0.030 SVOSOMG_06_06 + &
          #1.000 #0.030 SVOSOMG_07_06 + #1.000 #0.030 SVOSOMG_08_06 + #1.000 &
          #0.030 SVOSOMG_09_06 + #1.000 #0.030 SVOSOMG_10_06 + #1.000 #0.030 &
          SVOSOMG_11_06 + #1.000 #0.000 SVOSOMG_12_06 + #1.000 #0.030 SVOSOMG_04_07 + &
          #1.000 #0.030 SVOSOMG_05_07 + #1.000 #0.030 SVOSOMG_06_07 + #1.000 &
          #0.030 SVOSOMG_07_07 + #1.000 #0.030 SVOSOMG_08_07 + #1.000 #0.030 &
          SVOSOMG_09_07 + #1.000 #0.030 SVOSOMG_10_07 + #1.000 #0.030 SVOSOMG_11_07 + &
          #1.000 #0.030 SVOSOMG_12_07
 S426.7 SVOSOMG_14_01 + OH = #0.205 #0.361 SVOSOMG_14_02 + #0.205 #0.007 &
          SVOSOMG_14_03 + #0.205 #0.606 SVOSOMG_14_04 + #0.205 #0.026 SVOSOMG_14_05 + &
          #0.795 #0.077 SVOSOMG_01_01 + #0.795 #0.077 SVOSOMG_02_01 + #0.795 &
          #0.077 SVOSOMG_03_01 + #0.795 #0.077 SVOSOMG_04_01 + #0.795 #0.077 &
          SVOSOMG_05_01 + #0.795 #0.077 SVOSOMG_06_01 + #0.795 #0.077 SVOSOMG_07_01 + &
          #0.795 #0.077 SVOSOMG_08_01 + #0.795 #0.077 SVOSOMG_09_01 + #0.795 &
          #0.077 SVOSOMG_10_01 + #0.795 #0.077 SVOSOMG_11_01 + #0.795 #0.077 &
          SVOSOMG_12_01 + #0.795 #0.077 SVOSOMG_13_01 + #0.795 #0.077 SVOSOMG_01_02 + &
          #0.795 #0.077 SVOSOMG_02_02 + #0.795 #0.077 SVOSOMG_03_02 + #0.795 &
          #0.077 SVOSOMG_04_02 + #0.795 #0.077 SVOSOMG_05_02 + #0.795 #0.077 &
          SVOSOMG_06_02 + #0.795 #0.077 SVOSOMG_07_02 + #0.795 #0.077 SVOSOMG_08_02 + &
          #0.795 #0.077 SVOSOMG_09_02 + #0.795 #0.077 SVOSOMG_10_02 + #0.795 &
          #0.077 SVOSOMG_11_02 + #0.795 #0.077 SVOSOMG_12_02 + #0.795 #0.077 &
          SVOSOMG_13_02
 S427.7 SVOSOMG_14_02 + OH = #0.156 #0.361 SVOSOMG_14_03 + #0.156 #0.007 &
          SVOSOMG_14_04 + #0.156 #0.606 SVOSOMG_14_05 + #0.156 #0.026 SVOSOMG_14_06 + &
          #0.844 #0.054 SVOSOMG_01_01 + #0.844 #0.054 SVOSOMG_02_01 + #0.844 &
          #0.054 SVOSOMG_03_01 + #0.844 #0.054 SVOSOMG_04_01 + #0.844 #0.054 &
          SVOSOMG_05_01 + #0.844 #0.054 SVOSOMG_06_01 + #0.844 #0.054 SVOSOMG_07_01 + &
          #0.844 #0.054 SVOSOMG_08_01 + #0.844 #0.054 SVOSOMG_09_01 + #0.844 &
          #0.054 SVOSOMG_10_01 + #0.844 #0.054 SVOSOMG_11_01 + #0.844 #0.054 &
          SVOSOMG_12_01 + #0.844 #0.000 SVOSOMG_13_01 + #0.844 #0.054 SVOSOMG_01_02 + &
          #0.844 #0.054 SVOSOMG_02_02 + #0.844 #0.054 SVOSOMG_03_02 + #0.844 &
          #0.054 SVOSOMG_04_02 + #0.844 #0.054 SVOSOMG_05_02 + #0.844 #0.054 &
          SVOSOMG_06_02 + #0.844 #0.054 SVOSOMG_07_02 + #0.844 #0.054 SVOSOMG_08_02 + &
          #0.844 #0.054 SVOSOMG_09_02 + #0.844 #0.054 SVOSOMG_10_02 + #0.844 &
          #0.054 SVOSOMG_11_02 + #0.844 #0.054 SVOSOMG_12_02 + #0.844 #0.054 &
          SVOSOMG_13_02 + #0.844 #0.054 SVOSOMG_02_03 + #0.844 #0.054 SVOSOMG_03_03 + &
          #0.844 #0.054 SVOSOMG_04_03 + #0.844 #0.054 SVOSOMG_05_03 + #0.844 &
          #0.054 SVOSOMG_06_03 + #0.844 #0.054 SVOSOMG_07_03 + #0.844 #0.054 &
          SVOSOMG_08_03 + #0.844 #0.054 SVOSOMG_09_03 + #0.844 #0.054 SVOSOMG_10_03 + &
          #0.844 #0.054 SVOSOMG_11_03 + #0.844 #0.054 SVOSOMG_12_03 + #0.844 &
          #0.054 SVOSOMG_13_03
 S428.7 SVOSOMG_14_03 + OH = #0.125 #0.361 SVOSOMG_14_04 + #0.125 #0.007 &
          SVOSOMG_14_05 + #0.125 #0.606 SVOSOMG_14_06 + #0.125 #0.026 SVOSOMG_14_07 + &
          #0.875 #0.042 SVOSOMG_01_01 + #0.875 #0.042 SVOSOMG_02_01 + #0.875 &
          #0.042 SVOSOMG_03_01 + #0.875 #0.042 SVOSOMG_04_01 + #0.875 #0.042 &
          SVOSOMG_05_01 + #0.875 #0.042 SVOSOMG_06_01 + #0.875 #0.042 SVOSOMG_07_01 + &
          #0.875 #0.042 SVOSOMG_08_01 + #0.875 #0.042 SVOSOMG_09_01 + #0.875 &
          #0.042 SVOSOMG_10_01 + #0.875 #0.042 SVOSOMG_11_01 + #0.875 #0.042 &
          SVOSOMG_12_01 + #0.875 #0.000 SVOSOMG_13_01 + #0.875 #0.042 SVOSOMG_01_02 + &
          #0.875 #0.042 SVOSOMG_02_02 + #0.875 #0.042 SVOSOMG_03_02 + #0.875 &
          #0.042 SVOSOMG_04_02 + #0.875 #0.042 SVOSOMG_05_02 + #0.875 #0.042 &
          SVOSOMG_06_02 + #0.875 #0.042 SVOSOMG_07_02 + #0.875 #0.042 SVOSOMG_08_02 + &
          #0.875 #0.042 SVOSOMG_09_02 + #0.875 #0.042 SVOSOMG_10_02 + #0.875 &
          #0.042 SVOSOMG_11_02 + #0.875 #0.042 SVOSOMG_12_02 + #0.875 #0.000 &
          SVOSOMG_13_02 + #0.875 #0.042 SVOSOMG_02_03 + #0.875 #0.042 SVOSOMG_03_03 + &
          #0.875 #0.042 SVOSOMG_04_03 + #0.875 #0.042 SVOSOMG_05_03 + #0.875 &
          #0.042 SVOSOMG_06_03 + #0.875 #0.042 SVOSOMG_07_03 + #0.875 #0.042 &
          SVOSOMG_08_03 + #0.875 #0.042 SVOSOMG_09_03 + #0.875 #0.042 SVOSOMG_10_03 + &
          #0.875 #0.042 SVOSOMG_11_03 + #0.875 #0.042 SVOSOMG_12_03 + #0.875 &
          #0.042 SVOSOMG_13_03 + #0.875 #0.042 SVOSOMG_02_04 + #0.875 #0.042 &
          SVOSOMG_03_04 + #0.875 #0.042 SVOSOMG_04_04 + #0.875 #0.042 SVOSOMG_05_04 + &
          #0.875 #0.042 SVOSOMG_06_04 + #0.875 #0.042 SVOSOMG_07_04 + #0.875 &
          #0.042 SVOSOMG_08_04 + #0.875 #0.042 SVOSOMG_09_04 + #0.875 #0.042 &
          SVOSOMG_10_04 + #0.875 #0.042 SVOSOMG_11_04 + #0.875 #0.042 SVOSOMG_12_04 + &
          #0.875 #0.042 SVOSOMG_13_04
 S429.7 SVOSOMG_14_04 + OH = #0.103 #0.361 SVOSOMG_14_05 + #0.103 #0.007 &
          SVOSOMG_14_06 + #0.103 #0.632 SVOSOMG_14_07 + #0.897 #0.034 SVOSOMG_01_01 + &
          #0.897 #0.034 SVOSOMG_02_01 + #0.897 #0.034 SVOSOMG_03_01 + #0.897 &
          #0.034 SVOSOMG_04_01 + #0.897 #0.034 SVOSOMG_05_01 + #0.897 #0.034 &
          SVOSOMG_06_01 + #0.897 #0.034 SVOSOMG_07_01 + #0.897 #0.034 SVOSOMG_08_01 + &
          #0.897 #0.034 SVOSOMG_09_01 + #0.897 #0.034 SVOSOMG_10_01 + #0.897 &
          #0.034 SVOSOMG_11_01 + #0.897 #0.000 SVOSOMG_12_01 + #0.897 #0.000 &
          SVOSOMG_13_01 + #0.897 #0.034 SVOSOMG_01_02 + #0.897 #0.034 SVOSOMG_02_02 + &
          #0.897 #0.034 SVOSOMG_03_02 + #0.897 #0.034 SVOSOMG_04_02 + #0.897 &
          #0.034 SVOSOMG_05_02 + #0.897 #0.034 SVOSOMG_06_02 + #0.897 #0.034 &
          SVOSOMG_07_02 + #0.897 #0.034 SVOSOMG_08_02 + #0.897 #0.034 SVOSOMG_09_02 + &
          #0.897 #0.034 SVOSOMG_10_02 + #0.897 #0.034 SVOSOMG_11_02 + #0.897 &
          #0.034 SVOSOMG_12_02 + #0.897 #0.000 SVOSOMG_13_02 + #0.897 #0.034 &
          SVOSOMG_02_03 + #0.897 #0.034 SVOSOMG_03_03 + #0.897 #0.034 SVOSOMG_04_03 + &
          #0.897 #0.034 SVOSOMG_05_03 + #0.897 #0.034 SVOSOMG_06_03 + #0.897 &
          #0.034 SVOSOMG_07_03 + #0.897 #0.034 SVOSOMG_08_03 + #0.897 #0.034 &
          SVOSOMG_09_03 + #0.897 #0.034 SVOSOMG_10_03 + #0.897 #0.034 SVOSOMG_11_03 + &
          #0.897 #0.034 SVOSOMG_12_03 + #0.897 #0.000 SVOSOMG_13_03 + #0.897 &
          #0.034 SVOSOMG_02_04 + #0.897 #0.034 SVOSOMG_03_04 + #0.897 #0.034 &
          SVOSOMG_04_04 + #0.897 #0.034 SVOSOMG_05_04 + #0.897 #0.034 SVOSOMG_06_04 + &
          #0.897 #0.034 SVOSOMG_07_04 + #0.897 #0.034 SVOSOMG_08_04 + #0.897 &
          #0.034 SVOSOMG_09_04 + #0.897 #0.034 SVOSOMG_10_04 + #0.897 #0.034 &
          SVOSOMG_11_04 + #0.897 #0.034 SVOSOMG_12_04 + #0.897 #0.034 SVOSOMG_13_04 + &
          #0.897 #0.034 SVOSOMG_03_05 + #0.897 #0.034 SVOSOMG_04_05 + #0.897 &
          #0.034 SVOSOMG_05_05 + #0.897 #0.034 SVOSOMG_06_05 + #0.897 #0.034 &
          SVOSOMG_07_05 + #0.897 #0.034 SVOSOMG_08_05 + #0.897 #0.034 SVOSOMG_09_05 + &
          #0.897 #0.034 SVOSOMG_10_05 + #0.897 #0.034 SVOSOMG_11_05 + #0.897 &
          #0.034 SVOSOMG_12_05 + #0.897 #0.034 SVOSOMG_13_05
 S430.7 SVOSOMG_14_05 + OH = #0.086 #0.361 SVOSOMG_14_06 + #0.086 #0.639 &
          SVOSOMG_14_07 + #0.914 #0.030 SVOSOMG_01_01 + #0.914 #0.030 SVOSOMG_02_01 + &
          #0.914 #0.030 SVOSOMG_03_01 + #0.914 #0.030 SVOSOMG_04_01 + #0.914 &
          #0.030 SVOSOMG_05_01 + #0.914 #0.030 SVOSOMG_06_01 + #0.914 #0.030 &
          SVOSOMG_07_01 + #0.914 #0.030 SVOSOMG_08_01 + #0.914 #0.030 SVOSOMG_09_01 + &
          #0.914 #0.030 SVOSOMG_10_01 + #0.914 #0.030 SVOSOMG_11_01 + #0.914 &
          #0.000 SVOSOMG_12_01 + #0.914 #0.000 SVOSOMG_13_01 + #0.914 #0.030 &
          SVOSOMG_01_02 + #0.914 #0.030 SVOSOMG_02_02 + #0.914 #0.030 SVOSOMG_03_02 + &
          #0.914 #0.030 SVOSOMG_04_02 + #0.914 #0.030 SVOSOMG_05_02 + #0.914 &
          #0.030 SVOSOMG_06_02 + #0.914 #0.030 SVOSOMG_07_02 + #0.914 #0.030 &
          SVOSOMG_08_02 + #0.914 #0.030 SVOSOMG_09_02 + #0.914 #0.030 SVOSOMG_10_02 + &
          #0.914 #0.030 SVOSOMG_11_02 + #0.914 #0.000 SVOSOMG_12_02 + #0.914 &
          #0.000 SVOSOMG_13_02 + #0.914 #0.030 SVOSOMG_02_03 + #0.914 #0.030 &
          SVOSOMG_03_03 + #0.914 #0.030 SVOSOMG_04_03 + #0.914 #0.030 SVOSOMG_05_03 + &
          #0.914 #0.030 SVOSOMG_06_03 + #0.914 #0.030 SVOSOMG_07_03 + #0.914 &
          #0.030 SVOSOMG_08_03 + #0.914 #0.030 SVOSOMG_09_03 + #0.914 #0.030 &
          SVOSOMG_10_03 + #0.914 #0.030 SVOSOMG_11_03 + #0.914 #0.030 SVOSOMG_12_03 + &
          #0.914 #0.000 SVOSOMG_13_03 + #0.914 #0.030 SVOSOMG_02_04 + #0.914 &
          #0.030 SVOSOMG_03_04 + #0.914 #0.030 SVOSOMG_04_04 + #0.914 #0.030 &
          SVOSOMG_05_04 + #0.914 #0.030 SVOSOMG_06_04 + #0.914 #0.030 SVOSOMG_07_04 + &
          #0.914 #0.030 SVOSOMG_08_04 + #0.914 #0.030 SVOSOMG_09_04 + #0.914 &
          #0.030 SVOSOMG_10_04 + #0.914 #0.030 SVOSOMG_11_04 + #0.914 #0.030 &
          SVOSOMG_12_04 + #0.914 #0.000 SVOSOMG_13_04 + #0.914 #0.030 SVOSOMG_03_05 + &
          #0.914 #0.030 SVOSOMG_04_05 + #0.914 #0.030 SVOSOMG_05_05 + #0.914 &
          #0.030 SVOSOMG_06_05 + #0.914 #0.030 SVOSOMG_07_05 + #0.914 #0.030 &
          SVOSOMG_08_05 + #0.914 #0.030 SVOSOMG_09_05 + #0.914 #0.030 SVOSOMG_10_05 + &
          #0.914 #0.030 SVOSOMG_11_05 + #0.914 #0.030 SVOSOMG_12_05 + #0.914 &
          #0.030 SVOSOMG_13_05 + #0.914 #0.030 SVOSOMG_03_06 + #0.914 #0.030 &
          SVOSOMG_04_06 + #0.914 #0.030 SVOSOMG_05_06 + #0.914 #0.030 SVOSOMG_06_06 + &
          #0.914 #0.030 SVOSOMG_07_06 + #0.914 #0.030 SVOSOMG_08_06 + #0.914 &
          #0.030 SVOSOMG_09_06 + #0.914 #0.030 SVOSOMG_10_06 + #0.914 #0.030 &
          SVOSOMG_11_06 + #0.914 #0.030 SVOSOMG_12_06 + #0.914 #0.030 SVOSOMG_13_06
 S431.7 SVOSOMG_14_06 + OH = #0.071 #1.000 SVOSOMG_14_07 + #0.929 #0.027 &
          SVOSOMG_01_01 + #0.929 #0.027 SVOSOMG_02_01 + #0.929 #0.027 SVOSOMG_03_01 + &
          #0.929 #0.027 SVOSOMG_04_01 + #0.929 #0.027 SVOSOMG_05_01 + #0.929 &
          #0.027 SVOSOMG_06_01 + #0.929 #0.027 SVOSOMG_07_01 + #0.929 #0.027 &
          SVOSOMG_08_01 + #0.929 #0.027 SVOSOMG_09_01 + #0.929 #0.027 SVOSOMG_10_01 + &
          #0.929 #0.000 SVOSOMG_11_01 + #0.929 #0.000 SVOSOMG_12_01 + #0.929 &
          #0.000 SVOSOMG_13_01 + #0.929 #0.027 SVOSOMG_01_02 + #0.929 #0.027 &
          SVOSOMG_02_02 + #0.929 #0.027 SVOSOMG_03_02 + #0.929 #0.027 SVOSOMG_04_02 + &
          #0.929 #0.027 SVOSOMG_05_02 + #0.929 #0.027 SVOSOMG_06_02 + #0.929 &
          #0.027 SVOSOMG_07_02 + #0.929 #0.027 SVOSOMG_08_02 + #0.929 #0.027 &
          SVOSOMG_09_02 + #0.929 #0.027 SVOSOMG_10_02 + #0.929 #0.027 SVOSOMG_11_02 + &
          #0.929 #0.000 SVOSOMG_12_02 + #0.929 #0.000 SVOSOMG_13_02 + #0.929 &
          #0.027 SVOSOMG_02_03 + #0.929 #0.027 SVOSOMG_03_03 + #0.929 #0.027 &
          SVOSOMG_04_03 + #0.929 #0.027 SVOSOMG_05_03 + #0.929 #0.027 SVOSOMG_06_03 + &
          #0.929 #0.027 SVOSOMG_07_03 + #0.929 #0.027 SVOSOMG_08_03 + #0.929 &
          #0.027 SVOSOMG_09_03 + #0.929 #0.027 SVOSOMG_10_03 + #0.929 #0.027 &
          SVOSOMG_11_03 + #0.929 #0.000 SVOSOMG_12_03 + #0.929 #0.000 SVOSOMG_13_03 + &
          #0.929 #0.027 SVOSOMG_02_04 + #0.929 #0.027 SVOSOMG_03_04 + #0.929 &
          #0.027 SVOSOMG_04_04 + #0.929 #0.027 SVOSOMG_05_04 + #0.929 #0.027 &
          SVOSOMG_06_04 + #0.929 #0.027 SVOSOMG_07_04 + #0.929 #0.027 SVOSOMG_08_04 + &
          #0.929 #0.027 SVOSOMG_09_04 + #0.929 #0.027 SVOSOMG_10_04 + #0.929 &
          #0.027 SVOSOMG_11_04 + #0.929 #0.027 SVOSOMG_12_04 + #0.929 #0.000 &
          SVOSOMG_13_04 + #0.929 #0.027 SVOSOMG_03_05 + #0.929 #0.027 SVOSOMG_04_05 + &
          #0.929 #0.027 SVOSOMG_05_05 + #0.929 #0.027 SVOSOMG_06_05 + #0.929 &
          #0.027 SVOSOMG_07_05 + #0.929 #0.027 SVOSOMG_08_05 + #0.929 #0.027 &
          SVOSOMG_09_05 + #0.929 #0.027 SVOSOMG_10_05 + #0.929 #0.027 SVOSOMG_11_05 + &
          #0.929 #0.027 SVOSOMG_12_05 + #0.929 #0.000 SVOSOMG_13_05 + #0.929 &
          #0.027 SVOSOMG_03_06 + #0.929 #0.027 SVOSOMG_04_06 + #0.929 #0.027 &
          SVOSOMG_05_06 + #0.929 #0.027 SVOSOMG_06_06 + #0.929 #0.027 SVOSOMG_07_06 + &
          #0.929 #0.027 SVOSOMG_08_06 + #0.929 #0.027 SVOSOMG_09_06 + #0.929 &
          #0.027 SVOSOMG_10_06 + #0.929 #0.027 SVOSOMG_11_06 + #0.929 #0.027 &
          SVOSOMG_12_06 + #0.929 #0.027 SVOSOMG_13_06 + #0.929 #0.027 SVOSOMG_04_07 + &
          #0.929 #0.027 SVOSOMG_05_07 + #0.929 #0.027 SVOSOMG_06_07 + #0.929 &
          #0.027 SVOSOMG_07_07 + #0.929 #0.027 SVOSOMG_08_07 + #0.929 #0.027 &
          SVOSOMG_09_07 + #0.929 #0.027 SVOSOMG_10_07 + #0.929 #0.027 SVOSOMG_11_07 + &
          #0.929 #0.027 SVOSOMG_12_07 + #0.929 #0.027 SVOSOMG_13_07
 S432.7 SVOSOMG_14_07 + OH = #1.000 #0.028 SVOSOMG_01_01 + #1.000 #0.028 &
          SVOSOMG_02_01 + #1.000 #0.028 SVOSOMG_03_01 + #1.000 #0.028 SVOSOMG_04_01 + &
          #1.000 #0.028 SVOSOMG_05_01 + #1.000 #0.028 SVOSOMG_06_01 + #1.000 &
          #0.028 SVOSOMG_07_01 + #1.000 #0.028 SVOSOMG_08_01 + #1.000 #0.028 &
          SVOSOMG_09_01 + #1.000 #0.028 SVOSOMG_10_01 + #1.000 #0.000 SVOSOMG_11_01 + &
          #1.000 #0.000 SVOSOMG_12_01 + #1.000 #0.000 SVOSOMG_13_01 + #1.000 &
          #0.028 SVOSOMG_01_02 + #1.000 #0.028 SVOSOMG_02_02 + #1.000 #0.028 &
          SVOSOMG_03_02 + #1.000 #0.028 SVOSOMG_04_02 + #1.000 #0.028 SVOSOMG_05_02 + &
          #1.000 #0.028 SVOSOMG_06_02 + #1.000 #0.028 SVOSOMG_07_02 + #1.000 &
          #0.028 SVOSOMG_08_02 + #1.000 #0.028 SVOSOMG_09_02 + #1.000 #0.028 &
          SVOSOMG_10_02 + #1.000 #0.000 SVOSOMG_11_02 + #1.000 #0.000 SVOSOMG_12_02 + &
          #1.000 #0.000 SVOSOMG_13_02 + #1.000 #0.028 SVOSOMG_02_03 + #1.000 &
          #0.028 SVOSOMG_03_03 + #1.000 #0.028 SVOSOMG_04_03 + #1.000 #0.028 &
          SVOSOMG_05_03 + #1.000 #0.028 SVOSOMG_06_03 + #1.000 #0.028 SVOSOMG_07_03 + &
          #1.000 #0.028 SVOSOMG_08_03 + #1.000 #0.028 SVOSOMG_09_03 + #1.000 &
          #0.028 SVOSOMG_10_03 + #1.000 #0.028 SVOSOMG_11_03 + #1.000 #0.000 &
          SVOSOMG_12_03 + #1.000 #0.000 SVOSOMG_13_03 + #1.000 #0.028 SVOSOMG_02_04 + &
          #1.000 #0.028 SVOSOMG_03_04 + #1.000 #0.028 SVOSOMG_04_04 + #1.000 &
          #0.028 SVOSOMG_05_04 + #1.000 #0.028 SVOSOMG_06_04 + #1.000 #0.028 &
          SVOSOMG_07_04 + #1.000 #0.028 SVOSOMG_08_04 + #1.000 #0.028 SVOSOMG_09_04 + &
          #1.000 #0.028 SVOSOMG_10_04 + #1.000 #0.028 SVOSOMG_11_04 + #1.000 &
          #0.000 SVOSOMG_12_04 + #1.000 #0.000 SVOSOMG_13_04 + #1.000 #0.028 &
          SVOSOMG_03_05 + #1.000 #0.028 SVOSOMG_04_05 + #1.000 #0.028 SVOSOMG_05_05 + &
          #1.000 #0.028 SVOSOMG_06_05 + #1.000 #0.028 SVOSOMG_07_05 + #1.000 &
          #0.028 SVOSOMG_08_05 + #1.000 #0.028 SVOSOMG_09_05 + #1.000 #0.028 &
          SVOSOMG_10_05 + #1.000 #0.028 SVOSOMG_11_05 + #1.000 #0.028 SVOSOMG_12_05 + &
          #1.000 #0.000 SVOSOMG_13_05 + #1.000 #0.028 SVOSOMG_03_06 + #1.000 &
          #0.028 SVOSOMG_04_06 + #1.000 #0.028 SVOSOMG_05_06 + #1.000 #0.028 &
          SVOSOMG_06_06 + #1.000 #0.028 SVOSOMG_07_06 + #1.000 #0.028 SVOSOMG_08_06 + &
          #1.000 #0.028 SVOSOMG_09_06 + #1.000 #0.028 SVOSOMG_10_06 + #1.000 &
          #0.028 SVOSOMG_11_06 + #1.000 #0.028 SVOSOMG_12_06 + #1.000 #0.000 &
          SVOSOMG_13_06 + #1.000 #0.028 SVOSOMG_04_07 + #1.000 #0.028 SVOSOMG_05_07 + &
          #1.000 #0.028 SVOSOMG_06_07 + #1.000 #0.028 SVOSOMG_07_07 + #1.000 &
          #0.028 SVOSOMG_08_07 + #1.000 #0.028 SVOSOMG_09_07 + #1.000 #0.028 &
          SVOSOMG_10_07 + #1.000 #0.028 SVOSOMG_11_07 + #1.000 #0.028 SVOSOMG_12_07 + &
          #1.000 #0.028 SVOSOMG_13_07
 S433.7 SVOSOMG_15_01 + OH = #0.210 #0.361 SVOSOMG_15_02 + #0.210 #0.007 &
          SVOSOMG_15_03 + #0.210 #0.606 SVOSOMG_15_04 + #0.210 #0.026 SVOSOMG_15_05 + &
          #0.790 #0.071 SVOSOMG_01_01 + #0.790 #0.071 SVOSOMG_02_01 + #0.790 &
          #0.071 SVOSOMG_03_01 + #0.790 #0.071 SVOSOMG_04_01 + #0.790 #0.071 &
          SVOSOMG_05_01 + #0.790 #0.071 SVOSOMG_06_01 + #0.790 #0.071 SVOSOMG_07_01 + &
          #0.790 #0.071 SVOSOMG_08_01 + #0.790 #0.071 SVOSOMG_09_01 + #0.790 &
          #0.071 SVOSOMG_10_01 + #0.790 #0.071 SVOSOMG_11_01 + #0.790 #0.071 &
          SVOSOMG_12_01 + #0.790 #0.071 SVOSOMG_13_01 + #0.790 #0.071 SVOSOMG_14_01 + &
          #0.790 #0.071 SVOSOMG_01_02 + #0.790 #0.071 SVOSOMG_02_02 + #0.790 &
          #0.071 SVOSOMG_03_02 + #0.790 #0.071 SVOSOMG_04_02 + #0.790 #0.071 &
          SVOSOMG_05_02 + #0.790 #0.071 SVOSOMG_06_02 + #0.790 #0.071 SVOSOMG_07_02 + &
          #0.790 #0.071 SVOSOMG_08_02 + #0.790 #0.071 SVOSOMG_09_02 + #0.790 &
          #0.071 SVOSOMG_10_02 + #0.790 #0.071 SVOSOMG_11_02 + #0.790 #0.071 &
          SVOSOMG_12_02 + #0.790 #0.071 SVOSOMG_13_02 + #0.790 #0.071 SVOSOMG_14_02
 S434.7 SVOSOMG_15_02 + OH = #0.161 #0.361 SVOSOMG_15_03 + #0.161 #0.007 &
          SVOSOMG_15_04 + #0.161 #0.606 SVOSOMG_15_05 + #0.161 #0.026 SVOSOMG_15_06 + &
          #0.839 #0.050 SVOSOMG_01_01 + #0.839 #0.050 SVOSOMG_02_01 + #0.839 &
          #0.050 SVOSOMG_03_01 + #0.839 #0.050 SVOSOMG_04_01 + #0.839 #0.050 &
          SVOSOMG_05_01 + #0.839 #0.050 SVOSOMG_06_01 + #0.839 #0.050 SVOSOMG_07_01 + &
          #0.839 #0.050 SVOSOMG_08_01 + #0.839 #0.050 SVOSOMG_09_01 + #0.839 &
          #0.050 SVOSOMG_10_01 + #0.839 #0.050 SVOSOMG_11_01 + #0.839 #0.050 &
          SVOSOMG_12_01 + #0.839 #0.050 SVOSOMG_13_01 + #0.839 #0.000 SVOSOMG_14_01 + &
          #0.839 #0.050 SVOSOMG_01_02 + #0.839 #0.050 SVOSOMG_02_02 + #0.839 &
          #0.050 SVOSOMG_03_02 + #0.839 #0.050 SVOSOMG_04_02 + #0.839 #0.050 &
          SVOSOMG_05_02 + #0.839 #0.050 SVOSOMG_06_02 + #0.839 #0.050 SVOSOMG_07_02 + &
          #0.839 #0.050 SVOSOMG_08_02 + #0.839 #0.050 SVOSOMG_09_02 + #0.839 &
          #0.050 SVOSOMG_10_02 + #0.839 #0.050 SVOSOMG_11_02 + #0.839 #0.050 &
          SVOSOMG_12_02 + #0.839 #0.050 SVOSOMG_13_02 + #0.839 #0.050 SVOSOMG_14_02 + &
          #0.839 #0.050 SVOSOMG_02_03 + #0.839 #0.050 SVOSOMG_03_03 + #0.839 &
          #0.050 SVOSOMG_04_03 + #0.839 #0.050 SVOSOMG_05_03 + #0.839 #0.050 &
          SVOSOMG_06_03 + #0.839 #0.050 SVOSOMG_07_03 + #0.839 #0.050 SVOSOMG_08_03 + &
          #0.839 #0.050 SVOSOMG_09_03 + #0.839 #0.050 SVOSOMG_10_03 + #0.839 &
          #0.050 SVOSOMG_11_03 + #0.839 #0.050 SVOSOMG_12_03 + #0.839 #0.050 &
          SVOSOMG_13_03 + #0.839 #0.050 SVOSOMG_14_03
 S435.7 SVOSOMG_15_03 + OH = #0.131 #0.361 SVOSOMG_15_04 + #0.131 #0.007 &
          SVOSOMG_15_05 + #0.131 #0.606 SVOSOMG_15_06 + #0.131 #0.026 SVOSOMG_15_07 + &
          #0.869 #0.038 SVOSOMG_01_01 + #0.869 #0.038 SVOSOMG_02_01 + #0.869 &
          #0.038 SVOSOMG_03_01 + #0.869 #0.038 SVOSOMG_04_01 + #0.869 #0.038 &
          SVOSOMG_05_01 + #0.869 #0.038 SVOSOMG_06_01 + #0.869 #0.038 SVOSOMG_07_01 + &
          #0.869 #0.038 SVOSOMG_08_01 + #0.869 #0.038 SVOSOMG_09_01 + #0.869 &
          #0.038 SVOSOMG_10_01 + #0.869 #0.038 SVOSOMG_11_01 + #0.869 #0.038 &
          SVOSOMG_12_01 + #0.869 #0.038 SVOSOMG_13_01 + #0.869 #0.000 SVOSOMG_14_01 + &
          #0.869 #0.038 SVOSOMG_01_02 + #0.869 #0.038 SVOSOMG_02_02 + #0.869 &
          #0.038 SVOSOMG_03_02 + #0.869 #0.038 SVOSOMG_04_02 + #0.869 #0.038 &
          SVOSOMG_05_02 + #0.869 #0.038 SVOSOMG_06_02 + #0.869 #0.038 SVOSOMG_07_02 + &
          #0.869 #0.038 SVOSOMG_08_02 + #0.869 #0.038 SVOSOMG_09_02 + #0.869 &
          #0.038 SVOSOMG_10_02 + #0.869 #0.038 SVOSOMG_11_02 + #0.869 #0.038 &
          SVOSOMG_12_02 + #0.869 #0.038 SVOSOMG_13_02 + #0.869 #0.000 SVOSOMG_14_02 + &
          #0.869 #0.038 SVOSOMG_02_03 + #0.869 #0.038 SVOSOMG_03_03 + #0.869 &
          #0.038 SVOSOMG_04_03 + #0.869 #0.038 SVOSOMG_05_03 + #0.869 #0.038 &
          SVOSOMG_06_03 + #0.869 #0.038 SVOSOMG_07_03 + #0.869 #0.038 SVOSOMG_08_03 + &
          #0.869 #0.038 SVOSOMG_09_03 + #0.869 #0.038 SVOSOMG_10_03 + #0.869 &
          #0.038 SVOSOMG_11_03 + #0.869 #0.038 SVOSOMG_12_03 + #0.869 #0.038 &
          SVOSOMG_13_03 + #0.869 #0.038 SVOSOMG_14_03 + #0.869 #0.038 SVOSOMG_02_04 + &
          #0.869 #0.038 SVOSOMG_03_04 + #0.869 #0.038 SVOSOMG_04_04 + #0.869 &
          #0.038 SVOSOMG_05_04 + #0.869 #0.038 SVOSOMG_06_04 + #0.869 #0.038 &
          SVOSOMG_07_04 + #0.869 #0.038 SVOSOMG_08_04 + #0.869 #0.038 SVOSOMG_09_04 + &
          #0.869 #0.038 SVOSOMG_10_04 + #0.869 #0.038 SVOSOMG_11_04 + #0.869 &
          #0.038 SVOSOMG_12_04 + #0.869 #0.038 SVOSOMG_13_04 + #0.869 #0.038 &
          SVOSOMG_14_04
 S436.7 SVOSOMG_15_04 + OH = #0.109 #0.361 SVOSOMG_15_05 + #0.109 #0.007 &
          SVOSOMG_15_06 + #0.109 #0.632 SVOSOMG_15_07 + #0.891 #0.032 SVOSOMG_01_01 + &
          #0.891 #0.032 SVOSOMG_02_01 + #0.891 #0.032 SVOSOMG_03_01 + #0.891 &
          #0.032 SVOSOMG_04_01 + #0.891 #0.032 SVOSOMG_05_01 + #0.891 #0.032 &
          SVOSOMG_06_01 + #0.891 #0.032 SVOSOMG_07_01 + #0.891 #0.032 SVOSOMG_08_01 + &
          #0.891 #0.032 SVOSOMG_09_01 + #0.891 #0.032 SVOSOMG_10_01 + #0.891 &
          #0.032 SVOSOMG_11_01 + #0.891 #0.032 SVOSOMG_12_01 + #0.891 #0.000 &
          SVOSOMG_13_01 + #0.891 #0.000 SVOSOMG_14_01 + #0.891 #0.032 SVOSOMG_01_02 + &
          #0.891 #0.032 SVOSOMG_02_02 + #0.891 #0.032 SVOSOMG_03_02 + #0.891 &
          #0.032 SVOSOMG_04_02 + #0.891 #0.032 SVOSOMG_05_02 + #0.891 #0.032 &
          SVOSOMG_06_02 + #0.891 #0.032 SVOSOMG_07_02 + #0.891 #0.032 SVOSOMG_08_02 + &
          #0.891 #0.032 SVOSOMG_09_02 + #0.891 #0.032 SVOSOMG_10_02 + #0.891 &
          #0.032 SVOSOMG_11_02 + #0.891 #0.032 SVOSOMG_12_02 + #0.891 #0.032 &
          SVOSOMG_13_02 + #0.891 #0.000 SVOSOMG_14_02 + #0.891 #0.032 SVOSOMG_02_03 + &
          #0.891 #0.032 SVOSOMG_03_03 + #0.891 #0.032 SVOSOMG_04_03 + #0.891 &
          #0.032 SVOSOMG_05_03 + #0.891 #0.032 SVOSOMG_06_03 + #0.891 #0.032 &
          SVOSOMG_07_03 + #0.891 #0.032 SVOSOMG_08_03 + #0.891 #0.032 SVOSOMG_09_03 + &
          #0.891 #0.032 SVOSOMG_10_03 + #0.891 #0.032 SVOSOMG_11_03 + #0.891 &
          #0.032 SVOSOMG_12_03 + #0.891 #0.032 SVOSOMG_13_03 + #0.891 #0.000 &
          SVOSOMG_14_03 + #0.891 #0.032 SVOSOMG_02_04 + #0.891 #0.032 SVOSOMG_03_04 + &
          #0.891 #0.032 SVOSOMG_04_04 + #0.891 #0.032 SVOSOMG_05_04 + #0.891 &
          #0.032 SVOSOMG_06_04 + #0.891 #0.032 SVOSOMG_07_04 + #0.891 #0.032 &
          SVOSOMG_08_04 + #0.891 #0.032 SVOSOMG_09_04 + #0.891 #0.032 SVOSOMG_10_04 + &
          #0.891 #0.032 SVOSOMG_11_04 + #0.891 #0.032 SVOSOMG_12_04 + #0.891 &
          #0.032 SVOSOMG_13_04 + #0.891 #0.032 SVOSOMG_14_04 + #0.891 #0.032 &
          SVOSOMG_03_05 + #0.891 #0.032 SVOSOMG_04_05 + #0.891 #0.032 SVOSOMG_05_05 + &
          #0.891 #0.032 SVOSOMG_06_05 + #0.891 #0.032 SVOSOMG_07_05 + #0.891 &
          #0.032 SVOSOMG_08_05 + #0.891 #0.032 SVOSOMG_09_05 + #0.891 #0.032 &
          SVOSOMG_10_05 + #0.891 #0.032 SVOSOMG_11_05 + #0.891 #0.032 SVOSOMG_12_05 + &
          #0.891 #0.032 SVOSOMG_13_05 + #0.891 #0.032 SVOSOMG_14_05
 S437.7 SVOSOMG_15_05 + OH = #0.091 #0.361 SVOSOMG_15_06 + #0.091 #0.639 &
          SVOSOMG_15_07 + #0.909 #0.028 SVOSOMG_01_01 + #0.909 #0.028 SVOSOMG_02_01 + &
          #0.909 #0.028 SVOSOMG_03_01 + #0.909 #0.028 SVOSOMG_04_01 + #0.909 &
          #0.028 SVOSOMG_05_01 + #0.909 #0.028 SVOSOMG_06_01 + #0.909 #0.028 &
          SVOSOMG_07_01 + #0.909 #0.028 SVOSOMG_08_01 + #0.909 #0.028 SVOSOMG_09_01 + &
          #0.909 #0.028 SVOSOMG_10_01 + #0.909 #0.028 SVOSOMG_11_01 + #0.909 &
          #0.028 SVOSOMG_12_01 + #0.909 #0.000 SVOSOMG_13_01 + #0.909 #0.000 &
          SVOSOMG_14_01 + #0.909 #0.028 SVOSOMG_01_02 + #0.909 #0.028 SVOSOMG_02_02 + &
          #0.909 #0.028 SVOSOMG_03_02 + #0.909 #0.028 SVOSOMG_04_02 + #0.909 &
          #0.028 SVOSOMG_05_02 + #0.909 #0.028 SVOSOMG_06_02 + #0.909 #0.028 &
          SVOSOMG_07_02 + #0.909 #0.028 SVOSOMG_08_02 + #0.909 #0.028 SVOSOMG_09_02 + &
          #0.909 #0.028 SVOSOMG_10_02 + #0.909 #0.028 SVOSOMG_11_02 + #0.909 &
          #0.028 SVOSOMG_12_02 + #0.909 #0.000 SVOSOMG_13_02 + #0.909 #0.000 &
          SVOSOMG_14_02 + #0.909 #0.028 SVOSOMG_02_03 + #0.909 #0.028 SVOSOMG_03_03 + &
          #0.909 #0.028 SVOSOMG_04_03 + #0.909 #0.028 SVOSOMG_05_03 + #0.909 &
          #0.028 SVOSOMG_06_03 + #0.909 #0.028 SVOSOMG_07_03 + #0.909 #0.028 &
          SVOSOMG_08_03 + #0.909 #0.028 SVOSOMG_09_03 + #0.909 #0.028 SVOSOMG_10_03 + &
          #0.909 #0.028 SVOSOMG_11_03 + #0.909 #0.028 SVOSOMG_12_03 + #0.909 &
          #0.028 SVOSOMG_13_03 + #0.909 #0.000 SVOSOMG_14_03 + #0.909 #0.028 &
          SVOSOMG_02_04 + #0.909 #0.028 SVOSOMG_03_04 + #0.909 #0.028 SVOSOMG_04_04 + &
          #0.909 #0.028 SVOSOMG_05_04 + #0.909 #0.028 SVOSOMG_06_04 + #0.909 &
          #0.028 SVOSOMG_07_04 + #0.909 #0.028 SVOSOMG_08_04 + #0.909 #0.028 &
          SVOSOMG_09_04 + #0.909 #0.028 SVOSOMG_10_04 + #0.909 #0.028 SVOSOMG_11_04 + &
          #0.909 #0.028 SVOSOMG_12_04 + #0.909 #0.028 SVOSOMG_13_04 + #0.909 &
          #0.000 SVOSOMG_14_04 + #0.909 #0.028 SVOSOMG_03_05 + #0.909 #0.028 &
          SVOSOMG_04_05 + #0.909 #0.028 SVOSOMG_05_05 + #0.909 #0.028 SVOSOMG_06_05 + &
          #0.909 #0.028 SVOSOMG_07_05 + #0.909 #0.028 SVOSOMG_08_05 + #0.909 &
          #0.028 SVOSOMG_09_05 + #0.909 #0.028 SVOSOMG_10_05 + #0.909 #0.028 &
          SVOSOMG_11_05 + #0.909 #0.028 SVOSOMG_12_05 + #0.909 #0.028 SVOSOMG_13_05 + &
          #0.909 #0.028 SVOSOMG_14_05 + #0.909 #0.028 SVOSOMG_03_06 + #0.909 &
          #0.028 SVOSOMG_04_06 + #0.909 #0.028 SVOSOMG_05_06 + #0.909 #0.028 &
          SVOSOMG_06_06 + #0.909 #0.028 SVOSOMG_07_06 + #0.909 #0.028 SVOSOMG_08_06 + &
          #0.909 #0.028 SVOSOMG_09_06 + #0.909 #0.028 SVOSOMG_10_06 + #0.909 &
          #0.028 SVOSOMG_11_06 + #0.909 #0.028 SVOSOMG_12_06 + #0.909 #0.028 &
          SVOSOMG_13_06 + #0.909 #0.028 SVOSOMG_14_06
 S438.7 SVOSOMG_15_06 + OH = #0.077 #1.000 SVOSOMG_15_07 + #0.923 #0.025 &
          SVOSOMG_01_01 + #0.923 #0.025 SVOSOMG_02_01 + #0.923 #0.025 SVOSOMG_03_01 + &
          #0.923 #0.025 SVOSOMG_04_01 + #0.923 #0.025 SVOSOMG_05_01 + #0.923 &
          #0.025 SVOSOMG_06_01 + #0.923 #0.025 SVOSOMG_07_01 + #0.923 #0.025 &
          SVOSOMG_08_01 + #0.923 #0.025 SVOSOMG_09_01 + #0.923 #0.025 SVOSOMG_10_01 + &
          #0.923 #0.025 SVOSOMG_11_01 + #0.923 #0.000 SVOSOMG_12_01 + #0.923 &
          #0.000 SVOSOMG_13_01 + #0.923 #0.000 SVOSOMG_14_01 + #0.923 #0.025 &
          SVOSOMG_01_02 + #0.923 #0.025 SVOSOMG_02_02 + #0.923 #0.025 SVOSOMG_03_02 + &
          #0.923 #0.025 SVOSOMG_04_02 + #0.923 #0.025 SVOSOMG_05_02 + #0.923 &
          #0.025 SVOSOMG_06_02 + #0.923 #0.025 SVOSOMG_07_02 + #0.923 #0.025 &
          SVOSOMG_08_02 + #0.923 #0.025 SVOSOMG_09_02 + #0.923 #0.025 SVOSOMG_10_02 + &
          #0.923 #0.025 SVOSOMG_11_02 + #0.923 #0.025 SVOSOMG_12_02 + #0.923 &
          #0.000 SVOSOMG_13_02 + #0.923 #0.000 SVOSOMG_14_02 + #0.923 #0.025 &
          SVOSOMG_02_03 + #0.923 #0.025 SVOSOMG_03_03 + #0.923 #0.025 SVOSOMG_04_03 + &
          #0.923 #0.025 SVOSOMG_05_03 + #0.923 #0.025 SVOSOMG_06_03 + #0.923 &
          #0.025 SVOSOMG_07_03 + #0.923 #0.025 SVOSOMG_08_03 + #0.923 #0.025 &
          SVOSOMG_09_03 + #0.923 #0.025 SVOSOMG_10_03 + #0.923 #0.025 SVOSOMG_11_03 + &
          #0.923 #0.025 SVOSOMG_12_03 + #0.923 #0.000 SVOSOMG_13_03 + #0.923 &
          #0.000 SVOSOMG_14_03 + #0.923 #0.025 SVOSOMG_02_04 + #0.923 #0.025 &
          SVOSOMG_03_04 + #0.923 #0.025 SVOSOMG_04_04 + #0.923 #0.025 SVOSOMG_05_04 + &
          #0.923 #0.025 SVOSOMG_06_04 + #0.923 #0.025 SVOSOMG_07_04 + #0.923 &
          #0.025 SVOSOMG_08_04 + #0.923 #0.025 SVOSOMG_09_04 + #0.923 #0.025 &
          SVOSOMG_10_04 + #0.923 #0.025 SVOSOMG_11_04 + #0.923 #0.025 SVOSOMG_12_04 + &
          #0.923 #0.025 SVOSOMG_13_04 + #0.923 #0.000 SVOSOMG_14_04 + #0.923 &
          #0.025 SVOSOMG_03_05 + #0.923 #0.025 SVOSOMG_04_05 + #0.923 #0.025 &
          SVOSOMG_05_05 + #0.923 #0.025 SVOSOMG_06_05 + #0.923 #0.025 SVOSOMG_07_05 + &
          #0.923 #0.025 SVOSOMG_08_05 + #0.923 #0.025 SVOSOMG_09_05 + #0.923 &
          #0.025 SVOSOMG_10_05 + #0.923 #0.025 SVOSOMG_11_05 + #0.923 #0.025 &
          SVOSOMG_12_05 + #0.923 #0.025 SVOSOMG_13_05 + #0.923 #0.000 SVOSOMG_14_05 + &
          #0.923 #0.025 SVOSOMG_03_06 + #0.923 #0.025 SVOSOMG_04_06 + #0.923 &
          #0.025 SVOSOMG_05_06 + #0.923 #0.025 SVOSOMG_06_06 + #0.923 #0.025 &
          SVOSOMG_07_06 + #0.923 #0.025 SVOSOMG_08_06 + #0.923 #0.025 SVOSOMG_09_06 + &
          #0.923 #0.025 SVOSOMG_10_06 + #0.923 #0.025 SVOSOMG_11_06 + #0.923 &
          #0.025 SVOSOMG_12_06 + #0.923 #0.025 SVOSOMG_13_06 + #0.923 #0.025 &
          SVOSOMG_14_06 + #0.923 #0.025 SVOSOMG_04_07 + #0.923 #0.025 SVOSOMG_05_07 + &
          #0.923 #0.025 SVOSOMG_06_07 + #0.923 #0.025 SVOSOMG_07_07 + #0.923 &
          #0.025 SVOSOMG_08_07 + #0.923 #0.025 SVOSOMG_09_07 + #0.923 #0.025 &
          SVOSOMG_10_07 + #0.923 #0.025 SVOSOMG_11_07 + #0.923 #0.025 SVOSOMG_12_07 + &
          #0.923 #0.025 SVOSOMG_13_07 + #0.923 #0.025 SVOSOMG_14_07
 S439.7 SVOSOMG_15_07 + OH = #1.000 #0.026 SVOSOMG_01_01 + #1.000 #0.026 &
          SVOSOMG_02_01 + #1.000 #0.026 SVOSOMG_03_01 + #1.000 #0.026 SVOSOMG_04_01 + &
          #1.000 #0.026 SVOSOMG_05_01 + #1.000 #0.026 SVOSOMG_06_01 + #1.000 &
          #0.026 SVOSOMG_07_01 + #1.000 #0.026 SVOSOMG_08_01 + #1.000 #0.026 &
          SVOSOMG_09_01 + #1.000 #0.026 SVOSOMG_10_01 + #1.000 #0.026 SVOSOMG_11_01 + &
          #1.000 #0.000 SVOSOMG_12_01 + #1.000 #0.000 SVOSOMG_13_01 + #1.000 &
          #0.000 SVOSOMG_14_01 + #1.000 #0.026 SVOSOMG_01_02 + #1.000 #0.026 &
          SVOSOMG_02_02 + #1.000 #0.026 SVOSOMG_03_02 + #1.000 #0.026 SVOSOMG_04_02 + &
          #1.000 #0.026 SVOSOMG_05_02 + #1.000 #0.026 SVOSOMG_06_02 + #1.000 &
          #0.026 SVOSOMG_07_02 + #1.000 #0.026 SVOSOMG_08_02 + #1.000 #0.026 &
          SVOSOMG_09_02 + #1.000 #0.026 SVOSOMG_10_02 + #1.000 #0.026 SVOSOMG_11_02 + &
          #1.000 #0.000 SVOSOMG_12_02 + #1.000 #0.000 SVOSOMG_13_02 + #1.000 &
          #0.000 SVOSOMG_14_02 + #1.000 #0.026 SVOSOMG_02_03 + #1.000 #0.026 &
          SVOSOMG_03_03 + #1.000 #0.026 SVOSOMG_04_03 + #1.000 #0.026 SVOSOMG_05_03 + &
          #1.000 #0.026 SVOSOMG_06_03 + #1.000 #0.026 SVOSOMG_07_03 + #1.000 &
          #0.026 SVOSOMG_08_03 + #1.000 #0.026 SVOSOMG_09_03 + #1.000 #0.026 &
          SVOSOMG_10_03 + #1.000 #0.026 SVOSOMG_11_03 + #1.000 #0.026 SVOSOMG_12_03 + &
          #1.000 #0.000 SVOSOMG_13_03 + #1.000 #0.000 SVOSOMG_14_03 + #1.000 &
          #0.026 SVOSOMG_02_04 + #1.000 #0.026 SVOSOMG_03_04 + #1.000 #0.026 &
          SVOSOMG_04_04 + #1.000 #0.026 SVOSOMG_05_04 + #1.000 #0.026 SVOSOMG_06_04 + &
          #1.000 #0.026 SVOSOMG_07_04 + #1.000 #0.026 SVOSOMG_08_04 + #1.000 &
          #0.026 SVOSOMG_09_04 + #1.000 #0.026 SVOSOMG_10_04 + #1.000 #0.026 &
          SVOSOMG_11_04 + #1.000 #0.026 SVOSOMG_12_04 + #1.000 #0.000 SVOSOMG_13_04 + &
          #1.000 #0.000 SVOSOMG_14_04 + #1.000 #0.026 SVOSOMG_03_05 + #1.000 &
          #0.026 SVOSOMG_04_05 + #1.000 #0.026 SVOSOMG_05_05 + #1.000 #0.026 &
          SVOSOMG_06_05 + #1.000 #0.026 SVOSOMG_07_05 + #1.000 #0.026 SVOSOMG_08_05 + &
          #1.000 #0.026 SVOSOMG_09_05 + #1.000 #0.026 SVOSOMG_10_05 + #1.000 &
          #0.026 SVOSOMG_11_05 + #1.000 #0.026 SVOSOMG_12_05 + #1.000 #0.026 &
          SVOSOMG_13_05 + #1.000 #0.000 SVOSOMG_14_05 + #1.000 #0.026 SVOSOMG_03_06 + &
          #1.000 #0.026 SVOSOMG_04_06 + #1.000 #0.026 SVOSOMG_05_06 + #1.000 &
          #0.026 SVOSOMG_06_06 + #1.000 #0.026 SVOSOMG_07_06 + #1.000 #0.026 &
          SVOSOMG_08_06 + #1.000 #0.026 SVOSOMG_09_06 + #1.000 #0.026 SVOSOMG_10_06 + &
          #1.000 #0.026 SVOSOMG_11_06 + #1.000 #0.026 SVOSOMG_12_06 + #1.000 &
          #0.026 SVOSOMG_13_06 + #1.000 #0.000 SVOSOMG_14_06 + #1.000 #0.026 &
          SVOSOMG_04_07 + #1.000 #0.026 SVOSOMG_05_07 + #1.000 #0.026 SVOSOMG_06_07 + &
          #1.000 #0.026 SVOSOMG_07_07 + #1.000 #0.026 SVOSOMG_08_07 + #1.000 &
          #0.026 SVOSOMG_09_07 + #1.000 #0.026 SVOSOMG_10_07 + #1.000 #0.026 &
          SVOSOMG_11_07 + #1.000 #0.026 SVOSOMG_12_07 + #1.000 #0.026 SVOSOMG_13_07 + &
          #1.000 #0.026 SVOSOMG_14_07


These are the expanded type reactions

0LIST OF REACTANTS INPUT

         DECLARED TYPES:  1 = CONSTANT,  2 = BUILDUP,  3 = DUMMY, 4 = ACTIVE,  5 = ST. STATE,  6 = UNSPECIFIED

      1) [4] O3                   2) [4] NO                   3) [4] NO2                  4) [1] O2              
      5) [1] M                    6) [1] H2O                  7) [1] H2                   8) [1] HV              
      9) [1] OH                  10) [4] NO3                 11) [4] N2O5                12) [4] HNO3            
     13) [4] HONO                14) [4] HNO4                15) [4] SO2                 16) [4] SULF            
     17) [4] NPHE                18) [4] XN                  19) [4] RNO3                20) [4] PAN             
     21) [4] PAN2                22) [4] PBZN                23) [4] MAPAN               24) [4] CLNO            
     25) [4] CLONO               26) [4] CLNO2               27) [4] CLONO2              28) [4] ALK5            
     29) [4] BENZ                30) [4] BENZENE             31) [4] TRIMBNZ             32) [4] TOLUENE         
     33) [4] XYLENE              34) [4] ISOPRENE            35) [4] TERPENE             36) [4] SESQTRP         
     37) [4] ISOPRN              38) [4] TERP                39) [4] SESQ                40) [4] ARO1            
     41) [4] ARO2                42) [4] IVOC                43) [4] SVOC                44) [4] TOLSOMG         
     45) [4] XYLSOMG             46) [4] BNZSOMG             47) [4] ISPSOMG             48) [4] TRPSOMG         
     49) [4] IVOSOMG             50) [4] SVOSOMG             51) [5] O3P                 52) [5] O1D             
     53) [5] TBUO                54) [5] BZO                 55) [5] NO2EX               56) [5] xHO2            
     57) [5] xOH                 58) [5] xNO2                59) [5] xMEO2               60) [5] xMECO3          
     61) [5] xRCO3               62) [5] xMACO3              63) [5] xTBUO               64) [5] xCO             
     65) [5] xHCHO               66) [5] xCCHO               67) [5] xRCHO               68) [5] xACET           
     69) [5] xMEK                70) [5] xPROD2              71) [5] xGLY                72) [5] xMGLY           
     73) [5] xBACL               74) [5] xBALD               75) [5] xAFG1               76) [5] xAFG2           
     77) [5] xMACR               78) [5] xMVK                79) [5] xIPRD               80) [5] xRNO3           
     81) [5] zRNO3               82) [5] yROOH               83) [5] yR6OOH              84) [5] yRAOOH          
     85) [5] xCL                 86) [5] xCLCCHO             87) [5] xCLACET             88) [6] HO2             
     89) [6] CO                  90) [6] CO2                 91) [6] HO2H                92) [6] MEO2            
     93) [6] HCHO                94) [6] COOH                95) [6] MEOH                96) [6] RO2C            
     97) [6] RO2XC               98) [6] MECO3               99) [6] CCOOH              100) [6] RCO3            
    101) [6] N3                 102) [6] RCOOH              103) [6] BZCO3              104) [6] XC              
    105) [6] MACO3              106) [6] ACET               107) [6] CRES               108) [6] CCHO            
    109) [6] RCHO               110) [6] MEK                111) [6] HCOOH              112) [6] ROOH            
    113) [6] R6OOH              114) [6] PROD2              115) [6] RAOOH              116) [6] MGLY            
    117) [6] IPRD               118) [6] GLY                119) [6] AFG1               120) [6] AFG2            
    121) [6] HCOCO3             122) [6] BACL               123) [6] CATL               124) [6] AFG3            
    125) [6] xCNDPP             126) [6] yRAOOHP            127) [6] XYNL               128) [6] CNDPP           
    129) [6] BALD               130) [6] PHEN               131) [6] MACR               132) [6] MVK             
    133) [6] AFG4               134) [6] ETOH               135) [6] xAFG4              136) [6] CH4             
    137) [6] ETHENE             138) [6] ACETYLEN           139) [6] ALK1               140) [6] ALK2            
    141) [6] ALK3               142) [6] ALK4               143) [6] OLE1               144) [6] OLE2            
    145) [6] CL2                146) [6] CL                 147) [4] HCL                148) [6] CLO             
    149) [6] HOCL               150) [6] CLCCHO             151) [6] CLACET             152) [6] CLCHO           
    153) [4] NH3                154) [4] GDENS              155) [4] MMOLES             156) [4] XICBC           
    157) [4] YICBC              158) [4] ZICBC              159) [6] CHCL3              160) [4] BNZSOMG_01_00   
    161) [4] BNZSOMG_01_01      162) [4] BNZSOMG_01_02      163) [4] BNZSOMG_02_00      164) [4] BNZSOMG_02_01   
    165) [4] BNZSOMG_02_02      166) [4] BNZSOMG_02_03      167) [4] BNZSOMG_02_04      168) [4] BNZSOMG_03_00   
    169) [4] BNZSOMG_03_01      170) [4] BNZSOMG_03_02      171) [4] BNZSOMG_03_03      172) [4] BNZSOMG_03_04   
    173) [4] BNZSOMG_03_05      174) [4] BNZSOMG_03_06      175) [4] BNZSOMG_04_00      176) [4] BNZSOMG_04_01   
    177) [4] BNZSOMG_04_02      178) [4] BNZSOMG_04_03      179) [4] BNZSOMG_04_04      180) [4] BNZSOMG_04_05   
    181) [4] BNZSOMG_04_06      182) [4] BNZSOMG_04_07      183) [4] BNZSOMG_05_00      184) [4] BNZSOMG_05_01   
    185) [4] BNZSOMG_05_02      186) [4] BNZSOMG_05_03      187) [4] BNZSOMG_05_04      188) [4] BNZSOMG_05_05   
    189) [4] BNZSOMG_05_06      190) [4] BNZSOMG_05_07      191) [4] BNZSOMG_06_00      192) [4] BNZSOMG_06_01   
    193) [4] BNZSOMG_06_02      194) [4] BNZSOMG_06_03      195) [4] BNZSOMG_06_04      196) [4] BNZSOMG_06_05   
    197) [4] BNZSOMG_06_06      198) [4] BNZSOMG_06_07      199) [4] BNZSOMG_07_00      200) [4] BNZSOMG_07_01   
    201) [4] BNZSOMG_07_02      202) [4] BNZSOMG_07_03      203) [4] BNZSOMG_07_04      204) [4] BNZSOMG_07_05   
    205) [4] BNZSOMG_07_06      206) [4] BNZSOMG_07_07      207) [4] BNZSOMG_08_00      208) [4] BNZSOMG_08_01   
    209) [4] BNZSOMG_08_02      210) [4] BNZSOMG_08_03      211) [4] BNZSOMG_08_04      212) [4] BNZSOMG_08_05   
    213) [4] BNZSOMG_08_06      214) [4] BNZSOMG_08_07      215) [4] BNZSOMG_09_00      216) [4] BNZSOMG_09_01   
    217) [4] BNZSOMG_09_02      218) [4] BNZSOMG_09_03      219) [4] BNZSOMG_09_04      220) [4] BNZSOMG_09_05   
    221) [4] BNZSOMG_09_06      222) [4] BNZSOMG_09_07      223) [4] TOLSOMG_01_00      224) [4] TOLSOMG_01_01   
    225) [4] TOLSOMG_01_02      226) [4] TOLSOMG_02_00      227) [4] TOLSOMG_02_01      228) [4] TOLSOMG_02_02   
    229) [4] TOLSOMG_02_03      230) [4] TOLSOMG_02_04      231) [4] TOLSOMG_03_00      232) [4] TOLSOMG_03_01   
    233) [4] TOLSOMG_03_02      234) [4] TOLSOMG_03_03      235) [4] TOLSOMG_03_04      236) [4] TOLSOMG_03_05   
    237) [4] TOLSOMG_03_06      238) [4] TOLSOMG_04_00      239) [4] TOLSOMG_04_01      240) [4] TOLSOMG_04_02   
    241) [4] TOLSOMG_04_03      242) [4] TOLSOMG_04_04      243) [4] TOLSOMG_04_05      244) [4] TOLSOMG_04_06   
    245) [4] TOLSOMG_04_07      246) [4] TOLSOMG_05_00      247) [4] TOLSOMG_05_01      248) [4] TOLSOMG_05_02   
    249) [4] TOLSOMG_05_03      250) [4] TOLSOMG_05_04      251) [4] TOLSOMG_05_05      252) [4] TOLSOMG_05_06   
    253) [4] TOLSOMG_05_07      254) [4] TOLSOMG_06_00      255) [4] TOLSOMG_06_01      256) [4] TOLSOMG_06_02   
    257) [4] TOLSOMG_06_03      258) [4] TOLSOMG_06_04      259) [4] TOLSOMG_06_05      260) [4] TOLSOMG_06_06   
    261) [4] TOLSOMG_06_07      262) [4] TOLSOMG_07_00      263) [4] TOLSOMG_07_01      264) [4] TOLSOMG_07_02   
    265) [4] TOLSOMG_07_03      266) [4] TOLSOMG_07_04      267) [4] TOLSOMG_07_05      268) [4] TOLSOMG_07_06   
    269) [4] TOLSOMG_07_07      270) [4] XYLSOMG_01_00      271) [4] XYLSOMG_01_01      272) [4] XYLSOMG_01_02   
    273) [4] XYLSOMG_02_00      274) [4] XYLSOMG_02_01      275) [4] XYLSOMG_02_02      276) [4] XYLSOMG_02_03   
    277) [4] XYLSOMG_02_04      278) [4] XYLSOMG_03_00      279) [4] XYLSOMG_03_01      280) [4] XYLSOMG_03_02   
    281) [4] XYLSOMG_03_03      282) [4] XYLSOMG_03_04      283) [4] XYLSOMG_03_05      284) [4] XYLSOMG_03_06   
    285) [4] XYLSOMG_04_00      286) [4] XYLSOMG_04_01      287) [4] XYLSOMG_04_02      288) [4] XYLSOMG_04_03   
    289) [4] XYLSOMG_04_04      290) [4] XYLSOMG_04_05      291) [4] XYLSOMG_04_06      292) [4] XYLSOMG_04_07   
    293) [4] XYLSOMG_05_00      294) [4] XYLSOMG_05_01      295) [4] XYLSOMG_05_02      296) [4] XYLSOMG_05_03   
    297) [4] XYLSOMG_05_04      298) [4] XYLSOMG_05_05      299) [4] XYLSOMG_05_06      300) [4] XYLSOMG_05_07   
    301) [4] XYLSOMG_06_00      302) [4] XYLSOMG_06_01      303) [4] XYLSOMG_06_02      304) [4] XYLSOMG_06_03   
    305) [4] XYLSOMG_06_04      306) [4] XYLSOMG_06_05      307) [4] XYLSOMG_06_06      308) [4] XYLSOMG_06_07   
    309) [4] XYLSOMG_07_00      310) [4] XYLSOMG_07_01      311) [4] XYLSOMG_07_02      312) [4] XYLSOMG_07_03   
    313) [4] XYLSOMG_07_04      314) [4] XYLSOMG_07_05      315) [4] XYLSOMG_07_06      316) [4] XYLSOMG_07_07   
    317) [4] XYLSOMG_08_00      318) [4] XYLSOMG_08_01      319) [4] XYLSOMG_08_02      320) [4] XYLSOMG_08_03   
    321) [4] XYLSOMG_08_04      322) [4] XYLSOMG_08_05      323) [4] XYLSOMG_08_06      324) [4] XYLSOMG_08_07   
    325) [4] ISPSOMG_01_00      326) [4] ISPSOMG_01_01      327) [4] ISPSOMG_01_02      328) [4] ISPSOMG_02_00   
    329) [4] ISPSOMG_02_01      330) [4] ISPSOMG_02_02      331) [4] ISPSOMG_02_03      332) [4] ISPSOMG_02_04   
    333) [4] ISPSOMG_03_00      334) [4] ISPSOMG_03_01      335) [4] ISPSOMG_03_02      336) [4] ISPSOMG_03_03   
    337) [4] ISPSOMG_03_04      338) [4] ISPSOMG_03_05      339) [4] ISPSOMG_03_06      340) [4] ISPSOMG_04_00   
    341) [4] ISPSOMG_04_01      342) [4] ISPSOMG_04_02      343) [4] ISPSOMG_04_03      344) [4] ISPSOMG_04_04   
    345) [4] ISPSOMG_04_05      346) [4] ISPSOMG_04_06      347) [4] ISPSOMG_04_07      348) [4] ISPSOMG_05_00   
    349) [4] ISPSOMG_05_01      350) [4] ISPSOMG_05_02      351) [4] ISPSOMG_05_03      352) [4] ISPSOMG_05_04   
    353) [4] ISPSOMG_05_05      354) [4] ISPSOMG_05_06      355) [4] ISPSOMG_05_07      356) [4] TRPSOMG_01_00   
    357) [4] TRPSOMG_01_01      358) [4] TRPSOMG_01_02      359) [4] TRPSOMG_02_00      360) [4] TRPSOMG_02_01   
    361) [4] TRPSOMG_02_02      362) [4] TRPSOMG_02_03      363) [4] TRPSOMG_02_04      364) [4] TRPSOMG_03_00   
    365) [4] TRPSOMG_03_01      366) [4] TRPSOMG_03_02      367) [4] TRPSOMG_03_03      368) [4] TRPSOMG_03_04   
    369) [4] TRPSOMG_03_05      370) [4] TRPSOMG_03_06      371) [4] TRPSOMG_04_00      372) [4] TRPSOMG_04_01   
    373) [4] TRPSOMG_04_02      374) [4] TRPSOMG_04_03      375) [4] TRPSOMG_04_04      376) [4] TRPSOMG_04_05   
    377) [4] TRPSOMG_04_06      378) [4] TRPSOMG_04_07      379) [4] TRPSOMG_05_00      380) [4] TRPSOMG_05_01   
    381) [4] TRPSOMG_05_02      382) [4] TRPSOMG_05_03      383) [4] TRPSOMG_05_04      384) [4] TRPSOMG_05_05   
    385) [4] TRPSOMG_05_06      386) [4] TRPSOMG_05_07      387) [4] TRPSOMG_06_00      388) [4] TRPSOMG_06_01   
    389) [4] TRPSOMG_06_02      390) [4] TRPSOMG_06_03      391) [4] TRPSOMG_06_04      392) [4] TRPSOMG_06_05   
    393) [4] TRPSOMG_06_06      394) [4] TRPSOMG_06_07      395) [4] TRPSOMG_07_00      396) [4] TRPSOMG_07_01   
    397) [4] TRPSOMG_07_02      398) [4] TRPSOMG_07_03      399) [4] TRPSOMG_07_04      400) [4] TRPSOMG_07_05   
    401) [4] TRPSOMG_07_06      402) [4] TRPSOMG_07_07      403) [4] TRPSOMG_08_00      404) [4] TRPSOMG_08_01   
    405) [4] TRPSOMG_08_02      406) [4] TRPSOMG_08_03      407) [4] TRPSOMG_08_04      408) [4] TRPSOMG_08_05   
    409) [4] TRPSOMG_08_06      410) [4] TRPSOMG_08_07      411) [4] TRPSOMG_09_00      412) [4] TRPSOMG_09_01   
    413) [4] TRPSOMG_09_02      414) [4] TRPSOMG_09_03      415) [4] TRPSOMG_09_04      416) [4] TRPSOMG_09_05   
    417) [4] TRPSOMG_09_06      418) [4] TRPSOMG_09_07      419) [4] TRPSOMG_10_01      420) [4] TRPSOMG_10_02   
    421) [4] TRPSOMG_10_03      422) [4] TRPSOMG_10_04      423) [4] TRPSOMG_10_05      424) [4] TRPSOMG_10_06   
    425) [4] TRPSOMG_10_07      426) [4] TRPSOMG_11_00      427) [4] TRPSOMG_11_01      428) [4] TRPSOMG_11_02   
    429) [4] TRPSOMG_11_03      430) [4] TRPSOMG_11_04      431) [4] TRPSOMG_11_05      432) [4] TRPSOMG_11_06   
    433) [4] TRPSOMG_11_07      434) [4] TRPSOMG_12_00      435) [4] TRPSOMG_12_01      436) [4] TRPSOMG_12_02   
    437) [4] TRPSOMG_12_03      438) [4] TRPSOMG_12_04      439) [4] TRPSOMG_12_05      440) [4] TRPSOMG_12_06   
    441) [4] TRPSOMG_12_07      442) [4] TRPSOMG_13_00      443) [4] TRPSOMG_13_01      444) [4] TRPSOMG_13_02   
    445) [4] TRPSOMG_13_03      446) [4] TRPSOMG_13_04      447) [4] TRPSOMG_13_05      448) [4] TRPSOMG_13_06   
    449) [4] TRPSOMG_13_07      450) [4] TRPSOMG_14_00      451) [4] TRPSOMG_14_01      452) [4] TRPSOMG_14_02   
    453) [4] TRPSOMG_14_03      454) [4] TRPSOMG_14_04      455) [4] TRPSOMG_14_05      456) [4] TRPSOMG_14_06   
    457) [4] TRPSOMG_14_07      458) [4] TRPSOMG_15_00      459) [4] TRPSOMG_15_01      460) [4] TRPSOMG_15_02   
    461) [4] TRPSOMG_15_03      462) [4] TRPSOMG_15_04      463) [4] TRPSOMG_15_05      464) [4] TRPSOMG_15_06   
    465) [4] TRPSOMG_15_07      466) [4] IVOSOMG_01_00      467) [4] IVOSOMG_01_01      468) [4] IVOSOMG_01_02   
    469) [4] IVOSOMG_02_00      470) [4] IVOSOMG_02_01      471) [4] IVOSOMG_02_02      472) [4] IVOSOMG_02_03   
    473) [4] IVOSOMG_02_04      474) [4] IVOSOMG_03_00      475) [4] IVOSOMG_03_01      476) [4] IVOSOMG_03_02   
    477) [4] IVOSOMG_03_03      478) [4] IVOSOMG_03_04      479) [4] IVOSOMG_03_05      480) [4] IVOSOMG_03_06   
    481) [4] IVOSOMG_04_00      482) [4] IVOSOMG_04_01      483) [4] IVOSOMG_04_02      484) [4] IVOSOMG_04_03   
    485) [4] IVOSOMG_04_04      486) [4] IVOSOMG_04_05      487) [4] IVOSOMG_04_06      488) [4] IVOSOMG_04_07   
    489) [4] IVOSOMG_05_00      490) [4] IVOSOMG_05_01      491) [4] IVOSOMG_05_02      492) [4] IVOSOMG_05_03   
    493) [4] IVOSOMG_05_04      494) [4] IVOSOMG_05_05      495) [4] IVOSOMG_05_06      496) [4] IVOSOMG_05_07   
    497) [4] IVOSOMG_06_00      498) [4] IVOSOMG_06_01      499) [4] IVOSOMG_06_02      500) [4] IVOSOMG_06_03   
    501) [4] IVOSOMG_06_04      502) [4] IVOSOMG_06_05      503) [4] IVOSOMG_06_06      504) [4] IVOSOMG_06_07   
    505) [4] IVOSOMG_07_00      506) [4] IVOSOMG_07_01      507) [4] IVOSOMG_07_02      508) [4] IVOSOMG_07_03   
    509) [4] IVOSOMG_07_04      510) [4] IVOSOMG_07_05      511) [4] IVOSOMG_07_06      512) [4] IVOSOMG_07_07   
    513) [4] IVOSOMG_08_00      514) [4] IVOSOMG_08_01      515) [4] IVOSOMG_08_02      516) [4] IVOSOMG_08_03   
    517) [4] IVOSOMG_08_04      518) [4] IVOSOMG_08_05      519) [4] IVOSOMG_08_06      520) [4] IVOSOMG_08_07   
    521) [4] IVOSOMG_09_00      522) [4] IVOSOMG_09_01      523) [4] IVOSOMG_09_02      524) [4] IVOSOMG_09_03   
    525) [4] IVOSOMG_09_04      526) [4] IVOSOMG_09_05      527) [4] IVOSOMG_09_06      528) [4] IVOSOMG_09_07   
    529) [4] IVOSOMG_10_01      530) [4] IVOSOMG_10_02      531) [4] IVOSOMG_10_03      532) [4] IVOSOMG_10_04   
    533) [4] IVOSOMG_10_05      534) [4] IVOSOMG_10_06      535) [4] IVOSOMG_10_07      536) [4] IVOSOMG_11_00   
    537) [4] IVOSOMG_11_01      538) [4] IVOSOMG_11_02      539) [4] IVOSOMG_11_03      540) [4] IVOSOMG_11_04   
    541) [4] IVOSOMG_11_05      542) [4] IVOSOMG_11_06      543) [4] IVOSOMG_11_07      544) [4] IVOSOMG_12_00   
    545) [4] IVOSOMG_12_01      546) [4] IVOSOMG_12_02      547) [4] IVOSOMG_12_03      548) [4] IVOSOMG_12_04   
    549) [4] IVOSOMG_12_05      550) [4] IVOSOMG_12_06      551) [4] IVOSOMG_12_07      552) [4] IVOSOMG_13_00   
    553) [4] IVOSOMG_13_01      554) [4] IVOSOMG_13_02      555) [4] IVOSOMG_13_03      556) [4] IVOSOMG_13_04   
    557) [4] IVOSOMG_13_05      558) [4] IVOSOMG_13_06      559) [4] IVOSOMG_13_07      560) [4] IVOSOMG_14_00   
    561) [4] IVOSOMG_14_01      562) [4] IVOSOMG_14_02      563) [4] IVOSOMG_14_03      564) [4] IVOSOMG_14_04   
    565) [4] IVOSOMG_14_05      566) [4] IVOSOMG_14_06      567) [4] IVOSOMG_14_07      568) [4] IVOSOMG_15_00   
    569) [4] IVOSOMG_15_01      570) [4] IVOSOMG_15_02      571) [4] IVOSOMG_15_03      572) [4] IVOSOMG_15_04   
    573) [4] IVOSOMG_15_05      574) [4] IVOSOMG_15_06      575) [4] IVOSOMG_15_07      576) [4] SVOSOMG_01_00   
    577) [4] SVOSOMG_01_01      578) [4] SVOSOMG_01_02      579) [4] SVOSOMG_02_00      580) [4] SVOSOMG_02_01   
    581) [4] SVOSOMG_02_02      582) [4] SVOSOMG_02_03      583) [4] SVOSOMG_02_04      584) [4] SVOSOMG_03_00   
    585) [4] SVOSOMG_03_01      586) [4] SVOSOMG_03_02      587) [4] SVOSOMG_03_03      588) [4] SVOSOMG_03_04   
    589) [4] SVOSOMG_03_05      590) [4] SVOSOMG_03_06      591) [4] SVOSOMG_04_00      592) [4] SVOSOMG_04_01   
    593) [4] SVOSOMG_04_02      594) [4] SVOSOMG_04_03      595) [4] SVOSOMG_04_04      596) [4] SVOSOMG_04_05   
    597) [4] SVOSOMG_04_06      598) [4] SVOSOMG_04_07      599) [4] SVOSOMG_05_00      600) [4] SVOSOMG_05_01   
    601) [4] SVOSOMG_05_02      602) [4] SVOSOMG_05_03      603) [4] SVOSOMG_05_04      604) [4] SVOSOMG_05_05   
    605) [4] SVOSOMG_05_06      606) [4] SVOSOMG_05_07      607) [4] SVOSOMG_06_00      608) [4] SVOSOMG_06_01   
    609) [4] SVOSOMG_06_02      610) [4] SVOSOMG_06_03      611) [4] SVOSOMG_06_04      612) [4] SVOSOMG_06_05   
    613) [4] SVOSOMG_06_06      614) [4] SVOSOMG_06_07      615) [4] SVOSOMG_07_00      616) [4] SVOSOMG_07_01   
    617) [4] SVOSOMG_07_02      618) [4] SVOSOMG_07_03      619) [4] SVOSOMG_07_04      620) [4] SVOSOMG_07_05   
    621) [4] SVOSOMG_07_06      622) [4] SVOSOMG_07_07      623) [4] SVOSOMG_08_00      624) [4] SVOSOMG_08_01   
    625) [4] SVOSOMG_08_02      626) [4] SVOSOMG_08_03      627) [4] SVOSOMG_08_04      628) [4] SVOSOMG_08_05   
    629) [4] SVOSOMG_08_06      630) [4] SVOSOMG_08_07      631) [4] SVOSOMG_09_00      632) [4] SVOSOMG_09_01   
    633) [4] SVOSOMG_09_02      634) [4] SVOSOMG_09_03      635) [4] SVOSOMG_09_04      636) [4] SVOSOMG_09_05   
    637) [4] SVOSOMG_09_06      638) [4] SVOSOMG_09_07      639) [4] SVOSOMG_10_01      640) [4] SVOSOMG_10_02   
    641) [4] SVOSOMG_10_03      642) [4] SVOSOMG_10_04      643) [4] SVOSOMG_10_05      644) [4] SVOSOMG_10_06   
    645) [4] SVOSOMG_10_07      646) [4] SVOSOMG_11_00      647) [4] SVOSOMG_11_01      648) [4] SVOSOMG_11_02   
    649) [4] SVOSOMG_11_03      650) [4] SVOSOMG_11_04      651) [4] SVOSOMG_11_05      652) [4] SVOSOMG_11_06   
    653) [4] SVOSOMG_11_07      654) [4] SVOSOMG_12_00      655) [4] SVOSOMG_12_01      656) [4] SVOSOMG_12_02   
    657) [4] SVOSOMG_12_03      658) [4] SVOSOMG_12_04      659) [4] SVOSOMG_12_05      660) [4] SVOSOMG_12_06   
    661) [4] SVOSOMG_12_07      662) [4] SVOSOMG_13_00      663) [4] SVOSOMG_13_01      664) [4] SVOSOMG_13_02   
    665) [4] SVOSOMG_13_03      666) [4] SVOSOMG_13_04      667) [4] SVOSOMG_13_05      668) [4] SVOSOMG_13_06   
    669) [4] SVOSOMG_13_07      670) [4] SVOSOMG_14_00      671) [4] SVOSOMG_14_01      672) [4] SVOSOMG_14_02   
    673) [4] SVOSOMG_14_03      674) [4] SVOSOMG_14_04      675) [4] SVOSOMG_14_05      676) [4] SVOSOMG_14_06   
    677) [4] SVOSOMG_14_07      678) [4] SVOSOMG_15_00      679) [4] SVOSOMG_15_01      680) [4] SVOSOMG_15_02   
    681) [4] SVOSOMG_15_03      682) [4] SVOSOMG_15_04      683) [4] SVOSOMG_15_05      684) [4] SVOSOMG_15_06   
    685) [4] SVOSOMG_15_07   
0INITIAL DIFFUN CODE INSERTED:
       'RO2NO' = "K BR07"*"C NO"
       'RO2HO2' = "K BR08"*"C HO2"
       'RO2NO3' = "K BR09"*"C NO3"
       'RO2RO2' = ("K BR10" * "C MEO2") + ("K BR11" * ("C RO2C" + "C RO2XC"))
       'RO2RO3' = "K BR25" * ("C MECO3" + "C RCO3" + "C BZCO3" + "C MACO3")
       'RO2RO' = 'RO2NO' + 'RO2NO3' + 'RO2RO3' + 0.5*'RO2RO2'
       'RO2XRO' = 'RO2HO2' + (0.5*'RO2RO2')
       'RO2RO2M' = (0.5*'RO2RO2')
       'RO22NN' = 'RO2NO3' + 'RO2RO3' + (0.5*'RO2RO2')
           2  ITERATIONS REQUIRED TO ORDER STEADY STATE SPECIES.
0CODE PREPARED.  LENGTH=   8665 STATEMENTS
1          29-Nov-21 10:17 TEMP 300.                                                       
 
 
 LIST OF REACTANTS:  ( 556 ACTIVE,   6 CONSTANT,   7 BUILDUP,  36 STEADY STATE,   0 DUMMY, AND  80 DELETED
 
 TYPE           NO. NAME                 MWT    C'S   N'S  S'S     O'S   X'S   CONC(0)
 ------------   --- ----------------  -------  -----  ---  ---  -----  -----  -----------
 
 ACTIVE          14  O3                  48.00   0.00   0    0   3.00   0.00     0.00    
                 15  NO                  30.01   0.00   1    0   0.00   0.00     0.00    
                 16  NO2                 46.01   0.00   1    0   0.00   0.00     0.00    
                 17  NO3                 62.01   0.00   1    0   3.00   0.00     0.00    
                 18  N2O5               108.02   0.00   2    0   0.00   0.00     0.00    
                 19  HNO3                63.02   0.00   1    0   0.00   0.00     0.00    
                 20  HONO                47.02   0.00   1    0   0.00   0.00     0.00    
                 21  HNO4                79.02   0.00   1    0   0.00   0.00     0.00    
                 22  SO2                 64.06   0.00   0    1   0.00   0.00     0.00    
                 23  SULF                96.06   0.00   0    1   0.00   0.00     0.00    
                 24  NPHE               139.11   6.00   1    0   0.00   0.00     0.00    
                 25  XN                   0.00   0.00   1    0   0.00   0.00     0.00    
                 26  RNO3                 0.00   0.00   1    0   0.00   0.00     0.00    
                 27  PAN                121.05   2.00   1    0   0.00   0.00     0.00    
                 28  PAN2               135.08   3.00   1    0   0.00   0.00     0.00    
                 29  PBZN               183.13   7.00   1    0   0.00   0.00     0.00    
                 30  MAPAN              132.00   3.00   1    0   0.00   0.00     0.00    
                 31  CLNO                65.46   0.00   1    0   0.00   0.00     0.00    
                 32  CLONO               81.46   0.00   1    0   0.00   0.00     0.00    
                 33  CLNO2               81.46   0.00   1    0   0.00   0.00     0.00    
                 34  CLONO2              97.46   0.00   1    0   0.00   0.00     0.00    
                 35  ALK5               170.00  12.00   0    0   0.00   1.00     0.00    
                 36  BENZ                78.11   6.00   0    0   0.00   1.00     0.00    
                 37  BENZENE             78.11   6.00   0    0   0.00   1.00     0.00    
                 38  TRIMBNZ            120.19   9.00   0    0   0.00   1.00     0.00    
                 39  TOLUENE             92.14   7.00   0    0   0.00   1.00     0.00    
                 40  XYLENE             106.16   8.00   0    0   0.00   1.00     0.00    
                 41  ISOPRENE            68.12   5.00   0    0   0.00   1.00     0.00    
                 42  TERPENE              0.00  10.00   0    0   0.00   1.00     0.00    
                 43  SESQTRP              0.00  15.00   0    0   0.00   1.00     0.00    
                 44  ISOPRN              68.12   5.00   0    0   0.00   1.00     0.00    
                 45  TERP                 0.00  10.00   0    0   0.00   1.00     0.00    
                 46  ARO1                92.14   7.00   0    0   0.00   1.00     0.00    
                 47  ARO2               106.16   8.00   0    0   0.00   1.00     0.00    
                 48  IVOC                89.07   3.00   0    0   3.00   1.00     0.00    
                 49  SVOC               146.14   6.00   0    0   4.00   1.00     0.00    
                 50  HO2                  0.00   0.00   0    0   0.00   0.00     0.00    
                 51  CO                   0.00   0.00   0    0   0.00   0.00     0.00    
                 52  HO2H                 0.00   0.00   0    0   0.00   0.00     0.00    
                 53  MEO2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 54  HCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 55  COOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 56  MEOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 57  RO2C                 0.00   0.00   0    0   0.00   0.00     0.00    
                 58  RO2XC                0.00   0.00   0    0   0.00   0.00     0.00    
                 59  MECO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 60  CCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 61  RCO3                 0.00   0.00   0    0   0.00   0.00     0.00    
                 62  RCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 63  BZCO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 64  MACO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 65  ACET                 0.00   0.00   0    0   0.00   0.00     0.00    
                 66  CRES                 0.00   0.00   0    0   0.00   0.00     0.00    
                 67  CCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 68  RCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 69  MEK                  0.00   0.00   0    0   0.00   0.00     0.00    
                 70  HCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 71  ROOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 72  R6OOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 73  PROD2                0.00   0.00   0    0   0.00   0.00     0.00    
                 74  RAOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 75  MGLY                 0.00   0.00   0    0   0.00   0.00     0.00    
                 76  IPRD                 0.00   0.00   0    0   0.00   0.00     0.00    
                 77  GLY                  0.00   0.00   0    0   0.00   0.00     0.00    
                 78  AFG1                 0.00   0.00   0    0   0.00   0.00     0.00    
                 79  AFG2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 80  HCOCO3               0.00   0.00   0    0   0.00   0.00     0.00    
                 81  BACL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 82  CATL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 83  AFG3                 0.00   0.00   0    0   0.00   0.00     0.00    
                 84  XYNL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 85  BALD                 0.00   0.00   0    0   0.00   0.00     0.00    
                 86  PHEN                 0.00   0.00   0    0   0.00   0.00     0.00    
                 87  MACR                 0.00   0.00   0    0   0.00   0.00     0.00    
                 88  MVK                  0.00   0.00   0    0   0.00   0.00     0.00    
                 89  AFG4                 0.00   0.00   0    0   0.00   0.00     0.00    
                 90  ETOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 91  xAFG4                0.00   0.00   0    0   0.00   0.00     0.00    
                 92  CH4                  0.00   0.00   0    0   0.00   0.00     0.00    
                 93  ETHENE               0.00   0.00   0    0   0.00   0.00     0.00    
                 94  ACETYLEN             0.00   0.00   0    0   0.00   0.00     0.00    
                 95  ALK1                 0.00   0.00   0    0   0.00   0.00     0.00    
                 96  ALK2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 97  ALK3                 0.00   0.00   0    0   0.00   0.00     0.00    
                 98  ALK4                 0.00   0.00   0    0   0.00   0.00     0.00    
                 99  OLE1                 0.00   0.00   0    0   0.00   0.00     0.00    
                100  OLE2                 0.00   0.00   0    0   0.00   0.00     0.00    
                101  CL2                  0.00   0.00   0    0   0.00   0.00     0.00    
                102  CL                   0.00   0.00   0    0   0.00   0.00     0.00    
                103  HCL                 36.46   0.00   0    0   0.00   0.00     0.00    
                104  CLO                  0.00   0.00   0    0   0.00   0.00     0.00    
                105  HOCL                 0.00   0.00   0    0   0.00   0.00     0.00    
                106  CLCCHO               0.00   0.00   0    0   0.00   0.00     0.00    
                107  CLACET               0.00   0.00   0    0   0.00   0.00     0.00    
                108  NH3                 17.03   0.00   1    0   0.00   0.00     0.00    
                109  GDENS                0.00   0.00   0    0   0.00   0.00     0.00    
                110  MMOLES               0.00   0.00   0    0   0.00   0.00     0.00    
                111  XICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                112  YICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                113  ZICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                114  CHCL3                0.00   0.00   0    0   0.00   0.00     0.00    
                115  BNZSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                116  BNZSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                117  BNZSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                118  BNZSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                119  BNZSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                120  BNZSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                121  BNZSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                122  BNZSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                123  BNZSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                124  BNZSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                125  BNZSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                126  BNZSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                127  BNZSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                128  BNZSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                129  BNZSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                130  BNZSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                131  BNZSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                132  BNZSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                133  BNZSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                134  BNZSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                135  BNZSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                136  BNZSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                137  BNZSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                138  BNZSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                139  BNZSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                140  BNZSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                141  BNZSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                142  BNZSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                143  BNZSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                144  BNZSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                145  BNZSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                146  BNZSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                147  BNZSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                148  BNZSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                149  BNZSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                150  BNZSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                151  BNZSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                152  BNZSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                153  BNZSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                154  BNZSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                155  BNZSOMG_08_01      129.00   8.00   0    0   1.00   1.00     0.00    
                156  BNZSOMG_08_02      144.00   8.00   0    0   2.00   1.00     0.00    
                157  BNZSOMG_08_03      159.00   8.00   0    0   3.00   1.00     0.00    
                158  BNZSOMG_08_04      174.00   8.00   0    0   4.00   1.00     0.00    
                159  BNZSOMG_08_05      189.00   8.00   0    0   5.00   1.00     0.00    
                160  BNZSOMG_08_06      204.00   8.00   0    0   6.00   1.00     0.00    
                161  BNZSOMG_08_07      219.00   8.00   0    0   7.00   1.00     0.00    
                162  BNZSOMG_09_01      143.00   9.00   0    0   1.00   1.00     0.00    
                163  BNZSOMG_09_02      158.00   9.00   0    0   2.00   1.00     0.00    
                164  BNZSOMG_09_03      173.00   9.00   0    0   3.00   1.00     0.00    
                165  BNZSOMG_09_04      188.00   9.00   0    0   4.00   1.00     0.00    
                166  BNZSOMG_09_05      203.00   9.00   0    0   5.00   1.00     0.00    
                167  BNZSOMG_09_06      218.00   9.00   0    0   6.00   1.00     0.00    
                168  BNZSOMG_09_07      233.00   9.00   0    0   7.00   1.00     0.00    
                169  TOLSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                170  TOLSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                171  TOLSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                172  TOLSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                173  TOLSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                174  TOLSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                175  TOLSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                176  TOLSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                177  TOLSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                178  TOLSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                179  TOLSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                180  TOLSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                181  TOLSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                182  TOLSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                183  TOLSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                184  TOLSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                185  TOLSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                186  TOLSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                187  TOLSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                188  TOLSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                189  TOLSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                190  TOLSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                191  TOLSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                192  TOLSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                193  TOLSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                194  TOLSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                195  TOLSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                196  TOLSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                197  TOLSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                198  TOLSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                199  TOLSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                200  TOLSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                201  TOLSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                202  TOLSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                203  TOLSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                204  TOLSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                205  TOLSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                206  TOLSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                207  TOLSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                208  TOLSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                209  XYLSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                210  XYLSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                211  XYLSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                212  XYLSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                213  XYLSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                214  XYLSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                215  XYLSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                216  XYLSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                217  XYLSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                218  XYLSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                219  XYLSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                220  XYLSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                221  XYLSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                222  XYLSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                223  XYLSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                224  XYLSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                225  XYLSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                226  XYLSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                227  XYLSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                228  XYLSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                229  XYLSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                230  XYLSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                231  XYLSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                232  XYLSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                233  XYLSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                234  XYLSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                235  XYLSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                236  XYLSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                237  XYLSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                238  XYLSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                239  XYLSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                240  XYLSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                241  XYLSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                242  XYLSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                243  XYLSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                244  XYLSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                245  XYLSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                246  XYLSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                247  XYLSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                248  XYLSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                249  XYLSOMG_08_01      129.00   8.00   0    0   1.00   1.00     0.00    
                250  XYLSOMG_08_02      144.00   8.00   0    0   2.00   1.00     0.00    
                251  XYLSOMG_08_03      159.00   8.00   0    0   3.00   1.00     0.00    
                252  XYLSOMG_08_04      174.00   8.00   0    0   4.00   1.00     0.00    
                253  XYLSOMG_08_05      189.00   8.00   0    0   5.00   1.00     0.00    
                254  XYLSOMG_08_06      204.00   8.00   0    0   6.00   1.00     0.00    
                255  XYLSOMG_08_07      219.00   8.00   0    0   7.00   1.00     0.00    
                256  ISPSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                257  ISPSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                258  ISPSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                259  ISPSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                260  ISPSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                261  ISPSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                262  ISPSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                263  ISPSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                264  ISPSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                265  ISPSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                266  ISPSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                267  ISPSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                268  ISPSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                269  ISPSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                270  ISPSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                271  ISPSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                272  ISPSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                273  ISPSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                274  ISPSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                275  ISPSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                276  ISPSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                277  ISPSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                278  ISPSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                279  ISPSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                280  ISPSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                281  ISPSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                282  TRPSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                283  TRPSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                284  TRPSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                285  TRPSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                286  TRPSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                287  TRPSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                288  TRPSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                289  TRPSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                290  TRPSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                291  TRPSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                292  TRPSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                293  TRPSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                294  TRPSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                295  TRPSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                296  TRPSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                297  TRPSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                298  TRPSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                299  TRPSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                300  TRPSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                301  TRPSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                302  TRPSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                303  TRPSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                304  TRPSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                305  TRPSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                306  TRPSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                307  TRPSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                308  TRPSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                309  TRPSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                310  TRPSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                311  TRPSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                312  TRPSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                313  TRPSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                314  TRPSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                315  TRPSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                316  TRPSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                317  TRPSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                318  TRPSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                319  TRPSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                320  TRPSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                321  TRPSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                322  TRPSOMG_08_01      129.00   8.00   0    0   1.00   1.00     0.00    
                323  TRPSOMG_08_02      144.00   8.00   0    0   2.00   1.00     0.00    
                324  TRPSOMG_08_03      159.00   8.00   0    0   3.00   1.00     0.00    
                325  TRPSOMG_08_04      174.00   8.00   0    0   4.00   1.00     0.00    
                326  TRPSOMG_08_05      189.00   8.00   0    0   5.00   1.00     0.00    
                327  TRPSOMG_08_06      204.00   8.00   0    0   6.00   1.00     0.00    
                328  TRPSOMG_08_07      219.00   8.00   0    0   7.00   1.00     0.00    
                329  TRPSOMG_09_01      143.00   9.00   0    0   1.00   1.00     0.00    
                330  TRPSOMG_09_02      158.00   9.00   0    0   2.00   1.00     0.00    
                331  TRPSOMG_09_03      173.00   9.00   0    0   3.00   1.00     0.00    
                332  TRPSOMG_09_04      188.00   9.00   0    0   4.00   1.00     0.00    
                333  TRPSOMG_09_05      203.00   9.00   0    0   5.00   1.00     0.00    
                334  TRPSOMG_09_06      218.00   9.00   0    0   6.00   1.00     0.00    
                335  TRPSOMG_09_07      233.00   9.00   0    0   7.00   1.00     0.00    
                336  TRPSOMG_10_01      157.00  10.00   0    0   1.00   1.00     0.00    
                337  TRPSOMG_10_02      172.00  10.00   0    0   2.00   1.00     0.00    
                338  TRPSOMG_10_03      187.00  10.00   0    0   3.00   1.00     0.00    
                339  TRPSOMG_10_04      202.00  10.00   0    0   4.00   1.00     0.00    
                340  TRPSOMG_10_05      217.00  10.00   0    0   5.00   1.00     0.00    
                341  TRPSOMG_10_06      232.00  10.00   0    0   6.00   1.00     0.00    
                342  TRPSOMG_10_07      247.00  10.00   0    0   7.00   1.00     0.00    
                343  TRPSOMG_11_01      171.00  11.00   0    0   1.00   1.00     0.00    
                344  TRPSOMG_11_02      186.00  11.00   0    0   2.00   1.00     0.00    
                345  TRPSOMG_11_03      201.00  11.00   0    0   3.00   1.00     0.00    
                346  TRPSOMG_11_04      216.00  11.00   0    0   4.00   1.00     0.00    
                347  TRPSOMG_11_05      231.00  11.00   0    0   5.00   1.00     0.00    
                348  TRPSOMG_11_06      246.00  11.00   0    0   6.00   1.00     0.00    
                349  TRPSOMG_11_07      261.00  11.00   0    0   7.00   1.00     0.00    
                350  TRPSOMG_12_01      185.00  12.00   0    0   1.00   1.00     0.00    
                351  TRPSOMG_12_02      200.00  12.00   0    0   2.00   1.00     0.00    
                352  TRPSOMG_12_03      215.00  12.00   0    0   3.00   1.00     0.00    
                353  TRPSOMG_12_04      230.00  12.00   0    0   4.00   1.00     0.00    
                354  TRPSOMG_12_05      245.00  12.00   0    0   5.00   1.00     0.00    
                355  TRPSOMG_12_06      260.00  12.00   0    0   6.00   1.00     0.00    
                356  TRPSOMG_12_07      275.00  12.00   0    0   7.00   1.00     0.00    
                357  TRPSOMG_13_01      199.00  13.00   0    0   1.00   1.00     0.00    
                358  TRPSOMG_13_02      214.00  13.00   0    0   2.00   1.00     0.00    
                359  TRPSOMG_13_03      229.00  13.00   0    0   3.00   1.00     0.00    
                360  TRPSOMG_13_04      244.00  13.00   0    0   4.00   1.00     0.00    
                361  TRPSOMG_13_05      259.00  13.00   0    0   5.00   1.00     0.00    
                362  TRPSOMG_13_06      274.00  13.00   0    0   6.00   1.00     0.00    
                363  TRPSOMG_13_07      289.00  13.00   0    0   7.00   1.00     0.00    
                364  TRPSOMG_14_01      213.00  14.00   0    0   1.00   1.00     0.00    
                365  TRPSOMG_14_02      228.00  14.00   0    0   2.00   1.00     0.00    
                366  TRPSOMG_14_03      243.00  14.00   0    0   3.00   1.00     0.00    
                367  TRPSOMG_14_04      258.00  14.00   0    0   4.00   1.00     0.00    
                368  TRPSOMG_14_05      273.00  14.00   0    0   5.00   1.00     0.00    
                369  TRPSOMG_14_06      288.00  14.00   0    0   6.00   1.00     0.00    
                370  TRPSOMG_14_07      303.00  14.00   0    0   7.00   1.00     0.00    
                371  TRPSOMG_15_01      227.00  15.00   0    0   1.00   1.00     0.00    
                372  TRPSOMG_15_02      242.00  15.00   0    0   2.00   1.00     0.00    
                373  TRPSOMG_15_03      257.00  15.00   0    0   3.00   1.00     0.00    
                374  TRPSOMG_15_04      272.00  15.00   0    0   4.00   1.00     0.00    
                375  TRPSOMG_15_05      287.00  15.00   0    0   5.00   1.00     0.00    
                376  TRPSOMG_15_06      302.00  15.00   0    0   6.00   1.00     0.00    
                377  TRPSOMG_15_07      317.00  15.00   0    0   7.00   1.00     0.00    
                378  IVOSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                379  IVOSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                380  IVOSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                381  IVOSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                382  IVOSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                383  IVOSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                384  IVOSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                385  IVOSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                386  IVOSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                387  IVOSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                388  IVOSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                389  IVOSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                390  IVOSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                391  IVOSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                392  IVOSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                393  IVOSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                394  IVOSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                395  IVOSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                396  IVOSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                397  IVOSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                398  IVOSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                399  IVOSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                400  IVOSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                401  IVOSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                402  IVOSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                403  IVOSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                404  IVOSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                405  IVOSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                406  IVOSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                407  IVOSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                408  IVOSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                409  IVOSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                410  IVOSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                411  IVOSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                412  IVOSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                413  IVOSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                414  IVOSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                415  IVOSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                416  IVOSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                417  IVOSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                418  IVOSOMG_08_01      129.00   8.00   0    0   1.00   1.00     0.00    
                419  IVOSOMG_08_02      144.00   8.00   0    0   2.00   1.00     0.00    
                420  IVOSOMG_08_03      159.00   8.00   0    0   3.00   1.00     0.00    
                421  IVOSOMG_08_04      174.00   8.00   0    0   4.00   1.00     0.00    
                422  IVOSOMG_08_05      189.00   8.00   0    0   5.00   1.00     0.00    
                423  IVOSOMG_08_06      204.00   8.00   0    0   6.00   1.00     0.00    
                424  IVOSOMG_08_07      219.00   8.00   0    0   7.00   1.00     0.00    
                425  IVOSOMG_09_01      143.00   9.00   0    0   1.00   1.00     0.00    
                426  IVOSOMG_09_02      158.00   9.00   0    0   2.00   1.00     0.00    
                427  IVOSOMG_09_03      173.00   9.00   0    0   3.00   1.00     0.00    
                428  IVOSOMG_09_04      188.00   9.00   0    0   4.00   1.00     0.00    
                429  IVOSOMG_09_05      203.00   9.00   0    0   5.00   1.00     0.00    
                430  IVOSOMG_09_06      218.00   9.00   0    0   6.00   1.00     0.00    
                431  IVOSOMG_09_07      233.00   9.00   0    0   7.00   1.00     0.00    
                432  IVOSOMG_10_01      157.00  10.00   0    0   1.00   1.00     0.00    
                433  IVOSOMG_10_02      172.00  10.00   0    0   2.00   1.00     0.00    
                434  IVOSOMG_10_03      187.00  10.00   0    0   3.00   1.00     0.00    
                435  IVOSOMG_10_04      202.00  10.00   0    0   4.00   1.00     0.00    
                436  IVOSOMG_10_05      217.00  10.00   0    0   5.00   1.00     0.00    
                437  IVOSOMG_10_06      232.00  10.00   0    0   6.00   1.00     0.00    
                438  IVOSOMG_10_07      247.00  10.00   0    0   7.00   1.00     0.00    
                439  IVOSOMG_11_01      171.00  11.00   0    0   1.00   1.00     0.00    
                440  IVOSOMG_11_02      186.00  11.00   0    0   2.00   1.00     0.00    
                441  IVOSOMG_11_03      201.00  11.00   0    0   3.00   1.00     0.00    
                442  IVOSOMG_11_04      216.00  11.00   0    0   4.00   1.00     0.00    
                443  IVOSOMG_11_05      231.00  11.00   0    0   5.00   1.00     0.00    
                444  IVOSOMG_11_06      246.00  11.00   0    0   6.00   1.00     0.00    
                445  IVOSOMG_11_07      261.00  11.00   0    0   7.00   1.00     0.00    
                446  IVOSOMG_12_01      185.00  12.00   0    0   1.00   1.00     0.00    
                447  IVOSOMG_12_02      200.00  12.00   0    0   2.00   1.00     0.00    
                448  IVOSOMG_12_03      215.00  12.00   0    0   3.00   1.00     0.00    
                449  IVOSOMG_12_04      230.00  12.00   0    0   4.00   1.00     0.00    
                450  IVOSOMG_12_05      245.00  12.00   0    0   5.00   1.00     0.00    
                451  IVOSOMG_12_06      260.00  12.00   0    0   6.00   1.00     0.00    
                452  IVOSOMG_12_07      275.00  12.00   0    0   7.00   1.00     0.00    
                453  IVOSOMG_13_01      199.00  13.00   0    0   1.00   1.00     0.00    
                454  IVOSOMG_13_02      214.00  13.00   0    0   2.00   1.00     0.00    
                455  IVOSOMG_13_03      229.00  13.00   0    0   3.00   1.00     0.00    
                456  IVOSOMG_13_04      244.00  13.00   0    0   4.00   1.00     0.00    
                457  IVOSOMG_13_05      259.00  13.00   0    0   5.00   1.00     0.00    
                458  IVOSOMG_13_06      274.00  13.00   0    0   6.00   1.00     0.00    
                459  IVOSOMG_13_07      289.00  13.00   0    0   7.00   1.00     0.00    
                460  IVOSOMG_14_01      213.00  14.00   0    0   1.00   1.00     0.00    
                461  IVOSOMG_14_02      228.00  14.00   0    0   2.00   1.00     0.00    
                462  IVOSOMG_14_03      243.00  14.00   0    0   3.00   1.00     0.00    
                463  IVOSOMG_14_04      258.00  14.00   0    0   4.00   1.00     0.00    
                464  IVOSOMG_14_05      273.00  14.00   0    0   5.00   1.00     0.00    
                465  IVOSOMG_14_06      288.00  14.00   0    0   6.00   1.00     0.00    
                466  IVOSOMG_14_07      303.00  14.00   0    0   7.00   1.00     0.00    
                467  IVOSOMG_15_01      227.00  15.00   0    0   1.00   1.00     0.00    
                468  IVOSOMG_15_02      242.00  15.00   0    0   2.00   1.00     0.00    
                469  IVOSOMG_15_03      257.00  15.00   0    0   3.00   1.00     0.00    
                470  IVOSOMG_15_04      272.00  15.00   0    0   4.00   1.00     0.00    
                471  IVOSOMG_15_05      287.00  15.00   0    0   5.00   1.00     0.00    
                472  IVOSOMG_15_06      302.00  15.00   0    0   6.00   1.00     0.00    
                473  IVOSOMG_15_07      317.00  15.00   0    0   7.00   1.00     0.00    
                474  SVOSOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                475  SVOSOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                476  SVOSOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                477  SVOSOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                478  SVOSOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                479  SVOSOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                480  SVOSOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                481  SVOSOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                482  SVOSOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                483  SVOSOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                484  SVOSOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                485  SVOSOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                486  SVOSOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                487  SVOSOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                488  SVOSOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                489  SVOSOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                490  SVOSOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                491  SVOSOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                492  SVOSOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                493  SVOSOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                494  SVOSOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                495  SVOSOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                496  SVOSOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                497  SVOSOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                498  SVOSOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                499  SVOSOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                500  SVOSOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                501  SVOSOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                502  SVOSOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                503  SVOSOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                504  SVOSOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                505  SVOSOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                506  SVOSOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                507  SVOSOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                508  SVOSOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                509  SVOSOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                510  SVOSOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                511  SVOSOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                512  SVOSOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                513  SVOSOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
                514  SVOSOMG_08_01      129.00   8.00   0    0   1.00   1.00     0.00    
                515  SVOSOMG_08_02      144.00   8.00   0    0   2.00   1.00     0.00    
                516  SVOSOMG_08_03      159.00   8.00   0    0   3.00   1.00     0.00    
                517  SVOSOMG_08_04      174.00   8.00   0    0   4.00   1.00     0.00    
                518  SVOSOMG_08_05      189.00   8.00   0    0   5.00   1.00     0.00    
                519  SVOSOMG_08_06      204.00   8.00   0    0   6.00   1.00     0.00    
                520  SVOSOMG_08_07      219.00   8.00   0    0   7.00   1.00     0.00    
                521  SVOSOMG_09_01      143.00   9.00   0    0   1.00   1.00     0.00    
                522  SVOSOMG_09_02      158.00   9.00   0    0   2.00   1.00     0.00    
                523  SVOSOMG_09_03      173.00   9.00   0    0   3.00   1.00     0.00    
                524  SVOSOMG_09_04      188.00   9.00   0    0   4.00   1.00     0.00    
                525  SVOSOMG_09_05      203.00   9.00   0    0   5.00   1.00     0.00    
                526  SVOSOMG_09_06      218.00   9.00   0    0   6.00   1.00     0.00    
                527  SVOSOMG_09_07      233.00   9.00   0    0   7.00   1.00     0.00    
                528  SVOSOMG_10_01      157.00  10.00   0    0   1.00   1.00     0.00    
                529  SVOSOMG_10_02      172.00  10.00   0    0   2.00   1.00     0.00    
                530  SVOSOMG_10_03      187.00  10.00   0    0   3.00   1.00     0.00    
                531  SVOSOMG_10_04      202.00  10.00   0    0   4.00   1.00     0.00    
                532  SVOSOMG_10_05      217.00  10.00   0    0   5.00   1.00     0.00    
                533  SVOSOMG_10_06      232.00  10.00   0    0   6.00   1.00     0.00    
                534  SVOSOMG_10_07      247.00  10.00   0    0   7.00   1.00     0.00    
                535  SVOSOMG_11_01      171.00  11.00   0    0   1.00   1.00     0.00    
                536  SVOSOMG_11_02      186.00  11.00   0    0   2.00   1.00     0.00    
                537  SVOSOMG_11_03      201.00  11.00   0    0   3.00   1.00     0.00    
                538  SVOSOMG_11_04      216.00  11.00   0    0   4.00   1.00     0.00    
                539  SVOSOMG_11_05      231.00  11.00   0    0   5.00   1.00     0.00    
                540  SVOSOMG_11_06      246.00  11.00   0    0   6.00   1.00     0.00    
                541  SVOSOMG_11_07      261.00  11.00   0    0   7.00   1.00     0.00    
                542  SVOSOMG_12_01      185.00  12.00   0    0   1.00   1.00     0.00    
                543  SVOSOMG_12_02      200.00  12.00   0    0   2.00   1.00     0.00    
                544  SVOSOMG_12_03      215.00  12.00   0    0   3.00   1.00     0.00    
                545  SVOSOMG_12_04      230.00  12.00   0    0   4.00   1.00     0.00    
                546  SVOSOMG_12_05      245.00  12.00   0    0   5.00   1.00     0.00    
                547  SVOSOMG_12_06      260.00  12.00   0    0   6.00   1.00     0.00    
                548  SVOSOMG_12_07      275.00  12.00   0    0   7.00   1.00     0.00    
                549  SVOSOMG_13_01      199.00  13.00   0    0   1.00   1.00     0.00    
                550  SVOSOMG_13_02      214.00  13.00   0    0   2.00   1.00     0.00    
                551  SVOSOMG_13_03      229.00  13.00   0    0   3.00   1.00     0.00    
                552  SVOSOMG_13_04      244.00  13.00   0    0   4.00   1.00     0.00    
                553  SVOSOMG_13_05      259.00  13.00   0    0   5.00   1.00     0.00    
                554  SVOSOMG_13_06      274.00  13.00   0    0   6.00   1.00     0.00    
                555  SVOSOMG_13_07      289.00  13.00   0    0   7.00   1.00     0.00    
                556  SVOSOMG_14_01      213.00  14.00   0    0   1.00   1.00     0.00    
                557  SVOSOMG_14_02      228.00  14.00   0    0   2.00   1.00     0.00    
                558  SVOSOMG_14_03      243.00  14.00   0    0   3.00   1.00     0.00    
                559  SVOSOMG_14_04      258.00  14.00   0    0   4.00   1.00     0.00    
                560  SVOSOMG_14_05      273.00  14.00   0    0   5.00   1.00     0.00    
                561  SVOSOMG_14_06      288.00  14.00   0    0   6.00   1.00     0.00    
                562  SVOSOMG_14_07      303.00  14.00   0    0   7.00   1.00     0.00    
                563  SVOSOMG_15_01      227.00  15.00   0    0   1.00   1.00     0.00    
                564  SVOSOMG_15_02      242.00  15.00   0    0   2.00   1.00     0.00    
                565  SVOSOMG_15_03      257.00  15.00   0    0   3.00   1.00     0.00    
                566  SVOSOMG_15_04      272.00  15.00   0    0   4.00   1.00     0.00    
                567  SVOSOMG_15_05      287.00  15.00   0    0   5.00   1.00     0.00    
                568  SVOSOMG_15_06      302.00  15.00   0    0   6.00   1.00     0.00    
                569  SVOSOMG_15_07      317.00  15.00   0    0   7.00   1.00     0.00    
 
 CONSTANT         1  O2                   0.00   0.00   0    0   0.00   0.00    2.090E+05
                  2  M                    0.00   0.00   0    0   0.00   0.00    1.000E+06
                  3  H2O                  0.00   0.00   0    0   0.00   0.00    2.000E+04
                  4  H2                   0.00   0.00   0    0   0.00   0.00     0.00    
                  5  HV                   0.00   0.00   0    0   0.00   0.00     1.00    
                  6  OH                  17.00   0.00   0    0   1.00   0.00    1.015E-07
 
 BUILDUP-ONLY     7  CO2                  0.00   0.00   0    0   0.00   0.00     0.00    
                  8  N3                   0.00   0.00   0    0   0.00   0.00     0.00    
                  9  XC                   0.00   0.00   0    0   0.00   0.00     0.00    
                 10  xCNDPP               0.00   0.00   0    0   0.00   0.00     0.00    
                 11  yRAOOHP              0.00   0.00   0    0   0.00   0.00     0.00    
                 12  CNDPP                0.00   0.00   0    0   0.00   0.00     0.00    
                 13  CLCHO                0.00   0.00   0    0   0.00   0.00     0.00    
 
 STEADY STATE   570  O3P             
                571  O1D             
                572  TBUO            
                573  BZO             
                574  xHO2            
                575  xOH             
                576  xNO2            
                577  xMEO2           
                578  xMECO3          
                579  xRCO3           
                580  xMACO3          
                581  xTBUO           
                582  xCO             
                583  xHCHO           
                584  xCCHO           
                585  xRCHO           
                586  xACET           
                587  xMEK            
                588  xPROD2          
                589  xGLY            
                590  xMGLY           
                591  xBACL           
                592  xBALD           
                593  xAFG1           
                594  xAFG2           
                595  xMACR           
                596  xMVK            
                597  xIPRD           
                598  xRNO3           
                599  zRNO3           
                600  yROOH           
                601  yR6OOH          
                602  yRAOOH          
                603  xCL             
                604  xCLCCHO         
                605  xCLACET         
 
 DELETED        606  SESQ            
                607  TOLSOMG         
                608  XYLSOMG         
                609  BNZSOMG         
                610  ISPSOMG         
                611  TRPSOMG         
                612  IVOSOMG         
                613  SVOSOMG         
                614  NO2EX           
                615  BNZSOMG_01_00   
                616  BNZSOMG_02_00   
                617  BNZSOMG_03_00   
                618  BNZSOMG_04_00   
                619  BNZSOMG_05_00   
                620  BNZSOMG_06_00   
                621  BNZSOMG_07_00   
                622  BNZSOMG_08_00   
                623  BNZSOMG_09_00   
                624  TOLSOMG_01_00   
                625  TOLSOMG_02_00   
                626  TOLSOMG_03_00   
                627  TOLSOMG_04_00   
                628  TOLSOMG_05_00   
                629  TOLSOMG_06_00   
                630  TOLSOMG_07_00   
                631  XYLSOMG_01_00   
                632  XYLSOMG_02_00   
                633  XYLSOMG_03_00   
                634  XYLSOMG_04_00   
                635  XYLSOMG_05_00   
                636  XYLSOMG_06_00   
                637  XYLSOMG_07_00   
                638  XYLSOMG_08_00   
                639  ISPSOMG_01_00   
                640  ISPSOMG_02_00   
                641  ISPSOMG_03_00   
                642  ISPSOMG_04_00   
                643  ISPSOMG_05_00   
                644  TRPSOMG_01_00   
                645  TRPSOMG_02_00   
                646  TRPSOMG_03_00   
                647  TRPSOMG_04_00   
                648  TRPSOMG_05_00   
                649  TRPSOMG_06_00   
                650  TRPSOMG_07_00   
                651  TRPSOMG_08_00   
                652  TRPSOMG_09_00   
                653  TRPSOMG_11_00   
                654  TRPSOMG_12_00   
                655  TRPSOMG_13_00   
                656  TRPSOMG_14_00   
                657  TRPSOMG_15_00   
                658  IVOSOMG_01_00   
                659  IVOSOMG_02_00   
                660  IVOSOMG_03_00   
                661  IVOSOMG_04_00   
                662  IVOSOMG_05_00   
                663  IVOSOMG_06_00   
                664  IVOSOMG_07_00   
                665  IVOSOMG_08_00   
                666  IVOSOMG_09_00   
                667  IVOSOMG_11_00   
                668  IVOSOMG_12_00   
                669  IVOSOMG_13_00   
                670  IVOSOMG_14_00   
                671  IVOSOMG_15_00   
                672  SVOSOMG_01_00   
                673  SVOSOMG_02_00   
                674  SVOSOMG_03_00   
                675  SVOSOMG_04_00   
                676  SVOSOMG_05_00   
                677  SVOSOMG_06_00   
                678  SVOSOMG_07_00   
                679  SVOSOMG_08_00   
                680  SVOSOMG_09_00   
                681  SVOSOMG_11_00   
                682  SVOSOMG_12_00   
                683  SVOSOMG_13_00   
                684  SVOSOMG_14_00   
                685  SVOSOMG_15_00   

   20 VARIABLE COEFFICIENTS:
       1) LITTLE   1.000E-30    2) INERT     0.00        3) NROG      0.00        4) METHANE   0.00        5) RO2RO     0.00    
       6) RO2XRO    0.00        7) QY.ACET  5.000E-02    8) QY.MEK   0.175        9) QY.NO2A  1.500E-03   10) QY.NPHE  1.500E-02
      11) QY.BALD  6.000E-02   12) QY.PROD2 4.860E-03   13) RO2NO     0.00       14) RO22NN    0.00       15) RO2HO2    0.00    
      16) RO2RO2M   0.00       17) QY.CLAC  0.500       18) RO2NO3    0.00       19) RO2RO2    0.00       20) RO2RO3    0.00    

  736 CONSTANT COEFFICIENTS:
    2001)   2.00     2002)  0.610     2003)  0.390     2004)  0.800     2005)  0.200    
    2006)  0.500     2007)  0.750     2008)  0.250     2009)  0.600     2010)  0.400    
    2011)  0.440     2012)  0.410     2013)  0.150     2014)  0.100     2015)  0.900    
    2016)   2.24     2017)  0.560     2018)  -1.00     2019)   3.00     2020)   4.00    
    2021)  0.965     2022)  3.500E-02 2023)  0.620     2024)   1.38     2025)  0.380    
    2026)  0.967     2027)  3.900E-02 2028)  0.376     2029)  0.510     2030)  7.400E-02
    2031)  8.800E-02 2032)  0.504     2033)  0.300     2034)  0.509     2035)  0.491    
    2036) -0.491     2037)  0.143     2038)  0.142     2039)  0.457     2040) -0.455    
    2041)  0.700     2042)  0.744     2043)  0.251     2044)  4.000E-03 2045)  0.239    
    2046)  1.200E-02 2047)  0.205     2048)  3.400E-02 2049)  0.256     2050) -0.115    
    2051)  0.840     2052)  0.222     2053)  2.900E-02 2054)  9.000E-02 2055)  4.100E-02
    2056)  2.000E-02 2057)  7.500E-02 2058)  8.400E-02 2059)  0.160     2060)  0.782    
    2061)  7.700E-02 2062)  8.500E-02 2063)  2.600E-02 2064)  5.800E-02 2065)  0.698    
    2066)  0.858     2067)  1.700E-02 2068)  0.139     2069)  0.148     2070)  0.589    
    2071)  0.124     2072)  0.147     2073)  0.565     2074)  2.400E-02 2075)  0.448    
    2076)  3.000E-02 2077)  0.252     2078)  7.300E-02 2079)  0.713     2080)   2.67    
    2081)   1.40     2082)  0.170     2083)  0.130     2084)  4.000E-02 2085)  0.115    
    2086)  0.715     2087)  7.000E-02 2088)  0.230     2089)  0.110     2090)  0.210    
    2091)  0.995     2092)  0.930     2093)   1.65     2094)  0.330     2095) -8.000E-02
    2096)   1.90     2097)   6.00     2098)   7.00     2099)  9.500E-02 2100)  0.105    
    2101)  4.800E-02 2102)  6.500E-02 2103) -0.205     2104)  0.217     2105)  0.723    
    2106)  6.000E-02 2107)  0.521     2108)  0.201     2109)  0.334     2110)  0.407    
    2111)  0.129     2112)  0.107     2113)  0.267     2114)  0.783     2115)  0.284    
    2116)   1.02     2117)  0.173     2118)  0.305     2119)  0.695     2120)  0.195    
    2121)  0.206     2122)  0.733     2123)  0.117     2124)  0.561     2125)  0.114    
    2126)  0.274     2127)  0.153     2128)  1.900E-02 2129)  0.231     2130)  0.794    
    2131)  0.938     2132)  0.471     2133)  0.554     2134)  1.300E-02 2135)  0.258    
    2136)  7.000E-03 2137)  0.580     2138)  0.190     2139)  0.366     2140)  0.184    
    2141)  0.350     2142)  3.000E-03 2143)  0.163     2144)  0.264     2145) -0.575    
    2146)  0.902     2147)  9.800E-02 2148)  0.416     2149) -0.416     2150)  0.208    
    2151)  0.108     2152)  0.450     2153)  0.333     2154) -0.100     2155)   1.50    
    2156)  0.670     2157)  0.340     2158)  0.975     2159)  2.500E-02 2160)  0.675    
    2161) -0.725     2162)  0.164     2163)  6.400E-02 2164)  5.000E-02 2165)  0.475    
    2166)  0.950     2167)  0.351     2168) -5.000E-02 2169)  0.550     2170)  -2.20    
    2171)  0.289     2172)  0.336     2173)  5.500E-02 2174)  0.332     2175)  0.174    
    2176) -0.504     2177)  0.711     2178)  0.285     2179)  0.498     2180)  0.140    
    2181)  2.300E-02 2182)  0.742     2183)  0.372     2184)  4.700E-02 2185)  1.000E-03
    2186) -0.329     2187)  0.799     2188)  5.100E-02 2189)  0.572     2190)  0.227    
    2191)  0.218     2192)  8.000E-03 2193)  0.850     2194)  0.278     2195)  0.185    
    2196)   1.23     2197)  0.467     2198)  0.233     2199) -0.233     2200)  0.472    
    2201)  0.379     2202)  4.900E-02 2203)  0.473     2204)  7.100E-02 2205)  2.000E-03
    2206)  0.211     2207)  8.300E-02 2208)  0.402     2209)  0.329     2210)  0.528    
    2211)  0.877     2212)  0.913     2213)   1.59     2214)  8.700E-02 2215)  0.303    
    2216)  0.780     2217) -9.100E-02 2218)  0.189     2219)  0.313     2220)  0.976    
    2221)  0.175     2222)  1.100E-02 2223)  0.429     2224)  3.600E-02 2225)  1.000E-02
    2226)  3.100E-02 2227)  0.157     2228)  0.636     2229)  0.344     2230)  0.721    
    2231)  0.102     2232)  6.100E-02 2233)  0.214     2234)  6.300E-02 2235)  0.261    
    2236)  6.600E-02 2237)  0.591     2238)  0.396     2239)  8.100E-02 2240) -1.000E-03
    2241)   5.00     2242)  -3.00     2243)  -4.00     2244)  -6.00     2245)  -8.00    
    2246)   1.61     2247)  0.120     2248)  0.370     2249)  0.290     2250)  0.986    
    2251)  9.300E-02 2252)  0.907     2253)  0.624     2254)  0.320     2255)  0.357    
    2256)  0.833     2257)  0.266     2258)  0.192     2259)  0.275     2260)  0.122    
    2261)  0.204     2262)  0.401     2263)  0.936     2264)  0.749     2265)  0.187    
    2266)  0.813     2267)  0.240     2268)  0.116     2269)  0.570     2270)  0.314    
    2271)  0.704     2272) -0.105     2273)  0.236     2274)   1.25     2275)  0.445    
    2276)  0.983     2277) -4.600E-02 2278)  0.830     2279)   1.76     2280)  0.149    
    2281)  0.438     2282)  0.426     2283)  0.106     2284)  0.146     2285) -0.119    
    2286)  0.647     2287)   1.61     2288)  0.353     2289)  0.209     2290)  8.600E-02
    2291)   2.00     2292)  0.904     2293)   1.14     2294)  0.301     2295)  0.470    
    2296)  5.000E-03 2297)  6.000E-03 2298)  0.119     2299)  0.413     2300)  0.587    
    2301)  0.822     2302)  0.193     2303)  0.104     2304)  0.368     2305)  0.125    
    2306)  2.200E-02 2307)  0.112     2308)  3.700E-02 2309)  0.690     2310)  0.824    
    2311)   1.31     2312)  0.176     2313)  9.000E-03 2314)  0.546     2315)  0.454    
    2316)  0.437     2317)  0.113     2318)   1.22     2319)  0.914     2320)  0.966    
    2321)  0.788     2322)  0.481     2323)  0.136     2324)  7.600E-02 2325)  2.700E-02
    2326)  0.643     2327)  0.111     2328)  0.423     2329)  0.297     2330)  0.162    
    2331)  0.152     2332)  0.428     2333)  6.700E-02 2334)  0.315     2335)  1.800E-02
    2336)  3.300E-02 2337)  2.800E-02 2338)  2.100E-02 2339)  4.200E-02 2340)  6.900E-02
    2341)  0.128     2342)  0.409     2343)   1.18     2344)  0.154     2345)  0.322    
    2346)  0.269     2347)  1.400E-02 2348)  0.709     2349)  0.202     2350)  0.666    
    2351)  0.166     2352)  0.482     2353)  6.800E-02 2354)  0.138     2355)  0.403    
    2356)  0.286     2357)  8.900E-02 2358)  0.601     2359)   1.49     2360)  0.759    
    2361)   1.15     2362)  0.533     2363)  0.255     2364)   5.06     2365)  5.200E-02
    2366)  0.585     2367)  0.126     2368)  0.875     2369)  0.203     2370)  4.500E-02
    2371)  7.900E-02 2372)  0.220     2373)  0.165     2374)  4.300E-02 2375)  0.545    
    2376)   3.53     2377)  0.421     2378)   1.51     2379)  0.397     2380)   4.47    
    2381)  0.853     2382)   4.44     2383)   1.42     2384)  0.763     2385)  3.800E-02
    2386)   1.28     2387)  0.159     2388)  0.547     2389)  4.600E-02 2390)  0.908    
    2391)  0.851     2392)  0.680     2393)  0.198     2394)  0.541     2395)  0.237    
    2396)  0.109     2397)  0.686     2398)   1.26     2399)  0.630     2400)   1.26    
    2401) -0.370     2402)  0.414     2403)  0.588     2404)  0.197     2405)  0.586    
    2406) -0.287     2407)  0.145     2408)   1.08     2409)  0.502     2410)  0.186    
    2411)  0.676     2412)  0.280     2413)  0.855     2414)  0.348     2415)  0.404    
    2416)  0.802     2417)  8.200E-02 2418)  0.180     2419)  0.835     2420)   1.28    
    2421)  5.300E-02 2422)  0.625     2423)  0.947     2424)  0.538     2425)  0.730    
    2426)  0.712     2427)  0.420     2428)  0.762     2429)   1.17     2430)  0.738    
    2431)  0.177     2432)  0.671     2433)  0.970     2434)  0.488     2435) -9.000E-02
    2436)  9.400E-02 2437)   1.36     2438)  7.800E-02 2439)  0.343     2440)  0.253    
    2441)  0.827     2442)   1.74     2443)  0.287     2444)  0.412     2445)  0.247    
    2446)  0.327     2447)   1.54     2448)  0.352     2449)  8.000E-02 2450)  4.400E-02
    2451)  0.378     2452)   2.37     2453)  0.308     2454)  0.188     2455)  0.225    
    2456)  0.447     2457)   1.51     2458)  0.263     2459)  0.228     2460)  0.361    
    2461)  0.179     2462)  0.880     2463)  0.323     2464)  0.842     2465)  0.158    
    2466)  0.618     2467)  0.224     2468)   2.38     2469)  1.600E-02 2470)   2.26    
    2471)  0.582     2472)  0.548     2473)   3.54     2474)   1.00     2475)  0.514    
    2476)  0.207     2477)  0.192     2478)  6.929E-04 2479)  0.399     2480)  0.999    
    2481)   0.00     2482)  1.098E-03 2483)  0.999     2484)  4.054E-04 2485)   1.00    
    2486)  1.385E-03 2487)  0.999     2488)  0.333     2489)  0.286     2490)  2.877E-04
    2491)   1.00     2492)  0.667     2493)  1.608E-03 2494)  0.998     2495)  9.159E-04
    2496)  0.999     2497)  0.182     2498)  5.107E-04 2499)  0.999     2500)  0.154    
    2501)  2.231E-04 2502)   1.00     2503)  0.143     2504)  1.790E-03 2505)  0.998    
    2506)  0.111     2507)  1.823E-04 2508)  0.486     2509)   1.00     2510)  0.105    
    2511)  1.944E-03 2512)  0.998     2513)  0.167     2514)  1.252E-03 2515)  0.999    
    2516)  8.469E-04 2517)  0.999     2518)  5.595E-04 2519)  0.999     2520)  9.091E-02
    2521)  3.364E-04 2522)   1.00     2523)  8.333E-02 2524)  1.541E-04 2525)   1.00    
    2526)  2.077E-03 2527)  0.998     2528)  9.804E-04 2529)  0.999     2530)  7.407E-02
    2531)  4.699E-04 2532)   1.00     2533)  6.667E-02 2534)  6.452E-02 2535)  7.143E-02
    2536)  2.195E-03 2537)  0.998     2538)  1.503E-03 2539)  0.998     2540)  8.696E-02
    2541)  6.897E-02 2542)  8.106E-04 2543)  0.999     2544)  5.882E-02 2545)  5.876E-04
    2546)  0.999     2547)  5.263E-02 2548)  4.878E-02 2549)  0.694     2550)  0.273    
    2551)  0.275     2552)  0.306     2553)  4.140E-02 2554)  0.432     2555)  0.959    
    2556)  6.482E-02 2557)  0.935     2558)  2.443E-02 2559)  0.976     2560)  8.109E-02
    2561)  0.919     2562)  1.740E-02 2563)  0.983     2564)  9.351E-02 2565)  0.906    
    2566)  5.436E-02 2567)  0.946     2568)  3.068E-02 2569)  0.969     2570)  1.352E-02
    2571)  0.986     2572)  0.104     2573)  0.896     2574)  1.106E-02 2575)  0.989    
    2576)  0.112     2577)  0.888     2578)  7.357E-02 2579)  0.926     2580)  5.037E-02
    2581)  0.950     2582)  3.356E-02 2583)  0.966     2584)  2.032E-02 2585)  0.980    
    2586)  9.359E-03 2587)  0.991     2588)  0.119     2589)  0.881     2590)  5.808E-02
    2591)  0.942     2592)  2.826E-02 2593)  0.972     2594)  0.576     2595)  0.356    
    2596)  1.500E-02 2597)  1.034E-02 2598)  0.371     2599)  0.990     2600)  1.634E-02
    2601)  0.984     2602)  6.064E-03 2603)  0.994     2604)  2.058E-02 2605)  0.979    
    2606)  4.306E-03 2607)  0.996     2608)  2.385E-02 2609)  0.976     2610)  1.365E-02
    2611)  0.986     2612)  7.633E-03 2613)  0.992     2614)  3.342E-03 2615)  0.997    
    2616)  0.606     2617)  2.600E-02 2618)  5.852E-02 2619)  0.632     2620)  0.941    
    2621)  9.115E-02 2622)  0.909     2623)  3.466E-02 2624)  0.965     2625)  0.114    
    2626)  0.886     2627)  2.472E-02 2628)  0.975     2629)  0.131     2630)  0.869    
    2631)  7.662E-02 2632)  0.923     2633)  4.347E-02 2634)  0.957     2635)  1.923E-02
    2636)  0.981     2637)  0.144     2638)  0.856     2639)  1.574E-02 2640)  0.639    
    2641)  0.984     2642)  0.156     2643)  0.844     2644)  0.103     2645)  0.897    
    2646)  7.106E-02 2647)  0.929     2648)  4.752E-02 2649)  0.952     2650)  2.885E-02
    2651)  0.971     2652)  1.332E-02 2653)  0.987     2654)  0.165     2655)  0.835    
    2656)  8.179E-02 2657)  0.918     2658)  4.007E-02 2659)  0.960     2660)  0.174    
    2661)  0.826     2662)  0.123     2663)  0.877     2664)  6.812E-02 2665)  0.932    
    2666)  4.985E-02 2667)  0.950     2668)  0.182     2669)  0.818     2670)  9.945E-02
    2671)  0.901     2672)  6.250E-02 2673)  4.651E-02 2674)  4.255E-02 2675)  4.545E-02
    2676)  0.188     2677)  0.812     2678)  0.138     2679)  0.862     2680)  0.107    
    2681)  0.893     2682)  5.556E-02 2683)  8.425E-02 2684)  0.916     2685)  4.762E-02
    2686)  6.630E-02 2687)  0.934     2688)  4.167E-02 2689)  5.137E-02 2690)  0.949    
    2691)  3.774E-02 2692)  0.194     2693)  0.806     2694)  7.334E-02 2695)  0.927    
    2696)  3.704E-02 2697)  3.390E-02 2698)  3.571E-02 2699)  0.200     2700)  0.800    
    2701)  0.150     2702)  0.850     2703)  5.714E-02 2704)  0.120     2705)  0.880    
    2706)  4.444E-02 2707)  9.746E-02 2708)  0.903     2709)  7.977E-02 2710)  0.920    
    2711)  3.226E-02 2712)  6.505E-02 2713)  0.935     2714)  2.899E-02 2715)  3.030E-02
    2716)  0.205     2717)  0.795     2718)  7.692E-02 2719)  5.405E-02 2720)  0.125    
    2721)  0.875     2722)  3.448E-02 2723)  8.568E-02 2724)  0.914     2725)  2.985E-02
    2726)  2.667E-02 2727)  2.778E-02 2728)  0.210     2729)  0.790     2730)  0.161    
    2731)  0.839     2732)  3.846E-02 2733)  0.109     2734)  0.891     2735)  2.469E-02
    2736)  2.564E-02

 35 TYPES OF PHOTOLYSIS REACTIONS:

      FILE = NO2-06        LABEL =NO2 + HV = NO + O                                                               
                       FIRST REACTION = 1     
                     (WL,ABS*QY) =  (0.295 1.199E-19)  (0.300 1.433E-19)  (0.305 1.745E-19)  (0.310 2.020E-19)  (0.315 2.345E-19)
                                    (0.320 2.700E-19)  (0.325 3.020E-19)  (0.330 3.375E-19)  (0.335 3.800E-19)  (0.340 4.100E-19)
                                    (0.345 4.400E-19)  (0.350 4.913E-19)  (0.360 5.420E-19)  (0.370 5.835E-19)  (0.380 6.063E-19)
                                    (0.390 6.100E-19)  (0.400 2.563E-19)  (0.410 4.211E-20)  (0.420 2.373E-21)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = NO3NO-06      LABEL =NO3 + HV = NO + O2                                                              
                       FIRST REACTION = 16    
                     (WL,ABS*QY) =  (0.295 0.000E+00)  (0.300 0.000E+00)  (0.305 0.000E+00)  (0.310 0.000E+00)  (0.315 0.000E+00)
                                    (0.320 0.000E+00)  (0.325 0.000E+00)  (0.330 0.000E+00)  (0.335 0.000E+00)  (0.340 0.000E+00)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 7.901E-20)

      FILE = NO3NO2-6      LABEL =NO3 + HV = NO2 + O                                                              
                       FIRST REACTION = 17    
                     (WL,ABS*QY) =  (0.295 0.000E+00)  (0.300 0.000E+00)  (0.305 0.000E+00)  (0.310 0.000E+00)  (0.315 0.000E+00)
                                    (0.320 0.000E+00)  (0.325 0.000E+00)  (0.330 0.000E+00)  (0.335 0.000E+00)  (0.340 0.000E+00)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 1.050E-20)  (0.410 5.400E-20)  (0.420 1.175E-19)  (0.430 1.865E-19)
                                    (0.440 2.440E-19)  (0.450 3.795E-19)  (0.460 5.245E-19)  (0.470 7.320E-19)  (0.480 8.950E-19)
                                    (0.490 1.172E-18)  (0.500 5.652E-19)

      FILE = O3O1D-06      LABEL =O3 + HV = O1D + O2                                                              
                       FIRST REACTION = 18    
                     (WL,ABS*QY) =  (0.295 5.133E-19)  (0.300 2.589E-19)  (0.305 1.196E-19)  (0.310 2.552E-20)  (0.315 8.358E-21)
                                    (0.320 2.535E-21)  (0.325 8.010E-22)  (0.330 4.068E-22)  (0.335 1.864E-22)  (0.340 8.956E-23)
                                    (0.345 4.340E-23)  (0.350 1.295E-23)  (0.360 2.868E-24)  (0.370 7.144E-25)  (0.380 3.312E-25)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = O3O3P-06      LABEL =O3 + HV = O3P + O2                                                              
                       FIRST REACTION = 19    
                     (WL,ABS*QY) =  (0.295 5.575E-20)  (0.300 2.877E-20)  (0.305 2.892E-20)  (0.310 5.057E-20)  (0.315 3.062E-20)
                                    (0.320 1.764E-20)  (0.325 9.062E-21)  (0.330 4.678E-21)  (0.335 2.144E-21)  (0.340 1.030E-21)
                                    (0.345 4.991E-22)  (0.350 1.489E-22)  (0.360 3.298E-23)  (0.370 8.216E-24)  (0.380 5.479E-24)
                                    (0.390 9.325E-24)  (0.400 1.723E-23)  (0.410 3.103E-23)  (0.420 6.005E-23)  (0.430 9.578E-23)
                                    (0.440 1.630E-22)  (0.450 2.467E-22)  (0.460 3.775E-22)  (0.470 5.315E-22)  (0.480 7.990E-22)
                                    (0.490 9.605E-22)  (0.500 1.309E-21)

      FILE = HONO-06       LABEL =HONO + HV = HO. + NO                                                            
                       FIRST REACTION = 23    
                     (WL,ABS*QY) =  (0.295 4.716E-21)  (0.300 8.655E-21)  (0.305 1.441E-20)  (0.310 2.423E-20)  (0.315 3.220E-20)
                                    (0.320 4.963E-20)  (0.325 6.016E-20)  (0.330 9.115E-20)  (0.335 8.353E-20)  (0.340 1.478E-19)
                                    (0.345 7.675E-20)  (0.350 1.680E-19)  (0.360 1.696E-19)  (0.370 5.670E-20)  (0.380 1.028E-19)
                                    (0.390 4.414E-21)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = HNO3          LABEL =HNO3 + HV = products                                                            
                       FIRST REACTION = 28    
                     (WL,ABS*QY) =  (0.295 3.375E-21)  (0.300 2.065E-21)  (0.305 1.155E-21)  (0.310 6.100E-22)  (0.315 3.050E-22)
                                    (0.320 1.475E-22)  (0.325 6.900E-23)  (0.330 3.250E-23)  (0.335 1.600E-23)  (0.340 8.000E-24)
                                    (0.345 5.000E-24)  (0.350 1.000E-24)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = HNO4-06       LABEL =HO2NO2 + HV = PRODUCTS                                                          
                       FIRST REACTION = 34    
                     (WL,ABS*QY) =  (0.295 2.012E-20)  (0.300 1.226E-20)  (0.305 7.420E-21)  (0.310 4.449E-21)  (0.315 2.778E-21)
                                    (0.320 1.764E-21)  (0.325 1.149E-21)  (0.330 7.549E-22)  (0.335 4.783E-22)  (0.340 3.076E-22)
                                    (0.345 2.060E-22)  (0.350 8.580E-23)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = H2O2          LABEL =H2O2 + HV = 2 OH                                                                
                       FIRST REACTION = 41    
                     (WL,ABS*QY) =  (0.295 7.900E-21)  (0.300 5.950E-21)  (0.305 4.500E-21)  (0.310 3.400E-21)  (0.315 2.550E-21)
                                    (0.320 1.900E-21)  (0.325 1.450E-21)  (0.330 1.150E-21)  (0.335 8.500E-22)  (0.340 6.000E-22)
                                    (0.345 4.500E-22)  (0.350 1.000E-22)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = PAN           LABEL =PAN + HV = #.6 {MECO3 + NO2} + #.4 {MEO2 + CO2 + NO3}                           
                       FIRST REACTION = BR20  
                     (WL,ABS*QY) =  (0.295 2.582E-21)  (0.300 1.473E-21)  (0.305 8.714E-22)  (0.310 5.283E-22)  (0.315 3.206E-22)
                                    (0.320 1.951E-22)  (0.325 1.253E-22)  (0.330 8.221E-23)  (0.335 5.598E-23)  (0.340 3.403E-23)
                                    (0.345 2.247E-23)  (0.350 3.300E-24)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = HCHOR-06      LABEL =HCHO + HV = HCO. + H.                                                           
                       FIRST REACTION = BP01  
                     (WL,ABS*QY) =  (0.295 1.815E-20)  (0.300 2.347E-20)  (0.305 1.993E-20)  (0.310 1.675E-20)  (0.315 2.142E-20)
                                    (0.320 6.148E-21)  (0.325 1.751E-20)  (0.330 3.519E-21)  (0.335 6.380E-22)  (0.340 0.000E+00)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = HCHOM-06      LABEL =HCHO + HV = H2 + CO                                                             
                       FIRST REACTION = BP02  
                     (WL,ABS*QY) =  (0.295 8.994E-21)  (0.300 7.686E-21)  (0.305 7.979E-21)  (0.310 6.938E-21)  (0.315 1.030E-20)
                                    (0.320 3.590E-21)  (0.325 1.784E-20)  (0.330 3.899E-21)  (0.335 1.068E-20)  (0.340 6.273E-21)
                                    (0.345 3.477E-22)  (0.350 8.333E-22)  (0.360 2.505E-25)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CCHO_R        LABEL =CCHO + HV = CH3 + CHO                                                           
                       FIRST REACTION = BP09  
                     (WL,ABS*QY) =  (0.295 1.957E-20)  (0.300 1.511E-20)  (0.305 1.087E-20)  (0.310 5.940E-21)  (0.315 2.704E-21)
                                    (0.320 9.630E-22)  (0.325 1.905E-22)  (0.330 0.000E+00)  (0.335 0.000E+00)  (0.340 0.000E+00)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = C2CHO         LABEL =C2CHO + HV = C2H5. + CHO.                                                       
                       FIRST REACTION = BP12  
                     (WL,ABS*QY) =  (0.295 4.516E-20)  (0.300 4.037E-20)  (0.305 2.586E-20)  (0.310 1.738E-20)  (0.315 8.916E-21)
                                    (0.320 5.052E-21)  (0.325 2.141E-21)  (0.330 7.558E-22)  (0.335 2.154E-22)  (0.340 1.819E-23)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = ACET-06       LABEL =CH3-CO-CH3 + HV = Radical products                                              
                       FIRST REACTION = BP15  
                     (WL,ABS*QY) =  (0.295 9.281E-21)  (0.300 5.487E-21)  (0.305 2.106E-21)  (0.310 6.911E-22)  (0.315 2.217E-22)
                                    (0.320 6.786E-23)  (0.325 1.736E-23)  (0.330 3.841E-24)  (0.335 9.004E-25)  (0.340 2.283E-25)
                                    (0.345 4.355E-26)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = MEK-06        LABEL =MEK absorption cross sections                                                   
                       FIRST REACTION = BP17  
                     (WL,ABS*QY) =  (0.295 3.708E-20)  (0.300 2.799E-20)  (0.305 1.949E-20)  (0.310 1.216E-20)  (0.315 6.619E-21)
                                    (0.320 3.090E-21)  (0.325 1.206E-21)  (0.330 3.950E-22)  (0.335 1.050E-22)  (0.340 3.200E-23)
                                    (0.345 3.000E-24)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = COOH          LABEL =CH3OOH + HV = PRODUCTS                                                          
                       FIRST REACTION = BP23  
                     (WL,ABS*QY) =  (0.295 4.820E-21)  (0.300 3.630E-21)  (0.305 2.760E-21)  (0.310 2.105E-21)  (0.315 1.595E-21)
                                    (0.320 1.210E-21)  (0.325 9.200E-22)  (0.330 7.000E-22)  (0.335 5.400E-22)  (0.340 4.100E-22)
                                    (0.345 3.100E-22)  (0.350 2.125E-22)  (0.360 1.000E-22)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = GLY-07R       LABEL =HCOCHO + HV = HCO. + HCO.                                                       
                       FIRST REACTION = BP30  
                     (WL,ABS*QY) =  (0.295 1.335E-20)  (0.300 1.250E-20)  (0.305 1.118E-20)  (0.310 1.265E-20)  (0.315 1.010E-20)
                                    (0.320 8.504E-21)  (0.325 9.299E-21)  (0.330 6.068E-21)  (0.335 3.862E-21)  (0.340 3.774E-21)
                                    (0.345 4.134E-21)  (0.350 2.693E-21)  (0.360 2.940E-21)  (0.370 3.216E-21)  (0.380 3.382E-21)
                                    (0.390 3.340E-21)  (0.400 2.779E-21)  (0.410 2.823E-21)  (0.420 2.334E-21)  (0.430 1.226E-21)
                                    (0.440 8.665E-22)  (0.450 8.878E-22)  (0.460 8.288E-24)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = GLY-07M       LABEL =HCOCHO + HV = HCHO + H2                                                         
                       FIRST REACTION = BP31  
                     (WL,ABS*QY) =  (0.295 2.638E-20)  (0.300 2.399E-20)  (0.305 2.055E-20)  (0.310 1.964E-20)  (0.315 1.305E-20)
                                    (0.320 9.058E-21)  (0.325 8.183E-21)  (0.330 4.480E-21)  (0.335 2.331E-21)  (0.340 1.845E-21)
                                    (0.345 1.649E-21)  (0.350 8.151E-22)  (0.360 6.765E-22)  (0.370 5.717E-22)  (0.380 4.633E-22)
                                    (0.390 3.571E-22)  (0.400 2.288E-22)  (0.410 1.771E-22)  (0.420 1.127E-22)  (0.430 2.858E-23)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = MGLY-06       LABEL =CH3COCHO + HV = CH3CO. + HCO.                                                   
                       FIRST REACTION = BP34  
                     (WL,ABS*QY) =  (0.295 3.851E-20)  (0.300 3.316E-20)  (0.305 2.610E-20)  (0.310 1.880E-20)  (0.315 1.448E-20)
                                    (0.320 9.991E-21)  (0.325 6.127E-21)  (0.330 4.090E-21)  (0.335 2.696E-21)  (0.340 1.897E-21)
                                    (0.345 1.897E-21)  (0.350 2.530E-21)  (0.360 3.807E-21)  (0.370 5.165E-21)  (0.380 6.304E-21)
                                    (0.390 6.555E-21)  (0.400 6.689E-21)  (0.410 6.621E-21)  (0.420 5.153E-21)  (0.430 3.906E-21)
                                    (0.440 2.747E-21)  (0.450 8.781E-22)  (0.460 6.265E-23)  (0.470 5.915E-24)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = BACL-07       LABEL =CH3COCOCH3 + HV = 2 CH3CO.                                                      
                       FIRST REACTION = BP37  
                     (WL,ABS*QY) =  (0.295 2.190E-20)  (0.300 1.665E-20)  (0.305 1.135E-20)  (0.310 7.339E-21)  (0.315 5.499E-21)
                                    (0.320 4.919E-21)  (0.325 3.432E-21)  (0.330 3.557E-21)  (0.335 4.286E-21)  (0.340 5.093E-21)
                                    (0.345 5.906E-21)  (0.350 6.947E-21)  (0.360 9.491E-21)  (0.370 1.286E-20)  (0.380 1.528E-20)
                                    (0.390 1.467E-20)  (0.400 1.268E-20)  (0.410 1.151E-20)  (0.420 8.061E-21)  (0.430 5.528E-21)
                                    (0.440 3.909E-21)  (0.450 1.096E-21)  (0.460 1.111E-22)  (0.470 5.670E-24)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = BALD-06       LABEL =Benzaldehyde absorption cross sections recommended by Calvert et al (2002)      
                       FIRST REACTION = BP44  
                     (WL,ABS*QY) =  (0.295 1.338E-19)  (0.300 6.810E-20)  (0.305 7.310E-20)  (0.310 7.350E-20)  (0.315 9.010E-20)
                                    (0.320 8.810E-20)  (0.325 9.900E-20)  (0.330 8.780E-20)  (0.335 9.790E-20)  (0.340 8.580E-20)
                                    (0.345 6.510E-20)  (0.350 6.400E-20)  (0.360 2.825E-20)  (0.370 1.516E-20)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = AFG1          LABEL =Photolysis of lumped photoreactive aromatic unsaturated dicarbonyl product.     
                       FIRST REACTION = BP48  
                     (WL,ABS*QY) =  (0.295 1.650E-18)  (0.300 1.445E-18)  (0.305 1.215E-18)  (0.310 1.015E-18)  (0.315 8.595E-19)
                                    (0.320 7.205E-19)  (0.325 5.965E-19)  (0.330 4.995E-19)  (0.335 4.180E-19)  (0.340 3.525E-19)
                                    (0.345 3.015E-19)  (0.350 2.432E-19)  (0.360 1.835E-19)  (0.370 1.393E-19)  (0.380 1.036E-19)
                                    (0.390 7.650E-20)  (0.400 5.635E-20)  (0.410 3.395E-20)  (0.420 2.445E-20)  (0.430 1.398E-20)
                                    (0.440 4.902E-21)  (0.450 1.828E-21)  (0.460 4.563E-22)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = MACR-06       LABEL =Methacrolein total photolysis                                                   
                       FIRST REACTION = BP58  
                     (WL,ABS*QY) =  (0.295 1.264E-21)  (0.300 1.204E-21)  (0.305 1.107E-21)  (0.310 9.918E-22)  (0.315 8.359E-22)
                                    (0.320 6.669E-22)  (0.325 5.470E-22)  (0.330 4.222E-22)  (0.335 2.796E-22)  (0.340 2.208E-22)
                                    (0.345 1.479E-22)  (0.350 6.933E-23)  (0.360 2.524E-23)  (0.370 9.467E-24)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = MVK-06        LABEL =MVK total photolysis                                                            
                       FIRST REACTION = BP63  
                     (WL,ABS*QY) =  (0.295 5.184E-22)  (0.300 4.896E-22)  (0.305 4.330E-22)  (0.310 3.814E-22)  (0.315 3.223E-22)
                                    (0.320 2.525E-22)  (0.325 2.137E-22)  (0.330 1.508E-22)  (0.335 1.266E-22)  (0.340 6.351E-23)
                                    (0.345 5.835E-23)  (0.350 4.101E-23)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = IC3ONO2       LABEL =! I-C3H7ONO2 + HV = PRODUCTS                                                    
                       FIRST REACTION = BP71  
                     (WL,ABS*QY) =  (0.295 1.005E-20)  (0.300 6.650E-21)  (0.305 4.200E-21)  (0.310 2.550E-21)  (0.315 1.500E-21)
                                    (0.320 8.550E-22)  (0.325 4.900E-22)  (0.330 1.850E-22)  (0.335 0.000E+00)  (0.340 0.000E+00)
                                    (0.345 0.000E+00)  (0.350 0.000E+00)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CL2           LABEL =CL2 + HV = 2 CL 	                                                               
                       FIRST REACTION = CI01  
                     (WL,ABS*QY) =  (0.295 1.047E-19)  (0.300 1.355E-19)  (0.305 1.685E-19)  (0.310 1.980E-19)  (0.315 2.240E-19)
                                    (0.320 2.415E-19)  (0.325 2.505E-19)  (0.330 2.500E-19)  (0.335 2.400E-19)  (0.340 2.232E-19)
                                    (0.345 1.998E-19)  (0.350 1.600E-19)  (0.360 1.080E-19)  (0.370 6.700E-20)  (0.380 3.950E-20)
                                    (0.390 2.350E-20)  (0.400 1.550E-20)  (0.410 1.130E-20)  (0.420 8.450E-21)  (0.430 6.350E-21)
                                    (0.440 4.600E-21)  (0.450 3.200E-21)  (0.460 2.100E-21)  (0.470 8.000E-22)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLNO-06       LABEL =CLNO absorption cross sections recommended by IUPAC (2005)                      
                       FIRST REACTION = CI03  
                     (WL,ABS*QY) =  (0.295 9.839E-20)  (0.300 1.034E-19)  (0.305 1.105E-19)  (0.310 1.195E-19)  (0.315 1.286E-19)
                                    (0.320 1.373E-19)  (0.325 1.448E-19)  (0.330 1.495E-19)  (0.335 1.527E-19)  (0.340 1.520E-19)
                                    (0.345 1.490E-19)  (0.350 1.365E-19)  (0.360 1.197E-19)  (0.370 9.940E-20)  (0.380 7.840E-20)
                                    (0.390 5.983E-20)  (0.400 4.440E-20)  (0.410 3.370E-20)  (0.420 2.500E-20)  (0.430 2.202E-20)
                                    (0.440 2.052E-20)  (0.450 1.850E-20)  (0.460 2.237E-20)  (0.470 2.562E-20)  (0.480 2.315E-20)
                                    (0.490 1.782E-20)  (0.500 4.050E-22)

      FILE = CLONO         LABEL =CLONO + HV = CL + NO2	                                                          
                       FIRST REACTION = CI06  
                     (WL,ABS*QY) =  (0.295 1.355E-18)  (0.300 1.215E-18)  (0.305 1.095E-18)  (0.310 1.016E-18)  (0.315 8.920E-19)
                                    (0.320 7.785E-19)  (0.325 6.705E-19)  (0.330 5.820E-19)  (0.335 5.070E-19)  (0.340 3.970E-19)
                                    (0.345 3.130E-19)  (0.350 2.220E-19)  (0.360 1.192E-19)  (0.370 6.725E-20)  (0.380 3.225E-20)
                                    (0.390 1.450E-20)  (0.400 1.500E-21)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLNO2         LABEL =CLNO2 + HV = CL + NO2	                                                          
                       FIRST REACTION = CI07  
                     (WL,ABS*QY) =  (0.295 1.550E-19)  (0.300 1.420E-19)  (0.305 1.280E-19)  (0.310 1.129E-19)  (0.315 9.678E-20)
                                    (0.320 8.113E-20)  (0.325 6.598E-20)  (0.330 5.265E-20)  (0.335 4.115E-20)  (0.340 3.165E-20)
                                    (0.345 2.415E-20)  (0.350 1.595E-20)  (0.360 9.200E-21)  (0.370 3.450E-21)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLONO2-1      LABEL =CLONO2 + hv = CLO. + NO2                                                        
                       FIRST REACTION = CI14  
                     (WL,ABS*QY) =  (0.295 1.502E-20)  (0.300 1.074E-20)  (0.305 7.568E-21)  (0.310 5.083E-21)  (0.315 3.300E-21)
                                    (0.320 2.160E-21)  (0.325 1.421E-21)  (0.330 9.460E-22)  (0.335 6.381E-22)  (0.340 4.336E-22)
                                    (0.345 2.899E-22)  (0.350 1.383E-22)  (0.360 1.450E-23)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLONO2-2      LABEL =CLONO2 + hv = CL. + NO3                                                         
                       FIRST REACTION = CI15  
                     (WL,ABS*QY) =  (0.295 2.253E-20)  (0.300 1.611E-20)  (0.305 1.163E-20)  (0.310 8.617E-21)  (0.315 6.555E-21)
                                    (0.320 5.060E-21)  (0.325 3.974E-21)  (0.330 3.219E-21)  (0.335 2.707E-21)  (0.340 2.366E-21)
                                    (0.345 2.145E-21)  (0.350 1.974E-21)  (0.360 1.783E-21)  (0.370 1.405E-21)  (0.380 1.215E-21)
                                    (0.390 7.675E-22)  (0.400 5.400E-22)  (0.410 3.700E-22)  (0.420 2.425E-22)  (0.430 4.750E-23)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = HOCL-06       LABEL =HOCL absorption cross sections recommended by IUPAC (2005)                      
                       FIRST REACTION = CI19  
                     (WL,ABS*QY) =  (0.295 5.816E-20)  (0.300 6.078E-20)  (0.305 6.070E-20)  (0.310 5.784E-20)  (0.315 5.265E-20)
                                    (0.320 4.596E-20)  (0.325 3.867E-20)  (0.330 3.147E-20)  (0.335 2.505E-20)  (0.340 1.985E-20)
                                    (0.345 1.589E-20)  (0.350 1.219E-20)  (0.360 9.698E-21)  (0.370 8.018E-21)  (0.380 6.011E-21)
                                    (0.390 3.861E-21)  (0.400 2.108E-21)  (0.410 9.930E-22)  (0.420 6.300E-23)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLCCHO        LABEL =Chloroacetaldehyde absorption cross sections                                    
                       FIRST REACTION = CP19  
                     (WL,ABS*QY) =  (0.295 5.620E-20)  (0.300 5.150E-20)  (0.305 5.356E-20)  (0.310 4.093E-20)  (0.315 3.906E-20)
                                    (0.320 2.834E-20)  (0.325 1.828E-20)  (0.330 1.466E-20)  (0.335 5.964E-21)  (0.340 2.361E-21)
                                    (0.345 9.637E-22)  (0.350 3.370E-22)  (0.360 0.000E+00)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)

      FILE = CLACET        LABEL =Chloroacetone absorption cross sections                                         
                       FIRST REACTION = CP22  
                     (WL,ABS*QY) =  (0.295 9.713E-20)  (0.300 8.889E-20)  (0.305 7.495E-20)  (0.310 6.057E-20)  (0.315 4.406E-20)
                                    (0.320 2.939E-20)  (0.325 1.842E-20)  (0.330 9.244E-21)  (0.335 4.042E-21)  (0.340 1.661E-21)
                                    (0.345 7.154E-22)  (0.350 2.283E-22)  (0.360 3.200E-23)  (0.370 0.000E+00)  (0.380 0.000E+00)
                                    (0.390 0.000E+00)  (0.400 0.000E+00)  (0.410 0.000E+00)  (0.420 0.000E+00)  (0.430 0.000E+00)
                                    (0.440 0.000E+00)  (0.450 0.000E+00)  (0.460 0.000E+00)  (0.470 0.000E+00)  (0.480 0.000E+00)
                                    (0.490 0.000E+00)  (0.500 0.000E+00)
0DEFAULT TEMPERATURE =  298.12
 NO. REACTIONS =  805
 NO. SPECIES =  605
 NO. VARIABLE COEFFICIENTS =   20
 NO. CONSTANT COEFFICIENTS =  736
 KINETIC PARM BUFFER LENGTH = 1639
 PRODUCT ARRAY LENGTH =  41640
 
 MODEL = 'saprc14_rev1    ' STORED.

 **********************************************************************************************************************************

                                              E N D   P R E P   O U T P U T

 **********************************************************************************************************************************
