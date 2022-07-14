1**********************************************************************************************************************************

                              S A P R C   M O D E L   P R E P A R A T I O N   P R O G R A M

                                                     26-Sep-18 13:23 

                             TEMP 300.                                                       

 **********************************************************************************************************************************
0PARAMETERS SPECIFIED:
 

 REACTIONS INPUT:                                                              A        EA      B
                                                                               (IN PPM-MIN UNITS)
 1      NO2 + HV = NO + O3P                                                     PHOT. = NO2-06            3     8     0 =     2    45                            
 2      O3P + O2 + M = O3 + M                                              2.040E-05   0.00  -4.600      45     4     5 =     1     5                            
 3      O3P + O3 = #2 O2                                                   1.174E+04   4.09  -1.000      45     1     0 = -4001     4                            
 4      O3P + NO = NO2                                                    FALLOFF F= 0.600, N= 1.000     45     2     0 =     3                                  
                                                                      K0:  3.232E-03   0.00  -3.500
                                                                      KI:  4.404E+04   0.00  -1.000
 5      O3P + NO2 = NO + O2                                                8.073E+03  -0.37  -1.000      45     3     0 =     2     4                            
 6      O3P + NO2 = NO3                                                   FALLOFF F= 0.600, N= 1.000     45     3     0 =    10                                  
                                                                      K0:  8.978E-03   0.00  -3.800
                                                                      KI:  3.229E+04   0.00  -1.700
 7      O3 + NO = NO2 + O2                                                 4.404E+03   2.98  -1.000       1     2     0 =     3     4                            
 8      O3 + NO2 = O2 + NO3                                                2.055E+02   4.91  -1.000       1     3     0 =     4    10                            
 9      NO + NO3 = #2 NO2                                                  2.642E+04  -0.22  -1.000       2    10     0 = -4001     3                            
 10     NO + NO + O2 = #2 NO2                                              0.000E+00  -1.05  -2.000       2     2     4 = -4001     3                            
 11     NO2 + NO3 = N2O5                                                  FALLOFF F= 0.350, N= 1.330      3    10     0 =    11                                  
                                                                      K0:  1.293E-01   0.00  -6.100
                                                                      KI:  2.789E+03   0.00  -0.800
 12     N2O5 = NO2 + NO3                                                  FALLOFF F= 0.350, N= 1.330     11     0     0 =     3    10                            
                                                                      K0:  1.908E+12  21.86  -4.500
                                                                      KI:  5.820E+16  22.02   0.100
 13     N2O5 + H2O = #2 HNO3                                               3.670E-07   0.00  -1.000      11     6     0 = -4001    12                            
 14     N2O5 + H2O + H2O = #2 HNO3 + H2O                                   0.000E+00   0.00  -2.000      11     6     6 = -4001    12     6                      
 15     NO2 + NO3 = NO + NO2 + O2                                          6.606E+01   2.50  -1.000       3    10     0 =     2     3     4                      
 16     NO3 + HV = NO + O2                                                      PHOT. = NO3NO-06         10     8     0 =     2     4                            
 17     NO3 + HV = NO2 + O3P                                                    PHOT. = NO3NO2-6         10     8     0 =     3    45                            
 18     O3 + HV = O1D + O2                                                      PHOT. = O3O1D-06          1     8     0 =    46     4                            
 19     O3 + HV = O3P + O2                                                      PHOT. = O3O3P-06          1     8     0 =    45     4                            
 20     O1D + H2O = #2 OH                                                  2.393E+05  -0.12  -1.000      46     6     0 = -4001     9                            
 21     O1D + M = O3P + M                                                  3.494E+04  -0.19  -1.000      46     5     0 =    45     5                            
 22     OH + NO = HONO                                                    FALLOFF F= 0.600, N= 1.000      9     2     0 =    13                                  
                                                                      K0:  2.514E-02   0.00  -4.600
                                                                      KI:  5.284E+04   0.00  -1.100
 23     HONO + HV = OH + NO                                                     PHOT. = HONO-06          13     8     0 =     9     2                            
 24     OH + HONO = H2O + NO2                                              3.670E+03  -0.52  -1.000       9    13     0 =     6     3                            
 25     OH + NO2 = HNO3                                                   FALLOFF F= 0.600, N= 1.000      9     3     0 =    12                                  
                                                                      K0:  6.464E-02   0.00  -5.000
                                                                      KI:  4.110E+04   0.00  -1.000
 26     OH + NO3 = HO2 + NO2                                               2.936E+04   0.00  -1.000       9    10     0 =    82     3                            
 27     OH + HNO3 = H2O + NO3                                              3.523E+01  -0.91  -1.000       9    12     0 =     6    10                            
 28     HNO3 + HV = OH + NO2                                                    PHOT. = HNO3             12     8     0 =     9     3                            
 29     OH + CO = HO2 + CO2                                                2.114E+02   0.00  -1.000       9    83     0 =    82    84                            
 30     OH + O3 = HO2 + O2                                                 2.495E+03   1.87  -1.000       9     1     0 =    82     4                            
 31     HO2 + NO = OH + NO2                                                5.284E+03  -0.54  -1.000      82     2     0 =     9     3                            
 32     HO2 + NO2 = HNO4                                                  FALLOFF F= 0.600, N= 1.000     82     3     0 =    14                                  
                                                                      K0:  7.182E-03   0.00  -5.400
                                                                      KI:  4.257E+03   0.00  -2.100
 33     HNO4 = HO2 + NO2                                                  FALLOFF F= 0.600, N= 1.000     14     0     0 =    82     3                            
                                                                      K0:  5.461E+10  21.16  -3.400
                                                                      KI:  3.252E+17  22.20  -2.300
 34     HNO4 + HV = #.61 HO2 + #.61 NO2 + #.39 OH + #.39 NO3                    PHOT. = HNO4-06          14     8     0 = -4002    82 -4002     3                
                                                                                                                    -4003     9 -4003    10
 35     HNO4 + OH = H2O + NO2 + O2                                         1.908E+03  -0.75  -1.000      14     9     0 =     6     3     4                      
 36     HO2 + O3 = OH + #2 O2                                              2.980E-01  -1.38   3.570      82     1     0 =     9 -4001     4                      
 37     HO2 + HO2 = HO2H + O2                                              3.229E+02  -1.19  -1.000      82    82     0 =    85     4                            
 38     HO2 + HO2 + H2O = HO2H + O2 + H2O                                  1.106E-05  -5.56  -2.000      82    82     6 =    85     4     6                      
 39     NO3 + HO2 = #.8 OH + #.8 NO2 + #.8 O2 + #.2 HNO3 + #.2 O2          5.872E+03   0.00  -1.000      10    82     0 = -4004     9 -4004     3                
                                                                                                                    -4004     4 -4005    12
                                                                                                                    -4005     4
 40     NO3 + NO3 = #2 NO2 + O2                                            1.248E+03   4.87  -1.000      10    10     0 = -4001     3     4                      
 41     HO2H + HV = #2 OH                                                       PHOT. = H2O2             85     8     0 = -4001     9                            
 42     HO2H + OH = HO2 + H2O                                              2.642E+03   0.00  -1.000      85     9     0 =    82     6                            
 43     OH + HO2 = H2O + O2                                                7.046E+04  -0.50  -1.000       9    82     0 =     6     4                            
 44     OH + SO2 = HO2 + SULF                                             FALLOFF F= 0.600, N= 1.000      9    15     0 =    82    16                            
                                                                      K0:  1.185E-02   0.00  -6.300
                                                                      KI:  2.349E+03   0.00  -1.000
 45     OH + H2 = HO2 + H2O                                                1.130E+04   4.17  -1.000       9     7     0 =    82     6                            
 BR01   MEO2 + NO = NO2 + HCHO + HO2                                       3.376E+03  -0.71  -1.000      86     2     0 =     3    87    82                      
 BR02   MEO2 + HO2 = COOH + O2                                             5.079E+02  -1.55  -0.640      86    82     0 =    88     4                            
 BR03   MEO2 + HO2 = HCHO + O2 + H2O                                       4.903E+01  -1.55  -4.530      86    82     0 =    87     4     6                      
 BR04   MEO2 + NO3 = HCHO + HO2 + NO2                                      1.908E+03   0.00  -1.000      86    10     0 =    87    82     3                      
 BR05   MEO2 + MEO2 = MEOH + HCHO + O2                                     9.380E+01  -0.73  -2.800      86    86     0 =    89    87     4                      
 BR06   MEO2 + MEO2 = #2 HCHO + #2 HO2                                     1.086E+03   1.03  -1.000      86    86     0 = -4001    87 -4001    82                
 BR07   RO2C + NO = NO2                                                    3.817E+03  -0.75  -1.000      90     2     0 =     3                                  
 BR08   RO2C + HO2 =                                                       5.578E+02  -1.79  -1.000      90    82     0 =                                        
 BR09   RO2C + NO3 = NO2                                                   3.376E+03   0.00  -1.000      90    10     0 =     3                                  
 BR10   RO2C + MEO2 = #.5 HO2 + #.75 HCHO + #.25 MEOH                      2.936E+02   0.00  -1.000      90    86     0 = -4006    82 -4007    87                
                                                                                                                    -4008    89
 BR11   RO2C + RO2C =                                                      5.138E+01   0.00  -1.000      90    90     0 =                                        
 BR12   RO2XC + NO = XN                                                         SAME K AS BR07           91     2     0 =    18                                  
 BR13   RO2XC + HO2 =                                                           SAME K AS BR08           91    82     0 =                                        
 BR14   RO2XC + NO3 = NO2                                                       SAME K AS BR09           91    10     0 =     3                                  
 BR15   RO2XC + MEO2 = #.5 HO2 + #.75 HCHO + #.25 MEOH                          SAME K AS BR10           91    86     0 = -4006    82 -4007    87                
                                                                                                                    -4008    89
 BR16   RO2XC + RO2C =                                                          SAME K AS BR11           91    90     0 =                                        
 BR17   RO2XC + RO2XC =                                                         SAME K AS BR11           91    91     0 =                                        
 BR18   MECO3 + NO2 = PAN                                                 FALLOFF F= 0.300, N= 1.410     92     3     0 =    20                                  
                                                                      K0:  9.696E+00   0.00  -9.100
                                                                      KI:  1.776E+04   0.00  -1.900
 BR19   PAN = MECO3 + NO2                                                 FALLOFF F= 0.300, N= 1.410     20     0     0 =    92     3                            
                                                                      K0:  7.193E+12  24.05  -1.000
                                                                      KI:  2.400E+18  27.03   0.000
 BR20   PAN + HV = #.6 MECO3 + #.6 NO2 + #.4 MEO2 + &                            
        #.4 CO2 + #.4 NO3                                                       PHOT. = PAN              20     8     0 = -4009    92 -4009     3                
                                                                                                                    -4010    86 -4010    84
                                                                                                                    -4010    10
 BR21   MECO3 + NO = MEO2 + CO2 + NO2                                      1.101E+04  -0.58  -1.000      92     2     0 =    86    84     3                      
 BR22   MECO3 + HO2 = #.44 OH + #.44 MEO2 + #.44 CO2 + &                         
        #.41 CCOOH + #.15 O3 + #.15 CCOOH                                  7.633E+02  -1.95  -1.000      92    82     0 = -4011     9 -4011    86                
                                                                                                                    -4011    84 -4012    93
                                                                                                                    -4013     1 -4013    93
 BR23   MECO3 + NO3 = MEO2 + CO2 + NO2 + O2                                     SAME K AS BR09           92    10     0 =    86    84     3     4                
 BR24   MECO3 + MEO2 = #.1 CCOOH + #.1 HCHO + #.1 O2 + &                         
        #.9 HCHO + #.9 HO2 + #.9 MEO2 + #.9 CO2                            2.936E+03  -0.99  -1.000      92    86     0 = -4014    93 -4014    87                
                                                                                                                    -4014     4 -4015    87
                                                                                                                    -4015    82 -4015    86
                                                                                                                    -4015    84
 BR25   MECO3 + RO2C = MEO2 + CO2                                          6.459E+02  -2.13  -1.000      92    90     0 =    86    84                            
 BR26   MECO3 + RO2XC = MEO2 + CO2                                              SAME K AS BR25           92    91     0 =    86    84                            
 BR27   MECO3 + MECO3 = #2 MEO2 + #2 CO2 + O2                              4.257E+03  -0.99  -1.000      92    92     0 = -4001    86 -4001    84                
                                                                                                                        4
 BR28   RCO3 + NO2 = PAN2                                                  1.776E+04   0.00  -2.070      94     3     0 =    21                                  
 BR29   PAN2 = RCO3 + NO2                                                  4.980E+18  27.70   0.000      21     0     0 =    94     3                            
 BR30   PAN2 + HV = #.6 RCO3 + #.6 NO2 + #.4 RO2C + &                            
        #.4 xHO2 + #.4 yROOH + #.4 xCCHO + #.4 CO2 + &                           
        #.4 N3                                                                  PHOT. = PAN              21     8     0 = -4009    94 -4009     3                
                                                                                                                    -4010    90 -4010    50
                                                                                                                    -4010    76 -4010    60
                                                                                                                    -4010    84 -4010    95
 BR31   RCO3 + NO = NO2 + RO2C + xHO2 + yROOH + xCCHO + &                        
        CO2                                                                9.835E+03  -0.68  -1.000      94     2     0 =     3    90    50    76                
                                                                                                                       60    84
 BR32   RCO3 + HO2 = #.44 OH + #.44 RO2C + #.44 xHO2 + &                         
        #.44 xCCHO + #.44 yROOH + #.44 CO2 + &                                   
        #.41 RCOOH + #.15 O3 + #.15 RCOOH                                       SAME K AS BR22           94    82     0 = -4011     9 -4011    90                
                                                                                                                    -4011    50 -4011    60
                                                                                                                    -4011    76 -4011    84
                                                                                                                    -4012    96 -4013     1
                                                                                                                    -4013    96
 BR33   RCO3 + NO3 = NO2 + RO2C + xHO2 + yROOH + xCCHO + &                       
        CO2 + O2                                                                SAME K AS BR09           94    10     0 =     3    90    50    76                
                                                                                                                       60    84     4
 BR34   RCO3 + MEO2 = HCHO + HO2 + RO2C + xHO2 + xCCHO + &                       
        yROOH + CO2                                                             SAME K AS BR24           94    86     0 =    87    82    90    50                
                                                                                                                       60    76    84
 BR35   RCO3 + RO2C = RO2C + xHO2 + xCCHO + yROOH + CO2                         SAME K AS BR25           94    90     0 =    90    50    60    76                
                                                                                                                       84
 BR36   RCO3 + RO2XC = RO2C + xHO2 + xCCHO + yROOH + CO2                        SAME K AS BR25           94    91     0 =    90    50    60    76                
                                                                                                                       84
 BR37   RCO3 + MECO3 = #2 CO2 + MEO2 + RO2C + xHO2 + &                           
        yROOH + xCCHO + O2                                                      SAME K AS BR27           94    92     0 = -4001    84    86    90                
                                                                                                                       50    76    60     4
 BR38   RCO3 + RCO3 = #2 RO2C + #2 xHO2 + #2 xCCHO + &                           
        #2 yROOH + #2 CO2                                                       SAME K AS BR27           94    94     0 = -4001    90 -4001    50                
                                                                                                                    -4001    60 -4001    76
                                                                                                                    -4001    84
 BR39   BZCO3 + NO2 = PBZN                                                 2.011E+04   0.00  -1.000      97     3     0 =    22                                  
 BR40   PBZN = BZCO3 + NO2                                                 4.740E+18  27.82   0.000      22     0     0 =    97     3                            
 BR41   PBZN + HV = #.6 BZCO3 + #.6 NO2 + #.4 CO2 + &                            
        #.4 BZO + #.4 RO2C + #.4 NO3                                            PHOT. = PAN              22     8     0 = -4009    97 -4009     3                
                                                                                                                    -4010    84 -4010    48
                                                                                                                    -4010    90 -4010    10
 BR42   BZCO3 + NO = NO2 + CO2 + BZO + RO2C                                     SAME K AS BR31           97     2     0 =     3    84    48    90                
 BR43   BZCO3 + HO2 = #.44 OH + #.44 BZO + #.44 RO2C + &                         
        #.44 CO2 + #.41 RCOOH + &                                                
        #.15 O3 + #.44 RCOOH + #2.24 XC                                         SAME K AS BR22           97    82     0 = -4011     9 -4011    48                
                                                                                                                    -4011    90 -4011    84
                                                                                                                    -4012    96 -4013     1
                                                                                                                    -4011    96 -4016    98
 BR44   BZCO3 + NO3 = NO2 + CO2 + BZO + RO2C + O2                               SAME K AS BR09           97    10     0 =     3    84    48    90                
                                                                                                                        4
 BR45   BZCO3 + MEO2 = HCHO + HO2 + RO2C + BZO + CO2                            SAME K AS BR24           97    86     0 =    87    82    90    48                
                                                                                                                       84
 BR46   BZCO3 + RO2C = RO2C + BZO + CO2                                         SAME K AS BR25           97    90     0 =    90    48    84                      
 BR47   BZCO3 + RO2XC = RO2C + BZO + CO2                                        SAME K AS BR25           97    91     0 =    90    48    84                      
 BR48   BZCO3 + MECO3 = #2 CO2 + MEO2 + BZO + RO2C                              SAME K AS BR27           97    92     0 = -4001    84    86    48                
                                                                                                                       90
 BR49   BZCO3 + RCO3 = #2 CO2 + RO2C + xHO2 + yROOH + &                          
        xCCHO + BZO + RO2C                                                      SAME K AS BR27           97    94     0 = -4001    84    90    50                
                                                                                                                       76    60    48    90
 BR50   BZCO3 + BZCO3 = #2 BZO + #2 RO2C + #2 CO2                               SAME K AS BR27           97    97     0 = -4001    48 -4001    90                
                                                                                                                    -4001    84
 BR51   MACO3 + NO2 = MAPAN                                                     SAME K AS BR28           99     3     0 =    23                                  
 BR52   MAPAN = MACO3 + NO2                                                9.600E+17  26.80   0.000      23     0     0 =    99     3                            
 BR53   MAPAN + HV = #.6 MACO3 + #.6 NO2 + #.4 CO2 + &                           
        #.4 HCHO + #.4 MECO3 + #.4 NO3                                          PHOT. = PAN              23     8     0 = -4009    99 -4009     3                
                                                                                                                    -4010    84 -4010    87
                                                                                                                    -4010    92 -4010    10
 BR54   MACO3 + NO = NO2 + CO2 + HCHO + MECO3                                   SAME K AS BR31           99     2     0 =     3    84    87    92                
 BR55   MACO3 + HO2 = #.44 OH + #.44 HCHO + #.44 MECO3 + &                       
        #.44 CO2 + #.41 RCOOH + &                                                
        #.15 O3 + #.15 RCOOH + #.56 XC                                          SAME K AS BR22           99    82     0 = -4011     9 -4011    87                
                                                                                                                    -4011    92 -4011    84
                                                                                                                    -4012    96 -4013     1
                                                                                                                    -4013    96 -4017    98
 BR56   MACO3 + NO3 = NO2 + CO2 + HCHO + MECO3 + O2                             SAME K AS BR09           99    10     0 =     3    84    87    92                
                                                                                                                        4
 BR57   MACO3 + MEO2 = #2 HCHO + HO2 + CO2 + MECO3                              SAME K AS BR24           99    86     0 = -4001    87    82    84                
                                                                                                                       92
 BR58   MACO3 + RO2C = CO2 + HCHO + MECO3                                       SAME K AS BR25           99    90     0 =    84    87    92                      
 BR59   MACO3 + RO2XC = CO2 + HCHO + MECO3                                      SAME K AS BR25           99    91     0 =    84    87    92                      
 BR60   MACO3 + MECO3 = #2 CO2 + MEO2 + HCHO + MECO3 + O2                       SAME K AS BR27           99    92     0 = -4001    84    86    87                
                                                                                                                       92     4
 BR61   MACO3 + RCO3 = HCHO + MECO3 + RO2C + xHO2 + &                            
        yROOH + xCCHO + #2 CO2                                                  SAME K AS BR27           99    94     0 =    87    92    90    50                
                                                                                                                       76    60 -4001    84
 BR62   MACO3 + BZCO3 = HCHO + MECO3 + BZO + RO2C + &                            
        #2 CO2                                                                  SAME K AS BR27           99    97     0 =    87    92    48    90                
                                                                                                                    -4001    84
 BR63   MACO3 + MACO3 = #2 HCHO + #2 MECO3 + #2 CO2                             SAME K AS BR27           99    99     0 = -4001    87 -4001    92                
                                                                                                                    -4001    84
 BR64   TBUO + NO2 = RNO3 + #-1 XC + #-1 XC                                3.523E+04   0.00  -1.000      47     3     0 =    19 -4018    98 -4018                
                                                                                                                       98
 BR65   TBUO = ACET + MEO2                                                 4.500E+16  16.20   0.000      47     0     0 =   100    86                            
 BR66   BZO + NO2 = NPHE                                                   3.376E+04  -0.30  -1.000      48     3     0 =    17                                  
 BR67   BZO + HO2 = CRES + #-1 XC                                               SAME K AS BR08           48    82     0 =   101 -4018    98                      
 BR68   BZO = CRES + RO2C + xHO2 + #-1 XC                                  6.000E-02   0.00   0.000      48     0     0 =   101    90    50 -4018                
                                                                                                                       98
 RO01   xHO2 + #RO2RO = HO2                                                1.000E+00   0.00   0.000      50    -5     0 =    82                                  
 RO02   xHO2 + #RO2XRO =                                                   1.000E+00   0.00   0.000      50    -6     0 =                                        
 RO03   xOH + #RO2RO = OH                                                  1.000E+00   0.00   0.000      51    -5     0 =     9                                  
 RO04   xOH + #RO2XRO =                                                    1.000E+00   0.00   0.000      51    -6     0 =                                        
 RO05   xNO2 + #RO2RO = NO2                                                1.000E+00   0.00   0.000      52    -5     0 =     3                                  
 RO06   xNO2 + #RO2XRO = XN                                                1.000E+00   0.00   0.000      52    -6     0 =    18                                  
 RO07   xMEO2 + #RO2RO = MEO2                                              1.000E+00   0.00   0.000      53    -5     0 =    86                                  
 RO08   xMEO2 + #RO2XRO = XC                                               1.000E+00   0.00   0.000      53    -6     0 =    98                                  
 RO09   xMECO3 + #RO2RO = MECO3                                            1.000E+00   0.00   0.000      54    -5     0 =    92                                  
 RO10   xMECO3 + #RO2XRO = #2 XC                                           1.000E+00   0.00   0.000      54    -6     0 = -4001    98                            
 RO11   xRCO3 + #RO2RO = RCO3                                              1.000E+00   0.00   0.000      55    -5     0 =    94                                  
 RO12   xRCO3 + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      55    -6     0 = -4019    98                            
 RO13   xMACO3 + #RO2RO = MACO3                                            1.000E+00   0.00   0.000      56    -5     0 =    99                                  
 RO14   xMACO3 + #RO2XRO = #4 XC                                           1.000E+00   0.00   0.000      56    -6     0 = -4020    98                            
 RO15   xTBUO + #RO2RO = TBUO                                              1.000E+00   0.00   0.000      57    -5     0 =    47                                  
 RO16   xTBUO + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      57    -6     0 = -4020    98                            
 RO17   xCO + #RO2RO = CO                                                  1.000E+00   0.00   0.000      58    -5     0 =    83                                  
 RO18   xCO + #RO2XRO = XC                                                 1.000E+00   0.00   0.000      58    -6     0 =    98                                  
 BP01   HCHO + HV = #2 HO2 + CO                                                 PHOT. = HCHOR-06         87     8     0 = -4001    82    83                      
 BP02   HCHO + HV = H2 + CO                                                     PHOT. = HCHOM-06         87     8     0 =     7    83                            
 BP03   HCHO + OH = HO2 + CO + H2O                                         7.927E+03  -0.27  -1.000      87     9     0 =    82    83     6                      
 BP07   HCHO + NO3 = HNO3 + HO2 + CO                                       2.936E+03   4.83  -1.000      87    10     0 =    12    82    83                      
 BP08   CCHO + OH = MECO3 + H2O                                            6.459E+03  -0.73  -1.000     102     9     0 =    92     6                            
 BP09   CCHO + HV = CO + HO2 + MEO2                                             PHOT. = CCHO_R          102     8     0 =    83    82    86                      
 BP10   CCHO + NO3 = HNO3 + MECO3                                          2.055E+03   3.70  -1.000     102    10     0 =    12    92                            
 BP11   RCHO + OH = #.965 RCO3 + #.035 RO2C + &                                  
        #.035 xHO2 + #.035 xCO + #.035 xCCHO + #.035 yROOH                 7.486E+03  -0.81  -1.000     103     9     0 = -4021    94 -4022    90                
                                                                                                                    -4022    50 -4022    58
                                                                                                                    -4022    60 -4022    76
 BP12   RCHO + HV = RO2C + xHO2 + yROOH + xCCHO + CO + &                         
        HO2                                                                     PHOT. = C2CHO           103     8     0 =    90    50    76    60                
                                                                                                                       83    82
 BP13   RCHO + NO3 = HNO3 + RCO3                                           2.055E+03   3.18  -1.000     103    10     0 =    12    94                            
 BP14   ACET + OH = RO2C + xMECO3 + xHCHO + yROOH                          6.694E+01  -0.85   2.650     100     9     0 =    90    54    59    76                
 BP15   ACET + HV + #QY.ACET = #.62 MECO3 + #1.38 MEO2 + #.38 CO                PHOT. = ACET-06         100     8    -7 = -4023    92 -4024    86                
                                                                                                                    -4025    83
 BP16   MEK + OH = #.967 RO2C + #.039 RO2XC + #.039 zRNO3 + &                    
        #.376 xHO2 + #.51 xMECO3 + #.074 xRCO3 + &                               
        #.088 xHCHO + #.504 xCCHO + #.376 xRCHO + yROOH + &                      
        #.3 XC                                                             1.908E+03   0.05   1.000     104     9     0 = -4026    90 -4027    91                
                                                                                                                    -4027    75 -4028    50
                                                                                                                    -4029    54 -4030    55
                                                                                                                    -4031    59 -4032    60
                                                                                                                    -4028    61    76 -4033
                                                                                                                       98
 BP17   MEK + HV + #QY.MEK = MECO3 + RO2C + xHO2 + xCCHO + yROOH                PHOT. = MEK-06          104     8    -8 =    92    90    50    60                
                                                                                                                       76
 BP18   MEOH + OH = HCHO + HO2                                             4.184E+03   0.69  -1.000      89     9     0 =    87    82                            
 BP19   HCOOH + OH = HO2 + CO2                                             6.606E+02   0.00  -1.000     105     9     0 =    82    84                            
 BP20   CCOOH + OH = #.509 MEO2 + #.491 RO2C + #.509 CO2 + &                     
        #.491 xHO2 + #.491 xMGLY + #.491 yROOH + &                               
        #-0.491 XC                                                         6.165E+01  -1.70  -1.000      93     9     0 = -4034    86 -4035    90                
                                                                                                                    -4034    84 -4035    50
                                                                                                                    -4035    66 -4035    76
                                                                                                                    -4036    98
 BP21   RCOOH + OH = RO2C + xHO2 + #.143 CO2 + &                                 
        #.142 xCCHO + #.4 xRCHO + #.457 xBACL + yROOH + &                        
        #-0.455 XC                                                         1.761E+03   0.00  -1.000      96     9     0 =    90    50 -4037    84                
                                                                                                                    -4038    60 -4010    61
                                                                                                                    -4039    67    76 -4040
                                                                                                                       98
 BP22   COOH + OH = H2O + #.3 HCHO + #.3 OH + #.7 MEO2                     5.578E+03  -0.40  -1.000      88     9     0 =     6 -4033    87 -4033                
                                                                                                                        9 -4041    86
 BP23   COOH + HV = HCHO + HO2 + OH                                             PHOT. = COOH             88     8     0 =    87    82     9                      
 BP24   ROOH + OH = #.744 OH + #.251 RO2C + #.004 RO2XC + &                      
        #.004 zRNO3 + #.744 RCHO + #.239 xHO2 + &                                
        #.012 xOH + #.012 xHCHO + #.012 xCCHO + &                                
        #.205 xRCHO + #.034 xPROD2 + #.256 yROOH + &                             
        #-0.115 XC                                                         3.670E+04   0.00  -1.000     106     9     0 = -4042     9 -4043    90                
                                                                                                                    -4044    91 -4044    75
                                                                                                                    -4042   103 -4045    50
                                                                                                                    -4046    51 -4046    59
                                                                                                                    -4046    60 -4047    61
                                                                                                                    -4048    64 -4049    76
                                                                                                                    -4050    98
 BP25   ROOH + HV = RCHO + HO2 + OH                                             PHOT. = COOH            106     8     0 =   103    82     9                      
 BP26   R6OOH + OH = #.84 OH + #.222 RO2C + #.029 RO2XC + &                      
        #.029 zRNO3 + #.84 PROD2 + #.09 xHO2 + #.041 xOH + &                     
        #.02 xCCHO + #.075 xRCHO + #.084 xPROD2 + &                              
        #.16 yROOH + #.02 XC                                               8.220E+04   0.00  -1.000     107     9     0 = -4051     9 -4052    90                
                                                                                                                    -4053    91 -4053    75
                                                                                                                    -4051   108 -4054    50
                                                                                                                    -4055    51 -4056    60
                                                                                                                    -4057    61 -4058    64
                                                                                                                    -4059    76 -4056    98
 BP27   R6OOH + HV = OH + #.142 HO2 + #.782 RO2C + &                             
        #.077 RO2XC + #.077 zRNO3 + #.085 RCHO + &                               
        #.142 PROD2 + #.782 xHO2 + #.026 xCCHO + &                               
        #.058 xRCHO + #.698 xPROD2 + #.858 yR6OOH + &                            
        #.017 XC                                                                PHOT. = COOH            107     8     0 =     9 -4038    82 -4060                
                                                                                                                       90 -4061    91 -4061
                                                                                                                       75 -4062   103 -4038
                                                                                                                      108 -4060    50 -4063
                                                                                                                       60 -4064    61 -4065
                                                                                                                       64 -4066    77 -4067
                                                                                                                       98
 BP28   RAOOH + OH = #.139 OH + #.148 HO2 + #.589 RO2C + &                       
        #.124 RO2XC + #.124 zRNO3 + #.074 PROD2 + &                              
        #.147 MGLY + #.139 IPRD + #.565 xHO2 + #.024 xOH + &                     
        #.448 xRCHO + #.026 xGLY + #.030 xMEK + &                                
        #.252 xMGLY + #.073 xAFG1 + #.073 xAFG2 + &                              
        #.713 yR6OOH + #2.674 XC                                           2.070E+05   0.00  -1.000     109     9     0 = -4068     9 -4069    82                
                                                                                                                    -4070    90 -4071    91
                                                                                                                    -4071    75 -4030   108
                                                                                                                    -4072   110 -4068   111
                                                                                                                    -4073    50 -4074    51
                                                                                                                    -4075    61 -4063    65
                                                                                                                    -4076    63 -4077    66
                                                                                                                    -4078    69 -4078    70
                                                                                                                    -4079    77 -4080    98
 BP29   RAOOH + HV = OH + HO2 + #.5 GLY + #.5 MGLY + &                           
        #.5 AFG1 + #.5 AFG2 + &                                                  
        #.5 XC                                                                  PHOT. = COOH            109     8     0 =     9    82 -4006   112                
                                                                                                                    -4006   110 -4006   113
                                                                                                                    -4006   114 -4006    98
 BP30   GLY + HV = #2 CO + #2 HO2                                               PHOT. = GLY-07R         112     8     0 = -4001    83 -4001    82                
 BP31   GLY + HV = HCHO + CO                                                    PHOT. = GLY-07M         112     8     0 =    87    83                            
 BP32   GLY + OH = #.7 HO2 + #1.4 CO + #.3 HCOCO3                          4.550E+03  -0.68  -1.000     112     9     0 = -4041    82 -4081    83                
                                                                                                                    -4033   115
 BP33   GLY + NO3 = HNO3 + #.7 HO2 + #1.4 CO + #.3 HCOCO3                  4.110E+03   4.72  -1.000     112    10     0 =    12 -4041    82 -4081                
                                                                                                                       83 -4033   115
 BP80   HCOCO3 + NO = HO2 + CO + CO2 + NO2                                      SAME K AS BR31          115     2     0 =    82    83    84     3                
 BP81   HCOCO3 + NO2 = HO2 + CO + CO2 + NO3                                     SAME K AS BR28          115     3     0 =    82    83    84    10                
 BP82   HCOCO3 + HO2 =  #.44 HO2 + #.44 CO + &                                   
        #.44 CO2 + #.56 GLY + #.15 O3 + #.44 OH                                 SAME K AS BR22          115    82     0 = -4011    82 -4011    83                
                                                                                                                    -4011    84 -4017   112
                                                                                                                    -4013     1 -4011     9
 BP34   MGLY + HV = HO2 + CO + MECO3                                            PHOT. = MGLY-06         110     8     0 =    82    83    92                      
 BP35   MGLY + OH = CO + MECO3                                             2.202E+04   0.00  -1.000     110     9     0 =    83    92                            
 BP36   MGLY + NO3 = HNO3 + CO + MECO3                                     2.055E+03   3.77  -1.000     110    10     0 =    12    83    92                      
 BP37   BACL + HV = #2 MECO3                                                    PHOT. = BACL-07         116     8     0 = -4001    92                            
 BP38   CRES + OH = #.700 HO2 + #.100 BZO + #.170 xHO2 + &                       
        #.030 OH + #.170 RO2C + #.700 CATL + #.030 AFG3 + &                      
        #.085 xAFG1 + #.085 xAFG2 + #.085 xGLY + &                               
        #.085 xMGLY + #.100 xCNDPP + #.130 yRAOOH + &                            
        #.040 yRAOOHP + #.115 XC                                           2.349E+03  -1.93  -1.000     101     9     0 = -4041    82 -4014    48                
                                                                                                                    -4082    50 -4076     9
                                                                                                                    -4082    90 -4041   117
                                                                                                                    -4076   118 -4062    69
                                                                                                                    -4062    70 -4062    65
                                                                                                                    -4062    66 -4014   119
                                                                                                                    -4083    78 -4084   120
                                                                                                                    -4085    98
 BP39   CRES + NO3 = #.100 HNO3 + #.900 XN + #.700 HO2 + &                       
        #.100 BZO + #.170 xHO2 + #.030 OH + #.170 RO2C + &                       
        #.700 CATL + #.030 AFG3 + #.085 xAFG1 + &                                
        #.085 xAFG2 + #.085 xGLY + #.085 xMGLY + &                               
        #.130 yRAOOH + #.715 XC                                            2.055E+04   0.00  -1.000     101    10     0 = -4014    12 -4015    18                
                                                                                                                    -4041    82 -4014    48
                                                                                                                    -4082    50 -4076     9
                                                                                                                    -4082    90 -4041   117
                                                                                                                    -4076   118 -4062    69
                                                                                                                    -4062    70 -4062    65
                                                                                                                    -4062    66 -4083    78
                                                                                                                    -4086    98
 BP85   XYNL + OH = #.700 HO2 + #.070 BZO + #.230 xHO2 + &                       
        #.230 RO2C + #.700 CATL + #.115 xAFG1 + &                                
        #.115 xAFG2 + #.115 xGLY + #.115 xMGLY + &                               
        #.110 xCNDPP + #.020 yRAOOH + #.210 yRAOOHP + &                          
        #.995 XC                                                           1.083E+05   0.00  -1.000     121     9     0 = -4041    82 -4087    48                
                                                                                                                    -4088    50 -4088    90
                                                                                                                    -4041   117 -4085    69
                                                                                                                    -4085    70 -4085    65
                                                                                                                    -4085    66 -4089   119
                                                                                                                    -4056    78 -4090   120
                                                                                                                    -4091    98
 BP86   XYNL + NO3 = #.070 HNO3 + #.930 XN + #.700 HO2 + &                       
        #.070 BZO + #.230 xHO2 + #.230 RO2C + #.700 CATL + &                     
        #.115 xAFG1 + #.115 xAFG2 + #.115 xGLY + &                               
        #.115 xMGLY + #.020 yRAOOH + #1.655 XC                             4.492E+04   0.00  -1.000     121    10     0 = -4087    12 -4092    18                
                                                                                                                    -4041    82 -4087    48
                                                                                                                    -4088    50 -4088    90
                                                                                                                    -4041   117 -4085    69
                                                                                                                    -4085    70 -4085    65
                                                                                                                    -4085    66 -4056    78
                                                                                                                    -4093    98
 BP87   CATL + OH = #.400 HO2 + #.200 BZO + #.200 xHO2 + &                       
        #.200 OH + #.200 RO2C + #.200 AFG3 + #.100 xAFG1 + &                     
        #.100 xAFG2 + #.100 xGLY + #.100 xMGLY + &                               
        #.330 CNDPP + #.200 yRAOOH + #-.080 XC                             2.936E+05   0.00  -1.000     117     9     0 = -4010    82 -4005    48                
                                                                                                                    -4005    50 -4005     9
                                                                                                                    -4005    90 -4005   118
                                                                                                                    -4014    69 -4014    70
                                                                                                                    -4014    65 -4014    66
                                                                                                                    -4094   122 -4005    78
                                                                                                                    -4095    98
 BP88   CATL + NO3 = #.200 HNO3 + #.800 XN + #.400 HO2 + &                       
        #.200 BZO + #.200 xHO2 + #.200 OH + #.200 RO2C + &                       
        #.200 AFG3 + #.100 xAFG1 + #.100 xAFG2 + &                               
        #.100 xGLY + #.100 xMGLY + #.200 yRAOOH + &                              
        #1.900 XC                                                          2.495E+05   0.00  -1.000     117    10     0 = -4005    12 -4004    18                
                                                                                                                    -4010    82 -4005    48
                                                                                                                    -4005    50 -4005     9
                                                                                                                    -4005    90 -4005   118
                                                                                                                    -4014    69 -4014    70
                                                                                                                    -4014    65 -4014    66
                                                                                                                    -4005    78 -4096    98
 BP40   NPHE + OH = BZO + XN                                               5.138E+03   0.00  -1.000      17     9     0 =    48    18                            
 BP41   NPHE + HV + #QY.NO2A = HONO + #6 XC                                     PHOT. = NO2-06           17     8    -9 =    13 -4097    98                      
 BP42   NPHE + HV + #QY.NPHE = #6 XC + XN                                       PHOT. = NO2-06           17     8   -10 = -4097    98    18                      
 BP43   BALD + OH = BZCO3                                                  1.761E+04   0.00  -1.000     123     9     0 =    97                                  
 BP44   BALD + HV + #QY.BALD = #7 XC                                            PHOT. = BALD-06         123     8   -11 = -4098    98                            
 BP45   BALD + NO3 = HNO3 + BZCO3                                          1.967E+03   3.70  -1.000     123    10     0 =    12    97                            
 BP83   PHEN + OH = #.700 HO2 + #.100 BZO + #.095 xHO2 + &                       
        #.105 OH + #.095 RO2C + #.700 CATL + #.105 AFG3 + &                      
        #.048 xAFG1 + #.048 xAFG2 + #.095 xGLY + &                               
        #.065 yRAOOH + #.030 yRAOOHP + #-.205 XC                           6.899E+02  -2.42  -1.000     124     9     0 = -4041    82 -4014    48                
                                                                                                                    -4099    50 -4100     9
                                                                                                                    -4099    90 -4041   117
                                                                                                                    -4100   118 -4101    69
                                                                                                                    -4101    70 -4099    65
                                                                                                                    -4102    78 -4076   120
                                                                                                                    -4103    98
 BP84   PHEN + NO3 = #.100 HNO3 + #.900 XN + #.700 HO2 + &                       
        #.100 BZO + #.095 xHO2 + #.105 OH + #.095 RO2C + &                       
        #.700 CATL + #.105 AFG3 + #.048 xAFG1 + &                                
        #.048 xAFG2 + #.095 xGLY + #.065 yRAOOH + &                              
        #.030 yRAOOHP + #-.205 XC                                          5.578E+03   0.00  -1.000     124    10     0 = -4014    12 -4015    18                
                                                                                                                    -4041    82 -4014    48
                                                                                                                    -4099    50 -4100     9
                                                                                                                    -4099    90 -4041   117
                                                                                                                    -4100   118 -4101    69
                                                                                                                    -4101    70 -4099    65
                                                                                                                    -4102    78 -4076   120
                                                                                                                    -4103    98
 BP46   AFG1 + OH = #.217 MACO3 + #.723 RO2C + &                                 
        #.060 RO2XC + #.060 zRNO3 + #.521 xHO2 + &                               
        #.201 xMECO3 + #.334 xCO + #.407 xRCHO + &                               
        #.129 xMEK + #.107 xGLY + #.267 xMGLY + &                                
        #.783 yR6OOH + #.284 XC                                            1.086E+05   0.00  -1.000     113     9     0 = -4104    99 -4105    90                
                                                                                                                    -4106    91 -4106    75
                                                                                                                    -4107    50 -4108    54
                                                                                                                    -4109    58 -4110    61
                                                                                                                    -4111    63 -4112    65
                                                                                                                    -4113    66 -4114    77
                                                                                                                    -4115    98
 BP48   AFG1 + HV = #1.023 HO2 + #.173 MEO2 + &                                  
        #.305 MECO3 + #.500 MACO3 + #.695 CO + #.195 GLY + &                     
        #.305 MGLY + #.217 XC                                                   PHOT. = AFG1            113     8     0 = -4116    82 -4117    86                
                                                                                                                    -4118    92 -4006    99
                                                                                                                    -4119    83 -4120   112
                                                                                                                    -4118   110 -4104    98
 BP49   AFG2 + OH = #.217 MACO3 + #.723 RO2C + &                                 
        #.060 RO2XC + #.060 zRNO3 + #.521 xHO2 + &                               
        #.201 xMECO3 + #.334 xCO + #.407 xRCHO + &                               
        #.129 xMEK + #.107 xGLY + #.267 xMGLY + &                                
        #.783 yR6OOH + #.284 XC                                            1.086E+05   0.00  -1.000     114     9     0 = -4104    99 -4105    90                
                                                                                                                    -4106    91 -4106    75
                                                                                                                    -4107    50 -4108    54
                                                                                                                    -4109    58 -4110    61
                                                                                                                    -4111    63 -4112    65
                                                                                                                    -4113    66 -4114    77
                                                                                                                    -4115    98
 BP51   AFG2 + HV = PROD2 + #-1 XC                                              PHOT. = AFG1            114     8     0 =   108 -4018    98                      
 BP52   AFG3 + OH = #.206 MACO3 + #.733 RO2C + &                                 
        #.117 RO2XC + #.117 zRNO3 + #.561 xHO2 + &                               
        #.117 xMECO3 + #.114 xCO + #.274 xGLY + &                                
        #.153 xMGLY + #.019 xBACL + #.195 xAFG1 + &                              
        #.195 xAFG2 + #.231 xIPRD + #.794 yR6OOH + &                             
        #.938 XC                                                           1.372E+05   0.00  -1.000     118     9     0 = -4121    99 -4122    90                
                                                                                                                    -4123    91 -4123    75
                                                                                                                    -4124    50 -4123    54
                                                                                                                    -4125    58 -4126    65
                                                                                                                    -4127    66 -4128    67
                                                                                                                    -4120    69 -4120    70
                                                                                                                    -4129    73 -4130    77
                                                                                                                    -4131    98
 BP53   AFG3 + O3 = #.471 OH + #.554 HO2 + #.013 MECO3 + &                       
        #.258 RO2C + #.007 RO2XC + #.007 zRNO3 + #.580 CO + &                    
        #.190 CO2 + #.366 GLY + #.184 MGLY + #.350 AFG1 + &                      
        #.350 AFG2 + #.139 AFG3 + #.003 MACR + #.004 MVK + &                     
        #.003 IPRD + #.095 xHO2 + #.163 xRCO3 + &                                
        #.163 xHCHO + #.095 xMGLY + #.264 yR6OOH + &                             
        #-.575 XC                                                          2.099E-02   0.00  -1.000     118     1     0 = -4132     9 -4133    82                
                                                                                                                    -4134    92 -4135    90
                                                                                                                    -4136    91 -4136    75
                                                                                                                    -4137    83 -4138    84
                                                                                                                    -4139   112 -4140   110
                                                                                                                    -4141   113 -4141   114
                                                                                                                    -4068   118 -4142   125
                                                                                                                    -4044   126 -4142   111
                                                                                                                    -4099    50 -4143    55
                                                                                                                    -4143    59 -4099    66
                                                                                                                    -4144    77 -4145    98
 BP89   AFG4 + OH = #.902 RO2C + #.098 RO2XC + &                                 
        #.098 zRNO3 + #.902 xMECO3 + #.902 xRCHO + yROOH + &                     
        #.902 XC                                                           9.248E+04   0.00  -1.000     127     9     0 = -4146    90 -4147    91                
                                                                                                                    -4147    75 -4146    54
                                                                                                                    -4146    61    76 -4146
                                                                                                                       98
 BP54   MACR + OH = #.5 MACO3 + #.5 RO2C + #.5 xHO2 + &                          
        #.416 xCO + #.084 xHCHO + #.416 xMEK + &                                 
        #.084 xMGLY + #.5 yROOH + #-0.416 XC                               1.174E+04  -0.75  -1.000     125     9     0 = -4006    99 -4006    90                
                                                                                                                    -4006    50 -4148    58
                                                                                                                    -4058    59 -4148    63
                                                                                                                    -4058    66 -4006    76
                                                                                                                    -4149    98
 BP55   MACR + O3 = #.208 OH + #.108 HO2 + #.1 RO2C + &                          
        #.45 CO + #.117 CO2 + #.1 HCHO + #.9 MGLY + &                            
        #.333 HCOOH + #.1 xRCO3 + #.1 xHCHO + #.1 yROOH + &                      
        #-0.1 XC                                                           2.055E+00   4.17  -1.000     125     1     0 = -4150     9 -4151    82                
                                                                                                                    -4014    90 -4152    83
                                                                                                                    -4123    84 -4014    87
                                                                                                                    -4015   110 -4153   105
                                                                                                                    -4014    55 -4014    59
                                                                                                                    -4014    76 -4154    98
 BP56   MACR + NO3 = #.5 MACO3 + #.5 RO2C + #.5 HNO3 + &                         
        #.5 xHO2 + #.5 xCO + &                                                   
        #.5 yROOH + #1.5 XC + #.5 XN                                       2.202E+03   3.61  -1.000     125    10     0 = -4006    99 -4006    90                
                                                                                                                    -4006    12 -4006    50
                                                                                                                    -4006    58 -4006    76
                                                                                                                    -4155    98 -4006    18
 BP57   MACR + O3P = RCHO + XC                                             9.306E+03   0.00  -1.000     125    45     0 =   103    98                            
 BP58   MACR + HV = #.33 OH + #.67 HO2 + #.34 MECO3 + &                          
        #.33 MACO3 + #.33 RO2C + #.67 CO + #.34 HCHO + &                         
        #.33 xMECO3 + #.33 xHCHO + #.33 yROOH                                   PHOT. = MACR-06         125     8     0 = -4094     9 -4156    82                
                                                                                                                    -4157    92 -4094    99
                                                                                                                    -4094    90 -4156    83
                                                                                                                    -4157    87 -4094    54
                                                                                                                    -4094    59 -4094    76
 BP59   MVK + OH = #.975 RO2C + #.025 RO2XC + &                                  
        #.025 zRNO3 + &                                                          
        #.3 xHO2 + #.675 xMECO3 + #.3 xHCHO + &                                  
        #.675 xRCHO + #.3 xMGLY + yROOH + #-0.725 XC                       3.817E+03  -1.21  -1.000     126     9     0 = -4158    90 -4159    91                
                                                                                                                    -4159    75 -4033    50
                                                                                                                    -4160    54 -4033    59
                                                                                                                    -4160    61 -4033    66
                                                                                                                       76 -4161    98
 BP60   MVK + O3 = #.164 OH + #.064 HO2 + #.05 RO2C + &                          
        #.05 xHO2 + #.475 CO + #.124 CO2 + #.05 HCHO + &                         
        #.95 MGLY + #.351 HCOOH + #.05 xRCO3 + &                                 
        #.05 xHCHO + #.05 yROOH + #-0.05 XC                                1.248E+00   3.02  -1.000     126     1     0 = -4162     9 -4163    82                
                                                                                                                    -4164    90 -4164    50
                                                                                                                    -4165    83 -4071    84
                                                                                                                    -4164    87 -4166   110
                                                                                                                    -4167   105 -4164    55
                                                                                                                    -4164    59 -4164    76
                                                                                                                    -4168    98
 BP62   MVK + O3P = #.45 RCHO + #.55 MEK + #.45 XC                         6.341E+03   0.00  -1.000     126    45     0 = -4152   103 -4169   104                
                                                                                                                    -4152    98
 BP63   MVK + HV = #.4 MEO2 + #.6 CO + #.6 PROD2 + &                             
        #.4 MACO3 + #-2.2 XC                                                    PHOT. = MVK-06          126     8     0 = -4010    86 -4009    83                
                                                                                                                    -4009   108 -4010    99
                                                                                                                    -4170    98
 BP64   IPRD + OH = #.289 MACO3 + #.67 RO2C + #.67 xHO2 + &                      
        #.041 RO2XC + #.041 zRNO3 + #.336 xCO + #.055 xHCHO + &                  
        #.129 xCCHO + #.013 xRCHO + #.15 xMEK + &                                
        #.332 xPROD2 + #.15 xGLY + #.174 xMGLY + &                               
        #-0.504 XC + #.711 yR6OOH                                          9.086E+04   0.00  -1.000     111     9     0 = -4171    99 -4156    90                
                                                                                                                    -4156    50 -4055    91
                                                                                                                    -4055    75 -4172    58
                                                                                                                    -4173    59 -4111    60
                                                                                                                    -4134    61 -4013    63
                                                                                                                    -4174    64 -4013    65
                                                                                                                    -4175    66 -4176    98
                                                                                                                    -4177    77
 BP65   IPRD + O3 = #.285 OH + #.4 HO2 + #.048 RO2C + &                          
        #.048 xRCO3 + #.498 CO + #.14 CO2 + #.124 HCHO + &                       
        #.21 MEK + #.023 GLY + #.742 MGLY + #.1 HCOOH + &                        
        #.372 RCOOH + #.047 xCCHO + #.001 xHCHO + &                              
        #.048 yR6OOH + #-.329 XC                                           6.136E-03   0.00  -1.000     111     1     0 = -4178     9 -4010    82                
                                                                                                                    -4101    90 -4101    55
                                                                                                                    -4179    83 -4180    84
                                                                                                                    -4071    87 -4090   104
                                                                                                                    -4181   112 -4182   110
                                                                                                                    -4014   105 -4183    96
                                                                                                                    -4184    60 -4185    59
                                                                                                                    -4101    77 -4186    98
 BP66   IPRD + NO3 = #.15 MACO3 + #.15 HNO3 + #.799 RO2C + &                     
        #.799 xHO2 + #.051 RO2XC + #.051 zRNO3 + #.572 xCO + &                   
        #.227 xHCHO + #.218 xRCHO + #.008 xMGLY + &                              
        #.572 xRNO3 + #.85 yR6OOH + #.278 XN + #-1 XC + #.185 XC           1.468E+02   0.00  -1.000     111    10     0 = -4013    99 -4013    12                
                                                                                                                    -4187    90 -4187    50
                                                                                                                    -4188    91 -4188    75
                                                                                                                    -4189    58 -4190    59
                                                                                                                    -4191    61 -4192    66
                                                                                                                    -4189    74 -4193    77
                                                                                                                    -4194    18 -4018    98
                                                                                                                    -4195    98
 BP67   IPRD + HV = #1.233 HO2 + #.467 MECO3 + #.3 RCO3 + &                      
        #1.233 CO + #.3 HCHO + #.467 CCHO + #.233 MEK + &                        
        #-.233 XC                                                               PHOT. = MACR-06         111     8     0 = -4196    82 -4197    92                
                                                                                                                    -4033    94 -4196    83
                                                                                                                    -4033    87 -4197   102
                                                                                                                    -4198   104 -4199    98
 BP68   PROD2 + OH = #.472 HO2 + #.379 xHO2 + &                                  
        #.029 xMECO3 + #.049 xRCO3 + #.473 RO2C + &                              
        #.071 RO2XC + #.071 zRNO3 + #.002 HCHO + &                               
        #.211 xHCHO + #.001 CCHO + #.083 xCCHO + &                               
        #.143 RCHO + #.402 xRCHO + #.115 xMEK + &                                
        #.329 PROD2 + #.007 xPROD2 + #.528 yR6OOH + &                            
        #.877 XC                                                           2.275E+04   0.00  -1.000     108     9     0 = -4200    82 -4201    50                
                                                                                                                    -4053    54 -4202    55
                                                                                                                    -4203    90 -4204    91
                                                                                                                    -4204    75 -4205    87
                                                                                                                    -4206    59 -4185   102
                                                                                                                    -4207    60 -4037   103
                                                                                                                    -4208    61 -4085    63
                                                                                                                    -4209   108 -4136    64
                                                                                                                    -4210    77 -4211    98
 BP69   PROD2 + HV + #QY.PROD2 = #.913 xHO2 + #.4 MECO3 + #.6 RCO3 + &           
        #1.59 RO2C + #.087 RO2XC + #.087 zRNO3 + &                               
        #.303 xHCHO + #.163 xCCHO + #.78 xRCHO + yR6OOH + &                      
        #-.091 XC                                                               PHOT. = MEK-06          108     8   -12 = -4212    50 -4010    92                
                                                                                                                    -4009    94 -4213    90
                                                                                                                    -4214    91 -4214    75
                                                                                                                    -4215    59 -4143    60
                                                                                                                    -4216    61    77 -4217
                                                                                                                       98
 BP70   RNO3 + OH = #.189 HO2 + #.305 xHO2 + #.019 NO2 + &                       
        #.313 xNO2 + #.976 RO2C + #.175 RO2XC + &                                
        #.175 zRNO3 + #.011 xHCHO + #.429 xCCHO + &                              
        #.001 RCHO + #.036 xRCHO + #.004 xACET + &                               
        #.01 MEK + #.17 xMEK + #.008 PROD2 + &                                   
        #.031 xPROD2 + #.189 RNO3 + #.305 xRNO3 + &                              
        #.157 yROOH + #.636 yR6OOH + #.174 XN + #.04 XC                    1.057E+04   0.00  -1.000      19     9     0 = -4218    82 -4118    50                
                                                                                                                    -4128     3 -4219    52
                                                                                                                    -4220    90 -4221    91
                                                                                                                    -4221    75 -4222    59
                                                                                                                    -4223    60 -4185   103
                                                                                                                    -4224    61 -4044    62
                                                                                                                    -4225   104 -4082    63
                                                                                                                    -4192   108 -4226    64
                                                                                                                    -4218    19 -4118    74
                                                                                                                    -4227    76 -4228    77
                                                                                                                    -4175    18 -4084    98
 BP71   RNO3 + HV = #.344 HO2 + #.554 xHO2 + NO2 + &                             
        #.721 RO2C + #.102 RO2XC + #.102 zRNO3 + &                               
        #.074 HCHO + #.061 xHCHO + #.214 CCHO + &                                
        #.23 xCCHO + #.074 RCHO + #.063 xRCHO + &                                
        #.008 xACET + #.124 MEK + #.083 xMEK + &                                 
        #.19 PROD2 + #.261 xPROD2 + #.066 yROOH + &                              
        #.591 yR6OOH + #.396 XC                                                 PHOT. = IC3ONO2          19     8     0 = -4229    82 -4133    50                
                                                                                                                        3 -4230    90 -4231
                                                                                                                       91 -4231    75 -4030
                                                                                                                       87 -4232    59 -4233
                                                                                                                      102 -4088    60 -4030
                                                                                                                      103 -4234    61 -4192
                                                                                                                       62 -4071   104 -4207
                                                                                                                       63 -4138   108 -4235
                                                                                                                       64 -4236    76 -4237
                                                                                                                       77 -4238    98
 BP72   ETOH + OH = #.95 HO2 + #.05 xHO2 + #.05 RO2C + &                         
        #.081 xHCHO + #.95 CCHO + #.01 xCCHO + &                                 
        #.05 yROOH + #-.001 XC                                             8.059E+02  -1.05   1.000     128     9     0 = -4166    82 -4164    50                
                                                                                                                    -4164    90 -4239    59
                                                                                                                    -4166   102 -4225    60
                                                                                                                    -4164    76 -4240    98
 PO01   xHCHO + #RO2RO  = HCHO                                             1.000E+00   0.00   0.000      59    -5     0 =    87                                  
 PO02   xHCHO + #RO2XRO = XC                                               1.000E+00   0.00   0.000      59    -6     0 =    98                                  
 PO03   xCCHO + #RO2RO = CCHO                                              1.000E+00   0.00   0.000      60    -5     0 =   102                                  
 PO04   xCCHO + #RO2XRO = #2 XC                                            1.000E+00   0.00   0.000      60    -6     0 = -4001    98                            
 PO05   xRCHO + #RO2RO = RCHO                                              1.000E+00   0.00   0.000      61    -5     0 =   103                                  
 PO06   xRCHO + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      61    -6     0 = -4019    98                            
 PO07   xACET + #RO2RO = ACET                                              1.000E+00   0.00   0.000      62    -5     0 =   100                                  
 PO08   xACET + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      62    -6     0 = -4019    98                            
 PO09   xMEK + #RO2RO = MEK                                                1.000E+00   0.00   0.000      63    -5     0 =   104                                  
 PO10   xMEK + #RO2XRO = #4 XC                                             1.000E+00   0.00   0.000      63    -6     0 = -4020    98                            
 PO11   xPROD2 + #RO2RO = PROD2                                            1.000E+00   0.00   0.000      64    -5     0 =   108                                  
 PO12   xPROD2 + #RO2XRO = #6 XC                                           1.000E+00   0.00   0.000      64    -6     0 = -4097    98                            
 PO13   xGLY + #RO2RO = GLY                                                1.000E+00   0.00   0.000      65    -5     0 =   112                                  
 PO14   xGLY + #RO2XRO = #2 XC                                             1.000E+00   0.00   0.000      65    -6     0 = -4001    98                            
 PO15   xMGLY + #RO2RO = MGLY                                              1.000E+00   0.00   0.000      66    -5     0 =   110                                  
 PO16   xMGLY + #RO2XRO = #3 XC                                            1.000E+00   0.00   0.000      66    -6     0 = -4019    98                            
 PO17   xBACL + #RO2RO = BACL                                              1.000E+00   0.00   0.000      67    -5     0 =   116                                  
 PO18   xBACL + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      67    -6     0 = -4020    98                            
 PO19   xBALD + #RO2RO = BALD                                              1.000E+00   0.00   0.000      68    -5     0 =   123                                  
 PO20   xBALD + #RO2XRO = #7 XC                                            1.000E+00   0.00   0.000      68    -6     0 = -4098    98                            
 PO21   xAFG1 + #RO2RO = AFG1                                              1.000E+00   0.00   0.000      69    -5     0 =   113                                  
 PO22   xAFG1 + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      69    -6     0 = -4241    98                            
 PO23   xAFG2 + #RO2RO = AFG2                                              1.000E+00   0.00   0.000      70    -5     0 =   114                                  
 PO24   xAFG2 + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      70    -6     0 = -4241    98                            
 PO50   xAFG4 + #RO2XRO = #6 XC                                            1.000E+00   0.00   0.000     129    -6     0 = -4097    98                            
 PO51   xAFG4 + #RO2RO = AFG4                                              1.000E+00   0.00   0.000     129    -5     0 =   127                                  
 PO27   xMACR + #RO2RO = MACR                                              1.000E+00   0.00   0.000      71    -5     0 =   125                                  
 PO28   xMACR + #RO2XRO = #4 XC                                            1.000E+00   0.00   0.000      71    -6     0 = -4020    98                            
 PO29   xMVK + #RO2RO = MVK                                                1.000E+00   0.00   0.000      72    -5     0 =   126                                  
 PO30   xMVK + #RO2XRO = #4 XC                                             1.000E+00   0.00   0.000      72    -6     0 = -4020    98                            
 PO31   xIPRD + #RO2RO = IPRD                                              1.000E+00   0.00   0.000      73    -5     0 =   111                                  
 PO32   xIPRD + #RO2XRO = #5 XC                                            1.000E+00   0.00   0.000      73    -6     0 = -4241    98                            
 PO33   xRNO3 + #RO2RO = RNO3                                              1.000E+00   0.00   0.000      74    -5     0 =    19                                  
 PO34   xRNO3 + #RO2XRO = #6 XC + XN                                       1.000E+00   0.00   0.000      74    -6     0 = -4097    98    18                      
 PO35   zRNO3 + #RO2NO = RNO3 + #-1 XN                                     1.000E+00   0.00   0.000      75   -13     0 =    19 -4018    18                      
 PO36   zRNO3 + #RO22NN = PROD2 + HO2                                      1.000E+00   0.00   0.000      75   -14     0 =   108    82                            
 PO37   zRNO3 + #RO2XRO = #6 XC                                            1.000E+00   0.00   0.000      75    -6     0 = -4097    98                            
 PO38   yROOH + #RO2HO2 = ROOH + #-3 XC                                    1.000E+00   0.00   0.000      76   -15     0 =   106 -4242    98                      
 PO39   yROOH + #RO2RO2M = MEK + #-4 XC                                    1.000E+00   0.00   0.000      76   -16     0 =   104 -4243    98                      
 PO40   yROOH + #RO2RO =                                                   1.000E+00   0.00   0.000      76    -5     0 =                                        
 PO41   yR6OOH + #RO2HO2 = R6OOH + #-6 XC                                  1.000E+00   0.00   0.000      77   -15     0 =   107 -4244    98                      
 PO42   yR6OOH + #RO2RO2M  = PROD2 + #-6 XC                                1.000E+00   0.00   0.000      77   -16     0 =   108 -4244    98                      
 PO43   yR6OOH + #RO2RO =                                                  1.000E+00   0.00   0.000      77    -5     0 =                                        
 PO44   yRAOOH + #RO2HO2 = RAOOH + #-8 XC                                  1.000E+00   0.00   0.000      78   -15     0 =   109 -4245    98                      
 PO45   yRAOOH + #RO2RO2M  = PROD2 + #-6 XC                                1.000E+00   0.00   0.000      78   -16     0 =   108 -4244    98                      
 PO46   yRAOOH + #RO2RO =                                                  1.000E+00   0.00   0.000      78    -5     0 =                                        
 BE01   CH4 + OH = H2O + MEO2                                              2.716E+03   3.36  -1.000     130     9     0 =     6    86                            
 BE02   ETHENE + OH = RO2C + xHO2 + #1.61 xHCHO + &                              
        #.195 xCCHO + yROOH                                               FALLOFF F= 0.600, N= 1.000    131     9     0 =    90    50 -4246    59                
                                                                                                                    -4120    60    76
                                                                      K0:  3.591E+00   0.00  -6.500
                                                                      KI:  1.292E+04   0.00  -1.850
 BE03   ETHENE + O3 = #.16 OH + #.16 HO2 + #.51 CO + &                           
        #.12 CO2 + HCHO + #.37 HCOOH                                       1.342E+01   5.13  -1.000     131     1     0 = -4059     9 -4059    82                
                                                                                                                    -4029    83 -4247    84
                                                                                                                       87 -4248   105
 BE04   ETHENE + NO3 = RO2C + xHO2 + xRCHO + yROOH + &                           
        #-1 XC + XN                                                        4.844E+03   5.72   1.000     131    10     0 =    90    50    61    76                
                                                                                                                    -4018    98    18
 BE05   ETHENE + O3P = #.8 HO2 + #.51 MEO2 + #.29 RO2C + &                       
        #.51 CO + #.1 CCHO + #.29 xHO2 + #.278 xCO + &                           
        #.278 xHCHO + #.012 xGLY + #.29 yROOH + #.2 XC                     1.571E+04   1.59  -1.000     131    45     0 = -4004    82 -4029    86                
                                                                                                                    -4249    90 -4029    83
                                                                                                                    -4014   102 -4249    50
                                                                                                                    -4194    58 -4194    59
                                                                                                                    -4046    65 -4249    76
                                                                                                                    -4005    98
 BE06   ISOPRENE + OH = #.986 RO2C + #.093 RO2XC + &                             
        #.093 zRNO3 + #.907 xHO2 + #.624 xHCHO + #.23 xMACR + &                  
        #.32 xMVK + #.357 xIPRD + yR6OOH + &                                     
        #-1.000 XC + #0.833 XC + &                                               
        ISOPRXN                                                            3.728E+04  -0.81  -1.000      41     9     0 = -4250    90 -4251    91                
                                                                                                                    -4251    75 -4252    50
                                                                                                                    -4253    59 -4088    71
                                                                                                                    -4254    72 -4255    73
                                                                                                                       77 -4018    98 -4256
                                                                                                                       98    42
 BE07   ISOPRENE + O3 = #.266 OH + #.066 HO2 + &                                 
        #.192 RO2C + #.008 RO2XC + #.008 zRNO3 + #.275 CO + &                    
        #.122 CO2 + #.4 HCHO + #.1 PROD2 + #.39 MACR + &                         
        #.16 MVK + #.15 IPRD + #.204 HCOOH + &                                   
        #.192 xMACO3 + #.192 xHCHO + #.2 yR6OOH + &                              
        #-1.000 XC + #0.401 XC                                             1.154E+01   3.80  -1.000      41     1     0 = -4257     9 -4236    82                
                                                                                                                    -4258    90 -4192    91
                                                                                                                    -4192    75 -4259    83
                                                                                                                    -4260    84 -4010    87
                                                                                                                    -4014   108 -4003   125
                                                                                                                    -4059   126 -4013   111
                                                                                                                    -4261   105 -4258    56
                                                                                                                    -4258    59 -4005    77
                                                                                                                    -4018    98 -4262    98
 BE08   ISOPRENE + NO3 = #.936 RO2C + #.064 RO2XC + &                            
        #.064 zRNO3 + #.749 xHO2 + #.187 xNO2 + #.936 xIPRD + &                  
        yR6OOH + #-1 XC + #.936 XC + #.813 XN                              4.448E+03   0.89  -1.000      41    10     0 = -4263    90 -4163    91                
                                                                                                                    -4163    75 -4264    50
                                                                                                                    -4265    52 -4263    73
                                                                                                                       77 -4018    98 -4263
                                                                                                                       98 -4266    18
 BE09   ISOPRENE + O3P = #.25 MEO2 + #.24 RO2C + &                               
        #.01 RO2XC + #.01 zRNO3 + #.75 PROD2 + #.24 xMACO3 + &                   
        #.24 xHCHO + #.25 yR6OOH + &                                             
        #-1.0 XC                                                           5.138E+04   0.00  -1.000      41    45     0 = -4008    86 -4267    90                
                                                                                                                    -4225    91 -4225    75
                                                                                                                    -4007   108 -4267    56
                                                                                                                    -4267    59 -4008    77
                                                                                                                    -4018    98
 BE10   ACETYLEN + OH = #.7 OH + #.3 HO2 + #.3 CO + &                            
        #.7 GLY + #.3 HCOOH                                               FALLOFF F= 0.600, N= 1.000    132     9     0 = -4041     9 -4033    82                
                                                                                                                    -4033    83 -4041   112
                                                                                                                    -4033   105
                                                                      K0:  1.975E-01   0.00  -2.000
                                                                      KI:  1.218E+03   0.00   1.000
 BE11   ACETYLEN + O3 = #.5 OH + #1.5 HO2 + #1.5 CO + &                          
        #.5 CO2                                                            1.468E+01   8.15  -1.000     132     1     0 = -4006     9 -4155    82                
                                                                                                                    -4155    83 -4006    84
 BE12   BENZENE + OH = #.116 OH + #.29 RO2C + #.29 xHO2 + &                      
        #.024 RO2XC + #.024 zRNO3 + #.57 HO2 + #.57 CRES + &                     
        #.116 AFG3 + #.290 xGLY + #.029 xAFG1 + &                                
        #.261 xAFG2 + #.314 yRAOOH + &                                           
        #-1.0 XC + #0.024 XC + &                                                 
        #.764 BENZRO2                                                      3.420E+03   0.38  -1.000      35     9     0 = -4268     9 -4249    90                
                                                                                                                    -4249    50 -4074    91
                                                                                                                    -4074    75 -4269    82
                                                                                                                    -4269   101 -4268   118
                                                                                                                    -4249    65 -4053    69
                                                                                                                    -4235    70 -4270    78
                                                                                                                    -4018    98 -4074    98
                                                                                                                    -4271    36
 BE12a  BENZRO2 + NO = NO + BNZNRXN                                        3.963E+03   0.71  -1.000      36     2     0 =     2    37                            
 BE12b  BENZRO2 + HO2 = HO2 + BNZHRXN                                      2.789E+02   2.58  -1.000      36    82     0 =    82    38                            
 SSO3   SESQ + O3 = O3 + SESQRXN                                           1.703E+01   0.00  -1.000      43     1     0 =     1    44                            
 SSOH   SESQ + OH = OH + SESQRXN                                           2.892E+05   0.00  -1.000      43     9     0 =     9    44                            
 SSN3   SESQ + NO3 = NO3 + SESQRXN                                         2.789E+04   0.00  -1.000      43    10     0 =    10    44                            
 BL01   ALK1 + OH = xHO2 + RO2C + xCCHO + yROOH                            1.967E+03   0.99   1.000     133     9     0 =    50    90    60    76                
 BL02   ALK2 + OH = #.965 xHO2 + #.965 RO2C + &                                  
        #.035 RO2XC + #.035 zRNO3 + #.261 xRCHO + &                              
        #.704 xACET + yROOH + #-.105 XC                                    2.187E+03   0.17   1.000     134     9     0 = -4021    50 -4021    90                
                                                                                                                    -4022    91 -4022    75
                                                                                                                    -4235    61 -4272    62
                                                                                                                       76 -4273    98
 BL03   ALK3 + OH = #.695 xHO2 + #.236 xTBUO + &                                 
        #1.253 RO2C + #.07 RO2XC + #.07 zRNO3 + &                                
        #.026 xHCHO + #.445 xCCHO + #.122 xRCHO + &                              
        #.024 xACET + #.332 xMEK + #.983 yROOH + &                               
        #.017 yR6OOH + #-.046 XC                                           2.217E+03  -0.25  -1.000     135     9     0 = -4119    50 -4274    57                
                                                                                                                    -4275    90 -4087    91
                                                                                                                    -4087    75 -4063    59
                                                                                                                    -4276    60 -4260    61
                                                                                                                    -4074    62 -4174    63
                                                                                                                    -4277    76 -4067    77
                                                                                                                    -4278    98
 BL04   ALK4 + OH = #.83 xHO2 + #.01 xMEO2 + &                                   
        #.011 xMECO3 + #1.763 RO2C + #.149 RO2XC + &                             
        #.149 zRNO3 + #.002 xCO + #.029 xHCHO + &                                
        #.438 xCCHO + #.236 xRCHO + #.426 xACET + &                              
        #.106 xMEK + #.146 xPROD2 + yR6OOH + #-.119 XC                     5.505E+03  -0.09  -1.000     136     9     0 = -4279    50 -4225    53                
                                                                                                                    -4222    54 -4280    90
                                                                                                                    -4281    91 -4281    75
                                                                                                                    -4205    58 -4053    59
                                                                                                                    -4282    60 -4274    61
                                                                                                                    -4283    62 -4284    63
                                                                                                                    -4285    64    77 -4286
                                                                                                                       98
 BL05   ALK5 + OH = #.647 xHO2 + #1.605 RO2C + &                                 
        #.353 RO2XC + #.353 zRNO3 + #.04 xHCHO + &                               
        #.106 xCCHO + #.209 xRCHO + #.071 xACET + &                              
        #.086 xMEK + #.407 xPROD2 + yR6OOH + #2.004 XC + &                       
        ALKSOMG                                                            2.046E+04   0.00  -1.000      28     9     0 = -4287    50 -4288    90                
                                                                                                                    -4289    91 -4289    75
                                                                                                                    -4084    59 -4284    60
                                                                                                                    -4290    61 -4204    62
                                                                                                                    -4291    63 -4110    64
                                                                                                                       77 -4292    98   137
 BL06   OLE1 + OH = #.904 xHO2 + #.001 xMEO2 + &                                 
        #1.138 RO2C + #.095 RO2XC + #.095 zRNO3 + &                              
        #.7 xHCHO + #.301 xCCHO + #.47 xRCHO + &                                 
        #.005 xACET + #.026 xMACR + #.008 xMVK + &                               
        #.006 xIPRD + #.119 xPROD2 + #.413 yROOH + &                             
        #.587 yR6OOH + #.822 XC                                            9.072E+03  -1.00  -1.000     138     9     0 = -4293    50 -4185    53                
                                                                                                                    -4294    90 -4099    91
                                                                                                                    -4099    75 -4041    59
                                                                                                                    -4295    60 -4296    61
                                                                                                                    -4297    62 -4063    71
                                                                                                                    -4192    72 -4298    73
                                                                                                                    -4299    64 -4300    76
                                                                                                                    -4301    77 -4302    98
 BL07   OLE1 + O3 = #.116 HO2 + #.04 xHO2 + #.193 OH + &                         
        #.104 MEO2 + #.063 RO2C + #.004 RO2XC + &                                
        #.004 zRNO3 + #.368 CO + #.125 CO2 + #.5 HCHO + &                        
        #.147 CCHO + #.007 xCCHO + #.353 RCHO + &                                
        #.031 xRCHO + #.002 xACET + #.006 MEK + &                                
        #.185 HCOOH + #.022 CCOOH + #.112 RCOOH + &                              
        #.189 PROD2 + #.007 yROOH + #.037 yR6OOH + &                             
        #.69 XC                                                            4.624E+00   3.38  -1.000     138     1     0 = -4268    82 -4084    50                
                                                                                                                    -4303     9 -4304    86
                                                                                                                    -4234    90 -4044    91
                                                                                                                    -4044    75 -4305    83
                                                                                                                    -4306    84 -4006    87
                                                                                                                    -4072   102 -4136    60
                                                                                                                    -4289   103 -4226    61
                                                                                                                    -4205    62 -4298   104
                                                                                                                    -4195   105 -4307    93
                                                                                                                    -4308    96 -4218   108
                                                                                                                    -4136    76 -4309    77
                                                                                                                    -4310    98
 BL08   OLE1 + NO3 = #.824 xHO2 + #1.312 RO2C + &                                
        #.176 RO2XC + #.176 zRNO3 + #.009 xCCHO + &                              
        #.002 xRCHO + #.024 xACET + #.546 xRNO3 + &                              
        #.413 yROOH + #.587 yR6OOH + #.454 XN + #.572 XC                   6.943E+02   2.08  -1.000     138    10     0 = -4311    50 -4312    90                
                                                                                                                    -4313    91 -4313    75
                                                                                                                    -4314    60 -4205    61
                                                                                                                    -4074    62 -4315    74
                                                                                                                    -4300    76 -4301    77
                                                                                                                    -4316    18 -4189    98
 BL09   OLE1 + O3P = #.45 RCHO + #.437 MEK + #.113 PROD2 + &                     
        #1.224 XC                                                          2.187E+04   0.65  -1.000     138    45     0 = -4152   103 -4317   104                
                                                                                                                    -4318   108 -4319    98
 BL10   OLE2 + OH = #.914 xHO2 + #.966 RO2C + &                                  
        #.086 RO2XC + #.086 zRNO3 + #.209 xHCHO + &                              
        #.788 xCCHO + #.481 xRCHO + #.136 xACET + &                              
        #.076 xMEK + #.027 xMACR + #.002 xMVK + &                                
        #.037 xIPRD + #.022 xPROD2 + #.357 yROOH + &                             
        #.643 yR6OOH + #.111 XC                                            1.850E+04  -0.97  -1.000     139     9     0 = -4320    50 -4321    90                
                                                                                                                    -4291    91 -4291    75
                                                                                                                    -4290    59 -4322    60
                                                                                                                    -4323    61 -4324    62
                                                                                                                    -4325    63 -4326    71
                                                                                                                    -4205    72 -4309    73
                                                                                                                    -4307    64 -4255    76
                                                                                                                    -4327    77 -4328    98
 BL11   OLE2 + O3 = #.093 HO2 + #.039 xHO2 + #.423 OH + &                        
        #.29 MEO2 + #.147 xMECO3 + #.008 xRCO3 + &                               
        #.2 RO2C + #.003 RO2XC + #.003 zRNO3 + #.297 CO + &                      
        #.162 CO2 + #.152 HCHO + #.108 xHCHO + &                                 
        #.428 CCHO + #.067 xCCHO + #.315 RCHO + &                                
        #.018 xRCHO + #.048 ACET + #.031 MEK + &                                 
        #.001 xMEK + #.033 HCOOH + #.061 CCOOH + &                               
        #.222 RCOOH + #.028 MACR + #.021 MVK + &                                 
        #.042 PROD2 + #.069 yROOH + #.128 yR6OOH + &                             
        #.125 XC                                                           1.195E+01   2.49  -1.000     139     1     0 = -4251    82 -4027    50                
                                                                                                                    -4329     9 -4249    86
                                                                                                                    -4072    54 -4192    55
                                                                                                                    -4005    90 -4142    91
                                                                                                                    -4142    75 -4330    83
                                                                                                                    -4331    84 -4332    87
                                                                                                                    -4151    59 -4333   102
                                                                                                                    -4334    60 -4335   103
                                                                                                                    -4336    61 -4101   100
                                                                                                                    -4226   104 -4185    63
                                                                                                                    -4337   105 -4232    93
                                                                                                                    -4052    96 -4338   125
                                                                                                                    -4339   126 -4340   108
                                                                                                                    -4341    76 -4342    77
                                                                                                                    -4306    98
 BL12   OLE2 + NO3 = #.423 xHO2 + #.409 xNO2 + &                                 
        #.033 xMEO2 + #1.185 RO2C + #.136 RO2XC + &                              
        #.136 zRNO3 + #.074 xHCHO + #.546 xCCHO + &                              
        #.154 xRCHO + #.11 xACET + #.002 xMEK + &                                
        #.026 xMVK + #.007 xIPRD + #.322 xRNO3 + &                               
        #.357 yROOH + #.643 yR6OOH + #.269 XN + #.114 XC                   3.229E+02  -0.76  -1.000     139    10     0 = -4329    50 -4343    52                
                                                                                                                    -4337    53 -4344    90
                                                                                                                    -4324    91 -4324    75
                                                                                                                    -4030    59 -4315    60
                                                                                                                    -4345    61 -4089    62
                                                                                                                    -4205    63 -4063    72
                                                                                                                    -4136    73 -4346    74
                                                                                                                    -4255    76 -4327    77
                                                                                                                    -4347    18 -4125    98
 BL13   OLE2 + O3P = #.014 HO2 + #.007 xHO2 + &                                  
        #.007 xMACO3 + #.013 RO2C + #.001 RO2XC + &                              
        #.001 zRNO3 + #.006 xCO + #.074 RCHO + #.709 MEK + &                     
        #.006 xMACR + #.202 PROD2 + #.014 yROOH + &                              
        #.666 XC                                                           2.099E+04  -0.22  -1.000     139    45     0 = -4348    82 -4136    50                
                                                                                                                    -4136    56 -4134    90
                                                                                                                    -4185    91 -4185    75
                                                                                                                    -4298    58 -4030   103
                                                                                                                    -4349   104 -4298    71
                                                                                                                    -4350   108 -4348    76
                                                                                                                    -4351    98
 BL14   ARO1 + OH = #.166 HO2 + #.482 xHO2 + #.284 OH + &                        
        #.482 RO2C + #.068 RO2XC + #.068 zRNO3 + &                               
        #.218 xGLY + #.138 xMGLY + #.166 CRES + &                                
        #.049 xBALD + #.164 xAFG1 + #.193 xAFG2 + &                              
        #.284 AFG3 + #.077 xPROD2 + #.403 yRAOOH + &                             
        #.147 yR6OOH + #.002 XC + AR1SOMG                                  9.028E+03   0.00  -1.000      33     9     0 = -4352    82 -4353    50                
                                                                                                                    -4115     9 -4353    90
                                                                                                                    -4354    91 -4354    75
                                                                                                                    -4191    65 -4355    66
                                                                                                                    -4352   101 -4202    68
                                                                                                                    -4162    69 -4303    70
                                                                                                                    -4115   118 -4061    64
                                                                                                                    -4356    78 -4072    77
                                                                                                                    -4205    98    34
 BL15   ARO2 + OH = #.108 HO2 + #.58 xHO2 + #.202 OH + &                         
        #.58 RO2C + #.11 RO2XC + #.11 zRNO3 + #.116 xGLY + &                     
        #.286 xMGLY + #.104 xBACL + #.108 CRES + &                               
        #.039 xBALD + #.217 xAFG1 + #.21 xAFG2 + &                               
        #.202 AFG3 + #.035 xPROD2 + &                                            
        #.089 yR6OOH + #.601 yRAOOH + #1.486 XC + &                              
        #.804 ARO2RO2                                                      3.773E+04   0.00  -1.000      29     9     0 = -4151    82 -4137    50                
                                                                                                                    -4350     9 -4137    90
                                                                                                                    -4089    91 -4089    75
                                                                                                                    -4268    65 -4357    66
                                                                                                                    -4304    67 -4151   101
                                                                                                                    -4027    68 -4104    69
                                                                                                                    -4090    70 -4350   118
                                                                                                                    -4022    64 -4358    77
                                                                                                                    -4359    78 -4360    98
                                                                                                                    -4361    30
 BL15a  ARO2RO2 + NO = NO + XYLNRXN                                        3.963E+03   0.71  -1.000      30     2     0 =     2    31                            
 BL15b  ARO2RO2 + HO2 = HO2 + XYLHRXN                                      2.789E+02   2.58  -1.000      30    82     0 =    82    32                            
 BL16   TERP + OH = #.759 xHO2 + #.042 xRCO3 + &                                 
        #1.147 RO2C + #.2 RO2XC + #.2 zRNO3 + #.001 xCO + &                      
        #.264 xHCHO + #.533 xRCHO + #.036 xACET + &                              
        #.005 xMEK + #.009 xMGLY + #.014 xBACL + &                               
        #.002 xMVK + #.001 xIPRD + #.255 xPROD2 + yR6OOH + &                     
        #5.056 XC + TRPRXN                                                 2.745E+04  -0.86  -1.000      39     9     0 = -4362    50 -4340    55                
                                                                                                                    -4363    90 -4005    91
                                                                                                                    -4005    75 -4185    58
                                                                                                                    -4144    59 -4364    61
                                                                                                                    -4224    62 -4297    63
                                                                                                                    -4314    66 -4348    67
                                                                                                                    -4205    72 -4185    73
                                                                                                                    -4365    64    77 -4366
                                                                                                                       98    40
 BL17   TERP + O3 = #.052 HO2 + #.067 xHO2 + #.585 OH + &                        
        #.126 xMECO3 + #.149 xRCO3 + #.875 RO2C + &                              
        #.203 RO2XC + #.203 zRNO3 + #.166 CO + #.019 xCO + &                     
        #.045 CO2 + #.079 HCHO + #.15 xHCHO + #.22 xRCHO + &                     
        #.165 xACET + #.004 MEK + #.107 HCOOH + &                                
        #.043 RCOOH + #.001 xGLY + #.002 xMGLY + &                               
        #.055 xBACL + #.001 xMACR + #.001 xIPRD + &                              
        #.409 PROD2 + #.545 yR6OOH + #3.526 XC + TRPRXN                    1.405E+00   1.56  -1.000      39     1     0 = -4367    82 -4334    50                
                                                                                                                    -4368     9 -4369    54
                                                                                                                    -4281    55 -4370    90
                                                                                                                    -4371    91 -4371    75
                                                                                                                    -4352    83 -4128    58
                                                                                                                    -4372    84 -4373    87
                                                                                                                    -4013    59 -4374    61
                                                                                                                    -4375    62 -4044   104
                                                                                                                    -4112   105 -4376    96
                                                                                                                    -4185    65 -4205    66
                                                                                                                    -4173    67 -4185    71
                                                                                                                    -4185    73 -4343   108
                                                                                                                    -4377    77 -4378    98
                                                                                                                       40
 BL18   TERP + NO3 = #.162 xHO2 + #.421 xNO2 + &                                 
        #.019 xRCO3 + #1.509 RO2C + #.397 RO2XC + &                              
        #.397 zRNO3 + #.01 xCO + #.017 xHCHO + &                                 
        #.001 xCCHO + #.509 xRCHO + #.175 xACET + &                              
        #.001 xMGLY + #.003 xMACR + #.001 xMVK + &                               
        #.002 xIPRD + #.163 xRNO3 + yR6OOH + #.416 XN + &                        
        #4.473 XC + TRPRXN                                                 1.879E+03  -0.97  -1.000      39    10     0 = -4331    50 -4379    52                
                                                                                                                    -4128    55 -4380    90
                                                                                                                    -4381    91 -4381    75
                                                                                                                    -4225    58 -4067    59
                                                                                                                    -4185    60 -4034    61
                                                                                                                    -4221    62 -4185    66
                                                                                                                    -4142    71 -4185    72
                                                                                                                    -4205    73 -4143    74
                                                                                                                       77 -4148    18 -4382
                                                                                                                       98    40
 BL19   TERP + O3P = #.147 RCHO + #.853 PROD2 + #4.441 XC + TRPRXN         5.446E+04   0.00  -1.000      39    45     0 = -4072   103 -4383   108                
                                                                                                                    -4384    98    40
 CI01   CL2 + HV = #2 CL                                                        PHOT. = CL2             140     8     0 = -4001   141                            
 CI02   CL + NO + M = CLNO + M                                             2.729E-03   0.00  -3.800     141     2     5 =    24     5                            
 CI03   CLNO + HV = CL + NO                                                     PHOT. = CLNO-06          24     8     0 =   141     2                            
 CI04   CL + NO2 = CLONO                                                  FALLOFF F= 0.600, N= 1.000    141     3     0 =    25                                  
                                                                      K0:  4.669E-02   0.00  -4.000
                                                                      KI:  1.468E+05   0.00  -2.000
 CI05   CL + NO2 = CLNO2                                                  FALLOFF F= 0.600, N= 1.000    141     3     0 =    26                                  
                                                                      K0:  6.464E-03   0.00  -4.000
                                                                      KI:  1.468E+05   0.00  -2.000
 CI06   CLONO + HV = CL + NO2                                                   PHOT. = CLONO            25     8     0 =   141     3                            
 CI07   CLNO2 + HV = CL + NO2                                                   PHOT. = CLNO2            26     8     0 =   141     3                            
 CI08   CL + HO2 = HCL + O2                                                5.050E+04   0.00  -1.560     141    82     0 =   142     4                            
 CI09   CL + HO2 = CLO + OH                                                1.381E+04   0.00   1.100     141    82     0 =   143     9                            
 CI10   CL + O3 = CLO + O2                                                 4.110E+04   0.50  -1.000     141     1     0 =   143     4                            
 CI11   CL + NO3 = CLO + NO2                                               3.523E+04   0.00  -1.000     141    10     0 =   143     3                            
 CI12   CLO + NO = CL + NO2                                                9.101E+03  -0.59  -1.000     143     2     0 =   141     3                            
 CI13   CLO + NO2 = CLONO2                                                FALLOFF F= 0.600, N= 1.000    143     3     0 =    27                                  
                                                                      K0:  6.464E-03   0.00  -5.400
                                                                      KI:  2.202E+04   0.00  -2.900
 CI14   CLONO2 + HV = CLO + NO2                                                 PHOT. = CLONO2-1         27     8     0 =   143     3                            
 CI15   CLONO2 + HV = CL + NO3                                                  PHOT. = CLONO2-2         27     8     0 =   141    10                            
 CI16   CLONO2 = CLO + NO2                                                FALLOFF F= 0.600, N= 1.000     27     0     0 =   143     3                            
                                                                      K0:  6.576E+10  24.90  -2.000
                                                                      KI:  2.226E+17  24.90   3.500
 CI17   CL + CLONO2 = CL2 + NO3                                            9.101E+03  -0.29  -1.000     141    27     0 =   140    10                            
 CI18   CLO + HO2 = HOCL + O2                                              3.229E+03  -0.68  -1.000     143    82     0 =   144     4                            
 CI19   HOCL + HV = OH + CL                                                     PHOT. = HOCL-06         144     8     0 =     9   141                            
 CI20   CLO + CLO = #.29 CL2 + #1.42 CL + O2                               1.835E+04   3.89  -1.000     143   143     0 = -4249   140 -4385   141                
                                                                                                                        4
 CI21   OH + HCL = H2O + CL                                                2.495E+03   0.46  -1.000       9   142     0 =     6   141                            
 CI22   CL + H2 = HCL + HO2                                                5.725E+04   4.59  -1.000     141     7     0 =   142    82                            
 CP01   HCHO + CL = HCL + HO2 + CO                                         1.189E+05   0.06  -1.000      87   141     0 =   142    82    83                      
 CP02   CCHO + CL = HCL + MECO3                                            1.174E+05   0.00  -1.000     102   141     0 =   142    92                            
 CP03   MEOH + CL = HCL + HCHO + HO2                                       8.073E+04   0.00  -1.000      89   141     0 =   142    87    82                      
 CP04   RCHO + CL = HCL + #.9 RCO3 + #.1 RO2C + &                                
        #.1 xCCHO + &                                                            
        #.1 xCO + #.1 xHO2 + #.1 yROOH                                     1.806E+05   0.00  -1.000     103   141     0 =   142 -4015    94 -4014                
                                                                                                                       90 -4014    60 -4014
                                                                                                                       58 -4014    50 -4014
                                                                                                                       76
 CP05   ACET + CL = HCL + RO2C + xHCHO + xMECO3 + yROOH                    1.130E+05   1.99  -1.000     100   141     0 =   142    90    59    54                
                                                                                                                       76
 CP06   MEK + CL = HCL + #.975 RO2C + #.039 RO2XC + &                            
        #.039 zRNO3 + #.84 xHO2 + #.085 xMECO3 + &                               
        #.036 xRCO3 + #.065 xHCHO + #.07 xCCHO + &                               
        #.84 xRCHO + yROOH + #.763 XC                                      5.284E+04   0.00  -1.000     104   141     0 =   142 -4158    90 -4027                
                                                                                                                       91 -4027    75 -4051
                                                                                                                       50 -4062    54 -4224
                                                                                                                       55 -4102    59 -4087
                                                                                                                       60 -4051    61    76
                                                                                                                    -4386    98
 CP07   RNO3 + CL = HCL + #.038 NO2 + #.055 HO2 + &                              
        #1.282 RO2C + #.202 RO2XC + #.202 zRNO3 + &                              
        #.009 RCHO + #.018 MEK + #.012 PROD2 + &                                 
        #.055 RNO3 + #.159 xNO2 + #.547 xHO2 + &                                 
        #.045 xHCHO + #.300 xCCHO + #.020 xRCHO + &                              
        #.003 xACET + #.041 xMEK + #.046 xPROD2 + &                              
        #.547 xRNO3 + #.908 yR6OOH + #.201 XN + #-1 XC + #.851 XC          2.818E+05   0.00  -1.000      19   141     0 =   142 -4387     3 -4173                
                                                                                                                       82 -4388    90 -4350
                                                                                                                       91 -4350    75 -4314
                                                                                                                      103 -4336   104 -4046
                                                                                                                      108 -4173    19 -4389
                                                                                                                       52 -4390    50 -4372
                                                                                                                       59 -4033    60 -4056
                                                                                                                       61 -4142    62 -4055
                                                                                                                       63 -4391    64 -4390
                                                                                                                       74 -4392    77 -4108
                                                                                                                       18 -4018    98 -4393
                                                                                                                       98
 CP08   PROD2 + CL = HCL + #.314 HO2 + #.680 RO2C + &                            
        #.116 RO2XC + #.116 zRNO3 + #.198 RCHO + &                               
        #.116 PROD2 + #.541 xHO2 + #.007 xMECO3 + &                              
        #.022 xRCO3 + #.237 xHCHO + #.109 xCCHO + &                              
        #.591 xRCHO + #.051 xMEK + #.040 xPROD2 + &                              
        #.686 yR6OOH + #1.262 XC                                           2.936E+05   0.00  -1.000     108   141     0 =   142 -4270    82 -4394                
                                                                                                                       90 -4268    91 -4268
                                                                                                                       75 -4395   103 -4268
                                                                                                                      108 -4396    50 -4136
                                                                                                                       54 -4307    55 -4397
                                                                                                                       59 -4398    60 -4237
                                                                                                                       61 -4188    63 -4084
                                                                                                                       64 -4399    77 -4400
                                                                                                                       98
 CP09   GLY + CL = HCL + #.63 HO2 + #1.26 CO + #.37 RCO3 + &                     
        #-.37 XC                                                           1.189E+05   0.06  -1.000     112   141     0 =   142 -4401    82 -4402                
                                                                                                                       83 -4248    94 -4403
                                                                                                                       98
 CP10   MGLY + CL = HCL + CO + MECO3                                       1.174E+05   0.00  -1.000     110   141     0 =   142    83    92                      
 CP11   CRES + CL = HCL + xHO2 + xBALD + yR6OOH                            9.101E+04   0.00  -1.000     101   141     0 =   142    50    68    77                
 CP12   BALD + CL = HCL + BZCO3                                            1.174E+05   0.00  -1.000     123   141     0 =   142    97                            
 CP13   ROOH + CL = HCL + #.414 OH + #.588 RO2C + &                              
        #.414 RCHO + #.104 xOH + #.482 xHO2 + &                                  
        #.106 xHCHO + #.104 xCCHO + #.197 xRCHO + &                              
        #.285 xMEK + #.586 yROOH + #-0.287 XC                              2.437E+05   0.00  -1.000     106   141     0 =   142 -4404     9 -4405                
                                                                                                                       90 -4404   103 -4304
                                                                                                                       51 -4353    50 -4284
                                                                                                                       59 -4304    60 -4406
                                                                                                                       61 -4178    63 -4407
                                                                                                                       76 -4408    98
 CP14   R6OOH + CL = HCL + #.145 OH + #1.078 RO2C + &                            
        #.117 RO2XC + #.117 zRNO3 + #.145 PROD2 + #.502 xOH + &                  
        #.237 xHO2 + #.186 xCCHO + #.676 xRCHO + &                               
        #.28 xPROD2 + #.855 yR6OOH + #.348 XC                              4.404E+05   0.00  -1.000     107   141     0 =   142 -4409     9 -4410                
                                                                                                                       90 -4123    91 -4123
                                                                                                                       75 -4409   108 -4411
                                                                                                                       51 -4397    50 -4412
                                                                                                                       60 -4413    61 -4414
                                                                                                                       64 -4415    77 -4416
                                                                                                                       98
 CP15   RAOOH + CL = #.404 HCL + #.139 OH + #.148 HO2 + &                        
        #.589 RO2C + #.124 RO2XC + #.124 zRNO3 + &                               
        #.074 PROD2 + #.147 MGLY + #.139 IPRD + &                                
        #.565 xHO2 + #.024 xOH + #.448 xRCHO + &                                 
        #.026 xGLY + #.030 xMEK + #.252 xMGLY + &                                
        #.073 xAFG1 + #.073 xAFG2 + #.713 yR6OOH + &                             
        #2.674 XC                                                          6.297E+05   0.00  -1.000     109   141     0 = -4417   142 -4068     9                
                                                                                                                    -4069    82 -4070    90
                                                                                                                    -4071    91 -4071    75
                                                                                                                    -4030   108 -4072   110
                                                                                                                    -4068   111 -4073    50
                                                                                                                    -4074    51 -4075    61
                                                                                                                    -4063    65 -4076    63
                                                                                                                    -4077    66 -4078    69
                                                                                                                    -4078    70 -4079    77
                                                                                                                    -4080    98
 CP16   MACR + CL = #.25 HCL + #.165 MACO3 + #.802 RO2C + &                      
        #.033 RO2XC + #.033 zRNO3 + #.802 xHO2 + &                               
        #.541 xCO + #.082 xIPRD + #.18 xCLCCHO + &                               
        #.541 xCLACET + #.835 yROOH + #.208 XC                             5.651E+05   0.00  -1.000     125   141     0 = -4008   142 -4375    99                
                                                                                                                    -4418    90 -4337    91
                                                                                                                    -4337    75 -4418    50
                                                                                                                    -4396    58 -4419    73
                                                                                                                    -4420    80 -4396    81
                                                                                                                    -4421    76 -4150    98
 CP17   MVK + CL = #1.283 RO2C + #.053 RO2XC + #.053 zRNO3 + &                   
        #.322 xHO2 + #.625 xMECO3 + #.947 xCLCCHO + &                            
        yROOH + #.538 XC                                                   3.406E+05   0.00  -1.000     126   141     0 = -4422    90 -4423    91                
                                                                                                                    -4423    75 -4346    50
                                                                                                                    -4424    54 -4425    80
                                                                                                                       76 -4426    98
 CP18   IPRD + CL = #.401 HCL + #.084 HO2 + #.154 MACO3 + &                      
        #.73 RO2C + #.051 RO2XC + #.051 zRNO3 + &                                
        #.042 AFG1 + #.042 AFG2 + #.712 xHO2 + #.498 xCO + &                     
        #.195 xHCHO + #.017 xMGLY + #.009 xAFG1 + &                              
        #.009 xAFG2 + #.115 xIPRD + #.14 xCLCCHO + &                             
        #.42 xCLACET + #.762 yR6OOH + #.709 XC                             6.048E+05   0.00  -1.000     111   141     0 = -4262   142 -4058    82                
                                                                                                                    -4345    99 -4427    90
                                                                                                                    -4188    91 -4188    75
                                                                                                                    -4340   113 -4340   114
                                                                                                                    -4428    50 -4179    58
                                                                                                                    -4120    59 -4067    66
                                                                                                                    -4314    69 -4314    70
                                                                                                                    -4085    73 -4180    80
                                                                                                                    -4429    81 -4430    77
                                                                                                                    -4349    98
 CP19   CLCCHO + HV = HO2 + CO + RO2C + xCL + xHCHO + &                          
        yROOH                                                                   PHOT. = CLCCHO          145     8     0 =    82    83    90    79                
                                                                                                                       59    76
 CP20   CLCCHO + OH = RCO3 + #-1 XC                                        4.550E+03   0.00  -1.000     145     9     0 =    94 -4018    98                      
 CP21   CLCCHO + CL = HCL + RCO3 + #-1 XC                                  1.894E+04   0.00  -1.000     145   141     0 =   142    94 -4018    98                
 CP22   CLACET + HV + #QY.CLAC = MECO3 + RO2C + xCL + xHCHO + yROOH             PHOT. = CLACET          146     8   -17 =    92    90    79    59                
                                                                                                                       76
 CP23   xCL + #RO2RO = CL                                                  1.000E+00   0.00   0.000      79    -5     0 =   141                                  
 CP24   xCL + #RO2XRO =                                                    1.000E+00   0.00   0.000      79    -6     0 =                                        
 CP25   xCLCCHO + #RO2RO = CLCCHO                                          1.000E+00   0.00   0.000      80    -5     0 =   145                                  
 CP26   xCLCCHO + #RO2XRO = #2 XC                                          1.000E+00   0.00   0.000      80    -6     0 = -4001    98                            
 CP27   xCLACET + #RO2RO = CLACET                                          1.000E+00   0.00   0.000      81    -5     0 =   146                                  
 CP28   xCLACET + #RO2XRO = #3 XC                                          1.000E+00   0.00   0.000      81    -6     0 = -4019    98                            
 CE01   CH4 + CL = HCL + MEO2                                              1.072E+04   2.54  -1.000     130   141     0 =   142    86                            
 CE02   ETHENE + CL = #2 RO2C + xHO2 + xHCHO + CLCHO                      FALLOFF F= 0.600, N= 1.000    131   141     0 = -4001    90    50    59                
                                                                                                                      147
                                                                      K0:  5.746E-01   0.00  -5.300
                                                                      KI:  4.550E+05   0.00  -2.000
 CE03   ISOPRENE + CL = #.15 HCL + #1.168 RO2C + &                               
        #.085 RO2XC + #.085 zRNO3 + #.738 xHO2 + &                               
        #.177 xCL + #.275 xHCHO + #.177 xMVK + &                                 
        #.671 xIPRD + #.067 xCLCCHO + yR6OOH + #.018 XC                    7.046E+05   0.00  -1.000      41   141     0 = -4013   142 -4431    90                
                                                                                                                    -4062    91 -4062    75
                                                                                                                    -4432    50 -4433    79
                                                                                                                    -4259    59 -4433    72
                                                                                                                    -4434    73 -4334    80
                                                                                                                       77 -4336    98
 CE04   ACETYLEN + CL = HO2 + CO + XC                                     FALLOFF F= 0.600, N= 1.000    132   141     0 =    82    83    98                      
                                                                      K0:  1.867E-01   0.00  -4.400
                                                                      KI:  3.229E+05   0.00  -1.000
 BC01   ALK1 + CL = xHO2 + RO2C + HCL + xCCHO + yROOH                      1.218E+05   0.20  -1.000     133   141     0 =    50    90   142    60                
                                                                                                                       76
 BC02   ALK2 + CL = #.97 xHO2 + #.97 RO2C + #.03 RO2XC + &                       
        #.03 zRNO3 + HCL + #.482 xRCHO + #.488 xACET + &                         
        yROOH + #-.09 XC                                                   1.761E+05  -0.08  -1.000     134   141     0 = -4435    50 -4435    90                
                                                                                                                    -4076    91 -4076    75
                                                                                                                      142 -4353    61 -4436
                                                                                                                       62    76 -4437    98
 BC03   ALK3 + CL = #.835 xHO2 + #.094 xTBUO + &                                 
        #1.361 RO2C + #.07 RO2XC + #.07 zRNO3 + HCL + &                          
        #.078 xHCHO + #.34 xCCHO + #.343 xRCHO + &                               
        #.075 xACET + #.253 xMEK + #.983 yROOH + &                               
        #.017 yR6OOH + #.18 XC                                             2.730E+05   0.00  -1.000     135   141     0 = -4421    50 -4438    57                
                                                                                                                    -4439    90 -4087    91
                                                                                                                    -4087    75   142 -4440
                                                                                                                       59 -4157    60 -4441
                                                                                                                       61 -4057    62 -4442
                                                                                                                       63 -4277    76 -4067
                                                                                                                       77 -4420    98
 BC04   ALK4 + CL = #.827 xHO2 + #.003 xMEO2 + &                                 
        #.004 xMECO3 + #1.737 RO2C + #.165 RO2XC + &                             
        #.165 zRNO3 + HCL + #.003 xCO + #.034 xHCHO + &                          
        #.287 xCCHO + #.412 xRCHO + #.247 xACET + &                              
        #.076 xMEK + #.13 xPROD2 + yR6OOH + #.327 XC                       3.861E+05   0.00  -1.000     136   141     0 = -4443    50 -4142    53                
                                                                                                                    -4044    54 -4444    90
                                                                                                                    -4375    91 -4375    75
                                                                                                                      142 -4142    58 -4048
                                                                                                                       59 -4445    60 -4446
                                                                                                                       61 -4447    62 -4325
                                                                                                                       63 -4083    64    77
                                                                                                                    -4448    98
 BC05   ALK5 + CL = #.647 xHO2 + #1.541 RO2C + &                                 
        #.352 RO2XC + #.352 zRNO3 + HCL + #.022 xHCHO + &                        
        #.08 xCCHO + #.258 xRCHO + #.044 xACET + &                               
        #.041 xMEK + #.378 xPROD2 + yR6OOH + #2.368 XC                     6.180E+05   0.00  -1.000      28   141     0 = -4287    50 -4449    90                
                                                                                                                    -4450    91 -4450    75
                                                                                                                      142 -4307    59 -4451
                                                                                                                       60 -4135    61 -4452
                                                                                                                       62 -4055    63 -4453
                                                                                                                       64    77 -4454    98
 BC06   OLE1 + CL = #.902 xHO2 + #1.42 RO2C + &                                  
        #.098 RO2XC + #.098 zRNO3 + #.308 HCL + &                                
        #.025 xHCHO + #.146 xCCHO + #.051 xRCHO + &                              
        #.188 xMACR + #.014 xMVK + #.027 xIPRD + &                               
        #.225 xCLCCHO + #.396 xCLACET + #.413 yROOH + &                          
        #.587 yR6OOH + #1.361 XC                                           5.211E+05   0.00  -1.000     138   141     0 = -4146    50 -4385    90                
                                                                                                                    -4147    91 -4147    75
                                                                                                                    -4455   142 -4159    59
                                                                                                                    -4285    60 -4188    61
                                                                                                                    -4456    71 -4348    72
                                                                                                                    -4326    73 -4457    80
                                                                                                                    -4238    81 -4300    76
                                                                                                                    -4301    77 -4439    98
 BC07   OLE2 + CL = #.447 xHO2 + #.448 xCL + #.001 xMEO2 + &                     
        #1.514 RO2C + #.104 RO2XC + #.104 zRNO3 + &                              
        #.263 HCL + #.228 xHCHO + #.361 xCCHO + &                                
        #.3 xRCHO + #.081 xACET + #.04 xMEK + &                                  
        #.049 xMACR + #.055 xMVK + #.179 xIPRD + &                               
        #.012 xCLCCHO + #.18 xCLACET + #.357 yROOH + &                           
        #.643 yR6OOH + #.247 XC                                            5.622E+05   0.00  -1.000     139   141     0 = -4458    50 -4075    79                
                                                                                                                    -4185    53 -4459    90
                                                                                                                    -4304    91 -4304    75
                                                                                                                    -4460   142 -4461    59
                                                                                                                    -4462    60 -4033    61
                                                                                                                    -4239    62 -4084    63
                                                                                                                    -4202    71 -4173    72
                                                                                                                    -4463    73 -4046    80
                                                                                                                    -4420    81 -4255    76
                                                                                                                    -4327    77 -4447    98
 BC08   ARO1 + CL = #.88 xHO2 + #.88 RO2C + #.12 RO2XC + &                       
        #.12 zRNO3 + #.671 xBALD + #.21 xPROD2 + #.323 XC                  1.468E+05   0.00  -1.000      33   141     0 = -4464    50 -4464    90                
                                                                                                                    -4247    91 -4247    75
                                                                                                                    -4434    68 -4090    64
                                                                                                                    -4465    98
 BC09   ARO2 + CL = #.842 xHO2 + #.842 RO2C + &                                  
        #.158 RO2XC + #.158 zRNO3 + #.618 xBALD + &                              
        #.224 xPROD2 + #2.382 XC                                           3.200E+05   0.00  -1.000      29   141     0 = -4466    50 -4466    90                
                                                                                                                    -4467    91 -4467    75
                                                                                                                    -4468    68 -4469    64
                                                                                                                    -4470    98
 BC10   TERP + CL = #.252 xHO2 + #.068 xCL + &                                   
        #.034 xMECO3 + #.05 xRCO3 + #.016 xMACO3 + &                             
        #2.258 RO2C + #.582 RO2XC + #.582 zRNO3 + &                              
        #.548 HCL + #.035 xCO + #.158 xHCHO + &                                  
        #.185 xRCHO + #.274 xACET + #.007 xGLY + &                               
        #.003 xBACL + #.003 xMVK + #.158 xIPRD + &                               
        #.006 xAFG1 + #.006 xAFG2 + &                                            
        #.109 xCLCCHO + yR6OOH + #3.543 XC                                 8.015E+05   0.00  -1.000      39   141     0 = -4077    50 -4354    79                
                                                                                                                    -4048    54 -4164    55
                                                                                                                    -4471    56 -4472    90
                                                                                                                    -4473    91 -4473    75
                                                                                                                    -4474   142 -4022    58
                                                                                                                    -4467    59 -4195    61
                                                                                                                    -4126    62 -4136    65
                                                                                                                    -4142    67 -4142    72
                                                                                                                    -4467    73 -4298    69
                                                                                                                    -4298    70 -4398    80
                                                                                                                       77 -4475    98
 D1       HCL = HCL                                                        1.000E-25   0.00   0.000     142     0     0 =   142                                  
 D2       NH3 = NH3                                                        1.000E-25   0.00   0.000     148     0     0 =   148                                  
 D3       GDENS = GDENS                                                    1.000E-25   0.00   0.000     149     0     0 =   149                                  
 D4       MMOLES = MMOLES                                                  1.000E-25   0.00   0.000     150     0     0 =   150                                  
 D5       XICBC = XICBC                                                    1.000E-25   0.00   0.000     151     0     0 =   151                                  
 D6       YICBC = YICBC                                                    1.000E-25   0.00   0.000     152     0     0 =   152                                  
 D7       ZICBC = ZICBC                                                    1.000E-25   0.00   0.000     153     0     0 =   153                                  
 TS01   CHCL3 + OH = xCL + RO2C + yROOH + XC                               8.323E+02   1.00   1.000     154     9     0 =    79    90    76    98                
 TS02   CHCL3 + CL = xCL + RO2C + HCL + yROOH + XC                         4.844E+03   1.97  -1.000     154   141     0 =    79    90   142    76                
                                                                                                                       98
0 364 REACTIONS,  154 SPECIES,   17 VARIABLE COEFFICIENTS, 475 CONSTANT COEFFICIENTS, AND  35 PHOTOLYSIS FILES.


These are the som reactions

 BL14   ARO1 + OH = #0.166 HO2 + #0.482 xHO2 + #0.284 OH + #0.482 RO2C + &
          #0.068 RO2XC + #0.068 zRNO3 + #0.218 xGLY + #0.138 xMGLY + #0.166 &
          CRES + #0.049 xBALD + #0.164 xAFG1 + #0.193 xAFG2 + #0.284 AFG3 + &
          #0.077 xPROD2 + #0.403 yRAOOH + #0.147 yR6OOH + #0.002 XC + #1.000 &
          #0.123 AR1SOMG_07_01 + #1.000 #0.001 AR1SOMG_07_02 + #1.000 #0.002 &
          AR1SOMG_07_03 + #1.000 #0.874 AR1SOMG_07_04
 S1.1   AR1SOMG_02_01 + OH = #0.969 #0.123 AR1SOMG_02_02 + #0.969 #0.001 &
          AR1SOMG_02_03 + #0.969 #0.876 AR1SOMG_02_04 + #0.031 #1.000 AR1SOMG_01_01 + &
          #0.031 #1.000 AR1SOMG_01_02
 S2.1   AR1SOMG_02_02 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #2.000 &
          AR1SOMG_01_02
 S3.1   AR1SOMG_02_03 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_01_02
 S4.1   AR1SOMG_02_04 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_01_02
 S5.1   AR1SOMG_03_01 + OH = #0.996 #0.123 AR1SOMG_03_02 + #0.996 #0.001 &
          AR1SOMG_03_03 + #0.996 #0.002 AR1SOMG_03_04 + #0.996 #0.874 AR1SOMG_03_05 + &
          #0.004 #0.500 AR1SOMG_01_01 + #0.004 #0.500 AR1SOMG_02_01 + #0.004 &
          #0.500 AR1SOMG_01_02 + #0.004 #0.500 AR1SOMG_02_02
 S6.1   AR1SOMG_03_02 + OH = #0.868 #0.123 AR1SOMG_03_03 + #0.868 #0.001 &
          AR1SOMG_03_04 + #0.868 #0.002 AR1SOMG_03_05 + #0.868 #0.874 AR1SOMG_03_06 + &
          #0.132 #0.500 AR1SOMG_01_01 + #0.132 #0.000 AR1SOMG_02_01 + #0.132 &
          #0.500 AR1SOMG_01_02 + #0.132 #0.500 AR1SOMG_02_02 + #0.132 #0.500 &
          AR1SOMG_02_03
 S7.1   AR1SOMG_03_03 + OH = #1.000 #0.400 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.400 AR1SOMG_01_02 + #1.000 #0.000 AR1SOMG_02_02 + &
          #1.000 #0.400 AR1SOMG_02_03 + #1.000 #0.400 AR1SOMG_02_04
 S8.1   AR1SOMG_03_04 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.500 AR1SOMG_01_02 + #1.000 #0.000 AR1SOMG_02_02 + &
          #1.000 #0.000 AR1SOMG_02_03 + #1.000 #0.500 AR1SOMG_02_04
 S9.1   AR1SOMG_03_05 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_01_02 + #1.000 #0.000 AR1SOMG_02_02 + &
          #1.000 #0.000 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_02_04
 S10.1  AR1SOMG_03_06 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_01_02 + #1.000 #0.000 AR1SOMG_02_02 + &
          #1.000 #0.000 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_02_04
 S11.1  AR1SOMG_04_01 + OH = #0.999 #0.123 AR1SOMG_04_02 + #0.999 #0.001 &
          AR1SOMG_04_03 + #0.999 #0.002 AR1SOMG_04_04 + #0.999 #0.874 AR1SOMG_04_05 + &
          #0.001 #0.333 AR1SOMG_01_01 + #0.001 #0.333 AR1SOMG_02_01 + #0.001 &
          #0.333 AR1SOMG_03_01 + #0.001 #0.333 AR1SOMG_01_02 + #0.001 #0.333 &
          AR1SOMG_02_02 + #0.001 #0.333 AR1SOMG_03_02
 S12.1  AR1SOMG_04_02 + OH = #0.969 #0.123 AR1SOMG_04_03 + #0.969 #0.001 &
          AR1SOMG_04_04 + #0.969 #0.002 AR1SOMG_04_05 + #0.969 #0.874 AR1SOMG_04_06 + &
          #0.031 #0.286 AR1SOMG_01_01 + #0.031 #0.286 AR1SOMG_02_01 + #0.031 &
          #0.000 AR1SOMG_03_01 + #0.031 #0.286 AR1SOMG_01_02 + #0.031 #0.286 &
          AR1SOMG_02_02 + #0.031 #0.286 AR1SOMG_03_02 + #0.031 #0.286 AR1SOMG_02_03 + &
          #0.031 #0.286 AR1SOMG_03_03
 S13.1  AR1SOMG_04_03 + OH = #0.763 #0.123 AR1SOMG_04_04 + #0.763 #0.001 &
          AR1SOMG_04_05 + #0.763 #0.002 AR1SOMG_04_06 + #0.763 #0.874 AR1SOMG_04_07 + &
          #0.237 #0.222 AR1SOMG_01_01 + #0.237 #0.222 AR1SOMG_02_01 + #0.237 &
          #0.000 AR1SOMG_03_01 + #0.237 #0.222 AR1SOMG_01_02 + #0.237 #0.222 &
          AR1SOMG_02_02 + #0.237 #0.000 AR1SOMG_03_02 + #0.237 #0.222 AR1SOMG_02_03 + &
          #0.237 #0.222 AR1SOMG_03_03 + #0.237 #0.222 AR1SOMG_02_04 + #0.237 &
          #0.222 AR1SOMG_03_04
 S14.1  AR1SOMG_04_04 + OH = #1.000 #0.222 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.222 AR1SOMG_01_02 + &
          #1.000 #0.222 AR1SOMG_02_02 + #1.000 #0.000 AR1SOMG_03_02 + #1.000 &
          #0.222 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.222 &
          AR1SOMG_02_04 + #1.000 #0.222 AR1SOMG_03_04 + #1.000 #0.222 AR1SOMG_03_05
 S15.1  AR1SOMG_04_05 + OH = #1.000 #0.250 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.250 AR1SOMG_01_02 + &
          #1.000 #0.000 AR1SOMG_02_02 + #1.000 #0.000 AR1SOMG_03_02 + #1.000 &
          #0.250 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.250 &
          AR1SOMG_02_04 + #1.000 #0.000 AR1SOMG_03_04 + #1.000 #0.250 AR1SOMG_03_05 + &
          #1.000 #0.250 AR1SOMG_03_06
 S16.1  AR1SOMG_04_06 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.400 AR1SOMG_01_02 + &
          #1.000 #0.000 AR1SOMG_02_02 + #1.000 #0.000 AR1SOMG_03_02 + #1.000 &
          #0.000 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.400 &
          AR1SOMG_02_04 + #1.000 #0.000 AR1SOMG_03_04 + #1.000 #0.000 AR1SOMG_03_05 + &
          #1.000 #0.400 AR1SOMG_03_06
 S17.1  AR1SOMG_04_07 + OH = #1.000 #0.000 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_01_02 + &
          #1.000 #0.000 AR1SOMG_02_02 + #1.000 #0.000 AR1SOMG_03_02 + #1.000 &
          #0.000 AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.000 &
          AR1SOMG_02_04 + #1.000 #0.000 AR1SOMG_03_04 + #1.000 #0.000 AR1SOMG_03_05 + &
          #1.000 #0.000 AR1SOMG_03_06
 S18.1  AR1SOMG_05_01 + OH = #1.000 #0.123 AR1SOMG_05_02 + #1.000 #0.001 &
          AR1SOMG_05_03 + #1.000 #0.002 AR1SOMG_05_04 + #1.000 #0.874 AR1SOMG_05_05 + &
          #0.000 #0.250 AR1SOMG_01_01 + #0.000 #0.250 AR1SOMG_02_01 + #0.000 &
          #0.250 AR1SOMG_03_01 + #0.000 #0.250 AR1SOMG_04_01 + #0.000 #0.250 &
          AR1SOMG_01_02 + #0.000 #0.250 AR1SOMG_02_02 + #0.000 #0.250 AR1SOMG_03_02 + &
          #0.000 #0.250 AR1SOMG_04_02
 S19.1  AR1SOMG_05_02 + OH = #0.990 #0.123 AR1SOMG_05_03 + #0.990 #0.001 &
          AR1SOMG_05_04 + #0.990 #0.002 AR1SOMG_05_05 + #0.990 #0.874 AR1SOMG_05_06 + &
          #0.010 #0.200 AR1SOMG_01_01 + #0.010 #0.200 AR1SOMG_02_01 + #0.010 &
          #0.200 AR1SOMG_03_01 + #0.010 #0.000 AR1SOMG_04_01 + #0.010 #0.200 &
          AR1SOMG_01_02 + #0.010 #0.200 AR1SOMG_02_02 + #0.010 #0.200 AR1SOMG_03_02 + &
          #0.010 #0.200 AR1SOMG_04_02 + #0.010 #0.200 AR1SOMG_02_03 + #0.010 &
          #0.200 AR1SOMG_03_03 + #0.010 #0.200 AR1SOMG_04_03
 S20.1  AR1SOMG_05_03 + OH = #0.922 #0.123 AR1SOMG_05_04 + #0.922 #0.001 &
          AR1SOMG_05_05 + #0.922 #0.002 AR1SOMG_05_06 + #0.922 #0.874 AR1SOMG_05_07 + &
          #0.078 #0.154 AR1SOMG_01_01 + #0.078 #0.154 AR1SOMG_02_01 + #0.078 &
          #0.154 AR1SOMG_03_01 + #0.078 #0.000 AR1SOMG_04_01 + #0.078 #0.154 &
          AR1SOMG_01_02 + #0.078 #0.154 AR1SOMG_02_02 + #0.078 #0.154 AR1SOMG_03_02 + &
          #0.078 #0.000 AR1SOMG_04_02 + #0.078 #0.154 AR1SOMG_02_03 + #0.078 &
          #0.154 AR1SOMG_03_03 + #0.078 #0.154 AR1SOMG_04_03 + #0.078 #0.154 &
          AR1SOMG_02_04 + #0.078 #0.154 AR1SOMG_03_04 + #0.078 #0.154 AR1SOMG_04_04
 S21.1  AR1SOMG_05_04 + OH = #0.672 #0.123 AR1SOMG_05_05 + #0.672 #0.001 &
          AR1SOMG_05_06 + #0.672 #0.876 AR1SOMG_05_07 + #0.328 #0.143 AR1SOMG_01_01 + &
          #0.328 #0.143 AR1SOMG_02_01 + #0.328 #0.000 AR1SOMG_03_01 + #0.328 &
          #0.000 AR1SOMG_04_01 + #0.328 #0.143 AR1SOMG_01_02 + #0.328 #0.143 &
          AR1SOMG_02_02 + #0.328 #0.143 AR1SOMG_03_02 + #0.328 #0.000 AR1SOMG_04_02 + &
          #0.328 #0.143 AR1SOMG_02_03 + #0.328 #0.143 AR1SOMG_03_03 + #0.328 &
          #0.000 AR1SOMG_04_03 + #0.328 #0.143 AR1SOMG_02_04 + #0.328 #0.143 &
          AR1SOMG_03_04 + #0.328 #0.143 AR1SOMG_04_04 + #0.328 #0.143 AR1SOMG_03_05 + &
          #0.328 #0.143 AR1SOMG_04_05
 S22.1  AR1SOMG_05_05 + OH = #1.000 #0.143 AR1SOMG_01_01 + #1.000 #0.143 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.143 AR1SOMG_01_02 + #1.000 #0.143 AR1SOMG_02_02 + #1.000 &
          #0.000 AR1SOMG_03_02 + #1.000 #0.000 AR1SOMG_04_02 + #1.000 #0.143 &
          AR1SOMG_02_03 + #1.000 #0.143 AR1SOMG_03_03 + #1.000 #0.000 AR1SOMG_04_03 + &
          #1.000 #0.143 AR1SOMG_02_04 + #1.000 #0.143 AR1SOMG_03_04 + #1.000 &
          #0.000 AR1SOMG_04_04 + #1.000 #0.143 AR1SOMG_03_05 + #1.000 #0.143 &
          AR1SOMG_04_05 + #1.000 #0.143 AR1SOMG_03_06 + #1.000 #0.143 AR1SOMG_04_06
 S23.1  AR1SOMG_05_06 + OH = #1.000 #0.167 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.167 AR1SOMG_01_02 + #1.000 #0.167 AR1SOMG_02_02 + #1.000 &
          #0.000 AR1SOMG_03_02 + #1.000 #0.000 AR1SOMG_04_02 + #1.000 #0.167 &
          AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.000 AR1SOMG_04_03 + &
          #1.000 #0.167 AR1SOMG_02_04 + #1.000 #0.167 AR1SOMG_03_04 + #1.000 &
          #0.000 AR1SOMG_04_04 + #1.000 #0.167 AR1SOMG_03_05 + #1.000 #0.000 &
          AR1SOMG_04_05 + #1.000 #0.167 AR1SOMG_03_06 + #1.000 #0.167 AR1SOMG_04_06 + &
          #1.000 #0.167 AR1SOMG_04_07
 S24.1  AR1SOMG_05_07 + OH = #1.000 #0.222 AR1SOMG_01_01 + #1.000 #0.000 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.222 AR1SOMG_01_02 + #1.000 #0.000 AR1SOMG_02_02 + #1.000 &
          #0.000 AR1SOMG_03_02 + #1.000 #0.000 AR1SOMG_04_02 + #1.000 #0.222 &
          AR1SOMG_02_03 + #1.000 #0.000 AR1SOMG_03_03 + #1.000 #0.000 AR1SOMG_04_03 + &
          #1.000 #0.222 AR1SOMG_02_04 + #1.000 #0.000 AR1SOMG_03_04 + #1.000 &
          #0.000 AR1SOMG_04_04 + #1.000 #0.222 AR1SOMG_03_05 + #1.000 #0.000 &
          AR1SOMG_04_05 + #1.000 #0.222 AR1SOMG_03_06 + #1.000 #0.000 AR1SOMG_04_06 + &
          #1.000 #0.222 AR1SOMG_04_07
 S25.1  AR1SOMG_06_01 + OH = #1.000 #0.123 AR1SOMG_06_02 + #1.000 #0.001 &
          AR1SOMG_06_03 + #1.000 #0.002 AR1SOMG_06_04 + #1.000 #0.874 AR1SOMG_06_05 + &
          #0.000 #0.200 AR1SOMG_01_01 + #0.000 #0.200 AR1SOMG_02_01 + #0.000 &
          #0.200 AR1SOMG_03_01 + #0.000 #0.200 AR1SOMG_04_01 + #0.000 #0.200 &
          AR1SOMG_05_01 + #0.000 #0.200 AR1SOMG_01_02 + #0.000 #0.200 AR1SOMG_02_02 + &
          #0.000 #0.200 AR1SOMG_03_02 + #0.000 #0.200 AR1SOMG_04_02 + #0.000 &
          #0.200 AR1SOMG_05_02
 S26.1  AR1SOMG_06_02 + OH = #0.996 #0.123 AR1SOMG_06_03 + #0.996 #0.001 &
          AR1SOMG_06_04 + #0.996 #0.002 AR1SOMG_06_05 + #0.996 #0.874 AR1SOMG_06_06 + &
          #0.004 #0.154 AR1SOMG_01_01 + #0.004 #0.154 AR1SOMG_02_01 + #0.004 &
          #0.154 AR1SOMG_03_01 + #0.004 #0.154 AR1SOMG_04_01 + #0.004 #0.000 &
          AR1SOMG_05_01 + #0.004 #0.154 AR1SOMG_01_02 + #0.004 #0.154 AR1SOMG_02_02 + &
          #0.004 #0.154 AR1SOMG_03_02 + #0.004 #0.154 AR1SOMG_04_02 + #0.004 &
          #0.154 AR1SOMG_05_02 + #0.004 #0.154 AR1SOMG_02_03 + #0.004 #0.154 &
          AR1SOMG_03_03 + #0.004 #0.154 AR1SOMG_04_03 + #0.004 #0.154 AR1SOMG_05_03
 S27.1  AR1SOMG_06_03 + OH = #0.969 #0.123 AR1SOMG_06_04 + #0.969 #0.001 &
          AR1SOMG_06_05 + #0.969 #0.002 AR1SOMG_06_06 + #0.969 #0.874 AR1SOMG_06_07 + &
          #0.031 #0.118 AR1SOMG_01_01 + #0.031 #0.118 AR1SOMG_02_01 + #0.031 &
          #0.118 AR1SOMG_03_01 + #0.031 #0.118 AR1SOMG_04_01 + #0.031 #0.000 &
          AR1SOMG_05_01 + #0.031 #0.118 AR1SOMG_01_02 + #0.031 #0.118 AR1SOMG_02_02 + &
          #0.031 #0.118 AR1SOMG_03_02 + #0.031 #0.118 AR1SOMG_04_02 + #0.031 &
          #0.000 AR1SOMG_05_02 + #0.031 #0.118 AR1SOMG_02_03 + #0.031 #0.118 &
          AR1SOMG_03_03 + #0.031 #0.118 AR1SOMG_04_03 + #0.031 #0.118 AR1SOMG_05_03 + &
          #0.031 #0.118 AR1SOMG_02_04 + #0.031 #0.118 AR1SOMG_03_04 + #0.031 &
          #0.118 AR1SOMG_04_04 + #0.031 #0.118 AR1SOMG_05_04
 S28.1  AR1SOMG_06_04 + OH = #0.868 #0.123 AR1SOMG_06_05 + #0.868 #0.001 &
          AR1SOMG_06_06 + #0.868 #0.876 AR1SOMG_06_07 + #0.132 #0.105 AR1SOMG_01_01 + &
          #0.132 #0.105 AR1SOMG_02_01 + #0.132 #0.105 AR1SOMG_03_01 + #0.132 &
          #0.000 AR1SOMG_04_01 + #0.132 #0.000 AR1SOMG_05_01 + #0.132 #0.105 &
          AR1SOMG_01_02 + #0.132 #0.105 AR1SOMG_02_02 + #0.132 #0.105 AR1SOMG_03_02 + &
          #0.132 #0.105 AR1SOMG_04_02 + #0.132 #0.000 AR1SOMG_05_02 + #0.132 &
          #0.105 AR1SOMG_02_03 + #0.132 #0.105 AR1SOMG_03_03 + #0.132 #0.105 &
          AR1SOMG_04_03 + #0.132 #0.000 AR1SOMG_05_03 + #0.132 #0.105 AR1SOMG_02_04 + &
          #0.132 #0.105 AR1SOMG_03_04 + #0.132 #0.105 AR1SOMG_04_04 + #0.132 &
          #0.105 AR1SOMG_05_04 + #0.132 #0.105 AR1SOMG_03_05 + #0.132 #0.105 &
          AR1SOMG_04_05 + #0.132 #0.105 AR1SOMG_05_05
 S29.1  AR1SOMG_06_05 + OH = #0.598 #0.123 AR1SOMG_06_06 + #0.598 #0.877 &
          AR1SOMG_06_07 + #0.402 #0.100 AR1SOMG_01_01 + #0.402 #0.100 AR1SOMG_02_01 + &
          #0.402 #0.100 AR1SOMG_03_01 + #0.402 #0.000 AR1SOMG_04_01 + #0.402 &
          #0.000 AR1SOMG_05_01 + #0.402 #0.100 AR1SOMG_01_02 + #0.402 #0.100 &
          AR1SOMG_02_02 + #0.402 #0.100 AR1SOMG_03_02 + #0.402 #0.000 AR1SOMG_04_02 + &
          #0.402 #0.000 AR1SOMG_05_02 + #0.402 #0.100 AR1SOMG_02_03 + #0.402 &
          #0.100 AR1SOMG_03_03 + #0.402 #0.100 AR1SOMG_04_03 + #0.402 #0.000 &
          AR1SOMG_05_03 + #0.402 #0.100 AR1SOMG_02_04 + #0.402 #0.100 AR1SOMG_03_04 + &
          #0.402 #0.100 AR1SOMG_04_04 + #0.402 #0.000 AR1SOMG_05_04 + #0.402 &
          #0.100 AR1SOMG_03_05 + #0.402 #0.100 AR1SOMG_04_05 + #0.402 #0.100 &
          AR1SOMG_05_05 + #0.402 #0.100 AR1SOMG_03_06 + #0.402 #0.100 AR1SOMG_04_06 + &
          #0.402 #0.100 AR1SOMG_05_06
 S30.1  AR1SOMG_06_06 + OH = #1.000 #0.105 AR1SOMG_01_01 + #1.000 #0.105 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.000 AR1SOMG_05_01 + #1.000 #0.105 AR1SOMG_01_02 + #1.000 &
          #0.105 AR1SOMG_02_02 + #1.000 #0.105 AR1SOMG_03_02 + #1.000 #0.000 &
          AR1SOMG_04_02 + #1.000 #0.000 AR1SOMG_05_02 + #1.000 #0.105 AR1SOMG_02_03 + &
          #1.000 #0.105 AR1SOMG_03_03 + #1.000 #0.000 AR1SOMG_04_03 + #1.000 &
          #0.000 AR1SOMG_05_03 + #1.000 #0.105 AR1SOMG_02_04 + #1.000 #0.105 &
          AR1SOMG_03_04 + #1.000 #0.105 AR1SOMG_04_04 + #1.000 #0.000 AR1SOMG_05_04 + &
          #1.000 #0.105 AR1SOMG_03_05 + #1.000 #0.105 AR1SOMG_04_05 + #1.000 &
          #0.000 AR1SOMG_05_05 + #1.000 #0.105 AR1SOMG_03_06 + #1.000 #0.105 &
          AR1SOMG_04_06 + #1.000 #0.105 AR1SOMG_05_06 + #1.000 #0.105 AR1SOMG_04_07 + &
          #1.000 #0.105 AR1SOMG_05_07
 S31.1  AR1SOMG_06_07 + OH = #1.000 #0.125 AR1SOMG_01_01 + #1.000 #0.125 &
          AR1SOMG_02_01 + #1.000 #0.000 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.000 AR1SOMG_05_01 + #1.000 #0.125 AR1SOMG_01_02 + #1.000 &
          #0.125 AR1SOMG_02_02 + #1.000 #0.000 AR1SOMG_03_02 + #1.000 #0.000 &
          AR1SOMG_04_02 + #1.000 #0.000 AR1SOMG_05_02 + #1.000 #0.125 AR1SOMG_02_03 + &
          #1.000 #0.125 AR1SOMG_03_03 + #1.000 #0.000 AR1SOMG_04_03 + #1.000 &
          #0.000 AR1SOMG_05_03 + #1.000 #0.125 AR1SOMG_02_04 + #1.000 #0.125 &
          AR1SOMG_03_04 + #1.000 #0.000 AR1SOMG_04_04 + #1.000 #0.000 AR1SOMG_05_04 + &
          #1.000 #0.125 AR1SOMG_03_05 + #1.000 #0.125 AR1SOMG_04_05 + #1.000 &
          #0.000 AR1SOMG_05_05 + #1.000 #0.125 AR1SOMG_03_06 + #1.000 #0.125 &
          AR1SOMG_04_06 + #1.000 #0.000 AR1SOMG_05_06 + #1.000 #0.125 AR1SOMG_04_07 + &
          #1.000 #0.125 AR1SOMG_05_07
 S32.1  AR1SOMG_07_01 + OH = #1.000 #0.123 AR1SOMG_07_02 + #1.000 #0.001 &
          AR1SOMG_07_03 + #1.000 #0.002 AR1SOMG_07_04 + #1.000 #0.874 AR1SOMG_07_05 + &
          #0.000 #0.167 AR1SOMG_01_01 + #0.000 #0.167 AR1SOMG_02_01 + #0.000 &
          #0.167 AR1SOMG_03_01 + #0.000 #0.167 AR1SOMG_04_01 + #0.000 #0.167 &
          AR1SOMG_05_01 + #0.000 #0.167 AR1SOMG_06_01 + #0.000 #0.167 AR1SOMG_01_02 + &
          #0.000 #0.167 AR1SOMG_02_02 + #0.000 #0.167 AR1SOMG_03_02 + #0.000 &
          #0.167 AR1SOMG_04_02 + #0.000 #0.167 AR1SOMG_05_02 + #0.000 #0.167 &
          AR1SOMG_06_02
 S33.1  AR1SOMG_07_02 + OH = #0.998 #0.123 AR1SOMG_07_03 + #0.998 #0.001 &
          AR1SOMG_07_04 + #0.998 #0.002 AR1SOMG_07_05 + #0.998 #0.874 AR1SOMG_07_06 + &
          #0.002 #0.125 AR1SOMG_01_01 + #0.002 #0.125 AR1SOMG_02_01 + #0.002 &
          #0.125 AR1SOMG_03_01 + #0.002 #0.125 AR1SOMG_04_01 + #0.002 #0.125 &
          AR1SOMG_05_01 + #0.002 #0.000 AR1SOMG_06_01 + #0.002 #0.125 AR1SOMG_01_02 + &
          #0.002 #0.125 AR1SOMG_02_02 + #0.002 #0.125 AR1SOMG_03_02 + #0.002 &
          #0.125 AR1SOMG_04_02 + #0.002 #0.125 AR1SOMG_05_02 + #0.002 #0.125 &
          AR1SOMG_06_02 + #0.002 #0.125 AR1SOMG_02_03 + #0.002 #0.125 AR1SOMG_03_03 + &
          #0.002 #0.125 AR1SOMG_04_03 + #0.002 #0.125 AR1SOMG_05_03 + #0.002 &
          #0.125 AR1SOMG_06_03
 S34.1  AR1SOMG_07_03 + OH = #0.986 #0.123 AR1SOMG_07_04 + #0.986 #0.001 &
          AR1SOMG_07_05 + #0.986 #0.002 AR1SOMG_07_06 + #0.986 #0.874 AR1SOMG_07_07 + &
          #0.014 #0.095 AR1SOMG_01_01 + #0.014 #0.095 AR1SOMG_02_01 + #0.014 &
          #0.095 AR1SOMG_03_01 + #0.014 #0.095 AR1SOMG_04_01 + #0.014 #0.095 &
          AR1SOMG_05_01 + #0.014 #0.000 AR1SOMG_06_01 + #0.014 #0.095 AR1SOMG_01_02 + &
          #0.014 #0.095 AR1SOMG_02_02 + #0.014 #0.095 AR1SOMG_03_02 + #0.014 &
          #0.095 AR1SOMG_04_02 + #0.014 #0.095 AR1SOMG_05_02 + #0.014 #0.000 &
          AR1SOMG_06_02 + #0.014 #0.095 AR1SOMG_02_03 + #0.014 #0.095 AR1SOMG_03_03 + &
          #0.014 #0.095 AR1SOMG_04_03 + #0.014 #0.095 AR1SOMG_05_03 + #0.014 &
          #0.095 AR1SOMG_06_03 + #0.014 #0.095 AR1SOMG_02_04 + #0.014 #0.095 &
          AR1SOMG_03_04 + #0.014 #0.095 AR1SOMG_04_04 + #0.014 #0.095 AR1SOMG_05_04 + &
          #0.014 #0.095 AR1SOMG_06_04
 S35.1  AR1SOMG_07_04 + OH = #0.939 #0.123 AR1SOMG_07_05 + #0.939 #0.001 &
          AR1SOMG_07_06 + #0.939 #0.876 AR1SOMG_07_07 + #0.061 #0.083 AR1SOMG_01_01 + &
          #0.061 #0.083 AR1SOMG_02_01 + #0.061 #0.083 AR1SOMG_03_01 + #0.061 &
          #0.083 AR1SOMG_04_01 + #0.061 #0.000 AR1SOMG_05_01 + #0.061 #0.000 &
          AR1SOMG_06_01 + #0.061 #0.083 AR1SOMG_01_02 + #0.061 #0.083 AR1SOMG_02_02 + &
          #0.061 #0.083 AR1SOMG_03_02 + #0.061 #0.083 AR1SOMG_04_02 + #0.061 &
          #0.083 AR1SOMG_05_02 + #0.061 #0.000 AR1SOMG_06_02 + #0.061 #0.083 &
          AR1SOMG_02_03 + #0.061 #0.083 AR1SOMG_03_03 + #0.061 #0.083 AR1SOMG_04_03 + &
          #0.061 #0.083 AR1SOMG_05_03 + #0.061 #0.000 AR1SOMG_06_03 + #0.061 &
          #0.083 AR1SOMG_02_04 + #0.061 #0.083 AR1SOMG_03_04 + #0.061 #0.083 &
          AR1SOMG_04_04 + #0.061 #0.083 AR1SOMG_05_04 + #0.061 #0.083 AR1SOMG_06_04 + &
          #0.061 #0.083 AR1SOMG_03_05 + #0.061 #0.083 AR1SOMG_04_05 + #0.061 &
          #0.083 AR1SOMG_05_05 + #0.061 #0.083 AR1SOMG_06_05
 S36.1  AR1SOMG_07_05 + OH = #0.814 #0.123 AR1SOMG_07_06 + #0.814 #0.877 &
          AR1SOMG_07_07 + #0.186 #0.077 AR1SOMG_01_01 + #0.186 #0.077 AR1SOMG_02_01 + &
          #0.186 #0.077 AR1SOMG_03_01 + #0.186 #0.077 AR1SOMG_04_01 + #0.186 &
          #0.000 AR1SOMG_05_01 + #0.186 #0.000 AR1SOMG_06_01 + #0.186 #0.077 &
          AR1SOMG_01_02 + #0.186 #0.077 AR1SOMG_02_02 + #0.186 #0.077 AR1SOMG_03_02 + &
          #0.186 #0.077 AR1SOMG_04_02 + #0.186 #0.000 AR1SOMG_05_02 + #0.186 &
          #0.000 AR1SOMG_06_02 + #0.186 #0.077 AR1SOMG_02_03 + #0.186 #0.077 &
          AR1SOMG_03_03 + #0.186 #0.077 AR1SOMG_04_03 + #0.186 #0.077 AR1SOMG_05_03 + &
          #0.186 #0.000 AR1SOMG_06_03 + #0.186 #0.077 AR1SOMG_02_04 + #0.186 &
          #0.077 AR1SOMG_03_04 + #0.186 #0.077 AR1SOMG_04_04 + #0.186 #0.077 &
          AR1SOMG_05_04 + #0.186 #0.000 AR1SOMG_06_04 + #0.186 #0.077 AR1SOMG_03_05 + &
          #0.186 #0.077 AR1SOMG_04_05 + #0.186 #0.077 AR1SOMG_05_05 + #0.186 &
          #0.077 AR1SOMG_06_05 + #0.186 #0.077 AR1SOMG_03_06 + #0.186 #0.077 &
          AR1SOMG_04_06 + #0.186 #0.077 AR1SOMG_05_06 + #0.186 #0.077 AR1SOMG_06_06
 S37.1  AR1SOMG_07_06 + OH = #0.537 #1.000 AR1SOMG_07_07 + #0.463 #0.077 &
          AR1SOMG_01_01 + #0.463 #0.077 AR1SOMG_02_01 + #0.463 #0.077 AR1SOMG_03_01 + &
          #0.463 #0.000 AR1SOMG_04_01 + #0.463 #0.000 AR1SOMG_05_01 + #0.463 &
          #0.000 AR1SOMG_06_01 + #0.463 #0.077 AR1SOMG_01_02 + #0.463 #0.077 &
          AR1SOMG_02_02 + #0.463 #0.077 AR1SOMG_03_02 + #0.463 #0.077 AR1SOMG_04_02 + &
          #0.463 #0.000 AR1SOMG_05_02 + #0.463 #0.000 AR1SOMG_06_02 + #0.463 &
          #0.077 AR1SOMG_02_03 + #0.463 #0.077 AR1SOMG_03_03 + #0.463 #0.077 &
          AR1SOMG_04_03 + #0.463 #0.000 AR1SOMG_05_03 + #0.463 #0.000 AR1SOMG_06_03 + &
          #0.463 #0.077 AR1SOMG_02_04 + #0.463 #0.077 AR1SOMG_03_04 + #0.463 &
          #0.077 AR1SOMG_04_04 + #0.463 #0.077 AR1SOMG_05_04 + #0.463 #0.000 &
          AR1SOMG_06_04 + #0.463 #0.077 AR1SOMG_03_05 + #0.463 #0.077 AR1SOMG_04_05 + &
          #0.463 #0.077 AR1SOMG_05_05 + #0.463 #0.000 AR1SOMG_06_05 + #0.463 &
          #0.077 AR1SOMG_03_06 + #0.463 #0.077 AR1SOMG_04_06 + #0.463 #0.077 &
          AR1SOMG_05_06 + #0.463 #0.077 AR1SOMG_06_06 + #0.463 #0.077 AR1SOMG_04_07 + &
          #0.463 #0.077 AR1SOMG_05_07 + #0.463 #0.077 AR1SOMG_06_07
 S38.1  AR1SOMG_07_07 + OH = #1.000 #0.087 AR1SOMG_01_01 + #1.000 #0.087 &
          AR1SOMG_02_01 + #1.000 #0.087 AR1SOMG_03_01 + #1.000 #0.000 AR1SOMG_04_01 + &
          #1.000 #0.000 AR1SOMG_05_01 + #1.000 #0.000 AR1SOMG_06_01 + #1.000 &
          #0.087 AR1SOMG_01_02 + #1.000 #0.087 AR1SOMG_02_02 + #1.000 #0.087 &
          AR1SOMG_03_02 + #1.000 #0.000 AR1SOMG_04_02 + #1.000 #0.000 AR1SOMG_05_02 + &
          #1.000 #0.000 AR1SOMG_06_02 + #1.000 #0.087 AR1SOMG_02_03 + #1.000 &
          #0.087 AR1SOMG_03_03 + #1.000 #0.087 AR1SOMG_04_03 + #1.000 #0.000 &
          AR1SOMG_05_03 + #1.000 #0.000 AR1SOMG_06_03 + #1.000 #0.087 AR1SOMG_02_04 + &
          #1.000 #0.087 AR1SOMG_03_04 + #1.000 #0.087 AR1SOMG_04_04 + #1.000 &
          #0.000 AR1SOMG_05_04 + #1.000 #0.000 AR1SOMG_06_04 + #1.000 #0.087 &
          AR1SOMG_03_05 + #1.000 #0.087 AR1SOMG_04_05 + #1.000 #0.087 AR1SOMG_05_05 + &
          #1.000 #0.000 AR1SOMG_06_05 + #1.000 #0.087 AR1SOMG_03_06 + #1.000 &
          #0.087 AR1SOMG_04_06 + #1.000 #0.087 AR1SOMG_05_06 + #1.000 #0.000 &
          AR1SOMG_06_06 + #1.000 #0.087 AR1SOMG_04_07 + #1.000 #0.087 AR1SOMG_05_07 + &
          #1.000 #0.087 AR1SOMG_06_07


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
     29) [4] ARO2                30) [4] ARO2RO2             31) [4] XYLNRXN             32) [4] XYLHRXN         
     33) [4] ARO1                34) [4] AR1SOMG             35) [4] BENZENE             36) [4] BENZRO2         
     37) [4] BNZNRXN             38) [4] BNZHRXN             39) [4] TERP                40) [4] TRPRXN          
     41) [4] ISOPRENE            42) [4] ISOPRXN             43) [4] SESQ                44) [4] SESQRXN         
     45) [5] O3P                 46) [5] O1D                 47) [5] TBUO                48) [5] BZO             
     49) [5] NO2EX               50) [5] xHO2                51) [5] xOH                 52) [5] xNO2            
     53) [5] xMEO2               54) [5] xMECO3              55) [5] xRCO3               56) [5] xMACO3          
     57) [5] xTBUO               58) [5] xCO                 59) [5] xHCHO               60) [5] xCCHO           
     61) [5] xRCHO               62) [5] xACET               63) [5] xMEK                64) [5] xPROD2          
     65) [5] xGLY                66) [5] xMGLY               67) [5] xBACL               68) [5] xBALD           
     69) [5] xAFG1               70) [5] xAFG2               71) [5] xMACR               72) [5] xMVK            
     73) [5] xIPRD               74) [5] xRNO3               75) [5] zRNO3               76) [5] yROOH           
     77) [5] yR6OOH              78) [5] yRAOOH              79) [5] xCL                 80) [5] xCLCCHO         
     81) [5] xCLACET             82) [6] HO2                 83) [6] CO                  84) [6] CO2             
     85) [6] HO2H                86) [6] MEO2                87) [6] HCHO                88) [6] COOH            
     89) [6] MEOH                90) [6] RO2C                91) [6] RO2XC               92) [6] MECO3           
     93) [6] CCOOH               94) [6] RCO3                95) [6] N3                  96) [6] RCOOH           
     97) [6] BZCO3               98) [6] XC                  99) [6] MACO3              100) [6] ACET            
    101) [6] CRES               102) [6] CCHO               103) [6] RCHO               104) [6] MEK             
    105) [6] HCOOH              106) [6] ROOH               107) [6] R6OOH              108) [6] PROD2           
    109) [6] RAOOH              110) [6] MGLY               111) [6] IPRD               112) [6] GLY             
    113) [6] AFG1               114) [6] AFG2               115) [6] HCOCO3             116) [6] BACL            
    117) [6] CATL               118) [6] AFG3               119) [6] xCNDPP             120) [6] yRAOOHP         
    121) [6] XYNL               122) [6] CNDPP              123) [6] BALD               124) [6] PHEN            
    125) [6] MACR               126) [6] MVK                127) [6] AFG4               128) [6] ETOH            
    129) [6] xAFG4              130) [6] CH4                131) [6] ETHENE             132) [6] ACETYLEN        
    133) [6] ALK1               134) [6] ALK2               135) [6] ALK3               136) [6] ALK4            
    137) [6] ALKSOMG            138) [6] OLE1               139) [6] OLE2               140) [6] CL2             
    141) [6] CL                 142) [4] HCL                143) [6] CLO                144) [6] HOCL            
    145) [6] CLCCHO             146) [6] CLACET             147) [6] CLCHO              148) [4] NH3             
    149) [4] GDENS              150) [4] MMOLES             151) [4] XICBC              152) [4] YICBC           
    153) [4] ZICBC              154) [6] CHCL3              155) [4] AR1SOMG_01_00      156) [4] AR1SOMG_01_01   
    157) [4] AR1SOMG_01_02      158) [4] AR1SOMG_02_00      159) [4] AR1SOMG_02_01      160) [4] AR1SOMG_02_02   
    161) [4] AR1SOMG_02_03      162) [4] AR1SOMG_02_04      163) [4] AR1SOMG_03_00      164) [4] AR1SOMG_03_01   
    165) [4] AR1SOMG_03_02      166) [4] AR1SOMG_03_03      167) [4] AR1SOMG_03_04      168) [4] AR1SOMG_03_05   
    169) [4] AR1SOMG_03_06      170) [4] AR1SOMG_04_00      171) [4] AR1SOMG_04_01      172) [4] AR1SOMG_04_02   
    173) [4] AR1SOMG_04_03      174) [4] AR1SOMG_04_04      175) [4] AR1SOMG_04_05      176) [4] AR1SOMG_04_06   
    177) [4] AR1SOMG_04_07      178) [4] AR1SOMG_05_00      179) [4] AR1SOMG_05_01      180) [4] AR1SOMG_05_02   
    181) [4] AR1SOMG_05_03      182) [4] AR1SOMG_05_04      183) [4] AR1SOMG_05_05      184) [4] AR1SOMG_05_06   
    185) [4] AR1SOMG_05_07      186) [4] AR1SOMG_06_00      187) [4] AR1SOMG_06_01      188) [4] AR1SOMG_06_02   
    189) [4] AR1SOMG_06_03      190) [4] AR1SOMG_06_04      191) [4] AR1SOMG_06_05      192) [4] AR1SOMG_06_06   
    193) [4] AR1SOMG_06_07      194) [4] AR1SOMG_07_00      195) [4] AR1SOMG_07_01      196) [4] AR1SOMG_07_02   
    197) [4] AR1SOMG_07_03      198) [4] AR1SOMG_07_04      199) [4] AR1SOMG_07_05      200) [4] AR1SOMG_07_06   
    201) [4] AR1SOMG_07_07   
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
           8  ITERATIONS REQUIRED TO ORDER STEADY STATE SPECIES.
0CODE PREPARED.  LENGTH= 1631 STATEMENTS
1          26-Sep-18 13:23 TEMP 300.                                                       
 
 
 LIST OF REACTANTS:  ( 142 ACTIVE,   6 CONSTANT,   8 BUILDUP,  36 STEADY STATE,   0 DUMMY, AND   9 DELETED
 
 TYPE           NO. NAME                 MWT    C'S   N'S  S'S     O'S   X'S   CONC(0)
 ------------   --- ----------------  -------  -----  ---  ---  -----  -----  -----------
 
 ACTIVE          15  O3                  48.00   0.00   0    0   3.00   0.00     0.00    
                 16  NO                  30.01   0.00   1    0   0.00   0.00     0.00    
                 17  NO2                 46.01   0.00   1    0   0.00   0.00     0.00    
                 18  NO3                 62.01   0.00   1    0   0.00   0.00     0.00    
                 19  N2O5               108.02   0.00   2    0   0.00   0.00     0.00    
                 20  HNO3                63.02   0.00   1    0   0.00   0.00     0.00    
                 21  HONO                47.02   0.00   1    0   0.00   0.00     0.00    
                 22  HNO4                79.02   0.00   1    0   0.00   0.00     0.00    
                 23  SO2                 64.06   0.00   0    1   0.00   0.00     0.00    
                 24  SULF                96.06   0.00   0    1   0.00   0.00     0.00    
                 25  NPHE               139.11   6.00   1    0   0.00   0.00     0.00    
                 26  XN                   0.00   0.00   1    0   0.00   0.00     0.00    
                 27  RNO3                 0.00   0.00   1    0   0.00   0.00     0.00    
                 28  PAN                121.05   2.00   1    0   0.00   0.00     0.00    
                 29  PAN2               135.08   3.00   1    0   0.00   0.00     0.00    
                 30  PBZN               183.13   7.00   1    0   0.00   0.00     0.00    
                 31  MAPAN              132.00   3.00   1    0   0.00   0.00     0.00    
                 32  CLNO                65.46   0.00   1    0   0.00   0.00     0.00    
                 33  CLONO               81.46   0.00   1    0   0.00   0.00     0.00    
                 34  CLNO2               81.46   0.00   1    0   0.00   0.00     0.00    
                 35  CLONO2              97.46   0.00   1    0   0.00   0.00     0.00    
                 36  ALK5               170.00   0.00   0    0   0.00   1.00     0.00    
                 37  ARO2                 0.00   0.00   0    0   0.00   1.00     0.00    
                 38  ARO2RO2              0.00   0.00   0    0   0.00   1.00     0.00    
                 39  XYLNRXN            192.00   0.00   0    0   0.00   1.00     0.00    
                 40  XYLHRXN            192.00   0.00   0    0   0.00   1.00     0.00    
                 41  ARO1                92.14   7.00   0    0   0.00   1.00     0.00    
                 42  BENZENE              0.00   0.00   0    0   0.00   1.00     0.00    
                 43  BENZRO2              0.00   0.00   0    0   0.00   1.00     0.00    
                 44  BNZNRXN            144.00   0.00   0    0   0.00   1.00     0.00    
                 45  BNZHRXN            144.00   0.00   0    0   0.00   1.00     0.00    
                 46  TERP                 0.00   0.00   0    0   0.00   1.00     0.00    
                 47  TRPRXN             168.00   0.00   0    0   0.00   1.00     0.00    
                 48  ISOPRENE             0.00   0.00   0    0   0.00   1.00     0.00    
                 49  ISOPRXN             96.00   0.00   0    0   0.00   1.00     0.00    
                 50  SESQ                 0.00   0.00   0    0   0.00   1.00     0.00    
                 51  SESQRXN            378.00   0.00   0    0   0.00   1.00     0.00    
                 52  HO2                  0.00   0.00   0    0   0.00   0.00     0.00    
                 53  CO                   0.00   0.00   0    0   0.00   0.00     0.00    
                 54  HO2H                 0.00   0.00   0    0   0.00   0.00     0.00    
                 55  MEO2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 56  HCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 57  COOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 58  MEOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 59  RO2C                 0.00   0.00   0    0   0.00   0.00     0.00    
                 60  RO2XC                0.00   0.00   0    0   0.00   0.00     0.00    
                 61  MECO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 62  CCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 63  RCO3                 0.00   0.00   0    0   0.00   0.00     0.00    
                 64  RCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 65  BZCO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 66  MACO3                0.00   0.00   0    0   0.00   0.00     0.00    
                 67  ACET                 0.00   0.00   0    0   0.00   0.00     0.00    
                 68  CRES                 0.00   0.00   0    0   0.00   0.00     0.00    
                 69  CCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 70  RCHO                 0.00   0.00   0    0   0.00   0.00     0.00    
                 71  MEK                  0.00   0.00   0    0   0.00   0.00     0.00    
                 72  HCOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 73  ROOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 74  R6OOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 75  PROD2                0.00   0.00   0    0   0.00   0.00     0.00    
                 76  RAOOH                0.00   0.00   0    0   0.00   0.00     0.00    
                 77  MGLY                 0.00   0.00   0    0   0.00   0.00     0.00    
                 78  IPRD                 0.00   0.00   0    0   0.00   0.00     0.00    
                 79  GLY                  0.00   0.00   0    0   0.00   0.00     0.00    
                 80  AFG1                 0.00   0.00   0    0   0.00   0.00     0.00    
                 81  AFG2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 82  HCOCO3               0.00   0.00   0    0   0.00   0.00     0.00    
                 83  BACL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 84  CATL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 85  AFG3                 0.00   0.00   0    0   0.00   0.00     0.00    
                 86  XYNL                 0.00   0.00   0    0   0.00   0.00     0.00    
                 87  BALD                 0.00   0.00   0    0   0.00   0.00     0.00    
                 88  PHEN                 0.00   0.00   0    0   0.00   0.00     0.00    
                 89  MACR                 0.00   0.00   0    0   0.00   0.00     0.00    
                 90  MVK                  0.00   0.00   0    0   0.00   0.00     0.00    
                 91  AFG4                 0.00   0.00   0    0   0.00   0.00     0.00    
                 92  ETOH                 0.00   0.00   0    0   0.00   0.00     0.00    
                 93  xAFG4                0.00   0.00   0    0   0.00   0.00     0.00    
                 94  CH4                  0.00   0.00   0    0   0.00   0.00     0.00    
                 95  ETHENE               0.00   0.00   0    0   0.00   0.00     0.00    
                 96  ACETYLEN             0.00   0.00   0    0   0.00   0.00     0.00    
                 97  ALK1                 0.00   0.00   0    0   0.00   0.00     0.00    
                 98  ALK2                 0.00   0.00   0    0   0.00   0.00     0.00    
                 99  ALK3                 0.00   0.00   0    0   0.00   0.00     0.00    
                100  ALK4                 0.00   0.00   0    0   0.00   0.00     0.00    
                101  OLE1                 0.00   0.00   0    0   0.00   0.00     0.00    
                102  OLE2                 0.00   0.00   0    0   0.00   0.00     0.00    
                103  CL2                  0.00   0.00   0    0   0.00   0.00     0.00    
                104  CL                   0.00   0.00   0    0   0.00   0.00     0.00    
                105  HCL                 36.46   0.00   0    0   0.00   0.00     0.00    
                106  CLO                  0.00   0.00   0    0   0.00   0.00     0.00    
                107  HOCL                 0.00   0.00   0    0   0.00   0.00     0.00    
                108  CLCCHO               0.00   0.00   0    0   0.00   0.00     0.00    
                109  CLACET               0.00   0.00   0    0   0.00   0.00     0.00    
                110  NH3                 17.03   0.00   1    0   0.00   0.00     0.00    
                111  GDENS                0.00   0.00   0    0   0.00   0.00     0.00    
                112  MMOLES               0.00   0.00   0    0   0.00   0.00     0.00    
                113  XICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                114  YICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                115  ZICBC                0.00   0.00   0    0   0.00   0.00     0.00    
                116  CHCL3                0.00   0.00   0    0   0.00   0.00     0.00    
                117  AR1SOMG_01_01       31.00   1.00   0    0   1.00   1.00     0.00    
                118  AR1SOMG_01_02       46.00   1.00   0    0   2.00   1.00     0.00    
                119  AR1SOMG_02_01       45.00   2.00   0    0   1.00   1.00     0.00    
                120  AR1SOMG_02_02       60.00   2.00   0    0   2.00   1.00     0.00    
                121  AR1SOMG_02_03       75.00   2.00   0    0   3.00   1.00     0.00    
                122  AR1SOMG_02_04       90.00   2.00   0    0   4.00   1.00     0.00    
                123  AR1SOMG_03_01       59.00   3.00   0    0   1.00   1.00     0.00    
                124  AR1SOMG_03_02       74.00   3.00   0    0   2.00   1.00     0.00    
                125  AR1SOMG_03_03       89.00   3.00   0    0   3.00   1.00     0.00    
                126  AR1SOMG_03_04      104.00   3.00   0    0   4.00   1.00     0.00    
                127  AR1SOMG_03_05      119.00   3.00   0    0   5.00   1.00     0.00    
                128  AR1SOMG_03_06      134.00   3.00   0    0   6.00   1.00     0.00    
                129  AR1SOMG_04_01       73.00   4.00   0    0   1.00   1.00     0.00    
                130  AR1SOMG_04_02       88.00   4.00   0    0   2.00   1.00     0.00    
                131  AR1SOMG_04_03      103.00   4.00   0    0   3.00   1.00     0.00    
                132  AR1SOMG_04_04      118.00   4.00   0    0   4.00   1.00     0.00    
                133  AR1SOMG_04_05      133.00   4.00   0    0   5.00   1.00     0.00    
                134  AR1SOMG_04_06      148.00   4.00   0    0   6.00   1.00     0.00    
                135  AR1SOMG_04_07      163.00   4.00   0    0   7.00   1.00     0.00    
                136  AR1SOMG_05_01       87.00   5.00   0    0   1.00   1.00     0.00    
                137  AR1SOMG_05_02      102.00   5.00   0    0   2.00   1.00     0.00    
                138  AR1SOMG_05_03      117.00   5.00   0    0   3.00   1.00     0.00    
                139  AR1SOMG_05_04      132.00   5.00   0    0   4.00   1.00     0.00    
                140  AR1SOMG_05_05      147.00   5.00   0    0   5.00   1.00     0.00    
                141  AR1SOMG_05_06      162.00   5.00   0    0   6.00   1.00     0.00    
                142  AR1SOMG_05_07      177.00   5.00   0    0   7.00   1.00     0.00    
                143  AR1SOMG_06_01      101.00   6.00   0    0   1.00   1.00     0.00    
                144  AR1SOMG_06_02      116.00   6.00   0    0   2.00   1.00     0.00    
                145  AR1SOMG_06_03      131.00   6.00   0    0   3.00   1.00     0.00    
                146  AR1SOMG_06_04      146.00   6.00   0    0   4.00   1.00     0.00    
                147  AR1SOMG_06_05      161.00   6.00   0    0   5.00   1.00     0.00    
                148  AR1SOMG_06_06      176.00   6.00   0    0   6.00   1.00     0.00    
                149  AR1SOMG_06_07      191.00   6.00   0    0   7.00   1.00     0.00    
                150  AR1SOMG_07_01      115.00   7.00   0    0   1.00   1.00     0.00    
                151  AR1SOMG_07_02      130.00   7.00   0    0   2.00   1.00     0.00    
                152  AR1SOMG_07_03      145.00   7.00   0    0   3.00   1.00     0.00    
                153  AR1SOMG_07_04      160.00   7.00   0    0   4.00   1.00     0.00    
                154  AR1SOMG_07_05      175.00   7.00   0    0   5.00   1.00     0.00    
                155  AR1SOMG_07_06      190.00   7.00   0    0   6.00   1.00     0.00    
                156  AR1SOMG_07_07      205.00   7.00   0    0   7.00   1.00     0.00    
 
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
                 13  ALKSOMG              0.00   0.00   0    0   0.00   0.00     0.00    
                 14  CLCHO                0.00   0.00   0    0   0.00   0.00     0.00    
 
 STEADY STATE   157  O3P             
                158  O1D             
                159  TBUO            
                160  BZO             
                161  xHO2            
                162  xOH             
                163  xNO2            
                164  xMEO2           
                165  xMECO3          
                166  xRCO3           
                167  xMACO3          
                168  xTBUO           
                169  xCO             
                170  xHCHO           
                171  xCCHO           
                172  xRCHO           
                173  xACET           
                174  xMEK            
                175  xPROD2          
                176  xGLY            
                177  xMGLY           
                178  xBACL           
                179  xBALD           
                180  xAFG1           
                181  xAFG2           
                182  xMACR           
                183  xMVK            
                184  xIPRD           
                185  xRNO3           
                186  zRNO3           
                187  yROOH           
                188  yR6OOH          
                189  yRAOOH          
                190  xCL             
                191  xCLCCHO         
                192  xCLACET         
 
 DELETED        193  AR1SOMG         
                194  NO2EX           
                195  AR1SOMG_01_00   
                196  AR1SOMG_02_00   
                197  AR1SOMG_03_00   
                198  AR1SOMG_04_00   
                199  AR1SOMG_05_00   
                200  AR1SOMG_06_00   
                201  AR1SOMG_07_00   

   20 VARIABLE COEFFICIENTS:
       1) LITTLE   1.000E-30    2) INERT     0.00        3) NROG      0.00        4) METHANE   0.00        5) RO2RO     0.00    
       6) RO2XRO    0.00        7) QY.ACET  5.000E-02    8) QY.MEK   0.175        9) QY.NO2A  1.500E-03   10) QY.NPHE  1.500E-02
      11) QY.BALD  6.000E-02   12) QY.PROD2 4.860E-03   13) RO2NO     0.00       14) RO22NN    0.00       15) RO2HO2    0.00    
      16) RO2RO2M   0.00       17) QY.CLAC  0.500       18) RO2NO3    0.00       19) RO2RO2    0.00       20) RO2RO3    0.00    

  526 CONSTANT COEFFICIENTS:
    4001)   2.00     4002)  0.610     4003)  0.390     4004)  0.800     4005)  0.200    
    4006)  0.500     4007)  0.750     4008)  0.250     4009)  0.600     4010)  0.400    
    4011)  0.440     4012)  0.410     4013)  0.150     4014)  0.100     4015)  0.900    
    4016)   2.24     4017)  0.560     4018)  -1.00     4019)   3.00     4020)   4.00    
    4021)  0.965     4022)  3.500E-02 4023)  0.620     4024)   1.38     4025)  0.380    
    4026)  0.967     4027)  3.900E-02 4028)  0.376     4029)  0.510     4030)  7.400E-02
    4031)  8.800E-02 4032)  0.504     4033)  0.300     4034)  0.509     4035)  0.491    
    4036) -0.491     4037)  0.143     4038)  0.142     4039)  0.457     4040) -0.455    
    4041)  0.700     4042)  0.744     4043)  0.251     4044)  4.000E-03 4045)  0.239    
    4046)  1.200E-02 4047)  0.205     4048)  3.400E-02 4049)  0.256     4050) -0.115    
    4051)  0.840     4052)  0.222     4053)  2.900E-02 4054)  9.000E-02 4055)  4.100E-02
    4056)  2.000E-02 4057)  7.500E-02 4058)  8.400E-02 4059)  0.160     4060)  0.782    
    4061)  7.700E-02 4062)  8.500E-02 4063)  2.600E-02 4064)  5.800E-02 4065)  0.698    
    4066)  0.858     4067)  1.700E-02 4068)  0.139     4069)  0.148     4070)  0.589    
    4071)  0.124     4072)  0.147     4073)  0.565     4074)  2.400E-02 4075)  0.448    
    4076)  3.000E-02 4077)  0.252     4078)  7.300E-02 4079)  0.713     4080)   2.67    
    4081)   1.40     4082)  0.170     4083)  0.130     4084)  4.000E-02 4085)  0.115    
    4086)  0.715     4087)  7.000E-02 4088)  0.230     4089)  0.110     4090)  0.210    
    4091)  0.995     4092)  0.930     4093)   1.65     4094)  0.330     4095) -8.000E-02
    4096)   1.90     4097)   6.00     4098)   7.00     4099)  9.500E-02 4100)  0.105    
    4101)  4.800E-02 4102)  6.500E-02 4103) -0.205     4104)  0.217     4105)  0.723    
    4106)  6.000E-02 4107)  0.521     4108)  0.201     4109)  0.334     4110)  0.407    
    4111)  0.129     4112)  0.107     4113)  0.267     4114)  0.783     4115)  0.284    
    4116)   1.02     4117)  0.173     4118)  0.305     4119)  0.695     4120)  0.195    
    4121)  0.206     4122)  0.733     4123)  0.117     4124)  0.561     4125)  0.114    
    4126)  0.274     4127)  0.153     4128)  1.900E-02 4129)  0.231     4130)  0.794    
    4131)  0.938     4132)  0.471     4133)  0.554     4134)  1.300E-02 4135)  0.258    
    4136)  7.000E-03 4137)  0.580     4138)  0.190     4139)  0.366     4140)  0.184    
    4141)  0.350     4142)  3.000E-03 4143)  0.163     4144)  0.264     4145) -0.575    
    4146)  0.902     4147)  9.800E-02 4148)  0.416     4149) -0.416     4150)  0.208    
    4151)  0.108     4152)  0.450     4153)  0.333     4154) -0.100     4155)   1.50    
    4156)  0.670     4157)  0.340     4158)  0.975     4159)  2.500E-02 4160)  0.675    
    4161) -0.725     4162)  0.164     4163)  6.400E-02 4164)  5.000E-02 4165)  0.475    
    4166)  0.950     4167)  0.351     4168) -5.000E-02 4169)  0.550     4170)  -2.20    
    4171)  0.289     4172)  0.336     4173)  5.500E-02 4174)  0.332     4175)  0.174    
    4176) -0.504     4177)  0.711     4178)  0.285     4179)  0.498     4180)  0.140    
    4181)  2.300E-02 4182)  0.742     4183)  0.372     4184)  4.700E-02 4185)  1.000E-03
    4186) -0.329     4187)  0.799     4188)  5.100E-02 4189)  0.572     4190)  0.227    
    4191)  0.218     4192)  8.000E-03 4193)  0.850     4194)  0.278     4195)  0.185    
    4196)   1.23     4197)  0.467     4198)  0.233     4199) -0.233     4200)  0.472    
    4201)  0.379     4202)  4.900E-02 4203)  0.473     4204)  7.100E-02 4205)  2.000E-03
    4206)  0.211     4207)  8.300E-02 4208)  0.402     4209)  0.329     4210)  0.528    
    4211)  0.877     4212)  0.913     4213)   1.59     4214)  8.700E-02 4215)  0.303    
    4216)  0.780     4217) -9.100E-02 4218)  0.189     4219)  0.313     4220)  0.976    
    4221)  0.175     4222)  1.100E-02 4223)  0.429     4224)  3.600E-02 4225)  1.000E-02
    4226)  3.100E-02 4227)  0.157     4228)  0.636     4229)  0.344     4230)  0.721    
    4231)  0.102     4232)  6.100E-02 4233)  0.214     4234)  6.300E-02 4235)  0.261    
    4236)  6.600E-02 4237)  0.591     4238)  0.396     4239)  8.100E-02 4240) -1.000E-03
    4241)   5.00     4242)  -3.00     4243)  -4.00     4244)  -6.00     4245)  -8.00    
    4246)   1.61     4247)  0.120     4248)  0.370     4249)  0.290     4250)  0.986    
    4251)  9.300E-02 4252)  0.907     4253)  0.624     4254)  0.320     4255)  0.357    
    4256)  0.833     4257)  0.266     4258)  0.192     4259)  0.275     4260)  0.122    
    4261)  0.204     4262)  0.401     4263)  0.936     4264)  0.749     4265)  0.187    
    4266)  0.813     4267)  0.240     4268)  0.116     4269)  0.570     4270)  0.314    
    4271)  0.764     4272)  0.704     4273) -0.105     4274)  0.236     4275)   1.25    
    4276)  0.445     4277)  0.983     4278) -4.600E-02 4279)  0.830     4280)   1.76    
    4281)  0.149     4282)  0.438     4283)  0.426     4284)  0.106     4285)  0.146    
    4286) -0.119     4287)  0.647     4288)   1.61     4289)  0.353     4290)  0.209    
    4291)  8.600E-02 4292)   2.00     4293)  0.904     4294)   1.14     4295)  0.301    
    4296)  0.470     4297)  5.000E-03 4298)  6.000E-03 4299)  0.119     4300)  0.413    
    4301)  0.587     4302)  0.822     4303)  0.193     4304)  0.104     4305)  0.368    
    4306)  0.125     4307)  2.200E-02 4308)  0.112     4309)  3.700E-02 4310)  0.690    
    4311)  0.824     4312)   1.31     4313)  0.176     4314)  9.000E-03 4315)  0.546    
    4316)  0.454     4317)  0.437     4318)  0.113     4319)   1.22     4320)  0.914    
    4321)  0.966     4322)  0.788     4323)  0.481     4324)  0.136     4325)  7.600E-02
    4326)  2.700E-02 4327)  0.643     4328)  0.111     4329)  0.423     4330)  0.297    
    4331)  0.162     4332)  0.152     4333)  0.428     4334)  6.700E-02 4335)  0.315    
    4336)  1.800E-02 4337)  3.300E-02 4338)  2.800E-02 4339)  2.100E-02 4340)  4.200E-02
    4341)  6.900E-02 4342)  0.128     4343)  0.409     4344)   1.18     4345)  0.154    
    4346)  0.322     4347)  0.269     4348)  1.400E-02 4349)  0.709     4350)  0.202    
    4351)  0.666     4352)  0.166     4353)  0.482     4354)  6.800E-02 4355)  0.138    
    4356)  0.403     4357)  0.286     4358)  8.900E-02 4359)  0.601     4360)   1.49    
    4361)  0.804     4362)  0.759     4363)   1.15     4364)  0.533     4365)  0.255    
    4366)   5.06     4367)  5.200E-02 4368)  0.585     4369)  0.126     4370)  0.875    
    4371)  0.203     4372)  4.500E-02 4373)  7.900E-02 4374)  0.220     4375)  0.165    
    4376)  4.300E-02 4377)  0.545     4378)   3.53     4379)  0.421     4380)   1.51    
    4381)  0.397     4382)   4.47     4383)  0.853     4384)   4.44     4385)   1.42    
    4386)  0.763     4387)  3.800E-02 4388)   1.28     4389)  0.159     4390)  0.547    
    4391)  4.600E-02 4392)  0.908     4393)  0.851     4394)  0.680     4395)  0.198    
    4396)  0.541     4397)  0.237     4398)  0.109     4399)  0.686     4400)   1.26    
    4401)  0.630     4402)   1.26     4403) -0.370     4404)  0.414     4405)  0.588    
    4406)  0.197     4407)  0.586     4408) -0.287     4409)  0.145     4410)   1.08    
    4411)  0.502     4412)  0.186     4413)  0.676     4414)  0.280     4415)  0.855    
    4416)  0.348     4417)  0.404     4418)  0.802     4419)  8.200E-02 4420)  0.180    
    4421)  0.835     4422)   1.28     4423)  5.300E-02 4424)  0.625     4425)  0.947    
    4426)  0.538     4427)  0.730     4428)  0.712     4429)  0.420     4430)  0.762    
    4431)   1.17     4432)  0.738     4433)  0.177     4434)  0.671     4435)  0.970    
    4436)  0.488     4437) -9.000E-02 4438)  9.400E-02 4439)   1.36     4440)  7.800E-02
    4441)  0.343     4442)  0.253     4443)  0.827     4444)   1.74     4445)  0.287    
    4446)  0.412     4447)  0.247     4448)  0.327     4449)   1.54     4450)  0.352    
    4451)  8.000E-02 4452)  4.400E-02 4453)  0.378     4454)   2.37     4455)  0.308    
    4456)  0.188     4457)  0.225     4458)  0.447     4459)   1.51     4460)  0.263    
    4461)  0.228     4462)  0.361     4463)  0.179     4464)  0.880     4465)  0.323    
    4466)  0.842     4467)  0.158     4468)  0.618     4469)  0.224     4470)   2.38    
    4471)  1.600E-02 4472)   2.26     4473)  0.582     4474)  0.548     4475)   3.54    
    4476)   1.00     4477)  0.123     4478)  0.874     4479)  0.969     4480)  0.876    
    4481)  3.125E-02 4482)   0.00     4483)  0.996     4484)  4.115E-03 4485)  0.868    
    4486)  0.132     4487)  0.999     4488)  9.766E-04 4489)  0.333     4490)  0.286    
    4491)  0.763     4492)  0.237     4493)  0.222     4494)   1.00     4495)  3.200E-04
    4496)  0.990     4497)  1.024E-02 4498)  0.922     4499)  7.776E-02 4500)  0.154    
    4501)  0.672     4502)  0.328     4503)  0.143     4504)  0.167     4505)   1.00    
    4506)  1.286E-04 4507)  0.118     4508)  0.105     4509)  0.598     4510)  0.402    
    4511)   1.00     4512)  5.950E-05 4513)  0.998     4514)  1.904E-03 4515)  0.986    
    4516)  1.446E-02 4517)  9.524E-02 4518)  0.939     4519)  6.093E-02 4520)  8.333E-02
    4521)  0.814     4522)  0.186     4523)  7.692E-02 4524)  0.537     4525)  0.463    
    4526)  8.696E-02

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
 NO. REACTIONS =  402
 NO. SPECIES =  192
 NO. VARIABLE COEFFICIENTS =   20
 NO. CONSTANT COEFFICIENTS =  526
 KINETIC PARM BUFFER LENGTH = 1230
 PRODUCT ARRAY LENGTH = 4883
 
 MODEL = 'saprc14_rev1    ' STORED.

 **********************************************************************************************************************************

                                              E N D   P R E P   O U T P U T

 **********************************************************************************************************************************
