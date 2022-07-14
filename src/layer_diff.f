
C     ****************************************
C     * diffusion 
C     **************************************** 


C     WRITTEN BY Sam O'Donnell, Oct. 2020

C-----INPUTS--------------------------------------

C-----OUTPUTS-------------------------------------


      SUBROUTINE layer_diff(Kzl,numlay,temp_lay, pres_lay,
     & Nk_layo, Mk_layo, aadt, icomp, ibins)


      IMPLICIT NONE 

C-----INCLUDE FILES-------------------------------

C      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS-----------------------
     
      INTEGER v,n,I
      INTEGER :: numlay, icomp, ibins 
      DOUBLE PRECISION Kzl 
      DOUBLE PRECISION temp_lay
      DOUBLE PRECISION pres_lay
      DOUBLE PRECISION Nk_layer(ibins,numlay)
      DOUBLE PRECISION Mk_layer(ibins,icomp,numlay)
      DOUBLE PRECISION Nk_layo(ibins,numlay)
      DOUBLE PRECISION Mk_layo(ibins,icomp,numlay)
      DOUBLE PRECISION dz 
      DOUBLE PRECISION Jdiff
      DOUBLE PRECISION dCdt
      DOUBLE PRECISION aadt 
      DOUBLE PRECISION delmass
      DOUBLE PRECISION ratio

C      Nk_layo = Nk_layer
C      Mk_layo = Mk_layer


      DO v=1,numlay-1
        dz = (287.0*temp_lay(v)/9.8)*log(pres_lay(v)/pres_lay(v+1))
C        print*,'dz=',dz
        DO n=1,ibins
          Jdiff = Kzl*(Nk_layo(n,v)-Nk_layo(n,v+1))/dz
C          print*,'Kzl =',Kzl
C          print*,'DiffConc =',Nk_layo(n,v)-Nk_layo(n,v+1)
C          print*,'Jdiff =',Jdiff
          dCdt = Jdiff/dz
C          print*,'Bin number =',n 
C          print*,'Number of particles moving =',dCdt*aadt
          Nk_layo(n,v) = Nk_layo(n,v)-dCdt*aadt
          Nk_layo(n,v+1) = Nk_layo(n,v+1)+dCdt*aadt
          delmass = 0.d0
          IF (Jdiff.gt.0.d0) THEN
            ratio = ABS(dCdt*aadt/Nk_layo(n,v))
            DO I=1,icomp
              delmass = delmass+Mk_layo(n,I,v)*ratio
              Mk_layo(n,I,v) = Mk_layo(n,I,v)-Mk_layo(n,I,v)*ratio
              Mk_layo(n,I,v+1) = Mk_layo(n,I,v+1)+
     &          Mk_layo(n,I,v)*ratio
            ENDDO
          ELSEIF(Jdiff.lt.0.d0) THEN
            ratio = ABS(dCdt*aadt/Nk_layo(n,v+1))
            DO I=1,icomp
              delmass = delmass+Mk_layo(n,I,v+1)*ratio
              Mk_layo(n,I,v+1)=Mk_layo(n,I,v+1)-
     &          Mk_layo(n,I,v+1)*ratio
              Mk_layo(n,I,v) = Mk_layo(n,I,v)+
     &          Mk_layo(n,I,v+1)*ratio
            ENDDO
          ENDIF
C          print*,'Mass moving =',delmass
C          print*,'M/N moving =',delmass/(dCdt*aadt)
C          print*,'ratio =',ratio
C          ELSEIF(Jdiff.eq.0.d0) THEN
C            ratio = 0.d0
C          ENDIF
C          Nk_layo(n,v) = Nk_layo(n,v)-dCdt*aadt
C          Nk_layo(n,v+1) = Nk_layo(n,v+1)+dCdt*aadt
        ENDDO
      ENDDO

C      Mk_layer = Mk_layo
C      Nk_layer = Nk_layo
      RETURN 
      END
