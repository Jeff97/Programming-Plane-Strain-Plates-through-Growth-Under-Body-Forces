      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
      ! Surface name for a surface-based load definition (JLTYP=0). For a body force or an element-based surface load the surface name is passed in as blank.

      ! element C3D20 consists 27 Gauss points in each element
! TODO: one should modify the number of elements according to the mesh
      Parameter (NElem=3072, NGauss=27)
      REAL*8 X, Y, Z, RhoR, fZ, ReadX, ReadY, ReadZ, 
     &       TotalT, TargetF0, TargetF, Pi, ReadDetF, DetF

      COMMON /DLoadUSER/ ReadX(NElem, NGauss)
      COMMON /DLoadUSER/ ReadY(NElem, NGauss)
      COMMON /DLoadUSER/ ReadZ(NElem, NGauss)
      COMMON /DLoadUSER/ ReadDetF(NElem, NGauss)
     
        RhoR = 1000.0 ! kg/m^3
        fZ = -10 ! N/kg, downward

      IF((KSTEP == 1) .AND. (KINC == 1)) THEN
      ! KSPT---Step number, KINC---Increment number
      ! store the orginal coordinates

        ReadX(NOEL,NPT)=COORDS(1)
        ReadY(NOEL,NPT)=COORDS(2)
        ReadZ(NOEL,NPT)=COORDS(3)

        IF((NOEL == 25) .AND. (NPT == 1)) THEN
          WRITE(*,*) "NOEL", NOEL
          WRITE(*,*) "NPT", NPT
          WRITE(*,*) "TIME", TIME
          WRITE(*,*) "LAYER", LAYER
          WRITE(*,*) "KSTEP", KSTEP
          WRITE(*,*) "JLTYP", JLTYP
        END IF

        ! WRITE(*,*) "NOEL", NOEL
        ! WRITE(*,*) "NPT", NPT
        ! WRITE(*,*) "JLTYP", JLTYP
      END IF

      ! the required F is verified to be a scalar
      ! ABAQUS apply the body force in the deformed volume
      X = ReadX(NOEL,NPT)
      Y = ReadY(NOEL,NPT)
      Z = ReadZ(NOEL,NPT)

      TotalT = 1.0

      ! assign the det of F
      DetF = ReadDetF(NOEL, NPT)
      IF((KSTEP == 1) .AND. (KINC == 1)) THEN
        DetF = 1.0
      END IF

!     Body force in X-axis
      IF(JLTYP == 1) THEN
        TargetF = 0
      END IF

!     Body force in Y-axis
      IF(JLTYP == 2) THEN
        TargetF0 = RhoR*fZ
        TargetF = TargetF0/DetF ! Correct 
      END IF

!     Body force in Z-axis
      IF(JLTYP == 3) THEN
        TargetF = 0
      END IF

! TODO: for testing 
      ! TargetF = 0 
      F = TargetF*TIME(1)/TotalT

      RETURN
      END
! ------------------------- UMAT subroutine ---------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C LOCAL ARRAYS
C ----------------------------------------------------------------
C Ae - Elastic tensor at the end of the increment
C Some geometry quantities:
C       Er, Fr, Gr, Lr, Mr, Nr are related to reference surface
C       EE, FF, GG, LL, MM, NN are related to target surface
C ----------------------------------------------------------------
C
        Parameter (NElem=3072, NGauss=27)
        DOUBLE PRECISION :: X, Y, Z, RhoR, fZ, C0, h, Epsilon, 
     &  ReadX, ReadY, ReadZ, ReadDetF, DetF, Pi, TotalT, 
     &  Lambda1z0,Lambda1z1, Lambda2z0,Lambda2z1, LambdaSz0,LambdaSz1,
     &  Lambda1z2, EMOD, ENU, C10, D1, c1, c2, c3

      DOUBLE PRECISION :: MyStress(9), Ae(3,3), G(3,3)
      DOUBLE PRECISION :: F(3,3), b(3,3), Elasticity(9,9), bBar(3,3)
      DOUBLE PRECISION :: Ib, IbBar, DETAe

      INTEGER ::  ii, jj, kk, ll
      INTEGER ::  Seq(6), SeqEls(6)

      COMMON /DLoadUSER/ ReadX(NElem, NGauss)
      COMMON /DLoadUSER/ ReadY(NElem, NGauss)
      COMMON /DLoadUSER/ ReadZ(NElem, NGauss)
      COMMON /DLoadUSER/ ReadDetF(NElem, NGauss)
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C STATEV(1) - Initial coordinate COORDS(1)
C STATEV(2) - Initial coordinate COORDS(2)
C STATEV(3) - Initial coordinate COORDS(3)
C STATEV(4) - 
C STATEV(5) - 
C STATEV(9) - Norm of growth tensor
C ----------------------------------------------------------------
C ELASTIC PROPERTIES
        RhoR = 1000.0 ! kg/m^3
        fZ = -10.0 ! N/kg, downward
        C0 =  PROPS(1) ! Pa
        h = 0.01 ! m (total thickness = 2h)
        Epsilon = (RhoR*fZ)*1.0/(C0) ! 0.0001
        ! if Epsilon <0, beam bend downward- from 1 ->0.987
        ! if Epsilon >0, beam bend upward+  from 1 ->1.013
        Pi = 3.14159265359

        ENU = 0.4999
        C10 = C0
        D1 = (1.5*(1. - 2.*ENU))/(C10*(2. + 2.*ENU))

        Seq=(/ 1,5,9,2,3,6 /)
        SeqEls=(/ 1,5,9,4,7,8 /)
        ! C10 = EMOD/(4.0*(1.0+ENU))
        ! D1 = 6.0*(1.0-2.0*ENU)/EMOD
C get coordinate from COORDS variable
        IF (STATEV(1) .EQ. 0) THEN
                STATEV(1)=COORDS(1)
        END IF
        IF (STATEV(2) .EQ. 0) THEN
                STATEV(2)=COORDS(2)
        END IF
        IF (STATEV(3) .EQ. 0) THEN
                STATEV(3)=COORDS(3)
        END IF

C calculate the determinant of the deformation gradient in 3D

        F = DFGRD1
        call DetMatrix(F, DetF)

        ReadDetF(NOEL, NPT) = DetF
C store coordinate from STATEV
        ! X = ReadX(NOEL,NPT)
        ! Y = ReadY(NOEL,NPT)
        ! Z = ReadZ(NOEL,NPT)

        X = STATEV(1)
        Y = STATEV(2)
        Z = STATEV(3)
C assign growth functions Lambda

        Lambda1z0 = Sqrt(1 + X**2)
        Lambda2z0 = 1.0
        LambdaSz0 = 0.0

        Lambda1z1 = -(1/(1 + X**2))

        Lambda2z1 = 0.0
        LambdaSz1 = 0.0
!I cannot get total time from UMAT variable, please modify TotalT manually 
        TotalT=1.0

        DtltaG11 = (Lambda1z0 + Y*Lambda1z1-1.0)
     &             *(TIME(1)+DTIME)/TotalT
        DtltaG12 = 0.0
        DtltaG21 = DtltaG12
        DtltaG22 = (Lambda2z0 + Y*Lambda2z1-1.0)*(TIME(1)+DTIME)/TotalT

C G is growth tensor
        G11 = 1.0+DtltaG11
        G12 = DtltaG12
        G21 = DtltaG21
        G22 = 1.0+DtltaG22
        G33 = 1.0

        ! STATEV(1) = X
        ! STATEV(2) = Y
        ! STATEV(3) = Z
        STATEV(4) = DetF

        STATEV(5) = Lambda1z0
        STATEV(6) = Lambda1z1
        STATEV(7) = Lambda1z2
        STATEV(8) = G11
        STATEV(9) = Sqrt( (G11**2.0+G12**2.0+G21**2.0
     &                   +G22**2.0+G33**2.0)/3.0 )

C Elastic deformation tensor A=F.G**(-1)
C
        G = 0.0
        G(1,1) = G11
        G(2,2) = G22
        G(3,3) = G33
        G(1,2) = G12
        G(2,1) = G21

        call CalAe(Ae, F, G)
C
C Determinent of the elastic deformation tensor DETA=det(A)

        call DetMatrix(Ae, DETAe)

! left Cauchy-Green tensor b=F F^T
      b = 0.0
      DO ii=1, 3
        DO jj=1, 3
          DO KK=1, 3
            b(ii,jj) = b(ii,jj) + Ae(ii,KK)*Ae(jj,KK)
          END DO
        END DO
      END DO

! the first variant of b
      Ib = b(1,1) + b(2,2)+ b(3, 3)

! left Cauchy-Green tensor bBar = F F^T
      DO ii=1, 3
        DO jj=1, 3
          bBar(ii,jj) = b(ii,jj)*DETAe**(-2.0/3.0)
        END DO
      END DO

! the first variant of bBar
      IbBar = bBar(1,1) + bBar(2,2)+ bBar(3, 3)

! Calculate my Cauchy stress {11,12,13,21,22,23,31,32,33} ------------
      MyStress(1) = 2.0*C10*DETAe**(-5.0/3.0) 
     &  *(b(1,1) - Ib/3.0) + 2.0/D1*(DETAe - 1.0)
      MyStress(2) = 2.0*C10*DETAe**(-5.0/3.0) * b(1,2)
      MyStress(3) = 2.0*C10*DETAe**(-5.0/3.0) * b(1,3)

      MyStress(4) = 2.0*C10*DETAe**(-5.0/3.0) * b(2,1)
      MyStress(5) = 2.0*C10*DETAe**(-5.0/3.0) 
     &  *(b(2,2) - Ib/3.0) + 2.0/D1*(DETAe - 1.0)
      MyStress(6) = 2.0*C10*DETAe**(-5.0/3.0) * b(2,3)

      MyStress(7) = 2.0*C10*DETAe**(-5.0/3.0) * b(3,1)
      MyStress(8) = 2.0*C10*DETAe**(-5.0/3.0) * b(3,2)
      MyStress(9) = 2.0*C10*DETAe**(-5.0/3.0) 
     &  *(b(3,3) - Ib/3.0) + 2.0/D1*(DETAe - 1.0)
! Convert my stress to ABAQUS stress {11,22,33,12,13,23}
      DO ii=1,6
        STRESS(ii) = MyStress(Seq(ii))
      END DO

! Jacobian matrix of the constitutive model --------------------------
! Elasticity row = Sigma {11,21,31,12,22,32,13,23,33}
! Elasticity column = epsilon {11,21,31,12,22,32,13,23,33}
! here is the elasticity tensor that obtained from ABAQUS manual
! In ABAQUS UMAT, the elasticity tensor is 
! cc_ijkl = 1/J [(partial(J sigma_ij)/partial F_kM)F_lM 
!                 + (partial(J sigma_ij)/partial F_lM)F_kM]
      Elasticity = 0.0 

      ! delta_ik B_jl
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (ii==kk) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(0.5*bBar(jj,ll))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_jl B_ik
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (jj==ll) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(0.5*bBar(ii,kk))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_il B_jk
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (ii==ll) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(0.5*bBar(jj,kk))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_jk B_il
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (jj==kk) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(0.5*bBar(ii,ll))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_ij B_kl
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (ii==jj) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(-2.0/3.0*bBar(kk,ll))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_kl B_ij
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (kk==ll) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(-2.0/3.0*bBar(ii,jj))
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_ij delta_kl IB
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (ii==jj .AND. kk==ll) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/DETAe*C10*(2.0/9.0*IbBar)
             END IF
            END DO
          END DO
        END DO
      END DO

      ! delta_ij delta_kl 
      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            DO ll=1, 3
             IF (ii==jj .AND. kk==ll) THEN
              Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) = 
     1           Elasticity(3*(jj-1)+ii,3*(ll-1)+kk) 
     2           + 2.0/D1*(2.0*DETAe-1.0)
             END IF
            END DO
          END DO
        END DO
      END DO

! Convert my Elasticity to ABAQUS DDSDDE
! DDSDDE row = Sigma {11,22,33,12,13,23}
! DDSDDE column = epsilon {11,22,33,12,13,23}
      DO ii=1,6
        DO jj=1,6
          DDSDDE(ii,jj) = Elasticity(SeqEls(ii),SeqEls(jj))
        END DO
      END DO

! Abaqus use engineering shear strain,
! so DDSDDE(i,4), DDSDDE(i,5), DDSDDE(i,6) need to multiply 0.5
! is this still works in large deformation? I commented it
!       DO ii=4,6
!         DDSDDE(ii,4) = DDSDDE(ii,4)*0.5
!         DDSDDE(ii,5) = DDSDDE(ii,5)*0.5
!         DDSDDE(ii,6) = DDSDDE(ii,6)*0.5
!       END DO

      RETURN

      END SUBROUTINE UMAT
! ---------------------------- DetMatrix ------------------------------
      SUBROUTINE DetMatrix(M, DetM)

      DOUBLE PRECISION :: M(3,3), DetM

      DetM = -M(1,3)*M(2,2)*M(3,1) + M(1,2)*M(2,3)*M(3,1)
     &       +M(1,3)*M(2,1)*M(3,2) - M(1,1)*M(2,3)*M(3,2)
     &       -M(1,2)*M(2,1)*M(3,3) + M(1,1)*M(2,2)*M(3,3)

      END SUBROUTINE DetMatrix

! ---------------------------- InvMatrix ------------------------------
      SUBROUTINE InvMatrix(MInv, M, DetM)

      DOUBLE PRECISION :: MInv(3,3), M(3,3), DetM

      MInv(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/DetM
      MInv(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/DetM
      MInv(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/DetM
      MInv(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/DetM
      MInv(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/DetM
      MInv(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/DetM
      MInv(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/DetM
      MInv(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/DetM
      MInv(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/DetM

      END SUBROUTINE InvMatrix

! ------------------------------ CalAe --------------------------------
      SUBROUTINE CalAe(Ae, F, G)

      DOUBLE PRECISION :: Ae(3,3), F(3,3), G(3,3)

        Ae(1,1) = (F(1,3)*G(2,2)*G(3,1)-F(1,2)*G(2,3)*G(3,1)-F(1,3)*
     &    G(2,1)*G(3,2)+F(1,1)*G(2,3)*G(3,2)+F(1,2)*G(2,1)*G(3,3)-
     &    F(1,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(1,2) = (F(1,3)*G(1,2)*G(3,1)-F(1,2)*G(1,3)*G(3,1)-F(1,3)*
     &    G(1,1)*G(3,2)+F(1,1)*G(1,3)*G(3,2)+F(1,2)*G(1,1)*G(3,3)-
     &    F(1,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(1,3) = (F(1,3)*G(1,2)*G(2,1)-F(1,2)*G(1,3)*G(2,1)-F(1,3)*
     &    G(1,1)*G(2,2)+F(1,1)*G(1,3)*G(2,2)+F(1,2)*G(1,1)*G(2,3)-
     &    F(1,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

        Ae(2,1) = (F(2,3)*G(2,2)*G(3,1)-F(2,2)*G(2,3)*G(3,1)-F(2,3)*
     &    G(2,1)*G(3,2)+F(2,1)*G(2,3)*G(3,2)+F(2,2)*G(2,1)*G(3,3)-
     &    F(2,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(2,2) = (F(2,3)*G(1,2)*G(3,1)-F(2,2)*G(1,3)*G(3,1)-F(2,3)*
     &    G(1,1)*G(3,2)+F(2,1)*G(1,3)*G(3,2)+F(2,2)*G(1,1)*G(3,3)-
     &    F(2,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(2,3) = (F(2,3)*G(1,2)*G(2,1)-F(2,2)*G(1,3)*G(2,1)-F(2,3)*
     &    G(1,1)*G(2,2)+F(2,1)*G(1,3)*G(2,2)+F(2,2)*G(1,1)*G(2,3)-
     &    F(2,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

        Ae(3,1) = (F(3,3)*G(2,2)*G(3,1)-F(3,2)*G(2,3)*G(3,1)-F(3,3)*
     &    G(2,1)*G(3,2)+F(3,1)*G(2,3)*G(3,2)+F(3,2)*G(2,1)*G(3,3)-
     &    F(3,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(3,2) = (F(3,3)*G(1,2)*G(3,1)-F(3,2)*G(1,3)*G(3,1)-F(3,3)*
     &    G(1,1)*G(3,2)+F(3,1)*G(1,3)*G(3,2)+F(3,2)*G(1,1)*G(3,3)-
     &    F(3,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(3,3) = (F(3,3)*G(1,2)*G(2,1)-F(3,2)*G(1,3)*G(2,1)-F(3,3)*
     &    G(1,1)*G(2,2)+F(3,1)*G(1,3)*G(2,2)+F(3,2)*G(1,1)*G(2,3)-
     &    F(3,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

      END SUBROUTINE CalAe
! ---------------------------------------------------------------------
!       SUBROUTINE UTRACLOAD(ALPHA,T_USER,KSTEP,KINC,TIME,NOEL,NPT,
!      1 COORDS,DIRCOS,JLTYP,SNAME)
! C
!       INCLUDE 'ABA_PARAM.INC'
! C
!       DIMENSION T_USER(3), TIME(2), COORDS(3), DIRCOS(3,3)
!       CHARACTER*80 SNAME

!       Parameter (NNode=11603, NElem=1600, NGauss=27)
!       REAL*8 X, Y, Z, ParaC, ParaD, ReadX, ReadY, ReadZ, TotalT,
!      &        TargetF, Pi, Ce

!       COMMON /TractionUSER/ ReadX(NElem, NGauss)
!       COMMON /TractionUSER/ ReadY(NElem, NGauss)
!       COMMON /TractionUSER/ ReadZ(NElem, NGauss)

!       ! I use second order element C3D20 in the input file
!       ! there are 27 NPT in each element

!       ParaC = 161877.143
!       ParaD = 1.323755E-6
!       Pi = 3.14159265359
!       Ce = 2.7182818284

!       IF((KSTEP == 1) .AND. (KINC == 1)) THEN
!       ! during the first increment of first step
!       ! store the orginal coordinates

!         ReadX(NOEL,NPT)=COORDS(1)
!         ReadY(NOEL,NPT)=COORDS(2)
!         ReadZ(NOEL,NPT)=COORDS(3)

!       END IF

!       X = ReadX(NOEL,NPT)
!       Y = ReadY(NOEL,NPT)
!       Z = ReadZ(NOEL,NPT)

!       TotalT = 1.0

! !     Surface-based load, pressure
!       IF(SNAME == "ASSEMBLY_SURF-X0") THEN
!         TargetF = -((1133140.0*(-1.0+2.0*Y)*(3.0*(1.0+2.0*Y)
!      &      **(5.0/3.0)-(3.0*(3.0+2.0*Y))/(7.0*2.0**(1.0/3.0))))
!      &      /(9.0*(1.0+2.0*Y)**(2.0/3.0)))
!         ! WRITE(*,*) "SNAME", SNAME
!         T_USER(1) = 1.0
!         T_USER(2) = 0.0
!       END IF

!       IF(SNAME == "ASSEMBLY_SURF-X1") THEN
!         TargetF = ((1133140.0*(-1.0+2.0*Y)*(3.0*(1.0+2.0*Y)
!      &      **(5.0/3.0)-(3.0*(3.0+2.0*Y))/(7.0*2.0**(1.0/3.0))))
!      &      /(9.0*(1.0+2.0*Y)**(2.0/3.0)))
!         ! WRITE(*,*) "SNAME", SNAME
!         T_USER(1) = 1.0
!         T_USER(2) = 0.0
!       END IF
      
!       IF(SNAME == "ASSEMBLY_SURF-Y1") THEN
!         TargetF = 0
!         ! WRITE(*,*) "SNAME", SNAME
!         T_USER(1) = 0
!         T_USER(2) = 1
!       END IF

!       TotalT = 1.0
!       T_USER(3) = 0.0

! ! for testing, if body force is correct
!       ! TargetF = 0 
!       ALPHA = TargetF*TIME(1)/TotalT 

!       RETURN
!       END