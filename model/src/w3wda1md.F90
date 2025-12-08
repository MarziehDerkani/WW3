!> @file
!> @brief Contains module W3WDA1MD.
!>
!> @author M. Derkani (University of Western Australia) 
!> @date 03-Jan-2025
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Externally supplied data assimilation code.
!>
!> @details This module is intended to serve externally supplied 
!> data assimilation software to WAVEWATCH III interface W3WDASMD.
!>
!> @author M. Derkani (University of Western Australia) 
!> @date 02-Jan-2025
!>
MODULE W3WDA1MD
  PUBLIC :: W3WDA1
  PRIVATE :: SPCPAR, CALC_DA_DIANA, CALC_WEIGHT, &
             DA_HS_UPD2, INIT_GET_DATISEA, UPSPEC, &
             WSDUR, WSANA, ISCLOSE
  PRIVATE
  !/ Constants for duration limited energy growth (nondimensional)
  !/      E* = E0 tanh( A (t*)**B )
    REAL, PARAMETER :: E0 = 955.0
    REAL, PARAMETER :: A = 6.02E-5
    REAL, PARAMETER :: B = 0.695
  !/ Constants for nondimensional frequency-energy curve
  !/      E* = AF (f*)**BF
  !/REAL, PARAMETER :: AF = 1.68E-4  ! Linello et al. (1992)
  !/REAL, PARAMETER :: BF = -3.27
    REAL, PARAMETER :: AF = 5.054E-4 ! Toledano et al. (2022)
    REAL, PARAMETER :: BF = -2.959
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Data assimilation for mean wave parameters.
  !>
  !> @param[in] NDAT Number of records (rows) in data set.
  !> @param[in] MDAT Number of parameters in each data sets.
  !> @param[in] DATA0 Observations (mean parameters).
  !>
  !> @author M. Derkani (University of Western Australia)
  !> @date 02-Jan-2025
  !>
  SUBROUTINE W3WDA1 ( MDAT, NDAT, DATA0 )
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           H. L. Tolman            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         03-Jan-2025 |
    !/                  +-----------------------------------+
    !/
    !/    03-Jan-2025 : Origination.                        ( version 7.14 )
    !/
    !  1. Purpose :
    !
    !     WAVEWATCH III data assimilation routine for mean wave parameters.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       NDAT      I.A.   I   Records (rows) in data set.
    !       MDAT      I.A.   I   Number of parameters for each set.
    !       DATA0     R.A.   I   Observations (mean parameters).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing.
    !      EXTCDE    Subr. W3SERVMD Program abort.
    !      INIT_GET_ISEA Subr. W3PARALL Local JSEA to global ISEA
    !      INIT_GET_JSEA_ISPROC Subr. W3PARALL Global ISEA to JSEA/ISPROC
    !      W3DIST    Subr. W3GSRUMD Distance calculation
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3WDAS    Subr. W3WDASMD Wave data assimilation interface.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS, ONLY: RADIUS, TPI, DERA, RADE, UNDEF
    USE W3ADATMD, ONLY: CG, WN, U10, U10D, DW
#ifdef W3_MPI
    USE W3ADATMD, ONLY: MPI_COMM_WAVE
#endif
    USE W3WDATMD, ONLY: VA, ASF, UST
    USE W3GDATMD, ONLY: NSPEC, NK, NTH, NSEAL, SIG, MAPSF, &
                        MAPSTA, FLAGLL, XGRD, YGRD, ICLOSE, &
                        DA1METHOD
    USE W3PARALL, ONLY: INIT_GET_ISEA, INIT_GET_JSEA_ISPROC
    USE W3GSRUMD, ONLY: W3DIST
    USE W3ODATMD, ONLY: NDSO, NDSE, NDST, SCREEN, NAPROC, IAPROC, &
                        NAPLOG, NAPOUT, NAPERR, DIMP, &
                        WSCUT, PTMETH, FLCOMB
    USE W3PARTMD, ONLY: W3PART
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !
#ifdef W3_MPI
    USE MPI_F08
#endif
    !
    IMPLICIT NONE
    !
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    !/
    INTEGER, INTENT(IN)     :: MDAT, NDAT
    REAL, INTENT(IN)        :: DATA0(MDAT,NDAT)
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters :
    !/
    INTEGER                 :: IDAT, IX, IY, ISEA, JSEA, IK, ITH, &
                               MREC, DATPROC, DIMXP, NP, IP, W3PTMETH
    LOGICAL                 :: W3FLCOMB
    REAL, PARAMETER         :: HSMIN = 0.01
    REAL, PARAMETER         :: SFCUT = 0.50
    REAL                    :: HSDAT, TMDAT, S1DAT, S2DAT
    REAL                    :: HS, TM, S1, S2, A(NSPEC), E(NK,NTH), CG1(NK)
    REAL                    :: TMAN, HSAN, W, W2, WS, TSEA, USTAN
    REAL                    :: FACT, HSSEA, FMSEA, FMSEAAN, SEAFR, SEAAN
    REAL                    :: UABS, UDIR, DEPTH
    REAL, ALLOCATABLE       :: WP(:,:)
    REAL(KIND=8)            :: XDAT, YDAT, DKM, DIST2KM
#ifdef W3_MPI
    INTEGER                 :: IERR_MPI
#endif
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
    !/
    IF (MDAT.LT.4) THEN
      WRITE(NDSE,1000) MDAT
      CALL EXTCDE(99)
    END IF
    !
    !/ ------------------------------------------------------------------- /
    ! 1.  Initializations and test output
    ! 1.a Subroutine tracing
    !
#ifdef W3_S
    CALL STRACE (IENT, 'W3WDA1')
#endif
    !
    ! 2.  Data routine --------------------------------------------------- /
    !
    !    Store user-defined partition method (to be reset at end)
    W3PTMETH = PTMETH
    W3FLCOMB = FLCOMB
    IF (PTMETH.NE.1) PTMETH = 1              ! Standard partition method
    IF (FLCOMB.EQV..FALSE.) FLCOMB = .TRUE.  ! w/ combine wind sea
    DIMXP = ((NK+1)/2) * ((NTH-1)/2)
    ALLOCATE( WP(DIMP, 0:DIMXP) )
    !
    IF (FLAGLL) THEN
      DIST2KM = DBLE(RADIUS*DERA*1.000E-3)
    ELSE
      DIST2KM = DBLE(1.000E-3)
    ENDIF

    DO IDAT=1, NDAT
      TMDAT = UNDEF
      HSDAT = UNDEF
      XDAT = DBLE(DATA0(1,IDAT))
      YDAT = DBLE(DATA0(2,IDAT))
      IF (DATA0(3,IDAT) .LT. HSMIN) CYCLE

      CALL INIT_GET_DATISEA(XDAT, YDAT, ISEA)
      IF (ISEA.LE.0) CYCLE

      CALL INIT_GET_JSEA_ISPROC(ISEA, JSEA, DATPROC)
     !WRITE(*,'(2X,A8,3I6)') "MPI",IAPROC, NAPROC, DATPROC

      IF (IAPROC .EQ. DATPROC) THEN
        A(1:NSPEC) = VA(1:NSPEC,JSEA)
        CG1(1:NK) = CG(1:NK, ISEA)
        CALL SPCPAR(A, CG1, HSDAT, TMDAT, S1DAT, S2DAT)
      ! WRITE(*,'(2X,A8,4I6,4F8.3)') "DA-PROC", IAPROC, DATPROC, ISEA, JSEA, HSDAT, TMDAT, S1DAT, S2DAT
      ENDIF
#ifdef W3_MPI
      CALL MPI_BARRIER(MPI_COMM_WAVE, IERR_MPI)
      CALL MPI_BCAST(HSDAT, 1, MPI_REAL, DATPROC-1, MPI_COMM_WAVE, IERR_MPI)
      CALL MPI_BCAST(TMDAT, 1, MPI_REAL, DATPROC-1, MPI_COMM_WAVE, IERR_MPI)
      CALL MPI_BARRIER(MPI_COMM_WAVE, IERR_MPI)
#endif
     !WRITE(*,'(2X,A8,4I6,2F7.2)') "IAPROC ", IAPROC, DATPROC, ISEA, JSEA, HSDAT, TMDAT

      DO JSEA=1, NSEAL
        CALL INIT_GET_ISEA(ISEA, JSEA)
        IX = MAPSF(ISEA,1)
        IY = MAPSF(ISEA,2)
        IF (MAPSTA(IY,IX) .LT. 0) CYCLE

        DKM = W3DIST(FLAGLL,XDAT,YDAT,XGRD(IY,IX),YGRD(IY,IX))*DIST2KM
        IF (DKM .LT. 1000.0) THEN
          CG1(1:NK) = CG(1:NK, ISEA)
          A(1:NSPEC) = VA(1:NSPEC, JSEA)
          CALL SPCPAR(A, CG1, HS, TM, S1, S2)
         !WRITE(*,'(2X,A8,4I5,2F7.2,F10.1)') "    MOD", IAPROC, NAPROC, ISEA, JSEA, HS, TM, DKM
         !WRITE(*,'(2X,A,2F6.2,2F9.6)') "TEST [FC:HS,T01,S,SM]",HS, TM, S1, S2
    !
    !     a) Calculate correlation length scale
          SELECT CASE(DA1METHOD)
            CASE DEFAULT
              WRITE (NDSE,1010) DA1METHOD
              CALL EXTCDE(99)
            CASE(0) ! VOORRIPS ET AL (1997)
              CALL CALC_DA_VOORRIPS(DKM, WS)
            CASE(1) ! GREENSLADE & YOUNG (2004)
              CALL CALC_DA_DIANA(DKM, YGRD(IY,IX), YDAT, WS)
          END SELECT
    !
    !     b) Calculate weights for analysis
          CALL CALC_WEIGHT(WS, HS, DATA0(3,IDAT), HSDAT, &
                               TM, DATA0(4,IDAT), TMDAT, &
                            W, W2, TMAN, HSAN)
    !
    ! 3.  Actual data assimilation  -------------------------------------- /
    !     a) Convert spectrum A(k,th) to E(f,th) for partitioning
          DO IK=1, NK
            FACT = TPI * SIG(IK) / CG1(IK)
            DO ITH=1, NTH
              E(IK,ITH) = A(ITH+(IK-1)*NTH) * FACT
            END DO
          END DO
    !
          FMSEA = UNDEF   ! Mean frequency of wind sea (model guess)
          FMSEAAN = UNDEF ! Mean frequency of wind sea (analysis)
          SEAFR = UNDEF   ! Fraction of wind sea to total sig. wave height
          USTAN = UNDEF   ! Friction velocity (analysis)
          TSEA = UNDEF    ! Duration of wind sea (from growth curve)
    !
    !     b) Partition wave spectra and find wind sea partition
          UABS = U10(ISEA)*ASF(ISEA)
          UDIR = U10D(ISEA)*RADE
          DEPTH = DW(ISEA)
          IF (DEPTH.NE.DEPTH) THEN
            WRITE (NDSE,1020) ISEA,IX,IY,DW(ISEA)
            CALL EXTCDE(99)
          END IF
    !
          CALL W3PART(E, UABS, UDIR, DEPTH, WN(1:NK,ISEA), NP, WP, DIMXP)
          ! Array WP contains integral parameters describing partitions,
          ! where index IP=0 contains parameters for entire spectrum.
          DO IP=1, NP
            ! Scan for wind sea part (wind sea fraction >= threshold).
            IF (WP(6,IP).GE.WSCUT ) THEN
              HSSEA = WP(1, IP)
              FMSEA = 1.0/WP(13, IP) ! mean frequency
              SEAFR = MIN((HSSEA*HSSEA)/(WP(1,0)*WP(1,0)), 1.0)
            ! WRITE(*,"(A,1X,2I2,2F6.2,2F7.4)") "W3PART", PTMETH, IP,   &
            !      HSSEA, FMSEA, WP(6,IP), SEAFR
            END IF
          END DO
    !
    !     c) Estimate duration of wind sea from wind sea analysis
    !        (based on wind sea fraction calculated above)
          IF (SEAFR.GT.SFCUT .AND. HSSEA.GT.HSMIN) THEN
            ! Wind sea analysis based on wind sea fraction (from b)
            SEAAN = 4.0 * SQRT(SEAFR*(HSAN/4)*(HSAN/4))
            CALL WSDUR(UST(ISEA), HSSEA, TSEA)
    !
    !     d) Estimate analysis friction velocity and
    !        analysis mean frequency of wind sea
            CALL WSANA(UST(ISEA), TSEA, SEAAN, USTAN, FMSEAAN)
          END IF
    !
    !     e) Update spectrum by stretching and scaling
          CALL UPSPEC(A, CG1, HS, HSAN, FMSEA, FMSEAAN)
         !CALL SPCPAR(A, CG1, HS, TM, S1, S2)
         !WRITE(*,'(2X,A,2F6.2,2F9.6)') "TEST [AN:HS,T01,S,SM]", &
         !    HS, TM, S1, S2
    !
    ! 4.  Copy assimilated spectrum back to data structure --------------- /
          VA(1:NSPEC, JSEA) = A(1:NSPEC)
    !
         !IF (IAPROC .EQ. DATPROC) THEN
         !  CALL SPCPAR(VA(1:NSPEC, JSEA), CG1, HS, TM, S1, S2)
         !  WRITE(*,'(2X,A8,4I6,2F7.2,E10.3)') "DATPROC", IAPROC, DATPROC, ISEA, JSEA, HS, TM, WSTP
         !END IF !/ IAPROC == DATPROC
        END IF !/ (DKM .LT. 1000.0)
      END DO !/ JSEA..NSEAL
    END DO  !/ IDAT, NDAT
    !
    ! Restore user-defined partition method
    PTMETH = W3PTMETH
    FLCOMB = W3FLCOMB
    !
    RETURN
    !
    ! Formats
    !
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3WDA1 :'/             &
         '     DATA RECORD DIMENSION <4 : ',I8)
1010 FORMAT (/' *** WAVEWATCH III ERROR IN W3WDA1 :'/             &
         '     SCHEME W/ METHOD', I2,' NOT IMPLEMENTED.')
1020 FORMAT (/' *** WAVEWATCH III ERROR IN W3WDA1 :'/             &
         '     NaN found in depth at ISEA,IX,IY', 3I8)


  END SUBROUTINE W3WDA1
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Integrate mean paramters from wave spectra.
  !>
  !> Wave steepness S (deep water approximation,
  !> see Mendes and Oliveira (2021)):
  !>   S = Hs/L = Hs * 2pi/(g*T02**2)
  !>
  !> Spectral steepness SM
  !>   SM = M0 * OMEGA**4 / GRAV**2 Eq (4.69) in Young (1999)
  !> or expressed with wave number M0 * WNMEAN**2 / (4pi**2)
  !>
  !> Mendes and Oliveira (2021): Deep-water spectral
  !>     wave steepness offshore mainland Portugal
  !>     https://doi.org/10.1016/j.oceaneng.2021.109548
  !> Lionello et al. (1992): Assimilation of altimeter data
  !>     in a global third-generation wave model, JGR Ocean
  !> Young (1999): Wind generated ocean waves, Cambridge Press
  !>
  !> @param[in] A    Action density spectrum
  !> @param[in] CG   Group velocities
  !>
  !> @param[out] HSIG
  !> @param[out] T01
  !> @param[out] S
  !> @param[out] SM
  !>
  !> @author  @date
  !>
  SUBROUTINE SPCPAR ( A, CG, HSIG, T01, S, SM )
    !/
    USE CONSTANTS, ONLY: TPI, TPIINV, GRAV
    USE W3GDATMD, ONLY: NSPEC, NK, NTH, DDEN, SIG, FTE, DTH
    !
    IMPLICIT NONE
    REAL, INTENT(IN)  :: A(NSPEC), CG(NK)
    REAL, INTENT(OUT) :: HSIG, T01, S, SM
    REAL              :: M0, M1, M2, EB(NK), EBNK, FTE2
    INTEGER           :: IK, ITH
    !
    M0 = 0.0  ! 0th spectral moment (eqv to total sea surface variance)
    M1 = 0.0  ! 1st spectral moment
    M2 = 0.0  ! 2nd spectral moment
    T01 = 0.0
    HSIG = 0.0
    S = 0.0
    SM = 0.0
    !
    DO IK=1, NK
      EB(IK) = 0.0
      DO ITH=1, NTH
        EB(IK) = EB(IK) + A(ITH+(IK-1)*NTH)
      END DO
      EB(IK) = EB(IK) * DDEN(IK) / CG(IK)
      M0 = M0 + EB(IK)
      M1 = M1 + EB(IK) * SIG(IK)
      M2 = M2 + EB(IK) * SIG(IK) * SIG(IK)
    END DO
    !
    ! Add tail (beyond the discrete part of the spectrum)
    FTE2 = FTE/DTH/SIG(NK)
    EBNK = EB(NK) / DDEN(NK)
    M0 = M0 + EBNK * FTE2
    M1 = M1 + EBNK * SIG(NK) * FTE2 * (0.3333/0.25)
    M2 = M2 + EBNK * SIG(NK) * SIG(NK) * FTE2 * (0.5/0.25)

    M1 = M1 * TPIINV
    M2 = M2 * TPIINV*TPIINV

    T01 = M0 / MAX(M1, 1.0E-7)
    HSIG = 4.0 * SQRT( M0 )
    S = TPI * HSIG * M2 / (M0 * GRAV)
    SM = (TPI*M1)**4 * (M0**(-3)) / (GRAV*GRAV)
    !
    RETURN
    !
  END SUBROUTINE SPCPAR
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Calculate data assimilation weights based on
  !> Voorrips et al. (1997): "Assimilation of wve spectra
  !>   from pitch-and-roll buoys in a North Sea wave model",
  !>   JGR Oceans, 102(C3), 5829-5849.
  !>
  !> @param[in] RASMKM
  !> @param[in] LAT1
  !> @param[in] LAT2
  !>
  !> @param[out] WS
  !>
  !> @author  @date
  !>
  SUBROUTINE CALC_DA_VOORRIPS (RASMKM, WS)
    !/
    IMPLICIT NONE
    REAL(8), INTENT(IN)   :: RASMKM
    REAL, INTENT(OUT)     :: WS
    REAL(8)               :: L, W
    !
    L = 200.0
    !
    W = RASMKM / L
    WS = REAL(EXP(-W)**(3.0/2.0))
    !
  END SUBROUTINE CALC_DA_VOORRIPS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Calculate data assimilation weights based on
  !> Formula 7.11 in PhD thesis Greenslade (2003).
  !> Greenslade & Young (2004): "Background errors in a global
  !>   wave model determined from altimeter data", JGR Oceans,
  !>   109(C9), doi:10.1029/2004JC002324
  !>
  !> @param[in] RASMKM
  !> @param[in] LAT1
  !> @param[in] LAT2
  !>
  !> @param[out] WS
  !>
  !> @author  @date
  !>
  SUBROUTINE CALC_DA_DIANA (RASMKM, LAT1, LAT2, WS)
    !/
    IMPLICIT NONE
    REAL(8), INTENT(IN)   :: LAT1, LAT2, RASMKM
    REAL, INTENT(OUT)     :: WS
    REAL(8)               :: L1, L2, L12, W
    !
    L1 = 650.0 - 5.5 * ABS(LAT1)
    L2 = 650.0 - 5.5 * ABS(LAT2)
    L12 = SQRT(L1 * L2)
    !
    W = RASMKM / L12
    WS = REAL((1+W)*EXP(-W))
    !
  END SUBROUTINE CALC_DA_DIANA
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Calculate data assimilation weights based on
  !> structure no. 11 in Greenslade and Young (2004).
  !>
  !> @param[in] WS
  !> @param[in] HSWW3
  !> @param[in] HSBOUY
  !> @param[in] HSWW3b
  !> @param[in] TMWW3
  !> @param[in] TMBOUY
  !> @param[in] TMWW3b
  !>
  !> @param[out] W
  !> @param[out] W2
  !> @param[out] TMTRU
  !> @param[out] HSTRU
  !>
  !> @author  @date
  !>
  SUBROUTINE CALC_WEIGHT ( WS, HSWW3, HSBOUY, HSWW3b, &
                               TMWW3, TMBOUY, TMWW3b, &
                                  W, W2, TMTRU, HSTRU )
    !/
    IMPLICIT NONE
    REAL, INTENT(IN)   :: WS, HSWW3, HSBOUY, HSWW3b, TMBOUY, TMWW3b
    REAL, INTENT(OUT)  :: W, W2, TMTRU, HSTRU
    REAL               :: TMWW3
    !
    HSTRU = HSWW3 + WS * (HSBOUY - HSWW3b)
    TMTRU = TMWW3 + WS * (TMBOUY - TMWW3b)
    TMTRU = MIN(MAX(1.0, TMTRU), 25.0)
    TMWW3 = MIN(MAX(1.0, TMWW3), 25.0)
    !
    IF ( HSWW3 > 0.01 .AND. HSTRU > 0.01 ) THEN
      W2 = MIN(3.0, HSTRU / HSWW3)
      W = MIN(3.0, TMTRU / TMWW3)
    ELSE
      W2 = 1.0
      W = 1.0
    END IF
  !
  END SUBROUTINE CALC_WEIGHT
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Scale the spectrum for weight calculated from assimilation of Hs.
  !>
  !> @details Similar to UPD2 from ww3_uprstr
  !>
  !> @param[in] JSEA
  !> @param[in] W2
  !>
  !> @author  @date
  !>
  SUBROUTINE DA_HS_UPD2 ( JSEA, W2 )
    !/
    USE W3WDATMD, ONLY: VA
    USE W3GDATMD, ONLY: NSPEC
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: JSEA
    REAL, INTENT(IN)    :: W2
    INTEGER             :: I
    REAL                :: FACT
    !
    FACT = W2**2
    DO I=1, NSPEC
      VA(I,JSEA) = VA(I,JSEA)*FACT
    END DO
    !
  END SUBROUTINE DA_HS_UPD2
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Get ISEA for DA.
  !>
  !> @param[in] X longitude of observation point.
  !> @param[in] Y latitude of observation point.
  !>
  !> @param[out] ISEA index for DA location.
  !>
  !> @author  @date 
  !>
  SUBROUTINE INIT_GET_DATISEA ( X, Y, ISEA )
  !/
    USE W3GDATMD, ONLY: FLAGLL, NX, NY, MAPFS, MAPSTA, XGRD, YGRD
    USE W3GSRUMD, ONLY: W3DIST
  !/
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: X, Y
    INTEGER, INTENT(OUT)     :: ISEA
    REAL(KIND=8)             :: D, DIST, DMIN
    INTEGER                  :: IX, IY
  !/
    ISEA = -999
    D = HUGE(D)
    IF (FLAGLL) THEN
       DMIN = 1.0
    ELSE
       DMIN = 100.E3
    END IF
  !/
    DO IX=1, NX
      DO IY=1, NY
        IF ( MAPSTA(IY,IX) .LT. 0 ) CYCLE
        IF ((ABS(X-XGRD(IY,IX)).LT.DMIN).AND. &
            (ABS(Y-YGRD(IY,IX)).LT.DMIN)) THEN
          DIST = W3DIST(FLAGLL,X,Y,XGRD(IY,IX),YGRD(IY,IX))
          IF (DIST .LT. D) THEN
            D = DIST
            ISEA = MAPFS(IY,IX)
          END IF
        END IF
      END DO
    END DO
  !/
  END SUBROUTINE INIT_GET_DATISEA
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Modify spectrum by stretching and scaling following
  !>        method by Lionello et al. (1992). A new spectrum FN is
  !>        build from spectrum F in the form of
  !>             FN(sigma,theta) = A F(B*sigma,theta)
  !>
  !>        A swell dominant spectrum is updated using the
  !>        steepness criteria with a small correction applied:
  !>             DELTA=1-.006*(HSAN-HS)
  !>             A = DELTA*(ETAN/ETOT)**1.25
  !>             B = DELTA*(ETAN/ETOT)**.25
  !>        A wind sea dominant spectrum is computed with
  !>             B = FMSEA/FMSEAAN
  !>             A = (ETAN/ETOT)*B
  !>
  !> Lionello et al. (1992): "Assimilation of altimeter data in a global
  !>     third-generation wave model", JGR, 97(9), 14,453-14,474
  !>
  !> @param[inout] A    Action density spectrum A(k,theta)
  !> @param[in] CG      Group velocities
  !> @param[in] HS      Significant wave height (model guess)
  !> @param[in] HSAN    Analysis significant wave height
  !> @param[in] FM      Wind sea mean frequency (model guess)
  !> @param[in] FMAN    Analysis of wind sea mean frequency
  !>
  !> @author  @date
  !>
  SUBROUTINE UPSPEC ( A, CG, HS, HSAN, FM, FMAN )
  !/
    USE W3GDATMD, ONLY: NSPEC, NK, NTH, SIG, XFR
  !/
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: A(NSPEC)
    REAL, INTENT(IN)    :: HS, HSAN, FM, FMAN, CG(NK)
    INTEGER             :: IK, ITH, I1, I2, IKTH
    REAL                :: XHS, XR, XB, XL
    REAL                :: SPCUP(NSPEC), SPC1, SPC2
    REAL                :: SU, DSU, DELTA, XFAC, FAC1, FAC2
  !/
    XHS = HSAN/HS
  !/
    IF (FM.LT.1.0E-4) THEN
      ! The spectrum is mainly swell (wind sea mean
      ! frequency of the spectrum is negative)
      DELTA = 1.0 - 6.0E-3 * (HSAN-HS)
      XR = DELTA * (XHS**(2.5))
      XB = DELTA * SQRT(XHS)
    ELSE
      ! The spectrum is mainly wind sea
      XB = FM / FMAN
      XR = (XHS*XHS) * XB
    END IF
  !/
    SPCUP(1:NSPEC) = 0.0
    XL = ALOG(XFR)
  !/
    DO IK = 1, NK
      SU = SIG(IK) * XB
      I1 = INT(ALOG(SU/SIG(1))/XL) + 1
      I2 = I1 + 1
      IF (1.LE.I1 .AND. I2.LE.NK) THEN
        DSU = (SU - SIG(I1)) / (SIG(I2) - SIG(I1))
        ! Jacobians for spectral conversion (wn to sigma)
        FAC1 = SIG(I1)/CG(I1)
        FAC2 = SIG(I2)/CG(I2)
        XFAC = XR * CG(IK)/SIG(IK) ! scale and transform back
        DO ITH=1, NTH
          SPC1 = A(ITH + (I1-1)*NTH) * FAC1
          SPC2 = A(ITH + (I2-1)*NTH) * FAC2
          DELTA = DSU * (SPC2-SPC1)
          IKTH = ITH + (IK-1)*NTH
          SPCUP(IKTH) = MAX(0.0, SPC1 + DELTA) * XFAC
        END DO
      END IF
    END DO
  !/
    A(1:NSPEC) = SPCUP(1:NSPEC)
  !/
  END SUBROUTINE UPSPEC
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Find duration (time) form model estimate of the
  !>        wind sea partition (first guess) using nondimensional
  !>        growth curves (i.e., energy-duration)
  !>            E* = E0 tanh(A (t*)**B)
  !>        where nondimensional energy E* is given as
  !>            E* = E grav**2 / ust**4
  !>        nondimensional time t* is given as
  !>            t* = t grav/ust.
  !>
  !> @param[in]  UST    Friction velocity
  !> @param[in]  HS     Significant wave height (wind sea)
  !> @param[out] T      Duration (time in units of seconds)
  !>
  !> @author  @date
  !>
  SUBROUTINE WSDUR (UST, HS, T)
    !/
    USE CONSTANTS, ONLY: GRAV
    !/
    IMPLICIT NONE
    REAL, INTENT(IN)  :: UST, HS
    REAL, INTENT(OUT) :: T
    REAL              :: E, EST, TST
    !
    E = (HS/4)*(HS/4)
    EST = E * GRAV * GRAV / (UST**4) ! n.d. energy
    TST = MIN(EST / E0, 0.999999)    ! n.d. time
    TST = (1.0 + TST) / (1.0 - TST)  ! TANH
    TST = (0.5 * ALOG(TST) / A)**(1.0 / B)
    !
    T = TST * UST / GRAV
    !
    RETURN
    !
  END SUBROUTINE WSDUR
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Assuming that the estimate of the duration from first
  !>        guess friction velocity is correct, one can estimate
  !>        the dimensioness energt provided by the analysis and
  !>        find the corresponding analysis friction velocity
  !>        (ustan). The analysis value of the mean frequency
  !>        is derived using computed ustan and estan in
  !>        nondimensional frequency energy growth curve.
  !>        (Lionello et al., 1992)
  !>
  !> @param[in]  UST    Friction velocity
  !> @param[in]  T      Duration
  !> @param[in]  HSAN   Analysis significant wave height (wind sea)
  !> @param[out] USTAN  Analysis friction velocity
  !> @param[out] FMAN   Analysis mean frequency (wind sea)
  !>
  !> @author  @date
  !>
  SUBROUTINE WSANA (UST, T, HSAN, USTAN, FMAN)
    !/
    USE CONSTANTS, ONLY: GRAV
    !/
    IMPLICIT NONE
    REAL, INTENT(IN)  :: UST, T, HSAN
    REAL, INTENT(OUT) :: USTAN, FMAN
    REAL, PARAMETER   :: RTOL = 1.0E-5
    REAL, PARAMETER   :: ATOL = 1.0E-8
    INTEGER, PARAMETER :: ITERMAX = 12
    INTEGER           :: I
    REAL              :: E, EST, FST, T2, X2, RI, INC
    !
    FMAN = -999.9
    USTAN = UST
    E = (HSAN/4)*(HSAN/4)
    INC = 1.0
    !
    DO I=1, ITERMAX
      X2 = A*(GRAV * T / USTAN)**B
      IF (X2.GE.2.0) THEN
        USTAN = (E * GRAV * GRAV / E0)**(0.25)
      ELSE
        T2 = TANH(X2)
        RI = E * GRAV * GRAV / (USTAN**4) / E0
        INC = 1.0 - ((RI-T2)/(B*X2*(1.0-T2*T2)-4.0*T2))
      END IF
      IF (ISCLOSE(INC, 1.00000, RTOL, ATOL)) EXIT
      USTAN = USTAN * INC
    END DO
    !
    EST = E * GRAV * GRAV / (USTAN**4)
    FST = (EST / AF)**(1.0 / BF)
    FMAN = FST * GRAV / USTAN
    !
    RETURN
    !
  END SUBROUTINE WSANA
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Compare two floating-point numbers for approximate equality.
  !>
  !> @param[in]  A       floating point value 1
  !> @param[in]  B       floating point value 2
  !> @param[in]  REL_TOL relative tolerance
  !> @param[in]  ABS_TOL absolute tolerance
  !> @param[out] C       LOGICAL
  FUNCTION ISCLOSE(A, B, REL_TOL, ABS_TOL) RESULT(C)
    IMPLICIT NONE
    REAL, INTENT(IN) :: A, B
    REAL, INTENT(IN) :: REL_TOL, ABS_TOL
    LOGICAL :: C
    C = ABS(A - B) < MAX(ABS_TOL, REL_TOL * MAX(ABS(A), ABS(B)))
  END FUNCTION ISCLOSE
  !/ ------------------------------------------------------------------- /
  !/
  !/ End of module W3WDA1MD -------------------------------------------- /
  !/
END MODULE W3WDA1MD
