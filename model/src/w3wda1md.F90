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
  PRIVATE :: CALC_WAVE_PARAM, CALC_DA_DIANA, CALC_WEIGHT, &
             DA_HS_UPD2, INIT_GET_DATISEA
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
    USE CONSTANTS, ONLY: RADIUS, TPI, DERA, UNDEF
    USE W3ADATMD, ONLY: CG
#ifdef W3_MPI
    USE W3ADATMD, ONLY: MPI_COMM_WAVE
#endif
    USE W3WDATMD, ONLY: VA
    USE W3GDATMD, ONLY: NK, NTH, NSEAL, DDEN, SIG, MAPSF, &
                        MAPSTA, FLAGLL, XGRD, YGRD, ICLOSE
    USE W3PARALL, ONLY: INIT_GET_ISEA, INIT_GET_JSEA_ISPROC
    USE W3GSRUMD, ONLY: W3DIST
    USE W3ODATMD, ONLY: NDSO, NDSE, NDST, SCREEN, NAPROC, IAPROC, &
                        NAPLOG, NAPOUT, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
    !
    IMPLICIT NONE
    !
#ifdef W3_MPI
    INCLUDE "mpif.h"
#endif
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
    INTEGER                 :: IDAT, IX, IY, ISEA, JSEA, &
                               MREC, IK, ITH, DATPROC
    REAL, PARAMETER         :: HSMIN = 0.01
    REAL                    :: HSDAT, TMDAT
    REAL                    :: HS, TM
    REAL                    :: TMTRU, HSTRU, W, W2, WS
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
    ! 2.  Actual data assimilation routine ------------------------------- /
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
        CALL CALC_WAVE_PARAM(JSEA, ISEA, HSDAT, TMDAT)
      ! WRITE(*,'(2X,A8,4I6,2F7.2)') "DA-PROC", IAPROC, DATPROC, ISEA, JSEA, HSDAT, TMDAT
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
        IF (DKM .LT. 4000.0) THEN
          CALL CALC_WAVE_PARAM(JSEA, ISEA, HS, TM)
         !WRITE(*,'(2X,A8,4I5,2F7.2,F10.2)') "MOD", IAPROC, NAPROC, ISEA, JSEA, HS, TM, DKM
          CALL CALC_DA_DIANA(DKM, YGRD(IY,IX), YDAT, WS)
          CALL CALC_WEIGHT(WS, HS, DATA0(3,IDAT), HSDAT, &
                               TM, DATA0(4,IDAT), TMDAT, &
                            W, W2, TMTRU, HSTRU)
          CALL DA_HS_UPD2(JSEA, W2)
         !IF (IAPROC .EQ. DATPROC) THEN
         !  CALL CALC_WAVE_PARAM(JSEA, ISEA, HS, TM)
         !  WRITE(*,'(2X,A8,4I6,2F7.2)') "DATPROC", IAPROC, DATPROC, ISEA, JSEA, HS, TM
         !END IF !/ IAPROC == DATPROC
        END IF !/ (DKM .LT. 4000.0)
      END DO !/ JSEA..NSEAL
    END DO  !/ IDAT, NDAT
    !
    RETURN
    !
    ! Formats
    !
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3WDA1 :'/                &
              '     DATA RECORD DIMENSION <4 : ',I)
  END SUBROUTINE W3WDA1
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Get model parameters at observation point.
  !>
  !> @param[in] JSEA
  !> @param[in] ISEA
  !>
  !> @param[out] HSIG
  !> @param[out] TMEAN
  !>
  !> @author  @date
  !>
  SUBROUTINE CALC_WAVE_PARAM ( JSEA, ISEA, HSIG, TMEAN)
    !/
    USE CONSTANTS, ONLY: TPI
    USE W3ADATMD, ONLY: CG
    USE W3WDATMD, ONLY: VA
    USE W3GDATMD, ONLY: NK, NTH, DDEN, SIG
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: JSEA, ISEA
    REAL, INTENT(OUT)   :: HSIG, TMEAN
    REAL                :: EMEAN
    REAL, ALLOCATABLE   :: EB(:)
    INTEGER             :: IK, ITH
    !
    ALLOCATE(EB(NK))
    DO IK=1, NK
      EB(IK) = 0.0
      DO ITH=1, NTH
        EB(IK) = EB(IK) + VA(ITH+(IK-1)*NTH,JSEA)
      END DO
      EB(IK) = EB(IK) * DDEN(IK) / CG(IK,ISEA)
    END DO

    EMEAN = 0.0
    TMEAN = 0.0
    DO IK = 1, NK
        EMEAN = EMEAN + EB(IK)
        TMEAN = TMEAN + EB(IK) * SIG(IK)
    END DO
    DEALLOCATE(EB)

    TMEAN = TPI * EMEAN / MAX(TMEAN, 1.E-7)
    HSIG = 4.0 * SQRT( EMEAN )
    !
  END SUBROUTINE CALC_WAVE_PARAM
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief Calculate data assimilation weights based on
  !> Formula 7.11 in PhD thesis Greenslade (2003).
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
  !> structure no. 11 in Greenlase and Young (2004).
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
  !> @param[in] XBOUY Longitude of observation point.
  !> @param[in] YBOUY Latitude of observation point.
  !>
  !> @param[out] DASISEA ISEA for DA location.
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
  !/
  !/ End of module W3WDA1MD -------------------------------------------- /
  !/
END MODULE W3WDA1MD
