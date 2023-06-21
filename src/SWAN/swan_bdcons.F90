#include "swan_functions.h"
#define DEBUG
#undef DEBUG
!**********************************************************************
!*    
! Update may 2021: add a new boundary treatment to impose Mono or Multi
! Modal wave parameters along OBC nodes 
!                                                                     *
!**********************************************************************
      SUBROUTINE COMPUTE_IFILE_IT(JFILE, IT)

      USE VARS_WAVE
      USE TIMECOMM
      USE schism_glbl, ONLY: rkind
      USE schism_msgp
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: JFILE, IT
      REAL(rkind) :: DTMP
      INTEGER ITMP

      DTMP = TIMCO - BND_TIME_ALL_FILES(1,1)*DAY2SEC

      ITMP  = 0
      DO JFILE = 1, NUM_NETCDF_FILES_BND
        ITMP = ITMP + NDT_BND_FILE(JFILE)
        IF (ITMP .GT. INT(DTMP/WBDELT)) EXIT
      END DO

      ITMP = SUM(NDT_BND_FILE(1:JFILE-1))
      IT   = NINT(DTMP/WBDELT) - ITMP + 1

      IF (IT .GT. NDT_BND_FILE(JFILE)) THEN
        JFILE = JFILE + 1
        IT    = 1
      ENDIF

      IF(JFILE.GT.NUM_NETCDF_FILES_BND) &
         CALL parallel_abort("swan_bdcons, Error in COMPUTE_IFILE_IT")

      END SUBROUTINE

!**********************************************************************
!*    T. Guerin: Equivalent of COMPUTE_IFILE_IT for WW3 binary files  *
!**********************************************************************
      SUBROUTINE COMPUTE_IT(IT)

      USE VARS_WAVE
      USE TIMECOMM
      USE schism_glbl, ONLY: rkind
      USE schism_msgp
      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: IT
      REAL(rkind) :: DTMP
      INTEGER ITMP
      ! LOCAL should be global

      if(myrank.eq.0) print*,'COMPUTE_IT',WBTMJD,BND_TIME_ALL_FILES(1,1),DAY2SEC,WBDELT
      if(myrank.eq.0) print*,'COMPUTE_IT',LBINTER

      DTMP = TIMCO - BND_TIME_ALL_FILES(1,1)*DAY2SEC
      IT   = NINT(DTMP/WBDELT)+ 1
      IF(LBINTER) IT = IT + 1

      if(myrank.eq.0) print*,'COMPUTE_IT',IT,DTMP

      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!*  Called by SET_WAVE_BOUNDARY_CONDITION
!*
!**********************************************************************
      SUBROUTINE WAVE_BOUNDARY_CONDITION !(WBACOUT)

      USE SWCOMM3, ONLY:MDC,MSC,SPPARM
      USE OCPCOMM4
      USE M_GENARR
      USE VARS_WAVE
      USE schism_glbl, ONLY : rkind
      USE schism_msgp, ONLY : parallel_abort,myrank
      IMPLICIT NONE

!      REAL(rkind), INTENT(OUT)   :: WBACOUT(MDC,MSC,IOBCN_W)
      REAL(rkind),allocatable    :: WBAC_WAM(:,:) ! used only in WAM multi-mode forcing
      INTEGER                    :: IP, K
      INTEGER                    :: FSHAPE,DSHAPE
     
      INTEGER :: IP1, IP2,IPNEW,IPOLD ! quick fix to manage Nan

!AR: WAVE BOUNDARY

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)

!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
!
! Count number of active boundary points ...
!
!      IF(LWW3GLOBALOUT) THEN
!        IF (.NOT. ALLOCATED(WW3GLOBAL)) THEN
!          ALLOCATE(WW3GLOBAL(8,MNP), stat=istat)
!          IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 7')
!        END IF
!      END IF
      IF (TRIM(NESTING_TYPE_WAVE).EQ.'PARAM') THEN ! Parametric Wave Boundary is prescribed
!        WRITE(PRINTF,*) &
!       'Parametric Wave Boundary Condition is prescribed, only WW3'

        WBSPPARM = 0.

        CALL READ_NETCDF_WW3_PARAM
        CALL INTER_STRUCT_BOUNDARY(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,WBSPPARM)

!       IF (LWW3GLOBALOUT) &
!       CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,WW3GLOBAL)

        IPOLD = IOBCN_W
        DO IP = 1, IOBCN_W
!          CALL SPECTRAL_SHAPE(WBSPPARM(:,IP),WBACOUT(:,:,IP),.FALSE.,'CALL FROM WB 1', USE_OPTI_SPEC_SHAPE_BOUC)
           
           IF(WBSPPARM(2,IP) .GT. 0._rkind) IPOLD = IP

!          Quick Fix for NaN
!           IP1 = IP-1; 
!           IF(IP1.LE.0) IP1 = IOBCN_W
!           IP2 = IP+1; 
!           IF(IP2.GT.IOBCN_W) IP2 = 1

!           IF(WBSPPARM(1,IP1)/=WBSPPARM(1,IP1)) THEN
!             IP1 = IP1-1; 
!             IF(IP1.LE.0) IP1 = IOBCN_W
!           ENDIF
!           IF(WBSPPARM(1,IP2)/=WBSPPARM(1,IP2)) THEN
!             IP2 = IP2+1; 
!             IF(IP2.GT.IOBCN_W) IP2 = 1
!           ENDIF 

           IPNEW = IP           
      
           IF(WBSPPARM(1,IP) /= WBSPPARM(1,IP)) THEN ! NaN, take one neighbors node
            if(myrank.EQ.0) print*,'NAN value in OBC!'
!            IF(WBSPPARM(1,IP2).EQ.WBSPPARM(1,IP2)) THEN
!             IPNEW = IP2
!            ELSEIF(WBSPPARM(1,IP1).EQ.WBSPPARM(1,IP1)) THEN
!             IPNEW = IP1
!            ELSE
!             CALL parallel_abort('NaN in OBC')
!            ENDIF
            IPNEW = IPOLD   
           ELSE IF(WBSPPARM(2,IP) .EQ. 0._rkind) THEN 
            if(myrank.EQ.0) &
            print*,'swan_bdcons: Fake value in OBC!'
            IPNEW = IPOLD
           ENDIF

           IF(WBSPPARM(2,IP).EQ.0._rkind) THEN
            CALL parallel_abort('swan_bdcons: error')
           ENDIF

!          Use the SWAN spctral shape scheme, SSHAPE,
           SPPARM(1) = WBSPPARM(1,IPNEW)  !HSC1(IOB_NODE)
           SPPARM(2) = WBSPPARM(2,IPNEW)  !TPEAK(IOB_NODE)
           SPPARM(3) = WBSPPARM(3,IPNEW)  !DIRDEG1(IOB_NODE)
           SPPARM(4) = WBSPPARM(4,IPNEW)  !Directional Spreading, 20.0
           FSHAPE = INT(WBSPPARM(5,IPNEW))!spectral shape (1-4),
           FSHAPE = ABS(FSHAPE) ! A verifier 2 fois 
           DSHAPE = INT(WBSPPARM(6,IPNEW))!directional spreading in degree (1) or exponent (2)

           !print*,SPPARM(1),SPPARM(2) ,SPPARM(3),SPPARM(4),FSHAPE,DSHAPE
           CALL SSHAPE(WBAC(1,1,IP),SPCSIG,SPCDIR,FSHAPE,DSHAPE)

        END DO

      ELSE IF (TRIM(NESTING_TYPE_WAVE).EQ.'WAMPARAM') THEN ! Multi-Modal (WAM) 
               ! Parametric Wave Boundary is prescribed:
               ! Cumulate each wave partition and fill the wave Action Density 

        if(.not.allocated(WBAC_WAM)) allocate(WBAC_WAM(MDC,MSC))

        ! Reset ?
        WBAC(1:MDC,1:MSC,1:IOBCN_W) = 0._rkind  


        CALL GET_NETCDF_WAM_PARAM

        DO K=1,N_WAVE_WAM

           WBSPPARM = 0._rkind
           ! Get wave and Bulk spectrum function parameters
           ! for that wave partition
           CALL SPPARM_WAM(K,WBSPPARM)

           IPOLD = IOBCN_W
           DO IP = 1, IOBCN_W

             IF(WBSPPARM(2,IP) .GT. 0._rkind) IPOLD = IP
             IPNEW = IP

             IF(WBSPPARM(1,IP) /= WBSPPARM(1,IP)) THEN ! NaN, take one neighbors node
               IF(myrank.EQ.0) print*,'NAN value in OBC!'
               IPNEW = IPOLD
             ELSE IF(WBSPPARM(2,IP) .EQ. 0._rkind) THEN
               if(myrank.EQ.0) &
               print*,'swan_bdcons: Fake value in OBC!'
               IPNEW = IPOLD
             ENDIF

             IF(WBSPPARM(2,IP).EQ.0._rkind) &
              CALL parallel_abort('swan_bdcons: error')

!            Use the SWAN spctral shape scheme, SSHAPE,
             SPPARM(1) = WBSPPARM(1,IPNEW)  !Hs  (IOB_NODE)
             SPPARM(2) = WBSPPARM(2,IPNEW)  !Period TPEAK(IOB_NODE)
             SPPARM(3) = WBSPPARM(3,IPNEW)  !DIR (IOB_NODE)
             SPPARM(4) = WBSPPARM(4,IPNEW)  !Directional Spreading,
             FSHAPE = INT(WBSPPARM(5,IPNEW))!spectral shape (1-4),
             FSHAPE = ABS(FSHAPE) ! A verifier 2 fois
             DSHAPE = INT(WBSPPARM(6,IPNEW))!directional spreading in degree (1) or exponent (2)

             WBAC_WAM = 0._rkind

             CALL SSHAPE(WBAC_WAM,SPCSIG,SPCDIR,FSHAPE,DSHAPE)

!            Cumulate wave density from each wave system
             WBAC(1:MDC,1:MSC,IP) = WBAC(1:MDC,1:MSC,IP) + &
     &                              WBAC_WAM(1:MDC,1:MSC)
           ENDDO

        ENDDO
        
      ELSE IF (TRIM(NESTING_TYPE_WAVE).EQ.'SPEC') THEN ! Spectrum is prescribed

       CALL GET_BINARY_WW3_SPECTRA ! Read SPEC_WW3(MDC_WW3,MSC_WW3,IP_WW3), 
                                   ! Interpolate on SPEC_SWAN(MDC,MSC,IOBCN_W),
                                   ! Update WBAC(MDC,MSC,IOBCN_W)
      END IF

      !WRITE(PRINTF,*) 'END WAVE_BOUNDARY_CONDITION'

      END SUBROUTINE WAVE_BOUNDARY_CONDITION
!**********************************************************************
!* GET_BINARY_WW3_SPECTRA
!* Read a WAVEWATCHIII binary spectral file and do time and space
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITION
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
      SUBROUTINE GET_BINARY_WW3_SPECTRA !(WBAC)

      USE VARS_WAVE, ONLY : WBAC,NP_WW3,DR_WW3,DDIR_WW3,FQ_WW3
      USE VARS_WAVE, ONLY : MSC_WW3,MDC_WW3,XP_WW3,YP_WW3
      USE VARS_WAVE, ONLY : IOBCN_W,I_OBC_N_W
      USE SWCOMM3, ONLY : MDC,MSC,SPPARM,SLOW,SHIG,PI2_W,PI_W
      USE schism_glbl, ONLY : rkind
      USE swangriddata, ONLY:xcugrd,ycugrd
      USE schism_msgp, ONLY : parallel_abort,myrank

!      USE DATAPOOL, ONLY: NP_WW3, rkind, DR_WW3, DDIR_WW3, FQ_WW3, FRLOW, LNANINFCHK, DBG, WRITEDBGFLAG, FRHIGH
!      USE DATAPOOL, ONLY: LINHOM, IWBNDLC, XP, YP, XP_WW3, YP_WW3, STAT, WRITESTATFLAG, MSC, MDC, IWBMNP, MSC_WW3
!      USE DATAPOOL, ONLY: MDC_WW3
!# ifdef MPI_PARALL_GRID
!      USE DATAPOOL, ONLY: XLON, YLAT
!# endif
      IMPLICIT NONE
      !REAL(rkind), INTENT(OUT) :: WBAC(MDC,MSC,IOBCN_W) 
!      REAL(rkind) :: WBACOUT(MDC,MSC,IOBCN_W) !WW3 take care !
      INTEGER     :: IB,IPGL,IBWW3,TIME(2),IS
      REAL(rkind) :: SPEC_WW3(MDC_WW3,MSC_WW3,NP_WW3)
      REAL(rkind) :: SPEC_SWAN(MDC,MSC,NP_WW3)
      REAL(rkind) :: DIST(NP_WW3),TMP(NP_WW3),INDBWW3(NP_WW3)
      REAL(rkind) :: SPEC_WW3_TMP(MDC_WW3,MSC_WW3,NP_WW3)
      REAL(rkind) :: SPEC_WW3_UNSORT(MDC_WW3,MSC_WW3,NP_WW3)
      REAL(rkind) :: JUNK(MDC_WW3),DR_WW3_TMP(MDC_WW3)
      REAL(rkind) :: XPT,YPT

      ! SHIG   [CALCUL] =2*PI*FRHIG; highest spectral value of sigma
      ! SLOW   [CALCUL] =2*PI*FRLOW; lowest spectral value of sigma
      REAL(rkind) :: FRHIGH,FRLOW

      INTEGER     :: IFILE, IT
!     Local should be global ...
      INTEGER     :: WRITESTATFLAG,LINHOM
      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

      LINHOM = 1        ! Do spatial interp
      WRITESTATFLAG = 1 ! Print DBG info

!     STRACE           Tracing routine for debugging
      CALL STRACE (IENT,'GET_BINARY_WW3_SPECTRA')

      !CALL COMPUTE_IFILE_IT(IFILE, IT)
      CALL COMPUTE_IT(IT)

      FRHIGH = SHIG/PI2_W
      FRLOW  = SLOW/PI2_W
!
! Read spectra in file
!
      CALL READ_SPEC_WW3(IT,SPEC_WW3_UNSORT)
!
! Sort directions and carries spectra along (ww3 directions are not
! monotonic)
!
      DO IBWW3 = 1, NP_WW3
        DO IS = 1,MSC_WW3
          DR_WW3_TMP=DR_WW3
          SPEC_WW3_TMP(:,IS,IBWW3) = SPEC_WW3_UNSORT(:,IS,IBWW3)
          CALL SSORT2 (DR_WW3_TMP,SPEC_WW3_TMP(:,IS,IBWW3),JUNK,MDC_WW3,2)
          SPEC_WW3(:,IS,IBWW3) = SPEC_WW3_TMP(:,IS,IBWW3)
          DR_WW3 = DR_WW3_TMP
        ENDDO
        DDIR_WW3 = DR_WW3(2) - DR_WW3(1)

        if(myrank.EQ.0.and.WRITESTATFLAG.eq.1) &
          print*,'GET_BINARY_WW3_SPECTRA, AFTER SORTING', IBWW3, &
     &    SUM(SPEC_WW3(:,:,IBWW3))
      ENDDO
!
! Interpolate ww3 spectra on SWAN frequency grid
! GD: at the moment 360 deg spanning grids are mandatory
!
      IF((FQ_WW3(1).GT.FRLOW).OR.(FQ_WW3(MSC_WW3).LT.FRHIGH)) THEN
        IF (myrank.EQ.0.and.WRITESTATFLAG.eq.1) THEN
          print*,'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
          print*,'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
          print*,'WW3 spectra does not encompass the whole WWM spectra, &
     & please carefully check if this makes sense for your simulations'
        END IF
        CALL SPECTRALINT(SPEC_WW3,SPEC_SWAN)
      ELSE
        IF (myrank.EQ.0.and.WRITESTATFLAG.eq.1) THEN
          print*,'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
          print*,'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
        END IF
        CALL SPECTRALINT(SPEC_WW3,SPEC_SWAN)
      ENDIF

!      IF (LNANINFCHK .AND. WRITEDBGFLAG == 1) THEN
!        WRITE(DBG%FHNDL,*) 'SUMS AFTER INTERPOLATION', SUM(SPEC_WW3),SUM(SPEC_WWM)
!      ENDIF
!
! Interpolate ww3 spectra on swan boundary nodes
! GD: ww3 forcing works until here. Some more debugging is needed
!
      TMP = 0
      IF(LINHOM==1) THEN !nearest-neighbour interpolation
        DO IB=1,IOBCN_W
          IPGL = I_OBC_N_W(IB)
          XPT  = xcugrd(IPGL)
          YPT  = ycugrd(IPGL)

          IF (NP_WW3 .GT. 1) THEN
            DO IBWW3=1,NP_WW3
              DIST(IBWW3) = SQRT((XPT-XP_WW3(IBWW3))**2+(YPT-YP_WW3(IBWW3))**2)
              INDBWW3(IBWW3)=IBWW3
            ENDDO
            CALL SSORT2 (DIST, INDBWW3, TMP, NP_WW3, 2)
            CALL SHEPARDINT2D(2, 1./DIST(1:2),MSC,             & 
     &           MDC,SPEC_SWAN(:,:,INT(INDBWW3(1:2))),WBAC(:,:,IB),1)
!     &           MDC,SPEC_SWAN(:,:,INT(INDBWW3(1:2))),WBACOUT(:,:,IB),1)
            IF (myrank.EQ.0.and.WRITESTATFLAG.eq.1) THEN
              write(*,'(A20, 2F20.5,3F30.10)')        &
     &        'AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2),  &
     &        sum(SPEC_SWAN(:,:,INT(INDBWW3(1)))),              &
     &        sum(SPEC_SWAN(:,:,INT(INDBWW3(2)))), SUM(WBAC(:,:,IB))
!     &        sum(SPEC_SWAN(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
            END IF
          ELSE

!            WBACOUT(:,:,IB) = SPEC_SWAN(:,:,1)
            WBAC(:,:,IB) = SPEC_SWAN(:,:,1)


          ENDIF
        ENDDO ! IB

      ELSE ! LINHOM==0
        DO IB=1,IOBCN_W
!          WBACOUT(:,:,IB) = SPEC_WWM(:,:,1)
          WBAC(:,:,IB) = SPEC_SWAN(:,:,1)
        ENDDO
      ENDIF
      IF (myrank.EQ.0.and.WRITESTATFLAG.eq.1) THEN
        print*,'DONE GET_BINARY_WW3_SPECTRA'
      END IF

!      DO IB=1,IOBCN_W
!       DO IS =1,1
!        DO ID =1,1
!         WBAC = WBACOUT
!        ENDDO
!       ENDDO
!      ENDDO
       
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHEPARDINT2D(NP,WEIGHT,D1,D2,Z,ZINT,P)

      USE schism_glbl, ONLY : rkind
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: NP,P,D1,D2
      REAL(rkind), INTENT(IN)  :: WEIGHT(NP), Z(D1,D2,NP)
      REAL(rkind), INTENT(OUT) :: ZINT(D1,D2)
      INTEGER                  :: IP
      REAL                     :: SW
      SW=0
      ZINT=0
      DO IP=1,NP
        SW=SW+WEIGHT(IP)**P
        ZINT(:,:)=ZINT(:,:)+Z(:,:,IP)*(WEIGHT(IP)**P)
      ENDDO
      ZINT=ZINT/SW
      END SUBROUTINE
!**********************************************************************
!* SPECTRALINT
!* Interpolate spectrum on SWAN spectral grid.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)The spectral grid needs to be define over 360º and
!*             the min (resp. max) frequencies need to be smaller
!*             (resp. higher) than in WWM frequency grid.
!**********************************************************************
      SUBROUTINE SPECTRALINT(SPEC_WW3,SPEC_SWAN)

      USE VARS_WAVE
      USE schism_glbl, ONLY : rkind
      USE M_GENARR, ONLY : SPCDIR,SPCSIG
      USE SWCOMM3, ONLY : DDIR,PI2_W,PI_W,MDC,MSC

      IMPLICIT NONE
      REAL(rkind), INTENT(IN)  :: SPEC_WW3(MDC_WW3,MSC_WW3,NP_WW3)
      REAL(rkind), INTENT(OUT) :: SPEC_SWAN(MDC,MSC,NP_WW3)
      REAL(rkind) :: SPEC_WW3_TMP(MDC,MSC_WW3,NP_WW3)
      REAL(rkind) :: DF, M0_WW3, M1_WW3, M2_WW3, M0_SWAN, M1_SWAN, M2_SWAN
      INTEGER     :: IP,IS,ID
      REAL(rkind) :: JACOBIAN(MSC), AM, SM, FR(MSC)
      INTEGER     :: WRITESTATFLAG
      REAL(rkind), PARAMETER   :: ZERO     = 0._rkind
      REAL(rkind), PARAMETER   :: TWO      = 2._rkind
      REAL(rkind), PARAMETER   :: ONE      = 1._rkind
      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

      WRITESTATFLAG = 1 ! Print DBG info

!     STRACE           Tracing routine for debugging
      CALL STRACE (IENT,'SPECTRALINT')

      JACOBIAN = ONE/(SPCSIG*PI2_W)! ENERGY / HZ -> ACTION / RAD
      SPEC_WW3_TMP = ZERO
      SPEC_SWAN    = ZERO
      FR = SPCSIG/PI2_W

      IF (WRITESTATFLAG == 1) THEN
        WRITE(*,'(A20,I10,3F30.2)') 'BEFORE INTERPOLATION', SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_SWAN)
      ENDIF
      DO IP=1,NP_WW3
        DO IS=1,MSC_WW3
          CALL INTERLIND(MDC_WW3,MDC,DR_WW3,SPCDIR(:,1),SPEC_WW3(:,IS,IP),SPEC_WW3_TMP(:,IS,IP))
        ENDDO
        DO ID=1,MDC
          CALL INTERLIN (MSC_WW3,MSC,FQ_WW3,FR,SPEC_WW3_TMP(ID,:,IP),SPEC_SWAN(ID,:,IP))
        ENDDO
        M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
        DO ID = 1,MDC_WW3
          DO IS = 1,MSC_WW3-1
            DF = FQ_WW3(IS+1)-FQ_WW3(IS)
            AM = (SPEC_WW3(ID,IS+1,IP)+SPEC_WW3(ID,IS,IP))/TWO
            SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
            M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
            M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
            M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
          ENDDO
        ENDDO
        M0_SWAN = ZERO; M1_SWAN = ZERO; M2_SWAN = ZERO
        DO ID = 1,MDC
          DO IS = 1,MSC-1
            DF = FR(IS+1)-FR(IS)
            AM = (SPEC_SWAN(ID,IS+1,IP)+SPEC_SWAN(ID,IS,IP))/TWO
            SM = (FR(IS+1)+FR(IS))/TWO
            M0_SWAN =M0_SWAN+AM*DF*DDIR
            M1_SWAN =M1_SWAN+AM*SM*DF*DDIR
            M2_SWAN =M2_SWAN+AM*SM**2*DF*DDIR
          ENDDO
        ENDDO
        IF (WRITESTATFLAG == 1) THEN
          print*,'POINT NUMBER', IP
          WRITE(*,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_SWAN, &
                'M2 = ',M2_WW3, M2_SWAN
        END IF
      END DO

      IF (WRITESTATFLAG == 1) THEN
        WRITE(*,'(A20,I10,3F30.2)') 'AFTER INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_SWAN)
      ENDIF

! Do jacobian

      DO IP = 1, NP_WW3
        DO ID = 1, MDC
          SPEC_SWAN(ID,:,IP) = SPEC_SWAN(ID,:,IP) * JACOBIAN(:) ! convert to wave action in sigma,theta space
        END DO
      END DO

      IF (WRITESTATFLAG == 1) &
        print*,'CHECKING INTEGRATED PARAMETERS AFTER JACOBIAN'
      DO IP = 1, NP_WW3
        M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
        DO ID = 1,MDC_WW3
          DO IS = 1,MSC_WW3-1
            DF = FQ_WW3(IS+1)-FQ_WW3(IS)
            AM = (SPEC_WW3(ID,IS+1,IP)+SPEC_WW3(ID,IS,IP))/TWO
            SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
            M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
            M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
            M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
          ENDDO
        ENDDO
        M0_SWAN = ZERO; M1_SWAN = ZERO; M2_SWAN = ZERO
        DO ID = 1,MDC
          DO IS = 1,MSC-1
            DF = SPCSIG(IS+1)-SPCSIG(IS)
            SM = (SPCSIG(IS+1)+SPCSIG(IS))/TWO
            AM = (SPEC_SWAN(ID,IS+1,IP)+SPEC_SWAN(ID,IS,IP))/TWO * SM
            M0_SWAN =M0_SWAN+AM*DF*DDIR
            M1_SWAN =M1_SWAN+AM*SM*DF*DDIR
            M2_SWAN =M2_SWAN+AM*SM**2*DF*DDIR
          ENDDO
        ENDDO
        IF (WRITESTATFLAG == 1) THEN
          WRITE(*,*) 'POINT NUMBER', IP
          WRITE(*,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_SWAN, &
               'M2 = ',M2_WW3, M2_SWAN
        END IF
      END DO
      IF (WRITESTATFLAG == 1) THEN
        WRITE(*,'(A20,I10,3F30.2)') 'AFTER JACOBIAN', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_SWAN)
      ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTERLIN (NX1, NX2, X1, X2, Y1, Y2)

      USE schism_glbl, ONLY : rkind
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)                :: DX1(NX1-1)
      REAL(rkind), PARAMETER     :: THR = TINY(1.)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1) - X1(I)
      END DO

      DO I = 1, NX2
        DO J = 1, NX1 - 1
          IF (ABS(X2(I) - X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J) + (Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
          END IF
        END DO
      END DO

      END SUBROUTINE
!**********************************************************************
!* INTERLIND
!* Interpolate vector on a 2-pi periodic axis (directions in wwm grid).
!*
!* Called by SPECTRALINT
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)This routine should be togeher with INTERLIN either here
!*             or in wwm_aux.F90.
!**********************************************************************
      SUBROUTINE INTERLIND (NX1, NX2, X1, X2, Y1, Y2)
 
      USE schism_glbl, ONLY : rkind
      USE SWCOMM3, ONLY : PI2_W,PI_W
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), PARAMETER     :: THR = TINY(1.)

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)         :: DX1(NX1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1)-X1(I)
      END DO
      DX1(NX1) = X1(1)+PI2_W-X1(NX1)

      DO I = 1, NX2
        DO J = 1,NX1 - 1
          IF (ABS(X2(I)-X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
            EXIT
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J)+(Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
            EXIT
          ELSE IF (X2(I) .GT. X1(NX1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)-X1(NX1))
            EXIT
          ELSE IF (X2(I) .LT. X1(1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)+PI2_W-X1(NX1))
            EXIT
          END IF
        END DO
      END DO
      END SUBROUTINE


!**********************************************************************
!* READ_SPEC_WW3
!* (doing MPI exchanges if MULTIPLE_IN_BOUND = FALSE)
!* (otherwise, all node read the same data)
!*
!* Called by GET_BINARY_WW3_SPECTRA or GET_NC_WW3_SPECTRA
!*
!* Authors: Mathieu Dutour Sikiric
!**********************************************************************
      SUBROUTINE READ_SPEC_WW3(ISTEP,SPECOUT)

      USE VARS_WAVE, ONLY : MSC_WW3,MDC_WW3,NP_WW3,SERIAL
      USE schism_glbl, ONLY : rkind
      USE schism_msgp

      IMPLICIT NONE

      include 'mpif.h'

      INTEGER, INTENT(IN) :: ISTEP
!      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)
      REAL(rkind), INTENT(OUT) :: SPECOUT(MDC_WW3,MSC_WW3,NP_WW3)
      integer, allocatable :: send_rqst(:)
      integer, allocatable :: send_stat(:,:)
      integer :: siz, iProc, istat
      INTEGER :: IENT
      SAVE  IENT
      DATA  IENT /0/

!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'READ_SPEC_WW3')

      IF(.NOT.SERIAL) THEN

        siz=MDC_WW3*MSC_WW3*NP_WW3  ! SPECOUT is /MDC_WW3,MSC_WW3,NP_WW3

        IF (myrank.eq.0) THEN

          allocate(send_rqst(nproc-1),send_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) &
             CALL parallel_abort("READ_SPEC_WW3, allocate error 30")

          CALL READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)

          DO iProc=2,nproc
            CALL mpi_isend(SPECOUT, siz, rtype, iProc-1, 2034, comm, send_rqst(iProc-1), istat)
          END DO
          IF (nproc .gt. 1) THEN
            call mpi_waitall(nproc-1, send_rqst, send_stat,istat)
          END IF
          deallocate(send_rqst, send_stat)

        ELSE
          CALL MPI_RECV(SPECOUT, siz, rtype, 0, 2034, comm, istatus, istat)
        END IF

      ELSE ! SERIAL

       CALL READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)

      ENDIF ! SERIAL

      END SUBROUTINE
!**********************************************************************
!* READ_SPEC_NC_WW3_KERNEL, Called by READ_SPEC_WW3

!* Reads spectra in WAVEWATCHIII netcdf spectral file

!* Authors: Kevin MARTINS
!**********************************************************************
      SUBROUTINE READ_SPEC_NC_WW3_KERNEL(ISTEP,SPECOUT)

      USE VARS_WAVE
      USE schism_msgp, ONLY : parallel_abort,myrank
      USE NETCDF

      IMPLICIT NONE
      REAL(rkind), INTENT(OUT) :: SPECOUT(MDC_WW3,MSC_WW3,NP_WW3)
      INTEGER, INTENT(IN) :: ISTEP
      REAL :: SPEC_WW3(MDC_WW3,MSC_WW3)
      INTEGER :: IP, BND_NCID, ENERGY_VAR_ID, DD, FF, ISTAT
!     LOCAL but should be global
      INTEGER :: LBCSE  ! 1= Non-Stationnary LBC, 0 STATIONNARY

      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

      LBCSE = 1

!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'READ_SPEC_NC_WW3_KERNEL')

      IF (myrank .eq. 0) THEN

       ! Opening netcdf file
       ISTAT = NF90_OPEN(NESTING_FILE,NF90_NOWRITE,BND_NCID)
       IF (ISTAT /= 0) THEN
           Print *, 'Error while trying to open netcdf file'
           Print *, 'FILE=', TRIM(NESTING_FILE)
           Print *, 'One possible error is that the file is NC4 but'
           Print *, 'you linked to NC3'
           CALL parallel_abort("Error while trying to open netcdf file")
       END IF

       ! Loading the energy density variable ID
       ISTAT = NF90_INQ_VARID(BND_NCID,'efth',ENERGY_VAR_ID)
       if(ISTAT.ne.NF90_NOERR) &
          call parallel_abort("Netcdf Inquire varid error, efth")

!       ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
!          if(ISTAT.ne.NF90_NOERR) &
!            call parallel_abort("Netcdf Inquire varid error")

!          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
!          if(ISTAT.ne.NF90_NOERR) &
!            call parallel_abort("Netcdf Inquire var error")

!          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), &
!                  len = NDT_BND_FILE(JFILE))
!          if(ISTAT.ne.NF90_NOERR) &
!            call parallel_abort('Netcdf Inquire dim error')

!          WRITE(PRINTF,*) 'nb records IFILE=',JFILE, NDT_BND_FILE(JFILE)


       ! Reading spectra
       IF(LBCSE==1) THEN ! non-stationary ...
       DO IP = 1,NP_WW3
          !print*,IP,ISTEP,MDC_WW3,MSC_WW3
          ISTAT = NF90_GET_VAR(BND_NCID,ENERGY_VAR_ID,SPEC_WW3,start=(/1,1,IP,ISTEP/),count = (/MDC_WW3,MSC_WW3,1,1/))
          if(ISTAT.ne.NF90_NOERR) &
             call parallel_abort("ERROR WHILE READING SPECTRA, 3")
          !! Re-ordering spectrum (WW3 is [MDC_WW3,MSC_WW3] while WWM is [MSC_WW3,MDC_WW3])
          !!SPECOUT(:,:,IP) = TRANSPOSE(SPEC_WW3)
          SPECOUT(:,:,IP) = SPEC_WW3(:,:)
       ENDDO ! IP

       ELSE ! stationary ...
        DO IP = 1,NP_WW3
          ISTAT = NF90_GET_VAR(BND_NCID,ENERGY_VAR_ID,SPEC_WW3,start=(/1,1,IP,1/),count = (/MDC_WW3,MSC_WW3,1,1/))
          if(ISTAT.ne.NF90_NOERR) &
             call parallel_abort("ERROR WHILE READING SPECTRA, 4")
          ! Re-ordering spectrum (WW3 is [MDC_WW3,MSC_WW3] while WWM is [MSC_WW3,MDC_WW3])
          !SPECOUT(:,:,IP) = TRANSPOSE(SPEC_WW3)
          SPECOUT(:,:,IP) = SPEC_WW3(:,:)
        END DO
       ENDIF

      ISTAT = NF90_CLOSE(BND_NCID)
      if(ISTAT.ne.NF90_NOERR) &
         call parallel_abort("ERROR WHILE CLOSING NETCDF FILE, 5")

      ENDIF
      END SUBROUTINE READ_SPEC_NC_WW3_KERNEL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY_CONDITION

      USE OCPCOMM4
      USE VARS_WAVE
      USE TIMECOMM, ONLY:DTW,TIMCO
      USE schism_msgp, ONLY: myrank

      IMPLICIT NONE
      CHARACTER(len=29) :: CHR
      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/
!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'SET_WAVE_BOUNDARY_CONDITION')


      IF ( TIMCO*SEC2DAY > WBTMJD-1.E-8 .AND. &
           TIMCO*SEC2DAY < WBEMJD ) THEN 

         if(myrank.eq.0) print*,'Read next time step from boundary file'
!        Read next time step from boundary file ...
!        IF (LBINTER) THEN
          CALL WAVE_BOUNDARY_CONDITION !(WBACNEW)
          WBACNEW = WBAC
          DSPEC   = (WBACNEW-WBACOLD)/WBDELT*DTW
          WBAC    =  WBACOLD
          WBACOLD =  WBACNEW
!        ELSE ! .NOT. LBINTER
!          CALL WAVE_BOUNDARY_CONDITION(WBAC)
!        END IF ! LBINTER
          WBTMJD = WBTMJD + WBDELT*SEC2DAY
      ELSE ! Interpolate in time ... no need to read ...
         if(myrank.eq.0) print*,'Interpolate in time'
!        IF (LBINTER) THEN
          WBAC = WBAC + DSPEC
!        END IF
      END IF

!      CALL SET_WAVE_BOUNDARY

!      IF (LNANINFCHK) THEN
!        WRITE(DBG%FHNDL,*) ' FINISHED WITH BOUNDARY CONDITION ',SUM(AC2)
!        IF (SUM(AC2) .NE. SUM(AC2)) &
!     &      CALL WWM_ABORT('NAN IN BOUNDARY CONDTITION l.1959')
!      ENDIF
      END SUBROUTINE SET_WAVE_BOUNDARY_CONDITION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_BOUNDARY_CONDITION

      USE OCPCOMM4
      USE schism_msgp, ONLY : parallel_abort,myrank
      USE VARS_WAVE

      IMPLICIT NONE
      LOGICAL       :: DoAllocate
      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/
!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'INIT_WAVE_BOUNDARY_CONDITION')


!      WRITE(STAT%FHNDL, *) 'Begin of INIT_WAVE_BOUNDARY_CONDITION'
!TODO: Makes sure initial condition work also 
!when no wave boundary is set ...
!     IF (IBOUNDFORMAT == 2) THEN

!        IF( TRIM(NESTING_TYPE_WAVE) .EQ.'SPEC') THEN ! Spectrum is prescribed
!          CALL INIT_BINARY_WW3_SPECTRA
!        END IF

      IF (TRIM(NESTING_TYPE_WAVE) .EQ.'PARAM') THEN ! Parametric Wave Boundary is prescribed (NetCDF)

          CALL INIT_NETCDF_WW3_WAVEPARAMETER

      ELSEIF (TRIM(NESTING_TYPE_WAVE) .EQ.'WAMPARAM') THEN ! MultiModal Parametric Wave Boundary ala WAM is prescribed (NetCDF)
                                                          ! Netcdf. Format Vector with dim=IOB_NODE , not X,Y matrices
          CALL INIT_NETCDF_WAM_WAVEPARAMETER

      ELSEIF (TRIM(NESTING_TYPE_WAVE) .EQ.'SPEC') THEN ! Spectrum is prescribed (NetCDF)

          CALL INIT_NC_WW3_SPECTRA  

      ELSE

          CALL parallel_abort("ERROR IN INIT_WAVE_BOUNDARY_CONDITION")

      END IF

! Populate WBAC and WBACOLD  with first record

      IF( TRIM(NESTING_TYPE_WAVE) /= 'PARAM' .or.     &
          TRIM(NESTING_TYPE_WAVE) /= 'WAMPARAM' .or.  &
          TRIM(NESTING_TYPE_WAVE) /= 'SPEC') THEN

         LBINTER = .FALSE.

         CALL WAVE_BOUNDARY_CONDITION !(WBAC)
         !IF (LBINTER) WBACOLD = WBAC
         WBACOLD = WBAC

         LBINTER = .TRUE.
        
      END IF

      END SUBROUTINE INIT_WAVE_BOUNDARY_CONDITION


!**********************************************************************
!* INIT_NETCDF_WAM_WAVEPARAMETER                                      *
!*
!* Called by INIT_WAVE_BOUNDARY_CONDITION
!*
!* Read Multi-Modal wave paremeters from a WAM netcdf file like:
!* dimensions:
!*         time = UNLIMITED ; // (82 currently)
!*         nOpenBndNodes = 235 ;
!*         wave_partition = 3 ;
!* variables:
!*         double time(time) ;
!*                 time:base_date = 2019, 9, 21, 0 ;
!*                 time:units = "days since 2019-09-21" ;
!*                 time:standard_name = "time" ;
!*                 time:long_name = "time" ;
!*         float lat(nOpenBndNodes) ;
!*                 lat:long_name = "nodal latitude" ;
!*                 lat:units = "degrees_north" ;
!*         float lon(nOpenBndNodes) ;
!*                 lon:long_name = "nodal longitude" ;
!*                 lon:units = "degrees_east" ;
!*         float h(nOpenBndNodes) ;
!*                 h:long_name = "depth" ;
!*                 h:units = "m" ;
!*         int obc_nodes(nOpenBndNodes) ;
!*                 obc_nodes:long_name = "node node node\000\000long_na" ;
!*         float per_split(time, nOpenBndNodes, wave_partition) ;
!*                 per_split:long_name = "wave period in wave partitions" ;
!*                 per_split:units = "s" ;
!*         float dir_split(time, nOpenBndNodes, wave_partition) ;
!*                 dir_split:long_name = "wave direction in wave partitions" ;
!*                 dir_split:units = "degree" ;
!*         float hs_split(time, nOpenBndNodes, wave_partition) ;
!*                 hs_split:long_name = "wave significant height in wave partitions" ;
!*                 hs_split:units = "m" ;
!*         float dspr_split(time, nOpenBndNodes, wave_partition) ;
!*                 dspr_split:long_name = "one-sided directional spreading in wave partitions" ;
!*                 dspr_split:units = "degree" ;
!* 
!* // global attributes:
!*                 :bulk_spec_mode = "MEAN JONSWAP MODE" ;
!*                 :FSHAPE_WAM = 2 ;      # Wave spectra function :  FSHAPE_WAM   1:PM,2:JON,3:BIN,4:GAUS
!*                 :CHAR_WAM_PERIOD = 2 ; # 1:PEAK or 2:MEAN frequency
!*                 :CHAR_WAM_DSPR = 1 ;   # Directional distribution DEGREES or POWER ?   1:DEGR, 2:POWER
!*                 :data_source = "Wave parameters from MFWAM" ;
!*
!* Authors: Jerome Lefevre (05/2021)
!**********************************************************************
      SUBROUTINE INIT_NETCDF_WAM_WAVEPARAMETER

      USE TIMECOMM
      USE OCPCOMM4, ONLY: PRINTF
      USE VARS_WAVE, ONLY: FSHAPE_WAM, CHAR_WAM_DSPR,CHAR_WAM_PERIOD 
      USE VARS_WAVE, ONLY: NESTING_FILE,N_WAVE_WAM,N_OBN_WAM
      USE VARS_WAVE, ONLY: HS_WAM,DIR_WAM,PER_WAM,DSPR_WAM
      USE VARS_WAVE, ONLY: NDT_BND_FILE,BND_TIME_ALL_FILES,WBDELT
      USE VARS_WAVE, ONLY: SEC2DAY,DAY2SEC,WBBMJD,WBEMJD,SERIAL
      USE VARS_WAVE, ONLY: NUM_NETCDF_FILES_BND
      USE VARS_WAVE, ONLY : IOBCN_GL_W,I_OBC_GL_W
!      USE schism_glbl, ONLY: iplg,rkind
!      USE schism_glbl, ONLY: IOBCN_W,I_OBC_N_W
      USE schism_glbl, ONLY: iplg,rkind
      USE schism_msgp
      USE NETCDF

      IMPLICIT NONE

      include 'mpif.h'
      include 'netcdf.inc'

      INTEGER :: BND_NCID,ISTAT
      INTEGER :: OBN_DIM_ID, TIME_DIM_ID, WPART_DIM_ID
      INTEGER :: TIME_VAR_ID, OBN_VAR_ID
      INTEGER :: MAXSTEP_WAM
      character :: attr_name*50
      integer, allocatable, dimension(:) :: base_date
      integer, allocatable, dimension(:) :: I_OBN_WAM
      integer :: day, month, year, n_base_date, allocate_stat
      INTEGER :: I, JD_ORIGIN
      INTEGER :: iProc, eSize
      INTEGER :: iVect(7)
      LOGICAL :: lexist

      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/
!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'INIT_NETCDF_WAM_WAVEPARAMETER')

      NUM_NETCDF_FILES_BND = 1 ! Hard coded, always 1

      IF (myrank .eq. 0) THEN

        inquire(file=TRIM(NESTING_FILE),exist=lexist)
        IF(.NOT.lexist) &
           CALL parallel_abort("Missing WAM boundary file")

        ISTAT = NF90_OPEN(TRIM(NESTING_FILE), &
     &          NF90_NOWRITE, BND_NCID)
        IF (ISTAT /= 0) THEN
          Print *, 'Error while trying to open netcdf file'
          Print *, 'FILE=', TRIM(NESTING_FILE)
          Print *, 'One possible error is that the file is NC4 but'
          Print *, 'you linked to NC3'
          CALL parallel_abort("Error while trying to open netcdf file")
        END IF

        ISTAT = NF90_INQ_VARID(BND_NCID, 'time', TIME_VAR_ID)
        if(ISTAT.ne.NF90_NOERR) &
          call parallel_abort("Netcdf Inquire varid error")

        ! Reading time dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'time',TIME_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING TIME DIM ID")

        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,TIME_DIM_ID, &
     &          len = MAXSTEP_WAM)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING TIME DIM")

        ! Reading OBC node dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'nOpenBndNodes',OBN_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING OBC nodes DIM ID")

        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,OBN_DIM_ID, &
     &          len = N_OBN_WAM)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING OBC nodes DIM")

        ! OBC node check out ... 
        IF (N_OBN_WAM /= IOBCN_GL_W) &
        CALL parallel_abort("ERROR nOpenBndNodes /= IOBCN_GL_W")

        ! Reading Wave Partition dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'wave_partition',WPART_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING Wave Partition DIM ID")

        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,WPART_DIM_ID, &
     &          len = N_WAVE_WAM)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING Wave Partition DIM")

        ! Reading settings for Bulk Wave Spectrum Reconstruction
        !
        ! Wave spectra function :  FSHAPE_WAM   1:PM,2:JON,3:BIN,4:GAUS
        attr_name = 'FSHAPE_WAM'
        ISTAT = NF_GET_ATT_INT1(BND_NCID,NF90_GLOBAL,attr_name,& 
                FSHAPE_WAM)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR: GLB. ATTR. FSHAPE_WAM not found")
        IF((ISTAT /= 0).or.(FSHAPE_WAM.lt.1).or.(FSHAPE_WAM.gt.4)) &
        CALL parallel_abort("ERROR: ATTR. FSHAPE_WAM value error")

        ! Period: PEAK or MEAN frequency ?     1:PEAK or 2:MEAN frequency
        attr_name = 'CHAR_WAM_PERIOD'
        ISTAT = NF_GET_ATT_INT1(BND_NCID,NF90_GLOBAL,attr_name,&
                CHAR_WAM_PERIOD)
        IF (ISTAT /= 0) CALL parallel_abort(&
             "ERROR: GLB. ATTR. CHAR_WAM_PERIOD not found")
        IF((CHAR_WAM_PERIOD.lt.1).or.(CHAR_WAM_PERIOD.gt.2)) &
        CALL parallel_abort("ERROR: ATTR. CHAR_WAM_PERIOD value error")

        ! Directional distribution DEGREES or POWER ?   1:DEGR, 2:POWER
        attr_name = 'CHAR_WAM_DSPR'
        ISTAT = NF_GET_ATT_INT1(BND_NCID,NF90_GLOBAL,attr_name,&
                CHAR_WAM_DSPR)
        IF (ISTAT /= 0) CALL parallel_abort(&
             "ERROR: GLB. ATTR. CHAR_WAM_DSPR not found")
         IF((CHAR_WAM_DSPR.lt.1).or.(CHAR_WAM_DSPR.gt.2)) &
        CALL parallel_abort("ERROR: ATTR. CHAR_WAM_DSPR value error")


        ! Wave spectra function :  FSHAPE_WAM   1:PM,2:JON,3:BIN,4:GAUS
!        ISTAT = NF90_GET_ATT(BND_NCID, NF90_GLOBAL,'FSHAPE_WAM', &
!     &          FSHAPE_WAM)
!        IF((ISTAT /= 0).or.(FSHAPE_WAM.lt.1).or.(FSHAPE_WAM.gt.4)) &
!        CALL parallel_abort("ERROR: ATTR. FSHAPE_WAM value error")

        ! Period: PEAK or MEAN frequency ?     1:PEAK or 2:MEAN frequency
!        ISTAT = NF90_GET_ATT(BND_NCID, NF90_GLOBAL,'CHAR_WAM_PERIOD', &
!     &          CHAR_WAM_PERIOD)
!        IF((ISTAT /= 0).or.(CHAR_WAM_PERIOD.lt.1).or.(CHAR_WAM_PERIOD.gt.2)) &
!        CALL parallel_abort("ERROR: ATTR. CHAR_WAM_PERIOD value error") 

        ! Directional distribution DEGREES or POWER ?   1:DEGR, 2:POWER 
!        ISTAT = NF90_GET_ATT(BND_NCID, NF90_GLOBAL,'CHAR_WAM_DSPR', &
!     &          CHAR_WAM_DSPR)
!        IF((ISTAT /= 0).or.(CHAR_WAM_DSPR.lt.1).or.(CHAR_WAM_DSPR.gt.2)) &
!        CALL parallel_abort("ERROR: ATTR. CHAR_WAM_PERIOD value error")
 
        ! output some details
        WRITE(PRINTF,*) 'WAM NESTING FILE:',TRIM(NESTING_FILE)
        WRITE(PRINTF,*) '   nb records         = ',MAXSTEP_WAM
        WRITE(PRINTF,*) '   nb OBN nodes       = ',N_OBN_WAM
        WRITE(PRINTF,*) '   nb Wave Partitions = ',N_WAVE_WAM
        IF(FSHAPE_WAM.EQ.1) &
        WRITE(PRINTF,*) '   bulk spectrum = Pierson-Moskowitz'
        IF(FSHAPE_WAM.EQ.2) &
        WRITE(PRINTF,*) '   bulk spectrum = JONSWAP'
        IF(FSHAPE_WAM.EQ.3) &
        WRITE(PRINTF,*) '   bulk spectrum = ONE Frequency bin'
        IF(FSHAPE_WAM.EQ.4) &
        WRITE(PRINTF,*) '   bulk spectrum = Gaussian-shaped frequency'
        IF(CHAR_WAM_PERIOD.EQ.1) &
        WRITE(PRINTF,*) '   characteristic wave period = PEAK'
        IF(CHAR_WAM_PERIOD.EQ.2) &
        WRITE(PRINTF,*) '   characteristic wave period = MEAN'
        IF(CHAR_WAM_DSPR.EQ.1) &
        WRITE(PRINTF,*) '   directional width = DEGREES'
        IF(CHAR_WAM_DSPR.EQ.2) &
        WRITE(PRINTF,*) '   directional width = POWER'

        ! Read Times 
        ! ----------------------
        WRITE(PRINTF, *) 'INIT_NETCDF_WAM_WAVEPARAMETER ALLOCATE BEGIN'

        ALLOCATE(NDT_BND_FILE(1), stat=istat)
        IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.1")
        NDT_BND_FILE(1) = MAXSTEP_WAM

        ALLOCATE (BND_TIME_ALL_FILES(1,MAXSTEP_WAM), stat=istat)
        IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.2")
        BND_TIME_ALL_FILES = 0._rkind
        
        ! read all time steps in the proper format and transform in wwm time line
        !
        ! in WAM files, the time origine in time attributes is like:
        !        time:base_date = 2019, 9, 21, 0
        !
        ! get the base_date - the time and date which corresponds with time zero
        !       Note that only year,month,day of base_date are used (not 'hour').
        attr_name = 'base_date'

        ! get size of base_date (make sure it is at least 3)
        ISTAT =NF_INQ_ATTLEN(BND_NCID,TIME_DIM_ID,attr_name,n_base_date)
        if (n_base_date .lt. 3) &
          call parallel_abort('insufficient fields in base_date!')

        ! allocate space for base_date
        allocate(base_date(n_base_date), stat=istat)
        IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.6")
        base_date = 0

        ! read base_date
        ISTAT =NF_GET_ATT_INT(BND_NCID,TIME_DIM_ID,attr_name,base_date)

        ! convert base_date to integer Julian date
        year = base_date(1);  month = base_date(2);  day = base_date(3)

        JD_ORIGIN = jd(year,month,day)

        ! deallocate base_date
        deallocate(base_date)

        ! Boundary nodes, check out ...
        allocate(I_OBN_WAM(N_OBN_WAM))

        ISTAT = NF90_INQ_VARID(BND_NCID, 'obc_nodes',OBN_VAR_ID)
        if(ISTAT.ne.NF90_NOERR) &
        call parallel_abort("Netcdf Inquire obc_nodes error")

        ISTAT = NF90_GET_VAR(BND_NCID,OBN_VAR_ID,I_OBN_WAM,&
     &          start=(/1/),count = (/N_OBN_WAM/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING obc_nodes")

        DO I=1,N_OBN_WAM
           IF( I_OBN_WAM(I) /= I_OBC_GL_W(I) ) &
           CALL parallel_abort("ERROR in Boundary nodes index") 
        ENDDO

        deallocate(I_OBN_WAM)

        ! Reading time
!       ISTAT = NF90_INQ_VARID(BND_NCID,'time',TIME_VAR_ID)
!       IF (ISTAT /= 0) &
!       CALL parallel_abort("ERROR WHILE READING TIME VARIABLE")
        ISTAT = NF90_GET_VAR(BND_NCID,TIME_VAR_ID,BND_TIME_ALL_FILES,&
     &          start=(/1/),count = (/MAXSTEP_WAM/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING TIME")

       ! Time origine is : "days since yyyy, mm, dd    00:00:00"
       ! so add that date in jd
        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + JD_ORIGIN 

        ! Closing the netcdf file
        ISTAT = NF90_CLOSE(BND_NCID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE CLOSING NETCDF FILE")

      ENDIF ! myrank

      IF(.NOT.SERIAL) THEN
        IF (myrank.eq.0) THEN
          iVect(1)=MAXSTEP_WAM
          iVect(2)=N_OBN_WAM
          iVect(3)=N_WAVE_WAM
          iVect(4)=FSHAPE_WAM
          iVect(5)=CHAR_WAM_PERIOD
          iVect(6)=CHAR_WAM_DSPR
          DO IPROC=2,nproc
            CALL MPI_SEND(iVect,6,itype, iProc-1, 811, comm, ierr)
          END DO
          eSize=MAXSTEP_WAM
          DO IPROC=2,nproc
            CALL MPI_SEND(BND_TIME_ALL_FILES,eSize,rtype, iProc-1, 813, &
     &           comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(iVect,6,itype, 0, 811, comm, istatus, ierr)
          MAXSTEP_WAM= iVect(1)
          N_OBN_WAM  = iVect(2)
          N_WAVE_WAM = iVect(3)
          FSHAPE_WAM = iVect(4)
          CHAR_WAM_PERIOD= iVect(5)
          CHAR_WAM_DSPR  = iVect(6)
          ALLOCATE(BND_TIME_ALL_FILES(1,MAXSTEP_WAM), &
     &             NDT_BND_FILE(1), stat=istat)
          IF (istat/=0)   &
             CALL parallel_abort('WAM Nesting, allocate err. 7')
          NDT_BND_FILE(1) = MAXSTEP_WAM

          eSize=MAXSTEP_WAM
          CALL MPI_RECV(BND_TIME_ALL_FILES,eSize,rtype, 0, 813, &
     &         comm, istatus, ierr)
        END IF
      END IF


      WBDELT = (BND_TIME_ALL_FILES(1,2)-BND_TIME_ALL_FILES(1,1))*DAY2SEC
      IF(myrank.EQ.0) &
        WRITE(PRINTF,*) 'WAM Frequency (sec) =',REAL(WBDELT)
      ! Start Time in mjd (days)
      ! End Time in mjd (days)
      WBBMJD = BND_TIME_ALL_FILES(1,1)
      WBEMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WAM)   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Ultimate Time bound checking
      IF(myrank.EQ.0) THEN
       WRITE(PRINTF,*) 'SWAN Start (day):',TINIC*SEC2DAY
       WRITE(PRINTF,*) 'SWAN End   (day):',TFINC*SEC2DAY
       WRITE(PRINTF,*) 'First WAM record (day):',WBBMJD
       WRITE(PRINTF,*) 'Last  WAM record (day):',WBEMJD
      ENDIF

      IF( (WBBMJD.GT.TINIC*SEC2DAY) .OR.          &
          (WBEMJD.LT.TFINC*SEC2DAY)) THEN
         CALL parallel_abort( &
         'no appropriate time exists in WAM BC files')
      ENDIF

      ! Allocate Buckets 
      ALLOCATE (HS_WAM(N_WAVE_WAM,N_OBN_WAM),stat=istat)
      IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.3")
      ALLOCATE (PER_WAM(N_WAVE_WAM,N_OBN_WAM),stat=istat)
      IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.4")
      ALLOCATE (DIR_WAM(N_WAVE_WAM,N_OBN_WAM),stat=istat)
      IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.5")
      ALLOCATE (DSPR_WAM(N_WAVE_WAM,N_OBN_WAM),stat=istat)
      IF (istat/=0) CALL parallel_abort("WAM Nesting, allocate err.6")
      HS_WAM  = 0._rkind
      PER_WAM = 0._rkind
      DSPR_WAM= 0._rkind
      DIR_WAM = 0._rkind

      END SUBROUTINE

!**********************************************************************
!* GET_NETCDF_WAM_PARAM
!* Read a WAM multimode wave parameters netcdf file and do time and space
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITION
!*
!* Authors: Jerome Lefevre 05/2021
!**********************************************************************
      SUBROUTINE GET_NETCDF_WAM_PARAM 

      USE VARS_WAVE, ONLY : NDT_BND_FILE,IOBCN_W,I_OBC_N_W
!      USE SWCOMM3, ONLY : MDC,MSC,SPPARM,SLOW,SHIG,PI2_W,PI_W
      USE schism_glbl, ONLY : rkind
      USE schism_msgp, ONLY : parallel_abort,myrank

      INTEGER     :: IT
      INTEGER     :: IFILE ! not used 

      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

      WRITESTATFLAG = 1 ! Print DBG info
      IT =0
      IFILE=0

!     STRACE           Tracing routine for debugging
      CALL STRACE (IENT,'GET_NETCDF_WAM_PARAM')

      CALL COMPUTE_IFILE_IT(IFILE, IT)

      ! WAM: take the next record 
      IT = MIN(IT+1,NDT_BND_FILE(1))

      !CALL COMPUTE_IT(IT)
!      print*,'IT=',IT

      ! Read Wave parameters in file
      CALL READ_NETCDF_WAM_PARAM( IT )

!          th_time2(1,2)=th_time2(2,2)
!          th_time2(2,2)=th_time2(2,2)+th_dt2(2)
!        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for flux 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

!        rat=(time-th_time2(1,2))/th_dt2(2)
!        if(rat<-small1.or.rat>1+small1) then
!          write(errmsg,*) 'STEP: rat out in uv3D.th:',rat,time,th_time2(1:2,2)
!          call parallel_abort(errmsg)
!        endif
!        icount=0
!        icount2=0
!        do k=1,nope_global
!          if(iabs(ifltype(k))>=4) then
!            icount=icount+1
!            if(icount>nfltype2) call parallel_abort('Wrong counting 6')
!            do j=1,nond_global(k)
!              icount2=icount2+1
!              if(icount2>nnode_fl) call parallel_abort('Wrong counting vel')
!'
!              uthnd(1:nvrt,j,k)=(1-rat)*ath2(1,1:nvrt,icount2,1,2)+rat*ath2(1,1:nvrt,icount2,2,2) !ll frame if ics=2
!              vthnd(1:nvrt,j,k)=(1-rat)*ath2(2,1:nvrt,icount2,1,2)+rat*ath2(2,1:nvrt,icount2,2,2)
!            enddo !j
!          endif
!        enddo !k
!      endif !nfltype2

      END SUBROUTINE

!**********************************************************************
!* READ_NETCDF_WAM_PARAM
!* Authors: Jerome Lefevre 05/2021
!**********************************************************************

      SUBROUTINE READ_NETCDF_WAM_PARAM( IT )

      USE VARS_WAVE, ONLY: PER_WAM,DIR_WAM,HS_WAM,DSPR_WAM
      USE VARS_WAVE, ONLY: N_WAVE_WAM,N_OBN_WAM,NESTING_FILE
!      USE VARS_WAVE, ONLY: 
      USE schism_glbl, ONLY: skind,rkind
      USE schism_msgp, ONLY: parallel_abort,myrank,comm
      USE NETCDF

      IMPLICIT NONE

      include 'mpif.h'

      integer, intent(in) :: IT
!     CHARACTER(LEN=40)  :: EVAR
!     integer :: ITMP(NDX_BND,NDY_BND)
      REAL(rkind) :: scale_factor
      REAL(skind),allocatable :: RTMP(:,:)
!      REAL(rkind),allocatable :: buffer(:,:,:)
      REAL(4),allocatable :: buffer(:,:,:)   ! be coherent with 'mpi_real'
      integer :: BND_NCID, var_id, ISTAT
      LOGICAL :: lexist

      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

!     STRACE           Tracing routine for debugging
      CALL STRACE (IENT,'READ_NETCDF_WAM_PARAM')


!      if(time>th_time2(2,2)) then
!      ath2(:,:,:,1,2)=ath2(:,:,:,2,2)

!     PER_WAM(:,:,1)  = PER_WAM(:,:,2)
!     DIR_WAM(:,:,1)  = DIR_WAM(:,:,2)
!     HS_WAM(:,:,1)   = HS_WAM(:,:,2)
!     DSPR_WAM(:,:,1) = DSPR_WAM(:,:,2)

      allocate(buffer(4,N_WAVE_WAM, N_OBN_WAM),stat=istat)
      if(istat/= 0) call parallel_abort('buffer error')

      buffer = 0._rkind

      IF (myrank==0) THEN

         allocate(RTMP(N_WAVE_WAM, N_OBN_WAM),stat=istat)

         ISTAT = NF90_OPEN(TRIM(NESTING_FILE), &
     &          NF90_NOWRITE, BND_NCID)
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Open error")

         ! Read  'wave period in wave partitions'
         ISTAT = NF90_INQ_VARID(BND_NCID, 'per_split', var_id )
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

         scale_factor = 1._rkind

         ISTAT = NF90_GET_VAR(BND_NCID, var_id, RTMP ,  &
         start = (/ 1, 1, IT /), count = (/ N_WAVE_WAM, N_OBN_WAM, 1/))
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Get per_split error")
         buffer(1,1:N_WAVE_WAM,1:N_OBN_WAM) = MyREAL(RTMP)*scale_factor

         ! Read  'wave direction in wave partitions'
         ISTAT = NF90_INQ_VARID(BND_NCID, 'dir_split', var_id )
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

         scale_factor = 1._rkind

         ISTAT = NF90_GET_VAR(BND_NCID, var_id, RTMP ,  &
         start = (/ 1, 1, IT /), count = (/ N_WAVE_WAM, N_OBN_WAM, 1/))
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Get dir_split error")
         buffer(2,1:N_WAVE_WAM,1:N_OBN_WAM) = MyREAL(RTMP)*scale_factor

         ! Read  'wave significant height in wave partitions'
         ISTAT = NF90_INQ_VARID(BND_NCID, 'hs_split', var_id )
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

         scale_factor = 1._rkind

         ISTAT = NF90_GET_VAR(BND_NCID, var_id, RTMP ,  &
         start = (/ 1, 1, IT /), count = (/ N_WAVE_WAM, N_OBN_WAM, 1/))
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Get hs_split error")
         buffer(3,1:N_WAVE_WAM,1:N_OBN_WAM) = MyREAL(RTMP)*scale_factor

         ! Read  'one-sided directional spreading in wave partitions'
         ISTAT = NF90_INQ_VARID(BND_NCID, 'dspr_split', var_id )
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

         scale_factor = 1._rkind

         ISTAT = NF90_GET_VAR(BND_NCID, var_id, RTMP ,  &
         start = (/ 1, 1, IT /), count = (/ N_WAVE_WAM, N_OBN_WAM, 1/))
         if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Get dspr_split error")
         buffer(4,1:N_WAVE_WAM,1:N_OBN_WAM) = MyREAL(RTMP)*scale_factor

        ! Closing the netcdf file
        ISTAT = NF90_CLOSE(BND_NCID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE CLOSING NETCDF FILE")

        deallocate(RTMP)

      ENDIF !myrank

      CALL MPI_BCAST(buffer,4*N_WAVE_WAM*N_OBN_WAM,mpi_real,0,comm,istat)
      
      !PER_WAM(1:N_WAVE_WAM,1:N_OBN_WAM,2)  = buffer(1,1:N_WAVE_WAM,1:N_OBN_WAM)
      !DIR_WAM(1:N_WAVE_WAM,1:N_OBN_WAM,2)  = buffer(2,1:N_WAVE_WAM,1:N_OBN_WAM)
      !HS_WAM(1:N_WAVE_WAM,1:N_OBN_WAM,2)   = buffer(3,1:N_WAVE_WAM,1:N_OBN_WAM)
      !DSPR_WAM(1:N_WAVE_WAM,1:N_OBN_WAM,2) = buffer(4,1:N_WAVE_WAM,1:N_OBN_WAM)
      PER_WAM(1:N_WAVE_WAM,1:N_OBN_WAM)  = buffer(1,1:N_WAVE_WAM,1:N_OBN_WAM)
      DIR_WAM(1:N_WAVE_WAM,1:N_OBN_WAM)  = buffer(2,1:N_WAVE_WAM,1:N_OBN_WAM)
      HS_WAM(1:N_WAVE_WAM,1:N_OBN_WAM)   = buffer(3,1:N_WAVE_WAM,1:N_OBN_WAM)
      DSPR_WAM(1:N_WAVE_WAM,1:N_OBN_WAM) = buffer(4,1:N_WAVE_WAM,1:N_OBN_WAM)

      deallocate(buffer)

      END SUBROUTINE

!**********************************************************************
!* SPPARM_WAM
!* Authors: Jerome Lefevre 05/2021
!**********************************************************************
      SUBROUTINE SPPARM_WAM( K_WAV, VAL )

      USE VARS_WAVE, ONLY:IOBCN_W,IOBCN_GL_W,I_OBC_GL_W
      USE VARS_WAVE, ONLY:HS_WAM,PER_WAM,DIR_WAM,DSPR_WAM
      USE VARS_WAVE, ONLY:CHAR_WAM_PERIOD,CHAR_WAM_DSPR
      USE schism_glbl, ONLY:iplg,ipgl,ISONB_W,rkind
      USE schism_msgp, ONLY:parallel_abort,myrank

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: K_WAV  ! Wave partition number
      REAL(rkind), INTENT(OUT) :: VAL(8,IOBCN_W)
      INTEGER :: IP,ID
      INTEGER :: IVAR, I, J, icount

      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/

!     STRACE           Tracing routine for debugging
      CALL STRACE (IENT,'SPPARM_WAM')

!             print('(A,4I7)'),'rank,BND node,ISBND,ISONB_W=', &
!                    myrank,iplg(I_OBC_N_W(I)), &
!                    isbnd(1,I_OBC_N_W(I)),isonb_w(I_OBC_N_W(I))
!            enddo

      icount=0

      VAL = 0.0_rkind

      DO I=1,IOBCN_GL_W ! Global 

         ID = I_OBC_GL_W(I)
         if(ipgl(ID)%rank/=myrank)  cycle 

         IP = ipgl(ID)%id
         IF (ISONB_W(IP)/=2) & !cycle ! should never append...
         CALL parallel_abort("ERROR not on Wave Boundary?",ID)

         icount=icount+1
         ! eth(j,k)=(1-rat)*ath2(1,1,icount2,1,1)+rat*ath2(1,1,icount2,2,1)

         VAL(1,icount) = HS_WAM(K_WAV,I) 
         VAL(2,icount) = PER_WAM(K_WAV,I)
         VAL(3,icount) = DIR_WAM(K_WAV,I)
         VAL(4,icount) = DSPR_WAM(K_WAV,I)
         VAL(6,icount) = 1.
         VAL(7,icount) = 0.1
         VAL(8,icount) = 3.3

         ! FSHAPE_WAM: Wave spectra function :       1:PM,2:JON,3:BIN,4:GAUS
         ! CHAR_WAM_PERIOD: PEAK or MEAN frequency   1:PEAK or 2:MEAN frequency
         ! CHAR_WAM_DSPR Directional distribution DEGREES or POWER ?   1:DEGR, 2:POWER
         IF (CHAR_WAM_PERIOD==1) THEN ! PeaK 
            VAL(5,icount) = 2._rkind
         ELSE 
            VAL(5,icount) = -2._rkind
         ENDIF

         IF (CHAR_WAM_DSPR==1) THEN !Directional distribution DEGREES
            VAL(6,icount) = 1.
         ELSE
            VAL(6,icount) = 2.
         ENDIF
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
      ENDDO
      END SUBROUTINE

!**********************************************************************
!* INIT_NC_WW3_SPECTRA
!* Read the header of a WAVEWATCHIII netcdf spectral file to get
!* dimensions of the spectral grid and number of output locations
!* required for dynamic allocation.
!*
!* Called by INIT_WAVE_BOUNDARY_CONDITION
!*
!* Authors: Kevin Martins (07/2018)
!**********************************************************************
      SUBROUTINE INIT_NC_WW3_SPECTRA

      USE VARS_WAVE
      USE NETCDF
      USE schism_msgp, ONLY : parallel_abort,myrank
      USE TIMECOMM
      USE OCPCOMM4, ONLY: PRINTF
      USE SWCOMM3, ONLY : PI2_W,PI_W

      IMPLICIT NONE

      INTEGER :: BND_NCID,ISTAT
      INTEGER :: FREQ_DIM_ID, DIRS_DIM_ID, TIME_DIM_ID, STAT_DIM_ID
      INTEGER :: FREQ_VAR_ID, DIRS_VAR_ID, TIME_VAR_ID, STAT_LAT_ID, STAT_LON_ID
      REAL(rkind) :: eJD_WW3, eJD_WWM
      LOGICAL :: lexist
      INTEGER IENT
      SAVE  IENT
      DATA  IENT /0/
!     STRACE           Tracing routine for debugging
      CALL  STRACE (IENT,'INIT_NC_WW3_SPECTRA')

      !------------------------------------
      ! Reading grid and station dimensions
      !------------------------------------
      ! Opening netcdf file

!      IF (myrank .eq. 0) THEN

        inquire(file=TRIM(NESTING_FILE),exist=lexist)
        IF(.NOT.lexist) &
           CALL parallel_abort("Missing WW3 boundary file")

        ISTAT = NF90_OPEN(NESTING_FILE,NF90_NOWRITE,BND_NCID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("Error while trying to open netcdf file")

        ! Reading frequency dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'frequency',FREQ_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING FREQUENCY DIM ID")

        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,FREQ_DIM_ID,len = MSC_WW3)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING FREQUENCY DIM")

        ! Reading direction dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'direction',DIRS_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTION DIM ID")
        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,DIRS_DIM_ID,len = MDC_WW3)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTION DIM")

        ! Reading time dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'time',TIME_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTION DIM ID")
        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,TIME_DIM_ID,len = MAXSTEP_WW3)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTION DIM")

        ! Reading station dimension
        ISTAT = NF90_INQ_DIMID(BND_NCID,'station',STAT_DIM_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING STATION DIM ID")
 
        ISTAT = NF90_INQUIRE_DIMENSION(BND_NCID,STAT_DIM_ID,len = NP_WW3)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING STATION DIMENSION")

        !------------------------------------
        ! Reading variables
        !------------------------------------
        ! Arrays for the WW3 spectra
        print*,MSC_WW3,MDC_WW3,MAXSTEP_WW3,NP_WW3

        ALLOCATE(FQ_WW3(MSC_WW3),DR_WW3(MDC_WW3),         &
     &   BND_TIME_ALL_FILES(1,MAXSTEP_WW3),XP_WW3(NP_WW3),&
     &   YP_WW3(NP_WW3),stat=istat)
        IF (istat/=0) CALL parallel_abort("INIT_NC_WW3_SPECTRA, &
     &       WW3 freq/dirs allocate error")

        ! Reading longitude
        ISTAT = NF90_INQ_VARID(BND_NCID,'longitude',STAT_LON_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING LONGITUDE VARIABLE")

        ISTAT = NF90_GET_VAR(BND_NCID,STAT_LON_ID,XP_WW3,start=(/1,1/),count = (/NP_WW3,1/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING LONGITUDE")

        ! Reading latitude
        ISTAT = NF90_INQ_VARID(BND_NCID,'latitude',STAT_LAT_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING LATITUDE VARIABLE")

        ISTAT = NF90_GET_VAR(BND_NCID,STAT_LAT_ID,YP_WW3,start=(/1,1/),count = (/NP_WW3,1/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING LATITUDE")

        ! Convert coordinates to RAD
        !XP_WW3 = XP_WW3*DEGRAD
        !YP_WW3 = YP_WW3*DEGRAD

        ! Reading frequencies
        ISTAT = NF90_INQ_VARID(BND_NCID,'frequency',FREQ_VAR_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING FREQUENCY VARIABLE")
        ISTAT = NF90_GET_VAR(BND_NCID,FREQ_VAR_ID,FQ_WW3,start=(/1/),count = (/MSC_WW3/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING FREQUENCIES")

        ! Reading directions
        ISTAT = NF90_INQ_VARID(BND_NCID,'direction',DIRS_VAR_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTION VARIABLE")
        ISTAT = NF90_GET_VAR(BND_NCID,DIRS_VAR_ID,DR_WW3,start=(/1/),count = (/MDC_WW3/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING DIRECTIONS")

        ! Reading time
        ISTAT = NF90_INQ_VARID(BND_NCID,'time',TIME_VAR_ID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING TIME VARIABLE")
        ISTAT = NF90_GET_VAR(BND_NCID,TIME_VAR_ID,BND_TIME_ALL_FILES,start=(/1/),count = (/MAXSTEP_WW3/))
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE READING TIME")

        ! Closing the netcdf file
        ISTAT = NF90_CLOSE(BND_NCID)
        IF (ISTAT /= 0) &
        CALL parallel_abort("ERROR WHILE CLOSING NETCDF FILE")

!      ENDIF

        !------------------------------------
        ! Conventions (direction and time)
        !------------------------------------
        ! Convention for directions: care must be taken here as sometimes the directions in WW3
        ! are expressed with the origin taken towards the East. This corresponds to the last version of WW3.


        IF(MAXVAL(DR_WW3)>300) THEN
         DR_WW3 = MOD(DR_WW3*PI_W/180,PI2_W) ! Simple conversion to RAD
        ENDIF
        ! Converting time from WW3 convention to WWM convention
        ! In WW3, time is given in Julian days, with reference from [1990,1,1,0,0,0]
        !CALL DATE2JD(1858,11,17,0,0,0,eJD_WWM)
        !CALL DATE2JD(1990,1,1,0,0,0,eJD_WW3)
        ! in WW3 fils, the time origine is :
        ! Time origine is : "days since 1990-01-01 00:00:00"
        ! Compute time_origine in jd
        eJD_WW3 = jd(1990,1,1)

       ! Time origine is : "days since 1990-01-01 00:00:00"
       ! so add that date in jd
        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + eJD_WW3 !- eJD_WWM


        WBDELT =  (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC
        IF(myrank.EQ.0) &
        WRITE(PRINTF,*) 'WW3 Frequency (sec) =',REAL(WBDELT)
        ! Start Time in mjd (days)
        ! End Time in mjd (days)
        WBBMJD = BND_TIME_ALL_FILES(1,1)
        WBEMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WW3)

      ! Ultimate Time bound checking
      IF(myrank.EQ.0) THEN
       WRITE(PRINTF,*) 'SWAN Start (day):',TINIC*SEC2DAY
       WRITE(PRINTF,*) 'SWAN End   (day):',TFINC*SEC2DAY
       WRITE(PRINTF,*) 'First WW3 record (day):',WBBMJD
       WRITE(PRINTF,*) 'Last  WW3 record (day):',WBEMJD
      ENDIF

      IF( (WBBMJD.GT.TINIC*SEC2DAY) .OR.          &
          (WBEMJD.LT.TFINC*SEC2DAY)) THEN
         CALL parallel_abort( &
         'no appropriate time exists in WW3 BC files')
      ENDIF

      END SUBROUTINE

!**********************************************************************
!*    WW3 Like Forcing, Wave parametric                               *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_WW3_WAVEPARAMETER

      USE TIMECOMM
      USE OCPCOMM4, ONLY: PRINTF
      USE VARS_WAVE
      USE schism_glbl, ONLY: rkind
      USE schism_msgp
      USE NETCDF

      IMPLICIT NONE

      include 'mpif.h'

      INTEGER :: IT, JFILE, IVAR, BND_NCID,ISTAT
      INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
      REAL(rkind), ALLOCATABLE :: BND_TIME(:)
      character (len = *), parameter :: CallFct = &
      "INIT_NETCDF_WW3_WAVEPARAMETER"
      integer, dimension(nf90_max_var_dims) :: dimIDs
      INTEGER :: MAX_NDT_BND_FILE, iProc, eSize
      INTEGER :: iVect(4)
      character (len=200) :: BUFFER
      REAL(rkind) :: rVect(4)
      LOGICAL :: lexist
      INTEGER :: JD_ORIGIN
      integer, parameter :: wav_FHNDL = 20


      IF (myrank .eq. 0) THEN

        inquire(file=TRIM(NESTING_FILE),exist=lexist)
        IF(.NOT.lexist) &
           CALL parallel_abort("Missing WW3 boundary file")
 
        OPEN(wav_FHNDL,FILE=TRIM(NESTING_FILE),STATUS='OLD')
!!        WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME

        NUM_NETCDF_FILES_BND = 0
        DO
          READ( wav_FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND + 1
        END DO
        REWIND(wav_FHNDL)

!        WRITE(PRINTF,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND

        NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND/NUM_NETCDF_VAR_TYPES

!        WRITE(PRINTF,*) 'NUM_NETCDF_FILES_BND=', NUM_NETCDF_FILES_BND
!        WRITE(PRINTF,*) 'NUM_NETCDF_VAR_TYPES=', NUM_NETCDF_VAR_TYPES

        ALLOCATE(NETCDF_FILE_NAMES_BND(NUM_NETCDF_FILES_BND,&
     &                      NUM_NETCDF_VAR_TYPES), stat=istat)

        IF (istat/=0) &
            CALL parallel_abort("swan_bdcons, allocate error")

        DO IT = 1, NUM_NETCDF_FILES_BND
          DO IVAR = 1, NUM_NETCDF_VAR_TYPES
           READ( wav_FHNDL,'(a200)')BUFFER
           NETCDF_FILE_NAMES_BND(IT,IVAR)=trim(adjustL(BUFFER))
          END DO
        END DO
        CLOSE (wav_FHNDL)
!
! four files are read to set up the wave spectra Hs, Tm01, Dir, Sprd
!
        ALLOCATE(NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
        IF (istat/=0) &
            CALL parallel_abort("swan_bdcons, allocate error")
        NDT_BND_FILE = 0
!        DO JFILE = 1, NUM_NETCDF_FILES_BND
!          WRITE(PRINTF,'(I10,10X,5A200)') &
!          JFILE, NETCDF_FILE_NAMES_BND(JFILE,:)
!        END DO
!
! check number of time steps in netcdf file ... it is assumed that all 
! files have the same ammount of time steps ...
!
        DO JFILE = 1, NUM_NETCDF_FILES_BND
!          WRITE(PRINTF,*) &
!          'ifile=', jfile, 'file=', TRIM(NETCDF_FILE_NAMES_BND(JFILE,1))

          inquire(file=TRIM(NETCDF_FILE_NAMES_BND(JFILE,1)),exist=lexist)
          IF(.NOT.lexist) &
             CALL parallel_abort("Missing WW3 boundary condition file")

          ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES_BND(JFILE,1)), &
     &            NF90_NOWRITE, BND_NCID)
          IF (ISTAT /= 0) THEN
            Print *, 'Error while trying to open netcdf file'
            Print *, 'FILE=', TRIM(NETCDF_FILE_NAMES_BND(JFILE,1))
            Print *, 'One possible error is that the file is NC4 but'
            Print *, 'you linked to NC3'
            CALL parallel_abort("Error while trying to open netcdf file")
          END IF

          ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
          if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
          if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire var error")
          
          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), &
                  len = NDT_BND_FILE(JFILE))
          if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort('Netcdf Inquire dim error')
      
          WRITE(PRINTF,*) 'nb records IFILE=',JFILE, NDT_BND_FILE(JFILE)
        END DO
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(BND_NCID, 'longitude', ILON_ID)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILON_ID, dimids=dimIDs)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len=NDX_BND)

        ISTAT = nf90_inq_varid(BND_NCID, 'latitude', ILAT_ID)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILAT_ID, dimids=dimIDs)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len=NDY_BND)

        WRITE(PRINTF,*) 'Number of Gridpoints', NDX_BND, NDY_BND

        ALLOCATE (COORD_BND_X(NDX_BND), COORD_BND_Y(NDY_BND),stat=istat)
        IF (istat/=0)   &
           CALL parallel_abort('swan_bdcons, allocate error 11')
!
! read coordinates from files ....
!
        ISTAT = NF90_GET_VAR(BND_NCID, ILON_ID, COORD_BND_X)
        ISTAT = NF90_GET_VAR(BND_NCID, ILAT_ID, COORD_BND_Y)
!
! estimate offset ...
!
        OFFSET_X_BND = MINVAL(COORD_BND_X)
        OFFSET_Y_BND = MINVAL(COORD_BND_Y)
!
! resolution ...
!
        DX_BND  = ABS(MAXVAL(COORD_BND_X) - MINVAL(COORD_BND_X))/(NDX_BND-1)
        DY_BND  = ABS(MAXVAL(COORD_BND_Y) - MINVAL(COORD_BND_Y))/(NDY_BND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(BND_NCID)
!
! total number of time steps ... in all files
!
        NDT_BND_ALL_FILES = 0
!        write(PRINTF,*) NUM_NETCDF_FILES_BND
        DO IT = 1, NUM_NETCDF_FILES_BND
          NDT_BND_ALL_FILES = NDT_BND_ALL_FILES + NDT_BND_FILE(IT)
!          write(PRINTF,*) it, NDT_BND_FILE(it)
        END DO
!        WRITE(PRINTF,*) NDT_BND_ALL_FILES, NDT_BND_FILE

        WRITE(PRINTF, *) 'Reading WW3 Files DONE'

        WRITE(PRINTF, *) 'INIT_NETCDF_WW3_WAVEPARAMETER BEGIN'
        MAX_NDT_BND_FILE = MAXVAL(NDT_BND_FILE)
        ALLOCATE (BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), stat=istat)
        IF (istat/=0)    &
           CALL parallel_abort('swan_bdcons, allocate error 12')
        BND_TIME_ALL_FILES = 0._rkind
!
! read all time steps in the proper format and transform in wwm time line
!
        ! in WW3 fils, the time origine is :
        ! Time origine is : "days since 1990-01-01 00:00:00"
        ! Compute time_origine in jd
        JD_ORIGIN = jd(1990,1,1)

        !print*,"dbg correct 565"
        !JD_ORIGIN = JD_ORIGIN +564

        DO JFILE = 1, NUM_NETCDF_FILES_BND

          inquire(file=TRIM(NETCDF_FILE_NAMES_BND(JFILE,1)),exist=lexist)
          IF(.NOT.lexist) &
             CALL parallel_abort('Missing WW3 boundary condition file')

          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(JFILE,1),NF90_NOWRITE,BND_NCID)

          ALLOCATE (BND_TIME(NDT_BND_FILE(JFILE)), stat=istat)
          IF (istat/=0)   &
             CALL parallel_abort('swan_bdcons, allocate error 13')
          BND_TIME = 0._rkind

          ISTAT = NF90_GET_VAR(BND_NCID,ITIME_ID,BND_TIME)

          DO IT = 1, NDT_BND_FILE(JFILE)
            BND_TIME_ALL_FILES(JFILE,IT) = BND_TIME(IT)
!             CALL CT2MJD('19000101.000000',DTMP1)
!             CALL CT2MJD('19900101.000000',DTMP2)
!             CALL MJD2CT(DTMP1,chrdate)
!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
!             CALL MJD2CT(DTMP2,chrdate)
!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
!             CALL MJD2CT(0.0_rkind,chrdate)
!             WRITE(*,*) '00000000.000000', 0.0_rkind, chrdate
!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
!             IF (IT == 1 .AND. JFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
!             IF (IT == 1 .AND. JFILE ==1) WRITE(*,*) JFILE, IT, BND_TIME(IT), chrdate
          END DO
          DEALLOCATE(BND_TIME)
        END DO

       ! Time origine is : "days since 1990-01-01 00:00:00"
       ! so add that date in jd
        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + JD_ORIGIN !+ DT_DIFF_19901900

      END IF

      IF(.NOT.SERIAL) THEN
      !IF (.NOT. MULTIPLE_IN_BOUND) THEN
        IF (myrank.eq.0) THEN
          iVect(1)=NUM_NETCDF_FILES_BND
          iVect(2)=MAX_NDT_BND_FILE
          iVect(3)=NDX_BND
          iVect(4)=NDY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(iVect,4,itype, iProc-1, 811, comm, ierr)
          END DO
          rVect(1)=OFFSET_X_BND
          rVect(2)=OFFSET_Y_BND
          rVect(3)=DX_BND
          rVect(4)=DY_BND
          DO IPROC=2,nproc
            CALL MPI_SEND(rVect,4,rtype, iProc-1, 820, comm, ierr)
          END DO
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          DO IPROC=2,nproc
            CALL MPI_SEND(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(BND_TIME_ALL_FILES,eSize,rtype, iProc-1, 813, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(iVect,4,itype, 0, 811, comm, istatus, ierr)
          NUM_NETCDF_FILES_BND=iVect(1)
          MAX_NDT_BND_FILE=iVect(2)
          NDX_BND=iVect(3)
          NDY_BND=iVect(4)
          CALL MPI_RECV(rVect,4,rtype, 0, 820, comm, istatus, ierr)
          OFFSET_X_BND=rVect(1)
          OFFSET_Y_BND=rVect(2)
          DX_BND=rVect(3)
          DY_BND=rVect(4)
          ALLOCATE(BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAX_NDT_BND_FILE), NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
          IF (istat/=0)   & 
             CALL parallel_abort('swan_bdcons, allocate error 14')

          CALL MPI_RECV(NDT_BND_FILE,NUM_NETCDF_FILES_BND,itype, 0, 812, comm, istatus, ierr)
          eSize=NUM_NETCDF_FILES_BND*MAX_NDT_BND_FILE
          CALL MPI_RECV(BND_TIME_ALL_FILES,eSize,rtype, 0, 813, comm, istatus, ierr)
        END IF
      !END IF
      END IF

!      IF(myrank.EQ.0) WRITE(PRINTF,*) 'BND_TIME_ALL_FILES=',BND_TIME_ALL_FILES

      WBDELT =  (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC
      IF(myrank.EQ.0) &
      WRITE(PRINTF,*) 'WW3 Frequency (sec) =',REAL(WBDELT)

      ! Start Time in mjd (days)
      ! End Time in mjd (days) 
      WBBMJD = BND_TIME_ALL_FILES(1,1)
      WBEMJD = BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,              &
               MAX_NDT_BND_FILE)

      ! Ultimate Time bound checking
      IF(myrank.EQ.0) THEN
       WRITE(PRINTF,*) 'SWAN Start (day):',TINIC*SEC2DAY
       WRITE(PRINTF,*) 'SWAN End   (day):',TFINC*SEC2DAY
       WRITE(PRINTF,*) 'First WW3 record (day):',WBBMJD
       WRITE(PRINTF,*) 'Last  WW3 record (day):',WBEMJD
      ENDIF

      IF( (WBBMJD.GT.TINIC*SEC2DAY) .OR.          &
          (WBEMJD.LT.TFINC*SEC2DAY)) THEN
         CALL parallel_abort( &
         'no appropriate time exists in WW3 BC files')
      ENDIF

      ALLOCATE (HS_WW3(NDX_BND,NDY_BND),stat=istat); HS_WW3 = 0._rkind
      ALLOCATE (FP_WW3(NDX_BND,NDY_BND),stat=istat); FP_WW3 = 0._rkind
      ALLOCATE (T02_WW3(NDX_BND,NDY_BND),stat=istat); T02_WW3 = 0._rkind
      ALLOCATE (DSPR_WW3(NDX_BND,NDY_BND),stat=istat); DSPR_WW3 = 0._rkind
      ALLOCATE (DIR_WW3(NDX_BND,NDY_BND), stat=istat); DIR_WW3 = 0._rkind
      IF (istat/=0) CALL parallel_abort('swan_bdcons, allocate error 15')

      END SUBROUTINE INIT_NETCDF_WW3_WAVEPARAMETER

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3_IVAR(JFILE, IT, IVAR, VAR_READ)
      USE VARS_WAVE
      USE schism_glbl,ONLY: skind,rkind
      USE schism_msgp,ONLY: parallel_abort
      USE NETCDF
      IMPLICIT NONE

      integer, intent(in) :: JFILE, IT, IVAR
      character (len = *), parameter :: CallFct = "READ_NETCDF_WW3_IVAR"
      CHARACTER(LEN=40)  :: EVAR
      REAL(rkind), intent(out) :: VAR_READ(NDX_BND,NDY_BND)
      integer :: ITMP(NDX_BND,NDY_BND)
      REAL(skind) :: RTMP(NDX_BND,NDY_BND)
      REAL(rkind) :: scale_factor
      integer :: ncid, var_id, ISTAT
      LOGICAL :: lexist

      IF (IVAR .eq. 3) THEN
        EVAR='hs'
      ELSE IF (IVAR .eq. 2) THEN
        EVAR='fp'
      ELSE IF (IVAR .eq. 5) THEN
        EVAR='t02'
      ELSE IF (IVAR .eq. 4) THEN
        EVAR='spr'
      ELSE IF (IVAR .eq. 1) THEN
        EVAR='dir'
      ELSE
        Print *, 'IVAR=', IVAR
        CALL parallel_abort('Wrong IVAR')
      END IF

!      IF(myrank.EQ.0) THEN
!       print*,'READ_NETCDF_WW3_IVAR reading:',TRIM(NETCDF_FILE_NAMES_BND(JFILE,IVAR))
!       print*,'READ_NETCDF_WW3_IVAR reading vname:',TRIM(EVAR)
!       print*,'READ_NETCDF_WW3_IVAR reading IT',IT
!      ENDIF 

      inquire(file=TRIM(NETCDF_FILE_NAMES_BND(JFILE,IVAR)),exist=lexist)
      IF(.NOT.lexist) &
            CALL parallel_abort("Missing WW3 boundary condition file")

      ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(JFILE,IVAR),NF90_NOWRITE,ncid)
      if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Open error")

      ISTAT = nf90_inq_varid(ncid, TRIM(EVAR), var_id)
      if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Inquire varid error")

!      ISTAT = nf90_get_att(ncid, var_id, 'scale_factor', scale_factor)
!      if(ISTAT.ne.NF90_NOERR) &
!            call parallel_abort("Netcdf Get Att scale_factor error")
      scale_factor = 1._rkind

!      ISTAT = NF90_GET_VAR(ncid, var_id, ITMP, &
!       start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
!      if(ISTAT.ne.NF90_NOERR) &
!            call parallel_abort("Netcdf Get var error")
!      VAR_READ = MyREAL(ITMP) * scale_factor

      ISTAT = NF90_GET_VAR(ncid, var_id, RTMP ,  &
      start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
      if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Get var error")
      VAR_READ = MyREAL(RTMP) * scale_factor

      ISTAT = nf90_close(ncid)
      if(ISTAT.ne.NF90_NOERR) &
            call parallel_abort("Netcdf Close error")

      END SUBROUTINE READ_NETCDF_WW3_IVAR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3_SINGLE(JFILE,IT)
      USE VARS_WAVE
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JFILE, IT
      CALL READ_NETCDF_WW3_IVAR(JFILE, IT, 3, HS_WW3)
      CALL READ_NETCDF_WW3_IVAR(JFILE, IT, 2, FP_WW3)
      CALL READ_NETCDF_WW3_IVAR(JFILE, IT, 5, T02_WW3)
      CALL READ_NETCDF_WW3_IVAR(JFILE, IT, 4, DSPR_WW3)
      CALL READ_NETCDF_WW3_IVAR(JFILE, IT, 1, DIR_WW3)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3_PARAM

      USE VARS_WAVE
      USE OCPCOMM4
      USE schism_glbl, ONLY: rkind
      USE schism_msgp
      IMPLICIT NONE

      include 'mpif.h'

      INTEGER     :: JFILE, IT
      INTEGER     :: counter, ip, i, j
      REAL(rkind), ALLOCATABLE    :: U(:), V(:), H(:)
      REAL(rkind), SAVE           :: TIME, scale_factor
      INTEGER :: IX, IY, IPROC, ISTAT
      REAL(rkind), ALLOCATABLE :: ARR_send_recv(:)
      integer, allocatable :: bnd_rqst(:), bnd_stat(:,:)

      CALL COMPUTE_IFILE_IT(JFILE, IT)
!      WRITE(PRINTF,*)'READ_NETCDF_WW3_PARAM JFILE, IT=',JFILE, IT

      IF(.NOT.SERIAL) THEN 

      !IF (MULTIPLE_IN_BOUND) THEN
      !  CALL READ_NETCDF_WW3_SINGLE(JFILE,IT)
      !ELSE
        allocate(ARR_send_recv(5*NDX_BND*NDY_BND))
        IF (istat/=0) &
            CALL parallel_abort('swan_bdcons, allocate error 16')
        IF (myrank .eq. 0) THEN
          CALL READ_NETCDF_WW3_SINGLE(JFILE,IT)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              ARR_send_recv(J)=HS_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=FP_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=T02_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DSPR_WW3(IX,IY)
              J=J+1
              ARR_send_recv(J)=DIR_WW3(IX,IY)
            END DO
          END DO
          allocate(bnd_rqst(nproc-1), bnd_stat(MPI_STATUS_SIZE,nproc-1), stat=istat)
          IF (istat/=0) &
              CALL parallel_abort('swan_bdcons, allocate error 17')
          DO iProc=2,nproc
            CALL mpi_isend(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, iProc-1, 2032, comm, bnd_rqst(iProc-1), ierr)
          END DO
          IF (nproc > 1) THEN
            CALL MPI_WAITALL(nproc-1, bnd_rqst, bnd_stat, ierr)
          END IF
          deallocate(bnd_rqst, bnd_stat)
        ELSE
          CALL MPI_RECV(ARR_send_recv, 5*NDX_BND*NDY_BND, rtype, 0, 2032, comm, istatus, ierr)
          J=0
          DO IX=1,NDX_BND
            DO IY=1,NDY_BND
              J=J+1
              HS_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              FP_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              T02_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DSPR_WW3(IX,IY) = ARR_send_recv(J)
              J=J+1
              DIR_WW3(IX,IY) = ARR_send_recv(J)
            END DO
          END DO
        END IF
      !END IF

      ELSE ! SERIAL

       CALL READ_NETCDF_WW3_SINGLE(JFILE,IT)

      ENDIF ! SERIAL

#if 0
      IF (LWRITE_WW3_RESULTS) THEN
        OPEN(3012, FILE  = 'ergwiii.bin', FORM = 'UNFORMATTED')
        ALLOCATE(U(NDX_BND*NDY_BND), V(NDX_BND*NDY_BND), H(NDX_BND*NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons, allocate error 18')
        COUNTER = 1
        DO J = 1, NDY_BND
          DO I = 1, NDX_BND
            U(COUNTER) = HS_WW3(I,J)
            V(COUNTER) = DIR_WW3(I,J)
            H(COUNTER) = DSPR_WW3(I,J)
            COUNTER = COUNTER + 1
          END DO
        END DO
        TIME = TIME + 1.
        WRITE(3012) TIME
        WRITE(3012) (U(IP), V(IP), H(IP), IP = 1, NDX_BND*NDY_BND)
        DEALLOCATE(U,V,H)
      END IF
#endif

      END SUBROUTINE READ_NETCDF_WW3_PARAM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, &
     &                               MNPT, XPT, YPT, VAL, DOPEAK)

      USE VARS_WAVE, ONLY:IOBCN_W,HS_WW3,DIR_WW3,FP_WW3,T02_WW3,DSPR_WW3
      USE schism_glbl, ONLY:iplg,rkind
      USE schism_msgp, ONLY:parallel_abort

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      INTEGER, INTENT(IN)        :: MNPT
      REAL(rkind), INTENT(IN)    :: XPT(MNPT), YPT(MNPT)
      REAL(rkind), INTENT(OUT)   :: VAL(8,IOBCN_W)
      LOGICAL, INTENT(IN)        :: DOPEAK
      REAL(rkind) :: eVect(5)
      REAL(rkind) :: WX, WX1, WX2, WX3, WX4
      REAL(rkind) :: eVAL, HX1, HX2, LEN_X, LEN_Y
      REAL(rkind) :: DELTA_X, DELTA_Y
      INTEGER :: IVAR, I, J, IDX1, IDX2
      INTEGER :: IP, J_INT, I_INT
      DO IP = 1, MNPT
        LEN_X = XPT(IP) - OFFSET_X
        LEN_Y = YPT(IP) - OFFSET_Y
        I_INT = INT( LEN_X/DX ) + 1
        J_INT = INT( LEN_Y/DY ) + 1
        DELTA_X   = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
        DELTA_Y   = LEN_Y - (J_INT - 1) * DY !

        DO IVAR=1,5

          !IF(IVAR.EQ.1) THEN
          ! print*,'I_INT,J_INT',I_INT,J_INT
          !ENDIF

          DO I=0,1
            DO J=0,1
              IDX1=I_INT + I
              IDX2=J_INT + J
              IF (IVAR .eq. 1) THEN
                WX=HS_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 2) THEN
                WX=DIR_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 3) THEN
                WX=FP_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 4) THEN
                WX=T02_WW3(IDX1,IDX2)
              ELSE IF (IVAR .eq. 5) THEN
                WX=DSPR_WW3(IDX1,IDX2)
              ELSE
                CALL parallel_abort('Wrong IVAR')
              END IF
              IF ((I .eq. 0).and.(J .eq. 0)) THEN
                WX1=WX
              END IF
              IF ((I .eq. 0).and.(J .eq. 1)) THEN
                WX2=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 1)) THEN
                WX3=WX
              END IF
              IF ((I .eq. 1).and.(J .eq. 0)) THEN
                WX4=WX
              END IF
            END DO
          END DO
          HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
          HX2       = WX2 + (WX3-WX2)/DX * DELTA_X
          IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
            eVAL = 0.
          ELSE
            eVAL = HX1 + (HX2-HX1)/DY * DELTA_Y
          ENDIF
          eVect(IVAR)=eVAL
        END DO
        VAL(1,IP)=eVect(1)
        IF (DOPEAK) THEN
          IF (eVect(3) .gt. TINY(1.)) THEN
            VAL(2,IP)  = 1._rkind/eVect(3)
            VAL(5,IP)  = 2._rkind
          END IF
        ELSE
          VAL(2,IP) = eVect(4)
          VAL(5,IP) = -2.
        END IF
        VAL(3,IP)  = eVect(2)
        VAL(4,IP)  = eVect(5)
        VAL(6,IP)  = 1.
        VAL(7,IP)  = 0.1
        VAL(8,IP)  = 3.3

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!         1 - Pierson-Moskowitz
!         2 - JONSWAP
!         3 - BIN
!         4 - Gauss
!         positive peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTER_STRUCT_BOUNDARY(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)

      USE VARS_WAVE, ONLY: DOPEAK_BOUNDARY,IOBCN_W,I_OBC_N_W
      USE schism_glbl, ONLY:iplg,rkind
      USE swangriddata, ONLY:xcugrd,ycugrd
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NDX, NDY
      REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
      REAL(rkind), INTENT(OUT)   :: VAL(8,IOBCN_W)
      REAL(rkind)  :: XPT(IOBCN_W), YPT(IOBCN_W)
      INTEGER :: IP,ID
      DO IP=1,IOBCN_W
        ID = I_OBC_N_W(IP)
        XPT(IP)=xcugrd(ID)
        YPT(IP)=ycugrd(ID)
      END DO
      CALL SPPARM_INTER_STRUCT(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y, IOBCN_W, XPT, YPT, VAL, DOPEAK_BOUNDARY)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

      SUBROUTINE SSORT2 (X, Y, Z, N, KFLAG)

      USE schism_glbl, ONLY: rkind
      IMPLICIT NONE
!TS: Modified second array z(*) to carry along

!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
!     .. Scalar Arguments ..
      INTEGER KFLAG, N
!     .. Array Arguments ..
      REAL(rkind) X(*), Y(*), Z(*)
!     .. Local Scalars ..
      REAL(rkind) R, T, TT, TTY, TY, TZ, TTZ
      INTEGER I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!     None
!     .. Intrinsic Functions ..
      !INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SSORT
      NN = N
      IF (NN .LT. 1) THEN
        WRITE (*,*) 'The number of values to be sorted is not positive.'
         RETURN
      ENDIF
!
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         WRITE (*,*) 'The sort control parameter, &
     &                K, is not 2, 1, -1, or -2.'
         RETURN
      ENDIF
!
!     Alter array X to get decreasing order if needed
!
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF
!
      IF (KK .EQ. 2) GO TO 100
!
!     Sort X only
!
      M = 1
      I = 1
      J = NN
      R = 0.375E0
!
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
!
   30 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = X(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than than T, interchange with T
!
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50
!
!     Interchange these elements
!
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
!
!     Begin again on another portion of the unsorted array
!
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
!
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I
!
   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80
!
!     Sort X and carry Y along
!
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
!
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
!
  120 K = I
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)
      TZ = Z(IJ)
!
!     If first element of array is greater than T, interchange with T
!
      IF (X(I) .GT. T) THEN

         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)

         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)

         Z(IJ) = Z(I)
         Z(I) = TZ
         TZ = Z(IJ)
      ENDIF
      L = J
!
!     If last element of array is less than T, interchange with T
!
      IF (X(J) .LT. T) THEN

         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

         Z(IJ) = Z(J)
         Z(J) = TZ
         TZ = Z(IJ)
!
!        If first element of array is greater than T, interchange with T
!
         IF (X(I) .GT. T) THEN

            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)

            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)

            Z(IJ) = Z(I)
            Z(I) = TZ
            TZ = Z(IJ)

         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
!
!     Interchange these elements
!
      IF (K .LE. L) THEN

         TT = X(L)
         X(L) = X(K)
         X(K) = TT

         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY

         TTZ = Z(L)
         Z(L) = Z(K)
         Z(K) = TTZ

         GO TO 130
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
!
!     Begin again on another portion of the unsorted array
!
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
!
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
!
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      TZ = Z(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I
!
  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      Z(K+1) = Z(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      Z(K+1) = TZ
      GO TO 170
!
!     Clean up
!
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      END SUBROUTINE

