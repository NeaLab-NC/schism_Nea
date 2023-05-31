# if defined (USE_SWAN)
!==============================================================================|

MODULE VARS_WAVE

  USE schism_glbl, ONLY : rkind

  IMPLICIT NONE
  SAVE

  REAL(rkind),  PARAMETER      :: KDMAX = 10.D0
  REAL(rkind)                  :: ALPROL = 0.65

! Dirty Fix for "RADFLAG" (not implemented in param.nml yet)
! -----------------------
!  CHARACTER(LEN=3)             :: RADFLAG  = 'LON'
!  CHARACTER(LEN=3)             :: RADFLAG  = 'WON'  ! ala SWAN
  CHARACTER(LEN=3)             :: RADFLAG  = 'VOR'
! Dirty Fix for "IROLLER"  (not implemented in param.nml yet)
! -----------------------
! Some work todo for ROLLER in SWAN ... IROLLER = 0 
  INTEGER                      :: IROLLER = 0      

! Dirty Fix for "ZPROF_BREAK" (not implemented in param.nml yet)
!#ifdef USE_PLANARBEACH
!  INTEGER,  PARAMETER      :: ZPROF_BREAK = 2     ! The Type III from Uchiyama et al. / Ocean Modelling 34 (2010)
!#else
  INTEGER,  PARAMETER      :: ZPROF_BREAK = 6     ! Vertical distribution function of wave breaking source term, only used in 3D run
!#endif

! WWM like heritage
  REAL(rkind)    , ALLOCATABLE :: WK(:,:)
  REAL(rkind)    , ALLOCATABLE :: DS_BAND(:),DS_INCR(:)

  INTEGER :: INP_CUR_NTIME,INP_WI_NTIME,INP_FR_NTIME,INP_WLEV_NTIME
!   
  REAL(rkind),save,ALLOCATABLE :: AC2(:,:,:), COMPDA(:,:)
  REAL(rkind),save,ALLOCATABLE :: AC1(:,:,:)

  real(rkind),save,allocatable :: wwave_force(:,:,:)   ! wave (body) force 
  real(rkind),save,allocatable :: jpress(:)            ! wave-induced pressure
  real(rkind),save,allocatable :: sbr(:,:)             ! sbr(2,npa): momentum flux vector due to wave breaking 
                                                       ! (nearshore depth-induced breaking; see Bennis 2011)
  real(rkind),save,allocatable :: sbf(:,:)             ! sbf(2,npa): momentum lost by waves due to the bottom friction
  real(rkind),save,allocatable :: srol(:,:)            ! momentum flux vector due to roller
  real(rkind),save,allocatable :: sds(:,:)             ! momentum flux vector due to whitecapping

  real(rkind),save,allocatable :: stokes_hvel(:,:,:)            ! horizontal Stokes drift velocities (u,v) on node
  real(rkind),save,allocatable :: stokes_wvel(:,:)              ! vertical Stokes drift velocities (w) on node
  real(rkind),save,allocatable :: stokes_hvel_side(:,:,:)       ! horizontal Stokes drift velocities (u,v) on side
  real(rkind),save,allocatable :: stokes_wvel_side(:,:)         ! vertical  Stokes drift velocities (w) on side

  real(rkind),save,allocatable :: roller_stokes_hvel(:,:,:)     ! horizontal roller velocities (u,v) on node 
  real(rkind),save,allocatable :: roller_stokes_hvel_side(:,:,:)! horizontal roller velocities (u,v) on side
  real(rkind),save,allocatable :: eps_w(:),eps_r(:),eps_br(:)  ! Roller terms

  real(rkind),save,allocatable :: out_wwm(:,:)
  real(rkind),save,allocatable :: out_wwm_rol(:,:)
  real(rkind),save,allocatable :: out_wwm_windpar(:,:)

  real(rkind),save,allocatable :: tanbeta_x(:), tanbeta_y(:) ! MP from KM: bottom slope
  real(rkind),save,allocatable :: curx_wwm(:), cury_wwm(:)  !BM: coupling current
  real(rkind),save,allocatable :: wafo_opbnd_ramp(:) !BM: ramp on wave forces at open boundary
  real(rkind),save,allocatable :: taub_wc(:)
  real(rkind),save,allocatable :: delta_wbl(:)

  real(rkind),save,allocatable :: wave_sbrtot(:)
  real(rkind),save,allocatable :: wave_sbftot(:)
  real(rkind),save,allocatable :: wave_sintot(:)
  real(rkind),save,allocatable :: wave_sdstot(:)
  real(rkind),save,allocatable :: wave_svegtot(:)

  real(rkind),allocatable :: DROLP(:)           ! Mean wave direction [rad]: surface roller direction of propagation
  real(rkind),allocatable :: EROL1(:), EROL2(:) ! Bulk energy of surface rollers [m^3/s^2]
  real(rkind),allocatable :: CROLP(:)           ! Peak wave phase [m/s]: surface roller propagation speed, used for computing dissipation

! JL : same array names and same output fields as from WWM...
!  REAL(rkind), ALLOCATABLE :: out_wwm(:,:)
!  REAL(rkind), ALLOCATABLE :: out_wwm_windpar(:,:)
  REAL(rkind), ALLOCATABLE :: SPEC_DENSITY(:,:)
! REAL(rkind), ALLOCATABLE :: DISSURF(:),DISWCAP(:),DISBOT(:),QB(:) into out_wwm

! JL: see swanpre1.F90
  INTEGER,save                :: IOBCN_W       !!LOCAL NUMBER OF OPEN BOUNDARY NODES
  INTEGER,save                :: IOBCN_GL_W    !!GLOBAL NUMBER OF OPEN BOUNDARY NODES
  INTEGER,save,ALLOCATABLE    :: I_OBC_N_W(:)  !!OPEN BOUNDARY NODE LIST FOR SWAN (Local)
  INTEGER,save,ALLOCATABLE    :: I_OBC_GL_W(:) !!OPEN BOUNDARY NODE LIST FOR SWAN (Global)

! JL Experimental : load KN from kn.gr3 and pass to SWAN
  REAL(rkind),save, ALLOCATABLE, target :: KN(:) ! Roughness length (Nikuradse) from Madsen 1988

  LOGICAL :: SERIAL
  LOGICAL :: NESTING = .FALSE.

  CHARACTER(LEN=80) UGSWAN_VERSION !!STRING DESCRIBING VERSION
 
! NESTING : wve parameters or spectral density along OBC, wwm heritage
  CHARACTER(LEN=120) :: NESTING_FILE
  CHARACTER(LEN=80)  :: NESTING_TYPE_WAVE 
  CHARACTER(LEN=40)  :: NCDF_HS_NAME   = 'hs'
  CHARACTER(LEN=40)  :: NCDF_DIR_NAME  = 'dir'
  CHARACTER(LEN=40)  :: NCDF_SPR_NAME  = 'spr'
  CHARACTER(LEN=40)  :: NCDF_FP_NAME   = 'fp'
  CHARACTER(LEN=40)  :: NCDF_F02_NAME  = 't02'
  CHARACTER(LEN=200), ALLOCATABLE :: NETCDF_FILE_NAMES_BND(:,:)
  INTEGER                         :: NDX_BND, NDY_BND
  INTEGER                         :: NDT_BND_ALL_FILES
  INTEGER                         :: NUM_NETCDF_FILES_BND
  INTEGER, ALLOCATABLE            :: NDT_BND_FILE(:)
  INTEGER                         :: NUM_NETCDF_VAR_TYPES = 5
  REAL(rkind), ALLOCATABLE        :: BND_TIME_ALL_FILES(:,:)
  REAL(rkind),   ALLOCATABLE      :: COORD_BND_Y(:)
  REAL(rkind),   ALLOCATABLE      :: COORD_BND_X(:)
  REAL(rkind)                     :: DX_BND, DY_BND
  REAL(rkind)                     :: OFFSET_X_BND, OFFSET_Y_BND
!
! WAM MultiModal parametric wave forcing at boundary nodes (N_OBN_WAM) 
!     and splitted in N_WAVE_WAM partitions
  REAL(rkind),   ALLOCATABLE      :: HS_WAM(:,:)
  REAL(rkind),   ALLOCATABLE      :: PER_WAM(:,:)
  REAL(rkind),   ALLOCATABLE      :: DIR_WAM(:,:)
  REAL(rkind),   ALLOCATABLE      :: DSPR_WAM(:,:)
  INTEGER                         :: N_OBN_WAM, N_WAVE_WAM
! Bulk spectrum reconstruction settings  
  INTEGER                         :: FSHAPE_WAM      ! 1:PM,2:JON,3:BIN,4:GAUS 
  INTEGER                         :: CHAR_WAM_PERIOD ! 1:PEAK or 2:MEAN frequency
  INTEGER                         :: CHAR_WAM_DSPR   ! 1:DEGR, 2:POWER directional distribution 
! 
! WW3 parametric wave forcing
  REAL(rkind),   ALLOCATABLE      :: HS_WW3(:,:)
  REAL(rkind),   ALLOCATABLE      :: T02_WW3(:,:)
  REAL(rkind),   ALLOCATABLE      :: DIR_WW3(:,:)
  REAL(rkind),   ALLOCATABLE      :: FP_WW3(:,:)
  REAL(rkind),   ALLOCATABLE      :: DSPR_WW3(:,:)
! WW3 spectra 
  !REAL(rkind),   ALLOCATABLE             :: ALL_VAR_WW3(:,:,:)
  INTEGER            :: NP_WW3, MSC_WW3, MDC_WW3, MAXSTEP_WW3, TSTART_WW3(2)
  REAL(rkind)                     :: DTBOUND_WW3, DDIR_WW3
  REAL(rkind),   ALLOCATABLE      :: FQ_WW3(:)
  REAL(rkind),   ALLOCATABLE      :: DR_WW3(:)
  REAL(rkind),   ALLOCATABLE      :: XP_WW3(:), YP_WW3(:)
!
! ... wave boundary stuff, see swan_bdcons, wwm heritage
!
  REAL(rkind), ALLOCATABLE    :: WBAC   (:,:,:)
  REAL(rkind), ALLOCATABLE    :: WBAC_GL(:,:,:)
  REAL(rkind), ALLOCATABLE    :: WBACOLD(:,:,:)
  REAL(rkind), ALLOCATABLE    :: WBACNEW(:,:,:)
  REAL(rkind), ALLOCATABLE    :: DSPEC  (:,:,:)
  REAL(rkind), ALLOCATABLE    :: SPEG   (:,:,:)
  REAL(rkind)                 :: WBDELT
  REAL(rkind)                 :: WBTMJD
  REAL(rkind)                 :: WBBMJD
  REAL(rkind)                 :: WBEMJD
  REAL(rkind), ALLOCATABLE    :: WBSPPARM(:,:)
!  LOGICAL    :: DOPEAK_BOUNDARY = .TRUE. !if true error in  SPPARM_INTER_STRUCT
  LOGICAL    :: DOPEAK_BOUNDARY = .FALSE.
 
  INTEGER , ALLOCATABLE :: CROSS(:)
  INTEGER , ALLOCATABLE :: BGRIDP(:)                                   

  REAL(rkind)    , ALLOCATABLE :: BSPECS(:,:,:,:)                             
  REAL(rkind)    , ALLOCATABLE :: Sice(:,:,:)
  
  REAL(rkind)    , ALLOCATABLE :: BLKND(:), BLKNDC(:), OURQT(:)   
  INTEGER :: ITW,IT0            

  ! For Checkout (NaN)
  REAL(rkind)   , ALLOCATABLE :: AC2_SUM(:)
  REAL(rkind)   :: SUM_AC2

  ! Used in RADIATION_STRESS_SCHISM 

  LOGICAL       :: LETOT       = .FALSE.  ! JL ? WWM,see RADIATION_STRESS_SCHISM
  INTEGER                  :: IMET_DRY
  REAL(rkind), ALLOCATABLE :: RSXX(:), RSXY(:), RSYY(:) !FORCEXY(:,:)
  REAL(rkind), ALLOCATABLE :: SXX3D(:,:), SXY3D(:,:), SYY3D(:,:)

  ! Used in SwanComputeForce
  REAL(rkind), ALLOCATABLE :: WAVESTRX_2D(:),WAVESTRY_2D(:)
  REAL(rkind), ALLOCATABLE :: WAVESTRX_3D(:,:),WAVESTRY_3D(:,:)

  REAL(rkind) :: SOURCE_DTMAX,SOURCE_DTMIN

  CHARACTER(LEN=4) COMPUT

  REAL(rkind),  PARAMETER            :: DAY2SEC  = 86400.d0
  REAL(rkind),  PARAMETER            :: SEC2DAY  = 1.d0/DAY2SEC

END MODULE VARS_WAVE
# endif
