!******************************************************************
!
!   SUBROUTINE SWANCOMPUNSTRUC(ac2, ac1, compda, spcsig, spcdir, xytst, cross, it)
   SUBROUTINE SWANCOMPUNSTRUC( spcsig, spcdir, xytst, IT )
!
!******************************************************************
!subroutine SwanCompUnstruc ( ac2, ac1, COMPDA, spcsig, spcdir, xytst, cross, it )
!
!   Authors
!
!   40.80: Marcel Zijlema
!   40.85: Marcel Zijlema
!   40.95: Marcel Zijlema
!   41.02: Marcel Zijlema
!   41.07: Marcel Zijlema
!   41.10: Marcel Zijlema
!   41.20: Casey Dietrich
!   41.60: Marcel Zijlema
!   41.63: Marcel Zijlema
!
!   Updates
!
!   40.80,     July 2007: New subroutine
!   40.85,   August 2008: add propagation, generation and redistribution terms for output purposes
!   40.95,     June 2008: parallelization of unSWAN using MESSENGER of ADCIRC
!   41.02, February 2009: implementation of diffraction
!   41.07,   August 2009: bug fix: never-ending sweep is prevented
!   41.10,   August 2009: parallelization using OpenMP directives
!   41.20,     June 2010: extension to tightly coupled ADCIRC+SWAN model
!   41.60,     July 2015: more accurate computation of gradients of depth or wave number for turning rate
!   41.63,   August 2015: efficiency UnSWAN algorithm further improved; now one sweep per iteration
!
!   Purpose
!
!   Performs one time step for solution of wave action equation on unstructured grid
!
!   Method
!
!   A vertex-based algorithm is employed in which the variables are stored at the vertices of the mesh
!   The equation is solved in each vertex assuming a constant spectral grid resolution in all vertices
!   The propagation terms in both geographic and spectral spaces are integrated implicitly
!   Sources are treated explicitly and sinks implicitly
!   The calculation of the source terms is carried out in the original SWAN routines, e.g., SOURCE
!
!   The wave action equation is solved iteratively
!   A number of iterations are carried out until convergence is reached
!   In each iteration, a sweep through the vertices is executed
!   The solution of each vertex must be updated geographically before proceeding to the next one
!   The two upwave faces connecting the vertex to be updated enclose those wave directions that can be processed in the spectral space
!   A sweep is complete when all vertices are updated geographically (but not necessarily in whole spectral space)
!   An iteration is complete when all vertices are updated in both geographic and spectral spaces
!
!   Modules used
!
!    use ocpcomm4
!    use swcomm1
!    use swcomm2
!    use swcomm3
!    use swcomm4
!    use SwanGriddata
!    use SwanGridobjects
!    use SwanCompdata
!    use m_parall
!PUN    use SIZES, only: SZ, MNPROC
!PUN    use MESSENGER
!ADC    use couple2adcirc, only: MakeBoundariesReflective
!ADC    use NodalAttributes, only: FoundSwanWaveRefrac, LoadSwanWaveRefrac, SwanWaveRefrac
!

!   USE ALL_VARS, only:
!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!#  endif

   USE TIMECOMM
   USE OCPCOMM1
   USE OCPCOMM2
   USE OCPCOMM3
   USE OCPCOMM4
   USE SWCOMM1
   USE SWCOMM2
   USE SWCOMM3
   USE SWCOMM4
   USE SWANGRIDDATA
   USE SWANGRIDOBJECTS
   USE SWANCOMPDATA

!   USE M_GENARR, ONLY:XYTST,SPCSIG,SPCDIR

   USE VARS_WAVE
!  USE MOD_PREC

!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!!   USE MOD_NESTING, ONLY : NESTING_ON_WAVE, NESTING_TYPE_WAVE, SET_VAR_WAVE, SET_VAR_WAVE_AC2
!#  endif

#  if defined (WAVE_SETUP)
   USE MOD_WAVESETUP
#  endif

#  if defined (PLBC)
   USE MOD_PERIODIC_LBC
#  endif

!   USE ALL_VARS
   USE schism_glbl, ONLY : skind, rkind, dkind, nea, npa, iplg, ielg
!   USE schism_glbl, ONLY :, uu2, vv2, tr_nd, & !tnd, snd, &
!   USE schism_glbl, ONLY : kfp, idry, nvrt, ivcor,ipgl,fdb,lfdb,errmsg,     &
!   USE schism_glbl, ONLY : np_global,ne_global,                             &
!   USE schism_glbl, ONLY : np,ne,nvrt,i34,isbe,isbs,elside,elnode,ilnd,nlnd,nland,isbnd

!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!#  endif
   USE schism_msgp

!   USE TIMECOMM
!   USE OCPCOMM1
!   USE OCPCOMM2
!   USE OCPCOMM3
!   USE OCPCOMM4
!   USE SWCOMM1
!   USE SWCOMM2
!   USE SWCOMM3
!   USE SWCOMM4
!   USE M_GENARR
!   USE m_constants
!   USE m_xnldata
!   USE m_fileio

!#  if defined (EXPLICIT)
!   USE MOD_ACTION_EX
!#  else
!   USE MOD_ACTION_IM
!#  endif
!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!   USE MOD_NESTING, ONLY : NESTING_ON_WAVE, NESTING_TYPE_WAVE, SET_VAR_WAVE,SET_VAR_WAVE_AC2
!#  endif
!   USE ALL_VARS    !, ONLY : M,MT,NGL,SERIAL,PAR,MYID,NPROCS,MSR,     &
!                   !     IOBCN,I_OBC_N,AC2,COMPDA,ISBCE
!   USE VARS_WAVE
!#  if defined (WAVE_SETUP)
!   USE MOD_WAVESETUP
!#  endif

!#  if defined (PLBC)
!   USE MOD_PERIODIC_LBC
!#  endif


    IMPLICIT NONE


!#ifndef USE_MPIMODULE
   include 'mpif.h'
!#endif

!
!     ************************************************************************
!     *                                                                      *
!     *                  MAIN SUBROUTINE OF COMPUTATIONAL PART               *
!     *                                                                      *
!     *                               -- SWANCOMPUNSTRUC --                  *
!     *                                                                      *
!     ************************************************************************
!
!   Argument variables
!
!    integer, dimension(nfaces), intent(in)         :: cross  ! contains sequence number of obstacles for each face
                                                             ! where they crossing or zero if no crossing
    integer, intent(in)                            :: IT     ! counter in time space
    integer, dimension(NPTST), intent(in)          :: xytst  ! test points for output purposes
!
!    real(rkind), dimension(MDC,MSC,nverts), intent(in)    :: ac1    ! action density at previous time level
!    real(rkind), dimension(MDC,MSC,nverts), intent(inout) :: ac2    ! action density at current time level
!    real(rkind), dimension(nverts,MCMVAR), intent(inout)  :: COMPDA ! array containing space-dependent info (e.g. depth)
    real(rkind), dimension(MDC,6), intent(in)             :: spcdir ! (*,1): spectral direction bins (radians)
                                                             ! (*,2): cosine of spectral directions
                                                             ! (*,3): sine of spectral directions
                                                             ! (*,4): cosine^2 of spectral directions
                                                             ! (*,5): cosine*sine of spectral directions
                                                             ! (*,6): sine^2 of spectral directions
   real(rkind), dimension(MSC), intent(in)               :: spcsig ! relative frequency bins
!
!   Local variables
!
!   For Parallel run :
    REAL(dkind), DIMENSION(4) :: SBUF,RBUF1            ! For Sum/Max/Min of quantitie over partition
    REAL(rkind), ALLOCATABLE :: COMPDA_TMP1(:)         ! For Data exchange over partition
    INTEGER :: JDUM, IP, IDC, ISC

    CHARACTER*80 :: MSGSTR
    integer                               :: icell     ! cell index / loop counter
    integer                               :: id        ! loop counter over direction bins
    integer                               :: iddlow    ! minimum direction bin that is propagated within a sweep
    integer                               :: iddtop    ! maximum direction bin that is propagated within a sweep
    integer, parameter                    :: idebug=1  ! level of debug output:
                                                       ! 0 = no output
                                                       ! 1 = print extra output for debug purposes
    integer                               :: idtot     ! maximum number of bins in directional space for considered sweep
    integer                               :: idwmax    ! maximum counter for spectral wind direction
    integer                               :: idwmin    ! minimum counter for spectral wind direction
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer                               :: iface     ! face index
    integer                               :: inocnt    ! inocnv counter for calling thread
    integer                               :: inocnv    ! integer indicating number of vertices in which solver does not converged
    integer                               :: is        ! loop counter over frequency bins
    integer                               :: isslow    ! minimum frequency that is propagated within a sweep
    integer                               :: isstop    ! maximum frequency that is propagated within a sweep
    integer                               :: istat     ! indicate status of allocation
    integer                               :: istot     ! maximum number of bins in frequency space for considered sweep
    integer                               :: iter      ! iteration counter
    integer                               :: ivert     ! vertex index
    integer                               :: ivlow     ! lower index in range of vertices in calling thread
    integer                               :: ivup      ! upper index in range of vertices in calling thread
    integer                               :: j         ! loop counter
    integer                               :: jc        ! loop counter
    integer                               :: k         ! loop counter
    integer                               :: kvert     ! loop counter over vertices
    integer, dimension(2)                 :: link      ! local face connected to considered vertex where an obstacle crossed
    integer                               :: mnisl     ! minimum sigma-index occured in applying limiter
    integer                               :: mxnfl     ! maximum number of use of limiter in spectral space
    integer                               :: mxnfr     ! maximum number of use of rescaling in spectral space
    integer                               :: npfl      ! number of vertices in which limiter is used
    integer                               :: npfr      ! number of vertices in which rescaling is used
    integer                               :: swpnr     ! sweep number
    integer                               :: tid       ! thread number
    integer, dimension(3)                 :: vv        ! vertices in present cell
    integer                               :: vb        ! vertex of begin of present face
    integer                               :: ve        ! vertex of end of present face
    integer, dimension(2)                 :: vu        ! upwave vertices in present cell
    integer, dimension(24)                :: wwint     ! counters for 4 wave-wave interactions (see routine FAC4WW)
    !
    integer, dimension(:), allocatable    :: idcmax    ! maximum frequency-dependent counter in directional space
    integer, dimension(:), allocatable    :: idcmin    ! minimum frequency-dependent counter in directional space
    integer, dimension(:), allocatable    :: iscmax    ! maximum direction-dependent counter in frequency space
    integer, dimension(:), allocatable    :: iscmin    ! minimum direction-dependent counter in frequency space
    integer, dimension(:), allocatable    :: islmin    ! lowest sigma-index occured in applying limiter
    integer, dimension(:), allocatable    :: nflim     ! number of frequency use of limiter in each vertex
    integer, dimension(:), allocatable    :: nrscal    ! number of frequency use of rescaling in each vertex
!$  integer, dimension(:), allocatable    :: tlist     ! vertex list for calling thread
    !
    real(rkind)                                  :: abrbot    ! near bottom excursion
    real(rkind)                                  :: accur     ! percentage of active vertices in which required accuracy has been reached
    real(rkind)                                  :: acnrmo    ! norm of difference of previous iteration
    real(rkind), dimension(2)                    :: acnrms    ! array containing infinity norms
    real(rkind)                                  :: dal1      ! a coefficent for the 4 wave-wave interactions
    real(rkind)                                  :: dal2      ! another coefficent for the 4 wave-wave interactions
    real(rkind)                                  :: dal3      ! just another coefficent for the 4 wave-wave interactions
    real(rkind)                                  :: dhdx      ! derivative of depth in x-direction
    real(rkind)                                  :: dhdy      ! derivative of depth in y-direction
!PUN    real(SZ), dimension(1)                :: dum1      ! a dummy real(rkind) meant for UPDATER
!PUN    real(SZ), dimension(1)                :: dum2      ! a dummy real(rkind) meant for UPDATER
    real(rkind)                                  :: dummy     ! dummy variable (to be used in existing SWAN routine call)
    real(rkind)                                  :: etot      ! total wave energy density
    real(rkind)                                  :: fpm       ! Pierson Moskowitz frequency
    real(rkind)                                  :: frac      ! fraction of total active vertices
    real(rkind)                                  :: hm        ! maximum wave height
    real(rkind)                                  :: hs        ! significant wave height
    real(rkind)                                  :: kmespc    ! mean average wavenumber based on the WAM formulation
    real(rkind)                                  :: kteta     ! number of directional partitions
    real(rkind)                                  :: nwetp     ! total number of active vertices
    real(rkind)                                  :: qbloc     ! fraction of breaking waves
    real(rkind), dimension(2)                    :: rdx       ! first component of contravariant base vector rdx(b) = a^(b)_1
    real(rkind), dimension(2)                    :: rdy       ! second component of contravariant base vector rdy(b) = a^(b)_2
    real(rkind)                                  :: rhof      ! asymptotic convergence factor
    real(rkind)                                  :: rval1     ! a dummy value
    real(rkind)                                  :: rval2     ! a dummy value
    real(rkind)                                  :: smebrk    ! mean frequency based on the first order moment
    real(rkind)                                  :: snlc1     ! a coefficent for the 4 wave-wave interactions
    real(rkind)                                  :: stopcr    ! stopping criterion for stationary solution
    real(rkind)                                  :: thetaw    ! mean direction of the wind speed vector with respect to ambient current
    real(rkind)                                  :: ufric     ! wind friction velocity
    real(rkind), dimension(5)                    :: usrset    ! auxiliary array to store user-defined settings of 3rd generation mode
    real(rkind)                                  :: wind10    ! magnitude of the wind speed vector with respect to ambient current
    real(rkind), dimension(8)                    :: wwawg     ! weight coefficients for the 4 wave-wave interactions (see routine FAC4WW)
    real(rkind), dimension(8)                    :: wwswg     ! weight coefficients for the 4 wave-wave interactions semi-implicitly (see routine FAC4WW)
    real(rkind)                                  :: xis       ! difference between succeeding frequencies for computing 4 wave-wave interactions
    !
    real(rkind), dimension(:,:), allocatable     :: ac2old    ! array to store action density before solving system of equations
    real(rkind), dimension(:,:), allocatable     :: alimw     ! maximum energy by wind growth
                                                       ! this auxiliary array is used because the maximum value has to be checked
                                                       ! direct after solving the action balance equation
    real(rkind), dimension(:,:,:), allocatable   :: amat      ! coefficient matrix of system of equations in spectral space
    real(rkind), dimension(:,:), allocatable     :: rhs       ! right-hand side of system of equations in spectral space
    real(rkind), dimension(:,:), allocatable     :: cad       ! wave transport velocity in theta-direction
    real(rkind), dimension(:,:), allocatable     :: cas       ! wave transport velocity in sigma-direction
    real(rkind), dimension(:,:,:), allocatable   :: cax       ! wave transport velocity in x-direction
    real(rkind), dimension(:,:,:), allocatable   :: cay       ! wave transport velocity in y-direction
    real(rkind), dimension(:,:), allocatable     :: cgo       ! group velocity
    real(rkind), dimension(:,:), allocatable     :: da1c      ! implicit interaction contribution of first quadruplet, current bin (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: da1m      ! implicit interaction contribution of first quadruplet, current bin -1 (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: da1p      ! implicit interaction contribution of first quadruplet, current bin +1 (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: da2c      ! implicit interaction contribution of second quadruplet, current bin (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: da2m      ! implicit interaction contribution of second quadruplet, current bin -1 (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: da2p      ! implicit interaction contribution of second quadruplet, current bin +1 (unfolded space)
    real(rkind), dimension(:,:,:), allocatable   :: disc0     ! explicit part of dissipation in present vertex for output purposes
    real(rkind), dimension(:,:,:), allocatable   :: disc1     ! implicit part of dissipation in present vertex for output purposes
    real(rkind), dimension(:), allocatable       :: dkdx      ! derivative of wave number in x-direction
    real(rkind), dimension(:), allocatable       :: dkdy      ! derivative of wave number in y-direction
    real(rkind), dimension(:,:), allocatable     :: dmw       ! mud dissipation rate
    real(rkind), dimension(:,:), allocatable     :: dsnl      ! total interaction contribution of quadruplets to the main diagonal matrix
    real(rkind), dimension(:,:,:), allocatable   :: genc0     ! explicit part of generation in present vertex for output purposes
    real(rkind), dimension(:,:,:), allocatable   :: genc1     ! implicit part of generation in present vertex for output purposes
    real(rkind), dimension(:), allocatable       :: hscurr    ! wave height at current iteration level
    real(rkind), dimension(:), allocatable       :: hsdifc    ! difference in wave height of current and one before previous iteration
    real(rkind), dimension(:), allocatable       :: hsprev    ! wave height at previous iteration level
    real(rkind), dimension(:,:), allocatable     :: kwave     ! wave number
    real(rkind), dimension(:,:), allocatable     :: leakcf    ! leak coefficient in present vertex for output purposes
    real(rkind), dimension(:,:,:), allocatable   :: memnl4    ! auxiliary array to store results of 4 wave-wave interactions in full spectral space
    real(rkind), dimension(:,:,:), allocatable   :: obredf    ! action reduction coefficient based on transmission
    real(rkind), dimension(:,:,:), allocatable   :: redc0     ! explicit part of redistribution in present vertex for output purposes
    real(rkind), dimension(:,:,:), allocatable   :: redc1     ! implicit part of redistribution in present vertex for output purposes
    real(rkind), dimension(:,:), allocatable     :: reflso    ! contribution to the source term due to reflection
    real(rkind), dimension(:,:), allocatable     :: sa1       ! explicit interaction contribution of first quadruplet (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: sa2       ! explicit interaction contribution of second quadruplet (unfolded space)
    real(rkind), dimension(:,:), allocatable     :: sfnl      ! total interaction contribution of quadruplets to the right-hand side
    real(rkind), dimension(:,:,:,:), allocatable :: swtsda    ! several source terms computed at test points
!PUN    real(SZ), dimension(:), allocatable   :: temp      ! temporary array to store high precision data for UPDATER
    real(rkind), dimension(:), allocatable       :: tmcurr    ! mean period at current iteration level
    real(rkind), dimension(:), allocatable       :: tmdifc    ! difference in mean period of current and one before previous iteration
    real(rkind), dimension(:), allocatable       :: tmprev    ! mean period at previous iteration level
    real(rkind), dimension(:,:,:), allocatable   :: trac0     ! explicit part of propagation in present vertex for output purposes
    real(rkind), dimension(:,:,:), allocatable   :: trac1     ! implicit part of propagation in present vertex for output purposes
    real(rkind), dimension(:,:), allocatable     :: ue        ! energy density for computing 4 wave-wave interactions (unfolded space)
    !
    logical                               :: fguess    ! indicate whether first guess need to be applied or not
    logical                               :: lpredt    ! indicate whether action density in first iteration need to be estimated or not
    logical                               :: swpfull   ! indicate whether all necessary sweeps are done or not
    !
    logical, dimension(:,:), allocatable  :: anybin    ! true if bin is active in considered sweep
    logical, dimension(:,:), allocatable  :: anyblk    ! true if bin is blocked by a counter current based on a CFL criterion
    logical, dimension(:), allocatable    :: anywnd    ! true if wind input is active in considered bin
    logical, dimension(:,:), allocatable  :: groww     ! check for each frequency whether the waves are growing or not
                                                       ! in a spectral direction
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes

    integer ii,kk,jj
    !
!
!-----------------------------------------------------------------------
!                      End of variable definition
!-----------------------------------------------------------------------
!
    !
    ! print all the settings used in SWAN run
    !
     IF (LTRACE) CALL STRACE (IENT,'SWANCOMPUNSTRUCT')
    !PRINT*,'SWANCOMPUNSTRUCT 1'
    !PRINT *,'SPCDIR(:,2)', SPCDIR(:,2)
    !PRINT *,'SPCDIR(:,3)', SPCDIR(:,3)
    !PRINT *,'SPCDIR(1,2)', SPCDIR(1,2)
    !PRINT *,'SPCDIR(1,3)', SPCDIR(1,3)

    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! some initializations
    !
    ICMAX  = 3         ! stencil size
    PROPSL = PROPSC
    !
    IXCGRD(1) = -9999  ! to be used in routines SINTGRL and SOURCE so that Ursell number and
    IYCGRD(1) = -9999  ! quadruplets are calculated once in each vertex during an iteration
    !
    tid   = 0
    ivlow = 1
    ivup  = nverts
    !
    ! print all the settings used in SWAN run
    !
!    IF (IT .EQ. 1 .AND. ITEST.GE.1) CALL SWPRSET (SPCSIG)            ! 40.80
    IF( (myrank.eq.0) .and. (IT .EQ. 1) ) CALL SWPRSET (SPCSIG)
    !
    ! print test points
    !
!    if ( NPTST > 0 ) then
!       do j = 1, NPTST
!          write (PRINTF,107) j, xytst(j)
!       enddo
!    endif

!
    !
    ! allocation of shared arrays
    !WRITE(PRINTF,*) '! allocation of shared arrays'
    !
    allocate(islmin(nverts))
    allocate( nflim(nverts))
    allocate(nrscal(nverts))
    !
    allocate(hscurr(nverts))
    allocate(hsprev(nverts))
    allocate(hsdifc(nverts))
    allocate(tmcurr(nverts))
    allocate(tmprev(nverts))
    allocate(tmdifc(nverts))
    !
    ALLOCATE(SWTSDA(MDC,MSC,NPTSTA,MTSVAR))

!    *** quadruplets ***   !
    !WRITE(PRINTF,*) '! quadruplets'

    if ( IQUAD > 2 ) then
!      *** prior to every iteration full directional domain ***
       ALLOCATE(memnl4(MDC,MSC,nverts), stat = istat)
       if ( istat /= 0 ) then
          call msgerr ( 4, 'Allocation problem in SwanCompUnstruc: array memnl4 ' )
          return
       endif
    else
       ALLOCATE(memnl4(0,0,0))
    endif
    memnl4 = 0.
!
!  *** In case of SETUP expand array for setup data
!    IF(LSETUP > 0)THEN
!      MSTPDA = 23
!      ALLOCATE(SETPDA(nverts,MSTPDA))
!    ELSE
!      ALLOCATE(SETPDA(0,0))
!    END IF

!----------------------------------------------------------------------
!     Begin allocate shared arrays.
!----------------------------------------------------------------------
!
    !
    ! initialization of shared arrays
    !
    hscurr = 0.
    hsprev = 0.
    hsdifc = 0.
    tmcurr = 0.
    tmprev = 0.
    tmdifc = 0.
    !
    swtsda = 0.
    !
    tid = tid + 1

    !
    ! allocation of private arrays
    !
    allocate(   cad(MDC,MSC      ))
    allocate(   cas(MDC,MSC      ))
    allocate(   cax(MDC,MSC,ICMAX))
    allocate(   cay(MDC,MSC,ICMAX))
    allocate (  cgo(    MSC,ICMAX))
    allocate (kwave(    MSC,ICMAX))
    allocate (  dmw(    MSC,ICMAX))
    !
    allocate(idcmax(    MSC))
    allocate(idcmin(    MSC))
    allocate(iscmax(MDC    ))
    allocate(iscmin(MDC    ))
    allocate(anybin(MDC,MSC))
    !
    allocate(  amat(MDC,MSC,5))
    allocate(   rhs(MDC,MSC  ))
    allocate(ac2old(MDC,MSC  ))
    !
    allocate(anywnd(MDC))
    allocate(obredf(MDC,MSC,2))
    allocate(reflso(MDC,MSC))
    allocate( alimw(MDC,MSC))
    allocate( groww(MDC,MSC))
    allocate(anyblk(MDC,MSC))
    !
    allocate( disc0(MDC,MSC,MDISP))
    allocate( disc1(MDC,MSC,MDISP))
    allocate( genc0(MDC,MSC,MGENR))
    allocate( genc1(MDC,MSC,MGENR))
    allocate( redc0(MDC,MSC,MREDS))
    allocate( redc1(MDC,MSC,MREDS))
    allocate( trac0(MDC,MSC,MTRNP))
    allocate( trac1(MDC,MSC,MTRNP))
    allocate(leakcf(MDC,MSC      ))
    !
    allocate(dkdx(MSC))
    allocate(dkdy(MSC))

    !
    ! calculate ranges of spectral space for arrays related to 4 wave-wave interactions
    !
    !WRITE(PRINTF,*) '! ... to 4 wave-wave interactions'

    if ( IQUAD > 0 ) call FAC4WW ( xis, snlc1, dal1, dal2, dal3, spcsig, &
    &                              wwint, wwawg, wwswg )
    !
    if ( IQUAD > 0 ) then
       allocate(  ue(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate( sa1(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate( sa2(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       allocate(sfnl(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       if ( IQUAD == 1 ) then
          allocate(da1c(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da1p(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da1m(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2c(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2p(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(da2m(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
          allocate(dsnl(MSC4MI:MSC4MA,MDC4MI:MDC4MA))
       else
          allocate(da1c(0,0))
          allocate(da1p(0,0))
          allocate(da1m(0,0))
          allocate(da2c(0,0))
          allocate(da2p(0,0))
          allocate(da2m(0,0))
          allocate(dsnl(0,0))
       endif
    else
       allocate(  ue(0,0))
       allocate( sa1(0,0))
       allocate( sa2(0,0))
       allocate(sfnl(0,0))
       allocate(da1c(0,0))
       allocate(da1p(0,0))
       allocate(da1m(0,0))
       allocate(da2c(0,0))
       allocate(da2p(0,0))
       allocate(da2m(0,0))
       allocate(dsnl(0,0))
    endif
    !
    ! marks vertices active and non-active

    !WRITE(PRINTF,*)'SWANCOMPUNSTRUCT Mark Vertices'

    !
    nwetp = 0.
    do kvert = 1, nverts

       if ( COMPDA(kvert,JDP2) > DEPMIN ) then
          vert(kvert)%active = .true.
          nwetp = nwetp +1.
       else
          vert(kvert)%active = .false.
       endif
    enddo
    if ( it == 1 .and. ITEST > 0 ) write (PRINTF,108) nint(nwetp), nwetp*100./real(nverts)
    !
    ! First guess of action density will be applied if 3rd generation mode is employed 
    ! and wind is active (IWIND > 2)
    ! Note: this first guess is not used in nonstationary run (NSTATC > 0) or hotstart (ICOND = 4)
    !
    if ( IWIND > 2 .and. NSTATC == 0 .and. ICOND /= 4 ) then
       fguess = .true.
    else
       fguess = .false.
    endif
    !
    !
    !WRITE(PRINTF,*)'SWANCOMPUNSTRUCT ITER BEGIN,fguess ',fguess

    !IF(IT.EQ.1) THEN
    !  WRITE(PRINTF,*)'spcdir'
    !  WRITE(PRINTF,*)spcdir
    !  WRITE(PRINTF,*)'spcsig'
    !  WRITE(PRINTF,*)spcsig
    !ENDIF

    iterloop: do iter = 1, ITERMX

       !if(myrank.EQ.0) write (PRINTF,*) 'Iteration ',iter

       !
       ! some initializations
       !
       !
       !    *** IQUAD = 3: the nonlinear wave interactions are     ***
       !    *** calculated just once for an iteration. First,      ***
       !    *** set the auxiliary array equal zero before a        ***
       !    *** new iteration                                      *** 
       if ( IQUAD  > 2 ) memnl4(1:MDC,1:MSC,ivlow:ivup) = 0._rkind

       ! initialise Ursell number to 0 for each iteration
       if ( ITRIAD > 0 ) COMPDA(ivlow:ivup,JURSEL) = 0._rkind
       !
       ! initialise Dissipation and Leak to 0 for each iteration
       COMPDA(ivlow:ivup,JDISS) = 0._rkind
       COMPDA(ivlow:ivup,JLEAK) = 0._rkind
       ! and other stuff too
       COMPDA(ivlow:ivup,JDSXB) = 0._rkind
       COMPDA(ivlow:ivup,JDSXS) = 0._rkind
       COMPDA(ivlow:ivup,JDSXW) = 0._rkind
       COMPDA(ivlow:ivup,JDSXM) = 0._rkind
       COMPDA(ivlow:ivup,JDSXV) = 0._rkind
       COMPDA(ivlow:ivup,JDSXT) = 0._rkind
       COMPDA(ivlow:ivup,JGENR) = 0._rkind
       COMPDA(ivlow:ivup,JGSXW) = 0._rkind
       COMPDA(ivlow:ivup,JREDS) = 0._rkind
       COMPDA(ivlow:ivup,JRSXQ) = 0._rkind
       COMPDA(ivlow:ivup,JRSXT) = 0._rkind
       COMPDA(ivlow:ivup,JTRAN) = 0._rkind
       COMPDA(ivlow:ivup,JTSXG) = 0._rkind
       COMPDA(ivlow:ivup,JTSXT) = 0._rkind
       COMPDA(ivlow:ivup,JTSXS) = 0._rkind
       COMPDA(ivlow:ivup,JRADS) = 0._rkind
       COMPDA(ivlow:ivup,JQB  ) = 0._rkind

       ! Schism Coupling : Initialise SBR and SBF
       SBR = 0._rkind
       SBF = 0._rkind

       !
       ! initialise local (thread private) counter for SIP solver
       inocnt = 0
       !
       !    *** If a current is present and a penta-diagonal solver    ***
       !    *** is employed, it is possible that the solver does       ***
       !    *** not converged. For this, the counter INOCNV represents ***
       !    *** the number of geographical points in which the solver  ***
       !    *** did not converged    
       inocnv = 0
       !
       islmin = 9999
       nflim  = 0
       nrscal = 0
       !
       acnrms = -9999._rkind
       !
       ! During first iteration, first guess of action density is based on 2nd generation mode
       ! After first iteration, user-defined settings are re-used
       !
       if ( fguess ) then

          !
          if ( iter == 1 )then
             !
             ! save user-defined settings of 3rd generation mode
             ! Note: surf breaking, bottom friction and triads may be still active
             !
             usrset(1) = IWIND
             usrset(2) = IWCAP
             usrset(3) = IQUAD
             usrset(4) = PNUMS(20)
             usrset(5) = PNUMS(30)
             !
             ! first guess settings
             !
             IWIND     = 2        ! if first guess should be based on 1st generation mode, set IWIND = 1
             IWCAP     = 0
             IQUAD     = 0
             PNUMS(20) = 1.E22    ! no limiter
             PNUMS(30) = 0.       ! no under-relaxation
             !
!             write (PRINTF,101)
             !
          elseif ( iter == 2 ) then
             !
             ! re-activate user-defined settings of 3rd generation mode
             !
             IWIND     = usrset(1)
             IWCAP     = usrset(2)
             IQUAD     = usrset(3)
             PNUMS(20) = usrset(4)
             PNUMS(30) = usrset(5)
             !
!             write (PRINTF,102)
             !
          endif
          !
          ! JL, print info (DBG)
          !
!          if ( iter < 3 ) then
!             write (PRINTF,103) iter, PNUMS(20), PNUMS(30)
!             write (PRINTF,104) IWIND, IWCAP, IQUAD
!             write (PRINTF,105) ISURF, IBOT , ITRIAD
!             write (PRINTF,106) IVEG , ITURBV, IMUD
!          endif
          !
       endif
       !
       ! calculate diffraction parameter and its derivatives
       ! 
! JEROME : some additional job ToDo
!       if ( IDIFFR /= 0 ) call SWANDIFFPAR (ac2, COMPDA(1,JDP2), spcsig )
       if ( IDIFFR /= 0 ) call SWANDIFFPAR (    COMPDA(1,JDP2), spcsig )
       !
       ! all vertices are set untagged except non-active ones
       !
       do kvert = 1, nverts
          do jc = 1, vert(kvert)%noc         ! all cells around vertex
             vert(kvert)%updated(jc) = 0
          enddo
       enddo
       !

       do kvert = 1, nverts
          !
          ! in case of non-active vertex set action density equal to zero
          !
          if ( .not.vert(kvert)%active ) then
             !
             do jc = 1, vert(kvert)%noc
                icell = vert(kvert)%cell(jc)%atti(CELLID)
                vert(kvert)%updated(jc) = icell
             enddo
             !
             AC2(:,:,kvert) = 0._rkind
             !
          endif
          !
       enddo

!       if ( SCREEN /= PRINTF ) then
!          if ( NSTATC == 1 ) then
!             if ( MSR ) write (SCREEN,110) CHTIME, it, iter
!          else
!             write (PRINTF,120) iter
!             if ( MSR ) write (SCREEN,120) iter
!          endif
!       endif
       !
       ! loop over vertices in the grid
       !
       !WRITE(PRINTF,*)'loop over vertices in the grid'

!       print*,'First vert in domain ',myrank,vlist(ivlow),iplg(vlist(ivlow))
!       print*,'Last vert in domain ',myrank,vlist(ivup),iplg(vlist(ivup))

       vertloop: do kvert = ivlow, ivup
         !
         ivert = vlist(kvert)
!         print*,'ivert = ',ivert  
!         if(myrank.EQ.11) print*,'ivert=',iplg(ivert)
!         IF(ivert.EQ.78) CALL PSTOP

 !n        IF(kvert.GT.10) STOP
          
         !$ ivert = tlist(kvert)
!ADC         !
!ADC         ! allow SWAN to handle wave refraction as a nodal attribute
!ADC         if ( LoadSwanWaveRefrac .and. FoundSwanWaveRefrac ) then
!ADC            IREFR = nint(SwanWaveRefrac(ivert))
!ADC         endif
         !
         if ( vert(ivert)%active ) then   ! this active vertex needs to be updated
            !
            ! determine whether the present vertex is a test point
            !
            IPTST  = 0
            TESTFL = .false.
             if ( NPTST > 0 ) then
                do j = 1, NPTST
!                   if ( NGID(ivert) /= xytst(j) ) cycle
                   if ( ivert /= xytst(j) ) cycle
                   IPTST  = j
                   TESTFL = .true.
                enddo
             endif

            ! JL Print a resume
!            IF(TESTFL .AND. iter.EQ.3 .AND. IT.EQ.1 )THEN
!             DO j = 1, MNUMS
!              WRITE(PRINTF,'(A,I2,F10.2)')'PNUMS',j,PNUMS(j)
!             ENDDO
!             WRITE(PRINTF,'(A,I2)')'IWIND',IWIND
!             WRITE(PRINTF,'(A,F12.6)')'RDTIM',RDTIM
!            ENDIF
! JL Print a resume end

!            IF(TESTFL) THEN 
!              WRITE(PRINTF,*) 'TST Pt active,rk, depth',iplg(ivert),myrank,COMPDA(ivert,JDP2)
!            ENDIF
              

            !
            ! compute gradients of depth or wave number in present vertex meant for computing turning rate
            !
           call SWANGRADDEPTHORK ( COMPDA(1,JDP2), COMPDA(1,JMUDL2), &
                &            spcsig, dhdx, dhdy, dkdx, dkdy, ivert )

!           IF(TESTFL) THEN
!              WRITE(PRINTF,'(A,4I6,E12.4)')'rank,IT,iter,AC2 at ivert=',myrank,IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
!           ENDIF

!            WRITE(PRINTF,'(A,I8,4F10.6)')'AC2 at ivert=',iplg(ivert),AC2(1,1,ivert),AC2(MDC,1,ivert),AC2(1,MSC,ivert),AC2(MDC,MSC,ivert)
!            WRITE(PRINTF,'(A,I8,4F10.6)')'AC2 at ivert=',iplg(ivert),AC2(10,1,ivert),AC2(10,10,ivert),AC2(10,MSC,ivert),AC2(MDC,MSC,ivert)
!            WRITE(PRINTF,'(A,I8,4F10.6)')'AC2 at ivert=',iplg(ivert),AC2(15,1,ivert),AC2(15,15,ivert),AC2(15,MSC,ivert),AC2(MDC,MSC,ivert)
!            IF(IT.EQ.1.AND.iter.EQ.1) WRITE(PRINTF,'(A,I8,36E12.4)')'AC2 at ivert(ifrq=24)=',iplg(ivert),AC2(:,24,ivert)

!           IF(TESTFL.AND.IT.EQ.1.AND.ITER.EQ.1)THEN
!              print*,'rk,ivert,noc=',myrank,iplg(ivert),vert(ivert)%noc
!              print*,'loop over cells around considered vertex...'
!           ENDIF

!           IF(TESTFL) print*,''

            ! loop over cells around considered vertex
            celloop: do jc = 1, vert(ivert)%noc
               !
               icell = vert(ivert)%cell(jc)%atti(CELLID)
               !
               if ( vert(ivert)%updated(jc) == icell ) cycle celloop   ! this vertex in present cell is already updated
               !
               vv(1) = cell(icell)%atti(CELLV1)
               vv(2) = cell(icell)%atti(CELLV2)
               vv(3) = cell(icell)%atti(CELLV3)

               !
               ! pick up two upwave vertices
               !
               do k = 1, 3
                 if ( vv(k) == ivert ) then
                     vu(1) = vv(mod(k  ,3)+1)
                     vu(2) = vv(mod(k+1,3)+1)
                     exit
                 endif
               enddo

!               IF(TESTFL) THEN
!                print*,'rk,ivert, jc,icell=',myrank,iplg(ivert),jc,ielg(icell)
                !print*,'ivert, jc,icell (3vertex)=',iplg(ivert),jc,ielg(icell),iplg(vv(1)),iplg(vv(2)),iplg(vv(3))
                !print*,'ivert pick up two upwave vertices=',iplg(ivert),iplg(vu(1)),iplg(vu(2))
!               ENDIF
               !
               ! stores vertices of computational stencil
               !
               vs(1) = ivert
               vs(2) = vu(1)
               vs(3) = vu(2)

!               IF(TESTFL) THEN
!                 PRINT*,'rk,jc,STENCIL:',myrank,jc,iplg(vs(1)),iplg(vs(2)),iplg(vs(3))
!               ENDIF
!               PRINT *,'STENCIL =',vs(1),vs(2),vs(3)
!               PRINT *,'------------------------------'

!               PRINT*,'vertices of computational stencil',ivs(1),vs(2),vs(3)

               !
               KCGRD = vs    ! to be used in some original SWAN routines
               !
               swpnr = 0                                              ! this trick assures to calculate
                                                                      ! Ursell number and
               if ( all(mask=vert(ivert)%updated(:)==0) ) swpnr = 1   ! quadruplets only once in each vertex 
                                                                      ! during an iteration
               !
               ! compute wavenumber and group velocity in points of stencil
               !
               !PRINT*,'compute wavenumber and group velocity'

               call SWANDISPPARM ( kwave, cgo, dmw, COMPDA(1,JDP2), &
                  &                     COMPDA(1,JMUDL2), spcsig )
               !
               ! compute wave transport velocities in points of stencil for all directions
               !
               !PRINT*,'compute wave  transport velocities'

               !IF(TESTFL) THEN
               !  WRITE(PRINTF,*)'ICUR=',ICUR
               !  WRITE(PRINTF,*)ivert,'COMPDA(1,JVX2), COMPDA(1,JVY2)',COMPDA(ivert,JVX2), COMPDA(ivert,JVY2)
               !  WRITE(PRINTF,*)ivert,'COMPDA(1,JWX2), COMPDA(1,JWY2)',COMPDA(ivert,JWX2), COMPDA(ivert,JWY2)
               !ENDIF

               !IF(TESTFL) THEN
               ! WRITE(PRINTF,*)'CGO (min,mean,max) at:',iplg(ivert),'depth:',COMPDA(ivert,JDP2)
               !   do k = 1, ICMAX
               !      WRITE(PRINTF,*)minval(cgo(:,k)),sum(cgo(:,k))/MSC,maxval(cgo(:,k))
               !   enddo
               !ENDIF


               call SWANPROPVELX ( cax, cay, COMPDA(1,JVX2), COMPDA(1,JVY2), &
                  &                         cgo, spcdir(1,2), spcdir(1,3) )
               !
               ! compute local contravariant base vectors at present vertex
               !
               !PRINT*,'compute local contravariant base vectors'

               do k = 1, 3
                  if ( vv(k) == ivert ) then
                     rdx(1) = cell(icell)%geom(k)%rdx1
                     rdx(2) = cell(icell)%geom(k)%rdx2
                     rdy(1) = cell(icell)%geom(k)%rdy1
                     rdy(2) = cell(icell)%geom(k)%rdy2
                     exit
                  endif
               enddo

               !
               ! in case of spherical coordinates determine cosine of latitude (in degrees)
               ! and recalculate local contravariant base vectors
               !
               if(KSPHER.EQ.1) then
                  do k = 1, ICMAX
                     COSLAT(k) = cos(DEGRAD*(vert(vs(k))%attr(VERTY) + YOFFS))
                  enddo
                  do j = 1, 2
                     rdx(j) = rdx(j) / (COSLAT(1) * LENDEG)
                     rdy(j) = rdy(j) / LENDEG
                  enddo
               endif

               !IF(IT.EQ.1) print*,'IT,SWANPROPVELS',IT,LENDEG
                  !IF(TESTFL.AND.IT.EQ.1) THEN
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'COSLAT(1)',COSLAT(1)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'VERTY',vert(vs(1))%attr(VERTY)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'YOFFS',YOFFS
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'LENDEG',LENDEG
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'dep1',COMPDA(ivert,JDP1)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'dep2',COMPDA(ivert,JDP2)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'rdx1 =',rdx(1)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'rdx2 =',rdx(2)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'rdy1 =',rdy(1)
                  ! WRITE(PRINTF,*)'IT,SWANPROPVELS',IT,iplg(ivert),'rdy2 =',rdy(2)
                  !ENDIF



                 ! IF(TESTFL) THEN
                 !   WRITE(PRINTF,*)'SWANPROPVELS after LAT correc.',ivert,'rdx =',rdx,'rdy =',rdy
                 ! ENDIF

               !IF(TESTFL) THEN
               !WRITE(PRINTF,*)'spcsig',spcsig
               !WRITE(PRINTF,*)ivert,'cax =',cax
               !WRITE(PRINTF,*)ivert,'cay =',cay
               !ENDIF


               !
               ! compute spectral directions for the considered sweep in present vertex
               !
!               PRINT*,'compute spectral directions for the considered sweep'

               call SWANSWEEPSEL ( idcmin, idcmax, anybin, iscmin, iscmax, &
                                   iddlow, iddtop, idtot , isslow, isstop, &
                                   istot , cax   , cay   , rdx   , rdy   , &
                                   spcsig)


!                IF(TESTFL) THEN
!                  WRITE(PRINTF,'(A,3I5,A,2I5)')'Sweep',IT,iter,iplg(ivert),'iscmin, iscmax',iscmin, iscmax
!                  WRITE(PRINTF,'(A,3I5,A,2I5)')'Sweep',IT,iter,iplg(ivert),'isslow, isstop',isslow,isstop
!                  WRITE(PRINTF,'(A,3I5,A,2I5)')'Sweep',IT,iter,iplg(ivert),'dir bins',iddlow,iddtop
!                ENDIF

               !
               if ( idtot > 0 ) then
                  !
                  ! compute propagation velocities CAS and CAD in spectral space 
                  ! for the considered sweep in present vertex
                  !
!                  PRINT*,'compute propagation velocities in spectral space'

                  call SWANPROPVELS ( cad , cas , COMPDA(1,JVX2), COMPDA(1,JVY2), &
                                      COMPDA(1,JDP1) , COMPDA(1,JDP2), cax , cay, &
                                      kwave , cgo , spcsig , iddlow, iddtop     , &
                                      spcdir(1,2) , spcdir(1,3) , spcdir(1,4)   , &
                                      spcdir(1,5) , spcdir(1,6) , rdx , rdy     , &
                                      dhdx , dhdy , dkdx , dkdy )

                  !IF(TESTFL) THEN
                  ! WRITE(PRINTF,*)'SWANPROPVELS',IT,iter,iplg(ivert),'cad:'
                  ! WRITE(PRINTF,*) cad(1,1),cad(12,1),cad(MDC,1)
                  ! WRITE(PRINTF,*) cad(1,12),cad(12,12),cad(MDC,12)
                  ! WRITE(PRINTF,*) cad(1,MSC),cad(12,MSC),cad(MDC,MSC)
                  ! WRITE(PRINTF,*)'SWANPROPVELS',IT,iter,iplg(ivert),'cas:'
                  ! WRITE(PRINTF,*) cas(1,1),cas(12,1),cas(MDC,1)
                  ! WRITE(PRINTF,*) cas(1,12),cas(12,12),cas(MDC,12)
                  ! WRITE(PRINTF,*) cas(1,MSC),cas(12,MSC),cas(MDC,MSC)
                  !ENDIF

                  !
                  ! estimate action density in case of first iteration at cold start 
                  ! in stationary mode (since it is zero in first stationary run)
                  !
                  lpredt = .false.
                  if ( iter.EQ.1 .and. ICOND.NE.4 .and. NSTATC.EQ.0 ) then
                     lpredt = .true.
                     COMPDA(ivert,JHS) = 0.
                     goto 20
                  endif
 10               if ( lpredt ) then
!                    PRINT*,'STAT Mode : estimate action density SPREDT'

!                     call SPREDT (swpnr , ac2 , cax , cay , idcmin, idcmax,       &
                     call SPREDT (swpnr       , cax , cay , idcmin, idcmax,       &
                                  isstop, anybin, dummy, dummy, rdx , rdy , obredf)

                     lpredt = .false.
                  endif

!                  IF(TESTFL) THEN
!                      ! WRITE(PRINTF,*)'after source IMATRA',ivert,'rhs=',rhs(:,MSC)
!                       WRITE(PRINTF,*)'AC2 before SINTGRL ',ivert,' AC2=', maxval(AC2(:,:,ivert))
!                       WRITE(PRINTF,*)'AC2 before SINTGRL ',ivert,' AC2=', AC2(:,:,ivert)
!                  ENDIF

                  !
                  ! calculate various integral parameters for use in the source
                  ! terms
                  !
!                  call SINTGRL (spcdir , kwave , ac2 , COMPDA(1,JDP2), qbloc , &
                  call SINTGRL (spcdir , kwave ,            COMPDA(1,JDP2), qbloc ,  &
                                COMPDA(1,JURSEL),                                    &
                                rdx , rdy , dummy, etot , abrbot, COMPDA(1,JUBOT) ,  &
                                hs , COMPDA(1,JQB) , hm , kmespc , smebrk, kteta  ,  &
                                COMPDA(1,JPBOT) , COMPDA(1,JBOTLV), COMPDA(1,JGAMMA),&
                                COMPDA(1,JRESPL), swpnr ,                            &
                                iddlow , iddtop )

!                  IF(TESTFL .and. jc.EQ. vert(ivert)%noc) THEN
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'GAMBR',COMPDA(ivert,JGAMMA)
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'RESPL',COMPDA(ivert,JRESPL)
!                   !WRITE(PRINTF,*)'SINTGRL ',iplg(ivert),'QB, UBOT', COMPDA(ivert,JQB),COMPDA(ivert,JUBOT)
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'URSELL',COMPDA(ivert,JURSEL)
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'TMBOT',COMPDA(ivert,JPBOT)
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'ABRBOT', abrbot
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'ETOT',ETOT
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'HM',HM
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'HS',HS
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'QB_LOC',qbloc
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'AC2TOT',dummy
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'KMESPC',kmespc
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'SMEBRK',smebrk
!                   WRITE(PRINTF,'(A,3I7,A,E12.4)')'SINTGRL ',IT,iter,iplg(ivert),'KTETA',kteta
!                  ENDIF

                  !
                  !COMPDA(ivert,JHS) = hs

                  ! JL Add Fraction of breaking Wave in COMPDA (New) ??
                  ! COMPDA(ivert,JQB) = qbloc

 20               continue
                  !
                  ! compute transmission and/or reflection if obstacle is present 
                  ! in computational stencil
                  !
                  obredf = 1._rkind
                  reflso = 0._rkind
                  !

                  if ( NUMOBS > 0 ) then
                     !
                     ! determine obstacle for the link(s) in the computational stencil
                     !
                     do j = 1, cell(icell)%nof
                        !
                        ! determine vertices of the local face
                        !
                        vb = cell(icell)%face(j)%atti(FACEV1)
                        ve = cell(icell)%face(j)%atti(FACEV2)
                        !
                        if ( vb==ivert .or. ve==ivert ) then
                           !
                           if ( vb==vu(1) .or. ve==vu(1) ) then
                              !
                              iface = cell(icell)%face(j)%atti(FACEID)
                              link(1) = cross(iface)
                              !
                           elseif ( vb==vu(2) .or. ve==vu(2) ) then
                              !
                              iface = cell(icell)%face(j)%atti(FACEID)
                              link(2) = cross(iface)
                              !
                           endif
                           !
                        endif
                        !
                     enddo
                     !
                     if ( link(1)/=0 .or. link(2)/=0 ) then
                        !
                        call SWTRCF ( COMPDA(1,JWLV2), COMPDA(1,JHS), link, &
!                                      obredf, ac2, reflso, dummy , dummy  , &
                                      obredf,    reflso, dummy , dummy  , &
                                      dummy , cax , cay , rdx , rdy, anybin,&
                                      spcsig, spcdir )
                        !
                     endif
                     !
                  endif

!                  PRINT*,'swancompunstruc After HS OBSTACLE'


                  if (lpredt) goto 10
                  !
                  ! compute the transport part of the action balance equation
                  !
!                  PRINT*,'compute the transport part of the action balance eq.'

!                  IF(TESTFL) THEN
!                       WRITE(PRINTF,*)'IMATRA before SwanTranspAc',ivert,'rhs=',rhs(:,MSC)
!                       WRITE(PRINTF,*)'AC2 before SwanTranspAc ',ivert,' AC2=', maxval(AC2(:,:,ivert))
!!                       WRITE(PRINTF,*)'AC2 before SwanTranspAc ',ivert,' AC2=', AC2(:,:,ivert)
!                  ENDIF

!                  call SwanTranspAc ( amat  , rhs   , leakcf, ac2   , ac1   , &
                  call SwanTranspAc ( amat  , rhs   , leakcf,                 &
                                      cgo   , cax   , cay   , cad   , cas   , &
                                      anybin, rdx   , rdy   , spcsig, spcdir, &
                                      obredf, idcmin, idcmax, iscmin, iscmax, &
                                      iddlow, iddtop, isslow, isstop, anyblk, &
                                      trac0 , trac1 )

!                   DO ii=1,MSC
!                    DO jj=1,MDC 
!                     DO kk=1,5
!                     IF(amat(jj,ii,kk)/=amat(jj,ii,kk))THEN
!                      WRITE(PRINTF,*)'NAN in amat ID,IS,K,ivert',jj,ii,kk,ivert
!                      call parallel_abort('stop')
!                     ENDIF
!                     ENDDO
!                    ENDDO
!                   ENDDO

!                   IF(TESTFL) THEN
!                    WRITE(PRINTF,*)'IMATRA after SwanTranspAc',ivert,'rhs',rhs(:,MSC)
!                    WRITE(PRINTF,*)'AC2 after SwanTranspAc ',ivert,' AC2=', maxval(AC2(:,:,ivert))
!                   ENDIF

                  !
                  ! compute the source part of the action balance equation
                  !
!DBG
!DBG              OFFSRC = .TRUE.   ! NO SOURCE
!DBG              PNUMS(8) = 2.0    ! SOLMT1
!DBG              DYNDEP = .FALSE.  ! SOLMAT
!DBG              ICUR = 0          ! SOLMAT

                  if ( .not.OFFSRC ) then
                     !
                     ! initialize wind friction velocity and Pierson Moskowitz frequency
                     !
                     ufric = 1.e-15
                     fpm   = 1.e-15
                     !
                     ! compute the wind speed, mean wind direction, the PM frequency,
                     ! wind friction velocity and the minimum and maximum counters for
                     ! active wind input
                     !
                     if ( IWIND > 0 ) &
                     call WINDP1 ( wind10, thetaw, idwmin , idwmax , fpm, ufric,   &
                                   COMPDA(1,JWX2), COMPDA(1,JWY2), anywnd, spcdir, &
!                                  COMPDA(1,JVX2), COMPDA(1,JVY2), spcsig, ac2     &
                                   COMPDA(1,JVX2), COMPDA(1,JVY2), spcsig )
!ADC                               ,ivert )
                     !
                     ! compute the source terms
                     !
!                     PRINT *,'compute the source terms'

!        rhs = IMATRA  2D    Coefficients of right hand side of matrix         
!amat(:,:,1) = IMATDA  2D    Coefficients of diagonal of matrix
!amat(:,:,2) = IMAT5L  2D    Coefficients for implicit calculation in frequency space (lower diagonal)
!amat(:,:,3) = IMAT6U  2D    Coefficients for implicit calculation in frequency space (upper diagonal)
!amat(:,:,4) = IMATLA  2D    Coefficients of lower diagonal of matrix
!amat(:,:,5) = IMATUA  2D    Coefficients of upper diagonal of matrix

!                      IF(TESTFL) THEN
!                       WRITE(PRINTF,*)'before source',ivert,' rhs=',rhs(:,MSC)
!                       !WRITE(PRINTF,*)'before source',ivert,' amat(:,1:5,1)',amat(:,1:5,1)
!                      ENDIF

                     !IF(TESTFL) WRITE(PRINTF,*)'call SOURCE'


!JL: The Arguments ordering is very messy in 'swan41.10' ... Some changes are
!    applied
                     call SOURCE ( iter , IXCGRD(1) , IYCGRD(1) , swpnr , kwave, &
!                                  spcsig , spcdir(1,2) , spcdir(1,3), ac2 ,     &
                                   spcsig , spcdir,                              &
                                   COMPDA(1,JDP2) , amat(1,1,1) , rhs ,abrbot ,  &
                                   kmespc , dummy , COMPDA(1,JUBOT), ufric ,     &
                                   COMPDA(1,JVX2) , COMPDA(1,JVY2) , idcmin ,    & 
                                   idcmax , iddlow , iddtop , idwmin , idwmax ,  &
                                   isstop , hs , etot , qbloc ,                  &
                                   thetaw , hm , fpm , wind10 , dummy , groww,   &
                                   alimw , smebrk , kteta , snlc1 , dal1, dal2,  &
                                   dal3 , ue , sa1 , sa2 , da1c , da1p , da1m ,  &
                                   da2c , da2p , da2m , sfnl , dsnl , memnl4 ,   &
                                   wwint , wwawg , wwswg , cgo, COMPDA(1,JUSTAR),&
                                   COMPDA(1,JCHARN),&  ! JL add CHarnock
!                                   COMPDA(1,JZEL) , spcdir , anywnd , dmw ,disc0,&     
                                   COMPDA(1,JZEL), anywnd , dmw ,disc0,          &
                                   disc1 , genc0 , genc1 , redc0 , redc1 , xis,  &
                                   COMPDA(1,JFRC2) , it , COMPDA(1,JNPLA2) ,     &
                                   COMPDA(1,JTURB2) , COMPDA(1,JMUDL2),          &
                                   COMPDA(1,JURSEL) , anybin ,reflso ,           &
                                   swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB) ,  &
                                   swtsda(1,1,1,JPWCAP), swtsda(1,1,1,JPBTFR),   &
                                   swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S),     &
                                   swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPVEGT),   &
                                   swtsda(1,1,1,JPTURB), swtsda(1,1,1,JPMUD) ,   &
                                   swtsda(1,1,1,JPTRI) )

!                   IF(TESTFL) THEN
!                    WRITE(PRINTF,'(A,3I5,E12.4)')'RHS after SOURCE',IT,iter,iplg(ivert),SUM(rhs)
!                    WRITE(PRINTF,'(A,3I5,E12.4)')'AC2 after SOURCE',IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
!                   ENDIF


!                  DO ii=1,MSC
!                    DO jj=1,MDC
!                     DO kk=1,5
!                     IF(amat(jj,ii,kk)/=amat(jj,ii,kk))THEN
!                      WRITE(PRINTF,*)'after SOURCE NAN in amat ID,IS,K,ivert',jj,ii,kk,ivert
!                      call parallel_abort('stop')
!                     ENDIF
!                     ENDDO
!                    ENDDO
!                   ENDDO


                     !
                  endif
                  !
                  ! update action density by means of solving the action balance equation
                  !
                  !IF(TESTFL) THEN
                  ! WRITE(PRINTF,*)'update action density using action balance equation' 
                  ! WRITE(PRINTF,*)'AC2 after SOURCE ',ivert,' AC2=', maxval(AC2(:,:,ivert))
                  !ENDIF
!                  CALL PSTOP


                  if ( ACUPDA ) then
                     !
                     ! preparatory steps before solving system of equations
                     !

                   ! IF(TESTFL) WRITE(PRINTF,*)'call SOLPRE'

!                     PRINT *,'call SOLPRE ac2old',ac2old

!                     call SOLPRE(ac2        , ac2old     , rhs        , amat(1,1,4), &

                    call SOLPRE(              ac2old     , rhs        , amat(1,1,4), &
                                 amat(1,1,1), amat(1,1,5), amat(1,1,2), amat(1,1,3), &
                                 idcmin     , idcmax     , anybin     , idtot      , &
                                 istot      , iddlow     , iddtop     , isstop     , &
                                 spcsig     )

!                       IF(TESTFL) THEN
!                          WRITE(PRINTF,*)'IMATRA(rhs) after SOLPRE',iplg(ivert),'rhs=',sum(rhs(:,:))
!                          WRITE(PRINTF,*)'after solpre ',iplg(ivert),' AC2=',sum(AC2(:,:,ivert))
!                       ENDIF


!                  DO ii=1,MSC
!                    DO jj=1,MDC
!                     DO kk=1,5
!                     IF(amat(jj,ii,kk)/=amat(jj,ii,kk))THEN
!                      WRITE(PRINTF,*)'after SOLPRE NAN in amat ID,IS,K,ivert',jj,ii,kk,ivert
!                      call parallel_abort('stop')
!                     ENDIF
!                     ENDDO
!                    ENDDO
!                   ENDDO



                     !
                     if ( .not.DYNDEP .and. ICUR == 0 ) then
                        !
                        ! propagation in theta space only
                        ! solve tridiagonal system of equations using Thomas' algorithm
                        !
!                        PRINT *,'call SOLMAT rhs',rhs

!rhs = IMATRA
!amat(1,1,4) = IMATLA
!amat(1,1,1) = IMATDA
!amat(1,1,5) = IMATUA
!amat(1,1,2) = IMAT5L
!amat(1,1,3) = IMAT6U

                       ! IF(TESTFL) WRITE(PRINTF,*)'call SOLMAT'

!                        call SOLMAT ( idcmin     , idcmax     , ac2        , rhs, &
                        call SOLMAT ( idcmin     , idcmax     ,               rhs, &
                                      amat(1,1,1), amat(1,1,5), amat(1,1,4)     )

                        ! IF(TESTFL) THEN
                        !  WRITE(PRINTF,*)'after solmat IMATRA',ivert,'rhs=',rhs(:,MSC)
                        !  WRITE(PRINTF,*)'after solmat ',ivert,' AC2=',maxval(AC2(:,:,ivert))
                        ! ENDIF


!                        CALL PSTOP
                        !
                     else
                        !
                        ! propagation in both sigma and theta spaces
                        !
                        if ( int(PNUMS(8)) == 1 ) then
                           !
                           ! implicit scheme in sigma space
                           ! solve pentadiagonal system of equations using SIP solver
                           !
!                           PRINT *,'call SWSIP'
                          ! IF(TESTFL) WRITE(PRINTF,*)'call SWSIP'

!                          call SWSIP ( ac2        , amat(1,1,1)    , rhs            , amat(1,1,4), &
                           call SWSIP (              amat(1,1,1)    , rhs            , amat(1,1,4), &
                                        amat(1,1,5), amat(1,1,2)    , amat(1,1,3)    , ac2old     , &
                                        PNUMS(12)  , nint(PNUMS(14)), nint(PNUMS(13)), inocnt     , &
                                        iddlow     , iddtop         , isstop         , idcmin     , &
                                        idcmax     )

                           !
                        elseif (int(PNUMS(8)) == 2 ) then
                           !
                           ! explicit scheme in sigma space
                           ! solve tridiagonal system of equations using Thomas' algorithm
                           !
!                           PRINT *,'call SOLMT1'
                           !IF(TESTFL) WRITE(PRINTF,*)'call SOLMT1'
!                           call SOLMT1  ( idcmin , idcmax , ac2 , rhs    ,       &
                           call SOLMT1  ( idcmin , idcmax ,       rhs    ,       &
                                          amat(1,1,1), amat(1,1,5), amat(1,1,4), &
                                          isstop , anyblk , iddlow , iddtop )
                           !
                        endif
                        !
                     endif
                     !
                     ! if negative action density occur rescale with a factor
                     !
!                     PRINT *,'rescale with a factor'

                     if ( BRESCL ) &
!                     call RESCALE ( ac2, isstop, idcmin, idcmax, nrscal)
                     call RESCALE (       isstop, idcmin, idcmax, nrscal)

                     !
                     ! store propagation, generation, dissipation, redistribution, leak and radiation 
                     ! stress in present vertex
                     !
                     if ( LADDS )                                                     &
                     call ADDDIS ( COMPDA(1,JDISS), COMPDA(1,JLEAK),       anybin ,   &
                                   disc0 , disc1 , genc0 , genc1 , redc0 , redc1 ,    &
                                   trac0 , trac1 , amat(1,1,4) , amat(1,1,5) ,        &
                                   amat(1,1,2) , amat(1,1,3) ,                        &
                                   COMPDA(1,JDSXB), COMPDA(1,JDSXS), COMPDA(1,JDSXW), &
                                   COMPDA(1,JDSXV), COMPDA(1,JDSXT), COMPDA(1,JDSXM), &
                                   COMPDA(1,JGSXW), COMPDA(1,JGENR),                  &
                                   COMPDA(1,JRSXQ), COMPDA(1,JRSXT),                  &
                                   COMPDA(1,JREDS),                                   &
                                   COMPDA(1,JTSXG), COMPDA(1,JTSXT),                  &
                                   COMPDA(1,JTSXS), COMPDA(1,JTRAN),                  &
                                   leakcf         , COMPDA(1,JRADS), spcsig           )

                     !print*,'Bot, Cap.',COMPDA(ivert,JDSXB),COMPDA(1,JDSXW)
                     !
                     ! limit the change of the spectrum
                     !

                     if ( PNUMS(20) < 100. ) &
!                     call PHILIM ( ac2, ac2old, cgo, kwave, spcsig, anybin, &
                      call PHILIM (       ac2old, cgo, kwave, spcsig, anybin, &
                                   islmin, nflim, qbloc )
!                      IF(TESTFL) THEN
!                       WRITE(PRINTF,*)'kwave',kwave(1,1),kwave(12,1),kwave(MSC,1) 
!                       WRITE(PRINTF,*)'cgo',cgo(1,1),cgo(12,1),cgo(MSC,1)
!                      ENDIF
                     !
                     ! reduce the computed energy density if the value is larger then the limit value
                     ! as computed in SOURCE in case of first or second generation mode
                     !
                     if ( IWIND == 1 .or. IWIND == 2 ) &
!                     call WINDP3 ( isstop, alimw, ac2, groww, idcmin, idcmax )
                     call WINDP3 ( isstop, alimw,      groww, idcmin, idcmax )

                     !
                     ! store some infinity norms meant for convergence check
                     !
                     if ( PNUMS(21) == 2. ) &
!                     call SWACC ( ac2, ac2old, acnrms, isstop, idcmin, idcmax )
                     call SWACC (      ac2old, acnrms, isstop, idcmin, idcmax )
                     !

                  endif
                  !
               endif
               !
               ! tag updated vertex in present cell
               !
               vert(ivert)%updated(jc) = icell
               !

            enddo celloop

!            PRINT*,'End Cell Loop'

!ADC            !
!ADC            ! make boundaries reflective
!ADC            if ( .false. ) then
!ADC               !
!ADC               if ( vert(ivert)%atti(VMARKER)==1 .and. vert(ivert)%atti(VBC)==0 ) call MakeBoundariesReflective( ivert, ac2 )
!ADC               !
!ADC            endif
            !
         endif
         !
       enddo vertloop

      !CALL parallel_abort('STOP vertloop')

       !WRITE(PRINTF,*)'End Vert Loop'
       !PRINT*,'End Vert Loop'

       !
       ! global sum to inocnv counter
       !
       inocnv = inocnv + inocnt

       !
       ! exchange action densities with neighbours in parallel run
       ! # if defined (MULTIPROCESSOR)
!DBG
#if 0
       call parallel_barrier

       if ( NPTST > 0 ) then

           do kvert = 1, npa
            ivert = vlist(kvert)
            do j = 1, NPTST
              if ( ivert .EQ. xytst(j) ) then
               !WRITE(PRINTF,'(A,4I6,E12.4)')'rk,IT,iter,AC2 at ivert=',myrank,IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
               print'(A,4I6,E12.4)','bfr rk,IT,iter,AC2 at ivert=',myrank,IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
              endif
            enddo
           enddo
       endif
#endif
!DBG       

      IF(.NOT.SERIAL) THEN 
       CALL EXCHANGE_P4D_WWM(AC2)
      ENDIF

!DBG
#if 0
      if ( NPTST > 0 ) then
           do kvert = 1, npa
            ivert = vlist(kvert)
            do j = 1, NPTST
              if ( ivert .EQ. xytst(j) ) then
               !WRITE(PRINTF,'(A,4I6,E12.4)')'rk,IT,iter,AC2 at ivert=',myrank,IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
               print'(A,4I6,E12.4)','rk,IT,iter,AC2 at ivert=',myrank,IT,iter,iplg(ivert),SUM(AC2(:,:,ivert))
              endif
            enddo
           enddo

       endif
#endif
!DBG

       ! store the source terms assembled in test points per iteration in the files IFPAR, IFS1D and IFS2D
       !
!JL Skip that
# if 0
       if ( NPTST > 0 .and. NSTATM == 0 ) then
          if ( IFPAR > 0 ) write (IFPAR,151) iter
          if ( IFS1D > 0 ) write (IFS1D,151) iter
          if ( IFS2D > 0 ) write (IFS2D,151) iter
          call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), &
                        swtsda(1,1,1,JPBTFR), swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , &
                        swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , swtsda(1,1,1,JPVEGT), &
                        swtsda(1,1,1,JPTURB), swtsda(1,1,1,JPMUD) ,                       &
                        ac2 , spcsig , COMPDA(1,JDP2) , xytst , dummy )
       endif
# endif

       !
       ! Grab and Print Quantities useful for Diag / Debug purposes
       !
       if ( ITEST.GE.30 .or. idebug.EQ.1 ) then
          !
          SBUF = 0.0_dkind ; RBUF1 = 0.0_dkind

          ! Get the Total Number of active points
          nwetp = 0.
          do ivert = 1, nverts
             if ( vert(ivert)%active ) nwetp = nwetp +1.
          enddo

          RBUF1(1) = DBLE(nwetp)
          SBUF(1)  = DBLE(nwetp)
          !
          ! indicate number of vertices in which rescaling has taken place
          !
          npfr  = count(mask=nrscal>0)
          RBUF1(2) = DBLE(npfr )
          SBUF(2)  = DBLE(npfr )
          ! 
          ! indicate number of vertices in which limiter has taken place
          !
          npfl  = count(mask=nflim>0)
          RBUF1(3) = DBLE(npfl )
          SBUF(3)  = DBLE(npfl )
          !
          ! indicate number of vertices in which the SIP solver did not converged
          !
          RBUF1(4) = DBLE(inocnv) 
          SBUF(4)  = DBLE(inocnv)
          !
          ! Sum over partition
          !
!         # if defined (MULTIPROCESSOR)
          IF(.NOT.SERIAL) CALL MPI_REDUCE(SBUF,RBUF1,4,MPI_DOUBLE,MPI_SUM,0,comm,ierr)

          nwetp = RBUF1(1)
          npfr  = NINT(RBUF1(2))
          npfl  = NINT(RBUF1(3))
          inocnv = NINT(RBUF1(4))

          ! Get Quantities : 
          !  - Max number of frequency use of rescaling in each (nrscal)
          !  - Max number of frequency use of limiter in each vertex (nflim)
          SBUF = 0.0_dkind ; RBUF1 = 0.0_dkind
          ! 
          mxnfr = maxval(nrscal)
          mxnfl = maxval(nflim)
          RBUF1(1) = DBLE(mxnfr)
          SBUF(1)  = DBLE(mxnfr)
          RBUF1(2) = DBLE(mxnfl)
          SBUF(2)  = DBLE(mxnfl)
          !
          ! Get Max over partition
          !
!         # if defined (MULTIPROCESSOR)
          IF(.NOT.SERIAL) CALL MPI_REDUCE(SBUF,RBUF1,4,MPI_DOUBLE,MPI_MAX,0,comm,ierr)

          mxnfr = NINT(RBUF1(1))
          mxnfl = NINT(RBUF1(2))
          !
          ! indicate number of vertices in which rescaling has taken place
          frac = DBLE(npfr)*100./nwetp
          IF(myrank.EQ.0 .AND. npfr.GT.0) WRITE (PRINTF,130) 'rescaling', frac, mxnfr
!          IF(myrank.EQ.0 .AND. NSTATC.EQ.0 .AND. npfr.GT.0) WRITE(SCREEN,130) 'rescaling', frac, mxnfr
          !
          ! indicate number of vertices in which limiter has taken place
          frac = DBLE(npfl)*100./nwetp
          IF(myrank.EQ.0 .AND. npfl.GT.0) WRITE (PRINTF,130) 'limiter', frac, mxnfl
!          IF(myrank.EQ.0 .AND. NSTATC.EQ.0 .AND. npfl.GT.0) WRITE(SCREEN,130) 'limiter', frac, mxnfl
          !
          SBUF = 0.0_dkind ; RBUF1 = 0.0_dkind
          
          mnisl = minval(islmin)
          RBUF1(1) = DBLE(mnisl)
          SBUF(1)  = DBLE(mnisl)
!         # if defined (MULTIPROCESSOR)
          IF(.NOT.SERIAL) CALL MPI_REDUCE(SBUF,RBUF1,4,MPI_DOUBLE,MPI_MIN,0,comm,ierr)

          mnisl = NINT(RBUF1(1))
          !
! JL MUTE
!          IF(myrank.EQ.0 .AND. npfl.GT.0) WRITE(PRINTF,135) spcsig(mnisl)/PI2_W
!          IF(myrank.EQ.0 .AND. NSTATC.EQ.0 .AND. npfl.GT.0) WRITE(SCREEN,135) spcsig(mnisl)/PI2_W
          !
          ! indicate number of vertices in which the SIP solver did not converged
          IF(myrank.EQ.0 .AND. ((DYNDEP .or. ICUR.NE.0) .and. inocnv.NE.0)) WRITE(PRINTF,136) inocnv

       endif

       !
       ! exchange array COMPDA with neighbours in parallel run
       ! 
!      #      if defined (MULTIPROCESSOR)
       IF(.NOT.SERIAL) THEN
 
        IF(.NOT.ALLOCATED(COMPDA_TMP1)) ALLOCATE(COMPDA_TMP1(1:npa));   
        COMPDA_TMP1 = 0._rkind
        DO J = 1, MCMVAR
            COMPDA_TMP1(1:npa)=COMPDA(1:npa,J)
            CALL exchange_p2d(COMPDA_TMP1)
            COMPDA(1:npa,J) = COMPDA_TMP1(1:npa)
        ENDDO
        DEALLOCATE(COMPDA_TMP1)
 
       END IF

!#ifdef schism
       ! 
       ! exchange array sbr with neighbours in parallel run
       ! exchange array sbf with neighbours in parallel run
!      #      if defined (MULTIPROCESSOR)
       IF(.NOT.SERIAL) THEN

        IF(.NOT.ALLOCATED(COMPDA_TMP1)) ALLOCATE(COMPDA_TMP1(1:npa));

        COMPDA_TMP1 = 0._rkind
        DO J = 1, 2
            COMPDA_TMP1(1:npa) = sbr(J,1:npa)
            CALL exchange_p2d(COMPDA_TMP1)
            sbr(J,1:npa)  = COMPDA_TMP1(1:npa)

            COMPDA_TMP1(1:npa) = sbf(J,1:npa)
            CALL exchange_p2d(COMPDA_TMP1)
            sbf(J,1:npa)  = COMPDA_TMP1(1:npa)

        ENDDO

        DEALLOCATE(COMPDA_TMP1)

       END IF

!endif

       ! Computing Accuracy (SwanConvAccur) and check stop criteria (SwanConvStopc): 
       !
       ! info regarding the iteration process and the accuracy
       !
       if ( PNUMS(21) <= 1. ) then
          !
          if ( PNUMS(21) == 0. ) then
             !

             call SwanConvAccur ( accur, hscurr, tmcurr, COMPDA(1,JDHS), COMPDA(1,JDTM), &
!                                  xytst, spcsig, ac2, ivlow, ivup )
                                  xytst, spcsig,       ivlow, ivup )
             !
          else
             !
             call SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, tmdifc, &
!                                  COMPDA(1,JDHS), COMPDA(1,JDTM), xytst, spcsig, ac2, ivlow, ivup )
                                  COMPDA(1,JDHS), COMPDA(1,JDTM), xytst, spcsig,      ivlow, ivup )
             !
          endif

!         #if defined (MULTIPROCESSOR)
          IF(.NOT.SERIAL) CALL MPI_BCAST(accur,1,rtype,0,comm,ierr)
          

          if ( iter == 1 ) then
             accur = -9999._rkind
             if(myrank.EQ.0) write (PRINTF,142)
!             if ( NSTATC == 0 .and. MSR ) write (SCREEN,142)
          else
             if(myrank.EQ.0) write (PRINTF,140) accur, PNUMS(4)
!             if ( NSTATC == 0 .and. myrank.EQ.0 ) write (SCREEN,140) accur, PNUMS(4)
          endif

          !
          ! if accuracy has been reached then terminates iteration process
          !
          if ( accur >= PNUMS(4) ) exit iterloop
          !
       elseif ( PNUMS(21) == 2. ) then
          !
          write (PRINTF,141)
!          if ( NSTATC == 0 .and. myrank.EQ.0) write (SCREEN,141)
          if ( iter == 1 ) then
             stopcr = -9999.
             if(myrank.EQ.0) write (PRINTF,142)
!             if ( NSTATC == 0 .and. MSR ) write (SCREEN,142)
          else
             if ( acnrms(1) < 1.e-20 .or. acnrmo < 1.e-20 ) then
                if(myrank.EQ.0) write (PRINTF,143)
!                if ( NSTATC == 0 .and. myrank.EQ.0) write (SCREEN,143)
             else
                rhof   = acnrms(1)/acnrmo
                stopcr = PNUMS(1)*acnrms(2)*(1.-rhof)/rhof
                if(myrank.EQ.0) write (PRINTF,144) rhof, acnrms(1), stopcr
!                if ( NSTATC == 0 .and. myrank.EQ.0) write (SCREEN,144) rhof, acnrms(1), stopcr
             endif
          endif
          acnrmo = acnrms(1)

          !
          ! if accuracy has been reached then terminates iteration process
          !
          if ( acnrms(1) < stopcr ) exit iterloop
          !
       endif
       !
    enddo iterloop
    !
    ! deallocation of private arrays
    !
    deallocate(  cad)
    deallocate(  cas)
    deallocate(  cax)
    deallocate(  cay)
    deallocate(  cgo)
    deallocate(kwave)
    deallocate(  dmw)
    !
    deallocate(idcmax)
    deallocate(idcmin)
    deallocate(iscmax)
    deallocate(iscmin)
    deallocate(anybin)
    !
    deallocate(  amat)
    deallocate(   rhs)
    deallocate(ac2old)
    !
    deallocate(anywnd)
    deallocate(obredf)
    deallocate(reflso)
    deallocate( alimw)
    deallocate( groww)
    deallocate(anyblk)
    !
    deallocate( disc0)
    deallocate( disc1)
    deallocate( genc0)
    deallocate( genc1)
    deallocate( redc0)
    deallocate( redc1)
    deallocate( trac0)
    deallocate( trac1)
    deallocate(leakcf)
    !
    deallocate(dkdx)
    deallocate(dkdy)
    !
    deallocate(  ue)
    deallocate( sa1)
    deallocate( sa2)
    deallocate(sfnl)
    deallocate(da1c)
    deallocate(da1p)
    deallocate(da1m)
    deallocate(da2c)
    deallocate(da2p)
    deallocate(da2m)
    deallocate(dsnl)
    !
    ! end of parallel region
    !
    ! store the source terms assembled in test points per time step in the files IFPAR, IFS1D and IFS2D
    !
!JL SKip now
# if 0
    if ( NPTST > 0 .and. NSTATM == 1 ) then
       if ( IFPAR > 0 ) write (IFPAR,152) CHTIME
       if ( IFS1D > 0 ) write (IFS1D,152) CHTIME
       if ( IFS2D > 0 ) write (IFS2D,152) CHTIME
       call PLTSRC ( swtsda(1,1,1,JPWNDA), swtsda(1,1,1,JPWNDB), swtsda(1,1,1,JPWCAP), &
                     swtsda(1,1,1,JPBTFR), swtsda(1,1,1,JPWBRK), swtsda(1,1,1,JP4S)  , &
                     swtsda(1,1,1,JP4D)  , swtsda(1,1,1,JPTRI) , swtsda(1,1,1,JPVEGT), &
                     swtsda(1,1,1,JPTURB), swtsda(1,1,1,JPMUD) ,                       &
                     ac2 , spcsig , COMPDA(1,JDP2) , xytst , dummy               )
    endif
# endif
    !
    ! deallocation of shared arrays
    !
    deallocate(islmin)
    deallocate( nflim)
    deallocate(nrscal)
    !
    deallocate(hscurr)
    deallocate(hsprev)
    deallocate(hsdifc)
    deallocate(tmcurr)
    deallocate(tmprev)
    deallocate(tmdifc)
    !
    deallocate(swtsda)
    !
    deallocate(memnl4)
    !
    !
 101 format (// ' Settings of 2nd generation mode as first guess are used:')
 102 format (// ' User-defined settings of 3rd generation mode is re-used:')
 103 format (' ITER  ',i4,'    GRWMX ',e12.4,'    ALFA   ',e12.4)
 104 format (' IWIND ',i4,'    IWCAP ',i4   ,'            IQUAD  ',i4)
 105 format (' ISURF ',i4,'    IBOT  ',i4   ,'            ITRIAD ',i4)
 106 format (' IVEG  ',i4,'    ITURBV',i4   ,'            IMUD   ',i4,/)
 107 format (' Test points :',2i5)
 108 format (' Number of active points = ',i6,' (fillings-degree: ',f6.2,' %)')
 110 format ('+time ', a18, ', step ',i6, '; iteration ',i4)
 120 format (' iteration ', i4)
 130 format (1x,'use of ',a9,' in ',f6.2,' % of active vertices with maximum in spectral space = ',i4)
 135 format (1x,'lowest frequency occured above which limiter is applied = ',f7.4,' Hz')
 136 format (2x,'SIP solver: no convergence in ',i8,' vertices')
 140 format (' accuracy OK in ',f6.2,' % of wet grid points (',f6.2,' % required)',/ )
 141 format (7x,'Rho            ','Maxnorm   ','  Stoppingcrit.')
 142 format ('       not possible to compute accuracy, first iteration',/)
 143 format ('       norm less then 1e-20, no stopping criterion',/)
 144 format (1x,ss,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,/)
 151 format (i4, t41, 'iteration')
 152 format (a , t41, 'date-time')
    !
END SUBROUTINE SwanCompUnstruc
