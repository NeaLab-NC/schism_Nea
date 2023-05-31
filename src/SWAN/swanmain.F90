# if defined (USE_SWAN)
#include "swan_functions.h"
!***********************************************************************
!                                                                      *
!***********************************************************************
!
!    SWAN Main program
!
!***********************************************************************
!
!
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_SETUP(DT_SCHISM0,NSTEP_SWAN0)                                                   
!                                                                      *
!***********************************************************************
!
!     SWMAIN subroutine, calling SWINIT, SWREAD, SWCOMP 
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     Call SWINIT to initialize various common data
!                                                                   0
!         Call SWREAD to read and process user commands.
!              SWREAD populate several SWAN Unstruct. mesh objects (see
!              swangriddata,swangridobjects) by using informations from
!              the FVCOM mesh topology
!         -------------------------------------------------------------
!         Call SWPREP to check input and prepare computation
!         Call SWRBC to update COMPDA with wind, depth, current ... fields

!         If nonstationary computation is to be made                      40.00
!         Then start time step loop at IT=0 and                           40.00
!         Call SWINCO to calculate initial wave spectra                   40.00
!         -------------------------------------------------------------
!         Start Computation :
!               Call SWMAIN_LOOP from internal.F (SWAN in FVCOM)
!               Repeat :
!                    Call SNEXTI to update boundary conditions and input fields
!                    If IT>0                                              
!                      Call SWANCOMPUNSTRUC to calculate the wave field 
!                      Call SWANOUT to write results and create output
!     ----------------------------------------------------------------

!
!***********************************************************************
!
!   USE ALL_VARS
   USE schism_glbl, ONLY : skind,rkind,dkind,nea,npa,errmsg,iplg,it_main
   USE schism_glbl, only : start_year,start_month,start_day,start_hour,utc_start
   USE schism_glbl, only : RNDAY_SCHISM => rnday
   USE schism_glbl, only : ihot
   USE schism_glbl, only : nws,icou_elfe_wwm

   USE schism_msgp, ONLY : myrank,parallel_abort,nproc,parallel_barrier
   USE schism_msgp, ONLY : itype,comm

!   USE schism_glbl, ONLY :, uu2, vv2, tr_nd, & !tnd, snd, &
!   USE schism_glbl, ONLY : kfp, idry, nvrt, ivcor,ipgl,fdb,lfdb,errmsg,     &
!   USE schism_glbl, ONLY : np_global,ne_global,                             &
!   USE schism_glbl, ONLY : np,ne,nvrt,i34,isbe,isbs,elside,elnode,ilnd,nlnd,nland,isbnd

!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!#  endif      

   USE TIMECOMM                                                        
   USE OCPCOMM1, ONLY : REFDAY
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OUTP_DATA                                                       
   USE M_GENARR                                                        

   USE SWANGRIDDATA

!#  if defined (EXPLICIT)
!   USE MOD_ACTION_EX
!#  else
!   USE MOD_PETSC
!  USE MOD_ACTION_IM
!#  endif   
!   USE MOD_USGRID 
   USE VARS_WAVE
!#  if defined (SPHERICAL)   
!   USE MOD_SPHERICAL
!#  endif   
!#  if defined (NETCDF_IO)
!   USE MOD_NCDIO
!#  endif
!#  if defined (WAVE_SETUP)
!   USE MOD_WAVESETUP                                                     
!#  endif   
!
   IMPLICIT NONE
   include 'mpif.h'

!   INTEGER :: IUNIT, IOSTAT, IT0, ITW, SAVITE, ILEN, INERR, IERR

   INTEGER, INTENT(IN)     :: NSTEP_SWAN0
   REAL(rkind), INTENT(IN) :: DT_SCHISM0

   INTEGER :: IUNIT, IOSTAT,  SAVITE, ILEN, INERR, IERR
   INTEGER :: ISTAT, IF1, IL1, IDC,ISC,IP
!   CHARACTER :: PTYPE, PNAME *8, COMPUT *4, DTTIWR*18                     
   CHARACTER :: PTYPE, PNAME *8, DTTIWR*18                     
   CHARACTER*20 :: NUMSTR, CHARS(1)                                       
   CHARACTER*80 :: MSGSTR                                                 
   LOGICAL :: LOPEN     

! Now declared in module : mod_main_wave
!   INTEGER, ALLOCATABLE :: CROSS(:)                                  
!   INTEGER, ALLOCATABLE :: BGRIDP(:)                                   
!   REAL   , ALLOCATABLE :: BSPECS(:,:,:,:)                             
!   REAL   , ALLOCATABLE :: AC1(:,:,:), AC2_TMP(:,:,:),AC2_TMP1(:)                     
   REAL(rkind)   , ALLOCATABLE :: AC2_TMP(:,:,:)   !,AC2_TMP1(:)                     
!
!   REAL   , ALLOCATABLE :: BLKND(:), BLKNDC(:), OURQT(:)               

   REAL(rkind), ALLOCATABLE :: FTEMP(:)
   INTEGER :: I
   CHARACTER(LEN=100) :: NCFILE

   integer             :: start_jdate
   real(rkind) :: start_frac_jdate = -9999.0
   REAL(rkind), parameter :: ZERO = 0._rkind

!
!     --- initialize various data
   LEVERR=0                                                            
   MAXERR=1                                                            
   ITRACE=0                                                            
   INERR =0                                                            

   SERIAL = .TRUE.
   if(nproc.GT.1) SERIAL = .FALSE.

   CALL SWINIT (INERR)                                                 
   IF(INERR > 0) RETURN                                              
   
   COMPUT = '    '

!  --- read and process user commands
   IF(myrank==0) WRITE(PRINTF,*)'read and process user commands...'

! JL (SWREAD upgraded)
   CALL SWREAD !(COMPUT)   

! INPUTF unit : There is a conflict with open statement in function SCAN_FILE2,
! (SCAN_FILE2 use the same File Unit NuMBER)
! so close Do not forget to close it !
   CLOSE (UNIT=INPUTF)

!#  if !defined (EXPLICIT)
!   CALL PETSc_SET   
!#  endif     

! -----------------------------
! SCHISM < = > SWAN Interaction
! -----------------------------
! Note : some choice may replace those read befor with SWREAD
! from the INPUT File. So adapt for custom application ...

       IFLDYN(4) = 0 ! field 4: Friction, update ? .FALSE.
       IFLDYN(5) = 1 ! field 5 and 6: Wind update from schism
       IFLDYN(6) = 1
       VARWI     = 1

       IF (nws==0) THEN ! no atmos. forcing is applied
        IFLDYN(5) = 0
        IFLDYN(6) = 0
        VARWI     = 0
       ENDIF

       IF (icou_elfe_wwm == 1) THEN ! Full coupling
           !WLDEP       = DEP8
           !WATLEV      = ETA2
           !WATLEVOLD   = ETA1
           !DEP         = MAX(ZERO,WLDEP + WATLEV)
           !CURTXY(:,1) = UU2(NVRT,:)
           !CURTXY(:,2) = VV2(NVRT,:)
           !LSECU       = .TRUE.
           !LSEWL       = .TRUE.
           IFLDYN(2) = 1   ! field 2 and 3: current velocity, update  => LSECU = .TRUE.
           IFLDYN(3) = 1
           IFLDYN(4) = 1   ! field 4: Friction, update 
           IFLDYN(7) = 1   ! field 7: water level, update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
       ELSE IF (icou_elfe_wwm == 0) THEN ! No interaction at all, only wind from schism
           !WLDEP       = DEP8
           !WATLEV      = ZERO
           !WATLEVOLD   = ZERO
           !DEP         = WLDEP
           !CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           !CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           !LSECU       = .FALSE.
           !LSEWL       = .FALSE.
           COMPDA(:,JWLV2) = ZERO
           COMPDA(:,JWLV3) = ZERO
           COMPDA(:,JVX2)  = ZERO
           COMPDA(:,JVY2)  = ZERO
           COMPDA(:,JFRC2) = ZERO
           COMPDA(:,JFRC3) = ZERO
           IFLDYN(2) = 0   ! field 2 and 3: current velocity, no update  => LSECU = .FALSE.
           IFLDYN(3) = 0
           IFLDYN(4) = 0   ! field 4: Friction, no update
           IFLDYN(7) = 0   ! field 7: water level, no update => LSEWL = .FALSE.
           VARWLV = .FALSE.
           DYNDEP = .FALSE.
           VARFR  = .FALSE.
       ELSE IF (icou_elfe_wwm == 2) THEN ! Currents,friction and water levels in wwm but no radiation stress in SCHISM
           !WLDEP       = DEP8
           !WATLEV      = ETA2
           !WATLEVOLD   = ETA1
           !DEP         = MAX(ZERO, WLDEP + WATLEV)
           !CURTXY(:,1) = UU2(NVRT,:)
           !CURTXY(:,2) = VV2(NVRT,:)
           !LSECU       = .TRUE.
           !LSEWL       = .TRUE.
           IFLDYN(2) = 1   ! field 2 and 3: current velocity, update  => LSECU = .TRUE.
           IFLDYN(3) = 1
           IFLDYN(4) = 1   ! field 4: Friction, update
           IFLDYN(7) = 1   ! field 7: water level, update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
       ELSE IF (icou_elfe_wwm == 3) THEN ! No current, No Friction and no water levels in wwm but radiation stress in SCHISM
           !WLDEP       = DEP8
           !WATLEV      = ZERO
           !WATLEVOLD   = ZERO
           !DEP         = WLDEP
           !CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           !CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           !LSECU       = .FALSE.
           !LSEWL       = .FALSE.
           COMPDA(:,JWLV2) = ZERO
           COMPDA(:,JWLV3) = ZERO
           COMPDA(:,JVX2)  = ZERO
           COMPDA(:,JVY2)  = ZERO
           COMPDA(:,JFRC2) = ZERO
           COMPDA(:,JFRC3) = ZERO
           IFLDYN(2) = 0   ! field 2 and 3: current velocity, no update  => LSECU = .FALSE.
           IFLDYN(3) = 0   ! 
           IFLDYN(4) = 0   ! field 4: Friction, no update
           IFLDYN(7) = 0   ! field 7: water level, no update => LSEWL = .FALSE.
           VARWLV = .FALSE.
           DYNDEP = .FALSE.
           VARFR  = .FALSE. ! Variable Friction
       ELSE IF (icou_elfe_wwm == 4) THEN ! No current but water levels in wwm and radiation stresss in SCHISM
           !WLDEP       = DEP8
           !WATLEV      = ETA2
           !WATLEVOLD   = ETA1
           !DEP         = WLDEP
           !CURTXY(:,1) = 0.!UU2(NVRT,:)
           !CURTXY(:,2) = 0.!UU2(NVRT,:)
           !LSECU       = .FALSE.
           !LSEWL       = .TRUE.
           COMPDA(:,JVX2) = ZERO
           COMPDA(:,JVY2) = ZERO
           IFLDYN(2) = 0   ! field 2 and 3: current velocity, no update  => LSECU = .FALSE.
           IFLDYN(3) = 0
           IFLDYN(4) = 0   ! field 4: Friction, no update
           IFLDYN(7) = 1   ! field 7: water level, update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
           ICUR   = 0
       ELSE IF (icou_elfe_wwm == 5) THEN ! No current but water levels in SWAN and no radiation stress in SCHISM
           !COMPDA(:,JWLV2) = ZERO
           !COMPDA(:,JWLV3) = ZERO
           COMPDA(:,JVX2) = ZERO
           COMPDA(:,JVY2) = ZERO
           IFLDYN(2) = 0   ! field 2 and 3: current velocity, no update  => LSECU = .FALSE.
           IFLDYN(3) = 0
           IFLDYN(4) = 0   ! field 4: Friction, no update
           IFLDYN(7) = 1   ! field 7: water level, no update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
       ELSE IF (icou_elfe_wwm == 6) THEN ! Currents, friction but no water levels in wwm and radiation stress in SCHISM
           COMPDA(:,JWLV2) = ZERO
           COMPDA(:,JWLV3) = ZERO
           IFLDYN(2) = 1   ! field 2 and 3: current velocity, update  => LSECU = .TRUE.
           IFLDYN(3) = 1
           IFLDYN(4) = 1   ! field 4: Friction, update
           IFLDYN(7) = 0   ! field 7: water level, update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
       ELSE IF (icou_elfe_wwm == 7) THEN ! Currents, friction but no water levels in wwm and no radiation stress in SCHISM
           COMPDA(:,JWLV2) = ZERO
           COMPDA(:,JWLV3) = ZERO
           IFLDYN(2) = 1   ! field 2 and 3: current velocity, update  => LSECU = .TRUE.
           IFLDYN(3) = 1
           IFLDYN(4) = 1   ! field 4: Friction, update
           IFLDYN(7) = 0   ! field 7: water level, update => LSEWL = .TRUE.
           VARWLV = .TRUE.
           DYNDEP = .TRUE.
           VARFR  = .TRUE. ! Variable Friction
       ELSE
           call parallel_abort('see icou_elfe_wwm value')
       END IF

!
! --- Adding new pointers in array COMPDA
!
!#ifdef schism ! Define Mandatory wave parameters in SCHISM-SWAN

   ! IVTYPE = 54 'DISB' Sfric' 'Bottom friction dissipation' 'm2/s'
   MCMVAR = MCMVAR+1
   JDSXB  = MCMVAR
   ! IVTYPE = 55 'DISSU' 'Ssurf' 'Surf breaking dissipation' 'm2/s'
   MCMVAR = MCMVAR+1
   JDSXS  = MCMVAR
   ! IVTYPE = 56 'DISW' 'Swcap' 'Whitecapping dissipation' 'm2/s'
   MCMVAR = MCMVAR+1
   JDSXW  = MCMVAR
   ! IVTYPE = 36 'ZLEN' 'Zlen' 'Zero velocity thickness of boundary layer'
   MCMVAR = MCMVAR+1
   JZEL   = MCMVAR
   ! IVTYPE = xxx 'xxx'
   MCMVAR = MCMVAR+1
   JUSTAR = MCMVAR
   ! IVTYPE = xxx 'xxx'
   MCMVAR = MCMVAR+1
   JCHARN = MCMVAR
   ! IVTYPE = 38 'CDRAG' 'Cdrag' 'Drag coefficient'
   MCMVAR = MCMVAR+1
   JCDRAG = MCMVAR

!#endif

!
!  --- Replace previous definitions (DTW, RDTIM etc from SWREAD()
!
   DTW   = NSTEP_SWAN0*DT_SCHISM0  ! = SWAN DT, sec
   RDTIM = 1._rkind/DTW            ! RDTIM used in transport schemes (see swantransxxxx)
!
!  --- Set the SWAN Clock
!
   REFDAY = 0  ! Reference Day, set to 0
   !DTW = NSTEP_SWAN0*DT_SCHISM0     ! = SWAN DT, sec
   !ITW = INT(it_main/NSTEP_SWAN0) ! = SWAN IT,
   ITW = 0

! calculate the starting Julian date and starting fractional Julian date
   start_jdate = jd(start_year,start_month,start_day) ! days

   start_frac_jdate = real(start_jdate,rkind) &
     &             + (start_hour + utc_start) / 24._rkind
   IF(myrank==0) PRINT*,'start_frac_jdate = ', start_frac_jdate

   ! TIMCO     [    ] Time and date of the computation during the simulation (in
   !                  seconds since the reference day (REFDAY))
   !TIMCO =  Float(REFDAY) + (ITW*DTW)  ! sec
!   TIMCO = start_frac_jdate*DAY2SEC + (ITW*DTW)  ! sec

   ! TFINC     [    ] End time and date of the computation (in seconds since
   !                  the reference day (REFDAY))
   TFINC = start_frac_jdate*DAY2SEC + RNDAY_SCHISM*DAY2SEC  ! sec

   ! TINIC     [    ] Start time and date of the computation (in seconds since
   !                  the reference day (REFDAY))
   TINIC = start_frac_jdate*DAY2SEC  ! sec

   ! Hard coding for Non-Stationary case
   NSTATC = 1                                                    ! 40.00
!
!  --- SWAN cycles
!
   MTC = NINT ((TFINC - TINIC)/DTW)
   IF(MOD(TINIC-TFINC,DTW) > 0.01_rkind*DTW .AND.                     &
      MOD(TINIC-TFINC,DTW) < 0.99_rkind*DTW)                          &
      CALL parallel_abort &
      ('DTW is not a fraction of the computational period')

   TIMCO = TINIC   

   IF(myrank==0) PRINT*,'SWAN DTW (sec)   = ',DTW
   IF(myrank==0) PRINT*,'SWAN TIMCO (sec) = ',TIMCO
   IF(myrank==0) PRINT*,'SWAN Total Loops = ',MTC

!  --- allocate some arrays meant for computation                    

   IF(NUMOBS .GT. 0) THEN
!     IF(.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(2,M))            

      IF (OPTG.NE.5) THEN                                           ! 40.80
!        structured grid                                            ! 40.80
         MSGSTR = 'wrong Grid definition, OPTG should be 5'
         WRITE(errmsg,*) MSGSTR
         CALL parallel_abort(errmsg)

         ILEN = 2*MCGRD                                             ! 40.80
      ELSE                                                          ! 40.80
!        unstructured grid                                          ! 40.80
         ILEN = nfaces                                              ! 40.80
      ENDIF                                                         ! 40.80
      IF (.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(ILEN))

   ELSE
     IF(.NOT.ALLOCATED(CROSS)) ALLOCATE(CROSS(0))      ! since 40.80
   ENDIF      

! JL source of trouble
!   IF (NPTST.NE.NPTSTA) THEN                                        
!      MSGSTR = 'Something wrong NPTST should be = NPTSTA'
!      WRITE(*,*) MSGSTR
!      CALL PSTOP
!   ENDIF
                                                
   IF(.NOT.ALLOCATED(BSPECS)) ALLOCATE(BSPECS(MDC,MSC,NBSPEC,2))    
   IF(.NOT.ALLOCATED(BGRIDP)) ALLOCATE(BGRIDP(6*NBGRPT))            

!  --- do some preparations before computation                       

!   CALL SWPREP ( BSPECS, BGRIDP, CROSS , SPCDIR, SPCSIG )                            
!JL (SWPREP Upgraded) 
!   CALL SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID, YCGRID, KGRPNT,&    !40.31
!     &           KGRBND, SPCDIR, SPCSIG )                            !40.31
   CALL SWPREP (                        XCGRID, YCGRID, KGRPNT,&    !40.31
     &           KGRBND, SPCDIR, SPCSIG )                            !40.31

!   CALL SwanPrepComp ( CROSS )                                       !40.80
   CALL SWANPREPCOMP 

!  --- check all possible flags and if necessary change
!      if option is not correct

   CALL ERRCHK                                                       

!  --- initialisation of necessary grids for depth,
!      current, wind and friction

   IF (ALOCMP.AND.ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)              ! 40.97
   IF (.NOT.ALLOCATED(COMPDA)) THEN 
      ALLOCATE(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)                         ! 40.97
      ALOCMP = .FALSE.
   END IF

   IF(ISTAT .NE. 0)THEN             
     MSGSTR = 'Allocation problem: array COMPDA'                               
     WRITE(errmsg,*) MSGSTR
     CALL parallel_abort(errmsg)
   END IF        

!
!  --- Populate COMPDA with depth, current, wind, friction, mud,.. fields
!
!  JL: SWAN in SCHISM:  SWRBC update done
   CALL SWRBC 

!  SWAN in SCHISM : ALLOCATE ACTION DENSITY Now
   IF(.NOT.ALLOCATED(AC2)) ALLOCATE(AC2(MDC,MSC,MCGRD),STAT=ISTAT)
   IF(ISTAT.NE.0)THEN
     MSGSTR = 'Allocation problem: array AC2'
     WRITE(errmsg,*) MSGSTR
     CALL parallel_abort(errmsg)
   END IF
   AC2 = 0._rkind

   IF(ITEST.GE.2) THEN
    IF(.NOT.ALLOCATED(AC2_SUM)) ALLOCATE(AC2_SUM(MCGRD),STAT=ISTAT)
    AC2_SUM = 0._rkind
    SUM_AC2 = 0._rkind
   END IF

!  --- allocate AC1 in case of non-stationary situation or in case   
!      of using the S&L scheme                                       
! JL TEST sur AC1
   IF(NSTATM.EQ.1 .AND. MXITNS.GT.1 .OR. PROPSC.EQ.3)THEN        !40.31
     IF(.NOT.ALLOCATED(AC1))THEN                           
       ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)                     
     ELSE IF(SIZE(AC1).EQ.0)THEN                                  
       DEALLOCATE(AC1)                                             
       ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)                     
     END IF                                                         
     IF(ISTAT.NE.0)THEN                                         
       MSGSTR = 'Allocation problem: array AC1'                              
       WRITE(errmsg,*) MSGSTR
       CALL parallel_abort(errmsg)
     END IF                                                         
     AC1 = 0._rkind                                                   
   ELSE                                                              
     IF(.NOT.ALLOCATED(AC1)) ALLOCATE(AC1(0,0,0))                   
   ENDIF

!   IF(.NOT.ALLOCATED(AC1))THEN
!       ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
!   ELSE
!       DEALLOCATE(AC1)
!       ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
!   END IF
!   AC1 = AC2

   IF(LEVERR > MAXERR)THEN      
                                  
     WRITE (PRINTF, 6010) LEVERR
     IF(LEVERR < 4) WRITE (PRINTF, 6011)                           
6010 FORMAT(' ** No start of computation because of error level:',I3)
6011 FORMAT(' ** To ignore this error, change [maxerr] with the',     &
            ' SET command')                                          
   ELSE
!
     IF(ITEST >= 40)THEN                                           
       IF(NSTATC == 1)THEN                                         
         WRITE (PRINTF, '(" Type of computation: dynamic")')         
       ELSE                                                          
         IF(ONED)THEN                                              
           WRITE (PRINTF, '(" Type of computation: static 1-D")')    
         ELSE                                                        
           WRITE (PRINTF, '(" Type of computation: static 2-D")')    
         ENDIF                                                       
       ENDIF                                                         
     ENDIF
! 
!#ifdef schism 
! starting from a Hotstart or coldstart ?
     IF(ihot/=0) ICOND = 4
!#endif

     IF(NSTATC.EQ.1)THEN                                           
       IT0 = 0                                                       
       IF(ICOND.EQ.1)THEN  ! COLDSTART                             
!
!        --- compute default initial conditions (coldstart only)
!
         CALL SWINCO (                 XCGRID, YCGRID,      &           !40.31
     &                      KGRPNT, SPCDIR, SPCSIG, XYTST )             !40.31

!
!        --- reset ICOND to prevent second computation of initial condition
         ICOND = 0                                                   

       ELSE                !  ----  SWAN HOT-STARTING ----------

!#ifdef schism
         ! See Hydro/schism_init, last section in the tail:
         !  line 6677   USE_SWAN and ihot\=0 ...
         !  - last AC2 will be read from hotstart.nc and distribute to all ranks
!#endif
         
       ENDIF

!        --- check out for NaN
       IF(ITEST.GT.1) THEN
          DO I=1,MCGRD
           AC2_SUM(I) = sum(AC2(:,:,I))
          ENDDO
          SUM_AC2 = sum(AC2_SUM)

          IF(SUM_AC2/=SUM_AC2) THEN
           WRITE (PRINTF, '(" NaN in AC2 after SWINCO")')
           DO I=1,MCGRD
            WRITE (PRINTF,*) '"INDEX INDEX_GL AC2_SUM")',I,iplg(I),AC2_SUM(I)
           ENDDO
          write(errmsg,*)'NaN from SWAN after SWINCO'
          call parallel_abort(errmsg)
          ENDIF
       ENDIF

     ELSE
       IT0 = 1
     ENDIF
!
!    --- NESTING : Populate arrays with boundary forcing (WBAC** stuff)
!
     if(myrank.eq.0) print*,'SWMAIN_SETUP NESTING=',NESTING

     IF(NESTING) THEN

           IF(myrank==0) WRITE(PRINTF,*)&
          'SET THE INITIAL WAVE BOUNDARY CONDITION'

           !WBBMJD = TIMCO*SEC2DAY
           !WBTMJD = WBBMJD

           CALL INIT_WAVE_BOUNDARY_CONDITION
           ! reset WBBMJD
           WBBMJD = TIMCO*SEC2DAY
           WBTMJD = WBBMJD

     ENDIF

!         --- synchronize nodes                                           40.30
     IF(.NOT.SERIAL) call parallel_barrier

     IF(myrank==0) WRITE(PRINTF,*)'SET THE INITIAL WAVE BOUNDARY CONDITION END'

!     if(myrank==0) print*,"IT0  =",IT0
!     if(myrank==0) print*,"TIMCO=",TIMCO

!JL see Later, but may not work with unswan ...
#if 0
#    if defined(WAVE_SETUP)   
     IF(LSETUP > 0) CALL ALLOC_VARS_WSU   
#    endif
#endif
!    --- Ready to loop over time steps                                        
     
     ITW = IT0

   END IF ! NSTATC.EQ.1


   RETURN                                                              
   END SUBROUTINE SWMAIN_SETUP
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_LOOP(IT_SCHISM,DT_SCHISM0,NSTEP_SWAN0)                                                   
!                                                                      *
!***********************************************************************
!
!     SWMAIN subroutine, calling SWINIT, SWREAD, SWCOMP 
!
!***********************************************************************
!
!   USE ALL_VARS
!#  if defined (MULTIPROCESSOR)
!   USE MOD_PAR
!#  endif      
   USE schism_msgp, ONLY : myrank,parallel_abort,parallel_barrier
   USE schism_glbl, ONLY : errmsg,rkind,icou_elfe_wwm,dp,npa,nspool,nws
   USE schism_glbl, ONLY : iplg,RNDAY_SCHISM => rnday

   USE TIMECOMM                                                        
   USE OCPCOMM1, ONLY : REFDAY
   USE OCPCOMM2                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OUTP_DATA                                                       
   USE M_GENARR                                                        

!  USE MOD_INPUT, ONLY: NC_DAT

   USE VARS_WAVE

!#  if defined (SPHERICAL)   
!   USE MOD_SPHERICAL
!#  endif   
!#  if defined (NETCDF_IO)
!   USE MOD_NCDIO
!#  endif
!#  if defined (WAVE_SETUP)
!   USE MOD_WAVESETUP                                                     
!#  endif   
!
   IMPLICIT NONE

   INTEGER, INTENT(IN)     :: IT_SCHISM
   INTEGER, INTENT(IN)     :: NSTEP_SWAN0 !, icou_elfe_wwm
   REAL(rkind), INTENT(IN) :: DT_SCHISM0


!   INTEGER :: IUNIT, IOSTAT, IT0, ITW, SAVITE, ILEN, INERR, IERR
   INTEGER :: IUNIT, IOSTAT , SAVITE, ILEN, INERR, IERR
   INTEGER :: ISTAT, IF1, IL1, IDC,ISC,IP   
!   CHARACTER :: PTYPE, PNAME *8, COMPUT *4, DTTIWR*18                     
   CHARACTER :: PTYPE, PNAME *8, DTTIWR*18                     
   CHARACTER*20 :: NUMSTR, CHARS(1)                                       
   CHARACTER*80 :: MSGSTR                                                 
   LOGICAL :: LOPEN     

   REAL(rkind), ALLOCATABLE :: FTEMP(:)
   INTEGER :: I
   CHARACTER(LEN=100) :: NCFILE
   REAL(rkind), parameter :: ZERO = 0._rkind
#ifdef INCLUDE_TIMING
   REAL(rkind)        :: TIME1, TIME2, TIME3, TIME4
#endif


#ifdef INCLUDE_TIMING
   include 'mpif.h'
   TIME1 = mpi_wtime()
#endif

!   DTW = NSTEP_SWAN0*DT_SCHISM0     ! = SWAN DT, sec
   ITW = MAX(1,INT(IT_SCHISM/NSTEP_SWAN0)) ! = SWAN IT,

   ! TIMCO     [    ] Time and date of the computation during the simulation (in
   !                  seconds since the reference day (REFDAY))   
!   TIMCO =  TINIC + MyREAL(ITW)*DTW   ! sec

   ! TFINC     [    ] End time and date of the computation (in seconds since
   !                  the reference day (REFDAY))
!   TFINC =  MyREAL(REFDAY) + RNDAY_SCHISM*DAY2SEC  ! sec
   ! TINIC     [    ] Start time and date of the computation (in seconds since
   !                  the reference day (REFDAY))
!   TINIC = MyREAL(REFDAY)  ! sec

!   IF(myrank==0) PRINT*,'SWAN ITW = ',ITW                     
!   IF(myrank==0) PRINT*,'SWAN TIMCO (sec) = ',TIMCO
   ! Jerome TIMCO computation                      
   ! TIMCO     [    ] Time and date of the computation during the simulation (in
   !                  seconds since the reference day (REFDAY))
   !    TIMCO = DTW*IT0
!        TIMCO = DTW*ITW
   !    IT0 = ITW


   IF(myrank==0) PRINT*,'SWAN LOOP START, TIMCO (sec) = ',TIMCO
   IF(myrank==0) WRITE(PRINTF,*)'SWAN LOOP START, TIMCO (sec) = ',TIMCO

!
!print*,LEVERR,MAXERR
       IF(LEVERR > MAXERR)THEN                                    
         WRITE (PRINTF, 6030) LEVERR                                 
         IF(LEVERR < 4) WRITE (PRINTF, 6011)                       
6030     FORMAT(' ** No continuation of computation because ',       &
                  'of error level:',I3)                              
6011     FORMAT(' ** To ignore this error, change [maxerr] with the',&
            ' SET command') 
         MSGSTR = 'error level, see PRINT file'
         WRITE(errmsg,*) MSGSTR
         CALL parallel_abort(errmsg)                                         
       ENDIF    

!      --- synchronize nodes                                           40.30
       IF(.NOT.SERIAL) call parallel_barrier
!
!      --- update boundary conditions and input fields

!       IF(myrank==0)PRINT*,'Before SNEXTI'

       CALL SNEXTI (                                             &   ! 40.31
     &                    SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT,&   ! 40.31
     &                    XYTST , DEPTH , WLEVL , FRIC  , UXB   ,&   ! 40.31
     &                    UYB   , NPLAF , TURBF , MUDLF , WXI   ,&   ! 40.59 40.35 40.55 40.31
     &                    WYI   )
!      IF (STPNOW()) RETURN                                               !34.01

!
!         --- synchronize nodes                                           40.30
       IF(.NOT.SERIAL) call parallel_barrier

!       IF(myrank==0)PRINT*,'After SNEXTI'

!      --- Offline Mode, depth fixed 
!
       IF (icou_elfe_wwm == 0 .or. icou_elfe_wwm == 3) THEN
          COMPDA(1:npa,JDP1) = dp(1:npa)
          COMPDA(:,JDP2) = COMPDA(:,JDP1)
          IF (LSETUP.GT.0) THEN ! usefull for SETUP calc in SWAN
            COMPDA(:,JDPSAV) = COMPDA(:,JDP2)
           ENDIF
       END IF
!
!       IF(myrank==0) print*,'COMPUT=',COMPUT,ITW

       IF (COMPUT /= 'NOCO' .AND. ITW > 0) THEN                      

         SAVITE = ITEST                                              
         IF (ICOTES > ITEST) ITEST = ICOTES
!
!        --- compute action density for current time step
!
!         CALL SWCOMP( AC1, ITW, CROSS )                                      
!         CALL SWCOMP(ITW)                                      
         !PRINT*,'AFTER SWCOMP'
!         WRITE(IPT,*) 'CPU:',MYID,'AFTER SWCOMP'!

#ifdef INCLUDE_TIMING
         TIME2 = mpi_wtime()
#endif

! 03/2017 JLEFEVRE Update 
         IF(myrank==0) write(PRINTF,*)'! CALLING SWANCOMPUNSTRUC...',ITW

         CALL SWANCOMPUNSTRUC(spcsig, spcdir, xytst, ITW)

!        --- checkout for NAN
         IF(ITEST.GT.2) THEN
          DO I=1,MCGRD
           AC2_SUM(I) = sum(AC2(:,:,I))
          ENDDO
          SUM_AC2 = sum(AC2_SUM)

          IF(SUM_AC2/=SUM_AC2) THEN
           WRITE (PRINTF, '(" NaN in AC2 after SWANCOMPUNSTRUC")')
           DO I=1,MCGRD
            WRITE (PRINTF,*) '"INDEX INDEX_GL AC2_SUM")',I,iplg(I),AC2_SUM(I)
           ENDDO
          write(errmsg,*)'NaN from SWAN after SWANCOMPUNSTRUC'
          call parallel_abort(errmsg)
          ENDIF
         ENDIF


#ifdef INCLUDE_TIMING
         TIME3 = mpi_wtime()
#endif

!
!        --- set ICOND=4 for stationary computation, for next
!            (stationary) COMPUTE command                            
         ICOND = 4                                                   
!
!        --- check whether computed significant wave height at       
!            boundary differs from prescribed value given in         
!            boundary command values of incident Hs                  
!
!         if(mod(itw,CDF_INT) == 0)then
!         if(mod(itw,10) == 0)then
!         if(mod(itw,30) == 0)then
!           IF ( BNDCHK ) THEN                                          
!             CALL HSOBND (COMPDA(1,JHSIBC),DEPTH)         
!           ENDIF              
!	 end if  
         !PRINT*,'AFTER HSOBND'
!
!--NETCDF OUTPUT---------------------------------------------------------------!
!
! print*,'CDF_INT,ITW=',ITW
!#        if defined (NETCDF_IO)
         !IF(CDF_INT /= 0 .AND. CDF_OUT)THEN

!#       if defined(WAVE_ONLY)
!        IF (icou_elfe_wwm == 0 .OR. icou_elfe_wwm == 2 .OR. icou_elfe_wwm == 5 .OR. icou_elfe_wwm == 7) THEN
        ! No coupling with schism, call swanout for archive only
!          IF(myrank==0) print*,'IT_SCHISM,nspoo',IT_SCHISM,nspool,mod(IT_SCHISM,nspool)
!          IF(mod(IT_SCHISM,nspool)==0) then  ! save some cycle ...
!           IF(myrank==0) write(16,*)'! CALLING SWANOUT...'
!           CALL SWANOUT
!           IF(myrank==0) print*, 'AFTER SWANOUT'
!          ENDIF
!        ELSE
          IF(myrank==0) write(PRINTF,*)'! CALLING SWANOUT...',ITW
          CALL SWANOUT
          IF(myrank==0) print*, 'AFTER SWANOUT'
!        ENDIF

!--RESTART OUTPUT--------------------------------------------------------------!
!
#        if defined (NETCDF_IO)
!JQI         IF(CDF_RST /= 0)THEN
!JQI           IF(MOD(ITW,CDF_RST) == 0)THEN
!JQI             IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
!JQI             CALL OUT_NETCDF_RST(SPCDIR,SPCSIG,ITW)
!JQI           END IF
!JQI         END IF
#        endif

!!           --- write the SWAN hot-start file, if necessary               41.20
!            IF ( WriteSwanHotStart ) THEN
!               CALL BACKUP ( AC2,SPCSIG,SPCDIR,KGRPNT,XCGRID,YCGRID )
!               WriteSwanHotStart = .FALSE.
!            ENDIF

!         --- synchronize nodes                                           40.30
         call parallel_barrier
!
         ITEST = SAVITE                                              

#ifdef INCLUDE_TIMING
         TIME4 = mpi_wtime()
#endif

       ENDIF
!
       IF(ITW == IT0 .AND. .NOT.ALLOCATED(OURQT))THEN             
         ALLOCATE (OURQT(MAX_OUTP_REQ))                             
         OURQT = -9999.                                             
       ENDIF                                                         
!
       SAVITE = ITEST                                                
       IF(IOUTES > ITEST) ITEST = IOUTES

!
       IF(ERRPTS > 0) REWIND(ERRPTS)                               
       ITEST = SAVITE                                                

!      --- update time
       TIMCO = TIMCO + DTW
!       IF(myrank==0) PRINT*,'SWAN LOOP END, TIMCO (sec) = ',TIMCO

       IF(NSTATM == 1)THEN                                         
!         IF(NSTATC == 1 .AND. ITW < MTC) TIMCO = TIMCO + DTW         
!JQI           CHTIME = DTTIWR(ITMOPT, TIMCO)                              

!JL MUTE
!         IF(myrank==0 .AND. NSTATC == 1) WRITE (PRINTF, 222) CHTIME, TIMCO          
!222      FORMAT(' Time of computation ->  ',A,' in sec:', F16.0)      

         
       ENDIF                                                         

       ! WRITE(IPT,*) 'CPU:',MYID,'FINISHED SWANMAINLOOP'!
       !PRINT*,'FINISHED SWANMAINLOOP'


#ifdef INCLUDE_TIMING
       IF(myrank==0) THEN 
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') '-----TOTAL TIMINGS-----'
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') 'INTEGRATION        ', TIME3-TIME2
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') 'OUTPUT TO SCHISM   ', TIME4-TIME3
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') 'TOTAL TIME         ', TIME4-TIME1
           WRITE(PRINTF,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
!           CALL FLUSH(PRINTF)
       ENDIF
#endif



   RETURN                                                              
   END SUBROUTINE SWMAIN_LOOP
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWMAIN_CLOSE                                                   
!                                                                      *
   USE OCPCOMM4                                                        
   USE VARS_WAVE
   !USE MOD_PETSC

   IMPLICIT NONE

   INTEGER :: IUNIT
   LOGICAL :: LOPEN     

   DO IUNIT=1,HIOPEN                                                   
     INQUIRE(UNIT=IUNIT,OPENED=LOPEN)                                  
     IF(LOPEN .AND. IUNIT /= PRINTF) CLOSE(IUNIT)                       
   END DO                                                              
!
   INQUIRE(UNIT=PRINTF,OPENED=LOPEN)                                   
   IF(LOPEN) CLOSE(PRINTF)                                            
!
!  --- deallocate all allocated arrays                                 

   IF(ALLOCATED(AC1   )) DEALLOCATE(AC1   )                           
   IF(ALLOCATED(BGRIDP)) DEALLOCATE(BGRIDP)                           
   IF(ALLOCATED(BSPECS)) DEALLOCATE(BSPECS)                           
   IF(ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)                           
   IF(ALLOCATED(CROSS )) DEALLOCATE(CROSS )                           
   IF(ALLOCATED(OURQT )) DEALLOCATE(OURQT )                           
   IF(ALLOCATED(BLKND )) DEALLOCATE(BLKND )                           

   IF(ALLOCATED(AC2_SUM)) DEALLOCATE(AC2_SUM)

#  if !defined (EXPLICIT)
!   CALL PETSc_CLEANUP
   print*,'PETSc_CLEANUP ToDo'
#  endif
   
   RETURN                                                              
   END SUBROUTINE SWMAIN_CLOSE
 
!***********************************************************************
!                                                                      *
   SUBROUTINE SWINIT (INERR)                                           
!                                                                      *
!***********************************************************************
!
!     Initialize several variables and arrays
!
!***********************************************************************

   USE OCPCOMM1                                                        
   USE OCPCOMM2                                                        
   USE OCPCOMM3                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM1                                                         
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE TIMECOMM                                                        
   USE M_SNL4                                                          
   USE M_BNDSPEC                                                       
   USE SwanGriddata

   !USE ALL_VARS, ONLY: MGL,MT 
   USE schism_glbl, ONLY : np_global,npa,rkind
   USE schism_msgp, ONLY : myrank,parallel_abort


   IMPLICIT NONE                                       

   INTEGER :: INERR
   LOGICAL :: STPNOW                                                      

   INTEGER :: IGRID,IVT,IVTYPE,MXOUTAR
   
   VERTXT = BLANK                                                      
   VERNUM = 41.10
   WRITE (VERTXT, '(F5.2)') VERNUM                                     
!
   CALL OCPINI ('swaninit', .TRUE.,INERR) 
   IF(INERR > 0) RETURN                                              
   IF(STPNOW()) RETURN                                                
!
   WRITE (PRINTF, 6010) VERTXT                                         
6010 FORMAT (/,20X,'---------------------------------------',          &
             /,20X,'                 SWAN',                            &
             /,20X,'SIMULATION OF WAVES IN NEAR SHORE AREAS',          &
             /,20X,'         VERSION NUMBER ', A,                      &
             /,20X,'---------------------------------------',/)
!
!  ***** initial values for common variables *****
!  ***** names *****
   PROJID = 'SWAN'
   PROJNR = BLANK
   PROJT1 = BLANK
   PROJT2 = BLANK
   PROJT3 = BLANK
   FNEST  = BLANK
   FBCR   = BLANK
   FBCL   = BLANK
   UH     = 'm'
   UV     = 'm/s'
   UT     = 'sec'
   UL     = 'm'
   UET    = 'm3/s'
   UDI    = 'degr'
   UST    = 'm2/s2'
   UF     = 'N/m2'
   UP     = 'W/m'
   UAP    = 'W/m2'
   UDL    = 'm2/s'
!  ***** physical parameters *****
   GRAV_W   = 9.81_rkind
   WLEV   = 0._rkind
   CASTD  = 0._rkind             ! const. air-sea temp diff                
   CDCAP  = 99999._rkind
   PI_W     = 4._rkind*ATAN(1.)                                                
   PI2_W    = 2._rkind*PI_W
   UNDFLW = 1.E-15
   DNORTH = 90._rkind                                                        
   DEGRAD = PI_W/180._rkind
   RHO_W    = 1025._rkind
!  power of tail in spectrum, 1: E with f, 2: E with k,
!                             3: A with f, 4: A with k
   PWTAIL(1) = 4._rkind
   PWTAIL(2) = 2.5_rkind
   PWTAIL(3) = PWTAIL(1)+1._rkind                                            
   PWTAIL(4) = 3._rkind
!  ***** number of computational grid points ****                      
   MCGRD   = npa
!   MCGRDGL = 1           !JL Declared in M_PARALL            ! 40.30
   NGRBND  = 0                                                ! 40.00
!   NGRBGL  = 0                                               ! 40.30
   nverts  = 0                                                ! 40.80
   nvertsg = 0                                                        
!  time of computation                                                 
   TIMCO = -1.E10                                                      
   CHTIME = '    '                                                     
!  boundary conditions                                                 
   NBFILS = 0                                                          
   NBSPEC = 0                                                          
   NBGRPT = 0                                                          
!   NBGGL  = 0                                                          
   FSHAPE = 2                                                          
   DSHAPE = 2                                                          
   PSHAPE(1) = 3.3_rkind                                                     
   PSHAPE(2) = 0.1_rkind                                                     
!  ***** input grids *****
   DO IGRID = 1, NUMGRD
     XPG(IGRID)    = 0._rkind
     YPG(IGRID)    = 0._rkind
     ALPG(IGRID)   = 0._rkind
     COSPG(IGRID)  = 1._rkind                                                
     SINPG(IGRID)  = 0._rkind
!     DXG(IGRID)    = 0.
!     DYG(IGRID)    = 0.
!     MXG(IGRID)    = 0
!     MYG(IGRID)    = 0
     LEDS(IGRID)   = 0
     STAGX(IGRID)  = 0._rkind                                                
     STAGY(IGRID)  = 0._rkind                                                
     EXCFLD(IGRID) = -1.E20                                            
     IFLDYN(IGRID) = 0                                                 
     IFLTIM(IGRID) = -1.E20                                            
   END DO
!
   OPTG   = 5
   MXC    = 0
   MYC    = 0
!   MXCGL  = 0         !JL Declared in M_PARALL    
!   MYCGL  = 0         !JL Declared in M_PARALL
!   MXF    = 1                                                          
!   MXL    = 0                                                          
!   MYF    = 1                                                          
!   MYL    = 0                                                          

   MSC    = 0
   MDC    = 0
   MTC    = 1
   ICOMP  = 1
   ALPC   = 0._rkind
   FULCIR = .TRUE.
   SPDIR1 = 0._rkind
   asort  = -999._rkind                                                     ! 41.48
   CCURV  = .FALSE.
!  number of points needed in computational stencil:
   ICMAX  = 3
!  ***** numerical scheme *****
   NCOR   = 1
   NSTATM = -1                                                         
   NSTATC = -1                                                         
   NCOMPT = 0                                                          
!
!  initialise number of iterations stationary and nonstationary        
   MXITST = 15                                                         
   MXITNS = 1                                                          
   ITERMX = MXITST                                                     
   ICUR   = 0
   IDIF   = 0
   IINC   = 0
!
!  --- meaning IREFR:                                                  
!      IREFR = -1: limiter on Ctheta activated                         
!      IREFR =  1: No limiter on Ctheta                                
!      IREFR =  0: No refraction                                       
!
   IREFR  = 1 !1                                                          
   ITFRE  = 1
   IWIND  = 0
   IDRAG  = 2      ! 2nd Order Polynomial Fit, Default since  41.49 41.33
!     when coupled with ADCIRC, the default is to use ADCIRC drag formulation
!ADC      IDRAG  = 0
   IGEN   = 3                                                          
   IQUAD  = 2                                                          
   IWCAP  = 1                                                          
   ISURF  = 1
   IBOT   = 0
   ITRIAD = 0
   IMUD   = 0    ! 40.59
   IVEG   = 0    ! 40.55
   ITURBV = 0    ! 40.35
   VARWI  = .FALSE.
   VARFR  = .FALSE.
   VARWLV = .FALSE.                                                    
   VARAST = .FALSE. ! True means spatially variable air-sea t.d. 
   VARMUD = .FALSE.                                                  !  40.59
   VARNPL = .FALSE.                                                  !  40.55
   VARTUR = .FALSE.                                                  !  40.35
   U10    = 0._rkind
   WDIP   = 0._rkind
   INRHOG = 0      ! choice for output based on "variance" or "true energy" 
                   ! ie. wave dissipation 0 means unit = m2/s, 1 means unit= W/m2 
                   ! (by applying the factor rho*g) . Here m2/s
   DEPMIN = 0.08_rkind    !0.05
   SY0    = 3.3_rkind
   SIGMAG = 0.1_rkind
   XOFFS  = 0._rkind
   YOFFS  = 0._rkind
   LXOFFS = .FALSE.
   DYNDEP = .FALSE.   
!   LWDATE = 0                                                 
   NWAMN = 0
   MXOUTAR = 0
!
   FBS%NBS = -999                                                      
!
!  Set the defaults for the MDIA:                                      

   MDIA  = 6                                                           
   ALLOCATE(LAMBDA(MDIA),CNL4_1(MDIA),CNL4_2(MDIA))                   
   LAMBDA = (/0.08,0.09,0.11,0.15,0.16,0.29/)                         
   CNL4_1 = (/8.77,-13.82,10.02,-15.92,14.41,0.65/)                    
   CNL4_2 = CNL4_1                                                     
   CNL4_1 = CNL4_1 * ((2.*PI_W)**9)                                      
   CNL4_2 = CNL4_2 * ((2.*PI_W)**9)                                      
!
!  *** Initial conditions ***
   ICOND = 0                                                          
!
   BNAUT  = .FALSE.                                                    
   BNDCHK = .TRUE.                                                     
   BRESCL = .TRUE.                                                     
   ONED   = .FALSE.                                                    
   ACUPDA = .TRUE.                                                     
   OFFSRC = .FALSE.                                                 !   40.80
   LADDS  = .TRUE.  ! TrUE to get wave dissipation etc              !   40.85
   HSRERR = 0.1_rkind                                                        
!
!  higher order propagation and spherical coordinates                  
!
   PROJ_METHOD = 0                                                     
   PROPSS = 2                                                          
   PROPSN = 3                                                          
   PROPSC = 1                                                          
   PROPSL = 1                                                          
   PROPFL = 0                                                          
   WAVAGE = 0._rkind                                                         
   KSPHER = 0                                                          
   KREPTX = 0                                                          
   REARTH2 = 2.E7/PI_W    !Jianzhong Ge                                                  
   LENDEG = 2.E7/180._rkind                                                  
!
!  *** setup flag ***                                                  
   LSETUP = 0                                                          
!
!  *** flag for setup convergence                                      
!
   CSETUP = .TRUE.                                                     

!     flag for frequency dependent surf breaking                          41.06
!
   IFRSRF = 0                                                            !41.06
!
!     flag for wave directionality in surf breaking                       41.47
!
   IDISRF = 0   
!
!  PSETUP(1) is currently unused, but can be used as setup nesting flag
!  PSETUP(2) is the user defined correction for the level of the setup
!
   PSETUP(1) = 0.0_rkind                                                     
   PSETUP(2) = 0.0_rkind                                                     
!
!  *** ACCURACY criterion ***
!
!  *** relative error in significant wave height and mean period ***
!   PNUMS(1)  = 0.02              
   PNUMS(1)  = 0.02_rkind !JL Update
                                      
!  *** absolute error in significant wave heigth (m) ***
!   PNUMS(2)  = 0.03
   PNUMS(2)  = 0.005_rkind  ! JL update
!  *** absolute error in mean wave period (s) ***
!   PNUMS(3)  = 0.3
   PNUMS(3)  = 1000._rkind   !JL update
!  *** total number of wet gridpoints were accuracy has ***
!  *** been reached                                     ***
   PNUMS(4)  = 98.00_rkind
!
!  *** DIFFUSION schemes ***
!
!  *** Numerical diffusion over theta ***
   PNUMS(6)  = 0.5_rkind                                                     
!  *** Numerical diffusion over sigma ***
   PNUMS(7)  = 0.5_rkind                                                     
!  *** Explicit or implicit scheme in frequency space ***
!  *** default = implicit : PNUMS(8) = 1              ***
   PNUMS(8) = 1._rkind
!  *** diffusion coefficient for explicit scheme ***
   PNUMS(9) = 0.01_rkind
!
!  *** parameters for the SIP solver                        ***
!
!  *** Required accuracy to terminate the solver            ***
!  ***                                                      ***
!  ***  || Ax-b ||  <  eps2 * || b ||                       ***
!  ***                                                      ***
!  ***  eps2 = PNUMS(12)                                    ***
!
!  *** PNUMS(13) output for the solver. Possible values:    ***
!  ***     <0  : no output                                  ***
!  ***      0  : only fatal errors will be printed          ***
!  ***      1  : additional information about the iteration ***
!  ***           is printed                                 ***
!  ***      2  : gives a maximal amount of output           ***
!  ***           concerning the iteration process           ***
!
!  *** PNUMS(14) : maximum number of iterations             ***
!
   PNUMS(12) = 1.E-4
   PNUMS(13) = 0._rkind
   PNUMS(14) = 20._rkind
!
!  For the setup calculation, next parameters for the solver are used:
!
!  PNUMS(23) : required accuracy to terminate the solver
!  PNUMS(24) : output for the solver (see PNUMS(13) for meanings)
!  PNUMS(25) : maximum number of iterations
!
   PNUMS(23) = 1.E-6                                                   
   PNUMS(24) = 0._rkind                                                      
   PNUMS(25) = 1000._rkind                                                   
!
!  Maximum growth in spectral bin
!  The value is the default in the command GEN3 KOM
!
   PNUMS(20) = 0.1_rkind                                                     
!
!  Added coefficient for use with limiter on action (Qb switch)        
!
   PNUMS(28) = 1._rkind                                                      
!
!  *** set the values of PNUMS that are not used equal 0. ***
!
   PNUMS(5)  = 0.
!
!  The allowed global errors in the iteration procedure:               
!  PNUMS(15) for Hs and PNUMS(16) for Tm01                             
!
!      PNUMS(15) = 0.02                                                    30.82
!      PNUMS(16) = 0.02                                                    30.82
!     The next two values are meant for STOPC command
   PNUMS(15) = 0.005_rkind
   PNUMS(16) = 1000._rkind
!
!  coefficient for limitation of Ctheta                                
!  default no limitation on refraction                                 
!
   PNUMS(17) = -1._rkind                                                     
!
!  Limitation on Froude number; current velocity is reduced if greater
!  than Pnums(18)*Sqrt(grav*depth)
!
   !PNUMS(18) = 0.8                                                     
!  JEROME : limit ambient current over New caledonia reef, where the grid size
!  is small and current strong
   PNUMS(18) = 0.5_rkind

!
!  *** CFL criterion for explicit scheme in frequency space ***
!
   PNUMS(19) = 0.5_rkind * sqrt (2._rkind)
!
!  --- coefficient for type stopping criterion                         
!
!  PNUMS(21) = 0.          (ACCUR)                                            
   PNUMS(21) = 1._rkind  ! JL =1 (STOPC) in version 41.1
!
!  --- under-relaxation factor
!
   PNUMS(30) = 0.00_rkind                                                    
!
!     --- parameters for limiting Ctheta
!
   PNUMS(26) = 0.2_rkind                                                  !   41.06
   PNUMS(27) = 2.0_rkind                                                  !   41.06
   PNUMS(29) = 0.0_rkind                                                  !   41.06

!     --- computation of Ctheta based on wave number
!
   PNUMS(32) = 1._rkind                                                  !   41.07
!
!     --- parameters for limiting Csigma and Ctheta
!
   PNUMS(33) = 0.0_rkind                                                  !   41.35
   PNUMS(34) = 0.5_rkind                                                  !   41.35
   PNUMS(35) = 0.0_rkind                                                  !   41.35
   PNUMS(36) = 0.5_rkind                                                  !   41.35
!
!  *** (1) and (2): Komen et al. (1984) formulation ***
!
   PWCAP(1)  = 2.36E-5
   PWCAP(2)  = 3.02E-3
   PWCAP(9)  = 2._rkind                                                      
!     note that delta has been set to 1 since version 40.91A            !  41.41
   PWCAP(10) = 1._rkind                                                    !  41.41 34.00
   PWCAP(11) = 1._rkind                                                      
!
!  *** (3): Coefficient for Janssen(1989,1991) formulation ***
!  ** according to Komen et al. (1994) ***
!
   PWCAP(3)  = 4.5_rkind
   PWCAP(4)  = 0.5_rkind
!
!  *** (5): Coefficient for Longuet-Higgins ***
!
   PWCAP(5) = 1._rkind
!
!  *** (6): ALPHA in Battjes/Janssen ***
!
   PWCAP(6) = 0.88_rkind
   PWCAP(7) = 1._rkind
   PWCAP(8) = 0.75_rkind
!
!  *** (12): Proportionality coefficient in cumulative steepness       
!            method                                                    
!      (13): Power of cosine term in cumulative steepness method       
!
!   PWCAP(12) = 0.5       ! JL not used in 41.10                                        
!   PWCAP(13) = 2.0       ! JL not used in 41.10    
!
!  PBOT(1)   = 0.005          modified
   PBOT(1)   = 0.0_rkind                                                     
   PBOT(2)   = 0.015_rkind                                                   
   PBOT(3)   = 0.038_rkind
   PBOT(4)   = -0.08_rkind
   PBOT(5)   = 0.05_rkind
   PBOT(6)   = 2.65_rkind                                              !  41.51
   PBOT(7)   = 0.0001_rkind                                            !  41.51

!
   PSURF(1)  = 1.0_rkind
   PSURF(2)  = 0.73_rkind                                               

!
   PMUD(1)   = 0._rkind                                                  !  40.59
   PMUD(2)   = 1300._rkind                                               !  40.59
   PMUD(3)   = 0.0076_rkind                                              !  40.59
   PMUD(4)   = RHO_W                                               !  40.59
   PMUD(5)   = 1.3E-6                                              !  40.59
     
!
!  triad interactions
!   PTRIAD(1)  = 0.1                                                    
!   PTRIAD(2)  = 5.0                                                    
!   PTRIAD(3)  = 10.                                                    
!   PTRIAD(4)  = 0.2                                                    
!   PTRIAD(5)  = 0.01                                                   
   PTRIAD(1)  = 0.8_rkind                                                 !  41.44 40.61 30.82
   PTRIAD(2)  = 2.5_rkind                                                 !  40.61 40.56 30.82
   PTRIAD(3)  = 10._rkind                                                 !  40.23
   PTRIAD(4)  = 0.2_rkind                                                 !  40.13
   PTRIAD(5)  = 0.01_rkind                                                !  40.23
   PTRIAD(6)  = 0.95_rkind                                                !  41.46
!     original value of B=-0.75 in SPB appears to be reasonable
!     for unidirectional laboratory cases
!     however, this choice may not be appropriate for true
!     2D cases as it does not scale with the wave field
!     thus, B = 0 (for now)
!      PTRIAD(7)  = -0.75                                                  41.46
   PTRIAD(7)  = 0._rkind                                                  !  41.46
   PTRIAD(8)  = 0._rkind                                                  !  41.46
!
!  quadruplet interactions
   PQUAD(1) = 0.25_rkind                                                     
   PQUAD(2) = 3.E7                                                     
   PQUAD(3) = 5.5_rkind                                                      
   PQUAD(4) = 0.833_rkind                                                    
   PQUAD(5) = -1.25_rkind                                                    
!
   PWIND(1)  = 188.0_rkind
   PWIND(2)  = 0.59_rkind
   PWIND(3)  = 0.12_rkind
   PWIND(4)  = 250.0_rkind
   PWIND(5)  = 0.0023_rkind
   PWIND(6)  = -0.223_rkind
   PWIND(7)  = 0._rkind
   PWIND(8)  = -0.56_rkind
   PWIND(10) = 0.0036_rkind
   PWIND(11) = 0.00123_rkind
   PWIND(12) = 1.0_rkind
   PWIND(13) = 0.13_rkind
!  *** Janssen (1991) wave growth model ***
!  *** alpha ***
   PWIND(14) = 0.01_rkind
!   PWIND(14) = 0.0144
!  *** Charnock: Von Karman constant ***
   PWIND(15) = 0.41_rkind
!  *** rho air (density) ****
   PWIND(16) = 1.28_rkind
!  *** rho water (density) ***
   PWIND(17) = RHO_W
   PWIND(9)  = PWIND(16) / RHO_W
!  Coefficient in front of A term in 3d gen. growth term
!  default is 0; can be made non-zero in command GEN3 or GROWTH
   PWIND(31) = 0._rkind                                                      
!
!  --- coefficients for diffraction approximation                      
!
   IDIFFR    = 0                                                       
   PDIFFR(:) = 0._rkind                                                      
!
!  pointers in array COMPDA
   JDISS = 2
   JUBOT = 3
   JQB   = 4
   JSTP  = 5
   JDHS  = 6
   JDP1  = 7
   JDP2  = 8
   JVX1  = 9
   JVY1  = 10
   JVX2  = 11
   JVY2  = 12
   JVX3  = 13
   JVY3  = 14
   JDP3  = 15
   JWX2  = 16                                                         ! 40.00
   JWY2  = 17                                                         ! 40.00
   JWX3  = 18                                                         ! 40.00
   JWY3  = 19                                                         ! 40.00
   JDTM  = 20                                                         ! 40.00
   JLEAK = 21                                                         ! 40.00
   JWLV1 = 22                                                         ! 40.00
   JWLV3 = 23                                                         ! 40.00
   JWLV2 = 24                                                         ! 40.00
   JHSIBC = 25                                                        ! 40.00
   JHS    = 26                                                        ! 40.13
   JURSEL = 27                                                        ! 40.13
   JBOTLV = 28

   MCMVAR = 28                                   !41.38 40.65 40.61 40.51 40.13
!  subarray sequence number 1 is used only for unused subarrays        
   JFRC2 = 1                                                           
   JFRC3 = 1                                                           
   JUSTAR= 1
!  JL JCHARN Added ! Charnock and USTAR will be used for Air/sea coupling 
!  (need IWIND=4 and IWCAP = 2, following the JANSSEN et al. parameterization)
   JCHARN= 1
   JZEL  = 1                                                           
   JTAUW = 1                                                           
   JCDRAG= 1                                                           
!  added for air-sea temp. diff.:
   JASTD2= 1                                                           
   JASTD3= 1                                                         
!     added for fluid mud layer                                           40.59
   JMUDL1= 1                                                         ! 40.59
   JMUDL2= 1                                                         ! 40.59
   JMUDL3= 1                                                         ! 40.59
!     added for number of plants per square meter                         40.55
   JNPLA2= 1                                                         ! 40.55
   JNPLA3= 1                                                         ! 40.55
!     added for turbulent viscosity                                       40.35
   JTURB2= 1                                                         ! 40.35
   JTURB3= 1                                                         ! 40.35
   JGAMMA = 1                                                        ! 41.38
   JRESPL = 1                                                        ! 41.38
!
!  *** added for wave setup ***                                        
!
   JSETUP = 1                                                          
   JDPSAV = 1 
!
!     --- added for output purposes                                       40.65
!
   JPBOT  = 1                                                         ! 40.65 40.51
   JDSXB  = 1                                                         ! 40.65 40.61
   JDSXS  = 1                                                         ! 40.65 40.61
   JDSXW  = 1                                                         ! 40.65 40.61
   JDSXM  = 1                                                         ! 40.65 40.61
   JDSXV  = 1                                                         ! 40.65 40.61
   JDSXT  = 1                                                         ! 40.35
   JGENR  = 1                                                         ! 40.85
   JGSXW  = 1                                                         ! 40.85
   JREDS  = 1                                                         ! 40.85
   JRSXQ  = 1                                                         ! 40.85
   JRSXT  = 1                                                         ! 40.85
   JTRAN  = 1                                                         ! 40.85
   JTSXG  = 1                                                         ! 40.85
   JTSXT  = 1                                                         ! 40.85
   JTSXS  = 1                                                         ! 40.85
   JRADS  = 1                                                         ! 40.85
                                                         
!
!  Next pointers are for the plot of the source terms (SWTSDA)         
!
   JPWNDA = 1
   JPWNDB = 2
   JPWCAP = 3
   JPBTFR = 4
   JPWBRK = 5
   JP4S   = 6
   JP4D   = 7
   JPTRI  = 8
   JPVEGT = 9                                                         ! 40.55
   JPTURB = 10                                                        ! 40.35
   JPMUD  = 11                                                        ! 40.59
   MTSVAR = 11
!
!  ***** test output control *****
   ITEST  = 1                                                          
   INTES  = 0                                                          
   ICOTES = 0                                                          
   IOUTES = 0                                                          
   LTRACE = .FALSE.
   TESTFL = .FALSE.
!   NPTST  = 0
   NPTST = 1 ! JL should GE 1
   NPTSTA = 1
   LXDMP  = -1
   LYDMP  = 0
   NEGMES = 0
   MAXMES = 200
   IFPAR = 0                                                           
   IFS1D = 0                                                           
   IFS2D = 0                                                           
!  number of obstacles initialised at 0                                
   NUMOBS = 0                                                          
!  ***** output *****
   IUBOTR = 0
!  ***** plot output *****
!
   DO IVT = 1, NMOVAR
     OVKEYW(IVT) = 'XXXX'                                              
   ENDDO

   !     properties of output variables
!
      IVTYPE = 1
!     keyword used in SWAN command
      OVKEYW(IVTYPE) = 'XP'                                              !40.00
!     short name
      OVSNAM(IVTYPE) = 'Xp'
!     long name
      OVLNAM(IVTYPE) = 'X user coordinate'
!     unit name
      OVUNIT(IVTYPE) = UL
!     type (scalar/vector etc.)
      OVSVTY(IVTYPE) = 1
!     lower and upper limit
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
!     lowest and highest expected value
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
!     exception value
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 2
      OVKEYW(IVTYPE) = 'YP'                                              !40.00
      OVSNAM(IVTYPE) = 'Yp'
      OVLNAM(IVTYPE) = 'Y user coordinate'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 3
      OVKEYW(IVTYPE) = 'DIST'                                            !40.00
      OVSNAM(IVTYPE) = 'Dist'
      OVLNAM(IVTYPE) = 'distance along output curve'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 4
      OVKEYW(IVTYPE) = 'DEP'                                             !40.00
      OVSNAM(IVTYPE) = 'Depth'
      OVLNAM(IVTYPE) = 'Depth'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 5
      OVKEYW(IVTYPE) = 'VEL'                                             !40.00
      OVSNAM(IVTYPE) = 'Vel'
      OVLNAM(IVTYPE) = 'Current velocity'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -2.
      OVHEXP(IVTYPE) = 2.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 6
      OVKEYW(IVTYPE) = 'UBOT'                                            !40.00
      OVSNAM(IVTYPE) = 'Ubot'
      OVLNAM(IVTYPE)='RMS of maxima of orbital velocity near the bottom'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -10.
!
      IVTYPE = 7
      OVKEYW(IVTYPE) = 'DISS'                                            !40.00
      OVSNAM(IVTYPE) = 'Dissip'
      OVLNAM(IVTYPE) = 'Energy dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 8
      OVKEYW(IVTYPE) = 'QB'                                              !40.00
      OVSNAM(IVTYPE) = 'Qb'
      OVLNAM(IVTYPE) = 'Fraction breaking waves'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -1.
!
      IVTYPE = 9
      OVKEYW(IVTYPE) = 'LEA'                                             !40.00
      OVSNAM(IVTYPE) = 'Leak'
      OVLNAM(IVTYPE) = 'Energy leak over spectral boundaries'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 10
      OVKEYW(IVTYPE) = 'HS'                                              !40.00
      OVSNAM(IVTYPE) = 'Hsig'
      OVLNAM(IVTYPE) = 'Significant wave height'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = -9.
!                                                                         modified 10.09
      IVTYPE = 11
      OVKEYW(IVTYPE) = 'TM01'                                            !40.00
      OVSNAM(IVTYPE) = 'Tm01'                                            !20.81
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 12
      OVKEYW(IVTYPE) = 'RTP'                                             !40.00
      OVSNAM(IVTYPE) = 'RTpeak'                                          !40.41 20.81
      OVLNAM(IVTYPE) = 'Relative peak period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 13
      OVKEYW(IVTYPE) = 'DIR'                                             !40.00
      OVSNAM(IVTYPE) = 'Dir'
      OVLNAM(IVTYPE) = 'Average wave direction'
      OVUNIT(IVTYPE) = UDI                                               !30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 14
      OVKEYW(IVTYPE) = 'PDI'                                             !40.00
      OVSNAM(IVTYPE) = 'PkDir'
      OVLNAM(IVTYPE) = 'direction of the peak of the spectrum'
      OVUNIT(IVTYPE) = UDI                                               !30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 15
      OVKEYW(IVTYPE) = 'TDI'                                             !40.00
      OVSNAM(IVTYPE) = 'TDir'
      OVLNAM(IVTYPE) = 'direction of the energy transport'
      OVUNIT(IVTYPE) = UDI                                               !30.72
      OVSVTY(IVTYPE) = 2
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 360.
      OVEXCV(IVTYPE) = -999.
!
      IVTYPE = 16
      OVKEYW(IVTYPE) = 'DSPR'                                            !40.00
      OVSNAM(IVTYPE) = 'Dspr'
      OVLNAM(IVTYPE) = 'directional spreading'
      OVUNIT(IVTYPE) = UDI                                               !30.72
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 360.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 60.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 17
      OVKEYW(IVTYPE) = 'WLEN'                                            !40.00
      OVSNAM(IVTYPE) = 'Wlen'
      OVLNAM(IVTYPE) = 'Average wave length'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 200.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 18
      OVKEYW(IVTYPE) = 'STEE'                                            !40.00
      OVSNAM(IVTYPE) = 'Steepn'
      OVLNAM(IVTYPE) = 'Wave steepness'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 19
      OVKEYW(IVTYPE) = 'TRA'                                             !40.00
      OVSNAM(IVTYPE) = 'Transp'
      OVLNAM(IVTYPE) = 'Wave energy transport'
      OVUNIT(IVTYPE) = 'm3/s'                                            !40.00
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 20
      OVKEYW(IVTYPE) = 'FOR'                                             !40.00
      OVSNAM(IVTYPE) = 'WForce'
      OVLNAM(IVTYPE) = 'Wave driven force per unit surface'
      OVUNIT(IVTYPE) = UF                                                !30.72
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -1.E5
      OVULIM(IVTYPE) =  1.E5
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -9.                                               !40.80
!
      IVTYPE = 21
      OVKEYW(IVTYPE) = 'AAAA'                                            !40.00
      OVSNAM(IVTYPE) = 'AcDens'
      OVLNAM(IVTYPE) = 'spectral action density'
      OVUNIT(IVTYPE) = 'm2s'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 22
      OVKEYW(IVTYPE) = 'EEEE'                                            !40.00
      OVSNAM(IVTYPE) = 'EnDens'
      OVLNAM(IVTYPE) = 'spectral energy density'
      OVUNIT(IVTYPE) = 'm2'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 23
      OVKEYW(IVTYPE) = 'AAAA'                                            !40.00
      OVSNAM(IVTYPE) = 'Aux'
      OVLNAM(IVTYPE) = 'auxiliary variable'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E10
      OVULIM(IVTYPE) = 1.E10
      OVLEXP(IVTYPE) = -1.E10
      OVHEXP(IVTYPE) = 1.E10
      OVEXCV(IVTYPE) = -1.E10
!
      IVTYPE = 24
      OVKEYW(IVTYPE) = 'XC'                                              !40.00
      OVSNAM(IVTYPE) = 'Xc'
      OVLNAM(IVTYPE) = 'X computational grid coordinate'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 25
      OVKEYW(IVTYPE) = 'YC'                                              !40.00
      OVSNAM(IVTYPE) = 'Yc'
      OVLNAM(IVTYPE) = 'Y computational grid coordinate'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 26
      OVKEYW(IVTYPE) = 'WIND'                                            !40.00
      OVSNAM(IVTYPE) = 'Windv'
      OVLNAM(IVTYPE) = 'Wind velocity at 10 m above sea level'
      OVUNIT(IVTYPE) = UV                                                !30.72
      OVSVTY(IVTYPE) = 3
      OVLLIM(IVTYPE) = -100.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = -50.
      OVHEXP(IVTYPE) = 50.
      OVEXCV(IVTYPE) = 0.
!
      IVTYPE = 27
      OVKEYW(IVTYPE) = 'FRC'                                             !40.00
      OVSNAM(IVTYPE) = 'FrCoef'
      OVLNAM(IVTYPE) = 'Bottom friction coefficient'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!                                                                     new 10.09
      IVTYPE = 28
      OVKEYW(IVTYPE) = 'RTM01'                                           !40.00
      OVSNAM(IVTYPE) = 'RTm01'                                           !20.81
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 29                                                        !20.28
      OVKEYW(IVTYPE) = 'EEEE'                                            !40.00
      OVSNAM(IVTYPE) = 'EnDens'
      OVLNAM(IVTYPE) = 'energy density integrated over direction'        !20.28
      OVUNIT(IVTYPE) = 'm2'
      OVSVTY(IVTYPE) = 5
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 30                                                        !20.52
      OVKEYW(IVTYPE) = 'DHS'                                             !40.00
      OVSNAM(IVTYPE) = 'dHs'
      OVLNAM(IVTYPE) = 'difference in Hs between iterations'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 31                                                        !20.52
      OVKEYW(IVTYPE) = 'DRTM01'                                          !40.00
      OVSNAM(IVTYPE) = 'dTm'
      OVLNAM(IVTYPE) = 'difference in Tm between iterations'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 2.
      OVEXCV(IVTYPE) = -9.
!                                                                        !20.61
      IVTYPE = 32
      OVKEYW(IVTYPE) = 'TM02'                                            !40.00
      OVSNAM(IVTYPE) = 'Tm02'
      OVLNAM(IVTYPE) = 'Zero-crossing period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!                                                                        !20.61
      IVTYPE = 33
      OVKEYW(IVTYPE) = 'FSPR'                                            !40.00
      OVSNAM(IVTYPE) = 'FSpr'                                            !20.67
      OVLNAM(IVTYPE) = 'Frequency spectral width (Kappa)'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 34                                                        !20.67
      OVKEYW(IVTYPE) = 'URMS'                                            !40.00
      OVSNAM(IVTYPE) = 'Urms'
      OVLNAM(IVTYPE) = 'RMS of orbital velocity near the bottom'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 35                                                        !30.22
      OVKEYW(IVTYPE) = 'UFRI'                                            !40.00
      OVSNAM(IVTYPE) = 'Ufric'
      OVLNAM(IVTYPE) = 'Friction velocity'
      OVUNIT(IVTYPE) = UV
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 36                                                        !30.22
      OVKEYW(IVTYPE) = 'ZLEN'                                            !40.00
      OVSNAM(IVTYPE) = 'Zlen'
      OVLNAM(IVTYPE) = 'Zero velocity thickness of boundary layer'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 37                                                        !30.22
      OVKEYW(IVTYPE) = 'TAUW'                                            !40.00
      OVSNAM(IVTYPE) = 'TauW'
      OVLNAM(IVTYPE) = '    '
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 38                                                        !30.22
      OVKEYW(IVTYPE) = 'CDRAG'                                           !40.00
      OVSNAM(IVTYPE) = 'Cdrag'
      OVLNAM(IVTYPE) = 'Drag coefficient'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
!     *** wave-induced setup ***                                          32.02
!
      IVTYPE = 39                                                        !32.02
      OVKEYW(IVTYPE) = 'SETUP'                                           !40.00
      OVSNAM(IVTYPE) = 'Setup'                                           !32.02
      OVLNAM(IVTYPE) = 'Setup due to waves'                              !32.02
      OVUNIT(IVTYPE) = 'm'                                               !32.02
      OVSVTY(IVTYPE) = 1                                                 !32.02
      OVLLIM(IVTYPE) = -1.                                               !32.02
      OVULIM(IVTYPE) = 1.                                                !32.02
      OVLEXP(IVTYPE) = -1.                                               !32.02
      OVHEXP(IVTYPE) = 1.                                                !32.02
      OVEXCV(IVTYPE) = -9.                                               !32.02
!
      IVTYPE = 40                                                       !40.00
      OVKEYW(IVTYPE) = 'TIME'                                            !40.00
      OVSNAM(IVTYPE) = 'Time'
      OVLNAM(IVTYPE) = 'Date-time'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -99999.
!
      IVTYPE = 41                                                       !40.00
      OVKEYW(IVTYPE) = 'TSEC'                                            !40.00
      OVSNAM(IVTYPE) = 'Tsec'
      OVLNAM(IVTYPE) = 'Time in seconds from reference time'
      OVUNIT(IVTYPE) = 's'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100000.
      OVLEXP(IVTYPE) = -100000.
      OVHEXP(IVTYPE) = 1000000.
      OVEXCV(IVTYPE) = -99999.
!                                                        new             !40.00
      IVTYPE = 42
      OVKEYW(IVTYPE) = 'PER'                                             !40.00
      OVSNAM(IVTYPE) = 'Period'                                          !40.00
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!                                                        new             !40.00
      IVTYPE = 43
      OVKEYW(IVTYPE) = 'RPER'                                            !40.00
      OVSNAM(IVTYPE) = 'RPeriod'                                         !40.41!40.00
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 44                                                        !40.00
      OVKEYW(IVTYPE) = 'HSWE'                                            !40.00
      OVSNAM(IVTYPE) = 'Hswell'
      OVLNAM(IVTYPE) = 'Wave height of swell part'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 100.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 10.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 45
      OVKEYW(IVTYPE) = 'URSELL'                                          !40.03
      OVSNAM(IVTYPE) = 'Ursell'
      OVLNAM(IVTYPE) = 'Ursell number'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 46
      OVKEYW(IVTYPE) = 'ASTD'                                            !40.03
      OVSNAM(IVTYPE) = 'ASTD'
      OVLNAM(IVTYPE) = 'Air-Sea temperature difference'
      OVUNIT(IVTYPE) = 'K'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -50.
      OVULIM(IVTYPE) =  50.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 47                                                        !40.41
      OVKEYW(IVTYPE) = 'TMM10'
      OVSNAM(IVTYPE) = 'Tm_10'
      OVLNAM(IVTYPE) = 'Average absolute wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 48                                                        !40.41
      OVKEYW(IVTYPE) = 'RTMM10'
      OVSNAM(IVTYPE) = 'RTm_10'
      OVLNAM(IVTYPE) = 'Average relative wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 49
      OVKEYW(IVTYPE) = 'DIFPAR'                                          !40.21
      OVSNAM(IVTYPE) = 'DifPar'
      OVLNAM(IVTYPE) = 'Diffraction parameter'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -50.
      OVULIM(IVTYPE) =  50.
      OVLEXP(IVTYPE) = -10.
      OVHEXP(IVTYPE) =  10.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 50                                                        !40.51
      OVKEYW(IVTYPE) = 'TMBOT'
      OVSNAM(IVTYPE) = 'TmBot'
      OVLNAM(IVTYPE) = 'Near bottom wave period'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 51                                                        !40.51
      OVKEYW(IVTYPE) = 'WATL'
      OVSNAM(IVTYPE) = 'Watlev'
      OVLNAM(IVTYPE) = 'Water level'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 52                                                        !40.51
      OVKEYW(IVTYPE) = 'BOTL'
      OVSNAM(IVTYPE) = 'Botlev'
      OVLNAM(IVTYPE) = 'Bottom level'
      OVUNIT(IVTYPE) = UH
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = -1.E4
      OVULIM(IVTYPE) = 1.E4
      OVLEXP(IVTYPE) = -100.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 53                                                        !40.51
      OVKEYW(IVTYPE) = 'TPS'
      OVSNAM(IVTYPE) = 'TPsmoo'
      OVLNAM(IVTYPE) = 'Relative peak period (smooth)'
      OVUNIT(IVTYPE) = UT
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 54
      OVKEYW(IVTYPE) = 'DISB'                                            !40.61
      OVSNAM(IVTYPE) = 'Sfric'                                           !40.85
      OVLNAM(IVTYPE) = 'Bottom friction dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 55
      OVKEYW(IVTYPE) = 'DISSU'                                           !40.61
      OVSNAM(IVTYPE) = 'Ssurf'                                           !40.85
      OVLNAM(IVTYPE) = 'Surf breaking dissipation'                       !40.85
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 56
      OVKEYW(IVTYPE) = 'DISW'                                            !40.61
      OVSNAM(IVTYPE) = 'Swcap'                                           !40.85
      OVLNAM(IVTYPE) = 'Whitecapping dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 57
      OVKEYW(IVTYPE) = 'DISV'                                            !40.61
      OVSNAM(IVTYPE) = 'Sveg'                                            !40.85
      OVLNAM(IVTYPE) = 'Vegetation dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 58                                                        !40.64
      OVKEYW(IVTYPE) = 'QP'
      OVSNAM(IVTYPE) = 'Qp'
      OVLNAM(IVTYPE) = 'Peakedness'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 59                                                        !40.64
      OVKEYW(IVTYPE) = 'BFI'
      OVSNAM(IVTYPE) = 'BFI'
      OVLNAM(IVTYPE) = 'Benjamin-Feir index'
      OVUNIT(IVTYPE) = ' '
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 1000.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 60
      OVKEYW(IVTYPE) = 'GENE'                                            !40.85
      OVSNAM(IVTYPE) = 'Genera'
      OVLNAM(IVTYPE) = 'Energy generation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 61
      OVKEYW(IVTYPE) = 'GENW'                                            !40.85
      OVSNAM(IVTYPE) = 'Swind'
      OVLNAM(IVTYPE) = 'Wind source term'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 62
      OVKEYW(IVTYPE) = 'REDI'                                            !40.85
      OVSNAM(IVTYPE) = 'Redist'
      OVLNAM(IVTYPE) = 'Energy redistribution'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 63
      OVKEYW(IVTYPE) = 'REDQ'                                            !40.85
      OVSNAM(IVTYPE) = 'Snl4'
      OVLNAM(IVTYPE) = 'Total absolute 4-wave interaction'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 64
      OVKEYW(IVTYPE) = 'REDT'                                            !40.85
      OVSNAM(IVTYPE) = 'Snl3'
      OVLNAM(IVTYPE) = 'Total absolute 3-wave interaction'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 65
      OVKEYW(IVTYPE) = 'PROPA'                                           !40.85
      OVSNAM(IVTYPE) = 'Propag'
      OVLNAM(IVTYPE) = 'Energy propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 66
      OVKEYW(IVTYPE) = 'PROPX'                                           !40.85
      OVSNAM(IVTYPE) = 'Propxy'
      OVLNAM(IVTYPE) = 'xy-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 67
      OVKEYW(IVTYPE) = 'PROPT'                                           !40.85
      OVSNAM(IVTYPE) = 'Propth'
      OVLNAM(IVTYPE) = 'theta-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 68
      OVKEYW(IVTYPE) = 'PROPS'                                           !40.85
      OVSNAM(IVTYPE) = 'Propsi'
      OVLNAM(IVTYPE) = 'sigma-propagation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 69
      OVKEYW(IVTYPE) = 'RADS'                                            !40.85
      OVSNAM(IVTYPE) = 'Radstr'
      OVLNAM(IVTYPE) = 'Radiation stress'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 70
      OVKEYW(IVTYPE) = 'NPL'                                             !41.12
      OVSNAM(IVTYPE) = 'Nplant'
      OVLNAM(IVTYPE) = 'Plants per m2'
      OVUNIT(IVTYPE) = '1/m2'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 71
      OVKEYW(IVTYPE) = 'LWAVP'                                           !41.15
      OVSNAM(IVTYPE) = 'Lwavp'
      OVLNAM(IVTYPE) = 'Peak wave length'
      OVUNIT(IVTYPE) = UL
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 200.
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 72
      OVKEYW(IVTYPE) = 'DISTU'
      OVSNAM(IVTYPE) = 'Stur'
      OVLNAM(IVTYPE) = 'Turbulent dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 73
      OVKEYW(IVTYPE) = 'TURB'
      OVSNAM(IVTYPE) = 'Turb'
      OVLNAM(IVTYPE) = 'Turbulent viscosity'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 10.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 100.
      OVEXCV(IVTYPE) = -99.
!
      IVTYPE = 74
      OVKEYW(IVTYPE) = 'DISM'                                            !40.61
      OVSNAM(IVTYPE) = 'Smud'                                            !40.85
      OVLNAM(IVTYPE) = 'Fluid mud dissipation'
      OVUNIT(IVTYPE) = 'm2/s'
      OVSVTY(IVTYPE) = 1
      OVLLIM(IVTYPE) = 0.
      OVULIM(IVTYPE) = 1000.
      OVLEXP(IVTYPE) = 0.
      OVHEXP(IVTYPE) = 0.1
      OVEXCV(IVTYPE) = -9.
!
      IVTYPE = 100                                                       !41.62
      OVKEYW(IVTYPE) = 'PTHS'                                            !41.62
      OVSNAM(IVTYPE) = 'HsPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 01'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 101                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 02'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 102                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 03'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 103                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 04'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 104                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 05'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 105                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 06'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 106                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 07'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 107                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 08'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 108                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 09'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 109                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10HS'                                          !41.62
      OVSNAM(IVTYPE) = 'HsPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave height of partition 10'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 110                                                       !41.62  !make sure these defs are SWAN-consistent!
      OVKEYW(IVTYPE) = 'PTRTP'                                           !41.62
      OVSNAM(IVTYPE) = 'TpPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 01'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 111                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 02'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 112                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 03'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 113                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 04'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 114                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 05'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 115                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 06'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 116                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 07'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 117                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 08'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 118                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 09'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 119                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10RTP'                                         !41.62
      OVSNAM(IVTYPE) = 'TpPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Relative peak period of partition 10'            !41.62
      OVUNIT(IVTYPE) = UT                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 100.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 120                                                       !41.62
      OVKEYW(IVTYPE) = 'PTWLEN'                                          !41.62
      OVSNAM(IVTYPE) = 'WlPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 01'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 121                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 02'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 122                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 03'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 123                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 04'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 124                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 05'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 125                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 06'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 126                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 07'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 127                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 08'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 128                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 09'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 129                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10WLEN'                                        !41.62
      OVSNAM(IVTYPE) = 'WlPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave length of partition 10'             !41.62
      OVUNIT(IVTYPE) = UL                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1000.                                             !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 200.                                              !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 130                                                       !41.62  !make sure these defs are SWAN-consistent!
      OVKEYW(IVTYPE) = 'PTDIR'                                           !41.62
      OVSNAM(IVTYPE) = 'DrPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 01'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 131                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 02'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 132                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 03'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 133                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 04'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 134                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 05'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 135                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 06'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 136                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 07'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 137                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 08'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 138                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 09'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 139                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10DIR'                                         !41.62
      OVSNAM(IVTYPE) = 'DrPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Average wave direction of partition 10'          !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 2                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 360.                                              !41.62
      OVEXCV(IVTYPE) = -999.                                             !41.62
!
      IVTYPE = 140                                                       !41.62  !make sure these defs are SWAN-consistent!
      OVKEYW(IVTYPE) = 'PTDSPR'                                          !41.62
      OVSNAM(IVTYPE) = 'DsPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 01'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 141                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 02'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 142                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 03'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 143                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 04'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 144                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 05'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 145                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 06'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 146                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 07'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 147                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 08'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 148                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 09'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 149                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10DSPR'                                        !41.62
      OVSNAM(IVTYPE) = 'DsPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Directional spreading of partition 10'           !41.62
      OVUNIT(IVTYPE) = UDI                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 360.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 60.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 150                                                       !41.62
      OVKEYW(IVTYPE) = 'PTWFRAC'                                         !41.62
      OVSNAM(IVTYPE) = 'WfPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 01'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 151                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 02'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 152                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 03'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 153                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 04'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 154                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 05'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 155                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 06'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 156                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 07'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 157                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 08'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 158                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 09'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 159                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10WFRAC'                                       !41.62
      OVSNAM(IVTYPE) = 'WfPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Wind fraction of partition 10'                   !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 1.                                                !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 160                                                       !41.62
      OVKEYW(IVTYPE) = 'PTSTEE'                                          !41.62
      OVSNAM(IVTYPE) = 'StPT01'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 01'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 161                                                       !41.62
      OVKEYW(IVTYPE) = 'PT02STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT02'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 02'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 162                                                       !41.62
      OVKEYW(IVTYPE) = 'PT03STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT03'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 03'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 163                                                       !41.62
      OVKEYW(IVTYPE) = 'PT04STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT04'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 04'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 164                                                       !41.62
      OVKEYW(IVTYPE) = 'PT05STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT05'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 05'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 165                                                       !41.62
      OVKEYW(IVTYPE) = 'PT06STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT06'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 06'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 166                                                       !41.62
      OVKEYW(IVTYPE) = 'PT07STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT07'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 07'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 167                                                       !41.62
      OVKEYW(IVTYPE) = 'PT08STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT08'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 08'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 168                                                       !41.62
      OVKEYW(IVTYPE) = 'PT09STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT09'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 09'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 169                                                       !41.62
      OVKEYW(IVTYPE) = 'PT10STEE'                                        !41.62
      OVSNAM(IVTYPE) = 'StPT10'                                          !41.62
      OVLNAM(IVTYPE) = 'Wave steepness of partition 10'                  !41.62
      OVUNIT(IVTYPE) = ' '                                               !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 1.                                                !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 0.1                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 170                                                       !41.62
      OVKEYW(IVTYPE) = 'PARTIT'                                          !41.62
      OVSNAM(IVTYPE) = 'PARTIT'                                          !41.62
      OVLNAM(IVTYPE) = 'Spectral partions of all parameters'             !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
      IVTYPE = 171                                                       !41.62
      OVKEYW(IVTYPE) = 'NPART'                                           !41.62
      OVSNAM(IVTYPE) = 'Npart'                                           !41.62
      OVLNAM(IVTYPE) = 'Number of spectral partions'                     !41.62
      OVUNIT(IVTYPE) = UH                                                !41.62
      OVSVTY(IVTYPE) = 1                                                 !41.62
      OVLLIM(IVTYPE) = 0.                                                !41.62
      OVULIM(IVTYPE) = 100.                                              !41.62
      OVLEXP(IVTYPE) = 0.                                                !41.62
      OVHEXP(IVTYPE) = 10.                                               !41.62
      OVEXCV(IVTYPE) = -9.                                               !41.62
!
!
!  various parameters for computation of output quantities             
!
!  reference time for TSEC
   OUTPAR(1) = 0.
!  power in expression for PER and RPER
!  previous name: SPCPOW
   OUTPAR(2) = 1.
!  power in expression for WLEN
!  previous name: AKPOWR
   OUTPAR(3) = 1.
!  indicator for direction
!  =0: direction always w.r.t. user coordinates; =1: dir w.r.t. frame
   OUTPAR(4) = 0.
!  frequency limit for swell
   OUTPAR(5) = 0.1
!  0=integration over [0,inf], 1=integration over [fmin,fmax]
   OUTPAR(6:20) = 0.                                                  ! 40.87
!  lower bound of integration range for output parameters
   OUTPAR(21:35) = 0.                                                !  40.87
!  upper bound of integration range for output parameters
   OUTPAR(36:50) = 1000.                                             !  40.87

   RETURN
   END SUBROUTINE SWINIT
!
!***********************************************************************
!                                                                      *
!***********************************************************************
!                                                                      *
!      SUBROUTINE SWPREP ( BSPECS, BGRIDP, CROSS , XCGRID ,YCGRID , &      !40.31
      SUBROUTINE SWPREP (                         XCGRID ,YCGRID , &      !40.31
     &                    KGRPNT, KGRBND, SPCDIR, SPCSIG )                !40.31
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        !40.41
      USE SWCOMM1                                                         !40.41
      USE SWCOMM2                                                         !40.41
      USE SWCOMM3                                                         !40.41
      USE SWCOMM4                                                         !40.41
      USE M_BNDSPEC                                                       !40.31
!      USE M_PARALL                                                        !40.31
      USE M_DIFFR                                                         !40.21
      USE SwanGriddata                                                    !40.80

      USE VARS_WAVE, ONLY:  BSPECS, BGRIDP, CROSS
      USE schism_glbl, ONLY: rkind
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!

!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.09: Annette Kieftenburg
!     40.21: Agnieszka Herman
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     20.70, Jan. 96: new name, SPRCON is now called from this subr
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     30.70, Mar. 98: loop over grid points moved into subr SWOBST
!                     erroneous usage of SWPDIR removed
!                     close input files containing stationary input fields
!     30.82, Oct. 98: Added INTEGER declaration of array OBSTA(*)
!     40.00, Feb. 99: DYNDEP is made True, if depth or water level nonstationary
!     40.09, Aug. 00: If obstacle is on computational grid point it is moved a
!     bit
!     40.02, Oct. 00: Array KGRBND now has a dimension
!     40.02, Oct. 00: Initialisation of IERR
!     40.21, Aug. 01: allocation of arrays for diffraction
!     40.31, Nov. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Jun. 07: extension to unstructured grids
!
!  2. Purpose
!
!     do some preparations before computation is started

!
!  3. Method
!
!  4. Argument variables
!
      INTEGER, INTENT(INOUT) :: KGRBND(*)                                 !40.02
!
!     XCGRID: input  Coordinates of computational grid in x-direction     30.72
!     YCGRID: input  Coordinates of computational grid in y-direction     30.72
!
      REAL(rkind)    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                         !30.72
      REAL(rkind)    SPCDIR(MDC,6)  ,    SPCSIG(MSC)                             !40.31
!
!  5. SUBROUTINES CALLING
!
!     SWREAD
!
!  6. SUBROUTINES USED
!
!     OBSTMOVE
!     SWOBST                                                              40.09
!
!  7. ERROR MESSAGES
!
!     ---
!
!  8. REMARKS
!
!     ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Compute origin of computational grid in problem grid XPC, YPC
!       Compute origin of problem grid in computationl grid coordinates
!       Check input fields
!       Close files containing stationary input fields
!       Compute crossing of comp.grid lines with obstacles
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!      INTEGER   CROSS(2,MCGRD)
      INTEGER   KGRPNT(MXC,MYC)
!      INTEGER   BGRIDP(6*NBGRPT)                                          !40.31
!      REAL      BSPECS(MDC,MSC,NBSPEC,2)                                  !40.31
      INTEGER, ALLOCATABLE :: CROSSGL(:,:)                                !40.31
!
      TYPE(BSDAT) , POINTER :: CURRBS                                     !40.31
      TYPE(BGPDAT), POINTER :: CURBGP                                     !40.31
      SAVE IENT
      DATA IENT/0/
      CALL STRACE (IENT, 'SWPREP')
!
      IF ( ITEST.GE.100 ) THEN                                            !40.31
         WRITE (PRTEST,*) ' BNAUT after SWREAD = ',BNAUT                  !40.31
      END IF                 
!
!     coefficients for transformation from user coordinates to comp. coord.
!
      COSPC = COS(ALPC)
      SINPC = SIN(ALPC)
      XCP   = -XPC*COSPC - YPC*SINPC
      YCP   =  XPC*SINPC - YPC*COSPC
      ALCP  = -ALPC
!
!     wind direction w.r.t. computational grid
!
!     ALTMP = (WDIP - ALPC) / PI2              removed 30.50
      ALTMP = (WDIP) / PI2_W                                              !13/JAN
      WDIC  = PI2_W * (ALTMP - NINT(ALTMP))
!JL Skip Now
# if 0
      IF (MXC.LE.0 .AND. OPTG.NE.5) CALL MSGERR &                         !40.80
     & (3, 'no valid computational grid; check command CGRID')            !40.80
# endif
      IF (MCGRD.LE.1 .AND. nverts.LE.0) CALL MSGERR &                     !40.80
     & (3, 'no valid comp. grid; check command READ BOT or READ UNSTRU')  !40.80 30.50
!
      IF (LEDS(1).EQ.0) CALL MSGERR (3,'Bottom grid not defined')
      IF (LEDS(1).EQ.1) CALL MSGERR (3,'No bottom levels read')
      IF (IUBOTR.EQ.1 .AND. IBOT.EQ.0) &
     &    CALL MSGERR (1,'Bottom friction not on, UBOT not computed')
!
      IF (LEDS(2).EQ.2) THEN
        IF (LEDS(3).NE.2) &
     &  CALL MSGERR (3, 'VY not read, while VX is read')
!       ALBC  = ALPC - ALPG(2)
        ALBC  = - ALPG(2)
        COSVC = COS(ALBC)
        SINVC = SIN(ALBC)
      ENDIF
!
      IF (LEDS(4).EQ.2) VARFR = .TRUE.
!
      IF (VARFR .AND. IBOT.EQ.5) THEN                                     !41.51
         CALL MSGERR (1, &
     &           'Ripples model active, space varying friction ignored')
         VARFR = .FALSE.
      ENDIF
!
      IF (LEDS(5).EQ.2) THEN
        IF (LEDS(6).NE.2) &
     &  CALL MSGERR (3, 'WY not read, while WX is read')
        VARWI = .TRUE.
!       ALBC  = ALPC - ALPG(5)
        ALBC  = - ALPG(5)
        COSWC = COS(ALBC)
        SINWC = SIN(ALBC)
      ENDIF
!
      IF (LEDS(7).EQ.2) VARWLV = .TRUE.                                   !20.38
!
      IF (IFLDYN(1).EQ.1 .OR. IFLDYN(7).EQ.1) DYNDEP = .TRUE.             !40.00
!
      IF (OPTG.EQ.5) BNDCHK = .FALSE.                                     !40.80
!
      IF (IBOT.EQ.5) IDRAG = 1                                            !41.51

!
!     close input files containing stationary input fields                40.00
!
      DO IFLD = 1, NUMGRD
        IF (IFLDYN(IFLD).EQ.0 .AND. IFLNDS(IFLD).NE.0) THEN               !40.00
          CLOSE (IFLNDS(IFLD))
          IFLNDS(IFLD) = 0
        ENDIF
      ENDDO
!
!     computation of tail factors for moments of action spectrum
!     IP=0: action int. IP=1: energy int. IP=2: first moment of energy etc.
!
      DO IP = 0, 3
        PPTAIL = PWTAIL(1) - REAL(IP)
        PWTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))        !20.71
      ENDDO
!
!     check location of output areas
!JL Skip Now
# if 0
      CALL SPRCON ( XCGRID, YCGRID, KGRPNT, KGRBND )
# endif
!
!     *** find obstacles crossing the points in the stencil  ***
!     *** in structured grids                                ***          40.80
!
!JL Skip Now, some Job to do
# if 0
      IF (NUMOBS .GT. 0 .AND. OPTG.NE.5) THEN                             !40.80
        DO INDX = 1, MCGRD                                               !040697
          CROSS(1,INDX) = 0                                              !040697
          CROSS(2,INDX) = 0                                              !040697
        ENDDO                                                            !040697
!
        ITMP1  = MXC                                                      !40.31
        ITMP2  = MYC                                                      !40.31
        ITMP3  = MCGRD                                                    !40.31
        MXC    = MXCGL                                                    !40.31
        MYC    = MYCGL                                                    !40.31
        MCGRD  = MCGRDGL                                                  !40.31
        ALLOCATE(CROSSGL(2,MCGRDGL))                                      !40.31
        CROSSGL = 0                                                       !40.31
        CALL OBSTMOVE (XGRDGL, YGRDGL, KGRPGL)                            !40.31 40.09
        CALL SWOBST (XGRDGL, YGRDGL, KGRPGL, CROSSGL)                     !40.31 30.7N
        MXC    = ITMP1                                                    !40.31
        MYC    = ITMP2                                                    !40.31
        MCGRD  = ITMP3                                                    !40.31
!
        II = 1                                                            !40.31
        DO IX = MXF, MXL                                                  !40.31
           DO IY = MYF, MYL                                               !40.31
              INDX = KGRPGL(IX,IY)                                        !40.31
              IF ( INDX.NE.1 ) THEN                                       !40.31
                 II = II + 1                                              !40.31
                 CROSS(1:2,II) = CROSSGL(1:2,INDX)                        !40.31
              END IF                                                      !40.31
           END DO                                                         !40.31
        END DO                                                            !40.31
        DEALLOCATE(CROSSGL)   
        IF (ITEST .GE. 120) THEN                                          !060697
          WRITE(PRINTF,102)
 102      FORMAT('Links with obstacles crossing',/, &
     &    '     COMP COORD           LINK         VALUE')
          DO IIYY =2,MYC
            DO IIXX = 2,MXC
            I1 = KGRPNT(IIXX,IIYY)
              DO I2 = 1,2
                IF (CROSS(I2,I1) .NE. 0) &
     &          WRITE(PRINTF,101) IIXX,IIYY,I2,I1,CROSS(I2,I1)
 101            FORMAT(' POINT(',I4,',',I4,')', &
     &          '  CROSS(',I3,',',I5,')  = ',I5)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
# endif
!
      CURRBS => FBS                                                       !40.31
      DO                                                                  !40.31
         NBS = CURRBS%NBS                                                 !40.31
         IF (NBS.EQ.-999) EXIT                                            !40.31
         SPPARM(1:4) = CURRBS%SPPARM(1:4)                                 !40.31
         FSHAPE      = CURRBS%FSHAPE                                      !40.31
         DSHAPE      = CURRBS%DSHAPE                                      !40.31
         CALL SSHAPE (BSPECS(1,1,NBS,1), SPCSIG, SPCDIR, &                !40.31
     &                FSHAPE, DSHAPE)                                     !40.31
         IF (.NOT.ASSOCIATED(CURRBS%NEXTBS)) EXIT                         !40.31
         CURRBS => CURRBS%NEXTBS                                          !40.31
      END DO                                                              !40.31

      BGRIDP = 0                                                          !40.31
!JL Skip Now
# if 0
      IF (OPTG.NE.5) THEN                                                 !40.80
         IF (NBGGL.EQ.0) NBGGL  = NBGRPT                                  !40.51 40.31
         NBGRPT = 0                                                       !40.31
         DO IX = MXF, MXL                                                 !40.31
            DO IY = MYF, MYL                                              !40.31
               INDX = KGRPGL(IX,IY)                                       !40.31
               CURBGP => FBGP                                             !40.31
               DO II = 1, NBGGL                                           !40.31
                  INDXGR = CURBGP%BGP(1)                                  !40.31
                  IF ( INDXGR.EQ.INDX ) THEN                              !40.31
                     NBGRPT = NBGRPT + 1                                  !40.31
                     BGRIDP(6*NBGRPT-5) = KGRPNT(IX-MXF+1,IY-MYF+1)       !40.31
                     BGRIDP(6*NBGRPT-4) = CURBGP%BGP(2)                   !40.31
                     BGRIDP(6*NBGRPT-3) = CURBGP%BGP(3)                   !40.31
                     BGRIDP(6*NBGRPT-2) = CURBGP%BGP(4)                   !40.31
                     BGRIDP(6*NBGRPT-1) = CURBGP%BGP(5)                   !40.31
                     BGRIDP(6*NBGRPT  ) = CURBGP%BGP(6)                   !40.31
                  END IF                                                  !40.31

                  IF (.NOT.ASSOCIATED(CURBGP%NEXTBGP)) EXIT               !40.31
                  CURBGP => CURBGP%NEXTBGP                                !40.31
               END DO                                                     !40.31
            END DO                                                        !40.31
         END DO                                                           !40.31
      ELSE                                                                !40.80
# else
      IF (OPTG.EQ.5) THEN   ! JL Fix
# endif
         CURBGP => FBGP                                                   !40.80
         DO II = 1, NBGRPT                                                !40.80
            BGRIDP(6*II-5) = CURBGP%BGP(1)                                !40.80
            BGRIDP(6*II-4) = CURBGP%BGP(2)                                !40.80
            BGRIDP(6*II-3) = CURBGP%BGP(3)                                !40.80
            BGRIDP(6*II-2) = CURBGP%BGP(4)                                !40.80
            BGRIDP(6*II-1) = CURBGP%BGP(5)                                !40.80
            BGRIDP(6*II  ) = CURBGP%BGP(6)                                !40.80
            IF (.NOT.ASSOCIATED(CURBGP%NEXTBGP)) EXIT                     !40.80
            CURBGP => CURBGP%NEXTBGP                                      !40.80
         ENDDO                                                            !40.80
      ENDIF                                                               !40.80
!
!     --- allocate arrays for diffraction and set prop scheme to BSBT     40.41
!                                                                         40.21
      IF ( IDIFFR.EQ.1 ) THEN                                             !40.21
         IF (.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(1:MCGRD))        !40.21
         IF (.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(1:MCGRD))        !40.21
         IF (.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(1:MCGRD))        !40.21
         PROPSN = 1                                                       !40.41
         PROPSS = 1                                                       !40.41
      ELSE                                                                !40.41
         IF (.NOT.ALLOCATED(DIFPARAM)) ALLOCATE(DIFPARAM(0))              !40.21
         IF (.NOT.ALLOCATED(DIFPARDX)) ALLOCATE(DIFPARDX(0))              !40.21
         IF (.NOT.ALLOCATED(DIFPARDY)) ALLOCATE(DIFPARDY(0))              !40.21
      END IF                                                              !40.21
!
      RETURN
! * end of subroutine SWPREP *
      END

 
!**********************************************************************
!
   SUBROUTINE FLFILE_OLD (IGR1, IGR2,                                     &
                      ARR, ARR2, JX1, JX2, JX3, JY1, JY2, JY3,        &
                      COSFC, SINFC, COMPDA,                           &
                      XCGRID, YCGRID,                                 &
                      KGRPNT, IERR)
!  (This subroutine has not been used and tested yet)
!
!**********************************************************************

   USE schism_glbl, ONLY : skind,rkind, dkind

   USE TIMECOMM                                                        
   USE OCPCOMM4                                                        
   USE SWCOMM2                                                         
   USE SWCOMM3                                                         
   USE M_PARALL 
!   USE MOD_UNSTRUCT_GRID                                                       

   IMPLICIT NONE
!
!  2. PURPOSE
!
!     Update boundary conditions, update nonstationary input fields
!
!  4. Argument list
!
!     ARR      real  i  array holding values read from file (x-comp)
!     ARR2     real  i  array holding values read from file (y-comp)
!     TIMR2    real i/o time of last reading of input field
!     INTRV    real  i  time interval between input fields
!     TMENDR   real  i  end time of input field
!     IGR1     int   i  location in array COMPDA for interpolated input field data (x-comp)
!     IGR2     int   i  location in array COMPDA for interpolated input field data (y-comp)
!                       for a scalar field IGR2=0
!     JX1      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX2      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX3      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JY1      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY2      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY3      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     COSFC    real  i  cos of angle between input grid and computational grid
!     SINFC    real  i  sin of angle between input grid and computational grid
!     COMPDA   real i/o array holding values for computational grid points
!     XCGRID   real  i  x-coordinate of computational grid points
!     YCGRID   real  i  y-coordinate of computational grid points
!     KGRPNT   int   i  indirect addresses of computational grid points
!     NHDF     int   i  number of heading lines for a data file
!     NHDT     int   i  number of heading lines per time step
!     NHDC     int   i  number of heading lines before second component of vector field
!     IDLA     int   i  lay-out identifier for a data file
!     IDFM     int   i  format identifier for a data file
!     DFORM    char  i  format to read a data file
!     VFAC     real  i  multiplication factor applied to values from data file
!     IERR     int   o  error status: 0=no error, 9=end-of-file
!
!
   LOGICAL STPNOW                                                      
!
!  9. Structure
!
!     --------------------------------------------------------------
!     for all comp. grid points do
!         copy new values to old
!     --------------------------------------------------------------
!     repeat
!         if present time > time of last reading
!         then read new values from file
!              update time of last reading
!              interpolate values to computational grid
!         else exit from repeat
!     --------------------------------------------------------------
!     for all comp. grid points do
!         interpolate new values
!     --------------------------------------------------------------
!
! 10. SOURCE
!
!****************************************************************
!
   INTEGER    KGRPNT(MXC,MYC),                                     &
              IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
!
   REAL       COMPDA(MCGRD,MCMVAR),                                &
              XCGRID(MCGRD), YCGRID(MCGRD),                            &
              COSFC, SINFC
   REAL       ARR(*), ARR2(*)
!
!     local variables
!
   INTEGER    IENT, INDX, IX, IY
!     INDX       counter of comp. grid points
!     IX         index in x-dir of comput grid point
!     IY         index in y-dir of comput grid point
!
   REAL       SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
!JQI   REAL       TIMR1, XP, YP, UU, VV, VTOT, W1, W3,                &
   REAL       XP, YP, UU, VV, VTOT, W1, W3,                &
              SIZE1, SIZE2, SIZE3
   REAL(dkind)       TIMR1
!     TIMR1      time of one but last input field
!     XP         x-coord of one comput grid point
!     YP         y-coord of one comput grid point
!     UU         x-component of vector, or scalar value
!     VV         y-component of vector
!     VTOT       length of vector
!     W1         weighting coeff for interpolation in time
!     W3         weighting coeff for interpolation in time
!     DIRE       direction of interpolated vector
!     SIZE1      length of vector at time TIMR1
!     SIZE2      length of vector at time TIMCO
!     SIZE3      length of vector at time TIMR2
!
   INTEGER IG
   SAVE IENT
   DATA IENT /0/
   CALL STRACE (IENT, 'FLFILE')
!
   IERR = 0
!
   IF(JX1 > 1)THEN
     DO INDX = 2, MCGRD
       COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
     ENDDO
   ENDIF
   IF(IGR2 > 0 .AND. JY1 > 1)THEN
     DO INDX = 2, MCGRD
       COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
     ENDDO
   ENDIF
   TIMR1 = TIMCO - DTW
!
200 IF(TIMCO <= IFLTIM(IGR1)) GOTO 400
   TIMR1 = IFLTIM(IGR1)
   IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
   IF(IFLTIM(IGR1) > IFLEND(IGR1))THEN
     IFLTIM(IGR1) = 1.E10
     IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
     GOTO 400
   ENDIF
   IF(IFLNDS(IGR1) > 0)THEN                                         
     IF(INODE == MASTER)THEN
       CALL INAR2D( ARR, MXG(IGR1), MYG(IGR1),                       &
                    IFLNDF(IGR1),                                    &
		    IFLNDS(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),        &
		    IFLIDL(IGR1), IFLFAC(IGR1),                      &
		    IFLNHD(IGR1), IFLNHF(IGR1))
       IF(STPNOW()) RETURN                                           
     END IF
!JQI     CALL SWBROADC(IFLIDL(IGR1),1,SWINT)                              
     IF(IFLIDL(IGR1) < 0)THEN
!      end of file was encountered
       IFLTIM(IGR1) = 1.E10
       IF(IGR2 > 0) IFLTIM(IGR2) = IFLTIM(IGR1)
       GOTO 400
     ELSE
!JQI       CALL SWBROADC(ARR,MXG(IGR1)*MYG(IGR1),SWREAL)                   
     ENDIF
   ELSE                                                                
     IF(ITEST >= 20)THEN                                            
       CALL MSGERR (1,                                           &
                   'no read of input field because unit nr=0')                   
       WRITE (PRINTF, 208) IGR1
208    FORMAT (' field nr.', I2)
     ENDIF
   ENDIF                                                               
   IF(IGR2 > 0)THEN
     IFLTIM(IGR2) = IFLTIM(IGR1)
     IF(IFLNDS(IGR2) > 0)THEN                                       
       IF(INODE == MASTER)THEN                                      
         CALL INAR2D( ARR2, MXG(IGR2), MYG(IGR2),                   &
	              IFLNDF(IGR2),                                 &
		      IFLNDS(IGR2), IFLIFM(IGR2), IFLFRM(IGR2),     &
		      IFLIDL(IGR2), IFLFAC(IGR2), IFLNHD(IGR2), 0)
!        IF (STPNOW()) RETURN                                         
       END IF
!JQI       CALL SWBROADC(ARR2,MXG(IGR2)*MYG(IGR2),SWREAL)                  
     ENDIF                                                             
   ENDIF
!  Interpolation over the computational grid
   DO IG = 1, MCGRD
     IF(IGR2 == 0)THEN
!JQI       COMPDA(INDX,JX3) = ARR
     ELSE
!JQI       COMPDA(INDX,JX3) = ARR
!JQI       COMPDA(INDX,JY3) = ARR2
     ENDIF
   END DO
   
   GOTO 200
!
!   Interpolation in time
!
400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
     W1 = 1.-W3
     IF(ITEST >= 60) WRITE(PRTEST,402) IGR1,                          &
             TIMCO,IFLTIM(IGR1),W1,W3,JX1,JY1,JX2,JY2,JX3,JY3
402  FORMAT (' input field', I2, ' interp at ', 2F9.0, 2F8.3, 6I3)

   DO INDX = 2, MCGRD
     UU = W1 * COMPDA(INDX,JX2) + W3 * COMPDA(INDX,JX3)
     IF(IGR2 <= 0)THEN
       COMPDA(INDX,JX2) = UU
     ELSE
       VV = W1 * COMPDA(INDX,JY2) + W3 * COMPDA(INDX,JY3)
       VTOT = SQRT (UU*UU + VV*VV)
!
!      procedure to prevent loss of magnitude due to interpolation
!
       IF(VTOT > 0.)THEN
         SIZE1 = SQRT(COMPDA(INDX,JX2)**2 + COMPDA(INDX,JY2)**2)
         SIZE3 = SQRT(COMPDA(INDX,JX3)**2 + COMPDA(INDX,JY3)**2)
         SIZE2 = W1*SIZE1 + W3*SIZE3
!        SIZE2 is to be length of vector
         COMPDA(INDX,JX2) = SIZE2*UU/VTOT
         COMPDA(INDX,JY2) = SIZE2*VV/VTOT
       ELSE
         COMPDA(INDX,JX2) = UU
         COMPDA(INDX,JY2) = VV
       ENDIF
     ENDIF
   END DO	

   RETURN
!
   END SUBROUTINE FLFILE_OLD

!**********************************************************************
!
   SUBROUTINE FLWIND (IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR)
!
!**********************************************************************
!
!     Update nonstationary input fields, like Wind fields
!
!**********************************************************************

      USE TIMECOMM
      USE OCPCOMM4
      USE SWCOMM2
      USE SWCOMM3

      USE schism_glbl, ONLY : np_global,npa,np,skind,dkind,nws
      USE schism_glbl, ONLY : WINDX0=>WINDX, WINDY0=>WINDY
      USE schism_msgp, ONLY : myrank,parallel_abort,nproc

!# if defined (MULTIPROCESSOR)
!  USE MOD_PAR
!# endif
      USE VARS_WAVE, ONLY : INP_WI_NTIME,COMPDA

      IMPLICIT NONE

      LOGICAL STPNOW
      INTEGER :: IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
!
!     local variables
!
      INTEGER :: INDX,IENT
!
      CHARACTER*80 MSGSTR

      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'FLWIND')

      !REAL, ALLOCATABLE :: COMPDA_TMP1(:),COMPDA_TMP2(:)
!
      IERR = 0

      IF(JX1 > 1)THEN
       DO INDX = 1, npa
        COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
       ENDDO
      ENDIF

      IF(IGR2 > 0 .AND. JY1 > 1)THEN
       DO INDX = 1, npa
        COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
       ENDDO
      ENDIF

!      IF(nws.NE.2)THEN
! JL: quick Fix with PAHM 
      IF(nws.NE.2 .AND. nws.NE.-1 )THEN
         !MSGSTR = 'error in sbr FLWIND, nws!=2 not supported yet'
         MSGSTR = 'error in sbr FLWIND, nws!=2 or !=-1 not supported yet'
         CALL parallel_abort(MSGSTR)
      ENDIF

      IF(IGR1.NE.5 .or. IGR2.NE.6)THEN
         MSGSTR = 'error in sbr FLWIND, IGR1!=5 or IGR2!=6'
         CALL parallel_abort(MSGSTR)
      ENDIF
!
      COMPDA(1:npa,JX2) = WINDX0(1:npa)
      COMPDA(1:npa,JY2) = WINDY0(1:npa)

! No interpolation in time 
      COMPDA(1:npa,JX3) = COMPDA(1:npa,JX2)
      COMPDA(1:npa,JY3) = COMPDA(1:npa,JY2)

      !print*,'end  FLWIND'

      RETURN
!
  END SUBROUTINE FLWIND


!********************************************************************
!
      SUBROUTINE ERRCHK                                                   
!
!****************************************************************
!
      USE OCPCOMM4                                                        
      USE SWCOMM1                                                         
      USE SWCOMM2                                                         
      USE SWCOMM3                                                         
      USE SWCOMM4                                                         
      USE M_GENARR  
     
      USE schism_glbl, ONLY : ics
 
      IMPLICIT NONE                                                      
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Check all possible combinations of physical processes if
!     they are being activated and change value of settings if
!     necessary
!
!  3. Method
!
!     0      MESSAGE
!     1      WARNING
!     2      ERROR REPAIRABLE
!     3      SEVERE ERROR (calculation continues, however problems
!                          may arise)
!     4      TERMINATION ERROR (calculation is terminated )
!
!  4. Argument variables (updated 30.72)
!
!  6. Local variables
!
!     CHARS :     array to pass character info to MSGERR
!     MSGSTR:     string to pass message to call MSGERR
!
      CHARACTER*20 NUMSTR, CHARS(3)
      CHARACTER*80 MSGSTR
!
!  8. Subroutines used
!
!     MSGERR : Handles error messeges according to severity
!     NUMSTR : Converts integer/real to string
!     TXPBLA : Removes leading and trailing blanks in string
!
!  9. Subroutines calling
!
!     SWANCOM
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     JL : 03/2017 ERRCHK upgraded to version 41.10 
!
! 12. Structure
!
!     ------------------------------------------------------------
!     End of the subroutine ERRCHK
!     ------------------------------------------------------------
!
! 13. Source text
!
      LOGICAL  STPNOW                                                     
      LOGICAL  EQREAL  
      INTEGER :: GAMMA,II                                                   
      INTEGER, SAVE :: IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'ERRCHK')

!     --- choose scheme for stationary or nonstationary computation

!      PRINT*,'NSTATC = ',NSTATC
      IF (NSTATC.EQ.1) THEN                                               
        PROPSC = PROPSN                                                   
      ELSE                                                                
        PROPSC = PROPSS                                                   
      ENDIF                                                               
!
!     -----------------------------------------------------------------
!
!     *** WARNINGS AND ERROR MESSAGES ***
!
!     -----------------------------------------------------------------

!     *** Check formulation for whitecapping ***
!
      IF ( IWCAP.GT.8 ) THEN
         WRITE (MSGSTR, '(A,I2,A)') &
     &               'Unknown method for whitecapping (IWCAP=',IWCAP,')'
         CALL MSGERR( 1, TRIM(MSGSTR) )
         CALL MSGERR( 1, &
     &     'Whitecapping according to Komen et al. (1984) will be used')
         IWCAP     = 1
         PWCAP(1)  = 2.36E-5
         PWCAP(2)  = 3.02E-3
         PWCAP(9)  = 2.
         PWCAP(10) = 1.                                                  !41.41
         PWCAP(11) = 1.
      ENDIF
!
!     *** WAM cycle 3 physics ***
!
      IF((IWIND .EQ. 3 .OR. IWIND .EQ. 5) .AND. IWCAP .NE. 1             &
         .AND. IWCAP .NE. 7)THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Komen et al. (1984) for wind option G3/YAN')       
      ENDIF

      IF(IWCAP .EQ. 7 .AND. IWIND .NE. 5)THEN                            
        CALL MSGERR(1,'Activate wind option Yan (1987) in case of')       
        CALL MSGERR(1,'Alves & Banner (2003) white-capping method')       
      END IF                                                              
!
!     *** WAM cycle 4 physics ***
!
      IF(IWIND .EQ. 4 .AND. IWCAP .NE. 2)THEN
        CALL MSGERR(1,'Activate whitecapping mechanism according to')
        CALL MSGERR(1,'Janssen (1991) for wind option JANS        ')
      ENDIF
!
!     *** check numerical scheme in presence of a current ***
!
      IF(ICUR .EQ. 1)THEN
        IF(PNUMS(6) .EQ. 0.)THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in theta space ')
          CALL MSGERR(1,'-> set CDD = 1.')
          WRITE(PRINTF,*)
        ENDIF
        IF(PNUMS(7) .EQ. 0.)THEN
          CALL MSGERR(1,'In presence of a current it is recommended to')
          CALL MSGERR(1,'use an implicit upwind scheme in sigma space ')
          CALL MSGERR(1,'-> set CSS = 1.')
          WRITE(PRINTF,*)
        ENDIF
      END IF
!
!     check combination of REPeating option and grid type and dimension   
!
!JQI      IF(KREPTX > 0)THEN                                               
!       --- "GE" changed to "EQ" since OPTG.EQ.3 FOR CURVILINEAR          
!JQI        IF(OPTG == 3)                                             &
!JQI	  CALL MSGERR (3, 'Curvilinear grid cannot be REPeating')
!JQI        IF(PROPSC == 1 .AND. MXC < 1)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=1 for REPeating option')
!JQI        IF(PROPSC == 2 .AND. MXC < 2)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=2 for REPeating option')           
!JQI        IF(PROPSC == 3 .AND. MXC < 3)                             &
!JQI	  CALL MSGERR (3, 'MXC must be >=3 for REPeating option')           
!JQI      ENDIF
!
#     if defined (SPHERICAL)
      IF(OPTG /= 1 .AND. KSPHER /= 0)THEN                               
         CALL MSGERR (3, 'Spherical coordinates must be given in')        
         CALL MSGERR (3, 'uniform, recti-linear computational grid')      
      END IF                                                              
#     endif

      IF(ICS ==2 .AND. KSPHER /= 1)THEN
         CALL MSGERR (3, 'Spherical coordinates must be given ')
         CALL MSGERR (3, 'because ics==2 in param.in')
      END IF

      IF(ICS ==1 .AND. KSPHER /= 0)THEN
         CALL MSGERR (3, 'Carteseian coordinates must be given ')
         CALL MSGERR (3, 'because ics==1 in param.in')
      END IF
!
      IF(PROPSC == 2 .AND. NSTATC > 0)THEN                             
        CALL MSGERR (3, 'SORDUP scheme only in stationary run')           
      ENDIF
      IF(PROPSC == 3 .AND. NSTATC == 0)THEN                             
        CALL MSGERR (3, 'S&L scheme not in stationary run')               
      ENDIF
!
!     --- A warning about curvilinear and S&L scheme                      
      IF((PROPSC == 3).AND.(OPTG == 3))THEN                             
         CALL MSGERR(1,'the S&L scheme (higher order nonstationary')      
         CALL MSGERR(1,'IS NOT fully implemented for curvilinear')        
         CALL MSGERR(1,'coordinates. This may or may not be noticeable')  
         CALL MSGERR(1,'in simulations. DX and DY are approximated')      
         CALL MSGERR(1,'with DX and DY of two nearest cells.')            
         CALL MSGERR(1,'Note that SORDUP (higher order stationary)')      
         CALL MSGERR(1,'IS fully implemented for curvilinear coord.')     
         CALL MSGERR(1,'and differences are usually negligible.')         
      ENDIF                                                               
!
!     Here the various problems with quadruplets are checked
!     The combination of quadruplets and sectors is an error
!     in the calculation of quadruplets when the SECTOR option is
!     used in the CGRID command. This error should be corrected
!     in the future
!
      IF (ICUR.NE.0 .AND. (IQUAD.EQ.1 .OR. IQUAD.EQ.2)) THEN              !41.55
         CALL MSGERR(1,'Quadruplets will be updated per iteration')       !41.55
         CALL MSGERR(1,'instead of per sweep. This will increase')        !41.55
         CALL MSGERR(1,'amount of internal memory with a factor 2')       !41.55
         IQUAD = 3                                                        !41.55
      ENDIF

      IF(IWIND == 3 .OR. IWIND == 4)THEN                                
        IF(IQUAD == 0)THEN
          CALL MSGERR(2,'Quadruplets should be activated when SWAN  ')    
          CALL MSGERR(2,'is running in a third generation mode and  ')    
          CALL MSGERR(2,'wind is present                            ')    
        ENDIF
      ENDIF
!
      IF(IQUAD >= 1)THEN                                              
!
       IF(.NOT. FULCIR)THEN                                             
        IF((SPDIR2-SPDIR1) < (PI_W/12.))THEN                           
          CALL MSGERR(2,'A combination of using quadruplets with a'    )  
          CALL MSGERR(2,'sector of less than 30 degrees should be'     )  
          CALL MSGERR(2,'avoided at all times, it is likely to produce')  
          CALL MSGERR(2,'unreliable results and unexpected errors.'    )  
          CALL MSGERR(2,'Refer to the manual (CGRID) for details'      )  
        ELSE                                                              
          CALL MSGERR(1,'It is not recommended to use quadruplets'     )  
          CALL MSGERR(1,'in combination with calculations on a sector.')  
          CALL MSGERR(1,'Refer to the manual (CGRID) for details'      )  
        END IF                                                            
       END IF                                                             
!
       IF(IWIND == 0)THEN                                               
         CALL MSGERR(2,'It is not recommended to use quadruplets'     )   
         CALL MSGERR(2,'in combination with zero wind conditions.'    )   
       END IF                                                             
!
       IF(MSC == 3)THEN
         CALL MSGERR(4,'Do not activate quadruplets for boundary ')
         CALL MSGERR(4,'option BIN -> use other option           ')
         WRITE(PRINTF,*)
       END IF
!
      END IF                                                              

! JL Old Rule before version 41.10
#if 0
!
!     Check whether limiter should be de-activated                        
!
      IF(PNUMS(20) < 100.)THEN                                       
         IF(IGEN == 3)THEN
!        --- third generation mode
            IF(IQUAD == 0 .AND. ITRIAD == 0)THEN
               CALL MSGERR(1,'Limiter is de-activated in case of')
               CALL MSGERR(1,'no wave-wave interactions')
               PNUMS(20) = 1.E+20
            END IF
         ELSE IF(IGEN < 3)THEN
!        --- first or second generation mode
            IF(ITRIAD == 0)THEN
               CALL MSGERR(1,'Limiter is de-activated in case of')
               CALL MSGERR(1,'first or second generation mode')
               PNUMS(20) = 1.E+20
            END IF
         END IF
      END IF
# endif
!JL New Rule :

!     Check whether limiter should be de-activated
!
      IF (IQUAD.EQ.0 .AND. PNUMS(20).LT.100.) THEN                        !40.41
         CALL MSGERR(1,&
     &              'Limiter is de-activated in case of no quadruplets')  !40.41
         PNUMS(20) = 1.E+20
      END IF
!
!     Check resolution in frequency-space when DIA is used                
!
      IF(IQUAD > 0 .AND. IQUAD <= 3 .OR. IQUAD == 8)THEN               
         GAMMA = EXP(MyLOG(SHIG/SLOW)/REAL(MSC-1))                         
         IF(ABS(GAMMA-1.1) > 0.055)THEN                                
            CALL MSGERR(1,                                        &
	        'relative frequency resolution (df/f) deviates more')    
            CALL MSGERR(1,                                        &
	        'than 5% from 10%-resolution. This may be problematic')  
            CALL MSGERR(1,                                        &
	        'when quadruplets are approximated by means of DIA.')    
         END IF                                                           
      END IF                                                              
!
!     When using triads MSC must be less than 200!                        
!
      IF((ITRIAD > 0).AND.(MSC > 200))THEN                            
         CALL MSGERR(4,'When triads are active the number of     ')       
         CALL MSGERR(4,'directions chosen in the CGRID command   ')       
         CALL MSGERR(4,'must be less than 200                    ')       
      END IF                                                              
!
!JL Old Rule before version 41.10 
#if 0
!     When CSM is applied only DIA should be used                         
!
      IF(IWCAP == 6 .AND. IQUAD > 3 .AND. IQUAD /= 8)THEN              
         CALL MSGERR(3,                                            &
	     'In case of cumulative steepness method, only DIA '// &
	     'technique is supported')                                    
      END IF                                                              
# endif
!
!     When Alves and Banner and XNL are applied change the parameters     
!
      IF(IWCAP == 7 .AND. (IQUAD == 51 .OR. IQUAD == 52 .OR.        &
         IQUAD == 53))THEN                          
         IF(EQREAL(PWCAP( 1), 5.0E-5)) PWCAP( 1) = 5.0E-5                
         IF(EQREAL(PWCAP(12),1.75E-3)) PWCAP(12) = 1.95E-3               
         IF(EQREAL(PWCAP(10),    4.0)) PWCAP(10) = 4.                    
         IF(EQREAL(PWCAP( 9),    0.0)) PWCAP( 9) = 0.                    
         IF(EQREAL(PWCAP(11),    0.0)) PWCAP(11) = 0.                    
      END IF                                                              
!OMP!
!OMP      IF ( IQUAD.EQ.51 .OR. IQUAD.EQ.52 .OR. IQUAD.EQ.53 ) THEN           
!OMP         CALL MSGERR(4,'XNL not supported within OpenMP environment')     
!OMP      END IF                                                              
!
      IF(ITEST >= 120)THEN
        WRITE(PRINTF,3000) IWIND ,IQUAD, ICUR, IWCAP, MSC
3000    FORMAT(' ERRCHK : IWIND QUAD CUR WCAP MSC   : ',5I4)
        IF(IWIND > 0)THEN                                            
          DO II = 1, MWIND
            WRITE(PRINTF,30) II,PWIND(II)
 30         FORMAT(' PWIND(',I2,') = ',E11.4)
          ENDDO
        ENDIF
      END IF
!
      RETURN
      END SUBROUTINE ERRCHK

!***********************************************************************
!                                                                      *
  SUBROUTINE SWRBC !( COMPDA )                              !40.31 30.90
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                   !     40.41
      USE SWCOMM1                                                    !     40.41
      USE SWCOMM2                                                    !     40.41
      USE SWCOMM3                                                    !     40.41
      USE SWCOMM4                                                    !     40.41
      USE M_GENARR
      USE SwanGriddata                                               !     40.80
!
      USE schism_glbl, ONLY : np_global,npa,np,skind,dkind
      USE schism_glbl, ONLY : WINDX0=>WINDX, WINDY0=>WINDY
      USE VARS_WAVE, ONLY : COMPDA,KN 
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.60: Nico Booij
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     32.02: Roeland Ris & Cor van der Schelde (1D-version)
!     33.08: W. Erick Rogers
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     00.00, Nov. 86: existing version
!     00.01, Apr. 87: parameters added in call of subroutine SVALQI,
!                     some variable names changed (not common anymore)
!     30.60, Aug. 97: test on value of KGRPNT, skip part of code if
!                     KGRPNT(IX,IY)=1
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     32.01, Jan. 98: Initialise setup and saved depth for 1D-version
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     30.70, Mar. 98: proper water level stored in array COMPDA(*,JWLV2)
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Nov. 98: Now also interpolates curvilinear current input-fields
!                     (IGTYPE(2)=2)
!     40.00, Feb. 99: IDYNWI etc. replaced by IFLDYN(*)
!     33.08, July 98: minor changes related to the S&L scheme
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Dec. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!
!  2. Purpose
!
!     The depths and currents at a line in the computational grid are
!     determined and written to file with reference number NREF
!
!  3. Method
!
!     The depths and currents are computed by bilinear interpolation
!     and usually written to file INSTR.
!
!  4. Argument variables
!
!     COMPDA
!
!JL      REAL    COMPDA(MCGRD,MCMVAR)                                        !30.21
!
!  6. Local variables
!
!  8. Subroutines used
!
!     SVALQI (SWAN/SER)
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     ----------------------------------------------------------------
!     For every line IX of the computational grid do
!         For every point IY of this line do
!             Compute bottom grid coordinates as number of meshes
!             Call SVALQI to interpolate depth and current for the point
!             If current is on (ICUR = 1), then
!                 Compute current components relative to comp. grid
!             Else
!                 Current components are zero
!     ----------------------------------------------------------------
!
! 13. Source text
!
      REAL(rkind), parameter :: ZERO = 0._rkind

      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SWRBC')
!
!JL Skip that
# if 0
      IF (ITEST .GE. 100 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,*) '  ', ICUR, IGTYPE(2)
        WRITE(PRINTF,*) '******** In subroutine SWRBC *******'            !30.21
        IF (ICUR .EQ. 1 ) THEN
          WRITE(PRINTF,51)
        ELSE
          WRITE(PRINTF,50)
        ENDIF
      ENDIF
 50   FORMAT('P.index',5X,'coord.',14X,'depth')
 51   FORMAT('P.index',5X,'coord.',14X,'depth',13X,'UX',13X,'UY')
# endif
!
!     *** The arrays start to be filled in the second value ***
!     *** because in COMPDA(1,"variable"), is the default   ***
!     *** value for land points    version 30.21            ***
!
! JL: for debug purpose
#if 0
IF(myrank==0) THEN
PRINT *,'JDP1',JDP1
PRINT *,'JDP3',JDP3
PRINT *,'JWLV1',JWLV1
PRINT *,'JWLV3',JWLV3
PRINT *,'JVX1',JVX1
PRINT *,'JVX3',JVX3
PRINT *,'JFRC2',JFRC2
PRINT *,'JFRC3',JFRC3
PRINT *,'JWX2',JWX2
PRINT *,'JASTD2',JASTD2
PRINT *,'JASTD3',JASTD3
PRINT *,'JMUDL1',JMUDL1
PRINT *,'JMUDL2',JMUDL2
PRINT *,'JMUDL3',JMUDL3
PRINT *,'JNPLA2',JNPLA2
PRINT *,'JNPLA3',JNPLA3
PRINT *,'JTURB2',JTURB2
PRINT *,'JTURB3',JTURB3
PRINT *,'JBOTLV',JBOTLV
ENDIF
#endif

!     ***  Default values for land point  ***
      COMPDA(1,JDP2) = -1.                                                !40.00
      IF (JDP1.GT.1) COMPDA(1,JDP1) = -1.                                 !40.00
      IF (JDP3.GT.1) COMPDA(1,JDP3) = -1.                                 !40.00
      IF (VARWLV) THEN
        COMPDA(1,JWLV2) = ZERO                                            !40.00
!       next two lines only for nonstat water level
        IF (JWLV1.GT.1) COMPDA(1,JWLV1) = ZERO                            !40.00
        IF (JWLV3.GT.1) COMPDA(1,JWLV3) = ZERO                            !40.00
      ENDIF
      IF (ICUR.GT.0) THEN
        COMPDA(1,JVX2) = ZERO                                             !40.00
        COMPDA(1,JVY2) = ZERO                                             !40.00
        IF (JVX1.GT.1) COMPDA(1,JVX1) = ZERO                              !40.00
        IF (JVY1.GT.1) COMPDA(1,JVY1) = ZERO                              !40.00
        IF (JVX3.GT.1) COMPDA(1,JVX3) = ZERO                              !40.00
        IF (JVY3.GT.1) COMPDA(1,JVY3) = ZERO                              !40.00
      ENDIF
      IF (VARFR) THEN
        COMPDA(1,JFRC2) = ZERO                                            !40.00
        COMPDA(1,JFRC3) = ZERO                                            !40.00
      ELSEIF (JFRC2.GT.1) THEN                                            !41.51
        COMPDA(1,JFRC2) = ZERO                                            !41.51
      ENDIF
      IF (VARWI) THEN
        COMPDA(1,JWX2) = ZERO                                             !40.00
        COMPDA(1,JWY2) = ZERO                                             !40.00
        IF (JWX3.GT.1) COMPDA(1,JWX3) = ZERO                              !40.00
        IF (JWY3.GT.1) COMPDA(1,JWY3) = ZERO                              !40.00
      ENDIF
      IF (VARAST) THEN
        COMPDA(1,JASTD2) = ZERO                                           !40.03
        COMPDA(1,JASTD3) = ZERO                                           !40.03
      ENDIF
      IF (VARMUD) THEN                                                    !40.59
        COMPDA(1,JMUDL1) = ZERO                                           !40.59
        COMPDA(1,JMUDL2) = ZERO                                           !40.59
        COMPDA(1,JMUDL3) = ZERO                                           !40.59
      ENDIF
      IF (VARNPL) THEN                                                    !40.55
        COMPDA(1,JNPLA2) = ZERO                                           !40.55
        COMPDA(1,JNPLA3) = ZERO                                           !40.55
      ENDIF
      IF (VARTUR) THEN                                                    !40.35
        COMPDA(1,JTURB2) = ZERO                                           !40.35
        COMPDA(1,JTURB3) = ZERO                                           !40.35
      ENDIF
      COMPDA(1,JBOTLV) = ZERO                                             !40.65
!
!     --- unstructured grid                                               40.80
!
!     JL Notes: All quantities below are updated later in SNEXTI
!               Used for Init. purpose
!
      DO INDX = 1, nverts
!
         XP = xcugrd(INDX)                                                !40.80
         YP = ycugrd(INDX)                                                !40.80
!
!        ***** compute depth and water level *****
!
         DEP = DEPTH(INDX)
         COMPDA(INDX,JBOTLV) = DEP                                        !40.80
!
!        add constant water level
         DEPW = MAX(ZERO, DEP + WLEV)       !DEP + WLEV                   !40.80
         COMPDA(INDX,JDP2) = DEPW

!        ***In this step the water level at T+DT is copied to the  ***
!        ***water level at T (Only for first time computation)     ***
         IF (JDP1.GT.1) COMPDA(INDX,JDP1) = DEPW                          !40.80
         IF (JDP3.GT.1) COMPDA(INDX,JDP3) = DEPW                          !40.80
!
!        ***** compute current velocity *****
!
         IF (ICUR.EQ.1 .AND. IGTYPE(2) .GE. 1) THEN                       !40.80

            IF (DEPW.GT.0.) THEN
              UU = UXB(INDX,1)
              VV = UYB(INDX,1)
              VTOT = SQRT (UU*UU + VV*VV)
              CGMAX = PNUMS(18)*SQRT(GRAV_W*DEPW)
              IF (VTOT .GT. CGMAX) THEN
                CGFACT = CGMAX / VTOT
                UU = UU * CGFACT
                VV = VV * CGFACT
                IF (ERRPTS.GT.0) THEN
                   WRITE (ERRPTS, 212) INDX, 1
 212               FORMAT (I4, 1X, I2)
                ENDIF
              ENDIF
              COMPDA(INDX,JVX2) =  UU*COSVC + VV*SINVC
              COMPDA(INDX,JVY2) = -UU*SINVC + VV*COSVC

            ELSE
              COMPDA(INDX,JVX2) =  ZERO
              COMPDA(INDX,JVY2) =  ZERO
            ENDIF

            IF (JVX1.GT.1) COMPDA(INDX,JVX1) = COMPDA(INDX,JVX2)          !40.80
            IF (JVY1.GT.1) COMPDA(INDX,JVY1) = COMPDA(INDX,JVY2)          !40.80
            IF (JVX3.GT.1) COMPDA(INDX,JVX3) = COMPDA(INDX,JVX2)          !40.80
            IF (JVY3.GT.1) COMPDA(INDX,JVY3) = COMPDA(INDX,JVY2)          !40.80
         ENDIF
!
!         ***** compute variable friction coefficient *****
!
         IF (VARFR) THEN
!JQI            IF ( IGTYPE(4).EQ.3 ) THEN
!JQI               FRI = FRIC(INDX)
!JQI            ELSE
!JQI               FRI = SVALQI (XP, YP, 4, FRIC, 1, 0, 0)                    !40.80
!JQI            ENDIF
!JQI            COMPDA(INDX,JFRC2) = FRI                                      !40.80
!JQI            COMPDA(INDX,JFRC3) = FRI                                      !40.80
!            WRITE(6,*) 'IBOT,JFRC2,JFRC3',IBOT,JFRC2,JFRC3
            IF (IBOT.EQ.3) THEN
                !JL: MADSEN Scheme for Bottom Dissipation, use KN, Roughness Length
                FRI = KN(INDX)
                COMPDA(INDX,JFRC2) = FRI
                COMPDA(INDX,JFRC3) = FRI
!                WRITE (PRINTF, 214) INDX, FRI                         
! 214        FORMAT (' INDX, FRI: ', I7,2X,E11.4)                    

            ENDIF  
         ENDIF
!
!        ***** compute variable wind velocity *****
!
         IF (VARWI) THEN

          COMPDA(INDX,JWX2) = WINDX0(INDX)
          COMPDA(INDX,JWY2) = WINDY0(INDX)
          IF (JWX3.GT.1) COMPDA(INDX,JWX3) = COMPDA(INDX,JWX2)          !40.80
          IF (JWY3.GT.1) COMPDA(INDX,JWY3) = COMPDA(INDX,JWY2)          !40.80

         ENDIF
!
!     ***** compute variable air-sea temperature difference *****         40.80
!
         IF (VARAST) THEN
!JQI            IF ( IGTYPE(10).EQ.3 ) THEN
!JQI               ASTD = ASTDF(INDX)
!JQI            ELSE
!JQI               ASTD = SVALQI (XP, YP, 10, ASTDF, 1, 0, 0)                 !40.80
!JQI            ENDIF
!JQI            COMPDA(INDX,JASTD2) = ASTD                                    !40.80
!JQI            COMPDA(INDX,JASTD3) = ASTD                                    !40.80
         ENDIF
!
!     ***** compute number of plants per square meter *****               40.80
!
          IF (VARNPL) THEN
!JL             IF ( IGTYPE(11).EQ.3 ) THEN
                XNPL = NPLAF(INDX)
!JL             ELSE
!JL                XNPL = SVALQI (XP, YP, 11, NPLAF, 1 ,0, 0)               ! 40.80
!JL             ENDIF
             COMPDA(INDX,JNPLA2) = XNPL                                   !40.80
             COMPDA(INDX,JNPLA3) = XNPL                                   !40.80
          ENDIF
!
!     ***** compute turbulent viscosity *****                             40.35
!
          IF (VARTUR) THEN
             IF (PTURBV(2).LT.0.) THEN
!JL                IF ( IGTYPE(12).EQ.3 ) THEN
                   XTUR = TURBF(INDX)
!JL                ELSE
!JL                   XTUR = SVALQI (XP, YP, 12, TURBF, 0 ,0, 0)             !40.35
!JL                ENDIF
             ELSE
                UU = COMPDA(INDX,JVX2)
                VV = COMPDA(INDX,JVY2)
                XTUR = PTURBV(2) * SQRT(UU*UU+VV*VV) * COMPDA(INDX,JDP2)  !40.35
             ENDIF
             COMPDA(INDX,JTURB2) = XTUR                                   !40.35
             COMPDA(INDX,JTURB3) = XTUR                                   !40.35
          ENDIF
!
!     ***** compute fluid mud layer *****                                 40.80
!
          IF (VARMUD) THEN
!JL             IF ( IGTYPE(13).EQ.3 ) THEN
                XMUD = MUDLF(INDX)
!JL             ELSE
!JL                XMUD = SVALQI (XP, YP, 13, MUDLF, 1 ,0, 0)                !40.80
!JL             ENDIF
             COMPDA(INDX,JMUDL1) = XMUD                                   !40.80
             COMPDA(INDX,JMUDL2) = XMUD                                   !40.80
             COMPDA(INDX,JMUDL3) = XMUD                                   !40.80
          ENDIF
!
      ENDDO                                                               !40.80
!
!     *** initialise setup and saved depth ***                            32.02
!
      IF (LSETUP.GT.0) THEN                                               !32.02
        DO INDX = 1, MCGRD                                                !32.02
          COMPDA(INDX,JSETUP) = ZERO                                      !32.02
          COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)                         !32.02
        ENDDO                                                             !32.02
      ENDIF                                                               !32.02
!
!     --- initialize HSIBC                                                40.41
      COMPDA(:,JHSIBC) = ZERO                                             !40.41
!
!     --- initialize ZELEN and USTAR                                      40.41
      COMPDA(:,JZEL  ) = 2.E-33                                           !40.41
      COMPDA(:,JUSTAR) = 1.E-15                                           !40.41

! JL  --- initialize CHARN
      COMPDA(:,JCHARN)= ZERO
!
!     --- initialize UBOT and TMBOT                                       40.94
      COMPDA(:,JUBOT) = 0.                                                !40.94
      IF (JPBOT.GT.1) COMPDA(:,JPBOT) = ZERO                              !40.94
!
!     --- initialize GAMBR and RESPL                                      41.38
      IF (ISURF.EQ.6) THEN
         COMPDA(:,JGAMMA) = PSURF( 2)
         COMPDA(:,JRESPL) = PSURF(13)*COMPDA(:,JDP2)
      ENDIF
!
      RETURN
! * end of subroutine SWRBC *
      END


!*********************************************************************
!                                                                    *
      SUBROUTINE SNEXTI (                                       &  !40.31
     &                   SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT,&  !40.31
     &                   XYTST , DEPTH , WLEVL , FRIC  , UXB   ,&  !40.31
     &                   UYB   , NPLAF , TURBF , MUDLF , WXI   ,&  !40.59 40.35 40.55 40.31
     &                   WYI   )
!                                                                    *
!*********************************************************************
!
      USE TIMECOMM                                                       !40.41
      USE OCPCOMM4                                                       !40.41
      USE SWCOMM1                                                        !40.41
      USE SWCOMM2                                                        !40.41
      USE SWCOMM3                                                        !40.41
      USE SWCOMM4                                                        !40.41
      USE M_BNDSPEC                                                      !40.31
      USE SwanGriddata                                                   !40.80

      USE schism_glbl, ONLY : skind,rkind,dkind,iplg
      USE schism_glbl, ONLY : dp,eta2,uu2,vv2,nvrt,errmsg
      USE schism_msgp, ONLY : myrank,parallel_abort

      USE VARS_WAVE, ONLY : IOBCN_W,I_OBC_N_W,AC1,AC2,NESTING_TYPE_WAVE
      USE VARS_WAVE, ONLY : COMPDA,WBAC,NESTING,KN
      USE VARS_WAVE, ONLY : curx_wwm,cury_wwm
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.60: Nico Booij
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.74: IJsbrand Haagsma (Include version)
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     40.00, 40.13: Nico Booij
!     34.01: Jeroen Adema
!     40.14: Annette Kieftenburg
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     30.60, Jun. 97: condition for ATAN2 corrected
!     30.70, Sept 97: reduction of current only if depth is positive
!     30.72, Sept 97: INTEGER*4 replaced by INTEGER
!     30.72, Oct. 97: changed floating point comparison to avoid equality
!                     comparisons
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.70, Jan. 98: VNAM6 (nonstat current) corrected
!     30.70, Feb. 98: argument AUXW4 added in call of WAM nesting
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.70, Feb. 98: wind is no longer set to 0, id depth is negative
!     40.00, Nov. 97: complete revision of boundary value update
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description of several variables
!     34.01, Feb. 99: Introducing STPNOW
!     40.13, Mar. 01: Loop over INDX replaced by loop over IX, IY
!     40.14, Jun. 01: Waterlevel updated in case set-up is on
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Updates boundary conditions and input fields
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
! i   BGRIDP: data for interpolating to computational grid points
! i   KGRPNT: computational grid point addresses
! i   XYTST : test points
!
!JL   INTEGER  BGRIDP(*), XYTST(*), KGRPNT(MXC,MYC)
      INTEGER  XYTST(*), KGRPNT(MXC,MYC)
!
! i   AC1   : action density spectra on old time level
! i   AC2   : action density spectra on new time level
! i   BSPECS: boundary spectra
! i   COMPDA: values on computational grid
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
!
!      REAL     AC1(MDC,MSC,MCGRD)                                         !40.00 
!      REAL     AC2(MDC,MSC,MCGRD)                                         !30.21
!      REAL     BSPECS(MDC,MSC,NBSPEC,2)                                   !40.00
!      REAL     COMPDA(MCGRD,MCMVAR)                                       !30.72
      REAL(rkind)     SPCDIR(MDC,6)                                              !40.00
      REAL(rkind)     SPCSIG(MSC)                                                !30.01
      REAL(rkind)     XCGRID(MXC,MYC), YCGRID(MXC,MYC)                           !30.21
      REAL(rkind)     DEPTH(*), WLEVL(*), FRIC(*), UXB(*), UYB(*),&              !40.31
     &         NPLAF(*), TURBF(*),                         &              !40.35 40.55
     &         MUDLF(*),                                   &              !40.59
     &         WXI(*), WYI(*)                                             !40.31
!
!     TIMCO ..... Time (date) of computation
!     TFINC ..... Final time (date) of computation
!     DT    ..... Increment time for computation
!     TIMCU ..... Date to read the next current file.
!     TIMFR .....        "              friction
!     TIMWI .....        "              wind
!     TIMWL .....        "              water level
!     WEI??#..... Weights for linear interpolation for
!                 (??=) CUR, FRC, WIN, WLV,
!                 (#=) 2 for field at Ti and 1 for field at Ti+1
!     VARWE?..... Variation of the WEI??? in each DT for C, F, W , L
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     SWBROADC
!
      LOGICAL STPNOW                                                      !34.01
!
!  9. Subroutines calling
!
!     ---
!
! 10. Error messages
!
!     ---
!
! 11. Remarks
!
!     ---
!
! 12. Structure
!
!     -------------------------------------------------------------
!     ------------------------------------------------------------
!
! 13. Source text
!
      INTEGER     IERR
      REAL(rkind)        TSTVAL(10)                                              !40.00
      TYPE(BSPCDAT), POINTER :: CURBFL                                    !40.31
      LOGICAL LPB                                                         !40.80

!JL SWAN in FVCOM Local Vars
      INTEGER :: IOB,IOB_NODE,IBGRID,INDXGR,IS,ID,IG
      INTEGER :: K1,K2
      REAL(rkind)    :: W1,W2,ETOT,SIG2,HS,XP,YP,DEP,WLVL,DEPW,UU,VV
      REAL(rkind)    :: VTOT,CGMAX,CGFACT
      CHARACTER*80 :: MSGSTR
      REAL(skind), ALLOCATABLE :: AC2_TTMP(:,:,:)
!JL End SWAN in FVCOM Local Vars
      CHARACTER(LEN=80)  :: SIMTIME

      SAVE  IENT
      DATA  IENT/0/
      IF (LTRACE) CALL STRACE(IENT,'SNEXTI')

!     **   All action densities are shifted from array T+DT
!     **   to the array at time T
!
      IF (NSTATC.EQ.1) THEN                                               !30.70
        IF (ITERMX.GT.1                  &
     &       .OR. PROPSC.EQ.3            &                                !33.09
     &                                 ) THEN                             !30.70

         IF(.NOT.ALLOCATED(AC1))THEN
           MSGSTR = 'Allocation problem: array AC1'
           WRITE(errmsg,*) MSGSTR
           CALL parallel_abort(errmsg)

         ELSE IF(SIZE(AC1).EQ.0)THEN
           MSGSTR = 'Allocation problem: AC1 Size is Wrong'
           WRITE(errmsg,*) MSGSTR
           CALL parallel_abort(errmsg)
         END IF

         AC1(1:MDC,1:MSC,1:MCGRD) = AC2(1:MDC,1:MSC,1:MCGRD)

!          DO 50 IXY = 1, MCGRD                                            !30.21
!            IF(MOD(IXY,10).EQ.0) PRINT *,'IXY',IXY,'/',MCGRD
!            DO 55 ISS = 1, MSC
!              DO 60 IDD = 1, MDC
!                AC1(IDD,ISS,IXY) = AC2(IDD,ISS,IXY)                       !30.21
! 60           CONTINUE
! 55         CONTINUE
! 50       CONTINUE
!
        ENDIF                                                             !40.00
      ENDIF

! ------------------------------------------------
! --- update boundary conditions, 
! ------------------------------------------------

! --- Steady and Uniform wave spectral density ...

  IF(.NOT.NESTING) THEN

   DO IOB = 1, IOBCN_W
    IOB_NODE = I_OBC_N_W(IOB)
    !if(myrank==0) PRINT *,'IOB_NODE=',IOB_NODE

    !if(myrank==0) PRINT *,'SNEXTI updating AC2, node =',IOB_NODE
    !print*,'SNEXTI myrank=',myrank,'updating AC2, node =',IOB,IOB_NODE

    CALL SSHAPE(AC2(1,1,IOB_NODE),SPCSIG,SPCDIR,FSHAPE,DSHAPE)

    !if(myrank==0) PRINT *,'SNEXTI MAXVal AC2 =',maxval(AC2(:,:,IOB_NODE))

   END DO

  ELSE ! NESTING

   IF((TRIM(NESTING_TYPE_WAVE).EQ.'PARAM').OR. &
      (TRIM(NESTING_TYPE_WAVE).EQ.'WAMPARAM'))THEN

    CALL SET_WAVE_BOUNDARY_CONDITION

    DO IOB = 1, IOBCN_W
      IOB_NODE = I_OBC_N_W(IOB)
      AC2(:,:,IOB_NODE) = WBAC(:,:,IOB)
    END DO

    if(myrank==0) &
     WRITE(PRINTF,*)'SET_WAVE_BOUNDARY_CONDITION done'

   ELSE IF(TRIM(NESTING_TYPE_WAVE).EQ.'SPEC')THEN

!    CALL parallel_abort('SNEXTI SPEC ToDO')

    CALL SET_WAVE_BOUNDARY_CONDITION

    DO IOB = 1, IOBCN_W
      IOB_NODE = I_OBC_N_W(IOB)
      AC2(:,:,IOB_NODE) = WBAC(:,:,IOB)
    END DO

    if(myrank==0) &
     WRITE(PRINTF,*)'SET_WAVE_BOUNDARY_CONDITION done'

!    ALLOCATE(AC2_TTMP(MDC,MSC,0:MT));  AC2_TTMP = 0.0_skind
!    CALL SET_VAR_WAVE_AC2(WaveTime,AC2=AC2_TTMP)

!    DO IOB = 1, IOBCN_W
!      IOB_NODE = I_OBC_N_W(IOB)
!      AC2(:,:,IOB_NODE) = AC2_TTMP(:,:,IOB_NODE)
!    END DO
!    DEALLOCATE(AC2_TTMP)

   ELSE
    MSGSTR = 'THE PARAMETER NESTING_TYPE_WAVE SHOULD BE PARAM, &
             WAMPARAM or SPEC:'
    WRITE(errmsg,'(a,1x,a)') MSGSTR,TRIM(NESTING_TYPE_WAVE)
    CALL parallel_abort(errmsg)
   END IF

!   if(myrank==0) print*,'SNEXTI myrank=',myrank,' update AC2 along obc ends'

  END IF


! -----------------------------------------------
! --- update input fields (wind, water level etc.)
! -----------------------------------------------
! WRITE (PRTEST, *) '--- update input fields '

 122  FORMAT (A, 10(1X,E11.4))
!
!     fields 5 and 6: wind
!
      IF (IFLDYN(5) .EQ. 1) THEN

        if(myrank==0) WRITE(PRINTF,*) 'Using wind field from schism'
!JL : Call JL implementation, nws = 2 only !
        CALL FLWIND (5, 6, 0, JWX2, JWX3, 0, JWY2, JWY3, IERR)

!JL DBG
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' wind from file in test points'
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JWX3)                          !40.80
             ENDDO                                                        !40.80
          WRITE (PRTEST, 122) ' X-comp: ',              &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JWY3)                          !40.80
             ENDDO                                                        !40.80
          WRITE (PRTEST, 122) ' Y-comp: ',             &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF

      ELSE
!        if(myrank==0) print*,'SNEXTI IMDTW=',IMDTW
        WaveTime = WaveTime + IMDTW  
!        if(myrank==0) print*,'SNEXTI WaveTime updated=',WaveTime
      ENDIF

!# JL See later
#if 0
      IF (ITEST.GE.0 .AND. myrank==0) THEN
       !call PRINT_REAL_TIME(WaveTime,ipt,"Wave Time",TIMEZONE)
       !  Get the date and time of the current iterationr
       SIMTIME = Write_DateTime(WaveTime,6,"UTC")
       WRITE(IPT,101) SIMTIME
      ENDIF
  101 FORMAT(1X,"! unswan",3X,A26)
#endif

!
!  field 4: friction coeff.
!
   IF (IFLDYN(4) .EQ. 1) THEN

   !JL Use SCHISM Friction (swan online with schism)
    if(myrank==0) WRITE(PRINTF,*) 'Using Friction from schism'

      IF (IBOT.EQ.3) THEN
         ! MADSEN Scheme for Bottom Dissipation, use KN, Roughness Length (J.LEFEVRE)
         DO INDX = 1, nverts
          COMPDA(INDX,JFRC2) = KN(INDX)  
          COMPDA(INDX,JFRC3) = COMPDA(INDX,JFRC2) ! No interpolation in time
         ENDDO
      ELSE
         MSGSTR = 'Only the MADSEN Scheme using the roughness lengh KN is supported, &
             change (IBOT=3):'
         WRITE(errmsg,'(a,1x,a)') MSGSTR,TRIM(NESTING_TYPE_WAVE)
         CALL parallel_abort(errmsg)
      ENDIF

   ! DBG
      IF (NPTST.GT.0) THEN
       WRITE (PRTEST, '(A)') ' friction from schism at test points'
             DO IPTST = 1, MIN(10,NPTST)
                IVP = XYTST(IPTST)
                TSTVAL(IPTST) = COMPDA(IVP,JFRC3)
             ENDDO
          WRITE (PRTEST, 122) ' friction: ',              &
      &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
      ENDIF
   ! DBG end

   ENDIF ! IFLDYN(4) .EQ. 1

!
!  field 7: water level
!
   IF (IFLDYN(7) .EQ. 1) THEN
   !JL Use SCHISM water level (swan online with schism)
   !# if defined (WAVE_CURRENT_INTERACTION)
    if(myrank==0) WRITE(PRINTF,*) 'Using water levels from schism'

    DO INDX = 1, nverts
       COMPDA(INDX,JWLV2) = eta2(INDX)
       COMPDA(INDX,JWLV3) = eta2(INDX)
    ENDDO

   ! DBG
    IF (NPTST.GT.0) THEN
       WRITE (PRTEST, '(A)') ' water level from schism at test points'
             DO IPTST = 1, MIN(10,NPTST)                                  
                IVP = XYTST(IPTST)                                        
                TSTVAL(IPTST) = COMPDA(IVP,JWLV2)                         
             ENDDO                                                        
          WRITE (PRTEST, 122) ' eta: ',              &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
    ENDIF
   ! DBG end


   ENDIF ! IFLDYN(7) .EQ. 1
!
!  Depth
!
! #JL Use SCHISM total water column depth (swan online with schism)
!
!   if(myrank==0) WRITE(PRINTF,*) 'Using bathymetry from schism'
    DO INDX = 1, nverts
        DEPTH(INDX) = dp(INDX)
        COMPDA(INDX,JDP1) = COMPDA(INDX,JDP2)
        COMPDA(INDX,JDP2) = max(ZERO,DEPTH(INDX)+eta2(INDX))
!         PRINT *,' SNEXTI donot forget to set JDP2 back'
!         COMPDA(INDX,JDP2) = HG(INDX) ! wrong
        IF (LSETUP.GT.0) THEN ! usefull for SETUP calc in SWAN
            COMPDA(INDX,JDPSAV) = COMPDA(INDX,JDP2)
        ENDIF
    ENDDO
!
!   field 2 and 3: current velocity
!
    IF (IFLDYN(2) .EQ. 1) THEN

    if(myrank==0) WRITE(PRINTF,*) 'Using current velocities from schism'

      COMPDA(:,JVX2) = curx_wwm(:)
      COMPDA(:,JVY2) = cury_wwm(:)

      DO INDX = 1, nverts
!       reduce current velocity if Froude number is larger than PNUMS(18)
!       structured grid
        DEPW = COMPDA(INDX,JDP2)
        !WRITE(6,*) 'JVX2,JVY2,DEPW',COMPDA(INDX,JVX2),COMPDA(INDX,JVY2),DEPW
        IF (DEPW.GT.0.0_rkind) THEN
           UU =  COMPDA(INDX,JVX2)
           VV =  COMPDA(INDX,JVY2)
           VTOT = SQRT (UU*UU + VV*VV)
             
           CGMAX = PNUMS(18)*SQRT(GRAV_W*DEPW)
!          IF (DEPW.LE.5.0_skind .AND. VTOT.GE.0.8_skind) THEN
!              WRITE(6,*) 'CGMAX,VTOT',CGMAX,VTOT
!          ENDIF
           IF (VTOT .GT. CGMAX) THEN
              CGFACT = CGMAX / VTOT
              COMPDA(INDX,JVX2) = UU * CGFACT
              COMPDA(INDX,JVY2) = VV * CGFACT
           ENDIF

        ENDIF
      ENDDO

     ! DBG
      IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' current from schism at test points'
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JVX2)                          !40.80
             ENDDO                                                        !40.80
          WRITE (PRTEST, 122) ' X-comp: ',              &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JVY2)                          !40.80
             ENDDO                                                        !40.80
          WRITE (PRTEST, 122) ' Y-comp: ',             &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
      ENDIF

    ENDIF ! IFLDYN(2) .EQ. 1 
!
!     field 11: field containing number of plants per square meter
!
      IF (IFLDYN(11) .EQ. 1) THEN
!JL Skip now, to implement/check later ... 
#if 0
        CALL FLFILE (11, 0, NPLAF, (/0./),       &
     &               0, JNPLA2, JNPLA3, 0, 0, 0, &
     &               1., 0.,                     &
     &               COMPDA, XCGRID, YCGRID,     & 
     &               KGRPNT, IERR)
!        IF (STPNOW()) RETURN

        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' # plants/m2 from file in test points'
          IF (OPTG.NE.5) THEN                                             !40.80
             DO IPTST = 1, MIN(10,NPTST)
               IXP = XYTST(2*IPTST-1)
               IYP = XYTST(2*IPTST)
               TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JNPLA3)
             ENDDO
          ELSE                                                            !40.80
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JNPLA3)                        !40.80
             ENDDO                                                        !40.80
          ENDIF                                                           !40.80
          WRITE (PRTEST, 122) ' veg dens: ', &
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
# endif
      ENDIF
!
!     field 12: field containing turbulent viscosity
!
      IF (IFLDYN(12) .EQ. 1) THEN
!JL : Rationale? Skip now, to implement/check later ...
# if 0
        CALL FLFILE (12, 0, TURBF, (/0./),      &
     &               0, JTURB2, JTURB3, 0, 0, 0,&
     &               1., 0.,                    &
     &               COMPDA, XCGRID, YCGRID,    &
     &               KGRPNT, IERR)              
!       IF (STPNOW()) RETURN
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' turb visc from file in test points'
          IF (OPTG.NE.5) THEN                                             !40.80
             DO IPTST = 1, MIN(10,NPTST)
               IXP = XYTST(2*IPTST-1)
               IYP = XYTST(2*IPTST)
               TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JTURB3)
             ENDDO
          ELSE                                                            !40.80
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JTURB3)                        !40.80
             ENDDO                                                        !40.80
          ENDIF                                                           !40.80
          WRITE (PRTEST, 122) ' turb visc: ',&
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
# endif
      ENDIF
!
!     field 13: field containing fluid mud layer
!
      IF (IFLDYN(13) .EQ. 1) THEN
!JL : Skip now, to implement/check later ...
# if 0
        CALL FLFILE (13, 0, MUDLF, (/0./),            &
     &               JMUDL1, JMUDL2, JMUDL3, 0, 0, 0, &
     &               1., 0.,                          &
     &               COMPDA, XCGRID, YCGRID,          &
     &               KGRPNT, IERR)
!       IF (STPNOW()) RETURN
        IF (NPTST.GT.0) THEN
          WRITE (PRTEST, '(A)') ' mud layer from file in test points'
          IF (OPTG.NE.5) THEN                                             !40.80
             DO IPTST = 1, MIN(10,NPTST)
               IXP = XYTST(2*IPTST-1)
               IYP = XYTST(2*IPTST)
               TSTVAL(IPTST) = COMPDA(KGRPNT(IXP,IYP),JMUDL3)
             ENDDO
          ELSE                                                            !40.80
             DO IPTST = 1, MIN(10,NPTST)                                  !40.80
                IVP = XYTST(IPTST)                                        !40.80
                TSTVAL(IPTST) = COMPDA(IVP,JMUDL3)                        !40.80
             ENDDO                                                        !40.80
          ENDIF                                                           !40.80
          WRITE (PRTEST, 122) ' mud layer: ',
     &          (TSTVAL(IPTST), IPTST=1,MIN(10,NPTST))
        ENDIF
# endif
      ENDIF
!
!      IF (LTRACE) CALL STRACE(IENT,'SNEXTI END')
!     End of subroutine SNEXTI
!
      RETURN
      END
!
!****************************************************************
!JL SUBROUTINE RBFILE Not used, Skip
#if 0
!
      SUBROUTINE RBFILE (SPCSIG, SPCDIR, BFILED, BSPLOC, &
     &                   BSPDIR, BSPFRQ, BSPECS, XYTST )               !40.31 30.90
!
!****************************************************************
!
      USE TIMECOMM                                                        !40.41
      USE OCPCOMM2                                                        !40.41
      USE OCPCOMM4                                                        !40.41
      USE SWCOMM1                                                         !40.41
      USE SWCOMM2                                                         !40.41
      USE SWCOMM3                                                         !40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.82: IJsbrand Haagsma
!     30.90: IJsbrand Haagsma (Equivalence version)
!     40.00: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.05: Ekaterini E. Kriezi
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.61: Roop Lalbeharry
!     41.13: Nico Booij
!
!  1. Updates
!
!     40.00, Nov. 97: new subroutine
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     30.82, Oct. 98: Updated description of several variables
!     40.03, Nov. 99: after label 380 BFILED(15) is replaced by BFILED(14)
!            May  00: in calls of DTRETI now BFILED(6) is used as time option
!     40.05, Aug. 00: WW3 nesting and changes in the form of the code
!                     (use of f90 features), Revision of subroutine
!     40.02, Oct. 00: Avoided REWIND of uninitialised unit number NDSD
!     40.13, Apr. 01: GOTO 392 added for a single boundary file (case NDSL=0)
!     40.13, May  01: read heading lines in case of WAM free format file
!                     changed
!     40.31, Nov. 03: removing POOL-mechanism, reconsidering this subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.61, Nov. 06: variables USNEW, THWNEW no longer written in WAM4.5
!     41.13, Jul. 10: LWDATE introduced (length of date/time in WAM nest)
!
!  2. Purpose
!
!     read boundary spectra from one file and additional information of the
!     heading lines
!
!  3. Methode
!
!     read from boundary files, aditional information (like time),
!     form the head lines per time step, and the head lines per point spectrum,
!     read the spectrum of the boundary file.
!     Transform to spectral resolution used in SWAN to obtain boundary spectra.
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL,   INTENT(IN)     ::  SPCDIR(MDC,6)                            !30.82
      REAL,   INTENT(IN)     ::  SPCSIG(MSC)                              !30.82
      REAL,   INTENT(INOUT)  ::  BSPECS(MDC,MSC,NBSPEC,2)
!
!     BSPDIR: Spectral directions of input spectrum
!     BSPFRQ: Spectral frequencies of input spectrum
!
      REAL   , INTENT(INOUT)  :: BSPDIR(*)
      REAL   , INTENT(INOUT)  :: BSPFRQ(*)
!
!     BFILED  data concerning boundary condition files
!     BSPLOC  place in array BSPECS where to store interpolated spectra
!     XYTST   test points
!
      INTEGER, INTENT(INOUT)  ::  BFILED(*)
      INTEGER, INTENT(INOUT)  ::  BSPLOC(*)
      INTEGER, INTENT(IN)     ::  XYTST(*)
!
!  5. Parameter variables
!
!     --
!
!  6. Local variables
!
      INTEGER   DORDER, NDSL, NDSD, IHD, IBOUNC, IBSPEC, IERR
      INTEGER   ID, IS, JJ, NANG, NFRE, COUNT_IT, IENT, II               !40.05
      INTEGER   WWDATE, WWTIME                                           !40.05
!
!     DORDER    if <0, order of reading directions is reversed
!     NDSL      unit ref num for namelist file
!     NDSD      unit ref num for data file
!     IHD       counter for heading lines
!     IBOUNC    counter for boundary locations
!     IBSPEC    counter for spectra
!     ID        counter for directions
!     IS        counter for frequencies
!     JJ        counter
!     NANG      number of directions on file
!     NFRE      number of frequencies on file
!     COUNT_IT  counter of the time entering in the WW3 boundary file
!
!     WWDATE, WWTIME       time code in WaveWatch                         40.05
!
      REAL      TIMF1, TIMF2, UFAC, W1, BDEPTH, DUM_A, RFAC               !40.05
      REAL      XLON, XLAT, XDATE, EMEAN, THQ, FMEAN
!
!     TIMF1     time of reading old boundary condition
!     TIMF2     time of reading new boundary condition
!     UFAC      multiplication factor
!     W1        weighting coefficient used in interpolation
!     XLON      longitude
!     XLAT      latitude
!     XDATE     date-time read from WAM file
!     EMEAN     coefficient read from WAM file, ignored
!     THQ       coefficient read from WAM file, ignored
!     FMEAN     coefficient read from WAM file, ignored
!     DUM_A     dummy local variable (not used in any calculation)
!     BDEPTH    depth of the boundary points
!
      DOUBLE PRECISION DDATE
!     DDATE     date-time read from WAM file
!
      LOGICAL   NSTATF, UNFORM
!     NSTATF    if True time appears in bound. cond. file
!     UNFORM    if True reading is done unformatted
!
      CHARACTER BTYPE *4, HEDLIN *80, TIMSTR *20, &
     &          DATITM(5) *18, PTNME*12
!     BTYPE     type of boundary condition
!     HEDLIN    heading line
!     TIMSTR    time string
      CHARACTER (LEN=14) :: CDATE     ! date-time
!
      REAL, ALLOCATABLE :: SPAUX(:,:)
!
!  8. SUBROUTINES USED
!
!     COPYCH, DTRETI, LSPLIT, SSHAPE, RESPEC, MSGERR, SINTRP
!
!  9. SUBROUTINES CALLING
!
!       SNEXTI
!
!  10. ERROR MESSAGES
!
!        ---
!
!  11. REMARKS
!
!
!  12. STRUCTURE
!
!       ---------------------------------------------------------
!       If file contains stationary wave data
!       Then If time2 < 0
!            Then Read boundary values from file
!                 Transform to spectral resolution used in SWAN to
!                 obtain boundary spectra
!                 Make time2 = + Inf
!            ----------------------------------------------------
!       Else Make time1 = timco - DT
!            Repeat
!                 If timco > time2
!                 Then Exit from repeat
!                 -----------------------------------------------
!                 Make time1 = time2
!                 Make old field values = new values
!                 Read new values from file
!            ----------------------------------------------------
!            Interpolate in time between old and new values
!            to update old values
!            Transform to spectral resolution used in SWAN to
!            obtain boundary spectra
!       ---------------------------------------------------------
!
! 13. SOURCE
!
!****************************************************************
!
!
      SAVE COUNT_IT                                                       !40.31
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'RBFILE')
!
!
!     if data file is exhausted, return
      IF (BFILED(1).EQ.-1) GOTO 900
      NSTATF = (BFILED(1) .GT. 0)
      IERR   = 0
      CALL COPYCH (BTYPE, 'F', BFILED(7), 1, IERR)
!
      IF (BTYPE.EQ.'WAMW' .OR. BTYPE.EQ.'WAMC') THEN
        UNFORM = .TRUE.
      ELSE
        UNFORM = .FALSE.
      ENDIF
!
      IF (BFILED(2).LT.0) COUNT_IT = 0                                    !40.05
!
      DORDER = BFILED(9)
      TIMF1 = REAL(BFILED(2))
      TIMF2 = REAL(BFILED(3))
      IF (ITEST.GE.120) WRITE (PRINTF, 187) BFILED(1), BTYPE,&
     &      TIMF1, TIMF2, TIMCO
 187  FORMAT (' Boundary', I2, 2X, A, ' times: ', 3F10.1)
!
      ALLOCATE(SPAUX(MDC,MSC))                                           !40.31
!
!     if present time > time of last set of spectra, read new spectra
!
!     While Loop named GLOOP
!
      GLOOP : DO                                                         !40.05
!
        IF (TIMCO.GT.TIMF2) THEN                                         !40.05
!
!     then read from boundary nesting files all the information
!     and the spectral ones
!
!     COUNT_IT : counter of the time which enter in a WW3  boundary file
!     and read spectral - used to calculate NHED (number of heading lines
!     per time) in WW3N case.
!
          COUNT_IT = COUNT_IT+1                                          !40.05
          NDSL = BFILED(4)
          NDSD = BFILED(5)
!
!     in WW3 case rewind the boundary file
          IF(BTYPE.EQ.'WW3N') REWIND(NDSD)                               !40.05
!
          TIMF1 = TIMF2
!
!         move new spectra to old for all boundary points
          DO IBOUNC = 1, BFILED(8)
            IBSPEC = BSPLOC(IBOUNC)

            DO ID = 1, MDC
              DO IS = 1, MSC
                BSPECS(ID,IS,IBSPEC,1) = BSPECS(ID,IS,IBSPEC,2)
              ENDDO
            ENDDO

            IF (ITEST.GE.80) WRITE (PRTEST, *) ' spectrum moved ',&
     &         IBSPEC, TIMF1
          ENDDO
!
 210       CONTINUE
!
!     define the new  NHED in every time step for WW3 case
          IF (BTYPE.EQ.'WW3N') THEN
            BFILED(15) = BFILED(14)+(1+((BFILED(16)-1)+&
     &                 CEILING((BFILED(12)* BFILED(10))/7.))*&
     &                 BFILED(8))*(COUNT_IT-1)                            !40.05
          ENDIF                                                           !40.05
!
!     read heading lines per time step
!
!     HBFL is loop over the number of heading lines per time step
          HBFL : DO IHD = 1, BFILED(15)                                   !40.05
            IF (UNFORM) THEN
              READ (NDSD, END=380, ERR=920)
            ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN             !40.05
              READ (NDSD,*)                                               !40.05
            ELSE                                                          !40.05
              READ (NDSD, '(A)', END=380, ERR=920) HEDLIN
              IF (ITEST.GE.90) WRITE (PRINTF, 212) HEDLIN
 212          FORMAT (' heading line: ', A)
              IF (BTYPE.EQ.'SWNT') THEN
!               convert time string to time in seconds
                CALL DTRETI (HEDLIN(1:18), BFILED(6), TIMF2)              !40.03
              ENDIF
            ENDIF
          ENDDO HBFL
!
          IF (.NOT.NSTATF) TIMF2 = 0.
          IF (ITEST.GE.60) WRITE (PRINTF, 214)&
     &        TIMF1, TIMF2, TIMCO, BFILED(8), BFILED(15)
 214      FORMAT (' Boundary times ', 3F12.0, 2X, 4I4)
!
!         read additional information from the headers and spectrum from
!         the nesting boundary files (for all the cases)
!
!
!       BP_LOOP loop over the boundary nesting points
          BP_LOOP : DO  IBOUNC = 1, BFILED(8)                             !40.05
!
            IBSPEC = BSPLOC(IBOUNC)
!           division by 2*PI to account for difference in definition of freq.
!           Hz to rad/s
            NANG  = BFILED(10)                                            !40.00
!
!           calculate UFAC for the different nesting cases
            IF (BTYPE(1:3).EQ.'SWN' .AND. NANG.GT.0) THEN
!           in addition multiply by 180/PI to account for directions in degr
!           instead of radians (SWAN 2D spectral files)
              UFAC = 180./ (2.*PI_W**2)
            ELSE                                                          !40.05
              UFAC = 1./ (2.*PI_W)
            ENDIF
!           divide by Rho*Grav if quantity in file is energy density
            IF ((BFILED(17).EQ.1))  UFAC = UFAC / (RHO_W*GRAV_W)
!
!           read information from heading lines per spectrum
!
!          do loop over the numbers of the heading lines per spectrum
!
            DO IHD = 1, BFILED(16)
!
              IF (IBOUNC.EQ.1) THEN                                       !40.03
!             for the first spectrum (first point), read time from
!             heading line                                                40.03
                IF (BTYPE.EQ.'WAMW') THEN                                 !40.03
                  READ(NDSD, END=380, ERR=920) XLON, XLAT,&
     &               CDATE(1:LWDATE),&                                    !41.13
     &               EMEAN, THQ, FMEAN
                  IF (LEN_TRIM(CDATE) == 10 ) THEN
                     TIMSTR = TRIM(CDATE)
                  ELSEIF (LEN_TRIM(CDATE) == 12 ) THEN
                     TIMSTR = CDATE(1:10)
                  ELSEIF (LEN_TRIM(CDATE) == 14 ) THEN
                     TIMSTR = CDATE(3:12)
                  ENDIF
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  !40.03 40.05
                ELSE IF (BTYPE.EQ.'WAMC') THEN
                  READ(NDSD, END=380, ERR=920) XLON, XLAT, XDATE,&
     &                EMEAN, THQ, FMEAN
                  WRITE (TIMSTR,'(F11.0,9X)') XDATE
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  !40.03 40.05
                ELSE IF (BTYPE.EQ.'WAMF') THEN
                  READ(NDSD,*, END=380, ERR=920) XLON, XLAT, DDATE,&
     &              EMEAN, THQ, FMEAN
                  WRITE (TIMSTR,'(F11.0,9X)') DDATE
!                 convert time string to time in seconds
                  CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                  !40.03 40.05
                ELSE IF (BTYPE.EQ.'WW3N') THEN                            !40.05
!                  read from heading lines per spectrum , date ,time      40.05
                  IF(IHD.EQ.1) THEN                                       !40.05
                    READ(NDSD,*, END=380, ERR=920) WWDATE, WWTIME         !40.05
                    WRITE (TIMSTR, 118) WWDATE, WWTIME                    !40.05
 118                FORMAT (I8,'.',I6)                                    !40.05
!                    convert time string to time in seconds
                    CALL DTRETI (TIMSTR, BFILED(6), TIMF2)                !40.05
                  ELSE                                                    !40.05
!                   DUM_A dummy local variable used to read formatted     40.05
!                   files                                                 40.05
!                   this variable will be not used in any calculation     40.05
                    READ (NDSD,901) PTNME, DUM_A, DUM_A, BDEPTH,&         !40.05
     &                              DUM_A, DUM_A, DUM_A, DUM_A            !40.05
                  ENDIF                                                   !40.05
                ELSE                                                      !40.05
!                 SWAN files                                              40.05
                  READ (NDSD, '(A)', END=380, ERR=920) HEDLIN             !40.05
                  IF (ITEST.GE.100) WRITE (PRINTF, 212) HEDLIN            !40.05
                ENDIF
              ELSE
                IF (UNFORM) THEN
                  READ (NDSD, END=380, ERR=920)
                ELSE IF (BTYPE.EQ.'WW3N') THEN                            !40.13 40.05
!                 read from heading lines per spectrum depth of the       40.05
!                 boundary point                                          40.05
                  IF (IHD.GT.1) EXIT                                      !40.05
!                 exit is used because the header per spectrum in WW3     40.05
!                 for IBOUNC>1 is one line                                40.05
!                 DUM_A is a dummy local parameter used in formatted read 40.05
                  READ (NDSD,901) PTNME, DUM_A, DUM_A, BDEPTH,&
     &                            DUM_A, DUM_A, DUM_A, DUM_A              !40.05
                ELSE IF (BTYPE.EQ.'WAMF') THEN
!                 read HEDLIN replaced because data are sometimes written 40.13
!                 on two subsequent lines                                 40.13
                  READ(NDSD,*, END=380, ERR=920) XLON, XLAT, DDATE,&      !40.13
     &              EMEAN, THQ, FMEAN
                ELSE
                  READ (NDSD, '(A)', END=380, ERR=920) HEDLIN
                  IF (ITEST.GE.100) WRITE (PRINTF, 212) HEDLIN
                ENDIF
              ENDIF
!
              IF (BTYPE(1:3).EQ.'SWN') THEN
!             SWAN nesting: take proper action if heading line contains ZERO or NODATA
                IF (HEDLIN(1:6).EQ.'NODATA' .OR. HEDLIN(1:4).EQ.'ZERO')&
     &               THEN
                  DO IS = 1, MSC
                    DO ID = 1, MDC
                      BSPECS(ID,IS,IBSPEC,2) = 0.
                    ENDDO
                  ENDDO
!                 skip reading of values
                  CYCLE BP_LOOP                                           !40.05
                ELSE IF (HEDLIN(1:6).EQ.'FACTOR') THEN
                  READ (NDSD, *) RFAC
                  UFAC = UFAC * RFAC
!                 multiply factor read from file by UFAC (factor following from
!                 type of file)
                ELSE
!                 note: in case of 1D spectra heading line can be ignored
                  IF (NANG.GT.0) THEN                                     !40.00
                    CALL MSGERR (3,&
     &                'incorrect code in b.c. file: '//HEDLIN(1:20))      !40.00
                  ENDIF                                                   !40.00
                ENDIF
              ENDIF
!           end loop over heading lines per spectrum
            ENDDO
!
!           test output: which spectrum is processed                      40.00
!
            IF (ITEST.GE.60) THEN
              INQUIRE (UNIT=NDSD, NAME=FILENM)                            !40.00
              WRITE (PRTEST, 188) FILENM, CHTIME, IBOUNC, UFAC
 188          FORMAT&
     &   (' read spectrum ', A, '; time=', A, ' nr=', I3, F9.3)           !40.00
            ENDIF
!
!       start reading incoming wave data
!
            IF (BTYPE.EQ.'TPAR') THEN
              READ (NDSD, 222, END=380) HEDLIN
 222          FORMAT (A)
              CALL LSPLIT (HEDLIN, DATITM, 5)
              CALL DTRETI (DATITM(1), BFILED(6), TIMF2)                   !40.03
              DO II = 1, 4
                READ (DATITM(II+1), '(G12.0)') SPPARM(II)
              ENDDO
              IF (ITEST.GE.60) WRITE (PRTEST, *) ' TPAR boundary ',&
     &              TIMF2, (SPPARM(JJ), JJ=1,4)
              BFILED(3) = NINT(TIMF2)
              CALL SSHAPE (BSPECS(1,1,IBSPEC,2), SPCSIG, SPCDIR,&
     &                 FSHAPE, DSHAPE)
            ELSE
!         other (spectral) boundary conditions
              NANG  = BFILED(10)
              NFRE  = BFILED(12)
!
!         call RESPEC subroutine to read the spectum of the bound. files
!
             CALL RESPEC (BTYPE, NDSD, BFILED, UNFORM, DORDER,&           !40.00
     &           SPCSIG, SPCDIR, BSPFRQ, BSPDIR, BSPECS(1,1,IBSPEC,2),&   !40.31 30.90
     &           UFAC, IERR)
             IF (IERR.EQ.9) GOTO 380
            ENDIF
!
          END DO BP_LOOP                                                 !40.05
!
          IF (ITEST.GE.60) WRITE (PRINTF, 287) BTYPE, TIMF2
 287      FORMAT&
     &         (' Boundary data type ', A, ' processed, time: ', F10.1)
!
!        cycle back to will loop
          CYCLE GLOOP                                                    !40.05
!
!         if there are no more data on a boundary data file
!         close this file, and see if there is a next one
!
 380      CLOSE(NDSD)
!         read filename of next boundary file and open them
          IF (NDSL.GT.0) THEN                                            !40.05
            READ (NDSL, '(A)', END=390, ERR=930) FILENM
            IF (UNFORM) THEN
              OPEN (NDSD, FILE=FILENM, FORM='UNFORMATTED',&
     &               STATUS='OLD', ERR=930)
!              read heading lines
              DO IHD = 1, BFILED(14)                                     !40.03
                READ (NDSD, END=940, ERR=920)
              ENDDO
            ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN            !40.05
!             if it is WW3 open the new file only
              OPEN (NDSD, FILE=FILENM, FORM='FORMATTED',&                !40.05
     &               STATUS='OLD', ERR=930)                              !40.05
              COUNT_IT = 1                                               !40.05
            ELSE                                                         !40.05
              OPEN (NDSD, FILE=FILENM, FORM='FORMATTED',&
     &               STATUS='OLD', ERR=930)
              DO IHD = 1, BFILED(14)                                     !40.03
!               read heading lines
                READ (NDSD, '(A)', END=940, ERR=920) HEDLIN
                IF (ITEST.GE.80) WRITE (PRINTF, 212) HEDLIN
              ENDDO
            ENDIF
!
!      go back to statement 210 to start again the procedure of reading
!      info and spectrum from the new boundary file
!
            GOTO 210
!
          ELSE
!           boundary data are read from a single file
            GOTO 392                                                      !40.13
          ENDIF                                                           !40.05
!         close file containing filenames
 390      CLOSE (NDSL)
!
          BFILED(5) = 0
!         write message and close file containing spectra
 392      CALL MSGERR (1, 'data on boundary file exhausted')
!
          BFILED(4) = 0
          BFILED(1) = -1
          TIMF2 = 999999999.
          CYCLE GLOOP                                                     !40.05
!
!        (if necessary) data have been read from file, now interpolate in time
!
        ELSE                                                              !40.05
!       if present time <= time of last set of spectra then
!       transform to spectral resolution used in SWAN to
!       obtain boundary spectra
!
          IF (TIMF1.NE.TIMF2) THEN                                        !40.41
             W1 = (TIMF2-TIMCO) / (TIMF2-TIMF1)
          ELSE
             W1 = 0.
          END IF
          DO IBOUNC = 1, BFILED(8)
            IBSPEC = BSPLOC(IBOUNC)
            IF (IBOUNC.EQ.1 .AND. ITEST.GE.80) WRITE (PRTEST, 403)&
     &           TIMCO, W1, TIMF1, TIMF2, IBSPEC
 403        FORMAT (' interp in time ', F14.1, F8.3, 2F14.1, I4)
!
!       interpolate spectra in time; result has to be store in BSPECS(..,1)
!       first interpolate to auxiliary array
!
            CALL SINTRP (W1, 1.-W1, BSPECS(1,1,IBSPEC,1),&
     &               BSPECS(1,1,IBSPEC,2), SPAUX,        &                !40.31 30.90
     &               SPCDIR, SPCSIG)
!       use SINTRP to copy contents of aux. array to BSPECS(..,1)
!
            CALL SINTRP (1., 0., SPAUX,&                                  !40.31 30.90
     &               BSPECS(1,1,IBSPEC,2), BSPECS(1,1,IBSPEC,1),&
     &               SPCDIR, SPCSIG)
          ENDDO
          BFILED(2) = NINT(TIMCO)
          BFILED(3) = NINT(TIMF2)
          EXIT GLOOP
!       end of time comparison
        ENDIF                                                             !40.05

!     end of while loop
      END DO GLOOP                                                        !40.05
!
!
 901  FORMAT (A12,2F7.2,F10.1,2(F7.2,F6.1))
!
      DEALLOCATE(SPAUX)                                                   !40.31
!
 900  RETURN
!
 920  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    !40.00
      CALL MSGERR (4,&
     &     'error reading data from boundary file '//FILENM)              !40.00
      GOTO 900
 930  CALL MSGERR (4,&
     &     'error opening boundary file '//FILENM)                        !40.00
      GOTO 900
 940  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    !40.00
      CALL MSGERR (4,&
     &     'unexpected end of file on boundary file '//FILENM)            !40.00
      GOTO 900
!
!     End of subroutine RBFILE
!
      END
# endif
!****************************************************************
!JL SUBROUTINE RESPEC Not used, Skip
#if 0
!
      SUBROUTINE RESPEC (BTYPE, NDSD, BFILED, UNFORM, DORDER, &           !40.00
     &                   SPCSIG, SPCDIR, BSPFRQ, BSPDIR, LSPEC, UFAC,&
     &                   IERR)
!
!****************************************************************
!
      USE TIMECOMM                                                        !40.41
      USE OCPCOMM2                                                        !40.41
      USE OCPCOMM4                                                        !40.41
      USE SWCOMM3                                                         !40.41
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.81: Annette Kieftenburg
!     40.00, 40.13: Nico Booij
!     40.02: IJsbrand Haagsma
!     40.05: Ekaterini E. Kriezi
!     40.31: Tim Campbell and John Cazes
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.61: Roop Lalbeharry
!
!  1. Update
!
!     40.00, Nov. 98: new subroutine
!     30.81, Feb. 99: approximation for MS > 10 corrected
!     40.02, Feb. 00: initialisation of ISIGTA
!     40.05, Aug. 00: WW3 nesting and changes in the form of the code
!                     (use of f90 features), Revise version of subroutine
!     40.02, Sep. 00: Made BAUX0 allocatable
!     40.13, Apr. 01: message concerning ISIGTA removed
!     40.31, Jul. 03: bug fix
!     40.31, Nov. 03: removing POOL construction
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.41, Nov. 04: small corrections
!     40.61, Nov. 06: boundary spectra are read as written in WAM4.5
!
!  2. Purpose
!
!     read one 1D OR 2D boundary spectrum from file, and transform
!     to internal SWAN spectral resolution
!
!  3. Method
!
!  4. Argument variables
!
      INTEGER, INTENT(INOUT)  :: BFILED(*)                                !40.05
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL,    INTENT(IN)     :: SPCDIR(MDC,6)                            !30.82
      REAL,    INTENT(IN)     :: SPCSIG(MSC)                              !30.82
      REAL,    INTENT(IN)     :: BSPDIR(*)
      REAL,    INTENT(IN)     :: BSPFRQ(*)
      REAL,    INTENT(INOUT)  :: LSPEC(MDC,MSC)
      REAL,    INTENT(IN)     :: UFAC                                     !40.05
!
      INTEGER,  INTENT(IN) ::  NDSD                                       !40.00
      INTEGER,  INTENT(IN) ::  DORDER                                     !40.00
      INTEGER,  INTENT(INOUT) ::  IERR                                    !40.00
!
      LOGICAL   UNFORM
!
      CHARACTER BTYPE *4
!
!     BTYPE    char  inp   type of input
!     NDSD     int   inp   unit ref. number of input file
!     BFILED   int   inp   options for reading boundary condition file    40.00
!     UNFORM   log   inp   if True, unformatted reading is called for
!     DORDER   int   inp   if <0, order of directions has to be reversed
!     NANG     int   inp   num of spectral direction of input spectrum
!     NFRE     int   inp   num of spectral frequencies of input spectrum
!     BSPFRQ   real  inp   spectral frequencies of input spectrum
!     BSPDIR   real  inp   spectral directions of input spectrum
!     LSPEC    real  out   interpolated spectrum
!     UFAC     real  inp   factor used to multiply data
!     IERR     int   out   error status, 0: no error, 9: end of file
!
!  5. Parameter variables
!
!  6. Local variables
!
!     IANG      counter of directions
!     IFRE      counter of frequencies
!     ID        counter of directions
!     IS        counter of frequencies
!     ISIGTA    the last frequency which is determined by interpolation
!
      INTEGER   IANG, IFRE, ID, IS, ISIGTA,IENT,NFRE,NANG
!
      REAL, ALLOCATABLE :: BAUX0(:,:)                                     !40.05
      REAL, ALLOCATABLE :: BAUX1(:,:)                                     !40.31
      REAL, ALLOCATABLE :: BAUX2(:,:)                                     !40.31
      REAL, ALLOCATABLE :: BAUX3(:)                                       !40.31
      REAL, ALLOCATABLE :: BAUX4(:)                                       !40.31

      REAL      ETOT, ADEG, DD, ADIR, MS, CTOT, ACOS, CDIR
      REAL      DSUM, DSPR
!
!     ETOT      energy integrated over directions
!     ADEG      average direction in degr
!     DD        parameter for directional distribution
!     ADIR      average direction in rad
!     MS        power of Cos in directional distribution
!     CTOT      coefficient
!     DSPR      directional spread in rad
!     ACOS      cos of angle between direction and average dir.
!     CDIR      energy in one directional bin
!     BAUX0     auxiliary array for energy density
!     BAUX1     auxiliary array 1                                         40.31
!     BAUX2     auxiliary array 2                                         40.31
!     BAUX3     auxiliary array 3                                         40.31
!     BAUX4     auxiliary array 4                                         40.31
!
!  8. Subroutines used
!
!       GAMMAF (in SWANSER)
      REAL :: GAMMAF                                                      !40.03
      LOGICAL EQREAL                                                      !40.41
!
!  9. Subroutines calling
!
!     RBFILE
!
!  10. Error messages
!
!        ---
!
!  11. Remarks
!
!
!  12. Structure
!
!       ---------------------------------------------------------
!       for all frequencies of input spectrum do
!           if NANG = 0
!           then read energy density, av. direction and dir. spread
!                determine directional distribution
!           else read spectral energy densities (1 .. NANG)
!                redistribute to get densities for SWAN directions
!       ---------------------------------------------------------
!       for all spectral directions do
!           redistribute te get densities for SWAN frequencies
!       ---------------------------------------------------------
!
!  13. Source text
!
!****************************************************************
!
      SAVE      IENT
      DATA      IENT /0/
      CALL STRACE (IENT, 'RESPEC')
!
!     read the values from array BFILED for the number of
!     direction and frequencies
!
      NANG = BFILED(10)                                                   !40.05
      NFRE = BFILED(12)                                                   !40.05
      ALLOCATE (BAUX0(NFRE,NANG))                                         !40.02
      ALLOCATE (BAUX1(NANG,NFRE))                                         !40.31
      ALLOCATE (BAUX2(MDC ,NFRE))                                         !40.31
      ALLOCATE (BAUX3(NFRE))                                              !40.31
      ALLOCATE (BAUX4(MSC))                                               !40.31

      IF (ITEST.GE.60) THEN
        INQUIRE (UNIT=NDSD, NAME=FILENM)                                  !40.00
        WRITE (PRTEST, 7) FILENM, NANG, NFRE, UFAC, BFILED(18),&          !40.00
     &  BFILED(19), UNFORM
   7    FORMAT (' Entry RESPEC, reading file:', A, /, 6X,&
     &  I3, ' angles ', I3, ' freqs; factor ',&
     &  E12.4, ' dir.def:', I2, ' dir.spr.def:', I2, ' unform:', L2)      !40.00
      ENDIF
!
!     Initialisation
!
      ISIGTA = MSC                                                        !40.02
!
!     read spectral energy densities from b.c. file
      IF (NANG.EQ.0) THEN
!
!     1-D spectral input                                                  40.05
!
        DO IFRE=1,NFRE                                                    !40.05
!         1-D spectral input (only function of frequency)
          IF (IFRE.EQ.1) THEN
            READ (NDSD,*,END=940,ERR=920) ETOT, ADEG, DD
          ELSE
            READ (NDSD,*,END=930,ERR=920) ETOT, ADEG, DD
          ENDIF
!
          IF (EQREAL(ETOT,REAL(BFILED(11)))) THEN                        !40.41
!            in case of exception value, no energy
             ETOT = 0.
             ADEG = 0.
             DD   = 180./PI_W
          END IF
!
          IF (BFILED(18).EQ.1) THEN                                      !40.00
!           conversion from degrees to radians
            ADIR = PI_W * ADEG / 180.
          ELSE
!           conversion from Nautical to Cartesian conv.
            ADIR = PI_W * (180.+DNORTH-ADEG) / 180.
          ENDIF
!
          IF (BFILED(19).EQ.1) THEN                                       !40.00
!           DSPR is directional spread in radians
            DSPR = PI_W * DD / 180.
            IF (DSPR.NE.0.) THEN                                          !40.41
               MS = MAX (DSPR**(-2) - 2., 1.)
            ELSE
               MS = 1000.                                                 !40.41
            END IF
          ELSE
            MS = DD
          ENDIF
!
          IF (ITEST.GE.80) WRITE (PRTEST, 12) IFRE, ETOT,&                !40.00
     &          180.*ADIR/PI_W, MS
  12      FORMAT (' read freq ', I3, E10.3, '; Cart dir ',&               !40.00
     &          F7.1, '; Cos power ', F7.2)
!
!         generate distribution over directions
!
!         equations taken from Jahnke & Emde (chapter Factorial Function)
          IF (MS.GT.10.) THEN
            CTOT = SQRT(MS/(2.*PI_W)) * (1. + 0.25/MS)                      !30.81 40.00
          ELSE
            CTOT = 2.**MS * (GAMMAF(1.+0.5*MS))**2 / (PI_W*GAMMAF(1.+MS))   !40.00
          ENDIF
          DSUM = 0.
          DO ID = 1, MDC
            ACOS = COS(SPCDIR(ID,1) - ADIR)                               !30.82
            IF (ACOS .GT. 0.) THEN
              CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
            ELSE
              CDIR = 1.E-10
            ENDIF
            IF (ITEST.GE.20) DSUM = DSUM + CDIR * DDIR                    !40.00
            BAUX2(ID,IFRE) = CDIR * ETOT
          ENDDO
          IF (ITEST.GE.20) THEN
            IF (ABS(DSUM-1.).GT.0.1) WRITE (PRTEST, 138) DSUM, CTOT, MS   !40.00
 138        FORMAT (' integral over directions is ', F9.4,&
     &              ' with CTOT=', F10.3,'; power=', F8.2)                !40.00
          ENDIF
        END DO                                                            !40.05
!
      ELSE
!
!   (2D) fully spectral input
        IF (UNFORM) THEN
!          unformatted reading
           IF (DORDER.LT.0) THEN
              READ (NDSD,END=930,ERR=920)&
     &               ((BAUX1(IANG,IFRE),IANG=NANG,1,-1),IFRE=1,NFRE)
           ELSE
              READ (NDSD,END=930,ERR=920)&
     &               ((BAUX1(IANG,IFRE),IANG=1,NANG),IFRE=1,NFRE)
           ENDIF
        ELSEIF ((.NOT.UNFORM).AND.(BTYPE.EQ.'WW3N')) THEN                 !40.05
!
!         WW3 reading
!         BAUX0 local array to read the energy spectra from boundary files
!
          READ (NDSD,902,END=940,ERR=920)&                                !40.15
     &      ((BAUX0(IFRE,IANG),IFRE=1,NFRE),IANG=1,NANG,1)                !40.15
 902      FORMAT (7E11.3)                                                 !40.15
!
!         energy (variance) density   from E(FRQ,TH) to  E(TH,FRQ)
!
          DO IANG = 1,NANG                                                !40.05
            DO IFRE = 1,NFRE                                              !40.05
              BAUX1(IANG,IFRE) = BAUX0(IFRE,IANG)                         !40.05
            ENDDO                                                         !40.05
          ENDDO                                                           !40.05
        ELSE
!         format reading (except WW3)
          IF (DORDER.LT.0) THEN                                           !40.31
            READ (NDSD,*,END=930,ERR=920)  &                              !40.61
     &             ((BAUX1(IANG,IFRE),IANG=NANG,1,-1),IFRE=1,NFRE)        !40.61
          ELSE                                                            !40.31
            READ (NDSD,*,END=930,ERR=920)  &                              !40.61
     &             ((BAUX1(IANG,IFRE),IANG=1,NANG),IFRE=1,NFRE)           !40.61
          ENDIF                                                           !40.31
        ENDIF
!
        IF (ITEST.GE.120) THEN
          WRITE (PRINTF,*)' Spectra from file'
          DO IFRE = 1, NFRE
            WRITE (PRINTF,*) IFRE, (BAUX1(IANG,IFRE),IANG=1,NANG)
          ENDDO
        ENDIF

!       --- in case of exception value, no energy                         40.41

        DO IANG = 1,NANG
           DO IFRE = 1,NFRE
              IF (EQREAL(BAUX1(IANG,IFRE),REAL(BFILED(11)))) THEN
                 BAUX1(IANG,IFRE) = 0.
              END IF
           END DO
        END DO
!
!       --- transform to spectral directions used in SWAN
!           results appear in array BAUX2(MDC,NFRE)
!
        DO IFRE = 1, NFRE                                                 !40.05
          CALL CHGBAS (BSPDIR, SPCDIR, PI2_W, BAUX1(1,IFRE),   &
     &                 BAUX2(1,IFRE), NANG, MDC, ITEST, PRTEST)
        END DO                                                            !40.05
      ENDIF
!
!     interpolate energy densities to SWAN frequencies distribution
!
      IF (BSPFRQ(NFRE) .LT. SPCSIG(MSC)) THEN
        DO IS = MSC, 1, -1
          IF (SPCSIG(IS).LT.BSPFRQ(NFRE)) THEN
!           ISIGTA is the last frequency which is determined by interpolation
!           higher frequencies are determined by tail expression
            ISIGTA = IS
            EXIT                                                          !40.05
          ENDIF
        ENDDO
      ELSE
        ISIGTA = MSC
      ENDIF
!
!
!     UFAC is the product of the multiplication factor read from file
!     and the factor to transform from energy/Hz to energy/(rad/s)
!     and from energy/degr to energy/rad (latter only for 2d spectra)     40.00
!
      DO  ID = 1,MDC                                                      !40.05
        DO  IFRE = 1,NFRE                                                 !40.05
          BAUX3(IFRE) = UFAC * BAUX2(ID,IFRE)
        ENDDO                                                             !40.05

!       interpolate over frequency keeping energy constant, output BAUX4(MSC)
        CALL CHGBAS (BSPFRQ, SPCSIG, 0., BAUX3, BAUX4, NFRE, MSC, &
     &               ITEST, PRTEST)
        DO  IS=1,MSC                                                      !40.05
          IF (IS.LE.ISIGTA) THEN
!
!           to convert energy density to action density
!
            LSPEC(ID,IS) = BAUX4(IS)/SPCSIG(IS)
          ELSE
!
!           add a tail when IS > ISIGTA
!
            LSPEC(ID,IS) = LSPEC(ID,ISIGTA) * &
     &                 (SPCSIG(ISIGTA)/SPCSIG(IS))**(PWTAIL(1)+1)
          ENDIF
          IF (ITEST.GE.140) THEN
            WRITE (PRTEST, *) 'ID,IS,LSPEC(ID,IS)', ID,IS,LSPEC(ID,IS)
          ENDIF
        ENDDO                                                             !40.05
      ENDDO                                                               !40.05
!
      DEALLOCATE (BAUX0,BAUX1,BAUX2,BAUX3,BAUX4)                          !40.31 40.02

 900  IERR = 0
      RETURN
 920  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    !40.00
      CALL MSGERR (2, &
     &      'read error in boundary condition file '//FILENM)             !40.00


      RETURN
 930  INQUIRE (UNIT=NDSD, NAME=FILENM)                                    !40.00
      CALL MSGERR (2, &
     &      'insufficient data in boundary condition file '//FILENM)      !40.00

 940  IERR = 9

      RETURN

      END SUBROUTINE RESPEC
# endif

!**********************************************************************
!JL  Subroutine FLFILE Not used, Skip
# if 0

      SUBROUTINE FLFILE (IGR1, IGR2,                               &
     &                   ARR, ARR2, JX1, JX2, JX3, JY1, JY2, JY3,  &     !40.31
     &                   COSFC, SINFC, COMPDA,                     &     !40.31 30.90
     &                   XCGRID, YCGRID,                           &   
     &                   KGRPNT, IERR)
!
!**********************************************************************

      USE TIMECOMM                                                       !40.41
      USE OCPCOMM4                                                       !40.41
      USE SWCOMM2                                                        !40.41
      USE SWCOMM3                                                        !40.41
      USE M_PARALL                                                       !40.31
      USE SwanGriddata                                                   !40.80
!ADC      USE Couple2Swan, ONLY: ADCIRC_ETA2 => SWAN_ETA2,                    41.20
!ADC     &                       ADCIRC_UU2 => SWAN_UU2,
!ADC     &                       ADCIRC_VV2 => SWAN_VV2,
!ADC     &                       ADCIRC_WX2 => SWAN_WX2,
!ADC     &                       ADCIRC_WY2 => SWAN_WY2,
!ADC     &                       COUPCUR, COUPWIND, COUPWLV,
!ADC     &                       InterpoWeight
!ADC     &                      ,ADCIRC_Z0 => SWAN_Z0,
!ADC     &                       COUPFRIC

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.90: IJsbrand Haagsma (Equivalence version)
!     40.00: Nico Booij
!     34.01: Jeroen Adema
!     40.02: IJsbrand Haagsma
!     40.03, 40.13: Nico Booij
!     40.30: Marcel Zijlema
!     40.31: Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!     41.20: Casey Dietrich
!
!  1. Updates
!
!     40.00, Jan. 98: new subroutine replacing code in subr SNEXTI
!     30.90, Oct. 98: Introduced EQUIVALENCE POOL-arrays
!     34.01, Feb. 99: Introducing STPNOW
!     40.03, Aug. 00: condition added for calling INAR2D to prevent error
!                     in case command INP GRID is present and corresponding
!                     command READ is not.
!     40.02, Oct. 00: Avoided real/int conflict by replacing RPOOL for POOL in
!                     INAR2D
!     40.13, Mar. 01: misplaced error message moved to proper place
!     40.30, Mar. 03: introduction distributed-memory approach using MPI
!     40.31, Nov. 03: removing POOL-mechanism
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!     41.20, Mar. 10: extension to tightly coupled ADCIRC+SWAN model
!
!  2. PURPOSE
!
!     Update boundary conditions, update nonstationary input fields
!
!  3. METHOD
!
!
!  4. Argument list
!
!     ARR      real  i  array holding values read from file (x-comp)
!     ARR2     real  i  array holding values read from file (y-comp)
!     TIMR2    real i/o time of last reading of input field
!     INTRV    real  i  time interval between input fields
!     TMENDR   real  i  end time of input field
!     IGR1     int   i  location in array COMPDA for interpolated input field data (x-comp)
!     IGR2     int   i  location in array COMPDA for interpolated input field data (y-comp)
!                       for a scalar field IGR2=0
!     JX1      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX2      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JX3      int   i  location in array COMPDA for interpolated input field data (x-comp)
!     JY1      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY2      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     JY3      int   i  location in array COMPDA for interpolated input field data (y-comp)
!     COSFC    real  i  cos of angle between input grid and computational grid
!     SINFC    real  i  sin of angle between input grid and computational grid
!     COMPDA   real i/o array holding values for computational grid points
!     XCGRID   real  i  x-coordinate of computational grid points
!     YCGRID   real  i  y-coordinate of computational grid points
!     KGRPNT   int   i  indirect addresses of computational grid points
!     NHDF     int   i  number of heading lines for a data file
!     NHDT     int   i  number of heading lines per time step
!     NHDC     int   i  number of heading lines before second component of vector field
!     IDLA     int   i  lay-out identifier for a data file
!     IDFM     int   i  format identifier for a data file
!     DFORM    char  i  format to read a data file
!     VFAC     real  i  multiplication factor applied to values from data file
!     IERR     int   o  error status: 0=no error, 9=end-of-file
!
!
!  5. SUBROUTINES CALLING
!
!     SNEXTI
!
!  6. SUBROUTINES USED
!
!     INAR2D
!     MSGERR
!     STRACE
!     SWBROADC
!
      LOGICAL STPNOW                                                     !34.01
!
!  7. ERROR MESSAGES
!
!        ---
!
!  8. REMARKS
!
!
!  9. Structure
!
!     --------------------------------------------------------------
!     for all comp. grid points do
!         copy new values to old
!     --------------------------------------------------------------
!     repeat
!         if present time > time of last reading
!         then read new values from file
!              update time of last reading
!              interpolate values to computational grid
!         else exit from repeat
!     --------------------------------------------------------------
!     for all comp. grid points do
!         interpolate new values
!     --------------------------------------------------------------
!
! 10. SOURCE
!
!****************************************************************
!
      INTEGER    KGRPNT(MXC,MYC), &
     &           IGR1, IGR2, JX1, JX2, JX3, JY1, JY2, JY3, IERR
!
      REAL       COMPDA(MCGRD,MCMVAR),            &
     &           XCGRID(MXC,MYC), YCGRID(MXC,MYC),&
     &           COSFC, SINFC
      REAL       ARR(*), ARR2(*)
!
!     local variables
!
      INTEGER    IENT, INDX, IX, IY
!     INDX       counter of comp. grid points
!     IX         index in x-dir of comput grid point
!     IY         index in y-dir of comput grid point
!
      REAL       SVALQI
!     SVALQI     real function giving interpolated value of an input array
!
      REAL       TIMR1, XP, YP, UU, VV, VTOT, W1, W3,&
     &           SIZE1, SIZE2, SIZE3
!     TIMR1      time of one but last input field
!     XP         x-coord of one comput grid point
!     YP         y-coord of one comput grid point
!     UU         x-component of vector, or scalar value
!     VV         y-component of vector
!     VTOT       length of vector
!     W1         weighting coeff for interpolation in time
!     W3         weighting coeff for interpolation in time
!     DIRE       direction of interpolated vector
!     SIZE1      length of vector at time TIMR1
!     SIZE2      length of vector at time TIMCO
!     SIZE3      length of vector at time TIMR2
!
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'FLFILE')
!
      IERR = 0
!
      IF (JX1.GT.1) THEN
        DO INDX = 1, MCGRD
          COMPDA(INDX,JX1)=COMPDA(INDX,JX2)
        ENDDO
      ENDIF
      IF (IGR2.GT.0 .AND. JY1.GT.1) THEN
        DO INDX = 1, MCGRD
          COMPDA(INDX,JY1)=COMPDA(INDX,JY2)
        ENDDO
      ENDIF
      TIMR1 = TIMCO - DTW
!
 200  IF (TIMCO.LE.IFLTIM(IGR1)) GOTO 400
      TIMR1 = IFLTIM(IGR1)
      IFLTIM(IGR1) = IFLTIM(IGR1) + IFLINT(IGR1)
      IF (IFLTIM(IGR1) .GT. IFLEND(IGR1)) THEN
        IFLTIM(IGR1) = 1.E10
        IF (IGR2.GT.0) IFLTIM(IGR2) = IFLTIM(IGR1)
        GOTO 400
      ENDIF
      IF (IFLNDS(IGR1).GT.0) THEN                                         !40.03
        IF (INODE.EQ.MASTER) THEN
!ADC!          --- if we have coupled to the quantities from ADCIRC,          41.20
!ADC!              then grab them from memory instead of reading them
!ADC!              from the external file
!ADC           IF ( (IGR1.EQ.7).AND.COUPWLV ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_ETA2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_ETA2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.2).AND.COUPCUR ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_UU2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_UU2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.5).AND.COUPWIND ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_WX2(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_WX2(INDX,2))
!ADC              ENDDO
!ADC           ELSEIF ( (IGR1.EQ.4).AND.COUPFRIC ) THEN
!ADC              DO INDX = 1, nverts
!ADC                 ARR(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                  * REAL(ADCIRC_Z0(INDX,1))
!ADC     &                  + REAL(InterpoWeight)
!ADC     &                  * REAL(ADCIRC_Z0(INDX,2))
!ADC              ENDDO
!ADC           ELSE
           CALL INAR2D( ARR, MXG(IGR1), MYG(IGR1),               &        !40.31 40.02
     &                  IFLNDF(IGR1),                            &
     &                  IFLNDS(IGR1), IFLIFM(IGR1), IFLFRM(IGR1),&
     &                  IFLIDL(IGR1), IFLFAC(IGR1),              & 
     &                  IFLNHD(IGR1), IFLNHF(IGR1))
!           IF (STPNOW()) RETURN                                           !34.01
!ADC           ENDIF
        END IF
        CALL SWBROADC(IFLIDL(IGR1),1,SWINT)                               !40.30
        IF (IFLIDL(IGR1).LT.0) THEN
!         end of file was encountered
          IFLTIM(IGR1) = 1.E10
          IF (IGR2.GT.0) IFLTIM(IGR2) = IFLTIM(IGR1)
          GOTO 400
        ELSE
          CALL SWBROADC(ARR,MXG(IGR1)*MYG(IGR1),SWREAL)                   !40.31 40.30
        ENDIF
      ELSE                                                                !40.13
        IF (ITEST.GE.20) THEN                                             !40.13
          CALL MSGERR (1, &
     &    'no read of input field because unit nr=0')                     !40.13
          WRITE (PRINTF, 208) IGR1
 208      FORMAT (' field nr.', I2)
        ENDIF
      ENDIF                                                               !40.03
      IF (IGR2.GT.0) THEN
        IFLTIM(IGR2) = IFLTIM(IGR1)
        IF (IFLNDS(IGR2).GT.0) THEN                                       !40.03
          IF (INODE.EQ.MASTER) THEN                                       !40.30
!ADC!            added these lines to grab the y-components from memory       41.20
!ADC             IF ( (IGR2.EQ.3).AND.COUPCUR ) THEN
!ADC               DO INDX = 1, nverts
!ADC                 ARR2(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                   * REAL(ADCIRC_VV2(INDX,1))
!ADC     &                   + REAL(InterpoWeight)
!ADC     &                   * REAL(ADCIRC_VV2(INDX,2))
!ADC               ENDDO
!ADC             ELSEIF ( (IGR2.EQ.6).AND.COUPWIND ) THEN
!ADC               DO INDX = 1, nverts
!ADC                 ARR2(INDX) = (1.0 - REAL(InterpoWeight))
!ADC     &                   * REAL(ADCIRC_WY2(INDX,1))
!ADC     &                   + REAL(InterpoWeight)
!ADC     &                   * REAL(ADCIRC_WY2(INDX,2))
!ADC               ENDDO
!ADC             ELSE
             CALL INAR2D( ARR2, MXG(IGR2), MYG(IGR2),              &      !40.31 40.02
     &                    IFLNDF(IGR2),                            &
     &                    IFLNDS(IGR2), IFLIFM(IGR2), IFLFRM(IGR2),&
     &                    IFLIDL(IGR2), IFLFAC(IGR2), IFLNHD(IGR2), 0)
!            IF (STPNOW()) RETURN                                         !34.01
!ADC             ENDIF
          END IF
          CALL SWBROADC(ARR2,MXG(IGR2)*MYG(IGR2),SWREAL)                  !40.31 40.30
        ENDIF                                                             !40.03
      ENDIF
!     Interpolation over the computational grid
!JL Skip Now
# if 0
!     structured grid
      DO 230 IX = 1, MXC
        DO 240 IY = 1, MYC
          INDX = KGRPNT(IX,IY)
          IF (INDX.GT.1) THEN                                             !40.00
            XP = XCGRID(IX,IY)
            YP = YCGRID(IX,IY)
            UU = SVALQI (XP, YP, IGR1, ARR, 0, IX, IY)                    !40.31 30.90
            IF (IGR2.EQ.0) THEN
              COMPDA(INDX,JX3) = UU
            ELSE
              VV = SVALQI (XP, YP, IGR2, ARR2, 0, IX, IY)                 !40.31 30.90
              COMPDA(INDX,JX3) =  UU*COSFC + VV*SINFC
              COMPDA(INDX,JY3) = -UU*SINFC + VV*COSFC
            ENDIF
          ENDIF
 240    CONTINUE
 230  CONTINUE
# endif
!     unstructured grid
      DO INDX = 1, nverts                                                 !40.80
         XP = xcugrd(INDX)
         YP = ycugrd(INDX)
         IF ( IGTYPE(IGR1).EQ.3 ) THEN
            UU = ARR(INDX)
         ELSE
            UU = SVALQI (XP, YP, IGR1, ARR, 0, 0, 0)
         ENDIF
         IF (IGR2.EQ.0) THEN
            COMPDA(INDX,JX3) = UU
         ELSE
            IF ( IGTYPE(IGR2).EQ.3 ) THEN
               VV = ARR2(INDX)
            ELSE
               VV = SVALQI (XP, YP, IGR2, ARR2, 0, 0, 0)
            ENDIF
            COMPDA(INDX,JX3) =  UU*COSFC + VV*SINFC
            COMPDA(INDX,JY3) = -UU*SINFC + VV*COSFC
         ENDIF
      ENDDO                                                               !40.80
      GOTO 200
!
!         Interpolation in time
!
 400  W3 = (TIMCO-TIMR1) / (IFLTIM(IGR1)-TIMR1)
      W1 = 1.-W3
      IF (ITEST.GE.60) WRITE(PRTEST,402) IGR1,&
     &        TIMCO,IFLTIM(IGR1),W1,W3,JX1,JY1,JX2,JY2,JX3,JY3
 402  FORMAT (' input field', I2, ' interp at ', 2F9.0, 2F8.3, 6I3)
      DO 500 INDX = 1, MCGRD
        UU = W1 * COMPDA(INDX,JX2) + W3 * COMPDA(INDX,JX3)
        IF (IGR2.LE.0) THEN
          COMPDA(INDX,JX2) = UU
        ELSE
          VV = W1 * COMPDA(INDX,JY2) + W3 * COMPDA(INDX,JY3)
          VTOT = SQRT (UU*UU + VV*VV)
!
!         procedure to prevent loss of magnitude due to interpolation
!
          IF (VTOT.GT.0.) THEN
            SIZE1 = SQRT(COMPDA(INDX,JX2)**2 + COMPDA(INDX,JY2)**2)
            SIZE3 = SQRT(COMPDA(INDX,JX3)**2 + COMPDA(INDX,JY3)**2)
            SIZE2 = W1*SIZE1 + W3*SIZE3
!           SIZE2 is to be length of vector
            COMPDA(INDX,JX2) = SIZE2*UU/VTOT
            COMPDA(INDX,JY2) = SIZE2*VV/VTOT
          ELSE
            COMPDA(INDX,JX2) = UU
            COMPDA(INDX,JY2) = VV
          ENDIF
        ENDIF
 500  CONTINUE
      RETURN
!
!     End of subroutine FLFILE
      END
# endif
!
!***********************************************************************
!                                                                      *
!      SUBROUTINE SWINCO (AC2    ,COMPDA ,  &
      SUBROUTINE SWINCO (                  &
     &                   XCGRID ,YCGRID ,  &                             !30.72
     &                   KGRPNT ,SPCDIR ,  &
     &                   SPCSIG ,XYTST   )                               !30.72
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                       !40.41
      USE SWCOMM2                                                        !40.41
      USE SWCOMM3                                                        !40.41
      USE SWCOMM4                                                        !40.41
!      USE M_PARALL                                                       !40.31
      USE SwanGriddata                                                   !40.80

!      USE ALL_VARS
!      USE VARS_WAVE, ONLY : SERIAL,PAR,ART,   &
!                     IOBCN_W,I_OBC_N_W,AC2,COMPDA
      USE schism_glbl, ONLY : skind,rkind,dkind,ne_global,area,ipgl,iplg,ne,npa
      USE VARS_WAVE, ONLY : SERIAL,AC2,COMPDA

      !USE schism_msgp, ONLY : myrank,nproc,ierr,istatus,comm
      !USE schism_msgp, ONLY : rtype,itype
      USE schism_msgp

!# if defined (MULTIPROCESSOR)
!      USE MOD_PAR
!# endif

!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif


!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: The SWAN team                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     30.60: Nico Booij
!     30.70: Nico Booij
!     30.72: IJsbrand Haagsma
!     30.80, 40.13: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!            June 97: new for SWAN
!     30.60, Aug. 97: for zero wind velocity no change in action density
!                     modification to make procedure work for uniform wind
!                     maximum set to dim.less fetch loops over IS and ID
!                     swapped for efficiency
!     30.70, Sep. 97: output for test point added, argument XYTST added
!     30.72, Sept 97: Replaced DO-block with one CONTINUE to DO-block with
!                     two CONTINUE's
!     30.70, Feb. 98: computation of initial values revised argument list added
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!     30.80, Apr. 98: correction computation FDLSS
!     30.82, Apr. 98: Modified computation of FDLSS, FPDLSS, HSDLSS
!     30.82, Oct. 98: Updated description of several variables
!     40.13, Feb. 01: correction for 1D cases
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.80, Sep. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Imposing of wave initial conditions at a computational grid
!
!  3. Method
!
!     The initial conditions are given using the following equation       30.70
!     for dimensionless Hs as function of dimensionless fetch:            30.70
!
!     Hs = 0.00288 f**(0.45)                                              40.00
!     Tp = 0.46    f**(0.27)                                              40.00
!     average direction = wind direction
!     directional distribution: Cos**2
!
!     after computation of the integral parameters the subroutine SSHAPE  30.70
!     is used to compute the spectrum
!
!  4. Argument variables
!
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
! i   XCGRID: Coordinates of computational grid in x-direction            30.72
! i   YCGRID: Coordinates of computational grid in y-direction            30.72
!
      REAL(rkind)    SPCDIR(MDC,6)                                              !30.82
      REAL(rkind)    SPCSIG(MSC)                                                !30.72
      REAL(rkind)    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                        !30.72
!
!     AC2        real  i/o   action density spectra
!     COMPDA     real  inp   quantities in comp grid points
!     KGRPNT     real  inp   indirect addresses of comp grid points
!     XYTST      int   inp   test points
!
!  6. Local variables
!
!     FDLSS : Dimensionless fetch
!     TPDLSS: Dimensionless peak period
!     HSDLSS: Dimensionless significant wave height
!
      REAL(rkind)    FDLSS,  TPDLSS, HSDLSS                                     !30.80

      REAL(rkind) :: TLEN,FETCH,TDXY
!
!  7. Common blocks used
!
!
!  8. REMARKS
!
!  9. STRUCTURE
!
! 10. SOURCE TEXT
!
!JL      REAL     COMPDA(MCGRD,MCMVAR)
!
!JL      REAL     AC2(MDC,MSC,MCGRD)
!
      INTEGER  KGRPNT(MXC,MYC), XYTST(*)                                 !30.70
!
      LOGICAL  INTERN

! SWAN in FVCOM local Vars :
      REAL(rkind) :: ART_TMP 
!      REAL(skind), ALLOCATABLE :: ART_TMP(:)
!# if defined (MULTIPROCESSOR)
      REAL(rkind), ALLOCATABLE :: AC2_TMP(:)
      INTEGER :: ISC,IDC,IP
!# endif
      REAL(rkind) :: SPPARM_TMP(4)
! SWAN in FVCOM local Vars end 

      SAVE IENT
      DATA IENT/0/

      CALL STRACE(IENT,'SWINCO')

!     *** Fetch Computation 'the mean delta' ***

!JL Skip and Replace with Method below
# if 0
      IF (KSPHER.EQ.0) THEN
        TLEN  = (XCLEN + YCLEN)/2.
      ELSE
        COSYG = COS (DEGRAD * (YOFFS + 0.5*(YCGMIN+YCGMAX)))             !33.09
        TLEN  = LENDEG * (COSYG*XCLEN + YCLEN)/2.
      ENDIF
      TDXY  = FLOAT(MXCGL + MYCGL)/2.                                    !40.31
      IF ( nvertsg.NE.0 ) TDXY = FLOAT(nvertsg)                          !40.95 40.80
!
      FETCH = TLEN/TDXY
#endif


!JL, SWAN in SCHISM begins

     TLEN = 0.0
     TDXY = 0.0

     IF(SERIAL) THEN 
      ! AREA = Sum over the global grid
      DO IN = 1, ne_global
       TLEN = TLEN + area(IN) 
      END DO

     ENDIF

! # if defined (MULTIPROCESSOR)
     IF(.NOT.SERIAL) THEN
 
      ART_TMP = 0.0
      DO IN = 1, ne
       ART_TMP =  ART_TMP + area(IN)
      END DO

!      print*,'ART_TMP=',ART_TMP  
   
      call mpi_allreduce(ART_TMP,TLEN,1,rtype,MPI_SUM,comm,ierr)

!      if(myrank==0) print*,'TLEN=',TLEN

     END IF 

     TLEN = SQRT(TLEN)
     TDXY = SQRT(REAL(ne_global))
     FETCH = TLEN/TDXY

!# if defined (MULTIPROCESSOR)
     IF(.NOT.SERIAL) THEN
      CALL mpi_bcast(FETCH,1,rtype,0,comm,istat)
      CALL mpi_bcast(TLEN,1,rtype,0,comm,istat)
     END IF

      SPPARM_TMP(1:4) = SPPARM(1:4)
!JL, SWAN in SCHISM End

      SY0   = 3.3                                                        !30.70
      IF (ITEST.GE.60 .OR. NPTST.GT.0) WRITE (PRTEST, 62) FETCH          !30.70
  62  FORMAT (' test SWINCO, fetch:', E12.4)                             !30.70


!     --- unstructured grid
!
      DO INX = 1, nverts                                                 !40.80
!        internal vertices and ghost vertices only

! JL Add on
         AC2(:,:,INX) = 0.0_rkind

         IF ( vmark(INX) == 0 .or. vmark(INX) == 999 ) THEN              !40.95
            TESTFL = .FALSE.

            DO IPTST = 1, NPTST

!             IF(ipgl(INX)%rank.EQ.myrank) THEN
!               IF(ipgl(INX)%id.EQ.XYTST(IPTST)) TESTFL = .TRUE.
               !IF(iplg(INX).EQ.XYTST(IPTST)) TESTFL = .TRUE.
               IF(INX.EQ.XYTST(IPTST)) TESTFL = .TRUE.
             !ENDIF 

            ENDDO
!
            IF (VARWI) THEN
              WX = COMPDA(INX,JWX2)
              WY = COMPDA(INX,JWY2)
!
!             *** Local wind speed and direction ***
              WSLOC = SQRT(WX*WX + WY*WY)
              IF (WX .NE. 0. .OR. WY .NE. 0.) THEN
                WDLOC = ATAN2(WY,WX)
              ELSE
                WDLOC = 0.
              ENDIF
            ELSE
!             uniform wind field
              WSLOC = U10
              WDLOC = WDIP
            ENDIF
!
            IF (WSLOC .GT. 1.E-10) THEN

! Dimensionless Hs and Tp calculated according to K.K. Kahma & C.J. Calkoen,
! (JPO, 1992) and Pierson-Moskowitz for limit values.
!
!             calculate dimensionless fetch:
              FDLSS = GRAV_W * FETCH / (WSLOC*WSLOC)
!
!             calculate dimensionless significant wave height:
              HSDLSS = MIN (0.21_rkind, 0.00288_rkind*FDLSS**0.45_rkind)
              SPPARM(1) = HSDLSS * WSLOC**2 / GRAV_W
!             calculate dimensionless peak period:
              TPDLSS = MIN (1._rkind/0.13_rkind, 0.46_rkind*FDLSS**0.27_rkind)
              SPPARM(2) = WSLOC * TPDLSS / GRAV_W
              IF (SPPARM(1).LT.0.05) SPPARM(2) = 2._rkind
              SPPARM(3) = 180._rkind * WDLOC / PI_W
              SPPARM(4) = 2._rkind
              IF (TESTFL) WRITE (PRTEST, 65) xcugrd(INX), &
     &        ycugrd(INX), FDLSS, (SPPARM(JJ), JJ = 1, 3)
            ELSE
              SPPARM(1) = 0.02_rkind
              SPPARM(2) = 2._rkind
              SPPARM(3) = 0._rkind
              SPPARM(4) = 0._rkind
            ENDIF

            CALL SSHAPE (AC2(1,1,INX), SPCSIG, SPCDIR, 2, 2)
!
         ENDIF

         IF(TESTFL) print*,'SWINCO',myrank,iplg(INX),' AC2 sum=',SUM(AC2(:,:,INX))

      ENDDO

! Update at ghost nodes :
#if 0
!# if defined (MULTIPROCESSOR)
!      IF(.NOT.SERIAL) CALL EXCHANGE_P4D_WWM(AC2)
      IF(.NOT.SERIAL)THEN
        ALLOCATE(AC2_TMP(0:npa));   AC2_TMP = 0._rkind
        DO ISC = 1,MSC
          DO IDC = 1,MDC
           DO IP = 1,npa ! np ?
              AC2_TMP(IP)=AC2(IDC,ISC,IP)
           END DO
           CALL exchange_p2d(AC2_TMP)
           DO IP = 1,npa
              AC2(IDC,ISC,IP) = AC2_TMP(IP)
           END DO
          END DO
        END DO

        DEALLOCATE(AC2_TMP)

        call parallel_barrier  !synchronize
      END IF
#endif

      IF(.NOT.SERIAL) THEN
       CALL EXCHANGE_P4D_WWM(AC2)
      ENDIF

!DBG
      DO INX = 1, nverts                                                 !40.80

        DO IPTST = 1, NPTST
         IF(INX.EQ.XYTST(IPTST)) THEN
            print*,'SWINCO aft',myrank,iplg(INX),' AC2 sum=',SUM(AC2(:,:,INX))
         ENDIF
        ENDDO
      ENDDO
      


!      IF(PAR)THEN
!        ALLOCATE(AC2_TMP(0:MT));   AC2_TMP = 0.0_skind
!        DO ISC = 1,MSC
!          DO IDC = 1,MDC
!           DO IP = 1,MT
!              AC2_TMP(IP)=AC2(IDC,ISC,IP)
!           END DO
!           CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,AC2_TMP)
!           CALL AEXCHANGE(NC,MYID,NPROCS,AC2_TMP)
!           DO IP = 1,MT
!              AC2(IDC,ISC,IP) = AC2_TMP(IP)
!           END DO
!          END DO
!        END DO
!        DEALLOCATE(AC2_TMP)
!      END IF
!JQI  IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!# endif

      SPPARM(1:4) = SPPARM_TMP(1:4)


      !PRINT *,'SWINCO max AC2 Cell 1',maxval(AC2(:,:,1))
      !PRINT *,'SWINCO max AC2 Cell 78',maxval(AC2(:,:,78))
      !PRINT *,'SWINCO max AC2 Cell 79',maxval(AC2(:,:,79))


  65  FORMAT (' test point ',6(1X,E12.4))                        !30.70

      RETURN
! * end of subroutine SWINCO *
      END
!****************************************************************
!
      SUBROUTINE SWCLME
!
!****************************************************************
!
      USE M_WCAP
      USE M_SNL4
      USE M_GENARR
      USE M_PARALL
      USE M_DIFFR
      USE SwanGriddata
      USE SwanCompdata

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 1993-2016  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!  0. Authors
!
!     40.31: Marcel Zijlema
!     40.80: Marcel Zijlema
!
!  1. Updates
!
!     40.31, Oct. 03: New subroutine
!     40.80, Aug. 07: extension to unstructured grids
!
!  2. Purpose
!
!     Clean memory
!
!  3. Method
!
!     De-allocates several allocatable arrays
!
!  9. Subroutines calling
!
!     SWMAIN
!
! 13. Source text
!
      IF (ALLOCATED(SIGPOW))   DEALLOCATE(SIGPOW)
      IF (ALLOCATED(CNL4_1))   DEALLOCATE(CNL4_1)
      IF (ALLOCATED(CNL4_2))   DEALLOCATE(CNL4_2)
      IF (ALLOCATED(LAMBDA))   DEALLOCATE(LAMBDA)
      IF (ALLOCATED(KGRPNT))   DEALLOCATE(KGRPNT)
      IF (ALLOCATED(KGRBND))   DEALLOCATE(KGRBND)
      IF (ALLOCATED(XYTST ))   DEALLOCATE(XYTST )
!      IF (ALLOCATED(AC2   ))   DEALLOCATE(AC2   )
      IF (ALLOCATED(XCGRID))   DEALLOCATE(XCGRID)
      IF (ALLOCATED(YCGRID))   DEALLOCATE(YCGRID)
      IF (ALLOCATED(SPCSIG))   DEALLOCATE(SPCSIG)
      IF (ALLOCATED(SPCDIR))   DEALLOCATE(SPCDIR)
      IF (ALLOCATED(DEPTH ))   DEALLOCATE(DEPTH )
      IF (ALLOCATED(FRIC  ))   DEALLOCATE(FRIC  )
      IF (ALLOCATED(UXB   ))   DEALLOCATE(UXB   )
      IF (ALLOCATED(UYB   ))   DEALLOCATE(UYB   )
      IF (ALLOCATED(WXI   ))   DEALLOCATE(WXI   )
      IF (ALLOCATED(WYI   ))   DEALLOCATE(WYI   )
      IF (ALLOCATED(WLEVL ))   DEALLOCATE(WLEVL )
      IF (ALLOCATED(ASTDF ))   DEALLOCATE(ASTDF )
      IF (ALLOCATED(IBLKAD))   DEALLOCATE(IBLKAD)
      IF (ALLOCATED(XGRDGL))   DEALLOCATE(XGRDGL)
      IF (ALLOCATED(YGRDGL))   DEALLOCATE(YGRDGL)
      IF (ALLOCATED(KGRPGL))   DEALLOCATE(KGRPGL)
      IF (ALLOCATED(KGRBGL))   DEALLOCATE(KGRBGL)
      IF (ALLOCATED(DIFPARAM)) DEALLOCATE(DIFPARAM)
      IF (ALLOCATED(DIFPARDX)) DEALLOCATE(DIFPARDX)
      IF (ALLOCATED(DIFPARDY)) DEALLOCATE(DIFPARDY)
      IF (ALLOCATED(MUDLF ))   DEALLOCATE(MUDLF )
      IF (ALLOCATED(NPLAF ))   DEALLOCATE(NPLAF )
      IF (ALLOCATED(LAYH  ))   DEALLOCATE(LAYH  )
      IF (ALLOCATED(VEGDIL))   DEALLOCATE(VEGDIL)
      IF (ALLOCATED(VEGNSL))   DEALLOCATE(VEGNSL)
      IF (ALLOCATED(VEGDRL))   DEALLOCATE(VEGDRL)
      IF (ALLOCATED(TURBF ))   DEALLOCATE(TURBF )
!
      IF (ALLOCATED(xcugrd  )) DEALLOCATE(xcugrd  )
      IF (ALLOCATED(ycugrd  )) DEALLOCATE(ycugrd  )
!PUN      IF (ALLOCATED(xcugrdgl)) DEALLOCATE(xcugrdgl)
!PUN      IF (ALLOCATED(ycugrdgl)) DEALLOCATE(ycugrdgl)
      IF (ALLOCATED(ivertg  )) DEALLOCATE(ivertg  )
      IF (ALLOCATED( vmark  )) DEALLOCATE( vmark  )
      IF (ALLOCATED( vlist  )) DEALLOCATE( vlist  )
      IF (ALLOCATED( blist  )) DEALLOCATE( blist  )
!
      RETURN
      END
#endif
