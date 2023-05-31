#include "swan_functions.h"

MODULE MOD_WAVE_CURRENT

   USE schism_glbl,ONLY: skind,rkind
   USE schism_msgp
   USE VARS_WAVE
   
   IMPLICIT NONE

   CONTAINS

!====================================================================================|
   SUBROUTINE WAVE_CURRENT_SETUP

   USE schism_glbl,ONLY: rkind,skind,errmsg
   USE schism_glbl,ONLY: nvrt,msc2,mdc2,nsa,npa,nea
   USE schism_glbl,ONLY: fwvor_advxy_stokes,fwvor_advz_stokes
   USE schism_glbl,ONLY: fwvor_gradpress,fwvor_breaking
   USE schism_glbl,ONLY: fwvor_streaming

   IMPLICIT NONE
   integer :: istat

   ALLOCATE(WAVESTRX_2D(npa));                        WAVESTRX_2D = 0.0_rkind
   ALLOCATE(WAVESTRY_2D(npa));                        WAVESTRY_2D = 0.0_rkind
   ALLOCATE(WAVESTRX_3D(nvrt,npa));                   WAVESTRX_3D = 0.0_rkind
   ALLOCATE(WAVESTRY_3D(nvrt,npa));                   WAVESTRY_3D = 0.0_rkind

   ALLOCATE(KN(npa));                                 KN = 0.0_rkind

   ! Computed Fields using schemes below
   ALLOCATE(jpress(npa));                             jpress = 0.0_rkind
   !ALLOCATE(wwave_force(2,nvrt,nsa)); wwave_force = 0.0_rkind !already done in schism_init
   ALLOCATE(stokes_hvel(2,nvrt,npa));                 stokes_hvel = 0.0_rkind
   ALLOCATE(stokes_wvel(nvrt,npa));                   stokes_wvel = 0.0_rkind
   ALLOCATE(stokes_hvel_side(2,nvrt,nsa));            stokes_hvel_side = 0.0_rkind
   ALLOCATE(stokes_wvel_side(nvrt,nsa));              stokes_wvel_side = 0.0_rkind

   ALLOCATE(roller_stokes_hvel_side(2,nvrt,nsa));     roller_stokes_hvel_side = 0.0_rkind
   ALLOCATE(roller_stokes_hvel(2,nvrt,npa));          roller_stokes_hvel = 0.0_rkind

   ALLOCATE(eps_w(npa),eps_r(npa),eps_br(npa),delta_wbl(npa),taub_wc(npa),stat=istat)
   if(istat/=0) call parallel_abort('WAVE SETUP: allocation failed 1')
   eps_w=0.d0; eps_r=0.d0; eps_br=0.d0; delta_wbl=0.d0; taub_wc=0.d0

   ALLOCATE(wave_sbrtot(npa),wave_sbftot(npa),wave_sdstot(npa),wave_sintot(npa), stat=istat)
   if(istat/=0) call parallel_abort('WAVE SETUP: allocation failed 2')
   wave_sbrtot=0.0D0; wave_sbftot=0.0D0; wave_sintot=0.0D0; wave_sdstot=0.0D0

   IF (RADFLAG == 'VOR') THEN
     if(myrank==0) write(16,*)'Schism-Swan coupling with VF'
     ! BM:
     ! Flags for accounting (1) or not (0) for the different
     ! terms involved in the vortex force formalism (RADFLAG='VOR')
     fwvor_advxy_stokes = 1
     fwvor_advz_stokes  = 1
     fwvor_gradpress    = 1
     fwvor_breaking     = 1
     fwvor_streaming    = 1

     ! JEROME : ROLLER are not implemented in schim-swan coupling yet 
     IF (IROLLER == 1) THEN

       ALLOCATE(EROL1(npa), EROL2(npa), stat=istat)
       IF (istat/=0) call parallel_abort('WAVE_CURRENT_SETUP, allocate error 2')
       EROL1 =  0.0_rkind; EROL2 =  0.0_rkind;

       ALLOCATE(DROLP(npa),CROLP(npa), stat=istat)
       IF (istat/=0) call parallel_abort('WAVE_CURRENT_SETUP, allocate error 3')
       DROLP =  0.0_rkind; CROLP =  0.0_rkind;

     END IF

   ELSEIF (RADFLAG == 'WON') THEN
     if(myrank==0) write(16,*)'Schism-Swan coupling with LH from Adc-Swan'
   ELSEIF (RADFLAG == 'LON') THEN
     if(myrank==0) write(16,*)'Schism-Swan coupling with LH from schism-wwm'
   ELSE 
     write(errmsg,*)'WAVE SETUP: impossible ',TRIM(RADFLAG)
     call parallel_abort(errmsg)
   END IF

!   == SWAN out fields (see swanser->swanout for details) ==

   ALLOCATE(out_wwm(npa,35));                  out_wwm = 0.0_rkind
   ALLOCATE(out_wwm_windpar(npa,10));          out_wwm_windpar = 0.0_rkind
   ALLOCATE(out_wwm_rol(npa,35));              out_wwm_rol     = 0.0_rkind


   ALLOCATE(SPEC_DENSITY(npa,msc2));           SPEC_DENSITY = 0.0_rkind
   ! SWAN out fields from SSURF and SBOT
   ALLOCATE(sbr(2,npa));                       sbr  = 0.0_rkind
   ALLOCATE(sbf(2,npa));                       sbf  = 0.0_rkind
   ALLOCATE(sds(2,npa));                       sds  = 0.0_rkind
   ALLOCATE(srol(2,npa));                      srol = 0.0_rkind

   !BM/JL: ramp on wwave_force at open boundary
   ALLOCATE( wafo_opbnd_ramp(nsa) );           wafo_opbnd_ramp=1.0d0
   !BM/JL: coupling current for WWM
   ALLOCATE( curx_wwm(npa) );                  curx_wwm = 0.d0
   ALLOCATE( cury_wwm(npa) );                  cury_wwm = 0.d0

   ALLOCATE( WK(msc2,npa) )     ;              WK = 0.0_rkind
   ALLOCATE( DS_BAND(0:msc2+1) );              DS_BAND = 0._rkind
   ALLOCATE( DS_INCR(0:msc2+1) );              DS_INCR = 0._rkind

   ALLOCATE( SXX3D(NVRT,npa) );                SXX3D = 0._rkind
   ALLOCATE( SXY3D(NVRT,npa) );                SXY3D = 0._rkind
   ALLOCATE( SYY3D(NVRT,npa) );                SYY3D = 0._rkind

   ALLOCATE( RSXX(npa) );                      RSXX = 0._rkind
   ALLOCATE( RSXY(npa) );                      RSXY = 0._rkind
   ALLOCATE( RSYY(npa) );                      RSYY = 0._rkind

   RETURN
   END SUBROUTINE WAVE_CURRENT_SETUP 

!**********************************************************************
!*                                                                    *
!**********************************************************************
!* This routine reads a map of bottom roughness length 'KN', used to 
!  derive the bottom friction based on the Madsen et al. 1988 scheme.
!  KN is equivalent of Nikuradse sand grain roughness. See: O.S. Madsen, 
!  Y.K. Poon, H.C. Graber (1988). Spectral wave attenuation by bottom
!  friction: Theory.” Proc. 21st Int. Conf. Coastal Engineering, ASCE, 492-504.
!  Authors: J.Lefevre (05/2022)
!**********************************************************************
   SUBROUTINE READ_KN_MAP

   USE schism_glbl,ONLY: rkind,in_dir,len_in_dir,errmsg
   USE schism_glbl,ONLY: ne_global,np_global,ipgl
   USE VARS_WAVE,ONLY: kn
   USE schism_msgp

   IMPLICIT NONE
   logical :: lexist
   integer :: itmp1,itmp2,i,j,istat
   real(rkind) :: xtmp,ytmp,tmp 
   real(rkind),allocatable :: kn_tmp(:)

    ! Experimental : load Roughness Length 'KN' used later in SWAN
    ! recycle unit 32 (use for drag,albedo, etc ..)
    if(not(allocated(kn_tmp))) allocate(kn_tmp(np_global),stat=istat)
    if(istat/=0) call parallel_abort('READ_KN_MAP: allocation failed')

    kn_tmp = 0.05  ! default value for kn (Madsen 1988)
 
    inquire(file=in_dir(1:len_in_dir)//'kn.gr3',exist=lexist)

    if(lexist) then
      if(myrank==0) then
        open(32,file=in_dir(1:len_in_dir)//'kn.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
           &call parallel_abort('Check kn.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0) then
            write(errmsg,*)'Negative KN roughness length',tmp
            call parallel_abort(errmsg)
          endif
          kn_tmp(i) = tmp
          !if(ipgl(i)%rank==myrank) kn(ipgl(i)%id)=tmp
        enddo
        close(32)
        write(16,*)'done reading kn.gr3'
      endif
    else ! no kn.gr3 file, use default value
      if(myrank==0) write(16,*)'done initializing kn with uniform value (0.05)'
    endif !lexist

    call mpi_bcast(kn_tmp,np_global,rtype,0,comm,istat)
    
    do i=1,np_global
        if(ipgl(i)%rank==myrank) kn(ipgl(i)%id)=kn_tmp(i)
    enddo !i

    deallocate(kn_tmp)

    END SUBROUTINE READ_KN_MAP


!**********************************************************************
!*                                                                    *
!**********************************************************************
!* This routine reads a ramp (between 0 and 1) for wave forces starting
!  from the open boundary, as defined in wafo_ramp.gr3.
!  Activation: wafo_obcramp (0/1:off/on).
!  Subroutine called at the initialization only (in INITIALIZE_WWM).
!  Then, wafo_opbnd_ramp is applied to wave forces at each time step in
!  wwm_main (routine inside wwm_coupl_selfe.F90).
!  Authors: X. Bertin & B. Mengual (05/2020)
!**********************************************************************
    SUBROUTINE READ_WAFO_OPBND_RAMP
!Since the output from this routine is only used in the coupler
!(wwm_coupl_selfe.F90), geometric arrays from SCHISM are used to account
!for quads

    USE VARS_WAVE, ONLY: wafo_opbnd_ramp
    USE schism_glbl,ONLY: rkind,in_dir,len_in_dir,errmsg
    USE schism_glbl, ONLY : xcj,ycj,xnd,ynd,ns,nsa,np_global,isidenode,ipgl
    USE schism_msgp

    IMPLICIT NONE
    LOGICAL      :: lexist
    INTEGER      :: IS,n1,n2,icount,IB,IP,j,itmp,istat
    REAL(rkind)  :: dx,dy,xtmp,ytmp,tmp1
    REAL(rkind),allocatable :: rwafop_tmp(:),ramp_wafop(:)

    ! Init: no ramp
    wafo_opbnd_ramp=1.0d0 ! at sides

    if(not(allocated(rwafop_tmp))) allocate(rwafop_tmp(np_global),stat=istat)
    if(istat/=0) call parallel_abort('READ_WAFO: allocation failed 1')
    if(not(allocated(ramp_wafop))) allocate(ramp_wafop(np_global),stat=istat)
    if(istat/=0) call parallel_abort('READ_WAFO: allocation failed 2')

    rwafop_tmp=0.0d0
    ramp_wafop=0.0d0

    ! Open wafo_ramp.gr3 and read ramp at nodes
    inquire(file=in_dir(1:len_in_dir)//'wafo_ramp.gr3',exist=lexist)

    if(lexist) then
      if(myrank==0) then
        open(10,file=in_dir(1:len_in_dir)//'wafo_ramp.gr3',status='old')
        read(10,*)
        read(10,*) j,itmp
        if(itmp/=np_global) CALL parallel_abort('READ_WAFO: Check np_global in wafo_ramp.gr3')
        DO IP = 1,np_global
          read(10,*)itmp,xtmp,ytmp,tmp1
          if(tmp1<0.0d0.or.tmp1>1.0d0) then
            write(errmsg,*)'READ_WAFO: ramp_wafo <0 or >1!',tmp1
            call parallel_abort(errmsg)
          endif
          rwafop_tmp(IP) = tmp1
          !if(ipgl(i)%rank==myrank) kn(ipgl(i)%id)=tmp
        enddo
        close(10)
        write(16,*)'done reading wafo_ramp.gr3'
      endif
    else ! no kn.gr3 file, use default value
      if(myrank==0) CALL parallel_abort('READ_WAFO: wafo_obcramp=1 requires wafo_ramp.gr3')
    endif !lexist

    call mpi_bcast(rwafop_tmp,np_global,rtype,0,comm,istat)

    DO IP=1,np_global
       if(ipgl(IP)%rank==myrank) ramp_wafop(ipgl(IP)%id)=rwafop_tmp(IP)
    ENDDO !i

    deallocate(rwafop_tmp)

    ! Compute ramp at sides
    DO IS = 1,nsa
       n1 = isidenode(1,IS); n2 = isidenode(2,IS)
       wafo_opbnd_ramp(IS)=(ramp_wafop(n1)+ramp_wafop(n2))/2.d0
    ENDDO

    deallocate(ramp_wafop)

    END SUBROUTINE READ_WAFO_OPBND_RAMP


! ==============================================================
   SUBROUTINE INIT_SPECTRAL_GRID
   !USE VARS_WAVE
   USE SWCOMM3
   USE M_GENARR, only: SPCDIR,SPCSIG
   use schism_glbl, only: rkind
   USE schism_glbl, ONLY: msc2,mdc2,npa,dp

   IMPLICIT NONE
   INTEGER :: IS,IP
   REAL(rkind) :: DEPLOC,KWAVELOC,SPSIGLOC,CGLOC,NN,ND

    IF (msc2 .GE. 2) THEN
        DS_BAND(0)     = SPCSIG(2)- SPCSIG(1)
        DS_BAND(1)     = DS_BAND(0)
        DS_BAND(msc2)   = SPCSIG(msc2) - SPCSIG(msc2-1)
        DS_BAND(msc2+1) = DS_BAND(msc2)
        DS_INCR(0)     = DS_BAND(0)
        DS_INCR(1)     = DS_BAND(0)
        DS_INCR(msc2)   = DS_BAND(msc2)
        DS_INCR(msc2+1) = DS_INCR(msc2)
        DO IS = 2, msc2-1 ! Bandwith at gridpoints
          DS_BAND(IS) = (SPCSIG(IS)-SPCSIG(IS-1))/2._rkind + &
                        (SPCSIG(IS+1)-SPCSIG(IS))/2._rkind
        END DO
        DO IS = 2, msc2 ! Stepwidth between gridpoints K and K-1
          DS_INCR(IS) = SPCSIG(IS) - SPCSIG(IS-1)
        END DO
   END IF

! Compute Kwave

   DO IP = 1, npa
      DEPLOC = MAX(DEPMIN,dp(IP))
!           DEPLOC = DEP(IP)
!           WRITE(740+myrank,*) 'IP=', IP, ' DEPLOC=', DEPLOC
     DO IS = 1, msc2
      SPSIGLOC = SPCSIG(IS)
      CALL KSCIP1(1,SPSIGLOC,DEPLOC,KWAVELOC,CGLOC,NN,ND)
              !KWAVELOC: wave number K, CGLOC: group velocity,
              !NN:ratio of group and phase velocity and its derivative ND
      WK(IS,IP) = KWAVELOC
!     CG(IS,IP) = CGLOC
!     WC(IP,IS) = CGLOC/NN !WVC
     END DO
   END DO


   RETURN
   END SUBROUTINE INIT_SPECTRAL_GRID


#ifdef USE_SWAN
!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after
!Bennis et al., 2011)
!*  => Computation of the wave-induced pressure term at nodes (the
!gradient is computed directly
!*  at sides when calculating the forces) and the Stokes drift
!velocities. The latter are
!*  computed at all levels, at nodes and sides, and for both the wave
!and roller (kept separated).

!  Hacked from schism v5.9, LR : ToDo : Fix roller case, rename stokes_vel ...

!**********************************************************************
   SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM

        USE VARS_WAVE
        USE SWCOMM3, only: DMIN=>DEPMIN,DDIR,G9=>GRAV_W !=(SPDIR2-SPDIR1)/MDC, mesh size in theta-direction of computational grid
        USE M_GENARR, only: SPCDIR,SPCSIG
!     SPCDIR: (*,1); spectral directions (radians)
!             (*,2); cosine of spectral directions
!             (*,3); sine of spectral directions
!             (*,4); cosine^2 of spectral directions
!             (*,5); cosine*sine of spectral directions
!             (*,6); sine^2 of spectral directions
!     SPCSIG: Relative frequencies in computational domain in sigma-space
        USE schism_glbl, ONLY: errmsg, hmin_radstress, ns, kbs, kbe, nea, idry_e, &
                            &  isdel, indel, elnode, dldxy, zs, area,nsa,idry_s, &
                            &  isidenode,nne, idry, eta2, kbp, np, nvrt, msc2,mdc2
        use schism_glbl, only: MNP => npa,       & ! Nodes of the augmented domain
     &                         MNE => nea,       & ! Elements of the augmented domain
     &                         DEP8 => dp,       & ! depth in the augmented domain
     &                         ipgl,             & ! node global to local mapping
     &                         ielg,             & ! element local to global maping
     &                         ZETA => znl       ! Z-Levels of SCHISM

        USE schism_msgp
        IMPLICIT NONE

        INTEGER     :: IP, k, ID, IS, IL
        REAL(rkind) :: DEP_loc, D_loc, k_loc, kD_loc, z_loc, E_loc, Er_loc, JPress_loc, eMult
        REAL(rkind) :: Uint, Vint, Urint, Vrint
        REAL(rkind) :: USTOKES_loc(NVRT), VSTOKES_loc(NVRT),UrSTOKES_loc(NVRT), VrSTOKES_loc(NVRT)

        integer     :: IE, isd, j, l, n1, n2, n3, icount
        real(rkind) :: tmp0, tmp1, tmp2, ztmp, ubar, vbar, dhdx, dhdy
        real(rkind) :: STOKES_WVEL_ELEM(NVRT,MNE), ws_tmp1(NVRT,nsa),ws_tmp2(NVRT,nsa)

        !REAL(rkind) :: ACLOC(mdc2,msc2)

        real(rkind) :: dr_dxy_loc(2,NVRT,nsa)

        real(rkind) :: ustk2d(2),cff

!...    Computing Stokes drift horizontal velocities at nodes and
!       pressure term
        DO IP = 1, MNP
          IF(idry(IP) == 1) CYCLE

          ! Total water depth at the node
          DEP_loc = DEP8(IP)+ETA2(IP)
          D_loc = MAX( DEP_loc, hmin_radstress )

          ! Initialization of the local Stokes drift and J variables
          USTOKES_loc  = 0.D0; VSTOKES_loc  = 0.D0;
          UrSTOKES_loc = 0.D0; VrSTOKES_loc = 0.D0;
          JPress_loc   = 0.D0;

          ! Loop on the frequencies
          DO IS=1,msc2

            eMult = SPCSIG(IS)*DDIR*DS_INCR(IS)
            !CALL KSCIP1(1,SPCSIG(IS),DEP8(IP),KWAVELOC,CGLOC,NN,ND)
            !eWk=WK(IS,IP)
            !eWk=KWAVELOC
            k_loc  = MIN(KDMAX/DEP_loc,WK(IS,IP))
            kD_loc = MIN(KDMAX, WK(IS,IP)*D_loc)
            !eWkReal = kD/eDep
            !eSinh2kd= MySINH(2*kD)
            !eSinhkd = MySINH(kD)
            !eSinhkd2= eSinhkd**2
            !eSigma  = SPCSIG(IS)
            Uint   = 0.D0
            Vint   = 0.D0
            Urint  = 0.D0
            Vrint  = 0.D0
            ! Loop on the directions
            DO ID=1,mdc2
              E_loc  = AC2(ID,IS,IP)*eMult
              JPress_loc = JPress_loc + G9*k_loc*E_loc/DSINH(2.D0*kD_loc)
              Uint = Uint + SPCSIG(IS)*k_loc*SPCDIR(ID,2)*E_loc       !COSTH(ID)
              Vint = Vint + SPCSIG(IS)*k_loc*SPCDIR(ID,3)*E_loc       !SINTH(ID)
!              IF (IROLLER == 1) THEN
!                Er_loc   = RAC2(IS,ID,IP)*SPSIG(IS)*DDIR*DS_INCR(IS)
!                Urint    = Urint + SPSIG(IS)*k_loc*COSTH(ID)*Er_loc
!                Vrint    = Vrint + SPSIG(IS)*k_loc*SINTH(ID)*Er_loc
!              END IF
            END DO !ID

            ! Loop on the vertical nodes
             DO IL = KBP(IP), NVRT
              ! Here we need to compute z+h of Eq. C.1 of Bennis et al.
              ! (2011)
              ! In her framework, z varies from -h to eta, meaning that
              ! z+h corresponds to the distance to the bed
              ! -ZETA(KBP(IP),IP) corresponds to h, the depth at node IP
              ! (not the total water depth)
              ! Waves
              z_loc            = ZETA(IL,IP) - ZETA(KBP(IP),IP) !from bottom; 'z+D'
              USTOKES_loc(IL)  = USTOKES_loc(IL)  + Uint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              VSTOKES_loc(IL)  = VSTOKES_loc(IL)  + Vint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
              ! Surface rollers
!              IF (IROLLER == 1) THEN
!                UrSTOKES_loc(IL) = UrSTOKES_loc(IL) + Urint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
!                VrSTOKES_loc(IL) = VrSTOKES_loc(IL) + Vrint*DCOSH(2.D0*k_loc*z_loc)/DSINH(kD_loc)**2
!              END IF
             END DO !NVRT
          END DO !MSC

          ! Surface roller contribution to horizontal Stokes drift velocities
          ! NB: we do not just add the contribution and keep separated arrays.
          ! This is motivated by the fact that we do not want this contribution to
          ! influence Wst, which is computed from the continuity equation for waves only

          IF (IROLLER == 1) THEN
            Urint = 2.D0*COS(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)
            Vrint = 2.D0*SIN(DROLP(IP))*EROL2(IP)/(CROLP(IP)*D_loc)

            ! Homogeneous across depth
            UrSTOKES_loc = Urint
            VrSTOKES_loc = Vrint

            ! Making sure, the Stokes drift velocities do not blow up in very shallow water
            IF (D_loc < 2.D0*hmin_radstress) THEN
              UrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Urint)),Urint)
              VrSTOKES_loc = SIGN(MIN(0.1D0*SQRT(G9*D_loc),ABS(Vrint)),Vrint)
            END IF
          END IF

          ! Storing Stokes drift horizontal velocities 
          STOKES_HVEL(1,:,IP) = USTOKES_loc
          STOKES_HVEL(2,:,IP) = VSTOKES_loc
          ! Surface rollers
          IF (IROLLER == 1) THEN
            ! Smoothing the roller contribution to the Stokes drift velocity near the shoreline
            ! With this profile, U_st < 10% of computed U_st at h < DMIN, and U_st > 95% of computed U_st at h > 2.25*DMIN
            IF (D_loc < 1.5D0*DMIN) THEN
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc*(SINH(DEP_loc)/SINH(1.5D0))**2
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc*(SINH(DEP_loc)/SINH(1.5D0))**2
            ELSE
              ROLLER_STOKES_HVEL(1,:,IP) = UrSTOKES_loc
              ROLLER_STOKES_HVEL(2,:,IP) = VrSTOKES_loc
            END IF
          END IF

          ! Pressure term
          JPRESS(IP) = JPress_loc

          ! JEROME Replace WavPres from swanser with JPRESS in OUT_wwm(:,40)
          out_wwm(IP,40) = JPress_loc

! Jerome for Check out
#if 1
          ustk2d =  0.D0
          DO k=kbp(IP)+1,nvrt
           cff = (ZETA(k,IP) - ZETA(k-1,IP))/D_loc
           ustk2d(1) = ustk2d(1) + cff*USTOKES_loc(k)
           ustk2d(2) = ustk2d(2) + cff*VSTOKES_loc(k)
          ENDDO
          out_wwm(IP,41) = ustk2d(1)
          out_wwm(IP,42) = ustk2d(2)
#endif

        END DO !MNP

!...    Computing Stokes drift horizontal velocities at sides (in pframe if ics=2)
        ! The average of the values from vertically adjacent nodes is taken
        STOKES_HVEL_SIDE = 0.D0; ROLLER_STOKES_HVEL_SIDE = 0.D0
        DO IS = 1,nsa
          IF(idry_s(IS) == 1) CYCLE

          ! Indexes of surrounding nodes
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          DO k = kbs(IS),NVRT
            ! Waves
            STOKES_HVEL_SIDE(1,k,IS) = (STOKES_HVEL(1,k,n1) + STOKES_HVEL(1,k,n2))/2.D0
            STOKES_HVEL_SIDE(2,k,IS) = (STOKES_HVEL(2,k,n1) + STOKES_HVEL(2,k,n2))/2.D0
            ! Surface rollers
            IF (IROLLER == 1) THEN
              ROLLER_STOKES_HVEL_SIDE(1,k,IS) = (ROLLER_STOKES_HVEL(1,k,n1) + ROLLER_STOKES_HVEL(1,k,n2))/2.D0
              ROLLER_STOKES_HVEL_SIDE(2,k,IS) = (ROLLER_STOKES_HVEL(2,k,n1) + ROLLER_STOKES_HVEL(2,k,n2))/2.D0
            END IF
          END DO
        END DO !nsa

!...    Compute bottom Stokes drift z-vel. at elements
        STOKES_WVEL_ELEM = 0.D0
        DO IE = 1,nea
           IF(idry_e(IE) == 1) CYCLE

           ! Index of the surrounding nodes
           n1 = elnode(1,IE)
           n2 = elnode(2,IE)
           n3 = elnode(3,IE)
           IF(kbe(IE) == 0) THEN
             WRITE(errmsg,*)'Error: kbe(i) == 0'
             CALL parallel_abort(errmsg)
           END IF

           ubar = (STOKES_HVEL(1,max(kbp(n1),kbe(IE)),n1) + STOKES_HVEL(1,max(kbp(n2),kbe(IE)),n2) &
               & + STOKES_HVEL(1,max(kbp(n3),kbe(IE)),n3))/3.D0 !average bottom stokes-x-vel
           vbar = (STOKES_HVEL(2,max(kbp(n1),kbe(IE)),n1) + STOKES_HVEL(2,max(kbp(n2),kbe(IE)),n2) &
               & + STOKES_HVEL(2,max(kbp(n3),kbe(IE)),n3))/3.D0 !average bottom stokes-y-vel

           dhdx = DEP8(n1)*dldxy(1,1,IE) + DEP8(n2)*dldxy(2,1,IE) + DEP8(n3)*dldxy(3,1,IE) !eframe
           dhdy = DEP8(n1)*dldxy(1,2,IE) + DEP8(n2)*dldxy(2,2,IE) + DEP8(n3)*dldxy(3,2,IE)
           STOKES_WVEL_ELEM(kbe(IE),IE) = -dhdx*ubar - dhdy*vbar
        END DO !nea

!...    Compute bottom Stokes z-vel. at nodes
        STOKES_WVEL = 0.D0
        DO IP = 1,np !residents only
           IF(idry(IP) == 1) CYCLE

           !Bottom Stokes z-vel.
           tmp0 = 0.D0
           DO j = 1,nne(IP)
              ie = indel(j,IP)
              IF(idry_e(ie)==0) THEN
                STOKES_WVEL(kbp(IP),IP) = STOKES_WVEL(kbp(IP),IP) + STOKES_WVEL_ELEM(kbe(ie),ie)*area(ie)
              END IF
              tmp0 = tmp0 + area(ie)
           END DO !j
           STOKES_WVEL(kbp(IP),IP) = STOKES_WVEL(kbp(IP),IP)/tmp0
        END DO !np

!...    Compute horizontal gradient of Stokes x and y-vel. (to compute Stokes z-vel.)
        ws_tmp1 = 0.D0; ws_tmp2 = 0.D0
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,STOKES_HVEL(1,:,:),dr_dxy_loc)
        ws_tmp1(:,:) = dr_dxy_loc(1,:,:) !valid only in resident
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,STOKES_HVEL(2,:,:),dr_dxy_loc)
        ws_tmp2(:,:) = dr_dxy_loc(2,:,:)

!...    Compute Stokes z-vel. at all levels: STOKES_WVEL_SIDE(NVRT,nsa)
        STOKES_WVEL_SIDE = 0.D0
        DO IS = 1,ns !residents only
          IF(idry_s(IS) == 1) CYCLE
          n1 = isidenode(1,IS)
          n2 = isidenode(2,IS)

          !Bottom Stokes z-vel.
          STOKES_WVEL_SIDE(kbs(IS),IS) = (STOKES_WVEL(max(kbs(IS),kbp(n1)),n1) + STOKES_WVEL(max(kbs(IS),kbp(n2)),n2))/2.D0

          !Stokes z-vel. at all levels
          DO k = kbs(IS)+1, NVRT
            ztmp = zs(k,IS) - zs(k-1,IS)
            STOKES_WVEL_SIDE(k,IS) = STOKES_WVEL_SIDE(k-1,IS) &
                                   & -(ws_tmp1(k,IS)+ws_tmp1(k-1,IS))/2.D0*ztmp &
                                   & -(ws_tmp2(k,IS)+ws_tmp2(k-1,IS))/2.D0*ztmp
          END DO
        END DO !ns

      END SUBROUTINE STOKES_STRESS_INTEGRAL_SCHISM

!**********************************************************************
!  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!  => Computation of the conservative terms A1 and B1 from Eq. (11) and (12) respectively
!
! Bennis, A.-C., F. Ardhuin, and F. Dumas (2011), On the coupling of wave
! and three-dimensional circulation models: Choice of theoretical frame-
! work, practical implementation and adiabatic tests, Ocean Modell.,
! 40(3–4), 260–272, doi:10.1016/j.ocemod.2011.09.003.
!
! Guerin, Thomas & Bertin, X & Coulombier, Thibault & de Bakker, Anouk. (2018).
! Impacts of wave-induced circulation in the surf zone on wave setup.
! Ocean Modelling. 123. 10.1016/j.ocemod.2018.01.006.
!
! Nirnimesh Kumar, George Voulgaris, John C. Warner, Maitane Olabarrieta,
! Implementation of the vortex force formalism in the coupled
! ocean-atmosphere-wave-sediment transport (COAWST) modeling system
! for inner shelf and surf zone applications,
! Ocean Modelling, Volume 47, 2012, Pages 65-95
!**********************************************************************
!
!  References: See COMPUTE_CONSERVATIVE_VF_TERMS_SCHISM, WWMIII/wwm_coupl_selfe.F90
!
!**********************************************************************
      SUBROUTINE CONSERVATIVE_VF_TERMS_SCHISM
        USE VARS_WAVE
        USE schism_glbl, ONLY: kbs, ns, idry_e, isdel, elnode, dldxy, cori, zs, su2, sv2, nsa, idry_s
        USE schism_glbl, ONLY: uu2, vv2, kbp, nvrt, fwvor_advz_stokes,fwvor_advxy_stokes
        USE schism_glbl, ONLY: fwvor_gradpress
        use schism_glbl, only: MNP => npa        ! Nodes of the augmented domain
        use schism_glbl, only: dps,isidenode,eta2

        USE schism_msgp
!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        USE schism_glbl, only: DIAG_BHD_BAR,DIAG_SKC_BAR,DIAG_VOR_BAR
!#endif

        IMPLICIT NONE

        integer     :: IS, IE, k, l, icount
        real(rkind) :: dJ_dx_loc, dJ_dy_loc, du_loc, dv_loc, dz_loc, Ust_loc, Vst_loc
        real(rkind) :: du_dxy(2,NVRT,nsa), dv_dxy(2,NVRT,nsa)

        real(rkind) :: VorF_force(2,nvrt)  ! Bucket for Vortex Force
        real(rkind) :: StkC_force(2,nvrt)  ! Bucket for Stokes-Coriolis Force

!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        real(rkind),allocatable :: swild69(:,:) !used for exchange (deallocate immediately afterwards)
!        real(rkind) :: sum1,sum2,sum3,sum4
!        integer     :: istat

!        allocate(swild69(9,nsa),stat=istat)
!        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild69')
!#endif

!YJZ: init here for better readability
        WWAVE_FORCE = 0.d0
!#if defined USE_PLANARBEACH || defined USE_DUCK94
!!        diag_BOT_bar = 0.D0  ! Wave Bottom Dissipation
!        diag_BHD_bar = 0.D0  ! Bernouilli Head Pressure
!        diag_SKC_bar = 0.D0  ! Stokes-Coriolis
!        diag_VOR_bar = 0.D0  ! Vortex Force
!#endif

!...    Computing the spatial derivative of horizontal velocities
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,uu2,du_dxy)
        CALL hgrad_nodes(2,0,NVRT,MNP,nsa,vv2,dv_dxy)

!...    Main loop over the sides
        DO IS = 1,ns !resident
          IF(idry_s(IS) == 1) CYCLE

          !------------------------
          ! Pressure term (grad(J))
          icount = 0; dJ_dx_loc = 0.d0; dJ_dy_loc = 0.d0

          IF (fwvor_gradpress == 1) THEN ! BM

            DO l = 1,2 !elements
              IE = isdel(l,IS)
              IF(ie /= 0 .AND. idry_e(IE) == 0) THEN
                icount = icount + 1
                dJ_dx_loc = dJ_dx_loc + dot_product(JPRESS(elnode(1:3,IE)),dldxy(1:3,1,IE)) !in eframe
                dJ_dy_loc = dJ_dy_loc + dot_product(JPRESS(elnode(1:3,IE)),dldxy(1:3,2,IE))
              END IF
            END DO !l

            ! Averaging the values from the two surrounding elements
            IF(icount > 2) CALL parallel_abort('Pressure term:icount>2')
            IF(icount == 2) THEN
              dJ_dx_loc = dJ_dx_loc/2.D0
              dJ_dy_loc = dJ_dy_loc/2.D0
            END IF

          END IF

          !------------------------
          ! Saving wave forces: loop over the vertical
          ! Description of the terms:
          ! 1 - Terms with Coriolis force and the spatial derivative of
          ! horizontal velocities
          ! 2 - Term of the spatial variation de the wave-induced
          ! pressure (J)
          ! 3 - Term -w_s*(du/dz,dv/dz)
          du_loc = 0.D0; dv_loc = 0.D0; dz_loc = 1.D0

          ! init here for better 
          VorF_force = 0.d0
          StkC_force = 0.d0

          DO k = kbs(IS),NVRT

            IF (fwvor_advz_stokes == 1) THEN ! BM

              ! du/dz and dv/dz terms
              IF (k == kbs(IS)) THEN
                dz_loc = zs(k+1,IS) - zs(k,IS)
                du_loc = su2(k+1,IS) - su2(k,IS)
                dv_loc = sv2(k+1,IS) - sv2(k,IS)
              ELSE IF (k == NVRT) THEN
                dz_loc = zs(k,IS) - zs(k-1,IS)
                du_loc = su2(k,IS) - su2(k-1,IS)
                dv_loc = sv2(k,IS) - sv2(k-1,IS)
              ELSE
                dz_loc = zs(k+1,IS) - zs(k-1,IS)
                du_loc = su2(k+1,IS) - su2(k-1,IS)
                dv_loc = sv2(k+1,IS) - sv2(k-1,IS)
              END IF

            END IF

            ! Stokes drift velocity
            ! LRU team : switch off roller contribution, which is only accounted
            ! for within continuity equation. This is motivated by the fact that VF
            ! arises from the irrotational part of the wave motion as opposed
            ! to surface rollers.
            Ust_loc = 0.D0; Vst_loc = 0.D0
            IF (fwvor_advxy_stokes == 1) THEN
              !IF (IROLLER == 1) THEN
              !  Ust_loc = STOKES_HVEL_SIDE(1,k,IS) + ROLLER_STOKES_HVEL_SIDE(1,k,IS)
              !  Vst_loc = STOKES_HVEL_SIDE(2,k,IS) + ROLLER_STOKES_HVEL_SIDE(2,k,IS)
              !ELSE
              !  Ust_loc = STOKES_HVEL_SIDE(1,k,IS)
              !  Vst_loc = STOKES_HVEL_SIDE(2,k,IS)
              !END IF
              Ust_loc = STOKES_HVEL_SIDE(1,k,IS)
              Vst_loc = STOKES_HVEL_SIDE(2,k,IS)
            END IF

            ! -----------------------------------------------------------
            ! -- VORTEX FORCE 
            ! -----------------------------------------------------------

            ! -du/dy*v_s +dv/dx*v_s - W_s*du/dz
            ! ~~~~~~~~~~                                   -du/dy   *    v_s
            VorF_force(1,k) = VorF_force(1,k) - du_dxy(2,k,IS)*Vst_loc  
            ! +du/dy*u_s -dv/dx*u_s - W_s*dv/dz
            ! ~~~~~~~~~                                     du/dy   *    u_s
            VorF_force(2,k) = VorF_force(2,k) + du_dxy(2,k,IS)*Ust_loc 


            ! -du/dy*v_s +dv/dx*v_s - W_s*du/dz
            !            ~~~~~~~~~~                        dv/dx    *      v_s
            VorF_force(1,k) = VorF_force(1,k) + dv_dxy(1,k,IS)*Vst_loc
            ! +du/dy*u_s -dv/dx*u_s - W_s*dv/dz
            !            ~~~~~~~~~                        -dv/dx    *      u_s
            VorF_force(2,k) = VorF_force(2,k) - dv_dxy(1,k,IS)*Ust_loc


            ! -du/dy*v_s +dv/dx*v_s - W_s*du/dz
            !                         ~~~~~~~~~           - W_s        * du/dz
            VorF_force(1,k) = VorF_force(1,k) - STOKES_WVEL_SIDE(k,IS)*du_loc/dz_loc 
            ! +du/dy*u_s -dv/dx*u_s - W_s*dv/dz
            !                         ~~~~~~~~~           - W_s        * dv/dz
            VorF_force(2,k) = VorF_force(2,k) - STOKES_WVEL_SIDE(k,IS)*dv_loc/dz_loc                

            ! -----------------------------------------------------------
            ! -- STOKES-CORIOLIS FORCE 
            ! -----------------------------------------------------------

            ! +f*v_s
            StkC_force(1,k) = StkC_force(1,k) + (cori(IS)*Vst_loc)
            ! -f*u_s
            StkC_force(2,k) = StkC_force(2,k) - (cori(IS)*Ust_loc)

            ! Saving wave forces
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + StkC_force(1,k) + VorF_force(1,k) - dJ_dx_loc
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + StkC_force(2,k) + VorF_force(2,k) - dJ_dy_loc

          END DO  !k


!#if defined USE_PLANARBEACH || defined USE_DUCK94
!          sum1 = 0.d0 !integral; x-comp. for Vortex Force
!          sum2 = 0.d0 !integral; y-comp. for Vortex Force
!          sum3 = 0.d0 !integral; x-comp. for Stoke-Coriolis Force
!          sum4 = 0.d0 !integral; y-comp. for Stoke-Coriolis Force

!          DO k = kbs(IS)+1,NVRT 
!            ! --- Depth integ. the Vortex Force
!            sum1=sum1+(zs(k,IS)-zs(k-1,IS))*(VorF_force(1,k)+VorF_force(1,k-1))/2.d0
!            sum2=sum2+(zs(k,IS)-zs(k-1,IS))*(VorF_force(2,k)+VorF_force(2,k-1))/2.d0
!            ! --- Depth integ. the Stokes-Coriolis Force
!            sum3=sum3+(zs(k,IS)-zs(k-1,IS))*(StkC_force(1,k)+StkC_force(1,k-1))/2.d0
!            sum4=sum4+(zs(k,IS)-zs(k-1,IS))*(StkC_force(2,k)+StkC_force(2,k-1))/2.d0
!          ENDDO

!          ! --- Save the wave acceleration term due to the Vortex Force [m.s-2]
!          diag_VOR_bar(1,IS) = diag_VOR_bar(1,IS) + sum1/(dps(IS) + SUM(eta2(isidenode(:,IS)))/2.d0)
!          diag_VOR_bar(2,IS) = diag_VOR_bar(2,IS) + sum2/(dps(IS) + SUM(eta2(isidenode(:,IS)))/2.d0)
!          ! --- Save the wave acceleration term due to the Stokes-Coriolis Force [m.s-2]
!          diag_SKC_bar(1,IS) = diag_SKC_bar(1,IS) + sum3/(dps(IS) + SUM(eta2(isidenode(:,IS)))/2.d0)
!          diag_SKC_bar(2,IS) = diag_SKC_bar(2,IS) + sum4/(dps(IS) + SUM(eta2(isidenode(:,IS)))/2.d0)
!          ! --- Save the wave acceleration term due to the Bernouilli Head force
!          diag_BHD_bar(1,IS) = diag_BHD_bar(1,IS) -dJ_dx_loc
!          diag_BHD_bar(2,IS) = diag_BHD_bar(2,IS) -dJ_dy_loc
!#endif

        END DO !ns

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        swild69(1,:)=diag_VOR_bar(1,:); swild69(2,:)=diag_VOR_bar(2,:)
!        swild69(3,:)=diag_SKC_bar(1,:); swild69(4,:)=diag_SKC_bar(2,:)
!        swild69(5,:)=diag_BHD_bar(1,:); swild69(6,:)=diag_BHD_bar(2,:)
!        swild69(7:9,:) = 0.d0
!        CALL exchange_s2d_9(swild69)
!        diag_VOR_bar(1,:)=swild69(1,:); diag_VOR_bar(2,:)=swild69(2,:)
!        diag_SKC_bar(1,:)=swild69(3,:); diag_SKC_bar(2,:)=swild69(4,:)
!        diag_BHD_bar(1,:)=swild69(5,:); diag_BHD_bar(2,:)=swild69(6,:)
!        deallocate(swild69)
!#endif

      END SUBROUTINE CONSERVATIVE_VF_TERMS_SCHISM

!**********************************************************************
!*  This routine is used with RADFLAG=VOR (3D vortex formulation, after Bennis et al., 2011)
!*  => Computation of the non-conservative terms due to depth-induced breaking (term Fb from Eq. (11) and (12))
!**********************************************************************
      SUBROUTINE BREAKING_VF_TERMS_SCHISM
        USE VARS_WAVE  ! only: out_wwm
!        USE schism_glbl, only: errmsg,hmin_radstress,kbs,kbe,nsa,ns,np,ne,nea,idry_e, &
!                             & idry_s,isidenode,isbs,i34,dps,dldxy,dr_dxy,h0,zs,nvrt,eta2
        USE schism_glbl, ONLY: hmin_radstress, kbs, ns, isbs, dps, &
                          &    zs,nsa,idry_s,isidenode,eta2,nvrt
!        use schism_glbl, only: DMIN=>h0 ! Minimum water depth. THis must be same as h0 in SCHISM
        USE schism_glbl, ONLY: h0
        USE SWCOMM3, only: DEPMIN

!#JL: confusion between h0, et DMIN cf. wwm' datapool. I use DMIN=0.01 , as hardoded in wwm' datapool

        USE schism_msgp
!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        USE schism_glbl, only: diag_BRK_bar
!#endif
        IMPLICIT NONE

        INTEGER     :: IS, isd, k, j, l, n1, n2, n3, icount
        REAL(rkind) :: eta_tmp, tmp0, htot, sum_2D, sum_3D
        REAL(rkind) :: swild_2D(NVRT), swild_3D(NVRT)

        real(rkind) :: BBbrk_force(2,nvrt) ! Wave acceleration due to breaking
        real(rkind) :: BBrol_force(2,nvrt) ! Wave acceleration due to roller
        real(rkind) :: BBwhc_force(2,nvrt) ! acceleration due to whitecapping
!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        real(rkind) :: sum1, sum2
!#endif

        ! Apply lpp_filter : JL ToDo ??
!        IF (LPP_FILT_FLAG) CALL LPP_FILT(SBR(1,:))
!        IF (LPP_FILT_FLAG) CALL LPP_FILT(SBR(2,:))

        ! Init here
!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        diag_BRK_bar = 0.D0  ! Wave Breaking Acceleration
!#endif

        ! Compute sink of momentum due to wave breaking
        DO IS = 1, ns

          ! Check IF dry segment or open bnd segment
          IF(idry_s(IS) == 1 .or. isbs(IS) > 0) CYCLE

          ! Water depth at side
          n1 = isidenode(1,IS); n2 = isidenode(2,IS)
          eta_tmp = (eta2(n1) + eta2(n2))/2.D0
          !htot = max(h0,dps(IS)+eta_tmp,hmin_radstress) ! KM
          htot = max(h0,dps(IS)+eta_tmp)
          !htot = max(htot,hmin_radstress) ! je remets la version KM car seg fault
                                          ! ou bien 0.0 every et cycle cas htot==0.0 ?    
          ! init here for better
          BBbrk_force = 0.d0
          BBrol_force = 0.d0
          BBwhc_force = 0.d0

          IF(kbs(IS)+1 == NVRT) THEN !2D
            ! average between the two adjacent nodes

            ! Breaking acceleration + roller contribution
            IF (IROLLER == 1) THEN
              BBrol_force(1,:) = BBrol_force(1,:) -(SROL(1,n1) + SROL(1,n2))/2.D0/htot
              BBrol_force(2,:) = BBrol_force(2,:) -(SROL(2,n1) + SROL(2,n2))/2.D0/htot
              BBbrk_force(1,:) = BBbrk_force(1,:) -(1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2))/2.D0/htot 
              BBbrk_force(2,:) = BBbrk_force(2,:) -(1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2))/2.D0/htot
            ELSE
              BBbrk_force(1,:) = BBbrk_force(1,:) - (SBR(1,n1) + SBR(1,n2))/2.D0/htot
              BBbrk_force(2,:) = BBbrk_force(2,:) - (SBR(2,n1) + SBR(2,n2))/2.D0/htot
            ENDIF
            ! Whitecapping contribution ? SWAN, already in SBR ? 
            BBwhc_force(1,:) = BBwhc_force(1,:) -(SDS(1,n1) + SDS(1,n2))/2.D0/htot
            BBwhc_force(2,:) = BBwhc_force(2,:) -(SDS(2,n1) + SDS(2,n2))/2.D0/htot

            ! --- acceleration due to wave breaking and white capping(+wcap)
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) + BBbrk_force(1,:)
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) + BBbrk_force(2,:)

            ! --- acceleration due to wave roller
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) + BBrol_force(1,:)
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) + BBrol_force(2,:)

            ! --- acceleration due to whitecapping
            WWAVE_FORCE(1,:,IS) = WWAVE_FORCE(1,:,IS) + BBwhc_force(1,:)
            WWAVE_FORCE(2,:,IS) = WWAVE_FORCE(2,:,IS) + BBwhc_force(2,:)

          ELSE !3D

            ! Threshold on Hs
            tmp0 = (out_wwm(n1,1) + out_wwm(n2,1))/2.D0 !Hs

!           case PLANARBEACH, failed with "Infinity" value for sum_3D
!           Proposition : Keep the original tests
            IF(tmp0 <= 0.05D0) CYCLE
            IF(tmp0/htot < 0.1D0) CYCLE
!           end JL

            ! Vertical distribution function of qdm (due to wave breaking)
            !swild_2D = 0.D0;
            swild_3D = 0.D0
            DO k = kbs(IS), NVRT

              ! swild_2D(k) = 1.D0
              ! Homogeneous vertical distribution
              IF (ZPROF_BREAK == 1) swild_3D(k) = 1.D0
              IF (ZPROF_BREAK == 2) swild_3D(k) = cosh((dps(IS)+zs(k,IS))/(0.2D0*tmp0))
              IF (ZPROF_BREAK == 3) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**2.D0)
              IF (ZPROF_BREAK == 4) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**4.D0)
              IF (ZPROF_BREAK == 5) swild_3D(k) = 1.D0 - dtanh(((eta_tmp-zs(k,IS))/(0.5D0*tmp0))**8.D0)
              ! All in the two surface layers
              IF (ZPROF_BREAK == 6 .AND. k .GE. NVRT-1) swild_3D(k)=1.D0

            END DO !k

            ! Integral of the vertical distribution function
            !sum_2D = 0.0D0
            sum_3D = 0.0D0
            DO k = kbs(IS), NVRT-1
              !sum_2D = sum_2D + (swild_2D(k+1) + swild_2D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
              sum_3D = sum_3D + (swild_3D(k+1) + swild_3D(k))/2.D0*(zs(k+1,IS) - zs(k,IS))
            END DO !NVRT-1
            !IF(sum_2D*sum_3D == 0) CALL parallel_abort('Vertical profile in wave breaking-induced force: integral=0')
            IF(sum_3D == 0) CALL parallel_abort('Vertical profile in wave breaking-induced force: integral=0')
!'

            DO k = kbs(IS), NVRT
              ! Breaking acceleration + roller contribution
              IF (IROLLER == 1) THEN
                BBrol_force(1,k) = BBrol_force(1,k) -swild_3D(k)*(SROL(1,n1) + SROL(1,n2))/2.D0/sum_3D
                BBrol_force(2,k) = BBrol_force(2,k) -swild_3D(k)*(SROL(2,n1) + SROL(2,n2))/2.D0/sum_3D
                BBbrk_force(1,k) = BBbrk_force(1,k) -swild_3D(k)*(1.D0-ALPROL)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D
                BBbrk_force(2,k) = BBbrk_force(2,k) -swild_3D(k)*(1.D0-ALPROL)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D
              ELSE
                BBbrk_force(1,k) = BBbrk_force(1,k) -swild_3D(k)*(SBR(1,n1) + SBR(1,n2))/2.D0/sum_3D
                BBbrk_force(2,k) = BBbrk_force(2,k) -swild_3D(k)*(SBR(2,n1) + SBR(2,n2))/2.D0/sum_3D
              ENDIF
              ! Whitecapping contribution ? SWAN, already in SBR ?
              BBwhc_force(1,k) = BBwhc_force(1,k) -swild_3D(k)*(SDS(1,n1) + SDS(1,n2))/2.D0/sum_3D
              BBwhc_force(2,k) = BBwhc_force(2,k) -swild_3D(k)*(SDS(2,n1) + SDS(2,n2))/2.D0/sum_3D

              ! --- acceleration due to wave breaking and white capping(+wcap)
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + BBbrk_force(1,k)
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + BBbrk_force(2,k)
              ! --- acceleration due to wave roller
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + BBrol_force(1,k)
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + BBrol_force(2,k)
              ! --- acceleration due to whitecapping
              WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS) + BBwhc_force(1,k)
              WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS) + BBwhc_force(2,k)
            END DO
          END IF !2D/3D

!#JL ??          
          ! Smoothing wave forces near the shoreline
          ! With this profile, F < 10% of computed F at h < DMIN, and F > 95% of computed F at h > 2.25*DMIN
!          IF (htot < 3*DEPMIN) THEN 
!             WWAVE_FORCE(:,:,IS) = WWAVE_FORCE(:,:,IS)*tanh((0.5D0*htot/DEPMIN)**8.D0)
!             BBbrk_force(:,:) = BBbrk_force(:,:)*tanh((0.5D0*htot/DEPMIN)**8.D0)
!             BBrol_force(:,:) = BBrol_force(:,:)*tanh((0.5D0*htot/DEPMIN)**8.D0)
!          ENDIF

          ! --- Save the wave acceleration term due to wave surface breaking [m.s-2]

!#if defined USE_PLANARBEACH || defined USE_DUCK94
!          sum1 = 0.D0 !integral; x-comp.
!          sum2 = 0.D0 !integral
!          do k=kbs(IS)+1,nvrt !IS is wet
!            sum1=sum1+(zs(k,IS)-zs(k-1,IS))*(BBbrk_force(1,k)+BBbrk_force(1,k-1))/2.d0
!            sum1=sum1+(zs(k,IS)-zs(k-1,IS))*(BBrol_force(1,k)+BBrol_force(1,k-1))/2.d0
!!           sum1=sum1+(zs(k,IS)-zs(k-1,IS))*(BBwhc_force(1,k)+BBwhc_force(1,k-1))/2.d0

!            sum2=sum2+(zs(k,IS)-zs(k-1,IS))*(BBbrk_force(2,k)+BBbrk_force(2,k-1))/2.d0
!            sum2=sum2+(zs(k,IS)-zs(k-1,IS))*(BBrol_force(2,k)+BBrol_force(2,k-1))/2.d0
!!           sum2=sum2+(zs(k,IS)-zs(k-1,IS))*(BBwhc_force(2,k)+BBwhc_force(2,k-1))/2.d0
!          enddo !k
!          diag_BRK_bar(1,IS) = diag_BRK_bar(1,IS) + sum1/htot
!          diag_BRK_bar(2,IS) = diag_BRK_bar(2,IS) + sum2/htot
!#endif
        END DO !nsa

        ! Exchange between ghost regions
        CALL exchange_s3d_2(WWAVE_FORCE)

!#if defined USE_PLANARBEACH || defined USE_DUCK94
!        CALL exchange_s2d_2(diag_BRK_bar)
!#endif

      END SUBROUTINE BREAKING_VF_TERMS_SCHISM


!**********************************************************************
!*   This routine fill the bucket diag_WVE_BAR 
!*   #JL developpment for wave force diagnostics... (Not used #if 0) 
!*********************************************************************
#if 0
#if defined USE_PLANARBEACH || defined USE_DUCK94
      SUBROUTINE DIAG_WVE_FORCE
        USE VARS_WAVE  ! only:
        USE schism_glbl, ONLY: kbs, ns, dps, zs,nsa,idry_s,isidenode,eta2,nvrt,diag_WVE_bar
!        use schism_glbl, only: DMIN=>h0 ! Minimum water depth. THis must be same as h0 in SCHISM
        USE schism_glbl, ONLY: h0
        USE SWCOMM3, only: DEPMIN
        USE schism_msgp
        IMPLICIT NONE

        REAL(rkind)  :: tmp1,tmp2,htot
        INTEGER :: IS,k

        ! Init here
        diag_WVE_bar = 0.D0  ! Wave Force

!...  Loop on the resident sides
        do IS = 1,nsa
          if(idry_s(IS) == 1.or.nvrt==kbs(IS)+1) cycle
          tmp1=0.d0; tmp2=0.d0

          htot = dps(IS)+sum(eta2(isidenode(1:2,IS)))/2.d0
          htot = max(h0,htot)

          do k=kbs(IS)+1,nvrt
            !wwave_force at sidecenter and whole level, integrate along z
            tmp1 = tmp1+ (zs(k,IS)-zs(k-1,IS))*(wwave_force(1,k,IS)+wwave_force(1,k-1,IS))/2.d0
            tmp2 = tmp2+ (zs(k,IS)-zs(k-1,IS))*(wwave_force(2,k,IS)+wwave_force(2,k-1,IS))/2.d0
          enddo !k

          diag_WVE_bar(1,IS) = diag_WVE_bar(1,IS) + tmp1/htot
          diag_WVE_bar(2,IS) = diag_WVE_bar(2,IS) + tmp2/htot

        enddo !ns

        print*,'DIAG_WVE_FORCE'
        ! Exchange between ghost regions
        CALL exchange_s2d_2(diag_WVE_bar)
 
      END SUBROUTINE DIAG_WVE_FORCE
#endif
#endif

!**********************************************************************
!*   This routine fixes the wave forces to the barotropic gradient at the numerical shoreline (boundary between dry and wet elements)
!**********************************************************************
      SUBROUTINE SHORELINE_WAVE_FORCES
        USE VARS_WAVE  ! only:
        USE schism_glbl, only: errmsg,hmin_radstress,idry_e,idry_s,isidenode, &
                       & isdel,elnode,i34,dldxy,grav,thetai,ns,idry,isbnd,eta1,eta2
        USE schism_msgp
        IMPLICIT NONE

        INTEGER      :: IP,IS,INODE_1,INODE_2,IELEM
        REAL(rkind)  :: TMP_X,TMP_Y,BPGR(2)

!...  Loop on the resident sides
        do IS = 1,ns
          if(idry_s(IS) == 1) cycle

          ! Adjacent nodes index
          INODE_1 = isidenode(1,IS) ! Side node #1
          INODE_2 = isidenode(2,IS) ! Side node #2
          if(idry(INODE_1) == 1 .OR. idry(INODE_2) == 1) cycle

          ! Sides we are not interested in
          if(isdel(1,IS) == 0 .OR. isdel(2,IS) == 0) cycle ! Boundaries
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 0) cycle ! Case where both adjacent elements are wet
          if(idry_e(isdel(1,IS)) == 1 .AND. idry_e(isdel(2,IS)) == 1) cycle ! Case where both adjacent elements are dry (should never occur anyway)
          !if(isbnd(1,INODE_1) == 1 .OR. isbnd(1,INODE_2) == 1) cycle ! Case where the side touches open boundaries
          if(isbnd(1,INODE_1)>0.OR.isbnd(1,INODE_2)>0) cycle ! Case where the side touches open boundaries

          ! Reinitializing the wave force
          WWAVE_FORCE(:,:,IS) = 0.d0

          ! We are left with sides that belong to one dry element and one wet element
          ! We store the elements indexes for future use
          if(idry_e(isdel(1,IS)) == 0 .AND. idry_e(isdel(2,IS)) == 1) then
            IELEM = isdel(1,IS)
          elseif(idry_e(isdel(2,IS)) == 0 .AND. idry_e(isdel(1,IS)) == 1) then
            IELEM = isdel(2,IS)
          else
            cycle
          endif

          ! We compute the barotropic gradient at the shoreline (only the wet element is used)
          BPGR = 0
          do IP = 1,i34(IELEM)
            ! x and y-components of grad(eta)
            TMP_X = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,1,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,1,IELEM)
            TMP_Y = (1-thetai)*eta1(elnode(IP,IELEM))*dldxy(IP,2,IELEM) + thetai*eta2(elnode(IP,IELEM))*dldxy(IP,2,IELEM)
            ! Barotropic gradient = g*grad(eta)
            BPGR(1) = BPGR(1) + grav*TMP_X
            BPGR(2) = BPGR(2) + grav*TMP_Y
          enddo !IP

          ! Overwriting wave forces to balance out pressure gradient
          WWAVE_FORCE(1,:,IS) = BPGR(1)
          WWAVE_FORCE(2,:,IS) = BPGR(2)
        enddo !IS

        call exchange_s3d_2(WWAVE_FORCE)

      END SUBROUTINE SHORELINE_WAVE_FORCES

!**********************************************************************
!*                                                                    *
!**********************************************************************
!* This routine, called in main, applies a ramp to wave fores starting
!  from the open boundary.
!  This ramp is defined in the input file wafo_ramp.gr3, and is read
!  in wwm_initio  (subroutine READ_WAFO_OPBND_RAMP).
!  If wafo_obcramp==1, APPLY_WAFO_OPBND_RAMP is called in main at
!  each time step.
!  Authors: X. Bertin & B. Mengual (05/2020)
!**********************************************************************
      SUBROUTINE APPLY_WAFO_OPBND_RAMP

        USE VARS_WAVE, only: WWAVE_FORCE,wafo_opbnd_ramp
        USE schism_glbl, only : ns,kbs,idry_s,nsa,nvrt
        USE schism_msgp

        IMPLICIT NONE

        INTEGER      :: IS,k

        DO IS = 1,nsa
          IF(idry_s(IS) == 1) CYCLE
          DO k = kbs(IS),NVRT
            WWAVE_FORCE(1,k,IS) = WWAVE_FORCE(1,k,IS)*wafo_opbnd_ramp(IS)
            WWAVE_FORCE(2,k,IS) = WWAVE_FORCE(2,k,IS)*wafo_opbnd_ramp(IS)
          END DO ! NVRT
        END DO ! nsa

      END SUBROUTINE APPLY_WAFO_OPBND_RAMP
#endif /*USE_SWAN*/


!**********************************************************************
!*   This routine for RADFLAG=LON (Longuet-Higgin)
!*
!* Citation: Roland, A., Y. J. Zhang, H. V. Wang, Y. Meng, Y.-C. Teng, V. Maderich, 
!* I. Brovchenko, M. Dutour-Sikiric, and U. Zanke (2012), A fully coupled 3D 
!* wave-current interaction model on unstructured grids, J. Geophys. Res., 117, C00J33,
!*
!**********************************************************************
#ifdef USE_SWAN
      SUBROUTINE RADIATION_STRESS_SCHISM
        !USE DATAPOOL
!        use schism_glbl, only: iplg,errmsg,hmin_radstress,kbp,wwave_force,idry
        USE VARS_WAVE
        USE SWCOMM3, only: DDIR !=(SPDIR2-SPDIR1)/MDC, mesh size in theta-direction of computational grid
        USE SWCOMM3, only: dmin => DEPMIN, G9=> GRAV_W, PI_W
        USE SWCOMM3, only: SGHIGH => SHIG !=2*PI*FRHIG; highest spectral value of sigma
        USE SWCOMM3, only: SGLOW  => SLOW !=2*PI*FRLOW; lowest spectral value of sigma
        USE SWCOMM3, only: TAIL_ARR  => PWTAIL
        USE M_GENARR, only: SPCDIR,SPCSIG
        USE schism_msgp !, only : myrank,parallel_abort
        use schism_glbl, only: iplg,errmsg,hmin_radstress,msc2,mdc2,time_stamp
        use schism_glbl, only: idry_s,idry,kbp,isidenode,nsa,nvrt
        use schism_glbl, only: MNP => npa,       & ! Nodes of the augmented domain
     &                         MNE => nea,       & ! Elements of the augmented domain
     &                         DEP8 => dp,       & ! depth in the augmented domain
     &                         ipgl,             & ! node global to local mapping
     &                         ielg,             & ! element local to global maping
     &                         eta2,             &
     &                         ZETA => znl       ! Z-Levels of SCHISM

        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL(rkind)  :: ACLOC(mdc2,msc2)
        REAL(rkind)  :: COSE2, SINE2, COSI2
        REAL(rkind)  :: EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL(rkind)  :: m0, m0d, tmp, EHFR, ELOC, EFTAIL, ZZETA
        REAL(rkind)  :: DS, D, KW, KD, SINH2KD, SINHKW, COSH2KW, COSHKW, COSHKD, ETOTS, ETOTC, EWSIG, S11, S22
        REAL(rkind)  :: WNTMP,WKTMP,WCGTMP,WCTMP,WN,WKDEPTMP
        REAL(rkind)  :: WSTMP, DEPLOC
        REAL(rkind)  :: CGLOC,ND,KWAVELOC,NN 

        INTEGER      :: ND1,ND2
        REAL(rkind)  :: SINHKD,FSS(NVRT,MNP),FCS(NVRT,MNP),FSC(NVRT,MNP),FCC(NVRT,MNP)
        REAL(rkind)  :: dr_dxy(2,NVRT,nsa),HTOT,SXX3D0(NVRT,MNP),SYY3D0(NVRT,MNP),SXY3D0(NVRT,MNP)
        REAL(rkind)  :: WILD1(NVRT,MNP),WILD2(NVRT,MNP),WILD3(2,NVRT,nsa),WILD4(3,NVRT,MNP),DSPX,DSPY, WILD5(10)
        REAL(rkind), PARAMETER :: ZERO = 0._rkind
        REAL(rkind), PARAMETER :: ONE = 1._rkind
        REAL(rkind), PARAMETER :: TWO = 2._rkind
        REAL(rkind), PARAMETER :: WAVE_LENGTH_MIN = 0.01_rkind

!GD: imet_dry allows to choose between 2 different methods to compute the
!derivative at the sides between wet and dry elements:
!
! imet_dry=1 : only the values at the 2 nodes of the side are used to
! compute the derivative (this older method showed to provide inconsistent
! wave force at the wet/dry interface).
!
! imet_dry=2 : a 4-point stencil (the 3 wet nodes and an artificial
! node at the center of the side) is used to compute the derivative.
! This method is similar to using shape functions to compute the
! derivative at the center of the element and assigning this value to the
! the side center.

        IMET_DRY = 2

        SXX3D(:,:) = ZERO
        SYY3D(:,:) = ZERO
        SXY3D(:,:) = ZERO
        HTOT = dmin
        EFTAIL = ONE / (TAIL_ARR(1) - ONE)

        ETOT = ZERO
        MDIR = ZERO

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
!            IF (ABS(IOBP(IP)) .GT. 0) CYCLE
            !IF (DEP(IP) .LT. DMIN) CYCLE
            IF (IDRY(IP)==1) CYCLE

            DEPLOC = DEP8(IP)
!JL : in WWM, DEP8==dp not dp+eta2, why ?
            ACLOC(1:mdc2,1:msc2) = AC2(1:mdc2,1:msc2,IP)
            m0    = ZERO
            EWSIG  = ZERO
            ETOTS  = ZERO
            ETOTC  = ZERO

!            if(myrank.EQ.0.AND.modulo(IP,500).EQ.0) &
! &          print*,'RADIATION_STRESS_SCHISM,step 0: index_gl,depth,AC2(1:10,1:10)',&
! &          iplg(IP),DEPLOC,ACLOC(1:10,1:10)

            IF (msc2 .GE. 2) THEN
              DO ID = 1, mdc2
                m0d = ZERO
                DO IS = 2, msc2
                  tmp = 0.5_rkind*(SPCSIG(IS)*ACLOC(ID,IS)+SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPCSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (msc2 > 3) THEN
                  EHFR = ACLOC(ID,msc2) * SPCSIG(msc2)
                  m0 = m0 + DDIR * EHFR * SPCSIG(msc2) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * SPCDIR(ID,2) !COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SPCDIR(ID,3) !SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIGH - SGLOW
              DO ID = 1, mdc2
                m0d = ACLOC(ID,1) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. ZERO .and. .not. DEPLOC .lt. dmin) then
              EWS(IP) = EWSIG/m0
              WSTMP = EWSIG/m0
              !CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
!             Calculation of the wave number, group velocity, group number N
!             and the derivative of N w.r.t. depth (=ND)
!             CGLOC   output   group velocity
!             DEPLOC  input    local depth
!             WKTMP   output   wave number
!             WNTMP   output   ratio of group and phase velocity
!             ND      output   derivative of N with respect to D
!             WSTMP   input    rel. frequency for which wave parameters must be determined
              CALL KSCIP1(1,WSTMP,DEPLOC,WKTMP,CGLOC,WNTMP,ND)
              EWN(IP) = WNTMP      ! Group Number N
              EWK(IP) = WKTMP      ! Wave Number K
              MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
            ELSE
              EWS(IP)  = ZERO
              EWN(IP)  = ZERO
              EWK(IP)  = 10._rkind
              MDIR(IP) = ZERO
            END IF
#if 0
            if(myrank.EQ.0.AND.modulo(IP,500).EQ.0) THEN
             print*,'RADIATION_STRESS_SCHISM,step 0: index_gl,HS,DM,KLM,dep,wlv',&
 &           iplg(IP),out_wwm(IP,1),out_wwm(IP,9),out_wwm(IP,5),out_wwm(IP,34),out_wwm(IP,35)
             print*,'RADIATION_STRESS_SCHISM,step 1: index_gl,dep,m0,ews,ewn,ewk,mdir',&
 &           iplg(IP),DEPLOC,m0,EWS(IP),EWN(IP),EWK(IP),MDIR(IP)*180._rkind/PI_W 
             print*,'RADIATION_STRESS_SCHISM, check HS=4*sqrt(m0)', 4*sqrt(m0),out_wwm(IP,1)
             print*,'RADIATION_STRESS_SCHISM, check WV_number',EWK(IP),TWO*PI_W/MAX(out_wwm(IP,6),WAVE_LENGTH_MIN)
            endif
#endif

          END DO !IP -> MNP
        END IF !LETOT

!AR: Here comes the whole story ...
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2)
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono^X^O
! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!

        SXX3D = ZERO
        SXY3D = ZERO
        SYY3D = ZERO
        WWAVE_FORCE = ZERO

        IF (RADFLAG .EQ. 'LON') THEN
          RSXX = ZERO
          RSXY = ZERO
          RSYY = ZERO
          DO IP = 1, MNP
            DEPLOC = DEP8(IP)
            IF (IDRY(IP)==1) CYCLE

            IF (.NOT. LETOT) THEN
              ACLOC(1:mdc2,1:msc2) = AC2(1:mdc2,1:msc2,IP)
              DO ID = 1, mdc2
                DO IS = 2, msc2
                  ELOC  = 0.5_rkind * (SPCSIG(IS)*ACLOC(ID,IS)+SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS_INCR(IS)*DDIR
                  COSE2 = SPCDIR(ID,4) !COS(SPDIR(ID))**TWO  ! SPCDIR(*,4); cosine^2 of spectral directions
                  SINE2 = SPCDIR(ID,6) !SIN(SPDIR(ID))**TWO  ! SPCDIR(*,6); sine^2 of spectral directions
                  COSI2 = SPCDIR(ID,5) !COS(SPDIR(ID)) * SIN(SPDIR(ID)) ! SPCDIR(*,5); cosine*sine of spectral directions

                  CALL KSCIP1(1,SPCSIG(IS),DEPLOC,KWAVELOC,CGLOC,NN,ND)
                  !WN    = CG(IS,IP) / ( SPCSIG(IS)/WK(IS,IP) )
                  WN    = CGLOC / ( SPCSIG(IS)/KWAVELOC)
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - 0.5_rkind) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2          ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - 0.5_rkind) * ELOC
                ENDDO
              ENDDO
            ELSE IF (LETOT) THEN
              RSXX(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*SIN(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
              RSXY(IP) =  ETOT(IP) *  EWN(IP)* EWK(IP)*SIN(MDIR(IP))*EWK(IP)*COS(MDIR(IP))* ONE/EWK(IP)
              RSYY(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*COS(MDIR(IP)))**TWO/EWK(IP)**TWO+ONE)-0.5_rkind)
            END IF
#if 0
            if(myrank.EQ.0.AND.modulo(IP,500).EQ.0) &
 &          print'(A,f16.1,I7,3(1X,E12.4))','RADIATION_STRESS_SCHISM,step 2: index_gl,RSXX,RSXY,RSYY',&
 &          time_stamp,iplg(IP),RSXX(IP),RSXY(IP),RSYY(IP)
#endif
          END DO ! MNP

          DO IP = 1, MNP
            IF (IDRY(IP)==1) THEN
              SXX3D(:,IP) = ZERO
              SXY3D(:,IP) = ZERO
              SYY3D(:,IP) = ZERO
            ELSE
              SXX3D(:,IP) = RSXX(IP) / DEP8(IP) * G9
              SXY3D(:,IP) = RSXY(IP) / DEP8(IP) * G9
              SYY3D(:,IP) = RSYY(IP) / DEP8(IP) * G9
            END IF
          END DO
          !Store as double for force later
          SXX3D0 = SXX3D
          SXY3D0 = SXY3D
          SYY3D0 = SYY3D

        ELSE
          call parallel_abort('R.S.: unknown R.S. model')
        END IF !RADFLAG

!       Integrate over depth for checking
        RSXX = ZERO
        DO IP = 1, MNP
          IF(IDRY(IP)==1) CYCLE
          DO IL = KBP(IP)+1, NVRT
            RSXX(IP) = RSXX(IP) + 0.5_rkind*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) * ABS((ZETA(IL,IP) - ZETA(IL-1,IP)))/G9
          END DO !IL
        END DO !IP
!'
!       Computation in double precision here
!       SXX3D0() etc. should have dimension of m^2/s/s, defined at nodes and whole levels.
!       Use same arrays to temporarily store properly scaled Sxx etc
!       write(12,*)'Checking Sxx,Sxy,Syy:'
        do IP=1,MNP
          IF(IDRY(IP)==1) then
            SXX3D0(:,IP)=ZERO
            SYY3D0(:,IP)=ZERO
            SXY3D0(:,IP)=ZERO
            cycle
          endif

          do IL=KBP(IP),NVRT
!           D*(Sxx, Sxy, Syy)/rho in Xia et al. (2004)
!           After this the dimension should be m^3/s/s
            tmp=max(DEP8(IP)+ETA2(IP),hmin_radstress)
            !SXX3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXX3D0(IL,IP) !D*Sxx/rho
            !SXY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXY3D0(IL,IP) !D*Sxy/rho
            !SYY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SYY3D0(IL,IP) !D*Syy/rho
            SXX3D0(IL,IP)=tmp*SXX3D0(IL,IP) !D*Sxx/rho
            SXY3D0(IL,IP)=tmp*SXY3D0(IL,IP) !D*Sxy/rho
            SYY3D0(IL,IP)=tmp*SYY3D0(IL,IP) !D*Syy/rho
          enddo !k
        enddo !IP
!'
!       Compute radiation stress force
!       wwave_force(:,1:nsa,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!       and has a dimension of m/s/s
        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXX3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)

        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
!            write(*,*) sum(eta2), sum(dep8)
!            write(*,*) HTOT, eta2(isidenode(1,IS)), eta2(isidenode(2,IS)), DEP8(isidenode(1,IS)), DEP8(isidenode(1,IS)), isidenode(1,IS), isidenode(1,IS)
!            if (isidenode(1,IS) == 150 .AND. isidenode(2,IS) == 149) stop
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (999)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, nvrt
              WWAVE_FORCE(1,il,IS)=WWAVE_FORCE(1,il,IS)-dr_dxy(1,il,IS)/HTOT
            end do
          endif
        enddo !IS
        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SYY3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)

        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (998)')
            HTOT=max(HTOT,hmin_radstress)
            do il = 1, nvrt
              WWAVE_FORCE(2,il,IS)=WWAVE_FORCE(2,il,IS)-dr_dxy(2,il,IS)/HTOT
            end do
          endif
        enddo !IS

        call hgrad_nodes(IMET_DRY,0,NVRT,MNP,nsa,SXY3D0,dr_dxy)
        call exchange_s3d_2(dr_dxy)
!
!          write(12,*)'Checking R.S.'
!
        do IS=1,nsa
          if(idry_s(IS)==0) then
            HTOT=(eta2(isidenode(1,IS))+eta2(isidenode(2,IS))+DEP8(isidenode(1,IS))+DEP8(isidenode(2,IS)))/2
            if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (997)')
            HTOT=max(HTOT,hmin_radstress)
            WWAVE_FORCE(1,:,IS)=WWAVE_FORCE(1,:,IS)-dr_dxy(2,:,IS)/HTOT
            WWAVE_FORCE(2,:,IS)=WWAVE_FORCE(2,:,IS)-dr_dxy(1,:,IS)/HTOT
          endif

#if 0
            if(myrank.EQ.0.AND.modulo(IS,500).EQ.0) &
 &          print'(A,f16.1,I7,2(1X,E12.4))','WAVE_FORCE:',&
 &          time_stamp,IS,WWAVE_FORCE(1,nvrt,IS),WWAVE_FORCE(2,nvrt,IS)
#endif

        enddo !IS

      END SUBROUTINE RADIATION_STRESS_SCHISM
#endif /*USE_SWAN*/ 


#ifdef USE_SWAN
!*********************************************************************
!* This routine for RADFLAG=WON (Longuet-Higgin)
!* Second alternative : original UNSWAN-ADCIRC method and logic      *
!---------------------------------------------------------------------
      SUBROUTINE SwanComputeForce
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
!     Copyright (C) 1993-2017  Delft University of Technology
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
!   Authors
!
!   40.80: Marcel Zijlema
!   41.20: Casey Dietrich
!
!   Updates
!
!   40.80, February 2008: New subroutine
!   41.20, February 2008: extension to tightly coupled ADCIRC+SWAN model
!
!   Purpose
!
!   Computes wave-induced force in vertices
!
!   Method
!
!   First, compute radiation stresses in vertices
!   Next, compute gradients of radiation stresses in vertices
!   Finally, compute wave-induced force in vertices
!
!   Modules used
!
    USE VARS_WAVE
!    USE SWCOMM3, only: DDIR !=(SPDIR2-SPDIR1)/MDC, mesh size in theta-direction of computational grid
!    USE SWCOMM3, only: dmin => DEPMIN, G9=> GRAV_W, PI_W
!    USE SWCOMM3, only: SGHIGH => SHIG !=2*PI*FRHIG; highest spectral value of sigma
!    USE SWCOMM3, only: SGLOW  => SLOW !=2*PI*FRLOW; lowest spectral value of sigma
!    USE SWCOMM3, only: TAIL_ARR  => PWTAIL
    USE M_GENARR, only: SPCDIR,SPCSIG

    USE schism_msgp !, only : myrank,parallel_abort
    use schism_glbl, only: iplg,errmsg,hmin_radstress,kbp,msc2,mdc2,idry,nvrt, &
                           dp,eta2,znl,npa,nsa,idry_s,kbp,time_stamp,isidenode
    !use schism_glbl, only: nverts => npa,       & ! Nodes of the augmented domain
    ! &                        DEP8 => dp,       & ! depth in the augmented domain
    ! &                        ipgl,             & ! node global to local mapping
    ! &                        ielg,             & ! element local to global maping
    ! &                        ZETA => znl       ! Z-Levels of SCHISM
    use ocpcomm4 
    use swcomm2
    use swcomm3
    use swcomm4, ONLY : KSPHER, LENDEG!, COSLAT
    use SwanGriddata
    use SwanGridobjects
!ADC    use Couple2Adcirc, only: compda
!ADC    use GLOBAL,        only: rsnx2, rsny2
!
    implicit none
!NADC!
!NADC!   Argument variables
!NADC!
!NADC    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
!NADC    real, dimension(nverts), intent(in)         :: dep2   ! water depth at current time level
!NADC    real, dimension(nverts), intent(out)        :: fx     ! wave-induced force in x-direction
!NADC    real, dimension(nverts), intent(out)        :: fy     ! wave-induced force in y-direction
!NADC    real, dimension(MDC,6), intent(in)          :: spcdir ! (*,1): spectral direction bins (radians)
!NADC                                                          ! (*,2): cosine of spectral directions
!NADC                                                          ! (*,3): sine of spectral directions
!NADC                                                          ! (*,4): cosine^2 of spectral directions
!NADC                                                          ! (*,5): cosine*sine of spectral directions
!NADC                                                          ! (*,6): sine^2 of spectral directions
!NADC    real, dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell    ! index of present cell
    integer                               :: id       ! loop counter over direction bins
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: jc       ! loop counter
    integer                               :: jcell    ! index of next cell
    integer                               :: IL
    !
    integer, dimension(3)                 :: v        ! vertices in present cell
    !
    double precision                      :: area     ! twices the area of centroid dual around present vertex
    real(rkind)                           :: cslat    ! cosine of latitude
    real(rkind)                           :: deploc   ! local depth
    REAL(rkind)                           :: ACLOC(mdc2,msc2)
    real(rkind)                           :: dsxxdx   ! x-gradient of sxx
    real(rkind)                           :: dsxydx   ! x-gradient of sxy
    real(rkind)                           :: dsxydy   ! y-gradient of sxy
    real(rkind)                           :: dsyydy   ! y-gradient of syy
    real(rkind)                           :: sxx0     ! sxx in centroid of present cell
    real(rkind)                           :: sxx1     ! sxx in centroid of next cell
    real(rkind)                           :: sxy0     ! sxy in centroid of present cell
    real(rkind)                           :: sxy1     ! sxy in centroid of next cell
    real(rkind)                           :: syy0     ! syy in centroid of present cell
    real(rkind)                           :: syy1     ! syy in centroid of next cell
    real(rkind)                           :: sxxsumd  ! cumulated sxx over direction space
    real(rkind)                           :: sxysumd  ! cumulated sxy over direction space
    real(rkind)                           :: syysumd  ! cumulated syy over direction space
    real(rkind)                           :: sxxsums  ! cumulated sxx over frequency space
    real(rkind)                           :: sxysums  ! cumulated sxy over frequency space
    real(rkind)                           :: syysums  ! cumulated syy over frequency space
    double precision                      :: x0       ! x-coordinate of the centroid of present cell
    double precision                      :: x1       ! x-coordinate of the centroid of next cell
    double precision                      :: y0       ! y-coordinate of the centroid of present cell
    double precision                      :: y1       ! y-coordinate of the centroid of next cell
    !
    real(rkind), dimension(1)                    :: cg       ! group velocity
    real(rkind), dimension(1)                    :: k        ! wave number
    real(rkind), dimension(1)                    :: n        ! ratio of group and phase velocity
    real(rkind), dimension(1)                    :: nd       ! derivative of n with respect to depth
    real(rkind), dimension(1)                    :: sig      ! relative frequency
    real(rkind), dimension(nverts)               :: dep2     ! help array to store water depth
!    real(rkind), dimension(nverts)               :: fx       ! help array to store wave-induced force in x-direction
!    real(rkind), dimension(nverts)               :: fy       ! help array to store wave-induced force in y-direction
    !
!    real(rkind), dimension(:), allocatable       :: sxx      ! x-component of radiation stress in x-direction
!    real(rkind), dimension(:), allocatable       :: sxy      ! cross component of radiation stress in x/y-direction
!    real(rkind), dimension(:), allocatable       :: syy      ! y-component of radiation stress in y-direction
!    replaced by RSXX,RSXY,RSYY
    !
    character(80)                         :: msgstr   ! string to pass message
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes

    REAL(rkind) :: HTOT,tmp
    REAL(rkind), PARAMETER :: ZERO = 0._rkind
    REAL(rkind), PARAMETER :: ONE = 1._rkind
    REAL(rkind), PARAMETER :: TWO = 2._rkind
    REAL(rkind), PARAMETER :: WAVE_LENGTH_MIN = 0.01_rkind

    REAL(rkind) :: xmin,xmax

!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! allocation and initialization of radiation stresses
    !
    RSXX = ZERO
    RSXY = ZERO
    RSYY = ZERO
    dep2 = ZERO
    !
    ! compute radiation stresses in vertices
    !
    do ivert = 1, nverts ! nverts == npa
       !
       dep2(ivert) = dp(ivert) + eta2(ivert)
       deploc = dep2(ivert)
       acloc(1:mdc2,1:msc2) = ac2(1:mdc2,1:msc2,ivert)
       !
       if ( deploc <= DEPMIN ) cycle
       !
       sxxsums = ZERO
       sxysums = ZERO
       syysums = ZERO
       !
       ! calculate sum of each component contributed to the radiation stresses over all frequencies
       !
       do is = 1, msc2
          !
          ! compute group velocity over phase velocity (=n)
          !
          sig(1) = spcsig(is)
          call KSCIP1 (1,sig,deploc,k,cg,n,nd)
          !
          sxxsumd = ZERO
          sxysumd = ZERO
          syysumd = ZERO
          !
          ! calculate sum of each component contributed to the radiation stresses over all wave directions
          !
          do id = 1, mdc2
             !
             sxxsumd = sxxsumd + ( n(1)*(spcdir(id,4) + ONE) -0.5_rkind) * acloc(id,is)
             sxysumd = sxysumd +   n(1)* spcdir(id,5)                    * acloc(id,is)
             syysumd = syysumd + ( n(1)*(spcdir(id,6) + ONE) -0.5_rkind) * acloc(id,is)
             !
          enddo ! mdc2
          !
          ! integrate over frequency space
          !
          sxxsums = sxxsums + sxxsumd * spcsig(is)**2
          sxysums = sxysums + sxysumd * spcsig(is)**2
          syysums = syysums + syysumd * spcsig(is)**2
          !
       enddo ! msc2
       !
       ! compute the radiation stresses (divided by rho times gravitational acceleration)
       ! FRINTF is the frequency integration factor (=df/f)
       ! DDIR ise the band width in directional space = PI2_W / MDC
       RSXX(ivert) = sxxsums * DDIR * FRINTF
       RSXY(ivert) = sxysums * DDIR * FRINTF
       RSYY(ivert) = syysums * DDIR * FRINTF
       !

#if 0
            if(myrank.EQ.0.AND.modulo(ivert,500).EQ.0) &
 &          print'(A,f16.1,I7,3(1X,E12.4))','RADIATION_STRESS_SCHISM,step 2: index_gl,RSXX,RSXY,RSYY',&
 &          time_stamp,iplg(ivert),RSXX(ivert),RSXY(ivert),RSYY(ivert)
#endif

    enddo  !nverts

    ! compute wave-induced force in vertices
    !
    WAVESTRX_2D = ZERO ! fx
    WAVESTRY_2D = ZERO ! fy
    WAVESTRX_3D = ZERO
    WAVESTRY_3D = ZERO
    !
    vertexloop : do ivert = 1, nverts
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) cycle vertexloop    ! boundary vertex
       !
       ! first, compute gradients of the radiation stresses in vertices
       !
       area   = 0d0
       dsxxdx = ZERO
       dsxydx = ZERO
       dsxydy = ZERO
       dsyydy = ZERO
       xmax = 0._rkind
       xmin = huge(1)
       !
       ! loop over cells around considered vertex
       !
       do jc = 1, vert(ivert)%noc
          !
          ! get present cell and its vertices
          !
          icell = vert(ivert)%cell(jc)%atti(CELLID)
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop
          !
          ! determine centroid of present cell
          !
          x0 = MyREAL(cell(icell)%attr(CELLCX))
          y0 = MyREAL(cell(icell)%attr(CELLCY))
          if(x0<xmin) xmin = x0
          if(x0>xmax) xmax = x0
          !
          ! determine radiation stresses in centroid in present cell
          !
          sxx0 = ( RSXX(v(1)) + RSXX(v(2)) + RSXX(v(3)) )/ 3._rkind
          sxy0 = ( RSXY(v(1)) + RSXY(v(2)) + RSXY(v(3)) )/ 3._rkind
          syy0 = ( RSYY(v(1)) + RSYY(v(2)) + RSYY(v(3)) )/ 3._rkind
          !
          ! get next cell in counterclockwise direction
          !
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          !
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          !
          ! determine centroid of next cell
          !
          x1 = MyREAL(cell(jcell)%attr(CELLCX))
          y1 = MyREAL(cell(jcell)%attr(CELLCY))
          if(x1<xmin) xmin = x1
          if(x1>xmax) xmax = x1
          !
          ! determine radiation stresses in centroid of next cell
          !
          sxx1 = ( RSXX(v(1)) + RSXX(v(2)) + RSXX(v(3)) )/ 3._rkind
          sxy1 = ( RSXY(v(1)) + RSXY(v(2)) + RSXY(v(3)) )/ 3._rkind
          syy1 = ( RSYY(v(1)) + RSYY(v(2)) + RSYY(v(3)) )/ 3._rkind
          !
          ! compute contribution to area of centroid dual
          !
          area = area + x0*y1 - x1*y0
          !
          ! compute contribution to x-gradient of radiation stresses sxx and sxy
          !
          dsxxdx = dsxxdx + ( sxx0 + sxx1 ) * real( y1 - y0 )
          dsxydx = dsxydx + ( sxy0 + sxy1 ) * real( y1 - y0 )
          !
          ! compute contribution to y-gradient of radiation stresses sxy and syy
          !
          dsxydy = dsxydy + ( sxy0 + sxy1 ) * real( x1 - x0 )
          dsyydy = dsyydy + ( syy0 + syy1 ) * real( x1 - x0 )
          !
       enddo ! loop over cells around considered vertex

       !
       ! Case Crossing Dateline: in centroid dual, some elements are on both sides of the dateline
       !
       IF(KSPHER.EQ.1.and.(xmax-xmin).gt.180.0) THEN
        area   = 0d0
        dsxxdx = ZERO
        dsxydx = ZERO
        dsxydy = ZERO
        dsyydy = ZERO
        ! loop over cells around considered vertex
        do jc = 1, vert(ivert)%noc
          ! get present cell and its vertices
          icell = vert(ivert)%cell(jc)%atti(CELLID)
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          if ( dep2(v(1)) <= DEPMIN .or. dep2(v(2)) <= DEPMIN .or. dep2(v(3)) <= DEPMIN ) cycle vertexloop
          ! determine centroid of present cell
          x0 = MyREAL(cell(icell)%attr(CELLCX))
          y0 = MyREAL(cell(icell)%attr(CELLCY))
          x0 = mod(x0+360.0,360.0)
          ! determine radiation stresses in centroid in present cell
          sxx0 = ( RSXX(v(1)) + RSXX(v(2)) + RSXX(v(3)) )/ 3._rkind
          sxy0 = ( RSXY(v(1)) + RSXY(v(2)) + RSXY(v(3)) )/ 3._rkind
          syy0 = ( RSYY(v(1)) + RSYY(v(2)) + RSYY(v(3)) )/ 3._rkind
          ! get next cell in counterclockwise direction
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          ! determine centroid of next cell
          x1 = MyREAL(cell(jcell)%attr(CELLCX))
          y1 = MyREAL(cell(jcell)%attr(CELLCY))
          x1 = mod(x1+360.0,360.0)
          ! determine radiation stresses in centroid of next cell
          sxx1 = ( RSXX(v(1)) + RSXX(v(2)) + RSXX(v(3)) )/ 3._rkind
          sxy1 = ( RSXY(v(1)) + RSXY(v(2)) + RSXY(v(3)) )/ 3._rkind
          syy1 = ( RSYY(v(1)) + RSYY(v(2)) + RSYY(v(3)) )/ 3._rkind
          ! compute contribution to area of centroid dual
          area = area + x0*y1 - x1*y0
          ! compute contribution to x-gradient of radiation stresses sxx and sxy
          dsxxdx = dsxxdx + ( sxx0 + sxx1 ) * real( y1 - y0 )
          dsxydx = dsxydx + ( sxy0 + sxy1 ) * real( y1 - y0 )
          ! compute contribution to y-gradient of radiation stresses sxy and syy
          dsxydy = dsxydy + ( sxy0 + sxy1 ) * real( x1 - x0 )
          dsyydy = dsyydy + ( syy0 + syy1 ) * real( x1 - x0 )
        enddo ! loop over cells around considered vertex
       ENDIF

       !
       ! if area is non-positive, give error and go to next vertex
       !
       if ( .not. area > 0_dkind) then
!          write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
!          call msgerr( 2, trim(msgstr) )
!          return
           msgstr = 'Wrong cell ordering around centroid dual : CW instead CCW'
           WRITE(errmsg,*) MSGSTR
           CALL parallel_abort(errmsg)
       endif
       !
       dsxxdx =  dsxxdx/real(area)
       dsxydx =  dsxydx/real(area)
       dsxydy = -dsxydy/real(area)
       dsyydy = -dsyydy/real(area)
       !
       ! in case of spherical coordinates, transform back to Cartesian coordinates
       !
       if ( KSPHER.GT.0 ) then
          !
          cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
          !
          dsxxdx = dsxxdx/(cslat * LENDEG)
          dsxydx = dsxydx/(cslat * LENDEG)
          dsxydy = dsxydy/LENDEG
          dsyydy = dsyydy/LENDEG
          !
       endif
       !
       ! finally, compute wave-induced force
       !
       !WAVESTRX_2D(ivert) = -RHO_W * GRAV_W * ( dsxxdx  + dsxydy ) ! fx
       !WAVESTRY_2D(ivert) = -RHO_W * GRAV_W * ( dsxydx  + dsyydy ) ! fy
       WAVESTRX_2D(ivert) = -GRAV_W * ( dsxxdx  + dsxydy )
       WAVESTRY_2D(ivert) = -GRAV_W * ( dsxydx  + dsyydy )

#if 0
       if(myrank.EQ.0.AND.modulo(ivert,500).EQ.0) &
 &      print'(A,f16.1,I7,2(1X,E12.4))','WAVESTRX:',&
 &      time_stamp,iplg(ivert),WAVESTRX_2D(ivert),WAVESTRY_2D(ivert)
#endif

       HTOT = dp(ivert)
       WAVESTRX_2D(ivert) = WAVESTRX_2D(ivert) / HTOT
       WAVESTRY_2D(ivert) = WAVESTRY_2D(ivert) / HTOT

       HTOT = max(dep2(ivert),hmin_radstress)
       DO IL = kbp(ivert), nvrt
          WAVESTRX_3D(IL,ivert) = WAVESTRX_2D(ivert) * HTOT 
          WAVESTRY_3D(IL,ivert) = WAVESTRY_2D(ivert) * HTOT 
       END DO !IL

    enddo vertexloop

    if(.NOT.SERIAL) CALL exchange_p3dw(WAVESTRX_3D)  ! dims= (nvrt,npa)
    if(.NOT.SERIAL) CALL exchange_p3dw(WAVESTRY_3D)
!
! Compute wave force along side
!
    WWAVE_FORCE = ZERO

    do IS=1,nsa
       if(idry_s(IS)==0) then
         HTOT = MAX((dep2(isidenode(1,IS))+dep2(isidenode(2,IS)))/2.0D0,hmin_radstress)
         if(HTOT.LE.ZERO) call parallel_abort('RADIATION_STRESS: (999)')
         do IL = 1, nvrt
              WWAVE_FORCE(1,IL,IS) =  0.5_rkind*(WAVESTRX_3D(IL,isidenode(1,IS)) + WAVESTRX_3D(IL,isidenode(2,IS))) / HTOT
              WWAVE_FORCE(2,IL,IS) =  0.5_rkind*(WAVESTRY_3D(IL,isidenode(1,IS)) + WAVESTRY_3D(IL,isidenode(2,IS))) / HTOT
         end do
       endif
#if 0
            if(myrank.EQ.0.AND.modulo(IS,500).EQ.0) &
 &          print'(A,f16.1,I7,2(1X,E12.4))','WAVE_FORCE:',&
 &          time_stamp,IS,WWAVE_FORCE(1,nvrt,IS),WWAVE_FORCE(2,nvrt,IS)
#endif

    enddo !IS

    !
    ! deallocation of radiation stresses
    !
!    deallocate(sxx)
!    deallocate(sxy)
!    deallocate(syy)
    !
end subroutine SwanComputeForce
#endif /*USE_SWAN*/

!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL(rkind) FUNCTION DVEC2RAD(U,V)

         USE schism_glbl, ONLY: rkind
         USE SWCOMM3, ONLY: PI_W
         IMPLICIT NONE
         REAL(rkind)            :: U,V
!        REAL(rkind)  :: DVEC2RAD

         DVEC2RAD = MyATAN2(V,U) * 180.0_rkind/PI_W
         IF (DVEC2RAD < 0.0_rkind) DVEC2RAD = DVEC2RAD + 360.0_rkind
         DVEC2RAD = DVEC2RAD * PI_W/180.0_rkind

      END FUNCTION


END MODULE MOD_WAVE_CURRENT


!#JL EXPERIMENTAL
#ifdef USE_SWAN
      SUBROUTINE Manning2Madsen
!-----------------------------------------------------------------------
! Adapted from ADCIRC-SWAN :
!
! Casey Dietrich : This routine will convert the ADCIRC Manning's n values
!    into roughness lengths 'KN' that can be used with the Madsen friction 
!    formulation inside SWAN. 
!
! Jerome Lefevre: If nchi = -1 (Manning), update the hydrodynamic bottom 
!                 roughness lengh "rough_p"  
!
! see O.S. Madsen,
!  Y.K. Poon, H.C. Graber (1988). Spectral wave attenuation by bottom
!  friction: Theory.” Proc. 21st Int. Conf. Coastal Engineering, ASCE, 492-504.
! see Rational: https://ccht.ccee.ncsu.edu/integral-coupling-of-bottom-friction/
!
! Authors: J.Lefevre (01/2023)
! Note:  KN ~ z0 = bottom roughness length scale (m) in notations
!       - SWAN applies a lower limit of 0.05m for KN
!       - depending on nchi in param.nml, SCHISM may use Manning "N" or "rough_p", 
!         the hydrodynamic rougness length, actually imposed from manning.gr3 
!         or rough.gr3 respectively
!
!         - For example, in sed_friction.F90, rough_p is computed using : 
!
!             rough_p = MAX(z0s , MAX(z0cr , z0wr+z0bld)+z0sw)
!
!              with z0s Nikuradse roughness (=d50/12.)  ; z0cr Current ripple roughness length(m) ;
!                   z0sw Sand waves roughness length (m); z0wr  Wave ripple roughness length (m)
!
!         and rough_p is finaly the "max apparent roughness length"

!       - With nchi=-1 and SWAN = On, this routine compute at each time step:
!         - rough_p , the "max apparent roughness length" using the              

      use VARS_WAVE
      use schism_msgp !, only : myrank,parallel_abort
      use schism_glbl, only: rmanning,hmin_man,grav,idry,dp,eta2, &
                           npa,rough_p
      use schism_glbl, only: ipre2,xnd,ynd,out_dir,lfdb,fdb,iplg,np,len_out_dir

      IMPLICIT NONE

      INTEGER  :: i

      REAL(rkind) :: htot
      REAL(rkind) :: K = 0.4d0 ! von Karman constant
      REAL(rkind) :: N
      REAL(rkind) :: Z0

!$OMP parallel do default(shared) private(i)
      do i=1,npa
          if(idry(i)==1) cycle
!         Wet node
          htot = max(hmin_man,dp(i)+eta2(i)) !>0

          N  = rmanning(i)

! From Bretschneider et al. (1986) and by assuming that 
! hydraulic radius ~ Total depth, the following relationship between 
! Manning’s n and the friction length Z0 can be de deduced (Casey):
          Z0 = ( htot ) * EXP( -1.d0 * ( 1.d0 + K * htot**(1.d0/6.d0) &
     &      / ( N * SQRT(grav) ) ) )

! Update the total apparent roughness (m) used in hydro/schism_step to 
! compute the enhancement bottom friction due to wave, see wbl_Soulsby97
          rough_p(i) = Z0

!#Casey 110518: Enforce a lower limit on the Manning's n seen by SWAN.          
          if (N.LT.0.03d0) then
            Z0 = ( htot ) * EXP( -1.d0 * ( 1.d0 + K * htot**(1.d0/6.d0) &
     &      / ( 0.03d0 * SQRT(grav) ) ) ) 
          endif

! Update the total apparent roughness (m) used in SWAN 
!#Casey 091216: If we get a junk number, then use the default value.
          kn(i) = max(0.05d0,Z0)
      enddo !i

!     Dump Cdp for diagnostics
      if(ipre2/=0) then

        fdb='KN_000000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
        open(10,file=out_dir(1:len_out_dir)//fdb,status='replace')
        write(10,*)np,nproc
        do i=1,np
          write(10,'(i11,3(1x,e20.12))')iplg(i),xnd(i),ynd(i),kn(i)
        enddo !i
        close(10)
        if(myrank==0) write(16,*)'KN_ output done...'

        fdb='rough_000000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-5:lfdb),'(i6.6)') myrank
        open(10,file=out_dir(1:len_out_dir)//fdb,status='replace')
        write(10,*)np,nproc
        do i=1,np
          write(10,'(i11,3(1x,e20.12))')iplg(i),xnd(i),ynd(i),rough_p(i)
        enddo !i
        close(10)
        if(myrank==0) write(16,*)'rough_ output done...'

        call parallel_finalize
        stop
      endif

!$OMP end parallel do

!C-----------------------------------------------------------------------
      END SUBROUTINE Manning2Madsen
!C-----------------------------------------------------------------------
#endif
