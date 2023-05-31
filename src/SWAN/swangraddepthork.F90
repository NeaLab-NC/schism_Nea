#include "swan_functions.h"

subroutine SwanGradDepthorK ( dep2, mudl2, spcsig, dhdx, dhdy, dkdx, dkdy, ivert )
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
!   Authors
!
!   41.60: Marcel Zijlema
!
!   Updates
!
!   41.60, July 2015: New subroutine
!
!   Purpose
!
!   Computes gradients of depth or wave number in vertex
!   meant for computing turning rate
!
!   Method
!
!   Application of the Green-Gauss theorem with the assumption of
!   a constant gradient over the controle volume (centroid dual)
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects

    USE schism_msgp, ONLY : myrank,parallel_abort
    USE schism_glbl, ONLY : errmsg,iplg,rkind,dkind

!    USE VARS_WAVE, ONLY:MYID,MSR
!    USE MOD_PAR, ONLY: NGID,NGID_X,EGID,EGID_X
!    USE MOD_UTILS
!
    implicit none
!
!   Argument variables
!
    integer                , intent(in)   :: ivert  ! counter corresponding to current vertex
    !
    real(rkind), dimension(nverts), intent(in)   :: dep2   ! water depth at current time level
    real(rkind)                   , intent(out)  :: dhdx   ! derivative of depth in x-direction
    real(rkind)                   , intent(out)  :: dhdy   ! derivative of depth in y-direction
    real(rkind), dimension(MSC)   , intent(out)  :: dkdx   ! derivative of wave number in x-direction
    real(rkind), dimension(MSC)   , intent(out)  :: dkdy   ! derivative of wave number in y-direction
    real(rkind), dimension(nverts), intent(in)   :: mudl2  ! mud layer at current time level
    real(rkind), dimension(MSC)   , intent(in)   :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: icell    ! index of present cell
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: jc       ! loop counter
    integer                               :: jcell    ! index of next cell
    !
    integer, dimension(3)                 :: v        ! vertices in present cell
    !
    real(dkind)                                  :: area     ! twices the area of centroid dual around present vertex
    real(rkind), dimension(MSC)                  :: arr      ! auxiliary array
    real(rkind)                                  :: cslat    ! cosine of latitude
    real(rkind)                                  :: dpmax    ! maximum depth found in centroid dual
    real(rkind)                                  :: dpmin    ! minimum depth found in centroid dual
    real(rkind), dimension(3)                    :: dloc     ! local depth at vertices
    real(rkind), dimension(3)                    :: dm       ! local mud layer at vertices
    real(rkind)                                  :: dmaxc    ! maximum depth found per cell of centroid dual
    real(rkind)                                  :: dminc    ! minimum depth found per cell of centroid dual
    real(rkind), parameter                       :: drat= 5. ! ratio between maximum and minimum depths in centroid dual
    real(rkind)                                  :: h0       ! depth in centroid of present cell
    real(rkind)                                  :: h1       ! depth in centroid of next cell
    real(rkind), dimension(MSC)                  :: k0       ! wave number in centroid of present cell
    real(rkind), dimension(MSC)                  :: k1       ! wave number in centroid of next cell
    real(rkind), dimension(MSC,3)                :: kloc     ! local wave number at vertices
    real(dkind)                      :: x0       ! x-coordinate of the centroid of present cell
    real(dkind)                      :: x1       ! x-coordinate of the centroid of next cell
    real(dkind)                      :: y0       ! y-coordinate of the centroid of present cell
    real(dkind)                      :: y1       ! y-coordinate of the centroid of next cell
    !
    character(80)                         :: msgstr   ! string to pass message
    !
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes

    REAL(rkind) :: xmin,xmax

!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGradDepthorK')
    !
    dhdx = 0._rkind
    dhdy = 0._rkind
    !
    dkdx = 0._rkind
    dkdy = 0._rkind
    !
    ! if no frequency shift and no refraction, return
    !
    if ( (ITFRE == 0 .or. ICUR == 0) .and. IREFR == 0 ) return
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid

!    PRINT *,'SwanGradDepthorK BND Marker ?',vert(ivert)%atti(VMARKER)
!    PRINT *,'limiter on Ctheta',IREFR
!    PRINT *,'frequency shift  ',ITFRE
!    PRINT *,'PNUMS(32)        ',PNUMS(32)
!    PRINT *,'ICUR             ',ICUR

    !print*,'---- Grad Processing ivert ',ivert,'-------'

    !
    ! case gradients of depth
    !
    if ( (IREFR /= 0 .and. int(PNUMS(32)) == 0) .or. (ITFRE /= 0 .and. ICUR /= 0) ) then
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) return    ! boundary vertex
       !
       area  =  0_dkind
       dpmax = -99999.
       dpmin =  99999.

       xmax = 0._rkind
       xmin = huge(1)
       !
       ! loop over cells around considered vertex
!       PRINT *,'SwanGradDepthorK loop over cells around considered vertex'
       !
       do jc = 1, vert(ivert)%noc
          !
          ! get present cell and its vertices
          !
          icell = vert(ivert)%cell(jc)%atti(CELLID)
          !print*,'Grad Processing cell',icell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)

          !WRITE(PRTEST,*) 'CPU:',myrank,'Neigh Nodes :',iplg(v(1)),iplg(v(2)),iplg(v(3))
          !WRITE(PRTEST,*) 'CPU:',myrank,'Bath values :',dep2(v(1)),dep2(v(2)),dep2(v(3))
          !
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))

!          PRINT *,'Neigh Depths :',dloc(1),dloc(2),dloc(3)
          !
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
          !
          dminc = min ( min( dloc(1), dloc(2) ), dloc(3) )
          if ( dminc < dpmin ) dpmin = dminc
          !
          dmaxc = max ( max( dloc(1), dloc(2) ), dloc(3) )
          if ( dmaxc > dpmax ) dpmax = dmaxc
          !
          ! determine centroid of present cell
          !
          x0 = MyREAL(cell(icell)%attr(CELLCX))
          y0 = MyREAL(cell(icell)%attr(CELLCY))
          if(x0<xmin) xmin = x0
          if(x0>xmax) xmax = x0
          !
          ! determine depth in centroid in present cell
          !
          h0 = ( dloc(1) + dloc(2) + dloc(3) )/ 3.
          !
          ! get next cell in counterclockwise direction
          !
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          !
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          !
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          !
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
          !
          ! determine centroid of next cell
          !
          x1 = MyREAL(cell(jcell)%attr(CELLCX))
          y1 = MyREAL(cell(jcell)%attr(CELLCY))
          if(x1<xmin) xmin = x1          
          if(x1>xmax) xmax = x1

          ! Not crossing the dateline ?
          !IF(KSPHER.EQ.1.AND.abs(x1-x0)>180.0) THEN
          ! x0 = mod(x0+360.0,360.0)
          ! x1 = mod(x1+360.0,360.0)
          !ENDIF
          !IF(KSPHER.EQ.1) THEN !.and.(xmax-xmin).gt.180.0) THEN
          !  x0 = mod(x0+360.0,360.0)
          !  x1 = mod(x1+360.0,360.0)
          !ENDIF
          !
          ! determine depth in centroid of next cell
          !
          h1 = ( dloc(1) + dloc(2) + dloc(3) )/ 3.
          !
          ! compute contribution to area of centroid dual
          !
          area = area + x0*y1 - x1*y0
          !print*,'x0,y0,x1,y1,area',x0,y0,x1,y1,area
          !print*,'x0,y0',cell(icell)%attr(CELLCX),cell(icell)%attr(CELLCY)
          !
          ! compute contribution to x-gradient of depth
          !
          dhdx = dhdx + ( h0 + h1 ) * real( y1 - y0 )
          !
          ! compute contribution to y-gradient of depth
          !
          dhdy = dhdy + ( h0 + h1 ) * real( x1 - x0 )
          !
       enddo
#if 1
       ! Case Crossing Dateline: in centroid dual, some elements are on both sides of the dateline
       ! Redo the computation
       IF(KSPHER.GT.1.and.(xmax-xmin).gt.180.0) THEN
        area  =  0_dkind
        dpmax = -99999.
        dpmin =  99999.

        dhdx = 0._rkind
        dhdy = 0._rkind

        ! loop over cells around considered vertex
        do jc = 1, vert(ivert)%noc
          ! get present cell and its vertices
          icell = vert(ivert)%cell(jc)%atti(CELLID)
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
          dminc = min ( min( dloc(1), dloc(2) ), dloc(3) )
          if ( dminc < dpmin ) dpmin = dminc
          dmaxc = max ( max( dloc(1), dloc(2) ), dloc(3) )
          if ( dmaxc > dpmax ) dpmax = dmaxc
          ! determine centroid of present cell
          x0 = MyREAL(cell(icell)%attr(CELLCX))
          y0 = MyREAL(cell(icell)%attr(CELLCY))
          x0 = mod(x0+360.0,360.0)
          ! determine depth in centroid in present cell
          h0 = ( dloc(1) + dloc(2) + dloc(3) )/ 3.
          ! get next cell in counterclockwise direction
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 10
          ! determine centroid of next cell
          x1 = MyREAL(cell(jcell)%attr(CELLCX))
          y1 = MyREAL(cell(jcell)%attr(CELLCY))
          x1 = mod(x1+360.0,360.0)
          ! determine depth in centroid of next cell
          h1 = ( dloc(1) + dloc(2) + dloc(3) )/ 3.
          ! compute contribution to area of centroid dual
          area = area + x0*y1 - x1*y0
          ! compute contribution to x-gradient of depth
          dhdx = dhdx + ( h0 + h1 ) * real( y1 - y0 )
          ! compute contribution to y-gradient of depth
          dhdy = dhdy + ( h0 + h1 ) * real( x1 - x0 )
        enddo
       ENDIF
#endif
       !
       ! if ratio between max and min depth is too large, set gradients to zero and skip to next vertex
       !
       !if ( dpmax > drat * dpmin ) goto 10
       !
       ! if area is non-positive, give error and go to next vertex
       !
       if ( .not. area > 0_dkind) then
!          write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
!          call msgerr( 2, trim(msgstr) )
!          return

           !WRITE(12,*) 'AREA CPU:',myrank,'Nodes :',iplg(ivert),area
           write (msgstr, '(a,i6)') '1 Wrong cell ordering around centroid dual : CW instead CCW ',iplg(ivert)
           !msgstr = 'Wrong cell ordering around centroid dual : CW instead CCW',iplg(ivert)
           WRITE(errmsg,*) MSGSTR
           CALL parallel_abort(errmsg)
       endif
       !
       dhdx =  dhdx/real(area)
       dhdy = -dhdy/real(area)
       !
       ! in case of spherical coordinates, transform back to Cartesian coordinates
       !
       if ( KSPHER.GT.0 ) then
          !
          cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
          !
          dhdx = dhdx/(cslat * LENDEG)
          dhdy = dhdy/LENDEG
          !
       endif
       !
       return
 10    dhdx = 0.
       dhdy = 0.
       !
    endif
    !
    ! case gradients of wave number
    !
    if ( IREFR /= 0 .and. int(PNUMS(32)) == 1 ) then
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) return    ! boundary vertex
       !
       area  =  0d0
       dpmax = -99999.
       dpmin =  99999.

       xmax = 0._rkind
       xmin = huge(1)

       !
       ! loop over cells around considered vertex
!       PRINT *,'SwanGradDepthorK loop over cells around considered vertex'
!       PRINT *,'SwanGradDepthorK ivert', ivert
       !
       do jc = 1, vert(ivert)%noc
          !
          ! get present cell and its vertices
          !
          icell = vert(ivert)%cell(jc)%atti(CELLID)

!          PRINT *,'Neigh Cell ',jc,icell
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) THEN
!            write(12,*) IPLG(ivert),IPLG(v(1)),IPLG(v(2)),IPLG(v(3))
!          ENDIF

          ! dbg, case 2 mpi threads
          !IF(MYID.EQ.2) WRITE(PRINTF,*) 'CPU:',MYID,'Neigh Nodes :',NGID_X(v(1)),NGID_X(v(2)),NGID_X(v(3))
          !IF(MSR) WRITE(PRTEST,*) 'CPU:',MYID,'Neigh Nodes:',NGID_X(v(1)),NGID_X(v(2)),NGID_X(v(3))
          !
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          !
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
          !
          dminc = min ( min( dloc(1), dloc(2) ), dloc(3) )
          if ( dminc < dpmin ) dpmin = dminc
          !
          dmaxc = max ( max( dloc(1), dloc(2) ), dloc(3) )
          if ( dmaxc > dpmax ) dpmax = dmaxc
          !
          ! determine centroid of present cell
          !
          x0 = cell(icell)%attr(CELLCX)
          y0 = cell(icell)%attr(CELLCY)
          if(x0<xmin) xmin = x0
          if(x0>xmax) xmax = x0
          !
          ! compute wave numbers for all frequencies
          !
!          PRINT *,'compute wave numbers for all frequencies'

          call KSCIP1 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
          !
          if ( IMUD == 1 ) then
             !
             if (VARMUD) then
                dm(1) = mudl2(v(1))
                dm(2) = mudl2(v(2))
                dm(3) = mudl2(v(3))
             else
                dm(1) = PMUD(1)
                dm(2) = PMUD(1)
                dm(3) = PMUD(1)
             endif
             !
             call KSCIP2 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr, arr, dm(1))
             call KSCIP2 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr, arr, dm(2))
             call KSCIP2 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr, arr, dm(3))
             !
          endif
          !
          ! determine wave number in centroid in present cell
!          PRINT *,'determine wave number in centroid in present cell'
          !
          k0(:) = ( kloc(:,1) + kloc(:,2) + kloc(:,3) )/ 3.
          !
          ! get next cell in counterclockwise direction
!          PRINT *,'get next cell in counterclockwise direction'
          !
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)

          ! dbg, case 2 mpi threads
          !IF(MYID.EQ.2) WRITE(PRINTF,*) 'CPU:',MYID,'NEXTCELL',EGID_X(jcell)
          !IF(MSR) WRITE(PRTEST,*) 'CPU:',MYID,'NEXTCELL',EGID_X(jcell)

!          PRINT *,'NEXTCELL',jcell
          !
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)

!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) THEN
!            write(12,*) IPLG(ivert),IPLG(v(1)),IPLG(v(2)),IPLG(v(3))
!          ENDIF

          !
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          !
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
          !
          ! determine centroid of next cell
          !
          x1 = cell(jcell)%attr(CELLCX)
          y1 = cell(jcell)%attr(CELLCY)
          if(x1<xmin) xmin = x1
          if(x1>xmax) xmax = x1
          !
          ! compute wave numbers for all frequencies
          !
          call KSCIP1 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
          !
          if ( IMUD == 1 ) then
             !
             if (VARMUD) then
                dm(1) = mudl2(v(1))
                dm(2) = mudl2(v(2))
                dm(3) = mudl2(v(3))
             endif
             !
             call KSCIP2 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr, arr, dm(1))
             call KSCIP2 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr, arr, dm(2))
             call KSCIP2 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr, arr, dm(3))
             !
          endif
          !
          ! determine wave number in centroid of next cell
          !
          k1(:) = ( kloc(:,1) + kloc(:,2) + kloc(:,3) )/ 3.
          !
          ! Not crossing the dateline ?
          !IF(KSPHER.EQ.1) THEN !and.(xmax-xmin).gt.180.0) THEN
          ! x0 = mod(x0+360.0,360.0)
          ! x1 = mod(x1+360.0,360.0)
          !ENDIF
          !
          ! compute contribution to area of centroid dual
          !
          area = area + x0*y1 - x1*y0
!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) write(12,*) '0',IPLG(ivert),area
          !
          ! compute contribution to x-gradient of wave number
          !
          dkdx(:) = dkdx(:) + ( k0(:) + k1(:) ) * real( y1 - y0 )
          !
          ! compute contribution to y-gradient of wave number
          !
          dkdy(:) = dkdy(:) + ( k0(:) + k1(:) ) * real( x1 - x0 )
          !
!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) THEN 
!             write(12,*) '0a',SUM(dkdx(:)),SUM(dkdy(:))
!             write(12,*) '0b',x0,x1
!             IF(jc.eq.vert(ivert)%noc) write(12,*) '0a Done',IPLG(ivert),SUM(dkdx(:)),SUM(dkdy(:)),xmin,xmax
!          ENDIF 
       enddo

#if 1
       ! Case Crossing Dateline: in centroid dual, some elements are on both sides of the dateline
       ! Redo the computation 
       IF(KSPHER.GT.0.and.(xmax-xmin).gt.180.0) THEN
        area  =  0_dkind
        dpmax = -99999.
        dpmin =  99999.
        !
        dkdx = 0._rkind
        dkdy = 0._rkind

        ! loop over cells around considered vertex
        do jc = 1, vert(ivert)%noc
          ! get present cell and its vertices
          icell = vert(ivert)%cell(jc)%atti(CELLID)
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
          dminc = min ( min( dloc(1), dloc(2) ), dloc(3) )
          if ( dminc < dpmin ) dpmin = dminc
          dmaxc = max ( max( dloc(1), dloc(2) ), dloc(3) )
          if ( dmaxc > dpmax ) dpmax = dmaxc
          ! determine centroid of present cell
          x0 = cell(icell)%attr(CELLCX)
          y0 = cell(icell)%attr(CELLCY)
          x0 = mod(x0+360.0,360.0)
          ! compute wave numbers for all frequencies
          call KSCIP1 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
          if ( IMUD == 1 ) then
             if (VARMUD) then
                dm(1) = mudl2(v(1))
                dm(2) = mudl2(v(2))
                dm(3) = mudl2(v(3))
             else
                dm(1) = PMUD(1)
                dm(2) = PMUD(1)
                dm(3) = PMUD(1)
             endif
             !
             call KSCIP2 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr, arr, dm(1))
             call KSCIP2 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr, arr, dm(2))
             call KSCIP2 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr, arr, dm(3))
             !
          endif
          ! determine wave number in centroid in present cell
          k0(:) = ( kloc(:,1) + kloc(:,2) + kloc(:,3) )/ 3.
          ! get next cell in counterclockwise direction
          jcell = vert(ivert)%cell(jc)%atti(NEXTCELL)
          v(1) = cell(jcell)%atti(CELLV1)
          v(2) = cell(jcell)%atti(CELLV2)
          v(3) = cell(jcell)%atti(CELLV3)
          dloc(1) = dep2(v(1))
          dloc(2) = dep2(v(2))
          dloc(3) = dep2(v(3))
          if ( dloc(1) <= DEPMIN .or. dloc(2) <= DEPMIN .or. dloc(3) <= DEPMIN ) goto 20
          ! determine centroid of next cell
          x1 = cell(jcell)%attr(CELLCX)
          y1 = cell(jcell)%attr(CELLCY)
          x1 = mod(x1+360.0,360.0)
          ! compute wave numbers for all frequencies
          call KSCIP1 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr)
          call KSCIP1 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr)
          if ( IMUD == 1 ) then
             if (VARMUD) then
                dm(1) = mudl2(v(1))
                dm(2) = mudl2(v(2))
                dm(3) = mudl2(v(3))
             endif
             call KSCIP2 (MSC, spcsig, dloc(1), kloc(1,1), arr, arr, arr, arr, dm(1))
             call KSCIP2 (MSC, spcsig, dloc(2), kloc(1,2), arr, arr, arr, arr, dm(2))
             call KSCIP2 (MSC, spcsig, dloc(3), kloc(1,3), arr, arr, arr, arr, dm(3))
          endif
          ! determine wave number in centroid of next cell
          k1(:) = ( kloc(:,1) + kloc(:,2) + kloc(:,3) )/ 3.
          ! compute contribution to area of centroid dual
          area = area + x0*y1 - x1*y0
!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) write(12,*) '1',IPLG(ivert),area
          ! compute contribution to x-gradient of wave number
          dkdx(:) = dkdx(:) + ( k0(:) + k1(:) ) * real( y1 - y0 )
          ! compute contribution to y-gradient of wave number
          dkdy(:) = dkdy(:) + ( k0(:) + k1(:) ) * real( x1 - x0 )
!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) THEN
!             write(12,*) '1a',SUM(dkdx(:)),SUM(dkdy(:))
!             IF(jc.eq.vert(ivert)%noc) write(12,*) '1a Done',IPLG(ivert),SUM(dkdx(:)),SUM(dkdy(:))
!          ENDIF
       enddo
      ENDIF
#endif

       !
       ! if ratio between max and min depth is too large, set gradients to zero and skip to next vertex
       !
       !if ( dpmax > drat * dpmin ) goto 20
       !
       ! if area is non-positive, give error and go to next vertex

       !
       if ( .not. area > 0d0 ) then
!          write (msgstr, '(a,i5)') ' Area of centroid dual is negative or zero in vertex ', ivert
!          call msgerr( 2, trim(msgstr) )
!          return
           !msgstr = 'Wrong cell ordering around centroid dual : CW instead CCW'
           write (msgstr, '(a,i6)') '2 Wrong cell ordering around centroid dual : CW instead CCW ',iplg(ivert)
           WRITE(errmsg,*) MSGSTR
           CALL parallel_abort(errmsg)
       endif
       !
       dkdx(:) =  dkdx(:)/real(area)
       dkdy(:) = -dkdy(:)/real(area)
       !
       ! in case of spherical coordinates, transform back to Cartesian coordinates
       !
       if ( KSPHER.GT.0 ) then
          !
          cslat = cos(DEGRAD*(vert(ivert)%attr(VERTY) + YOFFS))
          !
          dkdx(:) = dkdx(:)/(cslat * LENDEG)
          dkdy(:) = dkdy(:)/LENDEG


!          IF(IPLG(ivert).eq.91.or.IPLG(ivert).eq.26048.or.IPLG(ivert).eq.52005.or.IPLG(ivert).eq.5646) THEN
!             write(12,*) '1a Done',IPLG(ivert),SUM(dkdx(:)),SUM(dkdy(:))
!          ENDIF

          !
       endif
       !
       return
 20    dkdx = 0.
       dkdy = 0.
       !
    endif
    !
end subroutine SwanGradDepthorK
