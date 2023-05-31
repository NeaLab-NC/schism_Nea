subroutine SwanVertlist
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
!   40.80: Marcel Zijlema
!   41.07: Casey Dietrich
!   41.48: Marcel Zijlema
!
!   Updates
!
!   40.80,  July 2007: New subroutine
!   41.07,  July 2009: small fix (assign ref.point to deepest point in case of no b.c.)
!   41.48, March 2013: including order along a user-given direction
!
!   Purpose
!
!   Makes vertex list with specific order
!
!   Method
!
!   Sorting based on distance with respect to a reference point
!   This reference point can be either a vertex with boundary condition or a deepest point
!
!   Modules used
!
    USE OCPCOMM4
    USE M_GENARR, only: DEPTH
    USE SWANGRIDDATA
    USE SWANGRIDOBJECTS
    USE SWANCOMPDATA

    USE schism_glbl,ONLY:rkind,iplg,np,npa
    USE schism_msgp,ONLY:myrank
    USE VARS_WAVE, ONLY:SERIAL

    USE SWCOMM4, ONLY:LENDEG,KSPHER
!
    implicit none
!
!   Local variables
!
    integer, save                   :: ient = 0   ! number of entries in this subroutine
    integer                         :: istat      ! indicate status of allocation
    integer                         :: itmp       ! temporary stored integer for swapping
    integer                         :: j          ! loop counter over vertices
    integer                         :: k          ! counter
    integer, dimension(1)           :: kd         ! location of minimum value in array dist
    integer, dimension(1)           :: kx         ! location of minimum value in array of x-coordinates of boundary vertices
    integer, dimension(1)           :: ky         ! location of minimum value in array of y-coordinates of boundary vertices
    !
    real(rkind)                            :: d1         ! distance of a point to origin
    real(rkind)                            :: d2         ! distance of another point to origin
    real(rkind)                            :: depmax     ! maximum depth found
    real(rkind)                            :: rtmp       ! temporary stored real(rkind) for swapping
    real(rkind)                            :: x0         ! x-coordinate of reference point
    real(rkind)                            :: y0         ! y-coordinate of reference point
    real(rkind)                            :: x1         ! x-coordinate of point
    real(rkind)                            :: y1         ! y-coordinate of point

    !
    real(rkind), dimension(:), allocatable :: dist       ! distance of each point with respect to reference point
    !
    type(verttype), dimension(:), pointer :: vert ! datastructure for vertices with their attributes

    character(len=80) :: FILEDBG
    integer :: ILPOS,I
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanVertlist')

    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    istat = 0
    if(.not.allocated(vlist)) allocate (vlist(nverts), stat = istat)
    if ( istat /= 0 ) then
       call msgerr ( 4, 'Allocation problem in SwanVertlist: array vlist ' )
       return
    endif
    !
    allocate (dist(nverts))

!    print*,'ASORT = ',asort

    !
    if ( asort /= -999. ) then
       !
       ! order direction asort
       !
       do j = 1, nverts
          dist(j) = vert(j)%attr(VERTX) * cos(asort) + vert(j)%attr(VERTY) * sin(asort)
       enddo
       !
    else
       !
       ! determine reference point
       !
       if ( all(mask=vert(:)%atti(VBC)==0) ) then
          !
          ! if no boundary condition is given then find the vertex with the maximum depth
          !
          depmax = -999.
          !
          do j = 1, nverts
             !
             if ( DEPTH(j) > depmax ) then
                !
                depmax = DEPTH(j)
                kx(1)  = j
                ky(1)  = j
                !
             endif
             !
          enddo
          !
       !JL: addon : start from open OBC unlike island ...
!       elseif(any(mask=vmark(:)==2)) then
!
!          kx = minloc(vert(:)%attr(VERTX), vmark(:)==2,KIND=rkind)
!          ky = minloc(vert(:)%attr(VERTY), vmark(:)==2,KIND=rkind)

       else
          !
          ! reference point is one of the vertices nearest to the origin where boundary condition is given
          !
          !kx = minloc(vert(:)%attr(VERTX), vert(:)%atti(VBC)/=0,KIND=rkind)
          !ky = minloc(vert(:)%attr(VERTY), vert(:)%atti(VBC)/=0,KIND=rkind)
          kx = minloc(vert(:)%attr(VERTX), mask=vert(:)%atti(VBC)/=0)
          ky = minloc(vert(:)%attr(VERTY), mask=vert(:)%atti(VBC)/=0)
          !
       endif
       !
       if ( kx(1) == ky(1) ) then
          x0 = vert(kx(1))%attr(VERTX)
          y0 = vert(ky(1))%attr(VERTY)
       else
          !
          d1 = sqrt((vert(kx(1))%attr(VERTX))**2+(vert(kx(1))%attr(VERTY))**2)
          d2 = sqrt((vert(ky(1))%attr(VERTX))**2+(vert(ky(1))%attr(VERTY))**2)
          !
          if ( d1 < d2 ) then
             x0 = vert(kx(1))%attr(VERTX)
             y0 = vert(kx(1))%attr(VERTY)
          else
             x0 = vert(ky(1))%attr(VERTX)
             y0 = vert(ky(1))%attr(VERTY)
          endif
          !
       endif
       !
       ! calculate distance of each point with respect to reference point
       !

       do j = 1, nverts
          dist(j) = sqrt((vert(j)%attr(VERTX)-x0)**2 + (vert(j)%attr(VERTY)-y0)**2)
       enddo

       IF ( KSPHER.GT.0 ) THEN
!           transform to local Cartesian coordinates
            x0 = x0*LENDEG*COS(PI_W*y0/180._rkind)
            y0 = y0*LENDEG
            do j = 1, nverts
             x1 = vert(j)%attr(VERTX)*LENDEG*COS(PI_W*vert(j)%attr(VERTY)/180._rkind)
             y1 = vert(j)%attr(VERTY)*LENDEG
             dist(j) = sqrt((x1-x0)**2 + (y1-y0)**2)
            enddo
       ELSE

        do j = 1, nverts
          dist(j) = sqrt((vert(j)%attr(VERTX)-x0)**2 + (vert(j)%attr(VERTY)-y0)**2)
        enddo

       ENDIF 



       !
    endif

!DBG
#if 0
      IF(.NOT.SERIAL) THEN
         FILEDBG = 'dist_details_after_0000'
         ILPOS = LEN_TRIM(FILEDBG)
         WRITE(FILEDBG(ILPOS-3:ILPOS),'(i4.4)') myrank
         print*,FILEDBG
         open(38,file='outputs/'//TRIM(FILEDBG),status='replace')
         write(38,*)'np,npa=',np,npa
         write(38,*)'INDEX_GL,vmark,dist'
         DO I = 1,nverts
          write(38,*)iplg(I),vmark(I),dist(I)
         ENDDO
         CLOSE(38)
      ENDIF
#endif
!END DBG


    !
    ! sort vertex list in order of increasing distance
    !
    vlist=(/ (j, j=1, nverts) /)
    !
    do j = 1, nverts-1
       !
       !kd = minloc(dist(j:nverts),KIND=rkind)
       kd = minloc(dist(j:nverts))
       k  = kd(1) + j-1
       !
       if ( k /= j ) then
          !
          rtmp     = dist(j)
          dist(j)  = dist(k)
          dist(k)  = rtmp
          !
          itmp     = vlist(j)
          vlist(j) = vlist(k)
          vlist(k) = itmp
          !
       endif
       !
    enddo
    !
!DBG
#if 0
      IF(.NOT.SERIAL) THEN
         FILEDBG = 'vlist_details_after_0000'
         ILPOS = LEN_TRIM(FILEDBG)
         WRITE(FILEDBG(ILPOS-3:ILPOS),'(i4.4)') myrank
         print*,FILEDBG
         open(38,file='outputs/'//TRIM(FILEDBG),status='replace')
         write(38,*)'np,npa=',np,npa
         write(38,*)'dist,sorted,vmark'
         DO I = 1,nverts
          write(38,'(f10.2,2I7)')dist(I),iplg(vlist(I)),vmark(vlist(I))
         ENDDO
         CLOSE(38)
      ENDIF
#endif
!END DBG

    deallocate(dist)

end subroutine SwanVertlist
