subroutine SwanGridCell ( ncells, nverts, xcugrd, ycugrd, kvertc )
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
!
!   Updates
!
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   Fills cell-based data structure
!
!   Method
!
!   Based on unstructured grid
!   Note: we restrict ourselves to triangles only!
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use swcomm4, ONLY : KSPHER
    use SwanGridobjects

    USE schism_glbl, ONLY : rkind
    USE schism_msgp, ONLY:parallel_abort

!    USE ALL_VARS, ONLY: NTVE,NBVE
!    USE MOD_SPHERICAL
!    USE MOD_UTILS
!    USE MOD_PREC
!
    implicit none
!
!   Argument variables
!
    integer, intent(in)                       :: ncells  ! number of cells in grid
    integer, intent(in)                       :: nverts  ! number of vertices in grid
    !
    integer, dimension(3, ncells), intent(in) :: kvertc  ! vertices of the cell
                                                         ! (must be filled by a gridgenerator!)
    !
    real(rkind), dimension(nverts), intent(in)       :: xcugrd  ! the x-coordinates of the grid vertices
    real(rkind), dimension(nverts), intent(in)       :: ycugrd  ! the y-coordinates of the grid vertices
!
!   Local variables
!
    integer                               :: I1,J1,J2
    integer                               :: icell    ! loop counter over cells / index of present cell
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: ivert    ! loop counter over vertices / index of present vertex
    integer                               :: j        ! loop counter
    integer                               :: jcell    ! index of next cell
    integer                               :: k        ! auxiliary integer / loop counter
    integer                               :: l        ! loop counter
    integer                               :: noc      ! number of cells around considered vertex
    integer, parameter                    :: nov = 3  ! number of vertices in present cell (triangles only)
    integer, dimension(3)                 :: v        ! vertices in present cell
    integer                               :: vc       ! considered vertex
    integer                               :: vn       ! first upwave vertex of next cell
    integer                               :: vp       ! last upwave vertex of present cell
    !
    integer, dimension(:), allocatable    :: ivlist   ! list of index vertices
    !
    real(rkind)                                  :: carea    ! area of the present cell (triangles only)
    real(rkind)                                  :: dx1      ! first component of covariant base vector a_(1)
    real(rkind)                                  :: dx2      ! second component of covariant base vector a_(1)
    real(rkind)                                  :: dy1      ! first component of covariant base vector a_(2)
    real(rkind)                                  :: dy2      ! second component of covariant base vector a_(2)
    real(rkind)                                  :: rdet     ! reciproke of determinant
    real(rkind)                                  :: th1      ! direction of one face pointing to present vertex
    real(rkind)                                  :: th2      ! direction of another face pointing to present vertex
    real(rkind)                                  :: thdiff   ! difference between th2 and th1
    real(rkind), dimension(2)                    :: vec12    ! translation vector of coordinates: vertex2 - vertex1
    real(rkind), dimension(2)                    :: vec13    ! translation vector of coordinates: vertex3 - vertex1
    real(rkind)                                  :: xc       ! x-coordinate of the cell-centroid
    real(rkind)                                  :: yc       ! y-coordinate of the cell-centroid
    !
    logical                               :: nxtcell  ! indicate whether there is next cell in counterclockwise direction
    !
!#  if defined (SPHERICAL)
!    REAL :: side1,side2,side3
!#  endif
    REAL(rkind) :: VX1,VX2,VX3, VY1,VY2,VY3, area1
    LOGICAL :: REVERSE = .FALSE.
    REAL(rkind) :: xmin,xmax


    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
    type(celltype), dimension(:), pointer :: cell     ! datastructure for cells with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanGridCell')

    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid

    ! 
! 03/2017 JLEFEVRE Update : In the original SWAN subroutine "swangridcell"
! assume both a CCW counting of vertices in triangle and a CCW counting of
! triangles around vertex forming the tracer control element
! (also called 'dual triangle' in SWAN)
! This CCW counting is also applied when swan reads an ADCIRC mesh
!
! But in FVCOM grid, surrounding triangles are counted CW and vertice ordering
! in triangle is CW !! SWAN want a CCW ordering : this grid topology is used
! for example in "swangraddepthork"
!
! We test the vertice ordering by computing the cell area in the first cell, 
! If area gt 0 , REVERSE = .FALSE. means vertice ordering is CCW, nothing to do
! If area lt 0 , REVERSE = .TRUE.  means vertice ordering is CW and 
!                            We reverse the vertex ordering
    !
    area1 = 0._rkind
#if 0
    icell = 1
    v(1) = kvertc(1,icell)
    v(2) = kvertc(2,icell)
    v(3) = kvertc(3,icell)
    VX1= xcugrd(v(1))
    VX2= xcugrd(v(2))
    VX3= xcugrd(v(3))
    VY1= ycugrd(v(1))
    VY2= ycugrd(v(2))
    VY3= ycugrd(v(3))
    !PRINT *,v(1),v(2),v(3)
    !PRINT *,'VX1=',VX1,'VX2=',VX2,'VX3=',VX3
    !PRINT *,'VY1=',VY1,'VY2=',VY2,'VY3=',VY3

    area1 = (VX2 - VX1) * (VY3 - VY1) - &
            (VX3 - VX1) * (VY2 - VY1)
    area1 = 0.5_rkind * area1
    !PRINT *,'area1=',area1
    IF(area1.GT.0._rkind) THEN
     !PRINT *,'POSITIVE'
     REVERSE = .FALSE.
    ELSE
     PRINT *,'area1=',area1
     PRINT *,'NEGATIVE or 0'
     CALL parallel_abort('swangridcell:NEGATIVE')
     REVERSE = .TRUE.
    ENDIF
    !PRINT *,'area1 ',area1
#endif

    !
    ! loop over all cells
    !
    do icell = 1, ncells

#if 1
       area1 = 0._rkind
       v(1) = kvertc(1,icell)
       v(2) = kvertc(2,icell)
       v(3) = kvertc(3,icell)

       VX1= xcugrd(v(1)); VX2= xcugrd(v(2)); VX3= xcugrd(v(3))
       VY1= ycugrd(v(1)); VY2= ycugrd(v(2)); VY3= ycugrd(v(3)) 

       IF(KSPHER.EQ.1)THEN
        IF(abs(xcugrd(v(1))-xcugrd(v(2)))>180.0.OR.&
           abs(xcugrd(v(3))-xcugrd(v(1)))>180.0.OR.&
           abs(xcugrd(v(3))-xcugrd(v(2)))>180.0)THEN
           VX1=mod(VX1+360.0,360.0) 
           VX2=mod(VX2+360.0,360.0)
           VX3=mod(VX3+360.0,360.0)
        ENDIF
       ENDIF

       area1 = (VX2 - VX1) * (VY3 - VY1) - &
               (VX3 - VX1) * (VY2 - VY1)
       area1 = 0.5_rkind * area1       

       !PRINT *,'area1,icell=',area1,icell

       IF(area1.GT.0._rkind) THEN
        !PRINT *,'POSITIVE'
        !REVERSE = .FALSE.
       ELSE
        !PRINT *,'VX1,VX2,VX3',VX1,VX2,VX3
        !PRINT *,'VY1,VY2,VY3',VY1,VY2,VY3
        PRINT *,'NEGATIVE or 0'
        CALL parallel_abort('swangridcell:NEGATIVE')
        !REVERSE = .TRUE.
       ENDIF
#endif

       !
       ! determine number of vertices and faces (triangles only!)
       !
       cell(icell)%nov = nov
       cell(icell)%nof = nov
       !
       ! identification number
       !
       cell(icell)%atti(CELLID) = icell
       !
       ! cell is triangle and initiatively active
       !
       cell(icell)%atti(CELLRECT) = 0
       cell(icell)%active         = .true.
       !
       ! store vertices of the cell
       !
!      Original
       v(1) = kvertc(1,icell)
       v(2) = kvertc(2,icell)
       v(3) = kvertc(3,icell)
       !  
       ! JL In case of wrong vertice ordering ...
       ! We reverse the vertex ordering now
       !
       IF(REVERSE) THEN
        v(1) = kvertc(1,icell)
        v(2) = kvertc(3,icell)
        v(3) = kvertc(2,icell)
       ENDIF
 
       cell(icell)%atti(CELLV1) = v(1)
       cell(icell)%atti(CELLV2) = v(2)
       cell(icell)%atti(CELLV3) = v(3)
       !
       IF(KSPHER.EQ.1)THEN

        IF(abs(xcugrd(v(1))-xcugrd(v(2)))>180.0.OR.&
           abs(xcugrd(v(3))-xcugrd(v(1)))>180.0.OR.&
           abs(xcugrd(v(3))-xcugrd(v(2)))>180.0)THEN
           ! JL : SWAN Fix, Crossing the date line ...
           vec12(1) = mod(xcugrd(v(2))+360.0,360.0) - mod(xcugrd(v(1))+360.0,360.0)
           vec13(1) = mod(xcugrd(v(3))+360.0,360.0) - mod(xcugrd(v(1))+360.0,360.0)
           vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
           vec13(2) = ycugrd(v(3)) - ycugrd(v(1))
        ELSE
           vec12(1) = xcugrd(v(2)) - xcugrd(v(1))
           vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
           vec13(1) = xcugrd(v(3)) - xcugrd(v(1))
           vec13(2) = ycugrd(v(3)) - ycugrd(v(1))
         ! print*,'vec12,vec13',vec12(1),vec12(2),vec13(1),vec13(2)
        ENDIF

       ELSE
        vec12(1) = xcugrd(v(2)) - xcugrd(v(1))
        vec12(2) = ycugrd(v(2)) - ycugrd(v(1))
        vec13(1) = xcugrd(v(3)) - xcugrd(v(1))
        vec13(2) = ycugrd(v(3)) - ycugrd(v(1))
       ENDIF

       !
       ! store area of the cell in carea
       !
       carea = 0.5_rkind*abs(vec12(1)*vec13(2) - vec13(1)*vec12(2))
       !
       ! store local covariant and contravariant base vectors at each vertex
       ! next, store directions of faces pointing to present vertex
       !
       dx1 = -vec12(1)
       dy1 = -vec12(2)
       dx2 = -vec13(1)
       dy2 = -vec13(2)

!       if(icell.EQ.2) THEN
!        print*,'vec12(1)=',vec12(1)
!        print*,'vec12(2)=',vec12(2)
!        print*,'vec13(1)=',vec13(1)
!        print*,'vec13(2)=',vec13(2)
!        print*,'carea=',carea
!        print*,'dx1=',dx1,'dy1=',dy1
!        print*,'dx2=',dx2,'dy2=',dy2
!       endif

       !
       do j = 1, nov
          !
          rdet = 1./(dy2*dx1 - dy1*dx2)
          !
          cell(icell)%geom(j)%det  =  dy2*dx1 - dy1*dx2
          cell(icell)%geom(j)%dx1  =  dx1
          cell(icell)%geom(j)%dy1  =  dy1
          cell(icell)%geom(j)%dx2  =  dx2
          cell(icell)%geom(j)%dy2  =  dy2
          cell(icell)%geom(j)%rdx1 =  dy2*rdet
          cell(icell)%geom(j)%rdy1 = -dx2*rdet
          cell(icell)%geom(j)%rdx2 = -dy1*rdet
          cell(icell)%geom(j)%rdy2 =  dx1*rdet
          !
          th1 = atan2(dy1,dx1)
          th2 = atan2(dy2,dx2)
          !
          thdiff = th1 - th2
          do
             if ( abs(thdiff) <= PI_W ) exit
             th1 = th1 - sign (2._rkind, thdiff) * PI_W
             thdiff = th1 - th2
          enddo
          !
          cell(icell)%geom(j)%th1 = th1
          cell(icell)%geom(j)%th2 = th2
          !
          dx1 = dx2 - dx1
          dy1 = dy2 - dy1
          dx2 = dx1 - dx2
          dy2 = dy1 - dy2
          !
       enddo
       !
       ! determine orientation of the mesh
       !
       if ( icell == 1 ) then
          !
          if ( dy2*dx1 > dy1*dx2 ) then
             CVLEFT = .false.             ! right-handed
          else
             CVLEFT = .true.              ! left-handed
          endif
          !
       endif
       !
       ! store coordinates of centroid
       !
       xc = 0._rkind
       yc = 0._rkind
       xmax = 0._rkind
       xmin = huge(1)
       do j = 1, nov
          if(xcugrd(kvertc(j,icell))<xmin) xmin = xcugrd(kvertc(j,icell))
          if(xcugrd(kvertc(j,icell))>xmax) xmax = xcugrd(kvertc(j,icell))
          xc = xc + xcugrd(kvertc(j,icell))
          yc = yc + ycugrd(kvertc(j,icell))
       enddo
       xc = xc / real(nov)
       yc = yc / real(nov)

       ! Case Crossing Dateline
      IF(KSPHER.EQ.1.and.(xmax-xmin).gt.180.0) THEN
        xc = 0._rkind
        yc = 0._rkind
        do j = 1, nov
          xc = xc + mod(xcugrd(kvertc(j,icell))+360.0,360.0)
          yc = yc + ycugrd(kvertc(j,icell))
        enddo
        xc = xc / real(nov)
        yc = yc / real(nov)
      ENDIF

       !
       cell(icell)%attr(CELLAREA) = carea
       cell(icell)%attr(CELLCX  ) = xc
       cell(icell)%attr(CELLCY  ) = yc

!       if(icell.EQ.2) THEN
!        print*,'CELLCX=',xc
!        print*,'CELLCY=',yc
!       endif

       !
       ! store vertices of each face of the cell
       !
       do j = 1, nov
          k = mod(j,nov)+1
          cell(icell)%face(j)%atti(FACEV1) = kvertc(j,icell)
          cell(icell)%face(j)%atti(FACEV2) = kvertc(k,icell)
       enddo
       !
    enddo
    !
    ! Build the list of Cells forming the Dual triangle
    !
    allocate(ivlist(nverts))
    !
    ! loop over all vertices
    !
    do ivert = 1, nverts
       !
       ! identify the considered vertex and store index
       !
       vc = vert(ivert)%atti(VERTID)
       ivlist(vc) = ivert
       !
       ! initialize number of cells around vertex
       !
       vert(ivert)%noc = 0
       !
    enddo
    !
    do icell = 1, ncells
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! add present cell to each of these vertices
       !
      do j = 1, nov
         !
         ivert = ivlist(v(j))
         noc   = vert(ivert)%noc +1
         !
         vert(ivert)%noc = noc
         vert(ivert)%cell(noc)%atti(CELLID) = icell
         !
      enddo
      !
    enddo
    !
    deallocate(ivlist)

    ! Here we look for NEXTCELL

    do ivert = 1, nverts
       !
       noc = vert(ivert)%noc
       !
       ! loop over cells around considered vertex
       !
       do j = 1, noc
          !
          ! get cell and its vertices
          !
          icell = vert(ivert)%cell(j)%atti(CELLID)
          !
          v(1) = cell(icell)%atti(CELLV1)
          v(2) = cell(icell)%atti(CELLV2)
          v(3) = cell(icell)%atti(CELLV3)
          !
          ! pick up last upwave vertex (counterclockwise counting of vertices is
          ! assumed)
          !
          do k = 1, nov
             if ( v(k) == ivert ) then
                vp = v(mod(k+1,nov)+1)
                exit
             endif
          enddo
          !
          ! search for next cell in counterclockwise direction
          !
          nxtcell = .false.
          !
          do l = 1, noc
     
             if(l==j) cycle
             !
             ! get a cell and its vertices
             !
             jcell = vert(ivert)%cell(l)%atti(CELLID)
             !
             v(1) = cell(jcell)%atti(CELLV1)
             v(2) = cell(jcell)%atti(CELLV2)
             v(3) = cell(jcell)%atti(CELLV3)
             !
             ! pick up first upwave vertex (counterclockwise counting of
             ! vertices is assumed)
             !
             do k = 1, nov
                if ( v(k) == ivert ) then
                   vn = v(mod(k,nov)+1)
                   exit
                endif
             enddo
             !
             ! check whether first upwave vertex of next cell equals last upwave
             ! vertex of present cell
             !
             if ( vn == vp ) then
                nxtcell = .true.
                exit
             endif
             !
          enddo
          !
          if ( .not.nxtcell ) jcell = 0
          vert(ivert)%cell(j)%atti(NEXTCELL) = jcell
          !
       enddo
       !
    enddo

end subroutine SwanGridCell
