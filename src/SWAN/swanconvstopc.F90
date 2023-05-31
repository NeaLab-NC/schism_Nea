#include "swan_functions.h"
!subroutine SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, &
!                    tmdifc, delhs, deltm, xytst, spcsig, ac2, ivlow, ivup )
subroutine SwanConvStopc ( accur, hscurr, hsprev, hsdifc, tmcurr, tmprev, &
                  tmdifc, delhs, deltm, xytst, spcsig,      ivlow, ivup )
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
!   40.93: Marcel Zijlema
!   41.10: Marcel Zijlema
!
!   Updates
!
!   40.80,   October 2007: New subroutine
!   40.93, September 2008: extended with curvature of mean period
!   41.10,    August 2009: parallelization using OpenMP directives
!
!   Purpose
!
!   Determine accuracy of wave height and period by means of curvature for convergence check
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use swcomm4
    use SwanGriddata
    use SwanGridobjects

    USE VARS_WAVE    !, ONLY : AC2
    USE schism_glbl, ONLY : rkind,dkind
    USE schism_msgp



!#  if defined (MULTIPROCESSOR)
!    USE MOD_PAR
!#  endif

!
    implicit none


!#ifndef USE_MPIMODULE
    include 'mpif.h'
!#endif

!
!   Argument variables
!
    integer, intent(in)                         :: ivlow  ! lower index in range of vertices in calling thread
    integer, intent(in)                         :: ivup   ! upper index in range of vertices in calling thread
    integer, dimension(NPTST), intent(in)       :: xytst  ! test points for output purposes
    !
    real(rkind), intent(inout)                         :: accur  ! percentage of active vertices in which required accuracy has been reached
!   real(rkind), dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time level
    real(rkind), dimension(nverts), intent(out)        :: delhs  ! difference in wave height between last 2 iterations in all vertices
    real(rkind), dimension(nverts), intent(out)        :: deltm  ! difference in mean period between last 2 iterations in all vertices
    real(rkind), dimension(nverts), intent(inout)      :: hscurr ! wave height at current iteration level
    real(rkind), dimension(nverts), intent(inout)      :: hsdifc ! difference in wave height of current and one before previous iteration
    real(rkind), dimension(nverts), intent(inout)      :: hsprev ! wave height at previous iteration level
    real(rkind), dimension(nverts), intent(inout)      :: tmcurr ! mean period at current iteration level
    real(rkind), dimension(nverts), intent(inout)      :: tmdifc ! difference in mean period of current and one before previous iteration
    real(rkind), dimension(nverts), intent(inout)      :: tmprev ! mean period at previous iteration level
    real(rkind), dimension(MSC), intent(in)            :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer                               :: id       ! loop counter over direction bins
    integer, save                         :: ient = 0 ! number of entries in this subroutine
    integer                               :: is       ! loop counter over frequency bins
    integer                               :: ivert    ! loop counter over vertices
    integer                               :: j        ! loop counter
    !
    real(rkind)                                  :: curvah   ! required accuracy with respect to curvature in wave height
    real(rkind)                                  :: curvat   ! required accuracy with respect to curvature in mean period
    real(rkind)                                  :: fact     ! auxiliary factor
    real(rkind)                                  :: hsabs    ! absolute difference in wave height between last 2 iterations
    real(rkind)                                  :: hscurv   ! curvature of iteration curve of wave height
    real(rkind)                                  :: hsdif0   ! value of hsdifc at previous iteration level
    real(rkind)                                  :: hsprev0  ! wave height at one before previous iteration level
    real(rkind)                                  :: hsrel    ! required accuracy with respect to relative error in wave height
    real(rkind)                                  :: m0       ! moment of zeroth order
    real(rkind)                                  :: m1       ! moment of first order
    integer                               :: npacc    ! number of vertices in which required accuracy has been reached
    integer                               :: npacct   ! npacc counter for calling thread
    integer                               :: nwetp    ! total number of active vertices
    integer                               :: nwetpt   ! nwetp counter for calling thread
    real(rkind)                                  :: tmabs    ! absolute difference in mean period between last 2 iterations
    real(rkind)                                  :: tmcurv   ! curvature of iteration curve of mean period
    real(rkind)                                  :: tmdif0   ! value of tmdifc at previous iteration level
    real(rkind)                                  :: tmprev0  ! mean period at one before previous iteration level
    real(rkind)                                  :: tmrel    ! required accuracy with respect to relative error in mean period
    !
!   For Parallel run :
    INTEGER, DIMENSION(2)  :: SBUF,RBUF1               ! For Sum/Max/Min of quantitie over partition

    logical                               :: lhead    ! logical indicating to write header
    logical                               :: tstfl    ! indicates whether vertex is a test point
    !
    type(verttype), dimension(:), pointer :: vert     ! datastructure for vertices with their attributes
!
!   Common variables
!
    common/convstopc/npacc,nwetp
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanConvStopc')

    !
    ! point to vertex object
    !
    vert => gridobject%vert_grid
    !
    npacc = 0
    nwetp = 0
    !
    npacct = 0
    nwetpt = 0
    !
    deltm = 0._rkind
    delhs = 0._rkind
    !
    lhead = .true.
    !
    ! calculate a set of accuracy parameters based on relative error and curvature for Hs and Tm
    !
    do ivert = ivlow, ivup
       !
       if ( vert(ivert)%active ) then
          !
          ! determine whether the present vertex is a test point
          !
          tstfl = .false.
          if ( NPTST > 0 ) then
             do j = 1, NPTST
                if ( ivert /= xytst(j) ) cycle
                tstfl = .true.
             enddo
          endif
          !
          ! count active points
          !
          nwetpt = nwetpt + 1
          !
          ! store wave height and mean period of previous iteration levels
          !
          hsprev0       = max( 1.e-20, hsprev(ivert) )
          hsprev(ivert) = max( 1.e-20, hscurr(ivert) )
          tmprev0       = max( 1.e-20, tmprev(ivert) )
          tmprev(ivert) = max( 1.e-20, tmcurr(ivert) )
          !
          ! compute wave height and mean period for present vertex
          !
          m0 = 0._rkind
          m1 = 0._rkind
          do is = 1, MSC
             do id = 1, MDC
                fact = spcsig(is)**2 * ac2(id,is,ivert)
                m0 = m0 + fact
                m1 = m1 + fact * spcsig(is)
             enddo
          enddo
          m0 = m0 * FRINTF * DDIR
          m1 = m1 * FRINTF * DDIR
          !
          if ( m0 > 0. ) then
             hscurr(ivert) = max ( 1.e-20, 4._rkind*sqrt(m0) )
          else
             hscurr(ivert) = 1.e-20
          endif
          if ( m1 > 0. ) then
             tmcurr(ivert) = max ( 1.e-20, PI2_W*(m0/m1) )
          else
             tmcurr(ivert) = 1.e-20
          endif
          !
          ! compute absolute differences in wave height and mean period between last 2 iterations
          !
          hsabs = abs ( hscurr(ivert) - hsprev(ivert) )
          tmabs = abs ( tmcurr(ivert) - tmprev(ivert) )
          !
          delhs(ivert) = hsabs
          deltm(ivert) = tmabs
          !
          ! compute curvature of wave height
          !
          hsdif0        = hsdifc(ivert)
          hsdifc(ivert) = 0.5*( hscurr(ivert) - hsprev0 )
          hscurv        = abs ( hsdifc(ivert) - hsdif0 )
          !
          ! compute curvature of mean period
          !
          tmdif0        = tmdifc(ivert)
          tmdifc(ivert) = 0.5*( tmcurr(ivert) - tmprev0 )
          tmcurv        = abs ( tmdifc(ivert) - tmdif0 )
          !
          ! compute required accuracies for wave height
          !
          hsrel  = PNUMS( 1) * hscurr(ivert)
          curvah = PNUMS(15) * hscurr(ivert)
          !
          ! compute required accuracies for mean period
          !
          tmrel  = PNUMS( 1) * tmcurr(ivert)
          curvat = PNUMS(16) * tmcurr(ivert)
          !
          ! count vertices where wave height and period have reached required accuracies
          !
!          if(myrank.EQ.0) THEN
!          WRITE(PRINTF,*) 'm0,m1',m0,m1 
!          WRITE(PRINTF,*) 'tmcurr(ivert),PNUMS(15),PNUMS(16)',tmcurr(ivert),PNUMS(15),PNUMS(16)
!
!          WRITE(PRINTF,*) 'hscurv,curvah,hsabs,max(hsrel,PNUMS(2))',&
!                          hscurv,curvah,hsabs,max(hsrel,PNUMS(2))
!          WRITE(PRINTF,*) 'tmcurv,curvat,tmabs,max(tmrel,PNUMS(3))',&
!                          tmcurv,curvat,tmabs,max(tmrel,PNUMS(3))
!          WRITE(PRINTF,*) 'npacct',npacct
!          WRITE(PRINTF,*) 'nwetpt',nwetpt
!          endif

          if ( (hscurv <= curvah .and. hsabs <= max(hsrel,PNUMS(2)) ) .and. &
               (tmcurv <= curvat .and. tmabs <= max(tmrel,PNUMS(3)) ) ) npacct = npacct + 1

          !if(myrank.EQ.0) WRITE(PRINTF,*) '2',npacct
          !
          if (tstfl) then
             if (lhead) write(PRTEST,11)
             write (PRTEST,12) ivert, hsabs, hsabs/hscurr(ivert), hscurv/hscurr(ivert), tmabs, tmabs/tmcurr(ivert), tmcurv/tmcurr(ivert)
             lhead = .false.
          endif
          !
       else
          !
          hscurr(ivert) = 1.e-20
          tmcurr(ivert) = 1.e-20
          !
       endif
       !
    enddo
    !
    ! global sum to npacc and nwetp
    !
    npacc = npacc + npacct
    nwetp = nwetp + nwetpt
    !
    ! perform global reductions in parallel run
    !
    SBUF = 0 ; RBUF1 = 0

    RBUF1(1) = nwetp
    SBUF(1)  = nwetp
    !
    RBUF1(2) = npacc
    SBUF(2)  = npacc

    !WRITE(PRTEST,*)'rank npacc',myrank,npacc

    !
    ! Sum over partition
    !
    ! # if defined (MULTIPROCESSOR)
    IF(.NOT.SERIAL) call MPI_REDUCE(SBUF,RBUF1,2,MPI_INTEGER,MPI_SUM,0,comm,ierr)
   
    nwetp = RBUF1(1)
    npacc = RBUF1(2)
!    print*,'nwetp,cpu=',nwetp,myrank
!    print*,'npacc,cpu=',npacc,myrank


    ! compute percentage of active vertices where required accuracy has been reached
    !
    !if ( npacc > 0. ) accur = ceiling(DBLE(npacc)*10000./DBLE(nwetp))/100.
    if( npacc.GT.0) accur = REAL(npacc)*100._rkind/REAL(nwetp)
    !
 11 format(13x,'dHabs          ','dHrel          ','Curvature H    ','dTabs          ','dTrel          ','Curvature T    ')
 12 format(1x,ss,'k=',i7,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2,'  ',1pe13.6e2)
    
    !
end subroutine SwanConvStopc
