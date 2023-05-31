!
!****************************************************************
! ! jerome not used in unswan
  SUBROUTINE SWAPAR(IG,NICMAX,DEP2,KWAVE,CGO)
!
!****************************************************************
!
!     computes the wave parameters K and CGO in the nearby
!     points, depending of the sweep direction.
!     The nearby points are indicated with the index IC (see
!     FUNCTION ICODE(_,_)
!
!  Method
!
!     The wave number K(IS,iC) is computed with the dispersion relation:
!
!     S = GRAV K(IS,IC)tanh(K(IS,IC)DEP(IX,IY))
!
!     where S = is logarithmic distributed via LOGSIG
!
!     The group velocity CGO in the case without current is equal to
!
!                    1       K(IS,IC)DEP(IX,IY)          S
!     CGO(IS,IC) = ( - + --------------------------) -----------
!                    2   2 sinh 2K(IS,IC)DEP(IX,IY)  |k(IS,IC)|
!
!
  USE M_GENARR
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4
!# if defined (EXPLICIT)
!  USE MOD_ACTION_EX
!# else         
!  USE MOD_ACTION_IM
!# endif    
!  USE ALL_VARS 
  USE schism_glbl, ONLY : npa,rkind,ic3,indel
   
  IMPLICIT NONE                                                     

  INTEGER :: IC,IS,ID,IG,NICMAX,INDX
  REAL(rkind)    :: NN(1:MSC), ND(1:MSC)                                    
  REAL(rkind)    :: DEP2(npa),KWAVE(MSC,NICMAX),CGO(MSC,NICMAX)
  REAL(rkind)    :: DEPLOC

  DO IC = 1, NICMAX
    IF(IC == 1)THEN
      INDX  = IG
      DEPLOC = DEP2(INDX)
      IF(DEPLOC <= DEPMIN)THEN
!     *** depth is negative ***
        DO IS = 1, MSC
          KWAVE(IS,IC) = -1.                                           
          CGO(IS,IC)   = 0.                                            
        END DO	 
      ELSE
!     *** call KSCIP1 to compute KWAVE and CGO ***
        CALL KSCIP1(MSC,SPCSIG,DEPLOC,KWAVE(1,IC),CGO(1,IC),NN,ND)                                 
      ENDIF
    ELSE                                                    
!JL Fil Later
#if 0
      INDX  = NBVE(IG,IC-1)
      DEPLOC = DEP2(NV(INDX,1))+DEP2(NV(INDX,2))+DEP2(NV(INDX,3))
      DEPLOC = DEPLOC/3.0
      IF(DEPLOC <= DEPMIN)THEN
!     *** depth is negative ***
        DO IS = 1, MSC
          KWAVE(IS,IC) = -1.                                           
          CGO(IS,IC)   = 0.                                            
        END DO	 
      ELSE
!     *** call KSCIP1 to compute KWAVE and CGO ***
        CALL KSCIP1(MSC,SPCSIG,DEPLOC,KWAVE(1,IC),CGO(1,IC),NN,ND)
      ENDIF
#endif

#if 1
      INDX = indel(IC-1,IG)
      DEPLOC = DEP2(ic3(1,INDX))+DEP2(ic3(2,INDX))+DEP2(ic3(3,INDX))
      DEPLOC = DEPLOC/3.0
      IF(DEPLOC <= DEPMIN)THEN
!     *** depth is negative ***
        DO IS = 1, MSC
          KWAVE(IS,IC) = -1.
          CGO(IS,IC)   = 0.
        END DO
      ELSE
!     *** call KSCIP1 to compute KWAVE and CGO ***
        CALL KSCIP1(MSC,SPCSIG,DEPLOC,KWAVE(1,IC),CGO(1,IC),NN,ND)
      ENDIF
#endif      
 
    END IF  
  ENDDO                                                              

  RETURN
  END SUBROUTINE SWAPAR
 
 !
!
!****************************************************************
! ! ! jerome not used in unswan
  SUBROUTINE SWAPAR1(I,IS,ID,DEP2,KWAVEL,CGOL)
!
!****************************************************************
!
!     computes the wave parameters K and CGO in the nearby
!     points, depending of the sweep direction.
!     The nearby points are indicated with the index IC (see
!     FUNCTION ICODE(_,_)
!
!  Method
!
!     The wave number K(IS,iC) is computed with the dispersion relation:
!
!     S = GRAV K(IS,IC)tanh(K(IS,IC)DEP(IX,IY))
!
!     where S = is logarithmic distributed via LOGSIG
!
!     The group velocity CGO in the case without current is equal to
!
!                    1       K(IS,IC)DEP(IX,IY)          S
!     CGO(IS,IC) = ( - + --------------------------) -----------
!                    2   2 sinh 2K(IS,IC)DEP(IX,IY)  |k(IS,IC)|
!
  USE M_GENARR
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4 
!# if defined (EXPLICIT)
!  USE MOD_ACTION_EX
!# else        
!  USE MOD_ACTION_IM
!# endif    
!  USE ALL_VARS, ONLY : NV,NTVE,MT
  USE schism_glbl, ONLY : npa,elnode,rkind
   
  IMPLICIT NONE                                                     
!
  INTEGER :: IC,IS,ID,I,INDX
  REAL(rkind)    :: NN(1:MSC), ND(1:MSC)                                    
  REAL(rkind)    :: DEP2(npa),KWAVEL,CGOL
  REAL(rkind)    :: DEPLOC,SPCSIGL

  DEPLOC = (DEP2(elnode(1,I))+DEP2(elnode(2,I))+DEP2(elnode(3,I)))/3.0
  IF(DEPLOC <= DEPMIN)THEN
!   *** depth is negative ***
    KWAVEL = -1.                                           
    CGOL   = 0.                                            
  ELSE
!   *** call KSCIP1 to compute KWAVE and CGO ***
    SPCSIGL = SPCSIG(IS)
    CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVEL,CGOL,NN,ND)                                 
  END IF

  RETURN
  END SUBROUTINE SWAPAR1
!    
!****************************************************************
! Jerome: not use in unswan, to remove 
   SUBROUTINE SPROXY (I1     ,IS     ,ID     ,CAXL   ,CAYL   ,  &
                      CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L   )
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3
   USE SWCOMM4
   USE OCPCOMM4
   USE M_DIFFR
!#  if defined (EXPLICIT)
!   USE MOD_ACTION_EX
!#  else      
!   USE MOD_ACTION_IM
!#  endif   

   IMPLICIT NONE

   INTEGER  IC,IS,ID,I1,NICMAX
   REAL     CAXL,CAYL,CG0L,ECOSL,ESINL,UX2L,UY2L

   CAXL = CG0L * ECOSL
   CAYL = CG0L * ESINL
!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN
     CAXL = CAXL + UX2L
     CAYL = CAYL + UY2L
   END IF

   RETURN
   END SUBROUTINE SPROXY
!
!

!****************************************************************
! ! Jerome: not use in unswan, to remove
   SUBROUTINE SPROXY2 (CAXL   ,CAYL   ,  &
                      CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L   )
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3                                                         
   USE SWCOMM4                                                         
   USE OCPCOMM4                                                        
   USE M_DIFFR 
!#  if defined (EXPLICIT)
!   USE MOD_ACTION_EX
!#  else      
!   USE MOD_ACTION_IM
!#  endif   
   
   IMPLICIT NONE                                              

   INTEGER  IC,IS,ID,I1,NICMAX
   REAL     CAXL(MDC,MSC),CAYL(MDC,MSC),CG0L(MDC,MSC),ECOSL(MDC),ESINL(MDC),UX2L,UY2L

   DO ID=1,MDC
     CAXL(ID,:) = CG0L(ID,:) * ECOSL(ID)
     CAYL(ID,:) = CG0L(ID,:) * ESINL(ID)
   END DO

!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN 
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN 
     CAXL = CAXL + UX2L
     CAYL = CAYL + UY2L
   END IF

   RETURN
   END SUBROUTINE SPROXY2
!
!
!****************************************************************
! ! Jerome: not use in unswan, to remove
   SUBROUTINE SPROXY3 (CAXL   ,CAYLA ,CAYLB,  &          !yzhang_w3
              CG0L   ,ECOSL  ,ESINL  ,UX2L   ,UY2L, DLTYETMPP, DLTXETMPP ,DLTXEA , DLTXEB)
!
!****************************************************************
!
!        computes the propagation velocities of energy in X-, Y-
!        -space, i.e., CAX, CAY, in the presence or absence of
!        currents, for the action balance equation.
!
!        The propagation velocities are computed for the fully 360
!        degrees sector.
!
!     METHOD
!
!        The next equation are calculated:
!
!              @X     _
!        CAX = -- = n C cos (id) + Ux  = CGO cos(id) + Ux
!              @T
!
!              @Y     _
!        CAY = -- = n C sin(id)  + Uy  = CGO sin(id) + Uy
!              @T
!                                                         _
!
!       ******************************************************************
!       *  attention! in the action balance equation the term            *
!       *  dx                                                            *
!       *  -- = CGO + U = CX  with x, CGO, U and CX vectors              *
!       *  dt                                                            *
!       *  is in the literature the term dx/dt often indicated           *
!       *  with CX and CY in the action balance equation.                *
!       *  In this program we use:    CAX = CGO + U                      *
!       ******************************************************************
!
!   ------------------------------------------------------------
!   If depth is negative ( DEP(IX,IY) <= 0), then,
!     For every point in S and D-direction do,
!       Give propagation velocities default values :
!       CAX(ID,IS,IC)     = 0.   {propagation velocity of energy in X-dir.}
!       CAY(ID,IS,IC)     = 0.   {propagation velocity of energy in Y-dir.}
!     ---------------------------------------------------------
!   Else if current is on (ICUR > 0) then,
!     For every point in S and D-direction do,  {using the output of SWAPAR}
!       S = logaritmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S cos(D)
!       CAX = ( - + ------------------------) --------- + UX2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)      S sin(D)
!       CAY = ( - + ------------------------) -------- + UY2(IX,IY)
!               2   sinh 2K(IS,IC)DEP2(IX,IY) |K(IS,IC)|
!
!       ------------------------------------------------------
!   Else if current is not on (ICUR = 0)
!     For every point in S and D-direction do
!       S = logarithmic distributed via LOGSIG
!       Compute propagation velocity in X-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S cos(D)
!       CAX = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!       ------------------------------------------------------
!       Compute propagation velocity in Y-direction:
!
!               1    K(IS,IC)DEP2(IX,IY)        S sin(D)
!       CAY = ( - + ------------------------) ----------
!               2   sinh 2K(IS,IC)DEP2(IX,IY)  |K(IS,IC)|
!
!     ----------------------------------------------------------
!   End IF
!   ------------------------------------------------------------
!****************************************************************
   USE SWCOMM3
   USE SWCOMM4
   USE OCPCOMM4
   USE M_DIFFR
!   USE MOD_ACTION_EX

   IMPLICIT NONE
   INTEGER  IC,IS,ID,I1,NICMAX
   REAL   CAXL(MDC,MSC),CAYLA(MDC,MSC),CAYLB(MDC,MSC),CG0L(MDC,MSC),ECOSL(MDC),ESINL(MDC),UX2L,UY2L
   REAL     DLTXETMPP,DLTXEA,DLTXEB,DLTYETMPP

   DO ID=1,MDC
     CAXL(ID,:) = CG0L(ID,:) * ECOSL(ID) * DLTYETMPP
     CAYLA(ID,:) = CG0L(ID,:) * ESINL(ID) * DLTXEA
     CAYLB(ID,:) = CG0L(ID,:) * ESINL(ID) * DLTXEB
   END DO

!
!    --- adapt the velocities in case of diffraction
!
   IF(IDIFFR == 1 .AND. PDIFFR(3) /= 0.)THEN
!     CAXL = CAXL*DIFPARAM(I1)      
!     CAYL = CAYL*DIFPARAM(I1)      
   END IF
!
!    --- ambient currents added
!
   IF(ICUR == 1)THEN
     CAXL = CAXL + UX2L*DLTYETMPP
     CAYLA = CAYLA + UY2L*DLTXETMPP
     CAYLB = CAYLB + UY2L*DLTXETMPP
   END IF

   RETURN
   END SUBROUTINE SPROXY3
!
!
!*******************************************************************
!
      SUBROUTINE ADDDIS (DISSXY     ,LEAKXY     ,  &
!     &                   AC2        ,ANYBIN     ,  &
     &                              ANYBIN     ,  &

     &                   DISC0      ,DISC1      ,  &
     &                   GENC0      ,GENC1      ,  &                      ! 40.85
     &                   REDC0      ,REDC1      ,  &                      ! 40.85
     &                   TRAC0      ,TRAC1      ,  &                      ! 40.85
     &                   IMATLA     ,IMATUA     ,  &                      ! 40.85
     &                   IMAT5L     ,IMAT6U     ,  &                      ! 40.85
     &                   DSXBOT     ,              &                      ! 40.67 40.61
     &                   DSXSRF     ,              &                      ! 40.67 40.61
     &                   DSXWCP     ,              &                      ! 40.67 40.61
     &                   DSXVEG     ,DSXTUR     ,  &                      ! 40.35 40.67 40.61
     &                   DSXMUD     ,              &                      ! 40.67 40.61
     &                   GSXWND     ,GENRXY     ,  &                      ! 40.85
     &                   RSXQUA     ,RSXTRI     ,  &                      ! 40.85
     &                   REDSXY     ,              &                      ! 40.85
     &                   TSXGEO     ,TSXSPT     ,  &                      ! 40.85
     &                   TSXSPS     ,TRANXY     ,  &                      ! 40.85
     &                   LEAKC1     ,RADSXY     ,SPCSIG     )             ! 40.85 30.72

!
!*******************************************************************
!
      USE SWCOMM3                                                        ! 40.41

      USE VARS_WAVE, ONLY : AC2

      USE M_GENARR, ONLY: SPCDIR
      USE VARS_WAVE, ONLY : WK,DS_INCR,SBR,SBF
      USE schism_glbl,ONLY:rkind,iplg
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
!     30.72: IJsbrand Haagsma
!     40.61: Marcel Zijlema
!     40.67: Nico Booij
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     20.53, Aug. 95: New subroutine
!     30.74, Nov. 97: Prepared for version with INCLUDE statements
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     40.61, Sep. 06: introduction of all separate dissipation coefficients
!     40.67, Jun. 07: more accurate computation of dissipation terms
!     40.85, Aug. 08: add also propagation, generation and redistribution terms
!                     and radiation stress
!
!  2. Purpose
!
!     Adds propagation, generation, dissipation, redistribution, leak and
!     radiation stress terms
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL(rkind)    SPCSIG(MSC)                                                ! 30.72
!
!     IX          Counter of gridpoints in x-direction
!     IY          Counter of gridpoints in y-direction
!     MXC         Maximum counter of gridppoints in x-direction
!     MYC         Maximum counter of gridppoints in y-direction
!     MSC         Maximum counter of relative frequency
!     MDC         Maximum counter of directional distribution
!
!     one and more dimensional arrays:
!     ---------------------------------
!     AC2       4D    Action density as function of D,S,X,Y and T
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
!     SWANCOMPUNSTRUC
!
! 11. Remarks
!
!     DISSXY and LEAKXY are dissipation and leak integrated over the
!     spectrum for each point in the computational grid
!     The same holds for DSXBOT, DSXSRF and DSXWCP for bottom friction-,
!     surf- and whitecapping dissipation, respectively
!     DISSC0 and DISSC1 give the dissipation distributed over the
!     spectral space in one point of the computational grid
!     The same holds for DISBOT, DISSRF and DISWCP for bottom friction-,
!     surf- and whitecapping dissipation, respectively
!
! 12. Structure
!
!     -------------------------------------------------------------
!     -------------------------------------------------------------
!
! 13. Source text
!
!
      INTEGER :: II               ! counter                               40.67
      REAL(rkind)    :: ADISSIP(1:MDISP)                                         !40.67
      REAL(rkind)    :: AGENERT(1:MGENR)                                         !40.85
      REAL(rkind)    :: AREDIST(1:MREDS)                                         !40.85
      REAL(rkind)    :: ATRANSP(1:MTRNP)                                         !40.85
      REAL(rkind)    :: DISSXY(MCGRD)    ,LEAKXY(MCGRD)      ,&
     &         DSXBOT(MCGRD)      ,                  &                     !40.67 40.61
     &         DSXSRF(MCGRD)      ,                  &                     !40.67 40.61
     &         DSXWCP(MCGRD)      ,                  &                     !40.67 40.61
     &         DSXMUD(MCGRD)      ,                  &                     !40.67 40.61
     &         DSXVEG(MCGRD)      ,                  &                     !40.67 40.61
     &         DSXTUR(MCGRD)      ,                  &                     !40.35
     &         GSXWND(MCGRD)      ,                  &                     !40.85
     &         RSXQUA(MCGRD)      ,                  &                     !40.85
     &         RSXTRI(MCGRD)      ,                  &                     !40.85
     &         TSXGEO(MCGRD)      ,                  &                     !40.85
     &         TSXSPT(MCGRD)      ,                  &                     !40.85
     &         TSXSPS(MCGRD)      ,                  &                     !40.85
     &         GENRXY(MCGRD)      ,REDSXY(MCGRD)    ,TRANXY(MCGRD), &      !40.85
     &         RADSXY(MCGRD)      ,                  &                     !40.85
     &         LEAKC1(MDC,MSC)  !,AC2(MDC,MSC,MCGRD)                       30.21
      REAL(rkind) :: DISC0(1:MDC,1:MSC,1:MDISP)       ! dissipation coeff.       40.67
      REAL(rkind) :: DISC1(1:MDC,1:MSC,1:MDISP)       ! dissipation coeff.       40.67
      REAL(rkind) :: GENC0(1:MDC,1:MSC,1:MGENR)       ! generation coeff.        40.85
      REAL(rkind) :: GENC1(1:MDC,1:MSC,1:MGENR)       ! generation coeff.        40.85
      REAL(rkind) :: REDC0(1:MDC,1:MSC,1:MREDS)       ! redistribution coeff.    40.85
      REAL(rkind) :: REDC1(1:MDC,1:MSC,1:MREDS)       ! redistribution coeff.    40.85
      REAL(rkind) :: TRAC0(1:MDC,1:MSC,1:MTRNP)       ! transport coeff.         40.85
      REAL(rkind) :: TRAC1(1:MDC,1:MSC,1:MTRNP)       ! transport coeff.         40.85
      REAL(rkind) :: IMATLA(MDC,MSC)           ,&                                 !40.85
     &               IMATUA(MDC,MSC)           ,&                                 !40.85
     &               IMAT5L(MDC,MSC)           ,&                                 !40.85
     &               IMAT6U(MDC,MSC)                                              !40.85
!
      REAL(rkind) :: ARADSTR, DSDD, SDSDD, S1, ACT1, S2
      REAL(rkind) :: ACT2, S3, ACT3, ACT4, ACT5, ACONTR
      INTEGER :: ISC, IDC, IDM, IDP
!
! JL add
      REAL(rkind) :: SINT,COST

      LOGICAL :: ANYBIN(MDC,MSC)
      INTEGER, SAVE :: IENT=0
      CALL STRACE (IENT, 'ADDDIS')


!      IF (IPLG(KCGRD(1))==9216)THEN
!        write(12,*)'Bfr SBR m2/s',KCGRD(1)
!        write(12,*) SBR(1,KCGRD(1)),SBR(2,KCGRD(1))
!      ENDIF

!
      ADISSIP(1:MDISP) = 0._rkind                                        ! 40.67
      AGENERT(1:MGENR) = 0._rkind                                        ! 40.85
      AREDIST(1:MREDS) = 0._rkind                                        ! 40.85
      ATRANSP(1:MTRNP) = 0._rkind                                        ! 40.85
      ARADSTR          = 0._rkind                                        ! 40.85
      DO 100 ISC = 1, MSC
        DSDD  = DDIR * FRINTF * SPCSIG(ISC)
        SDSDD = DSDD * SPCSIG(ISC)
        DO 90 IDC = 1, MDC
          IDM = MOD ( IDC - 2 + MDC , MDC ) + 1
          IDP = MOD ( IDC     + MDC , MDC ) + 1
!
          S1   = SPCSIG(ISC)
          ACT1 = AC2(IDC,ISC,KCGRD(1))
          IF (ISC.EQ.1) THEN
             S2   = 0._rkind
             ACT2 = 0._rkind
          ELSE
             S2   = SPCSIG(ISC-1)
             ACT2 = AC2(IDC,ISC-1,KCGRD(1))
          ENDIF
          IF (ISC.EQ.MSC) THEN
             S3   = 0._rkind
             ACT3 = 0._rkind
          ELSE
             S3   = SPCSIG(ISC+1)
             ACT3 = AC2(IDC,ISC+1,KCGRD(1))
          ENDIF
          IF (.NOT.FULCIR .AND. IDC.EQ.1) THEN
             ACT4 = 0._rkind
          ELSE
             ACT4 = AC2(IDM,ISC,KCGRD(1))
          ENDIF
          IF (.NOT.FULCIR .AND. IDC.EQ.MDC) THEN
             ACT5 = 0._rkind
          ELSE
             ACT5 = AC2(IDP,ISC,KCGRD(1))
          ENDIF
!
          IF (ANYBIN(IDC,ISC)) THEN
            LEAKXY(KCGRD(1)) = LEAKXY(KCGRD(1)) + SDSDD*            &
     &                      LEAKC1(IDC,ISC) * AC2(IDC,ISC,KCGRD(1))
!
!           --- compute for each dissipation term                         40.61
!
            DO II = 1, MDISP                                              !40.67
              ACONTR= SDSDD*(DISC0(IDC,ISC,II) + DISC1(IDC,ISC,II)*ACT1)  !40.85
              ADISSIP(II) = ADISSIP(II) + ACONTR                          !40.85 40.67
              ARADSTR     = ARADSTR     - ACONTR                          !40.85
            ENDDO                                                         !40.67
!
!           --- compute for each generation term                          40.85
!
            DO II = 1, MGENR                                              !40.85
              ACONTR= SDSDD*(GENC0(IDC,ISC,II) + GENC1(IDC,ISC,II)*ACT1)  !40.85
              AGENERT(II) = AGENERT(II) + ACONTR                          !40.85
              ARADSTR     = ARADSTR     + ACONTR                          !40.85
            ENDDO                                                         !40.85
!
!           --- compute for each redistribution term                      40.85
!
            DO II = 1, MREDS                                              !40.85
              ACONTR= SDSDD*(REDC0(IDC,ISC,II) + REDC1(IDC,ISC,II)*ACT1)  !40.85
              AREDIST(II) = AREDIST(II) + ABS(ACONTR)                     !40.85
              ARADSTR     = ARADSTR     + ACONTR                          !40.85
            ENDDO                                                         !40.85
!
!           --- compute for each propagation term
!
            ACONTR = SDSDD* (TRAC0(IDC,ISC,1) + TRAC1(IDC,ISC,1)*ACT1)    !40.85
            ATRANSP(1) = ATRANSP(1) + ABS(ACONTR)                         !40.85
            ARADSTR    = ARADSTR    - ACONTR                              !40.85
!
            ACONTR = SDSDD* (TRAC0(IDC,ISC,2)      +&                      !40.85
     &                       TRAC1(IDC,ISC,2)*ACT1 +&                      !40.85
     &                       IMATLA(IDC,ISC) *ACT4 +&                      !40.85
     &                       IMATUA(IDC,ISC) *ACT5 )                      !40.85
            ATRANSP(2) = ATRANSP(2) + ABS(ACONTR)                         !40.85
            ARADSTR    = ARADSTR    - ACONTR                              !40.85
!
            ACONTR = DSDD * (TRAC0(IDC,ISC,3)           +&                 !40.85
     &                       TRAC1(IDC,ISC,3)* S1 *ACT1 +&                 !40.85
     &                       IMAT5L(IDC,ISC) * S2 *ACT2 +&                 !40.85
     &                       IMAT6U(IDC,ISC) * S3 *ACT3 )                 !40.85
            ATRANSP(3) = ATRANSP(3) + ABS(ACONTR)                         !40.85
            ARADSTR    = ARADSTR    - ACONTR                              !40.85
!
            ARADSTR    = ABS(ARADSTR)                                     !40.85

#if 1
!#ifdef SCHISM
          ! JL: schism coupling

          ! Compute  sbr(2,npa): momentum flux vector due to nearshore
          ! depth-induced breaking; see Bennis 2011), see WWMIII/wwm_breaking.F90

            COST = SPCDIR(IDC,2)  !  SPCDIR(ID,2) !COSTH(ID)
            SINT = SPCDIR(IDC,3)  !  SPCDIR(ID,3) !SINTH(ID)
!  Brk
            ACONTR = (DISC0(IDC,ISC,2)+DISC1(IDC,ISC,2)*ACT1)
!+ WCAP ?
!           ACONTR = ACONTR +  (DISC0(IDC,ISC,1)+DISC1(IDC,ISC,1)*ACT1)
            ACONTR = ACONTR*DDIR*DS_INCR(ISC) !*SPCSIG(IS)
            SBR(1,KCGRD(1)) = SBR(1,KCGRD(1)) -GRAV_W*COST*WK(ISC,KCGRD(1))*ACONTR !/SPCSIG(ISC)
            SBR(2,KCGRD(1)) = SBR(2,KCGRD(1)) -GRAV_W*SINT*WK(ISC,KCGRD(1))*ACONTR  !/SPCSIG(ISC)

          ! Compute sbf(2,npa): momentum lost by waves due to the bottom friction (
          ! See method in WWMIII/wwm_friction.F90
          ! JL: Please, check again and again
            ACONTR = (DISC0(IDC,ISC,3)+DISC1(IDC,ISC,3)*ACT1)
            ACONTR = ACONTR*DDIR*DS_INCR(ISC) !*SPCSIG(ISC)
            SBF(1,KCGRD(1)) = SBF(1,KCGRD(1)) -COST*WK(ISC,KCGRD(1))*ACONTR !(probablement *G)
            SBF(2,KCGRD(1)) = SBF(2,KCGRD(1)) -SINT*WK(ISC,KCGRD(1))*ACONTR !(probablement *G)
!#endif

#endif

          ENDIF
  90    CONTINUE
 100  CONTINUE


!      IF (IPLG(KCGRD(1))==9216)THEN
!        write(12,*)'ADISSIP m2/s',KCGRD(1)
!        write(12,*) SBR(1,KCGRD(1)),SBR(2,KCGRD(1))
!      ENDIF


!
      DSXWCP(KCGRD(1)) = DSXWCP(KCGRD(1)) + ADISSIP(1)     ! whitecapping      40.67
      DSXSRF(KCGRD(1)) = DSXSRF(KCGRD(1)) + ADISSIP(2)     ! surf break        40.67
      DSXBOT(KCGRD(1)) = DSXBOT(KCGRD(1)) + ADISSIP(3)     ! bottom fric       40.67
      DSXVEG(KCGRD(1)) = DSXVEG(KCGRD(1)) + ADISSIP(5)     ! vegetation        40.67
      DSXTUR(KCGRD(1)) = DSXTUR(KCGRD(1)) + ADISSIP(6)     ! turbulence        40.35
      DSXMUD(KCGRD(1)) = DSXMUD(KCGRD(1)) + ADISSIP(7)     ! mud dissip        40.67
      DISSXY(KCGRD(1)) = DISSXY(KCGRD(1)) + SUM(ADISSIP)   ! total dissip      40.67
!
      GSXWND(KCGRD(1)) = GSXWND(KCGRD(1)) + AGENERT(1)     ! wind input        40.85
      GENRXY(KCGRD(1)) = GENRXY(KCGRD(1)) + SUM(AGENERT)   ! total generation  40.85
!
      RSXQUA(KCGRD(1)) = RSXQUA(KCGRD(1)) + AREDIST(1)     ! quadruplets       40.85
      RSXTRI(KCGRD(1)) = RSXTRI(KCGRD(1)) + AREDIST(2)     ! triads            40.85
      REDSXY(KCGRD(1)) = REDSXY(KCGRD(1)) + SUM(AREDIST)   ! total redistribution 40.85
!
      TSXGEO(KCGRD(1)) = TSXGEO(KCGRD(1)) + ATRANSP(1)     ! xy-propagation    40.85
      TSXSPT(KCGRD(1)) = TSXSPT(KCGRD(1)) + ATRANSP(2)     ! theta-propagation 40.85
      TSXSPS(KCGRD(1)) = TSXSPS(KCGRD(1)) + ATRANSP(3)     ! sigma-propagation 40.85
      TRANXY(KCGRD(1)) = TRANXY(KCGRD(1)) + SUM(ATRANSP)   ! total propagation 40.85
!
!       energy transfer between waves and currents due to radiation stress, see
!       page 439 of
!       the ICCE paper of Holthuijsen, L.H., Zijlema, M. and Van der Ham, P.J.
!       (2009)
!       Wave physics in a tidal inlet, in: J.M. Smith (Ed.), Proc. 31st Int.
!       Conf. on Coast. Engng.
!       ASCE, World Scientific Publishing, Singapore, 437-448
!
      RADSXY(KCGRD(1)) = RADSXY(KCGRD(1)) + ARADSTR                      ! 40.85
!
      IMATLA = 0.
      IMATUA = 0.
      IMAT5L = 0.
      IMAT6U = 0.
!
      RETURN
      END

!****************************************************************
!
!      SUBROUTINE SPREDT (SWPDIR     ,AC2        ,CAX       , &
      SUBROUTINE SPREDT (SWPDIR                 ,CAX       , &
     &                   CAY        ,IDCMIN     ,IDCMAX    , &
     &                   ISSTOP     ,ANYBIN     ,            &
     &                   XCGRID     ,YCGRID     ,            &            !41.53
     &                   RDX        ,RDY        ,OBREDF    )              !40.00
!
!****************************************************************
!
      USE SWCOMM2                                                         !41.53
      USE SWCOMM3                                                         !40.41
      USE SWCOMM4                                                         !40.41
      USE OCPCOMM4                                                        !40.41

      USE VARS_WAVE, ONLY : AC2
      USE schism_glbl,ONLY:rkind
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
!     0. Authors
!
!     40.00, 40.13: Nico Booij
!     40.41: Marcel Zijlema
!     41.53: Marcel Zijlema
!
!     1. UPDATE
!
!        40.00, Aug 98: introduction of obstacle reduction factor to
!                       obtain correct initialisation
!                       argument list changed, swcomm3 added
!        40.13, Aug 01: modification of action densities is skipped
!                       in case of Mode Noupdate
!        40.41, Oct 04: common blocks replaced by modules, include files removed
!        41.53, Oct 14: correction curvilinear grid
!
!     2. PURPOSE
!
!        to estimate the action density depending of the sweep
!        direction during the first iteration of a stationary
!        computation. The reason for this is that AC2 is zero
!        at first iteration and no initialisation is given in
!        case of stationarity (NSTATC=0). Action density should
!        be nonzero because of the computation of the source
!        terms. The estimate is based on solving the equation
!
!            dN       dN
!        CAX -- + CAY -- = 0
!            dx       dy
!
!        in an explicit manner. In the estimate, the transmission
!        through obstacles or reflection at obstacles is taken into
!        account
!
!     3. METHOD
!
!
!          [RDX1*CAX + RDY1*CAY]*N(i-1,j) + [RDX2*CAX + RDY2*CAY]*N(i,j-1)
! N(i,j) = ---------------------------------------------------------------
!                      (RDX1+RDX2) * CAX  +  (RDY1+RDY2) * CAY
!
!     4. PARAMETERLIST
!
!       INTEGERS:
!       ---------
!       IC           Dummy variable: ICode gridpoint:
!                    IC = 1  Top or Bottom gridpoint
!                    IC = 2  Left or Right gridpoint
!                    IC = 3  Central gridpoint
!                    Whether which value IC has, depends of the sweep
!                    If necessary ic can be enlarged by increasing
!                    the array size of ICMAX
!       IX           Counter of gridpoints in x-direction
!       IY           Counter of gridpoints in y-direction
!       IS           Counter of relative frequency band
!       ID           Counter of directional distribution
!       ICMAX        Maximum array size for the points of the molecule
!       MXC          Maximum counter of gridppoints in x-direction
!       MYC          Maximum counter of gridppoints in y-direction
!       MSC          Maximum counter of relative frequency
!       MDC          Maximum counter of directional distribution
!       KSX          Dummy variable to get the right sign in the
!                    numerical difference scheme in X-direction
!                    depending on the sweep direction, KSX = -1 or +1
!       KSY          Dummy variable to get the right sign in the
!                    numerical difference scheme in Y-direction
!                    depending on the sweep direction, KSY = -1 or +1
!       SWPDIR       Sweep direction (..) (identical at the description
!                    of the direction the wind is blowing)
!
!       REALS:
!       ------
!
!       DX           Length of spatial cell in X-direction
!       DY           Length of spatial cell in Y-direction
!       ALEN         Part of side length of an angle side
!       BLEN         Part of side length of an angle side
!       LDIAG        Length of the diagonal of grid cel
!       ALPHA        angle of propagation velocity
!       BETA         angle between DX end DY
!       GAMMA        PI - alpha - beta
!       PI           3,14.......
!       FAC_A        Factor representing the influence of the action-
!                    density depening of the propagation velocity
!       FAC_B        Factor representing the influence of the action-
!                    density depening of the propagation velocity
!
!       REAL arrays:
!       -------------
!
!       AC2    4D    Action density as function of D,S,X,Y at time T
!       CAX    3D    Wave transport velocity in x-direction, function of
!                    (ID,IS,IC)
!       CAY    3D    Wave transport velocity in y-direction, function of
!                    (ID,IS,IC)
!       IDCMIN 1D    frequency dependent counters in case of a current
!       IDCMAX 1D    frequency dependent counters in case of a current
!       ANYBIN 2D    Determines if a bin fall within a sweep
!       XCGRID       coordinates of computational grid in x-direction
!       YCGRID       coordinates of computational grid in y-direction
!
!     5. SUBROUTINES CALLING
!
!        SWANCOMPUNSTRUC
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every sweep direction do,
!     For every point in S and D direction in sweep direction do,
!       predict values for action density at new point from values
!       of neighbour gridpoints taking into account spectral propagation
!       direction (with currents !!) and the boundary conditions.
!       --------------------------------------------------------
!       If wave action AC2 is negative, then
!         Give wave action initial value 1.E-10
!     ---------------------------------------------------------
!   End of SPREDT
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
!
      INTEGER  IS    ,ID    , &
     &         SWPDIR,IDDUM ,ISSTOP                                       !40.00
!
      REAL(rkind)     FAC_A ,FAC_B, WEIG1, WEIG2, TCF1, TCF2, CDEN, CNUM
!
!      REAL(rkind)  :: AC2(MDC,MSC,MCGRD)                                         !30.21
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL(rkind)  :: CAX(MDC,MSC,MICMAX)                                        !40.22
      REAL(rkind)  :: CAY(MDC,MSC,MICMAX)                                        !40.22
      REAL(rkind)  :: RDX(MICMAX),  RDY(MICMAX), &                               !40.08 30.21
     &         OBREDF(MDC,MSC,2)                                          !40.00
      REAL(rkind)  :: XCGRID(MXC,MYC), YCGRID(MXC,MYC)                           !41.53
!
      INTEGER  IDCMIN(MSC)              ,&
     &         IDCMAX(MSC)
!
      LOGICAL  ANYBIN(MDC,MSC)
!
      INTEGER  IDIR, IDCUM, IDMIN, IDMAX, NCURID, NID(4)                  !41.53
      REAL(rkind)     IDX, IDY, CLAT                                             !41.53
!
      INTEGER IENT
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SPREDT')
!
      IF ( OPTG.EQ.3 .AND. CCURV ) THEN
!
!        --- curvilinear grid                                             41.53
!            consider the equation on a rectangular grid                  41.53
!            so that all bins will be filled                              41.53
!
        IDX=1./(XCGRID(IXCGRD(1),IYCGRD(1))-XCGRID(IXCGRD(2),IYCGRD(2)))
        IDY=1./(YCGRID(IXCGRD(1),IYCGRD(1))-YCGRID(IXCGRD(3),IYCGRD(3)))
!
        IF ( KSPHER.GT.0 ) THEN
           CLAT = COS(DEGRAD*(YCGRID(IXCGRD(1),IYCGRD(1))+YOFFS))
           IDX  = IDX / (CLAT * LENDEG)
           IDY  = IDY / LENDEG
        ENDIF
!
        IDCUM = 0
        DO IDIR = 1, 4
           NID(IDIR) = (MDC*IDIR)/4 - IDCUM
           IDCUM     = (MDC*IDIR)/4
        ENDDO
        !
        ! determine loop bounds in spectral space for current sweep
        !
        IDMIN = MDC+1
        IDMAX = 0
        !
        IDIR   = 1
        NCURID = 0
        !
        DO ID = 1, MDC
           !
           IF ( IDIR.EQ.SWPDIR ) THEN
              IDMIN = MIN(ID,IDMIN)
              IDMAX = MAX(ID,IDMAX)
           ENDIF
           NCURID = NCURID + 1
           !
           IF ( NCURID.GE.NID(IDIR) ) THEN
              IDIR   = IDIR + 1
              NCURID = 0
           ENDIF
           !
        ENDDO
!
        DO IS = 1, MSC
          DO ID = IDMIN, IDMAX
!
             IF ( NUMOBS.GT.0 ) THEN
                TCF1 = OBREDF(ID,IS,1)
                TCF2 = OBREDF(ID,IS,2)
             ELSE
                TCF1 = 1.
                TCF2 = 1.
             ENDIF
!
             CDEN = IDX * CAX(ID,IS,1) + IDY * CAY(ID,IS,1)
!
             CNUM = IDX * CAX(ID,IS,2) * TCF1 * AC2(ID,IS,KCGRD(2)) + &
     &              IDY * CAY(ID,IS,3) * TCF2 * AC2(ID,IS,KCGRD(3))
!
             IF (ACUPDA) AC2(ID,IS,KCGRD(1)) = CNUM / CDEN
!
          ENDDO
        ENDDO
!
      ELSE
        DO IS = 1, ISSTOP
          DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            IF ( ANYBIN(ID,IS) ) THEN
!
!             *** Computation of weighting coefs WEIG1 AND WEIG2 ***
!
              CDEN = RDX(1) * CAX(ID,IS,1) + RDY(1) * CAY(ID,IS,1)
              CNUM =  (RDX(1) + RDX(2)) * CAX(ID,IS,1) &
     &              + (RDY(1) + RDY(2)) * CAY(ID,IS,1)
              WEIG1 = CDEN/CNUM
              WEIG2 = 1. - WEIG1
!
              IF (NUMOBS .GT. 0) THEN
                TCF1 = OBREDF(ID,IS,1)                                    !40.00
                TCF2 = OBREDF(ID,IS,2)                                    !40.00
              ELSE
                TCF1 = 1.
                TCF2 = 1.
              ENDIF
              FAC_A = TCF1 * WEIG1 * AC2(ID,IS,KCGRD(2))                  !40.00
              FAC_B = TCF2 * WEIG2 * AC2(ID,IS,KCGRD(3))                  !40.00
!
              IF (ACUPDA) &                                               !40.13
     &           AC2(ID,IS,KCGRD(1)) = MAX ( 0. , (FAC_A + FAC_B))        !30.21
!
            END IF
          END DO
        END DO
      END IF
!
! JL Skip now
# if 0

      IF ( ITEST .GE. 140 .AND. TESTFL ) THEN
        WRITE(PRINTF,6019) KCGRD(1), SWPDIR
 6019   FORMAT(' PREDT : POINT INDX  SWPDIR         :',2I5)
        DO IS = 1, ISSTOP
          DO IDDUM = IDCMIN(IS)-1, IDCMAX(IS)+1
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
            WRITE (PRINTF,6020) IS, ID, AC2(ID,IS,KCGRD(1)),
     &                          AC2(ID,IS,KCGRD(2)),
     &                          AC2(ID,IS,KCGRD(3)),
     &                          ANYBIN(ID,IS)
 6020       FORMAT ('       : IS ID AC2 AC2(2) AC2(3) ANYBIN   :',
     &              2I5,3(E12.4),L4)
          END DO
        END DO
      END IF
# endif

!     End of the subroutine SPREDT
      RETURN
      END

!
!******************************************************************
! jerome not used in unswan ??
  SUBROUTINE SWPSEL(IDCMIN    ,IDCMAX    ,CAX       ,             &
                    CAY       ,ISCMIN    ,ISCMAX    ,             &
		    IDTOT     ,ISTOT     ,IDDLOW    ,             &
		    IDDTOP    ,ISSTOP    ,DEP2      ,             &
		    UX2       ,UY2       ,SPCDIR    )
!
!******************************************************************
!
!     compute the frequency dependent counters in directional space
!     in a situation with a current and without a current.
!     The counters are only computed for the gridpoint
!     considered. This means IC = 1 
!
!******************************************************************
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OCPCOMM4                                                        
! USE ALL_VARS, ONLY : MT                                                       
  USE schism_glbl, ONLY : npa,rkind

  IMPLICIT NONE
!
  REAL(rkind)    :: SPCDIR(MDC,6)                                               
  INTEGER :: IS ,ID ,IDSUM ,IDCLOW ,IDCHGH ,IDTOT ,ISTOT ,    & 
             IDDLOW ,IDDTOP ,ISSLOW ,ISSTOP ,IENT ,IDDUM ,ISCLOW ,    &
	     ISCHGH ,IX ,IY ,IC         
  REAL(rkind)    :: CAXMID ,CAYMID ,GROUP ,UABS ,THDIR    
  INTEGER :: IDCMIN(MSC) ,IDCMAX(MSC) ,ISCMIN(MDC) ,ISCMAX(MDC) ,SECTOR(MSC)
  REAL(rkind)    :: CAX(MDC,MSC,MICMAX) ,CAY(MDC,MSC,MICMAX) ,DEP2(npa) ,     &
             UX2(npa) ,UY2(npa) ,RDX(10) ,RDY(10)    
  LOGICAL :: LOWEST ,LOWBIN ,HGHBIN
!
! *** initialize arrays in frequency direction ***
!
  DO ID = 1, MDC
    ISCMIN(ID) = 1
    ISCMAX(ID) = 1
  END DO  
!
! *** initialize array's in theta direction ***
!
  DO IS = 1, MSC
    IDCMIN(IS) = 1
    IDCMAX(IS) = MDC
    SECTOR(IS) = 1
  END DO  
!
!     *** calculate minimum and maximum counters in frequency ***
!     *** space if a current is present: ISCMIN and ISCMAX    ***
!
  IDDLOW =  9999
  IDDTOP = -9999
  DO IS = 1 , MSC
    IF(SECTOR(IS) > 0)THEN
      IDDLOW = MIN ( IDDLOW , IDCMIN(IS) )
      IDDTOP = MAX ( IDDTOP , IDCMAX(IS) )
    END IF
  END DO
!
!     *** Determine counters ***
!
  DO IDDUM = IDDLOW, IDDTOP
    ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
    LOWEST = .TRUE.
    DO IS = 1, MSC
      IF(LOWEST)THEN
        ISCLOW = IS
        LOWEST = .FALSE.
      END IF
      ISCHGH = IS
    END DO
!
!       *** set the minimum and maximum counters in arrays ***
!
    IF(.NOT.LOWEST)THEN
      ISCMIN(ID) = ISCLOW
      ISCMAX(ID) = ISCHGH
    ELSE
    END IF
!
  END DO   !IDDUM
!
!     *** calculate the maximum number of counters in both ***
!     *** directional space and frequency space            ***
!
  IF(IDDLOW /= 9999)THEN
    IDTOT = ( IDDTOP - IDDLOW ) + 1
    IF(ICUR == 1)THEN
      IF(IDTOT < 3)THEN
        IDDTOP = IDDTOP + 1
        IF(IDTOT == 1) IDDLOW = IDDLOW - 1
        IDTOT = 3
      END IF
    END IF
  ELSE
    IDTOT = 0
  END IF
!
! *** set variables ***
!
  IDTOT  =     1
  ISTOT  =     1
  ISSLOW =  9999                                                     
  ISSTOP = -9999                                                     

  DO IS = 1, MSC
    IDCLOW  = 0
    IDCHGH  = 0
    IDSUM   = 0
    DO ID = 1, MDC
      IDSUM = IDSUM + 1
      ISSLOW = MIN(IS,ISSLOW)                                   
      ISSTOP = MAX(IS,ISSTOP)                                   
    END DO
  END DO  
!
  IF(ISSLOW /= 9999)THEN
    ISSLOW = 1
!   minimal value of ISSTOP is 4 (or MSC if MSC<4)                    
    IF(ICUR > 0) ISSTOP = MAX(MIN(4,MSC),ISSTOP)                    
    ISTOT = ( ISSTOP - ISSLOW ) + 1
  ELSE
    ISTOT = 0
  END IF
!
!     *** check if IDTOT is less then MDC ***
!
  IF(IDTOT > MDC)THEN
    IDDLOW = 1
    IDDTOP = MDC
    IDTOT  = MDC
  END IF
!
!     *** check if the lowest frequency is not blocked !    ***
!     *** this can occur in real cases if the depth is very ***
!     *** small and the current velocity is large           ***
!     *** the propagation velocity Cg = sqrt (gd) < U       ***
!
  IF(ICUR == 1 .AND. FULCIR .AND.                      &
     ISSLOW /= 1 .AND. ISSLOW /= 9999)THEN                          
!        CALL MSGERR (2,'The lowest freqency is blocked')                  
7002 FORMAT (I4, 1X, I4, 1X, I2)                                       
    IC = 1
    GROUP = SQRT ( GRAV_W * DEP2(KCGRD(IC)) )
  END IF

  RETURN
  END SUBROUTINE SWPSEL

!****************************************************************
!
      SUBROUTINE STRSD (DD      ,IDCMIN  ,&
     &                  IDCMAX  ,CAD     ,IMATLA  ,IMATDA  ,IMATUA  ,&
!     &                  IMATRA  ,AC2     ,ISSTOP  ,&
     &                  IMATRA  ,         ISSTOP  ,&
     &                  ANYBIN  ,LEAKC1  ,TRAC0   ,TRAC1            )    ! 40.85 40.41 30.21

!****************************************************************
!
      USE SWCOMM3                                                        ! 40.41
      USE SWCOMM4                                                        ! 40.41
      USE OCPCOMM4                                                       ! 40.41

      USE schism_glbl, ONLY :rkind,iplg
      USE schism_msgp, ONLY : myrank
      USE VARS_WAVE, ONLY: AC2

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
!     1. UPDATE
!
!        40.41, Oct. 04: common blocks replaced by modules, include files
!        removed
!        40.85, Aug. 08: store theta-propagation for output purposes
!
!     2. PURPOSE
!
!        comp. of @[CAD AC2]/@D initial & boundary
!
!     3. METHOD
!
!        Compute the derivative in D-direction only n the central
!        gridpoint considered:
!                                Central grid point     : IC = 1
!

!       Depending on the parameter PNUMS(6) either a central difference
!       scheme (PNUMS(6) = 0) or an upstream scheme (PNUMS(6) = 1) is
!       used. Points 1, 2 and 3 are three consecutive points on the
!       T-axis. 2 is the central point for which @(C*A)/@THETA and
!       @(C*W*A)/@THETA is computed.
!
!                 1       2       3
!              ---O-------O-------O--- > THETA
!
!
!        PNUMS() = 0.  central difference scheme
!        PNUMS() = 1.  upwind scheme
!
!        @[CAD AC2]
!        ----------  =
!           @D
!
!        CAD(ID+1,IS,1) AC2(ID+1,IS,IX,IY) - CAD(ID-1,IS,1) AC2(ID-1,IS,IX,IY)
!        --------------------------------------------------------------------
!                                         2*DD
!
!     4. PARAMETERLIST
!
!        IC          Dummy variable: ICode gridpoint:
!                      IC = 1  Top or Bottom gridpoint
!                      IC = 2  Left or Right gridpoint
!                      IC = 3  Central gridpoint
!                    Whether which value IC has, depends of the sweep
!                    If necessary IC can be enlarged by increasing
!                    the array size of ICMAX
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        ICMAX       Maximum counter for the points of the molecule
!        MXC         Maximum counter of gridpoints in x-direction
!        MYC         Maximum counter of gridpoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!        FULCIR      logical: if true, computation on a full circle
!
!        REALS:
!        ---------
!
!        DD          Width of spectral direction band
!        PNH         Equal to (1/2)*DD
!        PI          (3,14)
!
!        one and more dimensional arrays:
!        ---------------------------------
!
!        CAD     3D  Wave transport velocity in S-dirction, function of
!                    (ID,IS,IC)
!        IMATDA  2D  Coefficients of diagonal of matrix
!        IMATLA  2D  Coefficients of lower diagonal of matrix
!        IMATUA  2D  Coefficients of upper diagonal of matrix
!        IDCMIN  1D  frequency dependent counter
!        IDCMIN  1D  frequency dependent counter
!        ANYBIN  2D  see subr SWPSEL
!        LEAKC1  2D  leak coefficient
!
!     5. SUBROUTINES CALLING
!
!        swantranspac             ACTION
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!   -----------------------------------------------------------
!   For every S and D-direction in direction of sweep do
!     Compute the derivative in D-direction:
!     ---------------------------------------------------------
!     Store the results of the transport terms in the
!     arrays IMATDA, IMATLA, IMATUA
!   -------------------------------------------------------------
!   End of STRSD
!   ------------------------------------------------------------
!
!     10. SOURCE

!****************************************************************
!
      LOGICAL  BIN1, BIN2, BIN3
!
      INTEGER  IENT  ,IS    ,ID    ,IIDM  ,IIDP  ,&
     &         ISSTOP,IDDUM
!
      REAL(rkind)     DD    ,PNH   ,PN1   ,PN2
!
!      REAL(rkind)     AC2(MDC,MSC,MCGRD)         ,&                               !30.21
      REAL(rkind)     &
     &         CAD(MDC,MSC,ICMAX)         ,&
     &         IMATLA(MDC,MSC)            ,&
     &         IMATDA(MDC,MSC)            ,&
     &         IMATUA(MDC,MSC)            ,&
     &         IMATRA(MDC,MSC)            ,&
     &         LEAKC1(MDC,MSC)
      REAL(rkind)  :: TRAC0(MDC,MSC,MTRNP)                                       !40.85
      REAL(rkind)  :: TRAC1(MDC,MSC,MTRNP)                                       !40.85

      REAL(rkind)  :: A1,A3,C1,C2,C3,PCD1,PCD2,PCD3
      REAL(rkind)  :: RHS12,RHS23,DIAG12,DIAG23
!
      INTEGER  IDCMIN(MSC)                ,&
     &         IDCMAX(MSC)
!
      LOGICAL  ANYBIN(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'STRSD')
!
      PNH = 1._rkind / (2._rkind * DD)
      PN1 =  (1._rkind - PNUMS(6) ) * PNH
      PN2 =  (1._rkind + PNUMS(6) ) * PNH
!
      DO 200 IS = 1, ISSTOP
        DO 100 IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD (IDDUM-1+MDC, MDC) + 1
          C2 = CAD(ID,IS,1)
          BIN2 = ANYBIN(ID,IS)
          IF (BIN2) THEN
            IF (FULCIR .OR. ID.GT.1) THEN
              IIDM = MOD (IDDUM-2+MDC, MDC) + 1
              C1   = CAD(IIDM,IS,1)
              BIN1 = ANYBIN(IIDM,IS)

              !IF(AC2(IIDM,IS,KCGRD(1))/=AC2(IIDM,IS,KCGRD(1)))THEN
              !    print*,'Nan Found,AC2,IIDM,IS,KCGRD(1)',IIDM,IS,KCGRD(1)
              !ENDIF

              IF (.NOT.BIN1) A1 = AC2(IIDM,IS,KCGRD(1))
            ELSE
              IIDM = 0
              C1   = C2
              BIN1 = .FALSE.
              A1   = 0._rkind
            ENDIF
            IF (FULCIR .OR. ID.LT.MDC) THEN
              IIDP = MOD (IDDUM+MDC, MDC) + 1
              C3   = CAD(IIDP,IS,1)
              BIN3 = ANYBIN(IIDP,IS)

              !IF(AC2(IIDP,IS,KCGRD(1))/=AC2(IIDP,IS,KCGRD(1)))THEN
              !    print*,'Nan Found,AC2,IIDP,IS,KCGRD(1)',IIDP,IS,KCGRD(1)
              !ENDIF


              IF (.NOT.BIN3) A3 = AC2(IIDP,IS,KCGRD(1))                  ! 30.21
            ELSE
              IIDP = 0
              C3   = C2
              BIN3 = .FALSE.
              A3   = 0._rkind
            ENDIF
!
!           *** fill the lower diagonal and the diagonal ***
!
            IF ( C1 .GT. 1.E-8 .AND. C2 .GT. 1.E-8 ) THEN
              PCD1 = PN2 * C1
              PCD2 = PN1 * C2
            ELSE IF ( C1 .LT. -1.E-8 .AND. C2 .LT. -1.E-8 ) THEN
              PCD1 = PN1 * C1
              PCD2 = PN2 * C2
            ELSE
              PCD1 = PNH * C1
              PCD2 = PNH * C2
            END IF
!
            RHS12 = 0._rkind
            DIAG12 = 0._rkind ! Add JL
            IF (IIDM.EQ.0 .AND. C2.LT.0._rkind) THEN
!             fully upwind approximation at the boundary of the directional
!             sector
              DIAG12 = - PCD1 - PCD2
              LEAKC1(ID,IS) = -C2
            ELSE
              DIAG12 = - PCD2
              IF (BIN1) THEN
                IMATLA(ID,IS) = IMATLA(ID,IS) - PCD1
              ELSE
                RHS12 = PCD1 * A1
              ENDIF
            ENDIF
!
            IF ( C2 .GT. 1.E-8 .AND. C3 .GT. 1.E-8 ) THEN
              PCD2 = PN2 * C2
              PCD3 = PN1 * C3
            ELSE IF ( C2 .LT. -1.E-8 .AND. C3 .LT. -1.E-8 ) THEN
              PCD2 = PN1 * C2
              PCD3 = PN2 * C3
            ELSE
              PCD2 = PNH * C2
              PCD3 = PNH * C3
            END IF
!
            RHS23 = 0._rkind
            DIAG23 = 0._rkind ! Add JL
            IF (IIDP.EQ.0 .AND. C2.GT.0._rkind) THEN
!             full upwind approximation at the boundary
              DIAG23 = PCD2 + PCD3
              LEAKC1(ID,IS) = C2
            ELSE
              DIAG23 = PCD2
              IF (BIN3) THEN
                IMATUA(ID,IS) = IMATUA(ID,IS) + PCD3
              ELSE
                RHS23 = - PCD3 * A3
                !IF(RHS23/=RHS23)THEN
                !  print*,'Nan Found,PCD3,A3',PCD3,A3
                !ENDIF
              ENDIF
            ENDIF
            !IF(myrank.EQ.0) THEN
            !   print*,'RHS12,RHS23=',RHS12,RHS23
            !ENDIF
            IMATDA(ID,IS) = IMATDA(ID,IS) + DIAG12 + DIAG23
            IMATRA(ID,IS) = IMATRA(ID,IS) + RHS12 + RHS23
            TRAC0(ID,IS,2) = TRAC0(ID,IS,2) - RHS12 - RHS23              ! 40.85
            TRAC1(ID,IS,2) = TRAC1(ID,IS,2) + DIAG12 + DIAG23            ! 40.85
          ENDIF
!
 100    CONTINUE
 200  CONTINUE
!
!     *** test output
!
      IF ( ITEST .GE. 80 .AND. TESTFL ) THEN
        WRITE(PRINTF,5001) FULCIR
5001    FORMAT (' FULL CIRCLE                   ',L4)
        WRITE(PRINTF,6021) KCGRD(1), iplg(KCGRD(1)),ISSTOP, PNUMS(6)                    ! 30.21
6021    FORMAT (' STRSD :POINT INDEX_GL ISTOP CDD      :',3I7,E12.4)
        WRITE(PRINTF,5021) PN1, PN2, PNH ,DD
5021    FORMAT (' STRSD : PN1 PN2 PNH DD      :',4E12.4)
      END IF
!
!     End of subroutine STRSD
      RETURN
 END SUBROUTINE
!****************************************************************
!
     SUBROUTINE STRSSI(SPCSIG  ,&
     &                  CAS     ,IMAT5L  ,IMATDA  ,IMAT6U  ,ANYBIN  ,&
!     &                  IMATRA  ,AC2     ,ISCMIN  ,ISCMAX  ,IDDLOW  ,&
     &                  IMATRA  ,         ISCMIN  ,ISCMAX  ,IDDLOW  ,&
     &                  IDDTOP  ,TRAC0   ,TRAC1                     )     !40.85 0.41 30.21

!
!****************************************************************
!
      USE SWCOMM3                                                        ! 40.41
      USE SWCOMM4                                                        ! 40.41
      USE OCPCOMM4                                                       ! 40.41

      USE schism_glbl,ONLY:rkind,iplg
      USE VARS_WAVE, ONLY: AC2

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
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!     40.85: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.85, Aug. 08: store sigma-propagation for output purposes
!
!  2. Purpose
!
!     comp. of @[CAS AC2]/@S initial & boundary : IMPLICIT SCHEME
!
!  3. Method
!
!     Compute the derivative in S-direction only n the central
!     gridpoint considered:
!                             Central grid point     : IC = 1
!
!     Depending on the parameter PNUMS(7) either a central difference
!     scheme (PNUMS(7) = 0) or an upstream scheme (PNUMS(7) = 1) is
!     used. Points 1, 2 and 3 are three consecutive points on the
!     T-axis. 2 is the central point for which @(C*A)/@SIGMA and
!     @(C*W*A)/@SIGMA is computed.
!
!               1       2       3
!            ---O-------O-------O--- > SIGMA
!
!
!     PNUMS() = 0.  central difference scheme
!     PNUMS() = 1.  upwind scheme
!
!     @[CAS AC2]
!     ----------  =
!        @S
!
!     CAS(ID,IS+1,1) AC2(ID,IS+1,IX,IY) - CAS(ID,IS-1,1) AC2(ID,IS-1,IX,IY)
!     --------------------------------------------------------------------
!                                      2 DS
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL(rkind)    SPCSIG(MSC)                                                ! 30.72
!
!        IC          Dummy variable: ICode gridpoint:
!                      IC = 1  Top or Bottom gridpoint
!                      IC = 2  Left or Right gridpoint
!                      IC = 3  Central gridpoint
!                    Whether which value IC has, depends of the sweep
!                    If necessary IC can be enlarged by increasing
!                    the array size of ICMAX
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        ICMAX       Maximum counter for the points of the molecule
!        MXC         Maximum counter of gridpoints in x-direction
!        MYC         Maximum counter of gridpoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!                    one sweep
!
!
!        REALS:
!        ---------
!
!        DD          Width of spectral direction band
!        PNH         Equal to (1/2)*DD
!        PI          (3,14)
!
!        one and more dimensional arrays:
!        ---------------------------------
!
!        CAS     3D  Wave transport velocity in S-dirction, function of
!                    (ID,IS,IC)
!        IMATDA  2D  Coefficients of diagonal of matrix
!        IMAT5L  2D  Coefficients of lower diagonal of matrix
!        IMAT6U  2D  Coefficients of upper diagonal of matrix
!        ISCMIN  1D  Minimum counter in frequency space per direction
!        ISCMIN  1D  Maximum counter in frequency space per direction
!
!     5. SUBROUTINES CALLING
!
!        swantranspac     ACTION
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. Common blocks used
!
!
!        ---
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!   -----------------------------------------------------------
!   For every S and D-direction in direction of sweep do
!     Compute the derivative in S-direction:
!     ---------------------------------------------------------
!     Store the results of the transport terms in the
!     arrays IMATDA, IMAT5L, IMAT6U
!   -------------------------------------------------------------
!   End of STRSSI
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!****************************************************************
!
      INTEGER  IENT, IS      ,ID      ,IDDLOW  ,IDDTOP  ,IDDUM
!
      REAL(rkind)     DS      ,PNH     ,PN1     ,PN2  ,C1    ,C2      ,&
     &         C3      ,A1      ,A3      ,PCD1  ,PCD2 ,PCD3  ,RHS12   ,&
     &         RHS23   ,DIAG12  ,DIAG23
      REAL(rkind)     TC12    ,TC23    ,S1      ,S3
!
      LOGICAL  BIN1    ,BIN3
!
!      REAL(rkind)   AC2(MDC,MSC,MCGRD)         ,&
      REAL(rkind)   &
     &         CAS(MDC,MSC,ICMAX)         ,&
     &         IMAT5L(MDC,MSC)            ,&
     &         IMATDA(MDC,MSC)            ,&
     &         IMAT6U(MDC,MSC)            ,&
     &         IMATRA(MDC,MSC)
      REAL(rkind)  :: TRAC0(MDC,MSC,MTRNP)                                      ! 40.85
      REAL(rkind)  :: TRAC1(MDC,MSC,MTRNP)                                      ! 40.85
!
      INTEGER  ISCMIN(MDC)                ,&
     &         ISCMAX(MDC)
!
      LOGICAL  ANYBIN(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'STRSSI')
!
      DO 500 IDDUM = IDDLOW, IDDTOP
        ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
        IF (ISCMIN(ID).EQ.0) GOTO 500
        DO 400 IS = ISCMIN(ID), ISCMAX(ID)
          A1 = 0._rkind
          A3 = 0._rkind
          C2 = CAS(ID,IS,1)
          IF ( IS .EQ. 1 ) THEN
            C1   = 0._rkind
            A1   = 0._rkind
            BIN1 = .FALSE.
            C3   = CAS(ID,IS+1,1)
            BIN3 = ANYBIN(ID,IS+1)
            IF (.NOT.BIN3) A3 = AC2(ID,IS+1,KCGRD(1))                    ! 30.21
            DS   = SPCSIG(IS+1) - SPCSIG(IS)                             ! 30.72
            S1 = 0._rkind                                                      ! 40.85
            S3 = SPCSIG(IS+1)                                            ! 40.85
          ELSE IF ( IS .EQ. MSC ) THEN
            C1   = CAS(ID,IS-1,1)
            BIN1 = ANYBIN(ID,IS-1)
            IF (.NOT.BIN1) A1 = AC2(ID,IS-1,KCGRD(1))                    ! 30.21
            C3   = C2
            A3   = 0._rkind
            BIN3 = .FALSE.
            DS   = SPCSIG(IS) - SPCSIG(IS-1)                             ! 30.72
            S1 = SPCSIG(IS-1)                                            ! 40.85
            S3 = 0._rkind                                                      ! 40.85
          ELSE
            C1   = CAS(ID,IS-1,1)
            C3   = CAS(ID,IS+1,1)
            BIN1 = ANYBIN(ID,IS-1)
            BIN3 = ANYBIN(ID,IS+1)
            IF (.NOT.BIN1) A1 = AC2(ID,IS-1,KCGRD(1))                    ! 30.21
            IF (.NOT.BIN3) A3 = AC2(ID,IS+1,KCGRD(1))                    ! 30.21
            DS   = 0.5_rkind * ( SPCSIG(IS+1) - SPCSIG(IS-1) )                 ! 30.72

            S1 = SPCSIG(IS-1)                                            ! 40.85
            S3 = SPCSIG(IS+1)                                            ! 40.85
          END IF
!
          PNH = 1. / (2. * DS)
          PN1 =  (1. - PNUMS(7) ) * PNH
          PN2 =  (1. + PNUMS(7) ) * PNH
!
!         *** fill the lower diagonal and the diagonal ***
!
          IF ( C1 .GT. 1.E-8 .AND. C2 .GT. 1.E-8 ) THEN
            PCD1 = PN2 * C1
            PCD2 = PN1 * C2
          ELSE IF ( C1 .LT. -1.E-8 .AND. C2 .LT. -1.E-8 ) THEN
            PCD1 = PN1 * C1
            PCD2 = PN2 * C2
          ELSE
            PCD1 = PNH * C1
            PCD2 = PNH * C2
          END IF
!
          RHS12 = 0._rkind
          TC12  = 0._rkind
          IF ( IS .EQ. 1 .AND. C2.LT.0.) THEN
!           fully upwind approximation at the boundary of the frequency space
            DIAG12 = - PCD1 - PCD2
          ELSE
            DIAG12 = - PCD2
            IF (BIN1) THEN
              IMAT5L(ID,IS) = IMAT5L(ID,IS) - PCD1
            ELSE
              RHS12 = PCD1 * A1
              TC12  = RHS12* S1
            ENDIF
          ENDIF
!
          IF ( C2 .GT. 1.E-8 .AND. C3 .GT. 1.E-8 ) THEN
            PCD2 = PN2 * C2
            PCD3 = PN1 * C3
          ELSE IF ( C2 .LT. -1.E-8 .AND. C3 .LT. -1.E-8 ) THEN
            PCD2 = PN1 * C2
            PCD3 = PN2 * C3
          ELSE
            PCD2 = PNH * C2
            PCD3 = PNH * C3
          END IF
!
          RHS23 = 0._rkind
          TC23  = 0._rkind
          IF (IS .EQ. MSC .AND. C2.GT.0.) THEN
!           full upwind approximation at the boundary
            DIAG23 = PCD2 + PCD3
          ELSE
            DIAG23 = PCD2
            IF (BIN3) THEN
              IMAT6U(ID,IS) = IMAT6U(ID,IS) + PCD3
            ELSE
              RHS23 = - PCD3 * A3
              TC23  = RHS23 * S3
            ENDIF
          ENDIF
          IMATDA(ID,IS) = IMATDA(ID,IS) + DIAG12 + DIAG23

          IMATRA(ID,IS) = IMATRA(ID,IS) + RHS12 + RHS23
          TRAC0(ID,IS,3) = TRAC0(ID,IS,3) - TC12 - TC23                  ! 40.85
          TRAC1(ID,IS,3) = TRAC1(ID,IS,3) + DIAG12 + DIAG23              ! 40.85
400     CONTINUE
500   CONTINUE
!
!     *** test output ***
!
      IF ( TESTFL .AND. ITEST .GE. 35 ) THEN
        WRITE(PRINTF,111) KCGRD(1),iplg(KCGRD(1)), IDDLOW, IDDTOP
 111    FORMAT(' STRSSI: POINT INDEX_GL IDDLOW IDDTOP  :',4I7)
        WRITE(PRINTF,131) PNUMS(7)
 131    FORMAT(' STRSSI: CSS                  :',2E12.4)
        WRITE(PRINTF,*)
        WRITE(PRINTF,*) ' matrix coefficients in STRSSI'
        WRITE(PRINTF,*)
        WRITE(PRINTF,*)&
     & '   IS   ID    IMAT5L       IMATDA       IMAT6U    IMATRA    CAS'
        DO IDDUM = IDDLOW, IDDTOP
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          IF (ISCMIN(ID).GT.0) THEN
            DO IS = ISCMIN(ID), ISCMAX(ID)
              WRITE(PRINTF,2101) IS, ID, IMAT5L(ID,IS),IMATDA(ID,IS),&
     &                         IMAT6U(ID,IS),IMATRA(ID,IS),CAS(ID,IS,1)
2101          FORMAT(1X,2I4,4X,4E12.4,E10.2)
            ENDDO
          ENDIF
        ENDDO
      END IF
!
!     End of subroutine STRSSI
      RETURN
 END SUBROUTINE

!****************************************************************
!
      SUBROUTINE STRSSB (IDDLOW  ,IDDTOP  ,&
     &                   IDCMIN  ,IDCMAX  ,ISSTOP  ,CAX     ,CAY     ,&
!    &                   CAS     ,AC2     ,SPCSIG  ,IMATRA  ,&
    &                   CAS     ,          SPCSIG  ,IMATRA  ,&
     &                   ANYBLK  ,RDX     ,RDY     ,TRAC0            )   ! 41.07 40.41 30.21
!
!****************************************************************
!
      USE SWCOMM3                                                        ! 40.41
      USE SWCOMM4                                                        ! 40.41
      USE OCPCOMM4                                                       ! 40.41

      USE schism_glbl,ONLY:rkind,iplg
      USE VARS_WAVE, ONLY: AC2
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
!     30.72: IJsbrand Haagsma
!     40.41: Marcel Zijlema
!     41.07: Marcel Zijlema
!
!  1. Updates
!
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     41.07, Jul. 09: also central scheme blended with upwind scheme
!
!  2. Purpose
!
!     comp. of @[CAS AC2]/@S initial & boundary with an explicit
!     scheme. The energy near the blocking point is removed
!     from the spectrum based on a CFL criterion

!
!     The frequencies beyond ISSTOP are blocked in a 1-D situation
!     For a 2-D case the situation is somewhat more complicated (
!     see below)
!
!
!        ^  |                       1-D case
!     E()|  |          *            ========
!           |        *   *
!           |              *
!           |       *        *      / blocking frequency
!           |                .... /
!           |      *         ....| *
!           |        SWEEP 1 ....| o o *
!           |     *          ....| o o o o o*
!          0---------------------|-----------|---------
!           0                  ISSTOP       MSC   --> s
!
!                           -|---|-
!                              ^
!                              |---- CFL > 0.5sqrt(2) -> ANYBLK = true
!
!
!               ANYBIN = TRUE     ANYBIN = FALSE
!           |--------------------|-----------|
!
!
!  3. Method
!
!     Compute the derivative in s-direction:
!     The nearby points are indicated with the index IC (see
!     FUNCTION ICODE(_,_):
!     Central grid point     : IC = 1
!     Point in X-direction   : IC = 2
!     Point in Y-direction   : IC = 3
!
!     @[CAS AC2]
!     --------- =
!        @S
!
!     CAS*AC2(ID,IS) - CAS*AC2(ID,IS-1)     F(IS+0.5) - F(IS-0.5)
!     ---------------------------------- = -----------------------
!                   DS                                DS
!
!                  /  CAS(IS+0.5) * ( (1-0.5mu)*AC2(IS+1) + 0.5mu*AC2(IS) )
!                  IF CAS(IS+0.5) < 0
!     F(IS+0.5) =  |
!                  \  CAS(IS+0.5) * ( (1-0.5mu)*AC2(IS) + 0.5mu*AC2(IS+1) )
!                  IF CAS(IS+0.5) > 0
!
!                  /  CAS(IS-0.5) * ( (1-0.5mu)*AC2(IS-1) + 0.5mu*AC2(IS) )
!                  IF CAS(IS-0.5) > 0
!     F(IS-0.5) =  |
!                  \  CAS(IS-0.5) * ( (1-0.5mu)*AC2(IS) + 0.5mu*AC2(IS-1) )
!                  IF CAS(IS-0.5) < 0
!
!     with

!           0 <= mu <= 1 a blending factor
!
!           mu = 0 corresponds to 1st order upwind scheme
!           mu = 1 corresponds to 2nd order central scheme
!
!           default value, mu = 0.5
!
!    and
!
!           CAS(IS+0.5) = ( CAS(IS+1) + CAS(IS) ) / 2.
!
!           CAS(IS-0.5) = ( CAS(IS) + CAS(IS-1) ) / 2.
!
!
!     ------------------------------------------------------------
!     Courant-Friedlich-Levich criterion :
!
!                  | Cs |
!                  | -- |
!                  | ds |         <
!               ---------------   =  0.5 * sqrt(2.0)
!      CFL  =  | Cx |   | Cy |
!              | -- | + | -- |
!              | dx |   | dy |
!
!     For a bin in which the CFL criterion is larger two
!     ways are possible:
!
!            1)  Cs can be limited
!            2)  Action in bin can be set equal zero
!
!     --------------------------------------------------------------
!
!  4. Argument variables
!
!     SPCSIG: Relative frequencies in computational domain in sigma-space 30.72
!
      REAL(rkind)    SPCSIG(MSC)                                                ! 30.72
!
!        IC          Dummy variable: ICode gridpoint:
!                      IC = 1  Top or Bottom gridpoint
!                      IC = 2  Left or Right gridpoint
!                      IC = 3  Central gridpoint
!                    Whether which value IC has, depends of the sweep
!                    If necessary IC can be enlarged by increasing
!                    the array size of ICMAX
!        IX          Counter of gridpoints in x-direction
!        IY          Counter of gridpoints in y-direction
!        IS          Counter of relative frequency band
!        ID          Counter of directional distribution
!        ICMAX       Maximum counter for the points of the molecule
!        MXC         Maximum counter of gridpoints in x-direction
!        MYC         Maximum counter of gridpoints in y-direction
!        MSC         Maximum counter of relative frequency
!        MDC         Maximum counter of directional distribution
!        ISSTOP      Maximum frequency counter for wave components
!                    that are propagated within a sweep
!        IDDLOW      Minimum direction that is propagated within a
!                    sweep
!        IDDTOP      Idem maximum
!
!        REALS:
!        ---------
!
!        FSA_        Dummy variable
!
!        one and more dimensional arrays:
!        ---------------------------------
!
!        AC2     4D  Action density as function of D,S,X,Y at time T
!        CAS     3D  Wave transport velocity in S-dirction, function of
!                    (ID,IS,IC)
!        CAX, CAY    Propagation velocities in x-y space
!        IMATRA  2D  Coefficients of right hand side of matrix
!        ISCMIN  1D  Diractional dependent counter
!        ISCMIN  1D  Directional dependent counter
!        ANYBLK  2D  Determines if a bin is BLOCKED by a counter current
!                    based on a CFL criterion
!
!     5. SUBROUTINES CALLING
!
!        swantranspac   ACTION
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. Common blocks used
!
!
!     8. REMARKS
!
!        ---
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!   For every S and D-direction in direction of sweep do,
!     Determine if CFL criterion is satisfied
!     Compute the derivative in s-direction:
!     ---------------------------------------------------------
!     Compute transportation terms
!     Store the terms in the array IMATRA
!   -------------------------------------------------------------
!   End of STRSSB
!   -------------------------------------------------------------
!
!     10. SOURCE
!
!************************************************************************
      INTEGER  IENT, IS      ,ID      ,ISSTOP  ,&
     &         IDDLOW  ,IDDTOP  ,IDDUM
!
      REAL(rkind)     FSA     ,FLEFT   ,FRGHT   ,DS      ,CFLMAX  ,CFLCEN  ,&
     &         CAXCEN  ,CAYCEN  ,CASCEN  ,TX      ,TY      ,TS      ,&
     &         CASL    ,CASR    ,PN1     ,PN2
!
      REAL(rkind)     CAS(MDC,MSC,ICMAX)       ,&
     &         CAX(MDC,MSC,ICMAX)       ,&
     &         CAY(MDC,MSC,ICMAX)       ,&
!    &         AC2(MDC,MSC,MCGRD)       ,&
     &         IMATRA(MDC,MSC)          ,&
     &         RDX(MICMAX)              ,&                                ! 40.08
     &         RDY(MICMAX)                                                !40.08
      REAL(rkind)  :: TRAC0(MDC,MSC,MTRNP)                                       !41.07
!
      INTEGER  IDCMIN(MSC)              ,&
     &         IDCMAX(MSC)
!
      LOGICAL  ANYBLK(MDC,MSC)
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'STRSSB')
!
!     --- determine blending factor
!
      PN1 = 0.5_rkind*(1._rkind+PNUMS(7))                                            ! 41.07
      PN2 = 0.5_rkind*(1._rkind-PNUMS(7))                                            ! 41.07
!
!     *** initialization of ANYBLK and CFLMAX value ***
!
      DO IS = 1, MSC
        DO ID = 1, MDC
          ANYBLK(ID,IS) = .FALSE.
        ENDDO
      ENDDO
      CFLMAX = PNUMS(19)
!
      DO IS = 1, ISSTOP
        IF ( IS .EQ. 1 ) THEN
          DS = SPCSIG(IS+1) - SPCSIG(IS)                                 ! 30.72
        ELSE IF ( IS .EQ. MSC ) THEN
          DS = SPCSIG(IS) - SPCSIG(IS-1)                                 ! 30.72
        ELSE
          DS = 0.5_rkind * ( SPCSIG(IS+1) - SPCSIG(IS-1) )                     ! 30.72
        END IF
        DO IDDUM = IDCMIN(IS), IDCMAX(IS)
          ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
          CAXCEN = ABS ( CAX(ID,IS,1) )
          CAYCEN = ABS ( CAY(ID,IS,1) )
          CASCEN = ABS ( CAS(ID,IS,1) )
!
          TX     = RDX(1) * CAXCEN + RDX(2) * CAXCEN
          TY     = RDY(1) * CAYCEN + RDY(2) * CAYCEN
!
          TS     = CASCEN / DS
          CFLCEN = TS / MAX( 1.E-20 , ( TX + TY ) )
          FRGHT = 0._rkind
          FLEFT = 0._rkind
!
!         *** check if a bin can be propagated or if it is blocked ***
!
          IF ( CFLCEN .GT. CFLMAX ) THEN
!
!           *** de-activate bin in solver by ANYBLK ***
!
            ANYBLK(ID,IS) = .TRUE.
!
          ELSE
!
!           *** calculate transport in frequency space ***
!
            IF ( IS .EQ. 1 ) THEN
!             *** for first point an upwind scheme is used ***
              CASR  = 0.5_rkind * ( CAS(ID,IS,1) + CAS(ID,IS+1,1) )
              IF ( CASR .LT. 0._rkind ) THEN
                FRGHT = CASR * AC2(ID,IS+1,KCGRD(1))                     ! 30.21
              ELSE
                FRGHT = CASR * AC2(ID,IS  ,KCGRD(1))                     ! 30.21
              END IF
              FLEFT = 0._rkind
            ELSE IF ( IS .EQ. MSC ) THEN
!             *** for the last discrete point in frequency space ***
!             *** an upwind scheme is used                       ***
              CASL = CAS(ID,IS-1,1)
              CASR = CAS(ID,IS  ,1)
              IF ( CASL .LT. 0. ) THEN
                FLEFT = CASL * AC2(ID,IS  ,KCGRD(1))                     ! 30.21
              ELSE
                FLEFT = CASL * AC2(ID,IS-1,KCGRD(1))                     ! 30.21
              END IF
              IF ( CASR .LT. 0._rkind ) THEN
!               *** assumption has been made that the flux is ***
!               *** zero for the bin beyond MSC               ***
                FRGHT = 0._rkind
              ELSE
                FRGHT = CASR * AC2(ID,IS,KCGRD(1))                       ! 30.21
              END IF
            ELSE
!             *** point in frequency range ***
              CASL  = 0.5_rkind * ( CAS(ID,IS,1) + CAS(ID,IS-1,1) )
              CASR  = 0.5_rkind * ( CAS(ID,IS,1) + CAS(ID,IS+1,1) )
              IF ( CASL .LT. 0._rkind ) THEN
                FLEFT = CASL * ( PN1*AC2(ID,IS  ,KCGRD(1)) +&             ! 41.07
     &                           PN2*AC2(ID,IS-1,KCGRD(1)) )              !41.07 30.21
              ELSE
                FLEFT = CASL * ( PN1*AC2(ID,IS-1,KCGRD(1)) +&             ! 41.07
     &                           PN2*AC2(ID,IS  ,KCGRD(1)) )             ! 41.07 30.21
              END IF
              IF ( CASR .LT. 0. ) THEN
                FRGHT = CASR * ( PN1*AC2(ID,IS+1,KCGRD(1)) +&             ! 41.07
     &                           PN2*AC2(ID,IS  ,KCGRD(1)) )             ! 41.07 30.21
              ELSE
                FRGHT = CASR * ( PN1*AC2(ID,IS  ,KCGRD(1)) +&             ! 41.07
     &                           PN2*AC2(ID,IS+1,KCGRD(1)) )             ! 41.07 30.21
              END IF
            END IF
!
            FSA  = ( FRGHT - FLEFT ) / DS
!
!           *** all the terms are known, store in IMATRA ***

            IMATRA(ID,IS) = IMATRA(ID,IS) - FSA
            TRAC0(ID,IS,3) = TRAC0(ID,IS,3) + FSA                        !41.07
          ENDIF
!
!         *** test output ***

          IF ( ITEST .GE. 50 .AND. TESTFL ) THEN
            WRITE(PRINTF,670) IS,ID,FRGHT,FLEFT,CFLCEN,ANYBLK(ID,IS)
 670        FORMAT(' STRSSB: FR FL CFLC ANYBLK:',2I3,3E12.4,L3)
          END IF

!
        ENDDO
      ENDDO
!
!     *** test output ***
!
      IF ( ITEST .GE. 50 .AND. TESTFL ) THEN
        WRITE(PRINTF,200) MDC,MSC,MCGRD
 200    FORMAT(' BLOCKB : MDC MSC MCGRD    : ',3I5)
        WRITE(PRINTF,300) KCGRD(1), iplg(KCGRD(1)),ISSTOP, CFLMAX
 300    FORMAT(' BLOCKB : POINT INDEX_GL ISSTOP CFLMAX: ',3I7,F8.4)
        WRITE(PRINTF,400) IDDLOW, IDDTOP
 400    FORMAT (' Active bins within a sweep  -> ID: ',I3,' to ',I3)
        WRITE(PRINTF,*)
        WRITE(PRINTF,*)(' Propagation of bin if blocking can occur')
        WRITE(PRINTF,*)('   1) No blocking of bin -> ANYBLK = .F.')
        WRITE(PRINTF,*)('   2) Blocking of bin    -> ANYBLK = .T.')
        WRITE(PRINTF,*)
        DO IDDUM = IDDTOP+1, IDDLOW-1, -1
          ID = MOD ( IDDUM - 1 + MDC, MDC) + 1
            WRITE(PRINTF,500) ID, (ANYBLK(ID,IS),IS=1,MIN(ISSTOP,25))
 500        FORMAT(I4,25L3)
        ENDDO
        WRITE(PRINTF,600)(IS, IS=1+4, MIN(ISSTOP,25), 5 )
 600    FORMAT(6X,'1',9X,5(I3,12X))
        WRITE(PRINTF,*)
!
      ENDIF
!
!     End of subroutine STRSSB
      RETURN
  END SUBROUTINE
!
!****************************************************************

