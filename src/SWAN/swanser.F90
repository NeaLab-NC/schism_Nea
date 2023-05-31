#include "swan_functions.h"
!
!     SWAN - SERVICE ROUTINES
!
!  Contents of this file
!
!     READXY
!     REFIXY
!     INFRAM
!     DISTR
!     KSCIP1
!     AC2TST
!     CVCHEK                                                              30.60
!     CVMESH                                                              30.60
!     NEWTON                                                              30.60
!     EVALF                                                               30.60
!     SWOBST                                                              30.60
!     TCROSS                                                              40.04
!     SWTRCF
!     SSHAPE                                                              40.00
!     SINTRP                                                              40.00
!     HSOBND                                                              32.01
!     CHGBAS                                                              40.00
!     GAMMA                                                               40.00
!     WRSPEC                                                              40.00
!TIMG!     SWTSTA                                                              40.23
!TIMG!     SWTSTO                                                              40.23
!TIMG!     SWPRTI                                                              40.23
!     TXPBLA                                                              40.23
!     INTSTR                                                              40.23
!     NUMSTR                                                              40.23
!     SWCOPI                                                              40.23
!     SWCOPR                                                              40.23
!     SWI2B                                                               40.30
!     SWR2B                                                               40.30
!
!  functions:
!  ----------
!  DEGCNV  (converts from cartesian convention to nautical and            32.01
!           vice versa)                                                   32.01
!  ANGRAD  (converts radians to degrees)                                  32.01
!  ANGDEG  (converts degrees to radians)                                  32.01
!
!  subroutines:
!  ------------
!  HSOBND  (Hs is calculated after a SWAN computation at all sides.       32.01
!           The calculated wave height from SWAN is then compared with    32.01
!           the wave heigth as provided by the user                       32.01
!
!***********************************************************************
!                                                                      *
   SUBROUTINE READXY (NAMX, NAMY, XX, YY, KONT, XSTA, YSTA)
!                                                                      *
!***********************************************************************
!
   USE OCPCOMM1                                                        
   USE SWCOMM2                                                         
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!     40.22: John Cazes and Tim Campbell
!     40.13: Nico Booij
!     40.51: Marcel Zijlema
!
!  1. UPDATE
!
!       Nov. 1996               offset values are added to standard values
!                               because they will be subtracted later
!     40.13, Nov. 01: a valid value for YY is required if a valid value
!                     for XX has been given; ocpcomm1.inc reactivated
!     40.51, Feb. 05: correction to location points equal to offset values
!
!  2. PURPOSE
!
!       Read x and y, initialize offset values XOFFS and YOFFS
!
!  3. METHOD
!
!       ---
!
!  4. PARAMETERLIST
!
!       NAMX, NAMY   inp char    names of the two coordinates as given in
!                                the user manual
!       XX, YY       out real    values of x and y taking into account offset
!       KONT         inp char    what to be done if values are missing
!                                see doc. of INDBLE (Ocean Pack doc.)
!       XSTA, YSTA   inp real    standard values of x and y
!
!  5. SUBROUTINES CALLING
!
!
!
!  6. SUBROUTINES USED
!
!       INDBLE (Ocean Pack)
   LOGICAL EQREAL
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       ----------------------------------------------------------------
!       Read x and y in double prec.
!       If this is first couple of values
!       Then assign values to XOFFS and YOFFS
!            make LXOFFS True
!       ---------------------------------------------------------------
!       make XX and YY equal to x and y taking into account offset
!       ----------------------------------------------------------------
!
! 10. SOURCE TEXT
!
   DOUBLE PRECISION XTMP, YTMP
   CHARACTER  NAMX *(*), NAMY *(*), KONT *(*)
   SAVE  IENT
   DATA  IENT /0/
   CALL  STRACE (IENT,'READXY')
!
!JQI   CALL INDBLE (NAMX, XTMP, KONT, DBLE(XSTA)+DBLE(XOFFS))
   IF(CHGVAL)THEN                                                    
!       a valid value was given for XX                                    
!JQI     CALL INDBLE (NAMY, YTMP, 'REQ', DBLE(YSTA)+DBLE(YOFFS))           
   ELSE                                                                
!JQI     CALL INDBLE (NAMY, YTMP, KONT, DBLE(YSTA)+DBLE(YOFFS))
   ENDIF                                                               
   IF(.NOT.LXOFFS)THEN
     XOFFS = REAL(XTMP)
     YOFFS = REAL(YTMP)
     LXOFFS = .TRUE.
   ENDIF
!JQI   IF(.NOT.EQREAL(XOFFS,REAL(XTMP)))THEN                             
!JQI     XX = REAL(XTMP-DBLE(XOFFS))
!JQI   ELSE IF(OPTG == 3)THEN                                            
!JQI     XX = 1.E-5                                                       
!JQI   ELSE                                                                
!JQI     XX = 0.                                                          
!JQI   END IF                                                              
!JQI   IF(.NOT.EQREAL(YOFFS,REAL(YTMP)))THEN                             
!JQI     YY = REAL(YTMP-DBLE(YOFFS))
!JQI   ELSE IF(OPTG == 3)THEN                                            
!JQI     YY = 1.E-5                                                       
!JQI   ELSE                                                                
!JQI     YY = 0.                                                          
!JQI   END IF                                                              
!
   RETURN
   END SUBROUTINE READXY
 

!****************************************************************
!
   SUBROUTINE TXPBLA(TEXT,IF,IL)
!
!****************************************************************
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     determines the position of the first and the last non-blank
!     (or non-tabulator) character in the text-string
!
!  4. Argument variables
!
!     IF          position of the first non-blank character in TEXT
!     IL          position of the last non-blank character in TEXT
!     TEXT        text string
!
   INTEGER IF, IL
   CHARACTER*(*) TEXT
!
!  6. Local variables
!
!     FOUND :     TEXT is found or not
!     ITABVL:     integer value of tabulator character
!     LENTXT:     length of TEXT
!
   INTEGER LENTXT, ITABVL
   LOGICAL FOUND
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
!DOS      ITABVL = ICHAR('	')
!UNIX      ITABVL = ICHAR('	')
   LENTXT = LEN (TEXT)
   IF = 1
   FOUND = .FALSE.
100 IF(IF <= LENTXT .AND. .NOT. FOUND)THEN
     IF(.NOT. (TEXT(IF:IF) ==  ' ' .OR.                   &
       ICHAR(TEXT(IF:IF)) == ITABVL))THEN
       FOUND = .TRUE.
     ELSE
       IF = IF + 1
     ENDIF
     GOTO 100
   ENDIF
   IL = LENTXT + 1
   FOUND = .FALSE.
200 IF(IL > 1 .AND. .NOT. FOUND)THEN
    IL = IL - 1
     IF(.NOT. (TEXT(IL:IL) ==  ' ' .OR.                   &
       ICHAR(TEXT(IL:IL)) == ITABVL))THEN
       FOUND = .TRUE.
     ENDIF
     GOTO 200
   ENDIF

   RETURN
   END SUBROUTINE TXPBLA
 
!****************************************************************
!
    CHARACTER*20 FUNCTION NUMSTR ( IVAL, RVAL, FORM )
!
!****************************************************************
!
    USE OCPCOMM4                                                        
    USE schism_glbl, ONLY:rkind
!
    IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Convert integer or real to string with given format
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!     FORM        given format
!     RVAL        real to be converted
!
    INTEGER   IVAL
    REAL(rkind)      RVAL
    CHARACTER FORM*20
!
!  6. Local variables
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
    IF(IVAL /= INAN)THEN
      WRITE(NUMSTR,FORM) IVAL
    ELSE IF( RVAL /= RNAN)THEN
      WRITE (NUMSTR,FORM) RVAL
    ELSE
      NUMSTR = ''
    END IF

    RETURN
    END FUNCTION NUMSTR
!****************************************************************
!
   SUBROUTINE SWCOPI ( IARR1, IARR2, LENGTH )
!
!****************************************************************
!
   USE OCPCOMM4                                                        
!
   IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Copies integer array IARR1 to IARR2
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     IARR1       source array
!     IARR2       target array
!     LENGTH      array length
!
   INTEGER LENGTH
   INTEGER IARR1(LENGTH), IARR2(LENGTH)
!
!  6. Local variables
!
!     I     :     loop counter
!     IENT  :     number of entries
!
   INTEGER I, IENT
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     STRACE           Tracing routine for debugging
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
!     Trivial.
!
! 13. Source text
!
   SAVE IENT
   DATA IENT/0/
   IF(LTRACE) CALL STRACE (IENT,'SWCOPI')

!  --- check array length

   IF(LENGTH <= 0)THEN
     CALL MSGERR( 3, 'Array length should be positive' )
   END IF

!  --- copy elements of array IARR1 to IARR2

   DO I = 1, LENGTH
     IARR2(I) = IARR1(I)
   END DO

   RETURN
   END SUBROUTINE SWCOPI
 

!***********************************************************************
!                                                                      *
    SUBROUTINE KSCIP1(MMT,SIG,D,K,CG,N,ND)
!                                                                      *
!***********************************************************************
!
    USE OCPCOMM4                                                        
    USE SWCOMM3                                                         
    USE schism_glbl, ONLY:rkind
!
    IMPLICIT NONE
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Calculation of the wave number, group velocity, group number N
!     and the derivative of N w.r.t. depth (=ND)
!
!  3. Method
!
!     --
!
!  4. Argument variables
!
!     MMT     input    number of frequency points
!
      INTEGER   MMT
!
!     CG      output   group velocity
!     D       input    local depth
!     K       output   wave number
!     N       output   ratio of group and phase velocity
!     ND      output   derivative of N with respect to D
!     SIG     input    rel. frequency for which wave parameters
!                      must be determined
!
      REAL(rkind)      CG(MMT), D,                                    &
                       K(MMT), N(MMT), ND(MMT), SIG(MMT)
!
!  6. Local variables
!
!     C         phase velocity
!     FAC1      auxiliary factor
!     FAC2      auxiliary factor
!     FAC3      auxiliary factor
!     IENT      number of entries
!     IS        counter in frequency (sigma-space)
!     KND       dimensionless wave number
!     ROOTDG    square root of D/GRAV
!     WGD       square root of GRAV*D
!     SND       dimensionless frequency
!     SND2      = SND*SND
!
      INTEGER   IENT, IS
      REAL(rkind)      KND, ROOTDG, SND, WGD, SND2, C, FAC1, FAC2, FAC3
!
!  8. Subroutines used
!
!     --
!
!  9. Subroutines calling
!
!     SWOEXA, SWOEXF (Swan/Output)
!
! 10. Error messages
!
!     --
!
! 11. Remarks
!
!     --
!
! 12. Structure
!
!     -----------------------------------------------------------------
!      Compute non-dimensional frequency SND
!      IF SND >= 2.5, then
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND according to
!        deep water theory
!      ELSE IF SND =< 1.e-6
!        Compute wave number K, group velocity CGO, ratio of group
!        and phase velocity N and its derivative ND
!        according to extremely shallow water
!      ELSE
!        Compute wave number K, group velocity CGO and the ratio of
!        group and phase velocity N by Pade and other simple formulas.
!        Compute the derivative of N w.r.t. D = ND.
!     -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT /0/
      IF (LTRACE) CALL STRACE (IENT, 'KSCIP1')
!
      ROOTDG = SQRT(D/GRAV_W)                                               !30.81
      WGD    = ROOTDG*GRAV_W                                                !30.81
      DO 200 IS = 1, MMT
!       SND is dimensionless frequency
        SND = SIG(IS) * ROOTDG
        IF (SND .GE. 2.5) THEN
!       ******* deep water *******
          K(IS)  = SIG(IS) * SIG(IS) / GRAV_W                             !30.81
          CG(IS) = 0.5 * GRAV_W / SIG(IS)                                 !30.81
          N(IS)  = 0.5
          ND(IS) = 0.
        ELSE IF (SND.LT.1.E-6) THEN
!       *** very shallow water ***                                        30.81
          K(IS)  = SND/D                                                  !30.81
          CG(IS) = WGD
          N(IS)  = 1.
          ND(IS) = 0.
        ELSE
          SND2  = SND*SND                                                 !40.41
          C     = SQRT(GRAV_W*D/(SND2+1./(1.+0.666*SND2+0.445*SND2**2 &   !40.41
     &                                  -0.105*SND2**3+0.272*SND2**4)))   !40.41
          K(IS) = SIG(IS)/C                                               !40.41
          KND   = K(IS)*D                                                 !40.41
          FAC1  = 2.*KND/SINH(2.*KND)                                     !40.41
          N(IS) = 0.5*(1.+FAC1)                                           !40.41
          CG(IS)= N(IS)*C                                                 !40.41
          FAC2  = SND2/KND                                                !40.41
          FAC3  = 2.*FAC2/(1.+FAC2*FAC2)                                  !40.41
          FAC2  = -K(IS)*(2.*N(IS)-1.)/(2.*D*N(IS))                       !41.16
          ND(IS)= FAC1*(0.5/D - K(IS)/FAC3 + FAC2*(0.5/K(IS) - D/FAC3))   !41.16 40.41
        ENDIF
  200 CONTINUE

    RETURN

    END SUBROUTINE KSCIP1
!***********************************************************************
!                                                                      *
    SUBROUTINE KSCIP2(S,D,CG)
!                                                                      *
!***********************************************************************
!
    USE OCPCOMM4
    USE SWCOMM3
    USE schism_glbl,ONLY:rkind
!
    IMPLICIT NONE
    REAL(rkind)      CG(MSC), D, K,S(MSC)
    INTEGER   IENT, IS
    REAL(rkind)      KND, ROOTDG, SND, WGD, SND2, C, FAC1,N

    ROOTDG = SQRT(D/GRAV_W)
    WGD    = ROOTDG*GRAV_W
!     SND is dimensionless frequency
     DO IS=1,MSC
      SND = S(IS) * ROOTDG
      IF(SND >= 2.5)THEN
!     ******* deep water *******
        K = S(IS) * S(IS) / GRAV_W
        CG(IS) = 0.5 * GRAV_W / S(IS)
      ELSE IF(SND < 1.E-6)THEN
!     *** very shallow water ***  
!        print*,'IN VERY SHALLOW WATER'
!       stop                                      
        K  = SND/D
        CG(IS) = WGD
      ELSE
        SND2  = SND*SND
        C     = SQRT(GRAV_W*D/(SND2+1./(1.+0.666*SND2+0.445*SND2**2      &
                -0.105*SND2**3+0.272*SND2**4)))
        K = S(IS)/C
        KND   = K*D
        FAC1  = 2.*KND/SINH(2.*KND)
        N = 0.5*(1.+FAC1)
        CG(IS)= N*C
      ENDIF
    END DO
    RETURN
    END SUBROUTINE KSCIP2



!****************************************************************

!****************************************************************
!
      CHARACTER*20 FUNCTION INTSTR ( IVAL )
!
!****************************************************************
!
      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: M. Zijlema                                    |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!     40.23: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Feb. 03: New subroutine
!
!  2. Purpose
!
!     Convert integer to string
!
!  4. Argument variables
!
!     IVAL        integer to be converted
!
      INTEGER IVAL
!
!  6. Local variables
!
!     CVAL  :     character represented an integer of mantisse
!     I     :     counter
!     IPOS  :     position in mantisse
!     IQUO  :     whole quotient
!
      INTEGER I, IPOS, IQUO
      CHARACTER*1, ALLOCATABLE :: CVAL(:)
!
! 12. Structure
!
!     Trivial.
!
! 13. Source text
!
      IPOS = 1
 100  CONTINUE
      IF (IVAL/10**IPOS.GE.1.) THEN
         IPOS = IPOS + 1
         GO TO 100
      END IF
      ALLOCATE(CVAL(IPOS))

      DO I=IPOS,1,-1
         IQUO=IVAL/10**(I-1)
         CVAL(IPOS-I+1)=CHAR(INT(IQUO)+48)
         IVAL=IVAL-IQUO*10**(I-1)
      END DO

      WRITE (INTSTR,*) (CVAL(I), I=1,IPOS)

      RETURN
      END FUNCTION INTSTR

!
!********************************************************************
!                                                                   *
      REAL(rkind) FUNCTION GAMMA(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)
!
      USE schism_glbl,ONLY:rkind

      REAL(rkind) XX, YY, ABIG                                                   
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30._rkind/
      CALL STRACE (IENT, 'GAMMA')
      YY =  REAL( GAMMLN( XX ) )
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMA = MyEXP(YY)
      RETURN
      END FUNCTION GAMMA

!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION TCROSS (X1, X2, X3, X4, Y1, Y2, Y3, Y4, X1ONOBST)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        !40.41
!
      IMPLICIT NONE                                                       !40.04
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

!     40.00  Gerbrant van Vledder
!     40.04  Annette Kieftenburg
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!       30.70, Feb 98: argument list simplified
!                      subroutine changed into logical function
!       40.00, Aug 98: division by zero prevented
!       40.04, Aug 99: method corrected, IMPLICIT NONE added, XCONOBST added,
!                      introduced TINY and EPSILON (instead of comparing to 0)
!                      replaced 0 < LMBD,MIU by  0 <= LMBD,MIU
!                      XCONOBST added to argument list
!       40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!       Find if there is an obstacle crossing the stencil in used
!
!  3. Method
!
!     For the next situation (A, B and C are the points in the stencil,
!     D and E  are corners of the obstacle
!
!
!      obstacle --> D(X3,Y3)
!                    *
!                     *
!                      *
!        (X2,Y2)        * (XC,YC)
!            B-----------@--------------------------A (X1,Y1)
!                        ^*                         /
!                   _____| *                       /
!                  |        *                     /
!                  |         *                   /
!         crossing point      *                 /
!                              *               /
!                               E             /
!                              (X4,Y4)       /
!                                           C
!
!
!       The crossing point (@) should be found solving the next eqs.
!       for LMBD and MIU.
!
!       | XC |    | X1 |           | X2 - X1 |
!       |    | =  |    | +  LMBD * |         |
!       | YC |    | Y1 |           | Y2 - Y1 |
!
!
!       | XC |    | X3 |           | X4 - X3 |                            40.04
!       |    | =  |    | +  MIU  * |         |
!       | YC |    | Y3 |           | Y4 - Y3 |                            40.04
!
!
!     If solution exist and (0 <= LMBD <= 1 and 0 <= MIU <= 1)            40.04
!     there is an obstacle crossing the stencil
!
!  4. Argument variables
!
!     X1, Y1  inp    user coordinates of one end of grid link
!     X2, Y2  inp    user coordinates of other end of grid link
!     X3, Y3  inp    user coordinates of one end of obstacle side
!     X4, Y4  inp    user coordinates of other end of obstacle side
!     X1ONOBST outp   boolean which tells whether (X1,Y1) is on obstacle
!
      REAL(rkind)       EPS, X1, X2, X3, X4, Y1, Y2, Y3, Y4
      LOGICAL    X1ONOBST
!
!  5. Parameter variables
!
!  6. Local variables
!
!     A,B,C,D    dummy variables
!     DIV1       denominator of value of LMBD (or MIU)
!     E,F        dummy variables
!     IENT       number of entries
!     LMBD       coefficient in vector equation for stencil points (or obstacle)
!     MIU        coefficient in vector equation for obstacle (or stencil points)
!
      INTEGER    IENT
      REAL(rkind)       A, B, C, D, DIV1, E, F, LMBD, MIU
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOBST
!     SWTRCF
!     OBSTMOVE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     Calculate MIU and LMBD
!     If 0 <= MIU, LMBD <= 1                                              40.04
!     Then TCROSS is .True.
!     Else TCROSS is .False.
!
! 13. Source text
! ======================================================================
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'TCROSS')
!
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))            ! 40.04
      IF (EPS ==0.) EPS = TINY(X1)                                       ! 40.04
      A    = X2 - X1
!     A not equal to zero
      IF (ABS(A) .GT. TINY(X1)) THEN                                     ! 40.04
        B    = X4 - X3
        C    = X3 - X1
        D    = Y2 - Y1
        E    = Y4 - Y3
        F    = Y3 - Y1
      ELSE
!       exchange MIU and LMBD                                             40.04
        A    = X4 - X3
        B    = X2 - X1
        C    = X1 - X3
        D    = Y4 - Y3
        E    = Y2 - Y1
        F    = Y1 - Y3
      ENDIF
      DIV1 = ((A*E) - (D*B))                                             ! 40.00
!
!     DIV1 = 0 means that obstacle is parallel to line through            40.04
!     stencil points, or (X3,Y3) = (X4,Y4);                               40.04
!     A = 0 means trivial set of equations X4= X3 and X2 =X1              40.04
!
      IF ((ABS(DIV1).LE.TINY(X1)) .OR.&                                 ! 40.04
     &    (ABS(A).LE.TINY(X1))) THEN                                    ! 40.04
        MIU = -1.                                                       ! 40.00
        LMBD = -1.                                                      ! 40.04
      ELSE                                                              ! 40.00
        MIU  = ((D*C) - (A*F)) / DIV1
        LMBD = (C + (B*MIU)) / A
      END IF                                                            ! 40.00
!
      IF (MIU  .GE. 0. .AND. MIU  .LE. 1. .AND.&                        ! 40.04
     &    LMBD .GE. 0. .AND. LMBD .LE. 1.) THEN                         ! 40.04
!
!       Only (X1,Y1) is checked, because of otherwise possible double     40.04
!       counting                                                          40.04
        IF ((LMBD.LE.EPS .AND. ABS(X2-X1).GT.EPS).OR.&                   ! 40.04
     &      (MIU .LE.EPS .AND. ABS(X2-X1).LE.EPS))THEN                   ! 40.04
          X1ONOBST = .TRUE.                                              ! 40.04
        ELSE                                                             ! 40.04
          X1ONOBST = .FALSE.                                             ! 40.04
        ENDIF                                                            ! 40.04
!
!       *** test output ***
        IF (ITEST .GE. 120) THEN
          WRITE(PRINTF,70)X1,Y1,X2,Y2,X3,Y3,X4,Y4
  70      FORMAT(' Obstacle crossing  :',/,&
     &    ' Coordinates of comp grid points and corners of obstacle:',/,&
     &    ' P1(',F12.2,',',F12.2,')',' P2(',F12.2,',',F12.2,')',/,&
     &    ' P3(',F12.2,',',F12.2,')',' P4(',F12.2,',',F12.2,')')
        ENDIF
!
        TCROSS = .TRUE.
      ELSE
        TCROSS = .FALSE.
      ENDIF
!
!     End of subroutine TCROSS
      RETURN
      END


!***********************************************************************
!                                                                      *
      SUBROUTINE SWTRCF (WLEV2 , CHS   ,&                                 !40.31 40.00
     &                   LINK  , OBREDF,&                                 !40.03
!     &                   AC2   , REFLSO, KGRPNT, XCGRID,&                 !40.41 40.09
     &                            REFLSO, KGRPNT, XCGRID,&                 !40.41 40.09
     &                   YCGRID, CAX,    CAY,    RDX,    RDY,    ANYBIN,& !40.09
     &                   SPCSIG, SPCDIR)                                  !40.13 40.28
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        !40.41
      USE SWCOMM2                                                         !40.66
      USE SWCOMM3                                                         !40.41
      USE SWCOMM4                                                         !40.41
      USE M_OBSTA                                                         !40.31
!      USE M_PARALL                                                        !40.31
      USE SwanGriddata                                                    !40.80
      USE SwanCompdata                                                    !40.80

      USE VARS_WAVE, ONLY : AC2
      USE schism_glbl,ONLY:rkind

      IMPLICIT NONE                                                       !40.09
!
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
!     30.70
!     40.03  Nico Booij
!     40.08  Erick Rogers
!     40.09  Annette Kieftenburg
!     40.13  Nico Booij
!     40.14  Annette Kieftenburg
!     40.18  Annette Kieftenburg
!     40.28  Annette Kieftenburg
!     40.30  Marcel Zijlema
!     40.31  Marcel Zijlema
!     40.41: Marcel Zijlema
!     40.66: Marcel Zijlema
!     40.80: Marcel Zijlema
!     41.65: Marcel Zijlema
!
!  1. Updates
!
!     30.70, Feb. 98: water level (WLEV2) replaced depth
!                     incident wave height introduced using argument
!                     CHS (sign. wave height in whole comput. grid)
!     40.03, Jul. 00: LINK1 and LINK2 in argumentlist replaced by LINK
!     40.09, Nov. 99: IMPLICIT NONE added, Method corrected
!                     Reflection option for obstacle added
!     40.14, Dec. 00: Reflection call corrected: reduced to neighbouring
!                     linepiece of obstacle (bug fix 40.11D)
!            Jan. 01: Constant waterlevel taken into account as well (bug fix
!            40.11E)
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added
!     40.13, Aug. 02: subroutine restructured:
!                     loop in reflection procedure changed to avoid double
!                     reflection
!                     argument list of subr REFLECT revised
!                     argument SPCDIR added
!     40.30, Mar. 03: correcting indices of test point with offsets MXF, MYF
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent
!                     with other subroutines
!     40.31, Oct. 03: changes w.r.t. obstacles
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!     40.66, Mar. 07: extension with d'Angremond and Van der Meer transmission
!     40.80, Mar. 08: extension to unstructured grids
!     41.65, Jun. 16: extension frequency and direction dependent tranmission
!     coefficients
!
!  2. Purpose
!
!      take the value of transmission coefficient given
!      by the user in case obstacle TRANSMISSION
!
!      or
!
!      compute the transmision coeficient in case obstacle DAM
!      based on Goda (1967) [from Seelig (1979)]
!      or d'Angremond and Van der Meer formula's (1996)                   40.66
!
!      if reflections are switched on, calculate sourceterm in            40.09
!      subroutine REFLECT                                                 40.09
!
!  3. Method
!
!     Calculate transmission coefficient based on Goda (1967)             40.09
!     from Seelig (1979)                                                  40.09
!     Kt = 0.5*(1-sin {pi/(2*alpha)*(WATHIG/Hi +beta)})
!     where
!     Kt         transmission coefficient
!
!     alpha,beta coefficients dependent on structure of obstacle
!                and waves
!     WATHIG     = F = h-d is the freeboard of the dam, where h is the    40.09
!                crest level of the dam above the reference level and d   40.09
!                is the mean water level (relative to reference level)    40.09
!     Hi         incident (significant) wave height                       40.09
!                                                                         40.09
!     If reflection are switched on and obstacle is not exactly on line   40.09
!     of two neighbouring gridpoints, calculate reflections               40.09
!
!  4. Argument variables
!
!     AC2      input     Action density array                             40.09
!     ANYBIN   input     Set a particular bin TRUE or FALSE depending on  40.09
!                        SECTOR                                           40.09
!     CAX      input     Propagation velocity                             40.09
!     CAY      input     Propagation velocity                             40.09
!     CHS      input     Hs in all computational grid points
!     KCGRD    input     Grid address of points of computational stencil
!     LINK     input     indicates whether link in stencil                40.03
!                        crosses an obstacle                              40.03
!     OBREDF   output    Array of action density reduction coefficients
!                        (reduction at the obstacle)
!     REFLSO   inp/outp  contribution to the source term of action        40.41
!                        balance equation due to reflection               40.41
!     RDX,RDY  input     Array containing spatial derivative coefficients 40.09
!     WLEV2    input     Water level in grid points
!
      INTEGER  KGRPNT(MXC,MYC)
      INTEGER  LINK(2)                                                    !40.03
      REAL(rkind)     CHS(MCGRD), OBREDF(MDC,MSC,2), WLEV2(MCGRD)                !30.70
!      REAL     :: AC2(MDC,MSC,MCGRD)                                      !40.09 40.22
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   40.22
      REAL(rkind)     :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)                !40.09 40.22
      REAL(rkind)     :: REFLSO(MDC,MSC), RDX(MICMAX), RDY(MICMAX)               !40.41 40.09 40.22 40.08
      REAL(rkind)     :: SPCSIG(MSC), SPCDIR(MDC,6)                              !40.13 40.28
      LOGICAL  :: ANYBIN(MDC,MSC)                                         !40.09

!  5. Parameter variables
!
!  6. Local variables
!
!     ALOW     Lower limit for FVH
!     BK       crest width                                                40.66
!     BUPL     Upper limit for FVH
!     BVH      Bk/Hsin                                                    40.66
!     FD1      Coeff. for freq. dep. reflection: vertical displacement    40.28
!     FD2      Coeff. for freq. dep. reflection: shape parameter          40.28
!     FD3      Coeff. for freq. dep. reflection: directional coefficient  40.28
!     FD4      Coeff. for freq. dep. reflection: bending point of freq.   40.28
!     FVH      WATHIG/Hsin
!     HGT      elevation of top of obstacle above reference level
!     HSIN     incoming significant wave height
!     ID       counter in directional space
!     IENT     number of entries of this subroutine
!     ILINK    indicates which link is analyzed: 1 -> neighbour in x
!                                                2 -> neighbour in y
!     IS       counter in frequency space
!     ITRAS    indicates kind of obstacle: 0 -> constant transm
!                                          1 -> dam, Goda
!                                          2 -> dam, d'Angremond and      40.66
!                                                    Van der Meer         40.66
!     JP       counter for number of corner points of obstacles
!     L0P      wave length in deep water                                  40.66
!     LREFDIFF indicates whether reflected energy should be               40.18
!              scattered (1) or not (0)                                   40.18
!     LREFL    if LREFL=0, no reflection; if LREFL=1, constant            40.13
!              reflection coeff.                                          40.13
!     LRFRD    Indicates whether frequency dependent reflection is        40.28
!              active (#0.) or not (=0.)                                  40.28
!     NMPO     link number
!     NUMCOR   number of corner points of obstacle
!     OBET     user defined coefficient (beta) in formulation of
!              Goda/Seelig (1967/1979)
!     OBHKT    transmission coefficient in terms of wave height
!     OGAM     user defined coefficient (alpha) in formulation of
!              Goda/Seelig (1967/1979)
!     POWN     user defined power of redistribution function              40.28
!     REFLCOEF reflection coefficient in terms of action density
!     REFLTST  used to test Refl^2+Transm^2 <=1                           40.13
!     SLOPE    slope of obstacle                                          40.66
!     SQRTREF  dummy variable
!     TRCF     transmission coefficient in terms of action density
!              (user defined or calculated (in terms of waveheight))
!     X1, Y1   user coordinates of one end of grid link
!     X2, Y2   user coordinates of other end of grid link
!     X3, Y3   user coordinates of one end of obstacle side
!     X4, Y4   user coordinates of other end of obstacle side
!     XCGRID   Coordinates of computational grid in x-direction
!     XI0P     breaker parameter                                          40.66
!     XOBS     x-coordinate of obstacle point                             40.80
!     XONOBST  Indicates whether computational point (X1,Y1) is on        40.14
!              obstacle                                                   40.14
!     XV       x-coordinate of vertex of face                             40.80
!     YCGRID   Coordinates of computational grid in y-direction
!     YOBS     y-coordinate of obstacle point                             40.80
!     YV       y-coordinate of vertex of face                             40.80
!     WATHIG   freeboard of the dam (= HGT-waterlevel)
!
      INTEGER    ID, IENT, ILINK, ITRAS, IS, JP, ICGRD, LREFL,&
     &           NUMCOR, NMPO, ISIGM
      INTEGER    LREFDIFF, LRFRD                                          !40.31
      REAL(rkind)       ALOW, BUPL, FVH, HGT, HSIN, OBET, OBHKT,&
     &           SLOPE, BK, L0P, XI0P, BVH, EMAX, ETD, TP,&                !40.66
     &           FAC1, FAC2,&
     &           POWN, OGAM, REFLCOEF,&                                    !40.18 40.09
     &           X1, X2, X3, X4, Y1, Y2, Y3, Y4, WATHIG
      REAL(rkind)       TRCF(MSC,MDC)                                            !41.65
      REAL(rkind)       FD1, FD2, FD3, FD4                                       !40.31 40.28
      REAL(rkind)       SQRTREF                                                  !40.09
      LOGICAL    XONOBST                                                  !40.14
      LOGICAL :: REFLTST                                                  !40.13
      REAL(rkind)       XCGRID(MXC,MYC), YCGRID(MXC,MYC)
      INTEGER    ICC, JJ
      REAL(rkind)    :: XOBS(2), XV(2), YOBS(2), YV(2)                           !40.80
      LOGICAL :: SwanCrossObstacle                                        !40.80
      TYPE(OBSTDAT), POINTER :: COBST                                     !40.31
!
!  8. Subroutines used
!
!     MSGERR           Writes error message
!     REFLECT          Computes effect of reflection
!     TCROSS           Searches for crossing point if exist               40.14
!
      LOGICAL TCROSS                                                      !40.14
!
!  9. Subroutines calling
!
!     SWOMPU                                                              30.70
!
! 10. Error messages
!
! 11. Remarks
!
!     Here the formulation of the transmission coefficients concerns the  40.09
!     ratio of action densities!                                          40.09
!
! 12. Structure
!
!     ------------------------------------------------------------------
!     For both links from grid point (X1,Y1) do                           40.13
!         calculate transmission coefficients                             40.09
!         assign values to OBREDF                                         40.13
!         If there is reflection                                          40.13
!         Then select obstacle side which crosses the grid link           40.13
!              calculate reflection source terms                          40.09
!     ------------------------------------------------------------------  40.13
!
! 13. Source text
! ======================================================================
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWTRCF')
!
      REFLTST = .TRUE.                                                    !40.13
      DO 90 ILINK = 1 ,2
!       default transmission coefficient
        TRCF = 1.
        NMPO = LINK(ILINK)
        IF (NMPO .EQ. 0) THEN
           OBREDF(1:MDC,1:MSC,ILINK) = 1.                                 !40.13
           GOTO 90
        END IF
!       incoming wave height
        HSIN = CHS(KCGRD(ILINK+1))                                        !40.03
        IF (HSIN.LT.0.1E-4) HSIN = 0.1E-4
        COBST => FOBSTAC                                                  !40.31
        DO JJ = 1, NMPO-1                                                 !40.31
           IF (.NOT.ASSOCIATED(COBST%NEXTOBST)) EXIT                      !40.31
           COBST => COBST%NEXTOBST                                        !40.31
        END DO                                                            !40.31
        ITRAS  = COBST%TRTYPE                                             !40.31
        IF (ITRAS .EQ. 0) THEN
!       constant transmission coefficient
!         User defined transmission coefficient concerns ratio of         !40.09
!         wave heights, so                                                !40.09
          OBHKT   = COBST%TRCOEF(1)                                       !40.31
          TRCF(1,1) = OBHKT * OBHKT                                       !40.09
        ELSE IF (ITRAS .EQ. 1) THEN
!       transmission coefficient according to Goda and Seelig
          HGT    =  COBST%TRCOEF(1)                                       !40.31
          OGAM   =  COBST%TRCOEF(2)                                       !40.31
          OBET   =  COBST%TRCOEF(3)                                       !40.31
!         level of dam above the water surface (freeboard):
          WATHIG =  HGT - WLEV2(KCGRD(1)) - WLEV                          !40.14 30.70
!
!         *** Here the transmission coeff. is that of Goda and Seelig ***
          FVH  = WATHIG/HSIN
          ALOW = -OBET-OGAM                                               !40.09
          BUPL = OGAM-OBET
!
          IF (FVH.LT.ALOW) FVH = ALOW
          IF (FVH.GT.BUPL) FVH = BUPL
          OBHKT = 0.5*(1.0-SIN(PI_W*(FVH+OBET)/(2.0*OGAM)))
! JL Skip Now
# if 0
          IF (TESTFL) WRITE (PRTEST, 20) IXCGRD(1)+MXF-2,&                 !40.30
     &           IYCGRD(1)+MYF-2,&                                         !40.30
     &           ILINK, HGT, WATHIG, HSIN, OBHKT                          !40.01
  20      FORMAT (' test SWTRCF ', 2X, 3I5, ' dam level=', F6.2,&
     &            ' depth=', F6.2, ' Hs=', F6.2, ' transm=', F6.3)        !40.01
# endif
          IF (TESTFL .AND. ITEST.GE.140) WRITE (PRTEST, 22)&
     &           OGAM, OBET, ALOW, BUPL, FVH                              !40.01
  22      FORMAT (8X, 5E12.4)                                             !40.01
!
!         Formulation of Goda/Seelig concerns ratio of waveheights.       !40.09
!         Here we use action density, so                                  !40.09
          TRCF(1,1) = OBHKT * OBHKT
        ELSE IF (ITRAS.EQ.2) THEN                                         !40.66
!       d'Angremond and Van der Meer formulae (1996)                      !40.66
          HGT   = COBST%TRCOEF(1)
          SLOPE = COBST%TRCOEF(2)
          BK    = COBST%TRCOEF(3)

!         level of dam above the water surface:
          WATHIG =  HGT - WLEV2(KCGRD(1)) - WLEV
!
!         compute peak frequency of incoming wave
          EMAX = 0.
          ISIGM = -1
          DO IS = 1, MSC
             ETD = 0.
             DO ID = 1, MDC
                ETD = ETD + SPCSIG(IS)*AC2(ID,IS,KCGRD(ILINK+1))*DDIR
             END DO
             IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
             END IF
          END DO
          IF (ISIGM.LE.0) ISIGM=MSC
          TP=2.*PI_W/SPCSIG(ISIGM)
!
!         compute breaker parameter
          L0P  = MAX(1.E-8,1.5613*TP*TP)
          XI0P = TAN(SLOPE*PI_W/180.)/SQRT(HSIN/L0P)
!
!         compute transmission coefficient
          FVH = WATHIG/HSIN
          BVH = BK/HSIN
          IF (BVH.EQ.0.) THEN
             OBHKT= -0.40*FVH
             IF (OBHKT.LT.0.075) OBHKT = 0.075
             IF (OBHKT.GT.0.900) OBHKT = 0.9
          ELSE IF (BVH.LT.8.) THEN
             OBHKT= -0.40*FVH + 0.64*(BVH**(-0.31))*(1.-EXP(-0.50*XI0P))
             IF (OBHKT.LT.0.075) OBHKT = 0.075
             IF (OBHKT.GT.0.900) OBHKT = 0.9
          ELSE IF (BVH.GT.12.) THEN
             OBHKT= -0.35*FVH + 0.51*(BVH**(-0.65))*(1.-EXP(-0.41*XI0P))
             IF (OBHKT.GT.0.93-0.006*BVH) OBHKT = 0.93-0.006*BVH
             IF (OBHKT.LT.0.05          ) OBHKT = 0.05
          ELSE
!            linear interpolation
             FAC1 = -0.40*FVH + 0.64*( 8.**(-0.31))*(1.-EXP(-0.50*XI0P))
             IF (FAC1.LT.0.075) FAC1 = 0.075
             IF (FAC1.GT.0.900) FAC1 = 0.9
             FAC2 = -0.35*FVH + 0.51*(12.**(-0.65))*(1.-EXP(-0.41*XI0P))
             IF (FAC2.LT.0.050) FAC2 = 0.050
             IF (FAC2.GT.0.858) FAC2 = 0.858
             OBHKT = 3.*FAC1 - 2.*FAC2 + BVH*(FAC2 - FAC1)/4.
          END IF
          IF(OPTG.NE.5) THEN
            X1 = XCGRID(IXCGRD(1),IYCGRD(1))+XOFFS
            Y1 = YCGRID(IXCGRD(1),IYCGRD(1))+YOFFS
          ELSE
            X1 = xcugrd(vs(1))+XOFFS
            Y1 = ycugrd(vs(1))+YOFFS
          ENDIF
          WRITE (PRINTF, 23) X1, Y1,&
     &       ILINK, HGT, WATHIG, HSIN, SQRT(L0P/1.5613), XI0P, OBHKT      !40.66
  23      FORMAT (' Transmission: ', 2X, 2F12.4, I5, ' dam level=',F6.2,&
     &            ' board=', F6.2, ' Hs=', F6.2, ' Tp=', F6.2,&
     &            ' Xi0p=', F6.3, ' Kt=', F6.3)                           !40.66
          IF (TESTFL .AND. ITEST.GE.140) WRITE (PRTEST, 22)&
     &           TAN(SLOPE*PI_W/180.), FVH, BVH, L0P, XI0P
!
!         Formulation of d'Angremond concerns ratio of waveheights
!         Here we use action density, so
          TRCF(1,1) = OBHKT * OBHKT
        ELSE IF (ITRAS .EQ. 11) THEN                                      !41.65
!       frequency dependent transmission coefficients                     41.65
!       user defined values concerns ratio of wave heights, so            41.65
          DO IS = 1, MSC
             OBHKT      = COBST%TRCF1D(IS)                                !41.65
             TRCF(IS,1) = OBHKT * OBHKT                                   !41.65
          ENDDO
        ELSE IF (ITRAS .EQ. 12) THEN                                      !41.65
!       frequency and direction dependent transmission coefficients       41.65
!       user defined values concerns ratio of wave heights, so            41.65
          DO ID = 1, MDC
             DO IS = 1, MSC
                OBHKT       = COBST%TRCF2D(ID,IS)                         !41.65
                TRCF(IS,ID) = OBHKT * OBHKT                               !41.65
             ENDDO
          ENDDO
        ENDIF
!
!       assign values to array OBREDF
        IF (ITRAS.EQ.11) THEN                                             !41.65
           DO IS = 1, MSC
              DO ID = 1, MDC
                 OBREDF(ID,IS,ILINK) = TRCF(IS,1)
              ENDDO
           ENDDO
        ELSE IF (ITRAS.EQ.12) THEN
           DO IS = 1, MSC
              DO ID = 1, MDC
                 OBREDF(ID,IS,ILINK) = TRCF(IS,ID)
              ENDDO
           ENDDO
        ELSE
           DO IS = 1, MSC
              DO ID = 1, MDC
                 OBREDF(ID,IS,ILINK) = TRCF(1,1)
              ENDDO
           ENDDO
        ENDIF
!
!     *** REFLECTION ****
!     *** X1 X2 X3 ETC. are the coordinates of point according ***
!     *** with the scheme in the function TCROSS header        ***        !40.04
        LREFL = COBST%RFTYP1                                              !40.13
        IF ( LREFL.GT.0. ) THEN                                           !40.31
!         Reflections are activated                                       !40.09
          SQRTREF  = COBST%RFCOEF(1)                                      !40.31
          REFLCOEF = SQRTREF * SQRTREF                                    !40.09
          LREFDIFF = COBST%RFTYP2                                         !40.31
          POWN     = COBST%RFCOEF(2)                                      !40.31
          FD1      = COBST%RFCOEF(3)                                      !40.31
          FD2      = COBST%RFCOEF(4)                                      !40.31
          FD3      = COBST%RFCOEF(5)                                      !40.31
          FD4      = COBST%RFCOEF(6)                                      !40.31
          LRFRD    = COBST%RFTYP3                                         !40.31
!
          IF (OPTG.NE.5) THEN                                             !40.80
!            determine (X1,Y1) and (X2,Y2)
             ICC   = KCGRD(1)                                             !40.09
             ICGRD = 0                                                    !40.09
             IF ( ICC.GT.1 ) THEN                                         !40.09
               X1 = XCGRID(IXCGRD(1),IYCGRD(1))                           !40.09
               Y1 = YCGRID(IXCGRD(1),IYCGRD(1))                           !40.09
               IF (KGRPNT(IXCGRD(ILINK+1),IYCGRD(ILINK+1)).GT.1) THEN     !40.09
                 X2    = XCGRID(IXCGRD(ILINK+1),IYCGRD(ILINK+1))          !40.09
                 Y2    = YCGRID(IXCGRD(ILINK+1),IYCGRD(ILINK+1))          !40.09
                 ICGRD = KCGRD(ILINK+1)                                   !40.09
               ENDIF
             ENDIF
             IF (ICGRD.EQ.0) GOTO 90                                      !40.13
!            select obstacle side crossing the grid link
             X3 = COBST%XCRP(1)                                           !40.31
             Y3 = COBST%YCRP(1)                                           !40.31
             NUMCOR = COBST%NCRPTS                                        !40.31
             DO JP = 2, NUMCOR                                            !40.09
               X4 = COBST%XCRP(JP)                                        !40.31
               Y4 = COBST%YCRP(JP)                                        !40.31
               IF (TCROSS(X1,X2,X3,X4,Y1,Y2,Y3,Y4,XONOBST)) GOTO 70       !40.14
               X3 = X4                                                    !40.09
               Y3 = Y4                                                    !40.09
             ENDDO
          ELSE                                                            !40.80
!           determine begin and end points of link                        40.80
            X1 = xcugrd(vs(1))                                            !40.80
            Y1 = ycugrd(vs(1))                                            !40.80
            X2 = xcugrd(vs(ILINK+1))                                      !40.80
            Y2 = ycugrd(vs(ILINK+1))                                      !40.80
            XV(1) = X1                                                    !40.80
            YV(1) = Y1                                                    !40.80
            XV(2) = X2                                                    !40.80
            YV(2) = Y2                                                    !40.80
!           select obstacle side crossing the grid link                   !40.80
            X3 = COBST%XCRP(1)                                            !40.80
            Y3 = COBST%YCRP(1)                                            !40.80
            XOBS(1) = X3                                                  !40.80
            YOBS(1) = Y3                                                  !40.80
            DO JP = 2, COBST%NCRPTS                                       !40.80
               X4 = COBST%XCRP(JP)                                        !40.80
               Y4 = COBST%YCRP(JP)                                        !40.80
               XOBS(2) = X4                                               !40.80
               YOBS(2) = Y4                                               !40.80
               IF ( SwanCrossObstacle( XV, YV, XOBS, YOBS ) ) GOTO 70     !40.80
               X3 = X4                                                    !40.80
               Y3 = Y4                                                    !40.80
               XOBS(1) = X3                                               !40.80
               YOBS(1) = Y3                                               !40.80
            ENDDO                                                         !40.80
          ENDIF                                                           !40.80
!         no crossing found, skip procedure
          GOTO 90

! 70      CALL REFLECT(AC2, REFLSO, X1, Y1, X2, Y2,&                      ! 40.13
 70      CALL REFLECT(     REFLSO, X1, Y1, X2, Y2,&                      ! 40.13
     &                 X3, Y3, X4, Y4, CAX,&                              ! 40.13
     &                 CAY, RDX, RDY, ILINK,&                             ! 40.13
     &                 REFLCOEF, LREFDIFF, POWN, ANYBIN,&                 ! 40.13
     &                 LRFRD, SPCSIG, SPCDIR, FD1, FD2, FD3, FD4,&        ! 40.13
     &                 OBREDF, REFLTST)                                   !40.13
        END IF
!
        IF (ITEST .GE. 120)  WRITE (PRTEST,10)&
     &  IXCGRD(1)-1, IYCGRD(1)-1, NMPO, TRCF(1,1)
  10    FORMAT(' SWTRCF: Point=', 2I5, ' NMPO  = ', I5, ' transm ',&      !40.03
     &  F8.3)                                                             !40.03
  90  CONTINUE
      IF (.NOT.REFLTST) THEN                                              !40.13
        CALL MSGERR(3,'Kt^2 + Kr^2 > 1 ')                                 !40.13
        IF (ITEST.LT.50) THEN                                             !40.13
          WRITE (PRTEST, 74) IXCGRD(1)-1, IYCGRD(1)-1                     !40.13
  74      FORMAT (' Kt^2 + Kr^2 > 1 in grid point:', 2I4)                 !40.13
        ENDIF                                                             !40.13
      ENDIF                                                               !40.13
      RETURN
!     * end of SUBROUTINE SWTRCF
      END


!************************************************************************
!                                                                       *
!      SUBROUTINE REFLECT (AC2, REFLSO, X1, Y1, X2, Y2, X3, Y3,&            !40.13
      SUBROUTINE REFLECT (     REFLSO, X1, Y1, X2, Y2, X3, Y3,&            !40.13
     &                    X4, Y4, CAX, CAY, RDX, RDY,&                     !40.13
     &                    ILINK, REF0, LREFDIFF, POWN, ANYBIN,&            !40.13
     &                    LRFRD, SPCSIG, SPCDIR, FD1, FD2, FD3, FD4,&      !40.13
     &                    OBREDF, REFLTST)                                !40.13
!                                                                       *
!************************************************************************
!                                                                         !40.09
      USE OCPCOMM4                                                        !40.41
      USE SWCOMM3                                                         !40.41
      USE SWCOMM4                                                         !40.41

      USE VARS_WAVE, ONLY : AC2

      USE schism_glbl,ONLY:rkind
!                                                                         !40.09
      IMPLICIT NONE                                                       !40.09
!
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
!  0. Authors                                                             40.09
!                                                                         40.09
!     40.09  Annette Kieftenburg                                          40.09
!     40.13  Nico Booij
!     40.18  Annette Kieftenburg                                          40.18
!     40.28  Annette Kieftenburg                                          40.28
!     40.38  Annette Kieftenburg                                          40.38
!     40.41: Marcel Zijlema
!                                                                         40.09
!  1. Updates                                                             40.09
!                                                                         40.09
!     40.09, Nov. 99: Subroutine created                                  40.09
!     40.18, Apr. 01: Scattered reflection against obstacles added        40.18
!     40.28, Dec. 01: Frequency dependent reflection added                40.28
!     40.38, Feb. 02: Diffuse reflection against obstacles added          40.38
!     40.13, Sep. 02: Subroutine restructured
!                     assumptions changed: reflected energy can come
!                     from Th_norm-PI/2 to Th_norm+PI/2
!     40.08, Mar. 03: Dimensioning of RDX, RDX changed to be consistent   40.08
!                     with other subroutines                              40.08
!     40.13, Nov. 03: test on refl + transm added
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!                                                                         40.09
!  2. Purpose                                                             40.09
!                                                                         40.09
!     Computation of REFLECTIONS near obstacles                           40.09
!                                                                         40.09
!  3. Method                                                              40.09
!                                                                         40.09
!     Determine the angle of the obstacle,                                40.09
!     Determine the angles between which reflections should be taken      40.09
!     into account                                                        40.09
!     Determine redistribution function                                   40.18
!     determine expression of reflection coefficient for frequency
!     dependency, if appropriate                                          40.28
!     Determine reflected action density (corrected for angle obstacle    40.09
!     and if option is on: redistribute energy)                           40.18
!     Add reflected spectrum to contribution for the right hand side      40.41
!     40.09
!     of matrix equation                                                  40.09
!                                                                         40.09
!  4. Modules used                                                        40.18
!                                                                         40.18
!     --                                                                  40.18
!                                                                         40.18
!  5. Argument variables                                                  40.09
!                                                                         40.09
!     AC2      inp  action density
!     ANYBIN   inp  Determines whether a bin fall within a sweep
!     CAX      inp  Propagation velocity in x-direction                   40.09
!     CAY      inp  Propagation velocity in y-direction                   40.09
!     FD1      inp  Coeff. freq. dep. reflection: vertical displacement   40.28
!     FD2      inp  Coeff. freq. dep. reflection: shape parameter         40.28
!     FD3      inp  Coeff. freq. dep. reflection: directional coefficient 40.28
!     FD4      inp  Coeff. freq. dep. reflection: bending point of freq.  40.28
!     ILINK    inp  Indicates which link is analyzed: 1 -> neighbour in x 40.09
!                                                     2 -> neighbour in y 40.09
!     LREFDIFF inp  Indicates whether reflected energy should be          40.18
!                   scattered (1) or not (0)                              40.18
!     LRFRD    inp  Indicates whether frequency dependent reflection is   40.28
!                   active (#0.) or not (=0.)                             40.28
!     OBREDF   inp  transmission coefficients
!     POWN     inp  User defined power of redistribution function         40.18
!     REF0     inp  reflection coefficient in terms of action density
!     REFLSO   i/o  contribution to the source term due to reflection     40.41
!     REFLTST  i/o  used to test Refl^2+Transm^2 <=1                      40.13
!     RDX,RDY  inp  Array containing spatial derivative coefficients      40.09
!     SPCDIR(*,1)   spectral directions (radians)
!     SPCSIG   inp  Relative frequency (= 2*PI*Freq.)                     40.28
!     X1, Y1   inp  Coordinates of computational grid point under         40.09
!                   consideration                                         40.09
!     X2, Y2   inp  Coordinates of computational grid point neighbour     40.09
!     X3, Y3   inp  User coordinates of one end of obstacle side          40.09
!     X4, Y4   inp  User coordinates of other end of obstacle side        40.09
!                                                                         40.09
!      REAL       :: AC2(MDC,MSC,MCGRD)                                    !40.18
!     Changed ICMAX to MICMAX, since MICMAX doesn't vary over gridpoint   !40.22
      REAL(rkind)       :: CAX(MDC,MSC,MICMAX), CAY(MDC,MSC,MICMAX)              !40.18 40.22
      REAL(rkind)       :: REFLSO(MDC,MSC), OBREDF(MDC,MSC,2)                    !40.41 40.18
      REAL(rkind)       :: RDX(MICMAX), RDY(MICMAX)                              !40.18 40.08
      REAL(rkind)       :: FD1, FD2, FD3, FD4, SPCSIG(MSC), SPCDIR(MDC,6)        !40.41 40.28
      REAL(rkind)       :: REF0                                                  !40.18
      REAL(rkind)       :: X1, X2, X3, X4, Y1, Y2, Y3, Y4                        !40.18
      LOGICAL    :: ANYBIN(MDC,MSC)                                       !40.18
      INTEGER    :: ILINK                                                 !40.18
      REAL(rkind)       :: POWN                                                  !40.41 40.18
      INTEGER    :: LREFDIFF, LRFRD                                       !40.41
      LOGICAL    :: REFLTST                                               !40.13
!                                                                         !40.18
!      INTENT (IN)     AC2, CAX, CAY, OBREDF,&                              !40.18
      INTENT (IN)          CAX, CAY, OBREDF,&                              !40.18
     &                FD1, FD2, FD3, FD4, LRFRD,&                          !40.28
     &                RDX, RDY, SPCSIG, SPCDIR, X1, X2,&                   !40.18
     &                X3, X4, Y1, Y2, Y3, Y4, ANYBIN, ILINK,&              !40.18
     &                POWN, LREFDIFF                                      !40.18
      INTENT (IN OUT) REF0, REFLSO                                        !40.41 40.18
!                                                                         !40.09
!  6. Parameter variables                                                 !40.09
!                                                                         !40.18
!  7. Local variables                                                     !40.09
!                                                                         !40.09
      REAL(rkind) :: AC2REF     ! reflected action density of one spectral bin
      REAL(rkind) :: BETA         ! local angle of obstacle                      !40.09
      REAL(rkind) :: EPS                                                         !40.09
      REAL(rkind), ALLOCATABLE :: PRDIF(:)   ! scattering filter                 !40.13
      REAL(rkind)    :: TH_INC         ! direction of incident wave
      REAL(rkind)    :: TH_NORM        ! direction of normal to obstacle
      REAL(rkind)    :: TH_OUT         ! direction of outgoing wave
      REAL(rkind)    :: IANG           ! angle divided by DTheta (DDIR)
      REAL(rkind)    :: SUMRD          ! sum of PRDIF array
      REAL(rkind)    :: W1, W2         ! interpolation coefficients

      INTEGER :: ID             ! counter of directions
      INTEGER :: IS             ! counter of frequencies
      INTEGER :: MAXIDR         ! width of scattering filter
      INTEGER :: IDR            ! relative directional counter
      INTEGER :: ID_I1, ID_I2   ! counters of incoming directions
      INTEGER :: IDA, IDB       ! counters of incoming directions
      INTEGER :: IENT           ! number of entries
!                                                                   
!                                                                         40.18
!  8. Subroutines used                                                    40.09
!                                                                         40.09
!  9. Subroutines calling                                                 40.09
!                                                                         40.09
!     SWTRCF                                                              40.09
!                                                                         40.09
! 10. Error messages                                                      40.09
!                                                                         40.09
!     if obstacle linepiece is of length < EPS                            40.09
!                                                                         40.09
! 11. Remarks                                                             40.09
!                                                                         40.09
!    -In case the obstacle cuts exactly through computational grid point, 40.09
!     the obstacle should be moved a bit with subroutine OBSTMOVE.        40.09
!    -The length of the obstacle linepiece is assumed to be               40.09
!     'long enough' compared to grid resolution (> 0.5*sqrt(dx^2+dy^2))   40.09
!     (if this restriction is violated, the reflections due to an obsta-  40.09
!     cle of one straight line can be very different from a similar line  40.09
!     consisting of several pieces (because only the directions of the    40.09
!     spectrum that are directed towards the obstacle linepiece are       40.09
!     reflected).                                                         40.09
!    -There should be only one intersection per computational gridcell.   40.09
!     Therefore it is better to avoid sharp edges in obstacles.           40.18
!                                                                         40.09
! 12. Structure                                                           40.09
!                                                                         40.09
!     -----------------------------------------------------------------
!     Determine angle of obstacle, Beta                                   40.09
!     Determine angle of normal from (X1,Y1) to obstacle                  40.13
!     If there is constant diffuse reflection
!     Then determine scattering distribution
!     -----------------------------------------------------------------
!     For all frequencies do
!         If amount of reflection varies with frequency
!         Then determine reflection coefficient
!         -------------------------------------------------------------
!         If there is diffuse reflection
!         Then if scattering varies with frequency
!              Then determine power of cos
!              --------------------------------------------------------
!              Determine distribution
!         -------------------------------------------------------------
!         For all active directions do
!             Determine specular incoming direction
!             For directions of scattering filter do
!                 multiply incoming action with scattering coefficient
!                 add this to array AC2REF
!     -----------------------------------------------------------------
!                                                                         !40.09
!     Add reflected spectrum to right hand side of matrix equation        !40.09
!                                                                         !40.09
! 13. Source text                                                         !40.09
!                                                                         !40.18
      SAVE     IENT                                                       !40.09
      DATA     IENT /0/                                                   !40.09
      CALL STRACE (IENT, 'REFLECT')                                       !40.09
!                                                                         !40.09
      EPS = EPSILON(X1)*SQRT((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))             !40.09
      IF (EPS ==0) EPS = TINY(X1)                                         !40.09

      IF (LREFDIFF .EQ. 0) THEN
        ALLOCATE (PRDIF(0:0))
        MAXIDR = 0
        PRDIF(0) = 1.
      ELSE
        MAXIDR = MDC/2
        ALLOCATE (PRDIF(0:MDC/2))
      ENDIF
!                                                                         !40.09
!     Determine angle of obstacle BETA, and related                       !40.09
!                                                                         !40.09
      IF (.NOT. ((ABS(Y4-Y3).LE.EPS) .AND. (ABS(X4-X3).LE.EPS)) ) THEN    !40.09
        BETA = ATAN2((Y4-Y3),(X4-X3))                                     !40.09
      ELSE                                                                !40.09
        CALL MSGERR (3, 'Obstacle contains line piece of length 0!')      !40.09
      END IF                                                              !40.09
!     determine direction of normal                         (4)
!     this is the normal from (X1,Y1)                        |
!     towards the obstacle                                   |
!     see sketch to establish sign                    (2)----+------(1)
!                                                            |
!                                                           (3)
      IF ((X1-X3)*(Y4-Y3)-(Y1-Y3)*(X4-X3) .GT. 0.) THEN                   !40.13
        TH_NORM = BETA + 0.5*PI_W
      ELSE
        TH_NORM = BETA - 0.5*PI_W
      ENDIF                                                               !40.13

!     prepare directional filter in case of diffuse reflection
      IF (LREFDIFF.EQ.1) THEN
        PRDIF(1:MAXIDR) = 0.
        PRDIF(0) = 1.
        SUMRD = 1.
        DO ID = 1, MAXIDR
          PRDIF(ID) = (COS(ID*DDIR))**POWN
          IF (PRDIF(ID) .GT. 0.01) THEN
            SUMRD = SUMRD + 2.*PRDIF(ID)
          ELSE
            MAXIDR = ID-1
            GOTO 23
          ENDIF
        ENDDO
  23    DO ID = 0, MAXIDR
          PRDIF(ID) = PRDIF(ID) / SUMRD
        ENDDO
        IF (TESTFL .AND. ITEST.GE.50) THEN
          WRITE (PRTEST, 27) POWN, MAXIDR
  27      FORMAT (' power scattering filter:', F4.1, I3)
          IF (ITEST.GE.130) WRITE (PRTEST, 28)&
     &    (PRDIF(IDR), IDR=0, MAXIDR)
  28      FORMAT (10 F7.3)
        ENDIF
      ENDIF

      DO IS = 1, MSC
        IF (LRFRD.EQ.1) THEN
!         amount of reflection varies with wave frequency
          REF0 = FD1 + &                                                  !40.13
     &           FD2/PI_W * ATAN2(PI_W*FD3*(SPCSIG(IS)-FD4),FD2)              !40.13
          IF (REF0 > 1.) REF0 = 1.                                        !40.28
          IF (REF0 < 0.) REF0 = 0.                                        !40.28
!         >>> should REF0 not be squared? <<<
        ENDIF
!       check whether reflection + transmission <= 1                      !40.13
        DO ID = 1, MDC                                                    !40.13
          IF ((REF0 + OBREDF(ID,IS,ILINK)) .GT. 1.) THEN                  !40.13
            REFLTST = .FALSE.                                             !40.13
            IF (ITEST.GE.50) THEN                                         !40.13
              WRITE (PRTEST, 72) IXCGRD(1)-1, IYCGRD(1)-1, ILINK, &       !40.13
     &               IS, ID, REF0, OBREDF(ID,IS,ILINK)                    !40.13
  72          FORMAT (' Refl+Transm>1 in ', 2I4, 2X, 3I3, 2X, 2F6.2)      !40.13
            ENDIF                                                         !40.13
          ENDIF                                                           !40.13
        ENDDO                                                             !40.13

        IF (LREFDIFF.EQ.2) THEN
!         spreading varies with frequency; not yet implemented
        ENDIF
        DO ID = 1, MDC
          IF (ANYBIN(ID,IS)) THEN
            AC2REF = 0.
            TH_OUT = SPCDIR(ID,1)
!           corresponding incident direction (assuming specular reflection)
            TH_INC = 2.*BETA-TH_OUT
!           determine counter for which direction is TH_INC:
            IANG = MOD (TH_INC-SPCDIR(1,1), 2.*PI_W) / DDIR
!           incident angle is between ID_I1 and ID_I2
            ID_I1 = 1 + INT (IANG)
            ID_I2 = ID_I1+1
!           W1 and W2 are weighting coefficients for the above directions
!           by linear interpolation
            W2 = IANG + 1. - REAL(ID_I1)
            W1 = 1. - W2
            DO IDR = -MAXIDR, MAXIDR
              IDA = ID_I1 + IDR
              IDB = ID_I2 + IDR
              IF (FULCIR) THEN
                IDA = 1+MOD(2*MDC+IDA-1,MDC)
!               only outgoing reflected waves, i.e. not towards obstacle
                IF (COS(TH_NORM-SPCDIR(IDA,1)) .GT. 0.)&
     &          AC2REF = AC2REF +&
     &          REF0 * W1 * PRDIF(ABS(IDR)) * AC2(IDA,IS,KCGRD(1))
                IDB = 1+MOD(2*MDC+IDB-1,MDC)
                IF (COS(TH_NORM-SPCDIR(IDB,1)) .GT. 0.)&
     &          AC2REF = AC2REF +&
     &          REF0 * W2 * PRDIF(ABS(IDR)) * AC2(IDB,IS,KCGRD(1))
              ELSE
                IF (IDA.GE.1 .AND. IDA.LE.MDC) THEN
                  IF (COS(TH_NORM-SPCDIR(IDA,1)) .GT. 0.)&
     &            AC2REF = AC2REF +&
     &            REF0 * W1 * PRDIF(ABS(IDR)) * AC2(IDA,IS,KCGRD(1))
                ENDIF
                IF (IDB.GE.1 .AND. IDB.LE.MDC) THEN
                  IF (COS(TH_NORM-SPCDIR(IDB,1)) .GT. 0.)&
     &            AC2REF = AC2REF +&
     &            REF0 * W2 * PRDIF(ABS(IDR)) * AC2(IDB,IS,KCGRD(1))
                ENDIF
              ENDIF
            ENDDO
!           add reflected energy to right hand side of matrix
            REFLSO(ID,IS) = REFLSO(ID,IS) + AC2REF *&
     &      (RDX(ILINK)*CAX(ID,IS,1) + RDY(ILINK)*CAY(ID,IS,1))           !40.13
          END IF                                                          !40.09
        END DO                                                            !40.09
      END DO                                                              !40.09
!                                                                         !40.18
      IF (LREFDIFF .GT. 0) DEALLOCATE (PRDIF)
      RETURN                                                              !40.09
!     End of subroutine REFLECT                                           !40.09
      END                                                                 !40.09

!
!***********************************************************************
!                                                                      *
      SUBROUTINE SSHAPE (ACLOC, SPCSIG, SPCDIR, FSHAPL, DSHAPL)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                        
      USE SWCOMM3 

      USE schism_msgp, ONLY : parallel_abort                                                        
      USE schism_glbl, ONLY : skind,rkind
!      IMPLICIT NONE
!
!  0. Authors
!
!  1. Updates
!
!  2. Purpose
!
!     Calculating of energy density at boundary point (x,y,sigma,theta)
!
!  3. Method (updated...)
!
!     see: M. Yamaguchi: Approximate expressions for integral properties
!          of the JONSWAP spectrum; Proc. JSCE, No. 345/II-1, pp. 149-152,
!          1984.
!
!     computation of mean period: see Swan system documentation
!
!  4. Argument variables
!
!   o ACLOC : Energy density at a point in space
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
!
      REAL(rkind)    ACLOC(MDC,MSC)
      REAL(rkind)    SPCDIR(MDC,6)                                               
      REAL(rkind)    SPCSIG(MSC)                                                 
!
! i   DSHAPL: Directional distribution
! i   FSHAPL: Shape of spectrum:
!             =1; Pierson-Moskowitz spectrum
!             =2; Jonswap spectrum
!             =3; bin
!             =4; Gauss curve
!             (if >0: period is interpreted as peak per.
!              if <0: period is interpreted as mean per.)
!
      INTEGER FSHAPL, DSHAPL                                              
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!     LSHAPE   absolute value of FSHAPL
!
      INTEGER  ID, IS, LSHAPE
!
!     PKPER    peak period                                                30.80
!     APSHAP   aux. var. used in computation of spectrum
!     AUX1     auxiliary variable
!     AUX2     auxiliary variable
!     AUX3     auxiliary variable
!     COEFF    coefficient for behaviour around the peak (Jonswap)
!     CPSHAP   aux. var. used in computation of spectrum
!     CTOT     total energy
!     CTOTT    total energy (used for comparison)
!     DIFPER   auxiliary variable used to select bin closest
!              to given frequency
!     MPER
!     MS       power in directional distribution
!     RA       action density
!     SALPHA
!     SF       frequency (Hz)
!     SF4      SF**4
!     SF5      SF**5
!     FPK      frequency corresponding to peak period (1/PKPER)           30.80
!     FPK4     FPK**4
!     SYF      peakedness parameter
!
      REAL(rkind)     APSHAP, AUX1, AUX2, AUX3
      REAL(rkind)     COEFF ,SYF   ,MPER  ,CTOT  ,CTOTT,PKPER  ,DIFPER
      REAL(rkind)     MS
      REAL(rkind)     RA    ,SALPHA,SF   ,SF4   ,SF5   ,FPK   ,FPK4, FAC

!     JL Not declared vars, fixed now
      REAL(rkind)    PPSHAP,ESOM,ACOS,HSTMP,AM0,AM1,AS2,AS3,PPTAIL,APTAIL,EPTAIL,ADIR,CDIR,DSPR,CPSHAP
      INTEGER        ITPER,JJ,ISP

      REAL(rkind) DEGCNV        ! Function DEGCNV
!
!     LOGPM    indicates whether peak or mean frequency is used
!     DVERIF   logical used in verification of incident direction
!
      LOGICAL  LOGPM, DVERIF                                              
!
!     PSHAPE   coefficients of spectral distribution (see remarks)
!     SPPARM   array containing integral wave parameters (see remarks)
!
!  8. Subroutines used
!
!     ---
!
!  9. Subroutines calling
!
! 10. Error messages
!
! 11. Remarks
!
!     PSHAPE(1): SY0, peak enhancement factor (gamma) in Jonswap spectrum
!     PSHAPE(2): spectral width in case of Gauss spectrum in rad/s
!
!     SPPARM    real     input    incident wave parameters (Hs, Period,
!                                 direction, Ms (dir. spread))
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by the user (either peak or mean)      30.80
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!
!     ---------------------------------------------------------------------
!
!     In the case of a JONSWAP spectrum the initial conditions are given by
!                   _               _       _       _       _
!                  |       _   _ -4  |     |       | S - S   |
!             2    |      |  S  |    |     |       |      p  |
!          a g     |      |  _  |    |  exp|-1/2 * |________ |* 2/pi COS(T-T  )
! E(S,D )= ___  exp|-5/4 *|  S  |    | G   |       | e * S   |              wi
!      wa    5     |      |   p |    |     |_      |_     p _|
!           S      |      |_   _|    |
!                  |_               _|
!
!   where
!         S   : rel. frequency
!
!         D   : Dir. of wave component
!          wa
!
!         a   : equili. range const. (Phillips' constant)
!         g   : gravity acceleration
!
!         S   : Peak frequency
!          p
!
!         G   : Peak enhancement factor
!         e   : Peak width
!
!         T   : local wind direction
!          wi
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       case shape
!       =1:   calculate value of Pierson-Moskowitz spectrum
!       =2:   calculate value of Jonswap spectrum
!       =3:   calculate value of bin spectrum
!       =4:   calculate value of Gauss spectrum
!       else: Give error message because of wrong shape
!       ----------------------------------------------------------------
!       if LOGPM is True
!       then calculate average period
!            if it differs from given average period
!            then recalculate peak period
!                 restart procedure to compute spectral shape
!       ----------------------------------------------------------------
!       for all spectral bins do
!            multiply all action densities by directional distribution
!       ----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE(IENT,'SSHAPE')
!
      IF (ITEST >= 80) WRITE (PRTEST, 8) FSHAPL, DSHAPL,          &
                                         (SPPARM(JJ), JJ = 1,4)
   8  FORMAT (' entry SSHAPE ', 2I3, 4E12.4)
      IF(FSHAPL < 0)THEN
        LSHAPE = - FSHAPL
        LOGPM  = .FALSE.
      ELSE
        LSHAPE = FSHAPL
        LOGPM  = .TRUE.
      ENDIF

!
      IF(SPPARM(1) <= 0._rkind)                                         &
         CALL MSGERR(1,' sign. wave height at boundary is not positive')  
!
      PKPER = SPPARM(2)
      ITPER = 0
      IF(LSHAPE == 3)THEN
!       select bin closest to given period
        DIFPER = 1.E10
        DO IS = 1, MSC
          IF(MyABS(PKPER - PI2_W/SPCSIG(IS)) < DIFPER)THEN
            ISP = IS
            DIFPER = MyABS(PKPER - PI2_W/SPCSIG(IS))
          ENDIF
        ENDDO
      ENDIF
!
!     compute spectral shape using peak period PKPER                      
!
      FAC  = 1._rkind
 100  FPK  = (1._rkind/PKPER)                                                   
      FPK4 = FPK**4._rkind
      IF(LSHAPE == 1)THEN
        SALPHA = ((SPPARM(1) ** 2._rkind) * (FPK4)) * 5.0_rkind / 16.0_rkind
      ELSE IF(LSHAPE == 2)THEN
!       *** SALPHA = alpha*(grav**2)/(2.*pi)**4)
        SALPHA = (SPPARM(1)**2._rkind * FPK4) /                            &
	         ((0.06533_rkind*(PSHAPE(1)**0.8015_rkind)+0.13467_rkind)*16.0_rkind)
      ELSE IF(LSHAPE == 4)THEN
        AUX1 = SPPARM(1)**2 / ( 16.0_rkind* MySQRT (PI2_W) * PSHAPE(2))
        AUX3 = 2.0_rkind * PSHAPE(2)**2._rkind
      ENDIF

      IF(ITEST >= 80) &
      WRITE (PRTEST, *) 'SALPHA,HS,PKPER,FPK4,PSHAPE',SALPHA,SPPARM(1),PKPER,FPK4,PSHAPE
!
      CTOTT = 0.0_rkind
      DO IS = 1, MSC                                                  
!
        IF(LSHAPE == 1)THEN
!         *** LSHAPE = 1 : Pierson and Moskowitz ***
          SF = SPCSIG(IS) / PI2_W
          SF4 = SF**4._rkind
          SF5 = SF**5._rkind
          RA = (SALPHA/SF5)*MyEXP(-(5._rkind*FPK4)/(4._rkind*SF4))/(PI2_W*SPCSIG(IS))
          ACLOC(MDC,IS) = RA
! JL      ACLOC(MDC,IS) = MAX(0.0_skind,RA)

        ELSE IF(LSHAPE == 2)THEN
!         *** LSHAPE = 2 : JONSWAP ***
          SF = SPCSIG(IS)/(PI2_W)
          SF4 = SF**4._rkind
          SF5 = SF**5._rkind
          CPSHAP = 1.25_rkind * FPK4 / SF4
          IF(CPSHAP > 10.0_rkind)THEN                                         
            RA = 0.0_rkind
          ELSE
            RA = (SALPHA/SF5) * MyEXP(-CPSHAP)
          ENDIF
          IF(SF < FPK)THEN
            COEFF = 0.07_rkind
          ELSE
            COEFF = 0.09_rkind
          ENDIF
          APSHAP =  0.5_rkind * ((SF-FPK) / (COEFF*FPK)) **2
          IF(APSHAP > 10.0_rkind)THEN                                         
            SYF = 1.0_rkind
          ELSE
            PPSHAP = MyEXP(-APSHAP)
            SYF = PSHAPE(1)**PPSHAP
          ENDIF

          !IF(IS.GT.30) WRITE (PRTEST, *) 'SF,SF4,SF5,RA,SALPHA,COEFF,PPSHAP',SF,SF4,SF5,RA,SALPHA,COEFF,PPSHAP

          RA = SYF*RA/(SPCSIG(IS)*PI2_W)
          ACLOC(MDC,IS) = RA
! JL      ACLOC(MDC,IS) = MAX(0.0_skind,RA)

          IF(ITEST >= 120) WRITE (PRTEST, 112)                  &
	            SF, SALPHA, CPSHAP, APSHAP, SYF, RA
 112      FORMAT (' SSHAPE freq. ', 8E12.4)

        ELSE IF(LSHAPE == 3)THEN
!
!         *** all energy concentrated in one BIN ***
!
          IF(IS == ISP)THEN
            ACLOC(MDC,IS) = ( SPPARM(1)**2._rkind ) /                   &
	                    ( 16.0_rkind * SPCSIG(IS)**2._rkind * FRINTF )
! JL        ACLOC(MDC,IS) = MAX(0.0_skind,ACLOC(MDC,IS))
          ELSE
            ACLOC(MDC,IS) = 0.0_rkind
          ENDIF

        ELSE IF(LSHAPE == 4)THEN
!
!         *** energy Gaussian distributed (wave-current tests) ***
!
          AUX2 = ( SPCSIG(IS) - ( PI2_W / PKPER ) )**2._rkind
          RA = AUX1 * MyEXP ( -1.0_rkind * AUX2 / AUX3 ) / SPCSIG(IS)
          ACLOC(MDC,IS) = RA
! JL      ACLOC(MDC,IS) = MAX(0.0_skind,RA)
        ELSE
          IF(IS == 1)THEN
            CALL MSGERR (2,'Wrong type for frequency shape')
            WRITE (PRINTF, *) ' -> ', FSHAPL, LSHAPE
          ENDIF
        ENDIF
        IF (ITEST >= 10)                                         &
	        CTOTT = CTOTT + FRINTF * ACLOC(MDC,IS) * SPCSIG(IS)**2
      END DO   !END DO IS

      IF(ITEST >= 10)THEN
        IF(SPPARM(1) > 0.01_rkind)THEN
          HSTMP = 4.0_rkind * MySQRT(CTOTT)
          IF(MyABS(HSTMP-SPPARM(1)) > 0.1_rkind*SPPARM(1)) THEN         
           WRITE (PRINTF, 303) SPPARM(1), HSTMP
           CALL parallel_abort("SSHAPE, Unexpected Deviation in HS, skipping now")
          ENDIF
 303      FORMAT (' SSHAPE, deviation in Hs, should be ', F8.3,      &
                  ', calculated ', F8.3)
        ENDIF
      ENDIF
!
!     if mean frequency was given recalculate PKPER and restart
!
      IF (.NOT.LOGPM .AND. ITPER < 10)THEN
        ITPER = ITPER + 1
!       calculate average frequency
        AM0 = 0.0_rkind
        AM1 = 0.0_rkind
        DO IS = 1, MSC
          AS2 = ACLOC(MDC,IS) * (SPCSIG(IS))**2
          AS3 = AS2 * SPCSIG(IS)
          AM0 = AM0 + AS2
          AM1 = AM1 + AS3
        ENDDO
!       contribution of tail to total energy density
        PPTAIL = PWTAIL(1) - 1.0_rkind                                           
        APTAIL = 1.0_rkind / (PPTAIL * (1.0_rkind + PPTAIL * (FRINTH-1.)))              
        AM0 = AM0 * FRINTF + APTAIL * AS2                                 
        PPTAIL = PWTAIL(1) - 2.0_rkind                                           
        EPTAIL = 1.0_rkind / (PPTAIL * (1.0_rkind + PPTAIL * (FRINTH-1.)))              
        AM1 = AM1 * FRINTF + EPTAIL * AS3                                 
!       Mean period:
        IF( AM1 /= 0.0_rkind)THEN                                             
           MPER = PI2_W * AM0 / AM1
        ELSE                                                              
           CALL MSGERR(3, ' first moment is zero in calculating the')     
           CALL MSGERR(3, ' spectrum at boundary using param. bc.')       
        END IF                                                            

        IF(ITEST >= 80) WRITE (PRTEST, 72) ITPER, SPPARM(2), MPER,     &
	                                    PKPER
  72    FORMAT (' SSHAPE iter=', I2, '  period values:', 3F7.2)
        IF(MyABS(MPER-SPPARM(2)) > 0.01_rkind*SPPARM(2))THEN
!         modification suggested by Mauro Sclavo
          PKPER = (SPPARM(2) / MPER) * PKPER                              
          GOTO 100
        ENDIF
      ENDIF
!
      IF (ITPER >= 10)THEN
        CALL MSGERR(3, 'No convergence calculating the spectrum')         
        CALL MSGERR(3, 'at the boundary using parametric bound. cond.')   
      ENDIF
!
!     now introduce distribution over directions

      ADIR = PI_W * DEGCNV( SPPARM(3)) / 180.0_rkind
      IF(DSHAPL == 1)THEN
        DSPR = PI_W * SPPARM(4) / 180.0_rkind
        MS = MAX (DSPR**(-2._rkind) - 2.0_rkind, 1.0_rkind)
      ELSE
        MS = SPPARM(4)
      ENDIF
      IF(MS < 12.0_rkind)THEN
        CTOT = (2.0_rkind**MS) * (GAMMA(0.5_rkind*MS+1.0_rkind))**2_rkind / (PI_W * GAMMA(MS+1.))
      ELSE
        CTOT =  MySQRT (0.5_rkind*MS/PI_W) / (1.0_rkind - 0.25_rkind/MS)
      ENDIF
      IF(ITEST >= 100)THEN
        ESOM = 0.0_rkind
        DO IS = 1, MSC
          ESOM = ESOM + FRINTF * SPCSIG(IS)**2._rkind * ACLOC(MDC,IS)
        ENDDO
        WRITE (PRTEST, *) ' SSHAPE dir ', 4._rkind*MySQRT(MyABS(ESOM)),         &
	       SPPARM(1), CTOT, MS, GAMMA(0.5_rkind*MS+1._rkind), GAMMA(MS+1._rkind),   &
	       CTOT                                                     
      ENDIF
      DVERIF = .FALSE.
      CTOTT = 0.0_rkind
      DO ID = 1, MDC
        ACOS = MyCOS(SPCDIR(ID,1) - ADIR)
        IF(ACOS > 0.0_rkind)THEN
          CDIR = CTOT * MAX (ACOS**MS, 1.E-10)
          IF(.NOT.FULCIR)THEN
            IF(ACOS >= MyCOS(DDIR)) DVERIF = .TRUE.
          ENDIF
        ELSE
          CDIR = 0.0_rkind
        ENDIF
        IF(ITEST >= 10) CTOTT = CTOTT + CDIR * DDIR
        IF(ITEST >= 200) WRITE (PRTEST, 360) ID,SPCDIR(ID,1),CDIR
 360    FORMAT (' ID Spcdir Cdir: ',I3,3(1X,E11.4))
        DO IS = 1, MSC
          ACLOC(ID,IS) = CDIR * ACLOC(MDC,IS)
!JL             ACLOC(ID,IS) = MAX(0.0_skind,CDIR * ACLOC(MDC,IS))
        ENDDO
      ENDDO
      IF(ITEST >= 10)THEN
        IF(MyABS(CTOTT-1.) > 0.1) WRITE (PRINTF, 363) CTOTT
 363    FORMAT (' SSHAPE, integral of Cdir is not 1, but:', F6.3)
      ENDIF
      IF (.NOT.FULCIR .AND. .NOT.DVERIF)                   &
         CALL MSGERR (1, 'incident direction is outside sector')

      RETURN

      END SUBROUTINE SSHAPE

!***********************************************************************
!                                                                      *
      REAL(rkind) FUNCTION DEGCNV (DEGREE)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM3           
      USE schism_glbl, ONLY:rkind,dkind                                              

      IMPLICIT NONE
!
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering                              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmers: R.C. Ris, N. Booij,                          |
!     |              IJ.G. Haagsma, A.T.M.M. Kieftenburg,         |
!     |              M. Zijlema, E.E. Kriezi,                     |
!     |              R. Padilla-Hernandez, L.H. Holthuijsen       |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2004-2005  Delft University of Technology
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees from nautical to cartesian or vice versa.
!
!  3. METHOD
!
!       DEGCNV = 180 + dnorth - degree
!
!  4. PARAMETERLIST
!
!       DEGCNV      direction in cartesian or nautical degrees.
!       DEGREE      direction in nautical or cartesian degrees.
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!           Nautical convention           Cartesian convention
!
!                    0                             90
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!        270 --------+-------- 90       180 --------+-------- 0
!                    |                              |
!                    |                              |
!                    |                              |
!                    |                              |
!                   180                            270
!
!  9. STRUCTURE
!
!     ---------------------------------
!     IF (NAUTICAL DEGREES) THEN
!       CONVERT DEGREES
!     IF (DEGREES > 360 OR < 0) THEN
!       CORRECT DEGREES WITHIN 0 - 360
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      INTEGER :: IENT
      REAL(rkind) :: DEGREE
      
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'DEGCNV')
!
      IF ( BNAUT ) THEN
        DEGCNV = 180.0_rkind + DNORTH - DEGREE
      ELSE
        DEGCNV = DEGREE
      ENDIF
!
      IF (DEGCNV >= 360.0_rkind) THEN
        DEGCNV = MOD (DEGCNV, 360.0_rkind)
      ELSE IF (DEGCNV < 0.0_rkind) THEN
        DEGCNV = MOD (DEGCNV, 360.0_rkind) + 360.0_rkind
      ELSE
!       DEGCNV between 0 and 360; do nothing
      ENDIF
!
      RETURN
      END FUNCTION DEGCNV
!
!***********************************************************************
!                                                                      *
      REAL(rkind) FUNCTION ANGRAD (DEGREE)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM3                                                       !  40.41
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform degrees to radians
!
!  3. METHOD
!
!       ANGRAD = DEGREE * PI / 180
!
!  4. PARAMETERLIST
!
!       ANGRAD      radians
!       DEGREE      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[radian] = ANGLE[degrees} * PI / 180
!     ---------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGRAD')
!
      ANGRAD = DEGREE * PI_W / 180._rkind
!
!
!     *** end of subroutine ANGRAD ***
!
      RETURN
      END
!
!***********************************************************************
      REAL(rkind) FUNCTION ANGDEG (RADIAN)
!                                                                      *
!***********************************************************************
!
      USE SWCOMM3                                                        ! 40.41
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
!  1. UPDATE
!
!       SEP 1997: New for SWAN 32.01
!                 Cor van der Schelde - Delft Hydraulics
!       30.70, Feb. 98: test output suppressed (causes problem if subr
!                       is used during output
!
!  2. PURPOSE
!
!       Transform radians to degrees
!
!  3. METHOD
!
!       ANGDEG = RADIAN * 180 / PI
!
!  4. PARAMETERLIST
!
!       RADIAN      radians
!       ANGDEG      degrees
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       NONE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ---------------------------------
!     ANGLE[degrees] = ANGLE[radians} * 180 / PI
!     ---------------------------------
      SAVE IENT
      DATA IENT /0/
      CALL STRACE(IENT,'ANGDEG')
!
      ANGDEG = RADIAN * 180._rkind / PI_W
!
!
!     *** end of subroutine ANGDEG ***
!
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE HSOBND (AC2   ,SPCSIG,HSIBC ,KGRPNT)                    ! 40.00
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                       ! 40.41
      USE SWCOMM3                                                        ! 40.41
      USE M_PARALL                                                       ! 40.31
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
!     32.01: Roeland Ris
!     30.70: Nico Booij
!     40.00: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     32.01, Sep. 97: new for SWAN
!     30.72, Jan. 98: Changed number of elements for HSI to MCGRD
!     30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for
!     SWAN
!     30.70, Feb. 98: structure scheme corrected
!     40.00, Mar. 98: integration method changed (as in SNEXTI)
!                     structure corrected
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Compare computed significant wave height with the value of
!     the significant wave height as predescribed by the user. If
!     the values differ more than e.g. 10 % give an error message
!     and the gridpoints where the error has been located
!
!  3. Method
!
!     ---
!
!  4. Argument variables
!
!     SPCSIG: input  Relative frequencies in computational domain in      30.72
!                    sigma-space                                          30.72
!
      REAL    SPCSIG(MSC)                                                ! 30.72
!
!       REALS:
!       ------
!       AC2        action density
!       HSI        significant wave height at boundary (using SWAN
!                  resolution (has thus not to be equal to the WAVEC
!                  significant wave height )
!       ETOT       total energy in a gridpoint
!       DS         increment in frequency space
!       DDIR       increment in directional space
!       HSC        computed wave height after SWAN computation
!       EFTAIL     contribution of tail to spectrum
!
!       INTEGERS:
!       ---------
!       KGRPNT     values of grid indices
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       TRACE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ------------------------------------------------------------------
!     for all computational grid points do                                30.70
!         if HSI is non-zero
!         then compute Hs from action density array
!              if relative difference is large than HSRERR
!              then write error message
!    -------------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      REAL      AC2(MDC,MSC,MCGRD) ,HSIBC(MCGRD)                         ! 30.72
!
      REAL      ETOT  ,HSC                                               ! 40.00
!
      INTEGER   ID    ,IS     ,IX     ,IY    ,INDX
!
      LOGICAL   HSRR
!
      INTEGER   KGRPNT(MXC,MYC)
!
      SAVE IENT, HSRR
      DATA IENT/0/, HSRR/.TRUE./
      CALL STRACE (IENT, 'HSOBND')
!
!     *** initializing ***

      HSRR = .TRUE.                                                      ! 40.03
!
      DO IY = MYC, 1, -1
        DO IX = 1, MXC
          INDX = KGRPNT(IX,IY)
          IF ( HSIBC(INDX) .GT. 1.E-25 ) THEN
!           *** compute Hs for boundary point (without tail) ***
            ETOT  = 0.
            DO ID = 1, MDC
              DO IS = 1, MSC                                             ! 40.00
                ETOT = ETOT + SPCSIG(IS)**2 * AC2(ID,IS,INDX)            ! 40.00
              ENDDO
            ENDDO
            IF (ETOT .GT. 1.E-8) THEN
              HSC = 4. * SQRT(ETOT*FRINTF*DDIR)
            ELSE
              HSC = 0.
            ENDIF
            HSREL = ABS(HSIBC(INDX) - HSC) / HSIBC(INDX)                 ! 40.51 40.00
            IF (HSREL .GT. HSRERR) THEN
              IF ( HSRR ) THEN
                WRITE (PRINTF,*) ' ** WARNING : ',&
     &             'Differences in wave height at the boundary'
                WRITE (PRINTF,802) HSRERR
 802            FORMAT (' Relative difference between input and ',&
     &          'computation >= ', F6.2)
                WRITE (PRINTF,*) '                        Hs[m]',&
     &                           '      Hs[m]      Hs[-]'
                WRITE (PRINTF,*) '    ix    iy  index   (input)',&
     &                           ' (computed) (relative)'
                WRITE (PRINTF,*) ' ----------------------------',&
     &                           '----------------------'
                HSRR = .FALSE.
              ENDIF
              WRITE (PRINTF,'(2(1x,I5),I7,3(1x,F10.2))')&
     &                 IX+MXF-1, IY+MYF-1, INDX, HSIBC(INDX), HSC, HSREL
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      WRITE(PRINTF,*)
!
      IF ( ITEST .GE. 150 ) THEN
        WRITE(PRINTF,*) 'Values of wave height at boundary (HSOBND)'
        WRITE(PRINTF,*) '------------------------------------------'
        DO IY = MYC, 1, -1
          WRITE (PRINTF,'(13F8.3)') ( HSIBC(KGRPNT(IX,IY)), IX=1 , MXC)
        ENDDO
      ENDIF
!
!     *** end of subroutine HSOBND ***
!
      RETURN
      END
!*****************************************************************
!                                                                *
      SUBROUTINE CHGBAS (X1, X2, PERIOD, Y1, Y2, N1, N2,&
     &                   ITEST, PRTEST)
!                                                                *
!*****************************************************************
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  G. van Vledder, N. Booij                       |
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
!  0. Update history
!
!       ver 20.48: also accomodates periodic variables such as directions
!
!  1. Purpose
!
!       change x-basis of a discretized y-function
!
!  2. Method
!
!     A piecewise constant representation of the functions is assumed
!
!     first boundaries of a cell in X1 are determined
!     then it is determined whether there are overlaps with cells
!     in X2. if so Y1*common length is added to Y2
!     Finally Y2 values are divided by cell lengths
!
!  3. Parameter list
!
!     Name    I/O  Type  Description
!
!     X1       i    ra   x-coordinates of input grid
!     X2       i    ra   x-coordinates of output grid
!     PERIOD   i    r    period, i.e. x-axis is periodic if period>0
!                        e.g. spectral directions
!     Y1       i    ra   function values of input grid
!     Y2       o    ra   function values of output grid
!     N1       i    i    number of x-values of input grid
!     N2       i    i    number of x-values of output grid
!
!  4. Subroutines used
!
!     ---
!
!  5. Error messages
!
!  6. Remarks
!
!       Cell boundaries in X1 are: X1A and X1B
!       X2 is assumed to be monotonically increasing; this is checked
!       X1 is assumed to be monotonous but not necessarily increasing
!
!  7. Structure
!
!       ------------------------------------------------------------------
!       Make all values of Y2 = 0
!       For each cell in X1 do
!           determine boundaries of cell in X1
!           --------------------------------------------------------------
!           For each cell in X2 do
!               determine overlap with cell in X1; limits: RLOW and RUPP
!               add to Y2: Y1 * length of overlapping interval
!       ------------------------------------------------------------------
!       For each cell in X2 do
!           divide Y2 value by cell length
!       ------------------------------------------------------------------
!
!  8. Source text
!
      INTEGER  I1, I2, N1, N2, ITEST, PRTEST
      REAL     X1(N1), Y1(N1), X2(N2), Y2(N2), PERIOD
      REAL     X1A, X1B, X2A, X2B, RLOW, RUPP
      LOGICAL  TWICE
      SAVE IENT
      DATA IENT /0/
      CALL STRACE (IENT, 'CHGBAS')
!
!     initialize output data
!
      DO I2 = 1, N2
        Y2(I2) = 0.
      ENDDO
      DO I2 = 2, N2
        IF (X2(I2).LE.X2(I2-1))&
     &    CALL MSGERR (2, 'subr. CHGBAS: values of X2 not increasing')
      ENDDO
!     boundaries of the range in X2
      X2LO  = 1.5 * X2(1)  - 0.5 * X2(2)
      X2HI  = 1.5 * X2(N2) - 0.5 * X2(N2-1)
      TWICE = .FALSE.
!
!     loop over cells in X1
!
      DO 300 I1 = 1, N1
        IF (ABS(Y1(I1)) .LT. 1.E-20) GOTO 300
!
!       determine cell boundaries in X1
!
        IF (I1.EQ.1) THEN
          X1A = 1.5 * X1(1) - 0.5 * X1(2)
        ELSE
          X1A = 0.5 * (X1(I1) + X1(I1-1))
        ENDIF

        IF (I1.EQ.N1) THEN
          X1B = 1.5 * X1(N1) - 0.5 * X1(N1-1)
        ELSE
          X1B = 0.5 * (X1(I1) + X1(I1+1))
        ENDIF
!
!       swap X1A and X1B if X1A > X1B
!
        IF (X1A.GT.X1B) THEN
          RR  = X1A
          X1A = X1B
          X1B = RR
        ENDIF

        IF (PERIOD.LE.0.) THEN
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
        ELSE
!         X is periodic; move interval in X1 if necessary
          TWICE = .FALSE.
          IADD = 0
  60      IF (X1B.GT.X2HI) THEN
            X1A = X1A - PERIOD
            X1B = X1B - PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)&
     &         CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 60
          ENDIF
  70      IF (X1A.LT.X2LO) THEN
            X1A = X1A + PERIOD
            X1B = X1B + PERIOD
            IADD = IADD + 1
            IF (IADD.GT.99)&
     &           CALL MSGERR (2, 'endless loop in CHGBAS')
            GOTO 70
          ENDIF
          IF (X1A.GT.X2HI) GOTO 300
          IF (X1B.LT.X2LO) GOTO 300
          IF (X1A.LT.X2LO .AND. X1A+PERIOD.LT.X2HI) TWICE = .TRUE.
          IF (X1B.GT.X2HI .AND. X1B-PERIOD.GT.X2LO) TWICE = .TRUE.
        ENDIF
!
!       loop over cells in X2
!
 100    DO 200 I2 = 1, N2

          IF (I2.EQ.1) THEN
            X2A = X2LO
          ELSE
            X2A = 0.5 * (X2(I2) + X2(I2-1))
          ENDIF

          IF (I2.EQ.N2) THEN
            X2B = X2HI
          ELSE
            X2B = 0.5 * (X2(I2) + X2(I2+1))
          ENDIF
!
!         (RLOW,RUPP) is overlapping interval of (X1A,X1B) and (X2A,X2B)
!
          IF (X1A.LT.X2B) THEN
            RLOW = MAX (X1A, X2A)
          ELSE
            GOTO 200
          ENDIF

          IF (X1B.GT.X2A) THEN
            RUPP = MIN (X1B, X2B)
          ELSE
            GOTO 200
          ENDIF

          IF (RUPP.LT.RLOW) THEN
            CALL MSGERR (3, 'interpolation error')
            WRITE (PRTEST, 140) I1, X1A, X1B, I2, X2A, X2B
 140        FORMAT (' I, XA, XB ', 2(I3, 2(1X,E12.4)))
          ELSE
            Y2(I2) = Y2(I2) + Y1(I1) * (RUPP-RLOW)
          ENDIF
 200    CONTINUE
!
!       Cell in X1 covers both ends of sector boundary
        IF (TWICE) THEN
          IF (X1A.LT.X2LO) THEN
             X1A = X1A + PERIOD
             X1B = X1B + PERIOD
          ENDIF
          IF (X1B.GT.X2HI) THEN
             X1A = X1A - PERIOD
             X1B = X1B - PERIOD
          ENDIF
          TWICE = .FALSE.
          GOTO 100
        ENDIF
 300  CONTINUE
!
      DO I2 = 1, N2
        IF (I2.EQ.1) THEN
          CELLEN = X2(2) - X2(1)
        ELSE IF (I2.EQ.N2) THEN
          CELLEN = X2(N2) - X2(N2-1)
        ELSE
          CELLEN = 0.5 * (X2(I2+1) - X2(I2-1))
        ENDIF
!       divide Y2 by cell length
        Y2(I2) = Y2(I2) / CELLEN
      ENDDO
      IF (ITEST.GE.160) THEN
        WRITE (PRTEST, 84) N1, N2
  84    FORMAT (' test CHGBAS ', 2I5)
        WRITE (PRTEST, 85) (X1(II), II = 1, N1)
        WRITE (PRTEST, 85) (Y1(II), II = 1, N1)
        WRITE (PRTEST, 85) (X2(II), II = 1, N2)
        WRITE (PRTEST, 85) (Y2(II), II = 1, N2)
  85    FORMAT (10 (1X,E10.3))
      ENDIF
!
      RETURN
      END

!********************************************************************
!                                                                   *
      REAL(rkind) FUNCTION GAMMAF(XX)
!                                                                   *
!********************************************************************
!
!   Updates
!     ver 30.70, Oct 1997 by N.Booij: new subroutine
!
!   Purpose
!     Compute the transcendental function Gamma
!
!   Subroutines used
!     GAMMLN  (Numerical Recipes)

      USE schism_glbl,ONLY:rkind
!
      REAL(rkind) XX, YY, ABIG                                                   !40.00
      SAVE IENT, ABIG
      DATA IENT /0/, ABIG /30./
      CALL STRACE (IENT, 'GAMMAF')
      YY = REAL( GAMMLN(XX) ) 
      IF (YY.GT.ABIG) YY = ABIG
      IF (YY.LT.-ABIG) YY = -ABIG
      GAMMAF = MyEXP(YY)
      RETURN
      END
!********************************************************************
!                                                                   *
       REAL FUNCTION GAMMLN(XX)
!                                                                   *
!********************************************************************
!
!   Method:
!     function is copied from: Press et al., "Numerical Recipes"
!
      USE schism_glbl, ONLY:rkind,dkind

      !OUBLE PRECISION  COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      REAL(dkind) XX
      REAL(dkind) COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE WRSPEC (NREF, ACLOC)
!                                                                      *
!***********************************************************************

      USE OCPCOMM4                                                        !40.41
      USE SWCOMM3                                                         !40.41
      USE OUTP_DATA                                                       !40.13

      IMPLICIT NONE                                                       !40.13
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
!     40.00, 40.13: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. UPDATE
!
!     new subroutine, update 40.00
!     40.03, Mar. 00: precision increased; 2 decimals more in output table
!     40.13, July 01: variable format using module OUTP_DATA
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. PURPOSE
!
!     Writing of action density spectrum in Swan standard format
!
!  3. METHOD
!
!
!  4. Argument variables
!
!       NREF    int    input    unit ref. number of output file
!       ACLOC   real   local    2-D spectrum or source term at one
!                               output location

      INTEGER, INTENT(IN) :: NREF
      REAL, INTENT(IN)    :: ACLOC(1:MDC,1:MSC)

!  5. Parameter variables
!
!  6. Local variables
!
!       ID      counter of spectral directions
!       IS      counter of spectral frequencies
!
      INTEGER :: ID, IS
!
!       EFAC    multiplication factor written to file
!
      REAL    :: EFAC
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!     SWOUTP (SWAN/OUTP)
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!       ----------------------------------------------------------------
!       determine maximum value of ACLOC
!       if maximum = 0
!       then write 'ZERO' to file
!       else write 'FACTOR'
!            determine multiplication factor, write this to file
!            write values of ACLOC/factor to file
!       ----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0                                           !40.13
      IF (LTRACE) CALL STRACE (IENT, 'WRSPEC')
!
!     first determine maximum energy density
      EFAC = 0.
      DO ID = 1, MDC
        DO IS = 1, MSC
          IF (ACLOC(ID,IS).GE.0.) THEN
            EFAC = MAX (EFAC, ACLOC(ID,IS))
          ELSE
            EFAC = MAX (EFAC, 10.*ABS(ACLOC(ID,IS)))
          ENDIF
        ENDDO
      ENDDO
      IF (EFAC .LE. 1.E-10) THEN
        WRITE (NREF, 12) 'ZERO'                                           !40.00
  12    FORMAT (A4)
      ELSE
        EFAC = 1.01 * EFAC * 10.**(-DEC_SPEC)                             !40.13
!       factor PI/180 introduced to account for change from rad to degr
!       factor 2*PI to account for transition from rad/s to Hz
        WRITE (NREF, 95) EFAC * 2. * PI_W**2 / 180.
  95    FORMAT ('FACTOR', /, E18.8)                                       !40.13
        DO IS = 1, MSC
!         write spectral energy densities to file
          WRITE (NREF, FIX_SPEC) (NINT(ACLOC(ID,IS)/EFAC), ID=1,MDC)      !40.13
        ENDDO
      ENDIF
      RETURN
!     end of subroutine WRSPEC
      END
!****************************************************************
!
!      SUBROUTINE SWACC(AC2, AC2OLD, ACNRMS, ISSTOP, IDCMIN, IDCMAX)
      SUBROUTINE SWACC(      AC2OLD, ACNRMS, ISSTOP, IDCMIN, IDCMAX)
!
!****************************************************************
!
      USE SWCOMM3                                                         !40.41
      USE OCPCOMM4                                                        !40.41

      USE VARS_WAVE,ONLY:AC2
      USE schism_glbl,ONLY:skind,rkind,dkind,npa
!
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
!     40.23: Marcel Zijlema
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     40.23, Sep. 02: New subroutine
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Determine some infinity norms meant for stop criterion
!
!  4. Argument variables
!
!     AC2         action density
!     AC2OLD      action density at previous iteration
!     ACNRMS      array containing infinity norms
!     IDCMIN      integer array containing minimum counter of directions
!     IDCMAX      integer array containing maximum counter of directions
!     ISSTOP      maximum frequency counter in this sweep
!
      INTEGER IDCMIN(MSC), IDCMAX(MSC), ISSTOP
      !REAL    AC2(MDC,MSC,MCGRD), AC2OLD(MDC,MSC), ACNRMS(2)
      REAL(rkind)    AC2OLD(MDC,MSC), ACNRMS(2)
        
!
!  6. Local variables
!
!     DIFFAC:     difference between AC2 and AC2OLD
!     ID    :     counter of direction
!     IDDUM :     uncorrected counter of direction
!     IENT  :     number of entries
!     IS    :     counter of frequency
!
      INTEGER ID, IDDUM, IS, IENT
      REAL(rkind)    DIFFAC
!
!  7. Common blocks used
!
!
!  8. Subroutines used
!     STRACE           Tracing routine for debugging
!
!  9. Subroutines calling
!
!     SWOMPU (in SWANCOM1)
!
! 12. Structure
!
!     determine infinity norms |ac2 - ac2old| and |ac2|
!     in selected sweep
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'SWACC')

      DO IS = 1, ISSTOP
         DO IDDUM = IDCMIN(IS), IDCMAX(IS)
            ID = MOD ( IDDUM - 1 + MDC , MDC ) + 1
!
!           *** determine infinity norms |ac2 - ac2old| and |ac2|
!
            DIFFAC = ABS(AC2(ID,IS,KCGRD(1)) - AC2OLD(ID,IS))
            IF (DIFFAC.GT.ACNRMS(1)) ACNRMS(1) = DIFFAC
            IF (ABS(AC2(ID,IS,KCGRD(1))).GT.ACNRMS(2))&
     &                              ACNRMS(2) = ABS(AC2(ID,IS,KCGRD(1)))

         END DO
      END DO

      RETURN
      END

!***********************************************************************
!                                                                      *
!      SUBROUTINE SWANOUT(DEP2)    !AC2   ,SPCSIG,HSIBC,KGRPNT)                     
      SUBROUTINE SWANOUT
!#     if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
#if 1
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4
      USE SWCOMM1, ONLY : OVEXCV, OVLNAM,  OVSNAM, OUTPAR
      USE SWCOMM3                        
      USE M_WCAP , ONLY : SIGPOW                                
!      USE M_PARALL
!      USE MOD_USGRID  
      USE SWCOMM4, ONLY : TESTFL, NPTST
      USE M_GENARR, ONLY : XYTST,SPCSIG, SPCDIR

      USE VARS_WAVE ! REAL(rkind) HSC1, DIRDEG1,TPEAK,WLEN,QB1

!      USE MOD_STATION_TIMESERIES,ONLY:OUT_STATION_TIMESERIES_ON,OUT_WAVE_PARTITION,NSTA,NODE_STA
!      USE MOD_SPARSE_TIMESERIES,ONLY:OUT_WAVE_SPARSE_TIMESERIES_ON,OUT_WAVE_PARTITION_SPARSE 
!      USE W3PARTMD

      USE schism_glbl,ONLY:skind,rkind,dkind,npa,iplg
      USE schism_msgp, ONLY : myrank,parallel_abort

!  2. Purpose :
!    The following variables are send from SWAN to SCHISM, using the 2D array out_wwm.
!    Based on the same logic and computation method like in WWM. WWM Wave parameters are 
!    given in brackets, followed by SWAN names if they exist
!
!  3. Method : Compute wave parameters and 
!     put them in array out_wwm(npa,35) to use in schism_step.F
!     See details in swanmain to find equivalence between out_wwm index and IVTYPE index


!     Notes in schism <-> swan coupling, HS, TPP, PEAKDM, ORBITAL are mandatory  

!  out_wwm(:,1)  = 'HS'       = SWAN IVTYPE 10    ! Significant wave height           
!  out_wwm(:,2)  = 'TM01'     = SWAN IVTYPE 11    ! Mean average period              
!  out_wwm(:,3)  = 'TM02'     = SWAN IVTYPE 32    ! Zero down crossing period for comparison with buoy.
!  out_wwm(:,4)  = 'TM10'     = SWAN IVTYPE 47    ! Average period of wave runup/overtopping ...
!  out_wwm(:,5)  = 'KLM'                          ! Mean wave number
!  out_wwm(:,6)  = 'WLM'      = SWAN IVTYPE 17    ! Mean wave length
!  out_wwm(:,7)  = 'ETOTS'                        ! Etot energy in x-direction
!  out_wwm(:,8)  = 'ETOTC'                        ! Etot energy in y-direction
!  out_wwm(:,9)  = 'DM'       = SWAN IVTYPE 13    ! Mean average energy transport direction  
!  out_wwm(:,10) = 'DSPR'     = SWAN IVTYPE 16    ! Mean directional spreading     
!  out_wwm(:,11) = 'TPPD'                         ! Discrete peak period (sec)
!  out_wwm(:,12) = 'TPP'                          ! Continuous peak period based on higher order moments (sec)
!  out_wwm(:,13) = 'CPP'                          ! Peak phase vel. (m/s)
!  out_wwm(:,14) = 'WNPP'                         ! Peak n-factor
!  out_wwm(:,15) = 'CGPP'                         ! Peak group vel.  
!  out_wwm(:,16) = 'KPP'                          ! Peak wave number
!  out_wwm(:,17) = 'LPP'                          ! Peak wave length.
!  out_wwm(:,18) = 'PEAKDM'  = SWAN IVTYPE 14     ! Peak (dominant) direction (degr)
!  out_wwm(:,19) = 'PEAKDSPR'                     ! Peak directional spreading
!  out_wwm(:,20) = 'DPEAK'                        ! Discrete peak direction
!  out_wwm(:,21) = 'UBOT'     = SWAN IVTYPE 6     ! Orbital vel. (m/s)
!  out_wwm(:,22) = 'ORBITAL'  = SWAN IVTYPE 34 ?  ! RMS Orbital vel. (m/s)
!  out_wwm(:,23) = 'BOTEXPER'                     ! Bottom excursion (amplitude).
!  out_wwm(:,24) = 'TMBOT'    = SWAN IVTYPE 50    ! Bottom wave period (sec)
!  out_wwm(:,25) = 'URSELL'   = SWAN IVTYPE 45    ! Ursell number based on peak period ...
!  out_wwm(:,26) = 'USTAR'    = SWAN IVTYPE 35    ! Friction velocity          
!!  out_wwm(:,27) = 'Z0'       = SWAN IVTYPE 36    ! Roughness length        WRONG
!!  out_wwm(:,28) = 'ALPHA_CH'                     ! Charnock coefficient   WRONG
!  out_wwm(:,27) = 'ALPHA_CH'                     ! Charnock coefficient
!  out_wwm(:,28) = 'Z0'       = SWAN IVTYPE 36    ! Roughness length       

!  out_wwm(:,29) = !Roller energy dissipation rate (W/m²) @nodes {Drol} 2D    ! Not in SWAN, ToDo
!  out_wwm(:,30) = 'DISSURF' = SWAN IVTYPE 55 ! Total wave energy dissipation rate by depth-induced breaking (W/m²) {wave_sbrtot}
!  out_wwm(:,31) = 'DISBOT'  = SWAN IVTYPE 54 ! Total wave energy dissipation rate by bottom friction (W/m²) {wave_sbftot}
!  out_wwm(:,32) = 'DISWCAP' = SWAN IVTYPE 56 !Total wave energy dissipation rate by whitecapping (W/m²) {wave_sdstot}
!!  iof_wwm(31) = 0 !Total wave energy dissipation rate by vegetation (W/m²) {wave_svegtot} ! Not in SWAN, ToDo
!  out_wwm(:,33) = 0 !Total wave energy input rate from atmospheric forcing (W/m²) {wave_sintot} ! Not in SWAN, ToDo
!  out_wwm(:,34) = 1 !WWM_energy vector {waveEnergyDirX,Y}  2D vector

! Skip Them 

!!  out_wwm(:,29) = 'CD'       = SWAN IVTYPE 38    ! Drag coefficient      WRONG 
!!  out_wwm(:,30) = 'WIND-X'                       ! windx                 WRONG
!!  out_wwm(:,31) = 'WIND-Y'                       ! windy                 WRONG
!!  out_wwm(:,29) = 'WIND-X'                       ! windx               
!!  out_wwm(:,30) = 'WIND-Y'                       ! windy                
!!  out_wwm(:,31) = 'CD'       = SWAN IVTYPE 38    ! Drag coefficient
!!  out_wwm(:,32) = 'CURR-X'
!!  out_wwm(:,33) = 'CURR-Y'
!!  out_wwm(:,34) = 'DEPTH'    = SWAN IVTYPE 4
!!  out_wwm(:,35) = 'ELEVATION'= SWAN IVTYPE 51
!   ------------ New wave parameters ----------------------
!!  out_wwm(:,36) = 'DISSURF' = SWAN IVTYPE 55 ! Surface breaking dissipation 'm2/s' 
!!  out_wwm(:,37) = 'DISWCAP' = SWAN IVTYPE 56 ! Whitecapping dissipation' 'm2/s'
!!  out_wwm(:,38) = 'DISBOT'  = SWAN IVTYPE 54 ! Bottom friction dissipation 'm2/s'
!!  out_wwm(:,39) = 'QB'      = SWAN IVTYPE 8  ! fraction of breaking wave '-'
!!  out_wwm(:,40) = 'wavpres'                  ! Wave induced pressure [Pa]
!!  out_wwm(:,41) = 'ustx'    =                ! Depth averaged Stokes drift, x-comp 'm/s'
!!  out_wwm(:,42) = 'usty'    =                ! Depth averaged Stokes drift, y-comp 'm/s'
!  Other : wave stress determined from the wave dissipation
!!  out_wwm(:,43) = 'fsurx'    = DISSURF+DISWCAP, x comp 'm2/s'
!!  out_wwm(:,44) = 'fsury'    = DISSURF+DISWCAP, y comp 'm2/s'
!!  out_wwm(:,45) = 'fbotx'    = DISBOT, x comp 'm2/s'
!!  out_wwm(:,46) = 'fboty'    = DISBOT, y comp 'm2/s'

!  Notes : The following variables are send from SWAN to SCHISM, they are MANDATORY variables !!
!
!  HS , Significant Wave height (waveheight) [m]
!  TPP, Continuous peak period based on higher order moments (sec)
!  PEAKDM, Peak (dominant) direction (degr)
!  ORBITAL, RMS Orbital vel. (m/s)
!
!  4. Argument variables
!
!     SPCSIG: input  Relative frequencies in computational domain in      
!                    sigma-space                                          
!
!      REAL    SPCSIG(MSC)                                                 
!
!       REALS:
!       ------
!       AC2        action density
!       HSI        significant wave height at boundary (using SWAN
!                  resolution (has thus not to be equal to the WAVEC
!                  significant wave height )
!       ETOT       total energy in a gridpoint
!       DS         increment in frequency space
!       DDIR       increment in directional space
!       HSC        computed wave height after SWAN computation
!       EFTAIL     contribution of tail to spectrum
!
!       INTEGERS:
!       ---------
!       KGRPNT     values of grid indices
!
!  5. SUBROUTINES CALLING
!
!       ---
!
!  6. SUBROUTINES USED
!
!       TRACE
!
!  7. ERROR MESSAGES
!
!       NONE
!
!  8. REMARKS
!
!       NONE
!
!  9. STRUCTURE
!
!     ------------------------------------------------------------------
!     for all computational grid points do                                
!         if HSI is non-zero
!         then compute Hs from action density array
!              if relative difference is large than HSRERR
!              then write error message
!    -------------------------------------------------------------------
!
! 10. SOURCE TEXT
!
!***********************************************************************
!
      IMPLICIT NONE
!      REAL      AC2(MDC,MSC,MCGRD) ,HSIBC(MCGRD)                          
      REAL(rkind)      HSIBC(npa)                         
!
      REAL(rkind)      ETOT, EFTOT ,HSC                                                
!
      INTEGER   ID , IS ,IG, J, IVTYPE
!
      LOGICAL   HSRR
!
!      INTEGER   KGRPNT(MXC,MYC)
!
!      REAL(rkind) :: HSC1(0:MT),HSC1_TEMP(MGL),FTEMP1(MGL),FTEMP2(MGL)
!      REAL(rkind) :: DIRDEG1(0:MT),DIRDEG1_TEMP(MGL)
!      REAL(rkind) :: TPEAK(0:MT),TPEAK_TEMP(MGL)
!-----------------Jianzhong-----------------------------
!      REAL(rkind) :: HSC1_TEMP(MGL),FTEMP1(MGL),FTEMP2(MGL),WLEN_TEMP(MGL)
!      REAL(rkind) :: DIRDEG1_TEMP(MGL)
!      REAL(rkind) :: TPEAK_TEMP(MGL)
!      REAl :: QB1_TEMP(MGL)
!      REAL,ALLOCATABLE :: HSC1_TEMP(:),FTEMP1(:),FTEMP2(:),WLEN_TEMP(:)
!      REAL,ALLOCATABLE :: DIRDEG1_TEMP(:)
!      REAL,ALLOCATABLE :: TPEAK_TEMP(:)
!      REAl,ALLOCATABLE :: QB1_TEMP(:)
!      REAl,ALLOCATABLE :: Pwave_bot_TEMP(:), Ub_swan_TEMP(:) !lwu for output new vars               
!-------------------------------------------------------
 !     REAL(rkind) :: WK(MSC) see mod_wave_current
      REAL(rkind) :: WKLOC(MSC)
      REAL(rkind), PARAMETER  :: eps=1e-14_rkind
      REAL(rkind), PARAMETER  :: ONE=1._rkind
      REAL(rkind), PARAMETER  :: ZERO=0._rkind
      REAL(rkind), PARAMETER :: epsmax=50._rkind
      
      REAL(rkind) :: EEX,EEY,EAD,NN,ND
      REAL(rkind) :: FF1,FM, DIRDEG, DSPR, PEAKFF
      REAL(rkind) :: DS,ETAIL,EFTAIL,EHFR,EMAX,EDI,ETD,HSREL
      REAL(rkind) :: EKTOT,SIG2,SKK,PPTAIL,CETAIL,CKTAIL,OMEG,OMEG2,E1,E2
      REAL(rkind) :: DEPLOC,KWAVELOC,CGLOC,SPCSIGL
      REAL(rkind) :: THETA,UXLOC,UYLOC,UXD
      REAL(rkind) :: FRINTF_X_DDIR,HM,QBLOC,BRCOEF
      !REAL(rkind) :: DS_BAND(MSC+1) see mod_wave_current

! WWM like wave paramters
      REAL(rkind) :: TMP(4)
      REAL(rkind) :: FPP(1),KPP !KPP(1)
      REAL(rkind) :: CGPP,WNPP,CPP,TPP,LPP,PEAKDM,PEAKDSPR,DPEAK
      REAL(rkind) :: MAXAC,ETOTF3,ETOTF4,ETOTC4,ETOTS4
      REAL(rkind) :: KPPD,CGPD,CPPD,TPPD
      INTEGER :: IDIRM,ISIGMP
! WWM like wave paramters
! waves parameters: ORBITAL, ... parameters
      REAL(rkind) :: UB2,AB2
      REAL(rkind) :: UBOT,BOTEXPER,TMBOT,ORBITAL
      REAL(rkind) :: ACTOT_DSIG(MSC),ETOT_DSIG(MSC)
      REAL(rkind) :: SINH_K_X_DEP_2(MSC)
      REAL(rkind) :: ETOT_SIG2_DSIG(MSC)
      REAL(rkind) :: ETOT_SIG2_DSHKD2_DSIG(MSC)
      REAL(rkind) :: ETOT_DSHKD2_DSIG(MSC)

      REAL(rkind) :: HQUOT,HQUOTP
      REAL(rkind) :: ustx,usty,wavpres

      INTEGER :: ISIGM,ITP

      REAL(rkind) :: ACLOC(MDC,MSC)
      REAL(rkind) :: FMIN,FMAX,ECS,ALCQ
      REAL(rkind) :: WLM,KLM
      REAL(rkind) :: ETOT_BK

 ! FOR Breaking, wave pressure head and Stokes velocity
      REAL(rkind) :: DISSURF,DISWCAP,DISBOT,QB
      REAL(rkind) :: HSLOC, CFF, SME_P, SMS_x, SMS_y

 ! FOR DISSIPATION BASED FORCE
      REAL(rkind) :: SME_D,fsurx,fsury,fbotx,fboty
      REAL(rkind) :: WKMAXLOC,SPCMAXLOC ! wave number and freq. with maximum energy
      REAL(rkind) :: DIRMEANLOC(3) ! Mean wave direction + cos and sin function

!----------------------------lwu for output new vars------
!      REAL(rkind)              :: ETOT_DSIG(MSC)
!      REAL(rkind)              :: ACTOT_DSIG(MSC), ETOT_DSHKD2_DSIG(MSC)
!      REAL(rkind)              :: ETOT_DRK_DSIG(MSC), ETOT_SIG2_DSIG(MSC)
!      REAL(rkind)              :: ETOT_K_DSIG(MSC), ETOT_SIG_DSIG(MSC)
!      REAL(rkind)              :: ETOT_SIG2_DSHKD2_DSIG(MSC)
!      REAL(rkind)              :: ETOT_SIG4_DSIG(MSC)
!      REAL(rkind)              :: SINH_K_X_DEP_2(MSC)
!      REAL(rkind)              :: AB2, UB2
!      REAL(rkind)              :: EBSIN(MSC)
!      REAL(rkind)              :: EBCOS(MSC)
!      REAL(rkind)              :: EBSIN_INT(MSC)
!      REAL(rkind)              :: EBCOS_INT(MSC)
!      REAL(rkind)              :: EBCOSTOT, EBSINTOT
!-----------------------------for output wave partition vars----
      REAL(rkind)              :: SPEC(MSC,MDC)
      INTEGER                  :: DIMXP,NSPECC,ISC,KK
!----------------------------------------------------------
      REAL(rkind)              :: DEGCNV,SwanIntgratSpc

      REAL(rkind)              :: eQuot2, eProd2
      REAL(rkind)              :: eMult, eWkReal, eDep, eLoc
      REAL(rkind)              :: eSinh2kd, eSinhkd, eSinhkd2

      INTEGER :: IENT
      DATA IENT/0/, HSRR/.TRUE./

      CALL STRACE (IENT, 'SWANOUT')
!
!     *** initializing ***
!
      HSRR = .TRUE.                                                       

      out_wwm = ZERO
      out_wwm_windpar = ZERO
      SPEC_DENSITY = ZERO
!      DISSURF      = ZERO
!      DISWCAP      = ZERO
!      DISBOT       = ZERO
!      QB           = ZERO
      wavpres      = ZERO
      ustx         = ZERO
      usty         = ZERO
      DISSURF = ZERO
      DISWCAP = ZERO
      DISBOT  = ZERO
      QB = ZERO
      WKLOC = ZERO
      HSLOC = ZERO
      ACLOC = ZERO
      WKMAXLOC = ZERO
      SPCMAXLOC = ZERO
      FMIN = ZERO
      FMAX = ZERO
      ECS  = ZERO
      ALCQ = ZERO ! JL should be 0 ...

! JL : Change applied using swan 41.10, see swanout1.ftn
!
!     coefficient for high frequency tail
!
      EFTAIL = ONE / (PWTAIL(1) - ONE)

      FRINTF_X_DDIR = FRINTF * DDIR

      !DS_BAND(0:MSC+1) = ZERO
      !DS_BAND(0)     = SPCSIG(2)- SPCSIG(1)
      !DS_BAND(1)     = DS_BAND(0)
      !DS_BAND(MSC)   = SPCSIG(MSC) - SPCSIG(MSC-1)
      !DS_BAND(MSC+1) = DS_BAND(MSC)
      !DO IS = 2, MSC-1 ! Bandwith at gridpoints
      !    DS_BAND(IS) = (SPCSIG(IS)-SPCSIG(IS-1))/2._rkind + (SPCSIG(IS+1)-SPCSIG(IS))/2._rkind
      !END DO

      DO IG = 1, npa

      !  print*,IG,npa

        TESTFL = .false.
        if ( NPTST > 0 ) then
            do J = 1, NPTST
               IF( IG /= XYTST(J) ) cycle
               TESTFL = .true.
            enddo
        endif

         ACLOC(1:MDC,1:MSC) = AC2(1:MDC,1:MSC,IG)
         DEPLOC = COMPDA(IG,JDP2)
         UXLOC  = COMPDA(IG,JVX2)  ! Ambiant Velocity
         UYLOC  = COMPDA(IG,JVY2)

         ! Compute Kwave ! Already done in SWAN/mod_wave_current.F90
!        IF(DEPLOC <= DEPMIN)THEN
!          *** depth is negative ***
!            KWAVELOC = -1._rkind
!            WK(1:MSC) = KWAVELOC
!        ELSE
!        ! *** call KSCIP1 to compute KWAVE and CGO ***
!          DO IS=1, MSC
!             SPCSIGL = SPCSIG(IS)
!             CALL KSCIP1(1,SPCSIGL,DEPLOC,KWAVELOC,CGLOC,NN,ND)
!             !KWAVELOC: wave number K, CGLOC: group velocity, 
!             !NN:ratio of group and phase velocity and its derivative ND
!             WK(IS) = KWAVELOC
!          ENDDO             
!        ENDIF

         WKLOC(1:MSC) = WK(1:MSC,IG)
!
!        ----------------- HSC1 ---------------------
         IVTYPE = 10 ! significant wave height Hs (m)

!        significant wave height (Hs)
         !WRITE(PRINTF,*)' significant wave height (Hs)'
          IF (OUTPAR(6).EQ.0.) THEN                                       !40.87
!            integration over [0,inf]                                     !40.87
             !WRITE(PRINTF,*)' integration over [0,inf]'
!
             ETOT = ZERO
!            trapezoidal rule is applied
             DO ID=1, MDC
               DO IS=2,MSC
                 DS=SPCSIG(IS)-SPCSIG(IS-1)                               !30.72
                 EAD = 0.5_rkind*(SPCSIG(IS)*ACLOC(ID,IS)+  &             !30.72
     &                 SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS*DDIR               !30.72
                 ETOT = ETOT + EAD
               ENDDO
               IF (MSC .GT. 3) THEN                                       !10.20
!                contribution of tail to total energy density
                 EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                       !30.72
                 ETOT = ETOT + DDIR * EHFR * SPCSIG(MSC) * EFTAIL         !30.72
               ENDIF
             ENDDO
          ELSE                                                            !40.87
!            integration over [fmin,fmax]                                 !40.87
             !WRITE(PRINTF,*)' integration over [fmin,fmax]'
!
             FMIN = PI2_W*OUTPAR(21)                                      !40.87
             FMAX = PI2_W*OUTPAR(36)                                      !40.87
             ECS  = ONE                                                   !40.87
             ETOT = SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC , ECS , 0.  , 0.    , ACLOC   ,& !40.87
     &                             1  )                                   !40.87
          ENDIF 

          IF (ETOT .GE. ZERO) THEN                                        !30.00
             HSLOC = 4._rkind*SQRT(ETOT)
          ELSE
             HSLOC = ZERO
          ENDIF

          out_wwm(IG,1) = HSLOC

          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG,iplg(IG),OVSNAM(IVTYPE),out_wwm(IG,1) !40.00
          ENDIF

         ETOT_BK = ETOT

!        ----------------- TM01 ---------------------
         IVTYPE = 11 ! ! Mean average period (s)

           IF (OUTPAR(7).EQ.0) THEN                                      !40.87
!             integration over [0,inf]                                   !40.87
              ETOT = ZERO
              EAD  = ZERO
              EFTOT = ZERO
              PPTAIL = PWTAIL(1) - ONE                                   !40.00
              ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))     !20.61
              PPTAIL = PWTAIL(1) - 2._rkind                              !40.00
              EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))     !20.61
              DO ID=1, MDC
                 THETA = SPCDIR(ID,1) + ALCQ                             !20.43
                 UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 DO IS = 1, MSC
                   OMEG = SPCSIG(IS) + WKLOC(IS) * UXD                   !30.72
                   EAD = FRINTF * SPCSIG(IS)**2 * ACLOC(ID,IS)           !40.00
                   ETOT = ETOT + EAD
                   EFTOT = EFTOT + EAD * OMEG                            !20.66
                 ENDDO
                 IF (MSC .GT. 3) THEN                                    !10.20
!                  contribution of tail to total energy density
                   EAD   = SPCSIG(MSC)**2 * ACLOC(ID,MSC)                !40.00
                   ETOT  = ETOT + ETAIL * EAD
                   EFTOT = EFTOT + EFTAIL * OMEG * EAD
                 ENDIF
              ENDDO
           ELSE                                                          !40.87
!             integration over [fmin,fmax]                               !40.87
              FMIN = PI2_W*OUTPAR(22)                                    !40.87
              FMAX = PI2_W*OUTPAR(37)                                    !40.87
              ECS  = ONE                                                 !40.87
              ETOT =SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC , ECS , UXLOC, UYLOC, ACLOC   ,& !40.87
     &                             2  )                                  !40.87
              EFTOT=SwanIntgratSpc(1. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC , ECS , UXLOC, UYLOC, ACLOC   ,& !40.87
     &                             2  )                                  !40.87
           ENDIF                                                         !40.87
           IF (EFTOT.GT.ZERO) THEN
              out_wwm(IG,2) = 2._rkind*PI_W * ETOT / EFTOT
           ELSE
              out_wwm(IG,2) = ZERO
           ENDIF

          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG,iplg(IG),OVSNAM(IVTYPE),out_wwm(IG,2) !40.00
          ENDIF

!        ----------------- TM02 ---------------------
         IVTYPE = 32    ! Zero down crossing period for comparison with buoy.

           IF (OUTPAR(15).EQ.0) THEN                                      !40.87
!             integration over [0,inf]                                    !40.87
              ETOT  = ZERO
              EFTOT = ZERO
              EAD   = ZERO
              PPTAIL = PWTAIL(1) - ONE                                     !20.61
              ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))       !20.61
              PPTAIL = PWTAIL(1) - 3._rkind                                !20.61
              EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))       !20.61
              DO ID=1, MDC
                 IF (ICUR.GT.0) THEN
                   THETA = SPCDIR(ID,1) + ALCQ
                   UXD   = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 ENDIF
                 DO IS=1,MSC
                   EAD  = SPCSIG(IS)**2 * ACLOC(ID,IS) * FRINTF          !30.72
                   IF (ICUR.GT.0) THEN
                     OMEG  = SPCSIG(IS) + WKLOC(IS) * UXD                !30.72
                     OMEG2 = OMEG**2
                   ELSE
                     OMEG2 = SPCSIG(IS)**2                                !30.72
                   ENDIF
                   ETOT  = ETOT + EAD                                    !20.61
                   EFTOT = EFTOT + EAD * OMEG2                           !20.61
                 ENDDO
                 IF (MSC .GT. 3) THEN
!                  contribution of tail to total energy density
                   EAD  = SPCSIG(MSC)**2 * ACLOC(ID,MSC)                 !30.72
                   ETOT  = ETOT  + ETAIL * EAD                           !20.61
                   EFTOT = EFTOT + EFTAIL * OMEG2 * EAD                  !20.61
                 ENDIF
              ENDDO
           ELSE                                                           !40.87
!             integration over [fmin,fmax]                                !40.87
              FMIN = PI2_W*OUTPAR(30)                                       !40.87
              FMAX = PI2_W*OUTPAR(45)                                       !40.87
              ECS  = ONE                                                   !40.87
              IF (ICUR.GT.0) THEN                                         !40.87
                 ITP = 2                                                  !40.87
              ELSE                                                        !40.87
                 ITP = 1                                                  !40.87
              ENDIF                                                       !40.87
              ETOT  = SwanIntgratSpc(0.          , FMIN , FMAX, SPCSIG,&  !40.87
     &                               SPCDIR(1,1) , WKLOC, ECS , UXLOC ,&  !40.87
     &                               UYLOC       , ACLOC, ITP )           !40.87
              EFTOT = SwanIntgratSpc(2.          , FMIN , FMAX, SPCSIG,&  !40.87
     &                               SPCDIR(1,1) , WKLOC, ECS , UXLOC ,&  !40.87
     &                               UYLOC       , ACLOC, ITP )           !40.87
           ENDIF                                                          !40.87
           IF (EFTOT.GT.ZERO) THEN
              out_wwm(IG,3) = 2._rkind*PI_W * SQRT(ETOT/EFTOT)             !20.61
           ELSE
              out_wwm(IG,3) = ZERO
           ENDIF

          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG,iplg(IG),OVSNAM(IVTYPE),out_wwm(IG,3) !40.00
          ENDIF

!        ----------------- TM10 ---------------------
          IVTYPE = 47    ! Average period of wave runup/overtopping ...

           IF (OUTPAR(19).EQ.0.) THEN                                     !40.87
!             integration over [0,inf]                                    !40.87
              ETOT = ZERO
              EFTOT = ZERO
              OMEG2 = ZERO
              EAD = ZERO
              PPTAIL = PWTAIL(1)
              ETAIL  = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))
              PPTAIL = PWTAIL(1) - ONE
              EFTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-1.)))
              DO ID=1, MDC
                 THETA = SPCDIR(ID,1) + ALCQ
                 UXD = UXLOC*COS(THETA) + UYLOC*SIN(THETA)
                 DO IS = 1, MSC
                   OMEG = SPCSIG(IS) + WKLOC(IS) * UXD
                   OMEG2 = OMEG ** (-ONE)
                   EAD = OMEG2 * FRINTF * SPCSIG(IS)**2 * ACLOC(ID,IS)
                   ETOT = ETOT + EAD
                   EFTOT = EFTOT + EAD * OMEG
                 ENDDO
                 IF (MSC .GT. 3) THEN
!                  contribution of tail to total energy density
                   EAD = OMEG2 * SPCSIG(MSC)**2 * ACLOC(ID,MSC)
                   ETOT = ETOT + ETAIL * EAD
                   EFTOT = EFTOT + EFTAIL * OMEG * EAD
                 ENDIF
              ENDDO
           ELSE                                                           !40.87
!             integration over [fmin,fmax]                                !40.87
              FMIN = PI2_W*OUTPAR(34)                                     !40.87
              FMAX = PI2_W*OUTPAR(49)                                     !40.87
              ECS  = ONE                                                  !40.87
              ETOT =SwanIntgratSpc(-1., FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC, ECS , UXLOC, UYLOC, ACLOC    ,& !40.87
     &                             2  )                                   !40.87
              EFTOT=SwanIntgratSpc( 0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC, ECS , UXLOC, UYLOC, ACLOC    ,& !40.87
     &                             2  )                                   !40.87
           ENDIF                                                          !40.87
           IF (EFTOT.GT.ZERO) THEN
              out_wwm(IG,4) = 2._rkind*PI_W * ETOT / EFTOT
           ELSE
              out_wwm(IG,4) = ZERO
           ENDIF

!        -----------Mean wave length and Mean wave number------
        IVTYPE = 17 ! average wave length,(Wlen, m)

           IF (OUTPAR(11).EQ.0) THEN                                      !40.87
!             integration over [0,inf]                                    !40.87
              ETOT  = ZERO
              EKTOT = ZERO
!             new integration method involving FRINTF                     !20.59
              DO IS=1, MSC
                 SIG2 = (SPCSIG(IS))**2                                   !30.72
                 SKK  = SIG2 * (WKLOC(IS))**OUTPAR(3)                     !40.00
                 DO ID=1,MDC
                   ETOT  = ETOT + SIG2 * ACLOC(ID,IS)                     !20.59
                   EKTOT = EKTOT + SKK * ACLOC(ID,IS)                     !20.59
                 ENDDO
              ENDDO
              ETOT  = FRINTF * ETOT
              EKTOT = FRINTF * EKTOT

              IF (MSC .GT. 3) THEN                                        !10.20
!                contribution of tail to total energy density
                 PPTAIL = PWTAIL(1) - ONE                                 !40.00
                 CETAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))  !20.61
                 PPTAIL = PWTAIL(1) - 2._rkind                            !40.00
                 IF (PPTAIL.LE.ZERO) THEN
                   CALL MSGERR (2,'WLEN: error tail computation')
                   GOTO 480
                 ENDIF
                 CKTAIL = ONE / (PPTAIL * (ONE + PPTAIL * (FRINTH-ONE)))  !20.59
                 DO ID=1,MDC
                   ETOT   = ETOT + CETAIL * SIG2 * ACLOC(ID,MSC)          !20.59
                   EKTOT  = EKTOT + CKTAIL * SKK * ACLOC(ID,MSC)          !20.59
                 ENDDO
 480             CONTINUE
              ENDIF
              IF (EKTOT.GT.ZERO) THEN
                 WLM = PI2_W * (ETOT / EKTOT) ** (ONE/OUTPAR(3))          !40.00
                 KLM = PI2_W/WLM
              ELSE
                 WLM = ZERO
                 KLM = 10.0_rkind
              ENDIF
           ELSE                                                           !40.87
!             integration over [fmin,fmax]                                !40.87
              FMIN = PI2_W*OUTPAR(26)                                     !40.87
              FMAX = PI2_W*OUTPAR(41)                                     !40.87
              ECS  = ONE                                                  !40.87
              ETOT  = SwanIntgratSpc(OUTPAR(3)-1., FMIN , FMAX, SPCSIG,&  !40.87
     &                               SPCDIR(1,1) , WKLOC, ECS , 0.    ,&  !40.87
     &                               0.          , ACLOC, 3   )           !40.87
              EKTOT = SwanIntgratSpc(OUTPAR(3)   , FMIN , FMAX, SPCSIG,&  !40.87
     &                               SPCDIR(1,1) , WKLOC, ECS , 0.    ,&  !40.87
     &                               0.          , ACLOC, 3   )           !40.87
              IF (EKTOT.GT.ZERO) THEN                                       !40.87
                 WLM = PI2_W * ETOT / EKTOT                               !40.87
                 KLM = PI2_W/WLM
              ELSE                                                        !40.87
                 KLM = 10.0_rkind
                 WLM = ZERO                                               !40.87
              ENDIF                                                       !40.87
           ENDIF                                                          !40.87

           out_wwm(IG,5)  = KLM
           out_wwm(IG,6)  = WLM

          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), WLM         !40.00
          ENDIF

!  -----------------------ETOTS=EEY  ! Etot energy in y-direction
!  -----------------------ETOTC=EEX  ! Etot energy in x-direction
!  -----------------------DM=DIRDEG  ! Mean average energy transport direction
!  -----------------------DSPR   ! Mean directional spreading
!           IVTYPE = 13 ! DM, Average wave direction, Dir (degree)
!           IVTYPE = 16 ! DSPR, Mean directional spreading

           IF (OUTPAR(8).EQ.0.) THEN                                     !40.87
!             integration over [0,inf]                                   !40.87
              ETOT = 0.
              EEX  = 0.  ! = ETOTC
              EEY  = 0.  ! = ETOTS

              DO ID=1, MDC
                 EAD = 0.
                 DO IS=2,MSC
                   DS=SPCSIG(IS)-SPCSIG(IS-1)                             !30.72
                   EDI = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+       &            !30.72
     &                        SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS             !30.72
                   EAD = EAD + EDI
                 ENDDO
                 IF (MSC .GT. 3) THEN                                     !10.20
!                  contribution of tail to total energy density
                   EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                     !30.72
                   EAD = EAD + EHFR * SPCSIG(MSC) * EFTAIL                !30.72
                 ENDIF
                 EAD = EAD * DDIR
                 ETOT = ETOT + EAD
                 EEX  = EEX + EAD * SPCDIR(ID,2) ! = ETOTC
                 EEY  = EEY + EAD * SPCDIR(ID,3) ! = ETOTS
              ENDDO
           ELSE                                                           !40.87
!             integration over [fmin,fmax]                                !40.87
              FMIN = PI2_W*OUTPAR(23)                                     !40.87
              FMAX = PI2_W*OUTPAR(38)                                     !40.87
              ECS  = ONE                                                  !40.87
              ETOT= SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),& !40.87
     &                             WKLOC, ECS , 0.  , 0.    , ACLOC    ,& !40.87
     &                             1  )                                   !40.87
              EEX = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),&  !40.87
     &                             WKLOC, SPCDIR(1,2), 0., 0., ACLOC, 1)  !40.87
              EEY = SwanIntgratSpc(0., FMIN, FMAX, SPCSIG, SPCDIR(1,1),&  !40.87
     &                             WKLOC, SPCDIR(1,3), 0., 0., ACLOC, 1)  !40.87
           ENDIF                                                          

           IF (ETOT.GT.ZERO) THEN
              IF (BNAUT) THEN                                             !30.70
                 DIRDEG = ATAN2(EEY,EEX) * 180._rkind/PI_W                !10.15
              ELSE
                 DIRDEG = (ALCQ + ATAN2(EEY,EEX)) * 180._rkind/PI_W       !10.15
              ENDIF
              IF (DIRDEG.LT.ZERO) DIRDEG = DIRDEG + 360._rkind              !10.15
!
!             *** Convert (if necessary) from nautical degrees ***        32.01
!             *** to cartesian degrees                         ***        32.01
!
              DIRDEG = DEGCNV( DIRDEG )                                   !32.01
              FF1 = MIN (ONE, SQRT(MAX(ZERO,EEX*EEX+EEY*EEY))/ETOT)
              DSPR = SQRT(MAX(ZERO,2._rkind-2._rkind*FF1)) * 180._rkind/PI_W
           ELSE
              FM     = ZERO
              DIRDEG = ZERO
              DSPR   = ZERO
           ENDIF

           out_wwm(IG,7) = EEY
           out_wwm(IG,8) = EEX
           out_wwm(IG,9) = DIRDEG
           out_wwm(IG,10) = DSPR

          !Keep a trace for wave dissipation forces
           DIRMEANLOC(1) = DIRDEG*PI_W/180._rkind
           DIRMEANLOC(2) = cos(DIRMEANLOC(1))
           DIRMEANLOC(3) = sin(DIRMEANLOC(1))

          IVTYPE = 13 !Average wave direction, Dir (degree)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                    !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), DIRDEG !40.00
          ENDIF

          IVTYPE = 16 ! Mean directional spreading
          IF (TESTFL.AND.ITEST.GE.60) THEN                                    !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), DSPR   !40.00
          ENDIF

!  -----------------------TPPD  Discrete peak period (sec)
!  -----------------------TPP   = SWAN IVTYPE 12 ! Continuous peak period based on higher order moments 
!  -----------------------CPP   Peak phase vel. (m/s)
!  -----------------------WNPP  Peak n-factor
!  -----------------------CGPP  Peak group vel.
!  -----------------------KPP   Peak wave number
!  -----------------------LPP   Peak wave length.
!  -----------------------PEAKDM'  = SWAN IVTYPE 14     ! Peak (dominant) direction (degr)
!  -----------------------PEAKDSPR         ! Peak directional spreading
!  -----------------------DPEAK            ! Discrete peak direction

! Peak period, not the SWAN version but WWM version !
! (wwm_specparam.F90, SUBROUTINE PEAK_PARAMETER), see details in
! Thesis Jose-Henrique Alves ... correct citation is given there ... :)

       MAXAC = MAXVAL(ACLOC)
       IF (MAXAC .gt. ZERO .AND. DEPLOC .GT. DEPMIN) THEN
         ETOTF3 = ZERO
         ETOTF4 = ZERO
         ETOTC4 = ZERO
         ETOTS4 = ZERO
         HQUOT  = ZERO
         HQUOTP = ZERO
         DO IS = 1, MSC
           DO ID = 1, MDC
             HQUOT  = ACLOC(ID,IS)/MAXAC
             HQUOTP = HQUOT**4
             !DS = SPCSIG(IS)-SPCSIG(IS-1)                             ! 30.72
             DS = DS_BAND(IS)
             ETOTF3 = ETOTF3 + SPCSIG(IS) * HQUOTP * DS
             ETOTF4 = ETOTF4 +              HQUOTP * DS
             ETOTC4 = ETOTC4 + SPCDIR(ID,2) * HQUOTP * DS
             ETOTS4 = ETOTS4 + SPCDIR(ID,3) * HQUOTP * DS
           END DO
         END DO

         IF(ETOTF3 .GT. ZERO .AND. ETOTF4 .GT. ZERO) THEN
           FPP(1)    = ETOTF3/ETOTF4  !*PI2_W
           SPCMAXLOC = FPP(1) ! Keep a trace: Frequency with max energy
!           CALL KSCIP1(1,FPP(1),DEPLOC,KPP,CGPP,WNPP,ND)  ! KSCIP1(MMT,SIG,D,K,CG,N,ND)
           CALL KSCIP1(1,FPP(1),DEPLOC,TMP(1),TMP(2),TMP(3),TMP(4))  
           KPP = TMP(1)
           CGPP= TMP(2)
           WNPP= TMP(3)
           ND  = TMP(4)
           ! comparable to WWM function CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           ! KPP   = KWAVELOC ; WNPP  = NN; CGPP  = CGLOC; CPP   = CGLOC/NN
           WKMAXLOC = KPP ! Keep a trace: wave number with max energy
           TPP    = PI2_W/FPP(1)
           !LPP    = PI2_W/KPP(1)
           LPP    = PI2_W/KPP
           CPP    = CGPP/WNPP

!          call KSCIP1 (MSC, spcsig, deploc, kwave(1,ic), cgo(1,ic), n, nd)

           ! comparable to WWM function CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           ! KPP   = KWAVELOC ; WNPP  = NN; CGPP  = CGLOC; CPP   = CGLOC/NN
           !WKMAXLOC = KPP(1) ! Keep a trace: wave number with max energy
           !TPP    = PI2_W/FPP(1)
           !LPP    = PI2_W/KPP(1)
           !LPP    = PI2_W/KPP
           !CPP    = CGPP/WNPP

           IF (BNAUT) THEN                                             !30.70
              PEAKDM = ATAN2(ETOTS4,ETOTC4) * 180._rkind/PI_W                
           ELSE
              PEAKDM = (ALCQ + ATAN2(ETOTS4,ETOTC4)) * 180._rkind/PI_W       
           ENDIF
           IF (PEAKDM.LT.ZERO) PEAKDM = PEAKDM + 360._rkind             
!          *** Convert (if necessary) from nautical degrees ***        32.01
!          *** to cartesian degrees                         ***        32.01
           PEAKDM = DEGCNV( PEAKDM )                                   !32.01
           PEAKFF = MIN(ONE,SQRT(MAX(ZERO,ETOTC4*ETOTC4+ETOTS4*ETOTS4))/ETOTF4)
           PEAKDSPR = SQRT(MAX(ZERO,2._rkind-2._rkind*PEAKFF)) * 180._rkind/PI_W
         ELSE
           FPP = ZERO
           SPCMAXLOC = ZERO
           KPP = 10.0_rkind 
           WKMAXLOC = 10.0_rkind
           CGPP= ZERO   
           WNPP= ZERO   
           CPP = ZERO     
           TPP = ZERO
           LPP = ZERO
           PEAKDM = ZERO
           PEAKDSPR = ZERO
         END IF

         DPEAK = 1
         ETOT = ZERO
         IDIRM = -1
         DO ID = 1, MDC
            EAD = ZERO
            DO IS = 2, MSC
               DS = SPCSIG(IS)-SPCSIG(IS-1)
               E1 = SPCSIG(IS-1)*ACLOC(ID,IS-1)
               E2 = SPCSIG(IS)*ACLOC(ID,IS)
               EAD = EAD + DS*(E1+E2)
            END DO
            IF (EAD .GT. ETOT) THEN
               ETOT = EAD
               IDIRM = ID
            END IF
         END DO
         IF (IDIRM.GT.0) THEN
           DPEAK    = SPCDIR(IDIRM,1) * 180._rkind/PI_W
           DPEAK    = DEGCNV( DPEAK )
         ELSE
           DPEAK = ZERO
         END IF

       ELSE
         FPP = ZERO
         KPP = 10.0_rkind
         CGPP = ZERO
         WNPP = ZERO
         CPP = ZERO
         TPP = ZERO
         LPP = ZERO
         PEAKDM = ZERO
         PEAKDSPR = ZERO
         DPEAK = ZERO
       END IF

       ! TPPD ...Discrete Peak Period

       KPPD = ZERO
       CGPD = ZERO
       CPPD = ZERO
       TPPD = ZERO
       ETOT = ZERO
       ISIGMP = -1
       DO IS = 1, MSC
         EAD = ZERO
         DO ID = 1, MDC
            EAD = EAD + SPCSIG(IS)*ACLOC(ID,IS)*DDIR
         ENDDO
         IF (EAD > ETOT) THEN
           ETOT = EAD
           ISIGMP = IS
         END IF
       END DO
       IF (ISIGMP.GT.0) THEN
          TPPD = ONE/(SPCSIG(ISIGMP)/PI2_W)
           !CALL KSCIP1(1,SPCSIG(ISIGMP),DEPLOC,KPPD,CGPD,NN,ND)
           !CPPD = CGPD/NN
       ELSE
          TPPD = ZERO
           !CPPD  = ZERO
           !KPPD  = ZERO
           !CGPD  = ZERO
       END IF

       out_wwm(IG,11) = TPPD
       out_wwm(IG,12) = TPP     != SWAN IVTYPER 12 !Peak Period (s)
       out_wwm(IG,13) = CPP
       out_wwm(IG,14) = WNPP
       out_wwm(IG,15) = CGPP
       out_wwm(IG,16) = KPP !KPP(1)
       out_wwm(IG,17) = LPP
       out_wwm(IG,18) = PEAKDM  != SWAN IVTYPE 14  ! Peak (dominant) direction (degr)
       out_wwm(IG,19) = PEAKDSPR
       out_wwm(IG,20) = DPEAK

       IVTYPE = 12 !Peak Period (s)
       IF (TESTFL.AND.ITEST.GE.60) THEN                                    !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), TPP !40.00
       ENDIF

       IVTYPE = 14 ! Mean directional spreading
       IF (TESTFL.AND.ITEST.GE.60) THEN                                    !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), PEAKDM !40.00
       ENDIF

! JL to compare with
#if 0
!         --------------------- TPEAK (SWAN version---------------------
          IVTYPE = 12  ! Peak Period (S)

           EMAX = 0.
           ISIGM = -1
           DO IS = 1, MSC
              ETD = 0.
              DO ID = 1, MDC
                ETD = ETD + SPCSIG(IS)*ACLOC(ID,IS)*DDIR          ! 30.72
              ENDDO
              IF (ETD.GT.EMAX) THEN
                EMAX  = ETD
                ISIGM = IS
              ENDIF
           ENDDO
           IF (ISIGM.GT.0) THEN
             TPP = PI2_W/SPCSIG(ISIGM)                    ! 30.72
           ELSE
             TPP = ZERO
           ENDIF

           IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), TPP         !40.00
           ENDIF
#endif
! JL to compare with
#if 0
!         --------------------- PEAKDM (SWAN version)---------------------
          IVTYPE = 14 ! direction of the peak of the spectrum, PkDir (degree)

           EMAX = ZERO
           IDIRM = -1
           DO ID = 1, MDC
              EAD = ZERO
              DO IS = 2, MSC
                DS = SPCSIG(IS)-SPCSIG(IS-1)                             ! 30.72
                E1 = SPCSIG(IS-1)*ACLOC(ID,IS-1)                         ! 30.72
                E2 = SPCSIG(IS)*ACLOC(ID,IS)                             ! 30.72
                EAD = EAD + DS * (E1+E2)
              ENDDO
              IF (EAD.GT.EMAX) THEN
                EMAX  = EAD
                IDIRM = ID
              ENDIF
           ENDDO
           IF (IDIRM.GT.0) THEN
!
!            *** Convert (if necessary) from nautical degrees ***         32.01
!            *** to cartesian degrees                         ***         32.01
!
             IF (BNAUT) THEN                                              !30.70
               PEAKDM =  SPCDIR(IDIRM,1) *180._rkind/PI_W
             ELSE
               PEAKDM = (ALCQ + SPCDIR(IDIRM,1)) *180._rkind/PI_W
             ENDIF
             PEAKDM = DEGCNV( PEAKDM )        !32.01
!
           ELSE
             PEAKDM = OVEXCV(IVTYPE)
           ENDIF

           IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), PEAKDM         !40.00
           ENDIF
#endif

!  out_wwm(:,21) = 'UBOT'     = SWAN IVTYPE 6     ! Orbital vel. (m/s)
!  out_wwm(:,22) = 'ORBITAL'  = SWAN IVTYPE 34    ! RMS Orbital vel. (m/s)
!  out_wwm(:,23) = 'BOTEXPER'                     ! Bottom excursion amplitude.
!  out_wwm(:,24) = 'TMBOT'    = SWAN IVTYPE 50    ! Bottom wave period (sec)
!  SWAN Method : see SINTGRL
         UBOT    = ZERO
         BOTEXPER= 0.001_rkind
         TMBOT   = ZERO
         ORBITAL = ZERO

         IF(ETOT_BK.GT.ZERO)  THEN

           ETOT_DSIG             = ZERO
           ACTOT_DSIG            = ZERO
           SINH_K_X_DEP_2        = ZERO
           ETOT_SIG2_DSIG        = ZERO
           ETOT_SIG2_DSHKD2_DSIG = ZERO

           SINH_K_X_DEP_2(:)   = SINH(MIN(30._rkind,WKLOC(:)*DEPLOC))**2_rkind

           ACTOT_DSIG(:)       = SUM(ACLOC(:,:),DIM=1)*SIGPOW(:,1)*FRINTF_X_DDIR
           ETOT_DSIG(:)        = ACTOT_DSIG(:) * SIGPOW(:,1)
           ETOT_SIG2_DSIG(:)   = ACTOT_DSIG(:) / SIGPOW(:,3)                 
           ETOT_DSHKD2_DSIG(:) = ETOT_DSIG(:) / SINH_K_X_DEP_2(:)
           ETOT_SIG2_DSHKD2_DSIG(:) = ETOT_DSHKD2_DSIG(:)*SIGPOW(:,2)

           UB2           = SUM(ETOT_SIG2_DSHKD2_DSIG)
           AB2           = SUM(ETOT_DSHKD2_DSIG)
         ELSE
           UB2           = ZERO
           AB2           = ZERO
         ENDIF
      
!     Calculate the orbital velocity UBOT, orbital excursion ABRBOT and
!     near bottom wave period TMBOT                                       40.51
!
        IF ( UB2.GT.ZERO) UBOT     = SQRT ( UB2 )
        IF ( UB2.GT.ZERO) ORBITAL  = SQRT ( 2._rkind*UB2 )
        IF ( AB2.GT.ZERO) BOTEXPER = SQRT ( 2._rkind*AB2 )
        IF ( UB2.GT.ZERO .AND. AB2.GT.ZERO) TMBOT = PI2_W*SQRT(AB2/UB2)   !40.51
!
        out_wwm(IG,21) = UBOT
        out_wwm(IG,22) = ORBITAL
        out_wwm(IG,23) = BOTEXPER
        out_wwm(IG,24) = TMBOT

        IVTYPE = 6
        IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,21)
        ENDIF

        IVTYPE = 34
        IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,22)
        ENDIF

        IVTYPE = 50
        IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,24)
        ENDIF

! JL to compare with        
#if 0
!         ---------------  Ub_swan -------------------
          IVTYPE = 6 ! Orbital velocity at the bottom, Ubot)
! JL : We use COMPDA(:,JUBOT) [with JUBOT= 3] bottom orbital velocity within array COMPDA
          !Ub_swan(IG) = COMPDA(IG,JUBOT)
          out_wwm(IG,21) =  COMPDA(IG,JUBOT)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,21)
          ENDIF

          IVTYPE = 34 ! RMS Orbital vel. (m/s)
          out_wwm(IG,22) =  SQRT(2*COMPDA(IG,JUBOT)**2)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,22)
          ENDIF

          IVTYPE = 50 ! ! Bottom wave period (sec)
          out_wwm(IG,24) =  COMPDA(IG,JPBOT)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,24)
          ENDIF
#endif

!  out_wwm(:,25) = 'URSELL'   = SWAN IVTYPE 45    ! Ursell number based on peak period ...
!  out_wwm(:,26) = 'USTAR'    = SWAN IVTYPE 35    ! Friction velocity
!  out_wwm(:,27) = 'ALPHA_CH'                     ! Charnock coefficient
!  out_wwm(:,28) = 'Z0'       = SWAN IVTYPE 36    ! Roughness length

          IVTYPE = 45
          out_wwm(IG,25) =  COMPDA(IG,JURSEL)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,25)
          ENDIF

! JLefevre : USTAR / Zeff / CHARNOCK following the Janssen et al parameterization
!            To use later as input for Ocean-Atm coupling
          IF(IWIND.EQ.4 .AND. IWCAP.EQ.2) THEN
            out_wwm(IG,26) = COMPDA(IG,JUSTAR)
            out_wwm(IG,27) = COMPDA(IG,JCHARN)
            out_wwm(IG,28) = COMPDA(IG,JZEL)
          ENDIF

          IVTYPE = 36 ! 'ZLEN'
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         
            WRITE(PRINTF, 222) IG, iplg(IG), 'Ustar',  out_wwm(IG,26)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Charnock',  out_wwm(IG,27)
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,28)
          ENDIF

!  J.L. Replace Surface Z0 by bottom Z0, the apparent bottom roughness (kn)
          out_wwm(IG,28) = COMPDA(IG,JFRC2)
 
!  out_wwm(:,29) = 'CD'       = SWAN IVTYPE 38    ! Drag coefficient
!  out_wwm(:,30) = 'WIND-X'                       ! windx
!  out_wwm(:,31) = 'WIND-Y'                       ! windy

          IVTYPE = 38        
!          out_wwm(IG,29) =  COMPDA(IG,JCDRAG)
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), COMPDA(IG,JCDRAG)
          ENDIF

          out_wwm(IG,29) = COMPDA(IG,JCDRAG)
          out_wwm(IG,30) = COMPDA(IG,JWX2)
          out_wwm(IG,31) = COMPDA(IG,JWY2)

          IVTYPE = 38  ! CD, wind+wave surface Drag coefficient
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,29)
          ENDIF

!  out_wwm(:,32) = 'CURR-X'    = COMPDA(IG,JVX2)
!  out_wwm(:,33) = 'CURR-Y'    = COMPDA(IG,JVY2)
!  out_wwm(:,34) = 'DEPTH'    = SWAN IVTYPE 4   = DEPLOC
!  out_wwm(:,35) = 'ELEVATION'= SWAN IVTYPE 51  = COMPDA(IG,JWLV2)

          IF (TESTFL.AND.ITEST.GE.60) THEN
            WRITE(PRINTF, 222) IG, iplg(IG), 'U_Cur',  COMPDA(IG,JVX2)
            WRITE(PRINTF, 222) IG, iplg(IG), 'V_Cur',  COMPDA(IG,JVY2)
            WRITE(PRINTF, 222) IG, iplg(IG), 'DEPTH',  DEPLOC
            WRITE(PRINTF, 222) IG, iplg(IG), 'WLV',  COMPDA(IG,JWLV2)
          ENDIF

#if 1  
!            case(27)
!              out_name(ncount_2dnode)='rollerDissRate'  varout_2dnode(icount,:)=rho0*eps_r(1:np)
!            case(28)
!              out_name(ncount_2dnode)='dissRateDepBreaking'  varout_2dnode(icount,:)=wave_sbrtot(1:np) 
!            case(29)
!              out_name(ncount_2dnode)='dissRateBottFriction'   varout_2dnode(icount,:)=wave_sbftot(1:np)
!            case(30)
!              out_name(ncount_2dnode)='dissRateWhiteCapping'   varout_2dnode(icount,:)=wave_sdstot(1:np)
!            case(31)
!              out_name(ncount_2dnode)='energyInputAtmos'    varout_2dnode(icount,:)=wave_sintot(1:np)

          IVTYPE = 55
          DISSURF = COMPDA(IG,JDSXS)
          !out_wwm(IG,36) =  DISSURF
          wave_sbrtot(IG) = DISSURF * GRAV_W * RHO_W ! convert frm m2/s to W/m2
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), wave_sbrtot(IG)
          ENDIF

          IVTYPE = 56
          DISWCAP = COMPDA(IG,JDSXW)
          !out_wwm(IG,37) =  DISWCAP
          wave_sdstot(IG) = DISWCAP * GRAV_W * RHO_W  ! convert frm m2/s to W/m2
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), wave_sdstot(IG)
          ENDIF

          IVTYPE = 54
          DISBOT = COMPDA(IG,JDSXB)
          wave_sbftot(IG) =  DISBOT * GRAV_W * RHO_W  ! convert frm m2/s to W/m2
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), wave_sbftot(IG)
          ENDIF

!          IVTYPE = 8
!          QB = COMPDA(IG,JQB)
!          out_wwm(IG,39) =  QB
!          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
!            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,39)
!          ENDIF
#endif


!JL Add additional Diags
#if 0
!  out_wwm(:,36) = 'DISSURF' = SWAN IVTYPE 55 ! Surface breaking dissipation 'm2/s  COMPDA(IG,JDSXS)' 
!  out_wwm(:,37) = 'DISWCAP' = SWAN IVTYPE 56 ! Whitecapping dissipation' 'm2/s     COMPDA(IG,JDSXW)' 
!  out_wwm(:,38) = 'DISBOT'  = SWAN IVTYPE 54 ! Bottom friction dissipation 'm2/s   COMPDA(IG,JDSXB)'
!  out_wwm(:,39) = 'QB'      = SWAN IVTYPE 8  ! fraction of breaking wave '-'
! ToDo : ! vegetation Dissipation  COMPDA(IG,JDSXV)'
!        ! turbulence Dissipation  COMPDA(IG,JDSXT)'
!        ! mud dissip              COMPDA(IG,JDSXM)'
!        - total dissip (SUM)      COMPDA(IG,   ??  DISSXY) '

          IVTYPE = 55
          DISSURF = COMPDA(IG,JDSXS)
          out_wwm(IG,36) =  DISSURF
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,36)
          ENDIF

          IVTYPE = 56
          DISWCAP = COMPDA(IG,JDSXW)
          out_wwm(IG,37) =  DISWCAP
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,37)
          ENDIF

          IVTYPE = 54
          DISBOT = COMPDA(IG,JDSXB)
          out_wwm(IG,38) =  DISBOT
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,38)
          ENDIF

          IVTYPE = 8
          QB = COMPDA(IG,JQB)
          out_wwm(IG,39) =  QB
          IF (TESTFL.AND.ITEST.GE.60) THEN                                         !40.00
            WRITE(PRINTF, 222) IG, iplg(IG), OVSNAM(IVTYPE), out_wwm(IG,39)
          ENDIF

!  out_wwm(:,40) = 'wavpres'                  ! Wave induced pressure [Pa]
!  out_wwm(:,41) = 'ustx'    =                ! Depth averaged Stokes drift, x-comp 'm/s'
!  out_wwm(:,42) = 'usty'    =                ! Depth averaged Stokes drift, y-comp 'm/s'

          SMS_x=ZERO;   SMS_y=ZERO;   SME_P=ZERO
          ustx=ZERO;    usty=ZERO;  wavpres=ZERO
          eDep = max(DEPLOC,DEPMIN)
          DO IS=1,MSC
            eMult = SPCSIG(IS)*DDIR*DS_INCR(IS)
            CFF   = MIN(30._rkind,WKLOC(IS)*eDep) !=KD
            eWkReal = CFF/eDep
            eSinh2kd = MySINH(2*CFF)
            eSinhkd  = MySINH(CFF)
            eSinhkd2 = eSinhkd**2
            SMS_x  = ZERO
            SMS_y  = ZERO
            DO ID=1,MDC
              eLoc  = ACLOC(ID,IS)*eMult
              SME_P = SME_P + GRAV_W * (CFF/eSinh2kd)*(ONE/eDep) * eLoc
              SMS_x = SMS_x + eLoc * SPCDIR(ID,2)
              SMS_y = SMS_y + eLoc * SPCDIR(ID,3)
            END DO !ID
            eQuot2 = (eSinh2kd/(2*CFF))/eSinhkd2
            eProd2 = SPCSIG(IS)*eWkReal*eQuot2
            ustx = ustx + SMS_x*eProd2 !=STOKESBAROX, Depth averaged Stokes drift, x-comp 'm/s'
            usty = usty + SMS_y*eProd2 !=STOKESBAROY, Depth averaged Stokes drift, y-comp 'm/s'
          END DO ! IS
 
          wavpres = SME_P
          out_wwm(IG,40) = wavpres
          out_wwm(IG,41) = ustx
          out_wwm(IG,42) = usty

          IF (TESTFL.AND.ITEST.GE.60) THEN
            WRITE(PRINTF, 222) IG, iplg(IG), 'WavPres (Pa)', out_wwm(IG,40)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Ustokes (m/s)', out_wwm(IG,41)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Vstokes (m/s)', out_wwm(IG,42)
          ENDIF

!  Other : wave stress determined from the wave dissipation
!  out_wwm(:,43) = 'fsurx'    = DISSURF+DISWCAP,  X-component of the 2-D surface wave dissipation force [m2/s2]
!  out_wwm(:,44) = 'fsury'    = DISSURF+DISWCAP, y comp ..
!  out_wwm(:,45) = 'fbotx'    = DISBOT, X-component of the 2-D bed wave dissipation force [m2/s2]
!  out_wwm(:,46) = 'fboty'    = DISBOT, y comp ..'

  ! wave dissipation forcings
  !        fsurx = ZERO;  fsury = ZERO
  !        fbotx = ZERO;  fboty = ZERO
  !        SME_D = ZERO
  !        IF(DEPLOC.GT.DEPMIN.AND.SPCMAXLOC.GT.ZERO) THEN
  !               !Based on dissipation from the max of spectra, see coherensV2.11/swan41/wavecalculator.ftn90
  !               SME_D = GRAV_W*WKMAXLOC/SPCMAXLOC
  !               fsurx = (DISSURF+DISWCAP)*SME_D*DIRMEANLOC(2)
  !               fsury = (DISSURF+DISWCAP)*SME_D*DIRMEANLOC(3)
  !               fbotx = DISBOT*SME_D*DIRMEANLOC(2)
  !               fboty = DISBOT*SME_D*DIRMEANLOC(3)
  !        ENDIF

          out_wwm(IG,43) = SBR(1,IG) !fsurx
          out_wwm(IG,44) = SBR(2,IG) !fsury
          out_wwm(IG,45) = SBF(1,IG) !fbotx
          out_wwm(IG,46) = SBF(2,IG) !fboty

          IF (TESTFL.AND.ITEST.GE.60) THEN
            WRITE(PRINTF, 222) IG, iplg(IG), 'Wv Surf Diss (X, m2/s)', out_wwm(IG,43)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Wv Surf Diss (Y, m2/s)', out_wwm(IG,44)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Wv Bott Diss (X, m2/s)', out_wwm(IG,45)
            WRITE(PRINTF, 222) IG, iplg(IG), 'Wv Bott Diss (Y, m2/s)', out_wwm(IG,46)
          ENDIF
#endif

!*********************************************************************
!*       WIND Outputs                                                 *
!**********************************************************************
!       out_wwm_windpar(npa,10):
!         1) = WINDXY(IP,1) ! wind vector u10,x
!         2) = WINDXY(IP,2) ! wind vector u10,y
!         3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) ! wind magnitutde u10
!         4) = TAUW(IP)     ! wave stress from the discrete part of the spectra
!         5) = TAUHF(IP)    ! high freq. part of the waves.
!         6) = TAUTOT(IP)   ! total stress of the wave
!         7) = Z0(IP)       ! apparent rougnes lengths (m)
!         8) = UFRIC(IP)    ! ustar - frictional vel. (m/s)
!         9) = ALPHA_CH(IP) ! Charnock Parameter gz0/ustar**2
!        10) = CD(IP)       ! Drag Coefficient



         out_wwm_windpar(IG,1) = COMPDA(IG,JWX2) ! wind vector u10,x
         out_wwm_windpar(IG,2) = COMPDA(IG,JWY2) ! wind vector u10,y
         ! wind magnitutde u10
         out_wwm_windpar(IG,3) = SQRT(COMPDA(IG,JWX2)**2.+ COMPDA(IG,JWY2)**2.) 
         ! 4 wave stress from the discrete part of the spectra,  see SWAN/swancom3.F90 sdwind4
         ! 5 high freq. part of the waves.                       see SWAN/swancom3.F90 sdwind4
         ! 6 total stress of the wave                            see SWAN/swancom3.F90 sdwind4
         out_wwm_windpar(IG,4:6) = ZERO

! JLefevre : USTAR / Zeff / CHARNOCK following the Janssen et al parameterization
!            To use later as input for Ocean-Atm coupling
          IF(IWIND.EQ.4 .AND. IWCAP.EQ.2) THEN
           out_wwm_windpar(IG,7) = COMPDA(IG,JZEL)   ! apparent rougnes lengths [m]
           out_wwm_windpar(IG,8) = COMPDA(IG,JUSTAR) ! ustar [m/s] used in schism_step to replace prescribed wind surface stress
           out_wwm_windpar(IG,9) = COMPDA(IG,JCHARN) ! Charnock Parameter gz0/ustar**2 [-}
          ELSE
           out_wwm_windpar(IG,7:9) = ZERO
          ENDIF
         out_wwm_windpar(IG,10) = COMPDA(IG,JCDRAG) ! SWAN IVTYPE 38    ! Drag coefficient

!**********************************************************************
!*                                                                    *
!**********************************************************************

!----------------------- Energy Spectrum (m^2/hz)-----------------------
          DO ID = 1,MDC
             SPEC_DENSITY(IG,1:MSC) = SPEC_DENSITY(IG,1:MSC) + &
                      2._rkind*PI_W*SPCSIG(1:MSC)*ACLOC(ID,1:MSC)*DDIR
          END DO

       ENDDO  ! IG cycle
       !WRITE(PRINTF,*)
   
 222        FORMAT(' SWOEXA: POINT POINT_GL', 2I7, 2X, A, 1X, E12.4)

# endif
      RETURN
      END SUBROUTINE SWANOUT
!****************************************************************
!
      SUBROUTINE CVCHEK (KGRPNT, XCGRID, YCGRID)                       !30.72
!
!****************************************************************
!
      USE OCPCOMM4                                                     !40.41
      USE SWCOMM2                                                      !40.41
      USE SWCOMM3                                                      !40.41
      USE SWCOMM4                                                      !40.41
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
!  0. Authors
!
!  !30.72: IJsbrand Haagsma
!  !40.13: Nico Booij
!  !40.41: Marcel Zijlema
!
!  1. Updates
!
!            May  96: New subroutine
!  !30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!  !40.13, Mar. 01: messages corrected and extended
!  !40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Checks whether the given curvilinear grid is correct
!     also set the value of CVLEFT.
!
!  3. Method
!
!     Going around a mesh in the same direction the interior
!     of the mesh must be always in the same side if the
!     coordinates are correct
!
!  4. Argument variables
!
!     KGRPNT: input  Array of indirect addressing
!
      INTEGER KGRPNT(MXC,MYC)                                          !30.72
!
!     XCGRID: input  Coordinates of computational grid in x-direction  !30.72
!     YCGRID: input  Coordinates of computational grid in y-direction  !30.72
!
      REAL    XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                      !30.72
!
!
!     5. SUBROUTINES CALLING
!
!        SWRBC
!
!     6. SUBROUTINES USED
!
!        ---
!
!     7. ERROR MESSAGES
!
!        ---
!
!     8. REMARKS
!
!
!     9. STRUCTURE
!
!   ------------------------------------------------------------
!     FIRST = True
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy),
!                    K3 = KGRPNT(ix+1,iy+1)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1),
!                    K3 = KGRPNT(ix,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1),
!                    K3 = KGRPNT(ix,iy)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy),
!                    K3 = KGRPNT(ix+1,iy)
!                 ---------------------------------------------------
!                 If K1>1 and K2>1 and K3>1
!                 Then Det = (xpg(K3)-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (ypg(K3)-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If FIRST
!                      Then Make FIRST = False
!                           If Det>0
!                           Then Make CVleft = False
!                           Else Make CVleft = True
!                      ----------------------------------------------
!                      If ((CVleft and Det<0) or (not CVleft and Det>0))
!                      Then Write error message with IX, IY, ISIDE
!   ------------------------------------------------------------
!
!     10. SOURCE
!
!****************************************************************
!
!
      LOGICAL  FIRST
!
      SAVE IENT
      DATA IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'CVCHEK')
!
!     test output
!
      IF (ITEST .GE. 150 .OR. INTES .GE. 30) THEN
        WRITE(PRINTF,186)
 186    FORMAT(/,' ... Subroutine CVCHEK...',   &
     &  /,2X,'POINT( IX, IY),  INDEX,      COORDX,       COORDY')
        ICON = 0
        DO 5 IIY = 1, MYC
          DO 6 IIX = 1, MXC
            ICON = ICON + 1
            WRITE(PRINTF,7)IIX-1,IIY-1,KGRPNT(IIX,IIY),  &
     &      XCGRID(IIX,IIY)+XOFFS, YCGRID(IIX,IIY)+YOFFS               !30.72 40.13
 6        CONTINUE
 5      CONTINUE
      ENDIF
 7    FORMAT(4X,I5,1X,I5,3X,I4,5X,F10.2,4X,F10.2)
!
      FIRST = .TRUE.
!
      DO 10 IX = 1,MXC-1
        DO 15 IY = 1,MYC-1
          DO 20 ISIDE = 1,4
            IF (ISIDE .EQ. 1) THEN
              IX1 = IX                                                 !40.13
              IY1 = IY                                                 !40.13
              IX2 = IX+1                                               !40.13
              IY2 = IY                                                 !40.13
              IX3 = IX+1                                               !40.13
              IY3 = IY+1                                               !40.13
            ELSE IF (ISIDE .EQ. 2) THEN
              IX1 = IX+1                                               !40.13
              IY1 = IY                                                 !40.13
              IX2 = IX+1                                               !40.13
              IY2 = IY+1                                               !40.13
              IX3 = IX                                                 !40.13
              IY3 = IY+1                                               !40.13
            ELSE IF (ISIDE .EQ. 3) THEN
              IX1 = IX+1                                               !40.13
              IY1 = IY+1                                               !40.13
              IX2 = IX                                                 !40.13
              IY2 = IY+1                                               !40.13
              IX3 = IX                                                 !40.13
              IY3 = IY                                                 !40.13
            ELSE IF (ISIDE .EQ. 4) THEN
              IX1 = IX                                                 !40.13
              IY1 = IY+1                                               !40.13
              IX2 = IX                                                 !40.13
              IY2 = IY                                                 !40.13
              IX3 = IX+1                                               !40.13
              IY3 = IY                                                 !40.13
            ENDIF
            K1  = KGRPNT(IX1,IY1)                                      !40.13
            XC1 = XCGRID(IX1,IY1)                                      !40.13 30.72
            YC1 = YCGRID(IX1,IY1)                                      !40.13 30.72
            K2  = KGRPNT(IX2,IY2)                                      !40.13
            XC2 = XCGRID(IX2,IY2)                                      !40.13 30.72
            YC2 = YCGRID(IX2,IY2)                                      !40.13 30.72
            K3  = KGRPNT(IX3,IY3)                                      !40.13
            XC3 = XCGRID(IX3,IY3)                                      !30.72
            YC3 = YCGRID(IX3,IY3)                                      !30.72
            DET   = 0.
            IF (K1 .GE. 2 .AND. K2 .GE. 2 .AND. K3 .GE. 2) THEN
              DET = ((XC3 - XC1) * (YC2 - YC1)) - &
     &              ((YC3 - YC1) * (XC2 - XC1))
              IF (DET .EQ. 0.) THEN
!               three grid points on one line                          !40.13
                CALL MSGERR (2,'3 comp. grid points on one line')      !40.13
                WRITE (PRINTF, 112) &
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS,  &            !40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS,  &            !40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                !40.13
 112            FORMAT (3(1X, 2I3, 2(1X, F14.4)))                      !40.13
              ENDIF
!
              IF (FIRST) THEN
                FIRST = .FALSE.
                IF (DET .GT. 0.) THEN
                  CVLEFT = .FALSE.
                ELSE
                  CVLEFT = .TRUE.
                ENDIF
              ENDIF
              IF (     (      CVLEFT .AND. DET .GT. 0.)  &
     &            .OR. (.NOT. CVLEFT .AND. DET .LT. 0.)) THEN
!               crossing grid lines in a mesh                          !40.13
                CALL MSGERR (2,'Grid angle <0 or >180 degrees')        !40.13
                WRITE (PRINTF, 112)    &
     &               IX1-1, IY1-1, XC1+XOFFS, YC1+YOFFS, &             !40.13
     &               IX2-1, IY2-1, XC2+XOFFS, YC2+YOFFS, &             !40.13
     &               IX3-1, IY3-1, XC3+XOFFS, YC3+YOFFS                !40.13
              ENDIF
            ENDIF
 20       CONTINUE
 15     CONTINUE
 10   CONTINUE
      RETURN
!     *** end of subroutine CVCHEK ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE CVMESH (XP, YP, XC, YC, KGRPNT, XCGRID ,YCGRID, KGRBND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                     !40.41
      USE SWCOMM2                                                      !40.41
      USE SWCOMM3                                                      !40.41
      USE SWCOMM4                                                      !40.41
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
!  !30.72: IJsbrand Haagsma
!  !40.00, 40.13: Nico Booij
!  !40.02: IJsbrand Haagsma
!  !40.41: Marcel Zijlema
!
!  1. Updates
!
!  !30.21, Jun. 96: New for curvilinear version
!  !30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!  !40.00, May  98: procedure for points outside grid accelerated
!  !40.00, Feb  99: procedure extended for 1D case
!                     XOFFS and YOFFS added in write statements
!  !40.02, Mar. 00: Fixed bug that placed dry testpoints outside computational grid
!  !40.13, Mar. 01: message "CVMESH 2nd attempt .." suppressed
!  !40.41, Oct. 04: common blocks replaced by modules, include files removed
!  !40.41, Nov. 04: search for boundary points improved
!
!  2. Purpose
!
!     procedure to find location in curvilinear grid for a point
!     given in problem coordinates
!
!  3. Method
!
!     First attempt: use Newton-Raphson method to find XC and YC
!     (Note: in the program XC and YC indicate the mesh and position in
!     the mesh) in a few steps; this may be most efficient if a series of
!     points is processed, because the previous point provides a good
!     first estimate.
!     This procedure may fail if the number of iterations is larger than
!     a previously set limit (default=5).
!
!     If the first attempt fails then determine whether the points (XP,YP)
!     is inside the mesh. If so, then the Newton-Raphson procedure is used
!     again with the pivoting point like first guess. Otherwise, scan the
!     boundaries whether the point is on the boundaries. If this fails, it
!     may be concluded that the point (XP,YP) is outside the grid.
!
!  4. Argument variables
!
!     XCGRID  input  Coordinates of computational grid in x-direction  !30.72
!     YCGRID  input  Coordinates of computational grid in y-direction  !30.72
!     XP, YP  input  a point given in problem coordinates
!     XC, YC  outp   same point in computational grid coordinates
!
      REAL     XCGRID(MXC,MYC),    YCGRID(MXC,MYC)                     !30.72
      REAL     XP, YP, XC, YC
!
!     KGRPNT   input   array(MXC,MYC)  grid numbers
!                      if KGRPNT <= 1, point is not in comp. grid.
!     KGRBND   input   lists all boundary grid points consecutively
!
      INTEGER  KGRPNT(MXC,MYC), KGRBND(*)                              !40.00
!
!     Local variables
!
!     MXITNR   number of iterations in Newton-Raphson procedure
!     IX, IY   counter of computational grid point
!     K1       address of grid point
!     IXMIN    counter of grid point closest to (XP,YP)
!     IYMIN    counter of grid point closest to (XP,YP)
!     IBND     counter of boundary grid points
!
      INTEGER       :: IX, IY, K1, IXMIN, IYMIN, IBND
      INTEGER, SAVE :: MXITNR = 0
      INTEGER, SAVE :: IENT = 0
!
!     INMESH   if True, point (XP,YP) is inside the computational grid
!     FINDXY   if True, Newton-Raphson procedure succeeded
!     ONBND    if True, given point is on boundary                     !40.41
!
      LOGICAL  INMESH ,FINDXY, ONBND
!
!     DISMIN   minimal distance found                                  !40.41
!     XPC1     user coordinate of a computational grid point
!     YPC1     user coordinate of a computational grid point
!     XC0      grid coordinate of grid point closest to (XP,YP)
!     YC0      grid coordinate of grid point closest to (XP,YP)
!
      REAL       :: DISMIN                                             !40.41
      REAL       :: XPC1, YPC1, XC0, YC0
!
!  5. SUBROUTINES CALLING
!
!     SINCMP
!
!  6. SUBROUTINES USED
!
!       NEWTON
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!     --------------------------------------------------------------
!     Determine XC and YC from XP and YP using Newton-Raphson iteration
!     process
!     If (XC and YC were found) then
!       Procedure is ready; Return values of XC and YC
!       return
!     else
!     ---------------------------------------------------------------------
!     For ix=1 to MXC-1 do
!         For iy=1 to MYC-1 do
!             Inmesh = True
!             For iside=1 to 4 do
!                 Case iside=
!                 1: K1 = KGRPNT(ix,iy), K2 = KGRPNT(ix+1,iy)
!                 2: K1 = KGRPNT(ix+1,iy), K2 = KGRPNT(ix+1,iy+1)
!                 3: K1 = KGRPNT(ix+1,iy+1), K2 = KGRPNT(ix,iy+1)
!                 4: K1 = KGRPNT(ix,iy+1), K2 = KGRPNT(ix,iy)
!                 ----------------------------------------------------------
!                 If K1>0 and K2>0
!                 Then Det = (xp-xpg(K1))*(ypg(K2)-ypg(K1)) -
!                            (yp-ypg(K1))*(xpg(K2)-xpg(K1))
!                      If ((CVleft and Det>0) or (not CVleft and Det<0))
!                      Then Make Inmesh = False
!                      Else  Inmesh = true and XC = IX and YC = IY
!                 Else Make Inmesh = False
!             --------------------------------------------------------
!             If Inmesh
!             Then Determine XC and YC using Newton-Raphson iteration
!                  process
!                  Procedure is ready; Return values of XC and YC
!     ---------------------------------------------------------------------
!     No mesh is found: Make XC and YC = exception value
!     Return values of XC and YC
!     ---------------------------------------------------------------------
!
!****************************************************************
!
!
      IF (LTRACE) CALL STRACE (IENT,'CVMESH')
!
      IF (ONED) THEN
        CALL NEWT1D  (XP, YP, XCGRID, YCGRID, KGRPNT,   &              !40.00
     &                XC ,YC ,FINDXY)
        IF (.NOT.FINDXY) THEN
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                       !40.00
          ENDIF
        ENDIF
        GOTO 99
      ELSE
!       two-dimensional computation
        XC = 1.
        YC = 1.
!       --- First attempt, to find XC,YC with Newton-Raphson method
        MXITNR = 5
        CALL NEWTON  (XP, YP, XCGRID, YCGRID,    &                     !40.00
     &                MXITNR ,ITER, XC ,YC ,FINDXY)                    !40.41 40.02
        IF ((ITEST .GE. 150 .OR. INTES .GE. 20) .AND. FINDXY) THEN     !40.02
           WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC                 !40.03
        ENDIF
 25     FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4, &
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
        IF (FINDXY) GOTO 80
!
        IF (INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)) THEN
!         --- select grid point closest to (XP,YP)
          DISMIN = 1.E20
          DO 50 IX = 1,MXC
            DO 40 IY = 1,MYC
              K1  = KGRPNT(IX,IY)
              IF (K1.GT.1) THEN
                XPC1 = XCGRID(IX,IY)
                YPC1 = YCGRID(IX,IY)
                DISXY = SQRT ((XP-XPC1)**2 + (YP-YPC1)**2)
                IF (DISXY .LT. DISMIN) THEN
                  IXMIN  = IX
                  IYMIN  = IY
                  DISMIN = DISXY
                ENDIF
              ENDIF
  40        CONTINUE
  50      CONTINUE
!         second attempt using closest grid point as first guess
          MXITNR = 20
          XC0 = REAL(IXMIN)
          YC0 = REAL(IYMIN)
!         ITEST condition changed from 20 to 120                       !40.13
          IF (ITEST.GE.120) WRITE (PRTEST, 55) XP+XOFFS ,YP+YOFFS , &  !40.13
     &          XC0-1. ,YC0-1.
  55      FORMAT (' CVMESH 2nd attempt, (XP,YP)=','(',F12.4,',',F12.4,&
     &          '), (XC,YC)=','(',F9.2,',',F9.2,')')
          DO KORNER = 1, 4
            IF (KORNER.EQ.1) THEN
              XC = XC0 + 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.2) THEN
              XC = XC0 - 0.2
              YC = YC0 + 0.2
            ELSE IF (KORNER.EQ.3) THEN
              XC = XC0 - 0.2
              YC = YC0 - 0.2
            ELSE
              XC = XC0 + 0.2
              YC = YC0 - 0.2
            ENDIF
            CALL NEWTON  (XP, YP, XCGRID, YCGRID,     &                !40.00
     &                    MXITNR ,ITER, XC ,YC ,FINDXY)                !40.41 40.02
            IF (FINDXY) THEN
              IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF,25) XP+XOFFS ,YP+YOFFS ,XC ,YC            !40.00
              ENDIF
              GOTO 80
            ENDIF
          ENDDO
          IF (ITER.GE.MXITNR) THEN                                     !40.41
             WRITE (PRINTF, 75) XP+XOFFS, YP+YOFFS, MXITNR             !40.00
  75         FORMAT (' search for point with location ', 2F12.4,  &    !40.41
     &               ' fails in', I3, ' iterations')                   !40.41
          END IF                                                       !40.41
        ELSE
!         scan boundary to see whether the point is close to the boundary
          DISMIN=99999.                                                !40.41
          ONBND =.FALSE.                                               !40.41
          IX1 = 0                                                      !40.41
          IY1 = 0                                                      !40.51
          IX2 = 0
          DO IBND = 1, NGRBND
            IF (IX2.NE.0) THEN                                         !40.41
               IX1 = IX2
               IY1 = IY2
               XP1 = XP2
               YP1 = YP2
            END IF
            IX2 = KGRBND(2*IBND-1)
            IY2 = KGRBND(2*IBND)
            IF (IX2.NE.0 .AND. (ABS(IX2-IX1).GT.1 .OR.    &            !40.51
     &                          ABS(IY2-IY1).GT.1)) IX1 = 0            !40.51
            IF (IX2.GT.0) THEN
              XP2 = XCGRID(IX2,IY2)
              YP2 = YCGRID(IX2,IY2)
              IF (IBND.GT.1 .AND. IX1.GT.0) THEN                       !40.51
!               --- determine relative distance from boundary segment
!                   with respect to the length of that segment
                SLEN2  = (XP2-XP1)**2 + (YP2-YP1)**2
                RELDIS = ABS((XP-XP1)*(YP2-YP1)-(YP-YP1)*(XP2-XP1)) / &
     &                   SLEN2
                IF (RELDIS.LT.0.01) THEN                               !40.41
!                 --- determine location on the boundary section
                  IF (RELDIS-DISMIN.LE.0.01) THEN                      !40.41
                     DISMIN = RELDIS                                   !40.41
                     RELLOC = ((XP-XP1)*(XP2-XP1)+(YP-YP1)*(YP2-YP1))/ &
     &                        SLEN2
                     IF (RELLOC.GE.-0.001 .AND. RELLOC.LE.1.001) THEN  !40.41
                        RELLCM = RELLOC                                !40.41
                        IF (RELLCM.LT.0.01) RELLCM=0.                  !40.41
                        IF (RELLCM.GT.0.99) RELLCM=1.                  !40.41
                        IX1M  = IX1                                    !40.41
                        IX2M  = IX2                                    !40.41
                        IY1M  = IY1                                    !40.41
                        IY2M  = IY2                                    !40.41
                        ONBND = .TRUE.                                 !40.41
                     ENDIF                                             !40.41
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          IF (ONBND) THEN                                              !40.41
             XC = FLOAT(IX1M) + RELLCM * FLOAT(IX2M-IX1M) - 1.         !40.41
             YC = FLOAT(IY1M) + RELLCM * FLOAT(IY2M-IY1M) - 1.         !40.41
             IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
                WRITE(PRINTF, 65) XP+XOFFS, YP+YOFFS, XC, YC           !40.00
  65            FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4,&
     &                  ') is on the boundary, (XC,YC)=(',      &
     &                  F9.2,',',F9.2,')')
             ENDIF
             GOTO 80
          ENDIF
          XC = -99.
          YC = -99.
          IF (ITEST .GE. 150 .OR. INTES .GE. 20) THEN
            WRITE(PRINTF, 85) XP+XOFFS, YP+YOFFS                       !40.00
  85        FORMAT (' CVMESH: (XP,YP)=','(',F12.4,',',F12.4, &
     &              ') is outside grid')
          ENDIF
          GOTO 99
        ENDIF
      ENDIF                                                            !40.00
  80  IF (KGRPNT(INT(XC+3.001)-2,INT(YC+3.001)-2).LE.1) THEN
         WRITE (PRINTF, 90) XP+XOFFS, YP+YOFFS                         !40.41
  90     FORMAT (' point with location ',2F12.4,' is not active')      !40.41
         XC = -99.
         YC = -99.
      ENDIF
  99  RETURN
      END
!***********************************************************************
!                                                                      *
      LOGICAL FUNCTION INMESH (XP, YP, XCGRID ,YCGRID, KGRBND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                     !40.41
      USE SWCOMM2                                                      !40.41
      USE SWCOMM3                                                      !40.41
      USE M_PARALL                                                     !40.41
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
!     Nico Booij
!  !40.41: Marcel Zijlema
!  !40.51: Marcel Zijlema
!
!  1. Updates
!
!       New function for curvilinear version (ver. 40.00). May '98
!    !40.03, Dec 99: test output added; commons swcomm2 and ocpcomm4 added
!    !40.41, Oct. 04: common blocks replaced by modules, include files removed
!    !40.41, Nov. 04: search for points restricted to subdomain
!    !40.51, Feb. 05: determining number of crossing points improved
!
!  2. Purpose
!
!       procedure to find whether a given location is
!       in the (curvilinear) computational grid
!
!  3. Method  suggested by Gerbrant van Vledder
!
!       draw a line from the point (XP,YP) in vertical direction
!       determine the number of crossings with the boundary of the
!       grid; if this number is even the point is outside
!
!  4. Argument variables
!
!
!     KGRBND   int  input   array containing boundary grid points
!
      INTEGER  KGRBND(*)
!
!     XP, YP    real, input   a point given in problem coordinates
!     XCGRID    real, input   array(IX,IY) x-coordinate of a grid point
!     YCGRID    real, input   array(IX,IY) y-coordinate of a grid point
!
      REAL     XCGRID(MXC,MYC) ,YCGRID(MXC,MYC), &
     &         XP, YP
!
!  5. Parameter variables
!
!  6. Local variables
!
!     NUMCRS   number of crossings with boundary outline
!
      INTEGER  NUMCRS, IX1, IY1, IX2, IY2
      REAL     XP1, XP2, YP1, YP2, YPS, RELDIS, RELDO
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!       CVMESH
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!     --------------------------------------------------------------
!     numcros = 0
!     For all sections of the boundary do
!         determine coordinates of end points (XP1,YP1) and (XP2,YP2)
!         If (XP1<XP and XP2>XP) or (XP1>XP and XP2<XP)
!         then If not (YP1<YP and YP2<YP)
!                   if YPS>YP
!                   then numcros = numcros + 1
!     ---------------------------------------------------------------
!     If numcros is even
!     Then Inmesh = False
!     Else Inmesh = True
!     ---------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      CALL STRACE (IENT,'INMESH')
!
      IF (XP.LT.XCLMIN .OR. XP.GT.XCLMAX .OR.       &                  !40.41
     &    YP.LT.YCLMIN .OR. YP.GT.YCLMAX) THEN                         !40.41
        IF (ITEST.GE.70) WRITE (PRTEST, 22) XP+XOFFS, YP+YOFFS, &      !40.03
     &    XCLMIN+XOFFS, XCLMAX+XOFFS, YCLMIN+YOFFS, YCLMAX+YOFFS       !40.41
  22    FORMAT (1X, 2F12.4, ' is outside region ', 4F12.4)
        INMESH = .FALSE.
        GOTO 90
      ENDIF
!
      IF (NGRBND.LE.0) THEN
        CALL MSGERR (3, 'grid outline not yet determined')
        RETURN
      ENDIF
!
      NUMCRS = 0
      IX1    = 0
      IY1    = 0
      IX2    = 0
      RELDIS = -1.

!     loop over the boundary of the computational grid
      DO IBND = 1, NGRBND
        IF (IX2.NE.0) THEN
           IX1 = IX2
           IY1 = IY2
           XP1 = XP2
           YP1 = YP2
        END IF
        IX2 = KGRBND(2*IBND-1)
        IY2 = KGRBND(2*IBND)
        IF (IX2.NE.0 .AND. (ABS(IX2-IX1).GT.1 .OR.    &                !40.51
     &                      ABS(IY2-IY1).GT.1)) IX1 = 0                !40.51
        IF (IX2.GT.0) THEN
          XP2 = XCGRID(IX2,IY2)
          YP2 = YCGRID(IX2,IY2)
          IF (ITEST.GE.180) WRITE (PRTEST, 28) XP2+XOFFS, &            !40.03
     &    YP2+YOFFS
  28      FORMAT (' boundary point ', 2F12.4)
          IF (IBND.GT.1 .AND. IX1.GT.0) THEN                           !40.51
            IF (((XP1.GT.XP).AND.(XP2.LE.XP)).OR. &
     &          ((XP1.LE.XP).AND.(XP2.GT.XP))) THEN
              IF (YP1.GT.YP .OR. YP2.GT.YP) THEN
!               determine y-coordinate of crossing point
                YPS = YP1 + (XP-XP1) * (YP2-YP1) / (XP2-XP1)
!               determine relative distance from boundary segment      !40.51
!               with respect to the length of that segment             !40.51
                RELDO  = RELDIS                                        !40.51
                RELDIS = ABS(YP-YPS) / SQRT((XP2-XP1)**2 + (YP2-YP1)**2) !40.51
                IF (YPS.GT.YP.AND.ABS(RELDIS-RELDO).GT.0.1) THEN       !40.51
                  NUMCRS = NUMCRS + 1
                  IF (ITEST.GE.70) WRITE (PRTEST, 32) NUMCRS, &        !40.03
     &            XP+XOFFS, YP+YOFFS, YPS+YOFFS                        !40.03
  32              FORMAT (' crossing ', I1, ' point ', 3F12.4)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
!     point is inside the grid is number of crossings is odd
      IF (MOD(NUMCRS,2) .EQ. 1) THEN
        INMESH = .TRUE.
      ELSE
        INMESH = .FALSE.
      ENDIF
  90  RETURN
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWTON (XP, YP, XCGRID, YCGRID,   &                   !40.00
     &                   MXITNR, ITER, XC, YC, FIND)                   !40.41 40.02
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                     !40.41
      USE SWCOMM3                                                      !40.41
      USE SWCOMM4                                                      !40.41
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
!  0. Authors
!
!  !30.72: IJsbrand Haagsma
!  !30.80: Nico Booij
!  !30.82: IJsbrand Haagsma
!  !40.41: Marcel Zijlema
!
!  1. Updates
!
!  !30.21, Jun. 96: New for curvilinear version
!  !30.72, Feb. 98: Introduced generic names XCGRID, YCGRID and SPCSIG for SWAN
!  !30.82, Oct. 98: Updated description of several variables
!  !30.80, Oct. 98: computation of update of XC,YC modified to avoid
!                     division by 0
!  !40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Solve eqs. and find a point  (XC,YC) in a curvilinear grid (compt.
!     grid) for a given point (XP ,YP) in a cartesian grid (problem coord).
!
!  3. Method
!
!     In this subroutine the next equations are solved :
!
!                  @XP             @XP
!     XP(xc,yc) +  --- * @XC   +   --- * @YC  - XP(x,y) = 0
!                  @XC             @YC
!
!                  @YP             @YP
!     YP(xc,yc) +  --- * @XC   +    --- * @YC  - YP(x,y) = 0
!                  @XC             @YC
!
!     In the subroutine, next notation is used for the previous eqs.
!     XVC       + DXDXC * DXC   + DXDYC * DYC - XP  = 0.
!     YVC       + DYDXC * DXC   + DYDYC * DYC - YP  = 0.
!
!
!  4. Argument variables
!
! i   MXITNR: Maximum number of iterations                             !30.82
!
      INTEGER MXITNR, ITER                                             !40.41 40.02
!
!   o XC    : X-coordinate in computational coordinates                !30.82
! i   XCGRID: Coordinates of computational grid in x-direction         !30.72
! i   XP    : X-coordinate in problem coordinates                      !30.82
!   o YC    : Y-coordinate in computational coordinates                !30.82
! i   YCGRID: Coordinates of computational grid in y-direction         !30.72
! i   YP    : Y-coordinate in problem coordinates                      !30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                  !30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                  !30.82
!
!   o FIND  : Whether XC and YC are found                              !30.82
!
      LOGICAL FIND                                                     !30.82
!
!  6. SUBROUTINES USED
!
!     STRACE
!
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      INTEGER, SAVE :: IENT = 0
      IF (LTRACE) CALL STRACE (IENT,'NEWTON')
!
      DXC    = 1000.
      DYC    = 1000.
      TOLDC  = 0.001
      FIND   = .FALSE.
!
      IF (ITEST .GE. 200) THEN
        WRITE(PRINTF,*) ' Coordinates in subroutine NEWTON '
        DO J = 1, MYC
          DO I = 1, MXC
            WRITE(PRINTF,30) I ,J ,XCGRID(I,J) ,YCGRID(I,J)            !30.72
          ENDDO
        ENDDO
      ENDIF
 30   FORMAT(2(2X,I5),2(2X,E12.4))
!
      DO 14 K = 1 ,MXITNR
        ITER = K                                                       !40.41
        I1   = INT(XC)                                                 !40.00
        J1   = INT(YC)
        IF (I1 .EQ. MXC) I1 = I1 - 1
        IF (J1 .EQ. MYC) J1 = J1 - 1
        I2  = I1 + 1
        J2  = J1 + 1
        FJ1 = FLOAT(J1)
        FI1 = FLOAT(I1)
        FJ2 = FLOAT(J2)
        FI2 = FLOAT(I2)
!
        XVC   = (YC-FJ1)*((XC-FI1)*XCGRID(I2,J2)  + &
     &                    (FI2-XC)*XCGRID(I1,J2)) + &
     &          (FJ2-YC)*((XC-FI1)*XCGRID(I2,J1)  + &
     &                    (FI2-XC)*XCGRID(I1,J1))
        YVC   = (YC-FJ1)*((XC-FI1)*YCGRID(I2,J2)  + &
     &                    (FI2-XC)*YCGRID(I1,J2)) + &
     &          (FJ2-YC)*((XC-FI1)*YCGRID(I2,J1)  + &
     &                    (FI2-XC)*YCGRID(I1,J1))
        DXDXC = (YC -FJ1)*(XCGRID(I2,J2) - XCGRID(I1,J2)) + &
     &          (FJ2-YC )*(XCGRID(I2,J1) - XCGRID(I1,J1)) 
        DXDYC = (XC -FI1)*(XCGRID(I2,J2) - XCGRID(I2,J1)) + &
     &          (FI2-XC )*(XCGRID(I1,J2) - XCGRID(I1,J1))
        DYDXC = (YC -FJ1)*(YCGRID(I2,J2) - YCGRID(I1,J2)) + &
     &          (FJ2-YC )*(YCGRID(I2,J1) - YCGRID(I1,J1))
        DYDYC = (XC -FI1)*(YCGRID(I2,J2) - YCGRID(I2,J1)) + &
     &          (FI2-XC )*(YCGRID(I1,J2) - YCGRID(I1,J1))
!
        IF (ITEST .GE. 150) &
     &    WRITE(PRINTF,35) K, XC-1., YC-1., XP, YP, XVC, YVC           !40.00
 35     FORMAT(' NEWTON  iter=', I2, ' (XC,YC)=', 2(1X,F10.2),/, &     !40.00
     &         ' (XP,YP)=', 2(1X,F10.2),                         &
     &         '  X,Y(XC,YC) = ', 2(1X,F10.2))
        IF (ITEST .GE. 180) WRITE(PRINTF,36)                     &
     &     XCGRID(I1,J1), XCGRID(I1,J2), XCGRID(I2,J1), XCGRID(I2,J2), &
     &     YCGRID(I1,J1), YCGRID(I1,J2), YCGRID(I2,J1), YCGRID(I2,J2), &
     &                     DXDXC, DXDYC, DYDXC, DYDYC                  !40.00
 36     FORMAT(' NEWTON grid coord:', 8(1x, F10.0), /  &
     &         '        deriv=', 4(1X,F10.2))                          !40.00
!
!       *** the derivated terms of the eqs. are evaluated and  ***
!       *** the eqs. are solved                                ***
        DDEN = DXDXC*DYDYC - DYDXC*DXDYC                               !30.80
        DXP  = XP - XVC                                                !30.80
        DYP  = YP - YVC                                                !30.80
        IF ( DDEN.NE.0. ) THEN
           DXC = ( DYDYC*DXP - DXDYC*DYP) / DDEN                       !30.80
           DYC = (-DYDXC*DXP + DXDXC*DYP) / DDEN                       !30.80
        ENDIF
!
        XC = XC + DXC
        YC = YC + DYC
!
!       *** If the guess point (XC,YC) is outside of compt. ***
!       *** grid, put that point in the closest boundary    ***
        IF (XC .LT. 1. ) XC = 1.
        IF (YC .LT. 1. ) YC = 1.
        IF (XC .GT. MXC) XC = FLOAT(MXC)
        IF (YC .GT. MYC) YC = FLOAT(MYC)
!
        IF (ITEST .GE. 120 .OR. INTES .GE. 50 .OR. IOUTES .GE. 50) &
     &    WRITE(PRINTF,42) DXC, DYC, XC-1., YC-1.                      !40.00
 42     FORMAT(' (DXC,DYC)=', 2(1X,F10.2), ' (XC,YC)=', 2(1X,F10.2))   !40.00
!
!       *** If the accuracy is reached stop the iteration,  ***
        IF (ABS(DXC) .LE. TOLDC .AND. ABS(DYC) .LE. TOLDC) THEN
!
          FIND = .TRUE.
          XC = XC -1.
          YC = YC -1.
          RETURN
        ENDIF
!
 14   CONTINUE
      RETURN
!     *** end of subroutine NEWTON ***
      END
!***********************************************************************
!                                                                      *
      SUBROUTINE NEWT1D (XP, YP, XCGRID, YCGRID, KGRPNT,    &          !40.00
     &                   XC, YC, FIND)
!                                                                      *
!***********************************************************************
!
      USE OCPCOMM4                                                     !40.41
      USE SWCOMM2                                                      !40.41
      USE SWCOMM3                                                      !40.41
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
!  0. Authors
!
!  !40.00, 40.13: Nico Booij
!  !40.41: Marcel Zijlema
!
!  1. Updates
!
!  !40.00, Feb. 99: New (adaptation from subr NEWTON for 1D case)
!  !40.13, Feb. 01: DX and DY renamed to DELX and DELY (DX and DY are
!                     common var.); error in expression for RS corrected
!                     PRINTF replaced by PRTEST in test output
!  !40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     Finds broken coordinate XC for a given point XP in a rectilinear grid
!
!  3. Method
!
!     In this subroutine the step on the computational grid is selected
!     for which
!
!           (X-X1).(X2-X1)
!     0 <= --------------- <= 1
!          (X2-X1).(X2-X1)
!
!     where X, X1 and X2 are vectors; X corresponds to (Xp,Yp)
!     X1 and X2 are two neighbouring grid points
!
!  4. Argument variables
!
! i   KGRPNT: Grid adresses                                            !40.00
!
      INTEGER KGRPNT(MXC,MYC)                                          !40.00
!
!   o XC    : X-coordinate in computational coordinates                !30.82
! i   XCGRID: Coordinates of computational grid in x-direction         !30.72
! i   XP    : X-coordinate in problem coordinates                      !30.82
!   o YC    : Y-coordinate in computational coordinates                !30.82
! i   YCGRID: Coordinates of computational grid in y-direction         !30.72
! i   YP    : Y-coordinate in problem coordinates                      !30.82
!
      REAL    XC, XCGRID(MXC,MYC), XP                                  !30.82
      REAL    YC, YCGRID(MXC,MYC), YP                                  !30.82
!
!   o FIND  : Whether XC and YC are found                              !30.82
!
      LOGICAL FIND                                                     !30.82
!
!     Local variables:

      REAL :: DELX, DELY   ! grid line                                 !40.13
!
!  6. SUBROUTINES USED
!
!       ---
!
!  7. ERROR MESSAGES
!
!       ---
!
!  8. REMARKS
!
!       ---
!
!  9. STRUCTURE
!
!       -----------------------------------------------------------------
!       -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE     IENT
      DATA     IENT/0/
      IF (LTRACE) CALL STRACE (IENT,'NEWT1D')
!
      IF (ITEST .GE. 120) THEN                                         !40.13
        WRITE(PRTEST,*) ' Coordinates in subroutine NEWT1D '           !40.13
        DO I = 1, MXC
          WRITE(PRTEST,30) I, XCGRID(I,1)+XOFFS ,YCGRID(I,1)+YOFFS     !40.13
        ENDDO
      ENDIF
 30   FORMAT(2X,I5,2(2X,E12.4))
!
      FIND = .FALSE.
      DO 40 IX = 2 ,MXC
        IF (KGRPNT(IX-1,1).GT.1) THEN
          X1 = XCGRID(IX-1,1)
          Y1 = YCGRID(IX-1,1)
        ELSE
          GOTO 40
        ENDIF
        IF (KGRPNT(IX,1).GT.1) THEN
          X2 = XCGRID(IX,1)
          Y2 = YCGRID(IX,1)
        ELSE
          GOTO 40
        ENDIF
!       both ends of the step are valid grid points
!       now verify whether projection of (Xp,Yp) is within the step
        DELX = X2 - X1                                                 !40.13
        DELY = Y2 - Y1                                                 !40.13
        RS = ((XP - X1) * DELX + (YP - Y1) * DELY) /   &               !40.13
     &              (DELX * DELX + DELY * DELY)                        !40.13
        IF (RS.GE.0. .AND. RS.LE.1.) THEN
          FIND = .TRUE.
          XC = REAL(IX-2) + RS                                         !40.00
          YC = 0.
          GOTO 50
        ENDIF
  40  CONTINUE
  50  RETURN
!     *** end of subroutine NEWT1D ***
      END
!*******************************************************************
!                                                                  *
      SUBROUTINE SINTRP (W1, W2, FL1, FL2, FL, SPCDIR, SPCSIG)
!                                                                  *
!*******************************************************************
!
      USE OCPCOMM4                                                      !40.41
      USE SWCOMM3                                                       !40.41
!
!
!   --|-----------------------------------------------------------|--
!     |            Delft University of Technology                 |
!     | Faculty of Civil Engineering, Fluid Mechanics Group       |
!     | P.O. Box 5048,  2600 GA  Delft, the Netherlands           |
!     |                                                           |
!     | Authors :  Weimin Luo, Roeland Ris, Nico Booij            |
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
!  0. Authors
!
!     30.73: Nico Booij
!     30.82: IJsbrand Haagsma
!     40.00: Nico Booij
!     40.41: Marcel Zijlema
!
!  1. Updates
!
!     30.01, Jan. 96: New subroutine for SWAN Ver. 30.01
!     30.73, Nov. 97: revised
!     40.00, Apr. 98: procedure to maintain peakedness introduced
!     30.82, Oct. 98: Update description of several variables
!     30.82, Oct. 98: Made arguments in ATAN2 double precision to prevent
!                     underflows
!     40.41, Oct. 04: common blocks replaced by modules, include files removed
!
!  2. Purpose
!
!     interpolation of spectra
!
!  3. Method (updated...)
!
!     linear interpolation with peakedness maintained
!     interpolated average direction and frequency are determined
!     average direction and frequency of interpolated spectrum are determ.
!     shifts in frequency and direction are determined from spectrum 1 and
!     2 to the interpolated spectrum
!     bilinear interpolation in spectral space is used to calculate
!     contributions from spectrum 1 and 2.
!     in full circle cases interpolation crosses the boundary 0-360 degr.
!
!  4. Argument variables
!
!   o FL    : Interpolated spectrum.
! i   FL1   : Input spectrum 1.
! i   FL2   : Input spectrum 2.
! i   SPCDIR: (*,1); spectral directions (radians)                        30.82
!             (*,2); cosine of spectral directions                        30.82
!             (*,3); sine of spectral directions                          30.82
!             (*,4); cosine^2 of spectral directions                      30.82
!             (*,5); cosine*sine of spectral directions                   30.82
!             (*,6); sine^2 of spectral directions                        30.82
! i   SPCSIG: Relative frequencies in computational domain in sigma-space 30.82
! i   W1    : Weighting coefficient for spectrum 1.
! i   W2    : Weighting coefficient for spectrum 2.
!
      REAL    FL1(MDC,MSC), FL2(MDC,MSC), FL(MDC, MSC)
      REAL    SPCDIR(MDC,6)                                              !30.82
      REAL    SPCSIG(MSC)                                                !30.82
      REAL    W1, W2
!
!  5. Parameter variables
!
!  6. Local variables
!
!     ID       counter of directions
!     IS       counter of frequencies
!
      INTEGER  ID, IS
!
!     DOADD    indicates whether or not values have to be added
!
      LOGICAL  DOADD
!
!     ATOT1    integral over spectrum 1
!     ATOT2    integral over spectrum 2
!     AXTOT1   integral over x-component of spectrum 1
!     AXTOT2   integral over x-component of spectrum 2
!     AYTOT1   integral over y-component of spectrum 1
!     AYTOT2   integral over y-component of spectrum 2
!     ASTOT1   integral over Sigma * spectrum 1
!     ASTOT2   integral over Sigma * spectrum 2
!     ASIG1    average Sigma of spectrum 1
!     ASIG2    average Sigma of spectrum 2
!     DELD1    difference in direction between spectrum 1 and
!              the interpolated spectrum in number of directional steps
!     DELD2    same for spectrum 2
!     DELSG1   shift in frequency between spectrum 1 and interpolated
!              spectrum in number of frequency steps
!     DELSG2   same for spectrum 2
!
      REAL     ATOT1,  ATOT2,  AXTOT1, AXTOT2, AYTOT1, AYTOT2,&
     &         ASTOT1, ASTOT2
      REAL     ASIG1,  ASIG2
      REAL     DELD1,  DELD2,  DELSG1, DELSG2
!
!  8. Subroutines used
!
!  9. Subroutines calling
!
!      SNEXTI, RBFILE
!
! 10. Error messages
!
! 11. Remarks
!
! 12. Structure
!
!      -----------------------------------------------------------------
!      If W1 close to 1
!      Then copy FL from FL1
!      Else If W2 close to 1
!           Then copy FL from FL2
!           Else determine total energy in FL1 and FL2
!                If energy of FL1 = 0
!                Then make FL = W2 * FL2
!                Else If energy of FL2 = 0
!                     Then make FL = W1 * FL1
!                     Else determine average direction of FL1 and FL2
!                          make ADIR = W1 * ADIR1 + W2 * ADIR2
!                          determine average frequency of FL1 and FL2
!                          make ASIG = W1 * ASIG1 + W2 * ASIG2
!                          determine directional shift from FL1
!                          determine directional shift from FL2
!                          determine frequency shift from FL1
!                          determine frequency shift from FL2
!                          For all spectral components do
!                              compose FL from components of FL1 and FL2
!      -----------------------------------------------------------------
!
! 13. Source text
!
      SAVE IENT
      DATA IENT/0/
      CALL STRACE(IENT,'SINTRP')
!
!     interpolation of spectra
!     ------------------------
!
      IF (W1.GT.0.99) THEN
        DO 101 ID=1,MDC
          DO 102 IS=1,MSC
            FL(ID,IS) = FL1(ID,IS)
 102      CONTINUE
 101    CONTINUE
      ELSE IF (W1.LT.0.01) THEN
        DO 201 ID=1,MDC
          DO 202 IS=1,MSC
            FL(ID,IS) = FL2(ID,IS)
 202      CONTINUE
 201    CONTINUE
      ELSE
        ATOT1  = 0.
        ATOT2  = 0.
        AXTOT1 = 0.
        AXTOT2 = 0.
        AYTOT1 = 0.
        AYTOT2 = 0.
        ASTOT1 = 0.
        ASTOT2 = 0.
        DO 301 ID=1,MDC
          DO 302 IS=1,MSC
            ATOT1  = ATOT1  + FL1(ID,IS)
            AXTOT1 = AXTOT1 + FL1(ID,IS) * SPCDIR(ID,2)
            AYTOT1 = AYTOT1 + FL1(ID,IS) * SPCDIR(ID,3)
            ASTOT1 = ASTOT1 + FL1(ID,IS) * SPCSIG(IS)
            ATOT2  = ATOT2  + FL2(ID,IS)
            AXTOT2 = AXTOT2 + FL2(ID,IS) * SPCDIR(ID,2)
            AYTOT2 = AYTOT2 + FL2(ID,IS) * SPCDIR(ID,3)
            ASTOT2 = ASTOT2 + FL2(ID,IS) * SPCSIG(IS)
 302      CONTINUE
 301    CONTINUE
        IF (ATOT1.LT.1.E-9) THEN
          DO 401 ID=1,MDC
            DO 402 IS=1,MSC
              FL(ID,IS) = W2*FL2(ID,IS)
 402        CONTINUE
 401      CONTINUE
        ELSE IF (ATOT2.LT.1.E-9) THEN
          DO 501 ID=1,MDC
            DO 502 IS=1,MSC
              FL(ID,IS) = W1*FL1(ID,IS)
 502        CONTINUE
 501      CONTINUE
        ELSE
!         determine interpolation factors in Theta space
          AXTOT  = W1 * AXTOT1 + W2 * AXTOT2
          AYTOT  = W1 * AYTOT1 + W2 * AYTOT2
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 509)  ATOT1, ATOT2,&                           !40.02
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2
 509        FORMAT (' SINTRP factors ', 8E11.4, /, 15X, 4F7.3)
          ENDIF
!         DELD1 is the difference in direction between spectrum 1 and
!         the interpolated spectrum in number of directional steps
          DELD1  = REAL(ATAN2(DBLE(AXTOT*AYTOT1 - AYTOT*AXTOT1),&         !30.82
     &                        DBLE(AXTOT*AXTOT1 + AYTOT*AYTOT1))) / DDIR  !30.82
!         DELD2 is the difference between spectrum 2 and
!         the interpolated spectrum
          DELD2  = REAL(ATAN2(DBLE(AXTOT*AYTOT2 - AYTOT*AXTOT2),&         !30.82
     &                        DBLE(AXTOT*AXTOT2 + AYTOT*AYTOT2))) / DDIR  !30.82
          IDD1A  = NINT(DELD1)
          RDD1B  = DELD1 - REAL(IDD1A)
          IF (RDD1B .LT. 0.) THEN
            IDD1A = IDD1A - 1
            RDD1B = RDD1B + 1.
          ENDIF
          IDD1B  = IDD1A + 1
          RDD1B  = W1 * RDD1B
          RDD1A  = W1 - RDD1B
          IDD2A  = NINT(DELD2)
          RDD2B  = DELD2 - REAL(IDD2A)
          IF (RDD2B .LT. 0.) THEN
            IDD2A = IDD2A - 1
            RDD2B = RDD2B + 1.
          ENDIF
          IDD2B  = IDD2A + 1
          RDD2B  = W2 * RDD2B
          RDD2A  = W2 - RDD2B
!
!         determine interpolation factors in Sigma space
          ASIG1  = ASTOT1 / ATOT1
          ASIG2  = ASTOT2 / ATOT2
          ATOT   = W1 * ATOT1  + W2 * ATOT2
          ASTOT  = W1 * ASTOT1 + W2 * ASTOT2
          ASIG   = ASTOT / ATOT
!
!         DELSG1 is shift in frequency between spectrum 1 and interpolated
!         spectrum in number of frequency steps
          DELSG1 = ALOG (ASIG1 / ASIG) / FRINTF
          IDS1A  = NINT(DELSG1)
          RDS1B  = DELSG1 - REAL(IDS1A)
          IF (RDS1B .LT. 0.) THEN
            IDS1A = IDS1A - 1
            RDS1B = RDS1B + 1.
          ENDIF
          IDS1B  = IDS1A + 1
          RDS1A  = 1. - RDS1B
!
!         DELSG2 is shift in frequency between spectrum 2 and interpolated
!         spectrum in number of frequency steps
          DELSG2 = ALOG (ASIG2 / ASIG) / FRINTF
          IDS2A  = NINT(DELSG2)
          RDS2B  = DELSG2 - REAL(IDS2A)
          IF (RDS2B .LT. 0.) THEN
            IDS2A = IDS2A - 1
            RDS2B = RDS2B + 1.
          ENDIF
          IDS2B  = IDS2A + 1
          RDS2A  = 1. - RDS2B
!         test output
          IF (ITEST.GE.80) THEN
            WRITE (PRTEST, 510) ATOT, ATOT1, ATOT2,&
     &            AXTOT, AXTOT1, AXTOT2, AYTOT, AYTOT1, AYTOT2,&
     &            DELD1, DELD2, DELSG1, DELSG2
 510        FORMAT (' SINTRP factors ', 9E11.4, /, 15X, 4F7.3)
            WRITE (PRTEST, 512) IDS1A, RDS1A, IDS1B, RDS1B,&
     &            IDS2A, RDS2A, IDS2B, RDS2B,&
     &            IDD1A, RDD1A, IDD1B, RDD1B,&
     &            IDD2A, RDD2A, IDD2B, RDD2B
 512        FORMAT (' SINTRP ', 8(I2, F7.3))
          ENDIF
!
          DO 601 ID=1,MDC
            DO 602 IS=1,MSC
              FL(ID,IS) = 0.
 602        CONTINUE
 601      CONTINUE
          DO 611 ID=1,MDC
            DOADD = .TRUE.
            ID1A = ID + IDD1A
            IF (FULCIR) THEN
              IF (ID1A.LT.1)   ID1A = ID1A + MDC
              IF (ID1A.GT.MDC) ID1A = ID1A - MDC
            ELSE
              IF (ID1A.LT.1)   DOADD = .FALSE.
              IF (ID1A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 612 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD1A * RDS1A * FL1(ID1A,IS+IDS1A)
 612          CONTINUE
              DO 613 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD1A * RDS1B * FL1(ID1A,IS+IDS1B)
 613          CONTINUE
            ENDIF
 611      CONTINUE
          DO 621 ID=1,MDC
            DOADD = .TRUE.
            ID1B = ID + IDD1B
            IF (FULCIR) THEN
              IF (ID1B.LT.1)   ID1B = ID1B + MDC
              IF (ID1B.GT.MDC) ID1B = ID1B - MDC
            ELSE
              IF (ID1B.LT.1)   DOADD = .FALSE.
              IF (ID1B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 622 IS = MAX(1,1-IDS1A), MIN(MSC,MSC-IDS1A)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD1B * RDS1A * FL1(ID1B,IS+IDS1A)
 622          CONTINUE
              DO 623 IS = MAX(1,1-IDS1B), MIN(MSC,MSC-IDS1B)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD1B * RDS1B * FL1(ID1B,IS+IDS1B)
 623          CONTINUE
            ENDIF
 621      CONTINUE
          DO 631 ID=1,MDC
            DOADD = .TRUE.
            ID2A = ID + IDD2A
            IF (FULCIR) THEN
              IF (ID2A.LT.1)   ID2A = ID2A + MDC
              IF (ID2A.GT.MDC) ID2A = ID2A - MDC
            ELSE
              IF (ID2A.LT.1)   DOADD = .FALSE.
              IF (ID2A.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 632 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD2A * RDS2A * FL2(ID2A,IS+IDS2A)
 632          CONTINUE
              DO 633 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD2A * RDS2B * FL2(ID2A,IS+IDS2B)
 633          CONTINUE
            ENDIF
 631      CONTINUE
          DO 641 ID=1,MDC
            DOADD = .TRUE.
            ID2B = ID + IDD2B
            IF (FULCIR) THEN
              IF (ID2B.LT.1)   ID2B = ID2B + MDC
              IF (ID2B.GT.MDC) ID2B = ID2B - MDC
            ELSE
              IF (ID2B.LT.1)   DOADD = .FALSE.
              IF (ID2B.GT.MDC) DOADD = .FALSE.
            ENDIF
            IF (DOADD) THEN
              DO 642 IS = MAX(1,1-IDS2A), MIN(MSC,MSC-IDS2A)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD2B * RDS2A * FL2(ID2B,IS+IDS2A)
 642          CONTINUE
              DO 643 IS = MAX(1,1-IDS2B), MIN(MSC,MSC-IDS2B)
                FL(ID,IS) = FL(ID,IS) + &
     &                      RDD2B * RDS2B * FL2(ID2B,IS+IDS2B)
 643          CONTINUE
            ENDIF
 641      CONTINUE
        ENDIF
      ENDIF
!
!     Test output
      IF (ITEST.GE.80) THEN
        A1 = 0.
        A2 = 0.
        AA = 0.
        DO 801 ID=1,MDC
          DO 802 IS=1,MSC
            A1 = MAX(A1,FL1(ID,IS))
            A2 = MAX(A2,FL2(ID,IS))
            AA = MAX(AA,FL(ID,IS))
 802      CONTINUE
 801    CONTINUE
        WRITE (PRTEST, *) ' SINTRP, maxima ', A1, A2, AA
      ENDIF
!
      RETURN
!  end of subroutine of SINTRP
      END
