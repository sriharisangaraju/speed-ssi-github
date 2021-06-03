!    Copyright (C) 2012 The SPEED FOUNDATION
!    Author: Ilario Mazzieri
!
!    This file is part of SPEED.
!
!    SPEED is free software; you can redistribute it and/or modify it
!    under the terms of the GNU Affero General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    SPEED is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with SPEED.  If not, see <http://www.gnu.org/licenses/>.


!> @brief Makes not-honoring technique. Mechanical properties given node by node.
! Kumamoto 1D velocity model by Kobayashi et al 2017
! There is no interface for basin. Topography layer is just translated to a 
! depth depending on thickness of first 4 layers in 1D velocity model

     subroutine MAKE_MECH_PROP_CASE_046(rho, lambda, mu, gamma, qs, qp, & 
                                        xs, ys, zs, Depth, zs_all)
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all
     real*8              :: ni, VS, VP, xabs, yabs         
     real*8, dimension(1) :: val1              

     ! Properties of bottom most velocity layer in this block
     ! Also Same properties are used near absorbing boundaries
     ! as poisson's ratio of top 3 layers is high
     rho    = 2600.d0;
     lambda = rho * (4000.d0**2.d0 - 2.d0*2200.d0**2.d0);
     mu     = rho * 2200.d0**2.d0;
     qs     = 400.d0;
     qp     = 680.d0;
     gamma  = 4.d0*datan(1.d0)*5.d0/qs;
     
    !-------------------------------------------------------------------
    ! Checking if node falls near Absorbing Boundaries,
    ! element size ~ 50m @ top. So buffer to Absorbing boundary is
    ! 18 elements (900m ~ 450m*2 larger elements at bottom)
    xabs = MIN((xs - 13000.d0), (66100.d0 - xs));
    yabs = MIN((ys - 23000.d0), (68900.d0 - ys));

					
    ! Material properties of top 3 layers are assigned based on 
    ! depth of node from topography surface.
    if ((Depth .ge. 0.0d0) .and. (zs_all .ge. 0.0d0)) then     
      ! Node lies between topography surface (XYZ.out) and
      ! bottom surface (basin surface / ALL.out)

        if ((xabs .le. 900.d0) .or. (yabs .le. 900.d0) .or. (Depth .gt. 300.0d0)) then 
          return
        endif

       if ((Depth .ge. 0.0d0) .and. (Depth .le. 50.0d0)) then 
            !Depth from Topography Surface  <= 50m
            VS = 500.d0;
            VP  = 1900.d0;
            rho = 1800.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 100.d0;
            qp = 170.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs; 

       elseif ((Depth .gt. 50.0d0) .and. (Depth .le. 150.0d0)) then 
            !Depth from Topography Surface  50m to 150m
            VS = 900.d0;
            VP  = 2400.d0;
            rho = 2100.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 180.d0;
            qp = 306.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;       
             
       elseif ((Depth .gt. 150.0d0) .and. (Depth .le. 300.0d0)) then 
            !Depth from Topography Surface 150m to 300m
            VS = 1500.d0;
            VP  = 3400.d0;
            rho = 2500.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 300.d0;
            qp = 510.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;
       endif

     endif
     
     end subroutine MAKE_MECH_PROP_CASE_046          
