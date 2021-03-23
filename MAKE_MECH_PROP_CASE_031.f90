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



     subroutine MAKE_MECH_PROP_CASE_031(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
        vs30, thickness
     integer*4           :: sub_tag_all		        		
     real*8              :: ni, VS, VP, Depth_real    
     real*8, dimension(1) :: val1     
              
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
    !-------------------------------------------------------------------
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 1st Layer
    ! from ground surface to PE_B
                
    if (sub_tag_all.eq.1) then
            
            !val1(1) =  maxval((/ 200.d0, 90.d0 + 80.d0*dabs(zs)**0.30/))
            VS = 150 + 40*dabs(zs)**(0.4);
            VP  = VS*4.5d0;
            rho = 1800.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;   !CHECK f0 = 5Hz => FMAX in file.mate           

    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 2nd Layer
    ! from PE_B to NU_B
    elseif (sub_tag_all.eq.2) then
            val1(1) = minval((/ 600.d0,360.d0+0.68*dabs(zs) /));
            VS = val1(1);
            VP  = VS * 3.8d0;
            rho = 2050.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            


   ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 3rd Layer
   ! from NU_B to NS_B
    elseif (sub_tag_all.eq.3) then
            VS = 600.d0;
            VP  = VS * 3.2d0;
            rho = 2050.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;       
             

    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 4th Layer
    ! CK
    elseif (sub_tag_all.eq.4) then
            !VS = 365.d0 + 1.15*dabs(zs);
            VS = 2090.d0;
            if (dabs(zs) .le. 1200.d0) VS = 1515.d0 
            
            VP  = VS * 2.d0;
            rho = 2400.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 5th Layer
    ! ZE
    elseif (sub_tag_all.eq.5) then
            VS = 2850.d0;
            VP  = 5100.d0;
            rho = 2450.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 6th Layer
    ! RO (reservoir)
    elseif (sub_tag_all.eq.6) then
            VS = 2300.d0;
            VP  = 3900.d0;
            rho = 2450.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = qp/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    else ! LI (halfspace)
            VS = 2600.d0;
            VP  = 4500.d0;
            rho = 2650.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;
    endif

    !
    !-------------------------------------------------------------------
     
     
     
     end subroutine MAKE_MECH_PROP_CASE_031            
