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


     subroutine MAKE_MECH_PROP_CASE_005(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
                            vs30, thickness, sub_tag_all		
     real*8              :: ni, VS, VP, Depth_real         
              
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
     if (sub_tag_all.eq.1) then
            rho = 2100
            lambda = 4.557000E+09
            mu = 8.400000E+07
            gamma = 2.094395E-01
    
     ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 2nd Layer
     elseif (sub_tag_all.eq.2) then
            rho = 2100
            lambda = 6.289500E+09
            mu = 2.572500E+08
            gamma = 1.196797E-01 
            
     ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 3rd Layer
     elseif (sub_tag_all.eq.3) then
            rho = 2200
            lambda = 1.189100E+10
            mu = 9.295000E+08
            gamma = 6.444293E-02
    
     ! + MATERIAL INTO THE BEDROCK 
     elseif (sub_tag_all.eq.4) then
            if (zs.ge.-3000.0) then
                VS  = 0.4100*(-zs) + 2190.0000        !VS: S velocity in m/s
                VP  = 0.8100*(-zs) + 3690.0000        !VP: P velocity in m/s 
                rho = 0.0680*(-zs) + 2532.0000        !RHO: MASS DENSITY in kg/m^3
                lambda = rho * (VP**2 - 2*VS**2)
                mu = rho * VS**2
                gamma = 1.6111E-02
            else
                VS  = 0.0050*(-zs) + 3405.0000        !VS: S velocity in m/s
                VP  = 0.0050*(-zs) + 6105.0000        !VP: P velocity in m/s 
                rho = 0.0040*(-zs) + 2724.0000        !RHO: MASS DENSITY in kg/m^3
                lambda = rho * (VP**2 - 2*VS**2)
                mu = rho * VS**2
                gamma = 1.6111E-02
            endif
            
     ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 5th Layer (nu modified)
     elseif (sub_tag_all.eq.5) then
            rho = 2100
            lambda = 1.260000E+08
            mu = 8.400000E+07
            gamma = 2.094395E-01
    
     ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 6th Layer (nu modified)
     elseif (sub_tag_all.eq.6) then
            rho = 2100
            lambda = 3.858750E+08
            mu = 2.572500E+08
            gamma = 1.196797E-01 
            
     ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 7th Layer (nu modified)
     else 
            rho = 2200
            lambda = 1.394250E+09
            mu = 9.295000E+08
            gamma = 6.444293E-02
     endif

                
     end subroutine MAKE_MECH_PROP_CASE_005          
