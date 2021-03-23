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


     subroutine MAKE_MECH_PROP_CASE_006(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all) 
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
                            vs30, thickness
     integer*4           :: sub_tag_all		                            		
     real*8              :: ni, VS, VP, Depth_real         
              
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     


    if ((Depth .ge. 0.0d0) .and. (zs_all.ge.0.0d0)) then    
            ! Option #1
            ! VS  = 300.0d0 + 6.50d0*Depth  !VS: S velocity in m/s     
            ! VP  =  600.0d0 + 12.0d0*Depth         !VP: P velocity in m/s                  
            ! rho =  1900.0d0 + 1.25*Depth        !RHO: MASS DENSITY in kg/m^3            
            
            ! Option #2
             VS  = 300.0d0 + 30.0d0*(Depth)**0.67  !VS: S velocity in m/s     
             VP  =  600.0d0 + 30.0d0*(Depth)**0.77         !VP: P velocity in m/s                  
             rho =  1900.0d0 + 1.25*Depth        !RHO: MASS DENSITY in kg/m^3
             lambda = rho * (VP**2 - 2*VS**2)                   
             mu = rho * VS**2               
             gamma = 3.1416E-02 ! (Qs=100)   
    else            
            ! + MATERIAL INSIDE THE BEDROCK (Vs=1500m/s)            
             lambda = 1.0350E+10            
             mu = 5.1750E+09                
             rho = 2300.0d0                 
             gamma = 2.0944E-02             
    endif           
     
     end subroutine MAKE_MECH_PROP_CASE_006              
