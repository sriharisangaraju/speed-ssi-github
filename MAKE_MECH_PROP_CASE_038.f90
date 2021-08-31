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



     subroutine MAKE_MECH_PROP_CASE_038(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
                            vs30, thickness, sub_tag_all		
     real*8              :: ni, VS, VP, Depth_real
     real*8              :: VSini, VSfin, Zini, Zfin, Rini, Rfin
              
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
     if ((Depth .ge. 0.0d0).and.(zs_all .ge. 0.0d0)) then    
                 
     ! + MATERIAL INSIDE THE BASIN 
     
		  VS = 500					  
          VP  = 1.9 * VS
          rho = 1900   
                 
          lambda = rho * (VP**2 - 2*VS**2)            
          mu = rho * VS**2  
          qs = 50
          gamma = (3.1415*1.0)/qs
           			  
     else
     ! + MATERIAL INSIDE THE BEDROCK - FIRST LAYER OF CRUSTAL MODEL
     ! Depth_real = zs

         VS = 2100
         VP  = 3400
         rho =  2000         
         lambda = rho * (VP**2 - 2*VS**2)            
         mu = rho * VS**2  
         qs = 200
         gamma = (3.1415*1.0)/qs  
         
     endif
   
     
     
     
     end subroutine MAKE_MECH_PROP_CASE_038              
