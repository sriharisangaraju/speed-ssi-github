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


     subroutine MAKE_MECH_PROP_CASE_018(rho, lambda, mu, gamma, qs, qp, & !outputs
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
                   
     if ((Depth .ge. 0.0d0).and.(zs_all .ge. 0.0d0)) then                                    
         VS = 800 + 10*dsqrt(Depth);
         VP  = 1300.d0 + 15.d0*dsqrt(Depth);                           
         rho = 1800.d0 +  5.d0*dsqrt(Depth);                                 
         mu = rho * VS**2.d0;
         qs = 0.1d0*VS;
         gamma = 4.d0*datan(1.d0)/qs;
     else 
         VS = 2100 
         VP = 3500
         rho = 2200   
         lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);                           
         mu = rho * VS**2.d0;
         qs = 0.1d0*VS;
         gamma = 4.d0*datan(1.d0)/qs;
     endif     



 
     end subroutine MAKE_MECH_PROP_CASE_018   
