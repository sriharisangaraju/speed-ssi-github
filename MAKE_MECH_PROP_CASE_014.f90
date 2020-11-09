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


     subroutine MAKE_MECH_PROP_CASE_014(rho, lambda, mu, gamma, qs, qp, & !outputs
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
                                                
     if ((Depth .ge. 0.0d0) .and. (zs_all .ge. 0.0d0)) then    
                           
        ! + MATERIAL INSIDE THE ALLUVIAL BASIN 
        if (Depth .le. 100.0d0) then 
            VS = 300.d0                                                     
            VP = 520.d0                                        
            rho = 2200.d0                                                
            lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
            mu = rho * VS**2  
            gamma = (3.1415*(3.d0/3.d0))/(30.d0)                                                       
                            
        elseif (Depth .le. 200.0d0) then 
            VS = 400.d0
            VP = 700.d0                          
            rho = 2300.d0 
            lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
            mu = rho * VS**2
            gamma = (3.1415*(2.d0/3.d0))/(40.d0)  
                            
        elseif (Depth .le. 500.0d0) then                           
            VS = 500.d0     
            VP = 850.d0                          
            rho = 2400.d0                                  
            lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
            mu = rho * VS**2                          
            gamma = (3.1415*(2.d0/3.d0))/(50.d0)
                                                 
        else 
            VS = 1000.d0                        
            VP  = 1700.d0                          
            rho = 2400.d0       
            lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
            mu = rho * VS**2                          
            gamma = (3.1415*(2.d0/3.d0))/(100.d0)                                                 
        endif                
                                       
     else                                                    
        ! + MATERIAL INSIDE THE BEDROCK (Vs=1500m/s)            
        VS = 1500.d0                                                     
        VP = 2800.d0                                        
        rho = 2400.d0                                                
        lambda = rho * (VP**2 - 2.d0*VS**2)                                                    
        mu = rho * VS**2                                                                                 
        gamma = (3.1415*(3.d0/3.d0))/(150.d0)                                                               
     endif                                   
     
 
     end subroutine MAKE_MECH_PROP_CASE_014        
