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


     subroutine MAKE_MECH_PROP_CASE_019(rho, lambda, mu, gamma, qs, qp, & !outputs
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
     
     
     if (vs30 .le. 300.0d0) then                                                
	vp30 = 1800
     else  
	vp30 = 2000
     endif	
			
					  
     if ((Depth.ge.0.0d0).and.(zs_all.ge.0.0d0)) then  									 
		! + MATERIAL INSIDE THE BASIN 
                                          				 
		 VS = vs30 + (2000.0d0-vs30)*(Depth/1000.0d0)**(0.70d0)                                 	 
		 VP = vp30 + (4500.0d0-vp30)*(Depth/1000.0d0)**(0.70d0) 
		 rho = 2000.d0 + 0.40d0*Depth
                     
					 
		 if (Depth .le. 50.0d0) then                                                
                         QS = 20
	         elseif (Depth .le. 200.0d0) then  
			 QS = 50
		 elseif (Depth .le. 500.0d0) then
			 QS = 100
		 else
			 QS = 150
		 endif				 
					  
		lambda = rho * (VP**2 - 2*VS**2)                                                    
		mu = rho * VS**2                                                        			 
		gamma = (3.1415*(2/3))/QS       ! max freq of 2 Hz                                                 
                                                                      						                  						 
      else  
		! + MATERIAL INSIDE THE BEDROCK         			     
	                                              							 
	         VS =  2000.d0                                                        				                                                        				 
		 VP  = 4500.d0                          
		 rho = 2400.d0                                                        				 
                 QS = 200.d0
		 lambda = rho * (VP**2 - 2*VS**2)                                                    
		 mu = rho * VS**2                                                     			 
		 gamma = (3.1415*(2/3))/QS 				                                                       						         
      endif                                           
                                                      
                                                      
     end subroutine MAKE_MECH_PROP_CASE_019                                   
