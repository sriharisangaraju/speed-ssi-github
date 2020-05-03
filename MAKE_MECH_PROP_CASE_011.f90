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


     subroutine MAKE_MECH_PROP_CASE_011(rho, lambda, mu, gamma, qs, qp, & !outputs
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
     
     
     if ((Depth .ge. 0.0d0) .and. (zs_all .ge. 0.0d0)) then                  
                                                 
         if (Depth .lt. 15.0d0) then       
                 
             VS = 270.d0                                 
             ni = 0.45d0                                                           
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS 
             rho = 1700.d0      
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                
             gamma = (3.1415*(2/3))/(70.d0)       
                                                            
         elseif (Depth .lt. 50.0d0) then                                                       
         
             VS = 270.d0+11.5d0*(Depth-15.d0)                                                    
             ni = 0.45 - 0.0025*(Depth-15.d0)                                                    
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS
             rho = 1700.d0 + 5.d0*(Depth-15.d0)                                                  
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2  
             gamma = (3.1415*(2/3))/(70.d0+0.5d0*(Depth-15.d0))                                  
             
         else                                                                                
         
             VS= 270.d0 + 11.5d0*(50.d0-15.d0) + 0.7d0*(Depth-50)                                
             ni= 0.45d0 - 0.0025d0*(50.d0-15.d0) - 0.000075d0*(Depth-50)                         
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
             rho= 1700.d0 - 5.d0*(50.d0-15.d0) +0.5d0*(Depth-50)                                 
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                                       
             gamma = (3.1415*(2/3))/(70.d0 + 0.5d0*(50.d0-15.d0) + 0.0775d0*(Depth-50))          
             
         endif                                                                                      
                                                                                                                                                    
     else          
         Depth_real = abs(zs)
                                                                                                                            
         if (Depth_real .lt. 15.0d0) then               
         
             VS = 750.d0                                                                 
             ni = 0.30d0                                                        
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
             rho = 2000.d0                                             
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                
             gamma = (3.1415*(2/3))/(100.d0)
                                                                  
         elseif (Depth_real .lt. 50.0d0) then                                                       
         
             VS = 750.d0+14.d0*(Depth_real-15.d0)              
             ni = 0.30d0 - 0.0005d0*(Depth_real-15.d0)  
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
             rho = 2000.d0 + 6.5d0*(Depth_real-15.d0) 
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2       
             gamma = (3.1415*(2/3))/(100.d0+0.8d0*(Depth_real-15.d0))  
             
         else              
         
             VS= 750.d0 + 14.d0*(50.d0-15.d0) + 1.1d0*(Depth_real-50)         
             ni= 0.30d0 - 0.0005d0*(50.d0-15.d0) - 0.000022d0*(Depth_real-50) 
             VP  =  (((2.0d0*(1.0d0 - ni))/(1.0d0-2.0d0*ni))**0.5)*VS                            
             rho= 2000.d0 + 6.5d0*(50.d0-15.d0) + 0.26d0*(Depth_real-50)
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2  
             gamma = (3.1415*(2/3))/(100.d0 + 0.8d0*(50.d0-15.d0) + 0.049d0*(Depth_real-50))
             
         endif                                                  
     endif                                                                      
                                                      
                                                      
     end subroutine MAKE_MECH_PROP_CASE_011                                     
