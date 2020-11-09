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


     subroutine MAKE_MECH_PROP_CASE_012(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all, f_distance)
                                                      
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
                            vs30, thickness, f_distance
     integer*4           :: sub_tag_all		
     real*8              :: ni, VS, VP, Depth_real                                                 
                                                      
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
     

                         
     if ((Depth.ge.0.0d0).and.(zs_all.ge.0.0d0)) then                                    
     ! + MATERIAL INSIDE THE BASIN 
         if (Depth .le. 150.0d0) then       
             VS = 300.d0                   
             VP = 1500.d0                  
             rho = 1800.d0                  
             lambda = rho * (VP**2 - 2*VS**2)
             mu = rho * VS**2   
             qs = 0.1*VS;           
             gamma = (3.1415*(2.d0/3.d0))/qs !hy: fpeak = 2/3 Hz

         elseif(Depth .gt. 150.d0 .and. Depth .le. f_distance)  then          
             VS  =  300.d0 + 10.d0*(Depth-150.d0)**0.5 
             VP =  1500.d0 + 10.d0*(Depth-150.d0)**0.5
             rho = 1800.d0 +  6.d0*(Depth-150.d0)**0.5
             lambda = rho * (VP**2 - 2*VS**2)
             mu = rho * VS**2              
             qs = 0.1*vs;
             gamma = (3.1415*(2.d0/3.d0))/qs                

         else                                   
             VS =   800.d0 + 15.d0*(Depth-f_distance)**0.5                                 
             VP  = 2000.d0 + 15.d0*(Depth-f_distance)**0.5                           
             rho = 2100.d0 +  4.d0*(Depth-f_distance)**0.5                                 
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                          
             qs = 0.1*vs
             gamma = (3.1415*(2/3))/qs          

         endif                                                                        
                                                                                                                                                     
     else  
     ! + MATERIAL INSIDE THE BEDROCK         
         Depth_real = abs(zs)
         if (Depth_real .le. 1000.0d0) then                   
             VS = 1200.d0                                                                 
             VP  = 2300.d0                          
             rho = 2100.d0                                         
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                 
             gamma = (3.1415*(2/3))/(150.d0)                                                     
             !if(damping_type .eq. 2) then
             !   qs=120.d0; qp=230.d0;
             !endif
          elseif (Depth_real.le.3000.0d0) then                                                       
             VS = 2100.d0                                                                   
             VP  = 3500.d0                          
             rho = 2200.d0     
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                 
             gamma = (3.1415*(2/3))/(200.d0)                                
             !if(damping_type .eq. 2) then
             !   qs=200.d0; qp=400.d0;
             !endif

          elseif (Depth_real.le.6000.0d0) then
             VS = 2750.d0                                                                    
             VP  = 4750.d0                          
             rho = 2400.d0                                           
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                 
             gamma = (3.1415*(2/3))/(250.d0)                                  
             !if(damping_type .eq. 2) then
             !   qs=250.d0; qp=500.d0;
             !endif
           else
             VS = 3670.d0                                                                   
             VP  = 6340.d0                          
             rho = 2800.d0                   
             lambda = rho * (VP**2 - 2*VS**2)                                                    
             mu = rho * VS**2                                
             gamma = (3.1415*(2/3))/(350.d0)                                                
             !if(damping_type .eq. 2) then
             !   qs=350.d0; qp=700.d0;
             !endif
                                
           endif
     endif             

     
     
   
     end subroutine MAKE_MECH_PROP_CASE_012                            
