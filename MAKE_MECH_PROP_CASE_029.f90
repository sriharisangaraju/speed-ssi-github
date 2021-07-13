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


     subroutine MAKE_MECH_PROP_CASE_029(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
                                                      
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
                            vs30, thickness 
     integer*4           :: sub_tag_all		
     real*8              :: ni, VS, VP, Depth_real         
     real*8              :: VSini, VSfin, Zini, Zfin, Rini, Rfin

         
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

                ! generic bedrock outcrop with Vs30 = 1500 m/s - from Cotton et al. 2006
     if (Depth .le. 1) then
     VS = 1144  
     elseif(Depth .le. 30) then
     VSini = 1144
     VSfin = 1696
     Zini = 1
     Zfin = 30
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s
     elseif(Depth .le. 190) then
     VSini = 1696
     VSfin = 2381
     Zini = 30
     Zfin = 190
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s
     elseif(Depth .le. 4000) then
     VSini = 2381
     VSfin = 3454
     Zini = 190
     Zfin = 4000
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s 
     else 
     VS = 3440
     endif  
        
      VP  =  2.25*VS                                             !VP: P velocity in m/s - (Poisson  = 0.3 approx) 


     if (Depth .le. 100) then
      Rini = 2200;
      Rfin = 2400;
      Zini = 0;
      Zfin = 100;
      rho = Rini + (Rfin-Rini)*((Depth-Zini)/(Zfin-Zini))**0.50     !RHO: MASS DENSITY in kg/m^3       
     elseif(Depth .le. 1000) then
      Rini = 2400;
      Rfin = 2700;
      Zini = 100;
      Zfin = 1000;
      rho = Rini + (Rfin-Rini)*((Depth-Zini)/(Zfin-Zini))**0.50     !RHO: MASS DENSITY in kg/m^3   
     else
     rho = 2700.0d0
     endif
         


      qs = 0.1d0*VS                                              							 


!	         VS =  2000.d0                                                        				                                                        				 
!		 VP  = 4500.d0                          
!		 rho = 2400.d0                                                        				 
!                 qs = 200.d0

		 lambda = rho * (VP**2 - 2*VS**2)                                                    
 		 mu = rho * VS**2                                                     			 
		 gamma = (3.1415*(2/3))/QS 	

			                                                       						         
      endif                                           
                                                      
                                                      
     end subroutine MAKE_MECH_PROP_CASE_029                                   
