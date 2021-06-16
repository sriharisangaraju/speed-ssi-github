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



     subroutine MAKE_MECH_PROP_CASE_028(rho, lambda, mu, gamma, qs, qp, & !outputs
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
     !! NEWER VS-RULE
	      if (Depth .le. 150) then
	         VS = 281.64 + (2.0000*(548.33-281.64))/(1.0000+(15.0000/(Depth+0.1000))**1.2900)
		  else
		     VS = 975
		  endif	
		  					  
          VP  = 1.855 * VS
          rho = 1900 + VS / 1700 * 600          
          lambda = rho * (VP**2 - 2*VS**2)            
          mu = rho * VS**2  
          qs = 50
          gamma = (3.1415*(2/3))/qs
           
						  
     else
     ! + MATERIAL INSIDE THE BEDROCK - FIRST LAYER OF CRUSTAL MODEL
     ! Depth_real = zs
     
      ! generic bedrock outcrop with Vs_min = 1000 m/s and Vs30 = 1200 m/s - from Cotton et al. 2006
     if (Depth .le. 1) then
     VS = 1000  
     elseif(Depth .le. 30) then
     VSini = 1000
     VSfin = 1421
     Zini = 1
     Zfin = 30
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s
     elseif(Depth .le. 190) then
     VSini = 1421
     VSfin = 2217
     Zini = 30
     Zfin = 190
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s
     elseif(Depth .le. 1000) then
     VSini = 2217
     VSfin = 2600
     Zini = 190
     Zfin = 1000
     VS  = VSini + (VSfin-VSini)*((Depth-Zini)/(Zfin-Zini))**0.50     !VS: S velocity in m/s 
     else 
     VS = 2600
     endif  
        
      VP  =  sqrt(3.46)*VS                                             !VP: P velocity in m/s - (Poisson  = 0.3 approx) 

      if (Depth .le. 100) then
      Rini = 2200;
      Rfin = 2500;
      Zini = 0;
      Zfin = 100;
      rho = Rini + (Rfin-Rini)*((Depth-Zini)/(Zfin-Zini))**0.50     !RHO: MASS DENSITY in kg/m^3       
     elseif(Depth .le. 1000) then
      Rini = 2500;
      Rfin = 2840;
      Zini = 100;
      Zfin = 1000;
      rho = Rini + (Rfin-Rini)*((Depth-Zini)/(Zfin-Zini))**0.50     !RHO: MASS DENSITY in kg/m^3   
     else
     rho = 2840.0d0
     endif
         
     qs = 0.1d0*VS
     lambda = rho * (VP**2 - 2*VS**2)                   
     mu = rho * VS**2                                    
     gamma = (3.1415*(2/3))/qs   ! (Qs=Vs/10 f=0.67 Hz)
 
 
         ! generic bedrock outcrop with Vs_min = 800 m/s and Vs30 = 900 m/s - SUPERSEDED
 !        VS  = 800.0d0 + (1700.0d0-800.0d0)*(0.001d0*Depth)**0.50     !VS: S velocity in m/s
 !        VP  =  sqrt(3.46)*VS                                         !VP: P velocity in m/s - (Poisson  = 0.3 approx) 
 !        rho =  2200.0d0 + (2500.0d0-2200.0d0)*(0.001d0*Depth)**0.50  !RHO: MASS DENSITY in kg/m^3 
 !        qs = 0.1d0*VS
 !        lambda = rho * (VP**2 - 2*VS**2)                   
 !        mu = rho * VS**2                                    
 !        gamma = (3.1415*(2/3))/qs   ! (Qs=Vs/10 f=0.67 Hz)

     !    VS = 1700
     !    VP  = 3160
     !    rho =  2500         
     !    lambda = rho * (VP**2 - 2*VS**2)            
     !    mu = rho * VS**2  
     !    qs = 200
     !    gamma = (3.1415*(2/3))/qs  
     endif
   
     
     
     
     end subroutine MAKE_MECH_PROP_CASE_028              
