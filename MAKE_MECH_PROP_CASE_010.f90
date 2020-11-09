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


     subroutine MAKE_MECH_PROP_CASE_010(rho, lambda, mu, gamma, qs, qp, & !outputs
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
         if (Depth.le.300.0d0) then       ! 0 < z < 300         
             rho =  1700.0d0                                                             
             lambda = 2.9700E+08       !Vs = 300 m/s               
             mu = 1.5300E+08           !Vp = 596 m/s               
             gamma = 2.9920E-02        ! Qs = 70 (2 Hz)            
                                                 
         elseif ((Depth.gt.300.0d0).and.(Depth.le.700.0d0)) then  ! 300 < z < 700     
             rho =  2000.0d0                                                     
             lambda = 3.0000E+09      !Vs = 1000 m/s      
             mu = 2.0000E+09          !Vp = 1871 m/s      
             gamma = 2.0944E-02       ! Qs = 100 (2 Hz)          
                              
         elseif ((Depth.gt.700.0d0)) then ! 700 < z < 1500    
             rho =  2300.0d0                                                                 
             lambda = 1.1178E+10      !Vs = 1800 m/s                
             mu = 7.4520E+09          !Vp = 3368 m/s               
             gamma = 2.0944E-02 ! Qs = 100 (2 Hz)                              
         endif                                                                     
     else                            
         ! + MATERIAL INSIDE THE BEDROCK (Vs=3175 m/s)           
         lambda = 2.6217E+10                                
         mu = 2.6217E+10                                    
         rho = 2600.0d0                                     
         gamma = 1.0472E-02      
     endif             
 
                                                      
                                                      
     end subroutine MAKE_MECH_PROP_CASE_010                                    
