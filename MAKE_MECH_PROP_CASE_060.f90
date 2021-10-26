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
!    along with SPEED.  If not, see <https://urldefense.com/v3/__http://www.gnu.org/licenses/__;!!JTSHVUr6R1OOzg!aDbB_NDNNqzY5Rv71QLDFP6IQtf6PDbmxUO5-QKFq4zgkjrbIcf8DjlIuTdJhwQR9vU$ >.


!> @brief Makes not-honoring technique. Mechanical properties given node by node.


     subroutine MAKE_MECH_PROP_CASE_060(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
            
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all, vs30, thickness
     integer*4           :: sub_tag_all		                            		
     real*8              :: ni, VS, VP, Depth_real, xLeft, xRight, yUp, yDown 
!     real*8, parameter   :: bufferLayer  = 6000.0d0   
            
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
     !Left  Top    629160.875000  9354624.000000   -32.000000
     !Left  Down   628957.563000  9262946.000000   443.000000
     !Right Top    739448.625000  9354298.000000   -28.000000
     !Right Down   739071.438000  9262575.000000   406.000000 
     
     xLeft  = dabs(xs -  628957.563d0)
     xRight = dabs(xs -  739448.625d0)
     yUp    = dabs(ys - 9354624.000d0)
     yDown  = dabs(ys - 9262575.000d0)

                                
      if ((Depth .ge. 0.0d0).and.(zs_all .ge. 0.0d0)) then
     					if( (xLeft.le.6000) .or. (xRight.le.6000) .or. (yUp.le.6000) .or. (yDown.le.6000) ) then
      							! MATERIAL INSIDE THE BASIN MODIFIED AS OUTCROPPING BEDROCK 
                	       VS  = 2300.d0
				           VP  = 4100.d0
				           rho = 2000.d0
				           qs  = 150
				           qp  = 283 
     					else    
		                        ! MATERIAL INSIDE THE BASIN 
		                   VS  = 582.d0
		                   VP  = 1600.d0 
                           rho = 1200.d0            
                           qs  = 25
                           qp  = 44  
                        endif      
       else  
       ! OUTCROPPING BEDROCK 
              VS  = 2300.d0
              VP  = 4100.d0
              rho = 2000.d0
              qs  = 150
              qp  = 283            
      endif     
      
      lambda = rho * (VP**2 - 2*VS**2)
      mu = rho * VS**2 
      gamma = (3.1415*(2.d0/3.d0))/qs     

     end subroutine MAKE_MECH_PROP_CASE_060        
