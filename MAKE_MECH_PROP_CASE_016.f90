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


     subroutine MAKE_MECH_PROP_CASE_016(rho, lambda, mu, gamma, qs, qp, & !outputs
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
     
     
     if( vs30 .lt. 325.d0) then
         if ( Depth .le. 160.d0) then
            VS = 250.d0 + 43.d0*Depth**(0.5d0);
            VP  = 700.d0 + 45.d0*Depth**(0.5d0);
            rho = 1800.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1d0*VS;
            gamma = 4.d0*datan(1.d0)/qs;
         elseif(Depth .le. 2000.d0) then
            VS = 800.d0 + 37.13d0*(Depth-160d0)**(0.5)
            VP  = VS*1.6;
            rho = 1800 + 12.92d0*(Depth-160d0)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1d0*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        endif
     
     elseif (vs30 .lt. 500.d0) then
     
        if ( Depth .le. 80.d0) then
            VS = 325.d0 + 30.74*Depth**(0.5);
            VP  = 800.d0 + 42 *Depth**(0.5);;
            rho = 1850.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;               
        elseif(Depth .le. 120.d0) then
            VS = 600.d0 + 31.62*(Depth-80.d0)**(0.5);
            VP  = 1175.d0 + 26.72*(Depth-80.d0)**(0.5);
            rho =  1850.d0
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        elseif(Depth .le. 250.d0) then
            VS = 800.d0 + 40.75*(Depth-120.d0)**(0.5);
            VP  = VS*1.6;
            rho = 1850.d0 + 9.64*(Depth-120.d0)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;  
        elseif(Depth .le. 2000.d0) then
            VS = 700.d0 + 38.14*(Depth-30.d0)**(0.5);
            VP  = VS*1.6;
            rho = 1960.d0 + 9.43*(Depth-250.d0)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;  
            
        else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        endif          
     
     elseif (vs30 .lt. 700.d0) then
     
        if ( Depth .le. 50.d0) then
            VS = 500.d0 + 42.42*Depth**(0.5);
            VP  = 900.d0 + 42.42*Depth**(0.5);
            rho = 1900.d0;
            lambda = rho * (VP**2 - 2*VS**2);
            mu = rho * VS**2;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        elseif ( Depth .le. 250.d0) then
            VS = 800.d0 + 33.1*(Depth-50.d0)**(0.5);
            VP  = VS*1.6;
            rho = 1900.d0 + 4.89*(Depth-50.d0)**(0.5);
            lambda = rho * (VP**2 - 2*VS**2);
            mu = rho * VS**2;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        
       elseif(Depth .le. 2000.d0) then
            VS = 700.d0 + 38.14*(Depth-30.d0)**(0.5);
            VP  = VS*1.6;
            rho = 1960.d0 + 9.43*(Depth-250.d0)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        endif

     elseif (vs30 .lt. 900.d0) then
     
        if ( Depth .le. 4000.d0) then
            VS = 700.d0 + 37.9*(Depth)**(0.5)
            VP  = VS*1.6;
            rho = 1960.d0 + 8.885*(Depth)**(0.5) 
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        endif
                
     elseif (vs30 .lt. 1350.d0) then
     
         if ( Depth .le. 2000.d0) then
            VS = 900.d0 + 33.38 * (Depth)**(0.5);
            VP  = VS*1.6;
            rho = 2050.d0 + 215.1*(Depth*0.001)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*atan(1.d0)/qs;
         else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;
        endif
 
     else
            VS = 1350.d0 + 23.33*(Depth)**(0.5);
            VP = VS*1.6;
            rho =  2100 + 5.69*(Depth)**(0.5);
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = 0.1*VS;
            gamma = 4.d0*datan(1.d0)/qs;            
          
     endif
     
     ! specific for 2019.09.26 earthquake -- validation
     !stat_id1_x = 6.3269e+05; stat_id1_y = 4.5422e+06;
     !stat_id2_x = 6.5450e+05; stat_id2_y = 4.5376e+06;
     !     
     !if( ((xs(ic)-stat_id1_x)**2.d0 + (ys(ic)-stat_id1_y)**2.d0) .le. 3000 .or. &
     !    ((xs(ic)-stat_id2_x)**2.d0 + (ys(ic)-stat_id2_y)**2.d0) .le. 3000 ) then
     !   if ( Depth .le. 150.d0) then
     !       VS = 600.d0;
     !       VP  = 900.d0;
     !       rho = 1900.d0;
     !       lambda = rho * (VP**2 - 2*VS**2);
     !       mu = rho * VS**2;
     !       qs = 0.1*VS;
     !       gamma = 4.d0*datan(1.d0)/qs;
     !   elseif ( Depth .le. 250.d0) then
     !       VS = 800.d0 + 33.1*(Depth-50.d0)**(0.5);
     !       VP  = VS*1.6;
     !       rho = 1900.d0 + 4.89*(Depth-50.d0)**(0.5);
     !       lambda = rho * (VP**2 - 2*VS**2);
     !       mu = rho * VS**2;
     !       qs = 0.1*VS;
     !       gamma = 4.d0*datan(1.d0)/qs;
     !   
     !  elseif(Depth .le. 2000.d0) then
     !       VS = 700.d0 + 38.14*(Depth-30.d0)**(0.5);
     !       VP  = VS*1.6;
     !       rho = 1960.d0 + 9.43*(Depth-250.d0)**(0.5);
     !       lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
     !       mu = rho * VS**2.d0;
     !       qs = 0.1*VS;
     !       gamma = 4.d0*datan(1.d0)/qs;
     !   else
     !       VS = 1350.d0 + 23.33*(Depth)**(0.5);
     !       VP = VS*1.6;
     !       rho =  2100 + 5.69*(Depth)**(0.5);
     !       lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
     !       mu = rho * VS**2.d0;
     !       qs = 0.1*VS;
     !       gamma = 4.d0*datan(1.d0)/qs;
     !   endif
     !
     ! 
     !endif        

    
     !left
     if(dabs(xs - 576059.d0) .le. 2000.d0) then 
         VS = 3490.d0;
         VP = 5770.d0;
         rho = 2600.d0;
         lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
         mu = rho * VS**2.d0;
         qs = 0.1*VS;
         gamma = 4.d0*atan(1.d0)/qs;             
     !right
     elseif(dabs(xs - 740948.d0) .le. 2000.d0) then 
         VS = 3490.d0;
         VP = 5770.d0;
         rho = 2600.d0;
         lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
         mu = rho * VS**2.d0;
         qs = 0.1*VS;
         gamma = 4.d0*atan(1.d0)/qs;             
     !up
     elseif(dabs(ys - 4602206.d0) .le. 2000.d0) then 
         VS = 3490.d0;
         VP = 5770.d0;
         rho = 2600.d0;
         lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
         mu = rho * VS**2.d0;
         qs = 0.1*VS;
         gamma = 4.d0*atan(1.d0)/qs;             
     !down
     elseif(dabs(ys - 4502679.d0) .le. 2000.d0) then 
         VS = 3490.d0;
         VP = 5770.d0;
         rho = 2600.d0;
         lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
         mu = rho * VS**2.d0;
         qs = 0.1*VS;
         gamma = 4.d0*atan(1.d0)/qs;             
     endif      
         
              
              
     end subroutine MAKE_MECH_PROP_CASE_016                 
