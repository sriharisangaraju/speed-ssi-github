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



     subroutine MAKE_MECH_PROP_CASE_032(rho, lambda, mu, gamma, qs, qp, & !outputs
                                        xs, ys, zs, Depth, zs_all,&
                                        vs30, thickness, sub_tag_all)
              
     real*8, intent(out) :: rho, lambda, mu, gamma, qs, qp
     real*8, intent(in)  :: xs, ys, zs, Depth, zs_all,&
        vs30, thickness
     integer*4           :: sub_tag_all		   		
     real*8              :: ni, VS, VP, Depth_real         
     real*8, dimension(1) :: val1
                    
     rho    = 0.d0;
     lambda = 0.d0;
     mu     = 0.d0;
     gamma  = 0.d0;
     qs     = 0.d0;
     qp     = 0.d0
     
    !-------------------------------------------------------------------
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 1st Layer
    ! from ground surface to PE_B
    if (sub_tag_all.eq.1 .or. sub_tag_all.eq.2 ) then
            
            if (dabs(zs) .le. 15.d0) then
              VS = 150.d0;
              VP  = 1064.d0;
              rho = 1800.d0;
              if (check_case .eq. 1)  &
                  write(1000+mpi_id,*) xs,ys,zs, &
                   VS, VP, rho   
            elseif (dabs(zs) .le. 30.d0) then
              VS = 294.d0;
              VP  = 1321.d0;
              rho = 1800.d0;
              if (check_case .eq. 1)  &
                  write(1000+mpi_id,*) xs,ys,zs, &
                   VS, VP, rho   
            elseif (dabs(zs) .le. 50.d0) then
              VS = 332.d0;
              VP  = 1494.d0;
              rho = 1800.d0;
              if (check_case .eq. 1)  &
                  write(1000+mpi_id,*) xs,ys,zs, &
                   VS, VP, rho   
            elseif (dabs(zs) .le. 80.d0) then
              VS = 370.d0;
              VP  = 1664.d0;
              rho = 1800.d0;
            elseif (dabs(zs) .le. 120.d0) then
              VS = 408.d0;
              VP  = 1838.d0;
              rho = 1800.d0;
            elseif (dabs(zs) .le. 180.d0) then
              VS = 450.d0;
              VP  = 2024.d0;
              rho = 1800.d0;
            else ! NU_B 
              VS = 523.d0;
              VP  = 1988.d0;
              rho = 2050.d0;
            endif

            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs; 

   ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 3rd Layer
   ! from NU_B to NS_B
    elseif (sub_tag_all.eq.3) then
            VS = 600.d0;
            VP  = VS * 3.2d0;
            rho = 2050.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;       
             

    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 4th Layer
    ! CK
    elseif (sub_tag_all.eq.4) then

            VS = 2090.d0;
            if (dabs(zs) .le. 1200.d0) VS = 1515.d0 
            
            VP  = VS * 2.d0;
            rho = 2400.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 5th Layer
    ! ZE
    elseif (sub_tag_all.eq.5) then
            VS = 2850.d0;
            VP  = 5100.d0;
            rho = 2450.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    ! + MATERIAL INSIDE THE ALLUVIAL BASIN - 6th Layer
    ! RO (reservoir)
    elseif (sub_tag_all.eq.6) then
            VS = 2300.d0;
            VP  = 3900.d0;
            rho = 2450.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;            
    else ! LI (halfspace)
            VS = 2600.d0;
            VP  = 4500.d0;
            rho = 2650.d0;
            lambda = rho * (VP**2.d0 - 2.d0*VS**2.d0);
            mu = rho * VS**2.d0;
            qs = VS/10.d0;
            qp = VP/10.d0
            gamma = 4.d0*datan(1.d0)*5.d0/qs;
    endif

     
     end subroutine MAKE_MECH_PROP_CASE_032           
