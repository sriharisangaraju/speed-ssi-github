
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


!> @brief ...
!! @author Ilario Mazzieri
!> @date September, 2013 
!> @version 1.0
!> @param[in] ....
!> @param[out] ...

     subroutine   MAKE_RANDOM_PARAM(loc_n_num, nn_loc, &
                               nmat_rnd, rand_mat,  &
                               xs_loc, ys_loc, zs_loc, &
                               lambda_rnd, mu_rnd, rho_rnd, mpi_id)


     implicit none
    
  
     integer*4 :: nn_loc, mpi_id ,nmat_rnd
     integer*4, dimension(nmat_rnd) :: rand_mat
     integer*4, dimension(nn_loc) :: loc_n_num

     real*8, dimension(nn_loc) :: xs_loc, ys_loc, zs_loc 
     real*8, dimension(nn_loc), intent(inout) :: lambda_rnd, mu_rnd, rho_rnd
   
     character*70 :: file_rnd
     integer*4 :: size_values, unit_rnd, i, js, jr
     integer*4, dimension(nn_loc) :: nearest_el
     
     real*8 :: dist_sq, dist_sq_min
     real*8, dimension(nn_loc) :: dist

     real*8, dimension(:,:), allocatable :: values     
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     file_rnd = 'NH-FILES/NHCheck00000.dat'
     
     if (mpi_id .lt. 10) then                                        
           write(file_rnd(21:21),'(i1)') mpi_id;
     else if (mpi_id .lt. 100) then                                
           write(file_rnd(20:21),'(i2)') mpi_id;                       
     else if (mpi_id .lt. 1000) then                                
           write(file_rnd(19:21),'(i3)') mpi_id;                
     else if (mpi_id .lt. 10000) then                                
           write(file_rnd(18:21),'(i4)') mpi_id;                
     else                                                        
           write(file_rnd(17:21),'(i5)') mpi_id                
     endif                                                
    
     unit_rnd = 50 + mpi_id

     open(unit_rnd,file=file_rnd)
     read(unit_rnd,*) size_values

     allocate(values(size_values,6))
     values = 0.d0;
          
     do i = 1, size_values
        read(unit_rnd,*) values(i,1),values(i,2),values(i,3),values(i,4),values(i,5),values(i,6) 
     enddo

     close(unit_rnd)


     do js = 1, nn_loc

        dist_sq_min = 1.d30
        nearest_el(js) = -1

        do jr = 1, size_values

           dist_sq =  (xs_loc(js) - values(jr,1))**2 + (ys_loc(js) - values(jr,2))**2 &
                       + (zs_loc(js) - values(jr,3))**2  

           if ( dist_sq < dist_sq_min ) then
                dist_sq_min = dist_sq
                nearest_el(js) = jr
           endif

        enddo
        
        ! mu = vs^2 * ro
        mu_rnd(js) = values(nearest_el(js),6) * values(nearest_el(js),4)**2 
        ! lambda = vp^2 * ro - 2*mu
        lambda_rnd(js) = values(nearest_el(js),6) * values(nearest_el(js),5)**2 - 2.d0 * mu_rnd(js)
        rho_rnd(js) = values(nearest_el(js),6);
          
     enddo


     deallocate(values)     
     


     end subroutine MAKE_RANDOM_PARAM
     


