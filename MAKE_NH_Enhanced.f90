
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


!> @brief ...Not-Honoring Enhanced (NHE) Implementation
!! @author Srihari Sangaraju
!> @date August, 2020 
!> @version 1.0

!> @ Detailed Explanation
!> 1. User Can specify that a single block/multiple blocks in mesh 
!>    has to use material properties as per provided Tomography text file
!> 2. For the each GLL node in NHE-specified blocks, this subroutine  assigns 
!>    Mechanical properties of the corresponding nearest point it can find 
!>    in "Tomography text file"
!> 3. For the nodes in other blocks It assigns Mechanical properties
!>    specified in *.mate file
!> Note: In this version, This subroutine will be run, only if atleast 1 block is 
!>       specified to use NHE in *.mate file

!> @param[in] loc_n_num. Global node number of 'i'th local node is loc_n_num(i)
!> @param[in] nn_loc. No. of nodes in Local/Current Partition
!> @param[in] nmat_nhe No. of Blocks specified with NHE case
!> @param[in] nhe_mat Tag/Labels of Blocks where NHE has to be implemented
!> @param[in] xs_loc x-coordinate of spectral nodes
!> @param[in] ys_loc y-coordinate of spectral nodes
!> @param[in] zs_loc z-coordinate of spectral nodes

!> @param[out] lambda_nhe Lame coefficient lambda
!> @param[out] mu_nhe Lame coefficient mu
!> @param[out] rho_nhe material density 
!> @param[out] Qs_nhe quality factor for S-waves (Standard Linear Solid damping)
!> @param[out] Qp_nhe quality factor for P-waves (Standard Linear Solid damping)
!> @param[out] gamma_dmp_nhe damping coefficient gamma (Kosloff&Kosloff)


     subroutine   MAKE_NH_Enhanced(loc_n_num, nn_loc, nmat_nhe, nhe_mat, &
                               xs_loc, ys_loc, zs_loc, lambda_nhe, mu_nhe,&
                                rho_nhe, Qs_nhe, Qp_nhe, gamma_dmp_nhe, &
                                mpi_id)

     use kdtree2_module

     use speed_par

     implicit none

     include 'SPEED.MPI'
    
     integer*4 :: nn_loc, mpi_id ,nmat_rhe
     integer*4, dimension(nmat_rhe) :: nhe_mat
     integer*4, dimension(nn_loc) :: loc_n_num, nn_ind

     real*8, dimension(nn_loc) :: xs_loc, ys_loc, zs_loc 
     real*8, dimension(nn_loc), intent(inout) :: lambda_nhe, mu_nhe, rho_nhe
     real*8, dimension(nn_loc), intent(inout) :: Qs_nhe, mu_nhe, gamma_dmp_nhe
   
     character*70 :: file_tomo
     integer*4 :: npts_tomo, stat, error, n_neighbours
     real(kdkind) :: query_vec(3)
     real(kdkind), dimension(:,:), allocatable :: nodes_in_xyz, mech_prop

     type(kdtree2), pointer :: kd2_obj
     type(kdtree2_result) :: result_temp(1)

     real*8 :: t0, t1, time_elapsed
     integer*4 :: i, j, ipt, inode
     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Making GLL node list which falls inside NHE-specified Blocks
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !do ie = 1,nelem_loc
        !check material block. if is NHE block-go ahead, otherwise exit
        !loop over all nodes. prepare a list. no repetitions
        ! map them back to index in loc_n_num


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Reading Tomography Grid Points data
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! tomo_xyz_mech.in = text file with tomography data
     ! 1st Line = No. of Tomo Points (npts_tomo)
     ! 2nd to (npts_tomo+1)th line contains 8 columns:
     ! x_cord y_cord z_cord Rho Vs Vp Qs Qp    (Units must be same as given in *.mate file)
     file_tomo = 'tomo_xyz_mech.in'

     if(mpi_id .eq.0) then
      open(124,file=file_tomo)
      read(124,*) 
      read(124,*) npts_tomo
     endif
     call MPI_BCAST(npts_tomo,1,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)

     allocate(nodes_in_xyz(3,npts_tomo),stat=error)
     allocate(mech_prop(5,npts_tomo),stat=error)

     if(mpi_id.eq.0) then
      do ipt = 1, npts_in
        read(24,*)(nodes_in_xyz(i,ipt), i=1,3), (mech_prop(j,ipt), j=1,5)
      enddo
      close(24)
     endif

     call MPI_BCAST(nodes_in_xyz,,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
     call MPI_BCAST(mech_prop,,SPEED_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)

     call MPI_BARRIER(mpi_comm, mpi_ierr)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Defining Kd-tree Object and then search
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     kd2_obj => kdtree2_create(nodes_in_xyz,sort=.false.,rearrange=.true.)

     deallocate(nodes_in_xyz)

     n_neighbours = 1

     do inode = 1,nn_loc
        query_vec(1) = xs_loc(inode)
        query_vec(2) = ys_loc(inode)
        query_vec(3) = zs_loc(inode)
        call kdtree2_n_nearest(kd2_obj,query_vec,n_neighbours,result_temp)
        nn_ind(inode) = result_temp(1)%idx    ! Inder of Nearest Neighbor in Tomo Points
     enddo

     call kdtree2_destroy(kd2_obj)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Assigning Properties to Each node based on Their Nearest Neighbor
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     

     end subroutine MAKE_RANDOM_PARAM
     