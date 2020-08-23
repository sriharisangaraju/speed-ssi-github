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


     subroutine   MAKE_NH_Enhanced()

     	use speed_par

     	implicit none
 
      	include 'SPEED.MPI'

      	integer*4 :: count
      	integer*4, dimension(:), allocatable :: node_nhe_flag, NN_src_ind_loc



      	if (mpi_id.eq.0) write(*,'(A)')
        if (mpi_id.eq.0) write(*,'(A)')'---------------Setup Not-Honoring Enhanced ---------------' 

        mpi_comm = SPEED_COMM


        allocate(node_nhe_flag(nnod_loc))
        call MAKE_NH_Enhanced_initialise(nnod_loc, nmat_nhe, val_nhe, &
                               nmat, tag_mat, nelem_loc, con_nnz_loc, con_spx_loc, &
                               xx_spx_loc, yy_spx_loc, zz_spx_loc, &
                               count, &
                               node_nhe_flag, mpi_id, mpi_comm, mpi_file)

        ! count, node_nhe_flag


        allocate(NN_src_ind_loc(count))
        call MAKE_NH_Enhanced_NNSearch(nnod_loc, count, mpi_id, mpi_comm, &
                                      mpi_file, NN_src_ind_loc)


        ! NN_src_ind_loc


        allocate(lambda_nhe(nnod_loc), mu_nhe(nnod_loc), rho_nhe(nnod_loc))
        allocate(Qs_nhe_el(nelem_loc), Qp_nhe_el(nelem_loc)) ! gamma_nhe_el(nelem_loc))

        call MAKE_NH_Enhanced_ASSIGN_PROP(nnod_loc, nmat, prop_mat, sdeg_mat, &
        					   nelem_loc, con_nnz_loc, con_spx_loc, &
                               node_nhe_flag, count, NN_src_ind_loc, QS, QP
                               lambda_nhe, mu_nhe, rho_nhe, Qs_nhe_el, Qp_nhe_el, &
                               mpi_id, mpi_comm)


        deallocate(node_nhe_flag, NN_src_ind_loc)


      	if (mpi_id.eq.0) write(*,'(A)')
        if (mpi_id.eq.0) write(*,'(A)')'--------------- Completed ---------------' 

     end subroutine MAKE_NH_Enhanced