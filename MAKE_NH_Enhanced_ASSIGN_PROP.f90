
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


     subroutine   MAKE_NH_Enhanced_ASSIGN_PROP(nn_loc, nmat, prop_mat, sdeg_mat, &
                               ne_loc, cs_nnz_loc, cs_loc, &
                               node_nhe_flag, count, NN_src_ind_loc, QS, QP, &
                               lambda_nhe, mu_nhe, rho_nhe, &
                               Qs_nhe_el, Qp_nhe_el, mpi_id, mpi_comm) !gamma_dmp_nhe_el

     !use speed_par

     implicit none

     include 'SPEED.MPI'
    
     integer*4 :: nn_loc, mpi_id, nmat, count, mpi_comm, mpi_ierr
     integer*4, dimension(nmat) :: QS, QP, sdeg_mat
     integer*4, dimension(nn_loc) :: node_nhe_flag
     integer*4, dimension(count) ::  NN_src_ind_loc
     
     real*8, dimension(nn_loc), intent(inout) :: lambda_nhe, mu_nhe, rho_nhe
     real*4, dimension(nn_loc) :: Qs_nhe, Qp_nhe !gamma_dmp_nhe

     integer*4 :: ne_loc, cs_nnz_loc
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     real*4, dimension(ne_loc) :: Qs_nhe_el, Qp_nhe_el !gamma_dmp_nhe_el

     character*70 :: file_tomo
     integer*4 :: npts_tomo, stat, error
     real*4, dimension(:), allocatable :: tomo_rho, tomo_vs, tomo_vp, tomo_qs, tomo_qp
     real*8, dimension(nmat,4) :: prop_mat

     real*8 :: t0, t1, time_elapsed, dummy, vs_dum, vp_dum
     integer*4 :: i, j, ipt, inode, ie
     integer*4 :: im, istart, iend, nn


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
     call MPI_BCAST(npts_tomo,1,SPEED_INTEGER,0,mpi_comm,mpi_ierr)

     if(mpi_id.eq.0) then
        allocate(tomo_rho(npts_tomo),tomo_vs(npts_tomo),tomo_vp(npts_tomo),tomo_qs(npts_tomo),tomo_qp(npts_tomo))
        do ipt = 1, npts_in
          read(124,*)dummy, dummy, dummy, tomo_rho, tomo_vs, tomo_vp, tomo_qs, tomo_qp
        enddo
        close(124)
     endif


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Assigning Properties to Each node based on Their Nearest Neighbor
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!! Rho !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mpi_id.ne.0) then
        allocate(tomo_rho(npts_tomo))
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(tomo_rho, npts_tomo, SPEED_INTEGER, 0, mpi_comm, mpi_ierr)

      i = 0
      do inode=1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
            i = i + 1
            rho_nhe(inode) = tomo_rho(NN_src_ind_loc(i))
          else
            rho_nhe(inode) = prop_mat(node_nhe_flag(inode),1)
          endif
      enddo
      deallocate(tomo_rho)


      !!!!!!!!!!!!!!!!!!!!!!!!!!  Mu !!!!!!!!!!!!!!!!!!!!!!!
      if (mpi_id.ne.0) then
        allocate(tomo_vs(npts_tomo))
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(tomo_vs, npts_tomo, SPEED_INTEGER, 0, mpi_comm, mpi_ierr)

      i = 0
      do inode=1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
            i = i + 1
            vs_dum = tomo_vs(NN_src_ind_loc(i))
            mu_nhe(inode) = rho_nhe(inode) * vs_dum**2
          else
            mu_nhe(inode) = prop_mat(node_nhe_flag(inode),3)
          endif
      enddo
      deallocate(tomo_vs)

      !!!!!!!!!!!!!!!!!!!!!!!!!!  Lambda !!!!!!!!!!!!!!!!!!!!!!!
      if (mpi_id.ne.0) then
        allocate(tomo_vp(npts_tomo))
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(tomo_vp, npts_tomo, SPEED_INTEGER, 0, mpi_comm, mpi_ierr)

      i = 0
      do inode=1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
            i = i + 1
            vp_dum = tomo_vp(NN_src_ind_loc(i))
            vs_dum = mu_nhe(inode)/rho_nhe(inode)
            lambda_nhe(inode) = rho_nhe(inode) * (vp_dum**2 - 2*vs_dum) 
          else
            lambda_nhe(inode) = prop_mat(node_nhe_flag(inode),2)
          endif
      enddo
      deallocate(tomo_vp)


      !!!!!!!!!!!!!!!!!!!!!!!!!! Qs !!!!!!!!!!!!!!!!!!!!!!!
      if (mpi_id.ne.0) then
        allocate(tomo_qs(npts_tomo))
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(tomo_qs, npts_tomo, SPEED_INTEGER, 0, mpi_comm, mpi_ierr)

      i = 0
      do inode=1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
            i = i + 1
            Qs_nhe(inode) = tomo_qs(NN_src_ind_loc(i))
          else
            Qs_nhe(inode) = QS(node_nhe_flag(inode),1)
          endif
      enddo
      deallocate(tomo_qs)

      !!!!!!!!!!!!!!!!!!!!!!!!!! Qp !!!!!!!!!!!!!!!!!!!!!!!
      if (mpi_id.ne.0) then
        allocate(tomo_qp(npts_tomo))
      endif

      call MPI_BARRIER(mpi_comm, mpi_ierr)
      call MPI_BCAST(tomo_qp, npts_tomo, SPEED_INTEGER, 0, mpi_comm, mpi_ierr)

      i = 0
      do inode=1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
            i = i + 1
            Qp_nhe(inode) = tomo_qp(NN_src_ind_loc(i))
          else
            Qp_nhe(inode) = QP(node_nhe_flag(inode),1)
          endif
      enddo
      deallocate(tomo_qp)
     
     !!!!!!!!!!!!!!!!!!!!!!!!!! Gamma !!!!!!!!!!!!!!!!!!!!!!!!
     ! Only For Damping Type - 1 (Elastic case)
     ! For other Damping Types Its not needed or = 0

     !do im = 1, nmat
     !       if(QS(im) .eq. 0.d0) then 
     !          prop_mat(im,4) = 0.d0;
     !       else
     !          prop_mat(im,4) = 4.d0*datan(1.d0)*(fmax)/QS(im)
     !       endif   
     !enddo

     !! IF Gamma > 10^(-5) make_damp_yes_or_no = true.
     !! If Gamma < 10^(-5) damping is not assumed

     call MPI_BARRIER(mpi_comm, mpi_ierr)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!! Making Qs, Qp for Element instead of Each node

     Qs_nhe_el = 0;
     Qp_nhe_el = 0;

     do ie = 1,ne_loc
        im = cs_loc(cs_loc(ie -1))
        nn = sdeg_mat(im) + 1
        istart = cs_loc(ie-1) + 1
        iend = cs_loc(ie) - 1

        do j = istart,iend
            Qs_nhe_el(ie) = Qs_nhe_el(ie) + Qs_nhe(cs_loc(j))
            Qp_nhe_el(ie) = Qp_nhe_el(ie) + Qp_nhe(cs_loc(j))
        enddo
        Qs_nhe_el(ie) = Qs_nhe_el(ie)/(nn**3)
        Qs_nhe_el(ie) = Qs_nhe_el(ie)/(nn**3)
      enddo !ie=1,ne_loc

     end subroutine MAKE_NH_Enhanced_ASSIGN_PROP