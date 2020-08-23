
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


     subroutine   MAKE_NH_Enhanced_NNSearch(nn_loc, count, mpi_id, mpi_np, &
                                            mpi_comm, mpi_file, NN_src_ind_loc)

     use kdtree2_module

     use speed_par

     use speed_exit_codes

     implicit none

     include 'SPEED.MPI'
    
     integer*4 :: nn_loc, mpi_id,  mpi_np, count, mpi_comm, mpi_ierr
     integer*4, dimension(:), allocatable :: NN_src_ind_ip
     integer*4, dimension(count) :: NN_src_ind_loc
     integer*4, dimension(mpi_np) :: count_proc_nhe
   
     character*70 :: file_tomo, mpi_file, file_nhe_nnind, file_nhe_proc
     integer*4 :: npts_tomo, stat, error, n_neighbours
     real(kdkind) :: query_vec(3)
     real(kdkind), dimension(:), allocatable :: xs_ip, ys_ip, zs_ip
     real(kdkind), dimension(:,:), allocatable :: nodes_in_xyz

     type(kdtree2), pointer :: kd2_obj
     type(kdtree2_result) :: result_temp(1)

     real*8 :: t0, t1, time_elapsed
     real*8, dimension(5) :: dummy
     integer*4 :: i, j, ipt, inode, ip, ncount


     call MPI_BARRIER(mpi_comm, mpi_ierr)
     call MPI_ALLGATHER(count, 1, SPEED_INTEGER, count_proc_nhe, 1, SPEED_INTEGER, mpi_comm, mpi_ierr)


     if(mpi_id.eq.0) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Reading Tomography Grid Points data
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! tomo_xyz_mech.in = text file with tomography data
       ! 1st Line = No. of Tomo Points (npts_tomo)
       ! 2nd to (npts_tomo+1)th line contains 8 columns:
       ! x_cord y_cord z_cord Rho Vs Vp Qs Qp    (Units must be same as given in *.mate file)
       file_tomo = 'tomo_xyz_mech.in'

        open(124,file=file_tomo)
        read(124,*) 
        read(124,*) npts_tomo

       allocate(nodes_in_xyz(3,npts_tomo),stat=error)
       if (stat.ne.0) then
          write(*,*)'error: couldnt allocate memory for array,',&
          ' size=',npts_tomo
          call EXIT(EXIT_NORMAL)
       endif

       if(mpi_id.eq.0) then
        do ipt = 1, npts_in
          read(124,*)(nodes_in_xyz(i,ipt), i=1,3), (dummy(j), j=1,5)
        enddo
        close(124)
       endif


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Defining Kd-tree Object and then search
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       kd2_obj => kdtree2_create(nodes_in_xyz,sort=.false.,rearrange=.true.)
       n_neighbours = 1


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Getting NHE Nodes from Each MPI Process and Performing search
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       do ip = 0, mpi_np-1

          ! Reading coordinates of nodes where NN search have to be performed
          file_nhe_proc = 'nhexyz000000.mpi'
          file_nhe_nnind = 'nhenni000000.mpi'
          unit_mpi = 1500 + ip
          if (ip.lt. 10) then
              write(file_nhe_proc(12:12),'(i1)') ip
              write(file_nhe_nnind(12:12),'(i1)') ip
          else if (ip .lt. 100) then
              write(file_nhe_proc(11:12),'(i2)') ip
              write(file_nhe_nnind(11:12),'(i2)') ip
          else if (ip .lt. 1000) then 
              write(file_nhe_proc(10:12),'(i3)') ip
              write(file_nhe_nnind(10:12),'(i3)') ip
          else if (ip .lt. 10000) then
              write(file_nhe_proc(9:12),'(i4)') ip
              write(file_nhe_nnind(9:12),'(i4)') ip
          else if (ip .lt. 100000) then
              write(file_nhe_proc(8:12),'(i5)') ip
              write(file_nhe_nnind(8:12),'(i5)') ip
          else if (ip .lt. 1000000) then 
              write(file_nhe_proc(7:12),'(i6)') ip
              write(file_nhe_nnind(7:12),'(i6)') ip
          endif
            
          if(len_trim(mpi_file) .ne. 70) then                                         
             file_nhe_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_nhe_proc
             file_nhe_nnind_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_nhe_nnind
          else 
             file_nhe_new = file_nhe_proc     
             file_nhe_nnind_new = file_nhe_nnind
          endif
         
          open(unit_mpi,file=file_nhe_new)
          read(unit_mpi,*) ncount                
          if(ncount.ne.count_proc_nhe(ip+1) then
            write(*,*)'ncount in NHE_proc files are not consistent'
            call EXIT(EXIT_NO_NODES)
          endif

          if (ncount.gt.0) then
            allocate(xs_ip(ncount), ys_ip(ncount), zs_ip(ncount))
            allocate(NN_src_ind_ip(ncount))
            do inode = 1,ncount
                read(unit_mpi,*) xs_ip(inode), ys_ip(inode), zs_ip(inode)
            enddo
          endif
          close(unit_mpi)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Performing NN Search for nodes in each Partition
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (ncount.gt.0) then
            do inode = 1,ncount
              query_vec(1) = xs_ip(inode)
              query_vec(2) = xs_ip(inode)
              query_vec(3) = xs_ip(inode)
              call kdtree2_n_nearest(kd2_obj,query_vec,n_neighbours,result_temp)
              NN_src_ind_ip(inode) = result_temp(1)%idx    ! Inder of Nearest Neighbor in Tomo Points
           enddo
          endif

          open(unit_mpi,file=file_nhe_nnind_new) 
          if (ncount.gt.0) then
            do i = 1,ncount
                write(unit_mpi,*) NN_src_ind_ip(i)
            enddo
            deallocate(xs_ip, ys_ip, zs_ip, NN_src_ind_ip)
          endif
          close(unit_mpi)

       enddo ! ip = 0, mpi_np-1


       deallocate(nodes_in_xyz)
       call kdtree2_destroy(kd2_obj)


     endif !if mpi_id == 0

     call MPI_BARRIER(mpi_comm, mpi_ierr)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Reading Back the results of NN Search (in respective processor now)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      file_nhe_nnind = 'nhenni000000.mpi'
      unit_mpi = 1500 + ip
      if (ip.lt. 10) then
          write(file_nhe_nnind(12:12),'(i1)') ip
      else if (ip .lt. 100) then
          write(file_nhe_nnind(11:12),'(i2)') ip
      else if (ip .lt. 1000) then 
          write(file_nhe_nnind(10:12),'(i3)') ip
      else if (ip .lt. 10000) then
          write(file_nhe_nnind(9:12),'(i4)') ip
      else if (ip .lt. 100000) then
          write(file_nhe_nnind(8:12),'(i5)') ip
      else if (ip .lt. 1000000) then 
          write(file_nhe_nnind(7:12),'(i6)') ip
      endif
        
      if(len_trim(mpi_file) .ne. 70) then
         file_nhe_nnind_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_nhe_nnind
      else
         file_nhe_nnind_new = file_nhe_nnind
      endif
     
      open(unit_mpi,file=file_nhe_nnind_new)
      if (count.gt.0) then
        do inode = 1,count
            read(unit_mpi,*) NN_src_ind_loc(inode)
        enddo
      endif
      close(unit_mpi)

      call MPI_BARRIER(mpi_comm, mpi_ierr)

     end subroutine MAKE_NH_Enhanced_NNSearch
     