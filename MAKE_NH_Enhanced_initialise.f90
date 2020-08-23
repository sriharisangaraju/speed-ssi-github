
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

!> @param[out] count no. of Spectral nodes where Nearest Neighbor search has to be performed
!> @param[out] xs_loc_nhe x-coordinate of spectral nodes corresponding to NHE
!> @param[out] ys_loc_nhe
!> @param[out] ys_loc_nhe 
!> @param[out] node_nhe_flag contains material tags of blocks corresponding to each spectral node


     subroutine   MAKE_NH_Enhanced_initialise(nn_loc, nmat_nhe, val_nhe, &
                               nmat, tag_mat, ne_loc, cs_nnz_loc, cs_loc, &
                               xs_loc, ys_loc, zs_loc, &
                               count, &
                               node_nhe_flag, mpi_id, mpi_comm, mpi_file)

     implicit none
    
     use kdtree2_module

     integer*4 :: nn_loc, mpi_id ,nmat_rhe
     integer*4, dimension(nmat_rhe) :: val_case
     integer*4, dimension(nmat) :: tag_mat
     

     integer*4 :: ne_loc, cs_nnz_loc
     integer*4, dimension(0:cs_nnz_loc) :: cs_loc
     
     real*8, dimension(nn_loc) :: xs_loc, ys_loc, zs_loc
     integer*4, intent(inout) :: count
     integer*4, dimension(nn_loc), intent(inout) :: node_nhe_flag
     real(kdkind), dimension(:), allocatable :: xs_loc_nhe, ys_loc_nhe, zs_loc_nhe

     real*8 :: t0, t1, time_elapsed
     integer*4 :: i, j, ipt, inode, ie
     integer*4 :: im, istart, iend, mpi_ierr, mpi_comm
     
     character*70 :: file_nhe_proc, mpi_file
     

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Making GLL node list which falls inside NHE-specified Blocks
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     node_nhe_flag = 0

     do ie = 1,ne_loc
        im = cs_loc(cs_loc(ie -1))
        istart = cs_loc(ie-1) + 1
        iend = cs_loc(ie) - 1

        do j = istart,iend
          do i = 1,nmat_nhe

            if (val_nhe(i) .eq. tag_mat(im)) then
                ! Node Numbers of nodes where NHE has to be implemented
                !inode = cs_loc(j)  ! Local Node Number
               node_nhe_flag(cs_loc(j)) = 999
            else
                if (node_nhe_flag(cs_loc(j)).ne.999) node_nhe_flag(cs_loc(j)) = im
            endif

          enddo ! i = 1,nmat_nhe
        enddo 
      enddo !ie=1,ne_loc

      
      ! Making List of Local Nodes and their coordinates where
      ! Nearest Neighbor search has to be performed
      count = 0
      do inode =1,nn_loc
          if ((node_nhe_flag(inode).eq.999)) count = count + 1
      enddo

      if (count.gt.0) then
        allocate(xs_loc_nhe(count),ys_loc_nhe(count),zs_loc_nhe(count),loc_nod_num_nhe(count))
      endif

      count = 0
      do inode =1,nn_loc
          if ((node_nhe_flag(inode).eq.999))
              count = count + 1
              xs_loc_nhe(count) = xs_loc(inode)
              ys_loc_nhe(count) = ys_loc(inode)
              zs_loc_nhe(count) = zs_loc(inode)
          endif
      enddo

      call MPI_BARRIER(mpi_comm, mpi_ierr)

      file_nhe_proc = 'nhexyz000000.mpi'
      unit_mpi = 1500 + mpi_id                                 
      if (mpi_id.lt. 10) then                                        
          write(file_nhe_proc(12:12),'(i1)') mpi_id                
      else if (mpi_id .lt. 100) then                                
          write(file_nhe_proc(11:12),'(i2)') mpi_id                
      else if (mpi_id .lt. 1000) then                                
          write(file_nhe_proc(10:12),'(i3)') mpi_id                
      else if (mpi_id .lt. 10000) then                                
          write(file_nhe_proc(9:12),'(i4)') mpi_id                
      else if (mpi_id .lt. 100000) then        
          write(file_nhe_proc(8:12),'(i5)') mpi_id                
      else if (mpi_id .lt. 1000000) then                                
          write(file_nhe_proc(7:12),'(i6)') mpi_id                
      endif
        
      if(len_trim(mpi_file) .ne. 70) then                                         
         file_nhe_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_nhe_proc
      else 
         file_nhe_new = file_nhe_proc     
      endif
     
      open(unit_mpi,file=file_nhe_new)
      write(unit_mpi,*) count                
      if (count.gt.0) then
        do i = 1,count
                 write(unit_mpi,*) xs_loc_nhe(i), ys_loc_nhe(i), zs_loc_nhe(i)
        enddo
        deallocate(xs_loc_nhe, ys_loc_nhe, zs_loc_nhe)
      endif
      close(unit_mpi)        

      call MPI_BARRIER(mpi_comm, mpi_ierr)

     end subroutine MAKE_NH_Enhanced_initialise