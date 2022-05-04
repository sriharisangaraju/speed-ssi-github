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


!> @brief ...Writing VTK file to visualise in Paraview
!! @author Srihari Sangaraju
!> @date July, 2021 
!> @version 1.0

!> @param[in] loc_n_num. Global node number of 'i'th local node is loc_n_num(i)
!> @param[in] nn_loc. No. of nodes in Local/Current Partition
!> @param[in] nmat_nhe No. of Blocks specified with NHE case
!> @param[in] nhe_mat Tag/Labels of Blocks where NHE has to be implemented
!> @param[in] xs_loc x-coordinate of spectral nodes
!> @param[in] ys_loc y-coordinate of spectral nodes
!> @param[in] zs_loc z-coordinate of spectral nodes

    subroutine   WRITE_VTK_VOLUME_TIMESERIES(its, mpi_id, nn_loc, loc_n_num, nm, sdeg,&
                                            xx_loc, yy_loc, zz_loc, cs_nnz_loc, cs_loc, U2, stress)

      use speed_exit_codes

      implicit none

      include 'SPEED.MPI'

      integer*4 :: nm, nel_loc, nn_loc, cs_nnz_loc, mpi_id, its
      integer*4 :: ie, inode, nn, im
      integer*4 :: unit_mpi
      integer*4 :: ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8
    
      integer*4, dimension(nm) :: sdeg
      integer*4, dimension(nn_loc) :: loc_n_num
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc

      real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
      real*8, dimension(nn_loc*3) :: U2
      real*8, dimension(nn_loc*6) :: stress

      character*70 :: file_nhe_proc, file_nhe_new, mpi_file

      nel_loc = cs_loc(0) - 1

      mpi_file = 'MONITORS'
      unit_mpi = 2500 + mpi_id                   
      write(file_nhe_proc,'(A10,I5.5,A1,I8.8,A4)') 'DIS_X_PROC', mpi_id, '_', its, '.vtk'              

      if(len_trim(mpi_file) .ne. 70) then                                         
         file_nhe_new = mpi_file(1:len_trim(mpi_file)) // '/' // file_nhe_proc
      else 
         file_nhe_new = file_nhe_proc     
      endif
     
      !----------------------------------------------------------------------
      open(unit_mpi,file=file_nhe_new,status='replace')
      write(unit_mpi,'(a)') '# vtk DataFile Version 3.1'
      write(unit_mpi,'(a)') 'material model VTK file'
      write(unit_mpi,'(a)') 'ASCII'
      write(unit_mpi,'(a)') 'DATASET UNSTRUCTURED_GRID'
      write(unit_mpi, '(a,i12,a)') 'POINTS ', nn_loc, ' float'

      ! Node Coordinates
      do inode=1,nn_loc
          write(unit_mpi,'(3e20.12)') xx_loc(inode),yy_loc(inode),zz_loc(inode)
      enddo
      write(unit_mpi,*) ''

      ! Connectivity (note: node indices for vtk start at 0)
      write(unit_mpi,'(a,i12,i12)') "CELLS ",nel_loc,nel_loc*9
      do ie=1,nel_loc
        im = cs_loc(cs_loc(ie -1) + 0 )
        nn = sdeg(im) +1

        ic1 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(1 -1) + 1) - 1
        ic2 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(1 -1) + nn) - 1
        ic3 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(nn -1) + nn) - 1
        ic4 = cs_loc(cs_loc(ie -1) + nn*nn*(1 -1) + nn*(nn -1) + 1) - 1
        ic5 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(1 -1) + 1) - 1
        ic6 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(1 -1) + nn) - 1
        ic7 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(nn -1) + nn) - 1
        ic8 = cs_loc(cs_loc(ie -1) + nn*nn*(nn -1) + nn*(nn -1) + 1) - 1

        write(unit_mpi,'(9i12)') 8, ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8 
      enddo
      write(unit_mpi,*) ''

      ! vtkCellType: hexahedrons (ID = 12)
      write(unit_mpi,'(a,i12)') "CELL_TYPES ",nel_loc
      write(unit_mpi,'(6i12)') (12,ie=1,nel_loc)
      write(unit_mpi,*) ''

      ! Writing Scalar data
      write(unit_mpi,'(a,i12)') "POINT_DATA ",nn_loc

      ! UX
      write(unit_mpi,'(a)') "SCALARS stress_xy float"
      write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      do inode = 1,nn_loc
          !write(unit_mpi,*) U2(3*(inode-1)+1)
          write(unit_mpi,*) stress(6*(inode -1) +4)
      enddo
      write(unit_mpi,*) ''


      close(unit_mpi)   
      !------------------------------------------------------------------

    !   if (mpi_id.eq.0) write(*,'(A)')'Completed.' 

    end subroutine WRITE_VTK_VOLUME_TIMESERIES