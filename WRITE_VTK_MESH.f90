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

    subroutine   WRITE_VTK_MESH(nn_loc, loc_n_num, nm, tag_mat, prop_mat, &
                                  sdeg, xx_loc, yy_loc, zz_loc, cs_nnz_loc, &
                                  cs_loc, mpi_id)

      use speed_exit_codes

      implicit none

      include 'SPEED.MPI'

      integer*4 :: nm, nel_loc, nn_loc, cs_nnz_loc, mpi_id
      integer*4 :: ie, inode, nn, im
      integer*4 :: unit_mpi
      integer*4 :: ic1, ic2, ic3, ic4, ic5, ic6, ic7, ic8

      real*8 :: lambda, mu, qs, qp, dum, gamma

      integer*4, dimension(nm) :: sdeg, tag_mat
      integer*4, dimension(nn_loc) :: loc_n_num, mat_id
      integer*4, dimension(0:cs_nnz_loc) :: cs_loc

      real*8, dimension(nm,4) :: prop_mat
      real*8, dimension(nn_loc) :: xx_loc, yy_loc, zz_loc
      real*8, dimension(nn_loc) :: rho, VS, VP, poisn

      character*70 :: file_nhe_proc, file_nhe_new, mpi_file

      rho = 0.d0
      lambda = 0.d0
      mu = 0.d0
      qs = 0.d0
      qp = 0.d0

      VS = 0.d0
      VP = 0.d0

      nel_loc = cs_loc(0) - 1
      !Tagging Each node to blovkID in *.mate file
      do ie=1,nel_loc
        im = cs_loc(cs_loc(ie -1) + 0 )
        mat_id(cs_loc(cs_loc(ie-1) + 1 : cs_loc(ie)-1)) = im
      enddo

      do inode = 1,nn_loc

        rho(inode) = prop_mat(mat_id(inode),1);
        lambda = prop_mat(mat_id(inode),2);
        mu = prop_mat(mat_id(inode),3);
        gamma = prop_mat(mat_id(inode),4)

        VS(inode) = SQRT(mu/rho(inode))
        VP(inode) = SQRT((lambda/rho(inode)) + (2.d0*VS(inode)**2.d0))
        dum = (VP(inode)/VS(inode))**2.d0
        poisn(inode) = 0.5d0*(dum-2)/(dum-1)
      enddo


      if (mpi_id.eq.0) write(*,'(A)')
      if (mpi_id.eq.0) write(*,'(A)')'------Writing VTK file - Mechanical Properties----------' 

      mpi_file = 'MONITORS'
      unit_mpi = 2500 + mpi_id                   
      ! write(file_nhe_proc,'(A6,I6.6,A4)') 'NHMECH', mpi_id, '.vtk'      
      write(file_nhe_proc,'(A10,I5.5,A4)') 'DIS_X_PROC', mpi_id, '.vtk'        

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

      ! ! Writing Scalar data
      ! write(unit_mpi,'(a,i12)') "POINT_DATA ",nn_loc

      ! ! RHO
      ! write(unit_mpi,'(a)') "SCALARS Rho float"
      ! write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      ! do inode = 1,nn_loc
      !     write(unit_mpi,*) rho(inode)
      ! enddo
      ! write(unit_mpi,*) ''

      ! ! VS
      ! write(unit_mpi,'(a)') "SCALARS VS float"
      ! write(unit_mpi,'(a)') "LOOKUP_TABLE default"
      ! do inode = 1,nn_loc
      !     write(unit_mpi,*) VS(inode)
      ! enddo
      ! write(unit_mpi,*) ''

    !   ! VP
    !   ! write(unit_mpi,'(a)') "SCALARS VP float"
    !   write(unit_mpi,'(a)') "SCALARS Thick float"
    !   write(unit_mpi,'(a)') "LOOKUP_TABLE default"
    !   do inode = 1,nn_loc
    !       ! write(unit_mpi,*) VP(inode)
    !       write(unit_mpi,*) thick_nodes(inode)
    !   enddo
    !   write(unit_mpi,*) ''

    !   ! Poissons Ratio
    !   ! write(unit_mpi,'(a)') "SCALARS Poisson float"
    !   write(unit_mpi,'(a)') "SCALARS depth float"
    !   write(unit_mpi,'(a)') "LOOKUP_TABLE default"
    !   do inode = 1,nn_loc
    !       ! write(unit_mpi,*) poisn(inode)
    !       write(unit_mpi,*) zs_elev(inode)
    !   enddo
    !   write(unit_mpi,*) ''

      close(unit_mpi)   
      !------------------------------------------------------------------

      if (mpi_id.eq.0) write(*,'(A)')'Completed.' 

    end subroutine WRITE_VTK_MESH