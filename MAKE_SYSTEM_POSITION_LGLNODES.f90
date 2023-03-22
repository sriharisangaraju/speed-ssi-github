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

!> @brief Reads oscillator position from SYS.input, and
!>  finds to which local-LGL node it should be applied.
!! @author Ilario Mazzieri, Srihari Sangaraju
!> @date February, 2022
!> Outputs: (directly saved in module variables of speed_par)
!> locnode_buildID_map  - For each local LGL node, has information about SDOFs that have to be applied at that location.
!> node_counter_sdof    - Number of partitions that are sharing LGL node corresponding to each SDOF Oscillator.

subroutine MAKE_SYSTEM_POSITION_LGLNODES(nn_loc , local_node_num, cs_nnz_loc, cs_loc, &
                                          xs_loc, ys_loc, zs_loc, nsdof, sys_label, &
                                          x_system_lst, y_system_lst, z_system_lst, &
                                          mpi_id, mpi_comm, mpi_np)

  use speed_par, only : locnode_buildID_map, node_counter_sdof

  implicit none
  
  include 'SPEED.MPI'

  integer*4 :: nn_loc, nsdof, cs_nnz_loc
  integer*4 :: mpi_np, mpi_id, mpi_comm, mpi_ierr
  integer*4 :: i, j, max_nsdof_per_node
  integer*4, dimension(0:cs_nnz_loc) :: cs_loc
  integer*4, dimension(nn_loc) :: node_counter_loc
  integer*4, dimension(nsdof) :: node_counter_sdof_dum
  integer*4, dimension(nsdof) :: node_sys_loc, sys_label, nearestnode_globnum, nearestnode_locnum
  integer*4, dimension(nsdof*mpi_np) :: node_sys_glo
  integer*4, dimension(nn_loc) :: local_node_num

  !integer*4, dimension(:,:), allocatable  :: locnode_buildID_map   !INTENT(OUT)

  real*8, dimension(nn_loc) :: xs_loc, ys_loc, zs_loc
  real*8, dimension(nsdof) :: x_system_lst, y_system_lst, z_system_lst
  real*8, dimension(nsdof) :: dist_system_loc
  real*8, dimension(nsdof*mpi_np) :: dist_system_glo

  node_counter_loc = 0

  ! Find if LGL_gob_num is inside the local partition
  ! Finding Nearest Node to SDOF Oscillator across all partitions - in Global Numeration
  do i = 1,nsdof
    call GET_NEAREST_NODE(nn_loc, xs_loc, ys_loc, zs_loc, x_system_lst(i), y_system_lst(i), z_system_lst(i), &
                          node_sys_loc(i), dist_system_loc(i))
    node_sys_loc(i) = local_node_num(node_sys_loc(i))
  enddo

  call MPI_BARRIER(mpi_comm, mpi_ierr)

  call MPI_ALLGATHER(dist_system_loc, nsdof, SPEED_DOUBLE, dist_system_glo, nsdof, SPEED_DOUBLE, mpi_comm, mpi_ierr)
  call MPI_ALLGATHER(node_sys_loc, nsdof, SPEED_INTEGER, node_sys_glo, nsdof, SPEED_INTEGER, mpi_comm, mpi_ierr)

  node_sys_loc = 0

  call GET_MINVALUES(node_sys_glo, dist_system_glo, nsdof*mpi_np, node_sys_loc, nsdof, mpi_np)
  call MPI_BARRIER(mpi_comm, mpi_ierr)
  
  ! Finding if this nearest node is common in multiple partitions
  node_counter_sdof_dum = 0
  do i = 1, nsdof
    nearestnode_globnum(i) = node_sys_glo(node_sys_loc(i))

    call GET_INDLOC_FROM_INDGLO(local_node_num, nn_loc, nearestnode_globnum(i), nearestnode_locnum(i))

    if (nearestnode_locnum(i).ne.0) then
      node_counter_sdof_dum(i) = 1;
      node_counter_loc(nearestnode_locnum(i)) = node_counter_loc(nearestnode_locnum(i)) + 1
    endif
  enddo
  
  ! Finiding LGL nodes that are shared by different partitions
  allocate(node_counter_sdof(nsdof))
  call MPI_BARRIER(mpi_comm, mpi_ierr)
  call MPI_ALLREDUCE(node_counter_sdof_dum, node_counter_sdof, nsdof, SPEED_INTEGER, MPI_SUM, mpi_comm, mpi_ierr)

  max_nsdof_per_node = MAXVAL(node_counter_loc)  ! num on SDOFs allocated to a LGL node, id of 1st SDOF, id of 2nd SDOF etc...
  allocate(locnode_buildID_map(nn_loc, max_nsdof_per_node + 1))
  locnode_buildID_map = 0

  ! Mapping SDOFs to nearest LGL nodes; 
  do i = 1, nsdof
    nearestnode_globnum(i) = node_sys_glo(node_sys_loc(i))

    call GET_INDLOC_FROM_INDGLO(local_node_num, nn_loc, nearestnode_globnum(i), nearestnode_locnum(i))

    if (nearestnode_locnum(i).ne.0) then
      locnode_buildID_map(nearestnode_locnum(i), 1) = locnode_buildID_map(nearestnode_locnum(i), 1) + 1
      locnode_buildID_map(nearestnode_locnum(i), locnode_buildID_map(nearestnode_locnum(i), 1)+1 ) = sys_label(i)
      write(*,*) mpi_id, i, sys_label(i),'; node_counter = ', node_counter_sdof, '; coords = ', xs_loc(nearestnode_locnum(i)), &
      ys_loc(nearestnode_locnum(i)), zs_loc(nearestnode_locnum(i)) 
    endif
  enddo

end subroutine MAKE_SYSTEM_POSITION_LGLNODES
