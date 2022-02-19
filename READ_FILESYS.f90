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

!> @brief Reads SYS.input file
!! @author Aline Herlin
!> @date February, 2021
!> @version 1.0
!> @param[in] filec       file name
!> @param[in] num_nodes   number of oscillators
!> @param[out] label      system ID
!> @param[out] x, y, z    coord. of the structure

subroutine READ_FILESYS(filec, num_nodes, label, x, y, z)

  implicit none

  character*70 :: filec
  character*100000 :: input_line
  integer*4 :: i, num_nodes, ileft, iright, status
  real*8, dimension(num_nodes) :: x, y, z
  integer*4, dimension(num_nodes) :: label

  open(20,file=filec)
  read(20,'(A)',IOSTAT = status) input_line
  ileft = 1
  iright = len(input_line)
  read(input_line(ileft:iright),*) num_nodes

  do i = 1,num_nodes
    read(20,'(A)',IOSTAT = status) input_line
    if (status.ne.0) exit
    ileft = 1
    iright = len(input_line)
    read(input_line(ileft:iright),*) label(i), x(i), y(i), z(i)
  enddo

  close(20)

  return

end subroutine READ_FILESYS
