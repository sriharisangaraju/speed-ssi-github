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

!> @brief Reads inputs files with oscillator properties
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0

subroutine READ_SDOF_INPUT_FILES

  use SDOF_SYSTEM
  use speed_timeloop
  use speed_par, only : filename, sdof_file

  implicit none

  SDOFmon=701+mpi_id

  n_sdof = 0

  if ((SDOFnum.gt.0).and.(mpi_id.eq.0)) then

    filename="SDOFINFO.txt"

    if(len_trim(sdof_file) .ne. 70) then
      SDOFinfo = sdof_file(1:len_trim(sdof_file)) // '/' // filename
    else
      SDOFinfo = filename
    endif

    open(SDOFmon,file=SDOFinfo)
    read(SDOFmon,*) n_sdof
    close(SDOFmon)

    allocate(sys(n_sdof))			!!! SDOF system
    allocate(SDOFag(n_sdof,3),SDOFgd(n_sdof,3))
    SDOFag = 0; SDOFgd = 0
    call MAKE_SDOF_SYSTEM(SDOFinfo, deltat, mpi_id)

  endif

  return
end subroutine READ_SDOF_INPUT_FILES
