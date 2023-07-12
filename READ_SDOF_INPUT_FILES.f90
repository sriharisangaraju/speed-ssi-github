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

  use SPEED_SCI
  use speed_timeloop
  use speed_par, only : filename, sdof_file

  implicit none
  integer*4 :: unit_file

  bldinfo_fp=701+mpi_id

  n_bld = 0

  if ((SDOFnum.gt.0).and.(mpi_id.eq.0)) then
    
    ! Reading Config File - Damping Type, Unit Mass of building
    INQUIRE(FILE="Config.txt", EXIST=isConfigPresent)
    if (isConfigPresent) then
      open(unit_file,file='Config.txt')
      read(unit_file,*)
      read(unit_file,*) configtmp,MasspArea
      read(unit_file,*) kclose
      close(unit_file)
      if(kclose.ne.0) kclose=1.d0
    endif

    ! Reading BLDINFO.txt file
    filename="BLDINFO.txt"

    if(len_trim(sdof_file) .ne. 70) then
      BLDinfo = sdof_file(1:len_trim(sdof_file)) // '/' // filename
    else
      BLDinfo = filename
    endif

    open(bldinfo_fp,file=BLDinfo)
    read(bldinfo_fp,*) n_bld
    close(bldinfo_fp)

    allocate(sys(n_bld))  !!! SDOF system
    allocate(SDOFag(n_bld,3),SDOFgd(n_bld,3))
    allocate(ndt2_sys(n_bld))
    SDOFag = 0; SDOFgd = 0; ndt2_sys = 0;
    call MAKE_SDOF_SYSTEM(BLDinfo, deltat, mpi_id)

  endif

  return
end subroutine READ_SDOF_INPUT_FILES
