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

!> @brief Write system output files
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0
!> @param[in] tt1tmp  current time instant

subroutine WRITE_SDOF_OUTPUT_FILES(tt1tmp)

  use SDOF_SYSTEM
  use speed_timeloop
	implicit none

	integer*4 :: sfs, temp, i
  real*8 :: tt1tmp

  SDOFmon=10*(mpi_id+1)+7

  sfs = 0
  do i = 1, n_sdof
    temp = sys(i)%SFS
    if(temp.eq.1) sfs = temp
  enddo

  if (SDOFout(1).eq.1) then
    ! Relative Displacement
    if (sfs.eq.0) then
      open(SDOFmon,file=SDOFdisplX,position='append')			!!! x displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%SDOFIDR(1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFdisplY,position='append')			!!! y displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%SDOFIDR(2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFdisplZ,position='append')			!!! z displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%SDOFIDR(3)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      !!! soil displacement

      open(SDOFmon,file=SDOFgrdisplX,position='append')			!!! x displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFgrdisplY,position='append')			!!! y displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFgrdisplZ,position='append')			!!! z displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,3)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    elseif(sfs.eq.1) then

      open(SDOFmon,file=STRdisplX,position='append')			!!! structure x displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(1,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=STRdisplY,position='append')			!!! structure y displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(1,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDdisplX,position='append')			!!! soil x displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDdisplY,position='append')			!!! soil y displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDdisplZ,position='append')			!!! soil z displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFgd(i,3)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplX,position='append')			!!! foundation x displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(2,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplY,position='append')			!!! foundation y displacement
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(2,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplRX,position='append')			!!! foundation rotation from x motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(3,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplRY,position='append')			!!! foundation rotation from y motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(3,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplZX,position='append')			!!! foundation z displacement from x motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(4,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDdisplZY,position='append')			!!! foundation z displacement from y motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%u(4,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    endif
	endif

  if (SDOFout(2).eq.1) then

    if(sfs.eq.0) then
      !!! relative acceleration - SS

      open(SDOFmon,file=SDOFaccX,position='append')			!!! x acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%tempSDOFRA1(1)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFaccY,position='append')			!!! y acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%tempSDOFRA1(2)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFaccZ,position='append')			!!! z acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') sys(i)%tempSDOFRA1(3)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      !!! soil acceleration

      open(SDOFmon,file=SDOFgraccX,position='append')			!!! x ground acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,1)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFgraccY,position='append')			!!! y ground acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,2)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFgraccZ,position='append')			!!! y ground acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,3)
      enddo
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    elseif(sfs.eq.1) then

      open(SDOFmon,file=STRaccX,position='append')			!!! structure x acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(1,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=STRaccY,position='append')			!!! structure y acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(1,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDaccX,position='append')			!!! soil x acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDaccY,position='append')			!!! soil y acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=GRDaccZ,position='append')			!!! soil z acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') -SDOFag(i,3)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccX,position='append')			!!! foundation x acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(2,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccY,position='append')			!!! foundation y acceleration
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(2,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccRX,position='append')			!!! foundation rocking acceleration from x motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(3,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccRY,position='append')			!!! foundation rocking acceleration from y motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(3,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccZX,position='append')			!!! foundation z acceleration from x motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(4,1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDaccZY,position='append')			!!! foundation z acceleration from y motion
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%a(4,2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    endif
  endif

  if (SDOFout(3).eq.1) then

    if(sfs.eq.0) then
      open(SDOFmon,file=SDOFfX,position='append')			!!! x interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+1)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFfY,position='append')			!!! y interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+2)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=SDOFfZ,position='append')			!!! z interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+3)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    elseif(sfs.eq.1) then

      open(SDOFmon,file=STRfX,position='append')			!!! structure x force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%fs(1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=STRfY,position='append')			!!! structure y force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%fs(2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDfX,position='append')			!!! foundation x shear force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%fb(1)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=FNDfY,position='append')			!!! foundation y shear force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') sys(i)%fb(2)
      enddo

      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=INTfX,position='append')			!!! x interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+1)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=INTfY,position='append')			!!! y interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+2)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)

      open(SDOFmon,file=INTfZ,position='append')			!!! z interaction force
      write(SDOFmon,"(E16.7)",advance='NO') tt1tmp
      do i=1,n_sdof
        if(sys(i)%SFS.eq.1) write(SDOFmon,"(E16.7)",advance='NO') SDOFforceinput(3*(i-1)+3)
      end do
      write(SDOFmon,"(A1)") " "
      close(SDOFmon)
    endif
  endif

  return
end subroutine WRITE_SDOF_OUTPUT_FILES
