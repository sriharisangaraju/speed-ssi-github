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

!> @brief Creates files to store oscillator motion
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0

subroutine MAKE_SDOF_OUTPUT_FILES

  use SPEED_SCI
  use speed_timeloop
  use speed_par

  implicit none
  character*100 :: file
  integer*4 :: temp, sfs, ii, SDOFmon

  SDOFmon = 701 + mpi_id

  if (n_bld.gt.0) then

    sfs = 0
    do ii = 1, n_bld
      temp = sys(ii)%SFS
      if(temp.eq.1) sfs = temp
    enddo

    if (SDOFout(1).eq.1) then     ! displacement

      if(sfs.eq.0) then

        !!! structure relative displacement

        file="SDOF000000.DX"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFdisplX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFdisplX = file
        endif

        open(SDOFmon,file=SDOFdisplX,status='replace')
        close(SDOFmon)

        file="SDOF000000.DY"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFdisplY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFdisplY = file
        endif

        open(SDOFmon,file=SDOFdisplY,status='replace')
        close(SDOFmon)

        file="SDOF000000.DZ"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFdisplZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFdisplZ = file
        endif

        open(SDOFmon,file=SDOFdisplZ,status='replace')
        close(SDOFmon)

        !!! soil displacement

        file="SDOF000000.GDX"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgrdisplX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgrdisplX = file
        endif

        open(SDOFmon,file=SDOFgrdisplX,status='replace')
        close(SDOFmon)

        file="SDOF000000.GDY"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgrdisplY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgrdisplY = file
        endif

        open(SDOFmon,file=SDOFgrdisplY,status='replace')
        close(SDOFmon)

        file="SDOF000000.GDZ"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgrdisplZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgrdisplZ = file
        endif

        open(SDOFmon,file=SDOFgrdisplZ,status='replace')
        close(SDOFmon)

      elseif(sfs.eq.1) then

        !!! structure relative displacement

        file="STR000000.DX"

        if(len_trim(monitor_file) .ne. 70) then
          STRdisplX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRdisplX = file
        endif

        open(SDOFmon,file=STRdisplX,status='replace')
        close(SDOFmon)

        file="STR000000.DY"

        if(len_trim(monitor_file) .ne. 70) then
          STRdisplY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRdisplY = file
        endif

        open(SDOFmon,file=STRdisplY,status='replace')
        close(SDOFmon)

        !!! soil displacement

        file="GRD000000.DX"

        if(len_trim(monitor_file) .ne. 70) then
          GRDdisplX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDdisplX = file
        endif

        open(SDOFmon,file=GRDdisplX,status='replace')
        close(SDOFmon)

        file="GRD000000.DY"

        if(len_trim(monitor_file) .ne. 70) then
          GRDdisplY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDdisplY = file
        endif

        open(SDOFmon,file=GRDdisplY,status='replace')
        close(SDOFmon)

        file="GRD000000.DZ"

        if(len_trim(monitor_file) .ne. 70) then
          GRDdisplZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDdisplZ = file
        endif

        open(SDOFmon,file=GRDdisplZ,status='replace')
        close(SDOFmon)

        !!! foundation displacement

        file="FND000000.DX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplX = file
        endif

        open(SDOFmon,file=FNDdisplX,status='replace')
        close(SDOFmon)

        file="FND000000.DY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplY = file
        endif

        open(SDOFmon,file=FNDdisplY,status='replace')
        close(SDOFmon)

        file="FND000000.DRX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplRX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplRX = file
        endif

        open(SDOFmon,file=FNDdisplRX,status='replace')
        close(SDOFmon)

        file="FND000000.DRY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplRY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplRY = file
        endif

        open(SDOFmon,file=FNDdisplRY,status='replace')
        close(SDOFmon)

        file="FND000000.DZX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplZX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplZX = file
        endif

        open(SDOFmon,file=FNDdisplZX,status='replace')
        close(SDOFmon)

        file="FND000000.DZY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDdisplZY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDdisplZY = file
        endif

        open(SDOFmon,file=FNDdisplZY,status='replace')
        close(SDOFmon)
      endif
    endif

    if (SDOFout(2).eq.1) then       ! acceleration

      if(sfs.eq.0) then

        !!! structure total acceleration

        file="SDOF000000.AX"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFaccX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFaccX = file
        endif

        open(SDOFmon,file=SDOFaccX,status='replace')
        close(SDOFmon)

        file="SDOF000000.AY"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFaccY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFaccY = file
        endif

        open(SDOFmon,file=SDOFaccY,status='replace')
        close(SDOFmon)

        file="SDOF000000.AZ"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFaccZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFaccZ = file
        endif

        open(SDOFmon,file=SDOFaccZ,status='replace')
        close(SDOFmon)

        !!! soil acceleration

        file="SDOF000000.GAX"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgraccX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgraccX = file
        endif

        open(SDOFmon,file=SDOFgraccX,status='replace')
        close(SDOFmon)

        file="SDOF000000.GAY"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgraccY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgraccY = file
        endif

        open(SDOFmon,file=SDOFgraccY,status='replace')
        close(SDOFmon)

        file="SDOF000000.GAZ"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFgraccZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFgraccZ = file
        endif

        open(SDOFmon,file=SDOFgraccZ,status='replace')
        close(SDOFmon)
      elseif(sfs.eq.1) then

        !!! structure acceleration

        file="STR000000.AX"

        if(len_trim(monitor_file) .ne. 70) then
          STRaccX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRaccX = file
        endif

        open(SDOFmon,file=STRaccX,status='replace')
        close(SDOFmon)

        file="STR000000.AY"

        if(len_trim(monitor_file) .ne. 70) then
          STRaccY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRaccY = file
        endif

        open(SDOFmon,file=STRaccY,status='replace')
        close(SDOFmon)

        !!! soil acceleration

        file="GRD000000.AX"

        if(len_trim(monitor_file) .ne. 70) then
          GRDaccX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDaccX = file
        endif

        open(SDOFmon,file=GRDaccX,status='replace')
        close(SDOFmon)

        file="GRD000000.AY"

        if(len_trim(monitor_file) .ne. 70) then
          GRDaccY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDaccY = file
        endif

        open(SDOFmon,file=GRDaccY,status='replace')
        close(SDOFmon)

        file="GRD000000.AZ"

        if(len_trim(monitor_file) .ne. 70) then
          GRDaccZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          GRDaccZ = file
        endif

        open(SDOFmon,file=GRDaccZ,status='replace')
        close(SDOFmon)

        !!! foundation acceleration

        file="FND000000.AX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccX = file
        endif

        open(SDOFmon,file=FNDaccX,status='replace')
        close(SDOFmon)

        file="FND000000.AY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccY = file
        endif

        open(SDOFmon,file=FNDaccY,status='replace')
        close(SDOFmon)

        file="FND000000.ARX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccRX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccRX = file
        endif

        open(SDOFmon,file=FNDaccRX,status='replace')
        close(SDOFmon)

        file="FND000000.ARY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccRY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccRY = file
        endif

        open(SDOFmon,file=FNDaccRY,status='replace')
        close(SDOFmon)

        file="FND000000.AZX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccZX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccZX = file
        endif

        open(SDOFmon,file=FNDaccZX,status='replace')
        close(SDOFmon)

        file="FND000000.AZY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDaccZY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDaccZY = file
        endif

        open(SDOFmon,file=FNDaccZY,status='replace')
        close(SDOFmon)
      endif
    endif

    if (SDOFout(3).eq.1) then       ! interaction force

      if (sfs.eq.0) then

        !!! structure shear force

        file="SDOF000000.FX"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFfX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFfX = file
        endif

        open(SDOFmon,file=SDOFfX,status='replace')
        close(SDOFmon)

        file="SDOF000000.FY"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFfY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFfY = file
        endif

        open(SDOFmon,file=SDOFfY,status='replace')
        close(SDOFmon)

        file="SDOF000000.FZ"

        if(len_trim(monitor_file) .ne. 70) then
          SDOFfZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          SDOFfZ = file
        endif

        open(SDOFmon,file=SDOFfZ,status='replace')
        close(SDOFmon)
      elseif(sfs.eq.1) then

        !!! structure shear force

        file="STR000000.FX"

        if(len_trim(monitor_file) .ne. 70) then
          STRfX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRfX = file
        endif

        open(SDOFmon,file=STRfX,status='replace')
        close(SDOFmon)

        file="STR000000.FY"

        if(len_trim(monitor_file) .ne. 70) then
          STRfY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          STRfY = file
        endif

        open(SDOFmon,file=STRfY,status='replace')
        close(SDOFmon)

        !!! foundation shear force

        file="FND000000.FX"

        if(len_trim(monitor_file) .ne. 70) then
          FNDfX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDfX = file
        endif

        open(SDOFmon,file=FNDfX,status='replace')
        close(SDOFmon)

        file="FND000000.FY"

        if(len_trim(monitor_file) .ne. 70) then
          FNDfY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          FNDfY = file
        endif

        open(SDOFmon,file=FNDfY,status='replace')
        close(SDOFmon)

        !!! interaction force

        file="INT000000.FX"

        if(len_trim(monitor_file) .ne. 70) then
          INTfX = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          INTfX = file
        endif

        open(SDOFmon,file=INTfX,status='replace')
        close(SDOFmon)

        file="INT000000.FY"

        if(len_trim(monitor_file) .ne. 70) then
          INTfY = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          INTfY = file
        endif

        open(SDOFmon,file=INTfY,status='replace')
        close(SDOFmon)

        file="INT000000.FZ"

        if(len_trim(monitor_file) .ne. 70) then
          INTfZ = monitor_file(1:len_trim(monitor_file)) // '/' // file
        else
          INTfZ = file
        endif

        open(SDOFmon,file=INTfZ,status='replace')
        close(SDOFmon)
      endif
    endif
  endif

  return
end subroutine MAKE_SDOF_OUTPUT_FILES
