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

!> @brief Read input parameters for oscillators
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0
!> @param[in] file_sdof   file name (SDOF00000X.info)
!> @param[in] id          mpi process id
!> @param[in] dtsite      soil time step

subroutine MAKE_SDOF_SYSTEM(file_sdof,dtsite,id)

  use SDOF_SYSTEM

	implicit none

  integer*4 :: id, SDOFunit, i
  character*70 :: file_sdof
  real*8 :: dtsite
  real*8, dimension(4,4) :: mat_temp

  SDOFunit=1+(id+1)*10

  open(SDOFunit,file=file_sdof)			!!! open SDOFINFO.txt
  read(SDOFunit,*)      !!! here is the number of SDOF oscillators

  do i=1, n_sdof

    read (SDOFunit,"(2I10)") sys(i)%ID, sys(i)%const_law

    sys(i)%SFS = 0

    if (sys(i)%const_law.eq.1) then
      read (SDOFunit,"(1I10)") sys(i)%SFS
      if (sys(i)%SFS.eq.0) then
        read(SDOFunit,"(4E15.7)") sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%TN
      elseif (sys(i)%SFS.eq.1) then
        read(SDOFunit,"(13E15.7)") sys(i)%Ms, sys(i)%Mf, sys(i)%J, sys(i)%height, &
                                  sys(i)%Ks, sys(i)%K0, sys(i)%Kr, sys(i)%Kv, sys(i)%CSI, &
                                  sys(i)%C0, sys(i)%Cr, sys(i)%Cv, sys(i)%TN

        read(SDOFunit,"(2E15.7)") sys(i)%beta_newmark, sys(i)%gamma_newmark

        sys(i)%MAT_M = 0; sys(i)%MAT_KS = 0; sys(i)%MAT_KI = 0; sys(i)%MAT_C = 0

        sys(i)%MAT_M(1,1) = sys(i)%Ms; sys(i)%MAT_M(2,2) = sys(i)%Mf; sys(i)%MAT_M(3,3) = sys(i)%J; sys(i)%MAT_M(4,4) = sys(i)%Ms + sys(i)%Mf

        sys(i)%MAT_KS(1,1) = sys(i)%Ks; sys(i)%MAT_KS(1,2) = -sys(i)%Ks; sys(i)%MAT_KS(1,3) = -sys(i)%Ks*sys(i)%height
        sys(i)%MAT_KS(2,1) = -sys(i)%Ks; sys(i)%MAT_KS(2,2) = sys(i)%Ks; sys(i)%MAT_KS(2,3) = sys(i)%Ks*sys(i)%height
        sys(i)%MAT_KS(3,1) = -sys(i)%Ks*sys(i)%height; sys(i)%MAT_KS(3,2) = sys(i)%Ks*sys(i)%height; sys(i)%MAT_KS(3,3) = sys(i)%Ks*(sys(i)%height**2)

        sys(i)%MAT_KI(2,2) = sys(i)%K0; sys(i)%MAT_KI(3,3) = sys(i)%Kr; sys(i)%MAT_KI(4,4) = sys(i)%Kv

        sys(i)%Cs = datan(1.d0)*16.d0*sys(i)%Ms*sys(i)%CSI/sys(i)%TN

        sys(i)%MAT_C(1,1) = sys(i)%Cs; sys(i)%MAT_C(1,2) = -sys(i)%Cs; sys(i)%MAT_C(1,3) = -sys(i)%Cs*sys(i)%height
        sys(i)%MAT_C(2,1) = -sys(i)%Cs; sys(i)%MAT_C(2,2) = sys(i)%Cs + sys(i)%C0; sys(i)%MAT_C(2,3) = sys(i)%Cs*sys(i)%height
        sys(i)%MAT_C(3,1) = -sys(i)%Cs*sys(i)%height; sys(i)%MAT_C(3,2) = sys(i)%Cs*sys(i)%height; sys(i)%MAT_C(3,3) = sys(i)%Cs*(sys(i)%height**2) + sys(i)%Cr
        sys(i)%MAT_C(4,4) = sys(i)%Cv

        sys(i)%u = 0; sys(i)%v = 0; sys(i)%a = 0; sys(i)%f = 0
      endif

    elseif (sys(i)%const_law.eq.2) then
      read(SDOFunit,"(5E15.7)") sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%FY, sys(i)%TN

    elseif (sys(i)%const_law.eq.3) then
      read(SDOFunit,"(9E15.7)") sys(i)%Ms, sys(i)%Ks, sys(i)%Hs, sys(i)%Ss, sys(i)%CSI, sys(i)%FY, sys(i)%EH, sys(i)%EU, sys(i)%TN
      sys(i)%EY = sys(i)%FY/sys(i)%Ks
      sys(i)%FH = sys(i)%FY + sys(i)%Hs*(sys(i)%EH - sys(i)%EY)
      sys(i)%FU = sys(i)%FH + sys(i)%Ss*(sys(i)%EU - sys(i)%EH)
      sys(i)%branch = 1
      sys(i)%damage = 0

    endif

    sys(i)%SDOFIDR=0
    sys(i)%dSDOFIDR=0
    sys(i)%SDOFItF=0
    sys(i)%tempSDOFU0=0
    sys(i)%tempSDOFU1=0
    sys(i)%tempSDOFU=0

    sys(i)%dt=sys(i)%TN/(10.d0)  !(datan(1.d0)*4.d0)    - Stable time step for explicit central difference method for SDOF
    sys(i)%ndt=ceiling(dtsite/sys(i)%dt)
    sys(i)%dt=dtsite/sys(i)%ndt
    sys(i)%dt2=sys(i)%dt*sys(i)%dt

    sys(i)%Cs = datan(1.d0)*16.d0*sys(i)%Ms*sys(i)%CSI/sys(i)%TN      ! Damping Coefficient

    if (sys(i)%SFS.eq.1) then
      sys(i)%MAT_F = 1/(sys(i)%beta_newmark*sys(i)%dt2)*sys(i)%MAT_M + sys(i)%gamma_newmark/(sys(i)%beta_newmark*sys(i)%dt)*sys(i)%MAT_C + sys(i)%MAT_KS
      mat_temp = sys(i)%MAT_F + sys(i)%MAT_KI
      !call MAT_INVERSE(mat_temp, sys(i)%MAT_MCinv, 4) ! Issue with MPIf90 compiler for Apple silicon 
    endif
  enddo

	close(SDOFunit)
  return
end subroutine MAKE_SDOF_SYSTEM
