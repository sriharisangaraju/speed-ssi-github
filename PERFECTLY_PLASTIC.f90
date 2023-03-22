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

!---------------------Elastic - Perfecty Plastic---------------------------
!> @brief Defines stress-strain relation for EPP SDOF oscillator
!> @author Aline Herlin
!> @date January, 2021
!> @version 1.0
!> @param[out] s      force developed in the SDOF system
!> @param[in] de      variation of the displacement wrt previous time step
!> @param[in] E0		  stiffness (force-displ ratio)
!> @param[in] fy			yield strength

subroutine PERFECTLY_PLASTIC(E0, s, de, fy)

	implicit none

	real*8 :: E0, de, s, s0, fy, temp

	s0 = s
  if (s0.ge.-fy .and. s0.le.fy) then
    temp = s0 + E0 * de
    if ((temp.ge.-fy).and.(temp.le.fy)) then
      s = temp
    elseif (temp.lt.-fy) then
      s = -fy
    elseif (temp.gt.fy) then
      s = fy
    endif
  endif

	return
end subroutine PERFECTLY_PLASTIC
