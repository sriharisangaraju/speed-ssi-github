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

!> @brief Defines stress-strain relation for TRILINEAR SDOF oscillator
!> @author Aline Herlin
!> @date February, 2021
!> @version 1.0
!> @param[out] f      force
!> @param[in] e       drift
!> @param[in] de      variation of the displacement
!> @param[in] K, H, S, fy, fh, fu, ey, eh, eu, branch, damage 			parameters for constitutive law

subroutine TRILINEAR(K, H, S, f, de, e, fy, fh, fu, ey, eh, eu, branch, damage)

	implicit none

	real*8 :: K, H, S
	real*8 :: de, f, e0, temp, e, e_int, f0
	real*8 :: fy, fh, fu, ey, eh, eu, eh_base, finf, fsup
	integer*4 :: branch, damage, sgn


	f0 = f
	! e0 = e - de;
	if (e.ge.0.0) sgn = 1
	if (e.lt.0.0) sgn = -1

	if ((damage.eq.0).and.(abs(e).lt.eu)) then

		if (branch.eq.1) then			! LINEAR ELASTIC branch
			if (abs(e).lt.ey) then
				f = K*e
				branch = 1
			elseif ((abs(e).ge.ey).and.(abs(e).lt.eh)) then
				f = sgn*(fy + H*(abs(e)-ey))
				if (sgn.eq.1) branch = 2
				if (sgn.eq.-1) branch = 4
			elseif (abs(e).ge.eh) then
				f = sgn*(fh + S*(abs(e)-eh))
				if (sgn.eq.1) branch = 3
				if (sgn.eq.-1) branch = 5
			endif

		elseif (branch.eq.2) then			! HARDENING branch (positive)
			if (de.ge.0.0) then
				if (e.lt.eh) then
					f = fy + H*(e-ey)
					branch = 2
				elseif (e.ge.eh) then
					f = fh + S*(e-eh)
					branch = 3
				endif
			elseif (de.lt.0.0) then
				e_int = e0 - f0/K
				temp = f0 + K*de
				finf = K/(K-H) * (-fy + H*(e_int + ey))
				if (temp.gt.finf) then
					f = temp
					branch = 6
				elseif (temp.le.finf) then
					if (e.gt.-eh) then
						f = -fy + H*(e+ey)
						branch = 4
					elseif (e.le.-eh) then
						f = -fh + S*(e+eh)
						branch = 5
					endif
				endif
			endif

		elseif (branch.eq.3) then			! SOFTENING branch (positive)
			if (de.ge.0.0) then
				f = fh + S*(e-eh)
				branch = 3
			elseif (de.lt.0.0) then
				e_int = e0 - f0/K
				temp = f0 + K*de
				finf = K/(K-H) * (-fy + H*(e_int + ey))
				if (temp.gt.finf) then
					f = temp
					branch = 6
				elseif (temp.le.finf) then
					if (e.gt.-eh) then
						f= -fy + H*(e+ey)
						branch = 4
					elseif (e.le.-eh) then
						f = -fh + S*(e+eh)
						branch = 5
					endif
				endif
			endif

		elseif (branch.eq.4) then			! HARDENING branch (negative)
			if (de.le.0.0) then
				if (e.gt.-eh) then
					f = -fy + H*(e+ey)
					branch = 4
				elseif (e.le.-eh) then
					f = -fh + S*(e+eh)
					branch = 5
				endif
			elseif (de.gt.0.0) then
				e_int = e0 - f0/K
				temp = f0 + K*de
				fsup = K/(K-H) * (fy + H*(e_int - ey))
				if (temp.lt.fsup) then
					f = temp
					branch = 6
				elseif (temp.ge.fsup) then
					if (e.lt.eh) then
						f = fy + H*(e-ey)
						branch = 2
					elseif (e.ge.eh) then
						f = fh + S*(e-eh)
						branch = 3
					endif
				endif
			endif

		elseif (branch.eq.5) then			! SOFTENING branch (negative)
			if (de.le.0.0) then
				f = -fh + S*(e+eh)
				branch = 5
			elseif (de.gt.0.0) then
				e_int = e0 - f0/K
				temp = f0 + K*de
				fsup = K/(K-H) * (fy + H*(e_int - ey))
				if (temp.lt.fsup) then
					f = temp
					branch = 6
				elseif (temp.ge.fsup) then
					if (e.lt.eh) then
						f= fy + H*(e-ey)
						branch = 2
					elseif (e.ge.eh) then
						f = fh + S*(e-eh)
						branch = 3
					endif
				endif
			endif

		elseif (branch.eq.6) then			! transition branch
			e_int = e0 - f0/K
			temp = f0 + K*de
			eh_base = eh - fh/K
			if (abs(e_int).le.eh_base) then
				fsup = K/(K-H) * (fy + H*(e_int - ey))
				finf = K/(K-H) * (-fy + H*(e_int + ey))
			elseif (e_int.gt.eh_base) then
				fsup = K/(K-S) * (fh + S*(e_int - eh))
				finf = K/(K-H) * (-fy + H*(e_int + ey))
			elseif (e_int.lt.-eh_base) then
				finf = K/(K-S) * (-fh + H*(e_int + eh))
				fsup = K/(K-H) * (fy + H*(e_int - ey))
			endif
			if ((temp.gt.finf).and.(temp.lt.fsup)) then
				f = temp
				branch = 6
			elseif (temp.ge.fsup) then
				if (e.lt.eh) then
					f = fy + H*(e-ey)
					branch = 2
				elseif (e.ge.eh) then
					f = fh + S*(e-eh)
					branch = 3
				endif
			elseif (temp.le.finf) then
				if (e.gt.-eh) then
					f = -fy + H*(e+ey)
					branch = 4
				elseif (e.le.-eh) then
					f = -fh + S*(e+eh)
					branch = 5
				endif
			endif
		endif

	elseif ((damage.eq.1).or.(abs(e).ge.eu)) then
		damage = 1
		f = 0
	endif

return 
end subroutine TRILINEAR
