!    Copyright (c) 2012 The SPEED FOUNDATION
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

!> @brief Computes the displacement of the oscillator through central difference scheme
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0
!> @param[in] m           mass of the system
!> @param[in] c           damping of system
!> @param[in] dT          time step for system
!> @param[in] U1, U0      displacement at time instants n and n-1
!> @param[in] P           equivalent force [-Mu''-Fint]
!> @param[out] U          displacement at time instant n+1

subroutine CENTRAL_DIFFERENCE(m, c, dT, U, U1, U0, P)

	implicit none

	real*8 :: m, c, m1, m2, m3, m4
	real*8 :: U, U1, U0, P, dT, x, y, T

	m1=m/dT/dT + c/2./dT;			!!! coeff for u(n+1)
	m2=-2.*m/dT/dT;					!!! coeff for u(n)
	m3=m/dT/dT-c/2./dT;			!!! coeff for u(n-1)

	T=P-m2*U1-m3*U0

	U=T/m1

	return
end subroutine CENTRAL_DIFFERENCE
