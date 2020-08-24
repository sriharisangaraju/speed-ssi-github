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


!> @brief ...Not-Honoring Enhanced (NHE) Implementation
!! @author Srihari Sangaraju
!> @date August, 2020 
!> @version 1.0

    subroutine GET_MECH_PROP_NH_ENHANCED(ie, nn, nn_loc, cs_nnz_loc, cs_loc, &
    									rho_nhe, lambda_nhe, mu_nhe,  &
                                        rho_el, lambda_el, mu_el)!gamma_el

    implicit none

    integer*4 :: ie, cs_nnz_loc, nn, nn_loc
    integer*4 :: r, q, p, is, ic
    integer*4, dimension(cs_nnz_loc) :: cs_loc
    
    real*8, dimension(nn,nn,nn) :: rho_el, lambda_el, mu_el	!gamma_el
    real*8, dimension(nn_loc) :: rho_nhe, lambda_nhe, mu_nhe


    do r = 1,nn
	    do q = 1,nn
	       do p = 1,nn
	          is = nn*nn*(r -1) +nn*(q -1) +p
	          ic = cs_loc(cs_loc(ie -1) +is)

	          rho_el(p,q,r) = rho_nhe(ic)
	          lambda_el(p,q,r) = lambda_nhe(ic)
	          mu_el(p,q,r) = mu_nhe(ic)
	          !gamma_el(p,q,r) = gamma
	       enddo
	    enddo
	enddo


	end subroutine GET_MECH_PROP_NH_ENHANCED
