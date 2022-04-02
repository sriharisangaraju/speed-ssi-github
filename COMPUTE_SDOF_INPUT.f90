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

!> @brief Computes acceleration of sdof systems
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0
!> @param[in] sdof_num                    number of oscillators
!> @param[in] mpi_id                      process id
!> @param[in] elem_mlst                   list of monitored elements
!> @param[in] local_el_num                local element numbering
!> @param[in] ne_loc                      number of local elements
!> @param[in] cs_loc                      connectivity vector
!> @param[in] cs_nnz_loc                  length of cs_loc
!> @param[in] sdeg_mat                    spectral degree of materials
!> @param[in] nmat                        number of materials
!> @param[in] u2                          soil displacement at time n+1
!> @param[in] nnod_loc                    number of local nodes
!> @param[in] xr_mlst, yr_mlst, zr_mlst   oscillators coordinates
!> @param[in] dt2                         time step square
!> @param[in] mpi_np                      total number of mpi processes
!> @param[out] ub1,ub2,ub3                base displacement at time n+1,n,n-1
!> @param[out] SDOFinputab                base acceleration
!> @param[out] SDOFinputdispl             base displacement

subroutine COMPUTE_SDOF_INPUT(sdof_num, mpi_id, elem_mlst, local_el_num, ne_loc, &
                              cs_loc, cs_nnz_loc, sdeg_mat, nmat, &
                              u2, nnod_loc, &
                              xr_mlst, yr_mlst, zr_mlst, &
                              dt2, &
                              SDOFinputab, SDOFinputdispl, mpi_np, &
                              ub1, ub2, ub3)

  implicit none

  integer*4 :: imon, ielem, ie, im, nn, k, j, i, is, in, iaz, sdof_num, ishift
  integer*4 :: mpi_id, ne_loc, nmat, nnod_loc, cs_nnz_loc, mpi_np
  integer*4 :: nmpi
  integer*4, dimension(0:cs_nnz_loc) :: cs_loc
  integer*4, dimension(sdof_num) :: elem_mlst
  integer*4, dimension(nmat) :: sdeg_mat
  integer*4, dimension(ne_loc) :: local_el_num
  real*8 :: dt2
  real*8 :: uxm,uym,uzm
  real*8 :: axm,aym,azm
  real*8, dimension(:), allocatable :: ct, ww
  real*8, dimension(:,:), allocatable :: dd
  real*8, dimension(:,:,:), allocatable :: ux_el, uy_el, uz_el
  real*8, dimension(3*nnod_loc) :: u2
  real*8, dimension(sdof_num) ::xr_mlst, yr_mlst, zr_mlst
  real*8, dimension(3*sdof_num) :: SDOFinputab, SDOFinputdispl, ub1, ub2, ub3

  ub3(1:3*sdof_num)=ub2(1:3*sdof_num)
  ub2(1:3*sdof_num)=ub1(1:3*sdof_num)

  ishift = 0

  do imon = 1,sdof_num     !!! loop over the oscillators

    ielem = elem_mlst(imon)
    call GET_INDLOC_FROM_INDGLO(local_el_num, ne_loc, ielem, ie)
    if (ie .ne. 0) then
      im = cs_loc(cs_loc(ie -1) +0)     !!! material
      nn = sdeg_mat(im) +1              !!! n of quadrature nodes

      allocate(ct(nn),ww(nn),dd(nn,nn))
      allocate(ux_el(nn,nn,nn),uy_el(nn,nn,nn),uz_el(nn,nn,nn))

      call MAKE_LGL_NW(nn,ct,ww,dd)

      do k = 1,nn
        do j = 1,nn
          do i = 1,nn
            is = nn*nn*(k -1) +nn*(j -1) +i
            in = cs_loc(cs_loc(ie -1) + is)

            iaz = 3*(in -1) +1; ux_el(i,j,k) = u2(iaz)
            iaz = 3*(in -1) +2; uy_el(i,j,k) = u2(iaz)
            iaz = 3*(in -1) +3; uz_el(i,j,k) = u2(iaz)

          enddo
        enddo
      enddo

      !-------------------------------------------------------------
      ! ACCELERATION AND SOIL FORCES
      !-------------------------------------------------------------
      call GET_MONITOR_VALUE(nn,ct,ux_el,&
                             xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uxm)
      call GET_MONITOR_VALUE(nn,ct,uy_el,&
                             xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uym)
      call GET_MONITOR_VALUE(nn,ct,uz_el,&
                             xr_mlst(imon),yr_mlst(imon),zr_mlst(imon),uzm)

      if (dabs(uxm).lt.1.0e-30) uxm = 0.0e+00
      if (dabs(uym).lt.1.0e-30) uym = 0.0e+00
      if (dabs(uzm).lt.1.0e-30) uzm = 0.0e+00

      ub1(ishift+1) = uxm  !!! x displacement at the base of the structure at time n+1
      ub1(ishift+2) = uym  !!! y displacement at the base of the structure at time n+1
      ub1(ishift+3) = uzm  !!! z displacement at the base of the structure at time n+1

      ! Central Difference Scheme
      axm = (ub1(ishift+1) -2.0d0*ub2(ishift+1) +ub3(ishift+1)) / dt2
      aym = (ub1(ishift+2) -2.0d0*ub2(ishift+2) +ub3(ishift+2)) / dt2
      azm = (ub1(ishift+3) -2.0d0*ub2(ishift+3) +ub3(ishift+3)) / dt2

      if (dabs(axm).lt.1.0e-30) axm = 0.0e+00
      if (dabs(aym).lt.1.0e-30) aym = 0.0e+00
      if (dabs(azm).lt.1.0e-30) azm = 0.0e+00

      do nmpi=1,mpi_np
        if (mpi_id .eq. (nmpi-1)) then
          SDOFinputab(3*imon-2)=axm
          SDOFinputab(3*imon-1)=aym
          SDOFinputab(3*imon)=azm
          SDOFinputdispl(3*imon-2)=uxm
          SDOFinputdispl(3*imon-1)=uym
          SDOFinputdispl(3*imon)=uzm
        endif
      enddo

      deallocate(ct,ww,dd)
      deallocate(ux_el,uy_el,uz_el)

      ishift = ishift + 3
    endif !(ie .ne. 0)
  enddo

  return

end subroutine COMPUTE_SDOF_INPUT
