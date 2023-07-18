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

!> @brief Computation for SDOF shear model.
!> @author Aline Herlin
!> @date November, 2020
!> @version 1.0
!> @param[in] sID           oscillator ID
!> @param[in] gr_acc        base acceleration
!> @param[in] direction     direction of motion

subroutine SDOF_SHEAR_MODEL (sID, ndof, gr_acc, direction)

  use SPEED_SCI

  implicit none

  integer*4 :: sID, ndof, direction, Usign
  integer*4 :: idof, i, j
  real*8 :: gr_acc
  real*8, dimension(ndof) :: e, de, s, sysPeff0, sysPeff

  sysPeff0 = 0

  do idof=1,ndof
    sysPeff(idof) = gr_acc*sys(sID)%Ms(idof,idof)     !!!mass x (-u_g'') = inertial force from the ground acceleration
  enddo

  if (ndof.gt.1) then
    do idof=1,ndof-1
      sysPeff0(idof) = sysPeff(idof) - sys(sID)%IntForce(idof,direction) + sys(sID)%IntForce(idof+1,direction)     !!! SDOFItF = interaction force Ku(n)
    enddo
  endif
  sysPeff0(ndof) = sysPeff(ndof) - sys(sID)%IntForce(ndof,direction)

  ! returns sys(sID)%tempU = displacement at time n+1 = Relative displacement of structure w.r.t. ground
  if (sys(sID)%StructType.eq.1) then
    call CENTRAL_DIFFERENCE(ndof, sys(sID)%Ms(1,1), sys(sID)%Ms_inv(1,1), sys(sID)%Cs, sys(sID)%dt, &
                            sys(sID)%tempU(1,direction), sys(sID)%tempU1(1,direction), &
                            sys(sID)%tempU0(1,direction), sysPeff0, sys(sID)%flag_Minv )
  else
    call CENTRAL_DIFFERENCE(ndof, sys(sID)%Ms(1:ndof,1:ndof), sys(sID)%Ms_inv(1:ndof,1:ndof), sys(sID)%SysC(1:ndof,1:ndof), &
                            sys(sID)%dt, sys(sID)%tempU(1:ndof,direction), sys(sID)%tempU1(1:ndof,direction), &
                            sys(sID)%tempU0(1:ndof,direction), sysPeff0, sys(sID)%flag_Minv )
    ! if (direction.eq.1) then
    !   write(*,*) 'InputForces: Acc =', gr_acc, 'sysPeff = ', (sysPeff0(idof), idof = 1,sys(1)%NDOF)
    !   write(*,*) 'Input Un-1 = ', (sys(1)%tempU0(idof,1), idof=1,sys(1)%NDOF), 'Un = ', (sys(1)%tempU1(idof,1), idof=1,sys(1)%NDOF)
    !   write(*,*) 'Output Un+1 = ', (sys(1)%tempU(idof,1), idof=1,sys(1)%NDOF)
    ! endif
  endif

  ! Modify this for SSI - with MDOF (below code works for SSI-SDOF)
  if(sys(sID)%const_law.eq.3) then
    if (sys(sID)%tempU(1, direction).ge.0.0) Usign = 1
    if (sys(sID)%tempU(1, direction).lt.0.0) Usign = -1
    if ((abs(sys(sID)%tempU(1, direction)).ge.sys(sID)%EU).or.(sys(sID)%damage(direction).eq.1)) &
        sys(sID)%tempU(1,direction) = Usign*sys(sID)%EU
  endif

  ! Time Variation of Interstory Drift (IDR)
  sys(sID)%variIDR(1,direction) = sys(sID)%tempU(1,direction) - sys(sID)%IDR(1,direction)
  sys(sID)%IDR(1,direction)=sys(sID)%tempU(1,direction); 
  if (sys(sID)%StructType.eq.2) then
    do j=2,ndof
        sys(sID)%variIDR(j,direction) = (sys(sID)%tempU(j,direction) - sys(sID)%tempU(j-1,direction)) - sys(sID)%IDR(j,direction)
        sys(sID)%IDR(J,direction) = sys(sID)%tempU(j,direction) - sys(sID)%tempU(j-1,direction)
    enddo
    !Max IDR 
    do j=1,ndof
        sys(sID)%MaxIDR(j,direction) = max(sys(sID)%MaxIDR(j,direction), abs(sys(sID)%IDR(j,direction)/sys(sID)%Floor_h))
    enddo
  endif

  do idof= 1,ndof
      !!! Relative acceleration
      sys(sID)%tempRA1(idof,direction) = ( sys(sID)%tempU(idof,direction) + sys(sID)%tempU0(idof,direction) - &
                                          2.d0*sys(sID)%tempU1(idof, direction))/sys(sID)%dt2
      !!! Absolute acceleration AH
      sys(sID)%tempA1(idof,direction) = sys(sID)%tempRA1(idof,direction) - gr_acc   

      ! Updating displacement Values
      sys(sID)%tempU0(idof,direction)=sys(sID)%tempU1(idof,direction)
      sys(sID)%tempU1(idof,direction)=sys(sID)%tempU(idof,direction)
      sys(sID)%tempU(idof,direction)=0
  enddo

  do idof=1,ndof
      e(idof) = sys(sID)%IDR(idof,direction) - sys(sID)%variIDR(idof,direction)       ! IDR at (n)th time step
      de(idof)= sys(sID)%variIDR(idof,direction)    ! change in IDR at (n+1)th timestep
      s(idof) = sys(sID)%IntForce(idof,direction)   ! Shear Force at (n)th timestep

      if (sys(sID)%StructType.eq.1) then
          if (sys(sID)%const_law.eq.1) then
            call LINEAR_ELASTIC(sys(sID)%Ks, s, de)      !!! elastic constitutive law
          elseif (sys(sID)%const_law.eq.2) then
            call PERFECTLY_PLASTIC(sys(sID)%Ks, s, de, sys(sID)%FY)      !!! elastoplastic constitutive law
          elseif (sys(sID)%const_law.eq.3) then
            call TRILINEAR(sys(sID)%Ks, sys(sID)%Hs, sys(sID)%Ss, s, de, e, sys(sID)%FY, sys(sID)%FH, sys(sID)%FU, sys(sID)%EY, &
                sys(sID)%EH, sys(sID)%EU, sys(sID)%branch(direction), sys(sID)%damage(direction))      !!! trilinear constitutive law
          endif
      elseif (sys(sID)%StructType.eq.2) then

        if (sys(sID)%const_law.eq.1) then
          call LINEAR_ELASTIC(sys(sID)%props(idof,1,1), s(idof), de(idof))
        elseif (sys(sID)%const_law.eq.2) then
          if ((abs(sys(sID)%IDR(idof,direction))/sys(sID)%Floor_h).GT.sys(sID)%dval(idof,4)) then
              sys(sID)%MDOFIDeath(idof,direction) = 1
          endif
          call ksteel02(sys(sID)%props(idof,1:10,1), s(idof), e(idof), de(idof), sys(sID)%MDOFEt(idof,direction), &
                        sys(sID)%MDOFstatev(idof,1:11,direction), sys(sID)%MDOFspd(idof,direction), &
                        sys(sID)%MDOFyield(idof,direction), sys(sID)%MDOFIDeath(idof,direction) ) !sys(sID)%Ms(1:ndof,1:ndof), ndof
        endif
      endif

      sys(sID)%IntForce(idof,direction) = s(idof)

  enddo

  return
end subroutine SDOF_SHEAR_MODEL
