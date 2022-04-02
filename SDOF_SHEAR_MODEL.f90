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

subroutine SDOF_SHEAR_MODEL (sID, gr_acc, direction)

  use SDOF_SYSTEM

  implicit none

  integer*4 :: sID, direction, Usign
  real*8 :: e, de, s, SDOFPeff0, gr_acc, SDOFPeff

  SDOFPeff = gr_acc*sys(sID)%Ms     !!!mass x (-u_g'') = inertial force from the ground acceleration

  SDOFPeff0=SDOFPeff-sys(sID)%SDOFItF(direction)      !!! SDOFItF = interaction force Ku(n)

  call CENTRAL_DIFFERENCE(sys(sID)%Ms, sys(sID)%Cs, sys(sID)%dt, &
                          sys(sID)%tempSDOFU(direction), sys(sID)%tempSDOFU1(direction), &
                          sys(sID)%tempSDOFU0(direction), SDOFPeff0)

                          !!! returns sys(sID)%tempSDOFU = displacement at time n+1 = Relative displacement of structure w.r.t. ground

  if(sys(sID)%const_law.eq.3) then
    if (sys(sID)%tempSDOFU(direction).ge.0.0) Usign = 1
    if (sys(sID)%tempSDOFU(direction).lt.0.0) Usign = -1
    if ((abs(sys(sID)%tempSDOFU(direction)).ge.sys(sID)%EU).or.(sys(sID)%damage(direction).eq.1)) sys(sID)%tempSDOFU(direction) = Usign*sys(sID)%EU
  endif

  sys(sID)%dSDOFIDR(direction)=sys(sID)%tempSDOFU(direction)-sys(sID)%SDOFIDR(direction)      !!! variation of base drift

  sys(sID)%SDOFIDR(direction)=sys(sID)%tempSDOFU(direction)     !!! base drift

  sys(sID)%tempSDOFA1(direction)=(sys(sID)%tempSDOFU(direction) + &
     sys(sID)%tempSDOFU0(direction) - 2.d0*sys(sID)%tempSDOFU1(direction))/sys(sID)%dt2-gr_acc      !!! absolute acceleration AH

  sys(sID)%tempSDOFU0(direction)=sys(sID)%tempSDOFU1(direction)      !!! update values
  sys(sID)%tempSDOFU1(direction)=sys(sID)%tempSDOFU(direction)
  sys(sID)%tempSDOFU(direction)=0

  e=sys(sID)%SDOFIDR(direction)
  de=sys(sID)%dSDOFIDR(direction)
  s=sys(sID)%SDOFItF(direction)

  if (sys(sID)%const_law.eq.1) then
    call LINEAR_ELASTIC(sys(sID)%Ks, s, de)      !!! elastic constitutive law
  elseif (sys(sID)%const_law.eq.2) then
    call PERFECTLY_PLASTIC(sys(sID)%Ks, s, de, sys(sID)%FY)      !!! elastoplastic constitutive law
  elseif (sys(sID)%const_law.eq.3) then
    call TRILINEAR(sys(sID)%Ks, sys(sID)%Hs, sys(sID)%Ss, s, de, e, sys(sID)%FY, sys(sID)%FH, sys(sID)%FU, sys(sID)%EY, sys(sID)%EH, sys(sID)%EU, sys(sID)%branch(direction), sys(sID)%damage(direction))      !!! trilinear constitutive law
  endif

  sys(sID)%SDOFItF(direction)=s

  return
end subroutine SDOF_SHEAR_MODEL
