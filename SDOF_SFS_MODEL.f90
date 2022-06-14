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

!> @brief Response computation for 4DOF model.
!> @author Aline Herlin
!> @date February, 2021
!> @version 1.0
!> @param[in] sID           oscillator ID
!> @param[in] gr_acc        (negative) ground acceleration
!> @param[in] direction     motion direction

subroutine SDOF_SFS_MODEL(sID, gr_acc, direction)

  use SPEED_SCI

  implicit none

  integer*4 :: sID, direction
  real*8, dimension(3) :: gr_acc
  real*8, dimension(4) :: TN, TN1, TN2, TN3, TN4, temp, du, df, vec, a_old, v_old, u_old, f_old

  TN = 0; TN1 = 0; TN2 = 0; TN3 = 0; TN4 = 0; temp = 0; du = 0; df = 0; vec = 0
  a_old = sys(sID)%a(:,direction); v_old = sys(sID)%v(:,direction); u_old = sys(sID)%u(:,direction); f_old = sys(sID)%f(:,direction)

  TN1(1) = sys(sID)%Ms(1,1)*gr_acc(direction)
  TN1(2) = sys(sID)%Mf*gr_acc(direction)
  TN1(4) = (sys(sID)%Ms(1,1) + sys(sID)%Mf)*gr_acc(3)

  temp = (1 - 2*sys(sID)%beta_newmark)/(2*sys(sID)%beta_newmark)*a_old + 1/(sys(sID)%beta_newmark*sys(sID)%dt)*v_old
  TN2 = MATMUL(sys(sID)%MAT_M,temp)

  temp = (sys(sID)%gamma_newmark - 2*sys(sID)%beta_newmark)/(2*sys(sID)%beta_newmark)*sys(sID)%dt*a_old + &
          (sys(sID)%gamma_newmark - sys(sID)%beta_newmark)/(sys(sID)%beta_newmark)*v_old
  TN3 = MATMUL(sys(sID)%MAT_C,temp)

  TN4 = MATMUL(sys(sID)%MAT_KS,u_old)

  TN = TN1 + TN2 + TN3 - TN4 - f_old

  du = MATMUL(sys(sID)%MAT_MCinv,TN)

  df = MATMUL(sys(sID)%MAT_KI,du)

  sys(sID)%u(:,direction) = u_old + du

  sys(sID)%f(:,direction) = f_old + df

  vec = MATMUL(sys(sID)%MAT_KS,sys(sID)%u(:,direction))
  sys(sID)%fs(direction) = vec(1)
  sys(sID)%fb(direction) = sys(sID)%f(2,direction) + sys(sID)%Mf*(gr_acc(direction) - a_old(2))

  sys(sID)%v(:,direction) = v_old + (2*sys(sID)%beta_newmark - sys(sID)%gamma_newmark)/(2*sys(sID)%beta_newmark)*sys(sID)%dt*a_old + &
    sys(sID)%gamma_newmark/(sys(sID)%beta_newmark*sys(sID)%dt)*(du - v_old*sys(sID)%dt)

  sys(sID)%a(:,direction) = (2*sys(sID)%beta_newmark - 1)/(2*sys(sID)%beta_newmark)*a_old + &
    1/(sys(sID)%beta_newmark*sys(sID)%dt2)*(du - v_old*sys(sID)%dt)

  sys(sID)%IntForce(1,direction) = sys(sID)%f(2,direction)

  ! UnComment this after verifying
  !sys(sID)%SDOFItF(3) = sys(sID)%f(4,direction)

  return
end subroutine SDOF_SFS_MODEL
