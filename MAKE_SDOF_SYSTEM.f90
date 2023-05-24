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

subroutine MAKE_SDOF_SYSTEM(file_bldinfo,dtsite,id)

  use SPEED_SCI

  implicit none

  integer*4 :: id, SDOFunit, i, j, dummy
  character*70 :: file_bldinfo
  real*8 :: dtsite, tempint
  real*8, dimension(4,4) :: mat_temp

  SDOFunit=1+(id+1)*10

  open(SDOFunit,file=file_bldinfo)    !!! open SDOFINFO.txt
  read(SDOFunit,*) n_bld              !!! here is the number of SDOF oscillators

  do i=1, n_bld
    sys(i)%SFS = 0
    sys(i)%ForceApplicationtype = 1

    read (SDOFunit,*) sys(i)%ID, sys(i)%StructType, sys(i)%const_law

    if (sys(i)%StructType.eq.1) then

      sys(i)%NDOF = 1
      sys(i)%flag_Minv = 0;

      read (SDOFunit,*) sys(i)%SFS
      sys(i)%SFS = 0

      allocate(sys(i)%Ms(1,1), sys(i)%Ms_inv(1,1), sys(i)%Ks(1,1))
      allocate(sys(i)%IDR(1,3), sys(i)%variIDR(1,3))
      allocate(sys(i)%IntForce(1,3))
      allocate(sys(i)%tempU0(1,3), sys(i)%tempU1(1,3), sys(i)%tempU(1,3))
      allocate(sys(i)%tempA1(sys(i)%NDOF,3), sys(i)%tempRA1(sys(i)%NDOF,3))

      if (sys(i)%const_law.eq.1) then
        !read (SDOFunit,"(1I10)") sys(i)%SFS
        if (sys(i)%SFS.eq.0) then
          read(SDOFunit,*) sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%TN
          write(*,*) 'SDOF No: ', i, sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%TN
        elseif (sys(i)%SFS.eq.1) then
          read(SDOFunit,*) sys(i)%Ms, sys(i)%Mf, sys(i)%J, sys(i)%height, &
                                    sys(i)%Ks, sys(i)%K0, sys(i)%Kr, sys(i)%Kv, sys(i)%CSI, &
                                    sys(i)%C0, sys(i)%Cr, sys(i)%Cv, sys(i)%TN

          read(SDOFunit,*) sys(i)%beta_newmark, sys(i)%gamma_newmark

          sys(i)%MAT_M = 0; sys(i)%MAT_KS = 0; sys(i)%MAT_KI = 0; sys(i)%MAT_C = 0

          sys(i)%MAT_M(1,1) = sys(i)%Ms(1,1);     sys(i)%MAT_M(2,2) = sys(i)%Mf; 
          sys(i)%MAT_M(3,3) = sys(i)%J;           sys(i)%MAT_M(4,4) = sys(i)%Ms(1,1) + sys(i)%Mf

          sys(i)%MAT_KS(1,1) = sys(i)%Ks(1,1); 
          sys(i)%MAT_KS(1,2) = -sys(i)%Ks(1,1); 
          sys(i)%MAT_KS(1,3) = -sys(i)%Ks(1,1)*sys(i)%height;
          sys(i)%MAT_KS(2,1) = -sys(i)%Ks(1,1); 
          sys(i)%MAT_KS(2,2) = sys(i)%Ks(1,1); 
          sys(i)%MAT_KS(2,3) = sys(i)%Ks(1,1)*sys(i)%height;
          sys(i)%MAT_KS(3,1) = -sys(i)%Ks(1,1)*sys(i)%height; 
          sys(i)%MAT_KS(3,2) = sys(i)%Ks(1,1)*sys(i)%height; 
          sys(i)%MAT_KS(3,3) = sys(i)%Ks(1,1)*(sys(i)%height**2)

          sys(i)%MAT_KI(2,2) = sys(i)%K0; sys(i)%MAT_KI(3,3) = sys(i)%Kr; sys(i)%MAT_KI(4,4) = sys(i)%Kv

          sys(i)%Cs = datan(1.d0)*16.d0*sys(i)%Ms(1,1)*sys(i)%CSI/sys(i)%TN

          sys(i)%MAT_C(1,1) = sys(i)%Cs; sys(i)%MAT_C(1,2) = -sys(i)%Cs; sys(i)%MAT_C(1,3) = -sys(i)%Cs*sys(i)%height
          sys(i)%MAT_C(2,1) = -sys(i)%Cs; sys(i)%MAT_C(2,2) = sys(i)%Cs + sys(i)%C0; sys(i)%MAT_C(2,3) = sys(i)%Cs*sys(i)%height
          sys(i)%MAT_C(3,1) = -sys(i)%Cs*sys(i)%height; sys(i)%MAT_C(3,2) = sys(i)%Cs*sys(i)%height; 
          sys(i)%MAT_C(3,3) = sys(i)%Cs*(sys(i)%height**2) + sys(i)%Cr
          sys(i)%MAT_C(4,4) = sys(i)%Cv

          sys(i)%u = 0; sys(i)%v = 0; sys(i)%a = 0; sys(i)%f = 0
        endif

      elseif (sys(i)%const_law.eq.2) then
        read(SDOFunit,"(5E15.7)") sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%FY, sys(i)%TN

      elseif (sys(i)%const_law.eq.3) then
        ! Mass, Initial Stiffness, post-yield stiffness, Softening stiffness, dampingfactor, yield strength, displacement at peak strength, ultimate displacement, natural time period
        read(SDOFunit,"(9E15.7)") sys(i)%Ms, sys(i)%Ks, sys(i)%Hs, sys(i)%Ss, sys(i)%CSI, sys(i)%FY, sys(i)%EH, sys(i)%EU, sys(i)%TN
        sys(i)%EY = sys(i)%FY/sys(i)%Ks(1,1)          !Yield Displacement
        sys(i)%FH = sys(i)%FY + sys(i)%Hs*(sys(i)%EH - sys(i)%EY)         !Max. Peak strength
        sys(i)%FU = sys(i)%FH + sys(i)%Ss*(sys(i)%EU - sys(i)%EH)         !Ultimate failure strength after softening
        sys(i)%branch = 1
        sys(i)%damage = 0

        write(*,*) 'SDOF No: ', i, sys(i)%Ms, sys(i)%Ks, sys(i)%CSI, sys(i)%TN

      endif

      sys(i)%Cs = datan(1.d0)*16.d0*sys(i)%Ms(1,1)*sys(i)%CSI/sys(i)%TN      ! Damping Coefficient

    elseif (sys(i)%StructType.eq.2) then
      
      if(.not.isConfigPresent) then
        write(*,*) "*****For MDOF Structures, config.txt file is missing!*****"
        CALL EXIT()
      endif

      read(SDOFunit,*) sys(i)%NDOF, sys(i)%Floor_h, sys(i)%Area, sys(i)%ForceApplicationtype
      sys(i)%ForceApplicationtype = 1
      
      sys(i)%Height=sys(i)%NDOF*sys(i)%Floor_h        
            
      allocate(sys(i)%props(sys(i)%NDOF,10,1))
      allocate(sys(i)%dval(sys(i)%NDOF,4))
      allocate(sys(i)%Ms(sys(i)%NDOF,sys(i)%NDOF))
      allocate(sys(i)%Ks(sys(i)%NDOF,sys(i)%NDOF))
      allocate(sys(i)%SysC(sys(i)%NDOF,sys(i)%NDOF))
      allocate(sys(i)%Ms_Inv(sys(i)%NDOF,sys(i)%NDOF))
      allocate(sys(i)%MDOFEt(sys(i)%NDOF,2)) 
      allocate(sys(i)%MDOFstatev(sys(i)%NDOF,11,2))
      allocate(sys(i)%MDOFspd(sys(i)%NDOF,2))
      allocate(sys(i)%MDOFyield(sys(i)%NDOF,2))
      allocate(sys(i)%IDR(sys(i)%NDOF,3))
      allocate(sys(i)%variIDR(sys(i)%NDOF,3))
      allocate(sys(i)%IntForce(sys(i)%NDOF,3))
      allocate(sys(i)%tempU0(sys(i)%NDOF,3))
      allocate(sys(i)%tempU1(sys(i)%NDOF,3))
      allocate(sys(i)%tempU(sys(i)%NDOF,3))
      allocate(sys(i)%MaxIDR(sys(i)%NDOF,2))
      allocate(sys(i)%MDOFIDeath(sys(i)%NDOF,2))
      allocate(sys(i)%MDOFstate(sys(i)%NDOF,2))
      allocate(sys(i)%tempA1(sys(i)%NDOF,3))
      allocate(sys(i)%tempRA1(sys(i)%NDOF,3))
      sys(i)%Ms=0;       sys(i)%Ks=0;            sys(i)%SysC=0;      sys(i)%Ms_Inv=0; 
      sys(i)%MDOFEt=0;   sys(i)%MDOFstatev=0;    sys(i)%MDOFspd=0;   sys(i)%MDOFyield=0; 
      sys(i)%MaxIDR=0;   sys(i)%MDOFIDeath=0;    sys(i)%MDOFstate=0
      sys(i)%flag_Minv = 0;

      do j=1, sys(i)%NDOF
        read(SDOFunit,*)sys(i)%props(j,1:10,1)
      enddo

    endif

    sys(i)%IDR=0
    sys(i)%variIDR=0
    sys(i)%IntForce=0
    sys(i)%tempU0=0
    sys(i)%tempU1=0
    sys(i)%tempU=0

    if (sys(i)%SFS.eq.1) then
      sys(i)%MAT_F = 1/(sys(i)%beta_newmark*sys(i)%dt2) * sys(i)%MAT_M + &
                     sys(i)%gamma_newmark/(sys(i)%beta_newmark*sys(i)%dt)*sys(i)%MAT_C + sys(i)%MAT_KS
      mat_temp = sys(i)%MAT_F + sys(i)%MAT_KI
      call matinv(mat_temp, sys(i)%MAT_MCinv, 4)
    endif
  enddo

  MaxDOF_loc = 0
  do i=1, n_bld
    if (MaxDOF_loc.lt.sys(i)%NDOF) then
      MaxDOF_loc = sys(i)%NDOF
    endif

    if (sys(i)%StructType.eq.2) then
      read(SDOFunit,*) 
      do j=1, sys(i)%NDOF
        read(SDOFunit,*)sys(i)%dval(j,1:4)
      enddo 
    endif
  enddo    
  
  if (sys(1)%StructType.eq.2) read(SDOFunit,*) 
  do i=1, n_bld
    if (sys(i)%StructType.eq.2) then
      read(SDOFunit,*)dummy, sys(i)%T1, sys(i)%T2, sys(i)%Calpha, sys(i)%Cbeta, sys(i)%TN
    endif
  enddo

  close(SDOFunit)


  ! Making Mass, Damping Matrices
  do i=1, n_bld
    sys(i)%dt=sys(i)%TN/(datan(1.d0)*4.d0)    !- Stable time step for explicit central difference method for SDOF
    sys(i)%ndt=ceiling(dtsite/sys(i)%dt)
    sys(i)%dt=dtsite/sys(i)%ndt
    sys(i)%dt2=sys(i)%dt*sys(i)%dt

    !!!!!!!!!!!MassMatrix!!!!!!!!!!!!
    if (sys(i)%StructType.eq.2) then
      do j=1,sys(i)%NDOF
          sys(i)%Ms(j,j)=MasspArea*sys(i)%Area
      enddo
    endif

    !!!!!!!!!!!StiffMatrix!!!!!!!!!!!!
    if (sys(i)%StructType.eq.2) then
      do j=1,sys(i)%NDOF
          sys(i)%Ks(j,j)= sys(i)%Ks(j,j) + sys(i)%props(j,1,1) 
          if(j>1) then
              sys(i)%Ks(j-1,j-1)= sys(i)%Ks(j-1,j-1) + sys(i)%props(j,1,1)
              sys(i)%Ks(j,j-1)  = sys(i)%Ks(j,j-1)   - sys(i)%props(j,1,1)
              sys(i)%Ks(j-1,j)  = sys(i)%Ks(j-1,j)   - sys(i)%props(j,1,1)
          endif
      enddo
    endif

    !!!!!!!!!!! DampingMatrix - Rayleigh !!!!!!!!!!!!
    if (sys(i)%StructType.eq.2) then
      sys(i)%SysC(1:sys(i)%NDOF,1:sys(i)%NDOF)= sys(i)%Calpha * sys(i)%Ms(1:sys(i)%NDOF,1:sys(i)%NDOF) &
      + sys(i)%Cbeta * sys(i)%Ks(1:sys(i)%NDOF,1:sys(i)%NDOF)
      deallocate(sys(i)%Ks)
    endif
  enddo

  return
end subroutine MAKE_SDOF_SYSTEM
